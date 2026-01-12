#include <m_pd.h>
#include <string.h>
#include <math.h>

#define MINIMP3_IMPLEMENTATION
#include <minimp3_ex.h>

#define CLOWNRESAMPLER_IMPLEMENTATION
#include <clownresampler.h>

static t_class *mar_tilde_class;

typedef struct _mar_tilde {
    t_object x_obj;
    t_sample x_f;
    t_clock *clock;

    // MP3 state
    mp3dec_t mp3d;
    mp3dec_file_info_t info;
    int file_loaded;

    // Playback state
    size_t current_frame;
    size_t total_frames;
    int playing;
    int loop;
    int channels;
    int using_resampled;

    // resample
    float *resampled_buffer;
    size_t resampled_frames;
    int source_samplerate;
    // clownresampler
    ClownResampler_Precomputed cr_pre;
    ClownResampler_LowLevel_State cr_state;
    cc_s16l *cr_input;
    size_t cr_input_frames;

    // DSP
    int block_size;
    t_canvas *canvas;
    t_outlet *bang_out;
} t_mar_tilde;

// ─────────────────────────────────────
typedef struct {
    float *out;
    size_t frames_remaining;
    int channels;
} cr_out_cb_data;

// ─────────────────────────────────────
static cc_bool mar_cr_output_cb(void *user, const cc_s32f *frame, cc_u8f ch) {
    cr_out_cb_data *d = (cr_out_cb_data *)user;

    if (d->frames_remaining == 0) {
        return cc_false;
    }

    for (cc_u8f i = 0; i < ch; ++i) {
        *d->out++ = frame[i] * (1.0f / 32768.0f);
    }

    d->frames_remaining--;
    return cc_true;
}

// ─────────────────────────────────────
static int mar_tilde_resample_audio(t_mar_tilde *x) {
    int target_sr = sys_getsr();
    if (target_sr <= 0 || x->source_samplerate <= 0) {
        x->using_resampled = 0;
        return 0;
    }

    if (target_sr == x->source_samplerate) {
        x->using_resampled = 0;
        return 1;
    }

    size_t in_frames = x->info.samples / x->channels;

    /* Precompute kernel once */
    ClownResampler_Precompute(&x->cr_pre);

    if (!ClownResampler_LowLevel_Init(&x->cr_state, (cc_u8f)x->channels,
                                      (cc_u32f)x->source_samplerate, (cc_u32f)target_sr, 20000)) {
        pd_error(x, "[mar~] clownresampler init failed");
        return 0;
    }

    size_t pad = x->cr_state.lowest_level.integer_stretched_kernel_radius;
    x->cr_input_frames = in_frames;

    /* Allocate padded input buffer (int16) */
    x->cr_input = (cc_s16l *)malloc((in_frames + pad * 2) * x->channels * sizeof(cc_s16l));

    if (!x->cr_input) {
        return 0;
    }

    /* Zero padding */
    memset(x->cr_input, 0, pad * x->channels * sizeof(cc_s16l));
    memset(x->cr_input + (pad + in_frames) * x->channels, 0, pad * x->channels * sizeof(cc_s16l));

    /* Copy decoded MP3 samples */
    memcpy(x->cr_input + pad * x->channels, x->info.buffer,
           x->info.samples * sizeof(mp3d_sample_t));

    /* Estimate output size conservatively */
    x->resampled_frames = (size_t)ceil((double)in_frames * target_sr / x->source_samplerate) + 8;

    free(x->resampled_buffer);
    x->resampled_buffer = NULL;
    x->using_resampled = 0;

    x->resampled_buffer = (float *)malloc(x->resampled_frames * x->channels * sizeof(float));

    if (!x->resampled_buffer) {
        free(x->cr_input);
        x->cr_input = NULL;
        return 0;
    }

    cr_out_cb_data cb;
    cb.out = x->resampled_buffer;
    cb.frames_remaining = x->resampled_frames;
    cb.channels = x->channels;

    size_t frames_left = in_frames;

    ClownResampler_LowLevel_Resample(&x->cr_state, &x->cr_pre, x->cr_input + pad * x->channels,
                                     &frames_left, mar_cr_output_cb, &cb);

    x->resampled_frames -= cb.frames_remaining;
    x->total_frames = x->resampled_frames;
    x->using_resampled = 1;

    free(x->cr_input);
    x->cr_input = NULL;

    return 1;
}

// ─────────────────────────────────────
static void mar_clock_bang(t_mar_tilde *x) {
    outlet_bang(x->bang_out);
}

// ─────────────────────────────────────
static void mar_tilde_open(t_mar_tilde *x, t_symbol *s, int argc, t_atom *argv) {
    if (argc < 1 || argv[0].a_type != A_SYMBOL) {
        pd_error(x, "mar~ open: missing filename");
        return;
    }

    char dirbuf[MAXPDSTRING], *nameptr;
    char fullpath[MAXPDSTRING];

    int fd =
        canvas_open(x->canvas, atom_getsymbol(argv)->s_name, "", dirbuf, &nameptr, MAXPDSTRING, 1);

    if (fd < 0) {
        pd_error(x, "[mar~] cannot open file");
        return;
    }

    snprintf(fullpath, MAXPDSTRING, "%s/%s", dirbuf, nameptr);
    sys_close(fd);

    // Free previous resources
    if (x->file_loaded) {
        if (x->info.buffer) {
            free(x->info.buffer);
            x->info.buffer = NULL;
        }
        if (x->resampled_buffer) {
            free(x->resampled_buffer);
            x->resampled_buffer = NULL;
        }
    }

    x->using_resampled = 0;
    x->resampled_frames = 0;

    // Initialize decoder
    mp3dec_init(&x->mp3d);

    // Load MP3 file
    int error = mp3dec_load(&x->mp3d, fullpath, &x->info, NULL, NULL);
    if (error != 0 || x->info.buffer == NULL) {
        pd_error(x, "[mar~] MP3 decode failed");
        return;
    }

    x->file_loaded = 1;
    x->current_frame = 0;
    x->channels = x->info.channels;
    x->source_samplerate = x->info.hz;
    x->total_frames = x->info.samples / x->channels;
    x->using_resampled = 0;

    // Check if resampling is needed
    int target_sr = sys_getsr();
    if (x->source_samplerate != target_sr) {
        post("[mar~] Resampling from %d Hz to %d Hz", x->source_samplerate, target_sr);

        if (!mar_tilde_resample_audio(x)) {
            // Resampling failed, fall back to original
            pd_error(x, "[mar~] Resampling failed, using original sample rate");
            x->total_frames = x->info.samples / x->channels;
            x->using_resampled = 0;
        }
    } else {
        logpost(x, 3, "[mar~] Loaded %s - %d channels, %d Hz, %ld samples", nameptr, x->channels,
                x->source_samplerate, x->info.samples);
    }

    x->current_frame = 0;
    canvas_update_dsp();
}

// ─────────────────────────────────────
static void mar_tilde_float(t_mar_tilde *x, float f) {
    x->playing = f;
}

// ─────────────────────────────────────
static void mar_tilde_loop(t_mar_tilde *x, t_floatarg f) {
    x->loop = (f != 0);
}

// ─────────────────────────────────────
static t_int *mar_tilde_perform(t_int *w) {
    t_mar_tilde *x = (t_mar_tilde *)(w[1]);
    t_sample *out = (t_sample *)(w[2]);
    int n = (int)(w[3]);

    if (!x->file_loaded || !x->playing) {
        int ch = x->channels > 0 ? x->channels : 1;
        for (int i = 0; i < n; i++) {
            for (int c = 0; c < ch; c++) {
                out[c * n + i] = 0.0f;
            }
        }
        return w + 4;
    }

    int ch = x->channels;
    size_t total_frames = x->total_frames;
    mp3d_sample_t *buffer = x->info.buffer;
    float *resampled = x->resampled_buffer;
    int i = 0;

    while (i < n) {
        if (x->current_frame >= total_frames) {
            if (x->loop) {
                x->current_frame = 0;
            } else {
                x->playing = 0;
                clock_delay(x->clock, 0);
                while (i < n) {
                    for (int c = 0; c < ch; c++) {
                        out[c * n + i] = 0.0f;
                    }
                    i++;
                }
                break;
            }
        }

        size_t buffer_idx = x->current_frame * ch;
        if (x->using_resampled && resampled) {
            for (int c = 0; c < ch; c++) {
                out[c * n + i] = resampled[buffer_idx + c];
            }
        } else {
            for (int c = 0; c < ch; c++) {
                out[c * n + i] = (t_sample)buffer[buffer_idx + c] / 32768.0f;
            }
        }
        x->current_frame++;
        i++;
    }

    return w + 4;
}

// ─────────────────────────────────────
static void mar_tilde_dsp(t_mar_tilde *x, t_signal **sp) {
    signal_setmultiout(&sp[0], x->channels);
    x->block_size = sp[0]->s_n;
    dsp_add(mar_tilde_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
}

// ─────────────────────────────────────
static void mar_tilde_free(t_mar_tilde *x) {
    if (x->file_loaded && x->info.buffer) {
        free(x->info.buffer);
        x->info.buffer = NULL;
    }
    if (x->resampled_buffer) {
        free(x->resampled_buffer);
        x->resampled_buffer = NULL;
    }
}

// ─────────────────────────────────────
static void *mar_tilde_new(t_symbol *s, int argc, t_atom *argv) {
    t_mar_tilde *x = (t_mar_tilde *)pd_new(mar_tilde_class);

    // Initialize state
    x->file_loaded = 0;
    x->playing = 0;
    x->loop = 0;
    x->current_frame = 0;
    x->total_frames = 0;
    x->channels = 1;
    x->block_size = 64;
    x->canvas = canvas_getcurrent();
    x->clock = clock_new(x, (t_method)mar_clock_bang);
    x->resampled_buffer = NULL;
    x->resampled_frames = 0;
    x->using_resampled = 0;

    // Create outlets
    outlet_new(&x->x_obj, &s_signal);             // Audio outlet
    x->bang_out = outlet_new(&x->x_obj, &s_bang); // Bang when finished

    return (void *)x;
}

// ─────────────────────────────────────
extern "C" void mar_tilde_setup(void) {
    mar_tilde_class =
        class_new(gensym("mar~"), (t_newmethod)mar_tilde_new, (t_method)mar_tilde_free,
                  sizeof(t_mar_tilde), CLASS_DEFAULT, A_GIMME, 0);

    // DSP method
    class_addmethod(mar_tilde_class, (t_method)mar_tilde_dsp, gensym("dsp"), A_CANT, 0);

    // Control methods
    class_addmethod(mar_tilde_class, (t_method)mar_tilde_open, gensym("open"), A_GIMME, 0);
    class_addmethod(mar_tilde_class, (t_method)mar_tilde_loop, gensym("loop"), A_FLOAT, 0);
    class_addfloat(mar_tilde_class, (t_method)mar_tilde_float);
}
