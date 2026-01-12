#include <m_pd.h>
#include <string.h>
#include <math.h>

#define MINIMP3_IMPLEMENTATION
#include <minimp3_ex.h>

static t_class *mar_tilde_class;

typedef struct _mar_tilde {
    t_object x_obj;
    t_sample x_f;

    // MP3 state
    mp3dec_t mp3d;
    mp3dec_file_info_t info;
    int file_loaded;

    // Playback state
    size_t current_frame;
    int playing;
    int loop;
    int channels;

    // resample
    float *resampled_buffer;
    size_t resampled_frames;
    int source_samplerate;

    // DSP
    int block_size;
    t_canvas *canvas;
    t_outlet *bang_out;
} t_mar_tilde;

// ─────────────────────────────────────
static void hermite_resample(float *input, float *output, int channels, size_t input_frames,
                             size_t output_frames, double step) {
    double pos = 0.0;

    for (size_t i = 0; i < output_frames; i++) {
        int idx = (int)pos;
        float frac = (float)(pos - (double)idx);

        // Hermite interpolation coefficients
        float c0, c1, c2, c3;

        // Calculate coefficients for smooth interpolation
        {
            float frac2 = frac * frac;
            float frac3 = frac2 * frac;

            c0 = -0.5f * frac3 + frac2 - 0.5f * frac;
            c1 = 1.5f * frac3 - 2.5f * frac2 + 1.0f;
            c2 = -1.5f * frac3 + 2.0f * frac2 + 0.5f * frac;
            c3 = 0.5f * frac3 - 0.5f * frac2;
        }

        for (int ch = 0; ch < channels; ch++) {
            // Get 4 samples for interpolation with bounds checking
            int idx0 = idx - 1;
            int idx1 = idx;
            int idx2 = idx + 1;
            int idx3 = idx + 2;

            // Clamp indices to valid range
            idx0 = (idx0 < 0) ? 0 : idx0;
            idx1 = (idx1 < 0) ? 0 : idx1;
            idx2 = (idx2 < 0) ? 0 : idx2;
            idx3 = (idx3 < 0) ? 0 : idx3;

            idx0 = (idx0 >= input_frames) ? input_frames - 1 : idx0;
            idx1 = (idx1 >= input_frames) ? input_frames - 1 : idx1;
            idx2 = (idx2 >= input_frames) ? input_frames - 1 : idx2;
            idx3 = (idx3 >= input_frames) ? input_frames - 1 : idx3;

            float x0 = input[idx0 * channels + ch];
            float x1 = input[idx1 * channels + ch];
            float x2 = input[idx2 * channels + ch];
            float x3 = input[idx3 * channels + ch];

            // Apply Hermite interpolation
            output[i * channels + ch] = x0 * c0 + x1 * c1 + x2 * c2 + x3 * c3;
        }

        pos += step;
    }
}

// ─────────────────────────────────────
static int mar_tilde_resample_audio(t_mar_tilde *x) {
    int target_sr = sys_getsr();
    if (target_sr <= 0) {
        pd_error(x, "[mar~] invalid system samplerate (%d)", target_sr);
        target_sr = x->source_samplerate;
    }
    if (target_sr <= 0 || x->source_samplerate <= 0 || target_sr == x->source_samplerate) {
        return 1;
    }

    size_t input_frames = x->info.samples / x->channels;
    double ratio = (double)target_sr / (double)x->source_samplerate;
    if (ratio <= 0.0) {
        pd_error(x, "[mar~] bad resample ratio");
        return 0;
    }
    double step = 1.0 / ratio;
    x->resampled_frames = (size_t)ceil((double)input_frames * ratio) + 4;

    // Allocate float buffer for resampling
    float *float_buffer = (float *)malloc(x->info.samples * sizeof(float));
    if (!float_buffer) {
        pd_error(x, "[mar~] Failed to allocate float conversion buffer");
        return 0;
    }

    // Convert int16 to float
    for (size_t i = 0; i < x->info.samples; i++) {
        float_buffer[i] = (float)x->info.buffer[i] / 32768.0f;
    }

    // Allocate resampled buffer
    x->resampled_buffer = (float *)malloc(x->resampled_frames * x->channels * sizeof(float));
    if (!x->resampled_buffer) {
        pd_error(x, "[mar~] Failed to allocate resampled buffer");
        free(float_buffer);
        return 0;
    }

    // Perform resampling
    hermite_resample(float_buffer, x->resampled_buffer, x->channels, input_frames,
                     x->resampled_frames, step);

    // Clean up
    free(float_buffer);

    // Replace original buffer with resampled data
    free(x->info.buffer);

    // Convert back to int16 for compatibility with existing code
    x->info.buffer =
        (mp3d_sample_t *)malloc(x->resampled_frames * x->channels * sizeof(mp3d_sample_t));
    if (!x->info.buffer) {
        pd_error(x, "[mar~] Failed to allocate final buffer");
        free(x->resampled_buffer);
        x->resampled_buffer = NULL;
        return 0;
    }

    // Convert float back to int16
    for (size_t i = 0; i < x->resampled_frames * x->channels; i++) {
        // Clamp to prevent overflow
        float sample = x->resampled_buffer[i];
        if (sample > 1.0f) {
            sample = 1.0f;
        }
        if (sample < -1.0f) {
            sample = -1.0f;
        }
        x->info.buffer[i] = (mp3d_sample_t)(sample * 32767.0f);
    }

    // Update sample count
    x->info.samples = x->resampled_frames * x->channels;
    x->info.hz = target_sr;
    x->source_samplerate = target_sr;

    // Free the float buffer (we keep the int16 version)
    free(x->resampled_buffer);
    x->resampled_buffer = NULL;

    return 1;
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

    // Check if resampling is needed
    int target_sr = sys_getsr();
    if (x->source_samplerate != target_sr) {
        post("[mar~] Resampling from %d Hz to %d Hz", x->source_samplerate, target_sr);

        if (!mar_tilde_resample_audio(x)) {
            // Resampling failed, fall back to original
            pd_error(x, "[mar~] Resampling failed, using original sample rate");
            post("[mar~] Loaded %s - %d channels, %d Hz, %ld samples", nameptr, x->channels,
                 x->source_samplerate, x->info.samples);
        } else {
            post("[mar~] Resampling complete: now %d Hz, %ld samples", target_sr, x->info.samples);
        }
    } else {
        post("[mar~] Loaded %s - %d channels, %d Hz, %ld samples", nameptr, x->channels,
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
    size_t total_frames = x->info.samples / ch;
    mp3d_sample_t *buffer = x->info.buffer;
    int i = 0;

    while (i < n) {
        if (x->current_frame >= total_frames) {
            if (x->loop) {
                x->current_frame = 0;
            } else {
                x->playing = 0;
                outlet_bang(x->bang_out);
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
        for (int c = 0; c < ch; c++) {
            out[c * n + i] = (t_sample)buffer[buffer_idx + c] / 32768.0f;
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
}

// ─────────────────────────────────────
static void *mar_tilde_new(t_symbol *s, int argc, t_atom *argv) {
    t_mar_tilde *x = (t_mar_tilde *)pd_new(mar_tilde_class);

    // Initialize state
    x->file_loaded = 0;
    x->playing = 0;
    x->loop = 0;
    x->current_frame = 0;
    x->channels = 1;
    x->block_size = 64;
    x->canvas = canvas_getcurrent();

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
