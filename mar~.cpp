#include <m_pd.h>
#include <m_imp.h>
#include <g_canvas.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <memory>
#include <vector>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#define MINIMP3_IMPLEMENTATION
#include <minimp3_ex.h>

#define MINIFLAC_IMPLEMENTATION
#include <miniflac.h>

// Ogg Vorbis
#define STB_VORBIS_IMPLEMENTATION
#include <stb_vorbis.c>

// wav & aiff
#include <AudioFile.h>

// Resampler
#include <CDSPResampler.h>

static t_class *mar_tilde_class;

// ╭─────────────────────────────────────╮
// │        UNIFIED AUDIO BUFFER         │
// ╰─────────────────────────────────────╯
struct AudioBuffer {
    int channels;
    int samplerate;
    size_t frames;
    // Planar, contiguous: data[ch * frames + frame]
    std::vector<float> data;

    AudioBuffer() : channels(0), samplerate(0), frames(0) {}

    void allocate(int ch, size_t fr, int sr) {
        channels = ch;
        frames = fr;
        samplerate = sr;
        if (ch <= 0 || fr == 0) {
            clear();
            return;
        }
        data.assign((size_t)ch * fr, 0.0f);
    }

    void clear() {
        channels = 0;
        samplerate = 0;
        frames = 0;
        data.clear();
        data.shrink_to_fit();
    }

    bool empty() const {
        return frames == 0 || channels == 0 || data.empty();
    }

    float *channel_ptr(int c) {
        return data.data() + ((size_t)c * frames);
    }

    const float *channel_ptr(int c) const {
        return data.data() + ((size_t)c * frames);
    }
};

typedef struct _mar_tilde {
    t_object x_obj;
    t_sample x_f;
    t_clock *clock;

    // Unified audio buffers
    AudioBuffer audio;
    AudioBuffer resampled;
    int using_resampled;
    bool resample;

    // Playback state
    size_t current_frame;
    int playing;
    int loop;

    // DSP
    int block_size;
    int out_channels;
    t_canvas *canvas;
    t_outlet *bang_out;
} t_mar_tilde;

// ╭─────────────────────────────────────╮
// │              RESAMPLER              │
// ╰─────────────────────────────────────╯
static bool resample_audio(const AudioBuffer &input, AudioBuffer &output,
                           double target_samplerate) {
    const double source_sr = (double)input.samplerate;
    const double target_sr = target_samplerate;

    if (target_sr <= 0.0 || source_sr <= 0.0 || input.empty()) {
        return false;
    }

    if (fabs(target_sr - source_sr) < 0.1) {
        // No resampling needed
        return false;
    }

    const int ch = input.channels;
    const size_t in_frames = input.frames;

    const double ratio = target_sr / source_sr;
    size_t estimated_frames = (size_t)ceil((double)in_frames * ratio) + 64;

    std::vector<std::unique_ptr<r8b::CDSPResampler>> resamplers;
    resamplers.reserve(ch);
    try {
        for (int c = 0; c < ch; ++c) {
            resamplers.emplace_back(new r8b::CDSPResampler(source_sr, target_sr, 1024));
        }
    } catch (...) {
        return false;
    }

    std::vector<std::vector<double>> temp_output(ch);
    std::vector<size_t> write_pos((size_t)ch, 0);
    for (int c = 0; c < ch; ++c) {
        temp_output[c].resize(estimated_frames);
    }

    const size_t block = 1024;
    std::vector<double> in_block(block);
    double *out_block = nullptr;

    // Process input blocks
    for (size_t pos = 0; pos < in_frames; pos += block) {
        const size_t n = std::min(block, in_frames - pos);

        for (int c = 0; c < ch; ++c) {
            for (size_t i = 0; i < n; ++i) {
                in_block[i] = (double)input.channel_ptr(c)[pos + i];
            }

            int produced = resamplers[c]->process(in_block.data(), (int)n, out_block);
            if (produced > 0) {
                const size_t need = write_pos[(size_t)c] + (size_t)produced;
                if (need > temp_output[c].size()) {
                    temp_output[c].resize(need + 64);
                }
                std::copy(out_block, out_block + produced,
                          temp_output[c].begin() + (ptrdiff_t)write_pos[(size_t)c]);
                write_pos[(size_t)c] = need;
            }
        }
    }

    // Flush resamplers
    bool flushing = true;
    while (flushing) {
        int produced = -1;
        for (int c = 0; c < ch; ++c) {
            int p = resamplers[c]->process(nullptr, 0, out_block);
            if (c == 0) {
                produced = p;
            }
            if (p > 0) {
                const size_t need = write_pos[(size_t)c] + (size_t)p;
                if (need > temp_output[c].size()) {
                    temp_output[c].resize(need + 64);
                }
                std::copy(out_block, out_block + p,
                          temp_output[c].begin() + (ptrdiff_t)write_pos[(size_t)c]);
                write_pos[(size_t)c] = need;
            }
        }
        if (produced <= 0) {
            flushing = false;
        }
    }

    // Verify all channels have same length
    size_t out_frames = write_pos.empty() ? 0 : write_pos[0];
    for (int c = 1; c < ch; ++c) {
        if (write_pos[(size_t)c] != out_frames) {
            return false;
        }
    }

    // Allocate output buffer
    output.allocate(ch, out_frames, (int)target_sr);

    // Copy data
    for (int c = 0; c < ch; ++c) {
        float *dst = output.channel_ptr(c);
        const double *src = temp_output[c].data();
        std::transform(src, src + (ptrdiff_t)out_frames, dst, [](double v) { return (float)v; });
    }

    return true;
}

// ─────────────────────────────────────
static void mar_clock_bang(t_mar_tilde *x) {
    outlet_bang(x->bang_out);
}

// ╭─────────────────────────────────────╮
// │             MP3 LOADER              │
// ╰─────────────────────────────────────╯
static bool load_mp3(t_mar_tilde *x, const char *fullpath) {
    mp3dec_t mp3d;
    mp3dec_file_info_t info;

    mp3dec_init(&mp3d);
    int error = mp3dec_load(&mp3d, fullpath, &info, NULL, NULL);

    if (error != 0 || info.buffer == NULL) {
        pd_error(x, "[mar~] MP3 decode failed");
        return false;
    }

    const int channels = info.channels;
    const int samplerate = info.hz;
    const size_t total_samples = info.samples;
    const size_t frames = total_samples / channels;

    // Allocate unified buffer
    x->audio.allocate(channels, frames, samplerate);

    // Deinterleave and convert to float [-1, 1]
    const mp3d_sample_t *src = info.buffer;
    for (size_t f = 0; f < frames; ++f) {
        for (int c = 0; c < channels; ++c) {
            x->audio.channel_ptr(c)[f] = (float)src[f * channels + c] / 32768.0f;
        }
    }

    free(info.buffer);

    logpost(x, 3, "[mar~] Loaded MP3: %d channels, %d Hz, %zu frames", channels, samplerate,
            frames);

    return true;
}

// ╭─────────────────────────────────────╮
// │          OGG VORBIS LOADER          │
// ╰─────────────────────────────────────╯
static bool load_ogg(t_mar_tilde *x, const char *fullpath) {
    int channels = 0;
    int samplerate = 0;
    short *interleaved = nullptr;

    // Returns number of samples PER CHANNEL, and allocates interleaved output with malloc().
    const int samples_per_channel =
        stb_vorbis_decode_filename(fullpath, &channels, &samplerate, &interleaved);

    if (samples_per_channel <= 0 || !interleaved || channels <= 0 || samplerate <= 0) {
        if (interleaved) {
            free(interleaved);
        }
        pd_error(x, "[mar~] OGG (Vorbis) decode failed");
        return false;
    }

    const size_t frames = (size_t)samples_per_channel;
    x->audio.allocate(channels, frames, samplerate);

    constexpr float scale = 1.0f / 32768.0f;
    for (size_t f = 0; f < frames; ++f) {
        const size_t base = f * (size_t)channels;
        for (int c = 0; c < channels; ++c) {
            x->audio.channel_ptr(c)[f] = (float)interleaved[base + (size_t)c] * scale;
        }
    }

    free(interleaved);

    logpost(x, 3, "[mar~] Loaded OGG (Vorbis): %d channels, %d Hz, %zu frames", channels,
            samplerate, frames);
    return true;
}

// ╭─────────────────────────────────────╮
// │             FLAC LOADER             │
// ╰─────────────────────────────────────╯
static bool load_flac(t_mar_tilde *x, const char *fullpath) {
    FILE *file = fopen(fullpath, "rb");
    if (!file) {
        pd_error(x, "[mar~] Cannot open FLAC file");
        return false;
    }

    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    fseek(file, 0, SEEK_SET);

    uint8_t *flac_data = (uint8_t *)malloc(file_size);
    if (!flac_data) {
        fclose(file);
        pd_error(x, "[mar~] Memory allocation failed");
        return false;
    }

    size_t bytes_read = fread(flac_data, 1, file_size, file);
    fclose(file);

    if (bytes_read != (size_t)file_size) {
        free(flac_data);
        pd_error(x, "[mar~] Failed to read FLAC file");
        return false;
    }

    miniflac_t decoder;
    miniflac_init(&decoder, MINIFLAC_CONTAINER_UNKNOWN);

    int32_t *temp_samples[8] = {NULL};
    int channels = 0;
    int samplerate = 0;
    uint8_t bps = 0;
    size_t total_samples = 0;
    bool first_frame = true;

    uint32_t pos = 0;
    uint32_t length = file_size;
    uint32_t used = 0;
    MINIFLAC_RESULT res;

    // Allocate temporary sample buffers
    for (int i = 0; i < 8; i++) {
        temp_samples[i] = (int32_t *)malloc(sizeof(int32_t) * 65535);
        if (!temp_samples[i]) {
            pd_error(x, "[mar~] Memory allocation failed for sample buffers");
            for (int j = 0; j < i; j++) {
                free(temp_samples[j]);
            }
            free(flac_data);
            return false;
        }
    }

    while ((res = miniflac_decode(&decoder, &flac_data[pos], length, &used, temp_samples)) ==
           MINIFLAC_OK) {
        if (first_frame) {
            channels = miniflac_frame_channels(&decoder);
            samplerate = miniflac_frame_sample_rate(&decoder);
            bps = miniflac_frame_bps(&decoder);
            first_frame = false;
        }

        uint16_t block_size = miniflac_frame_block_size(&decoder);
        total_samples += block_size * channels;
        length -= used;
        pos += used;
        miniflac_sync(&decoder, &flac_data[pos], length, &used);
        pos += used;
        length -= used;
    }

    if (channels == 0 || samplerate == 0 || total_samples == 0) {
        pd_error(x, "[mar~] Invalid FLAC file or decode failed");
        for (int i = 0; i < 8; i++) {
            free(temp_samples[i]);
        }
        free(flac_data);
        return false;
    }

    size_t frames = total_samples / channels;
    x->audio.allocate(channels, frames, samplerate);
    pos = 0;
    length = file_size;
    size_t current_frame = 0;
    miniflac_init(&decoder, MINIFLAC_CONTAINER_UNKNOWN);

    while ((res = miniflac_decode(&decoder, &flac_data[pos], length, &used, temp_samples)) ==
           MINIFLAC_OK) {
        uint16_t block_size = miniflac_frame_block_size(&decoder);
        float normalization_factor = 1.0f / (float)(1 << (bps - 1));
        for (uint16_t i = 0; i < block_size && current_frame + i < frames; i++) {
            for (int c = 0; c < channels; c++) {
                float normalized = (float)temp_samples[c][i] * normalization_factor;
                x->audio.channel_ptr(c)[current_frame + i] = normalized;
            }
        }
        current_frame += block_size;
        length -= used;
        pos += used;
        miniflac_sync(&decoder, &flac_data[pos], length, &used);
        pos += used;
        length -= used;
    }

    // Clean up
    for (int i = 0; i < 8; i++) {
        free(temp_samples[i]);
    }
    free(flac_data);
    logpost(x, 3, "[mar~] Loaded FLAC: %d channels, %d Hz, %zu frames, %d-bit", channels,
            samplerate, frames, bps);
    return true;
}

// ╭─────────────────────────────────────╮
// │           WAV/AIFF LOADER           │
// ╰─────────────────────────────────────╯
static bool load_wav_aiff(t_mar_tilde *x, const char *fullpath) {
    AudioFile<float> a;

    if (!a.load(fullpath)) {
        pd_error(x, "[mar~] Failed to load WAV/AIFF file");
        return false;
    }

    const int channels = a.getNumChannels();
    const int samplerate = a.getSampleRate();
    const size_t frames = a.getNumSamplesPerChannel();

    if (channels <= 0 || frames == 0) {
        pd_error(x, "[mar~] Invalid audio file format");
        return false;
    }

    // Allocate unified buffer
    x->audio.allocate(channels, frames, samplerate);

    // Copy per-channel data (AudioFile already stores non-interleaved)
    for (int c = 0; c < channels; ++c) {
        float *dst = x->audio.channel_ptr(c);
        for (size_t f = 0; f < frames; ++f) {
            dst[f] = a.samples[c][f];
        }
    }

    logpost(x, 3, "[mar~] Loaded %s: %d channels, %d Hz, %zu frames, %d-bit", fullpath, channels,
            samplerate, frames, a.getBitDepth());

    return true;
}

// ─────────────────────────────────────
static void mar_tilde_open(t_mar_tilde *x, t_symbol *s, int argc, t_atom *argv) {
    if (argc < 1 || argv[0].a_type != A_SYMBOL) {
        pd_error(x, "[mar~] open: missing filename");
        canvas_update_dsp();
        return;
    }

    x->resampled.clear();
    x->audio.clear();

    char dirbuf[MAXPDSTRING], *nameptr;
    char fullpath[MAXPDSTRING];

    int fd =
        canvas_open(x->canvas, atom_getsymbol(argv)->s_name, "", dirbuf, &nameptr, MAXPDSTRING, 1);

    if (fd < 0) {
        pd_error(x, "[mar~]: %s: No such file or directory", atom_getsymbol(argv)->s_name);
        canvas_update_dsp();
        return;
    }

    snprintf(fullpath, MAXPDSTRING, "%s/%s", dirbuf, nameptr);
    logpost(x, 3, "Opening %s", fullpath);
    sys_close(fd);

    // Clear existing buffers
    x->audio.clear();
    x->resampled.clear();
    x->using_resampled = 0;
    x->current_frame = 0;
    x->playing = 0;

    const char *dot = strrchr(fullpath, '.');
    bool loaded = false;

    const auto t0 = std::chrono::steady_clock::now();

    auto report_time = [&]() {
        const auto t1 = std::chrono::steady_clock::now();
        const auto us = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
        logpost(x, 3, "[mar~] open took %lld µs", (long long)us);
    };

    if (dot && !strcasecmp(dot, ".mp3")) {
        loaded = load_mp3(x, fullpath);
    } else if (dot && (!strcasecmp(dot, ".wav") || !strcasecmp(dot, ".aiff") ||
                       !strcasecmp(dot, ".aif"))) {
        loaded = load_wav_aiff(x, fullpath);
    } else if (dot && !strcasecmp(dot, ".flac")) {
        loaded = load_flac(x, fullpath);
    } else if (dot && !strcasecmp(dot, ".ogg")) {
        loaded = load_ogg(x, fullpath);
    } else {
        x->playing = false;
        pd_error(x, "[mar~] Supported formats: .mp3, .wav, .aiff, .aif, .ogg, .flac");
        canvas_update_dsp();
        return;
    }
    report_time();

    if (!loaded) {
        x->playing = false;
        canvas_update_dsp();
        return;
    }

    // Resample if needed
    int target_sr = sys_getsr();
    if (x->audio.samplerate != target_sr && x->resample) {
        post("[mar~] Resampling from %d Hz to %d Hz", x->audio.samplerate, target_sr);

        if (resample_audio(x->audio, x->resampled, (double)target_sr)) {
            x->using_resampled = 1;
        } else {
            pd_error(x, "[mar~] Resampling failed, using original sample rate");
            x->using_resampled = 0;
        }
    }

    canvas_update_dsp();
}

// ─────────────────────────────────────
static void mar_tilde_bang(t_mar_tilde *x) {
    x->playing = 1;
    x->current_frame = 0;
}

// ─────────────────────────────────────
static void mar_tilde_float(t_mar_tilde *x, float f) {
    x->playing = (f != 0);
    x->current_frame = 0;
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

    const AudioBuffer &buf = x->using_resampled ? x->resampled : x->audio;
    const int out_ch = x->out_channels > 0 ? x->out_channels : 1;
    const int buf_ch = buf.channels > 0 ? buf.channels : 0;
    const int copy_ch = std::min(out_ch, buf_ch);

    if (buf.empty() || !x->playing) {
        std::fill(out, out + ((size_t)out_ch * (size_t)n), 0.0f);
        return w + 4;
    }

    const size_t total_frames = buf.frames;

    int out_pos = 0;
    while (out_pos < n) {
        if (x->current_frame >= total_frames) {
            clock_delay(x->clock, 0);
            if (x->loop) {
                x->current_frame = 0;
            } else {
                x->playing = 0;
                // Zero remainder
                for (int c = 0; c < out_ch; ++c) {
                    std::fill(out + ((size_t)c * (size_t)n) + (size_t)out_pos,
                              out + ((size_t)c * (size_t)n) + (size_t)n, 0.0f);
                }
                break;
            }
        }

        const size_t remaining = total_frames - x->current_frame;
        const int to_copy = (int)std::min<size_t>((size_t)(n - out_pos), remaining);

        for (int c = 0; c < copy_ch; ++c) {
            const float *src = buf.channel_ptr(c) + x->current_frame;
            float *dst = out + ((size_t)c * (size_t)n) + (size_t)out_pos;
            std::copy_n(src, to_copy, dst);
        }
        for (int c = copy_ch; c < out_ch; ++c) {
            float *dst = out + ((size_t)c * (size_t)n) + (size_t)out_pos;
            std::fill(dst, dst + (size_t)to_copy, 0.0f);
        }

        x->current_frame += (size_t)to_copy;
        out_pos += to_copy;
    }

    return w + 4;
}

// ─────────────────────────────────────
static void mar_tilde_dsp(t_mar_tilde *x, t_signal **sp) {
    const AudioBuffer &buf = x->using_resampled ? x->resampled : x->audio;
    int ch = buf.channels > 0 ? buf.channels : 1;

    signal_setmultiout(&sp[0], ch);
    x->out_channels = ch;
    x->block_size = sp[0]->s_n;
    dsp_add(mar_tilde_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
}

// ─────────────────────────────────────
static void mar_tilde_free(t_mar_tilde *x) {
    x->audio.clear();
    x->resampled.clear();
    if (x->clock) {
        clock_free(x->clock);
    }
}

// ─────────────────────────────────────
static void *mar_tilde_new(t_symbol *s, int argc, t_atom *argv) {
    t_mar_tilde *x = (t_mar_tilde *)pd_new(mar_tilde_class);

    x->playing = 0;
    x->loop = 0;
    x->current_frame = 0;
    x->using_resampled = 0;
    x->block_size = 64;
    x->out_channels = 1;
    x->canvas = canvas_getcurrent();
    x->clock = clock_new(x, (t_method)mar_clock_bang);
    x->resample = true;
    x->loop = false;

    int i = 0;
    while (i < argc) {
        if (argv[i].a_type == A_SYMBOL) {
            t_symbol *sym = atom_getsymbol(&argv[i]);
            if (strcmp("-loop", sym->s_name) == 0) {
                x->loop = true;
            } else if (strcmp("-nor", sym->s_name) == 0) {
                x->resample = false;
            } else {
                const char *dot = strrchr(sym->s_name, '.');
                bool isaudio = false;
                t_atom file[1];
                if (dot && !strcasecmp(dot, ".mp3")) {
                    isaudio = true;
                } else if (dot && !strcasecmp(dot, ".wav")) {
                    isaudio = true;
                } else if (dot && !strcasecmp(dot, ".aiff")) {
                    isaudio = true;
                } else if (dot && !strcasecmp(dot, ".aif")) {
                    isaudio = true;
                } else if (dot && !strcasecmp(dot, ".flac")) {
                    isaudio = true;
                } else if (dot && !strcasecmp(dot, ".ogg")) {
                    isaudio = true;
                }
                if (isaudio) {
                    SETSYMBOL(&file[0], sym);
                    mar_tilde_open(x, gensym("open"), 1, file);
                }
            }
        } else if (argv[i].a_type == A_FLOAT) {
            x->playing = atom_getfloat(&argv[i]) == 1;
        }

        i++;
    }

    outlet_new(&x->x_obj, &s_signal);
    x->bang_out = outlet_new(&x->x_obj, &s_bang);

    return (void *)x;
}

// ─────────────────────────────────────
extern "C" void mar_tilde_setup(void) {
    mar_tilde_class =
        class_new(gensym("mar~"), (t_newmethod)mar_tilde_new, (t_method)mar_tilde_free,
                  sizeof(t_mar_tilde), CLASS_DEFAULT, A_GIMME, 0);

    class_addmethod(mar_tilde_class, (t_method)mar_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(mar_tilde_class, (t_method)mar_tilde_open, gensym("open"), A_GIMME, 0);
    class_addmethod(mar_tilde_class, (t_method)mar_tilde_loop, gensym("loop"), A_FLOAT, 0);
    class_addfloat(mar_tilde_class, (t_method)mar_tilde_float);
    class_addbang(mar_tilde_class, (t_method)mar_tilde_bang);
    logpost(nullptr, 3, "[mar~] version 0.1.0");
}
