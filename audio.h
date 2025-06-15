#ifndef __AUDIO_H__
#define __AUDIO_H__

#define PCM_CAPTURE_DEVICE "hw:3,0"

// use pulse for bluetooth
#define PCM_PLAYBACK_DEVICE "pulse"

#define PREGAIN 3
#define BASS 10
#define MID 1
#define TREBLE 1

#define BASSLOWFREQ 50
#define BASSMIDFREQ 250
#define MIDTREBLEFREQ (2 * 1024)
#define TREBLEHIFREQ (6 * 1024)

#define PI 3.14159265358979323846
#define SCALE (1 << 15)

#define SAMPLERATE 48000
#define CHANNELS 2
#define FRAMES (1 * 1024)

#define MAXSECREC 10

#define FREQCONST SAMPLERATE / FRAMES

typedef struct {
  double real;
  double imag;
} ComplexNum;

typedef struct {
    char chunk_id[4];           // "RIFF"
    uint32_t chunk_size;
    char format[4];             // "WAVE"
    char subchunk1_id[4];       // "fmt "
    uint32_t subchunk1_size;
    uint16_t audio_format;      // PCM = 1
    uint16_t num_channels;
    uint32_t sample_rate;
    uint32_t byte_rate;
    uint16_t block_align;
    uint16_t bits_per_sample;
    char subchunk2_id[4];       // "data"
    uint32_t subchunk2_size;
} WAVHeader;

void keypress_handler(int signo);

void setupInterrupt();

void precompute_twiddle_factors(ComplexNum *twiddle_factors, int N);

void amplify_audio(int16_t *samples, int num_samples, int factor);

void eq_adjust(ComplexNum *complexBuf, int bass, int mid, int treble);

void distort_audio(int16_t *samples, int num_samples, int factor,
		   int16_t threshold);

void ifft_wrapper(ComplexNum *complexBuf, int16_t *buffer, ComplexNum *twiddle_factors);

void fft_wrapper(int16_t *buffer, ComplexNum *complexBuf,
		 ComplexNum *twiddle_factors);

void fft_fixed(ComplexNum *x, ComplexNum *twiddle_factors, int N);

int getHighest(ComplexNum *complexBuf);

void complexMul(ComplexNum *answer, ComplexNum a, ComplexNum b);

#endif
