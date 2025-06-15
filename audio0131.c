/*
TODO
Change to iterative fft like cooley tukey
it gets behind sometimes
get magnitude for getHighest instead of only real
increase size of buffer from int16* to 32 to amplitude higher
*/



#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <math.h>

#include "audio.h"

#define PCM_CAPTURE_DEVICE "hw:3,0" // Replace with your capture device
#define PCM_PLAYBACK_DEVICE "hw:3,0" // Replace with your playback device

#define AMPLIFICATION_FACTOR 10

#define PI 3.14159265358979323846
#define SCALE (1 << 15) // Fixed-point scale factor (2^15)

#define SAMPLERATE 48000
#define CHANNELS 2
#define FRAMES (1 * 1024)

int main() {
    int rc;
    snd_pcm_t *capture_handle, *playback_handle;
    snd_pcm_hw_params_t *capture_params, *playback_params;
    unsigned int rate = SAMPLERATE; // Sample rate
    int channels = CHANNELS;          // Stereo
    snd_pcm_uframes_t frames = FRAMES;
    char *buffer;
    FILE *output_file;

    // Open the output file
    output_file = fopen("captured_audio.raw", "wb");
    if (!output_file) {
        perror("Failed to open output file");
        return 1;
    }

    // Open the PCM device for capture
    rc = snd_pcm_open(&capture_handle, PCM_CAPTURE_DEVICE,
		      SND_PCM_STREAM_CAPTURE, 0);
    if (rc < 0) {
        fprintf(stderr, "Unable to open capture PCM device: %s\n",
		snd_strerror(rc));
        fclose(output_file);
        return 1;
    }

    // Open the PCM device for playback
    rc = snd_pcm_open(&playback_handle, PCM_PLAYBACK_DEVICE,
		      SND_PCM_STREAM_PLAYBACK, 0);
    if (rc < 0) {
        fprintf(stderr, "Unable to open playback PCM device: %s\n",
		snd_strerror(rc));
        fclose(output_file);
        snd_pcm_close(capture_handle);
        return 1;
    }

    // Set up capture parameters
    snd_pcm_hw_params_malloc(&capture_params);
    snd_pcm_hw_params_any(capture_handle, capture_params);
    snd_pcm_hw_params_set_access(capture_handle, capture_params,
				 SND_PCM_ACCESS_RW_INTERLEAVED);
    snd_pcm_hw_params_set_format(capture_handle, capture_params,
				 SND_PCM_FORMAT_S16_LE);
    snd_pcm_hw_params_set_channels(capture_handle, capture_params,
				   channels);
    snd_pcm_hw_params_set_rate_near(capture_handle, capture_params,
				    &rate, NULL);
    snd_pcm_hw_params_set_period_size_near(capture_handle,
					   capture_params, &frames, NULL);
    snd_pcm_hw_params(capture_handle, capture_params);

    // Set up playback parameters
    snd_pcm_hw_params_malloc(&playback_params);
    snd_pcm_hw_params_any(playback_handle, playback_params);
    snd_pcm_hw_params_set_access(playback_handle, playback_params,
				 SND_PCM_ACCESS_RW_INTERLEAVED);
    snd_pcm_hw_params_set_format(playback_handle, playback_params,
				 SND_PCM_FORMAT_S16_LE);
    snd_pcm_hw_params_set_channels(playback_handle, playback_params,
				   channels);
    snd_pcm_hw_params_set_rate_near(playback_handle, playback_params,
				    &rate, NULL);
    snd_pcm_hw_params_set_period_size_near(playback_handle,
					   playback_params, &frames, NULL);
    snd_pcm_hw_params(playback_handle, playback_params);

    // Allocate buffer
    int size = frames * channels;
    buffer = (char *)malloc(size * sizeof(int16_t));  // 2 bytes/sample
    if (buffer == NULL) {
      perror("malloc for buffer");
      exit(-1);
    }

    ComplexNum *complexBuf =
      (ComplexNum *) malloc(size * sizeof(ComplexNum));
    if (complexBuf == NULL) {
      perror("malloc for complex buf");
      exit(-1);
    }


    ComplexNum *twiddle_factors =
      (ComplexNum *) malloc(size * sizeof(ComplexNum));
    if (twiddle_factors == NULL) {
      perror("malloc for twiddles");
      exit(-1);
    }

    // create twiddle factors
    precompute_twiddle_factors(twiddle_factors, FRAMES);
    
    
    // Capture, save, and playback loop
    printf("Capturing, saving, and playing back audio...\n");
    while (1) {
        // Capture audio
        rc = snd_pcm_readi(capture_handle, buffer, frames);
        if (rc == -EPIPE) {
            fprintf(stderr, "Capture overrun occurred\n");
            snd_pcm_prepare(capture_handle);
        } else if (rc < 0) {
            fprintf(stderr, "Error capturing audio: %s\n",
		    snd_strerror(rc));
        }

	
	amplify_audio((int16_t *) buffer, frames * channels,
		      AMPLIFICATION_FACTOR);
	
	fft_wrapper((int16_t *) buffer, complexBuf, twiddle_factors);

	double freq = (double) getHighest(complexBuf);

	
	printf("Hz: %f\n", (freq * SAMPLERATE) / FRAMES);

	//ifft_wrapper(complexBuf, (int16_t *) buffer, twiddle_factors);
	
	/*
	int i;
	for (i = 0; i < size / 4; i++) {
	  printf("Frame %d:\n", i / channels);
	  printf("Channel 1: %d\n", buf2[i * channels]);
	  printf("Channel 2: %d\n", buf2[(i * channels) + channels]);
	}
	*/

	/*
	int16_t *buf2 = (int16_t *) buffer;
	printf("Channel 1: %d\n", buf2[0]);
        */
	
        // Save audio to file
        if (fwrite(buffer, size, 1, output_file) != 1) {
            perror("Failed to write to output file");
            break;
        }

        // Playback audio
        rc = snd_pcm_writei(playback_handle, buffer, frames);
        if (rc == -EPIPE) {
            fprintf(stderr, "Playback underrun occurred\n");
            snd_pcm_prepare(playback_handle);
        } else if (rc < 0) {
            fprintf(stderr, "Error playing audio: %s\n", snd_strerror(rc));
        }
    }

    // Cleanup
    fclose(output_file);
    snd_pcm_drain(capture_handle);
    snd_pcm_close(capture_handle);
    snd_pcm_drain(playback_handle);
    snd_pcm_close(playback_handle);
    free(buffer);

    return 0;
}



// Precompute sine and cosine tables (scaled by SCALE)
void precompute_twiddle_factors(ComplexNum *twiddle_factors, int N) {
    for (int k = 0; k < N / 2; k++) {
        double angle = -2.0 * PI * k / N;
        twiddle_factors[k].real = (int32_t)(cos(angle) * SCALE);
        twiddle_factors[k].imag = (int32_t)(sin(angle) * SCALE);
    }
}

void amplify_audio(int16_t *samples, int num_samples, int factor) {
    for (int i = 0; i < num_samples; i++) {
        int32_t amplified = (int32_t)(samples[i] * factor);

        // Clip to prevent overflow
        if (amplified > INT16_MAX) {
            amplified = INT16_MAX;
        } else if (amplified < INT16_MIN) {
            amplified = INT16_MIN;
        }

        samples[i] = (int16_t)amplified;
    }
}

void distort_audio(int16_t *samples, int num_samples, int factor,
		   int16_t threshold) {
  int negThreshold = threshold * -1;

  int i;
  for (i = 0; i < num_samples; i++) {
    int32_t amplified = (int32_t) (samples[i] * factor);

    if (amplified > threshold) {
      amplified = threshold;
    }
    else if (amplified < negThreshold) {
      amplified = negThreshold;
    }
    
    samples[i] = (int16_t) amplified;    
  }
}

void ifft_wrapper(ComplexNum *complexBuf, int16_t *buffer, ComplexNum *twiddle_factors) {
  // complex conjugate
  int i;
  for (i = 0; i < FRAMES; i++) {
    complexBuf[i].imag = complexBuf[i].imag * -1;
  }

  fft_fixed(complexBuf, twiddle_factors, FRAMES);

  // take complex conjugate and normalize and write to buffer
  for (i = 0; i < FRAMES; i++) {
    // dont need to take conjugate becuase we ignore complex buffer from here
    // complexBuf[i].imag = complexBuf[i] * -1;
    
    complexBuf[i].real = complexBuf[i].real / FRAMES;

    buffer[i] = (int16_t) complexBuf[i].real;
  }

  
}

void fft_wrapper(int16_t *buffer, ComplexNum *complexBuf,
		 ComplexNum *twiddle_factors) {
  int i;
  for (i = 0; i < FRAMES; i++) {
    complexBuf[i].real = buffer[i];
    complexBuf[i].imag = 0;
  }

  fft_fixed(complexBuf, twiddle_factors, FRAMES);
}

// Function to perform the FFT using fixed-point arithmetic
void fft_fixed(ComplexNum *x, ComplexNum *twiddle_factors,
	       int N) {
    // Base case: If N == 1, no further computation is needed
    if (N <= 1) {
        return;
    }

    // Allocate memory for even and odd parts
    ComplexNum *even = (ComplexNum *)malloc(N / 2 * sizeof(ComplexNum));
    ComplexNum *odd = (ComplexNum *)malloc(N / 2 * sizeof(ComplexNum));

    // Split the input into even and odd indices
    for (int i = 0; i < N / 2; i++) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    // Recursive FFT calls for even and odd parts
    fft_fixed(even, twiddle_factors, N / 2);
    fft_fixed(odd, twiddle_factors, N / 2);

    // Combine step
    for (int k = 0; k < N / 2; k++) {
        ComplexNum twiddle = twiddle_factors[k * (N / 2) / (N / 2)];

        // Multiply twiddle factor with odd[k]
        int64_t t_real = ((int64_t)twiddle.real * odd[k].real - (int64_t)twiddle.imag * odd[k].imag) / SCALE;
        int64_t t_imag = ((int64_t)twiddle.real * odd[k].imag + (int64_t)twiddle.imag * odd[k].real) / SCALE;

        // Combine results
        x[k].real = even[k].real + (int32_t)t_real;
        x[k].imag = even[k].imag + (int32_t)t_imag;

        x[k + N / 2].real = even[k].real - (int32_t)t_real;
        x[k + N / 2].imag = even[k].imag - (int32_t)t_imag;
    }

    // Free memory for temporary arrays
    free(even);
    free(odd);
}

int getHighest(ComplexNum *complexBuf) {
  int index = 0;
  int32_t highest = complexBuf[0].real;

  int i;
  for (i = 1; i < FRAMES; i++) {
    if (complexBuf[i].real > highest) {
      highest = complexBuf[i].real;
      index = i;
    }
  }

  // printf("%d\n", highest);
  
  return index;
}
