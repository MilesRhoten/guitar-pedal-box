/*
TODO
Change to iterative fft like cooley tukey
it gets behind sometimes
increase size of buffer from int16* to 32 to amplitude higher
*/



#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <math.h>
#include <signal.h>
#include <termios.h>

#include "audio.h"

volatile sig_atomic_t keyPressed = 0;
volatile int recording = 0;

int main() {
    int rc;
    snd_pcm_t *capture_handle, *playback_handle;
    snd_pcm_hw_params_t *capture_params, *playback_params;
    unsigned int rate = SAMPLERATE; // Sample rate
    int channels = CHANNELS;          // Stereo
    snd_pcm_uframes_t frames = FRAMES;
    char *buffer;
    char *loopBuf8Bit;
    int size = frames * channels;

    
    // Open the PCM device for capture
    rc = snd_pcm_open(&capture_handle, PCM_CAPTURE_DEVICE,
		      SND_PCM_STREAM_CAPTURE, 0);
    if (rc < 0) {
        fprintf(stderr, "Unable to open capture PCM device: %s\n",
		snd_strerror(rc));
        return 1;
    }

    // Open the PCM device for playback
    rc = snd_pcm_open(&playback_handle, PCM_PLAYBACK_DEVICE,
		      SND_PCM_STREAM_PLAYBACK, 0);
    if (rc < 0) {
        fprintf(stderr, "Unable to open playback PCM device: %s\n",
		snd_strerror(rc));
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
    buffer = (char *)malloc(size * sizeof(int16_t));  // 2 bytes/sample
    if (buffer == NULL) {
      perror("malloc for buffer");
      exit(-1);
    }

    loopBuf8Bit = (char *)malloc(size * sizeof(int16_t));  // 2 bytes/sample
    if (loopBuf8Bit == NULL) {
      perror("malloc for buffer");
      exit(-1);
    }

    int16_t *loopBuf = (int16_t*) loopBuf8Bit;

    
    ComplexNum *complexBuf = (ComplexNum *) malloc(size * sizeof(ComplexNum));
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
    

    // set up interrupt
    setupInterrupt();


    // wav file
    int wavFile = -1;
    ssize_t bytesWritten = 0;

    WAVHeader wavHeader = {
      .chunk_id = "RIFF",
      .chunk_size = 0,  // Will be updated later
      .format = "WAVE",
      .subchunk1_id = "fmt ",
      .subchunk1_size = 16,
      .audio_format = 1,  // PCM
      .num_channels = channels,
      .sample_rate = rate,
      .byte_rate = rate * channels * sizeof(int16_t),
      .block_align = channels * sizeof(int16_t),
      .bits_per_sample = 16,
      .subchunk2_id = "data",
      .subchunk2_size = 0  // Will be updated later
    };
    
    
    // Capture, save, and playback loop
    printf("Capturing, saving, and playing back audio...\n");
    while (1) {
      if (keyPressed == ' ') {
	if (recording == 0) {
	  printf("Now recording!\n");
	  recording = 1;
	  wavFile = open("recorded_audio.wav", O_RDWR | O_CREAT | O_TRUNC, 0644);
	  if (wavFile == -1) {
	    perror("Failed to open wav file");
	    exit(-1);
	  }
	  write(wavFile, &wavHeader, sizeof(WAVHeader));
	}
	else if (recording == 1) {
	  printf("Now looping!\n");
	  // write wav header
	  lseek(wavFile, 0, SEEK_SET);
	  wavHeader.chunk_size = bytesWritten + 36;
	  wavHeader.subchunk2_size = bytesWritten;
	  write(wavFile, &wavHeader, sizeof(WAVHeader));

	  // seek to end of header
	  lseek(wavFile, sizeof(WAVHeader), SEEK_SET);
	  recording = 2;
	}
	else if (recording == 2) {
	  printf("Now doing nothing!\n");
	  recording = 0;
	  close(wavFile);
	  snd_pcm_prepare(capture_handle);
	  snd_pcm_prepare(playback_handle);
	}
	else {
	  perror("recording");
	  exit(-1);
	}
	keyPressed = 0;
      }
      
      
      // Capture audio
      rc = snd_pcm_readi(capture_handle, buffer, frames);
      if (rc == -EPIPE) {
	fprintf(stderr, "Capture overrun occurred\n");
            snd_pcm_prepare(capture_handle);
	    snd_pcm_prepare(playback_handle);
      } else if (rc < 0) {
	fprintf(stderr, "Error capturing audio: %s\n",
		snd_strerror(rc));
      }

      int16_t *audioBuf = (int16_t *) buffer;
      
      amplify_audio(audioBuf, frames * channels,
		      PREGAIN);
      /*
      fft_wrapper(audioBuf, complexBuf, twiddle_factors);
      
      // double freq = (double) getHighest(complexBuf);	
      
      // printf("Hz: %f\n", (freq * SAMPLERATE) / FRAMES);
      
      // eq_adjust(complexBuf, BASS, MID, TREBLE);
      
      ifft_wrapper(complexBuf, audioBuf, twiddle_factors);
      */

	if (recording == 1) {
	  ssize_t bytes = write(wavFile, audioBuf, size * sizeof(int16_t));
	  if (bytes == -1) {
	    perror("error in wav write123");
	    exit(-1);
	  }
	  bytesWritten = bytesWritten + bytes;
	}

	if (recording == 2) {
	  ssize_t readSize = read(wavFile, loopBuf, size * sizeof(int16_t));
	  if (readSize == 0) {
	    lseek(wavFile, sizeof(WAVHeader), SEEK_SET);
	    readSize = read(wavFile, loopBuf, size);
	  }

	  // combine buffer
	  int i;
	  int newNum;
	  for (i = 0; i < size; i++) {
	    newNum = (audioBuf[i] / 2) + (loopBuf[i] / 2);
	    if (newNum > INT16_MAX) {
	      audioBuf[i] = (int16_t) INT16_MAX;
	    }
	    else if (newNum < INT16_MIN) {
	      audioBuf[i] = (int16_t) INT16_MIN;
	    }
	    else {
	      audioBuf[i] = newNum;
	    }	    	  
	  }
	}
	
	
        // Playback audio
        rc = snd_pcm_writei(playback_handle, buffer, frames);
        if (rc == -EPIPE) {
            fprintf(stderr, "Playback underrun occurred\n");
            snd_pcm_prepare(playback_handle);
	    snd_pcm_prepare(capture_handle);
        } else if (rc < 0) {
            fprintf(stderr, "Error playing audio: %s\n", snd_strerror(rc));
        }
    }

    // Cleanup
    snd_pcm_drain(capture_handle);
    snd_pcm_close(capture_handle);
    snd_pcm_drain(playback_handle);
    snd_pcm_close(playback_handle);
    free(buffer);

    return 0;
}

void keypress_handler(int signo) {
  char c;
  if (read(STDIN_FILENO, &c, 1) == 1) {
    keyPressed = c;
  }
}

void setupInterrupt() {
  struct termios term;
    
  // Set terminal to raw mode to capture keypresses immediately
  tcgetattr(STDIN_FILENO, &term);
  term.c_lflag &= ~(ICANON | ECHO);  // Disable line buffering & echo
  tcsetattr(STDIN_FILENO, TCSANOW, &term);
  
  // Enable asynchronous key detection
  int flags = fcntl(STDIN_FILENO, F_GETFL);
  fcntl(STDIN_FILENO, F_SETFL, flags | O_NONBLOCK | O_ASYNC);
  fcntl(STDIN_FILENO, F_SETOWN, getpid());
  signal(SIGIO, keypress_handler);
}


// Precompute sine and cosine tables (scaled by SCALE)
void precompute_twiddle_factors(ComplexNum *twiddle_factors, int N) {
    for (int k = 0; k < N; k++) {
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

void eq_adjust(ComplexNum *complexBuf, int bass, int mid, int treble) {
  int bassLow = BASSLOWFREQ / FREQCONST;
  int bassMid = BASSMIDFREQ / FREQCONST;
  int midTreble = MIDTREBLEFREQ / FREQCONST;
  int trebleHi = TREBLEHIFREQ / FREQCONST;


  int i;
  for (i = bassLow; i < bassMid; i++) {
    complexBuf[i].real = complexBuf[i].real * bass;
    complexBuf[i].imag = complexBuf[i].imag * bass;
  }
  
  for (i = bassMid; i < midTreble; i++) {
    complexBuf[i].real = complexBuf[i].real * mid;
    complexBuf[i].imag = complexBuf[i].imag * mid;
  }
  
  for (i = midTreble; i < trebleHi; i++) {
    complexBuf[i].real = complexBuf[i].real * treble;
    complexBuf[i].imag = complexBuf[i].imag * treble;
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
      
      // Multiply twiddle factor with odd[k]
      ComplexNum part;
      part.imag = sin((-2.0 * PI * k) / N);
      part.real = cos((-2.0 * PI * k) / N);
      ComplexNum total;
      complexMul(&total, part, odd[k]);
      
      // Combine results
      x[k].real = even[k].real + total.real;
      x[k].imag = even[k].imag + total.imag;
      
      x[k + (N / 2)].real = even[k].real - total.real;
      x[k + (N / 2)].imag = even[k].imag - total.imag;
    }
    
    // Free memory for temporary arrays
    free(even);
    free(odd);
}

int32_t getMagnitude(ComplexNum num) {
  //return sqrt((num.real << 1) + (num.imag << 1));
  return (num.real * num.real) + (num.imag * num.imag);
}

int getHighest(ComplexNum *complexBuf) {
  int index = 0;
  int32_t highest = getMagnitude(complexBuf[0]);
  
  int i;
  for (i = 1; i < FRAMES; i++) {
    if (getMagnitude(complexBuf[i]) > highest) {
      highest = getMagnitude(complexBuf[i]);
      index = i;
    }
  }

  // printf("%d\n", highest);
  
  return index;
}

void complexMul(ComplexNum *answer, ComplexNum a, ComplexNum b) {
  answer -> real = (a.real * b.real) - (a.imag * b.imag);
  answer -> imag = (a.real * b.imag) + (a.imag * b.real);
}
