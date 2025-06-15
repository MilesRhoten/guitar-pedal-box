CC = gcc

CFLAGS = -Wall -g

LDFLAGS = -lasound -lm

all : audio

audio : audio.o
	$(CC) -o audio audio.o $(LDFLAGS)

audio.o : audio.c audio.h
	$(CC) $(CFLAGS) -c -o audio.o audio.c

clean :
	rm *.o *~
