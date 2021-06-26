CC = gcc
CFLAGS = -std=gnu11 -Wall
LDFLAGS = -I.
LDLIBS = -lm

default: outA.txt read-stdin read-file
	cat outA.txt

outA.txt: read-cmdline
	./$< 1 2 3 4 5 > $@

clean:
	$(RM) read-cmdline read-stdin read-file
