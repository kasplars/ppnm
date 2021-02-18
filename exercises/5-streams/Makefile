CC = gcc
CFLAGS = -std=gnu11 -Wall
LDFLAGS = -I.
LDLIBS = -lm

default: outA.txt outB.txt outC.txt read-cmdline read-stdin read-file
	cat outA.txt

outA.txt: read-cmdline
	./$< 1 2 3 4 5 > $@

outB.txt:

outC.txt:

clean:
	$(RM) outA.txt outC.txt read-cmdline read-stdin read-file
