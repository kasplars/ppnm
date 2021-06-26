CC = gcc
CFLAGS = -std=gnu11 -O -fwrapv
CFLAGS += -Wall
LDFLAGS = -I.
LDLIBS = -lm

default: out.txt
	cat out.txt

epsilon: epsilon.o part3.o

part3.o:

epsilon.o:

out.txt: epsilon
	./$< > $@

clean:
	$(RM) epsilon *.o
