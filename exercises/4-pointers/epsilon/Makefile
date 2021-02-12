CC = gcc
CFLAGS = -std=gnu11
CFLAGS += -Wall
LDFLAGS = -I.
LDLIBS = -lm

out.txt: epsilon
	./$< > $@

clean:
	$(RM) epsilon *.o
