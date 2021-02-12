CC = gcc
CFLAGS = -O -std=gnu11
CFLAGS += -Wall
LDFLAGS = -I.
LDLIBS = -lm

out.txt: komplex
	./$< > $@

clean:
	$(RM) komplex *.o
