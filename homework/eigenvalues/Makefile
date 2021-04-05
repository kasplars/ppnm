CC = gcc
CFLAGS = -Wall -O1 -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt errors.txt
	./main 1> out.txt 2> errors.txt
	cat out.txt
	cat errors.txt

out.txt errors.txt: main

main: auxfile.o

.PHONEY: clean

clean:
	$(RM) main *.o errors.txt
