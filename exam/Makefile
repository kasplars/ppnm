CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt errors.txt
	cat out.txt

out.txt errors.txt: main
	./main 1> out.txt 2> errors.txt

.PHONEY: default clean

main: auxfile.o

clean: 
	$(RM) errors.txt *.o main
