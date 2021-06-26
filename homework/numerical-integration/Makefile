CC = gcc
CFLAGS = -Wall -O -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt errors.txt

out.txt errors.txt: main
	./main 1> out.txt 2> errors.txt
	
.PHONEY: default clean
	
clean: 
	$(RM) errors.txt main *.o
