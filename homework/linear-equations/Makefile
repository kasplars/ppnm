CC = gcc
CFLAGS = -Wall -O1 -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt errors.txt
	./qr-decomp 1> out.txt 2> errors.txt
	cat out.txt

out.txt: qr-decomp

errors.txt: qr-decomp

.PHONEY: clean

clean:
	$(RM) qr-decomp errors.txt
