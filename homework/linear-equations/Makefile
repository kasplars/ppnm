CC = gcc
CFLAGS = -Wall -O2 -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: time.pyxplot.png out.txt errors.txt
	./qr-decomp 1> out.txt 2> errors.txt
	cat out.txt

out.txt: qr-decomp

errors.txt: qr-decomp

time.pyxplot.png: outtime.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel texify("N");\
	set ylabel "Seconds";\
	set title "Time measuring of QR-decomposition";\
	plot \
	 "$<" using 1:3 with line title "GSL"\
	,"$<" using 1:2 with line title "Own implementation"\
	'| pyxplot

.PHONEY: clean

clean:
	$(RM) qr-decomp errors.txt
