CC = gcc
CFLAGS = -Wall -O1 -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt errors.txt olsplot.pyxplot.png
	./main 1> out.txt 2> errors.txt
	cat out.txt
	cat errors.txt

olsplot.pyxplot.png: plottingdata.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel ("Time in days");\
	set ylabel ("log(activity), rel. units");\
	set title "OLS";\
	plot \
	 "$<" using 1:2 with line title "OLS"\
	,"$<" using 3:4:5 with yerrorbars title "Data"\
	'| pyxplot

out.txt: main

errors.txt: main

.PHONEY: clean

clean:
	$(RM) main errors.txt plottingdata.txt
