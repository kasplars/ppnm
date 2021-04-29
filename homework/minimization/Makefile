CC = gcc
CFLAGS = -Wall -O -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt errors.txt figure.png

figure.png: plotdata.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "Energy $$E$$ [GeV]";\
	set ylabel "Cross section $$\sigma$$ [some units]";\
	set title "Higgs Boson data";\
	plot \
	 "$<" index 0 using 1:2:3 with yerrorbars title "Data"\
	,"$<" index 1 using 1:2 with line lt 2 title "Breit-Wigner fit"\
	'| pyxplot
	
plotdata.txt: main

out.txt errors.txt: main
	./main 1> out.txt 2> errors.txt

main: library.o

.PHONEY: default clean
	
clean: 
	$(RM) errors.txt *.o plotdata.txt
