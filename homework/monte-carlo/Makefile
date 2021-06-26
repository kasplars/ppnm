CC = gcc
CFLAGS = -Wall -O -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt errors.txt comparison.png

out.txt errors.txt: main
	./main 1> out.txt 2> errors.txt
	
comparison.png: comparison.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set tics out;\
	set xlabel "$$N$$";\
	set ylabel "error estimate";\
	set yrange [0:2.5];\
	set title "Error comparison";\
	plot \
	 "$<" using 1:3 with line lt 1 title "Pseudo: $$\sigma\\cdot\sqrt{N}\\cdot$$ scaling factor"\
	,"$<" using 1:5 with line lt 2 title "Quasi: $$\sigma \\cdot (\ln{N})^{d}/N\\cdot$$ scaling factor, d is dimensionality"\
	'| pyxplot


comparison.txt: main

.PHONEY: default clean
	
clean: 
	$(RM) errors.txt *.o main comparison.txt
