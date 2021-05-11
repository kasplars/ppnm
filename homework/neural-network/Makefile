CC = gcc
CFLAGS = -Wall -O -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt errors.txt neuralnetwork.png

out.txt errors.txt: main
	./main 1> out.txt 2> errors.txt
	
neuralnetwork.png: data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$f(x)$$";\
	set title "Neural network training";\
	plot \
	 "$<" index 0 using 1:2 w p ps 0.4 pt 2 title "data for training"\
	,"$<" index 1 using 1:2 with line lt 1 title "trained function"\
	,"$<" index 1 using 1:3 with line lt 2 title "true function"\
	'| pyxplot
	
data.txt: main

main: library.o ann.o

.PHONEY: default clean
	
clean: 
	$(RM) errors.txt *.o data.txt
