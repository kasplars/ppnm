CC = gcc
CFLAGS = -Wall -O -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt errors.txt neuralnetwork.png nndiffeq.png

out.txt errors.txt: main
	./main 1> out.txt 2> errors.txt
	
neuralnetwork.png: data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$f(x)$$";\
	set title "Neural network training on data set";\
	plot \
	 "$<" index 0 using 1:2 w p ps 0.4 pt 2 title "data for training"\
	,"$<" index 1 using 1:2 with line lt 1 title "trained function"\
	,"$<" index 1 using 1:3 with line lt 1 title "trained function differentiated"\
	,"$<" index 1 using 1:4 with line lt 1 title "trained function integrated"\
	,"$<" index 1 using 1:5 with line lt 2 title "true function"\
	'| pyxplot
	
nndiffeq.png: data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$f(x)$$";\
	set title "NN solution to the ODE $$y^{\prime\prime}+y=0$$ with boundary conditions at $$x=0$$: $$f(0)=0,f^{\prime}(0)=1$$";\
	plot \
	 "$<" index 2 using 1:2 with line lt 1 title "trained solution"\
	,"$<" index 2 using 1:3 with line lt 2 title "real solution"\
	'| pyxplot
	
data.txt: main

main: library.o ann.o

.PHONEY: default clean
	
clean: 
	$(RM) errors.txt *.o data.txt
