CC = gcc
CFLAGS = -Wall -O -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt errors.txt solution.png

out.txt errors.txt: main
	./main 1> out.txt 2> errors.txt

main: library.o

solution.png: odeoutput.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "$$r$$ [units of $$a_0$$]";\
	set ylabel "$$F_\\epsilon(r)$$";\
	set title "Radial wave function, Hartree and Bohr radius units";\
	plot \
	 "$<" index 0 using 1:2 with line lt 1 title "$$\\epsilon = -0.500315$$, $$r_{max}=8$$"\
	,"$<" index 1 using 1:2 with line lt 2 title "$$f_0(r)=r\\exp(-r)$$"\
	'| pyxplot
	
odeoutput.txt: main

.PHONEY: default clean
	
clean: 
	$(RM) errors.txt *.o *output.txt
