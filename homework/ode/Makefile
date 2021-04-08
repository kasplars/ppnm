CC = gcc
CFLAGS = -Wall -O -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt errors.txt ode.pdf sir.pdf

out.txt: main
	./$< 1> $@

errors.txt: main
	./$< 2> $@

main: library.o

ode.pdf: ode-out.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$t$$";\
	set ylabel "$$y$$";\
	set title "ODE";\
	plot \
	 "$<" using 1:2 with line lt 1 title "$$y$$"\
	,"$<" using 1:3 with line lt 2 title "$$y^\prime$$ "\
	'| pyxplot

sir.pdf: sir.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key outside;\
	set tics out;\
	set xlabel "$$t$$";\
	set ylabel "Number";\
	set title "ODE";\
	plot \
	 "$<" using 1:2 with line lt 1 title "$$S$$"\
	,"$<" using 1:3 with line lt 2 title "$$I$$ "\
	,"$<" using 1:4 with line lt 3 title "$$R$$ "\
	'| pyxplot

ode-out.txt: main

sir.txt: main

.PHONEY: default

clean:
	$(RM) main *.o errors.txt
