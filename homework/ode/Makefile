CC = gcc
CFLAGS = -Wall -O -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt errors.txt shm.pdf sir.pdf threebody.pdf

out.txt: main
	./$< 1> $@

errors.txt: main
	./$< 2> $@

main: library.o

shm.pdf: shm.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$t$$";\
	set ylabel "$$y$$";\
	set title "Simple harmonic motion";\
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
	set title "SIR-model";\
	plot \
	 "$<" using 1:2 with line lt 1 title "$$S$$"\
	,"$<" using 1:3 with line lt 2 title "$$I$$ "\
	,"$<" using 1:4 with line lt 3 title "$$R$$ "\
	'| pyxplot

threebody.pdf: threebody.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key outside;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$y$$";\
	set xrange [-1.5:1.5];\
	set yrange [-1.5:1.5];\
	set title "Three-body problem: 8-shape solution";\
	plot \
	 "$<" using 2:5 with line lt 1 title "Body 0"\
	,"$<" using 3:6 with line lt 1 title "Body 1"\
	,"$<" using 4:7 with line lt 1 title "Body 2"\
	'| pyxplot


shm.txt: main

sir.txt: main

threebody.txt: main

.PHONEY: default

clean:
	$(RM) main *.o errors.txt
