CC = gcc
CFLAGS = -Wall -O -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

all: out.txt errors.txt quantumenergy.png eigenfuncs.png outtime.png
	./main 1> out.txt 2> errors.txt
	cat out.txt
	cat errors.txt

out.txt errors.txt: main

quantumenergy.png: quantumenergyout.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$n$$";\
	set ylabel "$$E(n)$$ in units of $$2mL^2/\hbar^2$$";\
	set xrange [0:100];\
	set title "Quantum Particle in a Box";\
	plot \
	 "$<" using 1:2 with points pt 2 ps 0.3 title "calculated values"\
	,"$<" using 1:3 with points pt 1 ps 0.3 title "exact values"\
	'| pyxplot

eigenfuncs.png: eigenfuncs.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set tics out;\
	set xlabel "$$\xi$$";\
	set ylabel "$$u_n(\xi)$$";\
	set title "Eigenfunctions";\
	plot \
	 "$<" index 0 using 1:2 with line lt 2 color blue title "$$n=1$$ calculated eigenfunction"\
	,"$<" index 0 using 1:3 with line lt 1 lw 0.5 color blue title "$$n=1$$ exact eigenfunction"\
	,"$<" index 1 using 1:2 with line lt 2 color red title "$$n=2$$ calculated eigenfunction"\
	,"$<" index 1 using 1:3 with line lt 1 lw 0.5 color red title "$$n=2$$ exact eigenfunction"\
	,"$<" index 2 using 1:2 with line lt 2 color black title "$$n=3$$ calculated eigenfunction"\
	,"$<" index 2 using 1:3 with line lt 1 lw 0.5 color black title "$$n=3$$ exact eigenfunction"\
	'| pyxplot
	
outtime.png: outtime.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set tics out;\
	set xlabel "$$n$$";\
	set ylabel "seconds$$/n^{3}$$";\
	set title "Check if matrix diag. of $$n\\times n$$ matrix scales as $$O(n^{3})$$";\
	plot \
	 "$<" using 1:2 with points pt 3 ps 0.5 title "own implementation"\
	 ,"$<" using 1:3 with points pt 2 ps 0.5 title "gsl implementation"\
	 ,"$<" using 1:4 with points pt 1 ps 0.5 title "own optimized implementation"\
	'| pyxplot

quantumenergyout.txt: main

outtime.txt: main

eigenfuncs.txt: main

main: matrix.o

.PHONEY: clean all

clean:
	$(RM) main *.o errors.txt eigenfuncs.txt outtime.txt quantumenergyout.txt
