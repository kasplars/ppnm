CC = gcc
CFLAGS = -Wall -O2 -std=gnu11 -DNDEBUG
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
	set title "Error comparison";\
	a=1; b=-1; t(x)=a*x**b;\
	c=1; d=-1; r(x)=c*x**d;\
	fit t(x) withouterrors "$<" using 1:3 via a,b;\
	fit r(x) withouterrors "$<" using 1:5 via c,d;\
	plot \
	 "$<" using 1:3 with points pt 2 ps 0.4 title "Pseudo-random"\
	,"$<" using 1:5 with points pt 2 ps 0.4 title "Quasi-random"\
	, t(x) with lines title "$$%.2fN^{%.2f}$$"%(a,b)\
	, r(x) with lines title "$$%.2fN^{%.2f}$$"%(c,d)\
	'| pyxplot


comparison.txt: main

main: functions.o

.PHONEY: default clean
	
clean: 
	$(RM) errors.txt *.o main comparison.txt
