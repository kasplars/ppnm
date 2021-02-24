CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONEY:default
default: out.errfunc.png

outA.txt: firstfunc
	./$< > $@

outB.txt: errfunc
	./$< > $@

out.errfunc.png: outB.txt Makefile
	echo ' \
	set terminal png ;\
	set key top left;\
	set tics out;\
	set output "$@";\
	set xlabel "$$x$$";\
	set title "Error function";\
	plot \
	"$<" using 1:2 with line title "erf($$x$$)"\' | pyxplot

.PHONEY:clean
clean:
	$(RM) *.o firstfunc errfunc
