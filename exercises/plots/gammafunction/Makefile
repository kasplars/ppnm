#CFLAGS = `gsl-config --cflags`
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

gammaf.gnuplot.png: test.gpi
	cat test.gpi | sed 's/pyxplot/gnuplot/' | gnuplot

gammaf.pyxplot.png: data.txt Makefile tabulated-vals.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "gamma function from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-gamma"\
	,"$<" using 1:4 with line title "stirlings formula"\
	,"tabulated-vals.txt" using 1:2 with points title "factorial representation"\
	' | tee test.gpi | pyxplot

test.gpi:gammaf.pyxplot.png

data.txt: gammafunction
	./$< > $@

gammafunction: gammafunction.o stirling.o

.PHONEY:clean
clean:
	$(RM) *.o gammafunction data.txt test*
test:
	echo $(CFLAGS)
	echo $(LDLIBS)
