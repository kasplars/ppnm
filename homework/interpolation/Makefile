CC = gcc
CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt lspline.pyxplot.png outinteg.pyxplot.png
	cat out.txt

out.txt: linspline
	./$< > $@

lspline.pyxplot.png: outlspline.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel texify("z");\
	set ylabel texify("f(z)");\
	set title "Linear Spline";\
	plot \
	 "$<" using 1:3 with line title "GSL"\
	,"$<" using 1:2 with line title "Own implementation"\
	'| pyxplot

outinteg.pyxplot.png: outinteg.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel texify("z");\
	set ylabel texify("f(z)");\
	set title "Linear Spline Integrator";\
	plot \
	 "$<" using 1:3 with line title "GSL"\
	,"$<" using 1:2 with line title "Own implementation"\
	'| pyxplot

clean:
	$(RM) linspline
