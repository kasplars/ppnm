CC = gcc
CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt quadspline.pyxplot.png quadderiv.pyxplot.png quadinteg.pyxplot.png lspline.pyxplot.png linteg.pyxplot.png
	cat out.txt

out.txt: linspline
	./$< > $@

outquadspline.txt: quadspline
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

linteg.pyxplot.png: outlinteg.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel texify("z");\
	set ylabel texify("F(z)");\
	set title "Linear Spline Integral";\
	plot \
	 "$<" using 1:3 with line title "GSL"\
	,"$<" using 1:2 with line title "Own implementation"\
	'| pyxplot

quadspline.pyxplot.png: outquadspline.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel texify("z");\
	set ylabel texify("f(z)");\
	set title "Quadratic Spline";\
	plot \
	"$<" using 1:2 with line title "Own implementation"\
	'| pyxplot

quadderiv.pyxplot.png: outquadspline.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel texify("z");\
	set ylabel texify("Df(z)");\
	set title "Quadratic Spline Derivative";\
	plot \
	"$<" using 1:3 with line title "Own implementation"\
	'| pyxplot

quadinteg.pyxplot.png: outquadspline.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel texify("z");\
	set ylabel texify("F(z)");\
	set title "Quadratic Spline Integral";\
	plot \
	"$<" using 1:4 with line title "Own implementation"\
	'| pyxplot

clean:
	$(RM) linspline quadspline
