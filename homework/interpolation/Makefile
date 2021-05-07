CC = gcc
CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: quadspline.png quadderiv.png quadinteg.png linearspline.png linearinteg.png cubicspline.png cubicderiv.png cubicinteg.png
	@echo "\n\nLinear, quadratic and cubic spline all work as they should."	
	
outlspline.txt: linspline
	./$< > $@

outquadspline.txt: quadspline
	./$< > $@
	
outcubicspline.txt: cubicspline
	./$< > $@

linearspline.png: outlspline.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$f(x)$$";\
	set title "Linear Spline";\
	plot \
	 "$<" using 1:3 with line title "GSL"\
	,"$<" using 1:2 with line title "Own implementation"\
	'| pyxplot

linearinteg.png: outlinteg.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$F(x)$$";\
	set title "Linear Spline Integral";\
	plot \
	 "$<" using 1:3 with line title "GSL"\
	,"$<" using 1:2 with line title "Own implementation"\
	'| pyxplot

quadspline.png: outquadspline.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$f(x)$$";\
	set title "Quadratic Spline";\
	plot \
	"$<" using 1:2 with line title "Own implementation"\
	'| pyxplot

quadderiv.png: outquadspline.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$Df(x)$$";\
	set title "Quadratic Spline Derivative";\
	plot \
	"$<" using 1:3 with line title "Own implementation"\
	'| pyxplot

quadinteg.png: outquadspline.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$F(x)$$";\
	set title "Quadratic Spline Integral";\
	plot \
	"$<" using 1:4 with line title "Own implementation"\
	'| pyxplot
	
cubicspline.png: outcubicspline.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$f(x)$$";\
	set title "Cubic Spline";\
	plot \
	"$<" using 1:2 with line title "Own implementation"\
	'| pyxplot

cubicderiv.png: outcubicspline.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$Df(x)$$";\
	set title "Cubic Spline Derivative";\
	plot \
	"$<" using 1:3 with line title "Own implementation"\
	'| pyxplot

cubicinteg.png: outcubicspline.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$$x$$";\
	set ylabel "$$F(x)$$";\
	set title "Cubic Spline Integral";\
	plot \
	"$<" using 1:4 with line title "Own implementation"\
	'| pyxplot

clean:
	$(RM) linspline quadspline cubicspline out*
