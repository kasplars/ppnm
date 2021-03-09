CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

expfunc.pdf: expfunc.tex fig-pyxplot.pdf
	pdflatex $<
	pdflatex $<

out.txt: main
	./$< > $@

fig-pyxplot.pdf: fig-pyxplot.gpi out.txt
	pyxplot $<

clean:
	$(RM) *.log *.aux out* *.o fig*.pdf main
