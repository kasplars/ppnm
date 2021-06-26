CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
CFLAGS += -fopenmp
LDLIBS = -lm -lgomp

out.txt: montecarloA montecarloB
	./$< > $@
	./montecarloB >> $@
	
clean:
	$(RM) montecarloA montecarloB
