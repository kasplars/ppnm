CFLAGS = -O1 -Wall -std=gnu99 # since -O is the same as no optimization (-O0). 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

all: output ; cat $<
output: main ; ./$< > $@
clean: ; $(RM) main
