CFLAGS = -Wall -O1 -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: all
all: out.txt

out.txt: main
	./$< > out.txt

.PHONEY: clean
clean:
	$(RM) main.o main

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
