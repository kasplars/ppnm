CC = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

out.txt: math
	./math > out.txt

math: math.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

math.o: math.c
	$(CC) $(CFLAGS) -c math.c -o math.o

clean:
	$(RM) *.out *.o math
