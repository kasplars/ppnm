CC = gcc
CFLAGS = -O
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: hello
	./hello > out.txt

hello: hello.o
	$(CC) -o hello hello.o $(LDLIBS)

hello.o: hello.c
	$(CC) $(CFLAGS) -c hello.c

clean:
	$(RM) hello.o hello out.txt

test:
	echo $(LDLIBS)
	ech $(CC)
	echo $(RM)
