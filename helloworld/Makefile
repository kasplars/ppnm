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
	$(RM) hello.o hello 

upload: 
	git add --all
	git commit --all -m 'Initial commit'
	git push

test:
	echo $(LDLIBS)
	ech $(CC)
	echo $(RM)
