CC=mpicc

CFLAGS=-c -W -Wall -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wcast-align -Wformat-security -Wmissing-format-attribute  -Wno-long-long -Wcast-align -Winline -Werror -pedantic -pedantic-errors -std=gnu99 -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wold-style-definition -Wdeclaration-after-statement -Wbad-function-cast -Wnested-externs -Wunused -Wuninitialized -fPIC -g

all: run

run: main.o aux.o solve.o
	$(CC) main.o aux.o solve.o -o run -lm

main.o: main.c
	$(CC) $(CFLAGS) main.c

aux.o: aux.c
	$(CC) $(CFLAGS) aux.c

solve.o: solve.c
	$(CC) $(CFLAGS) solve.c

clean:
	rm -rf *.o run
