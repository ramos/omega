

CC=gcc
OPT=-Wall -ansi

INCL=-I/usr/include/
LIBS=-L/usr/lib/ -lgmp -lm

all:
	@make omega.x factoriales.x

omega.o: omega.c
	$(CC) $(OPT) $(INCL) -c -o omega.o omega.c

omega.x: omega.o
	$(CC) $(LIBS) -o omega.x omega.o

factoriales.o: factoriales.c
	$(CC) $(OPT) $(INCL) -c -o factoriales.o factoriales.c

factoriales.x: factoriales.o
	$(CC) $(LIBS) -o factoriales.x factoriales.o

