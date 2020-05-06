cc=gcc
CFLAGS=-Wall -g
LDLIBS=-lm

all: solver

clean:
	rm -rf solver
