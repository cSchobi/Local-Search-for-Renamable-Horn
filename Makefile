cc=gcc
CFLAGS=-Wall -g
LDLIBS=-lm

SOURCES=solver.c
EXECUTABLES=$(SOURCES:%.c=%)

all: solver

debug: CFLAGS += -DDEBUG
debug: all

clean:
	rm -rf $(EXECUTABLES)
