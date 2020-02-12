cc=gcc
CFLAGS=-Wall -g -DDEBUG

SOURCES=solver.c
EXECUTABLES=$(SOURCES:%.c=%)

all: $(EXECUTABLES)

clean:
	rm -rf $(EXECUTABLES)
