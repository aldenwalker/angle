UNAME = $(shell uname -s)
ifeq ($(UNAME),Darwin)
  CC=clang++
else
  CC=g++
endif
CFLAGS=-g -Wall -Wextra -pedantic
IFLAGS=-I/usr/X11R6/include
LFLAGS=-L/usr/X11R6/lib -lX11

all: entropy

entropy.o: entropy.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c entropy.cc

entropy: entropy.o
	$(CC) -o entropy entropy.o $(LFLAGS)

clean:
	rm *.o
	rm entropy

