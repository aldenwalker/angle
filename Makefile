UNAME = $(shell uname -s)
ifeq ($(UNAME),Darwin)
  CC=clang++
else
  CC=g++
endif
CFLAGS=-O3 #-g -Wall -Wextra -pedantic
IFLAGS=-I/usr/X11R6/include
LFLAGS=-L/usr/X11R6/lib -lX11

all: angle

lam.o: lam.cc lam.h
	$(CC) $(CFLAGS) $(IFLAGS) -c lam.cc 

angle.o: angle.cc angle_gui.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c angle.cc

angle_gui.o: angle_gui.cc angle_gui.h
	$(CC) $(CFLAGS) $(IFLAGS) -c angle_gui.cc

angle: angle.o angle_gui.o lam.o
	$(CC) -o angle angle.o angle_gui.o lam.o $(LFLAGS)

clean:
	rm *.o
	rm angle

