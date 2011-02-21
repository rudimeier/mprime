# Makefile for Linux 64-bit mprime
#

CC = gcc
CFLAGS = -I.. -I../gwnum -DX86_64 -O2 -Wno-unused-result

CPP = g++
CPPFLAGS = -I.. -I../gwnum -DX86_64 -O2

LFLAGS = -Wl,-M
LIBS = ../gwnum/gwnum.a ../gwnum/gwnum.ld -lm -lpthread -Wl,-Bstatic $(shell pkg-config --static --libs libcurl) -lstdc++ -Wl,-Bdynamic -ldl

FACTOROBJ = factor64.o
LINUXOBJS = prime.o menu.o
EXE      = mprime

#########################################################################

$(EXE): security $(LINUXOBJS) $(FACTOROBJ)
	$(CC) $(LFLAGS) -o $(EXE) $(LINUXOBJS) $(FACTOROBJ) $(LIBS)

security:
	[ ! -e ../security.h ] && touch ../security.h || true
	[ ! -e ../security.c ] && touch ../security.c || true
	[ ! -e ../secure5.c ] && touch ../secure5.c || true

clean:
	rm -f $(EXE) $(EXE2) $(LINUXOBJS) $(FACTOROBJ)

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
