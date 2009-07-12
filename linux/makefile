# Makefile for Linux 32-bit mprime
#

CC = gcc
CFLAGS = -I.. -I../gwnum -O2 -march=i486 -malign-double

CPP = g++
CPPFLAGS = -I.. -I../gwnum -O2 -march=i486 -malign-double

LFLAGS = -Wl,-M
LIBS = ../gwnum/gwnum.a ../gwnum/gwnum.ld -lm -lpthread -Wl,-Bstatic -lcurl -lstdc++ -Wl,-Bdynamic -ldl

FACTOROBJ = factor32.o
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
	rm -f $(EXE) $(EXE2) $(LINUXOBJS)

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
