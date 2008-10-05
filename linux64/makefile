# Makefile for Linux x86-64 mprime
#

CC = gcc
CFLAGS = -I.. -I../gwnum -DX86_64 -O2

CPP = g++
CPPFLAGS = -I.. -I../gwnum -DX86_64 -O2

LFLAGS = -Wl,-M
LIBS   = ../gwnum/gwnum.a ../gwnum/gwnum.ld -lm -lpthread -lcurl -lstdc++
LIBS2  = ../gwnum/gwnum.a ../gwnum/gwnum.ld -lm -lpthread -lcurl -lstdc++ -static

FACTOROBJ = factor64.o
LINUXOBJS = prime.o menu.o
EXE      = mprime
EXE2	 = sprime

#########################################################################

$(EXE): security $(LINUXOBJS) $(FACTOROBJ)
	$(CC) $(LFLAGS) -o $(EXE) $(LINUXOBJS) $(FACTOROBJ) $(LIBS)

$(EXE2): security $(LINUXOBJS) $(FACTOROBJ)
	$(CC) $(LFLAGS) -o $(EXE2) $(LINUXOBJS) $(FACTOROBJ) $(LIBS2)

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
