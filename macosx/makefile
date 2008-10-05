# Makefile for Mac OS X mprime
#

ENVP= MACOSX_DEPLOYMENT_TARGET=10.3

CC = gcc
CFLAGS = -I.. -I../gwnum -I../linux -O2 -g -march=i686 -malign-double

CPP = g++
CPPFLAGS = -I.. -I../gwnum -I../linux -O2 -march=i686 -malign-double

LFLAGS = -Wl,-M
LIBS   = ../gwnum/gwnum.a -lm -lpthread -lcurl -lstdc++

FACTOROBJ = factor32.o
LINUXOBJS = prime.o menu.o
EXE      = mprime

#########################################################################

$(EXE): security $(LINUXOBJS) $(FACTOROBJ)
	ar -x ../gwnum/gwnum.a multutil.o
	$(ENVP) $(CC) $(LFLAGS) -o $(EXE) $(LINUXOBJS) $(FACTOROBJ) multutil.o $(LIBS)

security:
	[ ! -e ../security.h ] && touch ../security.h || true
	[ ! -e ../security.c ] && touch ../security.c || true
	[ ! -e ../secure5.c ] && touch ../secure5.c || true

clean:
	rm -f $(EXE) $(EXE2) $(LINUXOBJS)

.c.o:
	$(ENVP) $(CC) $(CFLAGS) -c $<

.cpp.o:
	$(ENVP) $(CPP) $(CPPFLAGS) -c $<
