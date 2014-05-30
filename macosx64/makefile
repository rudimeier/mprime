# Makefile for Mac OS X 64-bit mprime
#

ENVP= MACOSX_DEPLOYMENT_TARGET=10.4

CC = gcc
CFLAGS = -I.. -I../gwnum -I../linux -O2 -DX86_64 -DCOMMAND_LINE_MPRIME -m64

CPP = g++
CPPFLAGS = -I.. -I../gwnum -I../linux -O2 -DX86_64 -DCOMMAND_LINE_MPRIME -m64

LFLAGS = -m64 -Wl,-no_pie
LIBS   = ../gwnum/amd64/release/gwnum.a -lm -lpthread -lcurl -framework IOKit -framework CoreFoundation -lstdc++

FACTOROBJ = ../prime95/macosx64/factor64.o
OBJS = prime.o menu.o
EXE = mprime

#########################################################################

$(EXE): security $(OBJS) $(FACTOROBJ)
	$(ENVP) $(CC) $(LFLAGS) -o $(EXE) $(OBJS) $(FACTOROBJ) $(LIBS)

security:
	[ ! -e ../security.h ] && touch ../security.h || true
	[ ! -e ../security.c ] && touch ../security.c || true
	[ ! -e ../secure5.c ] && touch ../secure5.c || true

clean:
	rm -f $(EXE) $(EXE2) $(OBJS)

.c.o:
	$(ENVP) $(CC) $(CFLAGS) -c $<

.cpp.o:
	$(ENVP) $(CPP) $(CPPFLAGS) -c $<
