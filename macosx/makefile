# Makefile for Mac OS X mprime
#

ENVP= MACOSX_DEPLOYMENT_TARGET=10.4

CC = gcc
CFLAGS = -I.. -I../gwnum -I../linux -O2 -m32 -march=i686 -DCOMMAND_LINE_MPRIME

CPP = g++
CPPFLAGS = -I.. -I../gwnum -I../linux -O2 -m32 -march=i686 -DCOMMAND_LINE_MPRIME

LFLAGS = -Wl,-no_pie
LIBS   = ../gwnum/release/gwnum.a -lm -lpthread -lcurl -framework IOKit -framework CoreFoundation -lstdc++

FACTOROBJ = ../prime95/macosx/factor32.o
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
