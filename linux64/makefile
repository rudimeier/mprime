# Makefile for Linux 64-bit mprime
#
#  Ugh, different linux variants require different makefiles.
#  The current makefile is for CentOS 5.10.  We prefer to link against
#  older Linux versions because linking on the latest, greatest version
#  will create an mprime executable that will not run on older
#  Linux versions because of glibc incompatibilites.
#
#  Some linux versions require some of the variations below:
#	"export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig"
#	CFLAGS = -I.. -I../gwnum -DX86_64 -O2 -Wno-unused-result
#	LFLAGS = -Wl,-M
#	LIBS = ../gwnum/gwnum.a ../gwnum/gwnum.ld -lm -lpthread -Wl,-Bstatic $(shell pkg-config --static --libs libcurl) -lstdc++ -Wl,-Bdynamic -ldl

CC = gcc
CFLAGS = -I.. -I../gwnum -I/usr/local/include -DX86_64 -O2 

CPP = g++
CPPFLAGS = -I.. -I../gwnum -DX86_64 -O2

LFLAGS = -Wl,-M -Wl,-L/usr/local/lib
LIBS = ../gwnum/gwnum.a ../gwnum/gwnum.ld -lm -lpthread -Wl,-Bstatic -lcurl -Wl,-Bdynamic -lrt -lstdc++ -ldl

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
