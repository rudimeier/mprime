# Makefile for Linux gwnum library
#

CC = gcc
CFLAGS = -I.. -O2 -march=i486 -malign-double

CPP = g++
CPPFLAGS = -I.. -I../qd -O2 -march=i486 -malign-double

LINUXOBJS = cpuid.o gwnum.o gwthread.o gwutil.o gwdbldbl.o giants.o ecmstag1.o

LIB = gwnum.a

#########################################################################

$(LIB): $(LINUXOBJS)
	ar -rs $(LIB) $(LINUXOBJS)

clean:
	rm -f $(LINUXOBJS)

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
