# Makefile for Linux and FreeBSD 32-bit gwnum library
#

CC = gcc
CFLAGS = -I.. -O2 -march=i486 -malign-double

CPP = g++
CPPFLAGS = -I.. -I../qd -O2 -march=i486 -malign-double

OBJS = cpuid.o gwnum.o gwtables.o gwthread.o gwutil.o gwdbldbl.o giants.o ecmstag1.o

LIB = gwnum.a

#########################################################################

$(LIB): $(OBJS)
	cp linux/gwnum.a .
	ar -rs $(LIB) $(OBJS)

clean:
	rm -f $(OBJS)

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
