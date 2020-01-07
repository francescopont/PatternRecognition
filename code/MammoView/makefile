###########################################################
# Description: The makefile for MammoView.
# Name: Michael Heath, University of South Florida
# Date: 6/21/2000
###########################################################

CC = gcc
CFLAGS = -O2
BASELINE = MammoView.o simple_edge.o howteck_header.o poly.o
INCLUDES = -I. -I../common -I/usr/openwin/include
LIBS = -L../common -lmammo -L/usr/openwin/lib -lX11 -lsocket -lm

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c 

all: MammoView

clean:
	\rm *.o MammoView

MammoView: $(BASELINE) ../common/libmammo.a
	$(CC) $(CFLAGS) $(INCLUDES) MammoView.o simple_edge.o howteck_header.o poly.o -o MammoView $(LIBS)

MammoView.o: MammoView.c ../common/breast.h ../common/mikesfileio.h ../common/virtual_image.h ../common/optical_density.h xwindow.h
simple_edge.o: simple_edge.c
howteck_header.o: howteck_header.c
poly.o: poly.c
