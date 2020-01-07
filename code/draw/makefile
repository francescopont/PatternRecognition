###########################################################
# Description: The makefile for drawimage.
# Name: Michael Heath, University of South Florida
# Date: 1/13/2000
###########################################################

CC = gcc
CFLAGS = -O2
BASELINE =  drawimage.o 
INCLUDES = -I. -I../common
LIBS = -L../common -lmammo -lm

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c 

all: drawimage

clean:
	\rm drawimage drawimage.o

drawimage: $(BASELINE) ../common/libmammo.a
	$(CC) $(CFLAGS) $(INCLUDES) drawimage.o -o drawimage $(LIBS)

drawimage.o: drawimage.c ../common/breast.h ../common/mikesfileio.h ../common/virtual_image.h
