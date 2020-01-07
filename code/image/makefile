###########################################################
# Description: The makefile for mkimage.
# Name: Michael Heath, University of South Florida
# Date: 1/25/2000
###########################################################

CC = gcc
CFLAGS = -O2
BASELINE =  mkimage.o 
INCLUDES = -I. -I../common
LIBS = -L../common -lmammo -lm

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c 

all: mkimage

clean:
	\rm mkimage mkimage.o

mkimage: $(BASELINE) ../common/libmammo.a
	$(CC) $(CFLAGS) $(INCLUDES) mkimage.o -o mkimage $(LIBS)

mkimage.o: mkimage.c ../common/breast.h ../common/mikesfileio.h ../common/virtual_image.h ../common/optical_density.h
