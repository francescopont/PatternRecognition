###########################################################
# Description: The makefile for afumfeature.
# Name: Michael Heath, University of South Florida
# Date: 1/13/2000
###########################################################

CC = gcc
CFLAGS = -O2
BASELINE =  afumfeature.o 
INCLUDES = -I. -I../common
LIBS = -L../common -lmammo -lm

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c 

all: afumfeature

clean:
	\rm afumfeature afumfeature.o

afumfeature: $(BASELINE) ../common/libmammo.a
	$(CC) $(CFLAGS) $(INCLUDES) afumfeature.o -o afumfeature $(LIBS)

afumfeature.o: afumfeature.c ../common/breast.h ../common/mikesfileio.h ../common/virtual_image.h ../common/optical_density.h
