###########################################################
# Description: The makefile for breastsegment.
# Name: Michael Heath, University of South Florida
# Date: 1/13/2000
###########################################################

CC = gcc
CFLAGS = -O2
BASELINE = breastsegment.o segment.o 
INCLUDES = -I. -I../common
LIBS = -L../common -lmammo -lm

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c 

all: breastsegment

clean:
	\rm breastsegment breastsegment.o segment.o

breastsegment: $(BASELINE) ../common/libmammo.a
	$(CC) $(CFLAGS) $(INCLUDES) breastsegment.o segment.o -o breastsegment $(LIBS)

breastsegment.o: breastsegment.c ../common/breast.h ../common/mikesfileio.h ../common/virtual_image.h
segment.o: segment.c ../common/breast.h ../common/virtual_image.h ../common/myalloc.h ../common/optical_density.h
