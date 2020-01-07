###########################################################
# Description: The makefile for combinesusp.
# Name: Mike Heath
# Date: 6/9/2000
###########################################################

CC = gcc
CFLAGS = -O2
BASELINE =  combinesusp.o detfeatures.o
INCLUDES = -I. -I../common
LIBS = -L../common -lmammo -lm

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c 

all: combinesusp

clean:
	\rm combinesusp combinesusp.o detfeatures.o

combinesusp: $(BASELINE) ../common/libmammo.a
	$(CC) $(CFLAGS) $(INCLUDES) combinesusp.o detfeatures.o -o combinesusp $(LIBS)

combinesusp.o: combinesusp.c ../common/breast.h ../common/mikesfileio.h ../common/virtual_image.h
detfeatures.o: detfeatures.c detfeatures.h ../common/mikesfileio.h
