###########################################################
# Description: The makefile for mktemplate.
# Name: Mike Heath
# Date: 1/13/2000
###########################################################

CC = gcc
CFLAGS = -O2
BASELINE =  mktemplate.o polyscan.o
INCLUDES = -I. -I../common
LIBS = -L../common -lmammo -lm

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c 

all: mktemplate

clean:
	\rm mktemplate mktemplate.o polyscan.o

mktemplate: $(BASELINE) ../common/libmammo.a
	$(CC) $(CFLAGS) $(INCLUDES) mktemplate.o polyscan.o -o mktemplate $(LIBS)

mktemplate.o: mktemplate.c ../common/breast.h ../common/mikesfileio.h ../common/virtual_image.h
polyscan.o: polyscan.c
