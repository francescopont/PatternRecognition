###########################################################
# Description: The makefile for DDSMeval.
# Name: Mike Heath
# Date: 1/13/2000
###########################################################

CC = gcc
CFLAGS = -O2
BASELINE =  DDSMeval.o polyscan.o
INCLUDES = -I. -I../common
LIBS = -L../common -lmammo -lm

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c 

all: DDSMeval

clean:
	\rm DDSMeval DDSMeval.o polyscan.o

DDSMeval: $(BASELINE) ../common/libmammo.a
	$(CC) $(CFLAGS) $(INCLUDES) DDSMeval.o polyscan.o -o DDSMeval $(LIBS)

DDSMeval.o: DDSMeval.c ../common/breast.h ../common/mikesfileio.h ../common/virtual_image.h
polyscan.o: polyscan.c
