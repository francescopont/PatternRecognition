###########################################################
# Description: The makefile for MIASeval.
# Name: Mike Heath
# Date: 6/9/2000
###########################################################

CC = gcc
CFLAGS = -O2
BASELINE =  MIASeval.o polyscan.o
INCLUDES = -I. -I../common
LIBS = -L../common -lmammo -lm

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c 

all: MIASeval

clean:
	\rm MIASeval MIASeval.o polyscan.o

MIASeval: $(BASELINE) ../common/libmammo.a
	$(CC) $(CFLAGS) $(INCLUDES) MIASeval.o polyscan.o -o MIASeval $(LIBS)

MIASeval.o: MIASeval.c ../common/breast.h ../common/mikesfileio.h ../common/virtual_image.h
polyscan.o: polyscan.c
