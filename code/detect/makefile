###########################################################
# Description: The makefile for detect.
# Name: Michael Heath, University of South Florida
# Date: 1/13/2000
###########################################################

CC = gcc
CFLAGS = -O2
BASELINE =  detect.o 
INCLUDES = -I.
LIBS = -lm

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c 

all: detect

clean:
	\rm detect detect.o

detect: $(BASELINE) 
	$(CC) $(CFLAGS) $(INCLUDES) detect.o -o detect $(LIBS)

detect.o: detect.c
