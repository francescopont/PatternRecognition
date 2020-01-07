/*******************************************************************************
* File: myalloc.h
* Overview: This file contains the function prototypes for the functions in
* mycalloc.c that we may want to call from other functions.
* Name: Michael Heath, University of South Florida
* Date: 7/29/99
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#ifndef _MYALLOC
#define _MYALLOC

#include <stdlib.h>

void *mycalloc(char *variable_name, size_t num, size_t size);
void *mymalloc(char *variable_name, size_t size);
void myfree(char *variable_name, void *ptr);
void mymemory_status();
void print_memory_inuse();
#endif
