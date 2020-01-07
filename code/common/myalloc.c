/*******************************************************************************
* File: mycalloc.c
* Overview: This code was written to help me in keeping track of memory.
* When the functions mycalloc() and mymalloc() are used, each block of
* allocated memory is kept track of with a linked list. When myfree() is used
* to free the memory, a check is done to be sure that memory was actually
* allocated at that address. A list of memory that is still allocated can be
* obtained by calling the function print_memory_inuse(). As you would expect,
* any memory not allocated by these functions is not kept track of at all.
* Name: Michael Heath, University of South Florida
* Date: 7/29/99
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include "myalloc.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

static unsigned long int bytes_allocated = 0;
static unsigned long int bytes_freed = 0;

/*******************************************************************************
* This structure is used locally (by this code) to keep track of memory that
* has been allocated.
*******************************************************************************/
struct MEMDATA{
   char *variable_name;
   unsigned long int address;
   unsigned long int bytes;
   int has_memory;
   struct MEMDATA *next;
};

static struct MEMDATA *start=NULL, *end=NULL;

void mymemory_status()
{
   printf("          Memory Status: total_bytes used = %lu bytes_freed = %lu bytes_currently_used = %lu\n",
      bytes_allocated, bytes_freed, bytes_allocated - bytes_freed);
   fflush(stdout);
}

void print_memory_inuse()
{
   struct MEMDATA *local=NULL;
   int i, numprintedchars;

   local = start;

   if(local == NULL) return;

   while(local != NULL){
      if(local->has_memory == 1){
          printf("          InUse: [%s]%n", local->variable_name, &numprintedchars);
          if(numprintedchars<50) for(i=0;i<(50-numprintedchars);i++) printf(" ");
          printf(" %10lu ", local->bytes);
          printf(" %p\n", (void *)(local->address));
      }
      local = local->next;
   }
}

void *mycalloc(char *variable_name, size_t num, size_t size)
{
   struct MEMDATA *newmemdata = NULL;
   void *newaddress = NULL;

   static void store_link(struct MEMDATA *memdata, struct MEMDATA **start, struct MEMDATA **end);

   newmemdata = (struct MEMDATA *) calloc(1, sizeof(struct MEMDATA));
   newmemdata->variable_name = (char *) calloc(strlen(variable_name)+1, sizeof(unsigned char));
   strcpy(newmemdata->variable_name, variable_name);
   newmemdata->bytes = (unsigned long int)num * (unsigned long int)size;

   bytes_allocated += (unsigned long int)num * (unsigned long int)size;

/*
   printf("          In mycalloc(%s, %lu): total_bytes used = %lu bytes_freed = %lu\n",
      variable_name, (unsigned long int)num * (unsigned long int)size,
      bytes_allocated, bytes_freed);
   fflush(stdout);
*/

   if((newaddress = (void *) calloc(num, size)) == NULL){
      fprintf(stderr, "Error callocing %s in mycalloc().\n", variable_name);
      exit(1);
   }

   newmemdata->address = (unsigned long) newaddress;
   newmemdata->has_memory = 1;
   store_link(newmemdata, &start, &end);

   return(newaddress);
}

void *mymalloc(char *variable_name, size_t size)
{
   struct MEMDATA *newmemdata = NULL;
   void *newaddress = NULL;

   static void store_link(struct MEMDATA *memdata, struct MEMDATA **start, struct MEMDATA **end);

   newmemdata = (struct MEMDATA *) calloc(1, sizeof(struct MEMDATA));
   newmemdata->variable_name = (char *) calloc(strlen(variable_name)+1, sizeof(unsigned char));
   strcpy(newmemdata->variable_name, variable_name);
   newmemdata->bytes = (unsigned long int)size;

   bytes_allocated += (unsigned long int)size;

/*
   printf("          In mymalloc(%s): total_bytes used = %lu bytes_freed = %lu\n",
      variable_name, bytes_allocated, bytes_freed);
   fflush(stdout);
*/

   newaddress = (void *)malloc(size);

   newmemdata->address = (unsigned long)newaddress;
   newmemdata->has_memory = 1;
   store_link(newmemdata, &start, &end);

   return(newaddress);
}

void myfree(char *variable_name, void *ptr)
{
   struct MEMDATA *memdata=NULL;

   static struct MEMDATA *find_link(struct MEMDATA *start, unsigned long mem_address);

   memdata = find_link(start, (unsigned long int)ptr);
   if(memdata == NULL){
      fprintf(stderr, "BIG TIME ERROR! TRYING TO FREE MEMORY FOR %s AT AN ADDRESS THAT WAS NOT ALLOCATED!!!\n",
         variable_name);
      exit(1);
   }
   fflush(stdout);

   free(ptr);
   memdata->has_memory = 0;
   bytes_freed += memdata->bytes;

/*
   printf("          In myfree(%s): total_bytes used = %lu bytes_freed = %lu\n",
      variable_name, bytes_allocated, bytes_freed);
   fflush(stdout);
*/
}

static void store_link(struct MEMDATA *memdata, struct MEMDATA **start, struct MEMDATA **end)
{
   if(*end == NULL){   /* first element in list */
      memdata->next = NULL;
      *end = memdata;
      *start = memdata;
      return;
   }

   (*end)->next = memdata;
   memdata->next = NULL;
   *end = memdata;
}

static struct MEMDATA *find_link(struct MEMDATA *start, unsigned long mem_address)
{
   if(start == NULL) return((struct MEMDATA *)NULL);

   while((start != NULL) && !((start->address == mem_address) && (start->has_memory == 1))) start = start->next;

   if(start != NULL) return(start);
   else return((struct MEMDATA *)NULL);
}
