/*******************************************************************************
* FILE: virtual_image.h
* Purpose: This code implements a read/write cache for images. The cache
* can support images in any of the common datatypes in 'C'.
* Name: Michael Heath, University of South Florida
* Date: 10/14/97
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#ifndef _VIRTUAL_IMAGE_
#define _VIRTUAL_IMAGE_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

/*******************************************************************************
* Definitions regarding datatypes.
*******************************************************************************/
#define UCHAR unsigned char
#define CHAR char
#define USHORT unsigned short int
#define SHORT short int
#define ULONG unsigned long int
#define LONG long int
#define FLOAT float
#define DOUBLE double

#define UCHARNUM 0
#define CHARNUM 1
#define USHORTNUM 2
#define SHORTNUM 3
#define ULONGNUM 4
#define LONGNUM 5
#define FLOATNUM 6
#define DOUBLENUM 7

static const int dt_size[8] = {1,1,2,2,4,4,4,8};
static const char *dt_name[] = {"UCHAR","CHAR","USHORT","SHORT","ULONG","LONG","FLOAT","DOUBLE"};

/*******************************************************************************
* Definitions for the cache.
*      r = rpage * NRB*RBS + rseg * RBS + rpix
*      c = cpage * NCB*CBS + cseg * CBS + cpix
*******************************************************************************/
#define INVALID_DATA 0
#define CLEAN_DATA 1
#define DIRTY_DATA 2

#define CACHE_READ_ONLY 0
#define CACHE_READ_WRITE 1

typedef struct{
   USHORT rpage, cpage;              /* Specifies which segment is in memeory. */
   USHORT rowlim, collim;            /* The number or colums and rows in this cache block. */
   USHORT dirtybit;                  /* INVALID_DATA, CLEAN_DATA, DIRTY_DATA */
   void *dataptr;                    /* A pointer to the data in this cach block. The memory is contiguous. */
}CACHESTRUCTURE;

typedef struct{
   FILE *fp;                         /* Pointer to the image data on disk. */
   char *filename;                   /* Name of the file containing the image being cached. */
   USHORT rbs, cbs,                  /* Number of rows and columns of pixels in a cache segment. */
          nrb, ncb;                  /* Number of rows and columns of cache segments in the cache. */
   USHORT rbs_bits, cbs_bits,        /* Number of bits of rows and columns of pixels in a cache segment. */
          nrb_bits, ncb_bits;        /* Number of bits of rows and columns of cache segments in the cache. */
   USHORT rbs_and_bits, cbs_and_bits,        /* Number of bits of rows and columns of pixels in a cache segment. */
          nrb_and_bits, ncb_and_bits;        /* Number of bits of rows and columns of cache segments in the cache. */
   USHORT depth;                     /* The depth of the cache. */
   int datatype;                     /* A number indicating the datatype of the data in the cache. */
   float resolution;                 /* Sample rate in microns of a (square) pixel. */
   USHORT swap_bytes;
   int headerbytes;                  /* The number of header bytes in an image file. */
   int rows, cols;                   /* Number of rows and columns in the image. */
   USHORT cachemode;                 /* The mode in which the cache is used (CACHE_READ_ONLY or CACHE_READ_WRITE). */
   USHORT filemode;                  /* The mode in which the cache file is used (CACHE_READ_ONLY or CACHE_READ_WRITE). */
   void *cachedata;                  /* The cache data itself. */
   CACHESTRUCTURE ***cachestructure; /* The structure that hold the cache tags, dirty bit etc. */
   double (*getpixel)();             /* A pointer to the function to get a pixel from the image. */
   void (*putpixel)();               /* A pointer to a function to set a pixel in the image. */
   void (*get_cache_segment)();      /* This is a pointer to the function for getting a segment of data into */
                                     /* the cache. This could be read from disk or could be computed. */
   void **sourceimlist;              /* A list of images that this procedure uses as input. */
   void *calculationdata;            /* Data used to get cache pages that are calculated. */

}CACHEIM;

/*******************************************************************************
* Function prototypes.
*******************************************************************************/
int allocate_cached_image(CACHEIM *image, int rows, int cols, int datatype,
    char *filename, char *mode, int max_depth, int max_cache_size, int headerbytes, USHORT swap_bytes,
    float resolution);

void deallocate_cached_image(CACHEIM *image);

void flush_cached_image(CACHEIM *image);

void read_cache_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg);
void write_cache_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg, int depth);
void print_virtual_image_information(CACHEIM image, char *varname);

/*
double get_cache_pixel_uchar(CACHEIM *image, int row, int col);
double get_cache_pixel_char(CACHEIM *image, int row, int col);
double get_cache_pixel_ushort(CACHEIM *image, int row, int col);
double get_cache_pixel_short(CACHEIM *image, int row, int col);
double get_cache_pixel_ulong(CACHEIM *image, int row, int col);
double get_cache_pixel_long(CACHEIM *image, int row, int col);
double get_cache_pixel_float(CACHEIM *image, int row, int col);
double get_cache_pixel_double(CACHEIM *image, int row, int col);

void put_cache_pixel_uchar(CACHEIM *image, int row, int col, double pixval);
void put_cache_pixel_char(CACHEIM *image, int row, int col, double pixval);
void put_cache_pixel_ushort(CACHEIM *image, int row, int col, double pixval);
void put_cache_pixel_short(CACHEIM *image, int row, int col, double pixval);
void put_cache_pixel_ulong(CACHEIM *image, int row, int col, double pixval);
void put_cache_pixel_long(CACHEIM *image, int row, int col, double pixval);
void put_cache_pixel_float(CACHEIM *image, int row, int col, double pixval);
void put_cache_pixel_double(CACHEIM *image, int row, int col, double pixval);
*/

#endif
