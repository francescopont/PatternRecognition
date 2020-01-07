/*******************************************************************************
* PROGRAM: virtual_image.c
* PURPOSE: This file contains the basic code for implementing a virtul image of
* any datatype.
* NAME: Michael Heath, University of South Florida
* Date: 10/14/97
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include "virtual_image.h"
#include "myalloc.h"
#include <string.h>

#define VERBOSEMODE 0

int get_csize(int x, int depth);
void find_cache_parameters(int rows, int cols, int max_depth, int max_size,
   int *a, int *b, int *depth, int bytes_per_pixel);

/*******************************************************************************
* Function prototypes for getting and putting pixels in virtual images.
* One set is for images that are completely in memory. The other set is
* for images that are only partially stored in the cache.
*******************************************************************************/
double get_fullyinmemory_cache_pixel_uchar(CACHEIM *image, int row, int col);
double get_fullyinmemory_cache_pixel_char(CACHEIM *image, int row, int col);
double get_fullyinmemory_cache_pixel_ushort(CACHEIM *image, int row, int col);
double get_fullyinmemory_cache_pixel_short(CACHEIM *image, int row, int col);
double get_fullyinmemory_cache_pixel_ulong(CACHEIM *image, int row, int col);
double get_fullyinmemory_cache_pixel_long(CACHEIM *image, int row, int col);
double get_fullyinmemory_cache_pixel_float(CACHEIM *image, int row, int col);
double get_fullyinmemory_cache_pixel_double(CACHEIM *image, int row, int col);
 
void put_fullyinmemory_cache_pixel_uchar(CACHEIM *image, int row, int col, double pixval);
void put_fullyinmemory_cache_pixel_char(CACHEIM *image, int row, int col, double pixval);
void put_fullyinmemory_cache_pixel_ushort(CACHEIM *image, int row, int col, double pixval);
void put_fullyinmemory_cache_pixel_short(CACHEIM *image, int row, int col, double pixval);
void put_fullyinmemory_cache_pixel_ulong(CACHEIM *image, int row, int col, double pixval);
void put_fullyinmemory_cache_pixel_long(CACHEIM *image, int row, int col, double pixval);
void put_fullyinmemory_cache_pixel_float(CACHEIM *image, int row, int col, double pixval);
void put_fullyinmemory_cache_pixel_double(CACHEIM *image, int row, int col, double pixval);
 
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
 
/*******************************************************************************
* A macro for getting a pixel of any datatype from an image that is held
* completely in the cache.
*******************************************************************************/
#define getpixel_fullyinmemory_macro(DTYPE){\
   return((double)(((DTYPE *)(image->cachedata))[row*(image->cols)+col]));\
}

/*******************************************************************************
* A macro for putting a pixel of any datatype from an image that is held
* completely in the cache.
*******************************************************************************/
#define putpixel_fullyinmemory_macro(DTYPE){\
   (((DTYPE *)(image->cachedata))[row*(image->cols)+col]) = (DTYPE)pixval;\
}

/*******************************************************************************
* A macro for getting a pixel of any datatype.
*******************************************************************************/
#define getpixel_macro(DTYPE){\
\
   int rpage, rseg, rpix, cpage, cseg, cpix;\
   CACHESTRUCTURE cstruct, cstruct2;\
   int d=0;\
\
   /****************************************************************************\
   * Translate the coordinates of the pixel into a cache location.\
   rseg = row/(image->rbs);\
   cseg = col/(image->cbs);\
   rpage = rseg / (image->nrb);\
   cpage = cseg / (image->ncb);\
   rseg = rseg % (image->nrb);\
   cseg = cseg % (image->ncb);\
   rpix = row % (image->rbs);\
   cpix = col % (image->cbs);\
   ****************************************************************************/\
   rseg = row >> (image->rbs_bits);\
   rpage = rseg >> (image->nrb_bits);\
   rseg = rseg & (image->nrb_and_bits);\
   rpix = row & (image->rbs_and_bits);\
\
   cseg = col >> (image->cbs_bits);\
   cpage = cseg >> (image->ncb_bits);\
   cseg = cseg & (image->ncb_and_bits);\
   cpix = col & (image->cbs_and_bits);\
\
   /****************************************************************************\
   * If the page we want is on the top of the stack we don't need to look any\
   * further for it. If this is true and the data is not INVALID_DATA then\
   * we can just get the pixel value and return it to the calling routine.\
   ****************************************************************************/\
   cstruct = image->cachestructure[rseg][cseg][0];\
   if((cstruct.cpage == cpage) && (cstruct.rpage == rpage) &&\
      (cstruct.dirtybit != INVALID_DATA)){\
      return((double)(((DTYPE *)(cstruct.dataptr))[(rpix<<(image->cbs_bits))+cpix]));\
   }\
\
   /****************************************************************************\
   * We need to look down the depth of the cache for this page. If we get to\
   * the bottom of the depth then we don't have the page in the cache and\
   * we will need to read it in from the disk file (or recompute it).\
   ****************************************************************************/\
   else{\
      for(d=1;d<(image->depth);d++){\
         cstruct2 = image->cachestructure[rseg][cseg][d];\
         image->cachestructure[rseg][cseg][d] = cstruct;\
         cstruct = cstruct2;\
\
         if((cstruct.cpage == cpage) && (cstruct.rpage == rpage)){\
            if(cstruct.dirtybit != INVALID_DATA){\
               /* printf("Found the desired information at depth = %d\n", d); */ \
               image->cachestructure[rseg][cseg][0] = cstruct;   /* put this page on top */\
               return((double)(((DTYPE *)(cstruct.dataptr))[(rpix<<(image->cbs_bits))+cpix]));\
            }\
            else d=image->depth;\
         }\
      }\
      image->cachestructure[rseg][cseg][0] = cstruct;   /* put this page on top */\
   }\
\
   /****************************************************************************\
   * We now know that the requested segment was not in the cache. Now we need to\
   * read in the segment from disk. If this is a read/write cache then we must\
   * first write out the page that was in the cache if it is dirty.\
   ****************************************************************************/\
   if((cstruct.dirtybit == DIRTY_DATA) && (image->filemode == CACHE_READ_WRITE))\
      write_cache_segment(image, cstruct.rpage, cstruct.cpage, rseg, cseg, 0);\
\
   /****************************************************************************\
   * Now we must read in the data that we want from disk.\
   ****************************************************************************/\
   (*(image->get_cache_segment))(image, rpage, cpage, rseg, cseg);\
   cstruct = image->cachestructure[rseg][cseg][0];\
\
   return((double)(((DTYPE *)(cstruct.dataptr))[(rpix<<(image->cbs_bits))+cpix]));\
}

/*******************************************************************************
* A macro for putting a pixel of any datatype.
*******************************************************************************/
#define putpixel_macro(DTYPE){\
\
   int rpage, rseg, rpix, cpage, cseg, cpix;\
   CACHESTRUCTURE cstruct, cstruct2;\
   int d=0;\
\
   /****************************************************************************\
   * Make sure that the cache is read/write. If not we have an error.\
   ****************************************************************************/\
   if(image->cachemode != CACHE_READ_WRITE){\
      fprintf(stderr, "ERROR! Trying to write to a read only cache.\n");\
      exit(1);\
   }\
\
   /****************************************************************************\
   * Translate the coordinates of the pixel into a cache location.\
   rseg = row/(image->rbs);\
   cseg = col/(image->cbs);\
   rpage = rseg / (image->nrb);\
   cpage = cseg / (image->ncb);\
   rseg = rseg % (image->nrb);\
   cseg = cseg % (image->ncb);\
   rpix = row % (image->rbs);\
   cpix = col % (image->cbs);\
   ****************************************************************************/\
   rseg = row >> (image->rbs_bits);\
   rpage = rseg >> (image->nrb_bits);\
   rseg = rseg & (image->nrb_and_bits);\
   rpix = row & (image->rbs_and_bits);\
\
   cseg = col >> (image->cbs_bits);\
   cpage = cseg >> (image->ncb_bits);\
   cseg = cseg & (image->ncb_and_bits);\
   cpix = col & (image->cbs_and_bits);\
\
   /****************************************************************************\
   * If the page we want is on the top of the stack we don't need to look any\
   * further for it. If this is true and the data is not INVALID_DATA then\
   * we can just get the pixel value and return it to the calling routine.\
   ****************************************************************************/\
   cstruct = image->cachestructure[rseg][cseg][0];\
   if((cstruct.cpage == cpage) && (cstruct.rpage == rpage) &&\
      (cstruct.dirtybit != INVALID_DATA)){\
      (((DTYPE *)(cstruct.dataptr))[(rpix<<(image->cbs_bits))+cpix]) = (DTYPE)pixval;\
      image->cachestructure[rseg][cseg][0].dirtybit = DIRTY_DATA;\
      return;\
   }\
\
   /****************************************************************************\
   * We need to look down the depth of the cache for this page. If we get to\
   * the bottom of the depth then we don't have the page in the cache and\
   * we will need to read it in from the disk file.\
   ****************************************************************************/\
   else{\
      for(d=1;d<(image->depth);d++){\
         cstruct2 = image->cachestructure[rseg][cseg][d];\
         image->cachestructure[rseg][cseg][d] = cstruct;\
         cstruct = cstruct2;\
\
         if((cstruct.cpage == cpage) && (cstruct.rpage == rpage)){\
            if(cstruct.dirtybit != INVALID_DATA){\
               image->cachestructure[rseg][cseg][0] = cstruct;   /* put this page on top */\
               (((DTYPE *)(cstruct.dataptr))[(rpix<<(image->cbs_bits))+cpix]) = (DTYPE)pixval;\
               image->cachestructure[rseg][cseg][0].dirtybit = DIRTY_DATA;\
               return;\
            }\
            else d=image->depth;\
         }\
      }\
      image->cachestructure[rseg][cseg][0] = cstruct;   /* put this page on top */\
   }\
\
   /****************************************************************************\
   * We now know that the requested segment was not in the cache. Now we need to\
   * read in the segment from disk. If this is a read/write cache then we must\
   * first write out the page that was in the cache if it is dirty.\
   ****************************************************************************/\
   if((cstruct.dirtybit == DIRTY_DATA) && (image->filemode == CACHE_READ_WRITE))\
      write_cache_segment(image, cstruct.rpage, cstruct.cpage, rseg, cseg, 0);\
\
   /****************************************************************************\
   * Now we must read in the data that we want from disk.\
   ****************************************************************************/\
   (*(image->get_cache_segment))(image, rpage, cpage, rseg, cseg);\
   cstruct = image->cachestructure[rseg][cseg][0];\
\
   /****************************************************************************\
   * Now we set the value of the cache pixel and mark the cache as dirty.\
   ****************************************************************************/\
   (((DTYPE *)(cstruct.dataptr))[(rpix<<(image->cbs_bits))+cpix]) = (DTYPE)pixval;\
   image->cachestructure[rseg][cseg][0].dirtybit = DIRTY_DATA;\
}

/*******************************************************************************
* Procedure: flush_cached_image
* Purpose: Because the image is cached using a write back cache, the image on
* disk may be inconsistent with the image data in the cache. This procedure
* goes through the entire cache and writes any dirty segments to the disk file.
* NAME: Michael Heath, University of South Florida
* Date: 10/14/97
*******************************************************************************/
void flush_cached_image(CACHEIM *image)
{
   CACHESTRUCTURE cstruct;
   int r, c, d;

   /****************************************************************************
   * If the cache mode is not CACHE_READ_WRITE, it makes no sense to flush it.
   * Just set every segment to clean.
   ****************************************************************************/
   if(image->filemode != CACHE_READ_WRITE){
      for(r=0;r<image->nrb;r++){
         for(c=0;c<image->ncb;c++){
            for(d=0;d<image->depth;d++){
               image->cachestructure[r][c][d].dirtybit = CLEAN_DATA;
            }
         }
      }
      return;
   }

   /****************************************************************************
   * If the file pointer is NULL we assume that the image is stored completely
   * in memory so we can not flush it to a file.
   ****************************************************************************/
   if(image->fp == NULL) return;

   for(r=0;r<image->nrb;r++){
      for(c=0;c<image->ncb;c++){
         for(d=0;d<image->depth;d++){
            cstruct = image->cachestructure[r][c][d];
            if(cstruct.dirtybit == DIRTY_DATA)
               write_cache_segment(image, cstruct.rpage, cstruct.cpage, r, c, d);
         }
      }
   }
}

/*******************************************************************************
* Procedure: allocate_cached_image
* Purpose: To allocate an image that can be cached.
* NAME: Michael Heath, University of South Florida
* Date: 10/13/97
*******************************************************************************/
int allocate_cached_image(CACHEIM *image, int rows, int cols, int datatype,
    char *filename, char *mode, int max_depth, int max_cache_size, int headerbytes,
    USHORT swap_bytes, float resolution)
{
   int r, c, d;
   char tmpfilename[100];
   void *imline=NULL;
   int cache_height=0, cache_width=0, cache_depth=0;
   UCHAR *tmp_ptr=NULL;
   int rbs=0, cbs=0, nrb=0, ncb=0;
   UCHAR *data_ptr=NULL;
   UCHAR thebyte;

   /****************************************************************************
   * Initialize all of the fields of the cache image datastructure to safe values.
   ****************************************************************************/
   image->fp = NULL;
   image->filename = NULL;
   image->rows = 0;
   image->cols = 0;
   image->rbs = 0;
   image->cbs = 0;
   image->nrb = 0;
   image->ncb = 0;
   image->rbs_bits = 0;
   image->cbs_bits = 0;
   image->nrb_bits = 0;
   image->ncb_bits = 0;
   image->rbs_and_bits = 0;
   image->cbs_and_bits = 0;
   image->nrb_and_bits = 0;
   image->ncb_and_bits = 0;
   image->swap_bytes = FALSE;
   image->depth = 0;
   image->resolution = 0.0;
   image->datatype = -1;
   image->cachemode = CACHE_READ_ONLY;
   image->filemode = CACHE_READ_ONLY;
   image->cachedata = NULL;
   image->cachestructure = NULL;
   image->headerbytes = 0;
   image->getpixel = NULL;
   image->putpixel = NULL;
   image->get_cache_segment = NULL;
   image->sourceimlist = NULL;
   image->calculationdata = NULL;

   /****************************************************************************
   * It is not worth it to cache very small images. Therefore, if an image is
   * smaller than 6MB we will read it fully into memory.
   ****************************************************************************/
   /* if(max_cache_size < (6*1024*1024)) max_cache_size = 6*1024*1024; */

   /****************************************************************************
   * If the user said that the cache could be as large as the image, then
   * make it the entire size of the image.
   ****************************************************************************/
   if(max_cache_size >= (rows * cols * dt_size[datatype])){
      rbs = (USHORT)rows;
      cbs = (USHORT)cols;
      nrb = 1;
      ncb = 1;
      cache_depth = 1;
   }

   /****************************************************************************
   * Calculate parameters regarding the cache size and how it is organized.
   * Note that cach_height and cache_width are both positive integer powers
   * of two.
   ****************************************************************************/
   else{
      find_cache_parameters(rows, cols, max_depth, max_cache_size, &cache_height,
         &cache_width, &cache_depth, (int)dt_size[datatype]);

      if(cache_height < 16) rbs = 2;
      else rbs = 4;
      nrb = cache_height / rbs;

      if(cache_width < 16) ncb = 2;
      else ncb = 4;
      cbs = cache_width / ncb;

   }

   if(VERBOSEMODE)
      printf("The cache will be: cache_depth=%d rbs=%d cbs=%d nrb=%d ncb=%d\n",
         cache_depth, rbs, cbs, nrb, ncb);

   /****************************************************************************
   * We know that rbs, cbs, nrb and ncb are all integer powers of two. We will
   * compute that power here so we can do bit shifting for multiplication
   * and division. Also we compute variables that can be used for very quick
   * modulus calculations because a%b == a&(b-1) when b is an integer power
   * of two.
   ****************************************************************************/
   image->rbs_bits = ceil(log((double)rbs)/log(2.0));
   image->cbs_bits = ceil(log((double)cbs)/log(2.0));
   image->nrb_bits = ceil(log((double)nrb)/log(2.0));
   image->ncb_bits = ceil(log((double)ncb)/log(2.0));

   image->rbs_and_bits = rbs - 1;
   image->cbs_and_bits = cbs - 1;
   image->nrb_and_bits = nrb - 1;
   image->ncb_and_bits = ncb - 1;

   image->swap_bytes = swap_bytes;

   if(VERBOSEMODE)
      printf("rbs_bits = %d cbs_bits = %d nrb_bits = %d ncb_bits = %d\n", 
         image->rbs_bits, image->cbs_bits, image->nrb_bits, image->ncb_bits);
  

   /****************************************************************************
   * Allocate the memory for the cache itself.
   ****************************************************************************/
   if((image->cachedata = (void *) mycalloc("image->cachedata", nrb*ncb*rbs*cbs*cache_depth,dt_size[datatype]))==NULL){
      fprintf(stderr, "Calloc error in allocate_cached_image().\n");
      return(0);
   }

   /****************************************************************************
   * Allocate memory for the CACHESTRUCTURE.
   ****************************************************************************/
   if((image->cachestructure = (CACHESTRUCTURE ***)mymalloc("image->cachestructure", nrb*sizeof(CACHESTRUCTURE **))) == NULL){
      fprintf(stderr, "Calloc error in allocate_cached_image().\n");
      return(0);
   }

   if((image->cachestructure[0] = (CACHESTRUCTURE **)mycalloc("image->cachestructure[0]",
      nrb*ncb,sizeof(CACHESTRUCTURE *))) == NULL){
      fprintf(stderr, "Calloc error in allocate_cached_image().\n");
      return(0);
   }
   for(r=1;r<nrb;r++) image->cachestructure[r] = image->cachestructure[0] + r * ncb;

   if((image->cachestructure[0][0] = (CACHESTRUCTURE *)mycalloc("image->cachestructure[0][0]", nrb*ncb*cache_depth,
      sizeof(CACHESTRUCTURE))) == NULL){
      fprintf(stderr, "Calloc error in allocate_cached_image().\n");
      return(0);
   }

   for(r=0;r<nrb;r++){
      for(c=0;c<ncb;c++){
         image->cachestructure[r][c] = image->cachestructure[0][0] + (r * ncb + c) * cache_depth;
      }
   }

   /****************************************************************************
   * Set up pointers for each cache segment to the beginning of the segment
   * of memory in the cache. These pointers will eliminate many calculations
   * of addresses in memory accesses.
   ****************************************************************************/
   tmp_ptr = ((UCHAR *)image->cachedata);
   for(r=0;r<nrb;r++){
      for(c=0;c<ncb;c++){
         for(d=0;d<cache_depth;d++){
            image->cachestructure[r][c][d].dataptr = tmp_ptr;
	    tmp_ptr += (rbs * cbs * dt_size[datatype]);
            image->cachestructure[r][c][d].dirtybit = INVALID_DATA;
         }
      }
   }
   tmp_ptr -= (nrb*ncb*rbs*cbs*cache_depth*dt_size[datatype]);

   /****************************************************************************
   * Now that we have successfully allocated memeory, assign values to some of
   * the values in the structure that do not depend on the type of the cache.
   ****************************************************************************/
   image->rows = rows;
   image->cols = cols;
   image->rbs = (USHORT)rbs;
   image->cbs = (USHORT)cbs;
   image->nrb = (USHORT)nrb;
   image->ncb = (USHORT)ncb;
   image->depth = (USHORT)cache_depth;
   image->datatype = datatype;
   image->headerbytes = headerbytes;
   image->resolution = resolution;

   /****************************************************************************
   * Set up pointers to the functions that get and put pixels into the cache.
   * This way we can call functions that work on any datatype without if/else
   * or switch statements.
   ****************************************************************************/
   if(((image->rbs) == rows) && ((image->cbs) == cols)){
      switch(datatype){
         case UCHARNUM:
            image->getpixel = get_fullyinmemory_cache_pixel_uchar;
            image->putpixel = put_fullyinmemory_cache_pixel_uchar;
            break;
         case CHARNUM:
            image->getpixel = get_fullyinmemory_cache_pixel_char;
            image->putpixel = put_fullyinmemory_cache_pixel_char;
            break;
         case USHORTNUM:
            image->getpixel = get_fullyinmemory_cache_pixel_ushort;
            image->putpixel = put_fullyinmemory_cache_pixel_ushort;
            break;
         case SHORTNUM:
            image->getpixel = get_fullyinmemory_cache_pixel_short;
            image->putpixel = put_fullyinmemory_cache_pixel_short;
            break;
         case ULONGNUM:
            image->getpixel = get_fullyinmemory_cache_pixel_ulong;
            image->putpixel = put_fullyinmemory_cache_pixel_ulong;
            break;
         case LONGNUM:
            image->getpixel = get_fullyinmemory_cache_pixel_long;
            image->putpixel = put_fullyinmemory_cache_pixel_long;
            break;
         case FLOATNUM:
            image->getpixel = get_fullyinmemory_cache_pixel_float;
            image->putpixel = put_fullyinmemory_cache_pixel_float;
            break;
         case DOUBLENUM:
            image->getpixel = get_fullyinmemory_cache_pixel_double;
            image->putpixel = put_fullyinmemory_cache_pixel_double;
            break;
      }
   }
   else{
      switch(datatype){
         case UCHARNUM:
            image->getpixel = get_cache_pixel_uchar;
            image->putpixel = put_cache_pixel_uchar;
            break;
         case CHARNUM:
            image->getpixel = get_cache_pixel_char;
            image->putpixel = put_cache_pixel_char;
            break;
         case USHORTNUM:
            image->getpixel = get_cache_pixel_ushort;
            image->putpixel = put_cache_pixel_ushort;
            break;
         case SHORTNUM:
            image->getpixel = get_cache_pixel_short;
            image->putpixel = put_cache_pixel_short;
            break;
         case ULONGNUM:
            image->getpixel = get_cache_pixel_ulong;
            image->putpixel = put_cache_pixel_ulong;
            break;
         case LONGNUM:
            image->getpixel = get_cache_pixel_long;
            image->putpixel = put_cache_pixel_long;
            break;
         case FLOATNUM:
            image->getpixel = get_cache_pixel_float;
            image->putpixel = put_cache_pixel_float;
            break;
         case DOUBLENUM:
            image->getpixel = get_cache_pixel_double;
            image->putpixel = put_cache_pixel_double;
            break;
      }
   }

   /****************************************************************************
   * Set a pointer to the function that is to be used to get a segment of the
   * image from the disk. Other procedures might replace this pointer with a
   * function to get processed data by actually doing the processing.
   ****************************************************************************/
   image->get_cache_segment = read_cache_segment;

   /****************************************************************************
   * If the image is to be created and it fits into memory then just set the
   * fields of the datastructure and return. We will not be using a file in this
   * case.
   ****************************************************************************/
   if((strcmp(mode,"w+b") == 0) && ((image->rbs) == rows) && ((image->cbs) == cols)){

      if(VERBOSEMODE)
         printf("Putting the entire image in memory.\n");

      image->fp = NULL;
      image->filename = NULL;
      image->cachemode = CACHE_READ_WRITE;
      image->filemode = CACHE_READ_ONLY;
      for(r=0;r<nrb;r++){
         for(c=0;c<ncb;c++){
            for(d=0;d<cache_depth;d++){
               image->cachestructure[r][c][d].dirtybit = CLEAN_DATA;
            }
         }
      }
      return(1);
   }

   /****************************************************************************
   * In this case we want to calculate the values as we need them. That is why
   * the "filename" has the value "tobecomputed".
   ****************************************************************************/
   if((strcmp(mode,"rb") == 0) && (strcmp(filename,"tobecomputed")==0)){

      if(VERBOSEMODE)
         printf("The image is to be computed as needed.\n");

      image->fp = NULL;
      image->filename = NULL;
      image->cachemode = CACHE_READ_WRITE;
      image->filemode = CACHE_READ_ONLY;

      return(1);
   }

   /****************************************************************************
   * If the image is to be opened rb or r+b it must already exist. Therefore
   * we open the file using the filename we were given.
   *     rb  - open a binary file for reading
   *     r+b - open a binary file for read/write
   ****************************************************************************/
   if((strcmp(mode,"rb") == 0) || (strcmp(mode,"r+b") == 0)){

      if(filename == NULL){
         fprintf(stderr, "No filename was supplied yet trying to open as r+b in allocate_cached_image().\n");
         return(0);
      }

      if((image->filename = (CHAR *) mycalloc("image->filename", strlen(filename)+1, sizeof(CHAR)))==NULL){
         fprintf(stderr, "Malloc error in allocate_cached_image().\n");
         return(0);
      }
      sprintf(image->filename, "%s", filename);

      if(strcmp(mode,"rb") == 0){

	 if(VERBOSEMODE)
            printf("The image is to be read only from a file.\n");

         if((image->fp = fopen(image->filename, mode)) == NULL){
            fprintf(stderr, "Error opening the image %s for reading in allocate_cached_image().\n",
               image->filename);
            return(0);
         }
         image->cachemode = CACHE_READ_ONLY;
         image->filemode = CACHE_READ_ONLY;
      }

      if(strcmp(mode,"r+b") == 0){

	 if(VERBOSEMODE)
            printf("The image is to be read/write from an existing file.\n");

         if((image->fp = fopen(image->filename, mode)) == NULL){
            fprintf(stderr, "Error opening the image %s for reading and writing in allocate_cached_image().\n",
               image->filename);
            return(0);
         }
         image->cachemode = CACHE_READ_WRITE;
         image->filemode = CACHE_READ_WRITE;
      }

      /*************************************************************************
      * If the image is to be fully in memory, then read in the whole image.
      * If the bytes need to be swapped, then go ahead and swap them.
      *************************************************************************/
      if((image->rbs == rows) && (image->cbs == cols)){

         if(((image->fp) = fopen((image->filename), mode)) == NULL){
            fprintf(stderr, "Error opening the file %s.\n", image->filename);
            exit(1);
         }

         fseek(image->fp, image->headerbytes, 0);
         if(fread(image->cachedata, dt_size[image->datatype], image->rows * image->cols,
	    image->fp) != (image->rows * image->cols)){
            fprintf(stderr, "Error reading from the file %s.\n", image->filename);
            exit(1);
         }

         if(image->swap_bytes == TRUE){

            /*************************************************************************
            * Since we must swap the bytes, swap them.
            *************************************************************************/
            data_ptr = (UCHAR *)(image->cachedata);
      
            if(dt_size[image->datatype] == 2){
               for(r=0;r<(image->rbs*image->cbs);r++,data_ptr+=2){
                  thebyte = *(data_ptr);
                  *data_ptr = *(data_ptr+1);
                  *(data_ptr+1) = thebyte;
               }
            }
            else if(((image->datatype == ULONGNUM) || (image->datatype == LONGNUM))){
               for(r=0;r<(image->rbs*image->cbs);r++,data_ptr+=dt_size[image->datatype]){
                  thebyte = *(data_ptr);
                  *data_ptr = *(data_ptr+3);
                  *(data_ptr+3) = thebyte;
      
                  thebyte = *(data_ptr+1);
                  *(data_ptr+1) = *(data_ptr+2);
                  *(data_ptr+2) = thebyte;
               }
            }

         }

         if(strcmp(mode,"rb") == 0){
            fclose(image->fp);
            image->fp = NULL;
         }
      }
      return(1);
   }

   /****************************************************************************
   * Because we have not returned yet, we need to create the file. The calling
   * function may have passed a filename to us, or may want us to create the
   * filename.
   ****************************************************************************/
   if(filename != NULL) sprintf(tmpfilename, "%s", filename);
   /*   else sprintf(tmpfilename, "tmp_%d_%d_%s.%lu", rows, cols, dt_name[datatype], (ULONG)rand()); */
   else sprintf(tmpfilename, "system_random_file");

   if((image->filename = (CHAR *) mycalloc("image->filename", strlen(tmpfilename)+1, sizeof(CHAR))) == NULL){
      fprintf(stderr, "Calloc error in allocate_cached_image().\n");
      return(0);
   }
   sprintf(image->filename, "%s", tmpfilename);

   /****************************************************************************
   * Open the file in the appropriate mode.
   *     w+b - create a binary file for read/write
   ****************************************************************************/
   if(strcmp(mode,"w+b") == 0){

      if(VERBOSEMODE)
         printf("The image is to be read/write with a file to be created.\n");

      if(filename != NULL){
         if((image->fp = fopen(image->filename, mode)) == NULL){
            fprintf(stderr, "Error creating the image %s for reading and writing in allocate_cached_image().\n",
               image->filename);
            return(0);
         }
      }
      else{
         if((image->fp = tmpfile()) == NULL){
            fprintf(stderr, "Error creating a temporary file for reading and writing in allocate_cached_image().\n");
            return(0);
         }
      }

      image->cachemode = CACHE_READ_WRITE;
      image->filemode = CACHE_READ_WRITE;

      /*************************************************************************
      * Write out space for the header bytes.
      *************************************************************************/
      for(r=0;r<image->headerbytes;r++) fputc((char)0, image->fp);

      if((imline = (void *) mycalloc("imline", cols,(size_t)(dt_size[datatype]))) == NULL){
         fprintf(stderr, "Calloc error in allocate_cached_image().\n");
         return(0);
      }

      for(r=0;r<rows;r++){
         if(fwrite(imline, dt_size[datatype], cols, image->fp) != cols){
            fprintf(stderr, "Error initializing the file %s in allocate_cached_image().\n", image->filename);
            return(0);
         }
      }
      fflush(image->fp);
      myfree("imline", imline);
      return(1);
   }

   /****************************************************************************
   * If we had a valid mode, we would have returned already.
   ****************************************************************************/
   fprintf(stderr, "Invalid mode in allocate_cached_image().\n");
   return(0);
}

/*******************************************************************************
* Procedure: deallocate_cached_image
* Purpose: To flush a cached image, close the associated file (if one is being
* used), deallocate all of the memory and reinitialize all of the fileds of the
* datastructure.
* NAME: Michael Heath, University of South Florida
* Date: 10/13/97
*******************************************************************************/
void deallocate_cached_image(CACHEIM *image)
{
   /****************************************************************************
   * Flush the image to make sure the disk file is updated with any changes
   * made to the cache.
   ****************************************************************************/
   flush_cached_image(image);

   /****************************************************************************
   * Close the disk file (if one exists).
   ****************************************************************************/
   if(image->fp != NULL) fclose(image->fp);

   /****************************************************************************
   * Free the memory associated with this image.
   ****************************************************************************/
   myfree("image->cachestructure[0][0]", image->cachestructure[0][0]);
   myfree("image->cachestructure[0]", image->cachestructure[0]);
   myfree("image->cachestructure", image->cachestructure);
   myfree("image->cachedata", image->cachedata);
   if(image->filename != NULL) myfree("image->filename", image->filename);
   if((image->sourceimlist) != NULL) myfree("image->sourceimlist", image->sourceimlist);
   if((image->calculationdata) != NULL) myfree("image->calculationdata", image->calculationdata);

   /****************************************************************************
   * Initialize all of the fields of the cache image datastructure to safe values.
   ****************************************************************************/
   image->fp = NULL;
   image->filename = NULL;
   image->rows = 0;
   image->cols = 0;
   image->rbs = 0;
   image->cbs = 0;
   image->nrb = 0;
   image->ncb = 0;
   image->rbs_bits = 0;
   image->cbs_bits = 0;
   image->nrb_bits = 0;
   image->ncb_bits = 0;
   image->rbs_and_bits = 0;
   image->cbs_and_bits = 0;
   image->nrb_and_bits = 0;
   image->ncb_and_bits = 0;
   image->depth = 0;
   image->resolution = 0.0;
   image->datatype = -1;
   image->cachemode = CACHE_READ_ONLY;
   image->filemode = CACHE_READ_ONLY;
   image->cachedata = NULL;
   image->cachestructure = NULL;
   image->headerbytes = 0;
   image->getpixel = NULL;
   image->putpixel = NULL;
   image->get_cache_segment = NULL;
   image->sourceimlist = NULL;
   image->calculationdata = NULL;
}

/*******************************************************************************
* Procedure: read_cache_segment
* Purpose: This procedure reads in a segment of the image from the file and
* places it in the cache. The cache segment that we are reading into will always
* be on top of the depth dimension.
* NAME: Michael Heath, University of South Florida
* Date: 10/14/97
*******************************************************************************/
void read_cache_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg)
{
   CACHESTRUCTURE cstruct;
   UCHAR *data_ptr=NULL;
   long seekposition=0;
   int r, base;
   UCHAR thebyte;

   if(VERBOSEMODE)
      printf("read_cache_segment()  (rpage=%3d,cpage=%3d,rseg=%3d,cseg=%3d)\n", rpage, cpage, rseg, cseg);

   if(image->cachestructure[rseg][cseg][0].dirtybit == DIRTY_DATA){
      fprintf(stderr, "Trying to read into dirty memory in read_cache_segment().\n");
      exit(1);
   }

   /****************************************************************************
   * Fill in the CACHESTRUCTURE information.
   ****************************************************************************/
   image->cachestructure[rseg][cseg][0].rpage = (USHORT)rpage;
   image->cachestructure[rseg][cseg][0].cpage = (USHORT)cpage;
   image->cachestructure[rseg][cseg][0].dirtybit = CLEAN_DATA;

   /****************************************************************************
   * Calculate the number of rows and the number of columns in this segment.
   * This will be image->rbs and image->cbs for most segments. Segments on the
   * right edge or on the bottom of the image may have less rows or columns in
   * the segment.
   ****************************************************************************/
   base = (rpage * image->nrb + rseg) * image->rbs;
   if((image->rows - base) >= image->rbs) image->cachestructure[rseg][cseg][0].rowlim = image->rbs;
   else image->cachestructure[rseg][cseg][0].rowlim = image->rows - base;

   base = (cpage * image->ncb + cseg) * image->cbs;
   if((image->cols - base) >= image->cbs) image->cachestructure[rseg][cseg][0].collim = image->cbs;
   else image->cachestructure[rseg][cseg][0].collim = image->cols - base;

   cstruct = image->cachestructure[rseg][cseg][0];

   /****************************************************************************
   * Read the segment from the file when the disk has valid data on it.
   * Otherwise just set the bytes to zero.
   ****************************************************************************/
   data_ptr = (UCHAR *)cstruct.dataptr;

   if(VERBOSEMODE)
      printf("      reading with fread (rpage=%d,cpage=%d,rseg=%d,cseg=%d)\n", rpage, cpage, rseg, cseg);

   seekposition = (((rpage * image->nrb + rseg) * image->rbs) * image->cols +
                   ((cpage * image->ncb + cseg) * image->cbs)) * dt_size[image->datatype] + image->headerbytes;

   for(r=0;r<cstruct.rowlim;r++){

/*
      printf("seekposition = %ld ", seekposition);
      fflush(stdout);

      printf("fileptr = %p ", image->fp);
      fflush(stdout);
*/

      if(fseek(image->fp, seekposition, 0) != 0){
         fprintf(stderr, "Fseek error in read_cache_segment().\n");
         exit(1);
      }
/*
      printf("%lu\n", ftell(image->fp));
      fflush(stdout);
*/

      if(fread((void *)data_ptr, (size_t)dt_size[image->datatype], cstruct.collim, image->fp) != cstruct.collim){
         fprintf(stderr, "Fread error in read_cache_segment().\n");
         exit(1);
      }
      data_ptr += image->cbs * dt_size[image->datatype];
      seekposition += image->cols * dt_size[image->datatype];
   }
   fflush(image->fp);

   /*************************************************************************
   * If we must swap the bytes, then swap them.
   *************************************************************************/
   if(image->swap_bytes == TRUE){

      data_ptr = (UCHAR *)cstruct.dataptr;

      if(dt_size[image->datatype] == 2){
         for(r=0;r<(image->rbs*image->cbs);r++,data_ptr+=2){
            thebyte = *(data_ptr);
            *data_ptr = *(data_ptr+1);
            *(data_ptr+1) = thebyte;
         }
      }
      else if(((image->datatype == ULONGNUM) || (image->datatype == LONGNUM))){
         for(r=0;r<(image->rbs*image->cbs);r++,data_ptr+=dt_size[image->datatype]){
            thebyte = *(data_ptr);
            *data_ptr = *(data_ptr+3);
            *(data_ptr+3) = thebyte;

            thebyte = *(data_ptr+1);
            *(data_ptr+1) = *(data_ptr+2);
            *(data_ptr+2) = thebyte;
         }
      }
   }
}

/*******************************************************************************
* Procedure: write_cache_segment
* Purpose: This procedure writes a segment of the image from the cache to a
* file. If the cache mode is CACHE_READ_ONLY an error occurs. If the dirtybit
* is INVALID_DATA or CLEAN_DATA, the data is not written to the file. The
* cache segment that we are writing could be anywhere in the depth dimension.
* Thus the parameter depth tells us where in depth this cache page is located.
* NAME: Michael Heath, University of South Florida
* Date: 10/14/97
*******************************************************************************/
void write_cache_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg, int depth)
{
   CACHESTRUCTURE cstruct;
   UCHAR *data_ptr=NULL;
   long seekposition=0;
   int r, base;
   UCHAR thebyte;

   /****************************************************************************
   * Make sure that we are allowed to write to the file.
   ****************************************************************************/
   if(image->filemode == CACHE_READ_ONLY){
      image->cachestructure[rseg][cseg][depth].dirtybit = CLEAN_DATA;
      return;
   }

   cstruct = image->cachestructure[rseg][cseg][depth];

   /****************************************************************************
   * If the cache segment is not dirty, don't write it out. Just return.
   ****************************************************************************/
   if(cstruct.dirtybit != DIRTY_DATA) return;

/*
   if(VERBOSEMODE)
      printf("write_cache_segment() (rpage=%3d,cpage=%3d,rseg=%3d,cseg=%3d)\n", rpage, cpage, rseg, cseg);
*/

   /****************************************************************************
   * Calculate the number of rows and the number of columns in this segment.
   * This will be image->rbs and image->cbs for most segments. Segments on the
   * right edge or on the bottom of the image may have less rows or columns in
   * the segment.
   ****************************************************************************/
   base = (rpage * image->nrb + rseg) * image->rbs;
   if((image->rows - base) >= image->rbs) image->cachestructure[rseg][cseg][depth].rowlim = image->rbs;
   else image->cachestructure[rseg][cseg][depth].rowlim = image->rows - base;

   base = (cpage * image->ncb + cseg) * image->cbs;
   if((image->cols - base) >= image->cbs) image->cachestructure[rseg][cseg][depth].collim = image->cbs;
   else image->cachestructure[rseg][cseg][depth].collim = image->cols - base;

   /****************************************************************************
   * Since we just changed something (rowlim and collim) in the cachestructure,
   * get a new compy of it to work with (cstruct).
   ****************************************************************************/
   cstruct = image->cachestructure[rseg][cseg][depth];

   /*************************************************************************
   * If we must swap the bytes, then swap them.
   *************************************************************************/
   if(image->swap_bytes == TRUE){

      data_ptr = (UCHAR *)cstruct.dataptr;

      if(dt_size[image->datatype] == 2){
         for(r=0;r<(image->rbs*image->cbs);r++,data_ptr+=2){
            thebyte = *(data_ptr);
            *data_ptr = *(data_ptr+1);
            *(data_ptr+1) = thebyte;
         }
      }
      else if(((image->datatype == ULONGNUM) || (image->datatype == LONGNUM))){
         for(r=0;r<(image->rbs*image->cbs);r++,data_ptr+=dt_size[image->datatype]){
            thebyte = *(data_ptr);
            *data_ptr = *(data_ptr+3);
            *(data_ptr+3) = thebyte;

            thebyte = *(data_ptr+1);
            *(data_ptr+1) = *(data_ptr+2);
            *(data_ptr+2) = thebyte;
         }
      }
   }

   /****************************************************************************
   * Write the data out to the file. Then set the dirty bit to CLEAN_DATA.
   ****************************************************************************/
   data_ptr = (UCHAR *)cstruct.dataptr;

   seekposition = (((rpage * image->nrb + rseg) * image->rbs) * image->cols +
                   ((cpage * image->ncb + cseg) * image->cbs)) * dt_size[image->datatype] + image->headerbytes;

   for(r=0;r<cstruct.rowlim;r++){
      if(fseek(image->fp, seekposition, 0) != 0){
         fprintf(stderr, "Fseek error in write_cache_segment().\n");
         exit(1);
      }
      if(fwrite((void *)data_ptr, (size_t)dt_size[image->datatype], cstruct.collim, image->fp) != cstruct.collim){
         fprintf(stderr, "Fwrite error in write_cache_segment().\n");
         exit(1);
      }
      data_ptr += image->cbs * dt_size[image->datatype];
      seekposition += image->cols * dt_size[image->datatype];
   }
   fflush(image->fp);
   image->cachestructure[rseg][cseg][depth].dirtybit = CLEAN_DATA;

   /*************************************************************************
   * If swapped the bytes, we must then swap them back. This is the price
   * that we pay for not copying the data into other memory before swapping
   * the bytes in th first place.
   *************************************************************************/
   if(image->swap_bytes == TRUE){

      data_ptr = (UCHAR *)cstruct.dataptr;

      if(dt_size[image->datatype] == 2){
         for(r=0;r<(image->rbs*image->cbs);r++,data_ptr+=2){
            thebyte = *(data_ptr);
            *data_ptr = *(data_ptr+1);
            *(data_ptr+1) = thebyte;
         }
      }
      else if(((image->datatype == ULONGNUM) || (image->datatype == LONGNUM))){
         for(r=0;r<(image->rbs*image->cbs);r++,data_ptr+=dt_size[image->datatype]){
            thebyte = *(data_ptr);
            *data_ptr = *(data_ptr+3);
            *(data_ptr+3) = thebyte;

            thebyte = *(data_ptr+1);
            *(data_ptr+1) = *(data_ptr+2);
            *(data_ptr+2) = thebyte;
         }
      }
   }
}

/*******************************************************************************
* Procedure: get_fullyinmemory_cache_pixel_uchar
* Procedure: get_fullyinmemory_cache_pixel_char
* Procedure: get_fullyinmemory_cache_pixel_ushort
* Procedure: get_fullyinmemory_cache_pixel_short
* Procedure: get_fullyinmemory_cache_pixel_ulong
* Procedure: get_fullyinmemory_cache_pixel_long
* Procedure: get_fullyinmemory_cache_pixel_float
* Procedure: get_fullyinmemory_cache_pixel_double
* Purpose: These procedures retrieve the value of a pixel in the cache image.
* It is only used for images that are stored completely in the cache.
* The actual work of the function is written in tha macro
* "getpixel_fullyinmemory_macro".
* NAME: Michael Heath, University of South Florida
* Date: 10/28/97
*******************************************************************************/
double get_fullyinmemory_cache_pixel_uchar(CACHEIM *image, int row, int col)
{
   getpixel_fullyinmemory_macro(UCHAR)
}
double get_fullyinmemory_cache_pixel_char(CACHEIM *image, int row, int col)
{
   getpixel_fullyinmemory_macro(CHAR)
}
double get_fullyinmemory_cache_pixel_ushort(CACHEIM *image, int row, int col)
{
   getpixel_fullyinmemory_macro(USHORT)
}
double get_fullyinmemory_cache_pixel_short(CACHEIM *image, int row, int col)
{
   getpixel_fullyinmemory_macro(SHORT)
}
double get_fullyinmemory_cache_pixel_ulong(CACHEIM *image, int row, int col)
{
   getpixel_fullyinmemory_macro(ULONG)
}
double get_fullyinmemory_cache_pixel_long(CACHEIM *image, int row, int col)
{
   getpixel_fullyinmemory_macro(LONG)
}
double get_fullyinmemory_cache_pixel_float(CACHEIM *image, int row, int col)
{
   getpixel_fullyinmemory_macro(FLOAT)
}
double get_fullyinmemory_cache_pixel_double(CACHEIM *image, int row, int col)
{
   getpixel_fullyinmemory_macro(DOUBLE)
}

/*******************************************************************************
* Procedure: put_fullyinmemory_cache_pixel_uchar
* Procedure: put_fullyinmemory_cache_pixel_char
* Procedure: put_fullyinmemory_cache_pixel_ushort
* Procedure: put_fullyinmemory_cache_pixel_short
* Procedure: put_fullyinmemory_cache_pixel_ulong
* Procedure: put_fullyinmemory_cache_pixel_long
* Procedure: put_fullyinmemory_cache_pixel_float
* Procedure: put_fullyinmemory_cache_pixel_double
* Purpose: This procedure sets the value of a pixel in the cache image.
* It is only used for images that are stored completely in the cache.
* The actual work of the function is written in tha macro
* "putpixel_fullyinmemory_macro".
* NAME: Michael Heath, University of South Florida
* Date: 10/14/97
*******************************************************************************/
void put_fullyinmemory_cache_pixel_uchar(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_fullyinmemory_macro(UCHAR)
}
void put_fullyinmemory_cache_pixel_char(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_fullyinmemory_macro(CHAR)
}
void put_fullyinmemory_cache_pixel_ushort(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_fullyinmemory_macro(USHORT)
}
void put_fullyinmemory_cache_pixel_short(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_fullyinmemory_macro(SHORT)
}
void put_fullyinmemory_cache_pixel_ulong(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_fullyinmemory_macro(ULONG)
}
void put_fullyinmemory_cache_pixel_long(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_fullyinmemory_macro(LONG)
}
void put_fullyinmemory_cache_pixel_float(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_fullyinmemory_macro(FLOAT)
}
void put_fullyinmemory_cache_pixel_double(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_fullyinmemory_macro(DOUBLE)
}

/*******************************************************************************
* Procedure: get_cache_pixel_uchar
* Procedure: get_cache_pixel_char
* Procedure: get_cache_pixel_ushort
* Procedure: get_cache_pixel_short
* Procedure: get_cache_pixel_ulong
* Procedure: get_cache_pixel_long
* Procedure: get_cache_pixel_float
* Procedure: get_cache_pixel_double
* Purpose: These procedures retrieve the value of a pixel in the cache image. If
* necessary, the a dirty segment is written to disk first and the appropriate
* segment is read in from disk. The actual work of the function is written in
* tha macro "getpixel_macro".
* NAME: Michael Heath, University of South Florida
* Date: 10/14/97
*******************************************************************************/
double get_cache_pixel_uchar(CACHEIM *image, int row, int col)
{
   getpixel_macro(UCHAR)
}
double get_cache_pixel_char(CACHEIM *image, int row, int col)
{
   getpixel_macro(CHAR)
}
double get_cache_pixel_ushort(CACHEIM *image, int row, int col)
{
   getpixel_macro(USHORT)
}
double get_cache_pixel_short(CACHEIM *image, int row, int col)
{
   getpixel_macro(SHORT)
}
double get_cache_pixel_ulong(CACHEIM *image, int row, int col)
{
   getpixel_macro(ULONG)
}
double get_cache_pixel_long(CACHEIM *image, int row, int col)
{
   getpixel_macro(LONG)
}
double get_cache_pixel_float(CACHEIM *image, int row, int col)
{
   getpixel_macro(FLOAT)
}
double get_cache_pixel_double(CACHEIM *image, int row, int col)
{
   getpixel_macro(DOUBLE)
}

/*******************************************************************************
* Procedure: put_cache_pixel_uchar
* Procedure: put_cache_pixel_char
* Procedure: put_cache_pixel_ushort
* Procedure: put_cache_pixel_short
* Procedure: put_cache_pixel_ulong
* Procedure: put_cache_pixel_long
* Procedure: put_cache_pixel_float
* Procedure: put_cache_pixel_double
* Purpose: This procedure sets the value of a pixel in the cache image. If
* necessary, the a dirty segment is written to disk first and the appropriate
* segment is read in from disk. The actual work of the function is written in
* tha macro "putpixel_macro".
* NAME: Michael Heath, University of South Florida
* Date: 10/14/97
*******************************************************************************/
void put_cache_pixel_uchar(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_macro(UCHAR)
}
void put_cache_pixel_char(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_macro(CHAR)
}
void put_cache_pixel_ushort(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_macro(USHORT)
}
void put_cache_pixel_short(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_macro(SHORT)
}
void put_cache_pixel_ulong(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_macro(ULONG)
}
void put_cache_pixel_long(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_macro(LONG)
}
void put_cache_pixel_float(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_macro(FLOAT)
}
void put_cache_pixel_double(CACHEIM *image, int row, int col, double pixval)
{
   putpixel_macro(DOUBLE)
}

/*******************************************************************************
* Procedure: find_cache_parameters
* Purpose: To calculate the parameters for the cache.
* NAME: Michael Heath, University of South Florida
* Date: 10/22/97
*******************************************************************************/
void find_cache_parameters(int rows, int cols, int max_depth, int max_size, int *a, int *b, int *depth, int bytes_per_pixel)
{
   int i=0, rsize=0, csize=0, cache_size=0;

/*   printf("\n"); */
   do{
      i += 1;
      rsize = get_csize(rows/4, i);
      csize = get_csize(cols, i);

      cache_size = i * rsize * csize * bytes_per_pixel;
/*
      printf("%2d (%10d %10d) (%10d %10d) %10d\n", i,
         rsize, csize, i*rsize, i*csize, cache_size);
*/
   }while((i < max_depth) && (cache_size > max_size));  
/*
   printf("%2d (%10d %10d) (%10d %10d) %10d\n", i,
      rsize, csize, i*rsize, i*csize, cache_size);
*/
   *a = rsize;
   *b = csize;
   *depth = i;
}

/*******************************************************************************
* Function: get_csize
* Purpose: To calculate 2^(ceil(log_base2(x/depth))).
* NAME: Michael Heath, University of South Florida
* Date: 10/22/97
*******************************************************************************/
int get_csize(int x, int depth)
{
   double tmp;
   int tmp_int;

   tmp_int = ceil( log( (double)x / (double)depth ) / log(2.0) );
   tmp = pow(2.0, (double)tmp_int);
   return((int)tmp);
}

void print_virtual_image_information(CACHEIM image, char *varname)
{
   printf("**********************************************************************\n");
   printf("********** Image %s ***************\n", varname);
   printf("fp = %p\n", (void *)image.fp);
   if(image.filename == NULL) printf("filename = NULL\n");
   else printf("filename = %s\n", image.filename);
   printf("rbs = %hu  cbs = %hu  nrb = %hu  ncb = %hu\n", image.rbs, image.cbs, image.nrb, image.ncb);
   printf("rbs_bits = %hu  cbs_bits = %hu  nrb_bits = %hu  ncb_bits = %hu\n",
      image.rbs_bits, image.cbs_bits, image.nrb_bits, image.ncb_bits);
   printf("rbs_and_bits = %hu  cbs_and_bits = %hu  nrb_and_bits = %hu  ncb_and_bits = %hu\n",
      image.rbs_and_bits, image.cbs_and_bits, image.nrb_and_bits, image.ncb_and_bits);
   printf("depth = %hu\n", image.depth);
   printf("datatype = %s\n", dt_name[image.datatype]);
   printf("resolution = %f\n", image.resolution);
   if(image.swap_bytes == TRUE) printf("swap_bytes = TRUE\n");
   else if(image.swap_bytes == FALSE) printf("swap_bytes = FALSE\n");
   printf("headerbytes = %d\n", image.headerbytes);
   printf("rows = %d   cols = %d\n", image.rows, image.cols);
   if(image.cachemode == CACHE_READ_ONLY) printf("cachemode = CACHE_READ_ONLY\n");
   else if(image.cachemode == CACHE_READ_WRITE) printf("cachemode = CACHE_READ_WRITE\n");
   if(image.filemode == CACHE_READ_ONLY) printf("filemode = CACHE_READ_ONLY\n");
   else if(image.filemode == CACHE_READ_WRITE) printf("filemode = CACHE_READ_WRITE\n");
   printf("get_cache_segment = %p\n", (void *)image.get_cache_segment);

   if(image.sourceimlist != NULL){
      printf("\n");
      /* print_virtual_image_information(*(((CACHEIM **)(image.sourceimlist))[0]), "SOURCE IMAGE"); */
   }
   printf("**********************************************************************\n");
}

