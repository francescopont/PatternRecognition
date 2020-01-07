/*******************************************************************************
* PROGRAM: convolve.c
* PURPOSE: This file contains the code for doing convolutions with virtual
* images.
* NAME: Michael Heath, University of South Florida
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include "virtual_image.h"
#include <string.h>
#include "myalloc.h"

#define XDIR 1
#define YDIR 2

typedef struct{
   int kernel_rows, kernel_cols;
   double *kernel_coefficients;
}CONVOLVE_DATA;

typedef struct{
   int kernel_size;
   int direction;
   double *kernel_coefficients;
}ONE_D_CONVOLVE_DATA;

void get_nonsep_convolved_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg);
void nonsep_convolve(CACHEIM *image, CACHEIM *source_image, int krows, int kcols,
    double *kernel, int max_depth, int max_cache_size, int datatype, char *mode);

void get_sep_convolved_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg);
void one_d_convolve(CACHEIM *image, CACHEIM *source_image, int ksize,
    double *kernel, int direction, int max_depth, int max_cache_size,
    int datatype, char *mode);

/*******************************************************************************
* Procedure: nonsep_convolve
* Name: Michael Heath, University of South Florida
* Date: 10/22/97
*******************************************************************************/
void nonsep_convolve(CACHEIM *image, CACHEIM *source_image, int krows, int kcols,
    double *kernel, int max_depth, int max_cache_size, int datatype, char *mode)
{
   CACHEIM *newim=NULL;
   CONVOLVE_DATA *convolve_data;
   int i, r, c, precomputeimage, rows, cols;

   /****************************************************************************
   * If the mode is "r+b" show an error because we can not already have the
   * result of the convolution in a file.
   ****************************************************************************/
   if(strcmp(mode, "r+b") == 0){
      fprintf(stderr, "Error. Tried to open a convolve image in mode \"r+b\".\n");
      exit(1);
   }

   rows = source_image->rows;
   cols = source_image->cols;

   if((rows*cols*dt_size[datatype]) <= max_cache_size) precomputeimage = 1;
   else precomputeimage = 0;

   /****************************************************************************
   * Allocate a virtual image.
   ****************************************************************************/
   if(precomputeimage){
      if(allocate_cached_image(image, rows, cols,
         datatype, "tobecomputed", "rb", max_depth, (rows*cols*dt_size[datatype])/4, 0, FALSE,
         source_image->resolution) == 0) exit(1);
   }
   else{
      if(allocate_cached_image(image, rows, cols,
         datatype, "tobecomputed", "rb", max_depth, max_cache_size, 0, FALSE,
         source_image->resolution) == 0) exit(1);
   }

   /****************************************************************************
   * Set up the data necessary for computing the convolution.
   ****************************************************************************/
   image->calculationdata = (void *) mymalloc("image->calculationdata", sizeof(CONVOLVE_DATA));

   convolve_data = (CONVOLVE_DATA *) image->calculationdata;
   convolve_data->kernel_rows = krows;
   convolve_data->kernel_cols = kcols;
   convolve_data->kernel_coefficients = (double *) calloc(krows*kcols, sizeof(double));
   for(i=0;i<(krows*kcols);i++) (convolve_data->kernel_coefficients)[i] = kernel[i];

   image->sourceimlist = (void **) mycalloc("image->sourceimlist", 1, sizeof(CACHEIM *));
   ((CACHEIM **)(image->sourceimlist))[0] = source_image;

   image->get_cache_segment = get_nonsep_convolved_segment;

   /****************************************************************************
   * If the mode is rb then we are finished. We have told the image how to
   * compute its values.
   ****************************************************************************/
   if((strcmp(mode, "rb") == 0) && (precomputeimage == 0)) return;

   /* printf("Precalculating the entire nonsep_convolved image.\n"); */

   /****************************************************************************
   * If the mode is "w+b" then we want to immediately calculate the result for
   * the entire image and cache the result in a file.
   ****************************************************************************/
   newim = (CACHEIM *) mycalloc("newim", 1, sizeof(CACHEIM));

   /****************************************************************************
   * Allocate a virtual image.
   ****************************************************************************/
   if(allocate_cached_image(newim, rows, cols,
      datatype, NULL, "w+b", max_depth, rows*cols*dt_size[datatype], 0, 0,
      source_image->resolution) == 0) exit(1);

   /****************************************************************************
   * Compute the result of the convolution for every pixel in the image.
   ****************************************************************************/
   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         (*(newim->putpixel))(newim, r, c, (*(image->getpixel))(image, r, c));
      }
   }

   deallocate_cached_image(image);
   *image = *newim;
   myfree("newim", newim);
}

/*******************************************************************************
* Procedure: get_nonsep_convolved_segment
* Purpose: This procedure calculates the output convolved values for a segment
* of the image. This function is called by get_pixel from virtual images
* of the type convolve.
* Name: Michael Heath, University of South Florida
* Date: 10/22/97
*******************************************************************************/
void get_nonsep_convolved_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg)
{
   int r=0, c=0, rr=0, cc=0, sr=0, sc=0, krows=0, kcols=0, start_row=0, start_col=0, num_rows=0, num_cols=0;
   double sum=0.0;
   USHORT hold_cachemode;
   CONVOLVE_DATA *convolve_data;
   CACHEIM *sourceimage=NULL;
   double *kernel=NULL;

   /****************************************************************************
   * Determine the portion of the image that this segment falls within. This
   * area contains the pixels we will be calculating.
   ****************************************************************************/
   image->cachestructure[rseg][cseg][0].rpage = (USHORT)rpage;
   image->cachestructure[rseg][cseg][0].cpage = (USHORT)cpage;
   image->cachestructure[rseg][cseg][0].dirtybit = CLEAN_DATA;
   /* cstruct = image->cachestructure[rseg][cseg][0]; */

   start_row = (rpage * image->nrb + rseg) * image->rbs;
   if((image->rows - start_row) >= image->rbs) num_rows = image->cachestructure[rseg][cseg][0].rowlim = image->rbs;
   else num_rows = image->cachestructure[rseg][cseg][0].rowlim = image->rows - start_row;

   start_col = (cpage * image->ncb + cseg) * image->cbs;
   if((image->cols - start_col) >= image->cbs) num_cols = image->cachestructure[rseg][cseg][0].collim = image->cbs;
   else num_cols = image->cachestructure[rseg][cseg][0].collim = image->cols - start_col;

   /****************************************************************************
   * We know that the cache segment that we want to update is in the cache.
   * Therefore we will temporarily allow ourselves to write into the cache.
   * It is very important that we replace the cache mode with its proper value
   * when we are finished with this function.
   ****************************************************************************/
   hold_cachemode = image->cachemode;
   image->cachemode = CACHE_READ_WRITE;

   /****************************************************************************
   * Get the data necessary to do the convolution from the image structure.
   ****************************************************************************/
   convolve_data = (CONVOLVE_DATA *)(image->calculationdata);
   krows = convolve_data->kernel_rows;
   kcols = convolve_data->kernel_cols;
   kernel = convolve_data->kernel_coefficients;
   sourceimage = (((CACHEIM **)image->sourceimlist)[0]);

   /*
   printf("get_nonsep_convolved_segment(image, %d, %d, %d, %d) %d %d %d %d\n",
   rpage, cpage, rseg, cseg, start_row, start_col, num_rows, num_cols);
   */

   /****************************************************************************
   * Do the convolution in this region of the image.
   ****************************************************************************/
   for(r=start_row;r<(start_row+num_rows);r++){
      for(c=start_col;c<(start_col+num_cols);c++){
         sum = 0;
         for(rr=0,sr=(r-krows/2);rr<krows;rr++,sr++){
            if((sr >= 0) && (sr < (image->rows))){
               for(cc=0,sc=(c-kcols/2);cc<kcols;cc++,sc++){
                  if((sc >= 0) && (sc < (image->cols))){
                     sum += kernel[rr*kcols+cc] * (*(sourceimage->getpixel))(sourceimage, sr, sc);
                  }
               }
            }
         }
         (*(image->putpixel))(image, r, c, sum);
      }
   }

   /****************************************************************************
   * Replace the cachemode with its previous value now that we are done updating
   * the cache.
   ****************************************************************************/
   image->cachemode = hold_cachemode;
}

/*******************************************************************************
* Procedure: sep_convolve
* Name: Michael Heath, University of South Florida
* Date: 7/14/98
*
typedef struct{
   int kernel_size;
   int direction;
   double *kernel_coefficients;
}ONE_D_CONVOLVE_DATA;
*
*******************************************************************************/
void one_d_convolve(CACHEIM *image, CACHEIM *source_image, int ksize,
    double *kernel, int direction, int max_depth, int max_cache_size,
    int datatype, char *mode)
{
   CACHEIM *newim=NULL;
   ONE_D_CONVOLVE_DATA *convolve_data;
   int i, r, c, precomputeimage, rows, cols;

   /* print_virtual_image_information(*source_image, "In one_d_convolve source_image."); */

   memset(image, 0, sizeof(CACHEIM));

   /****************************************************************************
   * If the mode is "r+b" show an error because we can not already have the
   * result of the convolution in a file.
   ****************************************************************************/
   if(strcmp(mode, "r+b") == 0){
      fprintf(stderr, "Error. Tried to open a convolve image in mode \"r+b\".\n");
      exit(1);
   }

   /****************************************************************************
   * Make sure that the user specified that the one-d kernel is either oriented
   * in the X or Y direction.
   ****************************************************************************/
   if(!((direction == XDIR) || (direction == YDIR))){
      printf("Invalid direction specification in one_d_convolve().\n");
      printf("   Specify %d for (1xN) or %d for (Nx1).\n", XDIR, YDIR);
      exit(1);
   }

   rows = source_image->rows;
   cols = source_image->cols;

   if((rows*cols*dt_size[datatype]) <= max_cache_size) precomputeimage = 1;
   else precomputeimage = 0;

   /****************************************************************************
   * Allocate a virtual image.
   ****************************************************************************/
   if((strcmp(mode, "rb") == 0) && (precomputeimage == 0)){
      if(allocate_cached_image(image, rows, cols, datatype, "tobecomputed",
         "rb", max_depth, (rows*cols*dt_size[datatype])/4, 0, FALSE,
         source_image->resolution) == 0) exit(1);
   }
   else{
      if(allocate_cached_image(image, rows, cols, datatype, "tobecomputed",
         "rb", max_depth, max_cache_size, 0, FALSE, source_image->resolution) == 0) exit(1);
   }

   /****************************************************************************
   * Set up the data necessary for computing the convolution.
   ****************************************************************************/
/*
   image->calculationdata = (void *) mymalloc("image->calculationdata", sizeof(ONE_D_CONVOLVE_DATA)+(ksize*sizeof(double)));

   convolve_data = (ONE_D_CONVOLVE_DATA *) image->calculationdata;
   convolve_data->kernel_size = ksize;
   convolve_data->direction = direction;
   convolve_data->kernel_coefficients = (double *)(((char *)(image->calculationdata)) + sizeof(ONE_D_CONVOLVE_DATA));
   for(i=0;i<ksize;i++){
      convolve_data->kernel_coefficients[i] = kernel[i];
   }
*/

   image->calculationdata = (void *) mymalloc("image->calculationdata", sizeof(ONE_D_CONVOLVE_DATA));

   convolve_data = (ONE_D_CONVOLVE_DATA *) image->calculationdata;
   convolve_data->kernel_size = ksize;
   convolve_data->direction = direction;
   convolve_data->kernel_coefficients = (double *) mycalloc("kernel_coeff", ksize, sizeof(double));
   for(i=0;i<ksize;i++){
      convolve_data->kernel_coefficients[i] = kernel[i];
   }

   image->sourceimlist = (void **) mycalloc("image->sourceimlist", 1, sizeof(CACHEIM *));
   ((CACHEIM **)(image->sourceimlist))[0] = source_image;

   image->get_cache_segment = get_sep_convolved_segment;

   /*
   print_virtual_image_information(*image, "In one_d_convolve image.");
   fflush(stdout);
   */

   /****************************************************************************
   * If the mode is rb then we are finished. We have told the image how to
   * compute its values.
   ****************************************************************************/
   if((strcmp(mode, "rb") == 0) && (precomputeimage == 0)) return;

   printf("Precalculating the entire sep_convolved image. (%d)\n", direction);

   /****************************************************************************
   * If the mode is "w+b" then we want to immediately calculate the result for
   * the entire image and cache the result in a file.
   ****************************************************************************/
   newim = (CACHEIM *) mycalloc("newim", 1, sizeof(CACHEIM));

   /****************************************************************************
   * Allocate a virtual image that is the entire size.
   ****************************************************************************/
   if(allocate_cached_image(newim, rows, cols, datatype, NULL, "w+b", max_depth,
      (rows*cols*dt_size[datatype]), 0, FALSE,
      source_image->resolution) == 0) exit(1);

   /****************************************************************************
   * Compute the result of the convolution for every pixel in the image.
   ****************************************************************************/
   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         (*(newim->putpixel))(newim, r, c, (*(image->getpixel))(image, r, c));
      }
   }

   deallocate_cached_image(image);
   *image = *newim;
   myfree("newim", newim);
}

/*******************************************************************************
* Procedure: get_sep_convolved_segment
* Purpose: This procedure calculates the output convolved values for a segment
* of the image. This function is called by get_pixel from virtual images
* of the type convolve.
* Name: Michael Heath, University of South Florida
* Date: 10/22/97
*******************************************************************************/
void get_sep_convolved_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg)
{
   int r=0, c=0, rr=0, cc=0, sr=0, sc=0, ksize=0, start_row=0, start_col=0, num_rows=0,
      num_cols=0, krows, kcols, direction;
   double sum=0.0;
   USHORT hold_cachemode;
   ONE_D_CONVOLVE_DATA *convolve_data;
   CACHEIM *sourceimage=NULL;
   double *kernel=NULL;

   /****************************************************************************
   * Determine the portion of the image that this segment falls within. This
   * area contains the pixels we will be calculating.
   ****************************************************************************/
   image->cachestructure[rseg][cseg][0].rpage = (USHORT)rpage;
   image->cachestructure[rseg][cseg][0].cpage = (USHORT)cpage;
   image->cachestructure[rseg][cseg][0].dirtybit = CLEAN_DATA;
   /* cstruct = image->cachestructure[rseg][cseg][0]; */

   start_row = (rpage * image->nrb + rseg) * image->rbs;
   if((image->rows - start_row) >= image->rbs) num_rows = image->cachestructure[rseg][cseg][0].rowlim = image->rbs;
   else num_rows = image->cachestructure[rseg][cseg][0].rowlim = image->rows - start_row;

   start_col = (cpage * image->ncb + cseg) * image->cbs;
   if((image->cols - start_col) >= image->cbs) num_cols = image->cachestructure[rseg][cseg][0].collim = image->cbs;
   else num_cols = image->cachestructure[rseg][cseg][0].collim = image->cols - start_col;

   /****************************************************************************
   * We know that the cache segment that we want to update is in the cache.
   * Therefore we will temporarily allow ourselves to write into the cache.
   * It is very important that we replace the cache mode with its proper value
   * when we are finished with this function.
   ****************************************************************************/
   hold_cachemode = image->cachemode;
   image->cachemode = CACHE_READ_WRITE;

   /****************************************************************************
   * Get the data necessary to do the convolution from the image structure.
   ****************************************************************************/
   convolve_data = (ONE_D_CONVOLVE_DATA *)(image->calculationdata);
   direction = convolve_data->direction;
   ksize = convolve_data->kernel_size;
   kernel = convolve_data->kernel_coefficients;
   sourceimage = (((CACHEIM **)image->sourceimlist)[0]);

   if(direction == XDIR){
      krows = 1;
      kcols = ksize;
   }
   if(direction == YDIR){
      krows = ksize;
      kcols = 1;
   }

   /*
   printf("get_sep_convolved_segment(image, %d, %d, %d, %d) %d %d %d %d\n",
   rpage, cpage, rseg, cseg, start_row, start_col, num_rows, num_cols);
   */

   /****************************************************************************
   * Do the convolution in this region of the image.
   ****************************************************************************/
   for(r=start_row;r<(start_row+num_rows);r++){
      for(c=start_col;c<(start_col+num_cols);c++){
         sum = 0;
         for(rr=0,sr=(r-krows/2);rr<krows;rr++,sr++){
            if((sr >= 0) && (sr < (image->rows))){
               for(cc=0,sc=(c-kcols/2);cc<kcols;cc++,sc++){
                  if((sc >= 0) && (sc < (image->cols))){
                     sum += kernel[rr*kcols+cc] * (*(sourceimage->getpixel))(sourceimage, sr, sc);
                  }
               }
            }
         }
         (*(image->putpixel))(image, r, c, sum);
      }
   }

   /****************************************************************************
   * Replace the cachemode with its previous value now that we are done updating
   * the cache.
   ****************************************************************************/
   image->cachemode = hold_cachemode;
}
