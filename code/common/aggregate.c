/*******************************************************************************
* PROGRAM: aggregate.c
* PURPOSE: This file contains the code for aggregating an image.
* NAME: Michael Heath, University of South Florida
* Note: One Numerical Recipes in C function is used in this code.
* DATE: 11/8/97
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include <string.h>
#include <stdlib.h>
#include "virtual_image.h"
#include "myalloc.h"

typedef struct{
   double aggfactor;
}AGGREGATION_DATA;

typedef struct{
   int aggfactor;
   double *localvals;
}AGGREGATION_MEDIAN_DATA;

void get_aggregate_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg);

void aggregate(CACHEIM *image, CACHEIM *source_image, double aggfactor,
    int max_depth, int max_cache_size, int datatype, char *mode);

int compare_double(const void *a, const void *b);
void aggregate_median(CACHEIM *image, CACHEIM *source_image, int aggfactor,
    int max_depth, int max_cache_size, int datatype, char *mode);
void get_aggregate_median_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg);

double select_kth(unsigned long k, unsigned long n, double arr[]);

void get_aggregate_double_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg);

/*******************************************************************************
* Procedure: aggregate
* Name: Michael Heath, University of South Florida
* Date: 11/8/97
*******************************************************************************/
void aggregate(CACHEIM *image, CACHEIM *source_image, double aggfactor,
    int max_depth, int max_cache_size, int datatype, char *mode)
{
   CACHEIM *newim=NULL;
   AGGREGATION_DATA *aggregation_data;
   int r, c;
   int agg_rows, agg_cols;
   int precomputeimage;

   /****************************************************************************
   * If the mode is "r+b" show an error because we can not already have the
   * result of the aggregation in a file.
   ****************************************************************************/
   if(strcmp(mode, "r+b") == 0){
      fprintf(stderr, "Error. Tried to open an aggregation image in mode \"r+b\".\n");
      exit(1);
   }

   /****************************************************************************
   * Determine how large the aggretated image will be.
   ****************************************************************************/
   agg_rows = ceil((double)(source_image->rows) / aggfactor);
   agg_cols = ceil((double)(source_image->cols) / aggfactor);

   if((agg_rows*agg_cols*dt_size[datatype]) <= max_cache_size) precomputeimage = 1;
   else precomputeimage = 0;

   /****************************************************************************
   * Allocate a virtual image. If the whole image is to be cached, then set
   * up a temporary image as 1/2 of the max_cache_size. This is necessary
   * because otherwise the computations will never be performed because
   * different getpixel and putpixel functions are used when the image is
   * completely contained in the image. Since those functions end up calling
   * the procedure that calculates the aggregates they would never be
   * called if we didn't do this.
   ****************************************************************************/
   if(precomputeimage){
      if(allocate_cached_image(image, agg_rows, agg_cols,
         datatype, "tobecomputed", "rb", max_depth,
	 (agg_rows * agg_cols * dt_size[datatype])/2, 0, FALSE,
         source_image->resolution*aggfactor) == 0) exit(1);
   }
   else{
      if(allocate_cached_image(image, agg_rows, agg_cols,
         datatype, "tobecomputed", "rb", max_depth, max_cache_size, 0, FALSE,
         source_image->resolution*aggfactor) == 0) exit(1);
   }

   /****************************************************************************
   * Set up the data necessary for computing the aggregation.
   ****************************************************************************/
   image->calculationdata = (void *) mycalloc("image->calculationdata", 1, sizeof(AGGREGATION_DATA));

   aggregation_data = (AGGREGATION_DATA *) image->calculationdata;
   aggregation_data->aggfactor = aggfactor;

   image->sourceimlist = (void **) mycalloc("image->sourceimlist", 1, sizeof(CACHEIM *));
   ((CACHEIM **)(image->sourceimlist))[0] = source_image;

   if(aggfactor == floor(aggfactor)) image->get_cache_segment = get_aggregate_segment;
   else image->get_cache_segment = get_aggregate_double_segment;

   /****************************************************************************
   * If the mode is rb and we are not going to precompute the image then we are
   * finished. We have told the image how to compute its values so we can return.
   ****************************************************************************/
   if((strcmp(mode, "rb") == 0) && (precomputeimage == 0)) return;

   /* printf("Precalculating the entire aggregation image.\n"); */

   /****************************************************************************
   * We must calculate the result for the entire image.
   ****************************************************************************/
   newim = (CACHEIM *) mycalloc("newim", 1, sizeof(CACHEIM));

   /****************************************************************************
   * Allocate a virtual image.
   ****************************************************************************/
   if(allocate_cached_image(newim, agg_rows, agg_cols,
      datatype, NULL, "w+b", max_depth, agg_rows*agg_cols*dt_size[datatype], 0, FALSE,
      source_image->resolution*aggfactor) == 0) exit(1);

   /****************************************************************************
   * Compute the result of the aggregation for every pixel in the image.
   ****************************************************************************/
   for(r=0;r<agg_rows;r++){
      for(c=0;c<agg_cols;c++){
         (*(newim->putpixel))(newim, r, c, (*(image->getpixel))(image, r, c));
      }
   }

   deallocate_cached_image(image);
   *image = *newim;
   myfree("newim", newim);
}

/*******************************************************************************
* Procedure: get_aggregate_segment
* Purpose: This procedure calculates the result of aggregation for a segment
* of the image. This function is called by get_pixel from virtual images
* of the type aggregation.
* Name: Michael Heath, University of South Florida
* Date: 11/8/97
*******************************************************************************/
void get_aggregate_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg)
{
   int r=0, c=0, rr=0, cc=0, sr=0, sc=0, aggfactor, start_row=0, start_col=0, num_rows=0, num_cols=0;
   double sum=0.0;
   USHORT hold_cachemode;
   AGGREGATION_DATA *aggregation_data;
   CACHEIM *sourceimage=NULL;
   int count=0;
   int source_rows, source_cols;

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
   * Get the data necessary to do the aggregation from the image structure.
   ****************************************************************************/
   aggregation_data = (AGGREGATION_DATA *)(image->calculationdata);
   aggfactor = (int)aggregation_data->aggfactor;
   sourceimage = (((CACHEIM **)image->sourceimlist)[0]);

   /*
   printf("get_aggregate_segment(image, %d, %d, %d, %d) %d %d %d %d\n", rpage, cpage,
   rseg, cseg, start_row, start_col, num_rows, num_cols);
   */

   source_rows = sourceimage->rows;
   source_cols = sourceimage->cols;

   /****************************************************************************
   * Do the aggregation in this region of the image.
   ****************************************************************************/
   for(r=start_row;r<(start_row+num_rows);r++){
      sr = r*aggfactor;
      for(c=start_col;c<(start_col+num_cols);c++){
	 sc = c*aggfactor;
	 count = 0;
         sum = 0.0;
         for(rr=sr;rr<(sr+aggfactor);rr++){
            if(rr < source_rows){
               for(cc=sc;cc<(sc+aggfactor);cc++){
                  if(cc < source_cols){
                     count++;
                     sum += (*(sourceimage->getpixel))(sourceimage, rr, cc);
                  }
               }
            }
         }
         (*(image->putpixel))(image, r, c, (sum/(double)count));
      }
   }

   /****************************************************************************
   * Replace the cachemode with its previous value now that we are done updating
   * the cache.
   ****************************************************************************/
   image->cachemode = hold_cachemode;
}

/*******************************************************************************
* Procedure: aggregate_median
* Purpose: This procedure does something like aggregate an image. In each
* aggfactor x aggfactor block of pixels, the median value (rather than the
* average value) is calculated. This median value is the value put in the
* output image for the block of pixels. A subsampling is basicall done by
* this procedure, but it is non-parametric in a sense. What I mean by this
* is that only the relative order of the pixels values (in the window) matters.
* Name: Michael Heath, University of South Florida
* Date: 7/28/98
*******************************************************************************/
void aggregate_median(CACHEIM *image, CACHEIM *source_image, int aggfactor,
    int max_depth, int max_cache_size, int datatype, char *mode)
{
   CACHEIM *newim=NULL;
   AGGREGATION_MEDIAN_DATA *aggregation_median_data;
   int r, c;
   int agg_rows, agg_cols;
   int precomputeimage;

   /****************************************************************************
   * If the mode is "r+b" show an error because we can not already have the
   * result of the aggregation in a file.
   ****************************************************************************/
   if(strcmp(mode, "r+b") == 0){
      fprintf(stderr, "Error. Tried to open an aggregation image in mode \"r+b\".\n");
      exit(1);
   }

   /****************************************************************************
   * Determine how large the aggretated image will be.
   ****************************************************************************/
   agg_rows = ceil((double)(source_image->rows) / (double)aggfactor);
   agg_cols = ceil((double)(source_image->cols) / (double)aggfactor);

   if((agg_rows*agg_cols*dt_size[datatype]) <= max_cache_size) precomputeimage = 1;
   else precomputeimage = 0;

   /****************************************************************************
   * Allocate a virtual image. If thw whole image is to be cached, then set
   * up a temporary image as 1/2 of the max_cache_size. This is necessary
   * because otherwise the computations will never be performed because
   * differt getpixel and putpixel functions are used when the image is
   * completely contained in the image. Since those functions end up calling
   * the procedure that calculates the aggregates they would never be
   * called if we didn't do this.
   ****************************************************************************/
   if(precomputeimage){
      if(allocate_cached_image(image, agg_rows, agg_cols,
         datatype, "tobecomputed", "rb", max_depth,
	 (agg_rows * agg_cols * dt_size[datatype])/2, 0, FALSE,
         source_image->resolution*aggfactor) == 0) exit(1);
   }
   else{
      if(allocate_cached_image(image, agg_rows, agg_cols,
         datatype, "tobecomputed", "rb", max_depth, max_cache_size, 0, FALSE,
         source_image->resolution*aggfactor) == 0) exit(1);
   }

   /****************************************************************************
   * Set up the data necessary for computing the aggregation.
   ****************************************************************************/
   image->calculationdata = (void *) mymalloc("image->calculationdata",
      aggfactor*aggfactor*sizeof(double) + sizeof(AGGREGATION_MEDIAN_DATA));

   memset(image->calculationdata, 0, aggfactor*aggfactor*sizeof(double) + sizeof(AGGREGATION_MEDIAN_DATA));

   aggregation_median_data = (AGGREGATION_MEDIAN_DATA *) image->calculationdata;
   aggregation_median_data->aggfactor = aggfactor;
   aggregation_median_data->localvals = (double *)((unsigned char *)image->calculationdata + sizeof(AGGREGATION_MEDIAN_DATA));

   image->sourceimlist = (void **) mycalloc("image->sourceimlist", 1, sizeof(CACHEIM *));
   ((CACHEIM **)(image->sourceimlist))[0] = source_image;

   image->get_cache_segment = get_aggregate_median_segment;

   /****************************************************************************
   * If the mode is rb and we are not going to precompute the image then we are
   * finished. We have told the image how to compute its values so we can return.
   ****************************************************************************/
   if((strcmp(mode, "rb") == 0) && (precomputeimage == 0)) return;

   printf("Precalculating the entire aggregation_median image.\n");

   /****************************************************************************
   * If the mode is "w+b" then we want to immediately calculate the result for
   * the entire image and cache the result in a file.
   ****************************************************************************/
   newim = (CACHEIM *) mycalloc("newim", 1, sizeof(CACHEIM));

   /****************************************************************************
   * Allocate a virtual image.
   ****************************************************************************/
   if(allocate_cached_image(newim, agg_rows, agg_cols,
      datatype, NULL, "w+b", max_depth, agg_rows*agg_cols*dt_size[datatype], 0, FALSE,
      source_image->resolution*aggfactor) == 0) exit(1);

   /****************************************************************************
   * Compute the result of the aggregation for every pixel in the image.
   ****************************************************************************/
   for(r=0;r<agg_rows;r++){
      for(c=0;c<agg_cols;c++){
         (*(newim->putpixel))(newim, r, c, (*(image->getpixel))(image, r, c));
      }
   }

   deallocate_cached_image(image);
   *image = *newim;
   myfree("newim", newim);
}

/*******************************************************************************
* Procedure: get_aggregate_median_segment
* Purpose: This procedure calculates the result of median aggregation for a
* segment of the image. This function is called by get_pixel from virtual images
* of the type aggregation.
* Name: Michael Heath, University of South Florida
* Date: 7/28/98
*******************************************************************************/
void get_aggregate_median_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg)
{
   int r=0, c=0, rr=0, cc=0, sr=0, sc=0, aggfactor, start_row=0, start_col=0, num_rows=0, num_cols=0;
   USHORT hold_cachemode;
   AGGREGATION_MEDIAN_DATA *aggregation_median_data;
   CACHEIM *sourceimage=NULL;
   int count=0;
   int source_rows, source_cols;
   double *localvals=NULL;

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
   * Get the data necessary to do the aggregation from the image structure.
   ****************************************************************************/
   aggregation_median_data = (AGGREGATION_MEDIAN_DATA *)(image->calculationdata);
   aggfactor = aggregation_median_data->aggfactor;
   localvals = aggregation_median_data->localvals;
   sourceimage = (((CACHEIM **)image->sourceimlist)[0]);

   /*
   printf("get_aggregate_median_segment(image, %d, %d, %d, %d) %d %d %d %d\n", rpage, cpage,
      rseg, cseg, start_row, start_col, num_rows, num_cols);
   */

   source_rows = sourceimage->rows;
   source_cols = sourceimage->cols;

   /****************************************************************************
   * Do the aggregation in this region of the image.
   ****************************************************************************/
   for(r=start_row;r<(start_row+num_rows);r++){
      sr = r*aggfactor;
      for(c=start_col;c<(start_col+num_cols);c++){
	 sc = c*aggfactor;
	 count = 0;
         for(rr=sr;rr<(sr+aggfactor);rr++){
            if(rr < source_rows){
               for(cc=sc;cc<(sc+aggfactor);cc++){
                  if(cc < source_cols){
                     localvals[count] = (*(sourceimage->getpixel))(sourceimage, rr, cc);
                     count++;
                  }
               }
            }
         }
         (*(image->putpixel))(image, r, c, select_kth((int)ceil(count/2.0), count, localvals-1));

         /*
         qsort(localvals, count, sizeof(double), compare_double);
         (*(image->putpixel))(image, r, c, (localvals[count/2]));
         */
      }
   }

   /****************************************************************************
   * Replace the cachemode with its previous value now that we are done updating
   * the cache.
   ****************************************************************************/
   image->cachemode = hold_cachemode;
}

int compare_double(const void *a, const void *b)
{
   if((*((double *)a)) < (*((double *)b))) return(-1);
   else if((*((double *)a)) > (*((double *)b))) return(1);
   else return(0);
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

double select_kth(unsigned long k, unsigned long n, double arr[])
{
	unsigned long i,ir,j,l,mid;
	double a,temp;

	l=1;
	ir=n;
	for (;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
				SWAP(arr[l],arr[ir])
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;
			SWAP(arr[mid],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
			}
			arr[l]=arr[j];
			arr[j]=a;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}
	}
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software -)+%R@(. */


/*******************************************************************************
* Procedure: get_aggregate_double_segment
* Purpose: This procedure calculates the result of aggregation for a segment
* of the image. This function is called by get_pixel from virtual images
* of the type aggregation. It is used when non integer aggregation values
* are used.
* Name: Michael Heath, University of South Florida
* Date: 9/3/98
*******************************************************************************/
void get_aggregate_double_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg)
{
   int r=0, c=0, rr=0, cc=0, start_row=0, start_col=0, num_rows=0, num_cols=0;
   USHORT hold_cachemode;
   AGGREGATION_DATA *aggregation_data;
   CACHEIM *sourceimage=NULL;
   int source_rows, source_cols;
   double rt, rb, cl, cr;
   double mt, mb, ml, mr;
   double sumpixel, sumweight, mult, val = 1.0;
   double aggfactor;

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
   * Get the data necessary to do the aggregation from the image structure.
   ****************************************************************************/
   aggregation_data = (AGGREGATION_DATA *)(image->calculationdata);
   aggfactor = aggregation_data->aggfactor;
   sourceimage = (((CACHEIM **)image->sourceimlist)[0]);

   /*
   printf("get_aggregate_double_segment(image, %d, %d, %d, %d) %d %d %d %d\n", rpage, cpage,
   rseg, cseg, start_row, start_col, num_rows, num_cols);
   */

   source_rows = sourceimage->rows;
   source_cols = sourceimage->cols;

   for(r=start_row;r<(start_row+num_rows);r++){
      rt = (double)r * aggfactor;
      rb = (double)(r + 1) * aggfactor;
      mt = ceil(rt) - rt;
      mb = rb - floor(rb);
      for(c=start_col;c<(start_col+num_cols);c++){

         cl = (double)c * aggfactor;
         cr = (double)(c + 1) * aggfactor;
         ml = ceil(cl) - cl;
         mr = cr - floor(cr);

         /* printf("[%3d][%3d] <= (%5.2f,%5.2f),(%5.2f,%5.2f)\n", big_r, big_c, rt, cl, rb, cr); */

         sumpixel = 0.0;
         sumweight = 0.0;

         for(rr=(int)floor(rt);rr<(int)ceil(rb);rr++){
            for(cc=(int)floor(cl);cc<(int)ceil(cr);cc++){

               if((rr>=0)&&(cc>=0)&&(rr<source_rows)&&(cc<source_cols)){

                  mult = 1.0;
                  val = (*(sourceimage->getpixel))(sourceimage, rr, cc);

                  if((double)rr < ceil(rt)) mult *= mt;
                  else if((double)rr >= floor(rb)) mult *= mb;
                  if((double)cc < ceil(cl)) mult *= ml;
                  else if((double)cc >= floor(cr)) mult *= mr;

                  sumweight += mult;
                  sumpixel += mult * val;
               }
            }
         }

         (*(image->putpixel))(image, r, c, (sumpixel / sumweight));
      }
   }

   /****************************************************************************
   * Replace the cachemode with its previous value now that we are done updating
   * the cache.
   ****************************************************************************/
   image->cachemode = hold_cachemode;
}
