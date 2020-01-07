/*******************************************************************************
* PROGRAM: cc_label.c
* PURPOSE: This file contains the code for connected component labeling an image.
* It works on gray-level images, so every pixel will be part of some component.
* This means that there is no background value like you have in binary-image
* connected component algorithms. Virtual images are used, so the image to process
* can be any datatype including floating point. Even though virtual images
* are used, the entire image is processed because it is the nature of the
* algorithm to process the entire connected component image at once. Other
* virtual image algorithms update only a block of the image as it is needed.
* NAME: Michael Heath, University of South Florida
* DATE: 2/15/98
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include "virtual_image.h"
#include <string.h>
#include <limits.h>
#include "myalloc.h"

#define VERBOSE 0

typedef struct{
   ULONG old_label, new_label;
   double value;
}LABELINFO;

void cc_label_0bkgd(CACHEIM *image, CACHEIM *source_image, char *filename,
    int max_depth, int max_cache_size, int datatype, char *mode, ULONG *numcomp);

void cc_label(CACHEIM *image, CACHEIM *source_image, char *filename,
    int max_depth, int max_cache_size, int datatype, char *mode, ULONG *numcomp);

void cc_relabel_by_intensity(CACHEIM *label_image, CACHEIM *intensity_image,
   ULONG num_components);

void cc_label_4connect(CACHEIM *image, CACHEIM *source_image, char *filename,
    int max_depth, int max_cache_size, int datatype, char *mode, ULONG *numcomp);

/*******************************************************************************
* Procedure: cc_label
* Name: Michael Heath, University of South Florida
* Date: 2/15/98
*******************************************************************************/
void cc_label(CACHEIM *image, CACHEIM *source_image, char *filename,
    int max_depth, int max_cache_size, int datatype, char *mode, ULONG *numcomp)
{
   int r, c, rr, pos, n, n_x, n_y, source_rows, source_cols;
   char necessary_mode[] = {'w', '+', 'b', '\0'};
   int time;
   ULONG comp=0;

   int fsa_x[4] = {-1,  1,  0, -1};
   int fsa_y[4] = { 0, -1, -1, -1};

   int rsa_x[4] = { 1, -1,  0,  1};
   int rsa_y[4] = { 0,  1,  1,  1};

   double value, neighbor_value;
   unsigned long int nextlabel, min_label, neighbor_label, max_label;

   unsigned long int maxval;

   ULONG *labels = NULL;

   source_rows = source_image->rows;
   source_cols = source_image->cols;

   /****************************************************************************
   * The conected component labeled image should be some integer data type. It's
   * datatype should not be float or double. The mode should be "w+b". This will
   * create a disk file image to use for cacheing and for storing the labeled
   * image if the file does not fit entirely in memory.
   ****************************************************************************/
   if((datatype == FLOATNUM) || (datatype == DOUBLENUM)) datatype = ULONGNUM;
   mode = necessary_mode;

   switch(datatype){
      case UCHARNUM:
         maxval = UCHAR_MAX;
         break;
      case CHARNUM:
         maxval = SCHAR_MAX;
         break;
      case USHORTNUM:
         maxval = USHRT_MAX;
         break;
      case SHORTNUM:
         maxval = SHRT_MAX;
         break;
      case ULONGNUM:
         maxval = ULONG_MAX;
         break;
      case LONGNUM:
         maxval = LONG_MAX;
         break;
      default:
         fprintf(stderr, "Error in cc_label(). Datatype not supported.\n");
         exit(1);
   }

   /****************************************************************************
   * Allocate a virtual image.
   ****************************************************************************/
   if(allocate_cached_image(image, source_rows, source_cols,
      datatype, filename, mode, max_depth, max_cache_size, 0, FALSE,
      source_image->resolution) == 0) exit(1);

   if(VERBOSE){
      printf("rbs = %d cbs = %d\n", image->rbs, image->cbs);
      printf("nrb = %d ncb = %d\n", image->nrb, image->ncb);
      printf("depth = %d\n", image->depth);
      printf("rbs_bits = %d cbs_bits = %d\n", image->rbs_bits, image->cbs_bits);
      printf("nrb_bits = %d ncb_bits = %d\n", image->nrb_bits, image->ncb_bits);
      printf("rbs_and_bits = %d cbs_and_bits = %d\n", image->rbs_and_bits, image->cbs_and_bits);
      printf("nrb_and_bits = %d ncb_and_bits = %d\n", image->nrb_and_bits, image->ncb_and_bits);
   }

   /****************************************************************************
   * Set up the data necessary for computing the connected compont labeled image.
   * Unlike some algorithms, this algorithm has no parameters.
   ****************************************************************************/
   image->calculationdata = (void *) NULL;

   /****************************************************************************
   * Set up a pointer to the source image.
   ****************************************************************************/
   image->sourceimlist = (void **) mycalloc("image->sourceimlist", 1, sizeof(CACHEIM *));
   ((CACHEIM **)(image->sourceimlist))[0] = source_image;

   if(VERBOSE) printf("Precalculating the entire cc-labeled image.\n");

   /****************************************************************************
   * Do the connected component labeling. We use a two pass procedure in which
   * the equivalence table is resolved and applied every time an equivalence is
   * found.
   ****************************************************************************/
   nextlabel = 0;

   /****************************************************************************
   * Do the foreward pass of the connected component labeling algorithm.
   ****************************************************************************/
   for(r=0,pos=0;r<source_rows;r++){
      for(c=0;c<source_cols;c++,pos++){

         value = (*(source_image->getpixel))(source_image, r, c);

         time = 0;
	 for(n=0;n<4;n++){

            n_x = c + fsa_x[n];
            n_y = r + fsa_y[n];

            if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

               neighbor_value = (*(source_image->getpixel))(source_image,
                                   n_y, n_x);

               if(value == neighbor_value){

                  neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                   n_y, n_x));

                  if(time == 0){
                     max_label = min_label = neighbor_label;
                  }
                  if(neighbor_label < min_label) min_label = neighbor_label;
                  if(neighbor_label > max_label) max_label = neighbor_label;
                  time++;
               }
            }
         }

         if(time == 0){
            (*(image->putpixel))(image, r, c, (double)nextlabel);

	    nextlabel++;
            if(nextlabel > maxval){
               fprintf(stderr, "Error in cc_label(). Nextlabel overflowed ");
               fprintf(stderr, "the size allowable by your datatype selection.\n");
               exit(1);
            }
         }
         else if(max_label == min_label){
            (*(image->putpixel))(image, r, c, (double)min_label);
         }
         else{

            (*(image->putpixel))(image, r, c, (double)min_label);

            for(rr=(pos-1);rr>=(pos-source_cols);rr--){

               n_x = rr % source_cols;
               n_y = rr / source_cols;

               if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

                  neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                     n_y, n_x));

                  if(neighbor_label == max_label){
                     (*(image->putpixel))(image, n_y, n_x, (double)min_label);
                  }
               }
            }
         }
      }
   }

   if(VERBOSE) printf("The maximum label is less than or equal to %lu\n", nextlabel-1);

   /****************************************************************************
   * Do the reverse pass of the connected component labeling algorithm.
   ****************************************************************************/
   for(r=(source_rows-1);r>=0;r--){
      for(c=(source_cols-1);c>=0;c--){

         pos = r * source_cols + c;

         value = (*(source_image->getpixel))(source_image, r, c);

         max_label = min_label = (unsigned long)((*(image->getpixel))(image, r, c));

	 for(n=0;n<4;n++){

            n_x = c + rsa_x[n];
            n_y = r + rsa_y[n];

            if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

               neighbor_value = (*(source_image->getpixel))(source_image,
                                 n_y, n_x);

               if(value == neighbor_value){

                  neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                   n_y, n_x));


                  if(neighbor_label < min_label) min_label = neighbor_label;
                  if(neighbor_label > max_label) max_label = neighbor_label;
               }
            }
         }

         if(max_label != min_label){

            (*(image->putpixel))(image, r, c, (double)min_label);

            for(rr=(pos+1);rr<=(pos+source_cols);rr++){

               n_x = rr % source_cols;
               n_y = rr / source_cols;

               if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

                  neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                   n_y, n_x));

                  if(neighbor_label == max_label){
                     (*(image->putpixel))(image, n_y, n_x, (double)min_label);
                  }
               }
            }
         }
      }
   }

   /****************************************************************************
   * Go through and relabel all of the components with values that are 
   * 0 -> (#unique_components - 1).
   ****************************************************************************/
   if((labels = (ULONG *) mycalloc("labels", nextlabel, sizeof(ULONG))) == NULL){
      fprintf(stderr, "Malloc error in cc_label().\n");
      exit(1);
   }

   for(r=0,pos=0;r<source_rows;r++){
      for(c=0;c<source_cols;c++,pos++){
         neighbor_label = (ULONG)((*(image->getpixel))(image, r, c));
         labels[neighbor_label]++;
      }
   }

   for(r=0,comp=0;r<nextlabel;r++){
      if(labels[r] != 0){
         labels[r] = comp;
         comp++;
      }
   }

   *numcomp = comp;
   if(VERBOSE) printf("There are actually %lu unique component labels in the image.\n", *numcomp);

   for(r=(source_rows-1);r>=0;r--){
      for(c=(source_cols-1);c>=0;c--){
         neighbor_label = (ULONG)((*(image->getpixel))(image, r, c));
         (*(image->putpixel))(image, r, c, (double)(labels[neighbor_label]));
      }
   }

   myfree("labels", labels);
}

/*******************************************************************************
* Procedure: cc_label_0bkgd
* Name: Michael Heath, University of South Florida
* Date: 2/15/98
*******************************************************************************/
void cc_label_0bkgd(CACHEIM *image, CACHEIM *source_image, char *filename,
    int max_depth, int max_cache_size, int datatype, char *mode, ULONG *numcomp)
{
   int r, c, rr, pos, n, n_x, n_y, source_rows, source_cols;
   char necessary_mode[] = {'w', '+', 'b', '\0'};
   int time;
   ULONG comp=0;

   int fsa_x[4] = {-1,  1,  0, -1};
   int fsa_y[4] = { 0, -1, -1, -1};

   int rsa_x[4] = { 1, -1,  0,  1};
   int rsa_y[4] = { 0,  1,  1,  1};

   double value, neighbor_value;
   unsigned long int nextlabel, min_label, neighbor_label, max_label;

   unsigned long int maxval;

   ULONG *labels = NULL;

   source_rows = source_image->rows;
   source_cols = source_image->cols;

   /****************************************************************************
   * The conected component labeled image should be some integer data type. It's
   * datatype should not be float or double. The mode should be "w+b". This will
   * create a disk file image to use for cacheing and for storing the labeled
   * image if the file does not fit entirely in memory.
   ****************************************************************************/
   if((datatype == FLOATNUM) || (datatype == DOUBLENUM)) datatype = ULONGNUM;
   mode = necessary_mode;

   switch(datatype){
      case UCHARNUM:
         maxval = UCHAR_MAX;
         break;
      case CHARNUM:
         maxval = SCHAR_MAX;
         break;
      case USHORTNUM:
         maxval = USHRT_MAX;
         break;
      case SHORTNUM:
         maxval = SHRT_MAX;
         break;
      case ULONGNUM:
         maxval = ULONG_MAX;
         break;
      case LONGNUM:
         maxval = LONG_MAX;
         break;
      default:
         fprintf(stderr, "Error in cc_label_0bkgd(). Datatype not supported.\n");
         exit(1);
   }

   /****************************************************************************
   * Allocate a virtual image.
   ****************************************************************************/
   if(allocate_cached_image(image, source_rows, source_cols,
      datatype, filename, mode, max_depth, max_cache_size, 0, FALSE,
      source_image->resolution) == 0) exit(1);

   if(VERBOSE){
      printf("rbs = %d cbs = %d\n", image->rbs, image->cbs);
      printf("nrb = %d ncb = %d\n", image->nrb, image->ncb);
      printf("depth = %d\n", image->depth);
      printf("rbs_bits = %d cbs_bits = %d\n", image->rbs_bits, image->cbs_bits);
      printf("nrb_bits = %d ncb_bits = %d\n", image->nrb_bits, image->ncb_bits);
      printf("rbs_and_bits = %d cbs_and_bits = %d\n", image->rbs_and_bits, image->cbs_and_bits);
      printf("nrb_and_bits = %d ncb_and_bits = %d\n", image->nrb_and_bits, image->ncb_and_bits);
   }

   /****************************************************************************
   * Set up the data necessary for computing the connected compont labeled image.
   * Unlike some algorithms, this algorithm has no parameters.
   ****************************************************************************/
   image->calculationdata = (void *) NULL;

   /****************************************************************************
   * Set up a pointer to the source image.
   ****************************************************************************/
   image->sourceimlist = (void **) mycalloc("image->sourceimlist", 1, sizeof(CACHEIM *));
   ((CACHEIM **)(image->sourceimlist))[0] = source_image;

   if(VERBOSE) printf("Precalculating the entire cc-labeled image.\n");

   /****************************************************************************
   * Do the connected component labeling. We use a two pass procedure in which
   * the equivalence table is resolved and applied every time an equivalence is
   * found.
   ****************************************************************************/
   nextlabel = 1;

   /****************************************************************************
   * Do the foreward pass of the connected component labeling algorithm.
   ****************************************************************************/
   for(r=0,pos=0;r<source_rows;r++){
      for(c=0;c<source_cols;c++,pos++){

         value = (*(source_image->getpixel))(source_image, r, c);

         if(value != 0.0){

            time = 0;
	    for(n=0;n<4;n++){

               n_x = c + fsa_x[n];
               n_y = r + fsa_y[n];

               if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

                  neighbor_value = (*(source_image->getpixel))(source_image,
                                      n_y, n_x);

                  if(value == neighbor_value){

                     neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                      n_y, n_x));

                     if(time == 0){
                        max_label = min_label = neighbor_label;
                     }
                     if(neighbor_label < min_label) min_label = neighbor_label;
                     if(neighbor_label > max_label) max_label = neighbor_label;
                     time++;
                  }
               }
            }

            if(time == 0){
               (*(image->putpixel))(image, r, c, (double)nextlabel);

	       nextlabel++;
               if(nextlabel > maxval){
                  fprintf(stderr, "Error in cc_label_0bkgd(). Nextlabel overflowed ");
                  fprintf(stderr, "the size allowable by your datatype selection.\n");
                  exit(1);
               }
            }
            else if(max_label == min_label){
               (*(image->putpixel))(image, r, c, (double)min_label);
            }
            else{

               (*(image->putpixel))(image, r, c, (double)min_label);

               for(rr=(pos-1);rr>=(pos-source_cols);rr--){

                  n_x = rr % source_cols;
                  n_y = rr / source_cols;

                  if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

                     neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                        n_y, n_x));

                     if(neighbor_label == max_label){
                        (*(image->putpixel))(image, n_y, n_x, (double)min_label);
                     }
                  }
               }
            }
         }
      }
   }

   if(VERBOSE) printf("The maximum label is less than or equal to %lu\n", nextlabel-1);

   /****************************************************************************
   * Do the reverse pass of the connected component labeling algorithm.
   ****************************************************************************/
   for(r=(source_rows-1);r>=0;r--){
      for(c=(source_cols-1);c>=0;c--){

         pos = r * source_cols + c;

         value = (*(source_image->getpixel))(source_image, r, c);

	 if(value != 0.0){

            max_label = min_label = (unsigned long)((*(image->getpixel))(image, r, c));

	    for(n=0;n<4;n++){

               n_x = c + rsa_x[n];
               n_y = r + rsa_y[n];

               if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

                  neighbor_value = (*(source_image->getpixel))(source_image,
                                    n_y, n_x);

                  if(value == neighbor_value){

                     neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                      n_y, n_x));


                     if(neighbor_label < min_label) min_label = neighbor_label;
                     if(neighbor_label > max_label) max_label = neighbor_label;
                  }
               }
            }

            if(max_label != min_label){

               (*(image->putpixel))(image, r, c, (double)min_label);

               for(rr=(pos+1);rr<=(pos+source_cols);rr++){

                  n_x = rr % source_cols;
                  n_y = rr / source_cols;

                  if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

                     neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                      n_y, n_x));

                     if(neighbor_label == max_label){
                        (*(image->putpixel))(image, n_y, n_x, (double)min_label);
                     }
                  }
               }
            }
         }
      }
   }

   /****************************************************************************
   * Go through and relabel all of the components with values that are 
   * 0 -> (#unique_components - 1).
   ****************************************************************************/
   if((labels = (ULONG *) mycalloc("labels", nextlabel, sizeof(ULONG))) == NULL){
      fprintf(stderr, "Malloc error in cc_label_0bkgd().\n");
      exit(1);
   }

   for(r=0,pos=0;r<source_rows;r++){
      for(c=0;c<source_cols;c++,pos++){
         neighbor_label = (ULONG)((*(image->getpixel))(image, r, c));
         labels[neighbor_label]++;
      }
   }

   for(r=0,comp=0;r<nextlabel;r++){
      if(labels[r] != 0){
         /* printf("labels[%d] = %lu  --> %lu\n", r, labels[r], comp); */
         labels[r] = comp;
         comp++;
      }
   }

   *numcomp = comp;
   if(VERBOSE) printf("There are actually %lu unique component labels in the image.\n", *numcomp);

   for(r=(source_rows-1);r>=0;r--){
      for(c=(source_cols-1);c>=0;c--){
         neighbor_label = (ULONG)((*(image->getpixel))(image, r, c));
         (*(image->putpixel))(image, r, c, (double)(labels[neighbor_label]));
      }
   }

   myfree("labels", labels);
}

/*******************************************************************************
* Procedure: cc_relabel_by_intensity
* Purpose: This procedure re-numbers the connected component labels in order of
* increasing image intensity values.
* Name: Michael Heath, University of South Florida
* Date: 3/2/98
*******************************************************************************/
void cc_relabel_by_intensity(CACHEIM *label_image, CACHEIM *intensity_image,
   ULONG num_components)
{
   LABELINFO *labelinfo = NULL;
   int r, c, rows, cols;
   ULONG this_label;
   double this_value;
   int compare_intensities(const void *, const void *);
   int compare_old_labels(const void *, const void *);

   rows = label_image->rows;
   cols = label_image->cols;

   if(VERBOSE) printf("There are %lu components in cc_relabel_by_intensity().\n", num_components);

   /****************************************************************************
   * Allocate an array of structures that will be used to relabel the cc_label
   * values. This array will serve as a LUT.
   ****************************************************************************/
   if((labelinfo = (LABELINFO *) mycalloc("labelinfo", num_components, sizeof(LABELINFO))) == NULL){
      fprintf(stderr, "Calloc error in cc_relabel_by_intensity().\n");
      exit(1);
   }

   /****************************************************************************
   * Go through the image and record the largest intensity value of each connected
   * component.
   ****************************************************************************/
   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         this_label = (ULONG)((*(label_image->getpixel))(label_image, r, c));
         if(this_label >= num_components) printf("%lu ", this_label);

         this_value = ((*(intensity_image->getpixel))(intensity_image, r, c));
         if(this_value > labelinfo[this_label].value) labelinfo[this_label].value = this_value;
      }
   }

   labelinfo[0].value = 0; /* Assume the background has a component value of 0. */

   /****************************************************************************
   * Record the value of each "old" connected component label. These are the
   * current values of the components. Keeping trach of them is necessary for
   * resorting the array of structures after we initially sort on the intensity
   * values.
   ****************************************************************************/
   for(r=0;r<num_components;r++) labelinfo[r].old_label = r;

   /****************************************************************************
   * Use the built in quick sort algorithm to sort the array based on increasing
   * intensity values.
   ****************************************************************************/
   if(VERBOSE) printf("First qsort().\n");
   qsort((void *)labelinfo, (size_t)(num_components), (size_t)sizeof(LABELINFO), compare_intensities);

   /****************************************************************************
   * Assign new labels with the lowest image intensity values recieving the
   * lowest values label numbers.
   ****************************************************************************/
   for(r=0;r<num_components;r++) labelinfo[r].new_label = r;

   /****************************************************************************
   * Use the built in quick sort algorithm to sort the array based on increasing
   * old_label values.
   ****************************************************************************/
   if(VERBOSE) printf("Second qsort().\n");
   qsort((void *)labelinfo, (size_t)(num_components), (size_t)sizeof(LABELINFO), compare_old_labels);

   if(VERBOSE) printf("Relabeling the image.\n");
   fflush(stdout);
   /****************************************************************************
   * Go through the image and apply our LUT to relabel each component with its
   * new component value.
   ****************************************************************************/
   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         this_label = (ULONG)((*(label_image->getpixel))(label_image, r, c));
         this_label = labelinfo[this_label].new_label;
         (*(label_image->putpixel))(label_image, r, c, (double)this_label);
      }
   }

   myfree("labelinfo", labelinfo);
}

int compare_intensities(const void *arg1, const void *arg2)
{
   double diff;
   diff = (((LABELINFO *)arg1)->value) - (((LABELINFO *)arg2)->value);

   if(diff < 0) return(-1);
   else if(diff > 0) return(1);
   else return(0);
}

int compare_old_labels(const void *arg1, const void *arg2)
{
   LONG diff;
   diff = (((LABELINFO *)arg1)->old_label) - (((LABELINFO *)arg2)->old_label);

   if(diff < 0) return(-1);
   else if(diff > 0) return(1);
   else return(0);
}

/*******************************************************************************
* Procedure: cc_label_4connect
* Name: Michael Heath, University of South Florida
* Date: 2/15/98
*******************************************************************************/
void cc_label_4connect(CACHEIM *image, CACHEIM *source_image, char *filename,
    int max_depth, int max_cache_size, int datatype, char *mode, ULONG *numcomp)
{
   int r, c, rr, pos, n, n_x, n_y, source_rows, source_cols;
   char necessary_mode[] = {'w', '+', 'b', '\0'};
   int time;
   ULONG comp=0;

   int fsa_x[4] = {-1,  1,  0, -1};
   int fsa_y[4] = { 0, -1, -1, -1};

   int rsa_x[4] = { 1, -1,  0,  1};
   int rsa_y[4] = { 0,  1,  1,  1};

   double value, neighbor_value;
   unsigned long int nextlabel, min_label, neighbor_label, max_label;

   unsigned long int maxval;

   ULONG *labels = NULL;

   source_rows = source_image->rows;
   source_cols = source_image->cols;

   /****************************************************************************
   * The conected component labeled image should be some integer data type. It's
   * datatype should not be float or double. The mode should be "w+b". This will
   * create a disk file image to use for cacheing and for storing the labeled
   * image if the file does not fit entirely in memory.
   ****************************************************************************/
   if((datatype == FLOATNUM) || (datatype == DOUBLENUM)) datatype = ULONGNUM;
   mode = necessary_mode;

   switch(datatype){
      case UCHARNUM:
         maxval = UCHAR_MAX;
         break;
      case CHARNUM:
         maxval = SCHAR_MAX;
         break;
      case USHORTNUM:
         maxval = USHRT_MAX;
         break;
      case SHORTNUM:
         maxval = SHRT_MAX;
         break;
      case ULONGNUM:
         maxval = ULONG_MAX;
         break;
      case LONGNUM:
         maxval = LONG_MAX;
         break;
      default:
         fprintf(stderr, "Error in cc_label_4connect(). Datatype not supported.\n");
         exit(1);
   }

   /****************************************************************************
   * Allocate a virtual image.
   ****************************************************************************/
   if(allocate_cached_image(image, source_rows, source_cols,
      datatype, filename, mode, max_depth, max_cache_size, 0, FALSE,
      source_image->resolution) == 0) exit(1);

   if(VERBOSE){
      printf("rbs = %d cbs = %d\n", image->rbs, image->cbs);
      printf("nrb = %d ncb = %d\n", image->nrb, image->ncb);
      printf("depth = %d\n", image->depth);
      printf("rbs_bits = %d cbs_bits = %d\n", image->rbs_bits, image->cbs_bits);
      printf("nrb_bits = %d ncb_bits = %d\n", image->nrb_bits, image->ncb_bits);
      printf("rbs_and_bits = %d cbs_and_bits = %d\n", image->rbs_and_bits, image->cbs_and_bits);
      printf("nrb_and_bits = %d ncb_and_bits = %d\n", image->nrb_and_bits, image->ncb_and_bits);
   }

   /****************************************************************************
   * Set up the data necessary for computing the connected compont labeled image.
   * Unlike some algorithms, this algorithm has no parameters.
   ****************************************************************************/
   image->calculationdata = (void *) NULL;

   /****************************************************************************
   * Set up a pointer to the source image.
   ****************************************************************************/
   image->sourceimlist = (void **) mycalloc("image->sourceimlist", 1, sizeof(CACHEIM *));
   ((CACHEIM **)(image->sourceimlist))[0] = source_image;

   if(VERBOSE) printf("Precalculating the entire cc-labeled image.\n");

   /****************************************************************************
   * Do the connected component labeling. We use a two pass procedure in which
   * the equivalence table is resolved and applied every time an equivalence is
   * found.
   ****************************************************************************/
   nextlabel = 0;

   /****************************************************************************
   * Do the foreward pass of the connected component labeling algorithm.
   ****************************************************************************/
   for(r=0,pos=0;r<source_rows;r++){
      for(c=0;c<source_cols;c++,pos++){

         value = (*(source_image->getpixel))(source_image, r, c);

         time = 0;
	 for(n=0;n<4;n+=2){

            n_x = c + fsa_x[n];
            n_y = r + fsa_y[n];

            if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

               neighbor_value = (*(source_image->getpixel))(source_image,
                                   n_y, n_x);

               if(value == neighbor_value){

                  neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                   n_y, n_x));

                  if(time == 0){
                     max_label = min_label = neighbor_label;
                  }
                  if(neighbor_label < min_label) min_label = neighbor_label;
                  if(neighbor_label > max_label) max_label = neighbor_label;
                  time++;
               }
            }
         }

         if(time == 0){
            (*(image->putpixel))(image, r, c, (double)nextlabel);

	    nextlabel++;
            if(nextlabel > maxval){
               fprintf(stderr, "Error in cc_label_4connect(). Nextlabel overflowed ");
               fprintf(stderr, "the size allowable by your datatype selection.\n");
               exit(1);
            }
         }
         else if(max_label == min_label){
            (*(image->putpixel))(image, r, c, (double)min_label);
         }
         else{

            (*(image->putpixel))(image, r, c, (double)min_label);

            for(rr=(pos-1);rr>=(pos-source_cols);rr--){

               n_x = rr % source_cols;
               n_y = rr / source_cols;

               if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

                  neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                     n_y, n_x));

                  if(neighbor_label == max_label){
                     (*(image->putpixel))(image, n_y, n_x, (double)min_label);
                  }
               }
            }
         }
      }
   }

   if(VERBOSE) printf("The maximum label is less than or equal to %lu\n", nextlabel-1);

   /****************************************************************************
   * Do the reverse pass of the connected component labeling algorithm.
   ****************************************************************************/
   for(r=(source_rows-1);r>=0;r--){
      for(c=(source_cols-1);c>=0;c--){

         pos = r * source_cols + c;

         value = (*(source_image->getpixel))(source_image, r, c);

         max_label = min_label = (unsigned long)((*(image->getpixel))(image, r, c));

	 for(n=0;n<4;n+=2){

            n_x = c + rsa_x[n];
            n_y = r + rsa_y[n];

            if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

               neighbor_value = (*(source_image->getpixel))(source_image,
                                 n_y, n_x);

               if(value == neighbor_value){

                  neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                   n_y, n_x));


                  if(neighbor_label < min_label) min_label = neighbor_label;
                  if(neighbor_label > max_label) max_label = neighbor_label;
               }
            }
         }

         if(max_label != min_label){

            (*(image->putpixel))(image, r, c, (double)min_label);

            for(rr=(pos+1);rr<=(pos+source_cols);rr++){

               n_x = rr % source_cols;
               n_y = rr / source_cols;

               if((n_x>=0)&&(n_y>=0)&&(n_x<source_cols)&&(n_y<source_rows)){

                  neighbor_label = (unsigned long)((*(image->getpixel))(image,
                                                   n_y, n_x));

                  if(neighbor_label == max_label){
                     (*(image->putpixel))(image, n_y, n_x, (double)min_label);
                  }
               }
            }
         }
      }
   }

   /****************************************************************************
   * Go through and relabel all of the components with values that are 
   * 0 -> (#unique_components - 1).
   ****************************************************************************/
   if((labels = (ULONG *) mycalloc("labels", nextlabel, sizeof(ULONG))) == NULL){
      fprintf(stderr, "Malloc error in cc_label_4connect().\n");
      exit(1);
   }

   for(r=0,pos=0;r<source_rows;r++){
      for(c=0;c<source_cols;c++,pos++){
         neighbor_label = (ULONG)((*(image->getpixel))(image, r, c));
         labels[neighbor_label]++;
      }
   }

   for(r=0,comp=0;r<nextlabel;r++){
      if(labels[r] != 0){
         labels[r] = comp;
         comp++;
      }
   }

   *numcomp = comp;
   if(VERBOSE) printf("There are actually %lu unique component labels in the image.\n", *numcomp);

   for(r=(source_rows-1);r>=0;r--){
      for(c=(source_cols-1);c>=0;c--){
         neighbor_label = (ULONG)((*(image->getpixel))(image, r, c));
         (*(image->putpixel))(image, r, c, (double)(labels[neighbor_label]));
      }
   }

   myfree("labels", labels);
}
