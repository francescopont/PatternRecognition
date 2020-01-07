/*******************************************************************************
* File: segment.c
* Overview: This code was written to segment the breast from a mammography
* image. Although it should work with different images from different scanners,
* it was developed using images from the Digital Database for Screening
* Mammography at the University of South Florida. The breast is detected and
* segmented and the boundary is chain-coded, smoothed and converted to a
* polygon. The code was written for images stored in a row-major format. The
* image can stored in the file with one or two bytes per pixel and can have
* header bytes. In the case of images with two bytes per pixel, the bytes
* can be swapped to allow for images written on computers with the opposite
* byte order (BIG-ENDIAN vs. SMALL-ENDIAN). An image cacheing mechanism is
* used to reduce the memory strain on a system because mammograms can be
* quite large.
* Note: This code was written as several separate files. It was then combined
* in this file to make the interface simpler. That is why there is somewhat
* more complexity in the code in this file than there needs to be.
* Name: Michael Heath, University of South Florida
* Date: 1/20/2000
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include "virtual_image.h"
#include "optical_density.h"
#include "breast.h"
#include "myalloc.h"
#include <string.h>

#define AIM_RESOLUTION 200.0

extern int VERBOSE;

void find_center(unsigned char *thresholded_image, int rows, int cols,
    int *center_row, int *center_col);
void estimate_tissue_image(unsigned char *thresholded_image, int rows, int cols,
   unsigned char **tissue_image, int center_row, int center_col, int basepos,
   float **rvals, float **cvals, int desired_points);
void bordercode(unsigned char *image, int rows, int cols, int **c_save,
   int row, int col, int *nn, int maxlen);
int find_threshold(CACHEIM *image);
int find_uthreshold(CACHEIM *image);

/*******************************************************************************
* Function: segment_breast
* Purpose: This function performs a segmentation of the breast and fills arrays
* of the coordinates of a polygon that outlines the segmented breast.
* Name: Michael Heath, University of South Florida
* Date: 7/22/99
*******************************************************************************/
void segment_breast(char *filename, int rows, int cols, int bpp, int hb,
       int swap, int projection, int apex, float resolution, int mode,
       double threshold_double, double uthreshold_double, int *numpoints,
       float **xcoord, float **ycoord, int mapmode, double A, double B)
{
   CACHEIM rawimage, aggregatedimage, *tempimageptr=NULL, imageptr;
   int aggfactor, r, c, center_row, center_col;
   unsigned char *thresholded_image=NULL, *tissue_image=NULL;
   int basepos;
   double od_offset, od_scale, min_density, max_density;
   int threshold, uthreshold;

   void aggregate_median(CACHEIM *image, CACHEIM *source_image, int aggfactor,
      int max_depth, int max_cache_size, int datatype, char *mode);

   memset(&rawimage, 0, sizeof(CACHEIM));
   memset(&aggregatedimage, 0, sizeof(CACHEIM));

   /****************************************************************************
   * Read in the image as a cached image.
   ****************************************************************************/
   if(bpp == 1){
      if(allocate_cached_image(&rawimage, rows, cols, UCHARNUM, filename,
         "rb", 12, (rows*cols*sizeof(UCHAR)/16), hb, (USHORT)swap,
         resolution) == 0) exit(1);
   }
   else if(bpp == 2){
      if(allocate_cached_image(&rawimage, rows, cols, USHORTNUM, filename,
         "rb", 12, (rows*cols*sizeof(USHORT)/16), hb, (USHORT)swap,
         resolution) == 0) exit(1);
   }
   else{
      fprintf(stderr, "Invalid bits per pixel (%d).\n", bpp);
      exit(1);
   }

   /****************************************************************************
   * Make a lower resolution copy of the image to process to find the breast
   * edge. This is done for computational speed. We could use the full size
   * image, but using a smaller resolution version is pretty much just as good,
   * and it will be processed faster. We reduce the resolution of the image
   * using a median/subsampling filter. The subsampling must be done by an
   * integer factor in each dimension using this approach.
   ****************************************************************************/
   if((rawimage.resolution == 0) || (rawimage.resolution >= 200.0)){
      aggfactor = 1;
      tempimageptr = &rawimage;
   }
   else{
      aggfactor = (int)ceil(AIM_RESOLUTION / rawimage.resolution);
      aggregate_median(&aggregatedimage, &rawimage, aggfactor, 12,
         (rawimage.rows * rawimage.cols * dt_size[USHORTNUM]) /
         (aggfactor*aggfactor*4), USHORTNUM, "rb");
      tempimageptr = &aggregatedimage;
   }

   if(VERBOSE)
      printf("The image is being aggregated by %d times for finding the breast boundary.\n", aggfactor);

   /****************************************************************************
   * If the user wants to work in scaled optical density space, map to optical
   * density.
   ****************************************************************************/
   if(mapmode == LINEARMAP){
      if(VERBOSE) printf("Mapping to optical density using linear mapping.\n");
      linear_optical_density(&imageptr, tempimageptr, 12,
         (rawimage.rows * rawimage.cols * dt_size[USHORTNUM]) / (aggfactor*aggfactor*4), "rb",
         &od_offset, &od_scale, &min_density, &max_density, A, B);
      threshold = od_to_scaled_od(threshold_double);
      uthreshold = od_to_scaled_od(uthreshold_double);
   }
   else if(mapmode == LOG10MAP){
      if(VERBOSE) printf("Mapping to optical density using log10 mapping.\n");
      log10_optical_density(&imageptr, tempimageptr, 12,
         (rawimage.rows * rawimage.cols * dt_size[USHORTNUM]) / (aggfactor*aggfactor*4), "rb",
         &od_offset, &od_scale, &min_density, &max_density, A, B);
      threshold = od_to_scaled_od(threshold_double);
      uthreshold = od_to_scaled_od(uthreshold_double);
   }
   else{
      if(VERBOSE) printf("Processing in gray level mode.\n");
      imageptr = *tempimageptr;
      threshold = (unsigned short int)floor(threshold_double);
      uthreshold = (unsigned short int)floor(uthreshold_double);
   }

   /****************************************************************************
   * Create a thresholded version of the reduced resolution image. All pixels
   * with a value less than the threshold value are set to 0. All pixels above
   * the upper threshold are set to 0 that 1) can be reached from the left or
   * right edge of the image by travelling along a row of pixels above the
   * threshold or 2) can be reached from the top or bottom of the image by
   * travelling along a column of pixels above the threshold. The rest of the
   * pixels are set to 255 to indicate that they might be in the breast.
   ****************************************************************************/
   thresholded_image = (unsigned char *) mycalloc("tmp_thresholded_image",
      imageptr.rows*imageptr.cols, sizeof(unsigned char));

   if(mode == AUTOMODE){
      threshold = find_threshold(&imageptr);
      if(VERBOSE) printf("The automatically found threshold value is %d.\n", threshold);
      uthreshold = find_uthreshold(&imageptr);
      if(VERBOSE) printf("The automatically found upper threshold value is %d.\n", uthreshold);
   }
   if((mode == MANUALMODE) && (threshold == 0.0)){
      threshold = find_threshold(&imageptr);
      if(VERBOSE) printf("The automatically found threshold value is %d.\n", threshold);
   }
   if((mode == MANUALMODE) && (uthreshold == 0.0)){
      uthreshold = find_uthreshold(&imageptr);
      if(VERBOSE) printf("The automatically found upper threshold value is %d.\n", uthreshold);
   }

   if(VERBOSE) printf("The threshold pixel value is %d\n", threshold);
   if(VERBOSE) printf("The upper threshold pixel value is %d\n", uthreshold);

   for(r=0;r<imageptr.rows;r++){

      for(c=0;c<imageptr.cols;c++){
         if(((*(imageptr.getpixel))(&imageptr, r, c)) < (double)threshold)
            thresholded_image[r*imageptr.cols+c] = 0;
         else thresholded_image[r*imageptr.cols+c] = 255;
      }
      c = 0;
      while((c<imageptr.cols) && (((*(imageptr.getpixel))(&imageptr, r, c)) > (double)uthreshold)){
         thresholded_image[r*imageptr.cols+c] = 0;
         c++;
      }
      c = imageptr.cols-1;
      while((c>=0) && (((*(imageptr.getpixel))(&imageptr, r, c)) > (double)uthreshold)){
         thresholded_image[r*imageptr.cols+c] = 0;
         c--;
      }
   }

   for(c=0;c<imageptr.cols;c++){
      r = 0;
      while((r<imageptr.rows) && (((*(imageptr.getpixel))(&imageptr, r, c)) > (double)uthreshold)){
         thresholded_image[r*imageptr.cols+c] = 0;
         r++;
      }
      r = imageptr.rows-1;
      while((r>=0) && (((*(imageptr.getpixel))(&imageptr, r, c)) > (double)uthreshold)){
         thresholded_image[r*imageptr.cols+c] = 0;
         r--;
      }
   }

   /****************************************************************************
   * Find the center of the largest blob in the image (of pixels above the
   * threshold determined above). This should definitly be inside of the breast
   * tissue. I am assuming that any labels and or border in the image is
   * smaller than the breast. This process basically finds the blob that
   * has the largest circle that can be inscribed in the blob that will not
   * touch a background pixel. This is a computationally poor implementation
   * of a method to do this. I am sure it could be made faster.
   ****************************************************************************/
   find_center(thresholded_image, imageptr.rows, imageptr.cols, &center_row, &center_col);

   /****************************************************************************
   * From the thresholded image, estimate the pixels that belong to breast
   * tissue from the thresholded image. Use the information about the
   * orientation of the breast in the image. The function estimate_tissue_image()
   * uses the thresholded image, the knowledge of which side of the image has
   * the base of the breast along it and a pixel position believed to be inside
   * the breast tissue to create a polygon of row and column coordinates of the
   * smoothed boundary of the breast. The breast border polygon is formed with
   * 100 vertices.
   ****************************************************************************/
   if(apex == APEXLEFT) basepos = 1; /* Indicates the base of the breast is on the right side of the image. */
   else basepos = 0;                 /* Indicates the base of the breast is on the left side of the image. */

   estimate_tissue_image(thresholded_image, imageptr.rows, imageptr.cols,
      &tissue_image, center_row, center_col, basepos, ycoord, xcoord, 100);
   *numpoints = 100;

   if(0){
      FILE *fp=NULL;
      char tissue_filename[200];

      sprintf(tissue_filename, "%s.tissue.pgm", filename);
      fp = fopen(tissue_filename, "wb");
      fprintf(fp, "P5\n%d %d\n255\n", imageptr.cols, imageptr.rows);
      fwrite(tissue_image, 1, imageptr.rows*imageptr.cols, fp);
      fclose(fp);
   }

   for(r=0;r<(*numpoints);r++){
      (*xcoord)[r] *= aggfactor;
      (*ycoord)[r] *= aggfactor;
   }
}

/*******************************************************************************
* Function: find_center
* Purpose: We believe that the largest region in a thresholded mammogram
* contains the breast tissue and that the center of this region will be a point
* inside the breast tissue. Calculating the central pixel is done to avoid
* bright regions that may be connected to the breast. We assume that any such
* regions are smaller in size than the breast tissue. A sort of medial axis
* transform is used to find the center of the breast. The code here could
* probably be rewritten to speed up this process a lot. Currently I just do
* repetative 3x3 erosions until no pixels are left in the binary image.
* Name: Michael Heath, University of South Florida
* Date: 5/23/98
*******************************************************************************/
void find_center(unsigned char *thresholded_image, int rows, int cols,
   int *center_row, int *center_col)
{
   unsigned short int *erodeim=NULL, *last_erodeim=NULL, *hold_ptr=NULL;
   int small_rows, small_cols;
   int r, c, rr, cc, pos, checkpos, k, itteration;
   int neighboroffset[8], sum, numchanged;

   /****************************************************************************
   * Reduce the resolution by 2 before finding the center. This is just done
   * for speed.
   ****************************************************************************/
   small_rows = rows/2;
   small_cols = cols/2;

   if((last_erodeim = (unsigned short int *) mycalloc("last_erodeim",
         small_rows*small_cols, sizeof(unsigned short int)))==NULL){
      fprintf(stderr, "Malloc error in find_center().\n");
      exit(1);
   }

   if((erodeim = (unsigned short int *) mycalloc("erodeim",
         small_rows*small_cols, sizeof(unsigned short int))) == NULL){
      fprintf(stderr, "Malloc error in find_center().\n");
      exit(1);
   }

   /****************************************************************************
   * Threshold the image.
   ****************************************************************************/
   for(r=0,pos=0;r<small_rows;r++){
      for(c=0;c<small_cols;c++,pos++){
         for(rr=(2*r);rr<(2*r+2);rr++){
            for(cc=(2*c);cc<(2*c+2);cc++){
               if(thresholded_image[rr*cols+cc] != 0){
                  last_erodeim[pos] = 1;
                  cc = 2*c+2;
                  rr = 2*r+2;
               }
            }
         }
      }
   }

   rows = small_rows;
   cols = small_cols;

   neighboroffset[0] = -cols - 1;
   neighboroffset[1] = -cols;
   neighboroffset[2] = -cols + 1;
   neighboroffset[3] = -1;
   neighboroffset[4] = +1;
   neighboroffset[5] = +cols - 1;
   neighboroffset[6] = +cols;
   neighboroffset[7] = +cols + 1;

   numchanged = 1;
   itteration = 1;
   while(numchanged != 0){

      numchanged = 0;

      memset(erodeim, 0, (size_t)(rows*cols*2));

      for(r=1;r<(rows-1);r++){
         for(c=1;c<(cols-1);c++){

            pos = r * cols + c;

            if(last_erodeim[pos] == itteration){

               sum = 0;
               for(k=0;k<8;k++){
                  checkpos = pos + neighboroffset[k];
                  if(last_erodeim[checkpos] == itteration) sum++;
               }

               if(sum == 0) erodeim[pos] = last_erodeim[pos];
               else if(sum < 8){
                  erodeim[pos] = last_erodeim[pos];
               }
               else if(sum == 8){
                  erodeim[pos] = last_erodeim[pos] + 1;
                  numchanged++;
               }
            }
            else erodeim[pos] = last_erodeim[pos];
         }
      } 

      hold_ptr = last_erodeim;
      last_erodeim = erodeim;
      erodeim = hold_ptr;

      itteration++;
   }

   checkpos = 0;
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         if(last_erodeim[pos] > last_erodeim[checkpos]) checkpos = pos;
      }
   }

   *center_row = checkpos / cols;
   *center_col = checkpos % cols;

   /****************************************************************************
   * Since we reduce the resolution by two, double the position of the center.
   ****************************************************************************/
   *center_row *= 2;
   *center_col *= 2;

   myfree("erodeim", erodeim);
   myfree("last_erodeim", last_erodeim);
   erodeim = NULL;
   last_erodeim = NULL;
   hold_ptr = NULL;
}

/*******************************************************************************
* Function: estimate_tissue_image
* Purpose: This function uses the thresholded image and a pixel location that is
* believed to be inside of the breast tissue to estimate which pixels contain
* breast tissue. This function is necessary, because there can be many pixels
* above the threshold that are not part of the breast. Connectivity is used to
* find the pixels that are part of the breast. A tissue image is allocated
* which identifies pixels that are believed to fall inside the breast region.
* A polygon is also constructed which is believed to contain the breast tissue.
* Name: Michael Heath, University of South Florida
* Date: 8/19/98
*******************************************************************************/
void estimate_tissue_image(unsigned char *thresholded_image, int rows, int cols,
   unsigned char **tissue_image, int center_row, int center_col, int basepos,
   float **rvals, float **cvals, int desired_points)
{
   unsigned char *tissue=NULL;
   unsigned char *tissue_edge=NULL;
   int i, j, r, c, pos, best_col, best_row, changed, crossed, chain_length, count=0;
   int *col_limit = NULL;
   int *chain_code = NULL;
   float thisdistance;
   float *sdat=NULL, *xdat=NULL, *ydat=NULL;
   float *sdat_tmp=NULL, *xdat_tmp=NULL, *ydat_tmp=NULL;
   float sumf, sumx, sumy;
   double *kernel = NULL;
   int ksize;
   double v0x, v0y, v1x, v1y, cross_product;

   void fill_in_holes(unsigned char *image, int rows, int cols, float processing_resolution);

   /****************************************************************************
   * Allocate the tissue image.
   ****************************************************************************/
   if((tissue = (unsigned char *) mycalloc("tissue", rows*cols,
         sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Malloc error in estimate_tissue_image().\n");
      exit(1);
   }
   *tissue_image = tissue;

   for(i=0;i<(rows*cols);i++) tissue[i] = thresholded_image[i];

   /****************************************************************************
   * Grow a region from the "center".
   ****************************************************************************/
   tissue[center_row * cols + center_col] = 254;

   for(r=center_row;r>=0;r--){
      for(c=center_col;c>=0;c--){
         pos = r * cols + c;
         if(tissue[pos] == 254){
            if(((c-1) >= 0) && (r >= 0) && (tissue[pos-1] == 255))
               tissue[pos-1] = 254;
            if(((c-1) >= 0) && ((r-1) >= 0) && (tissue[pos-cols-1] == 255))
               tissue[pos-cols-1] = 254;
            if((c >= 0) && ((r-1) >= 0) && (tissue[pos-cols] == 255))
               tissue[pos-cols] = 254;
         }
      }
      for(c=center_col;c<cols;c++){
         pos = r * cols + c;
         if(tissue[pos] == 254){
            if(((c+1) < cols) && (r >= 0) && (tissue[pos+1] == 255))
               tissue[pos+1] = 254;
            if(((c+1) < cols) && ((r-1) >= 0) && (tissue[pos-cols+1] == 255))
               tissue[pos-cols+1] = 254;
            if((c >= 0) && ((r-1) >= 0) && (tissue[pos-cols] == 255))
               tissue[pos-cols] = 254;
         }
      }
   }  

   for(r=center_row;r<rows;r++){
      for(c=center_col;c>=0;c--){
         pos = r * cols + c;
         if(tissue[pos] == 254){
            if(((c-1) >= 0) && (r < rows) && (tissue[pos-1] == 255))
               tissue[pos-1] = 254;
            if(((c-1) >= 0) && ((r+1) < rows) && (tissue[pos+cols-1] == 255))
               tissue[pos+cols-1] = 254;
            if((c >= 0) && ((r+1) < rows) && (tissue[pos+cols] == 255))
               tissue[pos+cols] = 254;
         }
      }
      for(c=center_col;c<cols;c++){
         pos = r * cols + c;
         if(tissue[pos] == 254){
            if(((c+1) < cols) && (r < rows) && (tissue[pos+1] == 255))
               tissue[pos+1] = 254;
            if(((c+1) < cols) && ((r+1) < rows) && (tissue[pos+cols+1] == 255))
               tissue[pos+cols+1] = 254;
            if((c >= 0) && ((r+1) < rows) && (tissue[pos+cols] == 255))
               tissue[pos+cols] = 254;
         }
      }
   }  

   /****************************************************************************
   * Remove any lines that can make it all the way across the image. Candidate
   * lines are those on the top and the bottom of the image, or those that touch
   * lines that cross the image on the top of the bottom of the image.
   ****************************************************************************/
   r = 0;
   do{
      c = 0;
      while((c < cols) && (tissue[r*cols+c] >= 254)) c++;
      if(c == cols){
         crossed = 1;
         for(c=0;c<cols;c++) tissue[r*cols+c] = 0;
      }
      else crossed = 0;
      r++;
   }while((r < rows) && (crossed == 1));

   r = rows - 1;
   do{
      c = 0;
      while((c < cols) && (tissue[r*cols+c] >= 254)) c++;
      if(c == cols){
         crossed = 1;
         for(c=0;c<cols;c++) tissue[r*cols+c] = 0;
      }
      else crossed = 0;
      r--;
   }while((r >= 0) && (crossed == 1));
 
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         if(tissue[pos] == 254) tissue[pos] = 255;
         else tissue[pos] = 0;
      }
   }

   fill_in_holes(tissue, rows, cols, AIM_RESOLUTION);  /* The value of AIM_RESOLUTION is irrelevent here. */

   /****************************************************************************
   * Create an image that will have all pixels inside the breast tissue, that
   * share a non breast tissue pixel in their immediate neighborhood.
   ****************************************************************************/
   if((tissue_edge = (unsigned char *) mycalloc("tissue_edge", rows*cols,
         sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Malloc error in estimate_tissue_image().\n");
      exit(1);
   }
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         if(tissue[pos] == 255){
            if(((r-1)<0) || ((c-1)<0) || ((r+1)==rows) || ((c+1)==cols)){
               tissue_edge[pos] = 255;
            }
            else{
               if(tissue[pos-cols-1] != 255) tissue_edge[pos] = 255;
               else if(tissue[pos-cols] != 255) tissue_edge[pos] = 255;
               else if(tissue[pos-cols+1] != 255) tissue_edge[pos] = 255;
               else if(tissue[pos-1] != 255) tissue_edge[pos] = 255;
               else if(tissue[pos+1] != 255) tissue_edge[pos] = 255;
               else if(tissue[pos+cols-1] != 255) tissue_edge[pos] = 255;
               else if(tissue[pos+cols] != 255) tissue_edge[pos] = 255;
               else if(tissue[pos+cols+1] != 255) tissue_edge[pos] = 255;
            }
         }
      }
   }

   /****************************************************************************
   * Form an array indicating the farthest column from the base of the breast
   * for each row of the image.
   ****************************************************************************/
   if((col_limit = (int *) mycalloc("col_limit", rows, sizeof(int))) == NULL){
      fprintf(stderr, "Calloc error in estimate_tissue_image().\n");
      exit(1);
   }
   if(basepos == 0){
      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++) if(tissue_edge[pos] == 255) col_limit[r] = c;
      }
   }
   else{
      for(r=0;r<rows;r++){
         for(c=cols-1,pos=((r+1)*cols-1);c>=0;c--,pos--)
            if(tissue_edge[pos] == 255) col_limit[r] = c;
      }
   }

   /****************************************************************************
   * Find the row in which the breast is the widest near the center of the
   * breast.
   ****************************************************************************/
   best_row = center_row;
   if(basepos == 0) best_col = 0;
   else best_col = cols-1;
   do{
      pos = 0;
      changed = 0;
      for(r=(center_row-10);r<=(center_row+10);r++){
         if((r>=0) && (r < rows)){
            pos += col_limit[r];
            count++;
         }
      }
      c = pos / count;
      if((basepos == 0) && (c > best_col)){
         best_col = c;
         best_row = r;
         changed = 1;
      } 
      else if((basepos == 1) && (c < best_col)){
         best_col = c;
         best_row = r;
         changed = 1;
      } 
   }while(changed == 1);

   best_col = col_limit[best_row];

   if(tissue_edge[best_row*cols+best_col] == 0){
      printf("Error the starting chain pixel has a value of %d\n", 
	 (int)tissue_edge[best_row*cols+best_col]);
   }

   /****************************************************************************
   * Chaincode the boundary of the breast.
   ****************************************************************************/
   bordercode(tissue, rows, cols, &chain_code, best_row, best_col,
      &chain_length, 50000);

   /****************************************************************************
   * Place the chain coded curve in arrays that express x and y as a function
   * of the path length.
   ****************************************************************************/
   sdat_tmp = (float *) mycalloc("sdat_tmp", chain_length, sizeof(float));
   xdat_tmp = (float *) mycalloc("xdat_tmp", chain_length, sizeof(float));
   ydat_tmp = (float *) mycalloc("ydat_tmp", chain_length, sizeof(float));

   if((sdat_tmp==NULL) || (xdat_tmp==NULL) || (ydat_tmp==NULL)){
      fprintf(stderr, "Malloc error of sdat_tmp, xdat_tmp or ydat_tmp.\n");
      exit(1);
   }

   r = best_row;
   c = best_col;

   sdat_tmp[0] = 0.0;
   xdat_tmp[0] = (float)c;
   ydat_tmp[0] = (float)r;

   for(i=0;i<(chain_length-1);i++){
      switch(chain_code[i]){
         case 0: r +=  0; c +=  1; thisdistance = 1.0;       break;
         case 1: r += -1; c +=  1; thisdistance = sqrt(2.0); break;
         case 2: r += -1; c +=  0; thisdistance = 1.0;       break;
         case 3: r += -1; c += -1; thisdistance = sqrt(2.0); break;
         case 4: r +=  0; c += -1; thisdistance = 1.0;       break;
         case 5: r +=  1; c += -1; thisdistance = sqrt(2.0); break;
         case 6: r +=  1; c +=  0; thisdistance = 1.0;       break;
         case 7: r +=  1; c +=  1; thisdistance = sqrt(2.0); break;
      }
      if((r>=0) && (c>=0) && (c<cols) && (r<rows)){
         sdat_tmp[i+1] = sdat_tmp[i] + thisdistance;
         xdat_tmp[i+1] = (float)c;
         ydat_tmp[i+1] = (float)r;
      }
   }

   /****************************************************************************
   * Smooth along the curve.
   ****************************************************************************/
   sdat = (float *) mycalloc("sdat", chain_length, sizeof(float));
   xdat = (float *) mycalloc("xdat", chain_length, sizeof(float));
   ydat = (float *) mycalloc("ydat", chain_length, sizeof(float));

   if((sdat==NULL) || (xdat==NULL) || (ydat==NULL)){
      fprintf(stderr, "Malloc error of sdat, xdat or ydat.\n");
      exit(1);
   }

   ksize = chain_length / 20;
   if((ksize % 2) == 0) ksize += 1;

   if((kernel = (double *) mycalloc("kernel", ksize, sizeof(double))) == NULL){
      fprintf(stderr, "Malloc error of kernel.\n");
      exit(1);
   }

   for(i=0;i<ksize;i++)
      kernel[i] = exp(-0.5*((double)i-(ksize/2))*((double)i-(ksize/2)) /
                  ((double)ksize/6.0 * (double)ksize/6.0));

   for(i=0;i<chain_length;i++){
      sumx = sumy = sumf = 0.0;
      for(j=0;j<ksize;j++){
         pos = (chain_length + (i + j - ksize/2)) % (chain_length);
         sumx += kernel[j] * xdat_tmp[pos];
         sumy += kernel[j] * ydat_tmp[pos];
         sumf += kernel[j];
      }
      xdat[i] = sumx / sumf;
      ydat[i] = sumy / sumf;
   }
   myfree("kernel", kernel);

   /****************************************************************************
   * Copy points along the boundary into arrays to pass back to the calling
   * function. Use the two points on the boundary to determine whether
   * the border is traversed clockwise or counterclockwise. Note if the cross
   * product is zero, we have an unhandled special case.
   ****************************************************************************/
   v0x = xdat[0] - center_col;
   v0y = ydat[0] - center_row;
   v1x = xdat[chain_length/10] - center_col;
   v1y = ydat[chain_length/10] - center_row;

   cross_product = v0x * v1y - v1x * v0y;

   if(((*rvals) = (float *) mycalloc("rvals", desired_points,
         sizeof(float))) == NULL){
      fprintf(stderr, "Malloc error of rvals.\n");
      exit(1);
   }
   if(((*cvals) = (float *) mycalloc("cvals", desired_points,
         sizeof(float))) == NULL){
      fprintf(stderr, "Malloc error of cvals.\n");
      exit(1);
   }

   if(((basepos == 0) && (cross_product < 0)) ||
      ((basepos == 1) && (cross_product > 0))){

      for(i=0;i<desired_points;i++){
         pos = (int)floor(((double)i /((double)desired_points - 0.0001)) *
                  chain_length);
         (*cvals)[i] = xdat[pos];
         (*rvals)[i] = ydat[pos];
      }
   }
   else{

      for(i=0;i<desired_points;i++){
         pos = (int)floor(((double)i /((double)desired_points - 0.0001)) *
                  chain_length);
         (*cvals)[desired_points-1-i] = xdat[pos];
         (*rvals)[desired_points-1-i] = ydat[pos];
      }
   }

   myfree("sdat_tmp", sdat_tmp);
   myfree("xdat_tmp", xdat_tmp);
   myfree("ydat_tmp", ydat_tmp);

   myfree("chain_code", chain_code);
   myfree("tissue_edge", tissue_edge);
   myfree("col_limit", col_limit);
   myfree("sdat", sdat);
   myfree("xdat", xdat);
   myfree("ydat", ydat);
}

/*******************************************************************************
* Function: bordercode
* Purpose: To chain code the border of a block of pixels that have the same
* value. The user should pass in the row and column position on the border
* (just inside) of the region they want to chain code the border.
* Name: Michael Heath, University of South Florida
* Date: 8/27/98
*******************************************************************************/
void bordercode(unsigned char *image, int rows, int cols, int **c_save,
   int row, int col, int *nn, int maxlen)
{
   int d, r, c, last_d, n=0, hit, k;
   int dx[8] = { 1, 1, 0,-1,-1,-1, 0, 1};
   int dy[8] = { 0,-1,-1,-1, 0, 1, 1, 1};
   int last_r, last_c, loop_time;
   int *chain = NULL;
   unsigned char **image2d = NULL, val;

   /****************************************************************************
   * Make a temporary 2d image by allocating an array of pointers and setting
   * them to the rows of the image.
   ****************************************************************************/
   if((image2d = (unsigned char **) mycalloc("image2d", rows,
         sizeof(unsigned char *))) == NULL){
      fprintf(stderr, "Error allocating an array in bordercode().\n");
      exit(1);
   }
   for(r=0;r<rows;r++) image2d[r] = image + r * cols;

   /****************************************************************************
   * Find a direction to point in a direction off of the blob.
   ****************************************************************************/
   val = image2d[row][col];
   d = 0;
   while((d <= 7) && ((r=row+dy[d]) >= 0) && ((c=col+dx[d]) >= 0) && (r<rows) &&
      (c<cols) && (image2d[r][c]==val)) d++;

   if(d == 8){
      printf("Point specified in bordercode is not on the edge of a region.\n");
      return;
   }

   /****************************************************************************
   * Allocate a temporary array to use to store the chain code with its maximum
   * possible length.
   ****************************************************************************/
   if((chain = (int *) mycalloc("chain", maxlen, sizeof(int))) == NULL){
      fprintf(stderr, "Malloc error in bordercode().\n");
      exit(1);
   }

   /****************************************************************************
   * Chaincode the border. We chain until we reach the initial point plus go
   * one more step to see that we not only approached the start point of the
   * chain, but also that we are chaining in the same direction.
   ****************************************************************************/
   last_d = d;
   r = row;
   c = col;
   loop_time = 0;
   last_r = -1;
   last_c = -1;
   do{

      hit = 0;

      if(loop_time != 0){
         last_r = r;
         last_c = c;
      }

      for(k=last_d+1;k<last_d+8;k++){
         d = k % 8;
         if(((r+dy[d]) >= 0) && ((r+dy[d]) < rows) && ((c+dx[d]) >= 0) &&
               ((c+dx[d]) < cols)){
            if(image2d[r+dy[d]][c+dx[d]] == val){
               hit = 1;
               r += dy[d];
               c += dx[d];
               break;
            }
         }
      }

      if(hit == 1){
         chain[n++] = d;
         last_d = (d + 5) % 8;
      }
      if(n == maxlen) break;

      loop_time++;

   }while(!((last_r == row) && (last_c == col) && (row+dy[chain[0]] == r) && (col+dx[chain[0]] == c)));

   n -= 1;
   /* if(VERBOSE) printf("The bordercode has length %d. [%d,%d] [%d,%d].\n", n, row, col, last_r, last_c); */
   *nn = n;

   /****************************************************************************
   * Allocate a new array to store the chain code with the exact length that
   * we need.
   ****************************************************************************/
   if(((*c_save) = (int *) mycalloc("c_save", n, sizeof(int))) == NULL){
      fprintf(stderr, "Malloc error in bordercode().\n");
      exit(1);
   }

   /****************************************************************************
   * Copy the array into the one we will save.
   ****************************************************************************/
   for(r=0;r<n;r++) (*c_save)[r] = chain[r];

   myfree("chain", chain);
   myfree("image2d", image2d);
}

/*******************************************************************************
* Function: find_threshold
* Purpose: This function finds the first local minimum in the median filtered
* histogram.
* Name: Michael Heath, University of South Florida
* Date: 4/30/98
*******************************************************************************/
int find_threshold(CACHEIM *image)
{
   int r, c, mode, p, q, k, divnumber;
   unsigned long int hist[65536], filteredhist[65536], m[7];
   unsigned short int max, value;

   /****************************************************************************
   * Calculate the histogram.
   ****************************************************************************/
   memset(hist, 0, 65536*sizeof(unsigned long int));

   for(r=0;r<image->rows;r++){
      for(c=0;c<image->cols;c++){
         value = (unsigned short int)(*(image->getpixel))(image, r, c);
         if(value != 0){   /* Skip because of digitally zeroed patient information. */
            hist[value]++;
         }
      }
   }

   /****************************************************************************
   * Find the maximum value in the image by examining the histogram.
   ****************************************************************************/
   max = 0;
   for(r=65535;r>=0;r--){
      if(hist[r] != 0){
         max = (unsigned short int)r;
         break;
      }
   }

   /****************************************************************************
   * Re-bin the histogram by summing up the counts in adjacent bins.
   ****************************************************************************/
   if(max > 256){
      divnumber = max / 256;

      /* if(VERBOSE) printf("Dividing by %d\n", divnumber); */

      for(r=0;r<65536;r++){
         c = r/divnumber;
         hist[c] += hist[r];
      }
      for(r=(c+1);r<65536;r++) hist[r] = 0;
   }
   else divnumber = 1;

   /****************************************************************************
   * Find the mode in the histogram. This should be a value in the background
   * of the mammogram.
   ****************************************************************************/
   mode = 0;
   for(r=0;r<65536;r++) if(hist[mode] < hist[r]) mode = r;

   /****************************************************************************
   * Run a seven element median filter over the histogram.
   ****************************************************************************/
   for(r=(0+3);r<(65536-3);r++){
      for(c=(r-3),k=0;c<=(r+3);c++,k++) m[k] = hist[c];
      for(p=0;p<(7-1);p++){
         for(q=(p+1);q<7;q++){
            if(m[p] > m[q]){
               value = m[p];
               m[p] = m[q];
               m[q] = value;
            }
         }
      }
      filteredhist[r] = m[3];
   }

   for(r=mode;r<(65536-1);r++){
      c = mode;
      while((c>=0) && (filteredhist[c] == filteredhist[r])) c--;
      if((c>=0) && (filteredhist[c] > filteredhist[r])){
         c = mode;
         while((c<65536) && (filteredhist[c] == filteredhist[r])) c++;
            if((c<65536) && (filteredhist[c] > filteredhist[r])){

               /* if((filteredhist[r] < filteredhist[r-1]) && (filteredhist[r] < filteredhist[r+1])) */
            return(r*divnumber);
         }
      }
   }

   /* if(VERBOSE) printf("Using the mode for the threshold.\n"); */
   return(mode*divnumber);
}

/*******************************************************************************
* Function: find_uthreshold
* Purpose: This function finds the first upper local minimum in the median filtered
* histogram.
* Name: Michael Heath, University of South Florida
* Date: 4/30/98
*******************************************************************************/
int find_uthreshold(CACHEIM *image)
{
   int r, c, p, q, k, divnumber;
   unsigned long int hist[65536], filteredhist[65536], m[7];
   unsigned short int max, value;

   /****************************************************************************
   * Calculate the histogram.
   ****************************************************************************/
   memset(hist, 0, 65536*sizeof(unsigned long int));

   for(r=0;r<image->rows;r++){
      for(c=0;c<image->cols;c++){
         value = (unsigned short int)(*(image->getpixel))(image, r, c);
         if(value != 0){   /* Skip because of digitally zeroed patient information. */
            hist[value]++;
         }
      }
   }

   /****************************************************************************
   * Find the maximum value in the image by examining the histogram.
   ****************************************************************************/
   max = 0;
   for(r=65535;r>=0;r--){
      if(hist[r] != 0){
         max = (unsigned short int)r;
         break;
      }
   }

   /****************************************************************************
   * Re-bin the histogram by summing up the counts in adjacent bins.
   ****************************************************************************/
   if(max > 256){
      divnumber = max / 256;

      for(r=0;r<65536;r++){
         c = r/divnumber;
         hist[c] += hist[r];
      }
      for(r=(c+1);r<65536;r++) hist[r] = 0;
   }
   else divnumber = 1;

   /****************************************************************************
   * Run a seven element median filter over the histogram.
   ****************************************************************************/
   for(r=(0+3);r<(65536-3);r++){
      for(c=(r-3),k=0;c<=(r+3);c++,k++) m[k] = hist[c];
      for(p=0;p<(7-1);p++){
         for(q=(p+1);q<7;q++){
            if(m[p] > m[q]){
               value = m[p];
               m[p] = m[q];
               m[q] = value;
            }
         }
      }
      filteredhist[r] = m[3];
   }

   r = 65535-3;
   while((r>0) && (filteredhist[r-1] >= filteredhist[r])) r--;
   return(r*divnumber);
}

/*******************************************************************************
* Function: fill_in_holes
* Purpose: In an image that has the value 0 and 255, this function finds all
* groups of zero valued pixels that are surrounded by 255 valued pixels and
* sets those zero valued pixels to 255.
*******************************************************************************/
void fill_in_holes(unsigned char *image, int rows, int cols, float processing_resolution)
{
   CACHEIM im, ccimage;
   int r, c, pos;
   unsigned long int numcomp=0, compid;
   unsigned char *neighbors=NULL;
   const unsigned char IMAGEEDGE=2;

   void cc_label_4connect(CACHEIM *image, CACHEIM *source_image, char *filename,
      int max_depth, int max_cache_size, int datatype, char *mode, ULONG *numcomp);

   memset(&im, 0, sizeof(CACHEIM));
   memset(&ccimage, 0, sizeof(CACHEIM));

   /****************************************************************************
   * Copy an image into a CACHEIM structure. We have a function to compute
   * the connected components of an image in this stired in this form.
   ****************************************************************************/
   if(allocate_cached_image(&im, rows, cols, UCHARNUM, NULL, "w+b", 12,
      rows*cols*sizeof(UCHAR), 0, 0, processing_resolution) == 0) exit(1);

   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         if(image[pos] == 0) (*(im.putpixel))(&im, r, c, 0.0);
         else (*(im.putpixel))(&im, r, c, 1.0);
      }
   }

   /****************************************************************************
   * Label the connected components in the image. Components are 8 connected.
   ****************************************************************************/
   cc_label_4connect(&ccimage, &im, NULL, 12, rows*cols*sizeof(ULONG), ULONGNUM, "rb", &numcomp);

   /* if(VERBOSE) printf("There are %lu connected components in the image.\n", numcomp); */

   /****************************************************************************
   * Find any component that has (intensity) value 0 that does not touch the
   * edge of the image.
   ****************************************************************************/
   neighbors = (UCHAR *) calloc(numcomp, sizeof(UCHAR));

   for(r=0;r<rows;r++){
      if(image[r*cols+0] == 0){
         compid = (*(ccimage.getpixel))(&ccimage, r, 0);
         neighbors[compid] = IMAGEEDGE;
      }
      if(image[r*cols+(cols-1)] == 0){
         compid = (*(ccimage.getpixel))(&ccimage, r, (cols-1));
         neighbors[compid] = IMAGEEDGE;
      }
   }

   for(c=0;c<cols;c++){
      if(image[c] == 0){
         compid = (*(ccimage.getpixel))(&ccimage, 0, c);
         neighbors[compid] = IMAGEEDGE;
      }
      if(image[(rows-1)*cols+c] == 0){
         compid = (*(ccimage.getpixel))(&ccimage, (rows-1), c);
         neighbors[compid] = IMAGEEDGE;
      }
   }

   /****************************************************************************
   * Relabel all pixels of any component that has (intensity) value 0 that does
   * not touch the edge of the image.
   ****************************************************************************/
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         if(image[pos] == 0){
            compid = (*(ccimage.getpixel))(&ccimage, r, c);
            if(neighbors[compid] != IMAGEEDGE) image[pos] = 255;
         }
      }
   }

   free(neighbors);
   deallocate_cached_image(&im);
   deallocate_cached_image(&ccimage);
}
