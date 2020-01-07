/*******************************************************************************
* PROGRAM: optical_density.c
* PURPOSE: This file contains the code for converting an image to an unsigned
* short integer image that is linearly related to optical density.
* NAME: Michael Heath, University of South Florida
* DATE: 11/9/97
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include "virtual_image.h"
#include "optical_density.h"
#include "myalloc.h"
#include <string.h>

/*******************************************************************************
* Procedure: get_optical_density_segment
* Purpose: This procedure gets a segment from an image and converts it to
* to optical density using the precalculated lut.
* NAME: Michael Heath, University of South Florida
* Date: 11/8/97
*******************************************************************************/
void get_optical_density_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg)
{
   int r=0, c=0, start_row=0, start_col=0, num_rows=0, num_cols=0;
   USHORT hold_cachemode, *map=NULL;
   CACHEIM *sourceimage=NULL;
   double pixel_value;
   OD_REMAPPER_DATA *od_remapper_data=NULL;

   /****************************************************************************
   * Determine the portion of the image that this segment falls within. This
   * area contains the pixels we will be calculating.
   ****************************************************************************/
   image->cachestructure[rseg][cseg][0].rpage = (USHORT)rpage;
   image->cachestructure[rseg][cseg][0].cpage = (USHORT)cpage;
   image->cachestructure[rseg][cseg][0].dirtybit = CLEAN_DATA;

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

   sourceimage = (((CACHEIM **)image->sourceimlist)[0]);

   /****************************************************************************
   * Get the remapping array.
   ****************************************************************************/
   if(image->calculationdata == NULL){
      printf("Bad calculation data!!!\n");
      exit(1);
   }

   od_remapper_data = (OD_REMAPPER_DATA *)image->calculationdata;
   map = od_remapper_data->map;

   /****************************************************************************
   * Do the transformation for this segment.
   ****************************************************************************/
   for(r=start_row;r<(start_row+num_rows);r++){
      for(c=start_col;c<(start_col+num_cols);c++){

         pixel_value = (*(sourceimage->getpixel))(sourceimage, r, c);
         pixel_value = (double)map[(int)pixel_value];

         (*(image->putpixel))(image, r, c, pixel_value);
      }
   }

   /****************************************************************************
   * Replace the cachemode with its previous value now that we are done updating
   * the cache.
   ****************************************************************************/
   image->cachemode = hold_cachemode;
}

double roundval(double x)
{
   if(x >= 0) x = floor(x + 0.5);
   else x = ceil(x - 0.5);

   return(x);
}

/*******************************************************************************
* Procedure: linear_optical_density
* Purpose: This function maps the image to optical density. It uses the equation
* optical_density = A + (B * gray_level). The real valued optical densities are
* then converted to integer values by scaling optical densities in the range
* MIN_DENSITY to MAX_DENSITY to fall in the range 65535->0.
* NAME: Michael Heath, University of South Florida
* Date: 9/2/99
*******************************************************************************/
void linear_optical_density(CACHEIM *image, CACHEIM *source_image, int max_depth,
   int max_cache_size, char *mode, double *od_offset, double *od_scale,
   double *min_density, double *max_density, double A, double B)
{
   CACHEIM *newim=NULL;
   int i, r, c;
   int datatype = USHORTNUM;
   OD_REMAPPER_DATA *od_remapper_data=NULL;
   double od;

   /****************************************************************************
   * If the mode is "r+b" show an error because we can not already have the
   * result of the transformation in a file.
   ****************************************************************************/
   if(strcmp(mode, "r+b") == 0){
      fprintf(stderr, "Error. Tried to open a linear_optical_density image in mode \"r+b\".\n");
      exit(1);
   }

   *od_offset = OD_OFFSET;
   *od_scale = OD_SCALE;
   *min_density = MIN_DENSITY;
   *max_density = MAX_DENSITY;

   /****************************************************************************
   * Allocate a virtual image. If thw whole image is to be cached, then set
   * up a temporary image as 1/2 of the max_cache_size. This is necessary
   * because otherwise the computations will never be performed because
   * differt getpixel and putpixel functions are used when the image is
   * completely contained in the image. Since those functions end up calling
   * the procedure that calculates the aggregates they would never be
   * called if we didn't do this.
   ****************************************************************************/
   if(((source_image->rows * source_image->cols)*dt_size[datatype]) <= max_cache_size){
      if(allocate_cached_image(image, source_image->rows, source_image->cols,
         datatype, "tobecomputed", "rb", max_depth,
	 (source_image->rows * source_image->cols * dt_size[datatype])/2, 0, FALSE,
         source_image->resolution) == 0) exit(1);
   }
   else{
      if(allocate_cached_image(image, source_image->rows, source_image->cols,
         datatype, "tobecomputed", "rb", max_depth, max_cache_size, 0, FALSE,
         source_image->resolution) == 0) exit(1);
   }

   /****************************************************************************
   * Set up the data necessary for computing the convolution.
   ****************************************************************************/
   if((image->calculationdata = (void *) mymalloc("image->calculationdata",
      sizeof(OD_REMAPPER_DATA) + (65536*sizeof(USHORT)))) == NULL){
      fprintf(stderr, "Malloc error of calculation data in linear_optical_density!\n");
      exit(1);
   }

   od_remapper_data = (OD_REMAPPER_DATA *)image->calculationdata;
   od_remapper_data->od_offset = *od_offset;
   od_remapper_data->od_scale = *od_scale;

   od_remapper_data->map = (USHORT *)(((char *)(image->calculationdata)) + sizeof(OD_REMAPPER_DATA));

   for(i=0;i<65536;i++){
      if(i == 0) i = 1;
      od = A + (B * (double)i);

      if(od > MAX_DENSITY) od = MAX_DENSITY;
      od -= MIN_DENSITY;

      if(od < 0.0) od = 0.0;

      od = roundval( OD_OFFSET - OD_SCALE * od );
      (od_remapper_data->map)[i] = (USHORT)od;
   }

   /****************************************************************************
   * Set up the data necessary for computing the transformation.
   ****************************************************************************/
   image->sourceimlist = (void **) mycalloc("image->sourceimlist", 1, sizeof(CACHEIM *));
   ((CACHEIM **)(image->sourceimlist))[0] = source_image;

   image->get_cache_segment = get_optical_density_segment;

   /****************************************************************************
   * If the mode is rb then we are finished. We have told the image how to
   * compute its values.
   ****************************************************************************/
   if((strcmp(mode, "rb") == 0) && ((source_image->rows * source_image->cols*dt_size[datatype]) > max_cache_size)) return;

   printf("Precalculating the entire linear_optical_density image.\n");

   /****************************************************************************
   * If the mode is "w+b" then we want to immediately calculate the result for
   * the entire image and cache the result in a file.
   ****************************************************************************/
   if((newim = (CACHEIM *) mycalloc("newim", 1, sizeof(CACHEIM))) == NULL){
      fprintf(stderr, "Malloc error of newim in linear_optical_density!\n");
      exit(1);
   }

   /****************************************************************************
   * Allocate a virtual image.
   ****************************************************************************/
   if(allocate_cached_image(newim, source_image->rows, source_image->cols,
      datatype, NULL, "w+b", max_depth, max_cache_size, 0, FALSE,
      source_image->resolution) == 0) exit(1);

   /****************************************************************************
   * Compute the result of the aggregation for every pixel in the image.
   ****************************************************************************/
   for(r=0;r<(source_image->rows);r++){
      for(c=0;c<(source_image->cols);c++){
         (*(newim->putpixel))(newim, r, c, (*(image->getpixel))(image, r, c));
      }
   }

   deallocate_cached_image(image);
   *image = *newim;
   myfree("newim", newim);
}

/*******************************************************************************
* Procedure: log10_optical_density
* Purpose: This function maps the image to optical density. It uses the equation
* optical_density = A + (B * log10(gray_level)). The real valued optical
* densities are then converted to integer values by scaling optical densities
* in the range MIN_DENSITY to MAX_DENSITY to fall in the range 65535->0.
* NAME: Michael Heath, University of South Florida
* Date: 9/2/99
*******************************************************************************/
void log10_optical_density(CACHEIM *image, CACHEIM *source_image, int max_depth,
   int max_cache_size, char *mode, double *od_offset, double *od_scale,
   double *min_density, double *max_density, double A, double B)
{
   CACHEIM *newim=NULL;
   int i, r, c;
   int datatype = USHORTNUM;
   OD_REMAPPER_DATA *od_remapper_data=NULL;
   double od;

   /****************************************************************************
   * If the mode is "r+b" show an error because we can not already have the
   * result of the transformation in a file.
   ****************************************************************************/
   if(strcmp(mode, "r+b") == 0){
      fprintf(stderr, "Error. Tried to open a log10_optical_density image in mode \"r+b\".\n");
      exit(1);
   }

   *od_offset = OD_OFFSET;
   *od_scale = OD_SCALE;
   *min_density = MIN_DENSITY;
   *max_density = MAX_DENSITY;

   /****************************************************************************
   * Allocate a virtual image. If thw whole image is to be cached, then set
   * up a temporary image as 1/2 of the max_cache_size. This is necessary
   * because otherwise the computations will never be performed because
   * differt getpixel and putpixel functions are used when the image is
   * completely contained in the image. Since those functions end up calling
   * the procedure that calculates the aggregates they would never be
   * called if we didn't do this.
   ****************************************************************************/
   if(((source_image->rows * source_image->cols)*dt_size[datatype]) <= max_cache_size){
      if(allocate_cached_image(image, source_image->rows, source_image->cols,
         datatype, "tobecomputed", "rb", max_depth,
	 (source_image->rows * source_image->cols * dt_size[datatype])/2, 0, FALSE,
         source_image->resolution) == 0) exit(1);
   }
   else{
      if(allocate_cached_image(image, source_image->rows, source_image->cols,
         datatype, "tobecomputed", "rb", max_depth, max_cache_size, 0, FALSE,
         source_image->resolution) == 0) exit(1);
   }

   /****************************************************************************
   * Set up the data necessary for computing the convolution.
   ****************************************************************************/
   if((image->calculationdata = (void *) mymalloc("image->calculationdata",
      sizeof(OD_REMAPPER_DATA) + (65536*sizeof(USHORT)))) == NULL){
      fprintf(stderr, "Malloc error of calculation data in log10_optical_density!\n");
      exit(1);
   }

   od_remapper_data = (OD_REMAPPER_DATA *)image->calculationdata;
   od_remapper_data->od_offset = *od_offset;
   od_remapper_data->od_scale = *od_scale;

   od_remapper_data->map = (USHORT *)(((char *)(image->calculationdata)) + sizeof(OD_REMAPPER_DATA));

   for(i=0;i<65536;i++){
      if(i == 0) i = 1;
      od = A + (B * log10((double)i));

      if(od > MAX_DENSITY) od = MAX_DENSITY;
      od -= MIN_DENSITY;

      if(od < 0.0) od = 0.0;

      od = roundval( OD_OFFSET - OD_SCALE * od );
      (od_remapper_data->map)[i] = (USHORT)od;
   }

   /****************************************************************************
   * Set up the data necessary for computing the transformation.
   ****************************************************************************/
   image->sourceimlist = (void **) mycalloc("image->sourceimlist", 1, sizeof(CACHEIM *));
   ((CACHEIM **)(image->sourceimlist))[0] = source_image;

   image->get_cache_segment = get_optical_density_segment;

   /****************************************************************************
   * If the mode is rb then we are finished. We have told the image how to
   * compute its values.
   ****************************************************************************/
   if((strcmp(mode, "rb") == 0) && ((source_image->rows * source_image->cols*dt_size[datatype]) > max_cache_size)) return;

   printf("Precalculating the entire log10_optical_density image.\n");

   /****************************************************************************
   * If the mode is "w+b" then we want to immediately calculate the result for
   * the entire image and cache the result in a file.
   ****************************************************************************/
   if((newim = (CACHEIM *) mycalloc("newim", 1, sizeof(CACHEIM))) == NULL){
      fprintf(stderr, "Malloc error of newim in log10_optical_density!\n");
      exit(1);
   }

   /****************************************************************************
   * Allocate a virtual image.
   ****************************************************************************/
   if(allocate_cached_image(newim, source_image->rows, source_image->cols,
      datatype, NULL, "w+b", max_depth, max_cache_size, 0, FALSE,
      source_image->resolution) == 0) exit(1);

   /****************************************************************************
   * Compute the result of the aggregation for every pixel in the image.
   ****************************************************************************/
   for(r=0;r<(source_image->rows);r++){
      for(c=0;c<(source_image->cols);c++){
         (*(newim->putpixel))(newim, r, c, (*(image->getpixel))(image, r, c));
      }
   }

   deallocate_cached_image(image);
   *image = *newim;
   myfree("newim", newim);
}

/*******************************************************************************
* Function:
* Purpose: The functions to remap to optical density stores scaled values of
* the optical densities. This function will give the rescaled value for the
* optical density that is passed to the function.
* Name: Michael Heath, University of South Florida
* Date: 1/28/2000
*******************************************************************************/
unsigned short int od_to_scaled_od(double od)
{
   /* printf("od = %f - %f * (%f - %f)\n", (float)OD_OFFSET, (float)OD_SCALE, od, (float)MIN_DENSITY); */

   if(od < MIN_DENSITY){
      printf("Warning in od_to_scaled_od().\n");
      printf("You are asking for the value for density=%f (using clipped density of 0.05).\n", od, (float)MIN_DENSITY);
      return((USHORT)OD_OFFSET);
   }

   if(od > MAX_DENSITY){
      printf("Warning in od_to_scaled_od().\n");
      printf("You are asking for the value for density=%f (using clipped density of %f).\n", od, (float)MAX_DENSITY);
      return((USHORT)(roundval(OD_OFFSET - OD_SCALE * (MAX_DENSITY-MIN_DENSITY))));
   }

   return((USHORT)(roundval(OD_OFFSET - OD_SCALE * (od-MIN_DENSITY))));
}

/*******************************************************************************
* Procedure: map_with_ushort_lut
* Purpose: This function maps the image with the supplied look-up-table.
* NAME: Michael Heath, University of South Florida
* Date: 9/2/99
*******************************************************************************/
void map_with_ushort_lut(CACHEIM *image, CACHEIM *source_image, int max_depth,
   int max_cache_size, int datatype, char *mode, unsigned short int *lut)
{
   CACHEIM *newim=NULL;
   int i, r, c;
   OD_REMAPPER_DATA *od_remapper_data=NULL;

   /****************************************************************************
   * If the mode is "r+b" show an error because we can not already have the
   * result of the transformation in a file.
   ****************************************************************************/
   if(strcmp(mode, "r+b") == 0){
      fprintf(stderr, "Error. Tried to open a map_with_ushort_lut image in mode \"r+b\".\n");
      exit(1);
   }

   /****************************************************************************
   * The datatype specifies the datatype of the output image from applying
   * this remapping function. We will only allow output image datatypes of
   * unsigned char and unsigned short int.
   ****************************************************************************/
   if(!((datatype == UCHARNUM) || (datatype == USHORTNUM))){
      fprintf(stderr, "The datatype must be UCHAR or USHORT in map_with_ushort_lut().\n");
      exit(1);
   }

   /****************************************************************************
   * Allocate a virtual image. If thw whole image is to be cached, then set
   * up a temporary image as 1/2 of the max_cache_size. This is necessary
   * because otherwise the computations will never be performed because
   * differt getpixel and putpixel functions are used when the image is
   * completely contained in the image. Since those functions end up calling
   * the procedure that calculates the aggregates they would never be
   * called if we didn't do this.
   ****************************************************************************/
   if(((source_image->rows * source_image->cols)*dt_size[datatype]) <= max_cache_size){
      if(allocate_cached_image(image, source_image->rows, source_image->cols,
         datatype, "tobecomputed", "rb", max_depth,
	 (source_image->rows * source_image->cols * dt_size[datatype])/2, 0, FALSE,
         source_image->resolution) == 0) exit(1);
   }
   else{
      if(allocate_cached_image(image, source_image->rows, source_image->cols,
         datatype, "tobecomputed", "rb", max_depth, max_cache_size, 0, FALSE,
         source_image->resolution) == 0) exit(1);
   }

   /****************************************************************************
   * Set up the data necessary for computing the lut.
   ****************************************************************************/
   if((image->calculationdata = (void *) mymalloc("image->calculationdata",
      sizeof(OD_REMAPPER_DATA) + (65536*sizeof(USHORT)))) == NULL){
      fprintf(stderr, "Malloc error of calculation data in map_with_ushort_lut!\n");
      exit(1);
   }

   od_remapper_data = (OD_REMAPPER_DATA *)image->calculationdata;
   od_remapper_data->od_offset = 0;
   od_remapper_data->od_scale = 0;

   od_remapper_data->map = (USHORT *)(((char *)(image->calculationdata)) + sizeof(OD_REMAPPER_DATA));

   for(i=0;i<65536;i++) (od_remapper_data->map)[i] = (USHORT)lut[i];

   /****************************************************************************
   * Set up the data necessary for computing the transformation.
   ****************************************************************************/
   image->sourceimlist = (void **) mycalloc("image->sourceimlist", 1, sizeof(CACHEIM *));
   ((CACHEIM **)(image->sourceimlist))[0] = source_image;

   image->get_cache_segment = get_optical_density_segment;

   /****************************************************************************
   * If the mode is rb then we are finished. We have told the image how to
   * compute its values.
   ****************************************************************************/
   if((strcmp(mode, "rb") == 0) && ((source_image->rows * source_image->cols*dt_size[datatype]) > max_cache_size)) return;

   printf("Precalculating the entire map_with_ushort_lut image.\n");

   /****************************************************************************
   * If the mode is "w+b" then we want to immediately calculate the result for
   * the entire image and cache the result in a file.
   ****************************************************************************/
   if((newim = (CACHEIM *) mycalloc("newim", 1, sizeof(CACHEIM))) == NULL){
      fprintf(stderr, "Malloc error of newim in map_with_ushort_lut!\n");
      exit(1);
   }

   /****************************************************************************
   * Allocate a virtual image.
   ****************************************************************************/
   if(allocate_cached_image(newim, source_image->rows, source_image->cols,
      datatype, NULL, "w+b", max_depth, max_cache_size, 0, FALSE,
      source_image->resolution) == 0) exit(1);

   /****************************************************************************
   * Compute the result of the aggregation for every pixel in the image.
   ****************************************************************************/
   for(r=0;r<(source_image->rows);r++){
      for(c=0;c<(source_image->cols);c++){
         (*(newim->putpixel))(newim, r, c, (*(image->getpixel))(image, r, c));
      }
   }

   deallocate_cached_image(image);
   *image = *newim;
   myfree("newim", newim);
}
