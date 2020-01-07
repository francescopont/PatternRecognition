/*******************************************************************************
* File: afumfeature.c
* Purpose: This file contains code to locate candidate mass regions in a
* mammogram using Average Fraction Under the Minimum (AFUM) filtering. This
* filter was developed by Michael Heath in the Computer Vision Lab and the
* University of South Florida.
* Name: Michael Heath, University of South Florida
* Date: 7/26/99
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "virtual_image.h"
#include "optical_density.h"
#include "breast.h"
#include "myalloc.h"
#include "mikesfileio.h"
#include "ICSio.h"
#include "overlay.h"

#define VERSION "1.0.0"
#define VERSIONDATE "January 21, 2000"

#define LCC 0
#define LMLO 1
#define RCC 2
#define RMLO 3

#define MAP_INTENSITIES 0 /* Make this 0 if you do not want to remap the intenisites. */
#define JUST_MAKE_IMAGE 0  /* Make this 0 to actually generate the feature images. */
#define AGGREGATE_BY_MEDIAN  /* Comment this out to aggregate by averageing rathern than by median. */

#define AIM_RESOLUTION 300.0

int VERBOSE = 0;

static void print_help();
static void get_command_line_arguments(int argc, char *argv[]);
void afum_process(char *filename, int full_rows, int full_cols, int bpp, int hb,
       int swap, int projection, int apex, float resolution, int numpoints,
       float *fxcoord, float *fycoord,
       int *feature_rows, int *feature_cols, int mapmode, double A, double B,
       unsigned char **aucfeatureim);
void aggregate(CACHEIM *image, CACHEIM *source_image, double aggfactor,
   int max_depth, int max_cache_size, int datatype, char *mode);
void aggregate_median(CACHEIM *image, CACHEIM *source_image, int aggfactor,
    int max_depth, int max_cache_size, int datatype, char *mode);
void polyscan_coords(int numpoints, int *x_coord, int *y_coord, unsigned char *image, int rows, int cols);

static int rows=0, cols=0, bpp=0, hb=0, swap=0, projection=0, apex=0, mapmode=NOMAP;
static float resolution;
static double A=0.0, B=0.0;
static char *image_filename = NULL, *segmentation_filename;
char *expanded[4] = {"LEFT_CC", "LEFT_MLO", "RIGHT_CC", "RIGHT_MLO"};

int main(int argc, char *argv[])
{
   FILE *fpauc=NULL;
   char spicfeature_filename[200], auc_filename[200];
   float ax_0, ax_1, ay_0, ay_1;
   int numpoints=0, feature_rows, feature_cols;
   float *xcoord=NULL, *ycoord = NULL, start_col, start_row, end_col, end_row;
   float startx, starty, endx, endy;
   unsigned char *aucfeatureim=NULL;
   int segmentation_rows, segmentation_cols;

   /****************************************************************************
   * Get the arguments from the command line.
   ****************************************************************************/
   get_command_line_arguments(argc, argv);

   /****************************************************************************
   * Read in the data from the breast segmentation file.
   ****************************************************************************/
   read_segmentation_file(segmentation_filename, &numpoints, &xcoord, &ycoord,
      &startx, &starty, &endx, &endy, &segmentation_rows, &segmentation_cols);

   /****************************************************************************
   * Convert the segmentation coordinates to our current working resolution.
   ****************************************************************************/
   convert_segmentation_coordinates(segmentation_rows, segmentation_cols,
      numpoints, xcoord, ycoord, &startx, &starty, &endx, &endy, rows, cols);

   /****************************************************************************
   * Compute the spiculated feature images.
   ****************************************************************************/
   afum_process(image_filename, rows, cols, bpp, hb,
       swap, projection, apex, resolution, numpoints,
       xcoord, ycoord, &feature_rows, &feature_cols,
       mapmode, A, B, &aucfeatureim);

   /****************************************************************************
   * Write the afum feature image out to a file.
   ****************************************************************************/
   sprintf(auc_filename, "%s.afum.pgm", image_filename);
   if((fpauc = fopen(auc_filename, "wb")) == NULL) return(0);
   fprintf(fpauc, "P5\n%d %d\n255\n", feature_cols, feature_rows);
   fwrite(aucfeatureim, 1, feature_rows*feature_cols, fpauc);
   fclose(fpauc);
}

/*******************************************************************************
* Function: get_command_line_arguments
* Purpose: This function is self explanitory. It just extracts the command line
* arguments.
* Name: Michael Heath, University of South Florida
* Date: 7/26/99
*******************************************************************************/
static void get_command_line_arguments(int argc, char *argv[])
{
   int i;
   char *ics_filename = NULL;
   int view;

   void get_ddsm_image_info(char *ics_filename, int view, char **filename, int *rows, int *cols,
      int *bpp, float *resolution, int *projection, int *apex);

   rows = 0;
   cols = 0;
   bpp = 0;
   image_filename = NULL;
   apex = -1;
   projection = -1;
   resolution = 0.0;

   resolution = ASSUMED_RESOLUTION;
   swap = 0;

   for(i=1;i<argc;i++){
      if(strcmp(argv[i], "-i") == 0){ image_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-inp") == 0){ cols = atoi(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-inl") == 0){ rows = atoi(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-bpp") == 0){ bpp = atoi(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-microns") == 0){ resolution = atof(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-version") == 0){
         printf("\n\n%s Version: %s %s\n\n", argv[0], VERSION, VERSIONDATE);
         exit(1);
      }
      else if(strcmp(argv[i], "-swap") == 0){ swap = 1; }
      else if(strcmp(argv[i], "-v") == 0){ VERBOSE = 1; }
      else if(strcmp(argv[i], "-projection") == 0){
         if(strcmp(argv[i+1], "MLO") == 0) projection = MLO;
         if(strcmp(argv[i+1], "CC") == 0) projection = CC;
         i++;
      }
      else if(strcmp(argv[i], "-apex") == 0){
         if(strcmp(argv[i+1], "LEFT") == 0) apex = APEXLEFT;
         if(strcmp(argv[i+1], "RIGHT") == 0) apex = APEXRIGHT;
      }
      else if(strcmp(argv[i], "-ics") == 0){ ics_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-view") == 0){
         if(strcmp(argv[i+1], "LEFT_MLO") == 0) view = LMLO;
         if(strcmp(argv[i+1], "RIGHT_MLO") == 0) view = RMLO;
         if(strcmp(argv[i+1], "LEFT_CC") == 0) view = LCC;
         if(strcmp(argv[i+1], "RIGHT_CC") == 0) view = RCC;
         i++;
      }
      else if(strcmp(argv[i], "-segment") == 0){ segmentation_filename = argv[i+1]; i++; }
      else if(strncmp(argv[i], "-linear", 7) == 0){
         mapmode = LINEARMAP;
         A = atof(argv[i+1]);
         B = atof(argv[i+2]);
         if(VERBOSE) printf("\nWe are going to use optical_density = %f + (%f * GL)\n", A, B);
         i+=2;
      }
      else if(strncmp(argv[i], "-log10", 6) == 0){
         mapmode = LOG10MAP;
         A = atof(argv[i+1]);
         B = atof(argv[i+2]);
         if(VERBOSE) printf("\nWe are going to use optical_density = %f + (%f * log10(GL))\n", A, B);
         i+=2;
      }
   }

   if(ics_filename != NULL){
      get_ddsm_image_info(ics_filename, view, &image_filename, &rows, &cols, &bpp, &resolution, &projection, &apex);
   }

   if((image_filename != NULL) && (segmentation_filename == NULL)){
      segmentation_filename = (char *) calloc(strlen(image_filename)+5, sizeof(char));
      sprintf(segmentation_filename, "%s.sgt", image_filename);
   }

   if((rows == 0) || (cols == 0) || (bpp == 0) || (image_filename == NULL) || (apex == -1) ||
      (projection == -1) || (resolution == 0.0) || (segmentation_filename == NULL)) print_help();

}

/*******************************************************************************
* Function: get_ddsm_image_info
* Purpose: This function retrieves some data about the image from the ics
* file for a DDSM case.
* Name: Michael Heath, University of South Florida
* Date: 1/30/2000
*******************************************************************************/
void get_ddsm_image_info(char *ics_filename, int view, char **filename, int *rows, int *cols,
   int *bpp, float *resolution, int *projection, int *apex)
{
   FILE *fp=NULL;
   char line[200], caseid[40];
   int bits, i;
   char tmpfilename[100] = {'\0'};

   /* C_0234_1.LEFT_CC.LJPEG.1 */

   /****************************************************************************
   * Open the .ics file.
   ****************************************************************************/
   if((fp = fopen(ics_filename, "r")) == NULL){
      fprintf(stderr, "Error opening the ics file named %s.\n", ics_filename);
      exit(1);
   }

   if((view == LMLO) || (view == RMLO)) *projection = MLO;
   if((view == LCC) || (view == RCC)) *projection = CC;

   while((fgets(line, 200, fp) != NULL) && !feof(fp)){
      if(strncmp(line, "filename", strlen("filename")) == 0){
         sscanf(line, "filename %s", caseid);

         if(caseid[0] == 'A'){
            if((view == LMLO) || (view == LCC)) *apex = APEXLEFT;
            if((view == RMLO) || (view == RCC)) *apex = APEXRIGHT;
         }
         if(caseid[0] == 'B'){
            if((view == LMLO) || (view == LCC)) *apex = APEXRIGHT;
            if((view == RMLO) || (view == RCC)) *apex = APEXLEFT;
         }
         if(caseid[0] == 'C'){
            if((view == LMLO) || (view == LCC)) *apex = APEXRIGHT;
            if((view == RMLO) || (view == RCC)) *apex = APEXLEFT;
         }
         if(caseid[0] == 'D'){
            if((view == LMLO) || (view == LCC)) *apex = APEXRIGHT;
            if((view == RMLO) || (view == RCC)) *apex = APEXLEFT;
         }
         *filename = (char *) calloc(50, sizeof(char));
         strncpy(tmpfilename, ics_filename, strlen(ics_filename) - 4);
         /* printf("tmpfilename = %s\n", tmpfilename); */
         for(i=(strlen(tmpfilename)-8);i<strlen(tmpfilename);i++)
            if(tmpfilename[i] == '-') tmpfilename[i] = '_';
         /* printf("tmpfilename = %s\n", tmpfilename); */
         sprintf(*filename, "%s.%s.LJPEG.1", tmpfilename, expanded[view]);
         /* printf("filename = %s\n", *filename); */
      }
      if(strncmp(line, expanded[view], strlen(expanded[view])) == 0){
         sscanf(line, "%*s %*s %d %*s %d %*s %d %*s %f", rows, cols, &bits, resolution);
         *bpp = (int)ceil((double)bits / 8.0);
      }
   }
}

/*******************************************************************************
* Function: print_help
* Purpose: This function is self explanitory. It just prints out command line help.
* Name: Michael Heath, University of South Florida
* Date: 7/26/99
*******************************************************************************/
static void print_help()
{

   printf("\n********************************************************************************\n");
   printf("This program computes a feature image using an Average Fraction Under the\n");
   printf("Minumum (AFUM) filtering process. Local maxima in this filtered output are\n");
   printf("likely positions of masslike objects. The larger the value of the local\n");
   printf("maximum the more susicious the position in the image.\n");
   printf("********************************************************************************\n\n");
   printf("<USAGE> afumfeature [-ics filename.ics -view VIEW] |\n");
   printf("                    [-i filename -inp # -inl # -bpp BPP -microns # -hb #\n");
   printf("                    -projection PROJECTION -apex APEX]\n");
   printf("                    [-swap] [-v] [-version]\n");
   printf("                    [-linear # #] | [-log10 # #]\n");
   printf("                    -segment filename.sgt\n\n");
   printf("   -ics        The name of an ics file from the DDSM database. This supplies the\n");
   printf("               image filename, view, inp, inl, bpp, micons, hb, projection\n");
   printf("               and apex automatically.\n");
   printf("   -view       The mammogram to process for this case. VIEW can be any one of\n");
   printf("               the following: LEFT_MLO,RIGHT_MLO,LEFT_CC or RIGHT_CC.\n");
   printf("   -segment    The breast segmentation file that was previously computed.\n");
   printf("   -i          The filename of the image to process.\n");
   printf("   -inp        The number of columns in the image.\n");
   printf("   -inl        The number of rows in the image.\n");
   printf("   -bpp        The number of bytes per pixel. BPP can be either 1 or 2.\n");
   printf("   -microns    The sampling resolution of the image in microns.\n");
   printf("   -projection The mammographic projection of this image. PROJECTION can be\n");
   printf("               either MLO or CC.\n");
   printf("   -apex       The side of the image that the nipple is facing. APEX can be\n");
   printf("               either LEFT or RIGHT.\n");
   printf("   -swap       If this flag is used the bytes of a 16-bit image will be swapped\n");
   printf("               when image data is read from a file.\n");
   printf("   -hb         The number of header bytes in the image file before image data\n");
   printf("               is found.\n");
   printf("   -linear     Rescale to optical density using a linear function of the form\n");
   printf("               optical_density = A + (B * GL). A and B are specified as double\n");
   printf("               precision values with a space in between them.\n");
   printf("   -log10      Rescale to optical density using a log function.\n");
   printf("               optical_density = A + (B * log10(GL)). A and B are specified as\n");
   printf("               double precision values with a space in between them.\n");
   printf("   -v          Run the program in verbose mode.\n");
   printf("   -version    Print the version and exit.\n\n");

   exit(1);
}

/*******************************************************************************
* Function: afum_process
* Purpose: This function generates the AFUM filtered image for a mammogram. The
* image is refuced in resolution and is filtered.
* Name: Michael Heath, University of South Florida
* Date: 7/26/99
*******************************************************************************/
void afum_process(char *filename, int full_rows, int full_cols, int bpp, int hb,
       int swap, int projection, int apex, float resolution, int numpoints,
       float *fxcoord, float *fycoord, int *feature_rows, int *feature_cols,
       int mapmode, double A, double B, unsigned char **aucfeatureim)
{
   CACHEIM rawimage, aggregatedimage, *imageptr=NULL, tempimage, tempimage0;
   int aggfactor, r, c, rows, cols, n;
   int basepos, pos;
   int *xcoord=NULL, *ycoord=NULL;
   float scale, processing_resolution=AIM_RESOLUTION;
   double od_offset, od_scale, min_density, max_density;
   UCHAR *pconverted_image=NULL;
   unsigned short int *lut=NULL;
   double xx;

   void newaucfeature_ushort(CACHEIM *image, unsigned char *featureim, int radius, int deltar,
       unsigned char *breastmask);

   memset(&rawimage, 0, sizeof(CACHEIM));
   memset(&aggregatedimage, 0, sizeof(CACHEIM));
   memset(&tempimage0, 0, sizeof(CACHEIM));
   memset(&tempimage, 0, sizeof(CACHEIM));

   /****************************************************************************
   * Read in the image as a cached image.
   ****************************************************************************/
   if(bpp == 1){
      if(allocate_cached_image(&rawimage, full_rows, full_cols, UCHARNUM, filename,
         "rb", 12, (full_rows*full_cols*sizeof(UCHAR)/16), hb, (USHORT)swap,
         resolution) == 0) exit(1);
   }
   else if(bpp == 2){
      if(allocate_cached_image(&rawimage, full_rows, full_cols, USHORTNUM, filename,
         "rb", 12, (full_rows*full_cols*sizeof(USHORT)/16), hb, (USHORT)swap,
         resolution) == 0) exit(1);
   }
   else{
      fprintf(stderr, "Invalid bits per pixel (%d).\n", bpp);
      exit(1);
   }

   /****************************************************************************
   * If the user wants to work in scaled optical density space, map to optical
   * density.
   ****************************************************************************/
   if(mapmode == LINEARMAP){
      if(VERBOSE) printf("Mapping to optical density using linear mapping.\n");
      linear_optical_density(&tempimage0, &rawimage, 12,
         (rawimage.rows * rawimage.cols * dt_size[USHORTNUM]) / 2, "rb",
         &od_offset, &od_scale, &min_density, &max_density, A, B);
   }
   else if(mapmode == LOG10MAP){
      if(VERBOSE) printf("Mapping to optical density using log10 mapping.\n");
      log10_optical_density(&tempimage0, &rawimage, 12,
         (rawimage.rows * rawimage.cols * dt_size[USHORTNUM]) / 16, "rb",
         &od_offset, &od_scale, &min_density, &max_density, A, B);
   }
   else{
      if(VERBOSE) printf("Processing in gray level mode.\n");
      tempimage0 = rawimage;
   }

   /****************************************************************************
   * Added this to be able to arbitrarily map the image intensities.
   ****************************************************************************/
   if(MAP_INTENSITIES){
      lut = (unsigned short int *) calloc(65536, sizeof(unsigned short int));
      for(r=0;r<65536;r++){
         xx = (double)r / 65535.0;
         xx = 65535.0 * (pow(10.0, xx) / 10.0);
         lut[r] = (unsigned short int)floor(xx);
      }
      map_with_ushort_lut(&tempimage, &tempimage0, 12,
      (rawimage.rows * rawimage.cols * dt_size[USHORTNUM]) / 16, USHORTNUM, "rb", lut);
      free(lut);
   }
   else tempimage = tempimage0;

   /****************************************************************************
   * Aggregate the image to effectivly change to a coarser resolution.
   ****************************************************************************/
   if((rawimage.resolution == 0) || (rawimage.resolution >= AIM_RESOLUTION)){
      aggfactor = 1;
      imageptr = &tempimage;
   }
   else{
      aggfactor = (int)ceil(AIM_RESOLUTION / rawimage.resolution);
#ifndef AGGREGATE_BY_MEDIAN
      if(VERBOSE) printf("Aggregating by averaging to resample the image.\n");
      aggregate(&aggregatedimage, &tempimage, (double)aggfactor, 12,
         (tempimage.rows * tempimage.cols * dt_size[USHORTNUM]) /
         (aggfactor*aggfactor*2), USHORTNUM, "rb");
#else
      if(VERBOSE) printf("Aggregating by finding the median to resample the image.\n");
      aggregate_median(&aggregatedimage, &tempimage, (int)aggfactor, 12,
         (tempimage.rows * tempimage.cols * dt_size[USHORTNUM]) /
         (aggfactor*aggfactor*2), USHORTNUM, "rb");
#endif
      imageptr = &aggregatedimage;
   }

   if(VERBOSE) printf("The aggfactor = %d.\n", aggfactor);

   /****************************************************************************
   * Set the basepos to indicate which side of the image has the APEX of the
   * breast.
   ****************************************************************************/
   if(apex == APEXLEFT) basepos = 1;
   else basepos = 0;

   *feature_rows = rows = imageptr->rows;
   *feature_cols = cols = imageptr->cols;

   /****************************************************************************
   * We will want to set all of the pixels outside of the breast to the
   * value zero. Since we already have a polygon of the breast area,
   * we will scan convert the polygon to an image. Recall that even though
   * the polygon may have been found in a subsampled version of the image,
   * the polygon verticies are saved in the full image resolution scale.
   ****************************************************************************/
   xcoord = (int *) mycalloc("xcoord", numpoints, sizeof(int));
   ycoord = (int *) mycalloc("ycoord", numpoints, sizeof(int));

   /* scale = (resolution / processing_resolution); */

   scale = (1.0 / (float)aggfactor);

   for(n=0;n<numpoints;n++){
      xcoord[n] = (int)floor(fxcoord[n] * scale + 0.5);
      ycoord[n] = (int)floor(fycoord[n] * scale + 0.5);
   }

   if((pconverted_image = (UCHAR *) mycalloc("pconverted_image", rows * cols, sizeof(UCHAR))) == NULL){
      fprintf(stderr, "Calloc error of image!\n");
      exit(1);
   }
   polyscan_coords(numpoints, xcoord, ycoord, pconverted_image, rows, cols);

   /****************************************************************************
   * Exclude some of the pixels from the side of the mammogram that the
   * patient is on.
   ****************************************************************************/
   if(apex == APEXLEFT){
      for(r=0;r<rows;r++){
         for(c=0;c<(int)floor(5000.0 / (double)(aggfactor * resolution));c++)
            pconverted_image[r*cols+c] = 0;
      }
   }
   else{
      for(r=0;r<rows;r++){
         for(c=(cols-1-(int)floor(5000.0 / (double)(aggfactor * resolution)));c<cols;c++)
            pconverted_image[r*cols+c] = 0;
      }
   }

   /****************************************************************************
   * Write out the scan converted image of the breast region. This will be
   * needed by the detect program.
   ****************************************************************************/
   if(1){
      FILE *fptmp;
      char tmpfilename[200];

      sprintf(tmpfilename, "%s.polyscan.pgm", filename);
      fptmp = fopen(tmpfilename, "wb");
      fprintf(fptmp, "P5\n%d %d\n255\n", cols, rows);
      fwrite(pconverted_image, 1, rows*cols, fptmp);
      fclose(fptmp);
   }

   /****************************************************************************
   * This allows the program to simply output the reduced size image. This is
   * not really needed by any other program.
   ****************************************************************************/
   if(1){
      FILE *fpim=NULL;
      char just_make_image_filename[200];
      unsigned short int *thisrow=NULL;

      sprintf(just_make_image_filename, "%s.reduced.pgm", filename);
      if((fpim = fopen(just_make_image_filename, "wb")) == NULL){
         fprintf(stderr, "Error opening the image %s for writing.\n", just_make_image_filename);
         exit(1);
      }
      fprintf(fpim, "P5\n%d %d\n65535\n", imageptr->cols, imageptr->rows);
      thisrow = (unsigned short int *) calloc(imageptr->cols, sizeof(unsigned short int));
      for(r=0;r<imageptr->rows;r++){
         for(c=0;c<imageptr->cols;c++)
            thisrow[c] = (unsigned short int)  (*(imageptr->getpixel))(imageptr, r, c);
         fwrite(thisrow, sizeof(unsigned short int), imageptr->cols, fpim);
      }
      fclose(fpim);
      free(thisrow);
      if(JUST_MAKE_IMAGE) exit(1);
   }

   /****************************************************************************
   * Compute my feature image. This is the Average Fraction Under the Minimum
   * filter that I developed.
   ****************************************************************************/
   if(((*aucfeatureim) = (unsigned char *) calloc(rows*cols, sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Calloc error in afum_process().\n");
      exit(1);
   }

   newaucfeature_ushort(imageptr, (*aucfeatureim), 30, 10, pconverted_image);

   /* if(imageptr != &rawimage) deallocate_cached_image(&aggregatedimage); */

   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         if(pconverted_image[pos] != 255){
            (*aucfeatureim)[pos] = 0;
         }
      }
   }

   myfree("xcoord", xcoord);
   myfree("ycoord", ycoord);
   imageptr = NULL;
   myfree("pconverted_image", pconverted_image);
}

/*******************************************************************************
* Function: newaucfeature_ushort
* Purpose: To compute a feature that indicates how likely each pixel is part
* of a local intensity hill.
* Name: Michael Heath, University of South Florida
* Date: 10/27/99
*******************************************************************************/
void newaucfeature_ushort(CACHEIM *image, unsigned char *featureim, int radius, int deltar,
    unsigned char *breastmask)
{
   int rows, cols, rpos, cpos, r, c, pos;
   unsigned short int *image16bit=NULL;

   void afum_filter_ushort(unsigned short int *image, int rows, int cols, unsigned char *afumimage,
       int radius, int deltar, unsigned char *mask);

   rows = image->rows;
   cols = image->cols;

   /****************************************************************************
   * Create a 16-bit image from the 16-bit CACHEIM image.
   ****************************************************************************/
   if((image16bit = (unsigned short int *) calloc(rows*cols, sizeof(unsigned short int))) == NULL){
      fprintf(stderr, "Calloc error in newaucfeature_ushort().\n");
      exit(1);
   }

   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         image16bit[pos] = (unsigned short int) (*(image->getpixel))(image, r, c);
      }
   }

   afum_filter_ushort(image16bit, rows, cols, featureim, radius, deltar, breastmask);

   free(image16bit);
}

/*******************************************************************************
* Function: afum_filter_ushort
* Purpose: To apply the AFUM filter to an image. This computes the average
* fraction of the pixels in circles with a range of radii that are lower valued
* than the minimum valued pixel in a disk with a smaller radius.
* Name: Michael Heath, University of South Florida
* Date: 11/29/99
*******************************************************************************/
void afum_filter_ushort(unsigned short int *image, int rows, int cols, unsigned char *afumimage,
    int radius, int deltar, unsigned char *mask)
{
   typedef struct{
      int radius, row, col;
   }AFUMFILTER;
   int r, c, d, pos, spos, rr, cc, r1, r2, n, lastr;
   static int last_radius=0, numpoints=0, *num_below_min=NULL, *num_at_radius=NULL;
   static AFUMFILTER tempafumfilter, *afumfilter=NULL;
   static unsigned short int *min_in_disk=NULL;
   double afum;

   /****************************************************************************
   * Set up the data structure we can use to make the filtering process efficient.
   ****************************************************************************/
   if(last_radius != radius){

      /*************************************************************************
      * Allocate memory for the structure.
      *************************************************************************/
      if(afumfilter != NULL) free(afumfilter);
      if((afumfilter = (AFUMFILTER *) calloc((radius*2+1)*(radius*2+1), sizeof(AFUMFILTER))) == NULL){
         fprintf(stderr, "Calloc error in afumfilter().\n");
         exit(1);
      }
      if(num_at_radius != NULL) free(num_at_radius);
      if((num_at_radius = (int *) calloc(radius+1, sizeof(int))) == NULL){
         fprintf(stderr, "Calloc error in afumfilter().\n");
         exit(1);
      }
      if(num_below_min != NULL) free(num_below_min);
      if((num_below_min = (int *) calloc(radius+1, sizeof(int))) == NULL){
         fprintf(stderr, "Calloc error in afumfilter().\n");
         exit(1);
      }
      if(min_in_disk != NULL) free(min_in_disk);
      if((min_in_disk = (unsigned short int *) calloc(radius+1, sizeof(unsigned short int))) == NULL){
         fprintf(stderr, "Calloc error in afumfilter().\n");
         exit(1);
      }

      /*************************************************************************
      * Compute the coordinates and distance for points inside the filter.
      *************************************************************************/
      numpoints = 0;
      last_radius = radius;
      for(r=(-radius);r<=radius;r++){
         for(c=(-radius);c<=radius;c++){
            d = (int)floor(sqrt((double)(r*r+c*c)));
            if(d <= radius){
               afumfilter[numpoints].radius = d;
               afumfilter[numpoints].row = r;
               afumfilter[numpoints].col = c;
               numpoints++;
               num_at_radius[d]++;
            }
         }
      }

      /*************************************************************************
      * Sort the elements by increasing radius. This sort is simple. It is not
      * very efficient but its computational burdon is very small compared to
      * applying the filter to the image.
      *************************************************************************/
      for(r=0;r<(numpoints-1);r++){
         for(c=(r+1);c<numpoints;c++){
            if(afumfilter[r].radius > afumfilter[c].radius){
               tempafumfilter = afumfilter[r];
               afumfilter[r] = afumfilter[c];
               afumfilter[c] = tempafumfilter;
            }
         }
         /* printf("Filter[%d] = %d,%d,%d\n", r, afumfilter[r].radius, afumfilter[r].row, afumfilter[r].col); */
      }
   }

   /****************************************************************************
   * Apply the filter to the image.
   ****************************************************************************/
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){

         if(mask[pos] != 0){
            for(n=0;n<=radius;n++){
               num_below_min[n] = 0;
               min_in_disk[n] = 65535;
            }
            min_in_disk[0] = image[pos];
            lastr = -1;
            for(n=1;n<numpoints;n++){
               if(lastr != afumfilter[n].radius){
                  if(afumfilter[n].radius != 0) min_in_disk[afumfilter[n].radius] = min_in_disk[afumfilter[n].radius-1];
                  lastr = afumfilter[n].radius;
               }
               rr = r + afumfilter[n].row;
               cc = c + afumfilter[n].col;
               if((rr>=0)&&(rr<rows)&&(cc>=0)&&(cc<cols)){
                  spos = rr*cols+cc;
                  if(mask[spos] != 0){
                     r2 = lastr;
                     if(image[spos] < min_in_disk[r2]) min_in_disk[r2] = image[spos];
                     r1 = r2 - deltar;
                     if((r1 >= 0) && (image[spos] < min_in_disk[r1])){
                        num_below_min[r2]++;
                     }
                  }
               }
            }
            afum = 0;
            for(n=deltar;n<=radius;n++){
               if(num_below_min[n] > num_at_radius[n]) fprintf(stderr, "Error! %d,%d\n", num_below_min[n], num_at_radius[n]);
               afum += (double)num_below_min[n] / (double)num_at_radius[n];
            }
            if((afum/(double)((radius-deltar) + 1)) < 1.0)
               afumimage[pos] = (unsigned char)floor(255.0 * afum / (double)((radius-deltar) + 1));
         }
      }
   }
}
