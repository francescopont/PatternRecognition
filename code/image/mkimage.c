/*******************************************************************************
* File: mkimage.c
* Purpose: This file contains code to create an image file from a decompressed
* LJPEG.1 image. The image can be mapped to opical density and can be
* resampled.
* Name: Michael Heath, University of South Florida
* Date: 1/25/2000
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
#define VERSIONDATE "January 25, 2000"

#define LCC 0
#define LMLO 1
#define RCC 2
#define RMLO 3

#define PGM 1
#define TIFF 2

int VERBOSE = 0;

static void print_help();
static void get_command_line_arguments(int argc, char *argv[]);
void image_convert_and_save(char *filename, int full_rows, int full_cols, int bpp, int hb,
       int swap, float resolution, int mapmode, double A, double B, float outmicrons,
       int scalebits, float scalelow, float scalehigh);
void aggregate(CACHEIM *image, CACHEIM *source_image, double aggfactor,
   int max_depth, int max_cache_size, int datatype, char *mode);
void aggregate_median(CACHEIM *image, CACHEIM *source_image, int aggfactor,
    int max_depth, int max_cache_size, int datatype, char *mode);

static int rows=0, cols=0, bpp=0, hb=0, swap=0, mapmode=NOMAP;
static float resolution;
static double A=0.0, B=0.0;
static char *image_filename = NULL;
static int filemode = PGM;
char *expanded[4] = {"LEFT_CC", "LEFT_MLO", "RIGHT_CC", "RIGHT_MLO"};
float outmicrons, scalelow, scalehigh;
int scalebits;

int main(int argc, char *argv[])
{
   /****************************************************************************
   * Get the arguments from the command line.
   ****************************************************************************/
   get_command_line_arguments(argc, argv);

   /****************************************************************************
   * Compute the spiculated feature images.
   ****************************************************************************/
   image_convert_and_save(image_filename, rows, cols, bpp, hb,
       swap, resolution, mapmode, A, B, outmicrons, scalebits, scalelow, scalehigh);
}

/*******************************************************************************
* Function: get_command_line_arguments
* Purpose: This function is self explanitory. It just extracts the command line
* arguments.
* Name: Michael Heath, University of South Florida
* Date: 1/25/2000
*******************************************************************************/
static void get_command_line_arguments(int argc, char *argv[])
{
   int i;
   char *ics_filename = NULL;
   int view, apex, projection;

   void get_ddsm_image_info(char *ics_filename, int view, char **filename, int *rows, int *cols,
      int *bpp, float *resolution, int *projection, int *apex);

   rows = 0;
   cols = 0;
   bpp = 0;
   image_filename = NULL;
   resolution = 0.0;
   outmicrons = 0;
   scalebits = 16;
   scalelow = scalehigh = -1.0;

   swap = 0;

   for(i=1;i<argc;i++){
      if(strcmp(argv[i], "-i") == 0){ image_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-inp") == 0){ cols = atoi(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-inl") == 0){ rows = atoi(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-bpp") == 0){ bpp = atoi(argv[i+1]); i++; }
      else if(strncmp(argv[i], "-tiff", 4) == 0) filemode = TIFF;
      else if(strcmp(argv[i], "-8bit") == 0) scalebits = 8;
      else if(strcmp(argv[i], "-16bit") == 0) scalebits = 16;
      else if(strcmp(argv[i], "-microns") == 0){ resolution = atof(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-outmicrons") == 0){ outmicrons = atof(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-version") == 0){
         printf("\n\n%s Version: %s %s\n\n", argv[0], VERSION, VERSIONDATE);
         exit(1);
      }
      else if(strcmp(argv[i], "-swap") == 0){ swap = 1; }
      else if(strcmp(argv[i], "-v") == 0){ VERBOSE = 1; }
      else if(strcmp(argv[i], "-ics") == 0){ ics_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-view") == 0){
         if(strcmp(argv[i+1], "LEFT_MLO") == 0) view = LMLO;
         if(strcmp(argv[i+1], "RIGHT_MLO") == 0) view = RMLO;
         if(strcmp(argv[i+1], "LEFT_CC") == 0) view = LCC;
         if(strcmp(argv[i+1], "RIGHT_CC") == 0) view = RCC;
         i++;
      }
      else if(strncmp(argv[i], "-scale", 6) == 0){
         scalelow = atof(argv[i+1]);
         scalehigh = atof(argv[i+2]);
         i+=2;
      }
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

   if((rows == 0) || (cols == 0) || (bpp == 0) || (image_filename == NULL) ||
      (resolution == 0.0)) print_help();

   if(outmicrons == 0.0) outmicrons = resolution;

   if((scalelow < 0.0) || (scalehigh < 0.0)){
      if(mapmode == NOMAP){
         if(bpp == 1){
            scalelow = 0.0;
            scalehigh = 255.0;
         }
         if(bpp == 2){
            scalelow = 0.0;
            scalehigh = 65535.0;
         }
      }
      else{
         scalelow = MAX_DENSITY;
         scalehigh = MIN_DENSITY;
      }
   }
}

/*******************************************************************************
* Function: get_ddsm_image_info
* Purpose: This function was written to simplify the commind line interface for
* images in the DDSM database. Since the ".ics" file for a case contains a lot
* of information about each image in the case, we can use it to supply a lot of
* the requred command line attributes. In particular the following attributes
* are extracted from the ".ics" file: rows, cols, resolution (in microns),
* projection (MLO or CC), the bits per pixel, and the orientation of the image
* [is the apex of the breast (i.e where the nipple * is) on the left or right
* side of the image]. The orientation may be left or right for either the right
* or left breast because some people scan mammograms upsidedown from others).
* Date: documented on 11/8/99
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
   printf("This program makes a copy of a decompressed LJPEG.1 image. The image can be\n");
   printf("reduced in size and can be remapped to be linearly related to optical density.\n");
   printf("********************************************************************************\n\n");
   printf("<USAGE> mkimage [-ics filename.ics -view VIEW] |\n");
   printf("                [-i filename -inp # -inl # -bpp BPP -microns # -hb #]\n");
   printf("                [-swap] [-v] [-version] [-linear # #] | [-log10 # #]\n");
   printf("                [-outmicrons #] [-tif] [-scale # #] [-8bit | -16bit]\n\n");
   printf("   -ics          The name of an ics file from the DDSM database. This supplies\n");
   printf("                 the image filename, view, inp, inl, bpp, micons, hb,\n");
   printf("                 automatically.\n");
   printf("   -view         The mammogram to process for this case. VIEW can be any one of\n");
   printf("                 the following: LEFT_MLO,RIGHT_MLO,LEFT_CC or RIGHT_CC.\n");
   printf("   -i            The filename of the image to process.\n");
   printf("   -inp          The number of columns in the image.\n");
   printf("   -inl          The number of rows in the image.\n");
   printf("   -bpp          The number of bytes per pixel. BPP can be either 1 or 2.\n");
   printf("   -microns      The sampling resolution of the image in microns.\n");
   printf("   -swap         If this flag is used the bytes of a 16-bit image will be\n");
   printf("                 swapped when image data is read from a file.\n");
   printf("   -hb           The number of header bytes in the image file before image data\n");
   printf("                 is found.\n");
   printf("   -linear       Rescale to optical density using a linear function of the form\n");
   printf("                 optical_density = A + (B * GL). A and B are specified as double\n");
   printf("                 precision values with a space in between them.\n");
   printf("   -log10        Rescale to optical density using a log function.\n");
   printf("                 optical_density = A + (B * log10(GL)). A and B are specified as\n");
   printf("                 double precision values with a space in between them.\n");
   printf("   -scale        Rescale and clip the image from the range (# #) to 0 and\n");
   printf("                 2^(8bits)-1 or 2^(16bits)-1, respectively. If you specify\n");
   printf("                 -linear or -log10, these numbers should be optical densities.\n");
   printf("   -8bit         Make the output image have 8-bits per pixel.\n");
   printf("   -16bit        Make the output image have 16-bits per pixel. This is the\n");
   printf("                 default.\n");
   printf("   -outmicrons   The desired resolution of the output image in microns.\n");
   printf("   -tif          Output the image in TIFF format rather than in 8bit PGM or\n");
   printf("                 pseudo 16bit PGM format.\n");
   printf("   -v            Run the program in verbose mode.\n");
   printf("   -version      Print the version and exit.\n\n");

   exit(1);
}

/*******************************************************************************
* Function: image_conver_and_save
* Purpose: This function can transform an image to optical density and can
* resize the image. The image is saved to disk.
* Name: Michael Heath, University of South Florida
* Date: 1/25/2000
*******************************************************************************/
void image_convert_and_save(char *filename, int full_rows, int full_cols, int bpp, int hb,
       int swap, float resolution, int mapmode, double A, double B, float outmicrons,
       int scalebits, float scallow, float scalehigh)
{
   FILE *fpim=NULL;
   CACHEIM rawimage, aggregatedimage, image, tempimage, tempimage0;
   int r, c, rows, cols, n;
   float aggfactor;
   float scale;
   double od_offset, od_scale, min_density, max_density;
   unsigned short int *lut=NULL;
   char image_filename[200];
   unsigned char *thisrow_uchar=NULL;
   unsigned short int *thisrow_ushort=NULL;
   int write_simple_tiff_CACHEIM(char *filename, CACHEIM *image);
   float lowval, highval;

   memset(&rawimage, 0, sizeof(CACHEIM));
   memset(&aggregatedimage, 0, sizeof(CACHEIM));
   memset(&tempimage0, 0, sizeof(CACHEIM));
   memset(&tempimage, 0, sizeof(CACHEIM));
   memset(&image, 0, sizeof(CACHEIM));

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
         (rawimage.rows * rawimage.cols * dt_size[USHORTNUM]) / 16, "rb",
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
   * Aggregate the image to effectivly change to a coarser resolution.
   ****************************************************************************/
   if((rawimage.resolution == 0) || (rawimage.resolution >= outmicrons)){
      aggfactor = 1.0;
      tempimage = tempimage0;
   }
   else{
      aggfactor = outmicrons / rawimage.resolution;
      if(VERBOSE) printf("Aggregating by averaging to resample the image.\n");
      aggregate(&aggregatedimage, &tempimage0, (double)aggfactor, 12,
         (tempimage0.rows * tempimage0.cols * dt_size[USHORTNUM]) /
         (aggfactor*aggfactor*2), USHORTNUM, "rb");
      tempimage = aggregatedimage;
   }

   if(VERBOSE) printf("The aggfactor = %f.\n", aggfactor);

   /****************************************************************************
   * Added this to be able to arbitrarily map the image intensities.
   ****************************************************************************/
   if(mapmode == NOMAP){
      if(!((scalebits == 16) && (scalelow == 0) && (scalehigh == 65535))){
         lut = (unsigned short int *) calloc(65536, sizeof(unsigned short int));
         if(scalelow > scalehigh){
            fprintf(stderr, "Error! The smaller valued graylevel must be specified first!\n");
            fprintf(stderr, "You used -scale %f %f\n", scalelow, scalehigh);
            fprintf(stderr, "This would invert the image and you are not allowed to do that. Sorry!\n");
            exit(1);
         }
         if(VERBOSE){
            if(scalebits == 8) printf("GrayLevel scaling: (%dbit) (%f, %f) -> (0,255)\n", scalebits, scalelow, scalehigh);
            else if(scalebits == 16) printf("GrayLevel scaling: (%dbit) (%f, %f) -> (0, 65535)\n", scalebits, scalelow, scalehigh);
         }
         for(r=0;r<65536;r++){
            if(r >= scalehigh){
               if(scalebits == 8) lut[r] = 255;
               else lut[r] = 65535;
            }
            else if(r < scalelow) lut[r] = 0;
            else{
               if(scalebits == 8) lut[r] = 255 * (r - scalelow)/ (scalehigh - scalelow);
               else lut[r] = 65535 * (r - scalelow)/ (scalehigh - scalelow);
            }
         }
      }
   }
   else{
      if(!((scalebits == 16) && (scalelow == 3.60) && (scalehigh == 0.05))){
         lut = (unsigned short int *) calloc(65536, sizeof(unsigned short int));
         if(scalelow < scalehigh){
            fprintf(stderr, "Error! The larger valued optical density must be specified first!\n");
            fprintf(stderr, "You used -scale %f %f\n", scalelow, scalehigh);
            fprintf(stderr, "This would invert the image and you are not allowed to do that. Sorry!\n");
            exit(1);
         }
         if(VERBOSE){
            if(scalebits == 8) printf("OpticalDensity scaling: (%dbit) (%f, %f) -> (0,255)\n", scalebits, scalelow, scalehigh);
            else if(scalebits == 16) printf("OpticalDensity scaling: (%dbit) (%f, %f) -> (0, 65535)\n", scalebits, scalelow, scalehigh);
         }
         lowval = (float)od_to_scaled_od((double)scalelow);
         highval = (float)od_to_scaled_od((double)scalehigh);
         for(r=0;r<65536;r++){
            if(r >= highval){
               if(scalebits == 8) lut[r] = 255;
               else lut[r] = 65535;
            }
            else if(r < lowval) lut[r] = 0;
            else{
               if(scalebits == 8) lut[r] = 255 * (r - lowval)/ (highval - lowval);
               else lut[r] = 65535 * (r - lowval)/ (highval - lowval);
            }
         }
      }
   }

   if(lut != NULL){
      if(scalebits == 8){
         map_with_ushort_lut(&image, &tempimage, 12, (tempimage.rows * tempimage.cols *
         dt_size[UCHARNUM]) / 16, UCHARNUM, "rb", lut);
      }
      if(scalebits == 16){
         map_with_ushort_lut(&image, &tempimage, 12, (tempimage.rows * tempimage.cols *
         dt_size[USHORTNUM]) / 16, USHORTNUM, "rb", lut);
      }
      free(lut);
   }
   else image = tempimage;

   if(filemode == PGM){
      sprintf(image_filename, "%s.image.pgm", filename);
      if((fpim = fopen(image_filename, "wb")) == NULL){
         fprintf(stderr, "Error opening the image %s for writing.\n", image_filename);
         exit(1);
      }
      if(scalebits == 8){
         fprintf(fpim, "P5\n%d %d\n# Aggfactor %f\n255\n", image.cols, image.rows, aggfactor);
         thisrow_uchar = (unsigned char *) calloc(image.cols, sizeof(unsigned char));
         for(r=0;r<image.rows;r++){
            for(c=0;c<image.cols;c++)
               thisrow_uchar[c] = (unsigned char)  (*(image.getpixel))(&image, r, c);
            fwrite(thisrow_uchar, sizeof(unsigned char), image.cols, fpim);
         }
         fclose(fpim);
         free(thisrow_uchar);
      }
      if(scalebits == 16){
         fprintf(fpim, "P5\n%d %d\n# Aggfactor %f\n65535\n", image.cols, image.rows, aggfactor);
         thisrow_ushort = (unsigned short int *) calloc(image.cols, sizeof(unsigned short int));
         for(r=0;r<image.rows;r++){
            for(c=0;c<image.cols;c++)
               thisrow_ushort[c] = (unsigned short int)  (*(image.getpixel))(&image, r, c);
            fwrite(thisrow_ushort, sizeof(unsigned short int), image.cols, fpim);
         }
         fclose(fpim);
         free(thisrow_ushort);
      }
   }
   else{
      sprintf(image_filename, "%s.image.tif", filename);
      write_simple_tiff_CACHEIM(image_filename, &image);
   }
}
