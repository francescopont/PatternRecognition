/*******************************************************************************
* FILE: breastsegment.c
* Purpose: This code segments the breast region. It runs on a single image.
* To Compile: There is a makefile for compiling this program.
* Name: Micahel Heath, University of South Florida
* Date: 1/20/2000
* Notes: 6/9/2000 - I added the function to find the major axis in the breast.
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "breast.h"
#include "mikesfileio.h"
#include "myalloc.h"
#include "virtual_image.h"

#define VERSION "1.1.0"
#define VERSIONDATE "June 9, 2000"

#define LCC 0
#define LMLO 1
#define RCC 2
#define RMLO 3

#define PGM 1
#define TIFF 2

int VERBOSE = 0;

static void print_help();
static void get_command_line_arguments(int argc, char *argv[]);

/*******************************************************************************
* Function is in segment.c.
*******************************************************************************/
void segment_breast(char *filename, int rows, int cols, int bpp, int hb,
       int swap, int projection, int apex, float resolution, int mode,
       double threshold_double, double uthreshold_double, int *numpoints, float **xcoord, float **ycoord,
       int mapmode, double A, double B);

static int rows=0, cols=0, bpp=0, hb=0, swap=0, projection=0, mode=0, apex=0, mapmode=NOMAP;
static double threshold = 0.0, uthreshold=0.0;
static float resolution;
static double A=0.0, B=0.0;
static char *image_filename = NULL;
char *expanded[4] = {"LEFT_CC", "LEFT_MLO", "RIGHT_CC", "RIGHT_MLO"};
int dosketch=0;
float sketchresolution = 1000.0;      /* This is expressed in microns. */
int fileformat = PGM;

int main(int argc, char *argv[])
{
   float ax_0=0, ax_1=0, ay_0=0, ay_1=0;
   int numpoints=0;
   float *xcoord=NULL, *ycoord = NULL;
   char segmentation_filename[200];
   char sketch_segmentation_filename[200];

   void render_segmentation_sketch(char *filename, int numpoints, float *xc,
       float *yc, float resolution, int rows, int cols, float sketchres);
   void find_major_axis(int numpoints, float *x, float *y, float resolution,
      int projection, int apex, float *ax_0, float *ay_0, float *ax_1, float *ay_1);

   /****************************************************************************
   * Get the command line parameters from the user.
   ****************************************************************************/
   get_command_line_arguments(argc, argv);
   sprintf(segmentation_filename, "%s.sgt", image_filename);
   if(fileformat == PGM)
      sprintf(sketch_segmentation_filename, "%s.sgt.sketch.pgm", image_filename);
   else sprintf(sketch_segmentation_filename, "%s.sgt.sketch.tif", image_filename);

   /****************************************************************************
   * Segment the breast from the image. This returns a polygon that includes
   * the breast region in the image.
   ****************************************************************************/
   segment_breast(image_filename, rows, cols, bpp, hb, swap, projection, apex, resolution, mode, threshold,
      uthreshold, &numpoints, &xcoord, &ycoord, mapmode, A, B);

   /****************************************************************************
   * Find the major axis of the breast. The axis is specified as two endpoints
   * of a line. The point 0 should be near the nipple and point 1 should be
   * near the chest wall.
   ****************************************************************************/
   find_major_axis(numpoints, xcoord, ycoord, resolution, projection, apex, &ax_0, &ay_0, &ax_1, &ay_1);

   /****************************************************************************
   * Write out the breast boundary to a file.
   ****************************************************************************/
   write_segmentation_file(segmentation_filename, numpoints, xcoord, ycoord, ax_0, ay_0, ax_1, ay_1, rows, cols);

   /****************************************************************************
   * Render a small image showing the segmentation result.
   ****************************************************************************/
   if(dosketch) render_segmentation_sketch(sketch_segmentation_filename, numpoints,
                    xcoord, ycoord, resolution, rows, cols, sketchresolution);
}

/*******************************************************************************
* Function: render_segmentation_sketch
* Purpose: This function can render a small image (low resolution) of the breast
* segmentation results. The use of this image is mostly to provide feedback to
* the user that the breast segmentation result is reasonable. They can not
* easily tell that by just looking at the ascii segmentation file.
* Name: Michael Heath, University of South Florida
* Date: 1/21/2000
*******************************************************************************/
void render_segmentation_sketch(char *filename, int numpoints, float *xc,
    float *yc, float resolution, int rows, int cols, float sketchres)
{
   FILE *fp=NULL;
   int *xcoord=NULL, *ycoord=NULL;
   unsigned char *sketch_image = NULL;
   int n, sketch_rows, sketch_cols;
   float scale;

   int write_simple_tiff_uchar(char *filename, unsigned char *image, int rows, int cols);

   /****************************************************************************
   * We will want to shrink the polygon to a lower resolution.
   ****************************************************************************/
   xcoord = (int *) mycalloc("xcoord", numpoints, sizeof(int));
   ycoord = (int *) mycalloc("ycoord", numpoints, sizeof(int));

   scale = (resolution / sketchres);

   for(n=0;n<numpoints;n++){
      xcoord[n] = (int)floor(xc[n] * scale + 0.5);
      ycoord[n] = (int)floor(yc[n] * scale + 0.5);
   }

   /****************************************************************************
   * Allocate memory in which we will render the segmentation region.
   ****************************************************************************/
   sketch_rows = (int)ceil(rows * scale);
   sketch_cols = (int)ceil(cols * scale);
   if((sketch_image = (UCHAR *) mycalloc("sketch_image", sketch_rows * sketch_cols, sizeof(UCHAR))) == NULL){
      fprintf(stderr, "Calloc error of image in render_segmentation_sketch()!\n");
      exit(1);
   }
   polyscan_coords(numpoints, xcoord, ycoord, sketch_image, sketch_rows, sketch_cols);

   /****************************************************************************
   * Write the image to a file in TIFF format.
   ****************************************************************************/
   if(fileformat == TIFF){
      write_simple_tiff_uchar(filename, sketch_image, sketch_rows, sketch_cols);
   }

   /****************************************************************************
   * Write the image to a file in PGM format.
   ****************************************************************************/
   else{
      if((fp = fopen(filename, "wb")) == NULL){
         fprintf(stderr, "Error opening the file %s for writing in render_segmentation_sketch().\n",
            filename);
         exit(1);
      }
      fprintf(fp, "P5\n%d %d\n255\n", sketch_cols, sketch_rows);
      fwrite(sketch_image, sizeof(unsigned char), sketch_rows * sketch_cols, fp);
      fclose(fp);
   }

   myfree("xcoord", xcoord);
   myfree("xcoord", ycoord);
   myfree("sketch_image", sketch_image);
}

/*******************************************************************************
* Function: get_command_line_arguments
* Purpose: This function is self explanitory. It just extracts the command line
* arguments.
* Name: Michael Heath, University of South Florida
* Date: 7/22/99
*******************************************************************************/
static void get_command_line_arguments(int argc, char *argv[])
{
   int i;
   char *ics_filename = NULL;
   int view;

   void get_ddsm_image_info(char *ics_filename, int view, char **filename,
      int *rows, int *cols, int *bpp, float *resolution, int *projection, int *apex);

   rows = 0;
   cols = 0;
   bpp = 0;
   image_filename = NULL;
   apex = -1;
   projection = -1;
   resolution = 0.0;
   hb = 0;

   mode = DEFAULTMODE;
   resolution = ASSUMED_RESOLUTION;
   swap = 0;

   for(i=1;i<argc;i++){
      if(strcmp(argv[i], "-i") == 0){ image_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-version") == 0){
         printf("\n\n%s Version: %s %s\n\n", argv[0], VERSION, VERSIONDATE);
         exit(1);
      }
      else if(strcmp(argv[i], "-inp") == 0){ cols = atoi(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-inl") == 0){ rows = atoi(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-hb") == 0){ hb = atoi(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-bpp") == 0){ bpp = atoi(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-microns") == 0){ resolution = atof(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-swap") == 0){ swap = 1; }
      else if(strcmp(argv[i], "-v") == 0){ VERBOSE = 1; }
      else if(strncmp(argv[i], "-tiff", 4) == 0){ fileformat = TIFF; }
      else if(strcmp(argv[i], "-sketchres") == 0){
         dosketch = 1;
         sketchresolution = atof(argv[i+1]);
         i++;
       }
      else if(strcmp(argv[i], "-projection") == 0){
         if(strcmp(argv[i+1], "MLO") == 0) projection = MLO;
         if(strcmp(argv[i+1], "CC") == 0) projection = CC;
         i++;
      }
      else if(strcmp(argv[i], "-apex") == 0){
         if(strcmp(argv[i+1], "LEFT") == 0) apex = APEXLEFT;
         if(strcmp(argv[i+1], "RIGHT") == 0) apex = APEXRIGHT;
      }
      else if(strcmp(argv[i], "-auto") == 0) mode = AUTOMODE;
      else if(strcmp(argv[i], "-thresh") == 0){
         mode = MANUALMODE;
         threshold = atof(argv[i+1]);
         i++;
      }
      else if(strcmp(argv[i], "-uthresh") == 0){
         mode = MANUALMODE;
         uthreshold = atof(argv[i+1]);
         i++;
      }
      else if(strcmp(argv[i], "-ics") == 0){ ics_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-view") == 0){
         if(strcmp(argv[i+1], "LEFT_MLO") == 0) view = LMLO;
         if(strcmp(argv[i+1], "RIGHT_MLO") == 0) view = RMLO;
         if(strcmp(argv[i+1], "LEFT_CC") == 0) view = LCC;
         if(strcmp(argv[i+1], "RIGHT_CC") == 0) view = RCC;
         i++;
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

   /****************************************************************************
   * If the user specified the auto mode they can not have specified either
   * threshold value.
   ****************************************************************************/
   if((mode == AUTOMODE) && (!((threshold == 0.0) && (uthreshold == 0.0)))){
      fprintf(stderr, "\nError! You can not specify both auto and a threshold or uthreshold value.\n");
      exit(1);
   }

   if((rows == 0) || (cols == 0) || (bpp == 0) || (image_filename == NULL) || (apex == -1) ||
      (projection == -1) || (resolution == 0.0)) print_help();
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
*
**********************************
* Example contents of an ics file.
**********************************
* ics_version 1.0
* filename C-0234-1
* DATE_OF_STUDY 22 12 1994
* PATIENT_AGE 40
* FILM 
* FILM_TYPE REGULAR
* DENSITY 2
* DATE_DIGITIZED 28 4 1998 
* DIGITIZER LUMISYS LASER
* SEQUENCE
* LEFT_CC LINES 4736 PIXELS_PER_LINE 2240 BITS_PER_PIXEL 12 RESOLUTION 50 NON_OVERLAY
* LEFT_MLO LINES 4760 PIXELS_PER_LINE 2680 BITS_PER_PIXEL 12 RESOLUTION 50 NON_OVERLAY
* RIGHT_CC LINES 4760 PIXELS_PER_LINE 2512 BITS_PER_PIXEL 12 RESOLUTION 50 OVERLAY
* RIGHT_MLO LINES 4688 PIXELS_PER_LINE 2608 BITS_PER_PIXEL 12 RESOLUTION 50 OVERLAY
*******************************************************************************/
void get_ddsm_image_info(char *ics_filename, int view, char **filename, int *rows, int *cols,
   int *bpp, float *resolution, int *projection, int *apex)
{
   FILE *fp=NULL;
   char line[200], caseid[20];
   int bits, i;
   char tmpfilename[100] = {'\0'};

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

      line[199] = '\0';

      if(strncmp(line, "filename", strlen("filename")) == 0){
         sscanf(line, "%*s %s", caseid);

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

         for(i=(strlen(tmpfilename)-8);i<strlen(tmpfilename);i++)
            if(tmpfilename[i] == '-') tmpfilename[i] = '_';

         sprintf(*filename, "%s.%s.LJPEG.1", tmpfilename, expanded[view]);

      }

      else if(strncmp(line, expanded[view], strlen(expanded[view])) == 0){
         sscanf(line, "%*s %*s %d %*s %d %*s %d %*s %f", rows, cols, &bits, resolution);
         if(VERBOSE){
            printf("The image file we will process has the following attributes:\n");
            printf("rows = %d\n", *rows);
            printf("cols = %d\n", *cols);
            printf("bits = %d\n", bits);
            printf("resolution = %f\n", *resolution);
         }
         *bpp = (int)ceil((double)bits / 8.0);
      }
   }
}

/*******************************************************************************
* Function: print_help
* Purpose: This function is self explanitory. It just prints out information.
* Name: Michael Heath, University of South Florida
* Date: 7/22/99
*******************************************************************************/
static void print_help()
{

   printf("\n********************************************************************************\n");
   printf("This program performs a segmentation of the breast from the background. The\n");
   printf("program will output a text file containing this information.\n");
   printf("********************************************************************************\n\n");
   printf("<USAGE> breastsegment [-ics filename.ics -view VIEW] |\n");
   printf("                      [-i filename -inp # -inl # -bpp BPP -microns # -hb #\n");
   printf("                      -projection PROJECTION -apex APEX] [-swap] [-v]\n");
   printf("                      [-linear # #] | [-log10 # #] [-sketchres #] [-version]\n");
   printf("                      [-auto] | [-thresh # -uthresh #] [-tif]\n\n");
   printf("   -ics        The name of an ics file from the DDSM database. This supplies the\n");
   printf("               image filename, view, inp, inl, bpp, micons, hb, projection and\n");
   printf("               apex automatically.\n");
   printf("   -view       The mammogram to process for this case. VIEW can be any one of\n");
   printf("               the following: LEFT_MLO,RIGHT_MLO,LEFT_CC or RIGHT_CC.\n");
   printf("   -i          The filename of the image to process.\n");
   printf("   -inp        The number of columns in the image.\n");
   printf("   -inl        The number of rows in the image.\n");
   printf("   -bpp        The number of bytes per pixel. BPP can be either 1 or 2.\n");
   printf("   -microns    The sampling resolution of the image in microns.\n");
   printf("   -swap       If this flag is used the bytes of a 16-bit image will be swapped\n");
   printf("               when image data is read from a file.\n");
   printf("   -hb         The number of header bytes in the image file before image data\n");
   printf("               is found.\n");
   printf("   -projection The mammographic projection of this image. PROJECTION can be\n");
   printf("               either MLO or CC.\n");
   printf("   -apex       The side of the image that the nipple is facing. APEX can be\n");
   printf("               either LEFT or RIGHT.\n");
   printf("   -auto       Find the threshold values automatically (not tested well!).\n");
   printf("   -thresh     A pixel value or optical density to threshold the breast region\n");
   printf("               from the background in the image.\n");
   printf("   -uthresh    A pixel value or optical density to threshold the breast from the\n");
   printf("               possible bright white border in the image (where there was no\n");
   printf("               screen or the film was blocked from exposure).\n");
   printf("   -linear     Rescale to optical density using a linear function of the form\n");
   printf("               optical_density = A + (B * GL). A and B are specified as double\n");
   printf("               precision values with a space in between them.\n");
   printf("   -log10      Rescale to optical density using a log function.\n");
   printf("               optical_density = A + (B * log10(GL)). A and B are specified as\n");
   printf("               double precision values with a space in between them.\n");
   printf("   -sketchres  Output an image that is a sketch of the segmented breast region\n");
   printf("               at the resolution specified in microns (0=background,255=breast).\n");
   printf("   -tif        If the sketch image is to be written, write it in TIFF format\n");
   printf("               rather than in PGM format.\n");
   printf("   -v          Run the program in verbose mode.\n");
   printf("   -version    Print the version and exit.\n\n");

   exit(1);
}

/*******************************************************************************
* Function: find_major_axis
* Purpose: The purpose of this function is to locate the major axis of the
* breast. Ideally this is a vector positioned at the base of the breast,
* perpendicular to the pectoral muscle that projects through the nipple.
* The local shape of the breast border is used to identify the major axis.
* The chest wall and nipple position are never explicitly located.
* Name: Mike Heath
* Date: 3/29/99
*******************************************************************************/
void find_major_axis(int numpoints, float *x, float *y, float resolution,
   int projection, int apex, float *ax_0, float *ay_0, float *ax_1, float *ay_1)
{
   int p, bestp, pixel_diameter, s1, s2, k, beginp, endp, begink, endk, i, j;
   double x0, y0, x1, y1, x2, y2, dist, mindist;
   double q0_x, q0_y, q1_x, q1_y, o_x, o_y, dx, dy, a, len;
   double p0_x, p0_y, p1_x, p1_y, p2_x, p2_y, p3_x, p3_y, p4_x, p4_y, ua, ub, min_ua;
   int p2, p3, p4, numintersections, pp = -1;
   double simm_a, simm_b, simm, bestsimm;
   double border_distance = 10.0;
   double surround_distance;
   int skipit, scp, ecp;

   /****************************************************************************
   * Compute the distance all the way around the boundary.
   ****************************************************************************/
   surround_distance = 0.0;

   for(p=0;p<numpoints;p++){
      if(p != 0) surround_distance += sqrt(((x[p]-x[p-1])*(x[p]-x[p-1]))+((y[p]-y[p-1])*(y[p]-y[p-1])));
   }
   surround_distance += sqrt(((x[p-1]-x[0])*(x[p-1]-x[0]))+((y[p-1]-y[0])*(y[p-1]-y[0])));

   if(VERBOSE) printf("The surround_distance is %f pixels.\n", surround_distance);
   if(VERBOSE) printf("This is %f cm.\n", surround_distance * resolution / 10000.0);

   /****************************************************************************
   * Set the distance that we will look (to the left and to the right) to
   * try an axis. The distance is specified here in cm and it is computed as
   * a fixed fraction of the distance around the segmented breast tissue.
   ****************************************************************************/
   if(projection == CC) border_distance = (surround_distance * resolution / 10000.0) / 4.0;
   if(projection == MLO) border_distance = (surround_distance * resolution / 10000.0) / 6.0;

   while((pp == -1) && (border_distance >= 3.0)){

      if(VERBOSE) printf("Using border_distance = %f cm.\n", border_distance);

      /*************************************************************************
      * Go through all of the points on the boundry testing each for the measure
      * of how good an axis would be through the point.
      *************************************************************************/
      for(bestp=0;bestp<numpoints;bestp++){

         /*************************************************************************
         * Find the beginning index of a border_distance cm portion of the chain.
         *************************************************************************/
         dist = 0.0;
         k = 0;
         s1 = bestp;
         x1 = x[s1] * resolution / 1000.0; /* gets us mm */
         y1 = y[s1] * resolution / 1000.0; /* gets us mm */
         do{
            s2 = (s1 - 1 + numpoints) % numpoints;
 
            x2 = x[s2] * resolution / 1000.0; /* gets us mm */
            y2 = y[s2] * resolution / 1000.0; /* gets us mm */
 
            dist += sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
 
            s1 = s2;
            x1 = x2;
            y1 = y2;
            k++;
         }while((dist < (10.0*border_distance)) && (k < 10000));  /* borderdistance cm */
         begink = -k;

         if(k < 10000) beginp = s1;

         /*************************************************************************
         * Find the ending index of a border_distance cm portion of the chain.
         *************************************************************************/
         dist = 0.0;
         k = 0;
         s1 = bestp;
         x1 = x[s1] * resolution / 1000.0; /* gets us mm */
         y1 = y[s1] * resolution / 1000.0; /* gets us mm */
         do{
            s2 = (s1 + 1 + numpoints) % numpoints;
 
            x2 = x[s2] * resolution / 1000.0; /* gets us mm */
            y2 = y[s2] * resolution / 1000.0; /* gets us mm */
 
            dist += sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
 
            s1 = s2;
            x1 = x2;
            y1 = y2;
            k++;
         }while((dist < (10.0*border_distance)) && (k < 10000));  /* border_distance cm */
         endk = k;

         if(k < 10000) endp = s1;

         /**********************************************************************
         * It was observed that sometimes the line connecting the points found
         * forwards and backwards along the chain (identified by beginp and endp)
         * was not completely inside the breast. We should limit the positions
         * that are found to those where this line is completely inside the
         * breast. Here we will check to make sure the line is completely
         * contained in the breast. This is accomplished by checking all of
         * the line segments in the range beginp to end p and to see that none
         * of them intersect the line from beginp to endp.
         **********************************************************************/

/* BEGINNING OF NEW SECTION */

         p1_x = x[bestp];
         p1_y = y[bestp];

         p2_x = (x[beginp] + x[endp]) / 2.0;
         p2_y = (y[beginp] + y[endp]) / 2.0;

         skipit = 0;
         for(i=begink;i<endk;i++){

            p3 = (bestp+numpoints+i)%numpoints;

            p3_x = x[p3];
            p3_y = y[p3];

            if(((p3_x - p1_x) * (p2_x - p1_x) + (p3_y - p1_y) * (p2_y - p1_y)) < 0.0){
               skipit = 1;
               break;
            }
         }

         if(skipit != 0) continue;

/*
         p1_x = x[beginp];
         p1_y = y[beginp];

         p2_x = x[endp];
         p2_y = y[endp];

         numintersections = 0;

         for(i=(begink+1);i<0;i++){
            p3 = (bestp+numpoints+i)%numpoints;
            p4 = (bestp+numpoints+i+1)%numpoints;

            p3_x = x[p3];
            p3_y = y[p3];

            p4_x = x[p4];
            p4_y = y[p4];

            ua = ((p4_x-p3_x)*(p1_y-p3_y) - (p4_y-p3_y)*(p1_x-p3_x)) /
                 ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

            ub = ((p2_x-p1_x)*(p1_y-p3_y) - (p2_y-p1_y)*(p1_x-p3_x)) /
                 ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

            if((ub>=0.0) && (ub<=1.0)){
               numintersections++;
	       if(numintersections == 1) min_ua = ua;
               else{
                  if(fabs(ua) < fabs(min_ua)) min_ua = ua;
               }
            }
         }

         for(i=0;i<(endk-1);i++){
            p3 = (bestp+numpoints+i)%numpoints;
            p4 = (bestp+numpoints+i+1)%numpoints;

            p3_x = x[p3];
            p3_y = y[p3];

            p4_x = x[p4];
            p4_y = y[p4];

            ua = ((p4_x-p3_x)*(p1_y-p3_y) - (p4_y-p3_y)*(p1_x-p3_x)) /
                 ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

            ub = ((p2_x-p1_x)*(p1_y-p3_y) - (p2_y-p1_y)*(p1_x-p3_x)) /
                 ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

            if((ub>=0.0) && (ub<=1.0)){
               numintersections++;
	       if(numintersections == 1) min_ua = ua;
               else{
                  if(fabs(ua) < fabs(min_ua)) min_ua = ua;
               }
            }
         }

         if(numintersections != 0) continue;
*/

         /**********************************************************************
         * Find a point (p2_x,p2_y) that is 0.5cm along the (possible) axis from
         * the origin to the midpoint of the distances travelled along the
         * chain.
         **********************************************************************/
/*
         p1_x = x[bestp];
         p1_y = y[bestp];

         p2_x = (x[beginp] + x[endp]) / 2.0;
         p2_y = (y[beginp] + y[endp]) / 2.0;

         dist = sqrt((p2_x-p1_x)*(p2_x-p1_x)+(p2_y-p1_y)*(p2_y-p1_y));

         dx = (p2_x - p1_x) / dist;
         dy = (p2_y - p1_y) / dist;

         p2_x = p1_x + 5000.0 / resolution;
         p2_y = p1_y + 5000.0 / resolution;

         p1_x = p2_x;
         p1_y = p2_y;

         numintersections = 0;

         for(i=begink;i<=endk;i++){
            p2 = (bestp+numpoints+i)%numpoints;

            p2_x = x[p2];
            p2_y = y[p2];

            if(i < 0){
               scp = i+1;
               ecp = 0;
            }
            else{
               scp = 0;
               ecp = i-1;
            }

            for(j=scp;j<ecp;j++){

               p3 = (bestp+numpoints+j)%numpoints;
               p4 = (bestp+numpoints+j+1)%numpoints;

               p3_x = x[p3];
               p3_y = y[p3];

               p4_x = x[p4];
               p4_y = y[p4];

               ua = ((p4_x-p3_x)*(p1_y-p3_y) - (p4_y-p3_y)*(p1_x-p3_x)) /
                    ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));
   
               ub = ((p2_x-p1_x)*(p1_y-p3_y) - (p2_y-p1_y)*(p1_x-p3_x)) /
                    ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

               if((ub>=0.0) && (ub<=1.0)) numintersections++;
            }
         }

         if(numintersections != 0) continue;
*/

/* */

/* END OF NEW SECTION */

         /*************************************************************************
         * Find the intersection of the stepped line with the border.
         *************************************************************************/
         q0_x = x[bestp];
         q0_y = y[bestp];

         q1_x = (x[beginp]+x[endp])/2.0;
         q1_y = (y[beginp]+y[endp])/2.0;

         if(apex == APEXRIGHT){
            if(q0_x < q1_x) continue;
         }
         else{
            if(q0_x > q1_x) continue;
         }
         if(projection == MLO){
            double normalizer;

	    if(q1_y > q0_y) continue;

            normalizer = sqrt((q1_x - q0_x) * (q1_x - q0_x) + (q1_y - q0_y) * (q1_y - q0_y));
            if(fabs((q1_y - q0_y)/normalizer) > fabs((q1_x - q0_x)/normalizer)) continue;
         }
         else{
            double normalizer;

            normalizer = sqrt((q1_x - q0_x) * (q1_x - q0_x) + (q1_y - q0_y) * (q1_y - q0_y));
            if(2*fabs((q1_y - q0_y)/normalizer) > fabs((q1_x - q0_x)/normalizer)) continue;

         }

         /* if(VERBOSE) printf("q0_x = %f q0_y = %f\n", q0_x, q0_y); */
         /* if(VERBOSE) printf("q1_x = %f q1_y = %f\n", q1_x, q1_y); */

         dx = x[endp] - x[beginp];
         dy = y[endp] - y[beginp];

         simm = 0.0;

         for(a=0.2;a<=1.0;a+=0.2){

            o_x = q0_x + a * (q1_x - q0_x);
            o_y = q0_y + a * (q1_y - q0_y);

            /* if(VERBOSE) printf("o_x = %f o_y = %f\n", o_x, o_y); */

            p1_x = o_x;
            p1_y = o_y;
   
            p2_x = p1_x + dx;   /* Any non-zero multiple of dx would work fine. */
            p2_y = p1_y + dy;   /* Any non-zero multiple of dy would work fine. */

            numintersections = 0;
            for(i=begink;i<0;i++){
               p3 = (bestp+numpoints+i)%numpoints;
               p4 = (bestp+numpoints+i+1)%numpoints;

               p3_x = x[p3];
               p3_y = y[p3];

               p4_x = x[p4];
               p4_y = y[p4];

               ua = ((p4_x-p3_x)*(p1_y-p3_y) - (p4_y-p3_y)*(p1_x-p3_x)) /
                    ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

               ub = ((p2_x-p1_x)*(p1_y-p3_y) - (p2_y-p1_y)*(p1_x-p3_x)) /
                    ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

               if((ub>=0.0) && (ub<=1.0)){
                  numintersections++;
	          if(numintersections == 1) min_ua = ua;
                  else{
                     if(fabs(ua) < fabs(min_ua)) min_ua = ua;
                  }
               }
            }

            if(numintersections != 0) simm_a = min_ua;

            numintersections = 0;
            for(i=0;i<endk;i++){
               p3 = (bestp+numpoints+i)%numpoints;
               p4 = (bestp+numpoints+i+1)%numpoints;

               p3_x = x[p3];
               p3_y = y[p3];

               p4_x = x[p4];
               p4_y = y[p4];

               ua = ((p4_x-p3_x)*(p1_y-p3_y) - (p4_y-p3_y)*(p1_x-p3_x)) /
                    ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

               ub = ((p2_x-p1_x)*(p1_y-p3_y) - (p2_y-p1_y)*(p1_x-p3_x)) /
                    ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

               if((ub>=0.0) && (ub<=1.0)){
                  numintersections++;
	          if(numintersections == 1) min_ua = ua;
                  else{
                     if(fabs(ua) < fabs(min_ua)) min_ua = ua;
                  }
               }
            }

            if(numintersections != 0) simm_b = min_ua;

            numintersections = 0;


            simm += (fabs(simm_a) - fabs(simm_b)) * (fabs(simm_a) - fabs(simm_b));

            /* simm += ((simm_a-simm_b)*(simm_a-simm_b)) / ((simm_a+simm_b)*(simm_a+simm_b)); */

         }

         if(sqrt((q0_x-q1_x)*(q0_x-q1_x)+(q0_y-q1_y)*(q0_y-q1_y)) >= 40){
            simm /= sqrt((q0_x - q1_x) * (q0_x - q1_x) + (q0_y - q1_y) * (q0_y - q1_y));
         }
         else simm = 1.0e6; 

         if(projection == CC){
            simm /= ((q0_x - q1_x)*(q0_x - q1_x)) / ((q0_x - q1_x) * (q0_x - q1_x) + (q0_y - q1_y) * (q0_y - q1_y));
         }

         /* if(VERBOSE) printf("The simm measure is %f\n", simm); */

         if(pp == -1){
            bestsimm = simm;
            pp = bestp;
         }
         else{
            if(simm < bestsimm){
               bestsimm = simm;
               pp = bestp;
            }
         }
      }

      if(pp == -1) border_distance -= 1.0;

   }

   if(pp == -1){
      if(VERBOSE) printf("No suitable axis was found!\n");
      myfree("x", x);
      myfree("y", y);
      return;
   }
   bestp = pp;

   /****************************************************************************
   * Find the beginning index of a border_distance cm portion of the chain.
   ****************************************************************************/
   dist = 0.0;
   k = 0;
   s1 = bestp;
   x1 = x[s1] * resolution / 1000.0; /* gets us mm */
   y1 = y[s1] * resolution / 1000.0; /* gets us mm */
   do{
      s2 = (s1 - 1 + numpoints) % numpoints;
 
      x2 = x[s2] * resolution / 1000.0; /* gets us mm */
      y2 = y[s2] * resolution / 1000.0; /* gets us mm */
 
      dist += sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
 
      s1 = s2;
      x1 = x2;
      y1 = y2;
      k++;
   }while((dist < (10.0*border_distance)) && (k < 10000));  /* border_distance cm */
   begink = -k;

   if(k < 10000) beginp = s1;

   /****************************************************************************
   * Find the ending index of a border_distance cm portion of the chain.
   ****************************************************************************/
   dist = 0.0;
   k = 0;
   s1 = bestp;
   x1 = x[s1] * resolution / 1000.0; /* gets us mm */
   y1 = y[s1] * resolution / 1000.0; /* gets us mm */
   do{
      s2 = (s1 + 1 + numpoints) % numpoints;
 
      x2 = x[s2] * resolution / 1000.0; /* gets us mm */
      y2 = y[s2] * resolution / 1000.0; /* gets us mm */
 
      dist += sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
 
      s1 = s2;
      x1 = x2;
      y1 = y2;
      k++;
   }while((dist < (10.0*border_distance)) && (k < 10000));  /* border_distance cm */
   endk = k;

   if(k < 10000) endp = s1;

   /****************************************************************************
   * Find the intersection of the stepped line with the border.
   ****************************************************************************/
   if(VERBOSE) printf("bestp = %d, begin_p = %d, end_p = %d\n", bestp, beginp, endp);
   if(VERBOSE) printf("begink = %d, endk = %d\n", begink, endk);

   if(VERBOSE){
      printf("bestp = %d, begin_p = %d, end_p = %d\n", bestp, beginp, endp);
      printf("The origin is at (%f, %f)\n", x[bestp]/7, y[bestp]/7);
      printf("The beginning is at (%f, %f)\n", x[beginp]/7, y[beginp]/7);
      printf("The end is at (%f, %f)\n", x[endp]/7, y[endp]/7);
   }

   simm = 0.0;

   q0_x = x[bestp];
   q0_y = y[bestp];

   q1_x = (x[beginp]+x[endp])/2.0;
   q1_y = (y[beginp]+y[endp])/2.0;

   if(VERBOSE) printf("q0_x = %f q0_y = %f\n", q0_x, q0_y);
   if(VERBOSE) printf("q1_x = %f q1_y = %f\n", q1_x, q1_y);

   dx = x[endp] - x[beginp];
   dy = y[endp] - y[beginp];

   for(a=0.2;a<=1.0;a+=0.2){

      o_x = q0_x + a * (q1_x - q0_x);
      o_y = q0_y + a * (q1_y - q0_y);

      if(VERBOSE) printf("o_x = %f o_y = %f\n", o_x, o_y);

      p1_x = o_x;
      p1_y = o_y;
   
      p2_x = p1_x + dx;   /* Any non-zero multiple of dx would work fine. */
      p2_y = p1_y + dy;   /* Any non-zero multiple of dy would work fine. */

      numintersections = 0;
      for(i=begink;i<0;i++){
         p3 = (bestp+numpoints+i)%numpoints;
         p4 = (bestp+numpoints+i+1)%numpoints;

         p3_x = x[p3];
         p3_y = y[p3];

         p4_x = x[p4];
         p4_y = y[p4];

         ua = ((p4_x-p3_x)*(p1_y-p3_y) - (p4_y-p3_y)*(p1_x-p3_x)) /
              ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

         ub = ((p2_x-p1_x)*(p1_y-p3_y) - (p2_y-p1_y)*(p1_x-p3_x)) /
              ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

         if((ub>=0.0) && (ub<=1.0)){
            numintersections++;
	    if(numintersections == 1) min_ua = ua;
            else{
               if(fabs(ua) < fabs(min_ua)) min_ua = ua;
            }
         }
      }

      if(numintersections != 0) simm_a = min_ua;

      numintersections = 0;
      for(i=0;i<endk;i++){
         p3 = (bestp+numpoints+i)%numpoints;
         p4 = (bestp+numpoints+i+1)%numpoints;

         p3_x = x[p3];
         p3_y = y[p3];

         p4_x = x[p4];
         p4_y = y[p4];

         ua = ((p4_x-p3_x)*(p1_y-p3_y) - (p4_y-p3_y)*(p1_x-p3_x)) /
              ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

         ub = ((p2_x-p1_x)*(p1_y-p3_y) - (p2_y-p1_y)*(p1_x-p3_x)) /
              ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

         if((ub>=0.0) && (ub<=1.0)){
            numintersections++;
	    if(numintersections == 1) min_ua = ua;
            else{
               if(fabs(ua) < fabs(min_ua)) min_ua = ua;
            }
         }
      }

      if(numintersections != 0) simm_b = min_ua;

      numintersections = 0;

      simm += (fabs(simm_a) - fabs(simm_b)) * (fabs(simm_a) - fabs(simm_b));
   }

   simm /= sqrt((q0_x - q1_x) * (q0_x - q1_x) + (q0_y - q1_y) * (q0_y - q1_y));

   /* if(VERBOSE) printf("The simm measure is %f\n", simm); */

   /****************************************************************************
   * Draw the axis on the image with a fixed scale. This will look like a ruler.
   ****************************************************************************/
   p1_x = q0_x;
   p1_y = q0_y;

   p2_x = q1_x;
   p2_y = q1_y;

   dx = p2_x - p1_x;
   dy = p2_y - p1_y;

   numintersections = 0;
   for(i=1;i<(numpoints-2);i++){
      p3 = (bestp+numpoints+i)%numpoints;
      p4 = (bestp+numpoints+i+1)%numpoints;

      p3_x = x[p3];
      p3_y = y[p3];

      p4_x = x[p4];
      p4_y = y[p4];

      ua = ((p4_x-p3_x)*(p1_y-p3_y) - (p4_y-p3_y)*(p1_x-p3_x)) /
           ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

      ub = ((p2_x-p1_x)*(p1_y-p3_y) - (p2_y-p1_y)*(p1_x-p3_x)) /
           ((p4_y-p3_y)*(p2_x-p1_x) - (p4_x-p3_x)*(p2_y-p1_y));

      if((ub>=0.0) && (ub<=1.0)){
         numintersections++;
         if(numintersections == 1) min_ua = ua;
         else{
            if(fabs(ua) < fabs(min_ua)) min_ua = ua;
         }
      }
   }

   q0_x = x[bestp]; /* This is the origin of the axis. */
   q0_y = y[bestp]; /* This is the origin of the axis. */

   q1_x = p1_x+min_ua*dx; /* This is the point of intersection near the base of the breast. */
   q1_y = p1_y+min_ua*dy; /* This is the point of intersection near the base of the breast. */

   if(VERBOSE) printf("q0_x = %f q0_y = %f\n", q0_x, q0_y);
   if(VERBOSE) printf("q1_x = %f q1_y = %f\n", q1_x, q1_y);

   dx = q1_x - q0_x;
   dy = q1_y - q0_y;
   len = sqrt(dx*dx+dy*dy);

   if(VERBOSE) printf("min_ua = %f len = %f\n", min_ua, len);

   a = (10 * 1000.0 / resolution) / 2.0;
   while(a < len){

      o_x = q0_x + a * dx/len;
      o_y = q0_y + a * dy/len;

      a += 10 * 1000.0 / resolution;
   }

   /****************************************************************************
   * Place the major axis that was located into the breast coord data structure.
   ****************************************************************************/
   *ax_0 = q0_x;
   *ay_0 = q0_y;
   *ax_1 = q1_x;
   *ay_1 = q1_y;
}
