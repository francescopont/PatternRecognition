/*******************************************************************************
* Program: drawimage
* Purpose: This program is a tool for creating an image showing the
* breast boundary from the segmentation file, the valid ground truth
* regions and the selected detections on an image.
*
* The program is run by specifying the needed files on the command line, and
* an optional option. The option "-l lesiontype" can be specified as either
* "-l CALCIFICATION" or "-l MASS". One of these will limit the ground truth
* regions to the specified type. This allows the rendering of regions of
* the specified type. The default is to render regions of both types.
*
* By default only regions that have a PATHOLOGY of MALIGNANT or BENIGN are
* used. All other ground truth regions are skipped. If one uses the optional
* -all flag, every region is scan converted and will be included in the
* template image that is created.
*
* Name: Michael Heath, University of South Florida
* Date: 1/21/2000
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ICSio.h"
#include "overlay.h"
#include "mikesfileio.h"

#define VERSION "1.0.0"
#define VERSIONDATE "January 26, 2000"

/*******************************************************************************
* Function prototypes.
*******************************************************************************/
static void print_help();
static int get_commandline_parameters(int argc, char *argv[], char **ics_filename,
   char **view, char *lesiontype, char **pathology, char **ms, char **mm,
   char **ct, char **cd, char **detection_filename, int *detmax, float *detthresh,
   int *use_all_regions, char **segmentation_filename, char **image_filename);
void drawimageline(unsigned char *image, int rows, int cols, int x0, int y0, int x1, int y1, int valtoplot);

char created_files[4][200];  /* This is used by read_ics_file(). */

int VERBOSE = 0;

int main(int argc, char *argv[])
{
   FILE *fp;
   char *ics_filename=NULL, *image_filename, *segmentation_filename, *detection_filename,
      *overlay_filename=NULL, draw_filename[200];
   char lesiontype[20], *view = NULL;
   int status;
   ICSDATA icsdata;
   OVERLAY_DATA overlay_data;
   int d, r, c, t, rr, cc, i;
   ABNORMALITY abnormality;
   REGION *truth_region=NULL;
   int truth_regions_valid = 0;
   char *ms = NULL, *mm = NULL, *ct = NULL, *cd = NULL;
   float resolution, desired_resolution, scale;
   int rows, cols, template_rows, template_cols, cp, gtpix, overlapping_regions=0;
   unsigned char *template=NULL;
   char template_filename[200] = {'\0'};
   int use_all_regions=0, use_labels=0, label;
   int *isvalid=NULL;
   int detmax;
   float detthresh;
   unsigned char *image = NULL;
   int numpoints;
   float *xcoord=NULL, *ycoord = NULL;
   float startx, starty, endx, endy;
   int segmentation_rows, segmentation_cols;
   int num_detections = 0;
   REGION *detections=NULL;
   int detection_rows, detection_cols;
   unsigned char *drawim=NULL, *redim=NULL, *grnim=NULL, *bluim=NULL;
   int overlay_rows, overlay_cols, pos;
   char *pathology=NULL;
   struct COLOR{
      int r, g, b;
   }color[256];

   int REDINDEX=1, GREENINDEX=2, BLUEINDEX=3, YELLOWINDEX=4;

   int read_pgm_image(char *infilename, unsigned char **image, int *rows, int *cols);
   int write_pgm_image(char *outfilename, unsigned char *image, int rows,
       int cols, char *comment, int maxval);
   void draw_segmentation(unsigned char *image, int rows, int cols, int numpoints,
       float *xcoord, float *ycoord, float startx, float starty,
       float endx, float endy, int grayval);
   void draw_detections(unsigned char *image, int rows, int cols, int num,
      REGION *detections, int grayval);
   int write_ppm_image(char *outfilename, unsigned char *image_red,
      unsigned char *image_grn, unsigned char *image_blu, int rows,
     int cols, char *comment, int maxval);
   void draw_valid_abnormalities(unsigned char *image, int rows, int cols,
       OVERLAY_DATA *overlay_data, int overlay_rows, int overlay_cols,
       int *isvalid, int grayval);

   memset(&overlay_data, 0, sizeof(OVERLAY_DATA));

   REDINDEX = 1;
   GREENINDEX = 2;
   BLUEINDEX = 3;
   YELLOWINDEX = 4;
   color[REDINDEX].r    = 255;  color[REDINDEX].g    =   0;  color[REDINDEX].b    =   0;
   color[GREENINDEX].r  =   0;  color[GREENINDEX].g  = 255;  color[GREENINDEX].b  =   0;
   color[BLUEINDEX].r   =   0;  color[BLUEINDEX].g   =   0;  color[BLUEINDEX].b   = 255;
   color[YELLOWINDEX].r = 255;  color[YELLOWINDEX].g = 255;  color[YELLOWINDEX].b =   0;

   /****************************************************************************
   * Get the command line parameters.
   ****************************************************************************/
   status = get_commandline_parameters(argc, argv, &ics_filename,
               &view, lesiontype, &pathology, &ms, &mm, &ct, &cd,
               &detection_filename, &detmax, &detthresh, &use_all_regions,
               &segmentation_filename, &image_filename);
   if(status == 0) exit(1);

   /****************************************************************************
   * If we are running in a verbose mode, print out the command line info.
   ****************************************************************************/
   if(VERBOSE){
      printf("\n\n************************************************************\n");
      printf(" The drawimage program is running in verbose mode.\n");
      printf("************************************************************\n\n");
      printf("   image.pgm   : %s\n", image_filename);
      if(ics_filename != NULL)
         printf("   file.ics    : %s\n", ics_filename);
      if(view != NULL)
         printf("   view        : %s\n", view);
      if(use_all_regions == 1){
         printf("   You selected the option to use all regions.\n");
      }
     
   }

   /****************************************************************************
   * Read the RAW 8-bit pgm image from the file.
   ****************************************************************************/
   if(read_pgm_image(image_filename, &image, &rows, &cols) == 0) exit(1);

   drawim = (unsigned char *)calloc(rows*cols, sizeof(unsigned char));
   /****************************************************************************
   * Read in the data from the breast segmentation file.
   ****************************************************************************/
   if(segmentation_filename != NULL){
      read_segmentation_file(segmentation_filename, &numpoints, &xcoord, &ycoord,
         &startx, &starty, &endx, &endy, &segmentation_rows, &segmentation_cols);

      /*************************************************************************
      * Convert the segmentation coordinates to our current working resolution.
      *************************************************************************/
      convert_segmentation_coordinates(segmentation_rows, segmentation_cols,
         numpoints, xcoord, ycoord, &startx, &starty, &endx, &endy, rows, cols);

      draw_segmentation(drawim, rows, cols, numpoints, xcoord, ycoord, startx, starty, endx, endy, YELLOWINDEX);
   }

   /****************************************************************************
   * Read in the detection file.
   ****************************************************************************/
   if(detection_filename != NULL){
      if(read_detection_file(detection_filename, &num_detections, &detections,
         &detection_rows, &detection_cols, detmax, detthresh) == 0)
         exit(1);

      /*************************************************************************
      * Rescale the dimensions to be in the same scale as the ground truth.
      *************************************************************************/
      convert_detection_coordinates(detection_rows, detection_cols, num_detections,
         detections, rows, cols);

      if(VERBOSE) printf("There are %d detections to plot.\n", num_detections);
      draw_detections(drawim, rows, cols, num_detections, detections, REDINDEX);
   }

   /****************************************************************************
   * Read in the overlay file.
   ****************************************************************************/
   if(ics_filename != NULL){
      if(read_ics_file(ics_filename, &icsdata, (char *)NULL) == 0) exit(1);

      overlay_filename = NULL;
      if((strcmp(view, "LEFT_CC") == 0) && (icsdata.left_cc.overlay_exists == 1)){
         overlay_filename = icsdata.left_cc.overlay_filename;
         overlay_rows = icsdata.left_cc.rows;
         overlay_cols = icsdata.left_cc.cols;
      }
      else if((strcmp(view, "LEFT_MLO") == 0) && (icsdata.left_mlo.overlay_exists == 1)){
         overlay_filename = icsdata.left_mlo.overlay_filename;
         overlay_rows = icsdata.left_mlo.rows;
         overlay_cols = icsdata.left_mlo.cols;
      }
      else if((strcmp(view, "RIGHT_CC") == 0) && (icsdata.right_cc.overlay_exists == 1)){
         overlay_filename = icsdata.right_cc.overlay_filename;
         overlay_rows = icsdata.right_cc.rows;
         overlay_cols = icsdata.right_cc.cols;
      }
      else if((strcmp(view, "RIGHT_MLO") == 0) && (icsdata.right_mlo.overlay_exists == 1)){
         overlay_filename = icsdata.right_mlo.overlay_filename;
         overlay_rows = icsdata.right_mlo.rows;
         overlay_cols = icsdata.right_mlo.cols;
      }

      if(overlay_filename != NULL){
         status = read_overlay_file(overlay_filename, &overlay_data);
         if(status == 0) exit(1);

         if(VERBOSE){
            printf("   The overlay file is %s.\n", overlay_filename);

            switch(overlay_data.total_abnormalities){
               case 0:  printf("   There are %d abnormalities in the overlay file.\n",
                           overlay_data.total_abnormalities);
                        break;
               case 1:  printf("   There is %d abnormality in the overlay file.\n",
                           overlay_data.total_abnormalities);
                        break;
               default: printf("   There are %d abnormalities in the overlay file.\n",
                           overlay_data.total_abnormalities);
            }
         }

         /**********************************************************************
         * Determine which abnormalities we will use. We are not using any
         * abnormalities that are not "MALIGNANT" of "BENIGN". Also we will only
         * use regions that match the lesiontype and the ms, mm, ct and cd.
         **********************************************************************/
         determine_valid_regions(overlay_data, lesiontype, pathology,
            ms, mm, ct, cd, &isvalid);

         draw_valid_abnormalities(drawim, rows, cols, &overlay_data, overlay_rows,
            overlay_cols, isvalid, GREENINDEX);

      }
   }

   /****************************************************************************
   * Allocate images for red, green and blue.
   ****************************************************************************/
   redim = (unsigned char *)calloc(rows*cols, sizeof(unsigned char));
   grnim = (unsigned char *)calloc(rows*cols, sizeof(unsigned char));
   bluim = (unsigned char *)calloc(rows*cols, sizeof(unsigned char));
   if((redim == NULL) || (grnim == NULL) || (bluim == NULL)){
      fprintf(stderr, "Calloc error!\n");
      exit(1);
   }

   /****************************************************************************
   * Convert the drawing to color and overlay it on the image.
   ****************************************************************************/
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         if(drawim[pos] == 0) redim[pos] = grnim[pos] = bluim[pos] = image[pos];
         else{
            redim[pos] = color[drawim[pos]].r;
            grnim[pos] = color[drawim[pos]].g;
            bluim[pos] = color[drawim[pos]].b;
         }
      }
   }

   /****************************************************************************
   * Write the PPM image to the file.
   ****************************************************************************/
   sprintf(draw_filename, "%s_draw.ppm", image_filename);
   if(write_ppm_image(draw_filename, redim, grnim, bluim, rows, cols, (char *)NULL, 255) == 0) exit(1);

   free(redim);
   free(grnim);
   free(bluim);

   return(1);
}

/*******************************************************************************
* Function: draw_valid_abnormalities
* Purpose: This function will draw the boundaries of the regions in the
* ground truth file that matched the user defined limits of which regions
* to use.
* Name: Michael Heath, University of South Florida
* Date: 1/28/2000
*******************************************************************************/
void draw_valid_abnormalities(unsigned char *image, int rows, int cols,
    OVERLAY_DATA *overlay_data, int overlay_rows, int overlay_cols,
    int *isvalid, int grayval)
{
   int i, v, x0, y0, x1, y1;
   OUTLINE *boundary;
   float scale_r, scale_c, scale;

   /****************************************************************************
   * Compute the scale factor for rescaling the coordinates from the original
   * image size to the size of the image we want to draw them in.
   ****************************************************************************/
   scale_r = (float)rows / (float) overlay_rows;
   scale_c = (float)cols / (float) overlay_cols;
   if(scale_r < scale_c) scale = scale_r;
   else scale = scale_c;

   /****************************************************************************
   * Draw each valid boundary on the image, rescling the coordinates as we go.
   ****************************************************************************/

   for(i=0;i<overlay_data->total_abnormalities;i++){

      if(isvalid[i] == 1){
         boundary = &(overlay_data->abnormalities[i].boundary);

         for(v=0;v<(boundary->length-1);v++){

            x0 = (int)floor(boundary->chain[v].c * scale);
            y0 = (int)floor(boundary->chain[v].r * scale);
            x1 = (int)floor(boundary->chain[v+1].c * scale);
            y1 = (int)floor(boundary->chain[v+1].r * scale);

            /* printf("Drawing (%4d,%4d)->(%4d,%4d)\n", x0, y0, x1, y1); */
            drawimageline(image, rows, cols, x0, y0, x1, y1, grayval);

         }

         x0 = x1;
         y0 = y1;
         x1 = (int)floor(boundary->chain[0].c * scale);
         y1 = (int)floor(boundary->chain[0].r * scale);

         /* printf("Drawing (%4d,%4d)->(%4d,%4d)\n", x0, y0, x1, y1); */
         drawimageline(image, rows, cols, x0, y0, x1, y1, grayval);
      }
   }
}

/*******************************************************************************
* Function: draw_detections
* Purpose: To draw a detection star at each detection site.
* Name: Michael Heath, University of South Florida
* Date: 1/26/2000
*******************************************************************************/
void draw_detections(unsigned char *image, int rows, int cols, int num,
   REGION *detections, int grayval)
{
   int i, x0, y0, hvdist = 6, diagdist = 4;

   /****************************************************************************
   * Draw the breast boundary on the image.
   ****************************************************************************/
   for(i=0;i<num;i++){
      x0 = (int)floor(detections[i].centroid_c);
      y0 = (int)floor(detections[i].centroid_r);

      /* printf("Drawing (%4d,%4d)\n", x0, y0); */
      drawimageline(image, rows, cols, x0-hvdist, y0, x0+hvdist, y0, grayval);
      drawimageline(image, rows, cols, x0, y0-hvdist, x0, y0+hvdist, grayval);
      drawimageline(image, rows, cols, x0-diagdist, y0-diagdist, x0+diagdist, y0+diagdist, grayval);
      drawimageline(image, rows, cols, x0-diagdist, y0+diagdist, x0+diagdist, y0-diagdist, grayval);
   }
}

/*******************************************************************************
* Function: draw_segmentation
* Purpose: To draw the breast boundary on the image.
* Name: Michael Heath, University of South Florida
* Date: 1/26/2000
*******************************************************************************/
void draw_segmentation(unsigned char *image, int rows, int cols, int numpoints,
    float *xcoord, float *ycoord, float startx, float starty,
    float endx, float endy, int grayval)
{
   int i, x0, x1, y0, y1;

   /****************************************************************************
   * Draw the breast boundary on the image.
   ****************************************************************************/
   for(i=0;i<numpoints;i++){
      x0 = (int)floor(xcoord[i]);
      y0 = (int)floor(ycoord[i]);
      x1 = (int)floor(xcoord[(i+1)%numpoints]);
      y1 = (int)floor(ycoord[(i+1)%numpoints]);

      /* printf("Drawing (%4d,%4d)->(%4d,%4d)\n", x0, y0, x1, y1); */
      drawimageline(image, rows, cols, x0, y0, x1, y1, grayval);
   }
}

/*******************************************************************************
* Function: get_commandline_parameters
* Purpose: To get the command line parameters and to check to make sure they
* were all specified. A default value of "BOTH" is used for lesiontype if the
* user did not specify this option.
* Name: Michael Heath, University of South Florida
* Date: 11/13/98
*******************************************************************************/
static int get_commandline_parameters(int argc, char *argv[], char **ics_filename,
   char **view, char *lesiontype, char **pathology, char **ms, char **mm,
   char **ct, char **cd, char **detection_filename, int *detmax, float *detthresh,
   int *use_all_regions, char **segmentation_filename, char **image_filename)
{
   int i;

   /****************************************************************************
   * Print help to the user if the program was run with no command line
   * parameters or if the user had a parameter beginning with "-h".
   ****************************************************************************/
   if(argc == 1){
      print_help();
      return(0);
   }
   for(i=1;i<argc;i++){
      if(strncmp(argv[i], "-h", 2) == 0){
         print_help();
         return(0);
      }
   }

   /****************************************************************************
   * Assign default parameters to each variable.
   ****************************************************************************/
   *ics_filename = (char *)NULL;
   *view = NULL;
   *use_all_regions = 0;
   strcpy(lesiontype, "NONE");
   *ms = NULL;
   *mm = NULL;
   *ct = NULL;
   *cd = NULL;
   *segmentation_filename = (char *)NULL;
   *image_filename = NULL;
   *detection_filename = NULL;
   *detmax = 100;
   *detthresh = 0.0;
   *pathology = NULL;

   /****************************************************************************
   * Extract the command line parameters.
   ****************************************************************************/
   for(i=1;i<argc;i++){
      if(strcmp(argv[i], "-ics") == 0){ *ics_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-view") == 0){ *view = argv[i+1]; i++; }
      if(strcmp(argv[i], "-segment") == 0){ *segmentation_filename = argv[i+1]; i++; }
      if(strcmp(argv[i], "-image") == 0){ *image_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-detection") == 0){
         if((i+3) < argc){
            *detection_filename = argv[i+1]; i++;
            *detmax = atoi(argv[i+1]); i++;
            *detthresh = atof(argv[i+1]); i++;
         }
      }
      else if(strcmp(argv[i], "-v") == 0) VERBOSE = 1;
      else if(strcmp(argv[i], "-version") == 0){
         printf("\n\n%s Version: %s %s\n\n", argv[0], VERSION, VERSIONDATE);
         exit(1);
      }
      else if(strcmp(argv[i], "-all") == 0) *use_all_regions = 1;
      else if(strcmp(argv[i], "-pathology") == 0){
         if((i+1) < argc){
            if(strcmp(argv[i+1], "MALIGNANT") == 0)
               *pathology = argv[i+1];
            else if(strcmp(argv[i+1], "BENIGN") == 0)
               *pathology = argv[i+1];
            else{
               fprintf(stderr, "mktemplete: Invalid pathology (MALIGNANT or BENIGN).\n");
               return(0);
            }
         }
         else{
            fprintf(stderr, "mktemplate: Error! No pathology was specified.\n");
            return(0);
         }
         i++;
      }
      else if(strcmp(argv[i], "-l") == 0){
         if((i+1) < argc){
            if(strcmp(argv[i+1], "CALCIFICATION") == 0)
               strcpy(lesiontype, argv[i+1]);
            else if(strcmp(argv[i+1], "MASS") == 0)
               strcpy(lesiontype, argv[i+1]);
            else if(strcmp(argv[i+1], "BOTH") == 0)
               strcpy(lesiontype, argv[i+1]);
            else{
               fprintf(stderr, "mktemplete: Invalid lesiontype (CALCIFICATION or MASS).\n");
               return(0);
            }
         }
         else{
            fprintf(stderr, "mktemplate: Error! No lesiontype was specified.\n");
            return(0);
         }
         i++;
      }
      else if(strcmp(argv[i], "-ms") == 0){
         if((i+1) < argc) *ms = argv[i+1];
         if((strcmp(lesiontype, "BOTH") == 0) || (strcmp(lesiontype, "CALCIFICATION") == 0)){
            fprintf(stderr, "\nError! You are trying to reference conflicting mass types.\n\n");
            exit(1);
         }
         strcpy(lesiontype, "MASS");
         i++;
      }
      else if(strcmp(argv[i], "-mm") == 0){
         if((i+1) < argc) *mm = argv[i+1];
         if((strcmp(lesiontype, "BOTH") == 0) || (strcmp(lesiontype, "CALCIFICATION") == 0)){
            fprintf(stderr, "\nError! You are trying to reference conflicting mass types.\n\n");
            exit(1);
         }
         strcpy(lesiontype, "MASS");
         i++;
      }
      else if(strcmp(argv[i], "-ct") == 0){
         if((i+1) < argc) *ct = argv[i+1];
         if((strcmp(lesiontype, "BOTH") == 0) || (strcmp(lesiontype, "MASS") == 0)){
            fprintf(stderr, "\nError! You are trying to reference conflicting mass types.\n\n");
            exit(1);
         }
         strcpy(lesiontype, "CALCIFICATION");
         i++;
      }
      else if(strcmp(argv[i], "-ms") == 0){
         if((i+1) < argc) *cd = argv[i+1];
         if((strcmp(lesiontype, "BOTH") == 0) || (strcmp(lesiontype, "MASS") == 0)){
            fprintf(stderr, "\nError! You are trying to reference conflicting mass types.\n\n");
            exit(1);
         }
         strcpy(lesiontype, "CALCIFICATION");
         i++;
      }
   }

   if(strcmp(lesiontype, "NONE") == 0) strcpy(lesiontype, "BOTH");

   /****************************************************************************
   * If the user specified the -all flag indicating they want to render all
   * regions, and they also specified a flag to limit the lesion type, warn
   * them they are doing something they probably don't mean to do.
   ****************************************************************************/
   if(*use_all_regions != 0){
      if(strcmp(lesiontype, "BOTH") != 0){
         fprintf(stderr, "Warning! You set -all to use all regions but you also are limiting the\n");
         fprintf(stderr, "lesion type. The -all flag overrides the limitations.\n");
      }
   }

   /****************************************************************************
   * Make sure the user entered all three required parameters.
   ****************************************************************************/
   if(*image_filename == NULL){
      fprintf(stderr, "drawimage: Error! You did not specify all needed parameters.\n");

      return(0);
   }

   /****************************************************************************
   * Make sure that the view is valid.
   ****************************************************************************/
   if(*ics_filename != NULL){
      if((*view == NULL) || (!((strcmp(*view, "LEFT_CC") == 0) || (strcmp(*view, "LEFT_MLO") == 0) ||
         (strcmp(*view, "RIGHT_CC") == 0) || (strcmp(*view, "RIGHT_MLO")) == 0))){
         fprintf(stderr, "drawimage: Error! Invalid view specified (%s).\n", *view);
         return(0);
      }
   }

   return(1);
}

/*******************************************************************************
* Function: print_help
* Purpose: To print out the command line help.
* Name: Michael Heath, University of South Florida
* Date: 1/26/2000
*******************************************************************************/
static void print_help()
{
   printf("\n\n********************************************************************************\n");
   printf("This program can be used to display the breast boundary, all or select regions\n");
   printf("marked in the ground truth (OVERLAY) file in a DDSM case and selected detected\n");
   printf("regions on an image. The input image must be in 8-bit PGM format. The output\n");
   printf("image is in PPM format and is 24bit color.\n");
   printf("********************************************************************************\n\n");

   printf("<USAGE> drawimage [-detection filename.det max# threshold#]\n");
   printf("   [-segment filename.sgt] [-v] [-version]\n");
   printf("   [-ics file.ics -view VIEW [-all][-l lesiontype ][-ms DESC][-mm DESC]\n");
   printf("   [-ct DESC][-cd DESC]] -image filename.pgm\n\n");

   printf("   -image       The input PGM image to draw on. The output is stored in a new\n");
   printf("                PPM file.\n");
   printf("   -detection   Overlay detections from the specified file. Only use at most\n");
   printf("                max# detections from the file. Use only detections with\n");
   printf("                suspiciousness values greater than or equal to threshold#.\n");
   printf("   -segment     Overlay the breast boundary in the specified file on the\n");
   printf("                image.\n");
   printf("   -ics         An ics file from the DDSM database. To overlay ground truth\n");
   printf("                regions from a DDSM case on the image you must specify both\n");
   printf("                -ics and -view. You can also specify other parameters to\n");
   printf("                select which ground truth\n");
   printf("                regions you want to have displayed.\n");
   printf("   -view        The view to process RIGHT_MLO, RIGHT_CC, LEFT_MLO or LEFT_CC.\n");
   printf("   -all         Use all of the regions (overrides -l, -ms, -mm, -ct and -cd).\n");
   printf("   -l           Skip regions that are not lesions you specify (CALCIFICATION or\n");
   printf("                MASS)\n");
   printf("   -ms          DESC-DESC-...-DESC    Skip masses without these shapes.\n");
   printf("                  \"ROUND\",\"OVAL\",\"LOBULATED\",\"IRREGULAR\"\n");
   printf("                  \"ARCHITECTURAL_DISTORTION\",\"TUBULAR\",\"LYMPH_NODE\"\n");
   printf("                  \"ASYMMETRIC_BREAST_TISSUE\",\"FOCAL_ASYMMETRIC_DENSITY\"\n");
   printf("   -mm DESC-DESC-...-DESC    Skip masses without these margins.\n");
   printf("                  \"CIRCUMSCRIBED\",\"MICROLOBULATED\",\"OBSCURED\",\"ILL_DEFINED\"\n");
   printf("                  \"SPICULATED\"\n");
   printf("   -ct DESC-DESC-...-DESC    Skip calcifications without this type.\n");
   printf("                  \"PUNCTATE\",\"AMORPHOUS\",\"PLEOMORPHIC\",\"ROUND_AND_REGULAR\"\n");
   printf("                  \"LUCENT_CENTER\",\"FINE_LINEAR_BRANCHING\",\"SKIN\",\"VASCULAR\"\n");
   printf("                  \"COARSE\",\"LARGE_RODLIKE\",\"EGGSHELL\",\"MILK_OF_CALCIUM\"\n");
   printf("                  \"SUTURE\",\"DYSTROPHIC\"\n");
   printf("   -cd DESC-DESC-...-DESC    Skip calcifications without this distribution.\n");
   printf("                  \"CLUSTERED\",\"LINEAR\",\"SEGMENTAL\",\"REGIONAL\"\n");
   printf("                  \"DIFFUSELY_SCATTERED\"\n");
   printf("   -version    Print the version of the software and exit.\n");
   printf("   -v          Run the program in verbose mode.\n\n");
}

/*******************************************************************************
* Function: drawimageline
* Purpose: This function draws a line on an image. This algorithm came from
* the book "Computer Graphics Principles and Practice", 2nd edition by
* Foley, vanDam, Feiner and Hughes. It looks like I butchered the code a
* bit, but I don't remember why. I might have changed the code so that
* we increment the variable (x or y) with the smaller magnitude of the slope
* so we won't draw a dotten line when the magnitude of the slope is high.
* Other than that, this is the midpoint line scan-conversion algorithm on page
* 78 of the previously mentioned book.
* Name: Mike Heath
* Date: 10/20/99
*******************************************************************************/
void drawimageline(unsigned char *image, int rows, int cols, int x0, int y0, int x1, int y1, int valtoplot)
{
   int dx, dy, incrE, incrNE, incrN, incrSE, incrNW, d, x, y, tmp;

   /****************************************************************************
   * If abs(dx) >= abs(dy) draw by incrementing x.
   ****************************************************************************/
   if(abs((x1-x0)) >= abs((y1-y0))){
      dx = x1 - x0;
      if(dx < 0){    /* Make sure we draw in positive x */
         tmp = x0;
         x0 = x1;
         x1 = tmp;
         tmp = y0;
         y0 = y1;
         y1 = tmp;
      }      
      dx = x1 - x0;
      dy = y1 - y0;

      /* Draw with E, NE */
      if(dy >= 0){
         d = 2*dy - dx;
         incrE = 2 * dy;
         incrNE = 2 * (dy - dx);
         x = x0;
         y = y0;
         image[y*cols+x] = (unsigned char)valtoplot;
         /* drawpoint(scanvas.image, y, x, wrows, wcols, 0, 0); */
         while(x < x1){
            if(d <= 0){
               d = d + incrE;
               x = x + 1;
            }
            else{
               d = d + incrNE;
               x = x + 1;
               y = y + 1;
            }
            image[y*cols+x] = (unsigned char)valtoplot;
            /* drawpoint(scanvas.image, y, x, wrows, wcols, 0, 0); */
         }
      }

      /* Draw with E, SE */
      else{
         d = 2 * dy + dx;
         incrE = 2 * dy;
         incrSE = 2 * (dy + dx);
         x = x0;
         y = y0;
         image[y*cols+x] = (unsigned char)valtoplot;
         /* drawpoint(scanvas.image, y, x, wrows, wcols, 0, 0); */
         while(x < x1){
            if(d <= 0){
               d = d + incrSE;
               x = x + 1;
               y = y - 1;
            }
            else{
               d = d + incrE;
               x = x + 1;
            }
            image[y*cols+x] = (unsigned char)valtoplot;
            /* drawpoint(scanvas.image, y, x, wrows, wcols, 0, 0); */
         }
      }
   }
   /****************************************************************************
   * Draw by incrementing y.
   ****************************************************************************/
   else{
      dy = y1 - y0;
      if(dy < 0){    /* Make sure we draw in positive y */
         tmp = x0;
         x0 = x1;
         x1 = tmp;
         tmp = y0;
         y0 = y1;
         y1 = tmp;
      }      
      dx = x1 - x0;
      dy = y1 - y0;

      /* Draw with N, NE */
      if(dx >= 0){
         d = dy - 2 * dx;
         incrN = -2 * dx;
         incrNE = 2 * (dy - dx);
         x = x0;
         y = y0;
         image[y*cols+x] = (unsigned char)valtoplot;
         /* drawpoint(scanvas.image, y, x, wrows, wcols, 0, 0); */
         while(y < y1){
            if(d <= 0){
               d = d + incrNE;
               y = y + 1;
               x = x + 1;
            }
            else{
               d = d + incrN;
               y = y + 1;
            }
            image[y*cols+x] = (unsigned char)valtoplot;
            /* drawpoint(scanvas.image, y, x, wrows, wcols, 0, 0); */
         }
      }

      /* Draw with N, NW */
      else{
         d = -dy - 2 * dx;
         incrN = -2 * dx;
         incrNW = 2 * (-dy - dx);
         x = x0;
         y = y0;
         image[y*cols+x] = (unsigned char)valtoplot;
         /* drawpoint(scanvas.image, y, x, wrows, wcols, 0, 0); */
         while(y < y1){
            if(d <= 0){
               d = d + incrN;
               y = y + 1;
            }
            else{
               d = d + incrNW;
               y = y + 1;
               x = x - 1;
            }
            image[y*cols+x] = (unsigned char)valtoplot;
            /* drawpoint(scanvas.image, y, x, wrows, wcols, 0, 0); */
         }

      }
   }
}
