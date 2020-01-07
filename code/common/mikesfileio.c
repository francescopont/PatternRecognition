/*******************************************************************************
* File: mikesfileio.c
* Purpose: This file contains code for reading and writing several types of
* files used in the detection of abnormailities in mammograms.
* Name: Michael Heath, University of South Florida
* Date: 7/26/99
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*
* Note: On 4/17/2000 the read_detection file was modified to allow the reading
* of the coordinates of the last detection in the file. Previously, the
* position was incorrect. It was actually the second to last position. This
* change fixed this bug. This probably has a little effect on the overall
* system because the last listed detection is the least suspicious one. Still,
* it was an error, and it was fixed. - Michael Heath
*******************************************************************************/
#include <string.h>
#include "mikesfileio.h"
#include "myalloc.h"

extern int VERBOSE;

/*******************************************************************************
* Function: convert_segmentation_coordinates
* Purpose: The segmentation coordinates can be rescaled to a different size
* by using this function.
* Name: Michael Heath, University of South Florida
* Date: 1/26/2000
*******************************************************************************/
void convert_segmentation_coordinates(int present_rows, int present_cols,
   int numpoints, float *xcoord, float *ycoord, float *startx, float *starty,
   float *endx, float *endy, int desired_rows, int desired_cols)
{
   int i;
   float scale_r, scale_c, scale;

   /****************************************************************************
   * Just return if no scaling needs to be done.
   ****************************************************************************/
   if((present_rows == desired_rows) && (present_cols == desired_cols)) return;

   /****************************************************************************
   * Compute the scale factor for rescaling the coordinates. Use the smaller
   * scale factor.
   ****************************************************************************/
   scale_r = (float)desired_rows / (float)present_rows;
   scale_c = (float)desired_cols / (float)present_cols;
   if(scale_r <= scale_c) scale = scale_r;
   else scale = scale_c;

   /****************************************************************************
   * Rescale the coordinates.
   ****************************************************************************/
   for(i=0;i<numpoints;i++){
      xcoord[i] *= scale;
      ycoord[i] *= scale;
   }
   *startx *= scale;
   *starty *= scale;
   *endx *= scale;
   *endy *= scale;
}

/*******************************************************************************
* Function: read_segmentation_file
* Purpose: This procedure reads the contents of a file that contains an outline
* of the breast. The variables startx, starty, endx and endy will be used by
* future versions of this software.
* Name: Michael Heath, University of South Florida
* Date: 7/26/99
*******************************************************************************/
int read_segmentation_file(char *filename, int *numpoints, float **xcoord, float **ycoord,
   float *startx, float *starty, float *endx, float *endy, int *rows, int *cols)
{
   FILE *fp=NULL;
   char aline[102];

   int i;

   *numpoints = 0;
   *xcoord = NULL;
   *ycoord = NULL;

   /****************************************************************************
   * Open the file for reading.
   ****************************************************************************/
   if((fp = fopen(filename, "r")) == NULL) return(0);

   /****************************************************************************
   * Skip any comments. These are lines that begin with the # character.
   ****************************************************************************/
   while((fgets(aline, 100, fp) != NULL) && !feof(fp)){
      if(strlen(aline) > 100){
         fprintf(stderr, "Error in segmentation file. Lines longer than 100 characters.\n");
         return(0);
      }
      if(aline[0] != '#') break;
   }

   /****************************************************************************
   * Read in the dimensions of the image in which the following coordinates
   * are specified. This allows all coordinates to be specified in any
   * resolution. We will know how we could rescale these coordinates to any
   * size.
   ****************************************************************************/
   sscanf(aline, "%d %d", cols, rows);

   /****************************************************************************
   * Read in the major axis of the breast. The position (startx,starty) is
   * the begining if the axis (near the nipple) and the location (endx,endy)
   * is the end of the axis which should lie near the edge of the film.
   ****************************************************************************/
   fscanf(fp, "%f%f%f%f", startx, starty, endx, endy);

   /****************************************************************************
   * Read in the chain of pixels that specify the breast border in the image.
   * There is a number specifying how many points there are. This is followed
   * by the x and y locations of each point.
   ****************************************************************************/
   fscanf(fp, "%d", numpoints);
   if((*xcoord = (float *) mycalloc("xcoord", *numpoints, sizeof(float))) == NULL){
      fprintf(stderr, "Calloc error in read_segmentation_file().\n");
      return(0);
   }
   if((*ycoord = (float *) mycalloc("ycoord", *numpoints, sizeof(float))) == NULL){
      fprintf(stderr, "Calloc error in read_segmentation_file().\n");
      return(0);
   }
   for(i=0;i<(*numpoints);i++)
      fscanf(fp, "%f%f", (*xcoord) + i, (*ycoord) + i);

   fclose(fp);
   return(1);
}

/*******************************************************************************
* Function: write_segmentation_file
* Purpose: This procedure reads the contents of a file that contains data about
* the image size and the outline of the breast.
* Name: Michael Heath, University of South Florida
* Date: 7/26/99
*******************************************************************************/
int write_segmentation_file(char *filename, int numpoints, float *xcoord, float *ycoord,
   float startx, float starty, float endx, float endy, int rows, int cols)
{
   FILE *fp=NULL;
   int i;

   /****************************************************************************
   * Open the file for writing.
   ****************************************************************************/
   if((fp = fopen(filename, "w")) == NULL){
      fprintf(stderr, "Error writing the breast coordinates file named %s.\n", filename);
      return(0);
   }

   /****************************************************************************
   * Write the size of the image in which the following coordinates are
   * specified. This would allow all coordinates to be at say 1/2 the size of
   * the DDSM image if the following rows and cols written here are DDSMcols/2
   * and DDSMrows/2. This way the coordinates can be scaled by a program that
   * might use them later to any size.
   ****************************************************************************/
   fprintf(fp, "%d %d\n", cols, rows);

   /****************************************************************************
   * Write out the major axis of the breast. The position (startx,starty) is
   * the begining if the axis (near the nipple) and the location (endx,endy)
   * is the end of the axis which should lie near the edge of the film.
   ****************************************************************************/
   fprintf(fp, "%f %f %f %f\n", startx, starty, endx, endy);

   /****************************************************************************
   * Write the chain of pixels that specify the breast border in the image.
   * There is a number specifying how many points there are. This is followed
   * by the x and y locations of each point.
   ****************************************************************************/
   fprintf(fp, "%d\n", numpoints);

   for(i=0;i<numpoints;i++){
      /* printf("In writing: [%3d] %f %f\n", i, xcoord[i], ycoord[i]); */
      fprintf(fp, "%f %f\n", xcoord[i], ycoord[i]);
   }

   fclose(fp);
   return(1);
}

/*******************************************************************************
* Function: read_detection_file
* Purpose: To read in the detected regions from a file. The file can contain
* detections in two different representations and they can be mixed in the
* same file. A file can have any number of detected regions including 0 regions.
* Detected regions can be of type POLYGON or type PROMPT. A PROMPT will have
* one point (pair of coordinates (x,y)) on the following line, and a POLYGON
* will have all vertices listed on sequential lines. Each point will be listed
* in integer form with the x-coordinate followed by the y coordinate. Only
* spaces can be between the coordinates. Comments can be listed on any line
* by themself and must have a # in the first column. No row of the file can have
* more than 100 characters in it. The first line after comment lines at the
* top of the file must cantain the number of columns and rows of the image
* in which the rest of the coordinates in the file are expressed.
*
* Please note that on 1/13/2000 I added the ability for a floating point
* number to be on the line with the keyword POLYGON or PROMPT. Placing a
* value on that line is optional. If it is there, it represents the
* suspiciousness of the detection. All of the detections should be listed
* in decreasing order of suspicion in the file. The presence of these
* numbers gives the user the option of running the program in a mode where
* all detections greater than or equal to a user input suspiciousness threshold
* to be used, rather than using all detections up to a user selected maximum
* number of detections. All suspiciousness values should be greater than or
* equal to zero.
*
* Name: Michael Heath, University of South Florida
* Date: 11/13/98
*
* # Example detection file.
* # Comment
* 100 175
* POLYGON 90.6
* 50 75
* 60 100
* 90 85
* 75 20
* 55 60
* PROMPT 20.3
* 100 150
*******************************************************************************/
int read_detection_file(char *detection_filename, int *num_detections, REGION **detections,
       int *detection_rows, int *detection_cols, int max, float suspiciousness_threshold)
{
   FILE *fp=NULL;
   int d=0, p, k;
   int num_polygons = 0, num_centroids = 0;
   int numpoints = 0, max_points = 0;
   char aline[110];
   const int polygon=1, centroid=2, unknown=3, comment=4;
   int last_type = 0, data;
   int *xcoord=NULL, *ycoord=NULL;
   float suspiciousness;

   *detection_rows = -1;
   *detection_cols = -1;

   /****************************************************************************
   * Make a first pass through the file to count the number of detections, and
   * the maximum number of points in any detection. Also make sure that no
   * line has more than 100 characters in it and that each PROMPT has only
   * one point and POLYGON detections have at least 3 points.
   ****************************************************************************/
   if((fp = fopen(detection_filename, "rb")) == NULL){
      fprintf(stderr, "Error opening the file %s.\n", detection_filename);
      return(0);
   }

   while((fgets(aline, 100, fp) != NULL) && !feof(fp)){
      if(strlen(aline) > 100){
         fprintf(stderr, "Error in the detection file. Lines longer than 100 characters.\n");
         return(0);
      }

      if(aline[0] == '#') continue;   /* Skip a comment line. */

      if((*detection_rows == -1) && (*detection_cols == -1)){
         sscanf(aline, "%d %d", detection_cols, detection_rows);
         continue;
      }

      if(strncmp("POLYGON", aline, strlen("POLYGON")) == 0){
         if((last_type == polygon) && (numpoints < 3)){
            fprintf(stderr, "There must be at least 3 points for a POLYGON.\n");
            return(0);
         }
         if((last_type == centroid) && (numpoints != 1)){
            fprintf(stderr, "There must be exactly one point for a PROMPT.\n");
            return(0);
         }
         if(numpoints > max_points) max_points = numpoints;
         numpoints = 0;
         num_polygons++;
         last_type = polygon;
      }

      if(strncmp("PROMPT", aline, strlen("PROMPT")) == 0){
         if((last_type == polygon) && (numpoints < 3)){
            fprintf(stderr, "There must be at least 3 points for a POLYGON.\n");
            return(0);
         }
         if((last_type == centroid) && (numpoints != 1)){
            fprintf(stderr, "There must be exactly one point for a PROMPT.\n");
            return(0);
         }
         if(numpoints > max_points) max_points = numpoints;
         numpoints = 0;
         num_centroids++;
         last_type = centroid;
      }

      else numpoints++;
   }
   if(numpoints > max_points) max_points = numpoints;

   fclose(fp);

   /*
   printf("There are %d POLYGON detectections in the file.\n", num_polygons);
   printf("There are %d PROMPT detectections in the file.\n", num_centroids);
   printf("The maximum number of points in any detected region is %d.\n", max_points);
   */

   *num_detections = (num_polygons + num_centroids);

   if((*num_detections) == 0) return(1);

   /****************************************************************************
   * Allocate an arrays for the x and y coordinates. We will read data into
   * these arrays before copying it into the OUTLINE structure or centroid
   * variables.
   ****************************************************************************/
   if(max_points < 16) max_points = 16;
   if((xcoord = (int *) calloc(max_points, sizeof(int))) == NULL){
      fprintf(stderr, "Calloc error!\n");
      return(0);
   }
   if((ycoord = (int *) calloc(max_points, sizeof(int))) == NULL){
      fprintf(stderr, "Calloc error!\n");
      return(0);
   }

   /****************************************************************************
   * Allocate an array of REGION structures.
   ****************************************************************************/
   if(((*detections) = (REGION *) calloc(num_polygons+num_centroids, sizeof(REGION))) == NULL){
      fprintf(stderr, "Calloc error!\n");
      return(0);
   }

   /****************************************************************************
   * Make a second pass through the file to read in all of the detections.
   ****************************************************************************/
   if((fp = fopen(detection_filename, "rb")) == NULL){
      fprintf(stderr, "Error opening the file %s.\n", detection_filename);
      return(0);
   }

   d = 0;         /* The detection region we are curently reading from the file. */
   last_type = 0; /* The previous detected region is neither POLYGON nor CALCIFICATION. */
   numpoints = 0;
   k = 0;

   /****************************************************************************
   * Skip the comment lines at the top of the file and skip the first line
   * after the comment. That first line contains the image dimensions that
   * we read in our first pass through the file.
   ****************************************************************************/
   while((fgets(aline, 110, fp) != NULL) && !feof(fp)){
      if(aline[0] != '#') break;
   }

   while((fgets(aline, 110, fp) != NULL)){
      data = unknown;
      if(aline[0] == '#') data = comment;
      if(strncmp("POLYGON", aline, strlen("POLYGON")) == 0){
         data = polygon;
         suspiciousness = 0.0;
         sscanf(aline+strlen("POLYGON"), "%f", &((*detections)[k].suspiciousness));
         k++;
      }
      if(strncmp("PROMPT", aline, strlen("PROMPT")) == 0){
         data = centroid;
         suspiciousness = 0.0;
         sscanf(aline+strlen("PROMPT"), "%f", &((*detections)[k].suspiciousness));
         k++;
      }

      if(sscanf(aline, "%d %d", &xcoord[numpoints], &ycoord[numpoints]) == 2){
         numpoints++;
      }

      /*************************************************************************
      * If we have a line with POLYGON or PROMPT on it and the last type
      * was either POLYGON or PROMPT, we must fill the coordinates into the
      * detection structure because we just finished reading the points for
      * either a PROMPT or a POLYGON. The value of last_type will tell us
      * whether we just finished reading a POLYGON or a PROMPT.
      *************************************************************************/
      if((data == polygon) || (data == centroid)){

         if(last_type == polygon){

            (*detections)[d].outline.chain = (CHAINPOINT *) calloc(numpoints, sizeof(CHAINPOINT));
            (*detections)[d].outline.length = numpoints;

            for(p=0;p<numpoints;p++){
               (*detections)[d].outline.chain[p].c = xcoord[p];
               (*detections)[d].outline.chain[p].r = ycoord[p];
            }
            d++;
         }

         if(last_type == centroid){
            (*detections)[d].outline.length = 0;
            (*detections)[d].centroid_c = xcoord[0];
            (*detections)[d].centroid_r = ycoord[0];

            d++;
         }

         numpoints = 0;
         last_type = data;  /* Store the type we will start to read in last_type. */
      }
   }

   /****************************************************************************
   * Finish the last detection region. There is no POLYGON or PROMPT tag
   * after the last detection region so the above code will not handle it.
   ****************************************************************************/
   if(last_type == polygon){

      (*detections)[d].outline.chain = (CHAINPOINT *) calloc(numpoints, sizeof(CHAINPOINT));
      (*detections)[d].outline.length = numpoints;

      for(p=0;p<numpoints;p++){
         (*detections)[d].outline.chain[p].c = xcoord[p];
         (*detections)[d].outline.chain[p].r = ycoord[p];
      }
      d++;
   }

   if(last_type == centroid){
      (*detections)[d].outline.length = 0;
      (*detections)[d].centroid_c = xcoord[0];
      (*detections)[d].centroid_r = ycoord[0];

      d++;
   }

   if(xcoord != NULL) free(xcoord);
   if(ycoord != NULL) free(ycoord);

   fclose(fp);

   /****************************************************************************
   * Cut the list of detections to the maximim number allowed.
   ****************************************************************************/
   if((*num_detections > max)){
      *num_detections = max;
      if(VERBOSE) printf("\n   Limiting the number of detections to %d.\n", max);
   }

   /****************************************************************************
   * Cut the list of detections using the suspiciousness threshold.
   ****************************************************************************/
   d = 0;
   while((d<(*num_detections)) && ((*detections)[d].suspiciousness >= suspiciousness_threshold)) d++;
   if(d < (*num_detections)){
      if(VERBOSE) printf("\n   Cutting %d detections using suspiciousness threshold of %f.\n",
	 (*num_detections)-d, suspiciousness_threshold);
      *num_detections = d;
   }

   /*
   for(d=0;d<(*num_detections);d++){
      printf("DET[%2d] = (%f,%f) (%f)\n", d,
            (*detections)[d].centroid_c,
            (*detections)[d].centroid_r,
            (*detections)[d].suspiciousness);
   }
   */

   return(1);
}

/*******************************************************************************
* Function: convert_detection_coordinates
* Purpose: The coordinates of the detections are specified in terms of some
* size image. We can use this function to rescale the points to a different
* size image.
* Name: Michael Heath, University of South Florida
* Date: 1/26/2000
*******************************************************************************/
void convert_detection_coordinates(int present_rows, int present_cols, int num_detections,
      REGION *detections, int desired_rows, int desired_cols)
{
   int numpoints, d, p;
   float scale_r, scale_c, scale;

   /****************************************************************************
   * If the present_rows and present_cols are equal to the desired_rows and
   * desired_cols just return because no rescaling needs to be done.
   ****************************************************************************/
   if((present_rows == desired_rows) && (present_cols == desired_cols)) return;

   /****************************************************************************
   * Compute the scale factor for rescaling the coordinates. Use the smaller
   * scale factor.
   ****************************************************************************/
   scale_r = (float)desired_rows / (float)present_rows;
   scale_c = (float)desired_cols / (float)present_cols;
   if(scale_r <= scale_c) scale = scale_r;
   else scale = scale_c;

   for(d=0;d<num_detections;d++){

      numpoints = detections[d].outline.length;

      detections[d].centroid_c *= scale;
      detections[d].centroid_r *= scale;

      for(p=0;p<numpoints;p++){
         detections[d].outline.chain[p].c *= scale;
         detections[d].outline.chain[p].r *= scale;
      }
   }
}

/*******************************************************************************
* Function: print_example_det_file
* Purpose: This function prints out an example of a detection file in the
* proper format.
* Name: Michael Heath, University of South Florida
* Date: 1/21/2000
*******************************************************************************/
void print_example_det_file()
{
   printf("# File example.det. You can have any number of lines (including zero) that begin\n");
   printf("# with a # character at the top of the file. These are comments. The first line\n");
   printf("# found after skipping all comment lines must contain the dimensions of the\n");
   printf("# image in which the rest of the coordinates in the file can be plotted. This\n");
   printf("# allows the coordinates to be easily rescaled to any new dimensions. The\n");
   printf("# dimensions must be expressed as two numbers, COLUMNS and ROWS, with a space\n");
   printf("# between them. In this example file there are 2611 COLUMNS and 5281 ROWS.\n");
   printf("# The remainder of the file is basically a list of detections ordered in\n");
   printf("# decreasing likelihood of suspicion. You can have any number of detections\n");
   printf("# between 0 and 100. Each detection can be either a PROMPT or a POLYGON. A\n");
   printf("# PROMPT is specified with two lines. The first line has the keyword PROMPT\n");
   printf("# on it and this can be followed by a non-negative suspiciousness value if\n");
   printf("# you want to place one there. The next line must have the coordinates of a\n");
   printf("# point in the image. Note that the coordinate 0,0 is in the upper left hand\n");
   printf("# corner of the image. This point is written as a COLUMN and a ROW, in that\n");
   printf("# order as whole numbers. A POLYGON is specified like a PROMPT, but there must\n");
   printf("# be at least three lines of coordinates to specify the vertices of the\n");
   printf("# polygon. When a polygon is specified, the program computes the centroid of\n");
   printf("# the polygon and treats this centroid like a prompt.\n");
   printf("2611 5281\n");
   printf("PROMPT 100.0\n");
   printf("1393 3283\n");
   printf("POLYGON 97.8\n");
   printf("1057 3200\n");
   printf("1040 3248\n");
   printf("1090 3230\n");
   printf("PROMPT 95.2\n");
   printf("686 3213\n");
   printf("PROMPT 30.3\n");
   printf("980 2436\n");
   printf("PROMPT 12.5\n");
   printf("1169 3262\n");
   printf("PROMPT 0.1\n");
   printf("448 3430\n");
   printf("\n\n");
}
