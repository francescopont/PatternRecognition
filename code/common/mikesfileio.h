/*******************************************************************************
* File: mikesfileio.h
* Purpose: This file contains prototypes for reading and writing several types of
* files used in the detection of abnormailities in mammograms.
* Name: Michael Heath, University of South Florida
* Date: 7/26/99
* Notes: 6/9/2000 - I added the macros to only allow this header file to be
*        included once.
*      
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#ifndef _MIKESFILEIO_
#define _MIKESFILEIO_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "overlay.h"

typedef struct{
   OUTLINE outline;
   double centroid_r, centroid_c;
   unsigned char *image;
   int rows, cols;
   int offset_r, offset_c;
   int valid;
   int matched;
   float suspiciousness;
}REGION;

void convert_segmentation_coordinates(int present_rows, int present_cols,
   int numpoints, float *xcoord, float *ycoord, float *startx, float *starty,
   float *endx, float *endy, int desired_rows, int desired_cols);
int read_segmentation_file(char *filename, int *numpoints, float **xcoord, float **ycoord,
   float *startx, float *starty, float *endx, float *endy, int *rows, int *cols);
int write_segmentation_file(char *filename, int numpoints, float *xcoord, float *ycoord,
   float startx, float starty, float endx, float endy, int rows, int cols);
int read_detection_file(char *detection_filename, int *num_detections, REGION **detections,
       int *detection_rows, int *detection_cols, int max, float suspiciousness_threshold);
void convert_detection_coordinates(int present_rows, int present_cols, int num_detections,
      REGION *detections, int desired_rows, int desired_cols);
void print_example_det_file();

#endif
