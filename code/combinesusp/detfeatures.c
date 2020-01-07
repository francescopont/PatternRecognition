/*******************************************************************************
* Program: detfeatures
* Purpose: This file contains the source code for computing some features
* for a detection site in a mammography image. The features are:
*
* 1) "radial" distance to the nipple position
* 2) FUM value for a range of radii
* 3) mean intensity for a circular region of a given diameter (in mm)
* 4) angle of the detection from the major axis (positive direction is
*    towards the axilla)
*
* Name: Michael Heath, University of South Florida
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
* Date: 5/15/2000
*
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "detfeatures.h"

#define VERSION "1.0.0"
#define VERSIONDATE "May 15, 2000"

extern int VERBOSE;

/*******************************************************************************
* Function: compute_detfeatures
* Purpose: To compute several features for a detection site in an image.
* Name: Michael Heath, University of South Florida.
* Date: 5/15/2000
*******************************************************************************/
void compute_detfeatures(unsigned short int *image, int rows, int cols, float microns,
   unsigned char *breastmask, float startx, float starty, float endx, float endy,
   int num_detections, REGION *detections, int radius, int deltar,
   DETFEATURE **detfeature)
{
   int i;
   float *nd=NULL, *angle=NULL, *pd=NULL;

   void compute_nipple_distance(float **nd, int num_detections, REGION *detections,
       float startx, float starty, float microns);
   void compute_detection_angle(float **angle, int num_detections, REGION *detections,
       float startx, float starty, float endx, float endy);
   void compute_median_in_circle(unsigned short int *image, int rows, int cols,
      unsigned char *breastmask, int row, int col, float diameter_mm, float microns, float *median);
   void compute_FUM(unsigned short int *image, int rows, int cols, unsigned char *mask,
       int radius, int deltar, int row, int col, int *num_fum, float **fum_vector);
   void compute_projected_distance(float **pd, int num_detections, REGION *detections,
       float startx, float starty, float endx, float endy, float microns);

   /****************************************************************************
   * Allocate an array of detection features.
   ****************************************************************************/
   if(((*detfeature) = calloc(num_detections, sizeof(DETFEATURE))) == NULL){
      fprintf(stderr, "Error allocating an array in compute_detfeatures().\n\n");
      exit(1);
   }

   /****************************************************************************
   * Compute the nipple distance and angle for all of the detections.
   ****************************************************************************/
   compute_nipple_distance(&nd, num_detections, detections, startx, starty, microns);
   compute_detection_angle(&angle, num_detections, detections, startx, starty, endx, endy);
   compute_projected_distance(&pd, num_detections, detections, startx, starty,
      endx, endy, microns);

   for(i=0;i<num_detections;i++){

      (*detfeature)[i].centroid_r = detections[i].centroid_r;
      (*detfeature)[i].centroid_c = detections[i].centroid_c;
      (*detfeature)[i].suspiciousness = detections[i].suspiciousness;

      (*detfeature)[i].nipple_distance = nd[i];
      (*detfeature)[i].angle = angle[i];

      (*detfeature)[i].projected_dist = pd[i];

      compute_FUM(image, rows, cols, breastmask, radius, deltar,
         (int)(detections[i].centroid_r+0.5), (int)(detections[i].centroid_c+0.5),
         &((*detfeature)[i].num_fum), &((*detfeature)[i].fum));

      compute_median_in_circle(image, rows, cols,
          breastmask, (int)(detections[i].centroid_r+0.5), (int)(detections[i].centroid_c+0.5),
          5.0, microns, &((*detfeature)[i].median));
   }

   free(nd);
   free(angle);
   free(pd);
}

/*******************************************************************************
* Function: compute_nipple_distance
* Purpose: This function computes the distance from the detection site to the
* nipple position for each detection.
* Name: Michael Heath, University of South Florida
* Date: 4/16/2000
*******************************************************************************/
void compute_nipple_distance(float **nd, int num_detections, REGION *detections,
    float startx, float starty, float microns)
{
   int i;
   double dx, dy;

   /****************************************************************************
   * Allocate memory for the array of distances.
   ****************************************************************************/
   if(((*nd) = (float *) calloc(num_detections, sizeof(float))) == NULL){
      fprintf(stderr, "Calloc error in compute_nipple_distance().\n");
      exit(1);
   }

   for(i=0;i<num_detections;i++){
      dx = detections[i].centroid_c - (double)startx;
      dy = detections[i].centroid_r - (double)starty;
      (*nd)[i] = sqrt(dx*dx+dy*dy) * microns;
   }
}

/*******************************************************************************
* Function: compute_projected_distance
* Purpose: This function computes the distance from the detection site 
* (projected onto the major axis) to the nipple position for each detection.
* Name: Michael Heath, University of South Florida
* Date: 5/20/2000
*******************************************************************************/
void compute_projected_distance(float **pd, int num_detections, REGION *detections,
    float startx, float starty, float endx, float endy, float microns)
{
   int i;
   double dx, dy, axis_length, ma_dx, ma_dy;

   /****************************************************************************
   * Allocate memory for the array of distances.
   ****************************************************************************/
   if(((*pd) = (float *) calloc(num_detections, sizeof(float))) == NULL){
      fprintf(stderr, "Calloc error in compute_nipple_distance().\n");
      exit(1);
   }

   axis_length = sqrt((double)((endx-startx)*(endx-startx)+(endy-starty)*(endy-starty)));

   ma_dx = (endx-startx) / axis_length;
   ma_dy = (endy-starty) / axis_length;

   for(i=0;i<num_detections;i++){
      dx = detections[i].centroid_c - (double)startx;
      dy = detections[i].centroid_r - (double)starty;
      (*pd)[i] = (dx*ma_dx+dy*ma_dy) * microns;
   }
}

/*******************************************************************************
* Function: compute_detection_angle
* Purpose: This function computes the angle between the detection and the major
* axis. Recall that the major axis is a line segment with the startx and starty
* being the estimated nipple position. Lower row numbers in the image are
* assumed to be oriented towards the axilla. This is the orientation of the
* images in the DDSM database at USF. The computed angle will be positive for
* detections towards the axilla from the major axis and will be negative for
* detections on the "other side" of the major axis. Angles are reported in
* degrees (not radians).
* Name: Michael Heath, University of South Florida
* Date: 4/16/2000
*******************************************************************************/
void compute_detection_angle(float **angle, int num_detections, REGION *detections,
    float startx, float starty, float endx, float endy)
{
   int i;
   double axis_dx, axis_dy, axis_length;
   double detection_dx, detection_dy, detection_length;
   double cos_angle;
   int cross_sign = 1, axis_sign = 1;

   /****************************************************************************
   * Allocate memory for the array of angles.
   ****************************************************************************/
   if(((*angle) = (float *) calloc(num_detections, sizeof(float))) == NULL){
      fprintf(stderr, "Calloc error in compute_nipple_distance().\n");
      exit(1);
   }

   axis_dx = (double)endx - (double)startx;
   axis_dy = (double)starty - (double)endy;
   axis_length = sqrt(axis_dx * axis_dx + axis_dy * axis_dy);

   for(i=0;i<num_detections;i++){
      detection_dx = detections[i].centroid_c - (double)startx;
      detection_dy = (double)starty - detections[i].centroid_r;
      detection_length = sqrt(detection_dx * detection_dx + detection_dy * detection_dy);

      cos_angle = (axis_dx * detection_dx + axis_dy * detection_dy) / (axis_length * detection_length);

      (*angle)[i] = acos(cos_angle) * 180.0 / 3.1415926535;

/*
      printf("(%f,%f) x (%f,%f) = %f\n", detection_dx/detection_length, detection_dy/detection_length,
         axis_dx/axis_length, axis_dy/axis_length,
         (detection_dx/detection_length) * (axis_dy/axis_length) - (detection_dy/detection_length) * (axis_dx/axis_length));
*/

      if(((detection_dx/detection_length)*(axis_dy/axis_length)-(detection_dy/detection_length)*(axis_dx/axis_length)) >= 0)
         cross_sign = +1;
      else cross_sign = -1;

      if(endx >= startx) axis_sign = -1;
      else axis_sign = 1;

      (*angle)[i] *= cross_sign * axis_sign;
   }
}

/*******************************************************************************
* Function: compute_median_in_circle
* Purpose: The find the median pixel value in a disk of a specified diameter.
* Name: Michael Heath, University of South Florida
* Date: 5/16/2000
*******************************************************************************/
void compute_median_in_circle(unsigned short int *image, int rows, int cols,
   unsigned char *breastmask, int row, int col, float diameter_mm, float microns, float *median)
{
   int r, c, d, n=0, i, j, dist_sq;
   int pixthreshdist;
   unsigned short int tempval;
   unsigned short int *vals;

   /****************************************************************************
   * Compute the radius in pixels.
   ****************************************************************************/
   pixthreshdist = (((double)diameter_mm/2.0) * 1000.0) / microns;

   /****************************************************************************
   * Fill an array with values of pixels in the disk that are inside the breast
   * tissue region of the image.
   ****************************************************************************/
   vals = (unsigned short int *) calloc(pixthreshdist*pixthreshdist, sizeof(unsigned short int));

   for(r=(row-pixthreshdist);r<=(row+pixthreshdist);r++){
      for(c=(col-pixthreshdist);c<=(col+pixthreshdist);c++){

         if((r>=0)&&(r<rows)&&(c>=0)&&(c<cols)&&(breastmask[r*cols+c]>0)){
            dist_sq = (r-row)*(r-row)+(c-col)*(c-col);
            d = (int)floor((double)dist_sq);

            if(d <= pixthreshdist){
               vals[n] = image[r*cols+c];
               n++;
            }
         }
      }
   }

   /****************************************************************************
   * Sort the first half of the array to find the median. This sort is slow,
   * but since this function will be called so few times and n is small, this
   * does not matter much.
   ****************************************************************************/
   for(i=0;i<(n/2);i++){
      for(j=0;j<n;j++){
         if(vals[i] > vals[j]){
            tempval = vals[i];
            vals[i] = vals[j];
            vals[j] = tempval;
         }
      }
   }

   *median = (float)vals[n/2];

   free(vals);
}

/*******************************************************************************
* Function: compute_FUM
* Purpose: To compute the AFUM for a pixel position. This computes the average
* fraction of the pixels in circles with a range of radii that are lower valued
* than the minimum valued pixel in a disk with a smaller radius. The AFUM value
* is calculated at the specified row and column. In addition to the AFUM value,
* all the the FUM values that were calculated are also returned.
* Name: Michael Heath, University of South Florida
* Date: 11/29/99
*******************************************************************************/
void compute_FUM(unsigned short int *image, int rows, int cols, unsigned char *mask,
    int radius, int deltar, int row, int col, int *num_fum, float **fum_vector)
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
   r = row;
   c = col;
   pos = r*cols + c;

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

      /*************************************************************************
      * Allocate memory to store the fum values.
      *************************************************************************/
      *num_fum = radius - deltar + 1;
      if(((*fum_vector) = (float *) calloc((*num_fum), sizeof(float))) == NULL){
         fprintf(stderr, "Error callocing an array in compute_FUM().\n\n");
         exit(1);

      }

      afum = 0;
      for(n=deltar;n<=radius;n++){
         if(num_below_min[n] > num_at_radius[n]) fprintf(stderr, "Error! %d,%d\n", num_below_min[n], num_at_radius[n]);
         (*fum_vector)[n-deltar] = (double)num_below_min[n] / (double)num_at_radius[n];
         afum += (double)num_below_min[n] / (double)num_at_radius[n];
         
      }
      /*
      if((afum/(double)((radius-deltar) + 1)) < 1.0)
         afumimage[pos] = (unsigned char)floor(255.0 * afum / (double)((radius-deltar) + 1));
      */
   }
   else{
      *num_fum = 0;
      *fum_vector = NULL;
   }
}
