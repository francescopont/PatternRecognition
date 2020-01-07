/*******************************************************************************
* File: detfeatures.h
*
* 1) "radial" distance to the nipple position
* 2) FUM value for a range of radii
* 3) mean intensity for a circular region of a given diameter (in mm)
* 4) angle of the detection from the major axis (positive direction is
*    towards the axilla)
*
* Name: Michael Heath, University of South Florida
* Date: 5/15/2000
*
*******************************************************************************/
#ifndef _DETFEATURES_
#define _DETFEATURES_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../common/mikesfileio.h"

/*******************************************************************************
* Create a data stucture to hold all of these features.
*******************************************************************************/
typedef struct{
   double centroid_r, centroid_c;   /* The coords of the detection site.      */
   float suspiciousness;            /* Suspicion value of the detection site. */

   float nipple_distance;           /* Detection to nupple distance in mm.    */
   float angle;                     /* Angle of det. with major axis.         */
   int num_fum;                     /* Number of sizes of the FUM.            */
   float *fum;                      /* FUM values at different radii of disk. */
   float median;                    /* Median intensity in disk of given diam */
   float projected_dist;            /* Detection projected onto the major axis. */

   float aug_suspiciousness;

}DETFEATURE;

void compute_detfeatures(unsigned short int *image, int rows, int cols, float microns,
   unsigned char *breastmask, float startx, float starty, float endx, float endy,
   int num_detections, REGION *detections, int radius, int deltar,
   DETFEATURE **detfeature);

#endif
