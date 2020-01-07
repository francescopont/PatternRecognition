/*******************************************************************************
* File: optical_density.h
* PURPOSE: This file contains the data structures and function prototypes for
* optical_density.c.
* NAME: Michael Heath, University of South Florida
* DATE: 11/9/97
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#ifndef _OPTICAL_DENSITY
#define _OPTICAL_DENSITY

#include "virtual_image.h"
#include "myalloc.h"
#include <string.h>

#define MIN_DENSITY 0.05
#define MAX_DENSITY 3.60

#define OD_OFFSET 65535.0
/* #define OD_SCALE (floor(( OD_OFFSET / (MAX_DENSITY - MIN_DENSITY) ))) */
#define OD_SCALE (( OD_OFFSET / (MAX_DENSITY - MIN_DENSITY) ))

typedef struct{
   double od_offset, od_scale;
   USHORT *map;
}OD_REMAPPER_DATA;

void get_optical_density_segment(CACHEIM *image, int rpage, int cpage, int rseg, int cseg);
double roundval(double x);

void linear_optical_density(CACHEIM *image, CACHEIM *source_image, int max_depth,
   int max_cache_size, char *mode, double *od_offset, double *od_scale,
   double *min_density, double *max_density, double A, double B);

void log10_optical_density(CACHEIM *image, CACHEIM *source_image, int max_depth,
   int max_cache_size, char *mode, double *od_offset, double *od_scale,
   double *min_density, double *max_density, double A, double B);

void map_with_ushort_lut(CACHEIM *image, CACHEIM *source_image, int max_depth,
   int max_cache_size, int datatype, char *mode, unsigned short int *lut);

unsigned short int od_to_scaled_od(double od);

#endif
