/*******************************************************************************
* PROGRAM: ICSio.h
* PURPOSE: This is the header file for ICSio.c
* NAME: Michael Heath, University of South Florida
* DATE: 9/12/97
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define JPEG_EXECUTABLE "/home/captiva1/heath/bin/jpeg"

/*******************************************************************************
* Definition of some data structures.
*******************************************************************************/
typedef struct{
   int day, month, year;
}DATE;

typedef struct{
   int image_exists;
   int rows;
   int cols;
   int bitsperpixel;
   float resolution;
   int overlay_exists;
   char compressed_filename[200];
   char uncompressed_filename[200];
   char overlay_filename[200];
   int did_i_decompress_it;
}ICSIMINFO;

typedef struct{
   float ics_version;
   char filename[200];
   DATE date_of_study;
   int patient_age;
   char media[40];
   char film_type[40];
   int density;
   DATE date_digitized;
   char digitizer_brand[80];
   char digitizer_model[80];
   char unknown[40];
   ICSIMINFO left_cc, left_mlo, right_cc, right_mlo;
}ICSDATA;

/*******************************************************************************
* Function prototypes.
*******************************************************************************/
int read_ics_file(char *filename, ICSDATA *icsdata, char *a_path);
int decompress_ics_images(ICSDATA *icsdata);
int decompress_ics_image(ICSDATA *icsdata, char *view);
int decompress_LJPEG_image(char *filename);
void ru_exit();
