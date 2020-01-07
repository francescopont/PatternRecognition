/*******************************************************************************
* File: overlay.h
* Purpose: This file contains header information for the file overlay.c.
* Name: Michael Heath, University of South Florida
* Date: 5/6/98
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#ifndef _OVERLAY
#define _OVERLAY

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_CORES 4         /* This is just here to limit memory consumption. */
#define MAX_DESCRIPTIONS 2  /* Each abnormality can has 2 LESION_TYPES. */

/*******************************************************************************
* Definitions of the data structures to use to hold the overlay data.
*******************************************************************************/
enum LESION_TYPE {MASS, CALCIFICATION, OTHER};

enum MASS_SHAPE {ROUND, OVAL, LOBULATED, IRREGULAR, ARCHITECTURAL_DISTORTION,
                 TUBULAR, LYMPH_NODE, ASYMMETRIC_BREAST_TISSUE, FOCAL_ASYMMETRIC_DENSITY};
enum MASS_MARGINS {CIRCUMSCRIBED, MICROLOBULATED, OBSCURED, ILL_DEFINED, SPICULATED};
enum CALCIFICATION_TYPE {PUNCTATE, AMORPHOUS, PLEOMORPHIC, ROUND_AND_REGULAR, LUCENT_CENTER, FINE_LINEAR_BRANCHING,
                         SKIN, VASCULAR, COARSE, LARGE_RODLIKE, EGGSHELL, MILK_OF_CALCIUM, SUTURE, DYSTROPHIC};
enum CALCIFICATION_DISTRIBUTION {CLUSTERED, LINEAR, SEGMENTAL, REGIONAL, DIFFUSELY_SCATTERED};

typedef struct{   /* The start point of a chain code and each additional   */
   int r, c;      /* point on the boundary are stored with this structure. */
}CHAINPOINT;      

typedef struct{       /* This data structure is used to store an outline.   */
   CHAINPOINT start;  /* It could be a boundary or a core. The chaincode   */
   int length;        /* is converted to points and the array of points is  */
   CHAINPOINT *chain; /* stored. The start point is stored in start and in the array. */
}OUTLINE;

union SHAPETYPE{      /* A MASS has a shape and a CALCIFICATION has a type. */
   short int shape, type;
};

union MARGINDISTRIBUTION{ /* A MASS has a margins and a CALCIFICATION has a distribution. */
   short int margin, distribution;
};

typedef struct{                  /* A lesion can be a MASS or a CALCIFICATION. The */
   enum LESION_TYPE lesion_type; /* lesion_type tells which type this is. Two      */
   union SHAPETYPE shapetype;    /* unions are used to store the (shape or type)   */
   union MARGINDISTRIBUTION      /* and the (margin or distribution).              */
	    margindistribution;
}DESCRIPTION;

typedef struct{
   int num_descriptions;
   DESCRIPTION description[MAX_DESCRIPTIONS];
   int assessment;
   int subtlety;
   char pathology[24];
   OUTLINE boundary;
   int num_cores;
   OUTLINE core[MAX_CORES];
}ABNORMALITY;

typedef struct{
   int total_abnormalities;
   ABNORMALITY *abnormalities;
}OVERLAY_DATA;

static char *MASS_SHAPE_string[] =
	{"ROUND","OVAL","LOBULATED","IRREGULAR","ARCHITECTURAL_DISTORTION", 
         "TUBULAR", "LYMPH_NODE", "ASYMMETRIC_BREAST_TISSUE", "FOCAL_ASYMMETRIC_DENSITY", NULL};
static char *MASS_MARGIN_string[] =
	{"CIRCUMSCRIBED","MICROLOBULATED","OBSCURED","ILL_DEFINED","SPICULATED", NULL};
static char *CALCIFICATION_TYPE_string[] =
        {"PUNCTATE","AMORPHOUS","PLEOMORPHIC","ROUND_AND_REGULAR","LUCENT_CENTER","FINE_LINEAR_BRANCHING",
         "SKIN", "VASCULAR", "COARSE", "LARGE_RODLIKE", "EGGSHELL", "MILK_OF_CALCIUM", "SUTURE", "DYSTROPHIC", NULL};
static char *CALCIFICATION_DISTRIBUTION_string[] =
	{"CLUSTERED","LINEAR","SEGMENTAL","REGIONAL","DIFFUSELY_SCATTERED", NULL};

/*******************************************************************************
* Function prototypes.
*******************************************************************************/

int read_overlay_file(char *filename, OVERLAY_DATA *overlay_data);
int read_abnormality(FILE *fp, ABNORMALITY *abnormality, int abnormality_number);

int determine_valid_regions(OVERLAY_DATA overlay_data, char *lesiontype,
    char *pathology, char *ms, char *mm, char *ct, char *cd, int **isvalid);

int get_description(FILE *fp, DESCRIPTION *description);
void read_chaincode(FILE *fp, OUTLINE *outline, char *type);

int inlist(char **stringptr, char *string_list[]);
void get_keyword(FILE *fp, char *keyword);
int check_keyword(char *keyword, const char *truth);
int get_intval(FILE *fp);
void get_string(FILE *fp, char *thestring);

int valtobit(int value);
int bittoval(int bit);

void overlay_exit();

int write_overlay_file(char *filename, OVERLAY_DATA overlay_data);
int write_abnormality(FILE *fp, ABNORMALITY *abnormality, int abnormality_number);
int put_description(FILE *fp, DESCRIPTION *description);
void write_keywords(FILE *fp, int bitmask, char *string_list[]);
void write_chaincode(FILE *fp, OUTLINE *outline, char *type);

#endif
