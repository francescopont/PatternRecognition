/*******************************************************************************
* Program: mktemplate
* Purpose: This program is a tool for creating a template of the
* ground truth regions from an OVERLAY file for a case in the DDSM database.
* An 8-bit pgm image is created.
* (http://marathon.csee.usf.edu/Mammography/Database.html) 
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
* By default the template is binary with the background assigned a value of
* zero and all pixels in all rendered regions set to 255.
*
* The option -label can be used to assign different numbers to each rendered
* region. The regions are numbered from 1 to n and the background is assigned
* a value of zero. If one wants to view this image they will want to rescale
* it for display.
*
* Name: Michael Heath, University of South Florida
* Date: 1/21/2000
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../common/ICSio.h"
#include "../common/overlay.h"
#include "../common/mikesfileio.h"

#define VERSION "1.0.0"
#define VERSIONDATE "January 21, 2000"

#define PGM 1
#define TIFF 2

/*******************************************************************************
* Function prototypes.
*******************************************************************************/
static void print_help();
static int get_commandline_parameters(int argc, char *argv[], char **ics_filename,
   char **view, char *lesiontype, char **pathology, char **ms, char **mm,
   char **ct, char **cd, float *desired_resolution, int *use_all_regions, int *use_labels, int *fileformat);

char created_files[4][200];  /* This is used by read_ics_file(). */

int VERBOSE = 0;

void polyscan_translated_outline(OUTLINE outline, unsigned char **image,
   int *rows, int *cols, int *offset_r, int *offset_c);

int main(int argc, char *argv[])
{
   FILE *fp;
   char *ics_filename=NULL, *overlay_filename=NULL;
   char lesiontype[20], *view = NULL;
   int status;
   ICSDATA icsdata;
   OVERLAY_DATA overlay_data;
   int d, r, c, t, rr, cc;
   ABNORMALITY abnormality;
   REGION *truth_region=NULL;
   int truth_regions_valid = 0;
   char *ms = NULL, *mm = NULL, *ct = NULL, *cd = NULL;
   float resolution, desired_resolution, scale;
   int rows, cols, template_rows, template_cols, cp, gtpix, overlapping_regions=0;
   unsigned char *template=NULL;
   char template_filename[200] = {'\0'};
   int use_all_regions=0, use_labels=0, label;
   int fileformat = PGM;
   int *isvalid=NULL;
   char *pathology = NULL;

   int write_simple_tiff_uchar(char *filename, unsigned char *image, int rows, int cols);

   memset(&overlay_data, 0, sizeof(OVERLAY_DATA));

   /****************************************************************************
   * Get the command line parameters.
   ****************************************************************************/
   status = get_commandline_parameters(argc, argv, &ics_filename,
               &view, lesiontype, &pathology, &ms, &mm, &ct, &cd, &desired_resolution,
               &use_all_regions, &use_labels, &fileformat);
   if(status == 0) exit(1);

   /****************************************************************************
   * If we are running in a verbose mode, print out the command line info.
   ****************************************************************************/
   if(VERBOSE){
      printf("\n\n************************************************************\n");
      printf(" The mktemplate program is running in verbose mode.\n");
      printf("************************************************************\n\n");
      if(ics_filename != NULL)
         printf("   file.ics    : %s\n", ics_filename);
      if(view != NULL)
         printf("   view        : %s\n", view);
      if(use_all_regions == 1){
         printf("   You selected the option to use all regions.\n");
      }
      if(use_labels == 1){
         printf("   The regions will be numbered 1 to n.\n");
      }
   }

   /****************************************************************************
   * Read in the ics file.
   ****************************************************************************/
   status = read_ics_file(ics_filename, &icsdata, (char *)NULL);
   if(status == 0) exit(1);

   /****************************************************************************
   * Read in the overlay file.
   ****************************************************************************/
   overlay_filename = NULL;
   if((strcmp(view, "LEFT_CC") == 0) && (icsdata.left_cc.overlay_exists == 1)){
      overlay_filename = icsdata.left_cc.overlay_filename;
      resolution = icsdata.left_cc.resolution;
      rows = icsdata.left_cc.rows;
      cols = icsdata.left_cc.cols;
      if(fileformat == PGM) sprintf(template_filename, "%s.template.pgm", icsdata.left_cc.uncompressed_filename);
      else sprintf(template_filename, "%s.template.tif", icsdata.left_cc.uncompressed_filename);
   }
   else if((strcmp(view, "LEFT_MLO") == 0) && (icsdata.left_mlo.overlay_exists == 1)){
      overlay_filename = icsdata.left_mlo.overlay_filename;
      resolution = icsdata.left_mlo.resolution;
      rows = icsdata.left_mlo.rows;
      cols = icsdata.left_mlo.cols;
      if(fileformat == PGM) sprintf(template_filename, "%s.template.pgm", icsdata.left_mlo.uncompressed_filename);
      else sprintf(template_filename, "%s.template.tif", icsdata.left_mlo.uncompressed_filename);
   }
   else if((strcmp(view, "RIGHT_CC") == 0) && (icsdata.right_cc.overlay_exists == 1)){
      overlay_filename = icsdata.right_cc.overlay_filename;
      resolution = icsdata.right_cc.resolution;
      rows = icsdata.right_cc.rows;
      cols = icsdata.right_cc.cols;
      if(fileformat == PGM) sprintf(template_filename, "%s.template.pgm", icsdata.right_cc.uncompressed_filename);
      else sprintf(template_filename, "%s.template.tif", icsdata.right_cc.uncompressed_filename);
   }
   else if((strcmp(view, "RIGHT_MLO") == 0) && (icsdata.right_mlo.overlay_exists == 1)){
      overlay_filename = icsdata.right_mlo.overlay_filename;
      resolution = icsdata.right_mlo.resolution;
      rows = icsdata.right_mlo.rows;
      cols = icsdata.right_mlo.cols;
      if(fileformat == PGM) sprintf(template_filename, "%s.template.pgm", icsdata.right_mlo.uncompressed_filename);
      else sprintf(template_filename, "%s.template.tif", icsdata.right_mlo.uncompressed_filename);
   }

   if(overlay_filename != NULL){
      status = read_overlay_file(overlay_filename, &overlay_data);
      if(status == 0) exit(1);

      /*************************************************************************
      * Determine the scale factor to resize the template of the ground truth regions.
      *************************************************************************/
      if(desired_resolution <= 0.0) desired_resolution = resolution;
      scale = resolution / desired_resolution;

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

      /****************************************************************************
      * If we are not going to use all of the regions, we must check through all
      * of the abnormalities to determine which ones we will use. We are not
      * using any abnormalities that are not "MALIGNANT" of "BENIGN". Also we
      * will only use regions that match the lesiontype and the ms, mm, ct and cd. 
      ****************************************************************************/
      if(use_all_regions != 1){
         determine_valid_regions(overlay_data, lesiontype, pathology, ms, mm,
            ct, cd, &isvalid);
      }
      else if(overlay_data.total_abnormalities != 0){
         isvalid = (int *) calloc(overlay_data.total_abnormalities, sizeof(int));
         for(d=0;d<overlay_data.total_abnormalities;d++) isvalid[d] = 1;
      }

      /*************************************************************************
      * Allocate an array of REGION structures.
      *************************************************************************/
      if((truth_region = (REGION *) calloc(overlay_data.total_abnormalities, sizeof(REGION))) == NULL){
         fprintf(stderr, "Calloc error!\n");
         exit(1);
      }

      /*************************************************************************
      * Scan convert each ground truth region that we will be checking for.
      * Do not use any ground truth outline that is not "MALIGNANT" or "BENIGN".
      * For example we will skip regions that are "BENIGN_WITHOUT_CALLBACK".
      * We will also skip "MASS" regions if the user specified "-l CALCIFICATION",
      * or "CALCIFICATION" regions if the user specified "-l MASS".
      *************************************************************************/
      for(d=0;d<overlay_data.total_abnormalities;d++){
         abnormality = overlay_data.abnormalities[d];

         truth_region[d].matched = 0;
         if(isvalid[d] == 1){
            truth_region[d].valid = 1;
            truth_regions_valid++;
         }

         /**********************************************************************
         * Scan convert the outline.
         **********************************************************************/
         if(truth_region[d].valid == 1){

            /*******************************************************************
            * Scale all of the coordinates in the boundary of the abnormality
            * to be in the new scale.
            *******************************************************************/
            truth_region[d].outline = abnormality.boundary; 
            truth_region[d].outline.start.r = (int)floor(scale * truth_region[d].outline.start.r);
            truth_region[d].outline.start.c = (int)floor(scale * truth_region[d].outline.start.c);
            
            for(cp=0;cp<truth_region[d].outline.length;cp++){
               truth_region[d].outline.chain[cp].r = (int)floor(scale * truth_region[d].outline.chain[cp].r);
               truth_region[d].outline.chain[cp].c = (int)floor(scale * truth_region[d].outline.chain[cp].c);
            }
            polyscan_translated_outline(truth_region[d].outline, &(truth_region[d].image),
               &(truth_region[d].rows), &(truth_region[d].cols),
               &(truth_region[d].offset_r), &(truth_region[d].offset_c));
         }

         if(VERBOSE){

            if(truth_region[d].valid == 1){
               printf("                   [rows = %d] [cols = %d]",
                  truth_region[d].rows, truth_region[d].cols);

               printf("[offset_r = %d] [offset_c = %d]\n",
                  truth_region[d].offset_r, truth_region[d].offset_c);
            }
         }
      }
   }
   else{
      if(VERBOSE){
         printf("   There is no overlay file for the %s image.\n", view);
         exit(1);
      }
   }

   /****************************************************************************
   * Determine the size of the image we need to use for rendering the regions.
   * Then we will allocate a template image with those dimensions.
   ****************************************************************************/
   template_rows = (int)ceil(scale * rows);
   template_cols = (int)ceil(scale * cols);

   if((template = (unsigned char *) calloc(template_rows * template_cols, sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Calloc error of the template image!\n");
      exit(1);
   }

   /****************************************************************************
   * If there are ground truth regions that are valid (i.e. are the type we
   * are interested in) then we must see how many of them the detections found.
   ****************************************************************************/
   label = 1;
   if(truth_regions_valid != 0){
      for(t=0;t<overlay_data.total_abnormalities;t++){

         if(truth_region[t].valid == 1){

            for(r=0;r<truth_region[t].rows;r++){
               rr = truth_region[t].offset_r + r;
               for(c=0;c<truth_region[t].cols;c++){
                  cc = truth_region[t].offset_c +c;

                  gtpix = truth_region[t].image[r * truth_region[t].cols + c];

                  if(gtpix != 0){   /* We are inside of this region. */
                     if(template[rr*template_cols+cc] != 0){
                        if(overlapping_regions != 0){
                           fprintf(stderr, "Caution! Some ground truth regions overlap!\n");
                           overlapping_regions = 1;
                        }
                     }
                     if(use_labels) template[rr*template_cols+cc] = label;
                     else template[rr*template_cols+cc] = 255;
                  }
               }
            }
            label++;
         }
      }
   }

   /****************************************************************************
   * Write out the template to a file.
   ****************************************************************************/
   if(template_filename[0] != '\0'){

      if(fileformat == PGM){
         if((fp = fopen(template_filename, "wb")) == NULL){
            fprintf(stderr, "Error opening the file %s for writing!\n\n", template_filename);
            exit(1);
         }
         fprintf(fp, "P5\n%d %d\n255\n", template_cols, template_rows);
         fwrite(template, sizeof(unsigned char), template_rows*template_cols, fp);
         fclose(fp);
      }
      else{
         write_simple_tiff_uchar(template_filename, template, template_rows, template_cols);
      }

   }

   return(1);
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
   char **ct, char **cd, float *desired_resolution, int *use_all_regions, int *use_labels, int *fileformat)
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
   *desired_resolution = 0.0;
   *use_all_regions = 0;
   *use_labels = 0;
   *pathology = NULL;
   strcpy(lesiontype, "NONE");

   *ms = NULL;
   *mm = NULL;
   *ct = NULL;
   *cd = NULL;

   /****************************************************************************
   * Extract the command line parameters.
   ****************************************************************************/
   for(i=1;i<argc;i++){
      if(strcmp(argv[i], "-ics") == 0){ *ics_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-view") == 0){ *view = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-resolution") == 0){ *desired_resolution = atof(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-v") == 0) VERBOSE = 1;
      else if(strncmp(argv[i], "-tiff", 4) == 0) *fileformat = TIFF;
      else if(strcmp(argv[i], "-label") == 0) *use_labels = 1;
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
   if((*ics_filename == NULL) || (*view==NULL)){
      fprintf(stderr, "mktemplate: Error! You did not specify all needed parameters.\n");

      /*************************************************************************
      * Make sure that the view is valid.
      *************************************************************************/
      if(!((strcmp(*view, "LEFT_CC") == 0) || (strcmp(*view, "LEFT_MLO") == 0) ||
         (strcmp(*view, "RIGHT_CC") == 0) || (strcmp(*view, "RIGHT_MLO") == 0))){
         fprintf(stderr, "mktemplate: Error! Invalid view specified (%s).\n", *view);
         *view = NULL;
      }

      if(*ics_filename != NULL)
         fprintf(stderr, "   file.ics = %s\n", *ics_filename);

      if(*view != NULL)
         fprintf(stderr, "   view =     %s\n", *view);

      return(0);
   }

   /****************************************************************************
   * Make sure that the view is valid.
   ****************************************************************************/
   if(!((strcmp(*view, "LEFT_CC") == 0) || (strcmp(*view, "LEFT_MLO") == 0) ||
      (strcmp(*view, "RIGHT_CC") == 0) || (strcmp(*view, "RIGHT_MLO") == 0))){
      fprintf(stderr, "mktemplate: Error! Invalid view specified (%s).\n", *view);
      return(0);
   }

   return(1);
}

/*******************************************************************************
* Function: print_help
* Purpose: To print out the command line help.
* Name: Michael Heath, University of South Florida
* Date: 11/13/98
*******************************************************************************/
static void print_help()
{
   printf("\n\n********************************************************************************\n");
   printf("This program can be used to create a template of all, or select regions marked\n");
   printf("in the ground truth (OVERLAY) files in a DDSM case. By default, all MALIGNANT\n");
   printf("and BENIGN regions are rendered in the template and the template is rendered\n");
   printf("with the resolution of the mammograms in the case.\n");
   printf("********************************************************************************\n\n");
   printf("<USAGE> mktemplate [-l lesiontype ] [-ms DESC] [-mm DESC] [-ct DESC] [-cd DESC]\n");
   printf("                   [-pathology PATHOLOGY] [-v] [-version] [-all] [-label]\n");
   printf("                   [-resolution #] [-tif] -ics file.ics -view VIEW\n\n");
   printf("   -ics        An ics file from the DDSM database.\n");
   printf("   -view       The view to process RIGHT_MLO or RIGHT_CC or LEFT_MLO or LEFT_CC.\n");
   printf("   -resolution The resolution of the template you create (in microns).\n");
   printf("   -version    Print the version of the software and exit.\n");
   printf("   -v          Run the program in verbose mode.\n");
   printf("   -l          Skip regions that are not lesions you specify (CALCIFICATION\n");
   printf("               or MASS)\n");
   printf("   -pathology  Skip regions without the specified PATHOLOGY. Valid values\n");
   printf("               are MALIGNANT and BENIGN.\n");
   printf("   -label      Assign labels 1 to n for all rendered regions.\n");
   printf("   -all        Use all of the regions (overrides -l, -ms, -mm, -ct and -cd).\n");
   printf("   -tif        Write the template out in TIFF format instead of PGM format.\n\n");
   printf("   -ms DESC-DESC-...-DESC    Skip masses without these shapes.\n");
   printf("         \"ROUND\",\"OVAL\",\"LOBULATED\",\"IRREGULAR\"\n");
   printf("         \"ARCHITECTURAL_DISTORTION\",\"TUBULAR\",\"LYMPH_NODE\"\n");
   printf("         \"ASYMMETRIC_BREAST_TISSUE\",\"FOCAL_ASYMMETRIC_DENSITY\"\n");
   printf("   -mm DESC-DESC-...-DESC    Skip masses without these margins.\n");
   printf("         \"CIRCUMSCRIBED\",\"MICROLOBULATED\",\"OBSCURED\",\"ILL_DEFINED\"\n");
   printf("         \"SPICULATED\"\n");
   printf("   -ct DESC-DESC-...-DESC    Skip calcifications without this type.\n");
   printf("         \"PUNCTATE\",\"AMORPHOUS\",\"PLEOMORPHIC\",\"ROUND_AND_REGULAR\"\n");
   printf("         \"LUCENT_CENTER\",\"FINE_LINEAR_BRANCHING\",\"SKIN\",\"VASCULAR\"\n");
   printf("         \"COARSE\",\"LARGE_RODLIKE\",\"EGGSHELL\",\"MILK_OF_CALCIUM\"\n");
   printf("         \"SUTURE\",\"DYSTROPHIC\"\n");
   printf("   -cd DESC-DESC-...-DESC    Skip calcifications without this distribution.\n");
   printf("         \"CLUSTERED\",\"LINEAR\",\"SEGMENTAL\",\"REGIONAL\"\n");
   printf("         \"DIFFUSELY_SCATTERED\"\n\n");
}
