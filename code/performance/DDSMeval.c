/*******************************************************************************
* Program: DDSMeval
* Purpose: This program is a tool for evaluating the results from a lesion
* detection algorithm that has been run on a case from the Digital Database
* for Screening Mammography (DDSM) at the University of South Florida.
* (http://marathon.csee.usf.edu/Mammography/Database.html) 
*
* The evaluation detection algorithm works by comparing each detected region
* to each ground truth region (specified in the OVERLAY file). A match is
* found when the centroid of the detected region falls within a ground truth
* region. Multiple detections can match the same ground truth region, and a
* single detection can match multiple ground truth regions.
*
* After the matching is done, all unmatched ground truth regions are counted
* as false negatives (misses), all matched regions are counted as true
* positives (hits) and all unmatched detected regions are counted as false
* alarms.
*
* The program is run by specifying the needed files on the command line, and
* an optional option. The option "-l lesiontype" can be specified as either
* "-l CALCIFICATION" or "-l MASS". One of these will limit the ground truth
* regions to the specified type. This allows the evaluation of an algorithm
* that detects only masses or calcifications. The default is to evaluate
* detections of both types of lesions.
*
* The command line option -max # was added to allow the user to only read the
* first so many (#) detections in the file. This can allow sampling for an FROC
* curve if the "threshold parameter" is the number of detections.
*
* The optional command line parameter -one was added to allow only one detected
* region to match with a ground truth region. The second, third, fourth etc.
* regions that match the ground truth region are counted as false alarms.
*
* The command line parameter -st was added to allow the selection of a
* subset of detections using a floating point number on the same line as
* the PROMPT or POLYGON. This allows an alternative approach of selecting
* some of the detections based on a suspiciousness threshold, rather than
* by specifying a maximum number of detections per image.
*
* Name: Michael Heath, University of South Florida
* Date: 1/20/2000
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
#define VERSIONDATE "January 21, 2000"

/*******************************************************************************
* Function prototypes.
*******************************************************************************/
static void print_help();
static int get_commandline_parameters(int argc, char *argv[], char **ics_filename,
   char **view, char **detection_filename, char *lesiontype, char **pathology,
   char **ms, char **mm, char **ct, char **cd, int *max, int *match_one_to_one,
   float *suspiciousness_threshold);

char created_files[4][200];  /* This is used by read_ics_file(). */
int VERBOSE = 0;

void polyscan_translated_outline(OUTLINE outline, unsigned char **image,
   int *rows, int *cols, int *offset_r, int *offset_c);

int main(int argc, char *argv[])
{
   char *ics_filename=NULL, *overlay_filename=NULL, *detection_filename=NULL;
   char lesiontype[20], *view = NULL;
   int status;
   ICSDATA icsdata;
   OVERLAY_DATA overlay_data;
   int num_detections = 0;
   REGION *detections=NULL;
   int d, p, r, c, i, t, k, match;
   unsigned long int sum_x, sum_y, count;
   ABNORMALITY abnormality;
   REGION *truth_region=NULL;
   int truth_regions_valid = 0;
   int false_positives = 0, false_negatives = 0, true_positives = 0;
   int target_r, target_c;
   char *ms = NULL, *mm = NULL, *ct = NULL, *cd = NULL;
   int max=100;
   int match_one_to_one = 0;
   float suspiciousness_threshold;
   int detection_rows, detection_cols;
   int rows, cols;
   int *isvalid = NULL;
   char *pathology = NULL;

   memset(&overlay_data, 0, sizeof(OVERLAY_DATA));

   /****************************************************************************
   * Get the command line parameters.
   ****************************************************************************/
   status = get_commandline_parameters(argc, argv, &ics_filename,
               &view, &detection_filename, lesiontype, &pathology, &ms, &mm,
               &ct, &cd, &max, &match_one_to_one, &suspiciousness_threshold);
   if(status == 0) exit(1);

   /****************************************************************************
   * If we are running in a verbose mode, print out the command line info.
   ****************************************************************************/
   if(VERBOSE){
      printf("\n\n************************************************************\n");
      printf(" The DDSMeval program is running in verbose mode.\n");
      printf("************************************************************\n\n");
      if(ics_filename != NULL)
         printf("   file.ics    : %s\n", ics_filename);
      if(view != NULL)
         printf("   view        : %s\n", view);
      if(detection_filename != NULL)
         printf("   file.det    : %s\n", detection_filename);
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
      rows = icsdata.left_cc.rows;
      cols = icsdata.left_cc.cols;
   }
   else if((strcmp(view, "LEFT_MLO") == 0) && (icsdata.left_mlo.overlay_exists == 1)){
      overlay_filename = icsdata.left_mlo.overlay_filename;
      rows = icsdata.left_mlo.rows;
      cols = icsdata.left_mlo.cols;
   }
   else if((strcmp(view, "RIGHT_CC") == 0) && (icsdata.right_cc.overlay_exists == 1)){
      overlay_filename = icsdata.right_cc.overlay_filename;
      rows = icsdata.right_cc.rows;
      cols = icsdata.right_cc.cols;
   }
   else if((strcmp(view, "RIGHT_MLO") == 0) && (icsdata.right_mlo.overlay_exists == 1)){
      overlay_filename = icsdata.right_mlo.overlay_filename;
      rows = icsdata.right_mlo.rows;
      cols = icsdata.right_mlo.cols;
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

      /****************************************************************************
      * Determine which abnormalities we will use. We are not using any
      * abnormalities that are not "MALIGNANT" of "BENIGN". Also we will only
      * use regions that match the lesiontype and the ms, mm, ct and cd.
      ****************************************************************************/
      determine_valid_regions(overlay_data, lesiontype, pathology, ms, mm, ct, cd, &isvalid);

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
            truth_region[d].outline = abnormality.boundary; 
            polyscan_translated_outline(truth_region[d].outline, &(truth_region[d].image),
               &(truth_region[d].rows), &(truth_region[d].cols),
               &(truth_region[d].offset_r), &(truth_region[d].offset_c));
         }

         if(VERBOSE){

            /*******************************************************************
            * If we wanted to, we could print out the scan converted ground
            * truth regions that are valid to files as images.
            ********************************************************************
            * FILE *fptmp;
            * char tmpfilename[30];
            *
            * if(truth_region[d].valid == 1){
            *    sprintf(tmpfilename, "%s_offsettruth_%02d.pgm", view, d+1);
            *    fptmp = fopen(tmpfilename, "wb");
            *    fprintf(fptmp, "P5\n%d %d\n255\n", truth_region[d].cols, truth_region[d].rows);
            *    fwrite(truth_region[d].image, 1, truth_region[d].cols * truth_region[d].rows, fptmp);
            *    fclose(fptmp);
            * }
            *******************************************************************/

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
      }
   }

   /****************************************************************************
   * Read in the detection file.
   ****************************************************************************/
   status = read_detection_file(detection_filename, &num_detections, &detections,
               &detection_rows, &detection_cols, max, suspiciousness_threshold);
   if(status == 0) exit(1);

   /****************************************************************************
   * Rescale the dimensions to be in the same scale as the ground truth.
   ****************************************************************************/
   convert_detection_coordinates(detection_rows, detection_cols, num_detections,
      detections, rows, cols);

   for(d=0;d<num_detections;d++){
      detections[d].valid = 1;
      detections[d].matched = 0;
   }

   if(VERBOSE){
      printf("\n   There are %d detected regions.\n\n", num_detections);
      for(d=0;d<num_detections;d++){
         if(detections[d].outline.length == 0){
            printf("      PROMPT %f\n", detections[d].suspiciousness);
            /* printf("%d %d\n", (int)detections[d].centroid_c, (int)detections[d].centroid_r); */
         }
         else{
            printf("      POLYGON(%d) %f\n", detections[d].outline.length, detections[d].suspiciousness);
            /*
            for(p=0;p<detections[d].outline.length;p++){
               printf("%d %d\n", detections[d].outline.chain[p].c, detections[d].outline.chain[p].r);
            }
            */
         }
      }
      printf("\n");
   }

   /****************************************************************************
   * Scan convert each detected region that is a polygon. Also compute the
   * centroid of each of these regions.
   ****************************************************************************/
   for(d=0;d<num_detections;d++){
      if(detections[d].outline.length != 0){

         /**********************************************************************
         * Scan convert the outline.
         **********************************************************************/
         polyscan_translated_outline(detections[d].outline, &(detections[d].image),
            &(detections[d].rows), &(detections[d].cols),
            &(detections[d].offset_r), &(detections[d].offset_c));

         /**********************************************************************
         * Compute the centroid of the outline.
         **********************************************************************/
         sum_x = sum_y = count = 0;
         for(r=0;r<detections[d].rows;r++){
            for(c=0;c<detections[d].cols;c++){
               if(detections[d].image[r*detections[d].cols+c] == 255){
                  sum_x += c;
                  sum_y += r;
                  count++;
               }
            }
         }
         detections[d].centroid_c = detections[d].offset_c + (double)(sum_x) / (double)(count);
         detections[d].centroid_r = detections[d].offset_r + (double)(sum_y) / (double)(count);

         if(VERBOSE){

            printf("      POLYGON_%02d:  [centroid_r = %.1f] [centroid_c = %.1f]\n",
               d+1, detections[d].centroid_r, detections[d].centroid_c);

            printf("                   [rows = %d] [cols = %d]",
               detections[d].rows, detections[d].cols);

            printf("[offset_r = %d] [offset_c = %d]\n",
               detections[d].offset_r, detections[d].offset_c);
         }
      }
      else{
         if(VERBOSE){
            printf("      PROMPT%02d: [centroid_r = %.1f] [centroid_c = %.1f]\n",
               d+1, detections[d].centroid_r, detections[d].centroid_c);
         }
      }
   }

   if(VERBOSE) printf("\n");

   /****************************************************************************
   * If there are ground truth regions that are valid (i.e. are the type we
   * are interested in) then we must see how many of them the detections found.
   ****************************************************************************/
   if(truth_regions_valid != 0){
      for(d=0;d<num_detections;d++){
         for(t=0;t<overlay_data.total_abnormalities;t++){

	    abnormality = overlay_data.abnormalities[t];
            match = 0;

            if(truth_region[t].valid == 1){

               target_r = detections[d].centroid_r - truth_region[t].offset_r;
               target_c = detections[d].centroid_c - truth_region[t].offset_c;

               /****************************************************************
               * Check to see if the centroid of the detected region falls
               * inside a valid ground truth region.
               ****************************************************************/
               if((target_r >= 0) && (target_r < truth_region[t].rows) &&
                  (target_c >= 0) && (target_c < truth_region[t].cols)){

                  if(truth_region[t].image[target_r * truth_region[t].cols + target_c] != 0){
                     match = 1;
                  }
               }

               if(match == 1){
                  detections[d].matched++;
                  truth_region[t].matched++;
               }
            }
         }
      }
   }

   /****************************************************************************
   * Now that we have matched up all of the detected regions and ground truth
   * regions, we must count up how many valid regions of each type had matches.
   ****************************************************************************/
   false_positives = false_negatives = true_positives = 0;
   for(d=0;d<num_detections;d++){
      if(detections[d].valid == 1){
         if(detections[d].matched == 0) false_positives++;
      }
   }
   for(d=0;d<overlay_data.total_abnormalities;d++){
      if(truth_region[d].valid == 1){
         if(truth_region[d].matched == 0) false_negatives++;
         else true_positives++;
      }
   }

   /****************************************************************************
   * Print out the results.
   ****************************************************************************/
   if(match_one_to_one == 0){
      printf("%2d %2d %2d #", true_positives, false_positives, false_negatives);
   }
   else printf("%2d %2d %2d #", true_positives, num_detections-true_positives, false_negatives);

   for(i=1;i<argc;i++) printf(" %s", argv[i]);
   printf("\n");

   if(VERBOSE) printf("\n");
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
   char **view, char **detection_filename, char *lesiontype, char **pathology,
   char **ms, char **mm, char **ct, char **cd, int *max, int *match_one_to_one,
   float *suspiciousness_threshold)
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
   *detection_filename = (char *)NULL;
   strcpy(lesiontype, "NONE");

   *pathology = NULL;
   *ms = NULL;
   *mm = NULL;
   *ct = NULL;
   *cd = NULL;
   *suspiciousness_threshold = 0.0;

   /****************************************************************************
   * Extract the command line parameters.
   ****************************************************************************/
   for(i=1;i<argc;i++){
      if(strcmp(argv[i], "-ics") == 0){ *ics_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-det") == 0){ *detection_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-view") == 0){ *view = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-version") == 0){
         printf("\n\n%s Version: %s %s\n\n", argv[0], VERSION, VERSIONDATE);
         exit(1);
      }
      else if(strcmp(argv[i], "-max") == 0){ *max = atoi(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-example") == 0){
         print_example_det_file();
         exit(1);
      }
      else if(strcmp(argv[i], "-st") == 0){ *suspiciousness_threshold = atof(argv[i+1]); i++; }
      else if(strcmp(argv[i], "-v") == 0) VERBOSE = 1;
      else if(strcmp(argv[i], "-one") == 0) *match_one_to_one = 1;
      else if(strcmp(argv[i], "-pathology") == 0){
         if((i+1) < argc){
            if(strcmp(argv[i+1], "MALIGNANT") == 0)
               *pathology = argv[i+1];
            else if(strcmp(argv[i+1], "BENIGN") == 0)
               *pathology = argv[i+1];
            else{
               fprintf(stderr, "DDSMeval: Invalid pathology (MALIGNANT or BENIGN).\n");
               return(0);
            }
         }
         else{
            fprintf(stderr, "DDSMeval: Error! No pathology was specified.\n");
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
               fprintf(stderr, "DDSMeval: Invalid lesiontype (CALCIFICATION or MASS).\n");
               return(0);
            }
         }
         else{
            fprintf(stderr, "DDSMeval: Error! No lesiontype was specified.\n");
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
      /*
      else{
         if(argv[i][0] == '-'){
            fprintf(stderr, "DDSMeval: Invalid option '%s' specified.\n", argv[i]);
            return(0);
         }
         filenames_specified++;
         switch(filenames_specified){
            case 1: *ics_filename = argv[i]; break;
            case 2: *view = argv[i]; break;
            case 3: *detection_filename = argv[i]; break;
         }
      }
      */
   }

   if(strcmp(lesiontype, "NONE") == 0) strcpy(lesiontype, "BOTH");

   /****************************************************************************
   * Make sure the user entered all three required parameters.
   ****************************************************************************/
   if((*ics_filename == NULL) || (*detection_filename == NULL) || (*view==NULL)){
      fprintf(stderr, "DDSMeval: Error! You did not specify all needed parameters.\n");

      /*************************************************************************
      * Make sure that the view is valid.
      *************************************************************************/
      if(!((strcmp(*view, "LEFT_CC") == 0) || (strcmp(*view, "LEFT_MLO") == 0) ||
         (strcmp(*view, "RIGHT_CC") == 0) || (strcmp(*view, "RIGHT_MLO") == 0))){
         fprintf(stderr, "DDSMeval: Error! Invalid view specified (%s).\n", *view);
         *view = NULL;
      }

      if(*ics_filename != NULL)
         fprintf(stderr, "   file.ics = %s\n", *ics_filename);

      if(*view != NULL)
         fprintf(stderr, "   view =     %s\n", *view);

      if(*detection_filename != NULL)
         fprintf(stderr, "   file.det = %s\n", *detection_filename);

      return(0);
   }

   /****************************************************************************
   * Make sure that the view is valid.
   ****************************************************************************/
   if(!((strcmp(*view, "LEFT_CC") == 0) || (strcmp(*view, "LEFT_MLO") == 0) ||
      (strcmp(*view, "RIGHT_CC") == 0) || (strcmp(*view, "RIGHT_MLO") == 0))){
      fprintf(stderr, "DDSMeval: Error! Invalid view specified (%s).\n", *view);
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
   printf("This program is a tool for evaluating the results from a lesion detection\n");
   printf("algorithm that has been run on a case from the Digital Database for Screening\n");
   printf("Mammography (DDSM) at the University of South Florida.\n");
   printf("(http://marathon.csee.usf.edu/Mammography/Database.html)\n");
   printf("The detection results must be placed in a file with a specified format.\n");
   printf("To see an example of the format of the file, run this program with the\n");
   printf("option -example. The output of this program is the number of true positives,\n");
   printf("false positives and false negatives in that order. Please note that only\n");
   printf("BENIGN and MALIGNANT lesions matching the provided lesion type and description\n");
   printf("are considered to be ground truth regions that should be detected.\n");
   printf("********************************************************************************\n\n");

   printf("<USAGE> DDSMeval [-l lesiontype ] [-ms DESC] [-mm DESC] [-ct DESC] [-cd DESC]\n");
   printf("                 [-pathology PATHOLOGY] [-v] [-version] [-example] [-max #]\n");
   printf("                 [-st #] [-one] -ics file.ics -view VIEW -det file.det\n\n");
   printf("   -ics        An ics file from the DDSM database.\n");
   printf("   -view       The view to process RIGHT_MLO, RIGHT_CC, LEFT_MLO or LEFT_CC.\n");
   printf("   -version    Print the version of the software and exit.\n");
   printf("   -v          Run the program in verbose mode.\n");
   printf("   -l          Skip regions that are not lesions you specify (CALCIFICATION\n");
   printf("               or MASS)\n");
   printf("   -det        The name of a detection file (If non-existent, 0 detections\n");
   printf("               are assumed).\n");
   printf("   -max        Use up to this many detections from the file, starting at the\n");
   printf("               top.\n");
   printf("   -pathology  Skip regions without the specified PATHOLOGY. Valid values\n");
   printf("               are MALIGNANT and BENIGN.\n");
   printf("   -st         Only use detections with a suspiciousness value>=this threshold.\n");
   printf("   -one        Only count the the first detection mapped to a ground truth\n");
   printf("               region as a true positive, all others mapped to the region are\n");
   printf("               false positives. When not using this option, multiple detections\n");
   printf("               mapping to the same ground truth region count as one true\n");
   printf("               positive detection and zero false positive detections.\n");
   printf("   -example    Prints out an example of the detection file format and exits.\n\n");
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
