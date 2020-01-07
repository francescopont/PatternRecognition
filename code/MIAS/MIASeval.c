/*******************************************************************************
* Program: MIASeval
* Purpose: This program is a tool for evaluating the results from a lesion
* detection algorithm that has been run on a case from the MIAS database.
*
* The evaluation detection algorithm works by comparing each detected region
* to each ground truth region (specified in a file that was already created).
* A match is
* found when the centroid of the detected region falls within a ground truth
* region. Multiple detections can match the same ground truth region, and a
* single detection can match multiple ground truth regions.
*
* After the matching is done, all unmatched ground truth regions are counted
* as false negatives (misses), all matched regions are counted as true
* positives (hits) and all unmatched detected regions are counted as false
* alarms.
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
* Date: 3/19/2000
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
#define VERSIONDATE "March 19, 2000"

#define MAXGTREGIONS 20

typedef struct{
   char image_filename[200], bgtissue[10], abnormality_class[10], severity[10];
   int rows, cols;
   int row, col, radius;
   int valid;
   int matched;
}MIASGT;

/*******************************************************************************
* Function prototypes.
*******************************************************************************/
static void print_help();
static int get_commandline_parameters(int argc, char *argv[], 
   char **detection_filename, char **gt_filename, int *max,
   int *match_one_to_one, float *suspiciousness_threshold);

void polyscan_translated_outline(OUTLINE outline, unsigned char **image,
   int *rows, int *cols, int *offset_r, int *offset_c);

int VERBOSE = 0;

int main(int argc, char *argv[])
{
   FILE *fpgt=NULL;
   char *detection_filename=NULL, *gt_filename=NULL;
   int status;
   int num_detections = 0;
   REGION *detections=NULL;
   int d, p, r, c, i, t, k, match;
   unsigned long int sum_x, sum_y, count;
   MIASGT truth_region[MAXGTREGIONS];
   int truth_regions_valid = 0;
   int false_positives = 0, false_negatives = 0, true_positives = 0;
   int target_r, target_c;
   int max=100;
   int match_one_to_one = 0;
   float suspiciousness_threshold;
   int detection_rows, detection_cols;
   int rows, cols;
   int *isvalid = NULL;
   char gtline[100];
   float distance;
   int num_truth_regions=0;

   memset(truth_region, 0, MAXGTREGIONS*sizeof(MIASGT));

   /****************************************************************************
   * Get the command line parameters.
   ****************************************************************************/
   status = get_commandline_parameters(argc, argv, &detection_filename, 
               &gt_filename, &max, &match_one_to_one, &suspiciousness_threshold);
   if(status == 0) exit(1);

   /****************************************************************************
   * If we are running in a verbose mode, print out the command line info.
   ****************************************************************************/
   if(VERBOSE){
      printf("\n\n************************************************************\n");
      printf(" The MIASeval program is running in verbose mode.\n");
      printf("************************************************************\n\n");
      if(detection_filename != NULL)
         printf("   file.det    : %s\n", detection_filename);
      if(gt_filename != NULL)
         printf("   file.gt     : %s\n", gt_filename);
   }

   /****************************************************************************
   * Read in the ground truth file.
   ****************************************************************************/
   num_truth_regions = 0;
   if((fpgt = fopen(gt_filename, "r")) != NULL){
      while((fgets(gtline, 100, fpgt) != NULL) && (!feof(fpgt))){
         if((gtline[0] != '#') && (strstr(gtline, "NORM") == NULL)){
            sscanf(gtline, "%s %d %d %s %s %s %d %d %d",
               truth_region[num_truth_regions].image_filename,
               &(truth_region[num_truth_regions].cols),
               &(truth_region[num_truth_regions].rows),
               truth_region[num_truth_regions].bgtissue,
               truth_region[num_truth_regions].abnormality_class,
               truth_region[num_truth_regions].severity,
               &(truth_region[num_truth_regions].col),
               &(truth_region[num_truth_regions].row),
               &(truth_region[num_truth_regions].radius));
            truth_region[num_truth_regions].valid = 1;
            num_truth_regions++;
         }
      }
      fclose(fpgt);
   }

   /****************************************************************************
   * Read in the detection file.
   ****************************************************************************/
   status = read_detection_file(detection_filename, &num_detections, &detections,
               &detection_rows, &detection_cols, max, suspiciousness_threshold);
   if(status == 0) exit(1);

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
            printf("      POLYGON(%d) %f\n", detections[d].outline.length,
               detections[d].suspiciousness);
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
   * Rescale the dimensions to be in the same scale as the ground truth.
   ****************************************************************************/
   if(num_truth_regions != 0){

      if(VERBOSE){
         for(d=0;d<num_truth_regions;d++){
            printf("%s %d %d %d\n", truth_region[d].abnormality_class,
               truth_region[d].row, truth_region[d].col, truth_region[d].radius);
         }

      }

      convert_detection_coordinates(detection_rows, detection_cols, num_detections,
         detections, truth_region[0].rows, truth_region[0].cols);

      /*************************************************************************
      * Scan convert each detected region that is a polygon. Also compute the
      * centroid of each of these regions.
      *************************************************************************/
      for(d=0;d<num_detections;d++){
         if(detections[d].outline.length != 0){

            /*******************************************************************
            * Scan convert the outline.
            *******************************************************************/
            polyscan_translated_outline(detections[d].outline, &(detections[d].image),
               &(detections[d].rows), &(detections[d].cols),
               &(detections[d].offset_r), &(detections[d].offset_c));

            /*******************************************************************
            * Compute the centroid of the outline.
            *******************************************************************/
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
   }

   /****************************************************************************
   * If there are ground truth regions that are valid (i.e. are the type we
   * are interested in) then we must see how many of them the detections found.
   ****************************************************************************/
   truth_regions_valid = num_truth_regions;
   if(truth_regions_valid != 0){
      for(d=0;d<num_detections;d++){
         for(t=0;t<num_truth_regions;t++){

            match = 0;

            if(truth_region[t].valid == 1){

               distance = sqrt((float)(truth_region[t].row - detections[d].centroid_r) *
                               (float)(truth_region[t].row - detections[d].centroid_r) +
                               (float)(truth_region[t].col - detections[d].centroid_c) *
                               (float)(truth_region[t].col - detections[d].centroid_c));

               /****************************************************************
               * Check to see if the centroid of the detected region falls
               * inside a valid ground truth region.
               ****************************************************************/
               if(distance < truth_region[t].radius){
                  match = 1;
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
   for(d=0;d<num_truth_regions;d++){
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
static int get_commandline_parameters(int argc, char *argv[], 
   char **detection_filename, char **gt_filename, int *max,
   int *match_one_to_one, float *suspiciousness_threshold)
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
   *detection_filename = (char *)NULL;
   *gt_filename = (char *)NULL;
   *suspiciousness_threshold = 0.0;

   /****************************************************************************
   * Extract the command line parameters.
   ****************************************************************************/
   for(i=1;i<argc;i++){
      if(strcmp(argv[i], "-det") == 0){ *detection_filename = argv[i+1]; i++; }
      else if(strcmp(argv[i], "-gt") == 0){ *gt_filename = argv[i+1]; i++; }
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
      /*
      else{
         if(argv[i][0] == '-'){
            fprintf(stderr, "MIASeval: Invalid option '%s' specified.\n", argv[i]);
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

   /****************************************************************************
   * Make sure the user entered all the required parameters.
   ****************************************************************************/
   if((*gt_filename == NULL) || (*detection_filename == NULL)){
      fprintf(stderr, "MIASeval: Error! You did not specify all needed parameters.\n");

      if(*detection_filename != NULL)
         fprintf(stderr, "   file.det = %s\n", *detection_filename);

      if(*gt_filename != NULL)
         fprintf(stderr, "   file.gt  = %s\n", *gt_filename);

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
   printf("algorithm that has been run on a case from the MIAS database.\n");
   printf("The detection results must be placed in a file with a specified format.\n");
   printf("To see an example of the format of the file, run this program with the\n");
   printf("option -example. The output of this program is the number of true positives,\n");
   printf("false positives and false negatives in that order.\n");
   printf("********************************************************************************\n\n");

   printf("<USAGE> MIASeval [-v] [-version] [-example] [-max #]\n");
   printf("                 [-st #] [-one] -gt file.gt -det file.det\n\n");
   printf("   -version    Print the version of the software and exit.\n");
   printf("   -v          Run the program in verbose mode.\n");
   printf("   -l          Skip regions that are not lesions you specify (CALCIFICATION\n");
   printf("               or MASS)\n");
   printf("   -det        The name of a detection file (If non-existent, 0 detections\n");
   printf("               are assumed).\n");
   printf("   -gt         The name of a ground truth file (If non-existent, 0 GT\n");
   printf("               regions are assumed).\n");
   printf("   -max        Use up to this many detections from the file, starting at the\n");
   printf("               top.\n");
   printf("   -st         Only use detections with a suspiciousness value>=this threshold.\n");
   printf("   -one        Only count the the first detection mapped to a ground truth\n");
   printf("               region as a true positive, all others mapped to the region are\n");
   printf("               false positives. When not using this option, multiple detections\n");
   printf("               mapping to the same ground truth region count as one true\n");
   printf("               positive detection and zero false positive detections.\n");
   printf("   -example    Prints out an example of the detection file format and exits.\n\n");
}
