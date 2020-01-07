/*******************************************************************************
* Program: combinesusp
* Purpose: This program will take input individual detection files for MLO and
* CC views of the same breast, and will output new detection files that
* incorporate the suspicion of the detections in the original files. We
* will use the breast segmentation files in this process, because they
* contain information we require, such as the nipple position that is stored
* as one endpoint of the major axis.
*
* The transformation of individual suspicion values to augmented suspiciousness
* values that combine the detection information between views, should have the
* following properties:
*
*  (1) A suspiciousness value should not be changed if we do not expect to
*      find a corresponding position in the other view. Reasons include...
*        (a) For some reason we do not believe the major axis is correct in
*            either of both views. This would mean that we would have little
*            to no confidence in the matching of detections between views.
*        (b) A detection site in one view is in a position in the breast such
*            that we do not expect the tissue was imaged in the second view.
*  (2) A suspiciousness value may increase (or decrease) based on the closeness
*      of the match in the other view and the original suspicion in both views.
*
* One transformation that seems to accomplish this is as follows:
*    GIVEN: A set of sites in each projection with associated suspicion values,
*           distances to the nipple position.
*    The output suspiciousness value is equal to
*       original_suspicion + avarage suspicion * F(delta_radial_position).
*
* The function F should be non-negative valued, have a value slightly higher
* than 1 at 0 and should decrease with increasing delta_radial_position.
* We could use the maximum output suspiciousness value for all pairs of
* corresponding positions.
*
* Name: Michael Heath, University of South Florida
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
* Date: 6/9/2000
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mikesfileio.h"
#include "detfeatures.h"

#define VERSION "1.1.0"
#define VERSIONDATE "June 9, 2000"

#define MAX_DETECTIONS 30
#define SUSPTHRESH 0.0

int VERBOSE = 0;

int main(int argc, char *argv[])
{
   char *MLO_segfile=NULL, *MLO_detfile=NULL, *CC_segfile=NULL, *CC_detfile=NULL;
   char *MLO_reducedfile=NULL, *CC_reducedfile=NULL;
   int MLO_num_detections = 0, CC_num_detections = 0;
   int MLO_detection_rows=0, MLO_detection_cols=0;
   int CC_detection_rows=0, CC_detection_cols=0;
   REGION *MLO_detections=NULL, *CC_detections=NULL;
   float microns;
   int status;
   int MLO_numpoints = 0, MLO_segmentation_rows=0, MLO_segmentation_cols=0;
   float *MLO_xcoord=NULL, *MLO_ycoord=NULL, MLO_startx=0.0, MLO_starty=0.0,
      MLO_endx=0.0, MLO_endy=0.0;
   int CC_numpoints = 0, CC_segmentation_rows=0, CC_segmentation_cols=0;
   float *CC_xcoord=NULL, *CC_ycoord=NULL, CC_startx=0.0, CC_starty=0.0,
      CC_endx=0.0, CC_endy=0.0;
   float length_MLO, length_CC;
   int i, r, c, m;
   double alpha=0.0, gamma=0.0;
   int num=0;
   char *MLO_templatefile = NULL, *CC_templatefile = NULL;
   unsigned char *MLO_template=NULL, *CC_template=NULL;
   int MLO_template_rows = 0, MLO_template_cols = 0;
   int CC_template_rows = 0, CC_template_cols = 0;
   int *MLO_orig_rank=NULL, *CC_orig_rank=NULL;
   char augment[10];
   DETFEATURE *MLO_detfeature=NULL, *CC_detfeature=NULL;
   unsigned short int *MLO_image=NULL, *CC_image=NULL;
   int MLO_rows, MLO_cols, CC_rows, CC_cols;
   unsigned char *MLO_breastmask=NULL, *CC_breastmask=NULL;

   int read_pgm_image_16bit(char *infilename, unsigned short int **image, int *rows,
      int *cols);
   void augment_suspicion(int view1_num_detections,
       DETFEATURE *view1_detfeature, int view2_num_detections,
       DETFEATURE *view2_detfeature, float view1_length,
       float view2_length, double alpha, double gamma, int num);
   void sort_detfeature(int num_detections, DETFEATURE *detfeature, int **orig_rank);
   void write_segmented_detections_file(char *old_detection_filename, int rows, int cols,
      int num_detections, DETFEATURE *detfeature);
   void compute_detfeatures(unsigned short int *image, int rows, int cols, float microns,
      unsigned char *breastmask, float startx, float starty, float endx, float endy,
      int num_detections, REGION *detections, int radius, int deltar,
      DETFEATURE **detfeature);
   void scan_convert_breast_region(int numpoints, float *xpoints, float *ypoints,
       int rows, int cols, unsigned char **breast_region, float microns,
       float startx, float endx, char *filename);

   /****************************************************************************
   * Get the command line parameters.
   ****************************************************************************/
   if(argc < 11){
      fprintf(stderr, "\n\n<USAGE> %s MLO.reduced.pgm MLO.sgt ", argv[0]);
      fprintf(stderr, "MLO.det CC.reduced.pgm CC.sgt CC.det microns alpha gamma num\n\n");
      exit(1);
   }
   MLO_reducedfile = argv[1];
   MLO_segfile = argv[2];
   MLO_detfile = argv[3];
   CC_reducedfile = argv[4];
   CC_segfile = argv[5];
   CC_detfile = argv[6];
   microns = atof(argv[7]);
   alpha = atof(argv[8]);
   gamma = atof(argv[9]);
   num = atoi(argv[10]);

   /****************************************************************************
   * Read in the 16-bit reduced PGM images.
   ****************************************************************************/
   if(read_pgm_image_16bit(MLO_reducedfile, &MLO_image, &MLO_rows, &MLO_cols) == 0) exit(1);
   if(read_pgm_image_16bit(CC_reducedfile, &CC_image, &CC_rows, &CC_cols) == 0) exit(1);

   /****************************************************************************
   * Read in the segmentation files.
   ****************************************************************************/
   status = read_segmentation_file(MLO_segfile, &MLO_numpoints, &MLO_xcoord, &MLO_ycoord,
               &MLO_startx, &MLO_starty, &MLO_endx, &MLO_endy,
               &MLO_segmentation_rows, &MLO_segmentation_cols);
   if(status == 0) exit(1);

   status = read_segmentation_file(CC_segfile, &CC_numpoints, &CC_xcoord, &CC_ycoord,
               &CC_startx, &CC_starty, &CC_endx, &CC_endy,
               &CC_segmentation_rows, &CC_segmentation_cols);
   if(status == 0) exit(1);

   /****************************************************************************
   * Read in the detection files.
   ****************************************************************************/
   status = read_detection_file(MLO_detfile, &MLO_num_detections, &MLO_detections,
               &MLO_detection_rows, &MLO_detection_cols, (int)MAX_DETECTIONS, (float)SUSPTHRESH);
   if(status == 0) exit(1);

   status = read_detection_file(CC_detfile, &CC_num_detections, &CC_detections,
               &CC_detection_rows, &CC_detection_cols, (int)MAX_DETECTIONS, (float)SUSPTHRESH);
   if(status == 0) exit(1);

   /****************************************************************************
   * Convert the segmentation coordinates to be in the same scale as the
   * detection coordinates.
   ****************************************************************************/
   convert_segmentation_coordinates(MLO_segmentation_rows, MLO_segmentation_cols,
      MLO_numpoints, MLO_xcoord, MLO_ycoord, &MLO_startx, &MLO_starty,
      &MLO_endx, &MLO_endy, MLO_detection_rows, MLO_detection_cols);

   convert_segmentation_coordinates(CC_segmentation_rows, CC_segmentation_cols,
      CC_numpoints, CC_xcoord, CC_ycoord, &CC_startx, &CC_starty,
      &CC_endx, &CC_endy, CC_detection_rows, CC_detection_cols);

   /****************************************************************************
   * Scan convert the breast region.
   ****************************************************************************/
   scan_convert_breast_region(MLO_numpoints, MLO_xcoord, MLO_ycoord,
       MLO_detection_rows, MLO_detection_cols, &MLO_breastmask, microns,
       MLO_startx, MLO_endx, MLO_segfile);
   scan_convert_breast_region(CC_numpoints, CC_xcoord, CC_ycoord,
       CC_detection_rows, CC_detection_cols, &CC_breastmask, microns,
       CC_startx, CC_endx, CC_segfile);

   /****************************************************************************
   * Compute the lengths of the axes.
   ****************************************************************************/
   length_MLO = sqrt((double)(MLO_startx - MLO_endx) * (double)(MLO_startx - MLO_endx) +
                     (double)(MLO_starty - MLO_endy) * (double)(MLO_starty - MLO_endy));
   length_MLO *= microns;
   length_CC = sqrt((double)(CC_startx - CC_endx) * (double)(CC_startx - CC_endx) +
                     (double)(CC_starty - CC_endy) * (double)(CC_starty - CC_endy));
   length_CC *= microns;

   /****************************************************************************
   * Compute the distance of each detection from the nipple position and store
   * the distance in microns.
   ****************************************************************************/
   compute_detfeatures(MLO_image, MLO_detection_rows, MLO_detection_cols, microns,
      MLO_breastmask, MLO_startx, MLO_starty, MLO_endx, MLO_endy,
      MLO_num_detections, MLO_detections, 30, 10, &MLO_detfeature);
   compute_detfeatures(CC_image, CC_detection_rows, CC_detection_cols, microns,
      CC_breastmask, CC_startx, CC_starty, CC_endx, CC_endy,
      CC_num_detections, CC_detections, 30, 10, &CC_detfeature);

   /****************************************************************************
   * Compute augmented suspiciousness values.
   ****************************************************************************/
   augment_suspicion(MLO_num_detections, MLO_detfeature,
      CC_num_detections, CC_detfeature, length_MLO, length_CC,
      alpha, gamma, num);
   augment_suspicion(CC_num_detections, CC_detfeature,
      MLO_num_detections, MLO_detfeature, length_CC, length_MLO,
      alpha, gamma, num);

   /****************************************************************************
   * If the axes lengths are not within 2.0cm of eachother in length, replace the
   * augmented suspiciousness values with the original suspiciousness values.
   * This may be done because we are uncertain on the correspondence when the
   * axes are not nearly the same length.
   ****************************************************************************/
   if(fabs(length_MLO - length_CC) > 20000.0){
      sprintf(augment, "NO ");
      /* printf("Skipping the augmentation on %s %s\n", MLO_detfile, CC_detfile); */
      for(i=0;i<MLO_num_detections;i++) MLO_detfeature[i].aug_suspiciousness = MLO_detfeature[i].suspiciousness;
      for(i=0;i<CC_num_detections;i++) CC_detfeature[i].aug_suspiciousness = CC_detfeature[i].suspiciousness;
   }
   else sprintf(augment, "YES");

   /****************************************************************************
   * Sort all of the detection data for the MLO and CC images using the augmented
   * suspicion as the key. Sort in decreasing order.
   ****************************************************************************/
   sort_detfeature(MLO_num_detections, MLO_detfeature, &MLO_orig_rank);
   sort_detfeature(CC_num_detections, CC_detfeature, &CC_orig_rank);

   /****************************************************************************
   * Write out new detection files.
   ****************************************************************************/
   write_segmented_detections_file(MLO_detfile, MLO_detection_rows, MLO_detection_cols,
      MLO_num_detections, MLO_detfeature);
   write_segmented_detections_file(CC_detfile, CC_detection_rows, CC_detection_cols,
      CC_num_detections, CC_detfeature);

   /****************************************************************************
   * Print out the information for each detection site.
   ****************************************************************************/
   if(0){
      int maxpoints, i, j;
      char MLO_string[100], CC_string[100];

      maxpoints = MLO_num_detections;
      if(maxpoints < CC_num_detections) maxpoints = CC_num_detections;

      for(i=0;i<maxpoints;i++){

         MLO_string[0] = '\0';
         CC_string[0] = '\0';

         if(i < MLO_num_detections){
            sprintf(MLO_string, "%6.1f %6.1f %6.1f %9.2f %6.1f %7.3f", MLO_detfeature[i].centroid_c,
               MLO_detfeature[i].centroid_r, MLO_detfeature[i].suspiciousness,
               MLO_detfeature[i].nipple_distance/10000.0, MLO_detfeature[i].aug_suspiciousness,
               MLO_detfeature[i].angle, MLO_detfeature[i].median);
         }
         if(i < CC_num_detections){
            sprintf(CC_string, "%6.1f %6.1f %6.1f %9.2f %6.1f %7.3f", CC_detfeature[i].centroid_c,
               CC_detfeature[i].centroid_r, CC_detfeature[i].suspiciousness,
               CC_detfeature[i].nipple_distance/10000.0, CC_detfeature[i].aug_suspiciousness,
               CC_detfeature[i].angle, CC_detfeature[i].median);
         }

         /* printf("%2d) %40s %40s\n", i+1, MLO_string, CC_string); */
      }
   }
}

/*******************************************************************************
* Function: write_segmented_detections_file
* Purpose: Write out the detections with the new, augmented suspiciousness
* values. Please note that the list of detections has been resorted in
* decreasing order of the augemented_suspiciousness values.
* Name: Michael Heath, University of South Florida
* Date: 4/17/2000
*******************************************************************************/
void write_segmented_detections_file(char *old_detection_filename, int rows, int cols,
   int num_detections, DETFEATURE *detfeature)
{
   FILE *fpdet=NULL;
   char detection_filename[200] = {'\0'};
   char nd_filename[200] = {'\0'};
   int p;

   /****************************************************************************
   * Write the detection coordinates out to a file.
   ****************************************************************************/
   sprintf(detection_filename, "%s.augdet", old_detection_filename);
   fpdet = fopen(detection_filename, "w");
   fprintf(fpdet, "# SOURCE FILENAME: %s\n", old_detection_filename);
   fprintf(fpdet, "%d %d", cols, rows);
   for(p=0;p<num_detections;p++){
      fprintf(fpdet, "\nPROMPT %f\n%d %d", detfeature[p].aug_suspiciousness,
         (int)detfeature[p].centroid_c, (int)detfeature[p].centroid_r); 
   }
   fclose(fpdet);
}

/*******************************************************************************
* Function: sort_detfeature
* Purpose: To sort the detections, augsuspicion and nd in decreasing order
* using augsuspicion as the key.
* Name: Michael Heath, University of South Florida
* Date: 4/17/2000
*******************************************************************************/
void sort_detfeature(int num_detections, DETFEATURE *detfeature, int **orig_rank)
{
   int i, j, tempint;
   DETFEATURE tempdetfeature;

   *orig_rank = (int *) calloc(num_detections, sizeof(int));

   for(i=0;i<num_detections;i++){
      (*orig_rank)[i] = i;
   }

   for(i=0;i<(num_detections-1);i++){
      for(j=(i+1);j<num_detections;j++){

         if(detfeature[i].aug_suspiciousness < detfeature[j].aug_suspiciousness){

            tempdetfeature = detfeature[i];
            detfeature[i] = detfeature[j];
            detfeature[j] = tempdetfeature;

            tempint = (*orig_rank)[i];
            (*orig_rank)[i] = (*orig_rank)[j];
            (*orig_rank)[j] = tempint;
         }

      }
   }
}

/*******************************************************************************
* Function: augment_suspicion
* Purpose: We have a list of detections in two views. Here we want to compute
* a new suspiciousness value for each detection in view one by combining the
* detections in the two views.
* Name: Michael Heath, University of South Florida
* Date: 4/17/2000
*******************************************************************************/
void augment_suspicion(int view1_num_detections,
    DETFEATURE *view1_detfeature, int view2_num_detections,
    DETFEATURE *view2_detfeature, float view1_length,
    float view2_length, double alpha, double gamma, int num)
{
   int i, j, view1_num, view2_num;
   double d, susp_view1, susp_view2, suspicion, max_suspicion;
   int num_augmented = 0, num_skipped = 0;

   double update_suspicion(double susp1, double susp2, double d,
      double alpha, double gamma);

   /****************************************************************************
   * Cycle through at most num detections in view 1 and in view 2.
   ****************************************************************************/
   view1_num = view1_num_detections;
   if(view1_num > num) view1_num = num;
   view2_num = view2_num_detections;
   if(view2_num > num) view2_num = num;

   /****************************************************************************
   * Compute the maximum augmented suspicion for each detection site.
   ****************************************************************************/
   for(i=0;i<view1_num;i++){

      suspicion = max_suspicion = 0.0;
      if((view1_detfeature[i].nipple_distance < view1_length) && (view1_detfeature[i].nipple_distance < view2_length)){
         susp_view1 = (double)view1_detfeature[i].suspiciousness;
         for(j=0;j<view2_num;j++){
            d = fabs((double)view1_detfeature[i].nipple_distance - (double)view2_detfeature[j].nipple_distance);
            susp_view2 = (double)view2_detfeature[j].suspiciousness;

            suspicion = update_suspicion(susp_view1, susp_view2, d, alpha, gamma);
            if(suspicion > max_suspicion) max_suspicion = suspicion;
         }
         view1_detfeature[i].aug_suspiciousness = max_suspicion;
         num_augmented++;
      }
      else{
         view1_detfeature[i].aug_suspiciousness = (double)view1_detfeature[i].suspiciousness;
         num_skipped++;
      }
   }
   for(i=view1_num;i<view1_num_detections;i++)
      view1_detfeature[i].aug_suspiciousness = (double)view1_detfeature[i].suspiciousness;
   /* printf("AUG = %d  SKIP=%d\n", num_augmented, num_skipped); */
}

/*******************************************************************************
* Function: update_suspicion
* Purpose: This function takes two suspiciousness values and a distance and
* calculates an updated suspiciousness value based on those three quantities.
* There are two other parameters in the equation that are constants.
*
*           G * (S1+S2)/2            S1 = suspicion in view 1
*      S1 + -------------            S2 = suspicion in view 2
*           1 + (a*d*a*d)            d  = distance between detections
*     -------------------            S  = stretch factor (constant)
*             2                      G  = gain factor (constant)
*
* Name: Michael Heath, University of South Florida
* Date: 4/17/2000
*******************************************************************************/
double update_suspicion(double susp1, double susp2, double d, double alpha,
   double gamma)
{
   double average_susp, val;

   d /= 1000.0;               /* Convert d from microns to mm. */

   average_susp = (susp1 + susp2) / 2.0;
   val = (susp1 +  gamma * average_susp / (1.0 + (alpha*d*alpha*d))) / 2.0;
   
   return(val);
}

/*******************************************************************************
* Function: scan_convert_breast_region
* Purpose: To create an image that is a mask of the breast region from the
* segmentation data.
* Name: Michael Heath, University of South Florida
* Date: 5/15/2000
*******************************************************************************/
void scan_convert_breast_region(int numpoints, float *xpoints, float *ypoints,
    int rows, int cols, unsigned char **breast_region, float microns,
    float startx, float endx, char *filename)
{
   int n, r, c;
   int *xcoord= NULL, *ycoord=NULL;

   void polyscan_coords(int numpoints, int *x_coord, int *y_coord, unsigned char *image, int rows, int cols);

   /****************************************************************************
   * We will want to set all of the pixels outside of the breast to the
   * value zero. Since we already have a polygon of the breast area,
   * we will scan convert the polygon to an image. Recall that even though
   * the polygon may have been found in a subsampled version of the image,
   * the polygon verticies are saved in the full image resolution scale.
   ****************************************************************************/
   xcoord = (int *) mycalloc("xcoord", numpoints, sizeof(int));
   ycoord = (int *) mycalloc("ycoord", numpoints, sizeof(int));

   for(n=0;n<numpoints;n++){
      xcoord[n] = (int)floor(xpoints[n] + 0.5);
      ycoord[n] = (int)floor(ypoints[n] + 0.5);
   }

   if(((*breast_region) = (unsigned char *) calloc(rows * cols, sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Calloc error of image!\n");
      exit(1);
   }
   polyscan_coords(numpoints, xcoord, ycoord, (*breast_region), rows, cols);

   /****************************************************************************
   * Exclude some of the pixels from the side of the mammogram that the
   * patient is on. (Is this backwards ?)
   ****************************************************************************/
   if(endx >= startx){    /* apex == APEXLEFT */
      for(r=0;r<rows;r++){
         for(c=0;c<(int)floor(5000.0 / (double)(microns));c++)
            (*breast_region)[r*cols+c] = 0;
      }
   }
   else{
      for(r=0;r<rows;r++){
         for(c=(cols-1-(int)floor(5000.0 / (double)(microns)));c<cols;c++)
            (*breast_region)[r*cols+c] = 0;
      }
   }

   /****************************************************************************
   * Write out the scan converted image of the breast region. This will be
   * needed by the detect program.
   ****************************************************************************/
   if(1){
      FILE *fptmp;
      char tmpfilename[200];

      sprintf(tmpfilename, "%s.polytest.pgm", filename);
      fptmp = fopen(tmpfilename, "wb");
      fprintf(fptmp, "P5\n%d %d\n255\n", cols, rows);
      fwrite((*breast_region), 1, rows*cols, fptmp);
      fclose(fptmp);
   }

   myfree("xcoord", xcoord);
   myfree("ycoord", ycoord);
}


/*******************************************************************************
* Function: read_pgm_image_16bit
* Purpose: To read an image from a file that stores it in (pseudo) 16-bit PGM
* format.
* Name: Michael Heath, University of South Florida
* Date: 5/15/2000
*******************************************************************************/
int read_pgm_image_16bit(char *infilename, unsigned short int **image, int *rows,
    int *cols)
{
   FILE *fp;
   char buf[71];
   int maxval = 0;

   /***************************************************************************
   * Open the input image file for reading if a filename was given. If no
   * filename was provided, set fp to read from standard input.
   ***************************************************************************/
   if(infilename == NULL) fp = stdin;
   else{
      if((fp = fopen(infilename, "r")) == NULL){
         fprintf(stderr, "Error reading the file %s in read_pgm_image_16bit().\n",
            infilename);
         return(0);
      }
   }

   /***************************************************************************
   * Verify that the image is in PGM format, read in the number of columns
   * and rows in the image and scan past all of the header information.
   ***************************************************************************/
   fgets(buf, 70, fp);
   if(strncmp(buf,"P5",2) != 0){
      fprintf(stderr, "The file %s is not in PGM format in ", infilename);
      fprintf(stderr, "read_pgm_image_16bit().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
   sscanf(buf, "%d %d", cols, rows);
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

   sscanf(buf, "%d", &maxval);

   if(maxval <= 255){
      fprintf(stderr, "Error reading the image in read_pgm_image_16bit().\n");
      exit(1);
   }

   /***************************************************************************
   * Allocate memory to store the image then read the image from the file.
   ***************************************************************************/
   if(((*image) = (unsigned short int *) malloc((*rows)*(*cols)*sizeof(unsigned short int))) == NULL){
      fprintf(stderr, "Memory allocation failure in read_pgm_image_16bit().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }
   fseek(fp, -(*rows)*(*cols)*sizeof(unsigned short int), 2);
   if((*rows) != fread((*image), sizeof(unsigned short int)*(*cols), (*rows), fp)){
      fprintf(stderr, "Error reading the image data in read_pgm_image_16bit().\n");
      if(fp != stdin) fclose(fp);
      free((*image));
      return(0);
   }

   if(fp != stdin) fclose(fp);
   return(1);
}

