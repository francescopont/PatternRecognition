/*******************************************************************************
* Program: detect.c
* Purpose: This program segments a filtered mammogram into regions. Some of
* the regions may be selected as being abnormal. Any such region will be written
* out in a "detection" file. The segmentation of the file is done by a watershed
* transformation. The algorithm was modified slightly to use a mask to specify
* that the watershed transform should only be applied to the pixels in the
* segmented breast region. Also, the watersheds are relabeled with the labels of
* the nearest catchment basin. It is afterall, only the catchment basins that
* we are interested in for those are the segmented regions. A reference to the
* source of the algorithm is provided in a comment for the function where it is
* implemented.
* To Compile:
* Name: Michael Heath, University of South Florida
* Date: 11/10/99
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define VERBOSE 0
#define MAXDETECTIONS 30
#define VERSION "1.0.0"
#define VERSIONDATE "January 22, 2000"

#define REORDER_CLOSE_BASINS 1  /* set to zero to not reorder to reduce the   */
                                /* suspiciousness of basins close to other    */
                                /* suspicios basins.                          */
#define CLOSEDIST 5000.0        /* Two basin floors can not be closer than    */
				/* this distance in microns (i.e. 5000.0=5mm) */

int error_case = 0;

int main(int argc, char *argv[])
{
   FILE *fpin=NULL, *fpdet=NULL;
   int rows, cols, p, minmatch, max_basin, *xcoord=NULL, *ycoord=NULL;
   unsigned char *image=NULL, *breastmaskim=NULL;
   unsigned short int *wshedim=NULL;
   char *feature_filename=NULL, *breastmask_filename=NULL,
      *gt_filename=NULL, filename[200];
   unsigned char *gray_image = NULL;
   unsigned char *min_feature=NULL;
   unsigned long int *xplace=NULL, *yplace = NULL;
   char detection_filename[200];
   int junk1, junk2, k=0, i;
   float resolution = 0.0, aggfactor = 1.0;

   void watershed_transform(unsigned char *image, int rows, int cols,
       unsigned short int **wshedim, unsigned char *breastmaskim);
   int read_pgm_image(char *infilename, unsigned char **image, int *rows, int *cols);
   void find_min_value_per_basin(unsigned short int *wshedim, unsigned char *image,
      int rows, int cols, unsigned char **min_feature, int *max_basin, int **xccord, int **ycoord);
   void reorder_close_basins(unsigned short int *wshed, int rows, int cols,
      int *xcoord, int *ycoord, int pixdist);

   /****************************************************************************
   * Get the command line parameters.
   ****************************************************************************/
   for(i=1,k=0;i<argc;i++){
      if(strcmp(argv[i], "-version") == 0){
         printf("\n\n%s Version: %s %s\n\n", argv[0], VERSION, VERSIONDATE);
         exit(1);
      }
      else{
         k++;
         if(k == 1) feature_filename = argv[i];
         else if(k == 2) breastmask_filename = argv[i];
         else if(k == 3) resolution = atof(argv[i]);
      }
   }
   if(k != 3){

      printf("\n********************************************************************************\n");
      printf("This program will produce a detection file containing a list of detection\n");
      printf("sites in decreasing order of suspicion. Three input files are required. The\n");
      printf("feature image is an 8-bit PGM image with larger values indicating a higher\n");
      printf("suspicion of cancer. The polyscan image is an 8-bit pgm image with pixel\n");
      printf("values of zero indicating the background of the image (i.e. outside the\n");
      printf("breast tissue). These files are produced by the program afumfeature. You\n");
      printf("could create an 8-bit PGM file of the probability of suspiciousness of\n");
      printf("each pixel with your own algorithm instead. You could also create a file\n");
      printf("to take the place of the filename.polyscan.pgm image using the -sketchres\n");
      printf("option with the segment program. The resolution argument specifies the\n");
      printf("resolution of these images in microns.\n");
      printf("********************************************************************************\n");

      printf("\n<USAGE> %s [-version] filename.feature.pgm filename.polyscan.pgm\n", argv[0]);
      printf("                      resolution\n\n");
      exit(1);
   }

   /****************************************************************************
   * Read in the feature image. Then grey level invert it.
   ****************************************************************************/
   if(read_pgm_image(feature_filename, &image, &rows, &cols) == 0) exit(1);
   for(p=0;p<(rows*cols);p++) image[p] = 255 - image[p];

   /****************************************************************************
   * Read in the breastmask image.
   ****************************************************************************/
   if(read_pgm_image(breastmask_filename, &breastmaskim, &rows, &cols) == 0) exit(1);

   /****************************************************************************
   * Apply the watershed transform to segment the image.
   ****************************************************************************/
   watershed_transform(image, rows, cols, &wshedim, breastmaskim);

   /****************************************************************************
   * The basins are labeled in order of increasing number of the minimum value
   * in each catchment basin. This function finds the minimum value in each
   * catchment basin and returns those minima in an array.
   ****************************************************************************/
   find_min_value_per_basin(wshedim, image, rows, cols, &min_feature, &max_basin,
      &xcoord, &ycoord);

   /****************************************************************************
   * It is possible that detections can be very close together. If the
   * macro REORDER_CLOSE_BASINS has the value 1, we will reduce the
   * suspiciousness of basins that are close to other, more suspicious basins.
   ****************************************************************************/
   if(REORDER_CLOSE_BASINS){
      reorder_close_basins(wshedim, rows, cols, xcoord, ycoord, (int)floor(CLOSEDIST / resolution));

      free(min_feature);
      min_feature = NULL;
      free(xcoord);
      xcoord = NULL;
      free(ycoord);
      ycoord = NULL;

      /*************************************************************************
      * The basins are labeled in order of increasing number of the minimum value
      * in each catchment basin. This function finds the minimum value in each
      * catchment basin and returns those minima in an array.
      *************************************************************************/
      find_min_value_per_basin(wshedim, image, rows, cols, &min_feature, &max_basin,
         &xcoord, &ycoord);
   }

   /****************************************************************************
   * Write the detection coordinates out to a file.
   ****************************************************************************/
   sprintf(detection_filename, "%s.det", feature_filename);
   fpdet = fopen(detection_filename, "w");
   fprintf(fpdet, "# FEATURE FILENAME: %s\n", feature_filename);
   fprintf(fpdet, "# BREASTMASK FILENAME: %s\n", breastmask_filename);
   fprintf(fpdet, "# RESOLUTION: %f\n", resolution);
   fprintf(fpdet, "%d %d", cols, rows);
   for(p=1;p<=MAXDETECTIONS;p++){
      if(p <= max_basin){
         fprintf(fpdet, "\nPROMPT %f\n%d %d",
            (float)(255-image[ycoord[p]*cols+xcoord[p]]), xcoord[p], ycoord[p]);
      }
      /*
      if(p <= max_basin) fprintf(fpdet, "\nPROMPT\n%d %d",
         (int)floor(aggfactor*xcoord[p]), (int)floor(aggfactor*ycoord[p]));
      */
   }
   fclose(fpdet);

   /****************************************************************************
   * Write the image out to a file.
   ****************************************************************************/
   if(0){
      FILE *fp=NULL;

      sprintf(filename, "%s_watershed.pgm", feature_filename);
      if((fp = fopen(filename, "wb")) == NULL){
         fprintf(stderr, "Error opening the file named %s for writing.\n", filename);
         exit(1);
      }
      fprintf(fp, "P5\n%d %d\n65535\n", cols, rows);
      fwrite(wshedim, 2, rows*cols, fp);
      fclose(fp);
   }
}

/*******************************************************************************
* Function: watershed_transform
* Purpose: To compute the watershed transform of an image using the algorithm
* described in the paper "Watersheds in Digital Spaces: An Efficient Algorithm
* Based on Immersion Simulations" by Luc Vincent and Pierre Soille that
* appeared in IEEE Transactions on Pattern Analysis and Machine Intelligence,
* Vol 13, No. 6, June 1991, pp. 583-598.
* Name: Michael Heath, University of South Florida
* Date: 11/3/99
*******************************************************************************/
void watershed_transform(unsigned char *image, int rows, int cols,
    unsigned short int **wshedim, unsigned char *breastmaskim)
{
   unsigned long int *chist=NULL, *sortedplace=NULL;
   int numgrayvals = 256, r, c, pos, n, h, hmin, hmax, p, pprime, pprimeprime,
      cpos, rpos;
   int rval[4] = {-1,0,0,1}, cval[4] = {0,-1,1,0};
   unsigned long int *pixelpos=NULL;
   const short int WSHED=0, INIT=-1, MASK=-2;
   short int *im0=NULL;
   unsigned short int *imd = NULL;
   unsigned short int current_dist;
   short int current_label = 0, ficticious_pixel = -1;
   long int *queue=NULL;
   int qspos, qrpos, max, rr, cc, thisdist;
   unsigned long int hoffset, hnum;
   int wsmin, wsmax, minlabel;

   void fifo_add(long int *queue, long int value, int max, int *spos, int rpos);
   int fifo_empty(int spos, int rpos);
   long int fifo_first(long int *queue, int max, int spos, int *rpos);
   void view_queue(long int *queue, int max, int spos, int rpos);

   /****************************************************************************
   * Allocate memory for a circular queue that can store up to 1/4 of the
   * size of the image.
   ****************************************************************************/
   max = rows*cols / 4;
   queue = (long int *) calloc(max, sizeof(long int));
   qspos = qrpos = 0;

   /****************************************************************************
   * Compute the cumulative histogram of the image.
   ****************************************************************************/
   if(VERBOSE) printf("Computing the cumulative histogram.\n");
   chist = (unsigned long int *) calloc(numgrayvals, sizeof(unsigned long int));
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         if(breastmaskim[pos] != 0) chist[image[pos]]++;
      }
   }

   hmin = 0;
   while((hmin < numgrayvals) && (chist[hmin] == 0)) hmin++;
   hmax = (numgrayvals-1);
   while((hmax > 0) && (chist[hmax] == 0)) hmax--;

   for(r=1;r<numgrayvals;r++) chist[r] += chist[r-1];

   /****************************************************************************
   * Set up an array of pointers to pixel locations in the image that is sorted
   * according to increasing gray level of the pixels.
   ****************************************************************************/
   if(VERBOSE) printf("Setting up pointers into the image.\n");
   pixelpos = (unsigned long int *) calloc(rows*cols, sizeof(unsigned long int));
   sortedplace = (unsigned long int *) calloc(numgrayvals, sizeof(unsigned long int));
   for(r=1;r<numgrayvals;r++) sortedplace[r] = chist[r-1];
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         if(breastmaskim[pos] != 0){
            pixelpos[sortedplace[image[pos]]] = pos;
            sortedplace[image[pos]]++;
         }
      }
   }
   free(sortedplace);

   /****************************************************************************
   * Allocate the output watershed image (im0) and a working distance image (imd).
   ****************************************************************************/
   im0 = (short int *) calloc(rows*cols, sizeof(short int));
   imd = (unsigned short int *) calloc(rows*cols, sizeof(unsigned short int));

   for(r=0,pos=0;r<rows;r++) for(c=0;c<cols;c++,pos++) im0[pos] = INIT;
   /* Note that imd is initialized to 0 by the calloc(). */

   /* printf("hmin = %d  hmax = %d\n", hmin, hmax); */

   if(VERBOSE) printf("Performing the watershed transform.\n");
   for(h=hmin;h<=hmax;h++){

      if(h == hmin){
         hoffset = 0;
         hnum = chist[h];
      }
      else{
         hoffset = chist[h-1];
         hnum = chist[h] - chist[h-1];
      }

      /*
      printf("Flooding gray value %d. There are %lu pixels with this value starting at %lu.\n", h, hnum, hoffset);
      */

      if(hnum == 0) continue;     /* Skip non-occuring gray values. */

      for(pos=hoffset;pos<(hoffset+hnum);pos++){
         p = pixelpos[pos];
         im0[p] = MASK;

         r = p / cols;
         c = p % cols;

         for(n=0;n<4;n++){
            rpos = r + rval[n];
            cpos = c + cval[n];
            if((rpos>=0)&&(cpos>=0)&&(rpos<rows)&&(cpos<cols)){
               pprime = rpos*cols+cpos;
               if((breastmaskim[pprime]!=0) && ((im0[pprime] > 0) || (im0[pprime] == WSHED))){
                  imd[p] = 1;
                  fifo_add(queue, p, max, &qspos, qrpos);
               }
            }
         }
      }

      current_dist = 1;
      fifo_add(queue, ficticious_pixel, max, &qspos, qrpos);

      while(1){
         error_case = 1;
         p = fifo_first(queue, max, qspos, &qrpos);
         if(p == ficticious_pixel){
            if(fifo_empty(qspos, qrpos) == 1) break;
            else{
               fifo_add(queue, ficticious_pixel, max, &qspos, qrpos);
               current_dist = current_dist + 1;
               error_case = 2;
               p = fifo_first(queue, max, qspos, &qrpos);
            }
         }

         r = p / cols;
         c = p % cols;

         for(n=0;n<4;n++){
            rpos = r + rval[n];
            cpos = c + cval[n];
            if((rpos>=0)&&(cpos>=0)&&(rpos<rows)&&(cpos<cols)){
               pprime = rpos*cols+cpos;
               if((breastmaskim[pprime]!=0) && ((imd[pprime] < current_dist) && ((im0[pprime] > 0) || (im0[pprime] == WSHED)))){
                  if(im0[pprime] > 0){
                     if((im0[p] == MASK) || (im0[p] == WSHED)) im0[p] = im0[pprime];
                     else if(im0[p] != im0[pprime]) im0[p] = WSHED;
                  }
                  else if(im0[p] == MASK) im0[p] = WSHED;
               }
               else if((im0[pprime] == MASK) && (imd[pprime] == 0)){
                  imd[pprime] = current_dist + 1;
                  fifo_add(queue, pprime, max, &qspos, qrpos);
               }
            }
         }
      }

      for(pos=hoffset;pos<(hoffset+hnum);pos++){

         p = pixelpos[pos];

         imd[p] = 0;

         if(im0[p] == MASK){
            current_label = current_label+1;
            fifo_add(queue, p, max, &qspos, qrpos);
            im0[p] = current_label;
            while(fifo_empty(qspos, qrpos) == 0){

	       error_case = 3;
               pprime = fifo_first(queue, max, qspos, &qrpos);

               r = pprime / cols;
               c = pprime % cols;

               for(n=0;n<4;n++){
                  rpos = r + rval[n];
                  cpos = c + cval[n];
                  if((rpos>=0)&&(cpos>=0)&&(rpos<rows)&&(cpos<cols)){
                     pprimeprime = rpos*cols+cpos;
                     if((breastmaskim[pprimeprime]!=0) && (im0[pprimeprime] == MASK)){
                        fifo_add(queue, pprimeprime, max, &qspos, qrpos);
                        im0[pprimeprime] = current_label;
                     }
                  }
               }
            }
         }
      } 
   }

   free(chist);
   free(pixelpos);
   free(imd);

   if(VERBOSE){
      printf("Removing the watersheds.\n");
      fflush(stdout);
   }

   /****************************************************************************
   * Make sure that the queue is empty.
   ****************************************************************************/
   if(fifo_empty(qspos, qrpos) != 1){
      printf("The queue should be empty and it isn't.\n");
      exit(1); 
   }

   /****************************************************************************
   * Get rid of all of the watersheds. To do this we must process pixels with
   * the value WSHED in increasing order of their distance from the nearest
   * catchment basin. Therefore we will find the minimum distance of each
   * WSHED pixel to the nearest catchment basin. Distance is in terms of
   * pixels (not Euclidean distance.)
   ****************************************************************************/
   fifo_add(queue, ficticious_pixel, max, &qspos, qrpos);
   (*wshedim) = (unsigned short int *) calloc(rows*cols, sizeof(unsigned short int));
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
	 if((breastmaskim[pos]!=0) && (im0[pos] == WSHED)){
            fifo_add(queue, pos, max, &qspos, qrpos);
            (*wshedim)[pos] = rows+cols;
         }
      }
   }

   thisdist = 0;

   while(fifo_empty(qspos, qrpos) == 0){

      pos = fifo_first(queue, max, qspos, &qrpos);

      if(pos == ficticious_pixel){
         thisdist++;
         if(fifo_empty(qspos, qrpos) == 1) break;
         /* printf("Computing pixels with distance %d.\n", thisdist); */
         fifo_add(queue, ficticious_pixel, max, &qspos, qrpos);
         pos = fifo_first(queue, max, qspos, &qrpos);
      }

      r = pos / cols;
      c = pos % cols;

      minlabel = rows+cols;

      for(rr=(r-1);rr<=(r+1);rr++){
         for(cc=(c-1);cc<=(c+1);cc++){
            if((rr>=0)&&(cc>=0)&&(rr<rows)&&(cc<cols)){
               if((breastmaskim[rr*cols+cc]!= 0) && ((*wshedim)[rr*cols+cc] < minlabel)){
                  minlabel = (*wshedim)[rr*cols+cc];
               }
            }
         }
      }

      if(minlabel < thisdist) (*wshedim)[pos] = thisdist+1;
      else fifo_add(queue, pos, max, &qspos, qrpos);
   }

   /****************************************************************************
   * Make sure that the queue is empty.
   ****************************************************************************/
   if(fifo_empty(qspos, qrpos) != 1){
      printf("The queue should be empty and it isn't.\n");
      exit(1); 
   }

   /****************************************************************************
   * Find the maximum label. Also place all pixels having the label WSHED in
   * the circular FIFO queue. These are pixels that will need to be relabeled
   * with the value assigned to a nearby catchment basin. We first put a
   * ficticious pixel value in the queue. This will allow us to later tell
   * when we have gone though the queue because we may have put values
   * into the queue in the interum.
   ****************************************************************************/
   fifo_add(queue, ficticious_pixel, max, &qspos, qrpos);
   wsmin = wsmax = im0[0];
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){

         if(breastmaskim[pos] != 0){

	    if(im0[pos] == WSHED) fifo_add(queue, pos, max, &qspos, qrpos);

            if(im0[pos] < wsmin) wsmin = im0[pos];
            else{
               if(im0[pos] > wsmax) wsmax = im0[pos];
            }

         }
      }
   }
   /* printf("wsmin = %d   wsmax = %d\n", wsmin, wsmax); */

   /****************************************************************************
   * Get rid of all of the watersheds. To do this we must process pixels with
   * the value WSHED in increasing order of their distance from the nearest
   * catchment basin. Therefore we will find the minimum distance of each
   * WSHED pixel to the nearest catchment basin. Distance is in terms of
   * pixels (not Euclidean distance.)
   ****************************************************************************/
   thisdist = 0;

   while(fifo_empty(qspos, qrpos) == 0){

      pos = fifo_first(queue, max, qspos, &qrpos);

      if(pos == ficticious_pixel){    /* We looped through the queue. */
         thisdist++;
         if(fifo_empty(qspos, qrpos) == 1) break;
         /* printf("Discovering pixels with distance %d.\n", thisdist); */
         fifo_add(queue, ficticious_pixel, max, &qspos, qrpos);
         pos = fifo_first(queue, max, qspos, &qrpos);
      }

      if((*wshedim)[pos] == thisdist){

         r = pos / cols;
         c = pos % cols;

         minlabel = wsmax+1;
         for(rr=(r-1);rr<=(r+1);rr++){
            for(cc=(c-1);cc<=(c+1);cc++){
               if((rr>=0)&&(cc>=0)&&(rr<rows)&&(cc<cols)&&(breastmaskim[rr*cols+cc]!=0)){
                  if((*wshedim)[rr*cols+cc] < thisdist){
                     if(im0[rr*cols+cc] < minlabel) minlabel = im0[rr*cols+cc];
                  }
               }
            }
         }

         im0[pos] = minlabel;
      }

      else fifo_add(queue, pos, max, &qspos, qrpos);
   }

   /****************************************************************************
   * Copy the signed watershed transformed image into an unsigned datatype. Also,
   * remove the watersheds to extend the catchment basins to touch one another.
   * Each catchment basis does of course have a unique value.
   ****************************************************************************/
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         if(im0[pos] < 0) (*wshedim)[pos] = 0;
         else (*wshedim)[pos] = im0[pos];
      }
   }

   free(queue);
   free(im0);
}

/*******************************************************************************
* Function: fifo_add
* Purpose: This function adds a long integer to a FIFO circular queue.
* This code was taken pretty much from "C The Complete Reference, 2nd Edition".
* Name: Michael Heath, University of South Florida
* Date: 11/3/99
*******************************************************************************/
void fifo_add(long int *queue, long int value, int max, int *spos, int rpos)
{
   /****************************************************************************
   * The queue is full if either spos is one less than rpos
   * or if spos is at the end of the queue array and rpos
   * is at the beginning.
   ****************************************************************************/
   if((((*spos)+1) == rpos) || ((((*spos)+1) == max) && (!rpos))){
      printf("The queue is full!\n");
      exit(1);
   }

   /* printf("Placing %ld in the queue.\n", value); */

   queue[(*spos)] = value;
   (*spos)++;
   if((*spos) == max) (*spos) = 0;  /* Loop back */
}

/*******************************************************************************
* Function: fifo_empty
* Purpose: This function checks to see if a FIFO circular queue is empty.
* This code was taken pretty much from "C The Complete Reference, 2nd Edition".
* Name: Michael Heath, University of South Florida
* Date: 11/3/99
*******************************************************************************/
int fifo_empty(int spos, int rpos)
{
   if(spos == rpos) return(1);
   else return(0);
}

/*******************************************************************************
* Function: fifo_first
* Purpose: This function gets a long integer from a FIFO circular queue.
* This code was taken pretty much from "C The Complete Reference, 2nd Edition".
* Name: Michael Heath, University of South Florida
* Date: 11/3/99
*******************************************************************************/
long int fifo_first(long int *queue, int max, int spos, int *rpos)
{
   long int valtoreturn = 0;

   if((*rpos) == max) *rpos = 0;  /* Loop back. */
   if((*rpos) == spos){
      printf("No events to perform (%d).\n", error_case);
      exit(1);
   }
   (*rpos)++;
   /* return(queue[(*rpos)-1]); */
   valtoreturn = queue[(*rpos)-1];
   if((*rpos) == max) *rpos = 0;  /* Loop back. */
   return(valtoreturn);
}

/*******************************************************************************
* Function: view_queue
* Purpose: This is a little utility function that is used in debugging the code.
* It simple prints out the contents of the queue.
* Name: Michael Heath, University of South Florida
* Date: 11/3/99
*******************************************************************************/
void view_queue(long int *queue, int max, int spos, int rpos)
{
   printf("\n");
   while(fifo_empty(spos, rpos) == 0){
      printf("Q[%ld] ", fifo_first(queue, max, spos, &rpos));
   }
   printf("\n");
}

/******************************************************************************
* Function: read_pgm_image
* Purpose: This function reads in an image in PGM format. The image can be
* read in from either a file or from standard input. The image is only read
* from standard input when infilename = NULL. Because the PGM format includes
* the number of columns and the number of rows in the image, these are read
* from the file. Memory to store the image is allocated in this function.
* All comments in the header are discarded in the process of reading the
* image. Upon failure, this function returns 0, upon sucess it returns 1.
******************************************************************************/
int read_pgm_image(char *infilename, unsigned char **image, int *rows,
    int *cols)
{
   FILE *fp;
   char buf[71];

   /***************************************************************************
   * Open the input image file for reading if a filename was given. If no
   * filename was provided, set fp to read from standard input.
   ***************************************************************************/
   if(infilename == NULL) fp = stdin;
   else{
      if((fp = fopen(infilename, "r")) == NULL){
         fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
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
      fprintf(stderr, "read_pgm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
   sscanf(buf, "%d %d", cols, rows);
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

   /***************************************************************************
   * Allocate memory to store the image then read the image from the file.
   ***************************************************************************/
   if(((*image) = (unsigned char *) malloc((*rows)*(*cols))) == NULL){
      fprintf(stderr, "Memory allocation failure in read_pgm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }
   if((*rows) != fread((*image), (*cols), (*rows), fp)){
      fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
      if(fp != stdin) fclose(fp);
      free((*image));
      return(0);
   }

   if(fp != stdin) fclose(fp);
   return(1);
}

/*******************************************************************************
* Function: unliken_excessive_basins
* Purpose: This function relabels the basins such that any basin that has pixels
* outside the breast is relabeled with a very large label. All basins that
* do not exceed the breast boundary are relabeled to fill in the missing labels.
* This is a preet adhoc routine.
* Name: Michael Heath, University of South Florida
* Date: 11/5/99
*******************************************************************************/
int unliken_excessive_basins(unsigned short int *wsimage, int rows, int cols, unsigned char *breastmaskim)
{
   unsigned short int max;
   int r, c, p, counter;
   unsigned short int *basin=NULL;

   /****************************************************************************
   * Find the maximum basin number.
   ****************************************************************************/
   max = wsimage[0];
   for(p=0;p<(rows*cols);p++) if(wsimage[p] > max) max = wsimage[p];

   basin = (unsigned short int *) calloc(max+1, sizeof(unsigned short int));
   for(r=0;r<=max;r++) basin[r] = r;

   for(r=0,p=0;r<rows;r++){
      for(c=0;c<cols;c++,p++){
         if(breastmaskim[p] == 0) basin[wsimage[p]] = 0;
      }
   }

   for(r=0;r<rows;r++){
      basin[wsimage[r*cols]] = 0;
      basin[wsimage[r*cols+cols-1]] = 0;
   }
   for(c=0;c<cols;c++){
      basin[wsimage[c]] = 0;
      basin[wsimage[(rows-1)*cols+c]] = 0;
   }

   counter = 0;
   for(r=1;r<=max;r++){
      if(basin[r] != 0){
         counter++;
         basin[r] = counter;
      }
   }
   for(r=1;r<=max;r++){
      if(basin[r] == 0){
         counter++;
         basin[r] = counter;
      }
   }

   for(r=0,p=0;r<rows;r++){
      for(c=0;c<cols;c++,p++){
         wsimage[p] = basin[wsimage[p]];
      }
   }

   free(basin);
}

/*******************************************************************************
* Function: find_minim_basin_value
* Purpose: This function computes the centroid (average row and average column)
* for each catchment basin (i.e. segmented region). It then finds the lowest
* numbered catchment basin that has its centroid in the ground truth region.
* Name: Michael Heath, University of South Florida
* Date: 11/10/99
*******************************************************************************/
int find_minim_basin_value(unsigned short int *wsimage, int rows, int cols,
    unsigned char *gtimage, unsigned long int **xplace, unsigned long int **yplace)
{
   int max, label, minlabel, rcoord, ccoord, p, r, c;
   unsigned long int *count=NULL, *sumx=NULL, *sumy=NULL;
   
   /****************************************************************************
   * Find the maximum basin number.
   ****************************************************************************/
   max = wsimage[0];
   for(p=0;p<(rows*cols);p++) if(wsimage[p] > max) max = wsimage[p];

   /* printf("There are %d bins.\n", max); */

   /****************************************************************************
   * Compute the average row and the average column of each region.
   ****************************************************************************/
   count = (unsigned long int *) calloc(max+1, sizeof(unsigned long int));
   sumx = (unsigned long int *) calloc(max+1, sizeof(unsigned long int));
   sumy = (unsigned long int *) calloc(max+1, sizeof(unsigned long int));

   for(r=0,p=0;r<rows;r++){
      for(c=0;c<cols;c++,p++){
         label = wsimage[p];
         count[label]++;
         sumx[label] += c;
         sumy[label] += r;
      }
   }

   if(gtimage == NULL){
      for(r=1;r<=max;r++){
         if(count[r] == 0) continue;
         sumy[r] = sumy[r] / count[r];
         sumx[r] = sumx[r] / count[r];
      }

      *xplace = sumx;
      *yplace = sumy;

      return(0);
   }

   /****************************************************************************
   * Find the minimum numbered region whose centroid falls in the ground truth
   * region.
   ****************************************************************************/
   minlabel = max;
   for(r=1;r<=max;r++){
      if(count[r] == 0) continue;
      rcoord = sumy[r] / count[r];
      ccoord = sumx[r] / count[r];
      /* printf("Trying %d %d\n", rcoord, ccoord); */
      if(gtimage[rcoord*cols+ccoord] != 0){
         minlabel = r;
         break;
      }
   }

   for(r=1;r<=max;r++){
      if(count[r] == 0) continue;
      sumy[r] = sumy[r] / count[r];
      sumx[r] = sumx[r] / count[r];
   }

   free(count);
/*
   free(sumx);
   free(sumy);
*/

   *xplace = sumx;
   *yplace = sumy;

   if(minlabel == max) return(0);
   return(minlabel);
}


/*******************************************************************************
* Function: find_min_value_per_basin
* Purpose: This function finds the minimum value of the feature image in each
* catchment basin.
* Name: Michael Heath, University of South Florida
* Date: 11/10/99
*******************************************************************************/
void find_min_value_per_basin(unsigned short int *wshedim, unsigned char *image,
    int rows, int cols, unsigned char **min_feature, int *max_basin,
    int **xcoord, int **ycoord)
{
   int max, label, p, r, c;
   int *xsum=NULL, *ysum=NULL, *count=NULL;

   /****************************************************************************
   * Find the maximum basin number.
   ****************************************************************************/
   max = wshedim[0];
   for(p=0;p<(rows*cols);p++) if(wshedim[p] > max) max = wshedim[p];
   *max_basin = max;

   /****************************************************************************
   * Compute the minumim in each nonzero basin.
   ****************************************************************************/
   (*min_feature) = (unsigned char *) calloc(max+1, sizeof(unsigned char));
   memset((*min_feature), 255, max+1);

   for(r=0,p=0;r<rows;r++){
      for(c=0;c<cols;c++,p++){
         label = wshedim[p];
         if(image[p] < (*min_feature)[label]) (*min_feature)[label] = image[p];
      }
   }

   /****************************************************************************
   * Find the position of the minimum in each basin.
   ****************************************************************************/
   xsum = (int *) calloc(max+1, sizeof(int));
   ysum = (int *) calloc(max+1, sizeof(int));
   count = (int *) calloc(max+1, sizeof(int));

   for(r=0,p=0;r<rows;r++){
      for(c=0;c<cols;c++,p++){
         label = wshedim[p];
         if(image[p] == (*min_feature)[label]){
            ysum[label] += r;
            xsum[label] += c;
            count[label]++;
         }
      }
   }

   for(r=1;r<=max;r++){
      if(count[r] != 0){
         xsum[r] /= count[r];
         ysum[r] /= count[r];
      }
   }
   
   *xcoord = xsum;
   *ycoord = ysum;

   free(count);
}


/*******************************************************************************
* Function: reorder_close_basins
* Purpose: We do not want multiple detections to lie too close to one another.
* To avoid this problem, we will reorder the bins to decrease the suspicion of
* bins that lie to close to other bins that will be declared as detections.
* Name: Michael Heath, University of South Florida
* Date: 12/13/99
*******************************************************************************/
void reorder_close_basins(unsigned short int *wshed, int rows, int cols,
   int *xcoord, int *ycoord, int pixdist)
{
   int r, c, pos, pixdist_sq, dist_sq;
   unsigned short int max;
   unsigned short int *lut=NULL;
   int *usebasin=NULL;

   /****************************************************************************
   * Find the maxim catchment basin label.
   ****************************************************************************/
   max = wshed[0];
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         if(wshed[pos] > max) max = wshed[pos];
      }
   }

   /****************************************************************************
   * Start with the lowest basin and remove any basins that are too close to it.
   * Then procede with basins in an increasing order with the same process.
   ****************************************************************************/
   usebasin = (int *) calloc(max+1, sizeof(int));
   for(pos=1;pos<=max;pos++) usebasin[pos] = 1;

   pixdist_sq = pixdist * pixdist;
   for(r=1;r<max;r++){
      if(usebasin[r] == 0) continue;
      for(c=(r+1);c<max;c++){
         if(usebasin[c] == 0) continue;
         dist_sq = (xcoord[c] - xcoord[r]) * (xcoord[c] - xcoord[r]) +
                   (ycoord[c] - ycoord[r]) * (ycoord[c] - ycoord[r]);
         if(dist_sq < pixdist_sq) usebasin[c] = 0;
      }
   }

   /****************************************************************************
   * Create a look-up-table to relabel the catchment basins.
   ****************************************************************************/
   lut = (unsigned short int *) calloc(max+1, sizeof(unsigned short int));
   c = 1;
   for(r=1;r<=max;r++){
      if(usebasin[r] == 1){
         lut[r] = c;
         c++;
      }
   }
   for(r=1;r<=max;r++){
      if(usebasin[r] == 0){
         lut[r] = c;
         c++;
      }
   }

   /****************************************************************************
   * Run the wshed image through the look-up-table.
   ****************************************************************************/
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         wshed[pos] = lut[wshed[pos]];
      }
   }

   free(usebasin);
   free(lut);
}
