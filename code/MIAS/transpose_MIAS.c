/*******************************************************************************
* Program: transpose_MIAS.c
* Purpose: To transpose an image file from the MIAS database. If the image is
* of the left breast, the image is flipped "upside-down" while it is being
* transposed.
* Name: Michael D. Heath, University of South Florida
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
* Date: 3/18/2000
*
*  Each entry has the prefix `mdb' followed by a serial number. The next
*  character, `l' or 'r', denotes left or right mammogram respectively from
*  the same examination. The final character, 's', 'm', 'l' or 'x' is the
*  size of the datafile:
*  
*                        Rows       Columns
*                        ----       -------
*  `s' (small)        1600pixels x 4320pixels
*  `m' (medium)       2048pixels x 4320pixels
*  `l' (large)        2600pixels x 4320pixels
*  `x' (extra large)  5200pixels x 4000pixels   ERROR!!!!!!
*  `x' (extra large)  4000pixels x 5200pixels
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#define PGM 0
#define LEFT 1
#define RIGHT 2

int main(int argc, char *argv[])
{
   FILE *fpin=NULL, *fpout=NULL;
   char *input_filename=NULL;
   char image_filename[200];
   char output_filename[200];
   int must_decompress = 0;
   int rows, cols, r, c;
   char systemline[200];
   unsigned char *rawimage=NULL;
   unsigned char *outline=NULL;
   int laterality=0;

   /****************************************************************************
   * Get the command line arguments.
   ****************************************************************************/
   if(argc != 2){
      fprintf(stderr, "\n\n<USAGE> %s MIASimage.Z\n\n", argv[0]);
      exit(1);
   }
   input_filename = argv[1];

   /****************************************************************************
   * Decompress the image.
   ****************************************************************************/
   if(input_filename[strlen(input_filename)-1] == 'Z'){
      memset(output_filename, 0, 200);
      strcpy(image_filename, input_filename);
      image_filename[strlen(input_filename)-2] = '\0';
      must_decompress = 1;
   }
   else{
      strcpy(image_filename, input_filename);
      must_decompress = 0;
   }
   if(must_decompress == 1){
      sprintf(systemline, "uncompress %s", input_filename);
      system(systemline);
   }

   /****************************************************************************
   * Detemine the size of the image file (rows,cols) by looking at the last
   * character in the filename and using it to determine the size.
   ****************************************************************************/
   switch(image_filename[strlen(image_filename)-1]){
      case 's': rows = 1600; cols = 4320; break;
      case 'm': rows = 2048; cols = 4320; break;
      case 'l': rows = 2600; cols = 4320; break;
      case 'x': rows = 4000; cols = 5200; break;
      default: fprintf(stderr, "Can not determine the size of the image for %c\n",
                  image_filename[strlen(image_filename)-1]);
   }

   /****************************************************************************
   * Detemine laterality of the breast from the filename.
   ****************************************************************************/
   switch(image_filename[strlen(image_filename)-2]){
      case 'r': laterality = RIGHT; break;
      case 'l': laterality = LEFT; break;
      default: fprintf(stderr, "Can not determine the laterality for %c\n",
                  image_filename[strlen(image_filename)-1]);
   }

/*
   printf("INPUT_FILENAME = %s\n", input_filename);
   printf("IMAGE_FILENAME = %s\n", image_filename);
   printf("MUST_DECOMPRESS = %d\n", must_decompress);
   printf("OUTPUT_FILENAME = %s\n", output_filename);
   printf("ROWS = %d   COLS = %d\n", rows, cols);
*/

   /****************************************************************************
   * Form the name of the output filename.
   ****************************************************************************/

   if(laterality == LEFT){
      if(PGM == 1){
         strcpy(output_filename, image_filename);
         strcat(output_filename, "_LMLO.pgm");
      }
      else{
         sprintf(output_filename, "%s_%d_%d_LMLO.raw", image_filename,
            rows, cols);
      }
   }
   else if(laterality == RIGHT){
      if(PGM == 1){
         strcpy(output_filename, image_filename);
         strcat(output_filename, "_RMLO.pgm");
      }
      else{
         sprintf(output_filename, "%s_%d_%d_RMLO.raw", image_filename,
            rows, cols);
      }
   }

   /****************************************************************************
   * Allocate memory for the entire image. Also allocate memory for one row
   * of the output image. This is an array with the length of the number of
   * rows in the input image because we are transposing it.
   ****************************************************************************/
   if((rawimage = (unsigned char *) calloc(rows*cols, sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Error allocating memory for the entire input image.\n");
      exit(1);
   }
   if((outline = (unsigned char *) calloc(rows, sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Error allocating memory for a row of the output image.\n");
      exit(1);
   }

   /****************************************************************************
   * Read in the input image.
   ****************************************************************************/
   if((fpin = fopen(image_filename, "rb")) == NULL){
      fprintf(stderr, "Error opening the file %s for reading!\n\n", image_filename);
      exit(1);
   }
   if(fread(rawimage, cols, rows, fpin) != rows){
      fprintf(stderr, "Error reading the input image.\n\n");
      exit(1);
   }
   fclose(fpin);

   /****************************************************************************
   * Write the output image to a file one line at a time.
   ****************************************************************************/
   if((fpout = fopen(output_filename, "wb")) == NULL){
      fprintf(stderr, "Error opening the file %s for writing!\n\n", output_filename);
      exit(1);
   }
   if(PGM == 1){
      fprintf(fpout, "P5\n%d %d\n255\n", rows, cols); /* transposing, remember? */
   }
   if(laterality == RIGHT){
      for(c=0;c<cols;c++){
         for(r=0;r<rows;r++){
            outline[r] = rawimage[r*cols+c];
         }
         fwrite(outline, 1, rows, fpout);
      }
   }
   else if(laterality == LEFT){
      for(c=(cols-1);c>=0;c--){
         for(r=0;r<rows;r++){
            outline[r] = rawimage[(rows-1-r)*cols+c];
         }
         fwrite(outline, 1, rows, fpout);
      }
   }

   fclose(fpout);

   /****************************************************************************
   * If we decompressed the image, recompress it.
   ****************************************************************************/
   if(must_decompress == 1){
      sprintf(systemline, "compress %s", image_filename);
      system(systemline);
   }

}
