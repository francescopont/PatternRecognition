/*******************************************************************************
* Program: transpose_MIASgt.c
* Purpose: To transpose the ground truth data for the MIAS images. This is
* done because the images were transposed and or flipped to get them into a
* an orientation that is similar to the orientation used in the DDSM databse.
* Since the images were transformed in this way, the ground truth is transformed
* in this way to (by this program). The program transpose_MIAS.c was used to
* transform the image files.
* If the image is of the left breast, the image is flipped "upside-down" while
* it is being transposed.
* Name: Michael D. Heath, University of South Florida
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
* Date: 3/19/2000
*
MIAS_TRUTH.txt
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
*  `x' (extra large)  5200pixels x 4000pixels   ERROR!!!!!
*  `x' (extra large)  4000pixels x 5200pixels

mdb001lm  G  CIRC  B  1815  1116  790
mdb002rl  G  CIRC  B  3091  1262  277

2nd column: Character of background tissue; 
                F - Fatty 
                G - Fatty-glandular
                D - Dense-glandular

3rd column: Class of abnormality present;
                CALC - Calcification
                CIRC - Well-defined/circumscribed masses
                SPIC - Spiculated masses
                MISC - Other, ill-defined masses
                ARCH - Architectural distortion
                ASYM - Asymmetry
                NORM - Normal

4th column: Severity of abnormality;
                B - Benign
                M - Malignant

*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LEFT 1
#define RIGHT 2

int main(int argc, char *argv[])
{
   FILE *fpin=NULL, *fpout=NULL;
   char output_filename[200];
   char image_filename[200], bgtissue[10], abnormality_class[10], severity[10];
   int rows, cols, r, c, ab_row, ab_col, ab_rad, pos, npos;
   int laterality=0;
   char fileline[200];
   int num_images = 0;
   int ab_col_trans, ab_row_trans;

   /****************************************************************************
   * Read in the input image.
   ****************************************************************************/
   if((fpin = fopen("MIAS_TRUTH.txt", "r")) == NULL){
      fprintf(stderr, "Error opening the file %s for reading!\n\n", "MIAS_TRUTH.txt");
      exit(1);
   }

   while((fgets(fileline, 200, fpin) != NULL) && (!feof(fpin))){

      if(strstr(fileline, "mdb") != NULL) num_images++;

      if(num_images == 0) continue;

      if(strstr(fileline, "mdb") != NULL){
         memset(image_filename, 0, 200);
         memset(bgtissue, 0, 10);
         pos = 0;
         sscanf(fileline, "%s%n", image_filename, &pos);
         sscanf(fileline+pos, "%s%n", bgtissue, &npos);
         pos += npos;
      }
      else if(fileline[0] == '\n'){
         pos = 0;
         continue;
      }

      memset(abnormality_class, 0, 10);
      memset(severity, 0, 10);
      ab_row = -1;
      ab_col = -1;
      ab_rad = -1;

      if(!((fileline[0] == ' ') || (fileline[0] == '\t'))){
         sscanf(fileline+pos, "%s%n", abnormality_class, &npos);
         pos += npos;
      }
      else{
         sscanf(fileline, "%s%n", abnormality_class, &pos);
      }

      if(strstr(abnormality_class, "NORM") == NULL){
         sscanf(fileline+pos, "%s%n", severity, &npos);
         pos += npos;

         if(strstr(fileline+pos, "*** see note 2") == NULL){
            sscanf(fileline+pos, "%d %d %d", &ab_col, &ab_row, &ab_rad);
         }
      }

      /*************************************************************************
      * Detemine the size of the image file (rows,cols) by looking at the last
      * character in the filename and using it to determine the size.
      *************************************************************************/
      switch(image_filename[strlen(image_filename)-1]){
         case 's': rows = 1600; cols = 4320; break;
         case 'm': rows = 2048; cols = 4320; break;
         case 'l': rows = 2600; cols = 4320; break;
         case 'x': rows = 4000; cols = 5200; break;
         default: fprintf(stderr, "Can not determine the size of the image for %c\n",
                     image_filename[strlen(image_filename)-1]);
      }

      /*************************************************************************
      * Detemine laterality of the breast from the filename.
      *************************************************************************/
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

      /*************************************************************************
      * Write the ground truth data out to a file. Append it to the file because
      * there may be more than one abnormality in a file.
      *************************************************************************/
      sprintf(output_filename, "%s.gt", image_filename);

      if((fpout = fopen(output_filename, "a")) == NULL){
         fprintf(stderr, "Error opening the file %s for writing.\n", output_filename);
         exit(1);
      }

      fprintf(fpout, "# %s", fileline);

      if(strstr(abnormality_class, "NORM") == NULL){
         if((ab_row != -1) && (ab_col != -1) && (ab_rad != -1)){
/*
            if(laterality == RIGHT){
               ab_col_trans = ab_row;
               ab_row_trans = ab_col;
            }
            else if(laterality == LEFT){
               ab_col_trans = rows-1-ab_row;
               ab_row_trans = cols-1-ab_col;
            }
*/
            if(laterality == RIGHT){
               ab_col_trans = rows-1-ab_row;
               ab_row_trans = ab_col;
            }
            else if(laterality == LEFT){
               ab_col_trans = ab_row;
               ab_row_trans = cols-1-ab_col;
            }
            fprintf(fpout, "%s %4d %4d %s %s %s %4d %4d %4d\n", image_filename, rows, cols,
               bgtissue, abnormality_class, severity, ab_col_trans, ab_row_trans, ab_rad);
         }
      }

      fclose(fpout);

/*
      if(strstr(abnormality_class, "NORM") != NULL){
         printf("[%s %s %s  ] %s", image_filename, bgtissue, abnormality_class, fileline);
      }
      else{
         if((ab_row != -1) && (ab_col != -1) && (ab_rad != -1)){
            if(laterality == RIGHT){
               ab_col_trans = ab_row;
               ab_row_trans = ab_col;
            }
            else if(laterality == LEFT){
               ab_col_trans = rows-1-ab_row;
               ab_row_trans = cols-1-ab_col;
            }
            printf("[%s %s %s %s %4d (%4d) %4d (%4d) %4d] %s", image_filename, bgtissue, abnormality_class, severity,
               ab_col, ab_col_trans, ab_row, ab_row_trans, ab_rad, fileline);
         }
         else{
            printf("[%s %s %s %s               ] %s", image_filename, bgtissue, abnormality_class, severity,
               fileline);
         }
      }
*/
   }

   fclose(fpin);
}
