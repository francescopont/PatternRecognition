/*******************************************************************************
* PROGRAM: ICSio.c
* PURPOSE: To handle Input and output for ICS image files.
* NAME: Michael Heath, University of South Florida
* DATE: 9/12/97
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/

#include "ICSio.h"

extern char created_files[4][200];

#define DEBUGMODE 0

/*******************************************************************************
* Procedure: read_ics_file
* Purpose: This procedure reads in the data from an ics file. This type of file
* contains information about the scanned mammographic images for a patient.
* Name: Michael Heath, University of South Florida
* Date: 9/9/97
* Note: I just wrote this function by looking at the contents of a few files.
* Here is a standard for this type of file and when I find documentation on the
* actual format, I will update this function. - Michael Heath
*
* ics_version 1.0
* filename B-9018-1
* DATE_OF_STUDY 00 00 0000
* PATIENT_AGE 00
* FILM 
* FILM_TYPE REGULAR
* DENSITY 0
* DATE_DIGITIZED 00 00 0000 
* DIGITIZER LUMISYS LASER_DENSOTOMETER
* SELECTED
* LEFT_CC LINES 2300 PIXELS_PER_LINE 1544 BITS_PER_PIXEL 12 RESOLUTION 100 OVERLAY
* LEFT_MLO LINES 2240 PIXELS_PER_LINE 1536 BITS_PER_PIXEL 12 RESOLUTION 100 OVERLAY
* RIGHT_CC LINES 2237 PIXELS_PER_LINE 1545 BITS_PER_PIXEL 12 RESOLUTION 100 OVERLAY
* RIGHT_MLO LINES 2225 PIXELS_PER_LINE 1547 BITS_PER_PIXEL 12 RESOLUTION 100 OVERLAY
*******************************************************************************/
int read_ics_file(char *filename, ICSDATA *icsdata, char *a_path)
{
   FILE *fp=NULL;
   char aline[100];
   int correct_ics_version_number = 0;
   char prefix;
   int case_number=0;
   char systemline[300];

   /**************************************************************************** 
   * Initialize the data structure for the ics data to all zeros.
   ****************************************************************************/
   memset(icsdata, 0, sizeof(ICSDATA));

   if((a_path != NULL) && (strcmp(a_path, "./") != 0)){
      sprintf(systemline, "ln -s %s%s", a_path, filename);
      system(systemline);
   }

   /**************************************************************************** 
   * Open the ics file for reading. Return if the file could not be opened.
   ****************************************************************************/
   if((fp = fopen(filename, "rb")) == NULL){
      fprintf(stderr, "Error opening the ics file named %s.\n", filename);
      return(0);
   }

   /****************************************************************************
   * Now extract the data from the file.
   ****************************************************************************/
   while(!feof(fp)){
      memset(aline, 0, 100);
      fgets(aline, 99, fp);

      /*************************************************************************
      * If we are in DEBUGMODE then print out the line of data.
      *************************************************************************/
      if(DEBUGMODE) printf("# %s", aline);

      if(strncmp(aline,"ics_version 1.0", strlen("ics_version 1.0")) == 0){
         correct_ics_version_number = 1;
         sscanf(aline, "%*s%f", &(icsdata->ics_version));
      }
      else if(strncmp(aline,"filename", strlen("filename")) == 0){
         sscanf(aline, "%*s%s", icsdata->filename);
         sscanf(icsdata->filename, "%c%*c%d", &prefix, &case_number);
      }
      else if(strncmp(aline,"DATE_OF_STUDY", strlen("DATE_OF_STUDY")) == 0){
         sscanf(aline, "%*s%d%d%d", &(icsdata->date_of_study.month),
            &(icsdata->date_of_study.day), &(icsdata->date_of_study.year));
      }
      else if(strncmp(aline,"PATIENT_AGE", strlen("PATIENT_AGE")) == 0){
         sscanf(aline, "%*s%d", &(icsdata->patient_age));
      }
      else if(strncmp(aline,"FILM ", strlen("FILM ")) == 0){
         sscanf(aline, "%s", icsdata->media);
      }
      else if(strncmp(aline,"FILM_TYPE", strlen("FILM_TYPE")) == 0){
         sscanf(aline, "%*s%s", icsdata->film_type);
      }
      else if(strncmp(aline,"DENSITY", strlen("DENSITY")) == 0){
         sscanf(aline, "%*s%d", &(icsdata->density));
      }
      else if(strncmp(aline,"DATE_DIGITIZED", strlen("DATE_DIGITIZED")) == 0){
         sscanf(aline, "%*s%d%d%d", &(icsdata->date_digitized.month),
            &(icsdata->date_digitized.day), &(icsdata->date_digitized.year));
      }
      else if(strncmp(aline,"DIGITIZER", strlen("DIGITIZER")) == 0){
         sscanf(aline, "%*s%s%s", icsdata->digitizer_brand, icsdata->digitizer_model);
      }
      else if(strncmp(aline,"SELECTED", strlen("SELECTED")) == 0){
         sscanf(aline, "%*s%s", icsdata->unknown);
      }

      /* LEFT CC */
      else if(strncmp(aline,"LEFT_CC", strlen("LEFT_CC")) == 0){
         if(strstr(aline, "LINES") != NULL){
            sscanf(strstr(aline, "LINES"), "%*s%d", &(icsdata->left_cc.rows));
         }
         if(strstr(aline, "PIXELS_PER_LINE") != NULL){
            sscanf(strstr(aline, "PIXELS_PER_LINE"), "%*s%d", &(icsdata->left_cc.cols));
         }
         if(strstr(aline, "BITS_PER_PIXEL") != NULL){
            sscanf(strstr(aline, "BITS_PER_PIXEL"), "%*s%d", &(icsdata->left_cc.bitsperpixel));
         }
         if(strstr(aline, "RESOLUTION") != NULL){
            sscanf(strstr(aline, "RESOLUTION"), "%*s%f", &(icsdata->left_cc.resolution));
         }
         if(strstr(aline, " OVERLAY") != NULL){
            icsdata->left_cc.overlay_exists = 1;
            sprintf(icsdata->left_cc.overlay_filename, "%c_%04d_1.LEFT_CC.OVERLAY", prefix, case_number);
            if((a_path != NULL) && (strcmp(a_path, "./") != 0)){
               sprintf(systemline, "ln -s %s%s", a_path, icsdata->left_cc.overlay_filename);
               system(systemline);
            }
         }
         if((icsdata->left_cc.rows==0) || (icsdata->left_cc.cols==0) ||
            (icsdata->left_cc.bitsperpixel==0) || (icsdata->left_cc.resolution==0)){
            fprintf(stderr, "Error in the data specification for image LEFT_CC. Assuming no image is available.\n");
            icsdata->left_cc.image_exists = 0;
         }
         else{
            icsdata->left_cc.image_exists = 1;
            sprintf(icsdata->left_cc.compressed_filename, "%c_%04d_1.LEFT_CC.LJPEG", prefix, case_number);
            if((a_path != NULL) && (strcmp(a_path, "./") != 0)){
               sprintf(systemline, "ln -s %s%s", a_path, icsdata->left_cc.compressed_filename);
               system(systemline);
            }
            sprintf(icsdata->left_cc.uncompressed_filename, "%c_%04d_1.LEFT_CC.LJPEG.1", prefix, case_number);
         }
      }

      /* LEFT MLO */
      else if(strncmp(aline,"LEFT_MLO", strlen("LEFT_MLO")) == 0){
         if(strstr(aline, "LINES") != NULL){
            sscanf(strstr(aline, "LINES"), "%*s%d", &(icsdata->left_mlo.rows));
         }
         if(strstr(aline, "PIXELS_PER_LINE") != NULL){
            sscanf(strstr(aline, "PIXELS_PER_LINE"), "%*s%d", &(icsdata->left_mlo.cols));
         }
         if(strstr(aline, "BITS_PER_PIXEL") != NULL){
            sscanf(strstr(aline, "BITS_PER_PIXEL"), "%*s%d", &(icsdata->left_mlo.bitsperpixel));
         }
         if(strstr(aline, "RESOLUTION") != NULL){
            sscanf(strstr(aline, "RESOLUTION"), "%*s%f", &(icsdata->left_mlo.resolution));
         }
         if(strstr(aline, " OVERLAY") != NULL){
            icsdata->left_mlo.overlay_exists = 1;
            sprintf(icsdata->left_mlo.overlay_filename, "%c_%04d_1.LEFT_MLO.OVERLAY", prefix, case_number);
            if((a_path != NULL) && (strcmp(a_path, "./") != 0)){
               sprintf(systemline, "ln -s %s%s", a_path, icsdata->left_mlo.overlay_filename);
               system(systemline);
            }
         }
         if((icsdata->left_mlo.rows==0) || (icsdata->left_mlo.cols==0) ||
            (icsdata->left_mlo.bitsperpixel==0) || (icsdata->left_mlo.resolution==0)){
            fprintf(stderr, "Error in the data specification for image LEFT_MLO. Assuming no image is available.\n");
            icsdata->left_mlo.image_exists = 0;
         }
         else{
            icsdata->left_mlo.image_exists = 1;
            sprintf(icsdata->left_mlo.compressed_filename, "%c_%04d_1.LEFT_MLO.LJPEG", prefix, case_number);
            if((a_path != NULL) && (strcmp(a_path, "./") != 0)){
               sprintf(systemline, "ln -s %s%s", a_path, icsdata->left_mlo.compressed_filename);
               system(systemline);
            }
            sprintf(icsdata->left_mlo.uncompressed_filename, "%c_%04d_1.LEFT_MLO.LJPEG.1", prefix, case_number);
         }
      }

      /* RIGHT CC */
      else if(strncmp(aline,"RIGHT_CC", strlen("RIGHT_CC")) == 0){
         if(strstr(aline, "LINES") != NULL){
            sscanf(strstr(aline, "LINES"), "%*s%d", &(icsdata->right_cc.rows));
         }
         if(strstr(aline, "PIXELS_PER_LINE") != NULL){
            sscanf(strstr(aline, "PIXELS_PER_LINE"), "%*s%d", &(icsdata->right_cc.cols));
         }
         if(strstr(aline, "BITS_PER_PIXEL") != NULL){
            sscanf(strstr(aline, "BITS_PER_PIXEL"), "%*s%d", &(icsdata->right_cc.bitsperpixel));
         }
         if(strstr(aline, "RESOLUTION") != NULL){
            sscanf(strstr(aline, "RESOLUTION"), "%*s%f", &(icsdata->right_cc.resolution));
         }
         if(strstr(aline, " OVERLAY") != NULL){
            icsdata->right_cc.overlay_exists = 1;
            sprintf(icsdata->right_cc.overlay_filename, "%c_%04d_1.RIGHT_CC.OVERLAY", prefix, case_number);
            if((a_path != NULL) && (strcmp(a_path, "./") != 0)){
               sprintf(systemline, "ln -s %s%s", a_path, icsdata->right_cc.overlay_filename);
               system(systemline);
            }
         }
         if((icsdata->right_cc.rows==0) || (icsdata->right_cc.cols==0) ||
            (icsdata->right_cc.bitsperpixel==0) || (icsdata->right_cc.resolution==0)){
            fprintf(stderr, "Error in the data specification for image RIGHT_CC. Assuming no image is available.\n");
            icsdata->right_cc.image_exists = 0;
         }
         else{
            icsdata->right_cc.image_exists = 1;
            sprintf(icsdata->right_cc.compressed_filename, "%c_%04d_1.RIGHT_CC.LJPEG", prefix, case_number);
            if((a_path != NULL) && (strcmp(a_path, "./") != 0)){
               sprintf(systemline, "ln -s %s%s", a_path, icsdata->right_cc.compressed_filename);
               system(systemline);
            }
            sprintf(icsdata->right_cc.uncompressed_filename, "%c_%04d_1.RIGHT_CC.LJPEG.1", prefix, case_number);
         }
      }

      /* RIGHT MLO */
      else if(strncmp(aline,"RIGHT_MLO", strlen("RIGHT_MLO")) == 0){
         if(strstr(aline, "LINES") != NULL){
            sscanf(strstr(aline, "LINES"), "%*s%d", &(icsdata->right_mlo.rows));
         }
         if(strstr(aline, "PIXELS_PER_LINE") != NULL){
            sscanf(strstr(aline, "PIXELS_PER_LINE"), "%*s%d", &(icsdata->right_mlo.cols));
         }
         if(strstr(aline, "BITS_PER_PIXEL") != NULL){
            sscanf(strstr(aline, "BITS_PER_PIXEL"), "%*s%d", &(icsdata->right_mlo.bitsperpixel));
         }
         if(strstr(aline, "RESOLUTION") != NULL){
            sscanf(strstr(aline, "RESOLUTION"), "%*s%f", &(icsdata->right_mlo.resolution));
         }
         if(strstr(aline, " OVERLAY") != NULL){
            icsdata->right_mlo.overlay_exists = 1;
            sprintf(icsdata->right_mlo.overlay_filename, "%c_%04d_1.RIGHT_MLO.OVERLAY", prefix, case_number);
            if((a_path != NULL) && (strcmp(a_path, "./") != 0)){
               sprintf(systemline, "ln -s %s%s", a_path, icsdata->right_mlo.overlay_filename);
               system(systemline);
            }
         }
         if((icsdata->right_mlo.rows==0) || (icsdata->right_mlo.cols==0) ||
            (icsdata->right_mlo.bitsperpixel==0) || (icsdata->right_mlo.resolution==0)){
            fprintf(stderr, "Error in the data specification for image RIGHT_MLO. Assuming no image is available.\n");
            icsdata->right_mlo.image_exists = 0;
         }
         else{
            icsdata->right_mlo.image_exists = 1;
            sprintf(icsdata->right_mlo.compressed_filename, "%c_%04d_1.RIGHT_MLO.LJPEG", prefix, case_number);
            if((a_path != NULL) && (strcmp(a_path, "./") != 0)){
               sprintf(systemline, "ln -s %s%s", a_path, icsdata->right_mlo.compressed_filename);
               system(systemline);
            }
            sprintf(icsdata->right_mlo.uncompressed_filename, "%c_%04d_1.RIGHT_MLO.LJPEG.1", prefix, case_number);
         }
      }

   }
   if(DEBUGMODE) printf("\n");

   /****************************************************************************
   * Make sure that the file is an ics file of the correct version.
   ****************************************************************************/
   if(correct_ics_version_number == 0){
      fprintf(stderr, "Error reading the \"ics\" file. File format incorrect.\n");
      fclose(fp);
      return(0);
   }

   /****************************************************************************
   * If we are in DEBUGMODE then print out the data from the datastructure.
   ****************************************************************************/
   if(DEBUGMODE){
      printf("ics_version %f\n", icsdata->ics_version);
      printf("filename %s\n", icsdata->filename);
      printf("DATE_OF_STUDY %d %d %d\n", icsdata->date_of_study.month,
         icsdata->date_of_study.day, icsdata->date_of_study.year);
      printf("PATIENT_AGE %d\n", icsdata->patient_age);
      printf("%s\n", icsdata->media);
      printf("FILM_TYPE %s\n", icsdata->film_type);
      printf("DENSITY %d\n", icsdata->density);
      printf("DATE_DIGITIZED %02d %02d %04d\n", icsdata->date_digitized.month,
         icsdata->date_digitized.day, icsdata->date_digitized.year);
      printf("DIGITIZER %s %s\n", icsdata->digitizer_brand, icsdata->digitizer_model);
      printf("%s\n", icsdata->unknown);
   
      if(icsdata->left_cc.image_exists){
         printf("LEFT_CC LINES %d PIXELS_PER_LINE %d BITS_PER_PIXEL %d RESOLUTION %f ",
	    icsdata->left_cc.rows, icsdata->left_cc.cols, icsdata->left_cc.bitsperpixel, icsdata->left_cc.resolution);
         if(icsdata->left_cc.overlay_exists) printf("OVERLAY\n");
         else printf("NO-OVERLAY\n");
         printf("      %s\n", icsdata->left_cc.compressed_filename);
         printf("      %s\n", icsdata->left_cc.uncompressed_filename);
         if(icsdata->left_cc.overlay_exists) printf("      %s\n", icsdata->left_cc.overlay_filename);
      }

      if(icsdata->left_mlo.image_exists){
         printf("LEFT_MLO LINES %d PIXELS_PER_LINE %d BITS_PER_PIXEL %d RESOLUTION %f ",
	    icsdata->left_mlo.rows, icsdata->left_mlo.cols, icsdata->left_mlo.bitsperpixel, icsdata->left_mlo.resolution);
         if(icsdata->left_mlo.overlay_exists) printf("OVERLAY\n");
         else printf("NO-OVERLAY\n");
         printf("      %s\n", icsdata->left_mlo.compressed_filename);
         printf("      %s\n", icsdata->left_mlo.uncompressed_filename);
         if(icsdata->left_mlo.overlay_exists) printf("      %s\n", icsdata->left_mlo.overlay_filename);
      }
   
      if(icsdata->right_cc.image_exists){
         printf("RIGHT_CC LINES %d PIXELS_PER_LINE %d BITS_PER_PIXEL %d RESOLUTION %f ",
	    icsdata->right_cc.rows, icsdata->right_cc.cols, icsdata->right_cc.bitsperpixel, icsdata->right_cc.resolution);
         if(icsdata->right_cc.overlay_exists) printf("OVERLAY\n");
         else printf("NO-OVERLAY\n");
         printf("      %s\n", icsdata->right_cc.compressed_filename);
         printf("      %s\n", icsdata->right_cc.uncompressed_filename);
         if(icsdata->right_cc.overlay_exists) printf("      %s\n", icsdata->right_cc.overlay_filename);
      }
   
      if(icsdata->right_mlo.image_exists){
         printf("RIGHT_MLO LINES %d PIXELS_PER_LINE %d BITS_PER_PIXEL %d RESOLUTION %f ",
	    icsdata->right_mlo.rows, icsdata->right_mlo.cols, icsdata->right_mlo.bitsperpixel, icsdata->right_mlo.resolution);
         if(icsdata->right_mlo.overlay_exists) printf("OVERLAY\n");
         else printf("NO-OVERLAY\n");
         printf("      %s\n", icsdata->right_mlo.compressed_filename);
         printf("      %s\n", icsdata->right_mlo.uncompressed_filename);
         if(icsdata->right_mlo.overlay_exists) printf("      %s\n", icsdata->right_mlo.overlay_filename);
      }
   }
   fclose(fp);
   return(1);
}

/*******************************************************************************
* Procedure: decompress_ics_images
* Purpose: This procedure decompresses the images that were specified in an
* ics file.
* Name: Michael Heath, University of South Florida
* Date: 9/9/97
*******************************************************************************/
int decompress_ics_images(ICSDATA *icsdata)
{
   FILE *fp;
   int bytesperpixel;
   int i;
   int number_of_created_files=0;

   int decompress_LJPEG_image(char *filename);

   /****************************************************************************
   * Initialize the list of created file to show no files.
   ****************************************************************************/
   for(i=0;i<4;i++) memset(created_files[i], 0, 200);

   /****************************************************************************
   * If there is a LEFT_CC image and it is not decompressed, then decompress it.
   * If the file is already exists in a decompressed form, and has the correct
   * file size, do not decomress it.
   ****************************************************************************/
   if(icsdata->left_cc.image_exists){

      if(icsdata->left_cc.bitsperpixel <= 8) bytesperpixel = 1;
      else bytesperpixel = 2;

      fp = NULL;
      if((fp = fopen(icsdata->left_cc.uncompressed_filename, "rb")) == NULL){
         if(decompress_LJPEG_image(icsdata->left_cc.compressed_filename) == 0) return(0);
         icsdata->left_cc.did_i_decompress_it = 1;
         sprintf(created_files[number_of_created_files], "%s", icsdata->left_cc.uncompressed_filename);
         number_of_created_files++;
      }
      else{
         fseek(fp, 0, 2);
         if(ftell(fp) != (bytesperpixel * icsdata->left_cc.rows * icsdata->left_cc.cols)){
            fclose(fp);
            if(decompress_LJPEG_image(icsdata->left_cc.compressed_filename) == 0) return(0);
            icsdata->left_cc.did_i_decompress_it = 1;
            sprintf(created_files[number_of_created_files], "%s", icsdata->left_cc.uncompressed_filename);
            number_of_created_files++;
         }
         else fclose(fp);
      }
   }

   /****************************************************************************
   * If there is a LEFT_MLO image and it is not decompressed, then decompress it.
   * If the file is already exists in a decompressed form, and has the correct
   * file size, do not decomress it.
   ****************************************************************************/
   if(icsdata->left_mlo.image_exists){

      if(icsdata->left_mlo.bitsperpixel <= 8) bytesperpixel = 1;
      else bytesperpixel = 2;

      fp = NULL;
      if((fp = fopen(icsdata->left_mlo.uncompressed_filename, "rb")) == NULL){
         if(decompress_LJPEG_image(icsdata->left_mlo.compressed_filename) == 0) return(0);
         icsdata->left_mlo.did_i_decompress_it = 1;
         sprintf(created_files[number_of_created_files], "%s", icsdata->left_mlo.uncompressed_filename);
         number_of_created_files++;
      }
      else{
         fseek(fp, 0, 2);
         if(ftell(fp) != (bytesperpixel * icsdata->left_mlo.rows * icsdata->left_mlo.cols)){
            fclose(fp);
            if(decompress_LJPEG_image(icsdata->left_mlo.compressed_filename) == 0) return(0);
            icsdata->left_mlo.did_i_decompress_it = 1;
            sprintf(created_files[number_of_created_files], "%s", icsdata->left_mlo.uncompressed_filename);
            number_of_created_files++;
         }
         else fclose(fp);
      }
   }

   /****************************************************************************
   * If there is a RIGHT_CC image and it is not decompressed, then decompress it.
   * If the file is already exists in a decompressed form, and has the correct
   * file size, do not decomress it.
   ****************************************************************************/
   if(icsdata->right_cc.image_exists){

      if(icsdata->right_cc.bitsperpixel <= 8) bytesperpixel = 1;
      else bytesperpixel = 2;

      fp = NULL;
      if((fp = fopen(icsdata->right_cc.uncompressed_filename, "rb")) == NULL){
         if(decompress_LJPEG_image(icsdata->right_cc.compressed_filename) == 0) return(0);
         icsdata->right_cc.did_i_decompress_it = 1;
         sprintf(created_files[number_of_created_files], "%s", icsdata->right_cc.uncompressed_filename);
         number_of_created_files++;
      }
      else{
         fseek(fp, 0, 2);
         if(ftell(fp) != (bytesperpixel * icsdata->right_cc.rows * icsdata->right_cc.cols)){
            fclose(fp);
            if(decompress_LJPEG_image(icsdata->right_cc.compressed_filename) == 0) return(0);
            icsdata->right_cc.did_i_decompress_it = 1;
            sprintf(created_files[number_of_created_files], "%s", icsdata->right_cc.uncompressed_filename);
            number_of_created_files++;
         }
         else fclose(fp);
      }
   }

   /****************************************************************************
   * If there is a RIGHT_MLO image and it is not decompressed, then decompress it.
   * If the file is already exists in a decompressed form, and has the correct
   * file size, do not decomress it.
   ****************************************************************************/
   if(icsdata->right_mlo.image_exists){

      if(icsdata->right_mlo.bitsperpixel <= 8) bytesperpixel = 1;
      else bytesperpixel = 2;

      fp = NULL;
      if((fp = fopen(icsdata->right_mlo.uncompressed_filename, "rb")) == NULL){
         if(decompress_LJPEG_image(icsdata->right_mlo.compressed_filename) == 0) return(0);
         icsdata->right_mlo.did_i_decompress_it = 1;
         sprintf(created_files[number_of_created_files], "%s", icsdata->right_mlo.uncompressed_filename);
         number_of_created_files++;
      }
      else{
         fseek(fp, 0, 2);
         if(ftell(fp) != (bytesperpixel * icsdata->right_mlo.rows * icsdata->right_mlo.cols)){
            fclose(fp);
            if(decompress_LJPEG_image(icsdata->right_mlo.compressed_filename) == 0) return(0);
            icsdata->right_mlo.did_i_decompress_it = 1;
            sprintf(created_files[number_of_created_files], "%s", icsdata->right_mlo.uncompressed_filename);
            number_of_created_files++;
         }
         else fclose(fp);
      }
   }

   return(1);
}

/*******************************************************************************
* Procedure: decompress_LJPEG_image
* Purpose: This procedure uses the system function to call the program that
* decompresses an image file that was compressed with a lossless JPEG
* compression program.
* Name: Michael Heath, University of South Florida
* Date: 9/9/97
*******************************************************************************/
int decompress_LJPEG_image(char *filename)
{
   FILE *fptest=NULL;
   char JPEG_decompress_line[200];
   char *jpeg_program = NULL;

   printf("Decompressing the file (%s)\n", filename);

   if((fptest = fopen(JPEG_EXECUTABLE, "r")) == NULL){
      jpeg_program = getenv("JPEG_PROGRAM");
      if(jpeg_program == NULL){
         fprintf(stderr, "\nError! You must define the environment variable for the jpeg program.\n");
         fprintf(stderr, "   (i.e. setenv JPEG_PROGRAM \"/home/captiva1/heath/bin/jpeg\")\n\n");
         exit(1);
      }
      if((fptest = fopen(jpeg_program, "r")) == NULL){
         fprintf(stderr, "\nError! Could not find the file %s\n\n", jpeg_program);
         exit(1);
      }
      else fclose(fptest);
   }
   else{
      jpeg_program = JPEG_EXECUTABLE;
      fclose(fptest);
   }

   sprintf(JPEG_decompress_line, "%s -d -s %s", jpeg_program, filename);
   if(system(JPEG_decompress_line) == (-1)){
      fprintf(stderr, "Error decompressing the file %s.\n", filename);
      return(0);
   }

   return(1);
}

/*******************************************************************************
* Procedure: decompress_ics_image
* Purpose: This procedure decompresses an image that was specified in an
* ics file.
* Name: Michael Heath, University of South Florida
* Date: 11/15/97
*******************************************************************************/
int decompress_ics_image(ICSDATA *icsdata, char *view)
{
   FILE *fp;
   int bytesperpixel;
   int i;
   int number_of_created_files=0;

   int decompress_LJPEG_image(char *filename);

   /* printf("decompress_ics_image called with view = %s\n", view); */

   /****************************************************************************
   * Initialize the list of created file to show no files.
   ****************************************************************************/
   for(i=0;i<4;i++) memset(created_files[i], 0, 200);

   /****************************************************************************
   * If there is a LEFT_CC image and it is not decompressed, then decompress it.
   * If the file is already exists in a decompressed form, and has the correct
   * file size, do not decomress it.
   ****************************************************************************/
   if((strcmp(view, "LEFT_CC") == 0) && (icsdata->left_cc.image_exists)){

      /* printf("%s %d\n", view, icsdata->left_cc.image_exists); */

      if(icsdata->left_cc.bitsperpixel <= 8) bytesperpixel = 1;
      else bytesperpixel = 2;

      fp = NULL;
      if((fp = fopen(icsdata->left_cc.uncompressed_filename, "rb")) == NULL){
         if(decompress_LJPEG_image(icsdata->left_cc.compressed_filename) == 0) return(0);
         icsdata->left_cc.did_i_decompress_it = 1;
         sprintf(created_files[number_of_created_files], "%s", icsdata->left_cc.uncompressed_filename);
         number_of_created_files++;
      }
      else{
         fseek(fp, 0, 2);
         if(ftell(fp) != (bytesperpixel * icsdata->left_cc.rows * icsdata->left_cc.cols)){
            fclose(fp);
            if(decompress_LJPEG_image(icsdata->left_cc.compressed_filename) == 0) return(0);
            icsdata->left_cc.did_i_decompress_it = 1;
            sprintf(created_files[number_of_created_files], "%s", icsdata->left_cc.uncompressed_filename);
            number_of_created_files++;
         }
         else fclose(fp);
      }
   }

   /****************************************************************************
   * If there is a LEFT_MLO image and it is not decompressed, then decompress it.
   * If the file is already exists in a decompressed form, and has the correct
   * file size, do not decomress it.
   ****************************************************************************/
   if((strcmp(view, "LEFT_MLO") == 0) && (icsdata->left_mlo.image_exists)){

      /* printf("%s %d\n", view, icsdata->left_mlo.image_exists); */

      if(icsdata->left_mlo.bitsperpixel <= 8) bytesperpixel = 1;
      else bytesperpixel = 2;

      fp = NULL;
      if((fp = fopen(icsdata->left_mlo.uncompressed_filename, "rb")) == NULL){
         if(decompress_LJPEG_image(icsdata->left_mlo.compressed_filename) == 0) return(0);
         icsdata->left_mlo.did_i_decompress_it = 1;
         sprintf(created_files[number_of_created_files], "%s", icsdata->left_mlo.uncompressed_filename);
         number_of_created_files++;
      }
      else{
         fseek(fp, 0, 2);
         if(ftell(fp) != (bytesperpixel * icsdata->left_mlo.rows * icsdata->left_mlo.cols)){
            fclose(fp);
            if(decompress_LJPEG_image(icsdata->left_mlo.compressed_filename) == 0) return(0);
            icsdata->left_mlo.did_i_decompress_it = 1;
            sprintf(created_files[number_of_created_files], "%s", icsdata->left_mlo.uncompressed_filename);
            number_of_created_files++;
         }
         else fclose(fp);
      }
   }

   /****************************************************************************
   * If there is a RIGHT_CC image and it is not decompressed, then decompress it.
   * If the file is already exists in a decompressed form, and has the correct
   * file size, do not decomress it.
   ****************************************************************************/
   if((strcmp(view, "RIGHT_CC") == 0) && (icsdata->right_cc.image_exists)){

      /* printf("%s %d\n", view, icsdata->right_cc.image_exists); */

      if(icsdata->right_cc.bitsperpixel <= 8) bytesperpixel = 1;
      else bytesperpixel = 2;

      fp = NULL;
      if((fp = fopen(icsdata->right_cc.uncompressed_filename, "rb")) == NULL){
         if(decompress_LJPEG_image(icsdata->right_cc.compressed_filename) == 0) return(0);
         icsdata->right_cc.did_i_decompress_it = 1;
         sprintf(created_files[number_of_created_files], "%s", icsdata->right_cc.uncompressed_filename);
         number_of_created_files++;
      }
      else{
         fseek(fp, 0, 2);
         if(ftell(fp) != (bytesperpixel * icsdata->right_cc.rows * icsdata->right_cc.cols)){
            fclose(fp);
            if(decompress_LJPEG_image(icsdata->right_cc.compressed_filename) == 0) return(0);
            icsdata->right_cc.did_i_decompress_it = 1;
            sprintf(created_files[number_of_created_files], "%s", icsdata->right_cc.uncompressed_filename);
            number_of_created_files++;
         }
         else fclose(fp);
      }
   }

   /****************************************************************************
   * If there is a RIGHT_MLO image and it is not decompressed, then decompress it.
   * If the file is already exists in a decompressed form, and has the correct
   * file size, do not decomress it.
   ****************************************************************************/
   if((strcmp(view, "RIGHT_MLO") == 0) && (icsdata->right_mlo.image_exists)){

      /* printf("%s %d\n", view, icsdata->right_mlo.image_exists); */

      if(icsdata->right_mlo.bitsperpixel <= 8) bytesperpixel = 1;
      else bytesperpixel = 2;

      fp = NULL;
      if((fp = fopen(icsdata->right_mlo.uncompressed_filename, "rb")) == NULL){
         if(decompress_LJPEG_image(icsdata->right_mlo.compressed_filename) == 0) return(0);
         icsdata->right_mlo.did_i_decompress_it = 1;
         sprintf(created_files[number_of_created_files], "%s", icsdata->right_mlo.uncompressed_filename);
         number_of_created_files++;
      }
      else{
         fseek(fp, 0, 2);
         if(ftell(fp) != (bytesperpixel * icsdata->right_mlo.rows * icsdata->right_mlo.cols)){
            fclose(fp);
            if(decompress_LJPEG_image(icsdata->right_mlo.compressed_filename) == 0) return(0);
            icsdata->right_mlo.did_i_decompress_it = 1;
            sprintf(created_files[number_of_created_files], "%s", icsdata->right_mlo.uncompressed_filename);
            number_of_created_files++;
         }
         else fclose(fp);
      }
   }

   return(1);
}

/*******************************************************************************
* Procedure: ru_exit
* Purpose: When we exit the program, we do not want to leave uncompressed image
* files laying around. This procedure removes the image files that we created.
* Name: Michael Heath, University of South Florida
* Date: 9/9/97
*******************************************************************************/
void ru_exit()
{
   int i;
   char system_line[200];

   for(i=0;i<4;i++){
      if(strlen(created_files[i]) != 0){
         sprintf(system_line, "\\rm %s", created_files[i]);
         system(system_line);
      }
   }

   exit(1);
}
