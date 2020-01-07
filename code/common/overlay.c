/*******************************************************************************
* File: overlay.c
* Purpose: This file contains code for reading the contents of an overlay file.
* Name: Michael Heath, University of South Florida
* Date: 5/6/98
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include "overlay.h"
#include <math.h>

/* #define TESTMODE   Comment out this line after testing the code. */

#define PRINTEVERYTHING 0

extern int VERBOSE;

int get_line_with_astring_in_it(FILE *fp, char *astring, char *line, char *filename);
int get_string_with_astring_in_it(FILE *fp, char *astring, char *line, char *filename);
int read_htmlthumbnail_file(char *filename);

/*******************************************************************************
* This is a segment of code that can be used to run the function that reads in
* the data from an overlay file.
*******************************************************************************/
#ifdef TESTMODE
void main(int argc, char *argv[])
{
   char *filename=NULL;
   OVERLAY_DATA overlay_data;
   int i;

   if(argc < 2){
      fprintf(stderr, "\n\n<USAGE> %s filename1.html [filename2.html ... filenameN.html]\n\n", argv[0]);
      exit(1);
   }

   for(i=1;i<argc;i++){
      filename = argv[i];
      if(read_htmlthumbnail_file(filename) == 0) exit(1);
   }

/*
   if(PRINTEVERYTHING) printf("Calling the function to read the file: %s\n", filename);
   if(read_overlay_file(argv[1], &overlay_data) == 0) exit(1);

   if(write_overlay_file(argv[2], overlay_data) == 0) exit(1);
*/

   printf("Successully finished processing the file %s.\n", argv[1]);
}
#endif

/*******************************************************************************
* Below is code for reading data from an OVERLAY file.
*******************************************************************************/

/*******************************************************************************
* Function: read_htmlthumbnail_file
* Purpose: This function reads in information from an HTML file and forms
* a line for each lesion listing relevent information about the case and the
* lesion.
*******************************************************************************/
int read_htmlthumbnail_file(char *filename)
{
   enum VOLTYPE {NORMALVOL=1, CANCERVOL=2, BENIGNVOL=4, BENIGNWCVOL=8};
   enum SCANNERTYPE {DBA=1, LUMISYS=2, HOWTEK=4};
   enum PATHOLOGYTYPE {MALIGNANT=1, BENIGN=2, BENIGN_WITHOUT_CALLBACK=4, UNPROVEN=8};

   FILE *fp=NULL;
   char line[200];
   int i, j, number_of_overlays, total_abnormalities, current_image;
   int density, left_cc_overlay=0, left_mlo_overlay=0, right_cc_overlay=0, right_mlo_overlay=0;
   char digitizer[20];
   int assessment, subtlety;
   char pathology[20];
   DESCRIPTION description1, description2;
   int num_descriptions;
   int voltype;
   int scannertype;
   int pathologytype;

   /****************************************************************************
   * Determine the type of the volume.
   ****************************************************************************/
   if(strstr(filename, "normal_") != NULL) voltype = NORMALVOL;
   else if(strstr(filename, "cancer_") != NULL) voltype = CANCERVOL;
   else if(strstr(filename, "benign_without_callback") != NULL) voltype = BENIGNWCVOL;
   else if(strstr(filename, "benign_") != NULL) voltype = BENIGNVOL;
   else{
      fprintf(stderr, "The volume type can not be determined from the filename %s.\n", filename);
      return(0);
   } 

   /****************************************************************************
   * Open the file.
   ****************************************************************************/
   if((fp = fopen(filename, "r")) == NULL){
      fprintf(stderr, "Error opening the HTML file %s for reading.\n", filename);
      return(0);
   }
   if(PRINTEVERYTHING) printf("Successfully opened the HTML file %s.\n", filename);

   /****************************************************************************
   * Make sure the file is an HTML file. Do this by looking for "HTML" on
   * some line of the file. It should be on the first line of the file because
   * the htmlthumbnail files were produced by a program that puts it there.
   ****************************************************************************/
   if(get_line_with_astring_in_it(fp, "HTML", line, filename) == 0){
      fclose(fp);
      return(0);
   }

   /****************************************************************************
   * Get the DENSITY, the DIGITIZER and which views have OVERLAYS.
   ics_version 1.0
   filename B-3027-1
   DATE_OF_STUDY 28 4 1995
   PATIENT_AGE 65
   FILM
   FILM_TYPE REGULAR
   DENSITY 2
   DATE_DIGITIZED 23 7 1997
   DIGITIZER LUMISYS LASER
   SEQUENCE
   LEFT_CC LINES 4720 PIXELS_PER_LINE 3128 BITS_PER_PIXEL 12 RESOLUTION 50 OVERLAY
   LEFT_MLO LINES 4696 PIXELS_PER_LINE 3184 BITS_PER_PIXEL 12 RESOLUTION 50 OVERLAY
   RIGHT_CC LINES 4680 PIXELS_PER_LINE 3112 BITS_PER_PIXEL 12 RESOLUTION 50 NON_OVERLAY
   RIGHT_MLO LINES 4720 PIXELS_PER_LINE 3120 BITS_PER_PIXEL 12 RESOLUTION 50 NON_OVERLAY
   ****************************************************************************/

   if(get_line_with_astring_in_it(fp, "DENSITY", line, filename) == 0){
      fclose(fp);
      return(0);
   }
   sscanf(line, "DENSITY %d", &density);

   if(get_line_with_astring_in_it(fp, "DIGITIZER", line, filename) == 0){
      fclose(fp);
      return(0);
   }

   sscanf(line, "DIGITIZER %s", digitizer);

   if(strstr(digitizer, "DBA") != NULL) scannertype = DBA;
   else if(strstr(digitizer, "LUMISYS") != NULL) scannertype = LUMISYS;
   else if(strstr(digitizer, "HOWTEK") != NULL) scannertype = HOWTEK;
   else{
      fprintf(stderr, "[[%s, %s]]", line, digitizer);
      fprintf(stderr, "The scanner type can not be determined from the html file %s.\n", filename);
      fclose(fp);
      return(0);
   } 

   if(get_line_with_astring_in_it(fp, "LEFT_CC", line, filename) == 0){
      fclose(fp);
      return(0);
   }
   if(strstr(line, "NON_OVERLAY") == NULL) left_cc_overlay = 1;

   if(get_line_with_astring_in_it(fp, "LEFT_MLO", line, filename) == 0){
      fclose(fp);
      return(0);
   }
   if(strstr(line, "NON_OVERLAY") == NULL) left_mlo_overlay = 1;

   if(get_line_with_astring_in_it(fp, "RIGHT_CC", line, filename) == 0){
      fclose(fp);
      return(0);
   }
   if(strstr(line, "NON_OVERLAY") == NULL) right_cc_overlay = 1;

   if(get_line_with_astring_in_it(fp, "RIGHT_MLO", line, filename) == 0){
      fclose(fp);
      return(0);
   }
   if(strstr(line, "NON_OVERLAY") == NULL) right_mlo_overlay = 1;

   /*
   if(left_cc_overlay) printf(" LEFT_CC");
   if(left_mlo_overlay) printf(" LEFT_MLO");
   if(right_cc_overlay) printf(" RIGHT_CC");
   if(right_mlo_overlay) printf(" RIGHT_MLO");
   printf("\n");
   */

   if(voltype == NORMALVOL){
      printf("%s", filename);
      printf(" %d", voltype);
      printf(" %d", (int)pow(2.0, (double)density-1));
      printf(" %d\n", scannertype);
      fclose(fp);
      return(1);
   }

   number_of_overlays = left_cc_overlay + left_mlo_overlay + right_cc_overlay + right_mlo_overlay;

   /****************************************************************************
   * Process the lines associated with each overlay file.
   ****************************************************************************/
   for(i=0;i<number_of_overlays;i++){
      if(get_line_with_astring_in_it(fp, "FILE", line, filename) == 0){
         fclose(fp);
         return(0);
      }
      /* printf("\n%s", line); */
      if(strstr(line, "LEFT_CC") != NULL){
         left_cc_overlay -= 1;
         current_image = 0;
      }
      if(strstr(line, "LEFT_MLO") != NULL){
         left_mlo_overlay -= 1;
         current_image = 2;
      }
      if(strstr(line, "RIGHT_CC") != NULL){
         right_cc_overlay -= 1;
         current_image = 1;
      }
      if(strstr(line, "RIGHT_MLO") != NULL){
         right_mlo_overlay -= 1;
         current_image = 3;
      }

      /*************************************************************************
      * Read in the number of abnormalities.
      *************************************************************************/
      if(get_line_with_astring_in_it(fp, "TOTAL_ABNORMALITIES", line, filename) == 0){
         fclose(fp);
         return(0);
      }
      sscanf(line, "TOTAL_ABNORMALITIES %d", &total_abnormalities);
      /* printf("TOTAL_ABNORMALITIES = %d\n", total_abnormalities); */

      /*************************************************************************
      * Read in the data for each abnormality.
      *************************************************************************/
      for(j=0;j<total_abnormalities;j++){
         if(get_line_with_astring_in_it(fp, "ABNORMALITY", line, filename) == 0){
            fclose(fp);
            return(0);
         }
         if(get_string_with_astring_in_it(fp, "LESION_TYPE", line, filename) == 0){
            fclose(fp);
            return(0);
         }

         /*  PROCESS LESION_TYPE */
         get_description(fp, &description1);
         num_descriptions = 1;

         fscanf(fp, "%s", line);
         if(strstr(line, "LESION_TYPE") != NULL){
            /*  PROCESS LESION_TYPE */
            get_description(fp, &description2);
            num_descriptions = 2;
         }

         if(strstr(line, "ASSESSMENT") == NULL){
            if(get_line_with_astring_in_it(fp, "ASSESSMENT", line, filename) == 0){
               fclose(fp);
               return(0);
            }
            sscanf(line, "ASSESSMENT %d", &assessment);
         }
         else{
            fgets(line, 200, fp);
            sscanf(line, "%d", &assessment);
         }

         if(get_line_with_astring_in_it(fp, "SUBTLETY", line, filename) == 0){
            fclose(fp);
            return(0);
         }
         sscanf(line, "SUBTLETY %d", &subtlety);

         if(get_line_with_astring_in_it(fp, "PATHOLOGY", line, filename) == 0){
            fclose(fp);
            return(0);
         }
         sscanf(line, "PATHOLOGY %s", pathology);

         if(strstr(pathology, "MALIGNANT") != NULL) pathologytype = MALIGNANT;
         else if(strstr(pathology, "BENIGN_WITHOUT_CALLBACK") != NULL) pathologytype = BENIGN_WITHOUT_CALLBACK;
         else if(strstr(pathology, "BENIGN") != NULL) pathologytype = BENIGN;
         else if(strstr(pathology, "UNPROVEN") != NULL) pathologytype = UNPROVEN;
         else{
            fprintf(stderr, "The pathology can not be determined from the html file %s.\n", filename);
            fclose(fp);
            return(0);
         } 

         /*
         printf("%s: DENSITY %d DIGITIZER %s ", filename, density, digitizer);
         printf("ASSESSMENT %d SUBTLETY %d PATHOLOGY %s ", assessment, subtlety, pathology);
         */

         printf("%s", filename);
         printf(" %d", voltype);
         printf(" %d", (int)pow(2.0, (double)density-1));
         printf(" %d", scannertype);
         printf(" %d", pathologytype);
         printf(" %3d", (int)pow(2.0, (double)assessment-1));
         printf(" %3d", (int)pow(2.0, (double)subtlety-1));

         if(num_descriptions == 1){
            if(description1.lesion_type == MASS)
               printf(" 1 %5d %5d     0     0", (int)description1.shapetype.shape, (int)description1.margindistribution.margin);
            if(description1.lesion_type == CALCIFICATION)
               printf(" 2     0     0 %5d %5d", (int)description1.shapetype.type, (int)description1.margindistribution.distribution);
            if(description1.lesion_type == OTHER)
               printf(" 4     0     0 %5d %5d", 0, 0);
         }
         if(num_descriptions == 2){
            if((description1.lesion_type == OTHER) || (description2.lesion_type == OTHER))
               fprintf(stderr, "\nError! A LESION_TYPE of OTHER must be the only LESION_TYPE in a description.\n");
            else if(description1.lesion_type == MASS)
               printf(" 3 %5d %5d %5d %5d", (int)description1.shapetype.shape, (int)description1.margindistribution.margin,
                                     (int)description2.shapetype.type, (int)description2.margindistribution.distribution);
            else if(description1.lesion_type == CALCIFICATION)
               printf(" 3 %5d %5d %5d %5d", (int)description2.shapetype.shape, (int)description2.margindistribution.margin,
                                     (int)description1.shapetype.type, (int)description1.margindistribution.distribution);
         }

         printf("\n");
      }
   }

   if((left_cc_overlay != 0) || (left_mlo_overlay != 0) || (left_cc_overlay != 0) || (left_mlo_overlay != 0)){
      fprintf(stderr, "Error reading the overlay part of the file %s\n", filename);
      fclose(fp);
      return(0);
   }

   fclose(fp);
   return(1);
}

int get_line_with_astring_in_it(FILE *fp, char *astring, char *line, char *filename)
{
   while(fgets(line, 200, fp) && !feof(fp) && (strstr(line, astring) == NULL)){
      if(ferror(fp)){
         fprintf(stderr, "An error occured while reading the file %s.\n", filename);
         fclose(fp);
         return(0);
      }
      if(feof(fp)){
         fprintf(stderr, "The file %s does not contain \"%s\".\n", filename, astring);
         fclose(fp);
         return(0);
      }
   }

   return(1);
}

int get_string_with_astring_in_it(FILE *fp, char *astring, char *line, char *filename)
{
   while(fscanf(fp, "%s", line) && !feof(fp) && (strstr(line, astring) == NULL)){
      if(ferror(fp)){
         fprintf(stderr, "An error occured while reading the file %s.\n", filename);
         fclose(fp);
         return(0);
      }
      if(feof(fp)){
         fprintf(stderr, "The file %s does not contain \"%s\".\n", filename, astring);
         fclose(fp);
         return(0);
      }
   }

   return(1);
}

/*******************************************************************************
* Function: read_overlay_file
* Purpose: This functions reads in all of the data from an overlay file.
* Name: Michael Heath, University of South Florida
* Date: 5/6/98
*******************************************************************************/
int read_overlay_file(char *filename, OVERLAY_DATA *overlay_data)
{
   FILE *fp=NULL;
   char keyword[100];
   int a;

   /****************************************************************************
   * Open the file.
   ****************************************************************************/
   if((fp = fopen(filename, "r")) == NULL){
      fprintf(stderr, "Error opening the OVERLAY file %s for reading.\n", filename);
      return(0);
   }
   if(PRINTEVERYTHING) printf("Successfully opened the file %s.\n", filename);

   /****************************************************************************
   * Read in the total number of abnormalities.
   ****************************************************************************/
   get_keyword(fp, keyword);
   if(check_keyword(keyword, "TOTAL_ABNORMALITIES") == 0) overlay_exit();
   overlay_data->total_abnormalities = get_intval(fp);
   if(PRINTEVERYTHING) printf("TOTAL_ABNORMALITIES %d\n", overlay_data->total_abnormalities);

   /****************************************************************************
   * Allocate memory for all of the abnormalities.
   ****************************************************************************/
   if((overlay_data->abnormalities = (ABNORMALITY *)
      calloc(overlay_data->total_abnormalities, sizeof(ABNORMALITY))) == NULL){
      fprintf(stderr, "Calloc error in read_overlay_file().\n");
      overlay_exit();
   }

   /****************************************************************************
   * Read in each of the abnormalities.
   ****************************************************************************/
   for(a=0;a<overlay_data->total_abnormalities;a++){
      read_abnormality(fp, (overlay_data->abnormalities) + a, a+1);
   }

   fclose(fp);
   return(1);
}

/*******************************************************************************
* Function: read_abnormality
* Purpose: To read in the data associated with one abnormality.
* Name: Michael Heath, University of South Florida
* Date: 5/6/98
*******************************************************************************/
int read_abnormality(FILE *fp, ABNORMALITY *abnormality, int abnormality_number)
{
   int number, description_number, b;
   char keyword[100];

   /****************************************************************************
   * Read in the number of the abnormality.
   ****************************************************************************/
   get_keyword(fp, keyword);
   if(check_keyword(keyword, "ABNORMALITY") == 0) overlay_exit();
   number = get_intval(fp);
   if(number != abnormality_number){
      fprintf(stderr, "Abnormality number different than expected.\n");
      return(0);
   }
   if(PRINTEVERYTHING) printf("ABNORMALITY %d\n", number);

   /****************************************************************************
   * Read in the description. There can be multiple LESION_TYPE, and each can
   * contain a description. Read in up to MAX_DESCRIPTIONS number of them.
   ****************************************************************************/
   get_keyword(fp, keyword);
   description_number = 0;
   while((description_number < MAX_DESCRIPTIONS) && (check_keyword(keyword, "LESION_TYPE"))){
      if(PRINTEVERYTHING) printf("LESION_TYPE ");
      get_description(fp, abnormality->description + description_number);
      description_number++;
      get_keyword(fp, keyword);
   }
   abnormality->num_descriptions = description_number;

   /****************************************************************************
   * Read in the assessment.
   ****************************************************************************/
   /* get_keyword(fp, keyword); Don't need this because it's already read in */
   if(check_keyword(keyword, "ASSESSMENT") == 0) overlay_exit();
   abnormality->assessment = get_intval(fp);
   if(PRINTEVERYTHING) printf("ASSESSMENT %d\n", abnormality->assessment);

   /****************************************************************************
   * Read in the subtlety.
   ****************************************************************************/
   get_keyword(fp, keyword);
   if(check_keyword(keyword, "SUBTLETY") == 0) overlay_exit();
   abnormality->subtlety = get_intval(fp);
   if(PRINTEVERYTHING) printf("SUBTLETY %d\n", abnormality->subtlety);

   /****************************************************************************
   * Read in the pathology..
   ****************************************************************************/
   get_keyword(fp, keyword);
   if(check_keyword(keyword, "PATHOLOGY") == 0) overlay_exit();
   get_string(fp, abnormality->pathology);
   if(PRINTEVERYTHING) printf("PATHOLOGY %s\n", abnormality->pathology);

   /****************************************************************************
   * Read in the total number of outlines.
   ****************************************************************************/
   get_keyword(fp, keyword);
   if(check_keyword(keyword, "TOTAL_OUTLINES") == 0) overlay_exit();
   abnormality->num_cores = get_intval(fp) - 1;
   if(PRINTEVERYTHING) printf("TOTAL_OUTLINES %d\n", abnormality->num_cores + 1);

   if(abnormality->num_cores > MAX_CORES){
      fprintf(stderr, "There are more codes than the program will allow.\n");
      overlay_exit();
   }

   /****************************************************************************
   * Read in each chain code. The first one will be the actual BOUNDARY and
   * additional chain codes will be a CORE.
   ****************************************************************************/
   read_chaincode(fp, &(abnormality->boundary), "BOUNDARY");
   for(b=0;b<abnormality->num_cores;b++){
      read_chaincode(fp, abnormality->core + b, "CORE");
   }

   return(1);
}

/*******************************************************************************
* Function: read_chaincode
* Purpose: This function reads in a chaincode from the file. The type indicates
* which type of outline we are expecting to read. It is either "BOUNDARY" or
* "CORE".
* Name: Michael Heath, University of South Florida
* Date: 5/7/98
*******************************************************************************/
void read_chaincode(FILE *fp, OUTLINE *outline, char *type)
{
   char aline[80000], *beginptr=NULL, *endptr=NULL;
   int numcharacters, length, a;

   /****************************************************************************
   * Read in anything until we either get to the end of the file or to a line
   * that begins with a tab or a B (for BOUNDARY).
   ****************************************************************************/
   do{
      fgets(aline, 80000, fp);
   }while(!feof(fp) && (!((aline[0] == 'B') || (aline[0] == '\t'))) );

   if(feof(fp)) exit(1);

   /****************************************************************************
   * Read the full line of data containing the chain code.
   ****************************************************************************/
   if(strstr(aline, type) == NULL){
      fprintf(stderr, "Found the wrong type of outline in read_chaincode().\n");
      overlay_exit();
   }

   if(PRINTEVERYTHING) printf("%s", type);

   /****************************************************************************
   * Extract the starting location of the chaincode.
   ****************************************************************************/
   fgets(aline, 80000, fp);
   sscanf(aline, "%d %d %n", &(outline->start.c), &(outline->start.r), &numcharacters);
   if(PRINTEVERYTHING) printf(" %d %d %d", outline->start.c, outline->start.r, numcharacters);

   /****************************************************************************
   ****************************************************************************/
   beginptr = aline + numcharacters;
   endptr = strchr(aline, '#');
   length = (int)(endptr - beginptr) / 2;

/*
   for(a=0;a<length;a++) printf("%c", beginptr[a*2]);

   CHAINPOINT start;  * It could be a boundary or a core. The chaincode  *
   int length;        * is converted to points and the array of points is *
   CHAINPOINT *chain; * stored. The start point is stored in start and in the array. *
*/

   /****************************************************************************
   * Allocate memory to store the outline. Fill the chain code into the this
   * memory. The chaincode is converted into an array of points. The first
   * point is included in this list of points.
   ****************************************************************************/
   outline->length = length + 1;  /* Allow for the start point. */
   if((outline->chain = (CHAINPOINT *) calloc(outline->length, sizeof(CHAINPOINT))) == NULL){
      fprintf(stderr, "Memory allocation failure in read_chaincode().\n");
      overlay_exit();
   }

   outline->chain[0].r = outline->start.r;
   outline->chain[0].c = outline->start.c;

   for(a=0;a<length;a++){
      switch(beginptr[a*2]){
         case '0':
            outline->chain[a+1].r = outline->chain[a].r - 1;
            outline->chain[a+1].c = outline->chain[a].c + 0;
            break;
         case '1':
            outline->chain[a+1].r = outline->chain[a].r - 1;
            outline->chain[a+1].c = outline->chain[a].c + 1;
            break;
         case '2':
            outline->chain[a+1].r = outline->chain[a].r + 0;
            outline->chain[a+1].c = outline->chain[a].c + 1;
            break;
         case '3':
            outline->chain[a+1].r = outline->chain[a].r + 1;
            outline->chain[a+1].c = outline->chain[a].c + 1;
            break;
         case '4':
            outline->chain[a+1].r = outline->chain[a].r + 1;
            outline->chain[a+1].c = outline->chain[a].c + 0;
            break;
         case '5':
            outline->chain[a+1].r = outline->chain[a].r + 1;
            outline->chain[a+1].c = outline->chain[a].c - 1;
            break;
         case '6':
            outline->chain[a+1].r = outline->chain[a].r + 0;
            outline->chain[a+1].c = outline->chain[a].c - 1;
            break;
         case '7':
            outline->chain[a+1].r = outline->chain[a].r - 1;
            outline->chain[a+1].c = outline->chain[a].c - 1;
            break;
         default:
            fprintf(stderr, "Error reading the chaincode (invald value).\n");
            overlay_exit();
      }
   }
   if(PRINTEVERYTHING) printf("(%d %d)\n", outline->chain[a].c, outline->chain[a].r);
}

/*******************************************************************************
* Function: get_description
* Purpose: This function reads in the description of an abnormaility from the
* file. The lesion should be either a MASS or a CALCIFICATION. If it is the
* former it will have SHAPE and MARGINS data oftherwise it will have a TYPE and
* DISTRIBUTION. Multiple values can be specified for these four fields. When
* multiple values are supplied, they will be separated with the '_' character.
* Name: Michael Heath, University of South Florida
* Date: 5/7/98
*******************************************************************************/
int get_description(FILE *fp, DESCRIPTION *description)
{
   char keyword[100];
   char line[200];
   char *next_shape_ptr = NULL, *next_margin_ptr=NULL;
   char *next_type_ptr=NULL, *next_distribution_ptr=NULL;
   int whichone;
   char *return_ptr=NULL;

   /****************************************************************************
   * Read in the type of the lesion.
   ****************************************************************************/
   get_keyword(fp, keyword);
   if(check_keyword(keyword, "MASS") == 1){
      description->lesion_type = MASS;
      if(PRINTEVERYTHING) printf("MASS");
      fgets(line, 200, fp);
      if((return_ptr = strchr(line, '\n')) != NULL) *return_ptr = '\0';
      next_shape_ptr = strstr(line, "SHAPE");
      next_margin_ptr = strstr(line, "MARGINS");
      if((next_shape_ptr == NULL) || (next_margin_ptr == NULL)){
         fprintf(stderr, "Error reading the description of a lesion type.\n");
         overlay_exit();
      }
      *(next_margin_ptr-1) = '\0';

      if(PRINTEVERYTHING) printf("\n[%s] ", next_shape_ptr);
      description->shapetype.shape = 0;
      while((whichone = inlist(&next_shape_ptr, MASS_SHAPE_string)) != -1){
         description->shapetype.shape += valtobit(whichone);
      }

      if(PRINTEVERYTHING) printf("\n[%s] ", next_margin_ptr);
      description->margindistribution.margin = 0;
      while((whichone = inlist(&next_margin_ptr, MASS_MARGIN_string)) != -1){
         description->margindistribution.margin += valtobit(whichone);
      }
      if(PRINTEVERYTHING){
         printf("\n");
         printf("SHAPE=%d\n", description->shapetype.shape);
         printf("MARGIN=%d\n", description->margindistribution.margin);
      }
   }
   else if(check_keyword(keyword, "CALCIFICATION") == 1){
      description->lesion_type = CALCIFICATION;
      if(PRINTEVERYTHING) printf("CALCIFICATION");
      fgets(line, 200, fp);
      if((return_ptr = strchr(line, '\n')) != NULL) *return_ptr = '\0';
      next_type_ptr = strstr(line, "TYPE");
      next_distribution_ptr = strstr(line, "DISTRIBUTION");
      if((next_type_ptr == NULL) || (next_distribution_ptr == NULL)){
         fprintf(stderr, "Error reading the description of a lesion type.\n");
         overlay_exit();
      }
      *(next_distribution_ptr-1) = '\0';

      if(PRINTEVERYTHING) printf("\n[%s] ", next_type_ptr);
      description->shapetype.type = 0;
      while((whichone = inlist(&next_type_ptr, CALCIFICATION_TYPE_string)) != -1){
         description->shapetype.type += valtobit(whichone);
      }

      if(PRINTEVERYTHING) printf("\n[%s] ", next_distribution_ptr);
      description->margindistribution.distribution = 0;
      while((whichone = inlist(&next_distribution_ptr, CALCIFICATION_DISTRIBUTION_string)) != -1){
         description->margindistribution.distribution += valtobit(whichone);
      }
      if(PRINTEVERYTHING){
         printf("\n");
         printf("TYPE=%d\n", description->shapetype.type);
         printf("DISTRIBUTION=%d\n", description->margindistribution.distribution);
      }
   }
   else if(check_keyword(keyword, "OTHER") == 1){
      description->lesion_type = OTHER;
      if(PRINTEVERYTHING) printf("OTHER");
      fgets(line, 200, fp);
      if((return_ptr = strchr(line, '\n')) != NULL) *return_ptr = '\0';

      /*
      next_type_ptr = strstr(line, "TYPE");
      next_distribution_ptr = strstr(line, "DISTRIBUTION");
      if((next_type_ptr == NULL) || (next_distribution_ptr == NULL)){
         fprintf(stderr, "Error reading the description of a lesion type.\n");
         overlay_exit();
      }
      *(next_distribution_ptr-1) = '\0';

      if(PRINTEVERYTHING) printf("\n[%s] ", next_type_ptr);
      description->shapetype.type = 0;
      while((whichone = inlist(&next_type_ptr, CALCIFICATION_TYPE_string)) != -1){
         description->shapetype.type += valtobit(whichone);
      }

      if(PRINTEVERYTHING) printf("\n[%s] ", next_distribution_ptr);
      description->margindistribution.distribution = 0;
      while((whichone = inlist(&next_distribution_ptr, CALCIFICATION_DISTRIBUTION_string)) != -1){
         description->margindistribution.distribution += valtobit(whichone);
      }
      */
      if(PRINTEVERYTHING){
         printf("\n");
         /*
         printf("TYPE=%d\n", description->shapetype.type);
         printf("DISTRIBUTION=%d\n", description->margindistribution.distribution);
         */
      }

   }
   else{
      fprintf(stderr, "Invalid LESION_TYPE.\n");
      return(0);
   }

   return(1);
}

/*******************************************************************************
* Function: inlist
* Purpose: To find the index of the string in the string list that first
* appears in the string pointed to by stringptr. If no match is found then -1
* is the value returned. Otherwise, the index of the string in the array of
* strings is returned. The pointer to stringptr is changed to NULL in the
* former case and to the next character after the end of the matching string in
* the latter case.
* Name: Michael Heath, University of South Florida
* Date: 5/7/98
*******************************************************************************/
int inlist(char **stringptr, char *string_list[])
{
   char *hold_ptr=NULL, *ptrval=NULL, *lowest_ptrval;
   int lowest_index=-1, index=0;
   
   hold_ptr = *stringptr;
   lowest_ptrval = *stringptr + strlen(*stringptr);

   index = 0;
   while(string_list[index] != NULL){
      /* if(PRINTEVERYTHING) printf("(TRYING %s)", string_list[index]); */
      ptrval = strstr(hold_ptr, string_list[index]);
      if(ptrval != NULL){
         if(ptrval < lowest_ptrval){
            lowest_ptrval = ptrval;
            lowest_index = index;
         }
      }
      index++;
   }

   if(lowest_index == -1){
      /* if(PRINTEVERYTHING) printf("  NOMATCH_RETURNING\n"); */
      *stringptr = NULL;
      return(-1);
   }

   if(PRINTEVERYTHING) printf("%s ", string_list[lowest_index]);

   *stringptr = lowest_ptrval + strlen(string_list[lowest_index]);

   return(lowest_index);
}

/*******************************************************************************
* Function: get_keyword
* Purpose: To scan in the next "word" from the file. It should be a keyword.
* Name: Michael Heath, University of South Florida
* Date: 5/7/98
*******************************************************************************/
void get_keyword(FILE *fp, char *keyword)
{
   fscanf(fp, "%s", keyword);
}

/*******************************************************************************
* Function:
* Purpose: To check if the string passed to the function is equal to the
* "truth" string.
* Name: Michael Heath, University of South Florida
* Date: 5/7/98
*******************************************************************************/
int check_keyword(char *keyword, const char *truth)
{
   if(strcmp(keyword, truth) == 0) return(1);
   else return(0);
}

/*******************************************************************************
* Function: get_intval
* Purpose: To read the next integer from the file.
* Name: Michael Heath, University of South Florida
* Date: 5/7/98
*******************************************************************************/
int get_intval(FILE *fp)
{
   int intval;
   fscanf(fp, "%d", &intval);
   return(intval);
}

/*******************************************************************************
* Function: get_string
* Purpose: To read in a string from the file.
* Name: Michael Heath, University of South Florida
* Date: 5/7/98
*******************************************************************************/
void get_string(FILE *fp, char *thestring)
{
   fscanf(fp, "%s", thestring);
}

/*******************************************************************************
* Function: overlay_exit
* Purpose: To print a message and exit.
* Name: Michael Heath, University of South Florida
* Date: 5/7/98
*******************************************************************************/
void overlay_exit()
{
   fprintf(stderr, "Error reading the overlay file.\n");
   exit(1);
}

/*******************************************************************************
* Function: valtobit
* Purpose: To convert a value to the value of an integer with that bit set.
* Name: Michael Heath, University of South Florida
* Date: 5/7/98
*******************************************************************************/
int valtobit(int value)
{
   switch(value){
      case 0: return(1);
      case 1: return(2);
      case 2: return(4);
      case 3: return(8);
      case 4: return(16);
      case 5: return(32);
      case 6: return(64);
      case 7: return(128);
      case 8: return(256);
      case 9: return(512);
      case 10: return(1024);
      case 11: return(2048);
      case 12: return(4096);
      case 13: return(8192);
      case 14: return(16384);
      case 15: return(32768);
      default:
         fprintf(stderr, "Error in valtobit().\n");
         exit(1);
   }
}

/*******************************************************************************
* Function: bittoval
* Purpose: To convert a value to the bit number in an integer that has that
* value.
* Name: Michael Heath, University of South Florida
* Date: 5/7/98
*******************************************************************************/
int bittoval(int bit)
{
   switch(bit){
      case 0: return(-1);
      case 1: return(0);
      case 2: return(1);
      case 4: return(2);
      case 8: return(3);
      case 16: return(4);
      case 32: return(5);
      case 64: return(6);
      case 128: return(7);
      case 256: return(8);
      case 512: return(9);
      case 1024: return(10);
      case 2048: return(11);
      case 4096: return(12);
      case 8192: return(13);
      case 16384: return(14);
      case 32768: return(15);
      default:
         fprintf(stderr, "Error in bittoval().\n");
         exit(1);
   }
}

/*******************************************************************************
* Below is code for writing data to an OVERLAY file.
*******************************************************************************/

/*******************************************************************************
* Function: write_overlay_file
* Purpose: To write the data from the OVERLAY_DATA structure out into a file.
* Name: Michael Heath, University of South Florida
* Date: 5/8/98
*******************************************************************************/
int write_overlay_file(char *filename, OVERLAY_DATA overlay_data)
{
   FILE *fp;
   int a;

   /****************************************************************************
   * Open the overlay file for writing.
   ****************************************************************************/
   if((fp = fopen(filename, "w")) == NULL){
      fprintf(stderr, "Error opening the file %s for writing.\n", filename);
      overlay_exit();
   }

   /****************************************************************************
   * Write the total number of abnormalities.
   ****************************************************************************/
   fprintf(fp, "TOTAL_ABNORMALITIES %d\n", overlay_data.total_abnormalities);

   /****************************************************************************
   * Write each of the abnormalities to the file.
   ****************************************************************************/
   for(a=0;a<overlay_data.total_abnormalities;a++){
      write_abnormality(fp, (overlay_data.abnormalities) + a, a+1);
   }

   fclose(fp);
   return(1);
}

/*******************************************************************************
* Function: write_abnormality
* Purpose: To write out the data associated with one abnormality.
* Name: Michael Heath, University of South Florida
* Date: 5/8/98
*******************************************************************************/
int write_abnormality(FILE *fp, ABNORMALITY *abnormality, int abnormality_number)
{
   int description_number, b;

   /****************************************************************************
   * Write the number of the abnormality.
   ****************************************************************************/
   fprintf(fp, "ABNORMALITY %d\n", abnormality_number);

   /****************************************************************************
   * Write the description. There can be multiple LESION_TYPE, and each can
   * contain a description.
   ****************************************************************************/
   for(description_number=0;description_number<abnormality->num_descriptions;description_number++){
      fprintf(fp, "LESION_TYPE ");
      put_description(fp, abnormality->description + description_number);
   }

   /****************************************************************************
   * Write the assessment.
   ****************************************************************************/
   fprintf(fp, "ASSESSMENT %d\n", abnormality->assessment);

   /****************************************************************************
   * Write the subtlety.
   ****************************************************************************/
   fprintf(fp, "SUBTLETY %d\n", abnormality->subtlety);

   /****************************************************************************
   * Write the pathology..
   ****************************************************************************/
   fprintf(fp, "PATHOLOGY %s\n", abnormality->pathology);

   /****************************************************************************
   * Write the total number of outlines.
   ****************************************************************************/
   fprintf(fp, "TOTAL_OUTLINES %d \n", abnormality->num_cores + 1);

   /****************************************************************************
   * Write out each chain code. The first one will be the actual BOUNDARY and
   * additional chain codes will be a CORE.
   ****************************************************************************/
   write_chaincode(fp, &(abnormality->boundary), "BOUNDARY");
   for(b=0;b<abnormality->num_cores;b++){
      write_chaincode(fp, abnormality->core + b, "CORE");
   }

   return(1);
}

/*******************************************************************************
* Function: put_description
* Purpose: This function writes the description of an abnormaility to a
* file. The lesion should be either a MASS or a CALCIFICATION. If it is the
* former it will have SHAPE and MARGINS data oftherwise it will have a TYPE and
* DISTRIBUTION. Multiple values can be specified for these four fields. When
* multiple values are supplied, they will be separated with the '_' character.
* Name: Michael Heath, University of South Florida
* Date: 5/7/98
*******************************************************************************/
int put_description(FILE *fp, DESCRIPTION *description)
{
   /****************************************************************************
   * Write out the description of a MASS.
   ****************************************************************************/
   if(description->lesion_type == MASS){

      fprintf(fp, "MASS ");

      fprintf(fp, "SHAPE ");
      write_keywords(fp, (int)(description->shapetype.shape), MASS_SHAPE_string);

      fprintf(fp, " MARGINS ");
      write_keywords(fp, (int)(description->margindistribution.margin), MASS_MARGIN_string);

      fprintf(fp, "\n");
   }
   /****************************************************************************
   * Write out the description of a CALCIFICATION.
   ****************************************************************************/
   else if(description->lesion_type == CALCIFICATION){

      fprintf(fp, "CALCIFICATION ");

      fprintf(fp, "TYPE ");
      write_keywords(fp, (int)(description->shapetype.type), CALCIFICATION_TYPE_string);

      fprintf(fp, " DISTRIBUTION ");
      write_keywords(fp, (int)(description->margindistribution.distribution), CALCIFICATION_DISTRIBUTION_string);

      fprintf(fp, "\n");
   }
   /****************************************************************************
   * Write out the description of a OTHER.
   ****************************************************************************/
   else if(description->lesion_type == OTHER){

      fprintf(fp, "OTHER");

      /*
      fprintf(fp, "TYPE ");
      write_keywords(fp, (int)(description->shapetype.type), CALCIFICATION_TYPE_string);

      fprintf(fp, " DISTRIBUTION ");
      write_keywords(fp, (int)(description->margindistribution.distribution), CALCIFICATION_DISTRIBUTION_string);
      */

      fprintf(fp, "\n");
   }
   else{
      fprintf(stderr, "Invalid LESION_TYPE.\n");
      return(0);
   }

   return(1);
}

/*******************************************************************************
* Function: write_keywords
* Purpose: This function writes out the keywords for the lesion.
* Name: Michael Heath, University of South Florida
* Date: 5/9/98
*******************************************************************************/
void write_keywords(FILE *fp, int bitmask, char *string_list[])
{
   int num_onbits=0, num, whichbit;

   /****************************************************************************
   * If no bits are on we just write N/A and return because there are no
   * keywords specified.
   ****************************************************************************/
   if(bitmask == 0){
      fprintf(fp, "N/A");
      return;
   }

   /****************************************************************************
   * Count the number of bits that are on.
   ****************************************************************************/
   num = bitmask;
   while(num != 0){
      if(num % 2) num_onbits++;
      num = num>>1;
   }

   /****************************************************************************
   * Write out the keywords.
   ****************************************************************************/
   num = bitmask;
   whichbit = 1;
   while(num != 0){
      if(num % 2){
         fprintf(fp, "%s", string_list[whichbit-1]);
         num_onbits--;
         if(num_onbits != 0) fprintf(fp, "-");
      }
      num = num>>1;
      whichbit++;
   }
}

/*******************************************************************************
* Function: write_chaincode
* Purpose: This function converts an outline to a chaincode and writes it 
* to a file. The type indicates which type of boundary it is. It is either
* "BOUNDARY" or "CORE".
* Name: Michael Heath, University of South Florida
* Date: 5/9/98
*******************************************************************************/
void write_chaincode(FILE *fp, OUTLINE *outline, char *type)
{
   int a, delta_r, delta_c;
   int directions[3][3] = {7, 0, 1, 6, 8, 2, 5, 4, 3};    /* 8 is never used. */

   /****************************************************************************
   * Write out the type of the boundary.
   ****************************************************************************/
   if(strcmp(type, "BOUNDARY")==0) fprintf(fp, "%s\n", type);
   if(strcmp(type, "CORE")==0) fprintf(fp, "\t%s\n", type);

   /****************************************************************************
   * Write out the starting location (column and then row).
   ****************************************************************************/
   fprintf(fp, "%d %d ", outline->start.c, outline->start.r);

   /****************************************************************************
   * Write out the chaincode.
   ****************************************************************************/
   for(a=1;a<outline->length;a++){
      delta_r = outline->chain[a].r - outline->chain[a-1].r;
      delta_c = outline->chain[a].c - outline->chain[a-1].c;

      fprintf(fp, "%d ", directions[1+delta_r][1+delta_c]);
   }
   fprintf(fp, "#\n");
}

/*******************************************************************************
* Function: determine_valid_regions
* Purpose: The overlay data is very complicated. We may only want some
* subset of the overlay regions so we can use this function to fill an array
* with ones or zeros indicating which regions we want to use. We can limit
* the lesion type, the pathology, the mass shape, the mass margins, the
* calcification type or the calcification distribution. In any case no regions
* without a pathology of BENIGN or MALIGNANT are to be used.
* Name: Michael Heath, University of South Florida
* Date: 1/28/2000
*******************************************************************************/
int determine_valid_regions(OVERLAY_DATA overlay_data, char *lesiontype,
    char *pathology, char *ms, char *mm, char *ct, char *cd, int **isvalid)
{
   ABNORMALITY abnormality;
   int d, c, isvalid_num = 0;
   char *next_desc_ptr = NULL;
   int description_int;
   int whichone;

   if(overlay_data.total_abnormalities == 0) return;

   /*************************************************************************
   * Allocate an array of integers to use for indicating which regions are
   * valid. That is they match the lesion_type and all relevent descriptions.
   *************************************************************************/
   if(((*isvalid) = (int *) calloc(overlay_data.total_abnormalities, sizeof(int))) == NULL){
      fprintf(stderr, "Calloc error in determine_valid_regions!\n");
      exit(1);
   }

   /*******************************************************************************
   * Do not use any ground truth outline that is not "MALIGNANT" or "BENIGN".
   * For example we will skip regions that are "BENIGN_WITHOUT_CALLBACK".
   * We will also skip "MASS" regions if the user specified "-l CALCIFICATION",
   * or "CALCIFICATION" regions if the user specified "-l MASS".
   ****************************************************************************/
   for(d=0;d<overlay_data.total_abnormalities;d++){
      abnormality = overlay_data.abnormalities[d];

      (*isvalid)[d] = 1;

      /*************************************************************************
      * If the ground truth region is neither "MALIGNANT" nor "BENIGN" we will
      * skip it.
      *************************************************************************/
      if(!((strcmp(abnormality.pathology, "MALIGNANT") == 0) || (strcmp(abnormality.pathology, "BENIGN") == 0))){
         (*isvalid)[d] = 0;
      }

      /*************************************************************************
      * If the pathology is not NULL the calling function asked this function
      * to limit the selection of abnormalities to either BENIGN or MALIGNANT.
      *************************************************************************/
      if(pathology != NULL){
         if(!(strcmp(abnormality.pathology, pathology) == 0)) (*isvalid)[d] = 0;
      }

      /*************************************************************************
      * If we are not yet skipping the region, check to see that the
      * lesion_type of the ground truth region is one that we are trying to find.
      * Remember that a ground truth region can have more than one lesion type.
      *************************************************************************/
      if((*isvalid)[d] == 1){
         if(strcmp(lesiontype, "CALCIFICATION") == 0){
            (*isvalid)[d] = 0;
            for(c=0;c<abnormality.num_descriptions;c++){
               if(abnormality.description[c].lesion_type == CALCIFICATION) (*isvalid)[d] = 1;
            }
         }
         if(strcmp(lesiontype, "MASS") == 0){
            (*isvalid)[d] = 0;
            for(c=0;c<abnormality.num_descriptions;c++){
               if(abnormality.description[c].lesion_type == MASS) (*isvalid)[d] = 1;
            }
         }
         if(strcmp(lesiontype, "BOTH") == 0){
            (*isvalid)[d] = 0;
            for(c=0;c<abnormality.num_descriptions;c++){
               if(abnormality.description[c].lesion_type == CALCIFICATION) (*isvalid)[d] = 1;
               if(abnormality.description[c].lesion_type == MASS) (*isvalid)[d] = 1;
            }
         }
         /* isvalid_num++; */
      }

      /**********************************************************************
      * If we are not yet skipping the region, check to see that the
      * description of the ground truth region is one that we are trying to find.
      * Remember that a ground truth region can have more than one lesion type.
      **********************************************************************/
      if((*isvalid)[d] == 1){
         if(strcmp(lesiontype, "CALCIFICATION") == 0){
            for(c=0;c<abnormality.num_descriptions;c++){
               if(abnormality.description[c].lesion_type == CALCIFICATION){

                  if(ct != NULL){
                     (*isvalid)[d] -= 1;
                     next_desc_ptr = ct;
                     description_int = 0;
                     while((whichone = inlist(&next_desc_ptr, CALCIFICATION_TYPE_string)) != -1){
                        description_int += valtobit(whichone);
                     }

                     if((description_int & abnormality.description[c].shapetype.type) != 0) (*isvalid)[d] += 1;
                  }

                  if(cd != NULL){
                     (*isvalid)[d] -= 1;
                     next_desc_ptr = cd;
                     description_int = 0;
                     while((whichone = inlist(&next_desc_ptr, CALCIFICATION_DISTRIBUTION_string)) != -1){
                        description_int += valtobit(whichone);
                     }

                     if((description_int & abnormality.description[c].margindistribution.distribution) != 0)
                        (*isvalid)[d] += 1;
                  }

               }
            }
         }

         if(strcmp(lesiontype, "MASS") == 0){
            for(c=0;c<abnormality.num_descriptions;c++){
               if(abnormality.description[c].lesion_type == MASS){

                  if(ms != NULL){
                     (*isvalid)[d] -= 1;
                     next_desc_ptr = ms;
                     description_int = 0;
                     while((whichone = inlist(&next_desc_ptr, MASS_SHAPE_string)) != -1){
                        description_int += valtobit(whichone);
                     }

                     if((description_int & abnormality.description[c].shapetype.shape) != 0) (*isvalid)[d] += 1;
                  }

                  if(mm != NULL){
                     (*isvalid)[d] -= 1;
                     next_desc_ptr = mm;
                     description_int = 0;
                     while((whichone = inlist(&next_desc_ptr, MASS_MARGIN_string)) != -1){
                        description_int += valtobit(whichone);
                     }

                     if((description_int & abnormality.description[c].margindistribution.margin) != 0)
                        (*isvalid)[d] += 1;
                  }

               }
            }
         }
   
         if((*isvalid)[d] == 1) isvalid_num++;
      }
   }

   /****************************************************************************
   * If we are in verbose mode, we can print out which abnormalities we are
   * skipping and which ones we are using.
   ****************************************************************************/
   if(VERBOSE){
      for(d=0;d<overlay_data.total_abnormalities;d++){
         abnormality = overlay_data.abnormalities[d];

         printf("   ABNORMALITY_%02d:", d+1);
         if((*isvalid)[d] == 0) printf(" SKIP");
         if((*isvalid)[d] == 1) printf(" USE ");
         printf(" %s ", abnormality.pathology);
         for(c=0;c<abnormality.num_descriptions;c++){
            /*
            if(abnormality.description[c].lesion_type == MASS) printf(" MASS");
            if(abnormality.description[c].lesion_type == CALCIFICATION) printf(" CALCIFICATION");
            */
            put_description(stdout, abnormality.description + c);
         }

         printf("\n");

      }
   }
}
