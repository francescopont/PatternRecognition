/*******************************************************************************
* Program: simpletiff
* Purpose: This code implements some functions for writing 8-bit and 16-bit
* images in the Tag Image File Format (TIFF). No read functions are provided.
* These functions simply help to write images out to files that can then be
* read by other programs. Please note that 16-bit TIFF images are not handled
* well by most programs. Many programs just use the upper or lower byte of
* data for each pixel. When we viewed 16-bit images in BigEndian (4d4d) format
* with different software this is what we experienced.
*   1) xv v3.10a displayed the most significant byte for each pixel
*   2) display (ImageMagick) displayed the least significant byte
*   3) matlab read the 16-bit data but displayed it very poorly - All 16bits
*      of data were there and were accessible.
*
* Name: Michael Heath, University of South Florida
*       Dr. Melanie Sutton, University of West Florida
* Date: 1/25/2000
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "virtual_image.h"
#include "myalloc.h"

struct TIFFHEADER{
   unsigned short int byteorder;
   unsigned short int version;
   unsigned long int ifdoffset;
};

struct TIFFTAG{
   unsigned short int tag;
   unsigned short int type;
   unsigned long int length;
   unsigned long int valoffset;
};

struct TIFFIFD{
   unsigned short int count;
   unsigned long int directory; /* struct TIFFTAG *directory; */
   unsigned long int nextifd;   /* struct TIFFIFD *nextifd;   */
};

/*
int main(int argc, char *argv[])
{
   unsigned char *image;
   int rows=256, cols=256, r, c;

   int write_simple_tiff_uchar(char *filename, unsigned char *image, int rows, int cols);

   image = (unsigned char *) calloc(rows*cols, sizeof(unsigned char));

   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         if(r > c) image[r*cols+c] = (unsigned char)r;
         else image[r*cols+c] = (unsigned char)c;
      }
   }

   write_simple_tiff_uchar("test.tif", image, rows, cols);
}
*/

/*******************************************************************************
* Function: write_simple_tiff_uchar
* Purpose: To write 8-bit image data to a TIFF file.
* Name: Michael Heath, University of South Florida
* Date: 1/25/2000
*******************************************************************************/
int write_simple_tiff_uchar(char *filename, unsigned char *image, int rows, int cols)
{
   FILE *fp=NULL;
   char *softwarename = "usfheath version 1.0.0";
   char *descriptionname = "This image was produced using software obtained from:\nhttp://marathon.csee.usf.edu/Mammography/software/heathusf.html";
   struct TIFFHEADER tiffheader = {0, 42, 0};
   struct TIFFIFD tiffifd = {13, 0, 0};
   struct TIFFTAG newsubfiletype = {254, 4, 1, 0};
   struct TIFFTAG imagewidth = {256, 3, 1, 0};
   struct TIFFTAG imagelength = {257, 3, 1, 0};
   struct TIFFTAG bitspersample = {258, 3, 1, 8};
   struct TIFFTAG compression = {259, 3, 1, 1};
   struct TIFFTAG photometric = {262, 3, 1, 1};
   struct TIFFTAG description = {270, 2, 0, 0};
   struct TIFFTAG stripoffsets = {273, 4, 0, 0};
   struct TIFFTAG rowsperstrip = {278, 3, 1, 1};
   struct TIFFTAG stripbytecounts = {279, 3, 1, 1};
   struct TIFFTAG planarconfiguration = {284, 3, 1, 1};
   struct TIFFTAG software = {305, 2, 0, 1};
   struct TIFFTAG sampleformat = {339, 3, 1, 1};
   unsigned long int *stripoffsetarr = NULL;
   unsigned short int *stripbytecountarr=NULL;
   unsigned long int ifdfilepos = 0;
   int r, rr, s, rowsperstripvalue, stripsperimage;
   char blankchar = '\0';

   unsigned short int get_byte_order();

   /****************************************************************************
   * Get the unsigned short int code for the byte order of this machine.
   ****************************************************************************/
   tiffheader.byteorder = get_byte_order();

   /****************************************************************************
   * Open the file for writing.
   ****************************************************************************/
   if((fp = fopen(filename, "w+b")) == NULL){
      fprintf(stderr, "Error opening the file %s for writing!\n", filename);
      return(0);
   }

   /****************************************************************************
   * Write the TIFF header to the file. Note that we will need to go back
   * later to write in the position of the Image File Directory because we do
   * not yet know where it will be in the file.
   ****************************************************************************/
   fwrite(&tiffheader, sizeof(struct TIFFHEADER), 1, fp);

   /****************************************************************************
   * Calculate the rows per strip to get nearly 8 kilobytes per strip. Also
   * compute the number of strips we need for the image.
   ****************************************************************************/
   r = 0;
   do{
      rowsperstripvalue = (8 * 1024 - r)/ (cols);
      if(rowsperstripvalue < 1) rowsperstripvalue = 1;
      stripsperimage = (rows + rowsperstripvalue - 1) / rowsperstripvalue;
      r++;
   }while(stripsperimage <= 2);

   /****************************************************************************
   * Fill in the image width value. Since it is less than 4 bytes it needs to
   * be left justified. If we are on a BigEndian machine shift the data 16bits.
   ****************************************************************************/
   imagewidth.valoffset = cols;
   if(tiffheader.byteorder == 0x4d4d)
      imagewidth.valoffset <<= 16;

   /****************************************************************************
   * Fill in the image length value. Since it is less than 4 bytes it needs to
   * be left justified. If we are on a BigEndian machine shift the data 16bits.
   ****************************************************************************/
   imagelength.valoffset = rows;
   if(tiffheader.byteorder == 0x4d4d)
      imagelength.valoffset <<= 16;

   /****************************************************************************
   * If we are on a BigEndian machine shift the data 16bits.
   ****************************************************************************/
   rowsperstrip.valoffset = rowsperstripvalue;
   if(tiffheader.byteorder == 0x4d4d){
      bitspersample.valoffset <<= 16;
      compression.valoffset <<= 16;
      photometric.valoffset <<= 16;
      planarconfiguration.valoffset <<= 16;
      sampleformat.valoffset <<= 16;
      rowsperstrip.valoffset <<= 16;
   }

   /****************************************************************************
   * Write the description to the file.
   ****************************************************************************/
   description.valoffset = ftell(fp);
   fwrite(descriptionname, sizeof(char), strlen(descriptionname), fp);
   description.length = strlen(descriptionname);
   if(description.length % 2){
      fputc((unsigned char)0, fp);
      description.length += 1;
   }
   else{
      fputc((unsigned char)0, fp);
      fputc((unsigned char)0, fp);
      description.length += 2;
   }

   /****************************************************************************
   * Write the name of the software and the version number to the file.
   ****************************************************************************/
   software.valoffset = ftell(fp);
   fwrite(softwarename, sizeof(char), strlen(softwarename), fp);
   software.length = strlen(softwarename);
   if(software.length % 2){
      fputc((unsigned char)0, fp);
      software.length += 1;
   }
   else{
      fputc((unsigned char)0, fp);
      fputc((unsigned char)0, fp);
      software.length += 2;
   }

   /****************************************************************************
   * Write the image data to the file and keep track of the strip offsets and
   * the length in bytes of each strip.
   ****************************************************************************/
   stripoffsets.length = stripsperimage;
   stripoffsetarr = (unsigned long int *) calloc(stripsperimage, sizeof(unsigned long int));
   stripbytecounts.length = stripsperimage;
   stripbytecountarr = (unsigned short int *) calloc(stripsperimage, sizeof(unsigned short int));
   for(s=0,r=0;s<stripsperimage;s++){ 
      stripoffsetarr[s] = ftell(fp);
      for(rr=0;rr<rowsperstripvalue;rr++,r++){
         if(r < rows){
            fwrite(&image[r*cols], sizeof(unsigned char), cols, fp);
            stripbytecountarr[s] += sizeof(unsigned char) * cols;
         }
      }
      if(stripbytecountarr[s] % 2){
         fputc((unsigned char)0, fp);
         stripbytecountarr[s] += 1;
      }
   }

   /****************************************************************************
   * Write the file offset for each strip.
   ****************************************************************************/
   stripoffsets.valoffset =  ftell(fp);
   fwrite(stripoffsetarr, sizeof(unsigned long int), stripsperimage, fp);

   /****************************************************************************
   * Write byte counts for each strip.
   ****************************************************************************/
   stripbytecounts.valoffset = ftell(fp);
   fwrite(stripbytecountarr, sizeof(unsigned short int), stripsperimage, fp);

   /****************************************************************************
   * Write the image file directory.
   ****************************************************************************/
   ifdfilepos = ftell(fp);
   fwrite(&(tiffifd.count), sizeof(unsigned short int), 1, fp);
   fwrite(&newsubfiletype, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&imagewidth, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&imagelength, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&bitspersample, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&compression, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&photometric, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&description, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&stripoffsets, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&rowsperstrip, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&stripbytecounts, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&planarconfiguration, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&software, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&sampleformat, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&(tiffifd.nextifd), sizeof(unsigned long int), 1, fp);

   /****************************************************************************
   * Go back and write the position of the image file directory in the header.
   ****************************************************************************/
   fseek(fp, 4, 0);
   fwrite(&ifdfilepos, sizeof(unsigned long int), 1, fp);

   free(stripoffsetarr);
   free(stripbytecountarr);
   fclose(fp);
   return(1);
}

/*******************************************************************************
* Function: write_simple_tiff_ushort
* Purpose: To write 16-bit image data to a TIFF file.
* Name: Michael Heath, University of South Florida
* Date: 1/25/2000
*******************************************************************************/
int write_simple_tiff_ushort(char *filename, unsigned short int *image, int rows, int cols)
{
   FILE *fp=NULL;
   char *softwarename = "usfheath version 1.0.0";
   char *descriptionname = "This image was produced using software obtained from:\nhttp://marathon.csee.usf.edu/Mammography/software/heathusf.html";
   struct TIFFHEADER tiffheader = {0, 42, 0};
   struct TIFFIFD tiffifd = {13, 0, 0};
   struct TIFFTAG newsubfiletype = {254, 4, 1, 0};
   struct TIFFTAG imagewidth = {256, 3, 1, 0};
   struct TIFFTAG imagelength = {257, 3, 1, 0};
   struct TIFFTAG bitspersample = {258, 3, 1, 16};
   struct TIFFTAG compression = {259, 3, 1, 1};
   struct TIFFTAG photometric = {262, 3, 1, 1};
   struct TIFFTAG description = {270, 2, 0, 0};
   struct TIFFTAG stripoffsets = {273, 4, 0, 0};
   struct TIFFTAG rowsperstrip = {278, 3, 1, 1};
   struct TIFFTAG stripbytecounts = {279, 3, 1, 1};
   struct TIFFTAG planarconfiguration = {284, 3, 1, 1};
   struct TIFFTAG software = {305, 2, 0, 1};
   struct TIFFTAG sampleformat = {339, 3, 1, 1};
   unsigned long int *stripoffsetarr = NULL;
   unsigned short int *stripbytecountarr=NULL;
   unsigned long int ifdfilepos = 0;
   int r, rr, s, rowsperstripvalue, stripsperimage;
   char blankchar = '\0';

   unsigned short int get_byte_order();

   /****************************************************************************
   * Get the unsigned short int code for the byte order of this machine.
   ****************************************************************************/
   tiffheader.byteorder = get_byte_order();

   /****************************************************************************
   * Open the file for writing.
   ****************************************************************************/
   if((fp = fopen(filename, "w+b")) == NULL){
      fprintf(stderr, "Error opening the file %s for writing!\n", filename);
      return(0);
   }

   /****************************************************************************
   * Write the TIFF header to the file. Note that we will need to go back
   * later to write in the position of the Image File Directory because we do
   * not yet know where it will be in the file.
   ****************************************************************************/
   fwrite(&tiffheader, sizeof(struct TIFFHEADER), 1, fp);

   /****************************************************************************
   * Calculate the rows per strip to get nearly 8 kilobytes per strip. Also
   * compute the number of strips we need for the image.
   ****************************************************************************/
   r = 0;
   do{
      rowsperstripvalue = (8 * 1024 - r)/ (2 * cols);
      if(rowsperstripvalue < 1) rowsperstripvalue = 1;
      stripsperimage = (rows + rowsperstripvalue - 1) / rowsperstripvalue;
      r++;
   }while(stripsperimage <= 2);

   /****************************************************************************
   * Fill in the image width value. Since it is less than 4 bytes it needs to
   * be left justified. If we are on a BigEndian machine shift the data 16bits.
   ****************************************************************************/
   imagewidth.valoffset = cols;
   if(tiffheader.byteorder == 0x4d4d)
      imagewidth.valoffset <<= 16;

   /****************************************************************************
   * Fill in the image length value. Since it is less than 4 bytes it needs to
   * be left justified. If we are on a BigEndian machine shift the data 16bits.
   ****************************************************************************/
   imagelength.valoffset = rows;
   if(tiffheader.byteorder == 0x4d4d)
      imagelength.valoffset <<= 16;

   /****************************************************************************
   * If we are on a BigEndian machine shift the data 16bits.
   ****************************************************************************/
   rowsperstrip.valoffset = rowsperstripvalue;
   if(tiffheader.byteorder == 0x4d4d){
      bitspersample.valoffset <<= 16;
      compression.valoffset <<= 16;
      photometric.valoffset <<= 16;
      planarconfiguration.valoffset <<= 16;
      sampleformat.valoffset <<= 16;
      rowsperstrip.valoffset <<= 16;
   }

   /****************************************************************************
   * Write the description to the file.
   ****************************************************************************/
   description.valoffset = ftell(fp);
   fwrite(descriptionname, sizeof(char), strlen(descriptionname), fp);
   description.length = strlen(descriptionname);
   if(description.length % 2){
      fputc((unsigned char)0, fp);
      description.length += 1;
   }
   else{
      fputc((unsigned char)0, fp);
      fputc((unsigned char)0, fp);
      description.length += 2;
   }

   /****************************************************************************
   * Write the name of the software and the version number to the file.
   ****************************************************************************/
   software.valoffset = ftell(fp);
   fwrite(softwarename, sizeof(char), strlen(softwarename), fp);
   software.length = strlen(softwarename);
   if(software.length % 2){
      fputc((unsigned char)0, fp);
      software.length += 1;
   }
   else{
      fputc((unsigned char)0, fp);
      fputc((unsigned char)0, fp);
      software.length += 2;
   }

   /****************************************************************************
   * Write the image data to the file and keep track of the strip offsets and
   * the length in bytes of each strip.
   ****************************************************************************/
   stripoffsets.length = stripsperimage;
   stripoffsetarr = (unsigned long int *) calloc(stripsperimage, sizeof(unsigned long int));
   stripbytecounts.length = stripsperimage;
   stripbytecountarr = (unsigned short int *) calloc(stripsperimage, sizeof(unsigned short int));
   for(s=0,r=0;s<stripsperimage;s++){ 
      stripoffsetarr[s] = ftell(fp);
      for(rr=0;rr<rowsperstripvalue;rr++,r++){
         if(r < rows){
            fwrite(&image[r*cols], sizeof(unsigned short int), cols, fp);
            stripbytecountarr[s] += sizeof(unsigned short int) * cols;
         }
      }
   }

   /****************************************************************************
   * Write the file offset for each strip.
   ****************************************************************************/
   stripoffsets.valoffset =  ftell(fp);
   fwrite(stripoffsetarr, sizeof(unsigned long int), stripsperimage, fp);

   /****************************************************************************
   * Write byte counts for each strip.
   ****************************************************************************/
   stripbytecounts.valoffset = ftell(fp);
   fwrite(stripbytecountarr, sizeof(unsigned short int), stripsperimage, fp);

   /****************************************************************************
   * Write the image file directory.
   ****************************************************************************/
   ifdfilepos = ftell(fp);
   fwrite(&(tiffifd.count), sizeof(unsigned short int), 1, fp);
   fwrite(&newsubfiletype, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&imagewidth, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&imagelength, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&bitspersample, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&compression, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&photometric, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&description, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&stripoffsets, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&rowsperstrip, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&stripbytecounts, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&planarconfiguration, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&software, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&sampleformat, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&(tiffifd.nextifd), sizeof(unsigned long int), 1, fp);

   /****************************************************************************
   * Go back and write the position of the image file directory in the header.
   ****************************************************************************/
   fseek(fp, 4, 0);
   fwrite(&ifdfilepos, sizeof(unsigned long int), 1, fp);

   free(stripoffsetarr);
   free(stripbytecountarr);
   fclose(fp);
   return(1);
}

/*******************************************************************************
* Function: get_byte_order
* Purpose: This function gets the byte order for the machine this code is
* running on. The value that is returned is the code used in the TIFF header
* to specify the byte order.
* Name: Michael Heath, University of South Florida
*       Dr. Melanie Sutton, University of West Florida
* Date: 1/25/2000
*******************************************************************************/
unsigned short int get_byte_order()
{
   unsigned short endian = 0x0001;
   unsigned char *littleEndian = (unsigned char *)&endian;

   if( *littleEndian ){
      /* printf( "\nIntel machine"); */
      return((unsigned short int)0x4949);
   }
   else{
      /* printf( "Not Intel machine"); */
      return((unsigned short int)0x4d4d);
   }
}


/*******************************************************************************
* Function: write_simple_tiff_CACHEIM
* Purpose: To write 16-bit image data to a TIFF file.
* Name: Michael Heath, University of South Florida
* Date: 1/25/2000
*******************************************************************************/
int write_simple_tiff_CACHEIM(char *filename, CACHEIM *image)
{
   FILE *fp=NULL;
   char *softwarename = "usfheath version 1.0.0";
   char *descriptionname = "This image was produced using software obtained from:\nhttp://marathon.csee.usf.edu/Mammography/software/heathusf.html";
   struct TIFFHEADER tiffheader = {0, 42, 0};
   struct TIFFIFD tiffifd = {13, 0, 0};
   struct TIFFTAG newsubfiletype = {254, 4, 1, 0};
   struct TIFFTAG imagewidth = {256, 3, 1, 0};
   struct TIFFTAG imagelength = {257, 3, 1, 0};
   struct TIFFTAG bitspersample = {258, 3, 1, 16};
   struct TIFFTAG compression = {259, 3, 1, 1};
   struct TIFFTAG photometric = {262, 3, 1, 1};
   struct TIFFTAG description = {270, 2, 0, 0};
   struct TIFFTAG stripoffsets = {273, 4, 0, 0};
   struct TIFFTAG rowsperstrip = {278, 3, 1, 1};
   struct TIFFTAG stripbytecounts = {279, 3, 1, 1};
   struct TIFFTAG planarconfiguration = {284, 3, 1, 1};
   struct TIFFTAG software = {305, 2, 0, 1};
   struct TIFFTAG sampleformat = {339, 3, 1, 1};
   unsigned long int *stripoffsetarr = NULL;
   unsigned short int *stripbytecountarr=NULL;
   unsigned char *thisrow_uchar=NULL;
   unsigned short int *thisrow_ushort=NULL;
   unsigned long int ifdfilepos = 0;
   int r, c, rr, s, rowsperstripvalue, stripsperimage, rows, cols;
   char blankchar = '\0';
   int bytesperpix = 1;

   unsigned short int get_byte_order();

   rows = image->rows;
   cols = image->cols;

   if(image->datatype == UCHARNUM) bytesperpix = 1;
   else if(image->datatype == USHORTNUM) bytesperpix = 2;
   else{
      fprintf(stderr, "Error in write_simple_tiff_CACHEIM()! Datatype is not UCHAR nor USHORTNUM.\n");
      exit(1);
   }
   bitspersample.valoffset = 8 * bytesperpix;

   /****************************************************************************
   * Get the unsigned short int code for the byte order of this machine.
   ****************************************************************************/
   tiffheader.byteorder = get_byte_order();

   /****************************************************************************
   * Open the file for writing.
   ****************************************************************************/
   if((fp = fopen(filename, "w+b")) == NULL){
      fprintf(stderr, "Error opening the file %s for writing!\n", filename);
      return(0);
   }

   /****************************************************************************
   * Write the TIFF header to the file. Note that we will need to go back
   * later to write in the position of the Image File Directory because we do
   * not yet know where it will be in the file.
   ****************************************************************************/
   fwrite(&tiffheader, sizeof(struct TIFFHEADER), 1, fp);

   /****************************************************************************
   * Calculate the rows per strip to get nearly 8 kilobytes per strip. Also
   * compute the number of strips we need for the image.
   ****************************************************************************/
   r = 0;
   do{
      rowsperstripvalue = (8 * 1024 - r)/ (cols * bytesperpix);
      if(rowsperstripvalue < 1) rowsperstripvalue = 1;
      stripsperimage = (rows + rowsperstripvalue - 1) / rowsperstripvalue;
      r++;
   }while(stripsperimage <= 2);

   /****************************************************************************
   * Fill in the image width value. Since it is less than 4 bytes it needs to
   * be left justified. If we are on a BigEndian machine shift the data 16bits.
   ****************************************************************************/
   imagewidth.valoffset = cols;
   if(tiffheader.byteorder == 0x4d4d)
      imagewidth.valoffset <<= 16;

   /****************************************************************************
   * Fill in the image length value. Since it is less than 4 bytes it needs to
   * be left justified. If we are on a BigEndian machine shift the data 16bits.
   ****************************************************************************/
   imagelength.valoffset = rows;
   if(tiffheader.byteorder == 0x4d4d)
      imagelength.valoffset <<= 16;

   /****************************************************************************
   * If we are on a BigEndian machine shift the data 16bits.
   ****************************************************************************/
   rowsperstrip.valoffset = rowsperstripvalue;
   if(tiffheader.byteorder == 0x4d4d){
      bitspersample.valoffset <<= 16;
      compression.valoffset <<= 16;
      photometric.valoffset <<= 16;
      planarconfiguration.valoffset <<= 16;
      sampleformat.valoffset <<= 16;
      rowsperstrip.valoffset <<= 16;
   }

   /****************************************************************************
   * Write the description to the file.
   ****************************************************************************/
   description.valoffset = ftell(fp);
   fwrite(descriptionname, sizeof(char), strlen(descriptionname), fp);
   description.length = strlen(descriptionname);
   if(description.length % 2){
      fputc((unsigned char)0, fp);
      description.length += 1;
   }
   else{
      fputc((unsigned char)0, fp);
      fputc((unsigned char)0, fp);
      description.length += 2;
   }

   /****************************************************************************
   * Write the name of the software and the version number to the file.
   ****************************************************************************/
   software.valoffset = ftell(fp);
   fwrite(softwarename, sizeof(char), strlen(softwarename), fp);
   software.length = strlen(softwarename);
   if(software.length % 2){
      fputc((unsigned char)0, fp);
      software.length += 1;
   }
   else{
      fputc((unsigned char)0, fp);
      fputc((unsigned char)0, fp);
      software.length += 2;
   }

   /****************************************************************************
   * Write the image data to the file and keep track of the strip offsets and
   * the length in bytes of each strip.
   ****************************************************************************/
   stripoffsets.length = stripsperimage;
   stripoffsetarr = (unsigned long int *) calloc(stripsperimage, sizeof(unsigned long int));
   stripbytecounts.length = stripsperimage;
   stripbytecountarr = (unsigned short int *) calloc(stripsperimage, sizeof(unsigned short int));

   if(bytesperpix == 1) thisrow_uchar = (unsigned char *) calloc(cols, sizeof(unsigned char));
   else if(bytesperpix == 2) thisrow_ushort = (unsigned short int *) calloc(cols, sizeof(unsigned short int));

   for(s=0,r=0;s<stripsperimage;s++){ 
      stripoffsetarr[s] = ftell(fp);
      for(rr=0;rr<rowsperstripvalue;rr++,r++){
         if(r < rows){
            if(bytesperpix == 1){
               for(c=0;c<cols;c++) thisrow_uchar[c] = (unsigned char)  (*(image->getpixel))(image, r, c);
               fwrite(thisrow_uchar, sizeof(unsigned char), cols, fp);
               stripbytecountarr[s] += sizeof(unsigned char) * cols;
            }
            else if(bytesperpix == 2){
               for(c=0;c<cols;c++) thisrow_ushort[c] = (unsigned short int)  (*(image->getpixel))(image, r, c);
               fwrite(thisrow_ushort, sizeof(unsigned short int), cols, fp);
               stripbytecountarr[s] += sizeof(unsigned short int) * cols;
            }
         }
      }
      if((bytesperpix == 1) && (stripbytecountarr[s] % 2)){
         fputc((unsigned char)0, fp);
         stripbytecountarr[s] += 1;
      }
   }
   if(bytesperpix == 1) free(thisrow_uchar);
   else if(bytesperpix == 2) free(thisrow_ushort);

   /****************************************************************************
   * Write the file offset for each strip.
   ****************************************************************************/
   stripoffsets.valoffset =  ftell(fp);
   fwrite(stripoffsetarr, sizeof(unsigned long int), stripsperimage, fp);

   /****************************************************************************
   * Write byte counts for each strip.
   ****************************************************************************/
   stripbytecounts.valoffset = ftell(fp);
   fwrite(stripbytecountarr, sizeof(unsigned short int), stripsperimage, fp);

   /****************************************************************************
   * Write the image file directory.
   ****************************************************************************/
   ifdfilepos = ftell(fp);
   fwrite(&(tiffifd.count), sizeof(unsigned short int), 1, fp);
   fwrite(&newsubfiletype, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&imagewidth, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&imagelength, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&bitspersample, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&compression, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&photometric, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&description, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&stripoffsets, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&rowsperstrip, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&stripbytecounts, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&planarconfiguration, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&software, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&sampleformat, sizeof(struct TIFFTAG), 1, fp);
   fwrite(&(tiffifd.nextifd), sizeof(unsigned long int), 1, fp);

   /****************************************************************************
   * Go back and write the position of the image file directory in the header.
   ****************************************************************************/
   fseek(fp, 4, 0);
   fwrite(&ifdfilepos, sizeof(unsigned long int), 1, fp);

   free(stripoffsetarr);
   free(stripbytecountarr);
   fclose(fp);
   return(1);
}
