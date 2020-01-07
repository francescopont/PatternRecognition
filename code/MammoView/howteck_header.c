/*******************************************************************************
* Program: This file includes functions for determining the vital header
* information (rows, cols, byte_order and header_bytes) from an image from the
* HowTeck scanner.
* Name: Mike Heath
* Date: 4/7/98
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>

/*******************************************************************************
* Function prototypes.
*******************************************************************************/
int read_howteck_image_header(char *filename,  int *rows, int *cols,
    int *headbytes, int *swapbytes);
short int chars_to_int(unsigned char a, unsigned char b);
short int chars_to_int_swapped(unsigned char a, unsigned char b);
short int read_short_int(FILE *fp, short int (*toint)(unsigned char, unsigned char));
short int find_rows(FILE *fp, short int (*toint)(unsigned char, unsigned char),
   short int (*read_sint)(FILE *, short int (*toint)(unsigned char, unsigned char)));
short int find_cols(FILE *fp, short int (*toint)(unsigned char, unsigned char),
   short int (*read_sint)(FILE *, short int (*toint)(unsigned char, unsigned char)));
long int find_headerbytes(FILE *fp, short int (*toint)(unsigned char, unsigned char),
   short int (*read_sint)(FILE *, short int (*toint)(unsigned char, unsigned char)));

/*
void main(int argc, char *argv[])
{
   char *filename=NULL;
   int rows=0, cols=0, headbytes=0, swapbytes=0;

   filename = argv[1];
   read_howteck_image_header(filename, &rows, &cols, &headbytes, &swapbytes);

   printf("rows = %d\n", rows);
   printf("cols = %d\n", cols);
   printf("headbytes = %d\n", headbytes);
   printf("swapbytes = %d\n", swapbytes);

}
*/

/*******************************************************************************
* Function: read_howteck_image_header
* Purpose: The HowTeck scanner used the DICOM image format. This is very
* complicated because tags are used so the header length is not fixed. In this
* function I only do 4 things. Therse are: 1) find the rows, 2) find the columns,
* 3) determine if the bytes need to be swapped, and 4) determine the number of
* feader bytes to skip over to get to the image data.
* Name: Mike Heath
* Date: 4/7/98
*******************************************************************************/
int read_howteck_image_header(char *filename,  int *rows, int *cols,
    int *headbytes, int *swapbytes)
{
   FILE *fp;
   short int intval1, intval2;
   short int (*toint)(unsigned char, unsigned char);

   /****************************************************************************
   * Open the file for reading.
   ****************************************************************************/
   if((fp = fopen(filename, "rb")) == NULL){
      fprintf(stderr, "Error opening the file %s for reading.\n", filename);
      return(0);
   }

   /****************************************************************************
   * Get the first two header bytes. They specify whether the file was written
   * in little endian or big endian format.
   ****************************************************************************/
   intval1 = fgetc(fp);
   intval2 = fgetc(fp);

   /* printf("%c%c\n", (char)intval1, (char)intval2); */

   intval1 = 'I';
   intval2 = 'I';

   if(((char)intval1 == 'I') && ((char)intval2 == 'I')){
      *swapbytes = 1;
      toint = chars_to_int_swapped;
   }
   else if(((char)intval1 == 'M') && ((char)intval2 == 'M')){
      *swapbytes = 0;
      toint = chars_to_int;
   }
   else{
      fprintf(stderr, "The file (%s) does not appear to be in HowTeck format.\n", filename);
      return(0);
   }

   /****************************************************************************
   * Get the next two characters from the file. These should make up the
   * number 42.
   ****************************************************************************/
/*
   intval1 = read_short_int(fp, toint);
   if(intval1 != 42){
      fprintf(stderr, "The magic number is %hd but it should be 42.\n", intval1);
      fprintf(stderr, "Error interpreting the image file.\n");
      return(0);
   }
*/

   /****************************************************************************
   * Read in the number of rows, column and the number of headerbytes from the
   * file.
   ****************************************************************************/
   *rows = find_rows(fp, toint, read_short_int);
   *cols = find_cols(fp, toint, read_short_int);
   *headbytes = find_headerbytes(fp, toint, read_short_int);

   fclose(fp);

   return(1);

}

/*******************************************************************************
* Function: find_rows
* Purpose: This function looks through the file for the tags indicating the
* rows in the image. It then reads in the number of rows and returns that
* value.
* Name: Mike Heath
* Date: 4/7/98
*******************************************************************************/
short int find_rows(FILE *fp, short int (*toint)(unsigned char, unsigned char),
   short int (*read_sint)(FILE *, short int (*toint)(unsigned char, unsigned char)))
{
   short int intval1, intval2;

   fseek(fp, 0, 0);

   intval1 = (*read_sint)(fp, toint);
   intval2 = (*read_sint)(fp, toint);

   while(!((intval1 == 40) && (intval2 == 16))){
      intval1 = intval2;
      intval2 = (*read_sint)(fp, toint);
   }

   /* printf("(%hx, %hx)\n", intval1, intval2); */

   intval1 = (*read_sint)(fp, toint);
   intval1 = (*read_sint)(fp, toint);
   intval1 = (*read_sint)(fp, toint);
   
   return(intval1);
}

/*******************************************************************************
* Function: find_cols
* Purpose: This function looks through the file for the tags indicating the
* columns in the image. It then reads in the number of columns and returns that
* value.
* Name: Mike Heath
* Date: 4/7/98
*******************************************************************************/
short int find_cols(FILE *fp, short int (*toint)(unsigned char, unsigned char),
   short int (*read_sint)(FILE *, short int (*toint)(unsigned char, unsigned char)))
{
   short int intval1, intval2;

   fseek(fp, 0, 0);

   intval1 = (*read_sint)(fp, toint);
   intval2 = (*read_sint)(fp, toint);

   while(!((intval1 == 40) && (intval2 == 17))){
      intval1 = intval2;
      intval2 = (*read_sint)(fp, toint);
   }

   /* printf("(%hx, %hx)\n", intval1, intval2); */

   intval1 = (*read_sint)(fp, toint);
   intval1 = (*read_sint)(fp, toint);
   intval1 = (*read_sint)(fp, toint);
   
   return(intval1);
}

/*******************************************************************************
* Function: find_headerbytes
* Purpose: This function determines the number of header bytes in the image.
* Name: Mike Heath
* Date: 4/7/98
*******************************************************************************/
long int find_headerbytes(FILE *fp, short int (*toint)(unsigned char, unsigned char),
   short int (*read_sint)(FILE *, short int (*toint)(unsigned char, unsigned char)))
{
   short int intval1, intval2;
   long int hb;

   fseek(fp, 0, 0);

   intval1 = (*read_sint)(fp, toint);
   intval2 = (*read_sint)(fp, toint);

   while(!((intval1 == 32736) && (intval2 == 16))){
      intval1 = intval2;
      intval2 = (*read_sint)(fp, toint);
      hb = ftell(fp);
      /* printf("%o (%x, %x)\n", hb, (long int)intval1, (long int)intval2); */
   }

   /* printf("(%x, %x)\n", (long int)intval1, (long int)intval2); */

   intval1 = (*read_sint)(fp, toint);
   intval1 = (*read_sint)(fp, toint);

   hb = ftell(fp);

   return(hb);
}

/*******************************************************************************
* Function: read_short_int
* Purpose: To read in a short integer from a file. The short int is swapped
* or not swapped depending on which function the function pointer toint
* points to.
* Name: Mike Heath
* Date: 4/7/98
*******************************************************************************/
short int read_short_int(FILE *fp, short int (*toint)(unsigned char, unsigned char))
{
   unsigned char intval1=0, intval2=0;
   short int intval;

   intval1 = fgetc(fp);
   intval2 = fgetc(fp);

   intval = (*toint)(intval1, intval2);

   /* printf("read_short() = %04x (%02x,%02x)\n", (long int)intval, (long int)intval1, (long int)intval2); */

   return(intval);
}

/*******************************************************************************
* Function: chars_to_int
* Purpose: To take two chars and form a short integer out of them.
* Name: Mike Heath
* Date: 4/7/98
*******************************************************************************/
short int chars_to_int(unsigned char a, unsigned char b)
{
   short int val = 0;

   val = a;
   val = val << 8;
   val += b;

   return(val);
}

/*******************************************************************************
* Function: chars_to_int
* Purpose: To take two chars and form a short integer out of them. The bytes are
* swapped in the process of forming the short integer.
* Name: Mike Heath
* Date: 4/7/98
*******************************************************************************/
short int chars_to_int_swapped(unsigned char a, unsigned char b)
{
   short int val = 0;

   val = b;
   val = val << 8;
   val += a;

   return(val);
}
