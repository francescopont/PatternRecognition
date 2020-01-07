/*******************************************************************************
* PROGRAM: MammoView
* PURPOSE: This program was written by Michael Heath at the University of
* South Florida. I developed this software as part of my dissertation
* research because I needed a program that would display images from the DDSM
* database without needing to reformat them. I also wanted enough
* control over the software to allow me to overlay graphic data such as the
* ground truth in the DDSM cases, breast segmentations, detections etc. on
* top of the images. All of the software was written in X without using any
* toolkits such as Motif. While this made the interface a little more crude,
* it did allow me to work on this software on other platforms, such as in
* Linux at home.
* 
* This program displays images stored in several formats. The program
* is named MammoView because a lot of digitized mammograms are stored in this
* format. Since these images are often very large, buffering is used for both
* the raw data stored in the images and intermediate images used in processing
* such as in aggregation.
*
* Presently the program can read images stored in:
* 1) compressed JPEG with an associated ".ics" file
* 2) Raw images with user specified image file parameters
* 3) both 8 and 16-bit PGM files
* 4) DBA 16-bit images (from Mass General Hospital)
* 5) LUMISYS images (from BOWMAN GRAY)
* 6) DICOM images from the HowTek scanner (From Mass General Hospital)
*
* Although the software was developed for my personal use, other students in
* the lab have found it useful for their research. For a couple of years I
* have made an executable version of the code available. Since I am finishing
* my degree, I wanted to make the source code available so anyone can pick it
* apart and modify it for their own purposes.
*
* I probably dont need to say this, but this software comes AS IS and has no
* warrenty of any kind. You may use it for research purposes. If you find the
* software useful, send me a short email. It will let me know that the time I
* spent to clean up the code and distribute it was not wasted. Good luck!
*
* NAME: Michael Heath, University of South Florida
* Copyright: Michael Heath and Dr. Kevin Bowyer 2000
* DATE: 6/21/2000
*
* NOTE: The header file xwindow.h must be in this directory for compiling the code.
*
* To Compile:
*Now there is a makefile for this program!
*******************************************************************************/
#include "virtual_image.h"
#include "xwindow.h"
#include "ICSio.h"
#include <string.h>
#include <math.h>
#include <limits.h>
#include "mikesfileio.h"
#include "optical_density.h"
#include "overlay.h"

static char codestring[20];

#define OverviewMaxRows (MaxScreenRows/2)
#define OverviewMaxCols (MaxScreenCols/2)

#define HIST_HEIGHT (MaxScreenRows/3)
#define HIST_WIDTH (MaxScreenCols/2)

#define OPTIONS_WIDTH 200

#define INTENSITY_MAPPING 0
#define DENSITY_MAPPING 1

#define PIXVALMODE 6
#define REGIONAVGMODE 7
#define ZOOMMODE 8

/*******************************************************************************
* Function prototypes.
*******************************************************************************/
void get_command_line_arguments(int argc, char *argv[], char *filename,
   int *rows, int *cols, int *bytesperpix, int *headerbytes, int *swap_bytes,
   int *fullresolution_height, int *fullresolution_width, char *scanner,
   int *startup_mode, char *view, char *direction, char **detection_filename,
   char **segmentation_filename,
   int *detmax, float *detthresh, char *lesiontype, char **ms, char **mm,
   char **ct, char **cd, char **pathology);

void setup_overview_image(CACHEIM *image, IMAGEXWIN *overview_image, int maxrows,
   int maxcols, unsigned short int **overview_data, char *window_title, int col_position);

void setup_fullresolution_window(CACHEIM *image, IMAGEXWIN *overview_image,
   IMAGEXWIN *fullresolution_image, IMAGEXWIN *histogram_image, int sub_rows,
   int sub_cols, USHORT **fullresolution_data, char *window_title);

void update_fullresolution_image(CACHEIM *image, IMAGEXWIN *overview_image,
   IMAGEXWIN *fullresolution_image, IMAGEXWIN *histogram_image, int x, int y,
   int scale_min, int scale_max, USHORT *fullresolution_data,
   int acquire_new_data, int num_detections, REGION *detections,
   int *isvalid, OVERLAY_DATA overlay_data);

void plot_histogram(unsigned long int *hist, int n, int lm, int rm, int tm,
   int bm, int *min_val, int *max_val, int *bin_width,
   IMAGEXWIN *histogram_image);

void rescale_overview_image(IMAGEXWIN *overview_image, unsigned short int *im,
   int num_detections, REGION *detections, int *isvalid, OVERLAY_DATA overlay_data);

void swap(unsigned char two_bytes[2]);

int overlay_polycurve(IMAGEXWIN *image, int rows_offset, int cols_offset,
   float subfactor, int *x_coord, int *y_coord, int numpoints, int color);

int overlay_polycurve_float(IMAGEXWIN *image, int rows_offset, int cols_offset,
   float subfactor, float *x_coord, float *y_coord, int numpoints, int color);

int readDBAImage_header(char *file_name, int *rows, int *cols,
   int *bytesperpixel, int *header_bytes, int *swapbytes, int *pixelskip,
   int *lineskip);

void aggregate(CACHEIM *image, CACHEIM *source_image, double aggfactor,
    int max_depth, int max_cache_size, int datatype, char *mode);

void aggregate_median(CACHEIM *image, CACHEIM *source_image, int aggfactor,
    int max_depth, int max_cache_size, int datatype, char *mode);

void setup_histogram_window(IMAGEXWIN *histogram_image, int hist_rows,
    int hist_cols, char *window_title);

void setup_buttons(IMAGEXWIN *imagemenu_image, int button_width, IMAGEXWIN *buttons,
   int num_buttons, char **button_labels, char *window_title, int row_position,
   int min_height, int min_width);

int overlay_groundtruth(IMAGEXWIN *image, int rows_offset, int cols_offset,
    float subfactor, int *isvalid, OVERLAY_DATA overlay_data);

char *find_which_button(XEvent theEvent, IMAGEXWIN *buttons, int num_buttons,
    char **button_labels);

void select_button(IMAGEXWIN *buttons, int num_buttons, char *button_labels[],
    char *alabel, unsigned long int color);

void unselect_button(IMAGEXWIN *buttons, int num_buttons, char *button_labels[],
    char *alabel, unsigned long int color);

void polygon_stats(CACHEIM *image, int *x_coord, int *y_coord, int num_points,
   int *count, double *mean, double *stdev);

int readIFSImage_header(char *file_name, int *rows, int *cols, int *bytesperpixel,
   int *header_bytes, int *swapbytes, int *pixelskip, int *lineskip);

void gaussian_smooth(unsigned short int *image, int rows, int cols, float sigma,
        float **smoothedim);

int read_howteck_image_header(char *filename,  int *rows, int *cols,
    int *headbytes, int *swapbytes);

int overlay_features(IMAGEXWIN *image, int rows_offset, int cols_offset,
   float subfactor);

int overlay_regions(IMAGEXWIN *image, int rows_offset, int cols_offset,
   float subfactor, int num_detections, REGION *detections);

void comensate_for_ahe(CACHEIM *image, CACHEIM *source_image, int max_depth,
   int max_cache_size, char *mode);

int overlay_nijmegen(IMAGEXWIN *image, int rows_offset, int cols_offset, float subfactor);

void update_zoom_window(CACHEIM *image, IMAGEXWIN *zoom_image, int x_buttonpos,
    int y_buttonpos, int zoom_factor);

/*******************************************************************************
* Some global variables.
*******************************************************************************/
XFontStruct *info_helvb14;
int verbosemode = 0;
Cursor arrow_cursor, crosshair_cursor, watch_cursor, question_cursor;

double od_offset = 0.0;
double od_scale = 1.0;
double min_density=0.0;
double max_density = 0.0;
double t_offset = 0.0;
double t_scale = 1.0;
double samplerate_in_microns = 0.0;

int current_aggfactor=1;
int scmin=0, scmax=65535;
int box_x, box_y, box_width, box_height;
int show_overlay = 0;
char created_files[4][200];
char overlay_filename[200];
float overview_subfactor=1.0;
int mapping = DENSITY_MAPPING;
int fullresolution_mode = PIXVALMODE;
int image_ul_x = 0, image_ul_y = 0;
int histogram_stretching=1;

int segmentation_toggle = 0;
int display_regions_toggle = 0;
int save_points_toggle = 0;

int test_x_coord[5] = {100, 100, 400, 300, 100};
int test_y_coord[5] = {100, 400, 500, 300, 50};
int test_numpoints = 5;

int last_min, last_max;

int VERBOSE = 0;

int seg_numpoints, seg_cols, seg_rows;
float seg_startx, seg_starty, seg_endx, seg_endy, *seg_xcoord=NULL, *seg_ycoord=NULL;
int other_seg_numpoints, other_seg_cols, other_seg_rows;
float other_seg_startx, other_seg_starty, other_seg_endx, other_seg_endy,
    *other_seg_xcoord=NULL, *other_seg_ycoord=NULL;
int other_num_detections;
REGION *other_detections;
float base_subfactor;

char *MIASgt_filename = NULL;

int ROWS;
int bothdet = 0;

main(int argc, char *argv[])
{
   CACHEIM image, image_raw, image_orig, image_full, image_agg2, image_agg4, image_agg8;
   char filename[200];
   int i, j, rows, cols, bytesperpix, headerbytes, pos;
   int swap_bytes=FALSE;
   char fontname[] = "-adobe-courier-medium-r-normal--14-140-75-75-m-90-iso8859-1";

   IMAGEXWIN overview_image, fullresolution_image, histogram_image, imagemenu_image, zoom_image;
   int subimage_rows=0, subimage_cols=0;
   int overview_x=0, overview_y=0;
   int scale_min, scale_max;
   USHORT *overview_data=NULL, *fullresolution_data=NULL;
   XEvent theEvent, theFutureEvent;
   char scanner[20] = {'\0'};
   char save_filename[140];
   char new_save_filename[140];
   FILE *fpout=NULL;
   XImage *tmpximage=NULL;
   char string1[100], string2[100], string3[100], string4[100];
   int cx, cy;
   int twidth=0;
   int basepos;
   int r, c;
   float *temp_float = NULL;
   long temp_long;
   int imagemenu_toggle = 0, histogram_toggle = 0;
   char *which_button_label = NULL;
   int overview_toggle = 0;
   int startup_mode = 0;
   char view[20] = {'\0'}, direction[20] = {'\0'};
   int zoom_factor = 1;
   int last_zoompos_x, last_zoompos_y;
   char *detection_filename = NULL;
   char *segmentation_filename = NULL;
   int detmax = 0;
   float detthresh = 0.0;
   int detection_rows, detection_cols, num_detections;
   int other_detection_rows, other_detection_cols, other_num_detections;
   REGION *detections=NULL;
   OVERLAY_DATA overlay_data;
   int *isvalid = NULL;
   char lesiontype[20];
   int status;
   char *ms = NULL, *mm = NULL, *ct = NULL, *cd = NULL;
   char *pathology = NULL;


   int zoommenu_toggle = 0;
   int num_zoombuttons=5;
   IMAGEXWIN zoommenu_image;
   IMAGEXWIN zoombuttons[5];
   char *zoombutton_labels[5] = {"5:1 ZOOM", "4:1 ZOOM", "3:1 ZOOM", "2:1 ZOOM","1:1 ZOOM"};

   int num_buttons=16;
   IMAGEXWIN buttons[16];
   char *button_labels[16];
   char label_0[] = "FULLRES";
   char label_1[] = "2X AGG";
   char label_2[] = "4X AGG";
   char label_3[] = "8X AGG";
   char label_4a[] = "NO OVERLAY";
   char label_4b[] = "OVERLAY";
   char label_5[] = "SAVE";
   char label_6[] = "SAVE SUBAREA";
   char label_7[] = "PIX VALUE";
   char label_8[] = "REGION AVG";
   char label_9[] = "DISTANCE";
   char label_10[] = "SHARPEN";
   char label_11[] = "OVERVIEW";
   char label_12[] = "HISTOGRAM";
   char label_13[] = "ZOOM";
   char label_14[] = "SEGMENTATION";
   char label_15[] = "QUIT";

   button_labels[0] = label_0;
   button_labels[1] = label_1;
   button_labels[2] = label_2;
   button_labels[3] = label_3;
   button_labels[4] = label_4a;
   button_labels[5] = label_5;
   button_labels[6] = label_6;
   button_labels[7] = label_7;
   button_labels[8] = label_8;
   button_labels[9] = label_9;
   button_labels[10] = label_10;
   button_labels[11] = label_11;
   button_labels[12] = label_12;
   button_labels[13] = label_13;
   button_labels[14] = label_14;
   button_labels[15] = label_15;

   overlay_filename[0] = '\0';

   /****************************************************************************
   * Zero our CACHEIM data structures.
   ****************************************************************************/
   memset((void *)&image, 0, sizeof(CACHEIM));
   memset((void *)&image_raw, 0, sizeof(CACHEIM));
   memset((void *)&image_full, 0, sizeof(CACHEIM));
   memset((void *)&image_agg2, 0, sizeof(CACHEIM));
   memset((void *)&image_agg4, 0, sizeof(CACHEIM));
   memset((void *)&image_agg8, 0, sizeof(CACHEIM));

   /****************************************************************************
   * Initialize the windows to tell the datastructures that they don't have
   * windows yet. Pixmaps will be allocated in display_image.
   ****************************************************************************/
   overview_image.have_a_pixmap = 0;
   fullresolution_image.have_a_pixmap = 0;
   histogram_image.have_a_pixmap = 0;
   imagemenu_image.have_a_pixmap = 0;
   zoommenu_image.have_a_pixmap = 0;
   zoom_image.have_a_pixmap = 0;

   /****************************************************************************
   * Get the command line arguments.
   ****************************************************************************/
   get_command_line_arguments(argc, argv, filename, &rows, &cols, &bytesperpix,
      &headerbytes, &swap_bytes, &subimage_rows, &subimage_cols, scanner,
      &startup_mode, view, direction, &detection_filename,
      &segmentation_filename, &detmax, &detthresh,
      lesiontype, &ms, &mm, &ct, &cd, &pathology);
   ROWS = rows;

   /****************************************************************************
   * In Linux, we would want to uncomment this line.
   ****************************************************************************/
   /* swap_bytes = !swap_bytes; */

   if(verbosemode){
      printf("****************************************");
      printf("****************************************\n");
      printf("The image is %d rows x %d cols\n", rows, cols);
      printf("The image has %d header bytes and %d bytes per pixel.\n",
         headerbytes, bytesperpix);

      if(swap_bytes == TRUE) printf("The bytes will be swapped.\n");
      else printf("The bytes will not be swapped.\n");

      printf("The sampling rate is %lf microns.\n", samplerate_in_microns);
      if(mapping == DENSITY_MAPPING) printf("You are working in density mode.\n");
      if(mapping == INTENSITY_MAPPING) printf("You are working in intensity mode.\n");
      printf("The scanner is: %s\n", scanner);
      printf("****************************************");
      printf("****************************************\n");
   }

   /****************************************************************************
   * If there is an overlay file then change the button label "NO OVERLAY" to
   * "OVERLAY".
   ****************************************************************************/
   if(overlay_filename[0] != '\0') button_labels[4] = label_4b;

   /****************************************************************************
   * If there is a MIAS ground truth filename then change the button label
   * "NO OVERLAY" to "OVERLAY".
   ****************************************************************************/
   if(MIASgt_filename != NULL) button_labels[4] = label_4b;

   /****************************************************************************
   * Read in the segmentation file.
   ****************************************************************************/
   if(segmentation_filename != NULL){
      char tmpfilename[200] = {'\0'}, *start=NULL;

      read_segmentation_file(segmentation_filename, &seg_numpoints, &seg_xcoord, &seg_ycoord,
         &seg_startx, &seg_starty, &seg_endx, &seg_endy, &seg_rows, &seg_cols);
      convert_segmentation_coordinates(seg_rows, seg_cols,
         seg_numpoints, seg_xcoord, seg_ycoord, &seg_startx, &seg_starty,
         &seg_endx, &seg_endy, rows, cols);

      /*************************************************************************
      * Read in the segmentation from the file for the other projection. Assume
      * The filename is the same except "CC" is replaced with "MLO" or "MLO" is
      * replaced with CC.
      *************************************************************************/
      if(bothdet == 1){
         if(strstr(segmentation_filename, "MLO") != NULL){
            strcpy(tmpfilename, segmentation_filename);

            start = strstr(tmpfilename, "MLO");
            strcpy(start+2, strstr(segmentation_filename, "MLO")+3);
            *start = 'C';
            *(start+1) = 'C';

            read_segmentation_file(tmpfilename, &other_seg_numpoints,
               &other_seg_xcoord, &other_seg_ycoord, &other_seg_startx, &other_seg_starty,
               &other_seg_endx, &other_seg_endy, &other_seg_rows, &other_seg_cols);

         }
         else if(strstr(segmentation_filename, "CC") != NULL){
            strcpy(tmpfilename, segmentation_filename);

            start = strstr(tmpfilename, "CC");
            strcpy(start+3, strstr(segmentation_filename, "CC")+2);
            *start = 'M';
            *(start+1) = 'L';
            *(start+2) = 'O';

            read_segmentation_file(tmpfilename, &other_seg_numpoints,
               &other_seg_xcoord, &other_seg_ycoord,
               &other_seg_startx, &other_seg_starty, &other_seg_endx, &other_seg_endy,
               &other_seg_rows, &other_seg_cols);
         }
      }
   }

   /****************************************************************************
   * Read in the detection file.
   ****************************************************************************/
   if(detection_filename != NULL){
      char tmpfilename[200] = {'\0'}, *start=NULL;

      read_detection_file(detection_filename, &num_detections, &detections,
               &detection_rows, &detection_cols, detmax, detthresh);
      convert_detection_coordinates(detection_rows, detection_cols, num_detections,
         detections, rows, cols);
      display_regions_toggle = 1;

      if(bothdet){
         /*************************************************************************
         * Read in the detections from the file for the other projection. Assume
         * The filename is the same except "CC" is replaced with "MLO" or "MLO" is
         * replaced with CC.
         *************************************************************************/
         if(strstr(detection_filename, "MLO") != NULL){
            strcpy(tmpfilename, detection_filename);

            start = strstr(tmpfilename, "MLO");
            strcpy(start+2, strstr(detection_filename, "MLO")+3);
            *start = 'C';
            *(start+1) = 'C';

            read_detection_file(tmpfilename, &other_num_detections, &other_detections,
                     &other_detection_rows, &other_detection_cols, detmax, detthresh);

            convert_detection_coordinates(other_detection_rows, other_detection_cols,
                     other_num_detections, other_detections, other_seg_rows,
                     other_seg_cols);

         }
         else if(strstr(detection_filename, "CC") != NULL){
            strcpy(tmpfilename, detection_filename);

            start = strstr(tmpfilename, "CC");
            strcpy(start+3, strstr(detection_filename, "CC")+2);
            *start = 'M';
            *(start+1) = 'L';
            *(start+2) = 'O';

            read_detection_file(tmpfilename, &other_num_detections, &other_detections,
                     &other_detection_rows, &other_detection_cols, detmax, detthresh);

            convert_detection_coordinates(other_detection_rows, other_detection_cols,
                     other_num_detections, other_detections, other_seg_rows,
                     other_seg_cols);
         }
      }
   }

   /****************************************************************************
   * Read in the overlay file.
   ****************************************************************************/
   if(overlay_filename[0] != '\0'){

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
   }

   /****************************************************************************
   * If a MIAS gt filename was specified, make the first character of the
   * overlay_filename something other than '\0'. This is simply done to fake
   * out the code into thinking there is an overlay file. This was necessary to
   * easily add the ability of this program to function with MIAS ground truth
   * files.
   ****************************************************************************/
   if(MIASgt_filename != NULL){
      overlay_filename[0] = 'A';
      overlay_filename[1] = '\0';
   }

   /****************************************************************************
   * Initialize the system for windows. These functions are found in the
   * xwindow.h header file. This was done for simplicity because there are
   * some global variables used in this program that are also used in the
   * xwindows.h file.
   ****************************************************************************/
   initX();
   initGrayScalecolormap();

   /****************************************************************************
   * Load a font.
   ****************************************************************************/
   if((info_helvb14 = XLoadQueryFont(theDisplay, fontname)) == NULL){
      fprintf(stderr, "Cannot load font: %s\n", fontname);
      exit(1);
   }

   /****************************************************************************
   * Set up a couple of cursors.
   ****************************************************************************/
   arrow_cursor = XCreateFontCursor(theDisplay, XC_top_left_arrow);
   crosshair_cursor = XCreateFontCursor(theDisplay, XC_crosshair);
   watch_cursor = XCreateFontCursor(theDisplay, XC_watch);
   question_cursor = XCreateFontCursor(theDisplay, XC_question_arrow);

   /****************************************************************************
   * Open an image that is on disk already.
   ****************************************************************************/
   if(bytesperpix == 1){
      if(allocate_cached_image(&image_raw, rows, cols, UCHARNUM, filename, "rb",
         12, (rows*cols*sizeof(UCHAR)/16), headerbytes, (USHORT)swap_bytes, samplerate_in_microns) == 0)
         exit(1);
   }
   else if(bytesperpix == 2){
      if(allocate_cached_image(&image_raw, rows, cols, USHORTNUM, filename, "rb",
         12, (rows*cols*sizeof(USHORT)/16), headerbytes, (USHORT)swap_bytes, samplerate_in_microns) == 0)
         exit(1);
   }

   image_orig = image_raw;

   /****************************************************************************
   * Add the processing step of mapping the image to optical density this is
   * required.
   *
   *   MGH_DBA_optical_density = 4.26700423014133 + (-0.90303289757264) * log10(pixel_value)
   *   MGH_HOWTEK_optical_density = 3.78928997845071 + (-0.00094568009377) * (pixel_value)
   *   LUMISYS_optical_density = 4.05977749300340 + (-0.00099080941710) * (pixel_value)
   *   ISMD_HOWTEK_optical_density = 3.96604095240593 + (-0.00099055807612) * (pixel_value)
   ****************************************************************************/
   if(mapping == DENSITY_MAPPING){
      if(strcmp(scanner, "DBA") == 0){
         if(verbosemode) printf("The image was scanned with the DBA scanner.\n");
         log10_optical_density(&image_full, &image_orig, 12,
            rows*cols*sizeof(USHORT)/8, "rb", &od_offset, &od_scale,
            &min_density, &max_density, 4.26700423014133, -0.90303289757264);
      }
      else if(strcmp(scanner, "LUMISYS") == 0){
         if(verbosemode) printf("The image was scanned with the LUMISYS scanner.\n");
         linear_optical_density(&image_full, &image_orig, 12,
            rows*cols*sizeof(USHORT)/8, "rb", &od_offset, &od_scale,
            &min_density, &max_density, 4.05977749300340, -0.00099080941710);
      }
      else if(strcmp(scanner, "HOWTEK") == 0){
         if(verbosemode) printf("The image was scanned with the HOWTEK scanner.\n");
         linear_optical_density(&image_full, &image_orig, 12,
            rows*cols*sizeof(USHORT)/8, "rb", &od_offset, &od_scale,
            &min_density, &max_density, 3.78928997845071, -0.00094568009377);
      }
      else if(strcmp(scanner, "HOWTEK_ISMD") == 0){
         if(verbosemode) printf("The image was scanned with the HOWTEK_ISMD scanner.\n");
         linear_optical_density(&image_full, &image_orig, 12,
            rows*cols*sizeof(USHORT)/8, "rb", &od_offset, &od_scale,
            &min_density, &max_density, 3.96604095240593, -0.00099055807612);
      }
      else{
         mapping = INTENSITY_MAPPING;
         image_full = image_orig;
      }
   }
   else image_full = image_orig;

   /****************************************************************************
   * Setup the buttons used in the image menu.
   ****************************************************************************/
   setup_buttons(&imagemenu_image, OPTIONS_WIDTH, buttons,
      num_buttons, button_labels, "MAIN MENU", 0, 0, 0);

   select_button(buttons, num_buttons, button_labels, "FULLRES", colors[BLUEINDEX].pixel);
   if((overlay_filename[0] != '\0') && (show_overlay == 1))
      select_button(buttons, num_buttons, button_labels, "OVERLAY", colors[BLUEINDEX].pixel);

   select_button(buttons, num_buttons, button_labels, "PIX VALUE", colors[BLUEINDEX].pixel);

   /****************************************************************************
   * Setup the buttons used in the zoom menu.
   ****************************************************************************/
   setup_buttons(&zoommenu_image, OPTIONS_WIDTH, zoombuttons,
      num_zoombuttons, zoombutton_labels, "ZOOM", imagemenu_image.rows+34+34, 512, OPTIONS_WIDTH+512);

   /****************************************************************************
   * Setup the buttons used in the zoom menu.
   ****************************************************************************/
   zoom_image = zoommenu_image;
   zoom_image.have_a_pixmap = 0;
   zoom_image.theDrawWindow = openWindow(OPTIONS_WIDTH, 0,
      512, 512, NORMAL_WINDOW, "ZOOM IMAGE", zoom_image.theBorderWindow,
      &(zoom_image.theGC));
   initEvents(zoom_image.theDrawWindow, NOT_IN_PALETTE);

   zoom_image.rows = 512;
   zoom_image.cols = 512;

   if((zoom_image.image = (unsigned char *) calloc(512*512, sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Malloc error.\n");
   }
   if((zoom_image.image_dsp = (unsigned char *) calloc(512*512, sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Malloc error.\n");
   }
   zoom_image.maskimage = NULL;
   zoom_image.theImage = NULL;

   display_image(&(zoom_image), 512, 512, 1);

   /****************************************************************************
   * Initially set up the display to be in full resolution mode.
   ****************************************************************************/
   current_aggfactor=1;
   image = image_full;

   /****************************************************************************
   * Set up the overview image in a window.
   ****************************************************************************/
   setup_overview_image(&image, &overview_image, OverviewMaxRows,
      OverviewMaxCols, &overview_data, filename, imagemenu_image.cols+10);

   overview_subfactor = (float)image.rows / (float)overview_image.rows;

   /****************************************************************************
   * If we are supposed to display the chaincodes, place them on the
   * overview image.
   ****************************************************************************/
   if((overlay_filename[0] != '\0') && (show_overlay == 1))
      overlay_groundtruth(&overview_image, 0, 0, overview_subfactor,
         isvalid, overlay_data);

   XDefineCursor(theDisplay, overview_image.theDrawWindow, watch_cursor);

   /****************************************************************************
   * Set up the fullresolution image in a window. The user could have specified
   * the size of the display window, but if they didn't, we calculate a size
   * to use.
   ****************************************************************************/
   if(subimage_rows == 0) subimage_rows = MaxScreenRows - 34;
   if(subimage_cols == 0) subimage_cols = MaxScreenCols - overview_image.cols - 18;

   setup_histogram_window(&histogram_image, HIST_HEIGHT, HIST_WIDTH, "Histogram Window");

   setup_fullresolution_window(&image, &overview_image, &fullresolution_image,
      &histogram_image, subimage_rows, subimage_cols, &fullresolution_data, filename);
   XSetLineAttributes(theDisplay, fullresolution_image.theGC, 2, LineSolid, CapRound, JoinMiter);
   XFlush(theDisplay);

   update_fullresolution_image(&image, &overview_image, &fullresolution_image,
      &histogram_image, 0, 0, -1, -1, fullresolution_data, 1, num_detections, detections,
      isvalid, overlay_data);

   XDefineCursor(theDisplay, overview_image.theDrawWindow, crosshair_cursor);
   XDefineCursor(theDisplay, fullresolution_image.theDrawWindow, crosshair_cursor);
   XFlush(theDisplay);

   /****************************************************************************
   * Map all of the windows so we can see them.
   ****************************************************************************/
   XMapWindow(theDisplay, fullresolution_image.theBorderWindow);
   XMapWindow(theDisplay, fullresolution_image.theDrawWindow);

   if(startup_mode == 1){

      XMapWindow(theDisplay, imagemenu_image.theBorderWindow);
      for(i=0;i<num_buttons;i++) XMapWindow(theDisplay, buttons[i].theDrawWindow);
      imagemenu_toggle = 1;

      XMapWindow(theDisplay, overview_image.theBorderWindow);
      XMapWindow(theDisplay, overview_image.theDrawWindow);
      select_button(buttons, num_buttons, button_labels, "OVERVIEW", colors[BLUEINDEX].pixel);
      overview_toggle = 1;

      XMapWindow(theDisplay, histogram_image.theBorderWindow);
      XMapWindow(theDisplay, histogram_image.theDrawWindow);
      select_button(buttons, num_buttons, button_labels, "HISTOGRAM", colors[BLUEINDEX].pixel);
      histogram_toggle = 1;
   }

   XFlush(theDisplay);

   /****************************************************************************
   * Enter the EventLoop. This is where all of the interaction is done until
   * the program is terminated. This is a really ugly mess.
   ****************************************************************************/
   do{

      BEGIN_LOOP:
 
      XNextEvent(theDisplay, &theEvent);
 
      /*************************************************************************
      * Handle expose events by repainting the screen.
      *************************************************************************/
      if(theEvent.type == Expose){
         if(theEvent.xbutton.window == overview_image.theDrawWindow)
            repaint(overview_image, theEvent.xexpose.x, theEvent.xexpose.y,
                    theEvent.xexpose.width, theEvent.xexpose.height,
                    theEvent.xexpose.x, theEvent.xexpose.y);
         if(theEvent.xbutton.window == fullresolution_image.theDrawWindow)
            repaint(fullresolution_image, theEvent.xexpose.x, theEvent.xexpose.y,
                    theEvent.xexpose.width, theEvent.xexpose.height,
                    theEvent.xexpose.x, theEvent.xexpose.y);
         if(theEvent.xbutton.window == histogram_image.theDrawWindow)
            repaint(histogram_image, theEvent.xexpose.x, theEvent.xexpose.y,
                    theEvent.xexpose.width, theEvent.xexpose.height,
                    theEvent.xexpose.x, theEvent.xexpose.y);
         if(theEvent.xbutton.window == zoom_image.theDrawWindow)
            repaint(zoom_image, theEvent.xexpose.x, theEvent.xexpose.y,
                    theEvent.xexpose.width, theEvent.xexpose.height,
                    theEvent.xexpose.x, theEvent.xexpose.y);
         for(i=0;i<num_buttons;i++){
            if(theEvent.xbutton.window == buttons[i].theDrawWindow)
               repaint(buttons[i], theEvent.xexpose.x, theEvent.xexpose.y,
                    theEvent.xexpose.width, theEvent.xexpose.height,
                    theEvent.xexpose.x, theEvent.xexpose.y);
         }
         for(i=0;i<num_zoombuttons;i++){
            if(theEvent.xbutton.window == zoombuttons[i].theDrawWindow)
               repaint(zoombuttons[i], theEvent.xexpose.x, theEvent.xexpose.y,
                    theEvent.xexpose.width, theEvent.xexpose.height,
                    theEvent.xexpose.x, theEvent.xexpose.y);
         }
      }

      /*************************************************************************
      * When we have a button press in a window, we handle it in this loop.
      * A button press on the quit button is not processed here.
      *************************************************************************/
      else if(theEvent.type == ButtonPress){

         /**********************************************************************
         * If the buttopress was in a "button" then process it here.
         **********************************************************************/
         if(((which_button_label = find_which_button(theEvent, buttons, num_buttons, button_labels)) != NULL) ||
            ((which_button_label = find_which_button(theEvent, zoombuttons, num_zoombuttons, zoombutton_labels)) != NULL)){

            XDefineCursor(theDisplay, overview_image.theDrawWindow, watch_cursor);
            XDefineCursor(theDisplay, fullresolution_image.theDrawWindow, watch_cursor);
            XDefineCursor(theDisplay, histogram_image.theDrawWindow, watch_cursor);
            XFlush(theDisplay);

            /*******************************************************************
            * This is the "FULLRES" button.
            *******************************************************************/
            if(strcmp("FULLRES", which_button_label) == 0){

               select_button(buttons, num_buttons, button_labels, "FULLRES", colors[REDINDEX].pixel);
               unselect_button(buttons, num_buttons, button_labels, "2X AGG", cell[0]);
               unselect_button(buttons, num_buttons, button_labels, "4X AGG", cell[0]);
               unselect_button(buttons, num_buttons, button_labels, "8X AGG", cell[0]);

               overview_x = (int)(((overview_x * (float)(overview_image.cols) /
                            (float)(image.cols)) * (float)image_full.cols) /
                            (float)(overview_image.cols));
               overview_y = (int)(((overview_y * (float)(overview_image.rows) /
                            (float)(image.rows)) * (float)image_full.rows) /
                            (float)(overview_image.rows));

               image = image_full;
               current_aggfactor = 1;

               scale_min = -1;
	       scale_max = -1;

               update_fullresolution_image(&image, &overview_image,
                  &fullresolution_image, &histogram_image, overview_x, overview_y,
                  scale_min, scale_max, fullresolution_data, 1, num_detections, detections,
                  isvalid, overlay_data);

               select_button(buttons, num_buttons, button_labels, "FULLRES", colors[BLUEINDEX].pixel);
            }
            /*******************************************************************
            * This is the "2X AGG" button.
            *******************************************************************/
            else if(strcmp("2X AGG", which_button_label) == 0){

               unselect_button(buttons, num_buttons, button_labels, "FULLRES", cell[0]);
               select_button(buttons, num_buttons, button_labels, "2X AGG", colors[REDINDEX].pixel);
               unselect_button(buttons, num_buttons, button_labels, "4X AGG", cell[0]);
               unselect_button(buttons, num_buttons, button_labels, "8X AGG", cell[0]);

               if(image_agg2.depth == 0){
                  aggregate(&image_agg2, &image_full, 2, 12, rows*cols*sizeof(USHORT)/16, USHORTNUM, "rb");
               }

	       if(verbosemode) printf("The 2X aggregated image resolution is %f microns.\n", image_agg2.resolution);

               overview_x = (int)(((overview_x * (float)(overview_image.cols) /
                            (float)(image.cols)) * (float)image_agg2.cols) /
                            (float)(overview_image.cols));
               overview_y = (int)(((overview_y * (float)(overview_image.rows) /
                            (float)(image.rows)) * (float)image_agg2.rows) /
                            (float)(overview_image.rows));

               image = image_agg2;
               current_aggfactor = 2;

               scale_min = -1;
	       scale_max = -1;

               update_fullresolution_image(&image, &overview_image,
                  &fullresolution_image, &histogram_image, overview_x, overview_y,
                  scale_min, scale_max, fullresolution_data, 1, num_detections, detections,
                  isvalid, overlay_data);

               select_button(buttons, num_buttons, button_labels, "2X AGG", colors[BLUEINDEX].pixel);
            }
            /*******************************************************************
            * This is the "4X AGG" button.
            *******************************************************************/
            else if(strcmp("4X AGG", which_button_label) == 0){

               unselect_button(buttons, num_buttons, button_labels, "FULLRES", cell[0]);
               unselect_button(buttons, num_buttons, button_labels, "2X AGG", cell[0]);
               select_button(buttons, num_buttons, button_labels, "4X AGG", colors[REDINDEX].pixel);
               unselect_button(buttons, num_buttons, button_labels, "8X AGG", cell[0]);

               if(image_agg4.depth == 0){
                  aggregate(&image_agg4, &image_full, 4, 12, rows*cols*sizeof(USHORT)/16, USHORTNUM, "rb");
               }

	       if(verbosemode) printf("The 4X aggregated image resolution is %f microns.\n", image_agg4.resolution);

               overview_x = (int)(((overview_x * (float)(overview_image.cols) /
                            (float)(image.cols)) * (float)image_agg4.cols) /
                            (float)(overview_image.cols));
               overview_y = (int)(((overview_y * (float)(overview_image.rows) /
                            (float)(image.rows)) * (float)image_agg4.rows) /
                            (float)(overview_image.rows));

               image = image_agg4;
               current_aggfactor = 4;

               scale_min = -1;
	       scale_max = -1;

               update_fullresolution_image(&image, &overview_image,
                  &fullresolution_image, &histogram_image, overview_x, overview_y,
                  scale_min, scale_max, fullresolution_data, 1, num_detections, detections,
                  isvalid, overlay_data);

               select_button(buttons, num_buttons, button_labels, "4X AGG", colors[BLUEINDEX].pixel);
            }
            /*******************************************************************
            * This is the "8X AGG" button.
            *******************************************************************/
            else if(strcmp("8X AGG", which_button_label) == 0){

               unselect_button(buttons, num_buttons, button_labels, "FULLRES", cell[0]);
               unselect_button(buttons, num_buttons, button_labels, "2X AGG", cell[0]);
               unselect_button(buttons, num_buttons, button_labels, "4X AGG", cell[0]);
               select_button(buttons, num_buttons, button_labels, "8X AGG", colors[REDINDEX].pixel);

               if(image_agg8.depth == 0){
                  aggregate(&image_agg8, &image_full, 8, 12, rows*cols*sizeof(USHORT)/60, USHORTNUM, "rb");
               }

	       if(verbosemode) printf("The 8X aggregated image resolution is %f microns.\n", image_agg8.resolution);

               overview_x = (int)(((overview_x * (float)(overview_image.cols) /
                            (float)(image.cols)) * (float)image_agg4.cols) /
                            (float)(overview_image.cols));
               overview_y = (int)(((overview_y * (float)(overview_image.rows) /
                            (float)(image.rows)) * (float)image_agg4.rows) /
                            (float)(overview_image.rows));

               image = image_agg8;
               current_aggfactor = 8;

               scale_min = -1;
	       scale_max = -1;

               update_fullresolution_image(&image, &overview_image,
                  &fullresolution_image, &histogram_image, overview_x, overview_y,
                  scale_min, scale_max, fullresolution_data, 1, num_detections, detections,
                  isvalid, overlay_data);

               select_button(buttons, num_buttons, button_labels, "8X AGG", colors[BLUEINDEX].pixel);
            }

            /*******************************************************************
            * This is the "OVERLAY" or "NO OVERLAY" button.
            *******************************************************************/
            else if((strcmp("OVERLAY", which_button_label) == 0) ||
                   (strcmp("NO OVERLAY", which_button_label) == 0)){

               if(overlay_filename[0] != '\0'){

                  if(show_overlay == 1){    /* Overlay is on, so turn it off. */
                     show_overlay = 0;
                  }
                  else{                     /* Overlay is off, so turn it on. */
                     show_overlay = 1;
                  }

                  select_button(buttons, num_buttons, button_labels, "OVERLAY", colors[REDINDEX].pixel);

                  scale_min = -1;
	          scale_max = -1;

                  update_fullresolution_image(&image, &overview_image,
                     &fullresolution_image, &histogram_image, overview_x, overview_y,
                     scale_min, scale_max, fullresolution_data, 0, num_detections, detections,
                     isvalid, overlay_data);

                  /*************************************************************
                  * Repaint the crop box on the positioning (overview) image.
                  *************************************************************/
                  display_image(&overview_image, overview_image.rows,
                     overview_image.cols, 1);

                  XSetForeground(theDisplay, overview_image.theGC,
                     colors[BLUEINDEX].pixel);  /* Blue */
            
                  XDrawRectangle(theDisplay, overview_image.thePixmap,
                     overview_image.theGC, box_x, box_y, box_width, box_height);
            
                  if((overlay_filename[0] != '\0') && (show_overlay == 1)){
                     overlay_groundtruth(&overview_image, 0, 0,
                        overview_subfactor, isvalid, overlay_data);
                  }

		  if(segmentation_toggle)
                     overlay_features(&overview_image, 0, 0, overview_subfactor);

                  if(display_regions_toggle)
                     overlay_regions(&overview_image, 0, 0, overview_subfactor,
                        num_detections, detections);


                  repaint(overview_image, 0, 0, overview_image.cols,
                     overview_image.rows, 0, 0);

                  if(show_overlay == 1){
                     select_button(buttons, num_buttons, button_labels, "OVERLAY", colors[BLUEINDEX].pixel);
                  }
                  else{
                     unselect_button(buttons, num_buttons, button_labels, "OVERLAY", cell[0]);
                  }
               }
            }

            /*******************************************************************
            * This is the "SAVE" button.
            *******************************************************************/
            else if(strcmp("SAVE", which_button_label) == 0){

               char user_response[100];

               select_button(buttons, num_buttons, button_labels, "SAVE", colors[REDINDEX].pixel);

               /****************************************************************
               * Ask the user for a filename for saving the 24-bit image.
               ****************************************************************/
               printf("Do you want to save the displayed image in color? (y/n) ");
               scanf("%s", user_response);

               if((user_response[0] == 'Y') || (user_response[0] == 'y')){

		  unsigned char apixel;
                  int working_rows, working_cols;

                  printf("   Enter a filename (\".ppm\" is automatically appended): ");
                  scanf("%s", save_filename);

                  sprintf(new_save_filename, "%s.ppm", save_filename);

                  working_rows = ((fullresolution_image.rows > image.rows) ? image.rows : fullresolution_image.rows);
                  working_cols = ((fullresolution_image.cols > image.cols) ? image.cols : fullresolution_image.cols);

                  if((tmpximage = XGetImage(theDisplay, fullresolution_image.thePixmap,
                     0, 0, fullresolution_image.cols, fullresolution_image.rows,
                     (ULONG)AllPlanes, (int)ZPixmap)) == NULL){
                     fprintf(stderr, "Error with XGetImage().\n");
                  }

                  if((fpout = fopen(new_save_filename, "wb")) != NULL){

                     fprintf(fpout, "P6\n%d %d\n", working_cols, working_rows);
                     fprintf(fpout, "# Source Image Filename: %s\n", filename);
                     fprintf(fpout, "# Upper Left Coordinate (col=%d,row=%d)\n", image_ul_x, image_ul_y);
                     fprintf(fpout, "# Aggregation Factor = %d\n", current_aggfactor);
                     if(mapping == DENSITY_MAPPING) fprintf(fpout, "# Mapping = DENSITY_MAPPING\n");
                     else if(mapping == INTENSITY_MAPPING) fprintf(fpout, "# Mapping = INTENSITY_MAPPING\n");
                     fprintf(fpout, "255\n");

                     for(r=0;r<working_rows;r++){
                        for(c=0;c<working_cols;c++){

			   apixel = tmpximage->data[r*(tmpximage->bytes_per_line) + tmpximage->xoffset + c];

                           /* fputc((UCHAR)(colors[unmapper[apixel]].red / 256), fpout); */
                           /* fputc((UCHAR)(colors[unmapper[apixel]].green / 256), fpout); */
                           /* fputc((UCHAR)(colors[unmapper[apixel]].blue / 256), fpout); */

                           fputc((UCHAR)(colors[apixel].red / 256), fpout);
                           fputc((UCHAR)(colors[apixel].green / 256), fpout);
                           fputc((UCHAR)(colors[apixel].blue / 256), fpout);
                        }
                     }

                     fclose(fpout);

		     XDestroyImage(tmpximage);
                  }
                  else fprintf(stderr, "Error opening the file %s to save the data to.\n", new_save_filename);
               }

               /****************************************************************
               * Ask the user for a filename for saving the unsigned short
               * integer data. Then save the data to a file.
               ****************************************************************/
               printf("Do you want to save the displayed image 16-bit greyscale? (y/n) ");
               scanf("%s", user_response);

               if((user_response[0] == 'Y') || (user_response[0] == 'y')){

                  int working_rows, working_cols;

                  printf("Enter a filename (\".16bit.pgm\" is automatically appended): ");
                  scanf("%s", save_filename);

                  working_rows = ((fullresolution_image.rows > image.rows) ? image.rows : fullresolution_image.rows);
                  working_cols = ((fullresolution_image.cols > image.cols) ? image.cols : fullresolution_image.cols);

                  sprintf(new_save_filename, "%s.16bit.pgm", save_filename);
                  if((fpout = fopen(new_save_filename, "wb")) != NULL){

                     fprintf(fpout, "P5\n%d %d\n", working_cols, working_rows);
                     fprintf(fpout, "# Source Image Filename: %s\n", filename);
                     fprintf(fpout, "# Upper Left Coordinate (col=%d,row=%d)\n", image_ul_x, image_ul_y);
                     fprintf(fpout, "# Aggregation Factor = %d\n", current_aggfactor);
                     if(mapping == DENSITY_MAPPING) fprintf(fpout, "# Mapping = DENSITY_MAPPING\n");
                     else if(mapping == INTENSITY_MAPPING) fprintf(fpout, "# Mapping = INTENSITY_MAPPING\n");
                     fprintf(fpout, "# RESOLUTION = %lf MICRONS\n", (double)current_aggfactor * samplerate_in_microns);
                     fprintf(fpout, "# SCANNER = %s\n", scanner);
                     fprintf(fpout, "65535\n\n");

                     for(r=0;r<working_rows;r++){
                        fwrite(fullresolution_data + r*fullresolution_image.cols, sizeof(USHORT), working_cols, fpout);
                     }
                     fclose(fpout);

                  }
                  else fprintf(stderr, "Error opening the file %s to save the data to.\n", new_save_filename);
               }

               unselect_button(buttons, num_buttons, button_labels, "SAVE", cell[0]);
            }

            /*******************************************************************
            * This is the "SAVE SUBAREA" button.
            *******************************************************************/
            else if(strcmp("SAVE SUBAREA", which_button_label) == 0){

               char user_response[100];
               int x_coord[200];
               int y_coord[200];
	       int num_points;

               select_button(buttons, num_buttons, button_labels, "SAVE SUBAREA", colors[REDINDEX].pixel);

               /****************************************************************
               * Ask the user for a filename for saving the 24-bit image.
               ****************************************************************/
               printf("\nEnter the name of a file in which to store the sub area (\".ppm\" is appended): ");
               scanf("%s", save_filename);

               printf("\nSelect two points in the fullresolution image by pressing the mouse.\n");
               printf("Button at an upper left and then a lower right point.\n\n");

               XDefineCursor(theDisplay, fullresolution_image.theDrawWindow, crosshair_cursor);
               XFlush(theDisplay);

	       do{
                  XNextEvent(theDisplay, &theEvent);
                  /*************************************************************
                  * Handle expose events by repainting the screen.
                  *************************************************************/
                  if(theEvent.type == Expose){
                     if(theEvent.xbutton.window == overview_image.theDrawWindow)
                        repaint(overview_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     if(theEvent.xbutton.window == fullresolution_image.theDrawWindow)
                        repaint(fullresolution_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     if(theEvent.xbutton.window == histogram_image.theDrawWindow)
                        repaint(histogram_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     for(i=0;i<num_buttons;i++){
                        if(theEvent.xbutton.window == buttons[i].theDrawWindow)
                           repaint(buttons[i], theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     }
                  }
               }while(theEvent.type != ButtonPress);

               /*************************************************************
               * Keep grabbing new coordinates for the polygon each time the
               * left mouse button is pressed.
               *************************************************************/
               num_points = 0;
               while((theEvent.type == ButtonPress) && (num_points < 2) &&
                     (theEvent.xbutton.window == fullresolution_image.theDrawWindow)){

                  /*
                  x_coord[num_points] = current_aggfactor * (image_ul_x + theEvent.xbutton.x);
                  if(x_coord[num_points] >= cols) x_coord[num_points] = cols - 1;
                  y_coord[num_points] = current_aggfactor * (image_ul_y + theEvent.xbutton.y);
                  if(y_coord[num_points] >= rows) y_coord[num_points] = rows - 1;
                  */

                  x_coord[num_points] = theEvent.xbutton.x;
                  y_coord[num_points] = theEvent.xbutton.y;

                  num_points++;

		  if(num_points != 2){
                     do{
                        XNextEvent(theDisplay, &theEvent);
                        /**********************************************************
                        * Handle expose events by repainting the screen.
                        **********************************************************/
                        if(theEvent.type == Expose){
                           if(theEvent.xbutton.window == overview_image.theDrawWindow)
                              repaint(overview_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                      theEvent.xexpose.width, theEvent.xexpose.height,
                                      theEvent.xexpose.x, theEvent.xexpose.y);
                           if(theEvent.xbutton.window == fullresolution_image.theDrawWindow)
                              repaint(fullresolution_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                      theEvent.xexpose.width, theEvent.xexpose.height,
                                      theEvent.xexpose.x, theEvent.xexpose.y);
                           if(theEvent.xbutton.window == histogram_image.theDrawWindow)
                              repaint(histogram_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                      theEvent.xexpose.width, theEvent.xexpose.height,
                                      theEvent.xexpose.x, theEvent.xexpose.y);
                           for(i=0;i<num_buttons;i++){
                              if(theEvent.xbutton.window == buttons[i].theDrawWindow)
                                 repaint(buttons[i], theEvent.xexpose.x, theEvent.xexpose.y,
                                      theEvent.xexpose.width, theEvent.xexpose.height,
                                      theEvent.xexpose.x, theEvent.xexpose.y);
                           }
                        }
                     }while(theEvent.type != ButtonPress);
                  }
               }

               if(num_points == 2){

		  unsigned char apixel;
                  int working_rows, working_cols, start_row, start_col, sub_ul_row, sub_ul_col;

                  sprintf(new_save_filename, "%s.ppm", save_filename);

                  if(y_coord[1] > y_coord[0]){
                     start_row = y_coord[0];
                     working_rows = y_coord[1] - y_coord[0];
                  }
                  else{
                     start_row = y_coord[1];
                     working_rows = y_coord[0] - y_coord[1];
                  }

                  if(x_coord[1] > x_coord[0]){
                     start_col = x_coord[0];
                     working_cols = x_coord[1] - x_coord[0];
                  }
                  else{
                     start_col = x_coord[1];
                     working_cols = x_coord[0] - x_coord[1];
                  }

                  sub_ul_row = current_aggfactor * (image_ul_y + start_row);
                  sub_ul_col = current_aggfactor * (image_ul_x + start_col);

                  if((tmpximage = XGetImage(theDisplay, fullresolution_image.thePixmap,
                     start_col, start_row, working_cols, working_rows,
                     (ULONG)AllPlanes, (int)ZPixmap)) == NULL){
                     fprintf(stderr, "Error with XGetImage().\n");
                  }

                  if((fpout = fopen(new_save_filename, "wb")) != NULL){

                     fprintf(fpout, "P6\n%d %d\n", working_cols, working_rows);
                     fprintf(fpout, "# Source Image Filename: %s\n", filename);
                     fprintf(fpout, "# Upper Left Coordinate (col=%d,row=%d)\n", sub_ul_col, sub_ul_row);
                     fprintf(fpout, "# Aggregation Factor = %d\n", current_aggfactor);
                     if(mapping == DENSITY_MAPPING) fprintf(fpout, "# Mapping = DENSITY_MAPPING\n");
                     else if(mapping == INTENSITY_MAPPING) fprintf(fpout, "# Mapping = INTENSITY_MAPPING\n");
                     fprintf(fpout, "255\n");

                     for(r=0;r<working_rows;r++){
                        for(c=0;c<working_cols;c++){

			   apixel = tmpximage->data[r*(tmpximage->bytes_per_line) + tmpximage->xoffset + c];

                           /* fputc((UCHAR)(colors[unmapper[apixel]].red / 256), fpout); */
                           /* fputc((UCHAR)(colors[unmapper[apixel]].green / 256), fpout); */
                           /* fputc((UCHAR)(colors[unmapper[apixel]].blue / 256), fpout); */

                           fputc((UCHAR)(colors[apixel].red / 256), fpout);
                           fputc((UCHAR)(colors[apixel].green / 256), fpout);
                           fputc((UCHAR)(colors[apixel].blue / 256), fpout);
                        }
                     }

                     fclose(fpout);

		     XDestroyImage(tmpximage);
                  }
                  else fprintf(stderr, "Error opening the file %s to save the data to.\n", new_save_filename);
               }

               unselect_button(buttons, num_buttons, button_labels, "SAVE SUBAREA", cell[0]);
            }

            /*******************************************************************
            * This is the "PIXVAL" button.
            *******************************************************************/
            else if(strcmp("PIX VALUE", which_button_label) == 0){

               if(fullresolution_mode == REGIONAVGMODE)
                  unselect_button(buttons, num_buttons, button_labels, "REGION AVG", cell[0]);
               if(zoommenu_toggle == 1){
                  select_button(buttons, num_buttons, button_labels, "ZOOM", colors[REDINDEX].pixel);
                  for(i=0;i<num_zoombuttons;i++) XUnmapWindow(theDisplay, zoombuttons[i].theDrawWindow);
                  XUnmapWindow(theDisplay, zoom_image.theDrawWindow);
                  XUnmapWindow(theDisplay, zoommenu_image.theBorderWindow);
                  unselect_button(buttons, num_buttons, button_labels, "ZOOM", cell[0]);
                  zoommenu_toggle = 0;
               }
               select_button(buttons, num_buttons, button_labels, "PIX VALUE", colors[BLUEINDEX].pixel);
               fullresolution_mode = PIXVALMODE;
            }

            /*******************************************************************
            * This is the "SEGMENTATION" button.
            *******************************************************************/
            else if(strcmp("SEGMENTATION", which_button_label) == 0){

               if(segmentation_toggle == 1){
                  unselect_button(buttons, num_buttons, button_labels, "SEGMENTATION", cell[0]);
                  segmentation_toggle = 0;

               }
               else if(segmentation_toggle == 0){
                  select_button(buttons, num_buttons, button_labels, "SEGMENTATION", colors[BLUEINDEX].pixel);

                  segmentation_toggle = 1;
               }

               scale_min = -1;
               scale_max = -1;

               update_fullresolution_image(&image, &overview_image,
                  &fullresolution_image, &histogram_image, overview_x, overview_y,
                  scale_min, scale_max, fullresolution_data, 0, num_detections, detections,
                  isvalid, overlay_data);

               /*************************************************************
               * Repaint the crop box on the positioning (overview) image.
               *************************************************************/
               display_image(&overview_image, overview_image.rows,
                  overview_image.cols, 1);

               XSetForeground(theDisplay, overview_image.theGC,
                  colors[BLUEINDEX].pixel);  /* Blue */

               XDrawRectangle(theDisplay, overview_image.thePixmap,
                  overview_image.theGC, box_x, box_y, box_width, box_height);

               if((overlay_filename[0] != '\0') && (show_overlay == 1)){
                  overlay_groundtruth(&overview_image, 0, 0,
                     overview_subfactor, isvalid, overlay_data);
               }

               if(segmentation_toggle)
                  overlay_features(&overview_image, 0, 0, overview_subfactor);

               if(display_regions_toggle)
                  overlay_regions(&overview_image, 0, 0, overview_subfactor,
                     num_detections, detections);

               repaint(overview_image, 0, 0, overview_image.cols,
                  overview_image.rows, 0, 0);
            }

            /*******************************************************************
            * This is the "REGION AVG" button.
            *******************************************************************/
            else if(strcmp("REGION AVG", which_button_label) == 0){

               if(fullresolution_mode == PIXVALMODE)
                  unselect_button(buttons, num_buttons, button_labels, "PIX VALUE", cell[0]);

               if(zoommenu_toggle == 1){
                  select_button(buttons, num_buttons, button_labels, "ZOOM", colors[REDINDEX].pixel);
                  for(i=0;i<num_zoombuttons;i++) XUnmapWindow(theDisplay, zoombuttons[i].theDrawWindow);
                  XUnmapWindow(theDisplay, zoom_image.theDrawWindow);
                  XUnmapWindow(theDisplay, zoommenu_image.theBorderWindow);
                  unselect_button(buttons, num_buttons, button_labels, "ZOOM", cell[0]);
                  zoommenu_toggle = 0;
               }
               select_button(buttons, num_buttons, button_labels, "REGION AVG", colors[BLUEINDEX].pixel);
               fullresolution_mode = REGIONAVGMODE;
            }
            /*******************************************************************
            * This is the "DISTANCE" button.
            *******************************************************************/
            else if(strcmp("DISTANCE", which_button_label) == 0){

               int x_coord[2];
               int y_coord[2];
	       int num_points;
               int count, p;
               double pixeldist;

               select_button(buttons, num_buttons, button_labels, "DISTANCE", colors[REDINDEX].pixel);

               XDefineCursor(theDisplay, fullresolution_image.theDrawWindow, crosshair_cursor);
               XFlush(theDisplay);

	       /****************************************************************
               * Copy the portion of the image that we will be drawing a block
               * of text on so we can put it back on the screen later.
	       ****************************************************************/
               if((tmpximage = XGetImage(theDisplay, fullresolution_image.thePixmap,
                  0, 0, 180, 7*info_helvb14->ascent, (ULONG)AllPlanes, (int)ZPixmap)) == NULL){
                  fprintf(stderr, "Error with XGetImage().\n");
               }

               XSetForeground(theDisplay, fullresolution_image.theGC, cell[32]);
               XFillRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                  0, 0, 180, 7*info_helvb14->ascent);
               XSetForeground(theDisplay, fullresolution_image.theGC, cell[0]);
               XDrawRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC, 0, 0,
                  180-1, (7*info_helvb14->ascent)-1);

               sprintf(string1, "Click on 2 ");
               twidth = XTextWidth(info_helvb14, string1, strlen(string1));
               sprintf(string2, "points to   ", count);
               sprintf(string3, "measure a   ");
               sprintf(string4, "distance.   ");
               XSetForeground(theDisplay, fullresolution_image.theGC, cell[253]);
   
               XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                  180/2 - twidth/2, 1.5*info_helvb14->ascent, string1, strlen(string1));
               XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                  180/2 - twidth/2, 3*info_helvb14->ascent, string2, strlen(string2));
               XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                  180/2 - twidth/2, 4.5*info_helvb14->ascent, string3, strlen(string3));
               XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                  180/2 - twidth/2, 6*info_helvb14->ascent, string4, strlen(string4));

               repaint(fullresolution_image, 0, 0, 180, 7*info_helvb14->ascent, 0, 0);

               num_points = 0;

               XSetForeground(theDisplay, fullresolution_image.theGC, colors[ORANGEINDEX].pixel);

	       do{
                  XNextEvent(theDisplay, &theEvent);
                  /*************************************************************
                  * Handle expose events by repainting the screen.
                  *************************************************************/
                  if(theEvent.type == Expose){
                     if(theEvent.xbutton.window == overview_image.theDrawWindow)
                        repaint(overview_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     if(theEvent.xbutton.window == fullresolution_image.theDrawWindow)
                        repaint(fullresolution_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     if(theEvent.xbutton.window == histogram_image.theDrawWindow)
                        repaint(histogram_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     for(i=0;i<num_buttons;i++){
                        if(theEvent.xbutton.window == buttons[i].theDrawWindow)
                           repaint(buttons[i], theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     }
                  }
               }while(!((theEvent.xbutton.window == fullresolution_image.theDrawWindow) && (theEvent.type == ButtonPress)));

               /*************************************************************
               * Keep grabbing new coordinates for the polygon each time the
               * left mouse button is pressed.
               *************************************************************/
               while((theEvent.type == ButtonPress) &&
                     (theEvent.xbutton.window == fullresolution_image.theDrawWindow) &&
                     (num_points < 2)){

                  x_coord[num_points] = current_aggfactor * (image_ul_x + theEvent.xbutton.x);
                  if(x_coord[num_points] >= cols) x_coord[num_points] = cols - 1;
                  y_coord[num_points] = current_aggfactor * (image_ul_y + theEvent.xbutton.y);
                  if(y_coord[num_points] >= rows) y_coord[num_points] = rows - 1;

                  XDrawArc(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     theEvent.xbutton.x - 4, theEvent.xbutton.y - 4,
                     8, 8, 0, 360*64);

                  if(num_points == 1){
                     XDrawLine(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                        (int)(x_coord[num_points-1]/current_aggfactor) - image_ul_x,
                        (int)(y_coord[num_points-1]/current_aggfactor) - image_ul_y,
                        (int)(x_coord[num_points]/current_aggfactor) - image_ul_x,
                        (int)(y_coord[num_points]/current_aggfactor) - image_ul_y);
                  }

                  /*
                  for(p=0;p<=num_points;p++){
                     printf("(%d %d) ", x_coord[p], y_coord[p]);
                  }
                  printf("\n");
                  */

                  num_points++;

                  repaint(fullresolution_image, 0, 0, fullresolution_image.cols,
                     fullresolution_image.rows, 0, 0);

                  if(num_points < 2){
                     do{
                        XNextEvent(theDisplay, &theEvent);
                        /**********************************************************
                        * Handle expose events by repainting the screen.
                        **********************************************************/
                        if(theEvent.type == Expose){
                           if(theEvent.xbutton.window == overview_image.theDrawWindow)
                              repaint(overview_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                      theEvent.xexpose.width, theEvent.xexpose.height,
                                      theEvent.xexpose.x, theEvent.xexpose.y);
                           if(theEvent.xbutton.window == fullresolution_image.theDrawWindow)
                              repaint(fullresolution_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                      theEvent.xexpose.width, theEvent.xexpose.height,
                                      theEvent.xexpose.x, theEvent.xexpose.y);
                           if(theEvent.xbutton.window == histogram_image.theDrawWindow)
                              repaint(histogram_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                      theEvent.xexpose.width, theEvent.xexpose.height,
                                      theEvent.xexpose.x, theEvent.xexpose.y);
                           for(i=0;i<num_buttons;i++){
                              if(theEvent.xbutton.window == buttons[i].theDrawWindow)
                                 repaint(buttons[i], theEvent.xexpose.x, theEvent.xexpose.y,
                                      theEvent.xexpose.width, theEvent.xexpose.height,
                                      theEvent.xexpose.x, theEvent.xexpose.y);
                           }
                        }
                     }while(theEvent.type != ButtonPress);
                  }
               }

               if(num_points == 2){

                  /*************************************************************
                  * Let the user know that the computer is working calculating
                  * the statistics.
                  *************************************************************/
                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[32]);
                  XFillRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     0, 0, 180, 7*info_helvb14->ascent);
                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[0]);
                  XDrawRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC, 0, 0,
                     180-1, (7*info_helvb14->ascent)-1);
   
                  sprintf(string1, "P1:%4d,%4d", x_coord[0], y_coord[0]);
                  twidth = XTextWidth(info_helvb14, string1, strlen(string1));
                  sprintf(string2, "P2:%4d,%4d", x_coord[1], y_coord[1]);
                  pixeldist = sqrt((double)(x_coord[0]-x_coord[1])*(double)(x_coord[0]-x_coord[1]) +
						 (double)(y_coord[0]-y_coord[1])*(double)(y_coord[0]-y_coord[1]));
                  sprintf(string3, "PDIST %6d", (int)pixeldist);
                  sprintf(string4, "%10.4fmm", (pixeldist * samplerate_in_microns / 1000.0));
                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[253]);
      
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 1.5*info_helvb14->ascent, string1, strlen(string1));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 3*info_helvb14->ascent, string2, strlen(string2));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 4.5*info_helvb14->ascent, string3, strlen(string3));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 6*info_helvb14->ascent, string4, strlen(string4));
   
                  repaint(fullresolution_image, 0, 0, 180, 7*info_helvb14->ascent, 0, 0);

                  /*
                  polygon_stats(&image_full, x_coord, y_coord, num_points, &count, &mean, &stdev);

                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[32]);
                  XFillRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     0, 0, 180, 7*info_helvb14->ascent);
                  if(mapping == DENSITY_MAPPING){
                     sprintf(string1, "REGION STATS (D)");
                     twidth = XTextWidth(info_helvb14, string1, strlen(string1));
                     sprintf(string2, "COUNT = %d", count);
                     sprintf(string3, "MEAN  = %0.3lf", (od_offset-(double)mean)/od_scale+min_density);
                     sprintf(string4, "STDEV = %0.3lf", stdev/od_scale);
                  }
                  else{
                     sprintf(string1, "REGION STATS (I)");
                     twidth = XTextWidth(info_helvb14, string1, strlen(string1));
                     sprintf(string2, "COUNT = %d", count);
                     sprintf(string3, "MEAN  = %0.3lf", mean);
                     sprintf(string4, "STDEV = %0.3lf", stdev/od_scale);
                  }
                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[253]);
   
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 1.5*info_helvb14->ascent, string1, strlen(string1));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 3*info_helvb14->ascent, string2, strlen(string2));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 4.5*info_helvb14->ascent, string3, strlen(string3));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 6*info_helvb14->ascent, string4, strlen(string4));

                  repaint(fullresolution_image, 0, 0, 180, 7*info_helvb14->ascent, 0, 0);
                  */
 
                  select_button(buttons, num_buttons, button_labels, "DISTANCE", colors[BLUEINDEX].pixel);

                  do{
                     XNextEvent(theDisplay, &theEvent);
                     /*************************************************************************
                     * Handle expose events by repainting the screen.
                     *************************************************************************/
                     if(theEvent.type == Expose){
                        if(theEvent.xbutton.window == overview_image.theDrawWindow)
                           repaint(overview_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height, theEvent.xexpose.x,
                                   theEvent.xexpose.y);
                        if(theEvent.xbutton.window == fullresolution_image.theDrawWindow)
                           repaint(fullresolution_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height, theEvent.xexpose.x,
                                   theEvent.xexpose.y);
                        if(theEvent.xbutton.window == histogram_image.theDrawWindow)
                           repaint(histogram_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height, theEvent.xexpose.x,
                                   theEvent.xexpose.y);
                        for(i=0;i<num_buttons;i++){
                           if(theEvent.xbutton.window == buttons[i].theDrawWindow)
                              repaint(buttons[i], theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height, theEvent.xexpose.x,
                                   theEvent.xexpose.y);
                        }
                     }
                  }while(theEvent.type != ButtonPress);

                  XPutImage(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC, tmpximage,
		     0, 0, 0, 0, 180, 7*info_helvb14->ascent);
                  XDestroyImage(tmpximage);
                  repaint(fullresolution_image, 0, 0, 180, 7*info_helvb14->ascent, 0, 0);
               }

               select_button(buttons, num_buttons, button_labels, "DISTANCE", colors[0].pixel);
            }

            /*******************************************************************
            * This is the "OVERVIEW" button.
            *******************************************************************/
            else if(strcmp("OVERVIEW", which_button_label) == 0){

	       if(overview_toggle == 0){
                  select_button(buttons, num_buttons, button_labels, "OVERVIEW", colors[REDINDEX].pixel);
                  XMapWindow(theDisplay, overview_image.theBorderWindow);
                  XMapWindow(theDisplay, overview_image.theDrawWindow);
                  select_button(buttons, num_buttons, button_labels, "OVERVIEW", colors[BLUEINDEX].pixel);
                  overview_toggle = 1;
               }
               else if(overview_toggle == 1){
                  select_button(buttons, num_buttons, button_labels, "OVERVIEW", colors[REDINDEX].pixel);
                  XUnmapWindow(theDisplay, overview_image.theBorderWindow);
                  XUnmapWindow(theDisplay, overview_image.theDrawWindow);
                  unselect_button(buttons, num_buttons, button_labels, "OVERVIEW", cell[0]);
                  overview_toggle = 0;
               }
            }
            /*******************************************************************
            * This is the "SHARPEN" button.
            *******************************************************************/
            else if(strcmp("SHARPEN", which_button_label) == 0){

               select_button(buttons, num_buttons, button_labels, "SHARPEN", colors[REDINDEX].pixel);

               
               gaussian_smooth(fullresolution_data, fullresolution_image.rows,
                  fullresolution_image.cols, 11.0, &temp_float); /* 4.5 */

               for(r=0,pos=0;r<(fullresolution_image.rows);r++){
                  for(c=0;c<(fullresolution_image.cols); c++, pos++){
                     temp_long = (long)fullresolution_data[pos] +
			1.05 * ((double)fullresolution_data[pos] - (double)temp_float[pos]);
                     if(temp_long < 0) fullresolution_data[pos] = 0;
                     else if(temp_long > 65535) fullresolution_data[pos] = 65535;
                     else fullresolution_data[pos] = (unsigned short int)temp_long;
                  }
               }

               free(temp_float);

               update_fullresolution_image(&image, &overview_image,
                  &fullresolution_image, &histogram_image, overview_x, overview_y,
                  scale_min, scale_max, fullresolution_data, 0, num_detections, detections,
                  isvalid, overlay_data);

               unselect_button(buttons, num_buttons, button_labels, "SHARPEN", cell[0]);

               select_button(buttons, num_buttons, button_labels, "SHARPEN", colors[BLUEINDEX].pixel);
            }

            /*******************************************************************
            * This is the "HISTOGRAM" button.
            *******************************************************************/
            else if(strcmp("HISTOGRAM", which_button_label) == 0){

	       if(histogram_toggle == 0){
                  select_button(buttons, num_buttons, button_labels, "HISTOGRAM", colors[REDINDEX].pixel);
                  XMapWindow(theDisplay, histogram_image.theBorderWindow);
                  XMapWindow(theDisplay, histogram_image.theDrawWindow);
                  select_button(buttons, num_buttons, button_labels, "HISTOGRAM", colors[BLUEINDEX].pixel);
                  histogram_toggle = 1;
               }
               else if(histogram_toggle == 1){
                  select_button(buttons, num_buttons, button_labels, "HISTOGRAM", colors[REDINDEX].pixel);
                  XUnmapWindow(theDisplay, histogram_image.theBorderWindow);
                  XUnmapWindow(theDisplay, histogram_image.theDrawWindow);
                  unselect_button(buttons, num_buttons, button_labels, "HISTOGRAM", cell[0]);
                  histogram_toggle = 0;
               }
            }

            /*******************************************************************
            * This is the "ZOOM" button.
            *******************************************************************/
            else if(strcmp("ZOOM", which_button_label) == 0){

               if(fullresolution_mode == PIXVALMODE)
                  unselect_button(buttons, num_buttons, button_labels, "PIX VALUE", cell[0]);
               if(fullresolution_mode == REGIONAVGMODE)
                  unselect_button(buttons, num_buttons, button_labels, "REGION AVG", cell[0]);

               if(zoommenu_toggle == 0){
                  select_button(buttons, num_buttons, button_labels, "ZOOM", colors[REDINDEX].pixel);
                  XMapWindow(theDisplay, zoommenu_image.theBorderWindow);
                  XMapWindow(theDisplay, zoom_image.theDrawWindow);
                  for(i=0;i<num_zoombuttons;i++) XMapWindow(theDisplay, zoombuttons[i].theDrawWindow);
                  select_button(buttons, num_buttons, button_labels, "ZOOM", colors[BLUEINDEX].pixel);
                  zoommenu_toggle = 1;
                  fullresolution_mode = ZOOMMODE;
                  select_button(zoombuttons, num_zoombuttons, zoombutton_labels,
                     zoombutton_labels[5-zoom_factor], colors[BLUEINDEX].pixel);
               }
               else if(zoommenu_toggle == 1){
                  select_button(buttons, num_buttons, button_labels, "ZOOM", colors[REDINDEX].pixel);
                  for(i=0;i<num_zoombuttons;i++) XUnmapWindow(theDisplay, zoombuttons[i].theDrawWindow);
                  XUnmapWindow(theDisplay, zoom_image.theDrawWindow);
                  XUnmapWindow(theDisplay, zoommenu_image.theBorderWindow);
                  unselect_button(buttons, num_buttons, button_labels, "ZOOM", cell[0]);
                  zoommenu_toggle = 0;
                  fullresolution_mode = PIXVALMODE;
               }
            }

            /*******************************************************************
            * This handles the zoom buttons.
            *******************************************************************/
            else{
               for(j=0;j<num_zoombuttons;j++){
                  if(strcmp(zoombutton_labels[j], which_button_label) == 0){
                     for(i=0;i<num_zoombuttons;i++)
                        unselect_button(zoombuttons, num_zoombuttons, zoombutton_labels, zoombutton_labels[i], cell[0]);
                     zoom_factor = 5 - j;
                     select_button(zoombuttons, num_zoombuttons, zoombutton_labels, zoombutton_labels[j],
                        colors[BLUEINDEX].pixel);

                     update_zoom_window(&image_full, &zoom_image, last_zoompos_x, last_zoompos_y, zoom_factor);

                     break;
                  }
               }
            }

            XDefineCursor(theDisplay, overview_image.theDrawWindow, crosshair_cursor);
            XDefineCursor(theDisplay, fullresolution_image.theDrawWindow, crosshair_cursor);
            XDefineCursor(theDisplay, histogram_image.theDrawWindow, arrow_cursor);
            XFlush(theDisplay);

         }

         /**********************************************************************
         * Mouse clicks in the fullresolution window are handled here. The
         * interpretation of the mouse click is context dependent. By this I
         * mean the the processing that needs to be done depends on the current
         * "mode".
         **********************************************************************/
         if(theEvent.xbutton.window == fullresolution_image.theDrawWindow){

            if(theEvent.xbutton.button == RIGHT_BUTTON){
               if(imagemenu_toggle == 0){
                  XMapWindow(theDisplay, imagemenu_image.theBorderWindow);
                  for(i=0;i<num_buttons;i++) XMapWindow(theDisplay, buttons[i].theDrawWindow);
                  imagemenu_toggle = 1;
               }
               else if(imagemenu_toggle == 1){
                  XUnmapWindow(theDisplay, imagemenu_image.theBorderWindow);
                  for(i=0;i<num_buttons;i++) XUnmapWindow(theDisplay, buttons[i].theDrawWindow);
                  imagemenu_toggle = 0;
               }
            }

            /*******************************************************************
            * When the user presses any mouse button in the fullresolution
            * window while in zoom mode, they get a zoomed up piece of the
            * image displayed.
            *******************************************************************/
            else if(fullresolution_mode == ZOOMMODE){

               update_zoom_window(&image_full, &zoom_image, theEvent.xbutton.x, theEvent.xbutton.y, zoom_factor);
               last_zoompos_x = theEvent.xbutton.x;
               last_zoompos_y = theEvent.xbutton.y;

/*
	       int est_file_x, est_file_y, est_rows, est_cols, er, ec, err, ecc, e_val;
               unsigned char e_val_uchar;
	       est_file_x = current_aggfactor * (image_ul_x + theEvent.xbutton.x);
               est_file_y = current_aggfactor * (image_ul_y + theEvent.xbutton.y);

	       est_rows = ceil(512.0 / (float)zoom_factor);
               est_cols = ceil(512.0 / (float)zoom_factor);

	       est_file_x -= (est_cols/2);
               if((est_file_x + est_cols) >= image_full.cols) est_file_x = image_full.cols-1-est_cols;
               if(est_file_x < 0) est_file_x = 0;

	       est_file_y -= (est_rows/2);
               if((est_file_y + est_rows) >= image_full.rows) est_file_y = image_full.rows-1-est_rows;
               if(est_file_y < 0) est_file_y = 0;

               for(er=0;er<est_rows;er++){
                  for(ec=0;ec<est_cols;ec++){
                     e_val =  (int)(*(image_full.getpixel))(&image_full, est_file_y+er, est_file_x+ec);

                     if((int)e_val < last_min) e_val_uchar= 0;
                     else if((int)e_val > last_max) e_val_uchar = 255;
                     else e_val_uchar = (unsigned char)((((float)e_val-(float)last_min)*255.0)/((float)last_max-(float)last_min));

                     for(err=0;err<zoom_factor;err++){
                        for(ecc=0;ecc<zoom_factor;ecc++){
                           zoom_image.image[(er*zoom_factor+err) * 512 + (ec*zoom_factor+ecc)] = e_val_uchar;
                        }
                     }
                  }
               }
	       XRaiseWindow(theDisplay, zoom_image.theBorderWindow);
               display_image(&(zoom_image), 512, 512, 1);
*/

            }

            /*******************************************************************
            * When the user presses any mouse button in the fullresolution
            * window while in pixel val mode, they get pixel values displayed
            * while they hold the button down.
            *******************************************************************/
            else if(fullresolution_mode == PIXVALMODE){

               if((tmpximage = XGetImage(theDisplay, fullresolution_image.thePixmap,
                  0, 0, 180, 7*info_helvb14->ascent, (ULONG)AllPlanes, (int)ZPixmap)) == NULL){
                  fprintf(stderr, "Error with XGetImage().\n");
               }

               do{

                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[32]);
                  XFillRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
		     0, 0, 180, 7*info_helvb14->ascent);
                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[0]);
                  XDrawRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC, 0, 0,
		     180-1, (7*info_helvb14->ascent)-1);

                  cx = image_ul_x + theEvent.xbutton.x;
                  if(cx == 0) cx = 1;
                  else if(cx == image.cols-1) cx = image.cols-2;
                  cy = image_ul_y + theEvent.xbutton.y;
                  if(cy == 0) cy = 1;
                  else if(cy == image.rows-1) cy = image.rows-2;

                  if(mapping == DENSITY_MAPPING){
                     sprintf(string1, "DENSITY (%4d,%4d)", current_aggfactor * cx, current_aggfactor * cy);
                     twidth = XTextWidth(info_helvb14, string1, strlen(string1));
                     sprintf(string2, "%0.3lf  %0.3lf  %0.3lf",
                        (od_offset-(double)(*(image.getpixel))(&image, cy-1, cx-1))/od_scale+min_density,
                        (od_offset-(double)(*(image.getpixel))(&image, cy-1, cx))/od_scale+min_density,
                        (od_offset-(double)(*(image.getpixel))(&image, cy-1, cx+1))/od_scale+min_density);
                     sprintf(string3, "%0.3lf  %0.3lf  %0.3lf",
                        (od_offset-(double)(*(image.getpixel))(&image, cy, cx-1))/od_scale+min_density,
                        (od_offset-(double)(*(image.getpixel))(&image, cy, cx))/od_scale+min_density,
                        (od_offset-(double)(*(image.getpixel))(&image, cy, cx+1))/od_scale+min_density);
                     sprintf(string4, "%0.3lf  %0.3lf  %0.3lf",
                        (od_offset-(double)(*(image.getpixel))(&image, cy+1, cx-1))/od_scale+min_density,
                        (od_offset-(double)(*(image.getpixel))(&image, cy+1, cx))/od_scale+min_density,
                        (od_offset-(double)(*(image.getpixel))(&image, cy+1, cx+1))/od_scale+min_density);
                  }
                  else{
                     sprintf(string1, "INTENSITY %4d,%4d", current_aggfactor * cx, current_aggfactor * cy);
                     twidth = XTextWidth(info_helvb14, string1, strlen(string1));
                     sprintf(string2, "  %5d %5d %5d ",
                        (int)((*(image.getpixel))(&image, cy-1, cx-1)),
                        (int)((*(image.getpixel))(&image, cy-1, cx)),
                        (int)((*(image.getpixel))(&image, cy-1, cx+1)));
                     sprintf(string3, "  %5d %5d %5d ",
                        (int)((*(image.getpixel))(&image, cy, cx-1)),
                        (int)((*(image.getpixel))(&image, cy, cx)),
                        (int)((*(image.getpixel))(&image, cy, cx+1)));
                     sprintf(string4, "  %5d %5d %5d ",
                        (int)((*(image.getpixel))(&image, cy+1, cx-1)),
                        (int)((*(image.getpixel))(&image, cy+1, cx)),
                        (int)((*(image.getpixel))(&image, cy+1, cx+1)));
                  }
                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[253]);

                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 1.5*info_helvb14->ascent, string1, strlen(string1));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 3*info_helvb14->ascent, string2, strlen(string2));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 4.5*info_helvb14->ascent, string3, strlen(string3));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 6*info_helvb14->ascent, string4, strlen(string4));

                  repaint(fullresolution_image, 0, 0, 180, 7*info_helvb14->ascent, 0, 0);

                  /*
                  XPeekEvent(theDisplay, &theFutureEvent);
		  if(theFutureEvent.type != ButtonRelease) XNextEvent(theDisplay, &theEvent);
                  */
               }while(0);      /* theFutureEvent.type != ButtonRelease); */

	       do{
                  XNextEvent(theDisplay, &theEvent);
                  /*************************************************************
                  * Handle expose events by repainting the screen.
                  *************************************************************/
                  if(theEvent.type == Expose){
                     if(theEvent.xbutton.window == overview_image.theDrawWindow)
                        repaint(overview_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     if(theEvent.xbutton.window == fullresolution_image.theDrawWindow)
                        repaint(fullresolution_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     if(theEvent.xbutton.window == histogram_image.theDrawWindow)
                        repaint(histogram_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     for(i=0;i<num_buttons;i++){
                        if(theEvent.xbutton.window == buttons[i].theDrawWindow)
                           repaint(buttons[i], theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     }
                  }
               }while(theEvent.type != ButtonRelease);

               XPutImage(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC, tmpximage,
                  0, 0, 0, 0, 180, 7*info_helvb14->ascent);
               XDestroyImage(tmpximage);

               repaint(fullresolution_image, 0, 0, 180, 7*info_helvb14->ascent, 0, 0);

               continue;
            }

            /*******************************************************************
            * When the user presses the left mouse button in the fullresolution
            * window while in region avg mode, they draw a polygon with
            * successive left mouse clicks. A right mouse button click defines
            * a new vertex and then completes the polygon.
            *******************************************************************/
            else if((fullresolution_mode == REGIONAVGMODE) && (theEvent.xbutton.button == LEFT_BUTTON)){

               FILE *fp_points=NULL;
               char points_filename[200];
               int x_coord[200];
               int y_coord[200];
	       int num_points;
               int count, p;
               double mean=0.0, stdev=0.0;

               select_button(buttons, num_buttons, button_labels, "REGION AVG", colors[REDINDEX].pixel);

               /****************************************************************
               * If the user already pressed the button to save each set of
               * points selected, open the file in append mode.
               ****************************************************************/
	       if(save_points_toggle == 1){
                  sprintf(points_filename, "%s.pts", filename);
                  if((fp_points = fopen(points_filename, "a")) == NULL){
                     fprintf(stderr, "Error opening the points file names %s.\n", points_filename);
                  }
               }

	       /****************************************************************
               * Copy the portion of the image that we will be drawing a block
               * of text on so we can put it back on the screen later.
	       ****************************************************************/
               if((tmpximage = XGetImage(theDisplay, fullresolution_image.thePixmap,
                  0, 0, 180, 7*info_helvb14->ascent, (ULONG)AllPlanes, (int)ZPixmap)) == NULL){
                  fprintf(stderr, "Error with XGetImage().\n");
               }

               XSetForeground(theDisplay, fullresolution_image.theGC, cell[32]);
               XFillRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                  0, 0, 180, 7*info_helvb14->ascent);
               XSetForeground(theDisplay, fullresolution_image.theGC, cell[0]);
               XDrawRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC, 0, 0,
                  180-1, (7*info_helvb14->ascent)-1);

               sprintf(string1, "Left click ");
               twidth = XTextWidth(info_helvb14, string1, strlen(string1));
               sprintf(string2, " adds points", count);
               sprintf(string3, "Right click ");
               sprintf(string4, " closes poly");
               XSetForeground(theDisplay, fullresolution_image.theGC, cell[253]);
   
               XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                  180/2 - twidth/2, 1.5*info_helvb14->ascent, string1, strlen(string1));
               XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                  180/2 - twidth/2, 3*info_helvb14->ascent, string2, strlen(string2));
               XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                  180/2 - twidth/2, 4.5*info_helvb14->ascent, string3, strlen(string3));
               XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                  180/2 - twidth/2, 6*info_helvb14->ascent, string4, strlen(string4));

               repaint(fullresolution_image, 0, 0, 180, 7*info_helvb14->ascent, 0, 0);

               /*************************************************************
               * Calculate the approximate coordiate in the fullresolution image.
               *************************************************************/
               num_points = 0;
               x_coord[num_points] = current_aggfactor * (image_ul_x + theEvent.xbutton.x);
               if(x_coord[num_points] >= cols) x_coord[num_points] = cols - 1;
               y_coord[num_points] = current_aggfactor * (image_ul_y + theEvent.xbutton.y);
               if(y_coord[num_points] >= rows) y_coord[num_points] = rows - 1;
               num_points++;

               XSetForeground(theDisplay, fullresolution_image.theGC, colors[ORANGEINDEX].pixel);

	       do{
                  XNextEvent(theDisplay, &theEvent);
                  /*************************************************************
                  * Handle expose events by repainting the screen.
                  *************************************************************/
                  if(theEvent.type == Expose){
                     if(theEvent.xbutton.window == overview_image.theDrawWindow)
                        repaint(overview_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     if(theEvent.xbutton.window == fullresolution_image.theDrawWindow)
                        repaint(fullresolution_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     if(theEvent.xbutton.window == histogram_image.theDrawWindow)
                        repaint(histogram_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     for(i=0;i<num_buttons;i++){
                        if(theEvent.xbutton.window == buttons[i].theDrawWindow)
                           repaint(buttons[i], theEvent.xexpose.x, theEvent.xexpose.y,
                                theEvent.xexpose.width, theEvent.xexpose.height,
                                theEvent.xexpose.x, theEvent.xexpose.y);
                     }
                  }
               }while(theEvent.type != ButtonPress);

               /*************************************************************
               * Keep grabbing new coordinates for the polygon each time the
               * left mouse button is pressed.
               *************************************************************/
               while((theEvent.type == ButtonPress) &&
                     (theEvent.xbutton.window == fullresolution_image.theDrawWindow)){

                  x_coord[num_points] = current_aggfactor * (image_ul_x + theEvent.xbutton.x);
                  if(x_coord[num_points] >= cols) x_coord[num_points] = cols - 1;
                  y_coord[num_points] = current_aggfactor * (image_ul_y + theEvent.xbutton.y);
                  if(y_coord[num_points] >= rows) y_coord[num_points] = rows - 1;

                  XDrawLine(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     (int)(x_coord[num_points-1]/current_aggfactor) - image_ul_x,
                     (int)(y_coord[num_points-1]/current_aggfactor) - image_ul_y,
                     (int)(x_coord[num_points]/current_aggfactor) - image_ul_x,
                     (int)(y_coord[num_points]/current_aggfactor) - image_ul_y);

                  repaint(fullresolution_image, 0, 0, fullresolution_image.cols,
                     fullresolution_image.rows, 0, 0);

                  if(theEvent.xbutton.button == RIGHT_BUTTON){

                     XDrawLine(theDisplay, fullresolution_image.thePixmap,
                        fullresolution_image.theGC,
                        (int)(x_coord[num_points]/current_aggfactor) - image_ul_x,
                        (int)(y_coord[num_points]/current_aggfactor) - image_ul_y,
                        (int)(x_coord[0]/current_aggfactor) - image_ul_x,
                        (int)(y_coord[0]/current_aggfactor) - image_ul_y);

                     repaint(fullresolution_image, 0, 0,
                        fullresolution_image.cols, fullresolution_image.rows, 0, 0);

                     num_points++;

                     break;
                  }

                  num_points++;

                  do{
                     XNextEvent(theDisplay, &theEvent);
                     /**********************************************************
                     * Handle expose events by repainting the screen.
                     **********************************************************/
                     if(theEvent.type == Expose){
                        if(theEvent.xbutton.window == overview_image.theDrawWindow)
                           repaint(overview_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height,
                                   theEvent.xexpose.x, theEvent.xexpose.y);
                        if(theEvent.xbutton.window == fullresolution_image.theDrawWindow)
                           repaint(fullresolution_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height,
                                   theEvent.xexpose.x, theEvent.xexpose.y);
                        if(theEvent.xbutton.window == histogram_image.theDrawWindow)
                           repaint(histogram_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height,
                                   theEvent.xexpose.x, theEvent.xexpose.y);
                        for(i=0;i<num_buttons;i++){
                           if(theEvent.xbutton.window == buttons[i].theDrawWindow)
                              repaint(buttons[i], theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height,
                                   theEvent.xexpose.x, theEvent.xexpose.y);
                        }
                     }
                  }while(theEvent.type != ButtonPress);

               }

               if(num_points >= 3){

                  /****************************************************************
                  * If the user already pressed the button to save each set of
                  * points selected, write the points to the file.
                  ****************************************************************/
	          if(save_points_toggle == 1){
                     fprintf(fp_points, "%s %d %d\n", filename, rows, cols);
                     fprintf(fp_points, "%d\n", num_points);
                     for(p=0;p<num_points;p++) fprintf(fp_points, "%d %d\n", y_coord[p], x_coord[p]);
                  }

                  /*************************************************************
                  * Let the user know that the computer is working calculating
                  * the statistics.
                  *************************************************************/
                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[32]);
                  XFillRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     0, 0, 180, 7*info_helvb14->ascent);
                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[0]);
                  XDrawRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC, 0, 0,
                     180-1, (7*info_helvb14->ascent)-1);
   
                  sprintf(string1, "Please Wait ");
                  twidth = XTextWidth(info_helvb14, string1, strlen(string1));
                  sprintf(string2, "------------", count);
                  sprintf(string3, "Calculating ");
                  sprintf(string4, "Statistics  ");
                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[253]);
      
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 1.5*info_helvb14->ascent, string1, strlen(string1));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 3*info_helvb14->ascent, string2, strlen(string2));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 4.5*info_helvb14->ascent, string3, strlen(string3));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 6*info_helvb14->ascent, string4, strlen(string4));
   
                  repaint(fullresolution_image, 0, 0, 180, 7*info_helvb14->ascent, 0, 0);

                  polygon_stats(&image_full, x_coord, y_coord, num_points, &count, &mean, &stdev);

                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[32]);
                  XFillRectangle(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     0, 0, 180, 7*info_helvb14->ascent);
                  if(mapping == DENSITY_MAPPING){
                     sprintf(string1, "REGION STATS (D)");
                     twidth = XTextWidth(info_helvb14, string1, strlen(string1));
                     sprintf(string2, "COUNT = %d", count);
                     sprintf(string3, "MEAN  = %0.3lf", (od_offset-(double)mean)/od_scale+min_density);
                     sprintf(string4, "STDEV = %0.3lf", stdev/od_scale);
                  }
                  else{
                     sprintf(string1, "REGION STATS (I)");
                     twidth = XTextWidth(info_helvb14, string1, strlen(string1));
                     sprintf(string2, "COUNT = %d", count);
                     sprintf(string3, "MEAN  = %0.3lf", mean);
                     sprintf(string4, "STDEV = %0.3lf", stdev/od_scale);
                  }
                  XSetForeground(theDisplay, fullresolution_image.theGC, cell[253]);
   
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 1.5*info_helvb14->ascent, string1, strlen(string1));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 3*info_helvb14->ascent, string2, strlen(string2));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 4.5*info_helvb14->ascent, string3, strlen(string3));
                  XDrawString(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC,
                     180/2 - twidth/2, 6*info_helvb14->ascent, string4, strlen(string4));

                  repaint(fullresolution_image, 0, 0, 180, 7*info_helvb14->ascent, 0, 0);
  
                  select_button(buttons, num_buttons, button_labels, "REGION AVG", colors[BLUEINDEX].pixel);

                  do{
                     XNextEvent(theDisplay, &theEvent);
                     /*************************************************************************
                     * Handle expose events by repainting the screen.
                     *************************************************************************/
                     if(theEvent.type == Expose){
                        if(theEvent.xbutton.window == overview_image.theDrawWindow)
                           repaint(overview_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height, theEvent.xexpose.x,
                                   theEvent.xexpose.y);
                        if(theEvent.xbutton.window == fullresolution_image.theDrawWindow)
                           repaint(fullresolution_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height, theEvent.xexpose.x,
                                   theEvent.xexpose.y);
                        if(theEvent.xbutton.window == histogram_image.theDrawWindow)
                           repaint(histogram_image, theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height, theEvent.xexpose.x,
                                   theEvent.xexpose.y);
                        for(i=0;i<num_buttons;i++){
                           if(theEvent.xbutton.window == buttons[i].theDrawWindow)
                              repaint(buttons[i], theEvent.xexpose.x, theEvent.xexpose.y,
                                   theEvent.xexpose.width, theEvent.xexpose.height, theEvent.xexpose.x,
                                   theEvent.xexpose.y);
                        }
                     }
                  }while(theEvent.type != ButtonPress);

                  XPutImage(theDisplay, fullresolution_image.thePixmap, fullresolution_image.theGC, tmpximage,
		     0, 0, 0, 0, 180, 7*info_helvb14->ascent);
                  XDestroyImage(tmpximage);
                  repaint(fullresolution_image, 0, 0, 180, 7*info_helvb14->ascent, 0, 0);
               }

               /****************************************************************
               * If the user already pressed the button to save each set of
               * points selected, open the file in append mode.
               ****************************************************************/
	       if(save_points_toggle == 1){
                  fclose(fp_points);
               }
            }
         }

         /**********************************************************************
         * Button presses in the histogram window are processed here. These
         * are used to remap the image displayed in the overview and the
         * fullresolution windows. Note that the fullresolution window could
         * contain fullresolution or aggregated image data.
         **********************************************************************/
         if(theEvent.xbutton.window == histogram_image.theDrawWindow){
	    if((theEvent.xbutton.button == LEFT_BUTTON) || (theEvent.xbutton.button == RIGHT_BUTTON)){

               XDefineCursor(theDisplay, overview_image.theDrawWindow, watch_cursor);
               XDefineCursor(theDisplay, fullresolution_image.theDrawWindow, watch_cursor);
               XDefineCursor(theDisplay, histogram_image.theDrawWindow, watch_cursor);
               XFlush(theDisplay);

               scale_min = -1;
	       scale_max = -1;

               if(theEvent.xbutton.button == LEFT_BUTTON) scale_min = (int)(theEvent.xbutton.x);
               if(theEvent.xbutton.button == RIGHT_BUTTON) scale_max = (int)(theEvent.xbutton.x);

               update_fullresolution_image(&image, &overview_image, &fullresolution_image, &histogram_image,
		  overview_x, overview_y, scale_min, scale_max, fullresolution_data, 0, num_detections, detections,
                  isvalid, overlay_data);

               XDefineCursor(theDisplay, overview_image.theDrawWindow, crosshair_cursor);
               XDefineCursor(theDisplay, fullresolution_image.theDrawWindow, crosshair_cursor);
               XDefineCursor(theDisplay, histogram_image.theDrawWindow, arrow_cursor);
               XFlush(theDisplay);
            }

	    else if(theEvent.xbutton.button == MIDDLE_BUTTON){

               XDefineCursor(theDisplay, overview_image.theDrawWindow, watch_cursor);
               XDefineCursor(theDisplay, fullresolution_image.theDrawWindow, watch_cursor);
               XDefineCursor(theDisplay, histogram_image.theDrawWindow, watch_cursor);
               XFlush(theDisplay);

               scale_min = -1;
	       scale_max = -1;

               rescale_overview_image(&overview_image, overview_data, num_detections, detections,
                  isvalid, overlay_data);

               XDefineCursor(theDisplay, overview_image.theDrawWindow, crosshair_cursor);
               XDefineCursor(theDisplay, fullresolution_image.theDrawWindow, crosshair_cursor);
               XDefineCursor(theDisplay, histogram_image.theDrawWindow, arrow_cursor);
               XFlush(theDisplay);
            }
	 }

         /**********************************************************************
         * Button presses in the overview window are processed here.
         **********************************************************************/
         if(theEvent.xbutton.window == overview_image.theDrawWindow){

            /*******************************************************************
            * Process a different region of the image and display the masks
            * to the screen.
            *******************************************************************/
            if(theEvent.xbutton.button == LEFT_BUTTON){

               XDefineCursor(theDisplay, overview_image.theDrawWindow, watch_cursor);
               XDefineCursor(theDisplay, fullresolution_image.theDrawWindow, watch_cursor);
               XDefineCursor(theDisplay, histogram_image.theDrawWindow, watch_cursor);
               XFlush(theDisplay);

               overview_x = (int)((float)theEvent.xbutton.x * (float)(image.cols) / (float)(overview_image.cols));
               overview_y = (int)((float)theEvent.xbutton.y * (float)(image.rows) / (float)(overview_image.rows));

               scale_min = -1;
	       scale_max = -1;

               update_fullresolution_image(&image, &overview_image, &fullresolution_image, &histogram_image,
		  overview_x, overview_y, scale_min, scale_max, fullresolution_data, 1, num_detections, detections,
                  isvalid, overlay_data);

               XDefineCursor(theDisplay, overview_image.theDrawWindow, crosshair_cursor);
               XDefineCursor(theDisplay, fullresolution_image.theDrawWindow, crosshair_cursor);
               XDefineCursor(theDisplay, histogram_image.theDrawWindow, arrow_cursor);
               XFlush(theDisplay);

            }
         }
      }
   }while( !((theEvent.xbutton.window == buttons[num_buttons-1].theDrawWindow) && (theEvent.type == ButtonPress)));
 
   XCloseDisplay(theDisplay);
   ru_exit();    /* Exit and delete the files that we uncompressed. */
}

/*******************************************************************************
* Procedure: find_which_button
* Purpose: This procedure returns the string in the button that was pressed.
* If no button was pressed then NUMM is returned.
* Name: Michael Heath, University of South Florida
*******************************************************************************/
char *find_which_button(XEvent theEvent, IMAGEXWIN *buttons, int num_buttons, char **button_labels)
{
   int i;
   char *ptr = NULL;

   for(i=0;i<num_buttons;i++){
      if(theEvent.xbutton.window == buttons[i].theDrawWindow){
         return(button_labels[i]);
      }
   }

   return(ptr);
}

/*******************************************************************************
* Procedure: select_button
* Purpose: To draw a blue border around a window to show it is selected.
* Name: Michael Heath, University of South Florida
* Date: 11/16/97
*******************************************************************************/
void select_button(IMAGEXWIN *buttons, int num_buttons, char *button_labels[], char *alabel, unsigned long int color)
{
   int which_button;

   for(which_button=0;which_button<num_buttons;which_button++){
      if(strcmp(alabel, button_labels[which_button]) == 0){

         XSetForeground(theDisplay, buttons[which_button].theGC, color); /* Blue */
         XDrawRectangle(theDisplay, buttons[which_button].thePixmap, buttons[which_button].theGC, 0, 0,
            buttons[which_button].cols-1, buttons[which_button].rows-1);

         repaint(buttons[which_button], 0, 0, buttons[which_button].cols, buttons[which_button].rows, 0, 0);
      }
   }
}

/*******************************************************************************
* Procedure: unselect_button
* Purpose: To draw a black border around a window to show it is unselected.
* Name: Michael Heath, University of South Florida
* Date: 11/16/97
*******************************************************************************/
void unselect_button(IMAGEXWIN *buttons, int num_buttons, char *button_labels[], char *alabel, unsigned long int color)
{
   int which_button;

   for(which_button=0;which_button<num_buttons;which_button++){
      if(strcmp(alabel, button_labels[which_button]) == 0){
         XSetForeground(theDisplay, buttons[which_button].theGC, color); /* Black */
         XDrawRectangle(theDisplay, buttons[which_button].thePixmap, buttons[which_button].theGC, 0, 0,
            buttons[which_button].cols-1, buttons[which_button].rows-1);
   
         repaint(buttons[which_button], 0, 0, buttons[which_button].cols, buttons[which_button].rows, 0, 0);
      }
   }
}

/*******************************************************************************
* Procedure: setup_buttons
*******************************************************************************/
void setup_buttons(IMAGEXWIN *imagemenu_image, int button_width,
   IMAGEXWIN *buttons, int num_buttons, char **button_labels, char *window_title,
   int row_position, int min_height, int min_width)
{
   int i, button_height = (1.5 * (info_helvb14->ascent + info_helvb14->descent));
   int actual_width, actual_height;

   actual_width = button_width;
   if(actual_width < min_width) actual_width = min_width;

   actual_height = num_buttons*button_height;
   if(actual_height < min_height) actual_height = min_height;

   /****************************************************************************
   * Set up the imagemenu image window.
   ****************************************************************************/
   imagemenu_image->theBorderWindow = openWindow(0, row_position,
      actual_width, actual_height, NORMAL_WINDOW, window_title, theRootWindow,
      &(imagemenu_image->theBorderGC));
   initEvents(imagemenu_image->theBorderWindow, IN_PALETTE);

   imagemenu_image->cols = actual_width;
   imagemenu_image->rows = actual_height;

   for(i=0;i<num_buttons;i++){

      buttons[i] = *imagemenu_image;

      buttons[i].have_a_pixmap = 0;

      buttons[i].theDrawWindow = openWindow(0, i*button_height,
         button_width, button_height, NORMAL_WINDOW, button_labels[i], buttons[i].theBorderWindow,
         &(buttons[i].theGC));
      initEvents(buttons[i].theDrawWindow, NOT_IN_PALETTE);

      buttons[i].rows = button_height;
      buttons[i].cols = button_width;

      if((buttons[i].image = (unsigned char *) calloc(3*button_height*button_width,
         sizeof(unsigned char))) == NULL){
         fprintf(stderr, "Malloc error.\n");
      }
      if((buttons[i].image_dsp = (unsigned char *) calloc(3*button_height*button_width,
         sizeof(unsigned char))) == NULL){
         fprintf(stderr, "Malloc error.\n");
      }
      buttons[i].maskimage = NULL;
      buttons[i].theImage = NULL;

      display_image(&(buttons[i]), button_height, button_width, 1);

      XSetForeground(theDisplay, buttons[i].theGC, cell[128]);
      XFillRectangle(theDisplay, buttons[i].thePixmap, buttons[i].theGC, 0, 0, button_width, button_height);

      unselect_button(buttons, num_buttons, button_labels, button_labels[i], cell[0]);

      XSetForeground(theDisplay, buttons[i].theGC, cell[196]); /* Lighter Gray */
      XDrawLine(theDisplay, buttons[i].thePixmap, buttons[i].theGC, 3, 3, buttons[i].cols-3, 3);
      XDrawLine(theDisplay, buttons[i].thePixmap, buttons[i].theGC, 3, 3, 3, buttons[i].rows-3);

      XSetForeground(theDisplay, buttons[i].theGC, cell[0]);
      XDrawLine(theDisplay, buttons[i].thePixmap, buttons[i].theGC,
	 buttons[i].cols-3, 3, buttons[i].cols-3, buttons[i].rows-3);
      XDrawLine(theDisplay, buttons[i].thePixmap, buttons[i].theGC,
	 3, buttons[i].rows-3, buttons[i].cols-3, buttons[i].rows-3);

      XSetFont(theDisplay, buttons[i].theGC, info_helvb14->fid);

      XDrawString(theDisplay, buttons[i].thePixmap, buttons[i].theGC,
         button_width/2 - XTextWidth(info_helvb14, button_labels[i], strlen(button_labels[i]))/2,
         buttons[i].rows - 2*info_helvb14->descent, button_labels[i], strlen(button_labels[i]));

      repaint(buttons[i], 0, 0, buttons[i].cols, buttons[i].rows, 0, 0);
   }
}

/*******************************************************************************
* Procedure: get_command_line_arguments
* Purpose: To extract the values from the command line and fill them into the
* correct variables. There are a lot of arguments so this procedure is lengthy.
* If the correct arguments are not supplied then help is printed to the user.
* Name: Michael Heath, University of South Florida
* Date: 7/17/97
*******************************************************************************/
void get_command_line_arguments(int argc, char *argv[], char *filename,
    int *rows, int *cols, int *bytesperpix, int *headerbytes, int *swap_bytes,
    int *fullresolution_height, int *fullresolution_width, char *scanner,
    int *startup_mode, char *view, char *direction, char **detection_filename,
    char **segmentation_filename,
    int *detmax, float *detthresh, char *lesiontype, char **ms, char **mm,
    char **ct, char **cd, char **pathology)
{
   FILE *fp=NULL;
   int i, b;
   int incomplete_input;
   int pixelskip, lineskip;
   ICSDATA icsdata;
   char thepath[100];
   char tempstring[200];
   char scanner_preference[100];
   char buf[100];
   int is_ISMD(char *filename);
   char *ics_filename = NULL;
   char tmpfilename[200];

   /****************************************************************************
   * Get the command line arguments. Make sure that the necessary arguments are
   * supplied on the command line. If any are not specified, inform the user
   * of the missing arguments.
   ****************************************************************************/
   if(argc == 1){
      printf("\n\n%s -i filename | -ics file.ics -view VIEW [-inp #pixels]\n", argv[0]);
      printf("             [-inl #lines] [-bpp #bytesperpixel] [-hb #headerbytes]\n");
      printf("             [-swap_bytes] [-sr sample_rate] [-w width] [-h height]\n");
      printf("             [-view whichview] [-direction DIRval] [-intensity]\n");
      printf("             [-density] [-verbose] [-ws #bits] [-manywin] [-version]\n");
      printf("             [-segmentation file.sgt] [-detection file.det max#det suspthresh]\n");
      printf("             [-l MASS | CALCIFICATION] [-pathology MALIGNANT | BENIGN]\n");
      printf("\n");
      printf("   Mandatory Arguments:\n");
      printf("   -------------------\n");
      printf("   -i filename         = Input image (filename.dba, filename.ics or filename.pgm).\n");
      printf("       --- or ---\n");
      printf("   -ics filename -view VIEW = List the ics filename and view may be any of\n");
      printf("                 (RIGHT_MLO, RIGHT_CC, LEFT_MLO or LEFT_CC)\n");
      printf("\n");
      printf("   Optional Arguments:\n");
      printf("   -------------------\n");
      printf("   -inp #pixels        = Number of pixels per line of the image.\n");
      printf("   -inl #lines         = Number of lines in the image.\n");
      printf("   -bpp #bytesperpixel = Number of bytes per pixel (1 or 2, default is 1). \n");
      printf("   -hb #headerbytes    = Number of header bytes in the image (default is 0).\n");
      printf("   -swap_bytes         = Swap the bytes when reading from or writing to the file.\n");
      printf("   -sr sample_rate     = Pixel sample rate in microns/pixel.\n");
      printf("   -w width            = Width of the fullresolution window.\n");
      printf("   -h height           = Height of the fullresolution window.\n");
      printf("   -view whichview     = LEFT_CC, LEFT_MLO, RIGHT_CC or RIGHT_MLO.\n");
      printf("   -direction DIRval   = Viewing film looking \"TOWARD_SOURCE\" or \"AWAY_FROM_SOURCE\"\n");
      printf("   -intensity          = Work in intensity.\n");
      printf("   -density            = Work in scaled density (DBA, LUMISYS, HOWTEK or HOWTEK_ISMD).\n");
      printf("   -verbose            = Print out information for debugging.\n");
      printf("   -ws #bits           = Workstation bits (8 or 24).\n");
      printf("   -manywin            = Start the program displaying many windows.\n");
      printf("   -version            = Print out the version number and exit.\n");
      printf("   -ics FILE -view VIEW = Allows the ground truth to be used when -i is used.\n");
      printf("   -segmentation file  = Allows the breast segmentation and axis to be displayed.\n");
      printf("   -detection file # # = Allows the detections to be displayed on the image. The\n");
      printf("                         first # is the max detections to use and the second is\n");
      printf("                         the suspiciousness threshold to apply.\n");
      printf("\n\n\n");
      exit(1);
   }

   /****************************************************************************
   * Get the command line arguments. Make sure that the necessary arguments are
   * supplied on the command line. If any are not specified, inform the user
   * of the missing arguments.
   ****************************************************************************/
   filename[0] = '\0';
   *rows = -1;
   *cols = -1;
   *bytesperpix = 1;
   *headerbytes = 0;
   *swap_bytes = FALSE;
   workstation_bits = 24;
   view[0] = '\0';
   direction[0] = '\0';
   sprintf(scanner_preference, "none");
   *detection_filename = NULL;
   *segmentation_filename = NULL;
   *detmax = 0;
   *detthresh = 0.0;
   strcpy(lesiontype, "NONE");
   *pathology = NULL;
   *ms = NULL;
   *mm = NULL;
   *ct = NULL;
   *cd = NULL;

   for(i=1;i<argc;i++){
      if(strcmp(argv[i],"-i")==0){ sprintf(filename, "%s", argv[i+1]); i++;}
      else if(strcmp(argv[i],"-miasgt")==0){ MIASgt_filename = argv[i+1]; i++;}
      else if(strcmp(argv[i],"-inl")==0){ *rows = atoi(argv[i+1]); i++;}
      else if(strcmp(argv[i],"-inp")==0){ *cols = atoi(argv[i+1]); i++;}
      else if(strcmp(argv[i],"-bpp")==0){ *bytesperpix = atoi(argv[i+1]); i++;}
      else if(strcmp(argv[i],"-hb")==0){ *headerbytes = atoi(argv[i+1]); i++;}
      else if(strcmp(argv[i],"-sr")==0){ samplerate_in_microns = atof(argv[i+1]); i++;}
      else if(strcmp(argv[i],"-verbose")==0){ verbosemode = 1;}
      else if(strcmp(argv[i],"-swap_bytes")==0){ *swap_bytes = TRUE;}
      else if(strcmp(argv[i],"-w")==0){ *fullresolution_width = atoi(argv[i+1]); i++;}
      else if(strcmp(argv[i],"-h")==0){ *fullresolution_height = atoi(argv[i+1]); i++;}
      else if(strcmp(argv[i],"-view")==0){ strcpy(view, argv[i+1]); i++;}
      else if(strcmp(argv[i],"-view")==0){ strcpy(direction, argv[i+1]); i++;}
      else if(strcmp(argv[i],"-version")==0){ printf("\n\n MammoView version: 6.21.2000\n"); exit(1); }
      else if(strcmp(argv[i],"-DBA")==0){ sprintf(scanner_preference, "DBA"); }
      else if(strcmp(argv[i],"-LUMISYS")==0){ sprintf(scanner_preference, "LUMISYS"); }
      else if(strcmp(argv[i],"-HOWTEK")==0){ sprintf(scanner_preference, "HOWTEK"); }
      else if(strcmp(argv[i],"-HOWTEK_ISMD")==0){ sprintf(scanner_preference, "HOWTEK_ISMD"); }
      else if(strcmp(argv[i], "-detection") == 0){
         if((i+3) < argc){
            *detection_filename = argv[i+1]; i++;
            *detmax = atoi(argv[i+1]); i++;
            *detthresh = atof(argv[i+1]); i++;
         }
      }
      else if(strcmp(argv[i], "-segmentation") == 0){
         if((i+1) < argc){
            *segmentation_filename = argv[i+1]; i++;
         }
      }
      else if(strcmp(argv[i], "-ics") == 0){
         if((i+1) < argc){
            ics_filename = argv[i+1]; i++;
         }
      }
      else if(strcmp(argv[i],"-intensity")==0){ mapping = INTENSITY_MAPPING;}
      else if(strcmp(argv[i],"-density")==0){ mapping = DENSITY_MAPPING;}
      else if(strcmp(argv[i],"-ws")==0){ workstation_bits = atoi(argv[i+1]); i++;}
      else if(strcmp(argv[i],"-manywin")==0){ *startup_mode = 1;}
      else if(strcmp(argv[i], "-pathology") == 0){
         if((i+1) < argc){
            if(strcmp(argv[i+1], "MALIGNANT") == 0)
               *pathology = argv[i+1];
            else if(strcmp(argv[i+1], "BENIGN") == 0)
               *pathology = argv[i+1];
            else{
               fprintf(stderr, "MammoView: Invalid pathology (MALIGNANT or BENIGN).\n");
               return;
            }
         }
         else{
            fprintf(stderr, "MammoView: Error! No pathology was specified.\n");
            return;
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
               fprintf(stderr, "MammoView: Invalid lesiontype (CALCIFICATION or MASS).\n");
               return;
            }
         }
         else{
            fprintf(stderr, "MammoView: Error! No lesiontype was specified.\n");
            return;
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
   }

   if(strcmp(lesiontype, "NONE") == 0) strcpy(lesiontype, "BOTH");

   /****************************************************************************
   * Make sure that the number of workstation bits is valid.
   ****************************************************************************/
   if(!((workstation_bits == 24) || (workstation_bits == 8))){
      fprintf(stderr, "Invalid value for ws.\n");
      exit(1);

   }

   /****************************************************************************
   * Try to interpret the image as a raw PGM image.
   ****************************************************************************/
   if(filename[0] != '\0'){
      if((fp = fopen(filename, "rb")) == NULL){
         fprintf(stderr, "Error opening the file %s\n", filename);
         exit(1);
      }

      fgets(buf, 90, fp);
      if(strncmp(buf,"P5",2) == 0){

         printf("Recognized a PGM image. Reading the data from the header...\n");

         if(strcmp(scanner_preference, "none") != 0){
            sprintf(scanner, "%s", scanner_preference);
         }

         do{
            fgets(buf, 90, fp);
            if(strncmp(buf,"# RESOLUTION = ", 15) == 0)
               sscanf(buf, "# RESOLUTION = %lf", &samplerate_in_microns);
            if(strncmp(buf,"# SCANNER = ", 12) == 0)
               sscanf(buf, "# SCANNER = %s", scanner);
         }while(buf[0] == '#');  /* skip all comment lines */

         sscanf(buf, "%d %d", cols, rows);

         do{
            fgets(buf, 90, fp);
            if(strncmp(buf,"# RESOLUTION = ", 15) == 0)
               sscanf(buf, "# RESOLUTION = %lf", &samplerate_in_microns);
            if(strncmp(buf,"# SCANNER = ", 12) == 0)
               sscanf(buf, "# SCANNER = %s", scanner);
         }while(buf[0] == '#');  /* skip all comment lines */

         if(strncmp(buf, "255", 3) == 0){
            *bytesperpix = 1;
            *swap_bytes = 0;
            pixelskip = 0;
            lineskip = 0;
            fseek(fp, 0, 2);
            *headerbytes = ftell(fp);
            *headerbytes -= (*rows)*(*cols);
         }
         else{
            *bytesperpix = 2;
            *swap_bytes = 0;
            pixelskip = 0;
            lineskip = 0;
            fseek(fp, 0, 2);
            *headerbytes = ftell(fp);
            *headerbytes -= (*rows)*(*cols)*2;
         }
      }

      if(fp != stdin) fclose(fp);

      /*************************************************************************
      * If the image filename ends in ".dba" then it is considered to be an
      * image from the DBA scanner and the image attributes are gotten from the
      * header.
      *************************************************************************/
      if(strcmp(&(filename[strlen(filename)-4]), ".dba") == 0){
         printf("Recognized a DBA image. Reading the data from the header...\n");

         sprintf(scanner, "DBA");

         readDBAImage_header(filename, rows, cols, bytesperpix,
            headerbytes, swap_bytes, &pixelskip, &lineskip);

         /**********************************************************************
         * We do not actually have the sample rate so we will make an educated
         * guess. We have received data scanned at 21, 42 and 210 microns from this
         * scanner. Images larger than 100 million bytes are 21 micron, images
         * less than 5 million bytes are 210 micron and all others are 42 micron.
         **********************************************************************/
         if((*rows)*(*cols)*2 >= 100000000) samplerate_in_microns = 21.0;
         else if((*rows)*(*cols)*2 >= 5000000) samplerate_in_microns = 42.0;
         else samplerate_in_microns = 210.0;
      }

      else if(strstr(filename, ".dcm") != NULL){

         if(read_howteck_image_header(filename, rows, cols,
            headerbytes, swap_bytes) == 0){
            fprintf(stderr, "Error reading the HowTeck image header.\n");
            exit(1);
         }

         *bytesperpix = 2;
         pixelskip = 0;
         lineskip = 0;

         /*************************************************************************
         * We again make an educated guess about the scanning rate. The HowTeck
         * scanner should be operated at 43.5 microns for clean films and at
         * perhaps 5 times that for marked films.
         *************************************************************************/
         if(is_ISMD(filename) == 1){
            sprintf(scanner, "HOWTEK_ISMD");
            if((*rows)*(*cols)*2 < 5000000) samplerate_in_microns = 174.0;
            else samplerate_in_microns = 43.5;  /* Assume 43.5 */
         }

         else{
            sprintf(scanner, "HOWTEK");
            if((*rows)*(*cols)*2 < 5000000) samplerate_in_microns = 217.5;
            else samplerate_in_microns = 43.5;  /* Assume 43.5 */
         }
      }

      else if((strstr(filename, "small") != NULL) || (strstr(filename, "full") != NULL)){

         sprintf(scanner, "LUMISYS");

         readIFSImage_header(filename, rows, cols, bytesperpix,
            headerbytes, swap_bytes, &pixelskip, &lineskip);

         if(strstr(filename, "full") != NULL) samplerate_in_microns = 50.0;
         else if(strstr(filename, "small") != NULL) samplerate_in_microns = 200.0;
      }
   }

   /****************************************************************************
   * If an ".ics" filename was specified then we get the data about the image
   * from the "ics" file.
   ****************************************************************************/
   if(ics_filename != NULL){
      printf("Recognized an ics data file. Reading the data from it...\n");

      strcpy(tmpfilename, ics_filename);

      /*************************************************************************
      * Separate out the filename into the path and the filename.
      *************************************************************************/
      b = strlen(tmpfilename) - 1;
      while((b > 0) && (tmpfilename[b] != '/')) b--;

      if(b == 0){
         sprintf(thepath, "./");

	 if(verbosemode){
            printf("No path was specified, the ics file should be located in this directory.\n");
            printf("PATH: %s\n", thepath);
            printf("ICS: %s\n", tmpfilename);
         }
      }
      else if(b < (strlen(tmpfilename) - 1)){
         sprintf(tempstring, "%s", tmpfilename);

         memset(thepath, 0, 100);
         strncpy(thepath, tmpfilename, b+1);

         sprintf(tmpfilename, "%s", tempstring+(b+1));

	 if(verbosemode){
            printf("A path was specified.\n");
            printf("PATH: %s\n", thepath);
            printf("ICS: %s\n", tmpfilename);
         }
      }

      if(read_ics_file(tmpfilename, &icsdata, thepath) == 0){
         exit(1);
      }

      if(view[0] == '\0'){
         fprintf(stderr, "ERROR! No view specified for the ICS file.\n");
         exit(1);
      }

      if(filename[0] == '\0'){

         sprintf(scanner, "%s", icsdata.digitizer_brand);

         if(strcmp(view, "LEFT_CC") == 0){
	    sprintf(filename, "%s", icsdata.left_cc.uncompressed_filename);
            *rows = icsdata.left_cc.rows;
            *cols = icsdata.left_cc.cols;
            samplerate_in_microns = icsdata.left_cc.resolution;
            if(icsdata.left_cc.overlay_exists) sprintf(overlay_filename, "%s", icsdata.left_cc.overlay_filename);
         }
         else if(strcmp(view, "LEFT_MLO") == 0){
	    sprintf(filename, "%s", icsdata.left_mlo.uncompressed_filename);
            *rows = icsdata.left_mlo.rows;
            *cols = icsdata.left_mlo.cols;
            samplerate_in_microns = icsdata.left_mlo.resolution;
            if(icsdata.left_mlo.overlay_exists) sprintf(overlay_filename, "%s", icsdata.left_mlo.overlay_filename);
         }
         else if(strcmp(view, "RIGHT_CC") == 0){
	    sprintf(filename, "%s", icsdata.right_cc.uncompressed_filename);
            *rows = icsdata.right_cc.rows;
            *cols = icsdata.right_cc.cols;
            samplerate_in_microns = icsdata.right_cc.resolution;
            if(icsdata.right_cc.overlay_exists) sprintf(overlay_filename, "%s", icsdata.right_cc.overlay_filename);
         }
         else if(strcmp(view, "RIGHT_MLO") == 0){
	    sprintf(filename, "%s", icsdata.right_mlo.uncompressed_filename);
            *rows = icsdata.right_mlo.rows;
            *cols = icsdata.right_mlo.cols;
            samplerate_in_microns = icsdata.right_mlo.resolution;
            if(icsdata.right_mlo.overlay_exists) sprintf(overlay_filename, "%s", icsdata.right_mlo.overlay_filename);
         }
         else{
            fprintf(stderr, "ERROR! Invalid view (%s) specified for the ICS file.\n", view);
            exit(1);
         }

         *bytesperpix = 2;
         *headerbytes = 0;
         *swap_bytes = FALSE;

         /* printf("VIEW = %s\n", view); */
         if(decompress_ics_image(&icsdata, view) == 0){
            fprintf(stderr, "Error decompressing the image.\n");
            exit(1);
         }
         base_subfactor = 1.0;
      }

      /*************************************************************************
      * Since a separate image filename was specified, we just want to find
      * if an "OVERLAY" file exists and what the resolution difference is
      * between the specified image and the image in which the OVERLAY files
      * ground truth information is specified. This will allow us to rescale
      * the overlay data to a new resoltution so we can overlay it on the
      * image we will be viewing.
      *************************************************************************/
      else{

         if(strcmp(view, "LEFT_CC") == 0){
            base_subfactor = (float)icsdata.left_cc.rows / (float)(*rows);
            samplerate_in_microns = icsdata.left_cc.resolution * base_subfactor;
            if(icsdata.left_cc.overlay_exists) sprintf(overlay_filename, "%s", icsdata.left_cc.overlay_filename);
         }
         else if(strcmp(view, "LEFT_MLO") == 0){
            base_subfactor = (float)icsdata.left_mlo.rows / (float)(*rows);
            samplerate_in_microns = icsdata.left_mlo.resolution * base_subfactor;
            if(icsdata.left_mlo.overlay_exists) sprintf(overlay_filename, "%s", icsdata.left_mlo.overlay_filename);
         }
         else if(strcmp(view, "RIGHT_CC") == 0){
            base_subfactor = (float)icsdata.right_cc.rows / (float)(*rows);
            samplerate_in_microns = icsdata.right_cc.resolution * base_subfactor;
            if(icsdata.right_cc.overlay_exists) sprintf(overlay_filename, "%s", icsdata.right_cc.overlay_filename);
         }
         else if(strcmp(view, "RIGHT_MLO") == 0){
            base_subfactor = (float)icsdata.right_mlo.rows / (float)(*rows);
            samplerate_in_microns = icsdata.right_mlo.resolution * base_subfactor;
            if(icsdata.right_mlo.overlay_exists) sprintf(overlay_filename, "%s", icsdata.right_mlo.overlay_filename);
         }
         else{
            fprintf(stderr, "ERROR! Invalid view (%s) specified for the ICS file.\n", view);
            exit(1);
         }
      }
   }

   /****************************************************************************
   * If the user left out any mandatory command line options, inform them of
   * the parameters they need to enter.
   ****************************************************************************/
   incomplete_input = 0;     /* Assume all user input was supplied. */

   if(filename[0] == '\0'){
      fprintf(stderr, "You must specify the input filename.\n");
      fprintf(stderr, "   -i filename\n");
      incomplete_input = 1;
   }
   if((*cols) == (-1)){
      fprintf(stderr, "You must specify the number of pixels per line.");
      fprintf(stderr, "   -inp #pixels\n");
      incomplete_input = 1;
   }
   if((*rows) == (-1)){
      fprintf(stderr, "You must specify the number of image lines.");
      fprintf(stderr, "   -inl #lines\n");
      incomplete_input = 1;
   }

   if(incomplete_input == 1){
      fprintf(stderr, "\n\n");
      exit(1);
   }

   if(!(((*bytesperpix)==1) || ((*bytesperpix)==2))){
      fprintf(stderr, "\n\nThe number of bytes per pixel must be 1 or 2.\n\n");
      exit(1);
   }
}

void setup_histogram_window(IMAGEXWIN *histogram_image, int hist_rows, int hist_cols, char *window_title)
{
   /****************************************************************************
   * Set up the histogram image window.
   ****************************************************************************/
   histogram_image->theBorderWindow = openWindow(0, MaxScreenRows-hist_rows-34,
      hist_cols, hist_rows, NORMAL_WINDOW, window_title, theRootWindow,
      &(histogram_image->theBorderGC));
   initEvents(histogram_image->theBorderWindow, IN_PALETTE);

   histogram_image->theDrawWindow = openWindow(0, 0, hist_cols, hist_rows,
      NORMAL_WINDOW, "Histogram", histogram_image->theBorderWindow,
      &(histogram_image->theGC));
   initEvents(histogram_image->theDrawWindow, NOT_IN_PALETTE );

   XSetFont(theDisplay, histogram_image->theGC, info_helvb14->fid);

   histogram_image->rows = hist_rows;
   histogram_image->cols = hist_cols;

   histogram_image->maskimage = NULL;
   histogram_image->theImage = NULL;
}

/*******************************************************************************
* Procedure: setup_fullresolution_window
* Purpose: This procedure sets up the fullresolution window.
* thresholded images.
* Name: Michael Heath, University of South Florida
* Date: 7/17/97
*******************************************************************************/
void setup_fullresolution_window(CACHEIM *image, IMAGEXWIN *overview_image,
   IMAGEXWIN *fullresolution_image, IMAGEXWIN *histogram_image,
   int sub_rows, int sub_cols, USHORT **fullresolution_data, char *window_title)
{
   /****************************************************************************
   * Since we are not zooming at all the fullresolution image should not be
   * larger than the image itself.
   ****************************************************************************/
   if(sub_rows > image->rows) sub_rows = image->rows;
   if(sub_cols > image->cols) sub_cols = image->cols;

   /****************************************************************************
   * Set up the fullresolution image window.
   ****************************************************************************/
   fullresolution_image->theBorderWindow = openWindow(MaxScreenCols-sub_cols-18, 0,
      sub_cols, sub_rows, NORMAL_WINDOW, window_title, theRootWindow,
      &(fullresolution_image->theBorderGC));
   initEvents(fullresolution_image->theBorderWindow, IN_PALETTE);

   fullresolution_image->theDrawWindow = openWindow(0, 0, sub_cols, sub_rows,
      NORMAL_WINDOW, "FullResolution", fullresolution_image->theBorderWindow,
      &(fullresolution_image->theGC));
   initEvents(fullresolution_image->theDrawWindow, NOT_IN_PALETTE );

   XSetFont(theDisplay, fullresolution_image->theGC, info_helvb14->fid);

   fullresolution_image->rows = sub_rows;
   fullresolution_image->cols = sub_cols;

/*
   *histogram_image = *fullresolution_image;

   histogram_image->theDrawWindow = openWindow(0, sub_rows, sub_cols, HIST_HEIGHT,
      NORMAL_WINDOW, "Histogram", histogram_image->theBorderWindow,
      &(histogram_image->theGC));
   initEvents(histogram_image->theDrawWindow, NOT_IN_PALETTE );

   XSetFont(theDisplay, histogram_image->theGC, info_helvb14->fid);

   histogram_image->rows = HIST_HEIGHT;
   histogram_image->cols = sub_cols;
*/

   /****************************************************************************
   * Allocate some memory for these images. This memory will be filled with
   * R,G,B values that are pixel interleaved.
   ****************************************************************************/
   if((fullresolution_image->image = (unsigned char *) calloc(sub_rows*sub_cols,
      sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Malloc error.\n");
   }
   if((fullresolution_image->image_dsp = (unsigned char *) calloc(sub_rows*sub_cols,
      sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Malloc error.\n");
   }
   fullresolution_image->maskimage = NULL;
   fullresolution_image->theImage = NULL;

   if((*fullresolution_data = (USHORT *) calloc(sub_rows*sub_cols, sizeof(USHORT))) == NULL){
      fprintf(stderr, "Malloc error.\n");
   }

   /****************************************************************************
   * Fill the fullresolution image with data.
   ****************************************************************************/
/*
   update_fullresolution_image(image, overview_image, fullresolution_image, histogram_image, 0, 0, -1, -1,
      *fullresolution_data, 1, num_detections, detections, isvalid, overlay_data);
*/
}

/*******************************************************************************
* Procedure: update_fullresolution_image
* Purpose: This procedure fills the fullresolution image with new data or just
* rescales the data that is already there.
* Name: Michael Heath, University of South Florida
* Date: 11/8/97
*******************************************************************************/
void update_fullresolution_image(CACHEIM *image, IMAGEXWIN *overview_image, IMAGEXWIN *fullresolution_image,
   IMAGEXWIN *histogram_image, int x, int y, int scale_min, int scale_max, USHORT *fullresolution_data,
   int acquire_new_data, int num_detections, REGION *detections, int *isvalid, OVERLAY_DATA overlay_data)
{
   int rows, cols;
   int val, r, c, pos;
   unsigned long int *hist = NULL;
   int bottom_of_range, top_of_range;
   int lm = 20, rm = 5, tm = 5, bm = 20;
   unsigned char val_uchar, *tmp_ptr=NULL;
   char thestring[100];
   static int text_top=0, text_left=0, text_width=0, text_height=0;
   float overview_scale, scale_width, hold_scale_width;
   int scale_width_pixels = 0;
   static int working_rows=0, working_cols=0;
   static int source_x, source_y;
   static int min=0, max=0;
   static int min_val, max_val;
   static int bin_width;
   double tmpval2square;

   unsigned long int background_color = cell[196];
   unsigned long int graph_color = colors[BLUEINDEX].pixel;
   unsigned long int number_color = cell[0];
   unsigned long int scalerange_color = cell[253];

   rows = fullresolution_image->rows;
   cols = fullresolution_image->cols;

   /****************************************************************************
   * The image might be smaller than the fullresolution image.
   ****************************************************************************/
   working_rows = ((rows > image->rows) ? image->rows : rows);
   working_cols = ((cols > image->cols) ? image->cols : cols);

   /****************************************************************************
   * If we are not using the same data as before, we must do a lot of work.
   ****************************************************************************/
   if(acquire_new_data){

      /*************************************************************************
      * Calculate the upper left corner of where we will grab data from the file.
      *************************************************************************/
      source_x = x - working_cols/2;
      if(source_x < 0) source_x = 0;
      if((source_x+working_cols) > (image->cols)) source_x = image->cols - working_cols;

      source_y = y - working_rows/2;
      if(source_y < 0) source_y = 0;
      if((source_y+working_rows) > (image->rows)) source_y = image->rows - working_rows;

      image_ul_x = source_x;
      image_ul_y = source_y;

      /*************************************************************************
      * Repaint the crop box on the positioning (overview) image.
      *************************************************************************/
      display_image(overview_image, overview_image->rows, overview_image->cols, 1);

      overview_scale = (float)overview_image->cols / (float)image->cols;
      box_x = (int)(source_x*overview_scale);
      box_y = (int)(source_y*overview_scale);
      box_width = (int)(cols*overview_scale);
      box_height = (int)(rows*overview_scale);

      XSetForeground(theDisplay, overview_image->theGC, colors[BLUEINDEX].pixel);  /* Blue */

      XDrawRectangle(theDisplay, overview_image->thePixmap, overview_image->theGC,
         box_x, box_y, box_width, box_height);

      if((overlay_filename[0] != '\0') && (show_overlay == 1))
         overlay_groundtruth(overview_image, 0, 0, overview_subfactor,
            isvalid, overlay_data);

      if(segmentation_toggle)
         overlay_features(overview_image, 0, 0, overview_subfactor);

      if(display_regions_toggle)
         overlay_regions(overview_image, 0, 0, overview_subfactor,
                        num_detections, detections);

      repaint(*overview_image, 0, 0, overview_image->cols, overview_image->rows, 0, 0);

      /*************************************************************************
      * Calculate statistics for scaling the image.
      *************************************************************************/
      hist = (unsigned long int *) calloc(65536, sizeof(unsigned long int));

      for(r=0;r<working_rows;r++){
	 pos = r*cols;
         for(c=0;c<working_cols;c++,pos++){
            tmpval2square = ((*(image->getpixel))(image, source_y+r, source_x+c));
            val = (int)floor(tmpval2square);
            hist[(int)val]++;
            fullresolution_data[pos] = (USHORT)val;
         }
      }

      r=0;
      c=0;
      while((r < 65535) && (c < (working_rows*working_cols/50))){
	 c += hist[r];
	 r++;
      }
      min = r;

      r=65535;
      c=0;
      while((r > 0) && (c < (working_rows*working_cols/400))){
	 c += hist[r];
	 r--;
      }
      max = r;

      /*************************************************************************
      * If the max and the min are the same, we will run into all kinds of
      * trouble with division calculations. Therefore I just add 1 to the min
      * to get the max. When the max is equal to the min we don't have any data
      * to display anyway so this should not matter.
      *************************************************************************/
      if(max == min) max = min+1;

      /*************************************************************************
      * Plot the histogram on the screen. Fit it into the histogram_image
      * Leave margin space specified by lm, rm, tm, bm. The histogram is
      * re-binned. The function returns the minimum and the maximum values in
      * the histogram and the computed bin width used in re-binning the histogram.
      *************************************************************************/
      plot_histogram(hist, 65536, lm, rm, tm, bm, &min_val, &max_val, &bin_width, histogram_image);

      free(hist);
   }
   else{
      if(scale_min >= 0) min = min_val + ((scale_min-lm)/histogram_stretching) * bin_width;
      if(scale_max >= 0) max = min_val + ((scale_max-lm)/histogram_stretching) * bin_width;
   }

   /*************************************************************************
   * Scale the image.
   *************************************************************************/
   if(min >= max){
      min = min_val;
      max = max_val;
   }
   else{
      if(min < min_val) min = min_val;
      if(max > max_val) max = max_val;
   }
   scmin = min;
   scmax = max;

   /*************************************************************************
   * Clear the scaling range from the window.
   *************************************************************************/
   XSetForeground(theDisplay, histogram_image->theGC, background_color);
   XFillRectangle(theDisplay, histogram_image->thePixmap, histogram_image->theGC, 0, histogram_image->rows-bm+2,
      histogram_image->cols, (bm/4 - 2));

   XFillRectangle(theDisplay, histogram_image->thePixmap, histogram_image->theGC,
      text_left, histogram_image->rows-bm+2, text_width, bm-2);

   /*************************************************************************
   * Color in part of the histogram window to show the scaling range.
   *************************************************************************/
   bottom_of_range = ((min - min_val) * histogram_stretching) / bin_width;
   top_of_range = ((max - min_val) * histogram_stretching) / bin_width;

   XSetForeground(theDisplay, histogram_image->theGC, scalerange_color);
   XFillRectangle(theDisplay, histogram_image->thePixmap, histogram_image->theGC,
      lm + bottom_of_range, histogram_image->rows-bm+2, top_of_range - bottom_of_range, (bm/4 - 2));

   repaint(*histogram_image, 0, 0, histogram_image->cols, histogram_image->rows, 0, 0);

   /****************************************************************************
   * Draw the numbers of the scaling range on the histogram image.
   ****************************************************************************/
   XSetForeground(theDisplay, histogram_image->theGC, number_color);  /* Black */

   if(mapping == DENSITY_MAPPING){
      sprintf(thestring, "XXXXXXXXXXXXXXXXXXXXX");
      text_width = XTextWidth(info_helvb14, thestring, strlen(thestring));

      sprintf(thestring, "(%0.3lf->0,%0.3lf->255)", (od_offset-(double)min)/od_scale+min_density,
         (od_offset-(double)max)/od_scale+min_density);
   }
   else{
      sprintf(thestring, "XXXXXXXXXXXXXXXXXXXXX");
      text_width = XTextWidth(info_helvb14, thestring, strlen(thestring));

      sprintf(thestring, "(%5d->0,%5d->255)", min, max);
   }

   text_height = info_helvb14->ascent;
   text_top = histogram_image->rows - text_height / 4;

   text_left = (histogram_image->cols/3) - text_width/2;

   XDrawString(theDisplay, histogram_image->thePixmap, histogram_image->theGC,
      text_left, histogram_image->rows - info_helvb14->ascent / 4, thestring, strlen(thestring));

   repaint(*histogram_image, 0, 0, histogram_image->cols, histogram_image->rows, 0, 0);

   /****************************************************************************
   * Scale the image.
   ****************************************************************************/
   if((rows != working_rows) || (cols != working_cols))
      memset((void *)(fullresolution_image->image), 0, rows*cols);

   for(r=0;r<working_rows;r++){
      pos = r*cols;
      tmp_ptr = fullresolution_image->image + pos;
      for(c=0;c<working_cols;c++,pos++){

         /* val = (*(image->getpixel))(image, source_y+r, source_x+c); */
	 val = (double)fullresolution_data[pos];

         /*******************************************************************    
         * Always put this value in the green layer.
         *******************************************************************/
         if((int)val < min) val_uchar= 0;
         else if((int)val > max) val_uchar = 255;
         else val_uchar = (unsigned char)((((float)val-(float)min)*255.0)/((float)max-(float)min));

	 *tmp_ptr = val_uchar;
	 tmp_ptr++;
      }  
   }
      
   display_image(fullresolution_image, rows, cols, 1);

   if((overlay_filename[0] != '\0') && (show_overlay == 1))
      overlay_groundtruth(fullresolution_image, current_aggfactor*source_y,
	 current_aggfactor*source_x, (float)current_aggfactor, isvalid, overlay_data);

   if(segmentation_toggle)
      overlay_features(fullresolution_image, current_aggfactor*source_y,
       current_aggfactor*source_x, (float)current_aggfactor);

   if(display_regions_toggle)
      overlay_regions(fullresolution_image, current_aggfactor*source_y,
       current_aggfactor*source_x, (float)current_aggfactor,
                        num_detections, detections);

   /****************************************************************************
   * Draw a scale on the full resolution image.
   ****************************************************************************/
   if(samplerate_in_microns != 0.0){
      hold_scale_width = 36.0;    /* In inches. */
      do{

         hold_scale_width *= 0.95;
         scale_width = hold_scale_width;

         if(scale_width >= 2.0) scale_width = floor(scale_width * 4.0) / 4.0;
         else scale_width = floor(scale_width * 10.0) / 10.0;

         scale_width_pixels = floor(scale_width * 2.54 / (1.0e-4 * (samplerate_in_microns * current_aggfactor)));

      }while(scale_width_pixels > (working_cols - 20));

      XSetForeground(theDisplay, fullresolution_image->theGC, colors[GREENINDEX].pixel);  /* Green */

      XDrawLine(theDisplay, fullresolution_image->thePixmap, fullresolution_image->theGC,
         10, working_rows-10, 10+scale_width_pixels, working_rows-10);

      XDrawLine(theDisplay, fullresolution_image->thePixmap, fullresolution_image->theGC,
         10, working_rows-10, 10, working_rows-20);
      XDrawLine(theDisplay, fullresolution_image->thePixmap, fullresolution_image->theGC,
         10+scale_width_pixels, working_rows-10, 10+scale_width_pixels, working_rows-20);
   
      sprintf(thestring, "%.2f in", scale_width);
      XDrawString(theDisplay, fullresolution_image->thePixmap, fullresolution_image->theGC,
         10 + scale_width_pixels/2 - XTextWidth(info_helvb14, thestring, strlen(thestring))/2,
         working_rows - 10 - info_helvb14->ascent / 4, thestring, strlen(thestring));

   }

   repaint(*fullresolution_image, 0, 0, fullresolution_image->cols, fullresolution_image->rows, 0, 0);

   last_min = min;
   last_max = max;
}

/*******************************************************************************
* Procedure: plot_function
* Purpose: Given a vector and a window with a size and margins, this procedure
* plots the function.
* Name: Michael Heath, University of South Florida
* Date: 11/7/97
*******************************************************************************/
void plot_histogram(unsigned long int *hist, int n, int lm, int rm, int tm, int bm,
   int *min_val, int *max_val, int *bin_width, IMAGEXWIN *histogram_image)
{
   unsigned long int *hist_rescaled = NULL;
   int histogram_width=0;
   int rows, cols, r, c, cc, min_hist_pos, max_hist_pos;
   int bar_height, max_bar_height;
   int line_left, line_right, line_bottom, line_top;
   char thestring[100];

   unsigned long int background_color = cell[196];
   unsigned long int graph_color = colors[BLUEINDEX].pixel;
   unsigned long int number_color = cell[0];

   rows = histogram_image->rows;
   cols = histogram_image->cols;
   
   /****************************************************************************
   * Find the minimum and maximum values in the image.
   ****************************************************************************/
   min_hist_pos = 0;
   while((min_hist_pos < (n-1)) && (hist[min_hist_pos] == 0)) min_hist_pos++;
   *min_val = min_hist_pos;

   max_hist_pos = n-1;
   while((max_hist_pos > 0) && (hist[max_hist_pos] == 0)) max_hist_pos--;
   *max_val = max_hist_pos;

   /****************************************************************************
   * Allocate memory that we will use to rescale the histogram into for plotting
   * on the screen.
   ****************************************************************************/

   if(max_hist_pos == min_hist_pos) max_hist_pos = min_hist_pos + 1;
   *max_val = max_hist_pos;

   *bin_width = (int)ceil(((float)(max_hist_pos-min_hist_pos) / (float)(cols - (lm + rm)))); 

   histogram_width = (int)ceil((float)(max_hist_pos - min_hist_pos) / (float)(*bin_width));

   if((hist_rescaled = (unsigned long int *) calloc(histogram_width, sizeof(unsigned long int))) == NULL){
      fprintf(stderr, "Malloc error! here\n");
      return;
   }

   /****************************************************************************
   * Re-bin the histogram to fit in the range we have on the screen.
   ****************************************************************************/
   r = 0;
   cc = min_hist_pos;
   while(r < (histogram_width)){
      c = 0;
      while((cc < max_hist_pos) && (c < (*bin_width))){
         hist_rescaled[r] += hist[cc];
         c++;
         cc++;
      }
      r++;
   }

   /****************************************************************************
   * If our histogram fills less than 1/2 out the allowable screen area then
   * we set the histogram_stretching variable to a scalar greater than one
   * to stretch out the histogram on the screen. This does not change the
   * numbering.
   ****************************************************************************/
   if(histogram_width < (cols - lm - tm))
      histogram_stretching = (int)floor((double) (cols - lm - tm) / (double) histogram_width);
   else histogram_stretching = 1;

   /****************************************************************************
   * Find the position of the maximum values of the histogram.
   * Changed to have r=1 on 1/13/98.
   ****************************************************************************/
   for(r=1,max_hist_pos=1;r<(histogram_width);r++){
      if(hist_rescaled[max_hist_pos] < hist_rescaled[r]) max_hist_pos = r;
   }

   /*************************************************************************
   * Create a pixmap if we don't already have one.
   *************************************************************************/
   if(histogram_image->have_a_pixmap == 0){
      histogram_image->thePixmap = XCreatePixmap(theDisplay, histogram_image->theDrawWindow, cols, rows, 8);
      histogram_image->have_a_pixmap = 1;
   }

   /*************************************************************************
   * Set the histogram window to gray.
   *************************************************************************/
   XSetForeground(theDisplay, histogram_image->theGC, background_color);  /* Black */
   XFillRectangle(theDisplay, histogram_image->thePixmap, histogram_image->theGC, 0, 0, cols, rows);

   /****************************************************************************
   * Draw a box to contain the plot that we are drawing.
   ****************************************************************************/
   line_top = tm - 1;
   line_bottom = rows - bm + 1;
   line_left = lm - 1;
   line_right = lm + histogram_width * histogram_stretching + 1;
   if(line_top < 0) line_top = 0;
   if(line_bottom >= rows) line_bottom = rows - 1;
   if(line_left < 0) line_left = 0;
   if(line_right >= cols) line_right = cols - 1;

   XSetForeground(theDisplay, histogram_image->theGC, number_color);
   XDrawRectangle(theDisplay, histogram_image->thePixmap, histogram_image->theGC,
      line_left, line_top, line_right-line_left, line_bottom-line_top);
   /* XDrawLine(theDisplay, histogram_image->thePixmap, histogram_image->theGC, x1, y1, x2, y2); */

   /****************************************************************************
   * Draw the minimum and maximum values on the x-axis in density.
   ****************************************************************************/
   XSetForeground(theDisplay, histogram_image->theGC, number_color);

   if(mapping == DENSITY_MAPPING) sprintf(thestring, "%s", "Density");
   else sprintf(thestring, "%s", "Intensity");

   XDrawString(theDisplay, histogram_image->thePixmap, histogram_image->theGC,
      (3*line_right)/4 - XTextWidth(info_helvb14, thestring, strlen(thestring))/2,
      rows - info_helvb14->ascent / 4, thestring, strlen(thestring));

   if(mapping == DENSITY_MAPPING)
      sprintf(thestring, "%.3lf", (od_offset-(double)(*min_val))/od_scale + min_density);
   else sprintf(thestring, "%d", (*min_val));

   XDrawString(theDisplay, histogram_image->thePixmap, histogram_image->theGC,
      line_left, rows - info_helvb14->ascent / 4, thestring, strlen(thestring));

   if(mapping == DENSITY_MAPPING)
      sprintf(thestring, "%.3lf", (od_offset-(double)(*max_val))/od_scale + min_density);
   else sprintf(thestring, "%d", (*max_val));

   XDrawString(theDisplay, histogram_image->thePixmap, histogram_image->theGC,
      line_right - XTextWidth(info_helvb14, thestring, strlen(thestring)),
      rows - info_helvb14->ascent/4, thestring, strlen(thestring));

   /*************************************************************************
   * Now, plot the histogram.
   *************************************************************************/
   max_bar_height = rows - (tm + bm);
   XSetForeground(theDisplay, histogram_image->theGC, graph_color);
   for(r=0;r<histogram_width;r++){
      bar_height = (int)(max_bar_height * (float)hist_rescaled[r] / (float)hist_rescaled[max_hist_pos]);
      if(bar_height > max_bar_height) bar_height = max_bar_height;

      XFillRectangle(theDisplay, histogram_image->thePixmap, histogram_image->theGC, lm+r*histogram_stretching,
         (rows - bm) - bar_height, histogram_stretching, bar_height);
   }
 
   free(hist_rescaled);
}

/*******************************************************************************
* Procedure: rescale_overview_image
* Purpose: To scale the overview image and display a box on it showing the
* cropping.
* Name: Michael Heath, University of South Florida
* Date: 11/13/97
*******************************************************************************/
void rescale_overview_image(IMAGEXWIN *overview_image, unsigned short int *im,
   int num_detections, REGION *detections, int *isvalid, OVERLAY_DATA overlay_data)
{
   int rows, cols, r, c;
   UCHAR *tmp_ptr=NULL;
   double val, scalenum;

   tmp_ptr = overview_image->image;
   rows = overview_image->rows;
   cols = overview_image->cols;

   scalenum = 256.0 / (double)(scmax-scmin+1);

   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++,im++){

         val = (double)(*im);

         if(val < scmin) val = 0.0;
         else if(val > scmax) val = 255.0;
         else val = floor( (val-scmin) * scalenum);
    
         (*tmp_ptr) = (unsigned char)val;
         tmp_ptr++;
      }  
   }

   display_image(overview_image, overview_image->rows, overview_image->cols, 1);

   /****************************************************************************
   * Redraw the overlay on the image.
   ****************************************************************************/
   if((overlay_filename[0] != '\0') && (show_overlay == 1))
      overlay_groundtruth(overview_image, 0, 0, overview_subfactor,
         isvalid, overlay_data);

   if(segmentation_toggle)
      overlay_features(overview_image, 0, 0, overview_subfactor);

   if(display_regions_toggle)
      overlay_regions(overview_image, 0, 0, overview_subfactor,
                        num_detections, detections);

   /****************************************************************************
   * Draw the cropping box on the image.
   ****************************************************************************/
   XSetForeground(theDisplay, overview_image->theGC, colors[BLUEINDEX].pixel);  /* Blue */

   XDrawRectangle(theDisplay, overview_image->thePixmap, overview_image->theGC,
      box_x, box_y, box_width, box_height);
   repaint(*overview_image, 0, 0, overview_image->cols, overview_image->rows, 0, 0);
}

/*******************************************************************************
* Procedure: setup_overview_image
* Purpose: This procedure reads in three bands of an image from a file and
* scales them to fit on the screen. It then puts the image up on the screen
* in a window.
* Name: Michael Heath, University of South Florida
* Date: 7/2/97
*******************************************************************************/
void setup_overview_image(CACHEIM *image, IMAGEXWIN *overview_image, int maxrows, int maxcols,
   unsigned short int **overview_data, char *window_title, int col_position)
{
   int rows, cols, bytesperpixel, r, c;
   int screenrows=0, screencols=0;
   float scalefactor=0;
   double val;
   unsigned char *tmp_ptr=NULL;
   int count=0;
   double sum=0.0, sum_sq=0.0, scalenum=0.0;
   double stdev=0.0;
   double min=0, max=0;
   unsigned short int *tmp_ushort_im=NULL;
   int pos;
   ULONG *hist=NULL;

   rows = image->rows;
   cols = image->cols;

   bytesperpixel = dt_size[image->datatype];

   /****************************************************************************
   * Determine the scale factor to use and the size to display the image on
   * the screen.
   ****************************************************************************/
   if((rows<=maxrows) && (cols<=maxcols)){
      scalefactor = 1.0;
      screenrows = rows;
      screencols = cols;
   }
   else{
      if(((float)rows/(float)cols) > ((float)maxrows/(float)maxcols)){  /* Rows are limiting. */
         if(rows < maxrows){
            scalefactor = 1.0;
            screenrows = rows;
            screencols = cols;
         }   
         else{
            scalefactor = (float)maxrows / (float)rows;
            screenrows = maxrows;
            screencols = scalefactor * cols;
         }   
      }
      else{                                                             /* Cols are limiting. */
         if(cols < maxcols){
            scalefactor = 1.0;
            screenrows = rows;
            screencols = cols;
         }   
         else{
            scalefactor = (float)maxcols / (float)cols;
            screencols = maxcols;
            screenrows = scalefactor * rows;
         }   
      }
   }     
      
   if(verbosemode) printf("The screenimage will be %d rows x %d cols.\n", screenrows, screencols);
 
   /****************************************************************************
   * Allocate memory for the resized image.
   ****************************************************************************/
   overview_image->rows = screenrows;
   overview_image->cols = screencols;
   if((overview_image->image = (unsigned char *) calloc(screenrows*screencols,
      sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Malloc error.\n");
      exit(1);
   }
   if((overview_image->image_dsp = (unsigned char *) calloc(screenrows*screencols,
      sizeof(unsigned char))) == NULL){
      fprintf(stderr, "Malloc error.\n");
      exit(1);
   }
   overview_image->maskimage = NULL;
   overview_image->theImage = NULL;

   *overview_data = tmp_ushort_im = (USHORT *) calloc(screenrows*screencols, sizeof(USHORT));

   /*************************************************************************    
   * Calculate statistics for the image.
   *************************************************************************/
   if(verbosemode) printf("Calculating the statistics for scaling the positioning image.\n");

   /****************************************************************************
   * Allocate a histogram array.
   ****************************************************************************/
   if((hist = (ULONG *) calloc(65536, sizeof(ULONG))) == NULL){
      fprintf(stderr, "Calloc error!\n");
      exit(1);
   }

   for(r=0,pos=0;r<screenrows;r++){
      if(verbosemode){
	 printf("\r   %d percent complete", (int)floor((100.0*r)/(double)(screenrows)));
	 fflush(stdout);
      }
      for(c=0;c<screencols;c++,pos++){
   
         val = (*(image->getpixel))(image, (int)(r/scalefactor), (int)(c/scalefactor));

	 tmp_ushort_im[pos] = (USHORT)val;
         hist[(USHORT)val]++; 
      }
   }
   if(verbosemode) printf("\r\n");

   r=0; c=0;
   while((r < 65535) && (c < (screenrows*screencols/50))){
      c += hist[r];
      r++;
   }
   min = r;

   r=65535; c=0;
   while((r > 0) && (c < (screenrows*screencols/400))){
      c += hist[r];
      r--;
   }
   max = r;

   /****************************************************************************
   * Scale the image for display.
   ****************************************************************************/
   tmp_ptr = overview_image->image;

   scalenum = 256.0 / (double)(max-min+1);

   for(r=0,pos=0;r<screenrows;r++){
      for(c=0;c<screencols;c++,pos++){

         /*val = (*(image->getpixel))(image, (int)(r/scalefactor), (int)(c/scalefactor)); */

         val = (double)(tmp_ushort_im[pos]);

         if(val < min) val = 0.0;
         else if(val > max) val = 255.0;
         else val = floor( (val-min) * scalenum);
    
         (*tmp_ptr) = (unsigned char)val;
         tmp_ptr++;
      }  
   }
    
   /****************************************************************************
   * Set up the image on the screen.
   ****************************************************************************/
   overview_image->theBorderWindow = openWindow(col_position, 0, screencols, screenrows,
                                     NORMAL_WINDOW, window_title, theRootWindow,
                                     &(overview_image->theBorderGC));
    
   initEvents(overview_image->theBorderWindow, IN_PALETTE);
 
   overview_image->theDrawWindow = openWindow(0, 0, screencols, screenrows,
                                NORMAL_WINDOW, image->filename,
                                overview_image->theBorderWindow,
                                &(overview_image->theGC));
 
   initEvents(overview_image->theDrawWindow, NOT_IN_PALETTE );
 
   display_image(overview_image, screenrows, screencols, 1);
}

/*******************************************************************************
* Procedure: swap
* Purpose: To swap two bytes.
* Name: Michael Heath, University of South Florida
* Date: 10/24/97
*******************************************************************************/
void swap(unsigned char two_bytes[2])
{
  unsigned char temp_char;

  temp_char = two_bytes[0];
  two_bytes[0] = two_bytes[1];
  two_bytes[1] = temp_char;
}

/*******************************************************************************
* Procedure: readDBAImage_header
* Purpose: To read in the header information from a DBA image.
* Name: Michael Heath, University of South Florida
* Date: 10/24/97
*******************************************************************************/
int readDBAImage_header(char *file_name, int *rows, int *cols, int *bytesperpixel,
   int *header_bytes, int *swapbytes, int *pixelskip, int *lineskip)
{

  FILE *file_ptr, *fopen();
  unsigned short temp_short;
  int orig_width;
  int orig_height;
  unsigned short pixel_skip;
  unsigned short line_skip;

  *bytesperpixel = 2;  /* We assume all dba images are two bytes per pixel. */
  *swapbytes = TRUE;

  file_ptr = fopen(file_name, "r");
  if(file_ptr == NULL)
  {
    printf("\nERROR: file %s not found\n", file_name);
    return(0);
  }

  fseek(file_ptr, 10, 0);

  fread(&temp_short, 2, 1, file_ptr);
  swap((unsigned char *)&temp_short);
  *cols = orig_width = (int)temp_short;

  fread(&temp_short, 2, 1, file_ptr);
  swap((unsigned char *)&temp_short);
  *pixelskip = pixel_skip = (int)temp_short;

  fread(&temp_short, 2, 1, file_ptr);
  swap((unsigned char *)&temp_short);
  *rows = orig_height = (int)temp_short ;

  fread(&temp_short, 2, 1, file_ptr);
  swap((unsigned char *)&temp_short);
  *lineskip = line_skip = (int)temp_short;

  fseek(file_ptr, 1024, 0);

  *header_bytes = ftell(file_ptr);

  fclose(file_ptr);

   return(1);
}

/*******************************************************************************
* Procedure: readIFSImage_header
* Purpose: To read in the header information from a IFS image.
* Name: Michael Heath, University of South Florida
* Date: 1/22/98
*******************************************************************************/
int readIFSImage_header(char *file_name, int *rows, int *cols, int *bytesperpixel,
   int *header_bytes, int *swapbytes, int *pixelskip, int *lineskip)
{
   FILE *fp=NULL;
   unsigned short temp_short;

   /****************************************************************************
   * Open the file for reading.
   ****************************************************************************/
   if((fp = fopen(file_name, "rb")) == NULL){
      fprintf(stderr, "Error opening the file %s in readIFSImage_header().\n", file_name);
      exit(1);
   }
  
   /****************************************************************************
   * The bytes do not need to be swapped in an IFS image.
   ****************************************************************************/
   *swapbytes = FALSE;

   /****************************************************************************
   * Determine the number of bytes per pixel by the filename. If it has "small"
   * in the filename then the image is 8 bits/pixel. Otherwise it is 16 bits
   * per pixel.
   ****************************************************************************/
   *bytesperpixel = 2;  /* We assume all IFS images are two bytes per pixel. */
   if(strstr(file_name, "small") != NULL) *bytesperpixel = 1;

   /****************************************************************************
   * Get the number of rows and columns in the image.
   ****************************************************************************/
   fseek(fp, 258, 0);
   fread(&temp_short, 2, 1, fp);
   *cols = (int)temp_short;

   fseek(fp, 322, 0);
   fread(&temp_short, 2, 1, fp);
   *rows = (int)temp_short;

   fclose(fp);

   /****************************************************************************
   * There are 512 header bytes in IFS images.
   ****************************************************************************/
   *header_bytes = 512;

   *pixelskip = 0;
   *lineskip = 0;

   return(1);
}

/*******************************************************************************
* Function: overlay_groundtruth
* Purpose: This function will draw the boundary of valid (selected) ground
* truth regions on the image.
* Name: Michael Heath, University of South Florida
* Date: 2/3/2000
*******************************************************************************/
int overlay_groundtruth(IMAGEXWIN *image, int rows_offset, int cols_offset,
    float subfactor, int *isvalid, OVERLAY_DATA overlay_data)
{
   FILE *fp=NULL;
   ABNORMALITY abnormality;
   int start_row=0, start_col=0;
   int a=0, d, p;
   int x_array[] = {0, 1, 1, 1, 0, -1, -1, -1};
   int y_array[] = {-1, -1, 0, 1, 1, 1, 0, -1};
   int x_pixel=0, y_pixel=0;
   int direction_to_go=0;
   char direction_to_go_string[4];
   char ch;
   int chaincolor;

   if(MIASgt_filename != NULL){

      FILE *fpgt=NULL;

      struct MIASGT{
         char image_filename[200], bgtissue[10], abnormality_class[10], severity[10];
         int rows, cols;
         int row, col, radius;
         int valid;
         int matched;
      } truth_region[40];
      char gtline[100];
      int num_truth_regions;
      int pixel_radius;

      /****************************************************************************
      * Read in the ground truth file.
      ****************************************************************************/
      XSetForeground(theDisplay, image->theGC, colors[REDINDEX].pixel); /* Red */
      XSetLineAttributes(theDisplay, image->theGC, 2, LineSolid, CapRound, JoinMiter);

      chaincolor = CHAINCOLOR1 - 1;
      num_truth_regions = 0;
      if((fpgt = fopen(MIASgt_filename, "r")) != NULL){
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

               chaincolor++;

               if(chaincolor > CHAINCOLOR9) chaincolor = CHAINCOLOR1;
               XSetForeground(theDisplay, image->theGC, colors[chaincolor].pixel);

               x_pixel = (int)((float)(truth_region[num_truth_regions].col * 50.0/samplerate_in_microns - cols_offset) / subfactor);
               y_pixel = (int)((float)(truth_region[num_truth_regions].row * 50.0/samplerate_in_microns - rows_offset) / subfactor);

               pixel_radius = (truth_region[num_truth_regions].radius * 50.0/samplerate_in_microns) / subfactor;

               XDrawArc(theDisplay, image->thePixmap, image->theGC,
                  x_pixel - pixel_radius, y_pixel - pixel_radius,
                  2*pixel_radius, 2*pixel_radius, 0, 360*64);

               num_truth_regions++;
            }
         }
         fclose(fpgt);
      }
   }

   /* Uncomment this to make the images used in the book chapter.
      XSetLineAttributes(theDisplay, image->theGC, 2, LineSolid, CapRound, JoinMiter);
      XSetForeground(theDisplay, image->theGC, cell[255]);
   */

   else{

      XSetForeground(theDisplay, image->theGC, colors[REDINDEX].pixel); /* Red */
      XSetLineAttributes(theDisplay, image->theGC, 2, LineSolid, CapRound, JoinMiter);

      chaincolor = CHAINCOLOR1;

      for(d=0;d<overlay_data.total_abnormalities;d++){
         abnormality = overlay_data.abnormalities[d];

         if(isvalid[d] == 1){
            if(chaincolor > CHAINCOLOR9) chaincolor = CHAINCOLOR1;
            XSetForeground(theDisplay, image->theGC, colors[chaincolor].pixel);
            XSetLineAttributes(theDisplay, image->theGC, 2, LineSolid, CapRound, JoinMiter);
            XFlush(theDisplay);

            /*******************************************************************
            * Draw the boundary on the image.
            *******************************************************************/
            for(p=0;p<abnormality.boundary.length;p++){

               x_pixel = (int)((float)(abnormality.boundary.chain[p].c/base_subfactor - cols_offset) / subfactor);
               y_pixel = (int)((float)(abnormality.boundary.chain[p].r/base_subfactor - rows_offset) / subfactor);

               if((x_pixel >= 0) && (x_pixel < image->cols) && (y_pixel >= 0) && (y_pixel < image->rows)){

                  /*
                  XDrawPoint(theDisplay, image->thePixmap, image->theGC, x_pixel, y_pixel);
                  */

                  /* Uncomment this to make the images used in the book chapter. */
                  XFillRectangle(theDisplay, image->thePixmap, image->theGC, x_pixel-1, y_pixel-1, 3, 3);

               }
            }
            /*******************************************************************
            * Draw the first core (if it exists) on the image.
            *******************************************************************/
            if(abnormality.num_cores != 0){
               for(p=0;p<abnormality.core[0].length;p++){

                  x_pixel = (int)((float)(abnormality.core[0].chain[p].c/base_subfactor - cols_offset) / subfactor);
                  y_pixel = (int)((float)(abnormality.core[0].chain[p].r/base_subfactor - rows_offset) / subfactor);

                  if((x_pixel >= 0) && (x_pixel < image->cols) && (y_pixel >= 0) && (y_pixel < image->rows)){

                     /*
                     XDrawPoint(theDisplay, image->thePixmap, image->theGC, x_pixel, y_pixel);
                     */

                     /* Uncomment this to make the images used in the book chapter. */
                     XFillRectangle(theDisplay, image->thePixmap, image->theGC, x_pixel-1, y_pixel-1, 3, 3);
                  }
               }
            }
            chaincolor++;
         }
      }   
   }
 
   /* repaint(*image, 0, 0, image->cols, image->rows, 0, 0); */

   return(1);
}

int overlay_nijmegen(IMAGEXWIN *image, int rows_offset, int cols_offset, float subfactor)
{
   FILE *fp=NULL;
   static char n_filename[100] = {'\0'};
   int total_abnormalities=0;
   char aline[100];
   int start_row=0, start_col=0;
   int a=0, radius;
   int x_pixel=0, y_pixel=0;
   int chaincolor, pixel_radius;

   /****************************************************************************
   * Open the overlay file.
   ****************************************************************************/
   if(n_filename[0] == '\0'){
      printf("Enter the name of the overlay file: ");
      scanf("%s", n_filename);
   }
   if((fp = fopen(n_filename, "r")) == NULL){
      fprintf(stderr, "Error opening the overlay file %s.\n", n_filename);
      return(0);
   }
 
   fscanf(fp, "%d", &total_abnormalities);
 
   XSetForeground(theDisplay, image->theGC, colors[REDINDEX].pixel); /* Red */

   chaincolor = CHAINCOLOR1 - 1;

   for(a=0;a<total_abnormalities;a++){

      chaincolor++;

      if(chaincolor > CHAINCOLOR9) chaincolor = CHAINCOLOR1;
      XSetForeground(theDisplay, image->theGC, colors[chaincolor].pixel);

      fscanf(fp, "%d %d %d", &start_col, &start_row, &radius);
 
      x_pixel = (int)((float)(start_col - cols_offset) / subfactor);
      y_pixel = (int)((float)(start_row - rows_offset) / subfactor);

      pixel_radius = radius / subfactor;

      XDrawArc(theDisplay, image->thePixmap, image->theGC,
         x_pixel - pixel_radius, y_pixel - pixel_radius,
         2*pixel_radius, 2*pixel_radius, 0, 360*64);

   }   
 
   fclose(fp);
 
   repaint(*image, 0, 0, image->cols, image->rows, 0, 0);

   return(1);
}

/*******************************************************************************
* Function: overlay_regions
* Purpose: This function draws regions stored in a file onto the image.
* is displayed.
* Name: Michael Heath, University of South Florida
* Date: 9/22/98
*******************************************************************************/
int overlay_regions(IMAGEXWIN *image, int rows_offset, int cols_offset,
   float subfactor, int num_detections, REGION *detections)
{
   int numregions = 0, numpoints = 0, p, r;
   int *x_coord=NULL, *y_coord=NULL;
   float x0, y0;
   int chaincolor, rpos, cpos;
   float pixel_diameter;
   int distance, x_pixel_start, y_pixel_start, x_pixel_end, y_pixel_end;
   int start_col, start_row;

   chaincolor = CHAINCOLOR1 - 1;

   for(r=0;r<num_detections;r++){

      chaincolor = CHAINCOLOR1 + r%10;

      /*************************************************************************
      * Overlay the border of the region.
      *************************************************************************/
      
      numpoints = detections[r].outline.length;

      if(numpoints == 0){

         cpos = detections[r].centroid_c;
         rpos = detections[r].centroid_r;

         pixel_diameter = (5000.0 / samplerate_in_microns) / subfactor;

         if(samplerate_in_microns == 0.0) pixel_diameter = 10;

         XSetForeground(theDisplay, image->theGC, colors[chaincolor].pixel);
         XDrawArc(theDisplay, image->thePixmap, image->theGC,
            (int)((cpos - cols_offset) / subfactor - pixel_diameter/2),
            (int)((rpos - rows_offset) / subfactor - pixel_diameter/2),
            pixel_diameter, pixel_diameter, 0, 360*64);

      }
      else{
         x_coord = (int *) calloc(numpoints, sizeof(int));
         y_coord = (int *) calloc(numpoints, sizeof(int));

         for(p=0;p<numpoints;p++){
            x_coord[p] = detections[r].outline.chain[p].c;
            y_coord[p] = detections[r].outline.chain[p].r;
         }

         overlay_polycurve(image, rows_offset, cols_offset,
            subfactor, x_coord, y_coord, numpoints, chaincolor);

         free(x_coord);
         free(y_coord);
      }

      if(bothdet == 1){

         start_col = seg_startx;
         start_row = seg_starty;

         x_pixel_start = (int)((float)(start_col - cols_offset) / subfactor);
         y_pixel_start = (int)((float)(start_row - rows_offset) / subfactor);

         cpos = other_detections[r].centroid_c;
         rpos = other_detections[r].centroid_r;

         distance = sqrt((cpos - other_seg_startx) * (cpos - other_seg_startx) +
                         (rpos - other_seg_starty) * (rpos - other_seg_starty));

         distance *= (float)ROWS / (float)seg_rows;

         distance /= subfactor;

         x_pixel_end = x_pixel_start - distance;
         y_pixel_end = y_pixel_start - distance;

         XSetForeground(theDisplay, image->theGC, colors[chaincolor].pixel);
         XDrawArc(theDisplay, image->thePixmap, image->theGC,
            x_pixel_end, y_pixel_end,
            2*distance, 2*distance, 0, 360*64);


         repaint(*image, 0, 0, image->cols, image->rows, 0, 0);
      }
   }
}

/*******************************************************************************
* Function: overlay_features
* Purpose: This function draws the breast coordinate features on an image that
* is displayed.
* Name: Michael Heath, University of South Florida
* Date: 8/10/98
*******************************************************************************/
int overlay_features(IMAGEXWIN *image, int rows_offset, int cols_offset,
   float subfactor)
{
   int start_row=0, start_col=0, end_col, end_row;
   int x_pixel_start=0, y_pixel_start=0, x_pixel_end=0, y_pixel_end=0;

   overlay_polycurve_float(image, rows_offset, cols_offset,
        subfactor, seg_xcoord, seg_ycoord, seg_numpoints, ORANGEINDEX);

   /****************************************************************************
   * Set the color to draw with.
   ****************************************************************************/
   XSetForeground(theDisplay, image->theGC, YELLOWINDEX);

   XSetLineAttributes(theDisplay, image->theGC, 2, LineSolid, CapRound, JoinMiter);

   start_col = seg_startx;
   start_row = seg_starty;
   end_col = seg_endx;
   end_row = seg_endy;

   x_pixel_start = (int)((float)(start_col - cols_offset) / subfactor);
   y_pixel_start = (int)((float)(start_row - rows_offset) / subfactor);

   x_pixel_end = (int)((float)(end_col - cols_offset) / subfactor);
   y_pixel_end = (int)((float)(end_row - rows_offset) / subfactor);

   /* Draw the major axis on the image. */
   XDrawLine(theDisplay, image->thePixmap, image->theGC, x_pixel_start, y_pixel_start,
      x_pixel_end, y_pixel_end);

}

/*******************************************************************************
* Function: overlay_polycurve
* Purpose: This function draws a curve of line segments in the window that is
* passed to the function. The idea is that this curve could be the breast
* boundary, the chest wall boundary the nipple location or whatever else the
* user wants drawn on the image.
* Name: Michael Heath, University of South Florida
* Date: 5/23/98
*******************************************************************************/
int overlay_polycurve(IMAGEXWIN *image, int rows_offset, int cols_offset,
   float subfactor, int *x_coord, int *y_coord, int numpoints, int color)
{
   int a=0;
   int start_row=0, start_col=0, end_col, end_row;
   int x_pixel_start=0, y_pixel_start=0, x_pixel_end=0, y_pixel_end=0;

   /****************************************************************************
   * If there is nothing to draw, just return.
   ****************************************************************************/
   if(numpoints < 1) return(1);

   /****************************************************************************
   * Set the color to draw with.
   ****************************************************************************/
   XSetForeground(theDisplay, image->theGC, colors[color].pixel); /* Violet */

   XSetLineAttributes(theDisplay, image->theGC, 2, LineSolid, CapRound, JoinMiter);

   end_col = x_coord[0];
   end_row = y_coord[0];

   for(a=1;a<numpoints;a++){

      start_col = end_col;
      start_row = end_row;

      end_col = x_coord[a];
      end_row = y_coord[a];
 
      x_pixel_start = (int)((float)(start_col - cols_offset) / subfactor);
      y_pixel_start = (int)((float)(start_row - rows_offset) / subfactor);

      x_pixel_end = (int)((float)(end_col - cols_offset) / subfactor);
      y_pixel_end = (int)((float)(end_row - rows_offset) / subfactor);

      XDrawLine(theDisplay, image->thePixmap, image->theGC, x_pixel_start, y_pixel_start,
         x_pixel_end, y_pixel_end);

   }   
 
   start_col = end_col;
   start_row = end_row;

   end_col = x_coord[0];
   end_row = y_coord[0];
 
   x_pixel_start = (int)((float)(start_col - cols_offset) / subfactor);
   y_pixel_start = (int)((float)(start_row - rows_offset) / subfactor);

   x_pixel_end = (int)((float)(end_col - cols_offset) / subfactor);
   y_pixel_end = (int)((float)(end_row - rows_offset) / subfactor);

   XDrawLine(theDisplay, image->thePixmap, image->theGC, x_pixel_start, y_pixel_start,
      x_pixel_end, y_pixel_end);

   /* repaint(*image, 0, 0, image->cols, image->rows, 0, 0); */

   return(1);
}

/*******************************************************************************
* Function: overlay_polycurve_float
* Purpose: This function draws a curve of line segments in the window that is
* passed to the function. The idea is that this curve could be the breast
* boundary, the chest wall boundary the nipple location or whatever else the
* user wants drawn on the image.
* Name: Michael Heath, University of South Florida
* Date: 5/23/98
*******************************************************************************/
int overlay_polycurve_float(IMAGEXWIN *image, int rows_offset, int cols_offset,
   float subfactor, float *x_coord, float *y_coord, int numpoints, int color)
{
   int a=0;
   int start_row=0, start_col=0, end_col, end_row;
   int x_pixel_start=0, y_pixel_start=0, x_pixel_end=0, y_pixel_end=0;

   /****************************************************************************
   * If there is nothing to draw, just return.
   ****************************************************************************/
   if(numpoints < 1) return(1);

   /****************************************************************************
   * Set the color to draw with.
   ****************************************************************************/
   XSetForeground(theDisplay, image->theGC, colors[color].pixel); /* Violet */

   end_col = x_coord[0];
   end_row = y_coord[0];

   for(a=1;a<numpoints;a++){

      start_col = end_col;
      start_row = end_row;

      end_col = x_coord[a];
      end_row = y_coord[a];
 
      x_pixel_start = (int)((float)(start_col - cols_offset) / subfactor);
      y_pixel_start = (int)((float)(start_row - rows_offset) / subfactor);

      x_pixel_end = (int)((float)(end_col - cols_offset) / subfactor);
      y_pixel_end = (int)((float)(end_row - rows_offset) / subfactor);

      XDrawLine(theDisplay, image->thePixmap, image->theGC, x_pixel_start, y_pixel_start,
         x_pixel_end, y_pixel_end);

   }   
 
   start_col = end_col;
   start_row = end_row;

   end_col = x_coord[0];
   end_row = y_coord[0];
 
   x_pixel_start = (int)((float)(start_col - cols_offset) / subfactor);
   y_pixel_start = (int)((float)(start_row - rows_offset) / subfactor);

   x_pixel_end = (int)((float)(end_col - cols_offset) / subfactor);
   y_pixel_end = (int)((float)(end_row - rows_offset) / subfactor);

   XDrawLine(theDisplay, image->thePixmap, image->theGC, x_pixel_start, y_pixel_start,
      x_pixel_end, y_pixel_end);

   /* repaint(*image, 0, 0, image->cols, image->rows, 0, 0); */

   return(1);
}

/*******************************************************************************
* Function: is_ISMD
* Purpose: This function checks the header of the .dcm image file for the
* words "Intelligent Systems". If this is found in the first 2000 bytes of the
* file, the image is assumed to come from the HOWTEK scanner at ISMD. This is
* a little crude (there are ways to do this by looking for certain tags), but
* it seems to work.
* Name: Michael Heath, University of South Florida
* Date: 3/3/99
*******************************************************************************/
int is_ISMD(char *filename)
{
   FILE *fp;
   char line[2000];
   int k;

   if((fp = fopen(filename, "rb")) == NULL){
      fprintf(stderr, "Error opening the file %s.\n", filename);
      exit(1);
   }

   fread(line, 1, 2000, fp);

   fclose(fp);

   for(k=0;k<2000;k++){
      if(line[k] == 0) line[k] = 32;
   }
   line[2000-1] = '\0';

   if(strstr(line, "Intelligent Systems") != NULL){
      return(1);
      /* printf("The image (%s) is from ISMD.\n", argv[1]); */
   }
   else{
      return(0);
      /* printf("The image (%s) is not from ISMD.\n", argv[1]); */
   }
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

void update_zoom_window(CACHEIM *image, IMAGEXWIN *zoom_image, int x_buttonpos, int y_buttonpos, int zoom_factor)
{
   int est_file_x, est_file_y, est_rows, est_cols, er, ec, err, ecc, e_val;
   unsigned char e_val_uchar;
   float left_fract, right_fract, top_fract, bottom_fract;
   int left_x, right_x, top_y, bottom_y;

   memset(zoom_image->image, 0, 512*512);

   est_file_x = current_aggfactor * (image_ul_x + x_buttonpos);
   est_file_y = current_aggfactor * (image_ul_y + y_buttonpos);

   if((image->rows * zoom_factor) < 512) est_rows = ceil((double)(image->rows) / (float)zoom_factor);
   else est_rows = ceil(512.0 / (float)zoom_factor);
   if((image->cols * zoom_factor) < 512) est_cols = ceil((double)(image->cols) / (float)zoom_factor);
   else est_cols = ceil(512.0 / (float)zoom_factor);

   est_file_x -= (est_cols/2);
   if((est_file_x + est_cols) >= image->cols) est_file_x = image->cols-1-est_cols;
   if(est_file_x < 0) est_file_x = 0;

    est_file_y -= (est_rows/2);
    if((est_file_y + est_rows) >= image->rows) est_file_y = image->rows-1-est_rows;
    if(est_file_y < 0) est_file_y = 0;

/*
    for(er=0;er<est_rows;er++){
       for(ec=0;ec<est_cols;ec++){
          e_val =  (int)(*(image->getpixel))(image, est_file_y+er, est_file_x+ec);

          if((int)e_val < last_min) e_val_uchar= 0;
          else if((int)e_val > last_max) e_val_uchar = 255;
          else e_val_uchar = (unsigned char)((((float)e_val-(float)last_min)*255.0)/((float)last_max-(float)last_min));

          for(err=0;err<zoom_factor;err++){
             for(ecc=0;ecc<zoom_factor;ecc++){
                zoom_image->image[(er*zoom_factor+err) * 512 + (ec*zoom_factor+ecc)] = e_val_uchar;
             }
          }
       }
    }
*/

   for(er=0;er<est_rows;er++){
      for(ec=0;ec<est_cols;ec++){
         e_val =  (int)(*(image->getpixel))(image, est_file_y+er, est_file_x+ec);

         if((int)e_val < last_min) e_val_uchar= 0;
         else if((int)e_val > last_max) e_val_uchar = 255;
         else e_val_uchar = (unsigned char)((((float)e_val-(float)last_min)*255.0)/((float)last_max-(float)last_min));

         zoom_image->image[(er*zoom_factor) * 512 + (ec*zoom_factor)] = e_val_uchar;
      }
   }

   for(er=0;er<est_rows;er++){
      for(ec=0;ec<est_cols;ec++){
         for(ecc=1;ecc<zoom_factor;ecc++){
            left_x = ec;
            right_x = ec+1;
            if(right_x == image->cols) right_x = image->cols - 1;
            left_fract = (float)ecc / (float)zoom_factor;
            right_fract = 1.0 - left_fract;
            zoom_image->image[(er*zoom_factor) * 512 + (ec*zoom_factor) + ecc] =
               left_fract * zoom_image->image[(er*zoom_factor) * 512 + (left_x*zoom_factor)] + 
               right_fract * zoom_image->image[(er*zoom_factor) * 512 + (right_x*zoom_factor)];
         }
      }
   }

   for(er=0;er<est_rows;er++){
      for(err=1;err<zoom_factor;err++){
         top_y = er;
         bottom_y = er+1;
         if(bottom_y == image->rows) bottom_y = image->rows-1;
         top_fract = (float)err / (float)zoom_factor;
         bottom_fract = 1.0 - top_fract;
         for(ecc=0;ecc<512;ecc++){
            zoom_image->image[(er*zoom_factor+err) * 512 + ecc] =
               top_fract * zoom_image->image[(top_y*zoom_factor) * 512 + ecc] + 
               bottom_fract * zoom_image->image[(bottom_y*zoom_factor) * 512 + ecc];
         }
      }
   }

   XRaiseWindow(theDisplay, zoom_image->theBorderWindow);
   display_image(zoom_image, 512, 512, 1);
}
