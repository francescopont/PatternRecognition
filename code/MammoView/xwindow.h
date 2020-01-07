/*******************************************************************************
* HEADER FILE: xwindow.h
* PURPOSE: This program contains the necessary functions to create windows in
* Xwindows and to display images in them. The images are displayed in palette
* color. This code has been run on 8-bit and 24-bit displays.
* AUTHOR: Mike Heath
* DATE: 4/10/95
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#define         POP_UP_WINDOW   1
#define         BORDER_WIDTH    0
#define         NORMAL_WINDOW   0
#define         IN_PALETTE      1
#define         NOT_IN_PALETTE  2
#define EV_MASK (ButtonPressMask   | \
		KeyPressMask       | \
		ExposureMask       | \
		StructureNotifyMask)

#define LEFT_BUTTON 1
#define MIDDLE_BUTTON 2
#define RIGHT_BUTTON 3

Display         *theDisplay;    /* -- Which display                  */
int             theScreen;      /* -- Which screen on the display    */
int             vis_depth;       /* -- Number of color planes         */
unsigned long   theBlackPixel;  /* -- System "Black" color           */
unsigned long   theWhitePixel;  /* -- System "White" color           */
Window          theRootWindow;  /* -- System-wide parent window      */
Visual          *theVisual;     /* -- System visual type             */
Colormap        theColormap;    /* -- System color map ID            */

int MaxScreenRows, MaxScreenCols;

unsigned long cell[256];
int colormap_size;
XVisualInfo xvis_info, *vis_list;
int         num_visuals;
int workstation_bits = 0;
short int unmapper[256];

/*******************************************************************************
* Define some colors to use for drawing in color.
*******************************************************************************/
#define REDINDEX 128
#define ORANGEINDEX 129
#define YELLOWINDEX 130
#define GREENINDEX 131
#define BLUEINDEX 132
#define INDIGOINDEX 133
#define VIOLETINDEX 134

#define CHAINCOLOR1 135
#define CHAINCOLOR2 136
#define CHAINCOLOR3 137
#define CHAINCOLOR4 138
#define CHAINCOLOR5 139
#define CHAINCOLOR6 140
#define CHAINCOLOR7 141
#define CHAINCOLOR8 142
#define CHAINCOLOR9 143

XColor *colors = NULL;

typedef struct{
   Window   theDrawWindow;
   Window   theBorderWindow;
   GC       theBorderGC;      /* -- Border window GC    */
   GC       theGC;            /* -- Draw the image  */
   unsigned char *image;
   unsigned char *image_dsp;
   unsigned char *maskimage;
   XImage *theImage;
   Pixmap  thePixmap;
   int rows, cols;
   int have_a_pixmap;
}IMAGEXWIN;

/****************************************************************************
* Function: pick_visual
* Purpose: To select a visual of appropriate "class" and "depth" from the
* list of visuals. Return 1 if successful, 0 if no matching visuals found.
* This was taken from viewpcx.
****************************************************************************/
static int pick_visual(int depth_wanted, int class_wanted)
{
    XVisualInfo *p_visinfo;
    int i, status = 0;

    for(i=0,p_visinfo=vis_list;i<num_visuals;i++, p_visinfo++){
       if((p_visinfo->class == class_wanted) && (p_visinfo->depth >= depth_wanted)){
          theVisual = p_visinfo->visual;
          vis_depth = p_visinfo->depth;
          xvis_info = *(p_visinfo);
          status = 1;
          break;
       }
    }
    return (status);
}

/****************************************************************************
* Function: print_visuals
* Purpose: List the available visuals that we can select from.
****************************************************************************/
void print_visuals()
{
   XVisualInfo *p_visinfo;
   int i;

   printf("Printing for %d Available Visuals\n", num_visuals);
   for(i=0,p_visinfo=vis_list;i<num_visuals;i++, p_visinfo++){
      printf("Class: %d\t Depth: %d\n", p_visinfo->class, p_visinfo->depth);
   }
}

/*******************************************************************************
* Establish a connection to the X-Server.
*******************************************************************************/
initX( )
{
   XVisualInfo vis_template;
   int bits_per_pixel = 8;

   /* Establish a connection to the X-Server */
   theDisplay = XOpenDisplay(NULL);

   /* Check if the connection was made */
   if(theDisplay == NULL){
      fprintf(stderr, "ERROR: Cannot establish a connection to the X Server ");
      fprintf(stderr, "%s\n", XDisplayName(NULL));
      exit(1);
   }

   theScreen = DefaultScreen( theDisplay);
   vis_template.screen = theScreen;

   /****************************************************************************
   * Get a list of the available visuals for this screen.
   ****************************************************************************/
   vis_list = XGetVisualInfo(theDisplay, VisualScreenMask, &vis_template,
                             &num_visuals);
   if(num_visuals == 0){
      fprintf(stderr, "No visuals found!\n");
      exit(0);
   }
   /* else print_visuals(); */

   /****************************************************************************
   * Search for a PseudoColor visual with depth >= image depth. You may want
   * other strategies here--for visuals such as "DirectColor" or for
   * StaticColor with 8 or more bit planes.
   ****************************************************************************/
   if(!pick_visual(bits_per_pixel, PseudoColor)){
      fprintf(stderr, "No appropriate visual...Exiting\n");
      exit(0);
   }

   theBlackPixel = BlackPixel(theDisplay, theScreen);
   theWhitePixel = WhitePixel(theDisplay, theScreen);
   theRootWindow = RootWindow(theDisplay, theScreen);

   MaxScreenRows = DisplayHeight(theDisplay, theScreen);
   MaxScreenCols = DisplayWidth(theDisplay, theScreen);

   /* printf("The screen size is %d x %d\n", MaxScreenCols, MaxScreenRows); */
}


initGrayScalecolormap()
{
   int i;
   if(workstation_bits == 8) initGrayScalecolormap_8bit();
   if(workstation_bits == 24) initGrayScalecolormap_24bit();

   for(i=0;i<135;i++){
      unmapper[colors[i].pixel] = i;
      /*
      printf("colors[%d] (%8d) (%8d) (%3d, %3d, %3d)\n", i, colors[i].pixel, (int)unmapper[colors[i].pixel],
	 colors[i].red/256, colors[i].green/256, colors[i].blue/256);
      */
   }
}

/*******************************************************************************
* Initialize the colormap.
*******************************************************************************/
initGrayScalecolormap_24bit()
{
   int i, j, theColormapSize, good;
   int warnings=0;
   int get_a_color(XColor *, int , int , int , int );
   int cindex;

   /****************************************************************************
   * Create the XColor entries for the colormap.
   ****************************************************************************/
   colormap_size = 256;
   if((colors = (XColor *) calloc(colormap_size, sizeof(XColor))) == NULL){
      fprintf(stderr, "No memory for setting up colormap\n");
      exit(1);
   }

   /****************************************************************************
   * Get 64 grey values in the colormap. This should produce an adequate but
   * not great grayscale for display.
   ****************************************************************************/
   for(i=0;i<128;i++){
      if(get_a_color(colors, i, 2*i, 2*i, 2*i) == 0){
         warnings++;
      }
      else{
         cell[2*i + 1] = colors[i].pixel;
         cell[2*i] = colors[i].pixel;
      }
   }

   if(get_a_color(colors, REDINDEX,    255,   0,  18) == 0) exit(1);
   if(get_a_color(colors, ORANGEINDEX, 255, 149,   0) == 0) exit(1);
   if(get_a_color(colors, YELLOWINDEX, 247, 255,   6) == 0) exit(1);
   if(get_a_color(colors, GREENINDEX,   14, 255,   0) == 0) exit(1);
   if(get_a_color(colors, BLUEINDEX,     0, 112, 238) == 0) exit(1);
   if(get_a_color(colors, INDIGOINDEX, 164,   0, 255) == 0) exit(1);
   if(get_a_color(colors, VIOLETINDEX, 230,   0, 255) == 0) exit(1);

   if(get_a_color(colors, CHAINCOLOR1, 255,   0,   0) == 0) exit(1);
   if(get_a_color(colors, CHAINCOLOR2,   0, 255,   0) == 0) exit(1);
   if(get_a_color(colors, CHAINCOLOR3,   0,   0, 255) == 0) exit(1);
   if(get_a_color(colors, CHAINCOLOR4,   0, 255, 255) == 0) exit(1);
   if(get_a_color(colors, CHAINCOLOR5, 255,   0, 255) == 0) exit(1);
   if(get_a_color(colors, CHAINCOLOR6, 164,   4, 255) == 0) exit(1);
   if(get_a_color(colors, CHAINCOLOR7,  60, 156, 255) == 0) exit(1);
   if(get_a_color(colors, CHAINCOLOR8, 255, 127,  26) == 0) exit(1);
   if(get_a_color(colors, CHAINCOLOR9, 255, 36,  132) == 0) exit(1);

   /****************************************************************************
   * Create the colormap.
   ****************************************************************************/
   theColormap = XCreateColormap(theDisplay, theRootWindow, theVisual, AllocAll);

   /****************************************************************************
   * Store the colors into the colormap.
   ****************************************************************************/
   XStoreColors(theDisplay, theColormap, colors, colormap_size);

   /* free(colors); */
}

/*******************************************************************************
* Procedure: get_a_color
* Purpose: To reserve a color in the colormap.
* Name: Mike Heath
* Date: 12/1/97
*******************************************************************************/
int get_a_color(XColor *map, int index, int redval, int grnval, int bluval)
{
   map[index].red    = redval * 65536L/256L;
   map[index].green  = grnval * 65536L/256L;
   map[index].blue   = bluval * 65536L/256L;
   map[index].pixel  = index;
   map[index].flags  = DoRed | DoGreen | DoBlue;

   return(1);
}

Window openWindow(int x, int y, int width, int height, int flag,
    char theTitle[], Window theParent, GC *theNewGC)
{
   XSetWindowAttributes    theWindowAttributes;
   XSizeHints              theSizeHints;
   unsigned        long    theWindowMask;
   Window                  theNewWindow;
   XWMHints                theWMHints;
   int screen;

   /*
   **      1) Set up the attributes desired for the window.
   **      Note that window managers may deny us some of
   **      these resources. 
   */

   theWindowAttributes.border_pixel      = theBlackPixel;
   theWindowAttributes.background_pixel  = cell[128];
   theWindowAttributes.colormap  = theColormap;

   if ( flag == POP_UP_WINDOW ) {
      theWindowAttributes.override_redirect = True;       
      theWindowAttributes.save_under        = True;       
      theWindowMask = ( CWColormap |
			CWBackPixel   |
                        CWBorderPixel |
                        CWSaveUnder   |
                        CWOverrideRedirect );
   }
   else{
      theWindowMask = CWColormap | CWBackPixel   | CWBorderPixel;
   }

   /* 2) Open a window on the display.  */
/*
   printf("In openWindow() we have:\n");
   printf("   visual.depth = %d\n", xvis_info.depth);
   printf("   visual.screen = %d\n", xvis_info.screen);
   printf("   visual.class = %d\n", xvis_info.class);
   printf("   visual.colormap_size = %d\n", xvis_info.colormap_size);
   printf("   visual.bits_per_rgb = %d\n", xvis_info.bits_per_rgb);

     typedef struct {
          Visual *visual;
          VisualID visualid;
          int screen;
          unsigned int depth;
          int class;
          unsigned long red_mask;
          unsigned long green_mask;
          unsigned long blue_mask;
          int colormap_size;
          int bits_per_rgb;
     } XVisualInfo;
*/

   theNewWindow = XCreateWindow(theDisplay, theParent, x, y,
                     (unsigned int)width, (unsigned int)height, 
                     (unsigned int)BORDER_WIDTH, xvis_info.depth, InputOutput,
                     theVisual, theWindowMask,
                     &theWindowAttributes);

   /*
   **      3) Send "Hints" to the Window Manager.
   **      Before this window will appear on the display,
   **      an X window manager may intecept the call and
   **      place the window where it wants to.  This next section
   **      tells the window manager "hints" as to where the
   **      window should go.
   */

   theWMHints.initial_state = NormalState;
   theWMHints.flags = StateHint;

   XSetWMHints(theDisplay, theNewWindow, &theWMHints);

   /* 4) Store the Window.  */
   XStoreName(theDisplay, theNewWindow, theTitle);

   /*
   **      5) Now tell the window manager about the size and location.
   **      we want for our windows.  USPosition means we are
   **      stating the User choose the position, same with the
   **      size.  PPosition and PSize would mean that the program
   **      choose the size.
   */

   theSizeHints.flags      = USPosition | USSize;
   theSizeHints.x          = x;
   theSizeHints.y          = y;
   theSizeHints.width      = width;
   theSizeHints.height     = height;

   XSetNormalHints(theDisplay, theNewWindow, &theSizeHints);

   /*
   **      6) Create a graphics context for the window.
   */

   if(createGC(theNewWindow, theNewGC) == 0){
      XDestroyWindow(theDisplay, theNewWindow);
      return((Window)0);
   }

   /*
   **      7) Ask X to place the window visibly on the screen.
   **      Up to now, the window has been created but has not
   **      appeared on the screen. Mapping the window places it
   **      visibly on the screen.
   */

   /* XMapWindow( theDisplay, theNewWindow ); */

   /*
   **      8) Flush out all the queued up X requests to the X server
   */
   /* XFlush( theDisplay ); */

   /*
   **      9) Return the window ID, which is needed to specify
   **      which window to draw to.
   */

   return(theNewWindow);
}

/*
**      createGC creates a graphics context for the given window.
**      A graphics context is necessary to draw into the window.
**
**      Returns 0 if there was an error, 1 if all is A-OK.
*/

createGC(Drawable theNewWindow, GC *theNewGC)
{
   XGCValues       theGCValues;

   *theNewGC = XCreateGC(theDisplay, theNewWindow, (unsigned long) 0, &theGCValues);

   if(*theNewGC == 0){         /* -- Unable to create a GC */
      return( 0 );             /* -- Error                 */
   }
   else{
      /* --  Set Foreground and Background defaults for the new GC */

      XSetForeground(theDisplay, *theNewGC, theBlackPixel);

      XSetBackground(theDisplay, *theNewGC, cell[128]);

      return(1);    /* -- A-OK                  */
   }
}

initEvents(Window theWindow, int inPalette)
{

   if ( inPalette == IN_PALETTE ) {
      XSelectInput(theDisplay, theWindow, EV_MASK);
   }
   else{
      XSelectInput(theDisplay, theWindow,
                   ( EV_MASK | ButtonMotionMask | ButtonReleaseMask ) );
                /* ( EV_MASK | PointerMotionMask ) ); */
   }
}

/***************************************************************************
* FUNCTION: display_image
* PURPOSE: This function places an image on the screen. All of the window
* information must already be in the structure. This function only pastes
* an image in an already existing window.
***************************************************************************/
void display_image(IMAGEXWIN *animage, int rows, int cols, int remap)
{
   XImage *theImage;
   int r, c;
   unsigned char thebyte;

   if(remap == 1){
      for(r=0;r<rows;r++)
         for(c=0;c<cols;c++) animage->image_dsp[r*cols+c] = cell[animage->image[r*cols+c]];
   }

   animage->theImage=XCreateImage(theDisplay, xvis_info.visual, xvis_info.depth, ZPixmap,
             0, animage->image_dsp, cols, rows, 8, 0);

   /* animage->theImage->byte_order = LSBFirst; */
   animage->theImage->bits_per_pixel = 8;
   animage->theImage->bytes_per_line = cols * (8 / 8);

   if(animage->have_a_pixmap == 1){
      XFreePixmap(theDisplay, animage->thePixmap);
      animage->have_a_pixmap = 0;
   }

   animage->thePixmap = XCreatePixmap(theDisplay, animage->theDrawWindow, cols, rows, 8);
   animage->have_a_pixmap = 1;

   XPutImage(theDisplay, animage->thePixmap, animage->theGC,
             animage->theImage, 0, 0, 0, 0, cols, rows);

   repaint(*animage, 0, 0, cols, rows, 0, 0);

}

repaint(IMAGEXWIN animage, int sx, int sy, int width, int height, int dx, int dy)
{
   XCopyArea(theDisplay, animage.thePixmap, animage.theDrawWindow, animage.theGC, sx,
      sy, width, height, dx, dy);
   XFlush(theDisplay);
}


/*******************************************************************************
* Initialize the colormap.
*******************************************************************************/
initGrayScalecolormap_8bit()
{
   int i, j, theColormapSize, good;
   int warnings=0;
   int get_a_single_color(XColor *, int , int , int , int );

   /****************************************************************************
   * Since we are running an a workstation that cannot support 24-bit mode,
   * we assume that there is a default colormap out that and that we can allocate
   * individual entries in it.
   ****************************************************************************/
   theColormap = DefaultColormap(theDisplay, theScreen);

   /****************************************************************************
   * Create the XColor entries for the colormap.
   ****************************************************************************/
   colormap_size = 256;
   if((colors = (XColor *) calloc(colormap_size, sizeof(XColor))) == NULL){
      fprintf(stderr, "No memory for setting up colormap\n");
      exit(1);
   }

   /****************************************************************************
   * Get come colors.
   ****************************************************************************/
   if(get_a_single_color(colors, REDINDEX,    255,   0,  18) == 0) exit(1);
   if(get_a_single_color(colors, ORANGEINDEX, 255, 149,   0) == 0) exit(1);
   if(get_a_single_color(colors, YELLOWINDEX, 247, 255,   6) == 0) exit(1);
   if(get_a_single_color(colors, GREENINDEX,   14, 255,   0) == 0) exit(1);
   if(get_a_single_color(colors, BLUEINDEX,     0, 112, 238) == 0) exit(1);
   if(get_a_single_color(colors, INDIGOINDEX, 164,   0, 255) == 0) exit(1);
   if(get_a_single_color(colors, VIOLETINDEX, 230,   0, 255) == 0) exit(1);

   if(get_a_single_color(colors, CHAINCOLOR1, 255,   0,   0) == 0) exit(1);
   if(get_a_single_color(colors, CHAINCOLOR2,   0, 255,   0) == 0) exit(1);
   if(get_a_single_color(colors, CHAINCOLOR3,   0,   0, 255) == 0) exit(1);
   if(get_a_single_color(colors, CHAINCOLOR4,   0, 255, 255) == 0) exit(1);
   if(get_a_single_color(colors, CHAINCOLOR5, 255,   0, 255) == 0) exit(1);
   if(get_a_single_color(colors, CHAINCOLOR6, 164,   4, 255) == 0) exit(1);
   if(get_a_single_color(colors, CHAINCOLOR7,  60, 156, 255) == 0) exit(1);
   if(get_a_single_color(colors, CHAINCOLOR8, 255, 127,  26) == 0) exit(1);
   if(get_a_single_color(colors, CHAINCOLOR9, 255, 36,  132) == 0) exit(1);

   /****************************************************************************
   * Get 64 grey values in the colormap. This should produce an adequate but
   * not great grayscale for display.
   ****************************************************************************/
   for(i=0;i<256;i+=4){
      if(get_a_single_color(colors, i/2, i, i, i) == 0){
         warnings++;
      }
      else{
         cell[i]   = colors[i/2].pixel;
         cell[i+1] = colors[i/2].pixel;
         cell[i+2] = colors[i/2].pixel;
         cell[i+3] = colors[i/2].pixel;
      }
   }
   if(warnings != 0){
      fprintf(stderr, "\nError allocating 64 grey levels for display. The display ");
      fprintf(stderr, "quality will likely\nbe very poor so the program is ");
      fprintf(stderr, "terminating. If you are running other\napplications, ");
      fprintf(stderr, "you may want to quit them and then re-run this program.\n\n");
      exit(1);
   }

   /****************************************************************************
   * Get 64 more grey values in the colormap. This should produce a better
   * quality grayscale for images.
   ****************************************************************************/
   warnings = 0;
   for(i=2;i<256;i+=4){
      if(get_a_single_color(colors, i/2, i, i, i) == 0){
         warnings++;
      }
      else{
         cell[i]   = colors[i/2].pixel;
         cell[i+1] = colors[i/2].pixel;
      }
   }

   if(warnings == 0) printf("You got 128 grey levels for displaying the image.\n");
   else printf("You got %d levels for displaying the image.\n", 128-warnings);
}

/*******************************************************************************
* Procedure: get_a_single_color
* Purpose: To reserve a color in the colormap.
* Name: Mike Heath
* Date: 12/1/97
*******************************************************************************/
int get_a_single_color(XColor *map, int index, int redval, int grnval, int bluval)
{
   map[index].red = redval * 65536L/256L;
   map[index].green = grnval * 65536L/256L;
   map[index].blue = bluval * 65536L/256L;
   if(XAllocColor(theDisplay, theColormap, &(map[index]))==0){
      printf("Warning: Cannot initialize the color (%d, %d, %d).\n", redval, grnval, bluval);
      return(0);
   }
   return(1);
}
