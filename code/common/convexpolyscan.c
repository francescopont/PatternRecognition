/*
 * Concave Polygon Scan Conversion
 * by Paul Heckbert
 * from "Graphics Gems", Academic Press, 1990
 */

/*******************************************************************************
 * concave: scan convert nvert-sided concave non-simple polygon with vertices at
 * (point[i].x, point[i].y) for i in [0..nvert-1] within the window win by
 * calling spanproc for each visible span of pixels.
 * Polygon can be clockwise or counterclockwise.
 * Algorithm does uniform point sampling at pixel centers.
 * Inside-outside test done by Jordan's rule: a point is considered inside if
 * an emanating ray intersects the polygon an odd number of times.
 * drawproc should fill in pixels from xl to xr inclusive on scanline y,
 * e.g:
 *	drawproc(y, xl, xr)
 *	int y, xl, xr;
 *	{
 *	    int x;
 *	    for (x=xl; x<=xr; x++)
 *		pixel_write(x, y, pixelvalue);
 *	}
 *
 *  Paul Heckbert	30 June 81, 18 Dec 89
 ******************************************************************************/

#include <stdio.h>
#include "graphicsgems.h"
#include "virtual_image.h"
#include "myalloc.h"

/* #define TESTCODE (Uncomment this command to compile a test program/) */

#define ALLOC(ptr, type, n)  ASSERT(ptr = (type *)calloc((n),sizeof(type)))

typedef struct {		/* window: a discrete 2-D rectangle */
    int x0, y0;			/* xmin and ymin */
    int x1, y1;			/* xmax and ymax (inclusive) */
} Window;

typedef struct {		/* a polygon edge */
    double x;	/* x coordinate of edge's intersection with current scanline */
    double dx;	/* change in x with respect to y */
    int i;	/* edge number: edge i goes from pt[i] to pt[i+1] */
} Edge;

/******************************************************************************/
/* Michael Heath's code.                                                               */
/******************************************************************************/
typedef struct{
   double x2d, y2d;                /* The 2D projected point */
}VERTEX;

typedef struct{
   int num;                        /* The number of vertices in the polygon */
   VERTEX *vertices;               /* A list of the vertices for the polygon */
}POLYGON;
/******************************************************************************/

static int n;			/* number of vertices */
static Point2 *pt;		/* vertices */

static int nact;		/* number of active edges */
static Edge *active;		/* active edge list:edges crossing scanline y */

static unsigned char **imageptr=NULL;
static CACHEIM *cache_image_ptr=NULL;

int compare_ind(), compare_active();
void polyscan_poly(POLYGON poly, unsigned char *image, int rows, int cols);
void polyscan_coords(int numpoints, int *x_coord, int *y_coord, unsigned char *image, int rows, int cols);
void polyscan_poly_cacheim(POLYGON poly, CACHEIM *cache_image);
void fill_line(int y, int xl, int xr);
void fill_cacheim_line(int y, int xl, int xr);

static cdelete(int); /* remove edge i from active list */
static cinsert(int, int);		/* append edge i to end of active list */


#ifdef TESTCODE
void main()
{
   FILE *fp=NULL;
   unsigned char *image=NULL;
   int i;
   int rows=512, cols=512;
   int numpoints=6, *x_coord, *y_coord;
   double fraction;

   x_coord = (int *) mycalloc("x_coord", numpoints, sizeof(int));
   y_coord = (int *) mycalloc("y_coord", numpoints, sizeof(int));

   image = (unsigned char *) mycalloc("image", rows*cols, 1);

   /*
   for(i=0;i<numpoints;i++){
      fraction = (double)i / (double)numpoints;
      x_coord[i] = fabs(cos(fraction * 6 * 3.1415926535)) * (cols / 4) * cos(fraction * 2 * 3.1415926535) + cols/2;
      y_coord[i] = fabs(cos(fraction * 6 * 3.1415926535)) * (rows / 4) * sin(fraction * 2 * 3.1415926535) + rows/2;
      printf("   %d  %d\n", x_coord[i], y_coord[i]);
   }
   */

   x_coord[0] =  100 ;   y_coord[0] =  100  ;
   x_coord[1] =  200 ;   y_coord[1] =  150  ;
   x_coord[2] =  300 ;   y_coord[2] =  100  ;
   x_coord[3] =  300 ;   y_coord[3] =  300  ;
   x_coord[4] =  200 ;   y_coord[4] =  100  ;
   x_coord[5] =  100 ;   y_coord[5] =  300  ;

   polyscan_coords(numpoints, x_coord, y_coord, image, rows, cols);

   fp = fopen("polytest.pgm", "wb");
   fprintf(fp, "P5\n%d %d\n255\n", cols, rows);
   fwrite(image, 1, rows*cols, fp);
   fclose(fp);

   myfree("image", image);
}
#endif

/*******************************************************************************
* Function: polyscan_poly
* Purpose: This code rasterizes a polygon, setting all pixels inside the
* polygon to have the value 255. The polycon can be concave (not all points
* connecting every ppossible pair of points inside the polygon need to be inside
* the polygon), can self intersect and the vertices can be specified in either
* clockwise or counterclockwise order. It does not matter. This routing just
* formats the data. The function concave does all of the work. That function
* came from the book "Graphics Gems" by Glassner.
* Name: Michael Heath, University of South Florida
* Date: 7/8/98
*******************************************************************************/
void polyscan_poly(POLYGON poly, unsigned char *image, int rows, int cols)
{
   int i;
   int nvert;			/* number of vertices */
   Point2 *point;			/* vertices of polygon */
   Window win;			/* screen clipping window */
   void (*spanproc)();		/* called for each span of pixels */

   nvert = poly.num;

   point = (Point2 *) calloc(nvert, sizeof(Point2));

   for(i=0;i<nvert;i++){
      point[i].x = poly.vertices[i].x2d;
      point[i].y = poly.vertices[i].y2d;
   }

   win.x0 = 0;
   win.y0 = 0;
   win.x1 = cols - 1;
   win.y1 = rows - 1;

   imageptr = (unsigned char **) calloc(rows, sizeof(unsigned char *));
   for(i=0;i<rows;i++) imageptr[i] = image + i*cols;

   concave(nvert, point, &win, fill_line);

   free(imageptr);
   free(point);
}

/*******************************************************************************
* Function: polyscan_poly_cacheim
* Purpose: This code rasterizes a polygon, setting all pixels inside the
* polygon to have the value 255. The polycon can be concave (not all points
* connecting every ppossible pair of points inside the polygon need to be inside
* the polygon), can self intersect and the vertices can be specified in either
* clockwise or counterclockwise order. It does not matter. This routing just
* formats the data. The function concave does all of the work. That function
* came from the book "Graphics Gems" by Glassner.
* Name: Michael Heath, University of South Florida
* Date: 7/8/98
*******************************************************************************/
void polyscan_poly_cacheim(POLYGON poly, CACHEIM *cache_image)
{
   int rows, cols;
   int i;
   int nvert;			/* number of vertices */
   Point2 *point;			/* vertices of polygon */
   Window win;			/* screen clipping window */
   void (*spanproc)();		/* called for each span of pixels */

   rows = cache_image->rows;
   cols = cache_image->cols;

   nvert = poly.num;

   point = (Point2 *) mycalloc("point", nvert, sizeof(Point2));

   for(i=0;i<nvert;i++){
      point[i].x = poly.vertices[i].x2d;
      point[i].y = poly.vertices[i].y2d;
   }

   win.x0 = 0;
   win.y0 = 0;
   win.x1 = cols - 1;
   win.y1 = rows - 1;

   cache_image_ptr = cache_image;

   concave(nvert, point, &win, fill_cacheim_line);

   cache_image_ptr = (CACHEIM *)NULL;
   myfree("point", point);
}

/*******************************************************************************
* Function: polyscan_coords
* Purpose: This code rasterizes a polygon, setting all pixels inside the
* polygon to have the value 255. The polycon can be concave (not all points
* connecting every possible pair of points inside the polygon need to be inside
* the polygon), can self intersect and the vertices can be specified in either
* clockwise or counterclockwise order. It does not matter. This routing just
* formats the data. The function concave does all of the work. That function
* came from the book "Graphics Gems" by Glassner.
* Name: Michael Heath, University of South Florida
* Date: 7/8/98
*******************************************************************************/
void polyscan_coords(int numpoints, int *x_coord, int *y_coord, unsigned char *image, int rows, int cols)
{
   int i;
   int nvert;			/* number of vertices */
   Point2 *point;			/* vertices of polygon */
   Window win;			/* screen clipping window */
   void (*spanproc)();		/* called for each span of pixels */

   nvert = numpoints;

   point = (Point2 *) mycalloc("point", nvert, sizeof(Point2));

   /* printf("rows = %d  cols = %d\n", rows, cols); */

   for(i=0;i<nvert;i++){
      point[i].x = x_coord[i];
      point[i].y = y_coord[i];

      /* printf("Points[%d] = (%f %f)\n", i, point[i].y, point[i].x); */
   }

   win.x0 = 0;
   win.y0 = 0;
   win.x1 = cols - 1;
   win.y1 = rows - 1;

   imageptr = (unsigned char **) mycalloc("imageptr", rows, sizeof(unsigned char *));
   for(i=0;i<rows;i++) imageptr[i] = image + i*cols;

   concave(nvert, point, &win, fill_line);

   myfree("imageptr", imageptr);
   myfree("point", point);
}

void fill_line(int y, int xl, int xr)
{
   int x;
   for(x=xl;x<=xr;x++) imageptr[y][x] = 255;
}

void fill_cacheim_line(int y, int xl, int xr)
{
   int x;
   for(x=xl;x<=xr;x++) (*(cache_image_ptr->putpixel))(cache_image_ptr, y, x, 255);
}

concave(nvert, point, win, spanproc)
int nvert;			/* number of vertices */
Point2 *point;			/* vertices of polygon */
Window *win;			/* screen clipping window */
void (*spanproc)();		/* called for each span of pixels */
{
    int k, y0, y1, y, i, j, xl, xr;
    int *ind;		/* list of vertex indices, sorted by pt[ind[j]].y */

    n = nvert;
    pt = point;
    if (n<=0) return;
    ALLOC(ind, int, n);
    ALLOC(active, Edge, n);

    /* create y-sorted array of indices ind[k] into vertex list */
    for (k=0; k<n; k++)
	ind[k] = k;
    qsort(ind, n, sizeof ind[0], compare_ind);	/* sort ind by pt[ind[k]].y */

    nact = 0;				/* start with empty active list */
    k = 0;				/* ind[k] is next vertex to process */
    y0 = MAX(win->y0, ceil(pt[ind[0]].y-.5));		/* ymin of polygon */
    y1 = MIN(win->y1, floor(pt[ind[n-1]].y-.5));	/* ymax of polygon */

    for (y=y0; y<=y1; y++) {		/* step through scanlines */
	/* scanline y is at y+.5 in continuous coordinates */

	/* check vertices between previous scanline and current one, if any */
	for (; k<n && pt[ind[k]].y<=y+.5; k++) {
	    /* to simplify, if pt.y=y+.5, pretend it's above */
	    /* invariant: y-.5 < pt[i].y <= y+.5 */
	    i = ind[k];	
	    /*
	     * insert or delete edges before and after vertex i (i-1 to i,
	     * and i to i+1) from active list if they cross scanline y
	     */
	    j = i>0 ? i-1 : n-1;	/* vertex previous to i */
	    if (pt[j].y <= y-.5)	/* old edge, remove from active list */
		cdelete(j);
	    else if (pt[j].y > y+.5)	/* new edge, add to active list */
		cinsert(j, y);
	    j = i<n-1 ? i+1 : 0;	/* vertex next after i */
	    if (pt[j].y <= y-.5)	/* old edge, remove from active list */
		cdelete(i);
	    else if (pt[j].y > y+.5)	/* new edge, add to active list */
		cinsert(i, y);
	}

	/* sort active edge list by active[j].x */
	qsort(active, nact, sizeof active[0], compare_active);

	/* draw horizontal segments for scanline y */
	for (j=0; j<nact; j+=2) {	/* draw horizontal segments */
	    /* span 'tween j & j+1 is inside, span tween j+1 & j+2 is outside */
	    xl = ceil(active[j].x-.5);		/* left end of span */
	    if (xl<win->x0) xl = win->x0;
	    xr = floor(active[j+1].x-.5);	/* right end of span */
	    if (xr>win->x1) xr = win->x1;
	    if (xl<=xr)
		(*spanproc)(y, xl, xr);		/* draw pixels in span */
	    active[j].x += active[j].dx;	/* increment edge coords */
	    active[j+1].x += active[j+1].dx;
	}
    }
}

static cdelete(i)		/* remove edge i from active list */
int i;
{
    int j;

    for (j=0; j<nact && active[j].i!=i; j++);
    if (j>=nact) return;	/* edge not in active list; happens at win->y0*/
    nact--;
    bcopy(&active[j+1], &active[j], (nact-j)*sizeof active[0]);
}

static cinsert(i, y)		/* append edge i to end of active list */
int i, y;
{
    int j;
    double dx;
    Point2 *p, *q;

    j = i<n-1 ? i+1 : 0;
    if (pt[i].y < pt[j].y) {p = &pt[i]; q = &pt[j];}
    else		   {p = &pt[j]; q = &pt[i];}
    /* initialize x position at intersection of edge with scanline y */
    active[nact].dx = dx = (q->x-p->x)/(q->y-p->y);
    active[nact].x = dx*(y+.5-p->y)+p->x;
    active[nact].i = i;
    nact++;
}

/* comparison routines for qsort */
compare_ind(u, v) int *u, *v; {return pt[*u].y <= pt[*v].y ? -1 : 1;}
compare_active(u, v) Edge *u, *v; {return u->x <= v->x ? -1 : 1;}
