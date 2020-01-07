/*******************************************************************************
* Name: Mike Heath
* File: Poly.c
* Purpose: This file contains a procedure for calculating the statistics of
* a polygonal region in a "virtual_image".
* Date: 1/8/98
*******************************************************************************/
#include "virtual_image.h"

#define VERBOSE 0

typedef struct{
   double x2d, y2d;                /* The 2D projected point */
}VERTEX;

typedef struct{
   int num;                        /* The number of vertices in the polygon */
   VERTEX *vertices;               /* A list of the vertices for the polygon */
}POLYGON;

int insidepoly(double x, double y, POLYGON poly);
double angle_radians(double x, double y);

/*******************************************************************************
* Function: polygon_stats
* Purpose: To compute the statistics for a polygonal area of an image.
*******************************************************************************/
void polygon_stats(CACHEIM *image, int *x_coord, int *y_coord, int num_points,
   int *count, double *mean, double *stdev)
{
   POLYGON poly;                     /* A structure used to hold the object */
   int min_x_coord, max_x_coord, min_y_coord, max_y_coord;
   int i, r, c;
   double pixvalue;

   /****************************************************************************
   * Find the bounding box of the polygon.
   ****************************************************************************/
   min_x_coord = max_x_coord = x_coord[0];
   for(i=0;i<num_points;i++){
      if(x_coord[i] < min_x_coord) min_x_coord = x_coord[i];
      if(x_coord[i] > max_x_coord) max_x_coord = x_coord[i];
   }
   min_y_coord = max_y_coord = y_coord[0];
   for(i=0;i<num_points;i++){
      if(y_coord[i] < min_y_coord) min_y_coord = y_coord[i];
      if(y_coord[i] > max_y_coord) max_y_coord = y_coord[i];
   }

   /****************************************************************************
   * Fill in the vertices into the POLYGON structure;
   ****************************************************************************/
   poly.num = num_points;
   if((poly.vertices = (VERTEX *) calloc(num_points, sizeof(VERTEX))) == NULL){
      fprintf(stderr, "Malloc error!\n");
      return;
   }
   if(VERBOSE) printf("The Polygon is: ");
   for(i=0;i<num_points;i++){
      poly.vertices[i].x2d = (double)x_coord[i];
      poly.vertices[i].y2d = (double)y_coord[i];
      if(VERBOSE) printf(" (%d, %d)", x_coord[i], y_coord[i]);
   }
   if(VERBOSE) printf("\n");

   /****************************************************************************
   * Go through each pixel in the bounding box and make it contribute to the
   * statistics if it is inside the polygon.
   ****************************************************************************/
   *count = 0;
   *mean = 0.0;
   *stdev = 0.0;
   for(r=min_y_coord;r<=max_y_coord;r++){
      for(c=min_x_coord;c<max_x_coord;c++){
         if(insidepoly((double)c, (double)r, poly)){

            *count += 1;

            pixvalue = (*(image->getpixel))(image, r, c);

            *mean += pixvalue;
            *stdev += pixvalue * pixvalue;
         }
      }
   }

   *stdev = (*stdev - (((*mean)*(*mean))/(double)(*count))) / (double)((*count)-1);
   *stdev = sqrt(*stdev);
   *mean /= (double)(*count);

   if(VERBOSE) printf("Count = %d, Mean = %lf, Stdev = %lf\n", *count, *mean, *stdev);

   free(poly.vertices);
}

/*******************************************************************************
* Function: insidepoly
* Purpose: Determine if the point (x,y) is inside the polygon.
*******************************************************************************/
int insidepoly(double x, double y, POLYGON poly)
{
   int v, numv;
   double v0x, v0y, v1x, v1y, length0, length1, dotproduct, theta = 0.0, jnk;
   double crossproduct, temptheta;

   numv = poly.num;

   v0x = poly.vertices[0].x2d - x;
   v0y = poly.vertices[0].y2d - y;
   length0 = sqrt(v0x*v0x + v0y*v0y);
   if(length0 == 0.0) return(0);

   for(v=1;v<numv;v++){

      v1x = poly.vertices[v].x2d - x;
      v1y = poly.vertices[v].y2d - y;
      length1 = sqrt(v1x*v1x + v1y*v1y);
      if(length1 == 0.0) return(0);

      dotproduct = (v0x * v1x + v0y * v1y)/(length0 * length1);
      crossproduct = (v0x*v1y - v1x*v0y)/(length0 * length1);
      temptheta = angle_radians(dotproduct, crossproduct);
      if(temptheta > M_PI) temptheta -= 2 * M_PI;

      theta += temptheta;

      v0x = v1x;
      v0y = v1y;
      length0 = length1;
   }

   /****************************************************************************
   * If the last point in the vertex list was not the first point in the vertex
   * list, then we must add in the last bit of the polygon. (i.e. the first
   * point).
   ****************************************************************************/
   if( !((poly.vertices[numv-1].x2d == poly.vertices[0].x2d) &&
         (poly.vertices[numv-1].y2d == poly.vertices[0].y2d)) ){

      /* printf("The polygon is not closed (the algorithm is closing it.).\n"); */

      v1x = poly.vertices[0].x2d - x;
      v1y = poly.vertices[0].y2d - y;
      length1 = sqrt(v1x*v1x + v1y*v1y);
      if(length1 == 0.0) return(0);

      dotproduct = (v0x * v1x + v0y * v1y)/(length0 * length1);
      crossproduct = (v0x*v1y - v1x*v0y)/(length0 * length1);
      temptheta = angle_radians(dotproduct, crossproduct);
      if(temptheta > M_PI) temptheta -= 2 * M_PI;

      theta += temptheta;
   }
   else{
      printf("The polygon was closed (i.e. the last_vertex = first_vertex).\n");
   }

   theta = fabs(theta);

   if(fabs(theta) < 1.0e-3) return(0);
   else return(1);

   /*
   if(fabs(2*M_PI - theta) < 1.0e-3) return(1);
   else return(0);
   */

   /*
   if(modf(fabs(theta), &jnk) < 1.0e-5) return(1);
   else return(0);
   */
}

/*******************************************************************************
* FUNCTION: angle_radians
* PURPOSE: This procedure computes the angle of a vector with components x and
* y. It returns this angle in radians with the answer being in the range
* 0 <= angle <2*PI.
*******************************************************************************/
/*
double angle_radians(double x, double y)
{
   double xu, yu, ang;

   xu = fabs(x);
   yu = fabs(y);

   if((xu == 0) && (yu == 0)) return(0);

   ang = atan(yu/xu);

   if(x >= 0){
      if(y >= 0) return(ang);
      else return(2*M_PI - ang);
   }
   else{
      if(y >= 0) return(M_PI - ang);
      else return(M_PI + ang);
   }
}
*/
