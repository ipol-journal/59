#pragma once

#include "io_png.h"
#include "mt.h"

/*----------------------------------------------------------------------------*/
/* Structures definition */
/*----------------------------------------------------------------------------*/

typedef struct cimage
{
  int nrow;        /* Number of rows (dy) */
  int ncol;        /* Number of columns (dx) */
  int nchannels;        /* Number of channels */

  unsigned char *red;     /* The red level plane (may be NULL) */
  unsigned char *green;   /* The green level plane (may be NULL) */
  unsigned char *blue;    /* The blue level plane (may be NULL) */
} * Cimage;

typedef struct gimage
{
  int nrow;        /* Number of rows (dy) */
  int ncol;        /* Number of columns (dx) */
  int nchannels;        /* Number of channels */

  unsigned char *val;     /* The graylevel plane (may be NULL) */
} * Gimage;

typedef struct pixel
{
  int x;
  int y;
  int nbneighbs;
} * Pixel;

typedef struct cand_dist
{
  float dist;
  long int patch;

} * Cand_dist;


void error(char * msg);
void write_png_image(Cimage v, char * filename);
void write_png_image_gray(Gimage image, char * filename);
void delete_Cimage(Cimage v);
void delete_Gimage(Gimage v);
void print_algo_info( char *in_file_name, int nrows, int ncols, int nchannels, int t, float tolerance, int out_sz);

Cimage new_Cimage(int nrow, int ncol);
Gimage new_Gimage(int nrow, int ncol);
Cimage create_image_map(int nrow, int ncol);
Cimage crop_image( Cimage v, int out_img_sz );
Gimage crop_image_gray( Gimage v, int out_img_sz );
Cimage add_border( Cimage v, int t );
Gimage add_border_gray( Gimage v, int t );
Cimage create_mask ( int nrow, int ncol, int t );
Gimage create_mask_gray( int nrow, int ncol, int t );
int is_grayscale(Cimage image, long num_pixels);
Cimage load_png_image(char *filename);
