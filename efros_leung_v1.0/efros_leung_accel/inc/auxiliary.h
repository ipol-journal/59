#pragma once

#include "io_png.h"
#include "mt.h"

/**
 * @file    auxiliary.h
 * @brief   Image structures definitions.
 * 
 * Definitions of common structures like images or coordinates.
 * 
 * @author  Cecilia Aguerrebere
 */

/**
 * @struct  cimage
 * @brief   A color image of integers.
 *
 * Values are integers in \f$0\dots 255\f$.
 * Can have up to three channels.
 */
typedef struct cimage
{
  int nrow;                 /**< @brief Number of rows (dy) */
  int ncol;                 /**< @brief Number of columns (dx) */
  int nchannels;            /**< @brief Number of channels */

  unsigned char *red;       /**< @brief The red level plane (may be NULL) */
  unsigned char *green;     /**< @brief The green level plane (may be NULL) */
  unsigned char *blue;      /**< @brief The blue level plane (may be NULL) */
} * Cimage;

/**
 * @struct  imageF
 * @brief   A gray-scale image of floats.
 */
typedef struct imageF
{
  int nrow;                 /**< @brief Number of rows (dy) */
  int ncol;                 /**< @brief Number of columns (dx) */

  float *val;               /**< @brief The red level plane (may be NULL) */

} * ImageF;

/**
 * @struct  pixel
 * @brief   Some pixel's properties.
 */
typedef struct pixel
{
  int x;                    /**< @brief X position in the image */
  int y;                    /**< @brief Y position in the image */
  int nbneighbs;            /**< @brief Number of known neighbors */
} * Pixel;

/**
 * @struct  cand_dist
 * @brief   Store a patch and a distance
 */
typedef struct cand_dist
{
  float dist;               /**< @brief Distance to another patch */
  long int patch;           /**< @brief Index of a patch in an array */    

} * Cand_dist;


void error(char * msg);
Cimage new_Cimage(int nrow, int ncol);
ImageF new_ImageF(int nrow, int ncol);
Cimage load_png_image(char *filename);
void write_png_image(Cimage v, char * filename);
void delete_Cimage(Cimage v);
void delete_imageF(ImageF v);
ImageF CImage_to_ImageF(Cimage image);
Cimage create_image_map(int nrow, int ncol);
void print_algo_info( char *in_file_name, int nrows, int ncols, int nchannels, int t, float tolerance, int dims_pca, int out_sz);
Cimage crop_image( Cimage v, int out_img_sz );
int is_grayscale(Cimage image, long num_pixels);
ImageF add_border( ImageF v, int t );
ImageF create_mask ( int nrow, int ncol, int t );

