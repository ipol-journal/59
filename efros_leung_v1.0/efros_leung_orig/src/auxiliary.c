#include "auxiliary.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * @file auxiliary.h
 */
 
 
/**
 * @brief   Fatal error function.
 * 
 * Print a message to standard-error output and exit.
 * 
 * @param   msg     String to be printed to standard-error output.
 */
void error(char * msg)
{
  fprintf(stderr,"Fatal Error: %s\n",msg);
  exit(EXIT_FAILURE);
}

/**
 * @brief   Create an empty color image.
 * 
 * Initialize the size of the image and allocate the memory.
 * 
 * @param   nrow    number of rows (dy) 
 * @param   ncol    number of columns (dx) 
 * 
 * @return  The created image
 * @warning Image is initialized to zero.
 */
Cimage new_Cimage(int nrow, int ncol)
{
  Cimage image;

  if( nrow == 0 || ncol == 0 ) error("new_Cimage: invalid image size.");

  image = (Cimage) malloc( sizeof(struct cimage) );
  if( image == NULL ) error("Not enough memory.");
  image->red = (unsigned char *) calloc( nrow*ncol, sizeof(unsigned char) );
  if( image->red == NULL ) error("not enough memory.");
  image->green = (unsigned char *) calloc( nrow*ncol, sizeof(unsigned char) );
  if( image->green == NULL ) error("not enough memory.");
  image->blue = (unsigned char *) calloc( nrow*ncol, sizeof(unsigned char) );
  if( image->blue == NULL ) error("not enough memory.");

  image->nrow = nrow;
  image->ncol = ncol; 
  image->nchannels = 3; 
  
  return (image);
}

/**
 * @brief   Create an empty gray level image.
 * 
 * Initialize the size of the image and allocate the memory.
 * 
 * @param   nrow    number of rows (dy) 
 * @param   ncol    number of columns (dx) 
 * 
 * @return  The created image
 * @warning Image initialized to zero.
 */
Gimage new_Gimage(int nrow, int ncol)
{
  Gimage image;

  if( nrow == 0 || ncol == 0 ) error("new_Cimage: invalid image size.");

  image = (Gimage) malloc( sizeof(struct gimage) );
  if( image == NULL ) error("Not enough memory.");
  image->val = (unsigned char *) calloc( nrow*ncol, sizeof(unsigned char) );
  if( image->val == NULL ) error("not enough memory.");
 
  image->nrow = nrow;
  image->ncol = ncol; 
  image->nchannels = 1; 
  
  return (image);
}


/**
 * @brief   Save a color image as a PNG file.
 * 
 * @param   image       color image to be saves
 * @param   filename    name (path) of the file
 */
void write_png_image(Cimage image, char * filename)
{
  
  int nrow = image->nrow;
  int ncol = image->ncol;
  int nchannel = 3;

  unsigned char *img = NULL;
  img = (unsigned char *) calloc( nrow*ncol*nchannel, sizeof(unsigned char) );

  /* copy each channel */
  memcpy(img, image->red, nrow * ncol * sizeof(unsigned char));
  memcpy(img + nrow * ncol, image->green, nrow * ncol * sizeof(unsigned char));
  memcpy(img + 2 * nrow * ncol, image->blue, nrow * ncol * sizeof(unsigned char)); 

  /* write and quit */  
  io_png_write_u8(filename, img, ncol, nrow, nchannel);
  free(img);
}

/**
 * @brief   Save a gray level image as a PNG file.
 * 
 * @param   image       gray level image to be saves
 * @param   filename    name (path) of the file
 */
void write_png_image_gray(Gimage image, char * filename)
{
  
  int nrow = image->nrow;
  int ncol = image->ncol;
  int nchannel = image->nchannels; 

  unsigned char *img = NULL;
  img = (unsigned char *) calloc( nrow*ncol*nchannel, sizeof(unsigned char) );

  /* copy each channel */
  memcpy(img, image->val, nrow * ncol * sizeof(unsigned char));

  /* write and quit */  
  io_png_write_u8(filename, img, ncol, nrow, nchannel);
  free(img);
}

/**
 * @brief   Delete a color image.
 * 
 * Free the memory used by the image.
 * 
 * @param   v   image to be deleted
 */
void delete_Cimage(Cimage v)
{
  if(v == NULL)
    return;

  if(v->red != NULL) free(v->red);
  if(v->green != NULL) free(v->green);
  if(v->blue != NULL) free(v->blue);
  
  free(v);
}

/**
 * @brief   Delete a gray level image.
 * 
 * Free the memory used by the image.
 * 
 * @param   v   image to be deleted
 */
void delete_Gimage(Gimage v)
{
  if(v == NULL)
    return;

  if(v->val != NULL) free(v->val);
  
  free(v);
}

/**
 * @brief   Create a color map.
 *  
 * @param   nrow    number of rows (dy) 
 * @param   ncol    number of columns (dx) 
 * 
 * @return  the created map
 */
Cimage create_image_map( int nrow, int ncol )
{
  int x,y,p = 0;
	
  Cimage map = new_Cimage(nrow, ncol);	
	
  for (y = 0; y < nrow; y++) 
    for (x = 0; x < ncol; x++) {
	     
      map->red[p] = (unsigned char) (x * 255 / ncol);
      map->green[p] = (unsigned char) (x * y * 255 / (ncol*nrow));
      map->blue[p] = (unsigned char) (y * 255 / nrow);
	     
      p++;
    }

  return map;
}
 
/**
 * @brief   Print information about the parameters.
 * 
 * @param t
 *      Half-size of the patches.
 * @param out_sz
 *      Size of the synthesized image.
 * @param tolerance
 *      Tolerance parameter \f$\varepsilon\f$.
 * @param dims_pca
 *      Number of PCA components used to compute distances
 *      between (the known parts of) patches.
 * @param weights
 *      If 'G', use Gaussian-weighted \f$l^2\f$ distance
 *      between patches instead of uniform \f$l^2\f$ distance.
 * @param in_file_name
 *      Example image name.
 * @param nrows
 *      Number of rows (dy) in example image.
 * @param ncols
 *      Number of columns (dx) in example image.
 * @param nchannels
 *      Number of channels in example image.
 */
void print_algo_info( char *in_file_name, int nrows, int ncols, int nchannels, int t, float tolerance, int out_sz)  

{
	
  printf("---------------------------------------------------\n");
  printf("-   Efros-Leung algorithm for texture synthesis   -\n");
  printf("---------------------------------------------------\n");
  printf("\n");
  printf("--------------------\n");
  printf("- Sample image\n");	
  printf("--------------------\n");
  printf("File = %s\n",in_file_name);
  printf("Size = %i x %i\n", nrows, ncols);
  printf("Channels = %i\n", nchannels);
  printf("\n");
  printf("--------------------\n");
  printf("- Input Parametres\n");
  printf("--------------------\n");
  printf("Patch size = %i x %i\n",2*t+1,2*t+1);
  printf("Tolerance = %.3f\n",tolerance);
  printf("Output image size = %i x %i\n",out_sz,out_sz);
  printf("\n");
  printf("--------------------\n");
  printf("- Warnings\n");
  printf("--------------------\n");
		
}

/**
 * @brief   Crop a color image.
 * 
 * @param   v           color image to be cropped
 * @param   out_img_sz  size of the new image
 * 
 * @return  the top-left part of the image
 */
Cimage crop_image( Cimage v, int out_img_sz )

{
  int i;
  Cimage cropped = new_Cimage(out_img_sz, out_img_sz);
	
  for (i = 0; i < out_img_sz; i++) {
		
    memcpy(cropped->red + i*out_img_sz, v->red + i*(out_img_sz + 1), out_img_sz * sizeof(unsigned char));
    memcpy(cropped->green + i*out_img_sz, v->green + i*(out_img_sz + 1), out_img_sz * sizeof(unsigned char));
    memcpy(cropped->blue + i*out_img_sz, v->blue + i*(out_img_sz + 1), out_img_sz * sizeof(unsigned char));
  }	
	
  return cropped;
		
}

/**
 * @brief   Crop a gray level image.
 * 
 * @param   v           gray level image to be cropped
 * @param   out_img_sz  size of the new image
 * 
 * @return  the top-left part of the image
 */
Gimage crop_image_gray( Gimage v, int out_img_sz )

{
  int i;
  Gimage cropped = new_Gimage(out_img_sz, out_img_sz);
	
  for (i = 0; i < out_img_sz; i++) {
		
    memcpy(cropped->val + i*out_img_sz, v->val + i*(out_img_sz + 1), out_img_sz * sizeof(unsigned char));
  }	
	
  return cropped;
		
}

/**
 * @brief   Add a border of zeros of width t to the input image.
 * 
 * @param   v       input image add the border.
 * @param   t       width of the border to be added.
 */
Cimage add_border( Cimage v, int t ) {

  int i,j,l,k;	
  int nrow_vb = v->nrow + 2*t;
  int ncol_vb = v->ncol + 2*t;
	
  Cimage vb = new_Cimage(nrow_vb, ncol_vb);

  for (i = 0; i < v->nrow; i++) {
	
    l = (i+t)*ncol_vb;
    k = i*v->ncol;
		
    for (j = 0; j < v->ncol; j++) {
					
      vb->red[l + j + t] = v->red[k + j];	
      vb->green[l + j + t] = v->green[k + j];	
      vb->blue[l + j + t] = v->blue[k + j];	
    }
			
  }

  return vb;
}

/**
 * @brief   Add a border of zeros of width t to the input image.
 * 
 * Idem add_border for gray level images. 
 * 
 * @see add_border
 */
Gimage add_border_gray( Gimage v, int t ) {

  int i,j,l,k;	
  int nrow_vb = v->nrow + 2*t;
  int ncol_vb = v->ncol + 2*t;
	
  Gimage vb = new_Gimage(nrow_vb, ncol_vb);

  for (i = 0; i < v->nrow; i++) {
	
    l = (i+t)*ncol_vb;
    k = i*v->ncol;
		
    for (j = 0; j < v->ncol; j++) {
					
      vb->val[l + j + t] = v->val[k + j];	
    }
			
  }

  return vb;

}

/**
 * @brief   Create an image of size (nrow+t) times (ncol+t). 
 * 
 * The image is filled with ones in the central region of size nrow x ncol and 
 * has a border of zeros of width t.
 * 
 * @param   nrow        number of rows of the region filled with ones.
 * @param   ncol        number of columns of the region filled with ones.
 * @param   t           width of the border filled with zeros.
 */
Cimage create_mask ( int nrow, int ncol, int t ) {

  int i,j,l;
  int nrow_mask = nrow + 2*t;
  int ncol_mask = ncol + 2*t;
	
  Cimage mask = new_Cimage(nrow_mask, ncol_mask);
	
  for (i = 0; i < nrow; i++) {
	
    l = (i+t)*ncol_mask;
		
    for (j = 0; j < ncol; j++) {
					
      mask->red[l + j + t] = 1;	
      mask->green[l + j + t] = 1;	
      mask->blue[l + j + t] = 1;	
    }
			
  }
	
	
  return mask;

}

/**
 * @brief   Create an image of size (nrow+t) times (ncol+t). 
 * 
 *  Idem create_mask for gray level images.
 * 
 * @see create_mask
 */
Gimage create_mask_gray( int nrow, int ncol, int t ) {

  int i,j,l;
  int nrow_mask = nrow + 2*t;
  int ncol_mask = ncol + 2*t;
	
  Gimage mask = new_Gimage(nrow_mask, ncol_mask);
	
  for (i = 0; i < nrow; i++) {
	
    l = (i+t)*ncol_mask;
		
    for (j = 0; j < ncol; j++) 					
      mask->val[l + j + t] = 1;						
  }	
	
  return mask;

}

/**
 * @brief   Check if an image is a real color image or a 3-channel gray image.
 * 
 * @param   image       input image to be tested.
 * @param   num_pixels  total number of pixels in the image (nrows * ncols)
 */
int is_grayscale(Cimage image, long num_pixels)
{
  long i;
	
  for( i = 0; i < num_pixels; i++)
    if(image->red[i] != image->green[i] || image->red[i] != image->blue[i])
      return 0;    /* Not a grayscale image */

  return 1;    /* This is a grayscale image */
}


/**
 * @brief   Load a color image from a PNG file.
 * 
 * @param   filename    name (path) of the image
 * 
 * @return  The loaded image
 */
Cimage load_png_image(char *filename)
{
  Cimage image;
  size_t nx, ny, nc;
  unsigned char *img = NULL;
  
  /* read the image */
  /* nx = number of rows */
  /* ny = number of columns */
  img = io_png_read_u8( filename, &ny, &nx, &nc);

  /* if error while reading */
  if (NULL == img) {
    fprintf(stderr, "Failed to read the image.\n");
    abort();
  }

  /* nx, ny and nc hols the image sizes */
  image = new_Cimage((int) nx, (int) ny);
  

  if(nc >= 3) {
    /* copy each channel */
    memcpy(image->red, img, nx * ny * sizeof(unsigned char));
    memcpy(image->green, img + nx * ny, nx * ny * sizeof(unsigned char));
    memcpy(image->blue, img + 2 * nx * ny, nx * ny * sizeof(unsigned char));
	
    image->nchannels = 3;
	
  } else {
    memcpy(image->red, img, nx * ny * sizeof(unsigned char));
    memcpy(image->green, img, nx * ny * sizeof(unsigned char));
    memcpy(image->blue, img, nx * ny * sizeof(unsigned char));
	
    image->nchannels = 1;
  }

  return(image);

}







  
