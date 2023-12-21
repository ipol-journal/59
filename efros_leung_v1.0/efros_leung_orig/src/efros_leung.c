#include "efros_leung.h"
#include <stdio.h>
#include <time.h>

#ifndef MIN_DICT_SIZE
#define MIN_DICT_SIZE   5
#endif 

#ifndef ALERT_DICT_SIZE
#define ALERT_DICT_SIZE  20 
#endif

/**
 * @file    efros_leung.h
 * @brief   Efros-Leung texture synthesis.
 * 
 * This file implements Efros & Leung's algorithm for texture synthesis for color images.
 * It is described in
 *      "Texture synthesis by non-parametric sampling"
 *      by A. A. Efros and T. K. Leung.
 *      in International Conference on Computer Vision, pages 1033–1038,
 *      September 1999.

 * @author  Cecilia Aguerrebere
 */

/**
 * @brief   Implements Efros-Leung texture synthesis algorithm for color images.
 *
 * @param t
 *      Half-size of the patches.
 * @param out_img_sz
 *      Size of the synthesized image.
 * @param tol
 *      Tolerance parameter \f$\varepsilon\f$.
 * @param init
 *      If non-zero, use a \em random patch of the image
 *      to initialize the algorithm.
 * @param w
 *      Sample image.
 * @param map_out
 *      Synthesis map of the sample image.
 * @param copy_map_out
 *      Synthesis map of the output image.
 * 
 * @todo
 *      File \em efros_leung.c, function \em efros_leung_synth:
 *      arguments \em map_out et \em copy_map_out ?    
 * 
 * @return  The synthesized image.
 */


Cimage efros_leung_synth(int t, int out_img_sz, float tol, int init, Cimage w, Cimage *map_out, Cimage *copy_map_out)
{
    
  /* -----------------------------------------------------------------*/
  /* Variable setting                                                 */
  /* -----------------------------------------------------------------*/
  int i, j;
  unsigned long int address_v, address_w;
  const int nrow_w = w->nrow;
  const int ncol_w = w->ncol;
  const int nrow_v = ((int)(out_img_sz/2))*2 + 1;  
  const int ncol_v = nrow_v;  
  const int half_ncol_v = (int) ((ncol_v - 1) / 2);
  const int half_nrow_v = (int) ((nrow_v - 1) / 2);
  const int half_ncol_init = 1;
  const int patchSize = 2 * t + 1;
  const int nrow_dict = nrow_w;
  const int ncol_dict = ncol_w;
  const long int size_dict = nrow_dict * ncol_dict;
  const float tolerance = (1.0 + tol);

  /* Initialize random generator */
  mt_init((unsigned long int) time(NULL));

  /* Verify that the input image is large enough to select the initial 3x3 seed */
  if (nrow_w < (2 * half_ncol_init + 1) || ncol_w < (2 * half_ncol_init + 1)) {
    error ("The dictionary is too too small.");
  }
  
  /* Verify that the patch size is not too big */ 
  if ( nrow_dict < 0 || ncol_dict < 0 ){
    error ("The patch size is too large.");
  }
  
  /* Verify that the input image is large enough to have a reasonable number of patch candidates */
  if ( size_dict < MIN_DICT_SIZE ) {	 
    error("The dictionary is too too small.");
  }  
    
  if ( size_dict < ALERT_DICT_SIZE ) {	
    printf("WARNING: The dictionary has only %li patches.\n",size_dict);
    printf("\n");
  }   
      
  printf("\n");
  printf("-------------------------------------------\n");
  printf("Pre-steps to textures synthesis          -- \n");
  printf("-------------------------------------------\n");
  printf("\n");
  printf("Dictionary size: %li patches.\n",size_dict);
  printf("\n");

  /* -----------------------------------------------------------------*/
  /* Pre-compute weights                                              */
  /* -----------------------------------------------------------------*/
  float *weights = NULL;  
  weights = gaussian_weights(patchSize, t);   

  /* -----------------------------------------------------------------*/
  /* Initialize output image                                          */
  /* -----------------------------------------------------------------*/
  Cimage v = new_Cimage(nrow_v, ncol_v);   /* Output color image */
  
  /* -----------------------------------------------------------------*/
  /* Allocate memory                                                  */
  /* -----------------------------------------------------------------*/
  unsigned char  * mask = (unsigned char *) calloc (ncol_v * nrow_v, sizeof (unsigned char));
  unsigned short * mask_neighb = (unsigned short *) calloc(patchSize*patchSize,sizeof(unsigned short));
    
  if (!mask || !mask_neighb )
    error ("Not enough memory");

  unsigned char * curr_patch = (unsigned char *) malloc(3 * patchSize * patchSize * sizeof(unsigned char));
  Pixel curr_px_list = (struct pixel *) malloc ((ncol_v * nrow_v) * sizeof (struct pixel));
  Cand_dist cand_list = (Cand_dist) malloc(size_dict*(sizeof(struct cand_dist)));  

  if(!curr_patch || !curr_px_list || !cand_list)
    error ("Not enough memory");
     
  /* -----------------------------------------------------------------*/
  /* Load patch list                                                  */
  /* -----------------------------------------------------------------*/
  Cimage wb = add_border( w, t );
  Cimage mask_patchs = create_mask ( w->nrow, w->ncol, t );
  
  unsigned char *patch_list = load_patch_list( wb, t );
  unsigned char *patch_list_mask = load_patch_list( mask_patchs, t );

  /* -----------------------------------------------------------------*/
  /* Create the copy/paste map to control copy level                  */
  /* -----------------------------------------------------------------*/
  Cimage map = create_image_map( nrow_w, ncol_w );
  Cimage copy_map = new_Cimage( nrow_v, ncol_v );
 
  /* -----------------------------------------------------------------*/
  /* Synthesized image initialization                                 */
  /* -----------------------------------------------------------------*/
  int initx = (int) (nrow_w / 2);
  int inity = (int) (ncol_w / 2);
  
  if (init != 0) {
    /* Get a random patch */
    initx = half_ncol_init  + (unsigned int) (mt_drand53() * (nrow_w - 2 * half_ncol_init));
    inity = half_ncol_init  + (unsigned int) (mt_drand53() * (ncol_w - 2 * half_ncol_init));
  }

  printf("\n");
  printf ("Seed (%d x %d) position in the sample image: (%d,%d)\n", 2 * half_ncol_init + 1, 2 * half_ncol_init + 1, initx, inity);

  for (i = -half_ncol_init; i < half_ncol_init + 1; i++)
    for (j = -half_ncol_init; j < half_ncol_init + 1; j++)
      {
        address_v = j + half_ncol_v + (i + half_nrow_v) * ncol_v;
        address_w = j + inity + (i + initx) * ncol_w;

        v->red[address_v] =  w->red[address_w];
        v->green[address_v] = w->green[address_w];
        v->blue[address_v] = w->blue[address_w];

        copy_map->red[address_v] =  map->red[address_w];
        copy_map->green[address_v] = map->green[address_w];
        copy_map->blue[address_v] = map->blue[address_w];

        mask[address_v] = 1;
      }

  /* -----------------------------------------------------------------*/
  /* START TEXTURE SYNTHESIS                                          */
  /* -----------------------------------------------------------------*/
  long int curr_list_sz;
  unsigned int neighb_updated;
  int total_cands;
  long int chosen;  
  int chosen_x, chosen_y;
  
  int corner = half_ncol_v - half_ncol_init - 1;
  int side = (2 * half_ncol_init + 1) + 2;
  unsigned long int filled_pixels = (2 * half_ncol_init + 1) * (2 * half_ncol_init + 1);

  printf("\n");
  printf("-------------------------------------------\n");
  printf("Starting texture synthesis             -- \n");
  printf("-------------------------------------------\n");
  
  /* Main loop */
  while (corner >= 0)
    {
      /* Generate the list of pixels to be filled. Ordered by number of neighbours (decreasing order). */
      curr_list_sz = generate_current_list( side, corner, curr_px_list, mask, t, ncol_v, nrow_v );

      /* For all the pixels in the list */
      for (i = 0; i < curr_list_sz; i++ ) {

		/* Find current filled neighbours */
		neighb_updated = known_neighbours( v, mask, curr_px_list, i, t, ncol_v, nrow_v, mask_neighb, curr_patch );      
	      
		/* Find candidates list */
		total_cands = find_candidates( patch_list, patch_list_mask, curr_patch, size_dict, t,  mask_neighb, neighb_updated, tolerance, cand_list, weights );
	
		/* Randomly draw a candidate number from the candidates list */
		chosen = random_choose( cand_list, total_cands );
	
		/* Retreive the coordinates on the sample image of the chosen candidate */
		retrieve_coords( chosen, ncol_w , &chosen_x, &chosen_y);
	
		/* From the previous coordinates (chosen_x, chosen_y)
		 * find the addresses in the sample image vector (address_w)
		 * and in the synthesized image vector (address_v) */
		address_v = curr_px_list[i].y + curr_px_list[i].x * ncol_v;
		address_w = chosen_y + chosen_x * ncol_w;
	
		/* Fill the synthesized image v with the chosen pixel value from the sample image w */
		v->red[address_v] = w->red[address_w];
		v->green[address_v] = w->green[address_w];
		v->blue[address_v] = w->blue[address_w];
	
		copy_map->red[address_v] =  map->red[address_w];
		copy_map->green[address_v] = map->green[address_w];
		copy_map->blue[address_v] = map->blue[address_w];
	
		/* Update the filled pixel in the mask */
		mask[address_v] = 1;
	      
      }

      /* ---------- Update the region to be filled ---------- */
      side = side + 2;
      corner = corner - 1;
      filled_pixels += curr_list_sz;
      printf ("\rSynthesis... %2d %%", (int) (100 * (float) filled_pixels) / (ncol_v * nrow_v));
      fflush (stdout);
    }

  *copy_map_out = copy_map;
  *map_out = map;
  
  /* -----------------------------------------------------------------*/
  /* THE WORK IS OVER :)                                              */
  /* -----------------------------------------------------------------*/
  free(mask);
  free(mask_neighb);
  free(curr_patch);
  free(curr_px_list);
  free(cand_list);
  free(patch_list);
  free(patch_list_mask);
  free(wb);
  free(mask_patchs);


  /* Exit */
  printf ("\rSynthesis... OK   \n");    
  
  /* The size of the synthesized image had been forced to be odd in order to make the filling by layers.  */
  /* Now, if the output image size was not odd, the image is cropped to fit the requested output size.    */
  if ( nrow_v != out_img_sz ) {
    Cimage v_crop = crop_image( v, out_img_sz );     
    return v_crop;
  } else
  
    return v;
}


/**
 * @brief   Implements Efros-Leung texture synthesis algorithm for gray level images.
 *
 * @param t
 *      Half-size of the patches.
 * @param out_img_sz
 *      Size of the synthesized image.
 * @param tol
 *      Tolerance parameter \f$\varepsilon\f$.
 * @param init
 *      If non-zero, use a \em random patch of the image
 *      to initialize the algorithm.
 * @param w
 *      Sample image.
 * @param map_out
 *      Synthesis map of the sample image.
 * @param copy_map_out
 *      Synthesis map of the output image.
 * 
 * @todo
 *      File \em efros_leung.c, function \em efros_leung_synth:
 *      arguments \em map_out et \em copy_map_out ?    
 * 
 * @return  The synthesized image.
 */

Gimage efros_leung_synth_gray(int t, int out_img_sz, float tol, int init, Gimage w, Cimage *map_out, Cimage *copy_map_out)
{
  
  /* -----------------------------------------------------------------*/
  /* Variable setting                                                 */
  /* -----------------------------------------------------------------*/
  int i, j;
  unsigned long int address_v, address_w;
  const int nrow_w = w->nrow;
  const int ncol_w = w->ncol;
  const int nrow_v = ((int)(out_img_sz/2))*2 + 1;  
  const int ncol_v = nrow_v;  
  const int half_ncol_v = (int) ((ncol_v - 1) / 2);
  const int half_nrow_v = (int) ((nrow_v - 1) / 2);
  const int half_ncol_init = 1;
  const int patchSize = 2 * t + 1;
  const int nrow_dict = nrow_w;
  const int ncol_dict = ncol_w;
  const long int size_dict = nrow_dict * ncol_dict;
  const float tolerance = (1.0 + tol);

  /* Initialize random generator */
  mt_init((unsigned long int) time(NULL));  
  
  /* Verify that the input image is large enough to select the initial 3x3 seed */
  if (nrow_w < (2 * half_ncol_init + 1) || ncol_w < (2 * half_ncol_init + 1)) {
    error ("The dictionary is too too small.");
  }
  
  /* Verify that the patch size is not too big */ 
  if ( nrow_dict < 0 || ncol_dict < 0 ){
	 error ("The patch size is too large.");
  }
  
  /* Verify that the input image is large enough to have a reasonable number of patch candidates */
  if ( size_dict < MIN_DICT_SIZE ) {	 
    error("The dictionary is too too small.");
  }  
    
  if ( size_dict < ALERT_DICT_SIZE ) {	
    printf("WARNING: The dictionary has only %li patches.\n",size_dict);
    printf("\n");
  }   
  
  printf("\n");
  printf("-------------------------------------------\n");
  printf("Pre-steps to textures synthesis          -- \n");
  printf("-------------------------------------------\n");
  printf("\n");
  printf("Dictionary size: %li patches.\n",size_dict);
  printf("\n");

  /* -----------------------------------------------------------------*/
  /* Pre-compute weights                                              */
  /* -----------------------------------------------------------------*/
  float *weights = NULL;  
  weights = gaussian_weights(patchSize, t);   

  /* -----------------------------------------------------------------*/
  /* Initialize output image                                          */
  /* -----------------------------------------------------------------*/
  Gimage v = new_Gimage(nrow_v, ncol_v);   /* Output color image */
  
  /* -----------------------------------------------------------------*/
  /* Allocate memory                                                  */
  /* -----------------------------------------------------------------*/
  unsigned char  * mask = (unsigned char *) calloc (ncol_v * nrow_v, sizeof (unsigned char));
  unsigned short * mask_neighb = (unsigned short *) calloc(patchSize*patchSize,sizeof(unsigned short));
    
  if (!mask || !mask_neighb )
    error ("Not enough memory");

  unsigned char * curr_patch = (unsigned char *) malloc(patchSize * patchSize * sizeof(unsigned char));
  Pixel curr_px_list = (struct pixel *) malloc ((ncol_v * nrow_v) * sizeof (struct pixel));
  Cand_dist cand_list = (Cand_dist) malloc(size_dict*(sizeof(struct cand_dist)));  

  if(!curr_patch || !curr_px_list || !cand_list)
    error ("Not enough memory");
     
  /* -----------------------------------------------------------------*/
  /* Load patch list                                                  */
  /* -----------------------------------------------------------------*/
  Gimage wb = add_border_gray( w, t );
  Gimage mask_patchs = create_mask_gray( w->nrow, w->ncol, t );
  
  unsigned char *patch_list = load_patch_list_gray( wb, t );
  unsigned char *patch_list_mask = load_patch_list_gray( mask_patchs, t );

  /* -----------------------------------------------------------------*/
  /* Create the copy/paste map to control copy level                  */
  /* -----------------------------------------------------------------*/
  Cimage map = create_image_map( nrow_w, ncol_w );
  Cimage copy_map = new_Cimage( nrow_v, ncol_v );
 
  /* -----------------------------------------------------------------*/
  /* Synthesized image initialization                                 */
  /* -----------------------------------------------------------------*/
  int initx = (int) (nrow_w / 2);
  int inity = (int) (ncol_w / 2);
  
  if (init != 0) {
    /* Get a random patch */
    initx = half_ncol_init  + (unsigned int) (mt_drand53() * (nrow_w - 2 * half_ncol_init));
    inity = half_ncol_init  + (unsigned int) (mt_drand53() * (ncol_w - 2 * half_ncol_init));
  }

  printf("\n");
  printf ("Seed (%d x %d) position in the sample image: (%d,%d)\n", 2 * half_ncol_init + 1, 2 * half_ncol_init + 1, initx, inity);

  for (i = -half_ncol_init; i < half_ncol_init + 1; i++)
    for (j = -half_ncol_init; j < half_ncol_init + 1; j++)
      {
        address_v = j + half_ncol_v + (i + half_nrow_v) * ncol_v;
        address_w = j + inity + (i + initx) * ncol_w;

        v->val[address_v] =  w->val[address_w];

        copy_map->red[address_v] =  map->red[address_w];
        copy_map->green[address_v] = map->green[address_w];
        copy_map->blue[address_v] = map->blue[address_w];

        mask[address_v] = 1;
      }

  /* -----------------------------------------------------------------*/
  /* START TEXTURE SYNTHESIS                                          */
  /* -----------------------------------------------------------------*/
  long int curr_list_sz;
  unsigned int neighb_updated;
  int total_cands;
  long int chosen;  
  int chosen_x, chosen_y;
  
  int corner = half_ncol_v - half_ncol_init - 1;
  int side = (2 * half_ncol_init + 1) + 2;
  unsigned long int filled_pixels = (2 * half_ncol_init + 1) * (2 * half_ncol_init + 1);

  printf("\n");
  printf("-------------------------------------------\n");
  printf("Starting textures synthesis             -- \n");
  printf("-------------------------------------------\n");
  
  /* Main loop */
  while (corner >= 0)
  {
    /* Generate the list of pixels to be filled. Ordered by number of neighbours (decreasing order). */
    curr_list_sz = generate_current_list( side, corner, curr_px_list, mask, t, ncol_v, nrow_v );

    /* For all the pixels in the list */
    for (i = 0; i < curr_list_sz; i++ ) {

      /* Find current filled neighbours */
      neighb_updated = known_neighbours_gray( v, mask, curr_px_list, i, t, ncol_v, nrow_v, mask_neighb, curr_patch );      
      
      /* Find candidates list */
      total_cands = find_candidates_gray( patch_list, patch_list_mask, curr_patch, size_dict, t,  mask_neighb, neighb_updated, tolerance, cand_list, weights );

      /* Randomly draw a candidate number from the candidates list */
      chosen = random_choose( cand_list, total_cands );

      /* Retreive the coordinates on the sample image of the chosen candidate */
      retrieve_coords( chosen, ncol_w , &chosen_x, &chosen_y);

      /* From the previous coordinates (chosen_x, chosen_y)
       * find the addresses in the sample image vector (address_w)
       * and in the synthesized image vector (address_v) */
      address_v = curr_px_list[i].y + curr_px_list[i].x * ncol_v;
      address_w = chosen_y + chosen_x * ncol_w;

      /* Fill the synthesized image v with the chosen pixel value from the sample image w */
      v->val[address_v] = w->val[address_w];

      copy_map->red[address_v] =  map->red[address_w];
      copy_map->green[address_v] = map->green[address_w];
      copy_map->blue[address_v] = map->blue[address_w];

      /* Update the filled pixel in the mask */
      mask[address_v] = 1;
      
    }

    /* ---------- Update the region to be filled ---------- */
    side = side + 2;
    corner = corner - 1;
    filled_pixels += curr_list_sz;
    printf ("\rSynthesis... %2d %%", (int) (100 * (float) filled_pixels) / (ncol_v * nrow_v));
    fflush (stdout);
  }

  *copy_map_out = copy_map;
  *map_out = map;
  
  /* -----------------------------------------------------------------*/
  /* THE WORK IS OVER :)                                              */
  /* -----------------------------------------------------------------*/
  free(mask);
  free(mask_neighb);
  free(curr_patch);
  free(curr_px_list);
  free(cand_list);
  free(patch_list);
  free(patch_list_mask);
  free(wb);
  free(mask_patchs);

  /* Exit */
  printf ("\rSynthesis... OK   \n");
  
  /* The size of the synthesized image had been forced to be odd in order to make the filling by layers.  */
  /* Now, if the output image size was not odd, the image is cropped to fit the requested output size.    */
  if ( nrow_v != out_img_sz ) {
    Gimage v_crop = crop_image_gray( v, out_img_sz );     
	return v_crop;
  } else
  
  return v;
}

	  
	  

	  
	  
