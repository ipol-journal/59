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
 * This file implements an accelerated version of the Efros & Leung's algorithm for texture synthesis.
 * It is described in
 *      "Exemplar-based texture synthesis : the Efros-Leung algorithm"
 *      by C. Aguerrebere, Y. Gousseau and G. Tartavel.
 *      in Image Processing Online (IPOL)
 *      December 2012.

 * @author  Cecilia Aguerrebere
 */

/**
 * @brief   Implements Efros-Leung texture synthesis algorithm.
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
 * @param dims_pca
 *      Number of PCA components used to compute distances
 *      between (the known parts of) patches.
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
Cimage efros_leung_synth(int t, int out_img_sz, float tol, int init, int dims_pca, Cimage w, Cimage *map_out, Cimage *copy_map_out)
{

  /* -----------------------------------------------------------------*/
  /* Variable setting                                                 */
  /* -----------------------------------------------------------------*/
  int i, j;
  unsigned long int address_v, address_w;
  const int ncol_w = w->ncol;
  const int nrow_w = w->nrow;  
  const int nrow_v = ((int)(out_img_sz/2))*2 + 1;
  const int ncol_v = nrow_v;
  const int half_ncol_v = (int) ((ncol_v - 1) / 2);
  const int half_nrow_v = (int) ((nrow_v - 1) / 2);
  const int half_ncol_init = 1;
  
  const int patchSize = 2 * t + 1;
  const int nrow_dict = nrow_w;
  const int ncol_dict = ncol_w;
  const long int size_dict = nrow_dict * ncol_dict;
  const long int size_dictUD = (ncol_w - 2*t)*(nrow_w - t);
  const long int size_dictRL = (ncol_w - t)*(nrow_w - 2*t);
  
  const int orig_dims = patchSize*t;
  const int patchSize_plus2 = patchSize + 2;  
  const float tolerance = 1 + tol;
  

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

  /* Verify that the dims_pca parameter is positive and smaller than the maximum number of PCA components */
  if ( dims_pca < 1 || dims_pca > patchSize*t ) {	 
    error("Incorrect number of PCA components.");
  }  
  

  printf("\n");
  printf("-------------------------------------------\n");
  printf("- Pre-steps to textures synthesis       -- \n");
  printf("-------------------------------------------\n");
  printf("\n");
  printf("Dictionary size: %li patches.\n",size_dict);
  printf("\n");
  
  /* -----------------------------------------------------------------*/
  /* Pre-compute weights                                              */
  /* -----------------------------------------------------------------*/
  float *weights = NULL;
  float *weights_pca = NULL;
  
  weights = gaussian_weights(patchSize, t);
  weights_pca = gaussian_weights_pca(patchSize, t);    

  /* -----------------------------------------------------------------*/
  /* RGB PCA: Reduce 3 dimensions (RGB) to an unique dimension.       */
  /* In the case of monochromatic images, convert CImage to ImageF.   */
  /* -----------------------------------------------------------------*/
  ImageF ww;
  
  if ( w->nchannels == 3 ) 
    ww = pca_rgb( w );
  else 
    ww = CImage_to_ImageF( w );

  /* -----------------------------------------------------------------*/
  /* Create combined center_patch - PCA dictionary                    */
  /* -----------------------------------------------------------------*/
  float *YU, *YD, *YR, *YL;
  float *WU, *WD, *WR, *WL;
  float *dictionaryU, *dictionaryD, *dictionaryR, *dictionaryL;
  float *meanU, *meanD, *meanR, *meanL;

  printf("Computing PCA base...");
  fflush(stdout);
  pca( ww, t, dims_pca, weights_pca, &WU, &WD, &WR, &WL, &YU, &YD, &YR, &YL, &meanU, &meanD, &meanR, &meanL );
  printf(" OK\n");

  printf("Building dictionary...");
  fflush(stdout);
  construct_dictionary_pca( ww, weights_pca, t, dims_pca, YU, YD, YR, YL, &dictionaryR, &dictionaryL, &dictionaryU, &dictionaryD );
  printf(" OK\n");

  free(YU);
  free(YD);
  free(YL);
  free(YR);

  /* -----------------------------------------------------------------*/
  /* Initialize output image                                          */
  /* -----------------------------------------------------------------*/
  ImageF vv = new_ImageF(nrow_v, ncol_v);  /* Float image for synthesis (1 component) */
  Cimage v = new_Cimage(nrow_v, ncol_v);   /* Output color image */

  /* -----------------------------------------------------------------*/
  /* Allocate memory for masks and patch list                         */
  /* -----------------------------------------------------------------*/
  unsigned char  * mask = (unsigned char *) calloc (ncol_v * nrow_v, sizeof (unsigned char));
  unsigned short * mask_neighb = (unsigned short *) calloc(patchSize,sizeof(unsigned short));
  unsigned short * mask_neighb_compl = (unsigned short *) calloc(patchSize*patchSize,sizeof(unsigned short));
  
  if (!mask || !mask_neighb || !mask_neighb_compl)
    error ("Not enough memory");

  float * curr_patch = (float *) malloc(patchSize * patchSize * sizeof(float));
  Pixel curr_px_list = (struct pixel *) malloc ((ncol_v * nrow_v) * sizeof (struct pixel));
  Cand_dist cand_list = (Cand_dist) malloc(size_dict*(sizeof(struct cand_dist)));  

  if(!curr_patch || !curr_px_list || !cand_list)
    error ("Not enough memory");

  /* -----------------------------------------------------------------*/
  /* Load patch list                                                  */
  /* -----------------------------------------------------------------*/
  ImageF wb = add_border( ww, t );
  ImageF mask_patchs = create_mask ( ww->nrow, ww->ncol, t );
  
  float *patch_list = load_patch_list( wb, t );
  float *patch_list_mask = load_patch_list( mask_patchs, t );


  /* -----------------------------------------------------------------*/
  /* Create the copy/paste map to control copy level                  */
  /* -----------------------------------------------------------------*/
  Cimage map = create_image_map( nrow_w, ncol_w );
  Cimage copy_map = new_Cimage( nrow_v, ncol_v );
 
  /* -----------------------------------------------------------------*/
  /* Synthesized image initialization from 3x3 seed                   */
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

        vv->val[address_v] = ww->val[address_w];
        mask[address_v] = 1;
      }

  /* -----------------------------------------------------------------*/
  /* START TEXTURE SYNTHESIS                                          */
  /* -----------------------------------------------------------------*/
  long int curr_list_sz;
  unsigned int neighb_updated;
  int total_cands;
  long int chosen;
  int pixels_pca, total_neighb;
  int chosen_x, chosen_y;
  
  int corner = half_ncol_v - half_ncol_init - 1;
  int side = (2 * half_ncol_init + 1) + 2;
  unsigned long int filled_pixels = (2 * half_ncol_init + 1) * (2 * half_ncol_init + 1);

  printf("\n");
  printf("-------------------------------------------\n");
  printf("- Starting texture synthesis             -- \n");
  printf("-------------------------------------------\n");
  

  /* Main loop */
  while (corner >= 0)
    {
      /* Generate the list of pixels to be filled. Ordered by number of neighbours (decreasing order). */
      curr_list_sz = generate_current_list( side, corner, curr_px_list, mask, t, ncol_v, nrow_v );

      /* ---------- Patches with PCA ---------- */
      if (side < patchSize_plus2 )
	pixels_pca = 0;
      else
	{
	  pixels_pca = 4*(side - 2 -2*t);
      
	  for (i = 0; i < pixels_pca; i++) {

	    /* Do the same as pixels outside the PCA part (below)
	     * but in the right dictionary (U: up, D: down, R: right, L: left) */
        
	    if ( curr_px_list[i].y == corner ) {
	      /* dictionary R */
	      total_neighb = known_neighboursRL( mask, curr_px_list, i, t, ncol_v, mask_neighb );
	      total_cands = find_candidatesR( vv, weights_pca, dictionaryR, curr_px_list, i, size_dictRL, mask_neighb, total_neighb, t, WR, meanR, orig_dims, dims_pca, tolerance, cand_list);
	      chosen = random_choose( cand_list, total_cands );          
	      retrieve_coords( chosen, ncol_w - t, t, 0, &chosen_x, &chosen_y);
          
	    } else if ( curr_px_list[i].y == corner + side - 1 ) {
	      /* dictionary L*/
	      total_neighb = known_neighboursRL( mask, curr_px_list, i, t, ncol_v, mask_neighb );
	      total_cands = find_candidatesL( vv, weights_pca, dictionaryL, curr_px_list, i, size_dictRL, mask_neighb, total_neighb, t, WL, meanL, orig_dims, dims_pca, tolerance, cand_list);
	      chosen = random_choose( cand_list, total_cands );          
	      retrieve_coords( chosen, ncol_w - t, t, t, &chosen_x, &chosen_y);
          
	    } else if ( curr_px_list[i].x == corner ) {
	      /* dictionary D */
	      total_neighb = known_neighboursUD( mask, curr_px_list, i, t, ncol_v, mask_neighb );
	      total_cands = find_candidatesD( vv, weights_pca, dictionaryD, curr_px_list, i, size_dictUD, mask_neighb, total_neighb, t, WD, meanD, orig_dims, dims_pca, tolerance, cand_list);
	      chosen = random_choose( cand_list, total_cands );          
	      retrieve_coords( chosen, ncol_w - 2*t, 0, t, &chosen_x, &chosen_y);
          
	    } else {
	      /* dictionary U */
	      total_neighb = known_neighboursUD( mask, curr_px_list, i, t, ncol_v, mask_neighb);
	      total_cands = find_candidatesU( vv, weights_pca, dictionaryU, curr_px_list, i, size_dictUD, mask_neighb, total_neighb, t, WU, meanU, orig_dims, dims_pca, tolerance, cand_list);
	      chosen = random_choose( cand_list, total_cands );          
	      retrieve_coords( chosen, ncol_w - 2*t, t, t, &chosen_x, &chosen_y);          
	    }

	    /* Find addres of the pixel to be filled (address_v) and the chosen pixel (address_w) */
	    address_v = curr_px_list[i].y + curr_px_list[i].x * ncol_v;
	    address_w = chosen_y + chosen_x * ncol_w;

	    /* Fill the pixel with the chosen candidate */
	    v->red[address_v] = w->red[address_w];
	    v->green[address_v] = w->green[address_w];
	    v->blue[address_v] = w->blue[address_w];
        
	    /* Update the synthesis map */
	    copy_map->red[address_v] =  map->red[address_w];
	    copy_map->green[address_v] = map->green[address_w];
	    copy_map->blue[address_v] = map->blue[address_w];
		
	    /* Fill the 1-component image vv */
	    vv->val[address_v] = ww->val[address_w];
        
	    /* Update the mask of already filled pixels*/
	    mask[address_v] = 1;

	  }
      
	}

      /* ---------- Patches without PCA ---------- */
      for (i = pixels_pca; i < curr_list_sz; i++ ) {
      
	/* Find current filled neighbours */
	neighb_updated = known_neighbours( vv, mask, curr_px_list, i, t, ncol_v, nrow_v, mask_neighb_compl, curr_patch );
      
	/* Find candidates list */
	total_cands = find_candidates( patch_list,patch_list_mask, curr_patch, size_dict, t,  mask_neighb_compl, neighb_updated, tolerance, cand_list, weights );

	/* Randomly draw a candidate number from the candidates list */
	chosen = random_choose( cand_list, total_cands );

	/* Retreive the coordinates on the sample image of the chosen candidate */
	retrieve_coords( chosen, ncol_w , 0, 0, &chosen_x, &chosen_y);

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

	/* Fill the 1-component image vv */
	vv->val[address_v] = ww->val[address_w];

	/* Update the mask of already filled pixels*/
	mask[address_v] = 1;
      
      }

      /* ---------- Update the region to be filled ---------- */
      side = side + 2;
      corner = corner - 1;
      filled_pixels += curr_list_sz;

    }
  
  /* -----------------------------------------------------------------*/
  /* THE WORK IS OVER :)                                              */
  /* -----------------------------------------------------------------*/
  free(mask);
  free(mask_neighb);
  free(mask_neighb_compl);
  free(curr_patch);
  free(curr_px_list);
  free(cand_list);
  free(patch_list);

  /* Exit */
  printf ("\rSynthesis... OK   \n");
 
  *copy_map_out = copy_map;
  *map_out = map;
  
  /* The size of the synthesized image had been forced to be odd in order to make the filling by layers.  */
  /* Now, if the output image size was not odd, the image is cropped to fit the requested output size.    */
  if ( nrow_v != out_img_sz ) {
	  
    Cimage v_crop = crop_image( v, out_img_sz );     
 	
    return v_crop;
	
  } else 
    return v;
  

}
