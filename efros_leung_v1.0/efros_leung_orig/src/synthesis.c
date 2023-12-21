#include "synthesis.h"

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/**
 * @file    synthesis.h
 * @brief   Sub-functions for texture synthesis.
 * @author  Cecilia Aguerrebere
*/

#ifndef MAX_DIST
#define MAX_DIST   100000000.
#endif 

#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif /* !M_PI */

#ifndef NEAR_ZERO
#define NEAR_ZERO   .000001
#endif 


/**
 * @brief   Extract the neighbors' mask and values at a given pixel.
 * 
 * @param   v       the image (color image)
 * @param   mask    the mask for the partially synthesized image \em v
 * @param   current list of pixels
 * @param   i       index of the pixel in the \em current list
 * @param   t       size of the half-patch
 * @param   ncols   number of columns (dx) \em v
 * @param   nrows   number of rows (dy) in \em v
 * @param   mask_neighb     output array for the neighbors' mask
 * @param   curr_patch      output array for the neighbors' values
 * 
 * @return  number of known neighbors
 */
int known_neighbours( Cimage v, unsigned char* mask, Pixel current, int i, int t, int ncols_v, int nrows_v, unsigned short *mask_neighb, unsigned char *curr_patch )

{
  int k, l;
  int counter = 0;
  int counter2 = 0;
  unsigned short neighbs = 0;
  
  /* For all the pixels in the patch centered at pixel i (current pixel being filled), find which pixels are already filled. */
  for (k = -t; k < t + 1; k++)   
    for (l = -t; l < t + 1; l++) {
      
      /* Border check: only consider the pixels on the patch that are included in image. */
      if (((current[i].y + l) >= 0) && ((current[i].y + l) < ncols_v) && ((current[i].x + k) >= 0) && ((current[i].x + k) < nrows_v))
        
        /* Check if pixel (k,l) is filled */  
        if ( mask[current[i].y + l + (current[i].x + k) * ncols_v] == 1 ) {
          
          /* If pixel (k,l) is filled, save in mask_neighg the number identifying the known neighbour. */
          mask_neighb[counter] = neighbs;     
          
          /* Save also the current patch for posterior computations */   
	  curr_patch[counter2] = v->red[current[i].y + l + (current[i].x + k) * ncols_v];   /* Aprovecho para guardar el patch actual */
	  curr_patch[counter2 + 1] = v->green[current[i].y + l + (current[i].x + k) * ncols_v];
	  curr_patch[counter2 + 2] = v->blue[current[i].y + l + (current[i].x + k) * ncols_v];         
          
          counter++;
	  counter2 = counter2 + 3;
        }
      neighbs++;
    }
  
  return counter;
  
}

/**
 * @brief   Extract the neighbors' mask and values at a given pixel.
 * 
 *  Idem known_neighbours for gray level images.
 *  
 * @see known_neighbours
 */
 int known_neighbours_gray( Gimage v, unsigned char* mask, Pixel current, int i, int t, int ncols_v, int nrows_v, unsigned short *mask_neighb, unsigned char *curr_patch )

{
  int k, l;
  int counter = 0;
  int counter2 = 0;
  unsigned short neighbs = 0;
  
  /* For all the pixels in the patch centered at pixel i (current pixel being filled), find which pixels are already filled. */
  for (k = -t; k < t + 1; k++)   
    for (l = -t; l < t + 1; l++) {
      
      /* Border check: only consider the pixels on the patch that are included in image. */
      if (((current[i].y + l) >= 0) && ((current[i].y + l) < ncols_v) && ((current[i].x + k) >= 0) && ((current[i].x + k) < nrows_v))
        
        /* Check if pixel (k,l) is filled */  
        if ( mask[current[i].y + l + (current[i].x + k) * ncols_v] == 1 ) {
          
          /* If pixel (k,l) is filled, save in mask_neighg the number identifying the known neighbour. */
          mask_neighb[counter] = neighbs;     
          
          /* Save also the current patch for posterior computations */   
	  curr_patch[counter2] = v->val[current[i].y + l + (current[i].x + k) * ncols_v];   /* Aprovecho para guardar el patch actual */
          
          counter++;
	  counter2++;
        }
      neighbs++;
    }
  
  return counter;
  
}

/**
 * @brief   Randomly choose an element from a list.
 * 
 * The choice is uniform.
 * 
 * @param   cand_list   a list of candidates
 * @param   size_list   the number of candidates
 * 
 * @return  the index of a random candidate
 */
long int random_choose ( Cand_dist cand_list, long int size_list )
{

  /* Randomly choose a candidate. */   
  long int  alea = (long int) (mt_drand53() * size_list);
  
  /* Retrieve chosen patch number */
  long int patch = cand_list[alea].patch;  
  
  return patch;

}


/**
 * @brief   Generate the current list of pixels to be filled.
 * 
 * @param   side
 * @param   corner
 * @param   current     output list of pixels
 * @param   mask        mask of already filled pixels
 * @param   t           half-size of patches
 * @param   ncols       number of columns (dx) for mask
 * @param   nrows       number of rows (dy) for mask
 * 
 * @todo    argument side and corner: still undocumented
 * @return  size of the list
 */
 long int generate_current_list ( int side, int corner, Pixel current, unsigned char *mask, int t, int ncols, int nrows )
{
  int i, j, k, l;
  long int list_sz, counter, counter2, counter3, lenght, alea;
  struct pixel temp;
  int cornerPlusSide = corner + side;
  int cornerPlusSideMinusOne = cornerPlusSide - 1;
  int sideMinusCorner = side - corner;
  int twoSideMinusCornerMinusOne = 2*side - corner - 1;
  int twoSideMinusCornerMinusOne2 = twoSideMinusCornerMinusOne + side - 2;
  
  list_sz = 4 * (side - 1);

  /* Save the coordinates of the pixels in the outer layer */
#pragma omp parallel for private(counter,counter2)
  
  for (i = corner; i < cornerPlusSide; i++)
    {
      counter = i - corner;      
      current[counter].x = i;
      current[counter].y = corner;

      counter2 = i + sideMinusCorner;
      current[counter2].x = i;
      current[counter2].y = cornerPlusSideMinusOne;

    }

#pragma omp parallel for private(counter2,counter3)
  
  for (j = corner + 1; j < cornerPlusSide - 1; j++)
    {
      counter2 = twoSideMinusCornerMinusOne + j; 
      current[counter2].x = cornerPlusSideMinusOne;
      current[counter2].y = j;

      counter3 = twoSideMinusCornerMinusOne2 + j;
      current[counter3].x = corner;
      current[counter3].y = j;

    }
     
      
  /* Count the number of known neighbours for each pixel in the list */
#pragma omp parallel for private(k,l) 
 
  for (i = 0; i < list_sz; i++)
    {
      current[i].nbneighbs = 0;
      for (k = -t; k < t + 1; k++)
        for (l = -t; l < t + 1; l++)
          if ( ((current + i)->y + l) >= 0 && ((current + i)->y + l) < ncols && ((current + i)->x + k) >= 0 && ((current + i)->x + k) < nrows )
            if (mask[(current + i)->y + l + ((current + i)->x + k) * ncols] == 1)
              current[i].nbneighbs++;
    }

  /* Randomly permute the list elements (Knuth shuffle) */
  lenght = list_sz-1; 
      
  while (lenght > 1) 
    { 
      alea = (int) ((lenght+1) * mt_drand53()); 
      lenght -= 1; 
      temp = current[lenght]; 
      current[lenght] = current[alea]; 
      current[alea] = temp; 
    } 

  /* Order the list in decreasing order of known neighbours  */
  qsort (current, list_sz, sizeof (struct pixel), compare_pixel);
  
  return list_sz;

}

/**
 * @brief   Load the patch of an image into an array.
 * 
 * Facilitate posterior access.
 * 
 * @param   w   image  (color image)
 * @param   t   half-size of patches
 * 
 * @return  the patches array extracted from the image
 */
 unsigned char *load_patch_list( Cimage w, int t )
{
  int long addr1, addrU;
  int k, i, j, p;
  int patchSize = (2*t+1)*(2*t+1)*3;
  int patchSizeCh = (2*t+1)*(2*t+1);  
  int T = 2*t + 1;
  int nrows = w->nrow;  
  int ncols = w->ncol;
  int Kcols = ncols - 2*t;
  int Krows = nrows - 2*t;
  int totalPatchs = Kcols * Krows;
  int init = t * ncols + t;

  int index[totalPatchs];
  unsigned char *patch_list = (unsigned char *)malloc(totalPatchs*patchSize*sizeof(unsigned char));

  k = 0;
  for ( i=0; i < Krows ; i++) 
    for ( j=0; j < Kcols ; j++) {   
      index[k] = init + j + i * ncols;
      
      k++;
    }

  /* Fill dictionary of known patches */
  for ( p=0; p < totalPatchs ; p++) {

    addr1 = p*patchSize;
    addrU = index[p] - init;
		
    for ( k=0; k < T; k++ ){      
      memcpy( patch_list + addr1 + k*T,  w->red + addrU + k*ncols, T );           
    } 
    for ( k=0; k < T; k++ ){      
      memcpy( patch_list + addr1 + k*T + patchSizeCh,  w->green + addrU + k*ncols, T );           
    } 
    for ( k=0; k < T; k++ ){      
      memcpy( patch_list + addr1 + k*T + patchSizeCh*2,  w->blue + addrU + k*ncols, T );           
    } 

  }
  
  return patch_list;

}

/**
 * @brief   Load the patch of an image into an array.
 * 
 * Idem load_patch_list for gray level images.
 * 
 * @see load_patch_list
 */
 unsigned char *load_patch_list_gray( Gimage w, int t )
{
  int long addr1, addrU;
  int k, i, j, p;
  int patchSize = (2*t+1)*(2*t+1);
  int T = 2*t + 1;
  int nrows = w->nrow;  
  int ncols = w->ncol;
  int Kcols = ncols - 2*t;
  int Krows = nrows - 2*t;
  int totalPatchs = Kcols * Krows;
  int init = t * ncols + t;

  int index[totalPatchs];
  unsigned char *patch_list = (unsigned char *)malloc(totalPatchs*patchSize*sizeof(unsigned char));

  k = 0;
  for ( i=0; i < Krows ; i++) 
    for ( j=0; j < Kcols ; j++) {   
      index[k] = init + j + i * ncols;
      
      k++;
    }

  /* Fill dictionary of known patches */
  for ( p=0; p < totalPatchs ; p++) {

    addr1 = p*patchSize;
    addrU = index[p] - init;
		
    for ( k=0; k < T; k++ )
      memcpy( patch_list + addr1 + k*T,  w->val + addrU + k*ncols, T );           	 
  }
  
  return patch_list;

}

/**
 * @brief   Get the set of candidates for a given pixel to be filled.
 * 
 * Given the dictionary and the current pixel to be filled, find all the possible  
 * candidates, i.e. all the patches p verifying:
 * 
 *  dist( p, current ) <= dist_min * (1 + tolerance)   
 * 
 * with dist_min = minimum distance between the current patch and the dictionary elements.
 * 
 * The list is given in cand_list, and the number of elements on the list is given
 * by size_cand_list. The cand_list is sorted in increasing order of distance value.
 *
 * @param   patch_list      list of patches extracted from the input image
 * @param   patch_list_mask list of mask of patches extracted from the mask of the input image
 * @param   curr_patch      current patch
 * @param   totalPatchs     number of patches in the patch_list
 * @param   t               half-size of the patch
 * @param   mask_neighb     mask of the known pixels in the current patch
 * @param   neighbs_updated number of known neighbors
 * @param   tolerance       value of \f$(1+\varepsilon)^2\f$
 * @param   cand_list       utput for the candidates (sorted from closest) 
 * @param   weights         weights used for the distance between patches
 *
 * @return  the number of candidates
 */
int find_candidates( unsigned char *patch_list,unsigned char *patch_list_mask, unsigned char *curr_patch, int totalPatchs, int t, unsigned short *mask_neighb, int neighbs_updated, float tolerance, Cand_dist cand_list, float *weights ) 
{

  float radius;
  
  int h = compute_distances( totalPatchs, t, neighbs_updated, mask_neighb, patch_list,patch_list_mask, curr_patch, tolerance, cand_list, &radius, weights );
 
  qsort (cand_list, (size_t) h, (size_t) sizeof (struct cand_dist),*compare_cand_dist);
  
  int size_cand_list = 1;
  while ( size_cand_list < h && cand_list[size_cand_list].dist <= radius ) {
   
    size_cand_list++;
  }
 
  return size_cand_list;

}

/**
 * @brief   Get the set of candidates for a given pixel to be filled.
 * 
 *  Idem find_candidates for gray level images.
 * 
 * @see find_candidates
 */
int find_candidates_gray( unsigned char *patch_list,unsigned char *patch_list_mask, unsigned char *curr_patch, int totalPatchs, int t, unsigned short *mask_neighb, int neighbs_updated, float tolerance, Cand_dist cand_list, float *weights ) 
{

  float radius;
  
  int h = compute_distances_gray( totalPatchs, t, neighbs_updated, mask_neighb, patch_list,patch_list_mask, curr_patch, tolerance, cand_list, &radius, weights );
 
  qsort (cand_list, (size_t) h, (size_t) sizeof (struct cand_dist),*compare_cand_dist);
  
  int size_cand_list = 1;
  while ( size_cand_list < h && cand_list[size_cand_list].dist <= radius ) {
   
    size_cand_list++;
  }
 
  return size_cand_list;

}

/**
 * @brief Retrieve the coordinates in the example image from a patch index.
 * 
 * @param   p       patch index
 * @param   ncols   number of columns (dx) in the image
 * @param   x       output for \em x position in the image
 * @param   y       output for \em y position in the image
 */
void retrieve_coords( int p, int ncols, int* x, int* y) {

  int x_aux;

  x_aux = (int)(p / ncols);
  *y = p - ncols * x_aux;
  *x = x_aux;  

}

/**
 * @brief   Compute the distances between a patch and all the others ones.
 * 
 * @param   totalPatchs    number of patches in the dictionary
 * @param   t               half-size of the patch
 * @param   neighbs_updated number of known neighbors
 * @param   mask_neighb     mask of known neighbors
 * @param   patch_list      the list of patches
 * @param   curr_patch      the current patch
 * @param   tolerance       value of \f$(1+\varepsilon)^2\f$
 * @param   cand_list       output for the distances
 * @param   radius_out      threshold for distance (computation stops if reached)
 * @param   weights         weights used for the distance between patches
 */
 int compute_distances( int totalPatchs, int t, int neighbs_updated, unsigned short *mask_neighb, unsigned char *patch_list,unsigned char *patch_list_mask, unsigned char *curr_patch, float tolerance, Cand_dist cand_list, float *radius_out, float *weights) 
{	
  int l,k, neighb;
  long int p;
  int patchSize = (2*t+1)*(2*t+1)*3;
  int patchSizeCh = (2*t+1)*(2*t+1);
  float dist_min = MAX_DIST;
  float radius;
  float diffR, diffG, diffB;
  int h = 0;  

  radius = dist_min*tolerance;  

#pragma omp parallel for private(p,l,k,neighb,diffR,diffG,diffB) shared(h,cand_list,radius,dist_min)  
  for ( p=0; p < totalPatchs ; p++) {

    float dist = 0.;
    float total_weight = 0.;
    int known_pixels = 0;
    l = 0;
    k = 0;    
   
    
    while ( dist <= radius && l < neighbs_updated) {

      neighb = *(mask_neighb + l);
      diffR = *(patch_list + p*patchSize + neighb) - *(curr_patch + k);  /* Esta bien que sea (curr_patch + l), mirar known_neighbours */      
      diffG = *(patch_list + p*patchSize + patchSizeCh + neighb) - *(curr_patch + k + 1);
      diffB = *(patch_list + p*patchSize + 2*patchSizeCh + neighb) - *(curr_patch + k + 2);
      
      float mask_w = (float) *(patch_list_mask + p*patchSize + neighb);
      
      total_weight += mask_w * weights[neighb];                 
      
      dist = dist + (diffR*diffR + diffG*diffG + diffB*diffB) * weights[neighb] * mask_w;      
      
      l++;
      k = k + 3;
            
      known_pixels += mask_w;

    }
  
    dist = dist / total_weight;		
    
    	
    if ( dist < dist_min && known_pixels == neighbs_updated && l == neighbs_updated ) {			
#pragma omp critical      
      if ( dist < dist_min ) {        
		dist_min = dist;
		radius = dist_min * tolerance; 
      }      
    }      
    
    
    if ( dist <= radius && l == neighbs_updated && known_pixels == neighbs_updated) { /* If completed the previous two calculs of distances is because dist < radius, so include it in the candidates list */

#pragma omp critical
      if ( dist <= radius ) {
		  	
		cand_list[h].dist = dist;
		cand_list[h].patch = p;
		h++;
		
      }      
    }   
	
  }
    
  *radius_out = radius;
  
  return h;
	
}

/**
 * @brief   Compute the distances between a patch and all the others ones.
 * 
 * Idem compute_distances for gray images.
 * 
 * @see compute_distances
 */
 int compute_distances_gray( int totalPatchs, int t, int neighbs_updated, unsigned short *mask_neighb, unsigned char *patch_list,unsigned char *patch_list_mask, unsigned char *curr_patch, float tolerance, Cand_dist cand_list, float *radius_out, float *weights) 
{	
  int l,k, neighb;
  long int p;
  int patchSize = (2*t+1)*(2*t+1);
  float dist_min = MAX_DIST;
  float radius;
  float diffVal;
  int h = 0;  

  radius = dist_min*tolerance;  

#pragma omp parallel for private(p,l,k,neighb,diffVal) shared(h,cand_list,radius,dist_min)  
  for ( p=0; p < totalPatchs ; p++) {

    float dist = 0.;
    float total_weight = 0.;
    int known_pixels = 0;
    l = 0;
    k = 0;    
   
    
    while ( dist <= radius && l < neighbs_updated) {

      neighb = *(mask_neighb + l);
      diffVal = *(patch_list + p*patchSize + neighb) - *(curr_patch + k);  /* Esta bien que sea (curr_patch + l), mirar known_neighbours */      
      
      float mask_w = (float) *(patch_list_mask + p*patchSize + neighb);
      
      total_weight += mask_w * weights[neighb];
      
      dist = dist + diffVal*diffVal * weights[neighb] * mask_w;
      
      l++;
      k++;
            
      known_pixels += mask_w;

    }
  
    dist = dist / total_weight;		
    
    
    if ( dist < dist_min && known_pixels == neighbs_updated && l == neighbs_updated ) {			
#pragma omp critical
      
      if ( dist < dist_min ) {        
	dist_min = dist;
	radius = dist_min * tolerance; 
      }      
    }      
    
    
    if ( dist <= radius && l == neighbs_updated && known_pixels == neighbs_updated) { /* If completed the previous two calculs of distances is because dist < radius, so include it in the candidates list */

#pragma omp critical
      if ( dist <= radius ) {
		
	cand_list[h].dist = dist;
	cand_list[h].patch = p;
	h++;
      }      
    }   
	
  }
    
  *radius_out = radius;
  
  return h;
	
}

/**
 * @brief   Pre-compute Gaussian weights.
 * 
 * @param   patchSize   size of the patch
 * @param   t           half-size of the patch (should have patchSize = 2t+1)
 * 
 * @return  array of Gaussian weights, length = patchSize x patchSize
 */
float * gaussian_weights( int patchSize, int t ) 
{
  int i,j;      
  float* gaussian = (float *) malloc (patchSize * patchSize * sizeof (float));
  if (!gaussian)
    error ("Not enough memory");    
    
  float s2 = (patchSize / 6.4) * (patchSize / 6.4);
  
  for (i = 0; i < patchSize; i++)
    for (j = 0; j < patchSize; j++)
      {
        gaussian[j + i * patchSize] = (1 / (2 * M_PI * s2)) * exp(-((i - t) * (i - t) + (j - t) * (j - t)) / (2 * s2));        
      }

  return gaussian;

}



/**
 * @brief   Pre-compute uniform weights.
 * 
 * @param   patchSize   size of the patch
 * 
 * @return  array of uniform weights, length = patchSize x patchSize
 */
 float * uniform_weights( int patchSize ) 
{
  int i;      
  float* uniform = (float *) malloc (patchSize * patchSize * sizeof (float));
  if (!uniform)
    error ("Not enough memory");      
  
  for (i = 0; i < patchSize*patchSize; i++) {
      
    uniform[i] = 1.0;   
         
  }

  return uniform;

}

/**
 * @brief   Comparison between pixel's number of neighbors.
 * 
 * This function is used by \em qsort to sort a list of pixel
 * in decreasing order of known neighbors.
 * 
 * @param   a   pointer to a pixel
 * @param   b   pointer to another pixel
 * 
 * @return  a negative, zero, or positive value respectively if the
 *      number of known neighbors is more, the same, or less
 *      for \em a than for \em b.
 */
int compare_pixel(const void * a, const void * b)
{
  Pixel aa = (Pixel) a;
  Pixel bb = (Pixel) b;
  /* Inverse order */
  return ((bb->nbneighbs) - (aa->nbneighbs));
}

/**
 * @brief   Comparison between candidates' distances.
 * 
 * This function is used by \em qsort to sort a list of candidates.
 * 
 * @param   a   pointer to a candidate
 * @param   b   pointer to another candidate
 * 
 * @return  a negative, zero, or positive value respectively if
 *      the distance of \em a is smaller than, equal to, or greater than
 *      the distance of \em b.
 */
int compare_cand_dist(const void * a, const void * b)
{
  Cand_dist aa = (Cand_dist) a;
  Cand_dist bb = (Cand_dist) b;
  return ((int) (aa->dist - bb->dist));
}
