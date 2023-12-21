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
  
  /* Decreasing order */
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

/**
 * @brief   Extract the neighbors' mask and values at a given pixel.
 * 
 * @param   v       the image
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
int known_neighbours( ImageF v, unsigned char* mask, Pixel current, int i,
		      int t, int ncols, int nrows, unsigned short *mask_neighb, float *curr_patch )
{
  int k, l;
  int counter = 0;
  unsigned short neighb = 0;
  
  /* For all the pixels in the patch centered at pixel i (current pixel being filled), find which pixels are already filled. */
  for (k = -t; k < t + 1; k++)   
    for (l = -t; l < t + 1; l++) {
      
      /* Border check: only consider the pixels on the patch that are included in image. */
      if (((current[i].y + l) >= 0) && ((current[i].y + l) < ncols) && ((current[i].x + k) >= 0) && ((current[i].x + k) < nrows))
        
        /* Check if pixel (k,l) is filled */  
        if ( mask[current[i].y + l + (current[i].x + k) * ncols] == 1 ) {
          
          /* If pixel (k,l) is filled, save in mask_neighb the number identifying the known neighbour. */
          mask_neighb[counter] = neighb;     
          
          /* Save also the current patch for posterior computations */   
          curr_patch[counter] = v->val[current[i].y + l + (current[i].x + k) * ncols]; 
          
          counter++;
        }
      neighb++;
    }
  
  return counter;
  
}

/**
 * @brief   Extract neighbors on the centered column of a patch.
 * @see     known_neighbours
 */
int known_neighboursRL( unsigned char* mask, Pixel current, int i, int t, int ncols, unsigned short *mask_neighb )
{
  int k;
  int counter = 0;
  unsigned short neighb = 0;
  
  /* For all the pixels in the central column of the patch centered at pixel i (current pixel being filled), */
  /* find which pixels are already filled. */
  for (k = -t; k <= t ; k++)  { 
    
    /* If pixel (x + k,y) is filled, save in mask_neighg the number identifying the known neighbour. */
    if (mask[current[i].y + (current[i].x + k) * ncols] == 1) {
      
      mask_neighb[counter] = neighb;      
      counter++;
      
    }
    neighb++;
  }
    
  return counter;

}

/**
 * @brief   Extract neighbors on the centered row of a patch.
 * @see     known_neighbours
 */
int known_neighboursUD( unsigned char* mask, Pixel current, int i, int t, int ncols, unsigned short *mask_neighb )
{
  int k;
  int counter = 0;
  unsigned short neighb = 0;

  /* For all the pixels in the row part of the patch centered at pixel i (current pixel being filled), */
  /* find which pixels are already filled. */
  for (k = -t; k <= t ; k++)  {
          
    /* If pixel (x,y + k) is filled, save in mask_neighg the number identifying the known neighbour. */
    if (mask[current[i].y + k + current[i].x * ncols] == 1) {     

      mask_neighb[counter] = neighb; 
      counter++;
    }
    neighb++;
  }

  return counter;
}

/**
 * @brief   Get the set of candidates for a given pixel to be filled.
 * 
 * Only for a patch with a known upper-part.
 * 
 * The candidate pixels \f$q\f$ for a pixel \f$p\f$ is the set:
 * \f[\Omega_p=\left\{q\;\middle|\;d(p,q)\leq (1+\varepsilon)\cdot d_{min}(p)\right\}\f]
 * where
 * \f[d_{min}(p)=\min_q\left\{d(p,q)\right\}\f]
 * and the \f$l^2\f$ distance \f$d(p,q)\f$ is computed on the known
 * part of the patches centered at \f$p\f$ and \f$q\f$.
 * 
 * @param   v               the example image
 * @param   weights_pca     weights used for the distance between patches
 * @param   dictionary      dictionary of patches extracted from \em v
 * @param   current         list of pixels
 * @param   i               index of the current pixel \f$p\f$ in \em current list
 * @param   totalPatches    number of patches in the dictionary
 * @param   mask_neighb     mask of the known pixels in the current patch
 * @param   total_neighb    number of known neighbors
 * @param   t               half-size of the patch
 * @param   mean            mean vector for patches dictionary
 * @param   WU              projection axis for the patches
 * @param   m               dimension of the patch 
 * @param   n               number of dimension kept by PCA
 * @param   tolerance       value of \f$(1+\varepsilon)^2\f$
 * @param   cand_list       utput for the candidates (sorted from closest)
 * 
 * @return  the number of candidates
 */
int find_candidatesU( ImageF v, float* weights_pca, float *dictionary,
		      Pixel current, int i, int totalPatches, unsigned short *mask_neighb, int total_neighb, int t,
		      float* WU, float* mean, int m, int n, float tolerance, Cand_dist cand_list)
{
  /* WU matrix of size m x n */

  int k, l;  
  int ncols = v->ncol;
  int T = 2*t + 1; 
  int Ttimest = T*t; 
  int Tt = 0;
  float radius;
  
  float *patch = (float *)malloc(T*sizeof(float));
  float *patch_aux = (float *)malloc(T*t*sizeof(float));
  
  float *patchPlus_t = patch + t;
    
  /* Store central part of patch */
  for (l = -t; l <= t ; l++)    
    patchPlus_t[l] = v->val[current[i].y + l + current[i].x * ncols] * weights_pca[t + l + Ttimest];
  
  /* Store upper part of the patch to be projected on the PCA base */
  for (k = -t; k < 0; k++) {
    Tt = (k + t)*T;
    for (l = -t; l <= t ; l++){
      patch_aux[Tt + l + t] = v->val[current[i].y + l + (current[i].x + k) * ncols] * weights_pca[Tt + l + t];                      
    }
  }       
  
  /* Find the projection of patch_aux on the PCA base */
  float *patch_pca = pca_project_vector( patch_aux, WU, mean, m, n );  

  /* Compute the distances from the current patch to all the patches in the dictionary. 
   * Distances are given in cands_list.     
   *               
   * radius = (1 + tolerance) * dist_min                                                        
   */
  int h = compute_distances_pca( totalPatches, t, tolerance, total_neighb, mask_neighb, patch, patch_pca, dictionary, n, cand_list, &radius );
  
  /* Order the candidates list in increasing order of distance */
  qsort (cand_list, h, sizeof (struct cand_dist),*compare_cand_dist);
  
  /* Find all the candidates verifying  dist( p, current ) <= radius */
  int size_cand_list = 1; 
  while ( size_cand_list < h && (cand_list + size_cand_list)->dist <= radius){ 
    size_cand_list++;
  }

  return size_cand_list;

}

/**
 * @brief   Get the set of candidates for a given pixel to be filled.
 * 
 * Only for a patch with a known down part.
 * 
 * @see     find_candidatesU
 */
int find_candidatesD( ImageF v, float* weights_pca, float *dictionary,
		      Pixel current, int i, int totalPatches, unsigned short *mask_neighb, int total_neighb, int t,
		      float* WU, float* mean, int m, int n, float tolerance, Cand_dist cand_list)
{
  /* WU matrix of size m x n */

  int k, l;  
  int ncols = v->ncol;
  int T = 2*t + 1;  
  int Ttimest = T*t;
  int Tt = 0;
  int gauss_off = (t+1)*T + t;
  float radius;
  
  float *patch = (float *)malloc(T*sizeof(float));
  float *patch_aux = (float *)malloc(T*t*sizeof(float));

  float *patchPlus_t = patch + t;
  
  /* Store central part of patch */
  for (l = -t; l <= t ; l++)    
    patchPlus_t[l] = v->val[current[i].y + l + current[i].x * ncols] * weights_pca[t + l + Ttimest];

  /* Store bottom part of patch to be projected on the PCA base */
  for (k = 1; k <= t; k++) {
    Tt = (k - 1)*T;
    for (l = -t; l <= t; l++) {
      patch_aux[Tt + l + t] = v->val[current[i].y + l + (current[i].x + k) * ncols] * weights_pca[gauss_off + Tt + l];     
    }
  }   
  
  float *patch_pca = pca_project_vector( patch_aux, WU, mean, m, n );  

  int h = compute_distances_pca( totalPatches, t, tolerance, total_neighb, mask_neighb, patch, patch_pca, dictionary, n, cand_list, &radius );
  
  qsort (cand_list, h, sizeof (struct cand_dist),*compare_cand_dist);
  
  int size_cand_list = 1;
 
  while ( size_cand_list < h && (cand_list + size_cand_list)->dist <= radius){ 
    size_cand_list++;
  }

  return size_cand_list;

}


/**
 * @brief   Get the set of candidates for a given pixel to be filled.
 * 
 * Only for a patch with a known right part.
 * 
 * @see     find_candidatesU
 */
int find_candidatesR( ImageF v, float *weights_pca, float *dictionary,
		      Pixel current, int i, int totalPatches, unsigned short *mask_neighb, int total_neighb, int t,
		      float* WU, float* mean, int m, int n, float tolerance, Cand_dist cand_list)
{
  /* WU matrix of size m x n */
  
  int l,  k;  
  int ncols = v->ncol;
  int T = 2*t + 1;  
  int Tt = 0;
  int Ttimest = T*t;
  int gauss_off = t*(1+T);
  float radius;  

  float *patch = (float *)malloc(T*sizeof(float));
  float *patch_aux = (float *)malloc(T*t*sizeof(float)); 
  
  float *patchPlus_t = patch + t;
  
  /* Store central part of patch */
  for (l = -t; l < t + 1; l++)
    
    patchPlus_t[l] = v->val[current[i].y + (current[i].x + l) * ncols] * weights_pca[t + l*T + Ttimest];
         
  /* Store upper part of patch to be projected on the PCA base */
  for (k = -t; k <= t; k++) {
    Tt = (k + t)*t;
    for (l = 1; l <= t; l++){
      patch_aux[Tt + l -1] = v->val[current[i].y + l + (current[i].x + k) * ncols] * weights_pca[gauss_off + l + k*T];                        
    }
  }  
    
  float *patch_pca = pca_project_vector( patch_aux, WU, mean, m, n );  

  int h = compute_distances_pca( totalPatches, t, tolerance, total_neighb, mask_neighb, patch, patch_pca, dictionary, n, cand_list, &radius );
  
  qsort (cand_list, h, sizeof (struct cand_dist),*compare_cand_dist);
  
  int size_cand_list = 1;
 
  while ( size_cand_list < h && (cand_list + size_cand_list)->dist <= radius){ 
    size_cand_list++;
  }   

  return size_cand_list;

}

/**
 * @brief   Get the set of candidates for a given pixel to be filled.
 * 
 * Only for a patch with a known left part.
 * 
 * @see     find_candidatesU
 */
int find_candidatesL( ImageF v, float* weights_pca, float *dictionary,
		      Pixel current, int i, int totalPatches, unsigned short *mask_neighb, int total_neighb, int t,
		      float* WU, float* mean, int m, int n, float tolerance, Cand_dist cand_list)
{
  /* WU matrix of size m x n */

  int l,  k;  
  int ncols = v->ncol;
  int T = 2*t + 1;  
  int Tt = 0;
  int Ttimest = T*t;
  int gauss_off = t*(1+T);
  float radius;  

  float *patch = (float *)malloc(T*sizeof(float));
  float *patch_aux = (float *)malloc(T*t*sizeof(float)); 
  
  float *patchPlus_t = patch + t;
  
  /* Store central part of patch */
  for (l = -t; l < t + 1; l++)
    
    patchPlus_t[l] = v->val[current[i].y + (current[i].x + l) * ncols] * weights_pca[t + T*l + Ttimest];
        
  /* Store left part of patch to be projected on the PCA base */
  for (k = -t; k <= t; k++) {
    Tt = (k + t)*t;
    for (l = -t; l < 0; l++) {
      patch_aux[Tt + l + t] = v->val[current[i].y + l + (current[i].x + k) * ncols] * weights_pca[gauss_off + l + k*T];             
    }
  }
    
  float *patch_pca = pca_project_vector( patch_aux, WU, mean, m, n );  

  int h = compute_distances_pca( totalPatches, t, tolerance, total_neighb, mask_neighb, patch, patch_pca, dictionary, n, cand_list, &radius );
  
  qsort (cand_list, h, sizeof (struct cand_dist),*compare_cand_dist);
  
  int size_cand_list = 1;
 
  while ( size_cand_list < h && (cand_list + size_cand_list)->dist <= radius){ 
    size_cand_list++;
  }
  
  return size_cand_list;

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
long int random_choose ( Cand_dist cand_list, long int size_list ){

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
long int generate_current_list ( int side, int corner, Pixel current,
				 unsigned char *mask, int t, int ncols, int nrows )
{
  int i, j, k, l;
  long int list_sz, counter, counter2, counter3, length, alea;
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
  length = list_sz-1; 
      
  while (length > 1) 
    { 
      alea = (int) ((length+1) * mt_drand53()); 
      length -= 1; 
      temp = current[length]; 
      current[length] = current[alea]; 
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
 * @param   w   image
 * @param   t   half-size of patches
 * 
 * @return  the patches array extracted from the image
 */
float *load_patch_list( ImageF w, int t )
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
  float *patch_list = (float *)malloc(totalPatchs*patchSize*sizeof(float));

  k = 0;
  for ( i=0; i < Krows ; i++) 
    for ( j=0; j < Kcols ; j++) {   
      index[k] = init + j + i * ncols;      
      k++;
    }

  /* Fill dictionary of patches */
  for ( p=0; p < totalPatchs ; p++) {

    addr1 = p*patchSize;
    addrU = index[p] - init;
        
    for ( k=0; k < T; k++ ){      
      memcpy( patch_list + addr1 + k*T,  w->val + addrU + k*ncols, T*sizeof(float) );           
    }

  }
  
  return patch_list;

}

/**
 * @brief   Get the set of candidates for a given pixel to be filled.
 * 
 * Don't use the PCA. The dictionary is simply the patch list.
 * 
 * @see find_candidatesU
 */
int find_candidates( float *patch_list,float *patch_list_mask, float *curr_patch, int totalPatchs, int t,
		     unsigned short *mask_neighb, int neighbs_updated, float tolerance,
		     Cand_dist cand_list, float* weights ) 
{

  float radius;
  
  int h = compute_distances( totalPatchs, t, neighbs_updated, mask_neighb, patch_list, patch_list_mask, curr_patch, tolerance, cand_list, &radius, weights );
 
  qsort (cand_list, h, sizeof (struct cand_dist),*compare_cand_dist);
  
  int size_cand_list = 1;
  while ( size_cand_list < h && (cand_list + size_cand_list)->dist <= radius ) {
    size_cand_list++;
  }
 
  return size_cand_list;

}

/**
 * @brief Retrieve the coordinates in the example image from a patch index.
 * 
 * @param   p       patch index
 * @param   ncols   number of columns (dx) in the image
 * @param   t       half-size of the patch
 * @param   x       output for \em x position in the image
 * @param   y       output for \em y position in the image
 */
void retrieve_coords( int p, int ncols, int tx, int ty, int* x, int* y) {

  int x_aux;

  x_aux = (int)(p / ncols);
  *y = p - ncols * x_aux + ty;
  *x = x_aux + tx;  

}

/**
 * @brief   Compute (using the PCA) the distances between a patch and all the others ones.
 * 
 * Each considered element has two parts: the first part corresponds to the central 
 * part of the patch (central column or central row) of size (2t + 1); 
 * the second part corresponds to the PCA projection of one of four patch subregions 
 * (up, down, right or left), it is of size equal to the PCA chosen dimensions \em n.
 * 
 * @param   totalPatches    number of patches in the dictionary
 * @param   t               half-size of the patch
 * @param   tolerance       value of \f$(1+\varepsilon)^2\f$
 * @param   total_neighb    number of known neighbors
 * @param   mask_neighb     mask of known neighbors
 * @param   patch
 * @param   patch_pca       the current patch coordinates in the PCA basis
 * @param   dictionary      dictionary of patches in the PCA basis
 * @param   n               number of dimension kept by PCA
 * @param   cand_list       output for the distances
 * @param   radius_out      threshold for distance (computation stops if reached)
 * 
 * @todo    uncommented parameter \em patch
 */
int compute_distances_pca( int totalPatches, int t, float tolerance,
			   int total_neighb, unsigned short *mask_neighb, float *patch, float *patch_pca,
			   float *dictionary, int n, Cand_dist cand_list, float *radius_out )
{       
  int l, k;
  long int p;   
  int T = 2*t + 1;      
  float dist_min = MAX_DIST;
  float radius = dist_min * tolerance;
  int neighb;
  int h = 0;
  float dist, diff_cent, diff_pca;
  int tdicc = T + n; /* size of a patch in the dicctionary */ 
  
  
#pragma omp parallel for private(dist,l,neighb,diff_cent,k,diff_pca) shared(h,radius,dist_min)

  /*For each patch compute distance and update candidates list if it corresponds*/
  for ( p=0; p<totalPatches; p++) {   
           
    /* Distance of central part of the patch*/
    dist = 0.;    
    l = 0;
    
    if ( total_neighb > 0 ) { /* If there are other pixels known in the central part of the current pixel. */
      while ( l < total_neighb && dist <= radius ) {
      
        neighb = mask_neighb[l];
        diff_cent = patch[neighb] - dictionary[p*tdicc + neighb];
        dist = dist + diff_cent * diff_cent;
        l++;
      }
    
    } 
    
    k = 0;
    if ( l == total_neighb ) {  /* If all central known pixels where considered in the previous calculation */

      /* Distance of the PCA part of the patch*/       
      while ( k < n && dist <= radius ) {
        
        diff_pca = patch_pca[k] - dictionary[p*tdicc + T + k];
        dist = dist + diff_pca * diff_pca;
        
        k++;
        
      }      
    }    
    
    if ( dist < dist_min ) {
#pragma omp critical
      if ( dist < dist_min ) {                
        dist_min = dist;
        radius = dist_min * tolerance;
      }      
    }      

  
    /* If completed the previous two calculs of distances is because dist < radius, so include it in the candidates list */ 
    if ( l == total_neighb && k == n &&  dist <= radius) { 
#pragma omp critical
      if ( l == total_neighb && k == n &&  dist <= radius) {
        cand_list[h].dist = dist;
        cand_list[h].patch = p;
        h++;      
      }
          
    } 

  }  

  *radius_out = radius;  /* radius_out = dist_min * tolerance */
  
  return h;   /* h = number of candidates in the candidate list */

}

/**
 * @brief   Compute (without PCA) the distances between a patch and all the others ones.
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
int compute_distances( int totalPatchs, int t, int neighbs_updated, unsigned short *mask_neighb,
		       float *patch_list, float *patch_list_mask,float *curr_patch, float tolerance,
		       Cand_dist cand_list, float *radius_out, float *weights ) 
{       
  int l, neighb;
  long int p;
  int patchSize = (2*t+1)*(2*t+1);
  float dist_min = MAX_DIST;
  float radius;
  float dist, diff;
  int h = 0; 

  radius = dist_min*tolerance;  
  
#pragma omp parallel for private(p,dist,l,neighb,diff) shared(h,cand_list,radius,dist_min)  
  for ( p=0; p < totalPatchs ; p++) {

    float total_weight = 0.;
    int known_pixels = 0;
    dist = 0.;
    l = 0;

    while ( dist <= radius && l < neighbs_updated){

      neighb = mask_neighb[l];
      diff = patch_list[p*patchSize + neighb] - curr_patch[l];    
      dist = dist + diff*diff * weights[neighb] * patch_list_mask[p*patchSize + neighb];
      total_weight += weights[neighb] * patch_list_mask[p*patchSize + neighb];
      
      known_pixels += patch_list_mask[p*patchSize + neighb];
      
      l++;
    }  
    
    dist = dist / total_weight;
    
    if ( dist < dist_min && known_pixels == neighbs_updated ) {           
#pragma omp critical
      if ( dist < dist_min ) {
        dist_min = dist;
        radius = dist_min * tolerance; 
      }      
    }      
    
    /* If completed the previous two calculs of distances is because dist < radius, so include it in the candidates list */    
    if ( dist <= radius && l == neighbs_updated && known_pixels == neighbs_updated ) { 

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
 * @brief   Pre-compute the root square of Gaussian weights.
 * 
 * @param   patchSize   size of the patch
 * @param   t           half-size of the patch (should have patchSize = 2t+1)
 * 
 * @return  array of root square of Gaussian weights, length = patchSize x patchSize
 */
float * gaussian_weights_pca( int patchSize, int t ) 
{
  int i,j;      
  float* gaussian = (float *) malloc (patchSize * patchSize * sizeof (float));
  if (!gaussian)
    error ("Not enough memory");

    
  float s2 = (patchSize / 6.4) * (patchSize / 6.4);
  
  for (i = 0; i < patchSize; i++)
    for (j = 0; j < patchSize; j++)
      {
        gaussian[j + i * patchSize] = (1 / sqrt(2 * M_PI * s2)) * exp(-((i - t) * (i - t) + (j - t) * (j - t)) / (4 * s2));
        
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
float * uniform_weights_pca( int patchSize ) 
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
