#include "pca.h"

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

/**
 * @file pca.h
 * @brief   PCA (Principal Components Analysis) functions.
 */

/**
 * @brief   Generate the matrix of \f$k\f$ PCA axis.
 * 
 * Compute the PCA axis--i.e. the eigenvectors of \f$XX^t\f$
 * and also the \f$W\f$ matrix of size \f$m\times m\f$
 * in the SVD of \f$X=WSV^t\f$.
 * 
 * Select the \f$k\f$ eigenvectors with the largest eigenvalues.
 * 
 * @param   X   data matrix, size \f$m\times n\f$, centered
 * @param   m   dimension of the data
 * @param   n   number of data
 * @param   k   number of components to be kept    
 * 
 * @warning \f$X\f$ is assumed centered
 * @return  Matrix of size \f$m\times k\f$:
 *      the \f$k\f$ PVA axis, normalized.
 */
float *pca_generate_base( float *X, int m, int n, int k )
{
  float *W, *eigenvals, *W_k;
    
  W = svd( X, m, n, &eigenvals );  

  W_k = (float *) malloc(m * k * sizeof(float));
  if (!W_k)
    error ("Not enough memory");

  memcpy(W_k, W, k*m*sizeof(float));  
  
  free(W);

  return W_k;
}

/**
 * @brief   Project the data on some axis
 * 
 * @param   X   data matrix, size \f$m\times n\f$, centered
 * @param   W   projection axis, size \f$m\times k\f$, normalized columns
 * @param   k   number of dimension kept by PCA
 * @param   m   dimension of the data
 * @param   n   number of data
 * 
 * @warning \f$X\f$ is assumed centered
 * @return  Projection \f$Y=W^tX\f$, size \f$m\times n\f$
 */
float *pca_project( float *X, float *W, int k, int m, int n )
{
  float *Y;

  Y = matrix_prod( W, X, k, n, m, 'T', 'N' );  

  return Y;
}

/**
 * @brief   Project a vector on some axis
 * 
 * @param   x   data vector, size \f$m\f$
 * @param   mean    mean vector to use for \f$x\f$, size \f$m\f$
 * @param   W   projection axis, size \f$m\times k\f$, normalized columns
 * @param   k   number of dimension kept by PCA
 * @param   m   dimension of the data
 * 
 * @return  Projection \f$y=W^tx\f$, size \f$k\f$
 */
float *pca_project_vector( float *x, float *W, float* mean, int m, int k )
{
  int i;
  float *y;
  float alpha = 1.;

  float *x_centered = (float *) malloc(m*sizeof(float));
  if (!x_centered)
    printf("Not enough memory");

  /* Center data vector before projection */
#pragma omp parallel for 
  for (i = 0; i < m; i++) {
    x_centered[i] = x[i] - mean[i];
  }

  y = matrix_vector_prod( W, x_centered, m, k, 'T', alpha ); 

  return y;
}

/**
 * @brief   Extract all half-patches from the image.
 * 
 * Built dictionaries will be matrices of size \f$m\cdot n\f$
 * containing \f$n\f$ patches of dimensions \f$m=t\cdot T\f$ .
 * 
 * Consecutive values are all the components of a same patch
 * (the columns of the dictionaries) with respect to Lapack
 * column major order.
 * 
 * @param   w               the image
 * @param   t               half-size of the patch
 *      (full size is \f$T=2\cdot t+1\f$)
 * @param   gauss_wg        Gaussian weights
 * @param   dictionaryU_out pointer to the output dictionary,
 *      dictionary will be extracted from all the upper half-patches
 * @param   dictionaryD_out the same with the down part of patches
 * @param   dictionaryR_out the same with the right part of patches
 * @param   dictionaryL_out the same with the left part of patches
 * @param   meanU           pointer to the output mean,
 *      computed on all the upper half-patches
 * @param   meanD           the same with the down part of patches
 * @param   meanR           the same with the right part of patches
 * @param   meanL           the same with the left part of patches
 */
void pca_generate_matrix( ImageF w, int t, float *gauss_wg,
			  float **dictionaryU_out, float **dictionaryD_out,
			  float **dictionaryR_out, float **dictionaryL_out,
			  float **meanU, float **meanD, float **meanR, float **meanL )
{
 
  int i, j, k, l, p; 
  int nrows = w->nrow;  
  int ncols = w->ncol;
  int KcolsUD = ncols - 2*t;
  int KrowsUD = nrows - t;
  int totalPatchUD = KcolsUD * KrowsUD;
  int KcolsRL = ncols - t;
  int KrowsRL = nrows - 2*t;
  int totalPatchRL = KcolsRL * KrowsRL;

  int T = 2*t + 1;
  int Tt = T*(t+1);
  int tPlusOne = t + 1;
  int addr1, addrU, addrD, addrR, addrL;

  float *dictionaryR, *dictionaryL, *dictionaryU, *dictionaryD;

  int indexU[totalPatchUD];
  int indexD[totalPatchUD];
  int indexR[totalPatchRL];
  int indexL[totalPatchRL];
  
  dictionaryU = (float *) malloc( t*T * totalPatchUD * sizeof(float));
  dictionaryD = (float *) malloc( t*T * totalPatchUD * sizeof(float));
  dictionaryR = (float *) malloc( t*T * totalPatchRL * sizeof(float));
  dictionaryL = (float *) malloc( t*T * totalPatchRL * sizeof(float));


  k = 0;
  for ( i=0; i < nrows - t ; i++) 
    for ( j=0; j <= ncols-T ; j++) { 	  
      indexU[k] = j + i * ncols;
      k++;
    }
  k = 0;
  for ( i=1; i < nrows - t +1 ; i++) 
    for ( j=0; j <= ncols-T ; j++) { 	  
      indexD[k] = j + i * ncols;
      k++;
    }
  k = 0;
  for ( i=0; i <= nrows-T ; i++) 
    for ( j=1; j < ncols-t+1 ; j++) { 	  
      indexR[k] = j + i * ncols;
      k++;
    }
  k = 0;
  for ( i=0; i <= nrows-T ; i++) 
    for ( j=0; j < ncols-t; j++) { 	  
      indexL[k] = j + i * ncols;
      k++;
    }

  /* Fill dictionary with known patches */
  for ( p=0; p < totalPatchUD ; p++) {

    addr1 = p*T*t;
    addrU = indexU[p];
    addrD = indexD[p];
 
    for ( k=0; k < t; k++ )      
      for ( l=0; l < T ; l++) {

        dictionaryU[addr1 + k*T + l] = w->val[addrU + k*ncols + l] * gauss_wg[ k*T + l];
        dictionaryD[addr1 + k*T + l] = w->val[addrD + k*ncols + l] * gauss_wg[ Tt + k*T + l];
      
      }  
  }
  /* Fill dictionary with known patches */
  for ( p=0; p < totalPatchRL ; p++) {

    addr1 = p*T*t;
    addrR = indexR[p];
    addrL = indexL[p];

    for ( k=0; k < T; k++) 
      for ( l=0; l < t ; l++) {

        dictionaryR[addr1 + k*t + l] = w->val[addrR + k*ncols + l] * gauss_wg[tPlusOne + k*T + l] ;
        dictionaryL[addr1 + k*t + l] = w->val[addrL + k*ncols + l] * gauss_wg[k*T + l]; 
                
      }
  }

  /* Center data to perform PCA */

  /* Compute the mean */
  /* Create unitary vector */
  float *unitary = (float *)malloc(totalPatchUD*sizeof(float));
  for (p = 0; p < totalPatchUD; p++)
    unitary[p] = 1.0;
  
  float alpha = 1.0 / totalPatchUD;
  int m = T*t;

  *meanU = matrix_vector_prod( dictionaryU, unitary, m, totalPatchUD, 'N', alpha );
  *meanD = matrix_vector_prod( dictionaryD, unitary, m, totalPatchUD, 'N', alpha );

  float diffU, diffD;
  
  /* Substract the mean */
#pragma omp parallel for private(i,diffU,diffD)

  for ( p = 0; p < totalPatchUD; p++ )
    for ( i = 0; i < m; i++)
      {
        diffU = *(dictionaryU + p*m + i) - *(*meanU + i);
        *(dictionaryU + p*m + i) = diffU;

        diffD = *(dictionaryD + p*m + i) - *(*meanD + i);
        *(dictionaryD + p*m + i) = diffD;
      }
      
  free(unitary);
  
  /*-----------------------------------------------------------*/
  unitary = (float *)malloc(totalPatchRL*sizeof(float));
  for (p = 0; p < totalPatchRL; p++)
    unitary[p] = 1.0;
  
  alpha = 1.0 / totalPatchRL;
  m = T*t;

  *meanR = matrix_vector_prod( dictionaryR, unitary, m, totalPatchRL, 'N', alpha );
  *meanL = matrix_vector_prod( dictionaryL, unitary, m, totalPatchRL, 'N', alpha );

  float diffR, diffL;
  
  /* Substract the mean */
#pragma omp parallel for private(i,diffR,diffL)

  for ( p = 0; p < totalPatchRL; p++ )
    for ( i = 0; i < m; i++)
      {
        diffR = *(dictionaryR + p*m + i) - *(*meanR + i);
        *(dictionaryR + p*m + i) = diffR;

        diffL = *(dictionaryL + p*m + i) - *(*meanL + i);
        *(dictionaryL + p*m + i) = diffL;
      }

  
  *dictionaryU_out = dictionaryU;
  *dictionaryD_out = dictionaryD;
  *dictionaryR_out = dictionaryR;
  *dictionaryL_out = dictionaryL;

}

/**
 * @brief   Decompose all half-patches of the image in a PCA basis.
 * 
 * @param   w           the image
 * @param   t           half-size of the patch
 *      (full size is \f$T=2\cdot t+1\f$)
 * @param   k           number of PCA components to be kept
 * @param   gauss_wg    Gaussian weights
 * @param   WU          pointer to the output data, extracted from
 *      all the upper half-patches
 * @param   WD          the same with the down part of patches
 * @param   WR          the same with the right part of patches
 * @param   WL          the same with the left part of patches
 * @param   YU          pointer to the output projected data from
 *      all the upper half-patches
 * @param   YD          the same with the down part of patches
 * @param   YR          the same with the right part of patches
 * @param   YL          the same with the left part of patches
 * @param   meanU       pointer to the output mean for the upper parts
 * @param   meanD       pointer to the output mean for the down parts
 * @param   meanR       pointer to the output mean for the right parts
 * @param   meanL       pointer to the output mean for the left parts
 * 
 * @see     pca_generate_base, pca_project
 */
void pca( ImageF w, int t, int k, float* gauss_wg,
	  float **WU, float **WD, float **WR, float **WL,
	  float **YU, float **YD, float **YR, float **YL,
	  float **meanU, float **meanD, float **meanR,float **meanL ) {

  /* k = number of dimensions kept by PCA */  
  float *dictionaryU, *dictionaryD, *dictionaryR, *dictionaryL;
  float *dictionaryU_aux, *dictionaryD_aux, *dictionaryR_aux, *dictionaryL_aux;
  int T = 2*t + 1;
  int m = t * T;
  int nUD = (w->ncol - 2*t)*(w->nrow - t);
  int nRL = (w->ncol - t)*(w->nrow - 2*t);
  
  dictionaryU_aux = (float *) malloc( t*T * nUD * sizeof(float));
  dictionaryD_aux = (float *) malloc( t*T * nUD * sizeof(float));
  dictionaryR_aux = (float *) malloc( t*T * nRL * sizeof(float));
  dictionaryL_aux = (float *) malloc( t*T * nRL * sizeof(float));

  /* Generate matrix to compute PCA (in column major order as required by Lapack). */  
  pca_generate_matrix( w, t, gauss_wg, &dictionaryU, &dictionaryD, &dictionaryR, &dictionaryL, &(*meanU), &(*meanD), &(*meanR), &(*meanL) );  
  
  /* Generate the auxiliary dictionaries to compute the bases since Lapack overwrites memory. */  
  memcpy( dictionaryU_aux, dictionaryU, t*T * nUD * sizeof(float) );
  memcpy( dictionaryD_aux, dictionaryD, t*T * nUD * sizeof(float) );
  memcpy( dictionaryR_aux, dictionaryR, t*T * nRL * sizeof(float) );
  memcpy( dictionaryL_aux, dictionaryL, t*T * nRL * sizeof(float) );

  /* Generate the corresponding PCA base for each dictionary (U, D, R, L). */  
  *WU = pca_generate_base( dictionaryU_aux, m, nUD, k );
  *WD = pca_generate_base( dictionaryD_aux, m, nUD, k );
  *WR = pca_generate_base( dictionaryR_aux, m, nRL, k );
  *WL = pca_generate_base( dictionaryL_aux, m, nRL, k );  
  
  /* Project each dictionary over its corresponding base keeping k dimensions. */  
  *YU = pca_project( dictionaryU, *WU, k, m, nUD );
  *YD = pca_project( dictionaryD, *WD, k, m, nUD );
  *YR = pca_project( dictionaryR, *WR, k, m, nRL );
  *YL = pca_project( dictionaryL, *WL, k, m, nRL );
  
  /* Free memory*/
  free(dictionaryU_aux);
  free(dictionaryD_aux);
  free(dictionaryR_aux);
  free(dictionaryL_aux);
  free(dictionaryU);
  free(dictionaryD);
  free(dictionaryR);
  free(dictionaryL);
  
}

/**
 * @brief   Extract dictionaries and PCA projections from the image.

 * @param   w               the image
 * @param   t               half-size of the patch
 *      (full size is \f$T=2\cdot t+1\f$)
 * @param   gauss_wg_pca    Gaussian weights
 * @param   k               number of PCA components to be kept
 * @param   YU              pointer to the output projected data from
 *      all the upper half-patches
 * @param   YD              the same with the down part of patches
 * @param   YR              the same with the right part of patches
 * @param   YL              the same with the left part of patches
 * @param   dictionaryU     pointer to the output dictionary,
 *      dictionary will be extracted from all the upper half-patches
 * @param   dictionaryD     the same with the down part of patches
 * @param   dictionaryR     the same with the right part of patches
 * @param   dictionaryL     the same with the left part of patches
 */
void construct_dictionary_pca( ImageF w, float* gauss_wg_pca, int t, int k,
			       float *YU, float *YD, float *YR, float *YL,
			       float **dictionaryR, float **dictionaryL,
			       float **dictionaryU, float **dictionaryD )
{
  int i, j, p, h; 
  int nrows = w->nrow;  
  int ncols = w->ncol;
  int KcolsUD = ncols - 2*t;
  int KrowsUD = nrows - t;
  int totalPatchUD = KcolsUD * KrowsUD;
  int KcolsRL = ncols - t;
  int KrowsRL = nrows - 2*t;
  int totalPatchRL = KcolsRL * KrowsRL;
  int T = 2*t + 1;
  int Tt = T*t;
  int dicSize = k + T;
  int addr1;

  int indexU[totalPatchUD];
  int indexD[totalPatchUD];
  int indexR[totalPatchRL];
  int indexL[totalPatchRL];

  /* Keep k PCA coordinates and the T of the middle patch */
  *dictionaryU = (float *) malloc( (k + T) * totalPatchUD * sizeof(float));
  *dictionaryD = (float *) malloc( (k + T) * totalPatchUD * sizeof(float));
  *dictionaryR = (float *) malloc( (k + T) * totalPatchRL * sizeof(float));
  *dictionaryL = (float *) malloc( (k + T) * totalPatchRL * sizeof(float));

  /* index[h] keeps the index of each patch. This makes it easier to recover the image. */
  h = 0;
  for ( i=t; i < nrows ; i++) 
    for ( j=0; j <= ncols-T ; j++) { 	  
      indexU[h] = j + i * ncols;
      h++;
    }
  h = 0;
  for ( i=0; i < nrows-t ; i++) 
    for ( j=0; j <= ncols-T ; j++) { 	  
      indexD[h] = j + i * ncols;
      h++;
    }
  h = 0;
  for ( i=0; i <= nrows-T ; i++) 
    for ( j=0; j < ncols-t ; j++) { 	  
      indexR[h] = j + i * ncols;
      h++;
    }
  h = 0;
  for ( i=0; i <= nrows-T ; i++) 
    for ( j=t; j < ncols; j++) { 	  
      indexL[h] = j + i * ncols;
      h++;
    }

  /* Fill dictionary of patches known on the upper/down half */
  for ( p=0; p < totalPatchUD ; p++) {

    addr1 = p*dicSize;
   
    for (i=0;i<T;i++) {
      *( *dictionaryU + addr1 + i) =  w->val[indexU[p] + i] * gauss_wg_pca[Tt + i];     
      *(*dictionaryD + addr1 + i) = w->val[indexD[p] + i] * gauss_wg_pca[Tt + i];  
    }

    memcpy( *dictionaryU + addr1 + T, YU + p*k, k*sizeof(float) );   
    memcpy( *dictionaryD + addr1 + T, YD + p*k, k*sizeof(float) );   
  }
  
  /* Fill dictionary of patches known on the right/left half */
  for ( p=0; p < totalPatchRL ; p++) {

    addr1 = p*dicSize;
 
    for ( i=-t; i <= t ; i++) {

      *(*dictionaryR + addr1 + i + t) =  w->val[indexR[p] + (i+t)*ncols] * gauss_wg_pca[t + (t + i)*T];
      *(*dictionaryL + addr1 + i + t) =  w->val[indexL[p] + (i+t)*ncols] * gauss_wg_pca[t + (t + i)*T];    
    }

    memcpy( *dictionaryR + addr1 + T, YR + p*k, k*sizeof(float) );
    memcpy( *dictionaryL + addr1 + T, YL + p*k, k*sizeof(float) );
    
  }

}

/**
 * @brief   Convert a color image into a gray-scale image using a PCA on RGB channels.
 * 
 * @param   w   color image
 * 
 * @return  gray-scale image
 */
ImageF pca_rgb( Cimage w )
{
  int m = 3;
  int totalPixels = (w->nrow)*(w->ncol);
  int n = totalPixels;
  int k = 1;
  float diff;

  /* Copy matrix information in column major order as required by Lapack */
        
  float *w_rgb = (float *)malloc(totalPixels*m*sizeof(float));
  float *w_rgb_aux = (float *)malloc(totalPixels*m*sizeof(float));
    
  int i, p, h = 0;    
  for (i=0; i<totalPixels; i++) {               
    w_rgb[h] = (float)(w->red[i]);    
    h++;
    w_rgb[h] = (float)(w->green[i]);       
    h++;
    w_rgb[h] = (float)(w->blue[i]);       
    h++;      
  }     
    
  /* Compute the mean */
  /* Create unitary vector */
  float *unitary = (float *)malloc(totalPixels*sizeof(float));
  for (p = 0; p < totalPixels; p++)
    *(unitary + p) = 1;
  
  float alpha = 1. / totalPixels;
  float *mean = matrix_vector_prod( w_rgb, unitary, m, totalPixels, 'N', alpha );    
   
  /* Subtract the mean to center data*/
  for (p = 0; p < totalPixels; p++)
    for (i = 0; i<m; i++) {
      diff = w_rgb[i + p*m] - mean[i];
      w_rgb[i + p*m] = diff;
    }   
    
  memcpy(w_rgb_aux,w_rgb,totalPixels*3*sizeof(float));
                
  /* Generate PCA base */    
  float *rgb_base = pca_generate_base( w_rgb_aux, m, n, k );
        
  /* Project w into the generated PCA subspace */
  float *rgb_projection = pca_project( w_rgb, rgb_base, k, m, n );
        
  /* The output must be of type ImageF to interface with posterior processing */
  ImageF pca_image = new_ImageF(w->nrow, w->ncol);  
        
  /* Must write data in row major order as it was before for image data */
  memcpy(pca_image->val, rgb_projection , totalPixels*sizeof(float));   
        
  return pca_image;
        
}
