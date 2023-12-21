#pragma once

#include "auxiliary.h"
#include "linear_algebra.h"
#include "svd.h"

float *pca_generate_base( float *X, int m, int n, int k );
float *pca_project( float *X, float *W, int k, int m, int n );
float *pca_project_vector( float *x, float *W, float *mean, int m, int k );
void pca_generate_matrix( ImageF w, int t, float *gauss_wg, float **dictionaryU_out,  float **dictionaryD_out, float **dictionaryR_out, float **dictionaryL_out, float **meanU, float **meanD, float **meanR, float **meanL  );
void pca( ImageF w, int t, int k, float* gauss_wg, float **WU, float **WD, float **WR, float **WL, float **YU, float **YD, float **YR, float **YL, float **meanU, float **meanD, float **meanR, float **meanL );
void construct_dictionary_pca( ImageF w, float* gauss_wg_pca, int t, int k, float *YU, float *YD, float *YR, float *YL, float **dictionaryR, float **dictionaryL, float **dictionaryU, float **dictionaryD );
ImageF pca_rgb( Cimage w );
