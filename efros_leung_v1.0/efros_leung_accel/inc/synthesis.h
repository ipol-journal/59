#pragma once

#include <stdlib.h>
#include "pca.h"
#include "mt.h"

/*int minimum( float *array, int nthreads );*/
int compare_pixel(const void * a, const void * b);
/*int compare_int(const void * a, const void * b);*/
int compare_cand_dist(const void * a, const void * b);

int known_neighbours( ImageF w, unsigned char* mask, Pixel current, int i, int t, int larg, int haut, unsigned short *mask_neighb, float *curr_patch );
int known_neighboursRL( unsigned char* mask, Pixel current, int i, int t, int larg, unsigned short *mask_neighb );
int known_neighboursUD( unsigned char* mask, Pixel current, int i, int t, int larg, unsigned short *mask_neighb );

int find_candidatesU( ImageF v, float* weights_pca, float *dictionary, Pixel current, int i, int totalPatches, unsigned short *mask_neighb, int total_neighb, int t, float* WU, float *mean, int m, int n, float tolerance, Cand_dist cand_list);
int find_candidatesD( ImageF v, float* weights_pca, float *dictionary, Pixel current, int i, int totalPatches, unsigned short *mask_neighb, int total_neighb, int t, float* WU, float *mean, int m, int n, float tolerance, Cand_dist cand_list);
int find_candidatesR( ImageF v, float* weights_pca, float *dictionary, Pixel current, int i, int totalPatches, unsigned short *mask_neighb, int total_neighb, int t, float* WU, float* mean, int m, int n, float tolerance, Cand_dist cand_list);
int find_candidatesL( ImageF v, float* weights_pca, float *dictionary, Pixel current, int i, int totalPatches, unsigned short *mask_neighb, int total_neighb, int t, float* WU, float *mean, int m, int n, float tolerance, Cand_dist cand_list);
int find_candidates( float *patch_list, float *patch_list_mask, float *curr_patch, int totalPatchs, int t, unsigned short *mask_neighb, int voisins_updated, float tolerance, Cand_dist cand_list, float* weights );
float * load_patch_list( ImageF w, int t );

long int generate_current_list ( int cote, int coin, Pixel current, unsigned char *mask, int t, int larg, int haut );

long int random_choose ( Cand_dist cand_list, long int size_list );
void retrieve_coords( int p, int larg, int tx, int ty, int* x, int* y);

int compute_distances_pca( int totalPatches, int t, float tolerance, int total_neighb, unsigned short *mask_neighb, float *patch, float *patch_pca, float *dictionary, int n, Cand_dist cand_list, float *radius_out );
int compute_distances( int totalPatchs, int t, int voisins_updated, unsigned short *mask_neighb, float *patch_list, float *patch_list_mask, float *curr_patch, float tolerance, Cand_dist cand_list, float *radius_out, float *weights );

float * gaussian_weights( int patchSize, int t );
float * gaussian_weights_pca( int patchSize, int t );
float * uniform_weights_pca( int patchSize );




