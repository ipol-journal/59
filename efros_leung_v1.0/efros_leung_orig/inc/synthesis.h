#pragma once

#include <stdlib.h>
#include "auxiliary.h"
#include "mt.h"

int compare_pixel(const void * a, const void * b);
int compare_cand_dist(const void * a, const void * b);

unsigned char *load_patch_list( Cimage w, int t );
unsigned char *load_patch_list_gray( Gimage w, int t );

float * gaussian_weights( int patchSize, int t );
float * uniform_weights( int patchSize );

long int generate_current_list ( int cote, int coin, Pixel current, unsigned char *mask, int t, int larg, int haut );
int known_neighbours( Cimage v, unsigned char* mask, Pixel current, int i, int t, int larg, int haut, unsigned short *mask_neighb, unsigned char *curr_patch );
int known_neighbours_gray( Gimage v, unsigned char* mask, Pixel current, int i, int t, int ncols_v, int nrows_v, unsigned short *mask_neighb, unsigned char *curr_patch );

int find_candidates( unsigned char *patch_list,unsigned char *patch_list_mask, unsigned char *curr_patch, int totalPatchs, int t, unsigned short *mask_neighb, int voisins_updated, float tolerance, Cand_dist cand_list, float *weights );
int find_candidates_gray( unsigned char *patch_list,unsigned char *patch_list_mask, unsigned char *curr_patch, int totalPatchs, int t, unsigned short *mask_neighb, int voisins_updated, float tolerance, Cand_dist cand_list, float *weights );

int compute_distances( int totalPatchs, int t, int voisins_updated, unsigned short *mask_neighb, unsigned char *patch_list,unsigned char *patch_list_mask, unsigned char *curr_patch, float tolerance, Cand_dist cand_list, float *radius_out, float *weights);
int compute_distances_gray( int totalPatchs, int t, int neighbs_updated, unsigned short *mask_neighb, unsigned char *patch_list,unsigned char *patch_list_mask, unsigned char *curr_patch, float tolerance, Cand_dist cand_list, float *radius_out, float *weights);

long int random_choose ( Cand_dist cand_list, long int size_list );
void retrieve_coords( int p, int ncols, int* x, int* y); 

