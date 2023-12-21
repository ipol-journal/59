#pragma once

#include "synthesis.h"
#include "auxiliary.h"

Cimage efros_leung_synth(int t, int out_img_sz, float tol, int init, int dims_pca, Cimage w, Cimage *map, Cimage *copy_map);
