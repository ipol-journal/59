#pragma once

#include "synthesis.h"
#include "auxiliary.h"

Cimage efros_leung_synth(int t, int out_img_sz, float tol, int init, Cimage w, Cimage *map_out, Cimage *copy_map_out);

Gimage efros_leung_synth_gray(int t, int out_img_sz, float tol, int init, Gimage w, Cimage *map_out, Cimage *copy_map_out);
