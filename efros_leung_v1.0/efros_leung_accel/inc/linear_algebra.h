#pragma once


float * matrix_prod( float *A, float *B, int m, int n, int k, char transaA, char transaB );
float * matrix_vector_prod( float *A, float *x, int m, int n, char trans, float alpha );
