#include <stdlib.h>
#include <cblas.h>
#include <float.h>
#include "svd.h"
#include "auxiliary.h"


static int sgesdd( char jobz, int m, int n, float *A, int lda, float *s, float *U, int ldu, float *VT, int ldvt, float *work, int lwork, int *iwork )

{
  extern void sgesdd_( char *jobz, int *m, int *n, float *A, int *lda, float *s, float *U, int *ldu, float *VT, int *ldvt, float *work, int *lwork, int *iwork, int *info );
  
  int info;
  sgesdd_( &jobz, &m, &n, A, &lda, s, U, &ldu, VT, &ldvt, work, &lwork, iwork, &info );
  
  return info;
}

/* Input: Matrix A of size m x n. m = dimensions, n = number of samples (patches) */
/* The code assumes m < n */

float *svd( float *A, int m, int n, float **s_out )
{
  char jobz = 'S';   /* The first min(M,N) columns of U and the first  min(M,N) rows of V**T are returned in the arrays U and VT */
  
  int lda = m; /* The leading dimension of the array A.  LDA >= max(1,M). */

  float* s;  /* Output: The singular values of A, sorted so that S(i) >= S(i+1) */

  float* U; /* Output: If JOBZ = 'S', U contains the first min(M,N) columns of U (the left singular vectors, stored columnwise). */

  int ldu = m;  /* The leading dimension of the array U. */

  float *VT;  /* Output: if JOBZ = 'S', VT contains the first min(M,N) rows of V**T (the right singular vectors, stored rowwise). */

  /* The leading dimension of the array VT.  */
  int ldvt = m;   

  float *work;   /* Output */
  float lwork;     /* If LWORK = -1, then a workspace query is assumed */
  int *iwork;    /* Integer array, dimension (8*min(M,N))*/ 
  
  /*int info;*/
        
  s = (float *) malloc(m * sizeof(float));
  if (!s)
    error ("Not enough memory");

  U = (float *) malloc(m * m * sizeof(float));
  if (!U)
    error ("Not enough memory");

  VT = (float *) malloc(m * n * sizeof(float));
  if (!VT)
    error ("Not enough memory for VT");

  iwork = (int *) malloc(8 * m * sizeof(int));
  if (!iwork)
    error ("Not enough memory");
 
  /*Do a query to know the optimum BufferSize */
  sgesdd( jobz, m, n, A, lda, s, U, ldu, VT, ldvt, &lwork, -1, iwork );
        
  work = (float *) malloc(lwork * sizeof(float));
  if (!work)
    error ("Not enough memory");

  sgesdd( jobz, m, n, A, lda, s, U, ldu, VT, ldvt, work, (int) lwork, iwork );
        
  *s_out = s;

  free(VT);
  free(iwork);
  free(work);

  return U;
        
}
