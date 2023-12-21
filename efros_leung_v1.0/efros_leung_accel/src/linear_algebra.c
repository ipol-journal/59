#include <stdio.h> 
#include <stdlib.h>
#include <cblas.h>
#include <float.h>
#include "linear_algebra.h"

static void sgemm( char transaA, char transaB, int m, int n, int k, float alpha, float *A, int lda, float *B, int ldb, float beta, float *C, int ldc)
/* *  Purpose */
/* *  ======= */
/* * */
/* *  SGEMM  performs one of the matrix-matrix operations */
/* * */
/* *     C := alpha*op( A )*op( B ) + beta*C, */
/* * */
/* *  where  op( X ) is one of */
/* * */
/* *     op( X ) = X   or   op( X ) = X**T, */
/* * */
/* *  alpha and beta are scalars, and A, B and C are matrices, with op( A ) */
/* *  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */


/*   TRANSA - CHARACTER*1. */
/*          On entry, TRANSA specifies the form of op( A ) to be used in */
/*           the matrix multiplication as follows: */

/*               TRANSA = 'N' or 'n',  op( A ) = A. */

/*               TRANSA = 'T' or 't',  op( A ) = A**T. */

/*               TRANSA = 'C' or 'c',  op( A ) = A**T. */

/* *  M      - INTEGER. */
/* *           On entry,  M  specifies  the number  of rows  of the  matrix */
/* *           op( A )  and of the  matrix  C.  M  must  be at least  zero. */
/* *           Unchanged on exit. */
/* * */
/* *  N      - INTEGER. */
/* *           On entry,  N  specifies the number  of columns of the matrix */
/* *           op( B ) and the number of columns of the matrix C. N must be */
/* *           at least zero. */
/* *           Unchanged on exit. */
/* * */
/* *  K      - INTEGER. */
/* *           On entry,  K  specifies  the number of columns of the matrix */
/* *           op( A ) and the number of rows of the matrix op( B ). K must */
/* *           be at least  zero. */
/* *           Unchanged on exit. */
{
  extern void sgemm_( char *transaA, char *transaB, int *m, int *n, int *k, float *alpha, float *A, int *lda, float *B, int *ldb, float *beta, float *C, int *ldc);

    
  sgemm_( &transaA, &transaB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
            
}


static void sgemv( char trans, int m, int n, float alpha, float *A, int lda, float *x, int incx, float beta, float *y, int incy )
/* *  Purpose */
/* *  ======= */
/* * */
/* *  SGEMV  performs one of the matrix-vector operations */
/* * */
/* *     y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y, */
/* * */
/* *  where alpha and beta are scalars, x and y are vectors and A is an */
/* *  m by n matrix. */
/* * */
/* *  Arguments */
/* *  ========== */
/* * */
/* *  TRANS  - CHARACTER*1. */
/* *           On entry, TRANS specifies the operation to be performed as */
/* *           follows: */
/* * */
/* *              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */
/* * */
/* *              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y. */
/* * */
/* *              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y. */
/* * */
/* *           Unchanged on exit. */
/* * */
/* *  M      - INTEGER. */
/* *           On entry, M specifies the number of rows of the matrix A. */
/* *           M must be at least zero. */
/* *           Unchanged on exit. */
/* * */
/* *  N      - INTEGER. */
/* *           On entry, N specifies the number of columns of the matrix A. */
/* *           N must be at least zero. */
/* *           Unchanged on exit. */
/* * */
/* *  ALPHA  - REAL            . */
/* *           On entry, ALPHA specifies the scalar alpha. */
/* *           Unchanged on exit. */
/* * */
/* *  A      - REAL             array of DIMENSION ( LDA, n ). */
/* *           Before entry, the leading m by n part of the array A must */
/* *           contain the matrix of coefficients. */
/* *           Unchanged on exit. */
/* * */
/* *  LDA    - INTEGER. */
/* *           On entry, LDA specifies the first dimension of A as declared */
/* *           in the calling (sub) program. LDA must be at least */
/* *           max( 1, m ). */
/* *           Unchanged on exit. */
/* * */
/* *  X      - REAL             array of DIMENSION at least */
/* *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/* *           and at least */
/* *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/* *           Before entry, the incremented array X must contain the */
/* *           vector x. */
/* *           Unchanged on exit. */
/* * */
/* *  INCX   - INTEGER. */
/* *           On entry, INCX specifies the increment for the elements of */
/* *           X. INCX must not be zero. */
/* *           Unchanged on exit. */
/* * */
/* *  BETA   - REAL            . */
/* *           On entry, BETA specifies the scalar beta. When BETA is */
/* *           supplied as zero then Y need not be set on input. */
/* *           Unchanged on exit. */
/* * */
/* *  Y      - REAL             array of DIMENSION at least */
/* *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/* *           and at least */
/* *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/* *           Before entry with BETA non-zero, the incremented array Y */
/* *           must contain the vector y. On exit, Y is overwritten by the */
/* *           updated vector y. */
/* * */
/* *  INCY   - INTEGER. */
/* *           On entry, INCY specifies the increment for the elements of */
/* *           Y. INCY must not be zero. */
/* *           Unchanged on exit. */
{
  extern void sgemv_( char *trans, int *m, int *n, float *alpha, float *A, int *lda, float *x, int *incx, float *beta, float *y, int *incy );
    
  sgemv_( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
            
}


float * matrix_prod( float *A, float *B, int m, int n, int k, char transaA, char transaB )
{

  /* op(A): m x k*/
  /* op(B): k x n*/
  /* C: m x n*/

  float alpha = 1;
  float beta = 0;

  int lda = ( transaA == 'N' ) ? m : k; 
  int ldb = ( transaB == 'N' ) ? k : n;
  int ldc = m;

  float *C;


  C = (float *) malloc(m * n * sizeof(float));
 
  sgemm( transaA, transaB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc); 
        
  return C;
        
}


float * matrix_vector_prod( float *A, float *x, int m, int n, char trans, float alpha )
{
  /* A: matrix of size m x n */

  float beta = 0;

  int lda = m;
  int incx = 1;
  int incy = 1;

  float *y;
  int sizeY;

  sizeY = ( trans == 'N' ) ? m : n;
  
  y = (float *) malloc(sizeY * sizeof(float));

  sgemv( trans, m, n, alpha, A, lda, x, incx, beta, y, incy );

  return y;
}



