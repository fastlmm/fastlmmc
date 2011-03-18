#if !defined( BlasLapack_h )
#define BlasLapack_h

/*
 * This file sets up the appropriate math library to be used by FastLmmC
 *  Currently supported libraries are:
 *    Intel's Math Kernel Library (MKL) http://software.intel.com/en-us/articles/intel-mkl/
 *    AMD's Core Math Library (ACML) http://www.amd.com/acml
 *
 *  You must select one of the libraries by defining USE_MKL or USE_ACML when 
 *    compiling FastLmmC
 */
#if defined( USE_MKL )
#include "mkl.h"
#elif defined( USE_ACML )
#include "acml.h"
#else
#error USE_MKL or USE_ACML must be defined to set the appropriate math library usage
#endif

class BlasLapack
   {
public:
   /*
    * BLAS Level 1:
    */

   // vector scaling by a constant
   // x[::incx] := alpha * x[::incx]
   // x : [1+(n-1)*abs(incx)] vector
   // alpha double constant
   static void scal( const int *n, const double *alpha, double *x, const int *incx );

   // vector copy
   // Y[::incY] := X[::incX]
   // X : [1+(n-1)*abs(incX)] vector
   // Y : [1+(n-1)*abs(incY)] vector
   static void copy( const int *n, const double *X, const int *incX, double *Y, const int *incY );

   /*
    * BLAS Level 2:
    */
   
   // general matrix vector product
   // y[::incy] := alpha * ( trans(A[:,:]) * x[::incx] ) + beta * y[::incy]
   // trans :     'N' normal, 'T' transposed
   // A        : [ldA * n] matrix,  ldA = max( 1, m )
   // x for 'N': [1+(n-1)*abs(incX)] vector
   //   for 'T': [1+(m-1)*abs(incX)] vector
   // y for 'N': [1+(m-1)*abs(incX)] vector
   //   for 'T': [1+(n-1)*abs(incX)] vector
   static void gemv( const char *trans, const int *m, const int *n, const double *alpha, const double *A, const int *ldA, const double *x, const int *incx, const double *beta, double *y, const int *incy );

   /*
    * BLAS Level 3:
    */

   // general matrix matrix product
   // C[:,:] := alpha * ( transA(A[:,:]) * transB(B[:,:]) ) + beta * C[:,:]
   // transA, transB:    'N' normal, 'T' transposed
   // transA(A) : [m * k] matrix,    m >= 0, k >= 0
   // transB(B) : [k * n] matrix,    n >= 0
   // A  for 'N': [ldA * k] matrix,   ldA >= max( 1, m )
   //    for 'T': [ldA * m] matrix,   ldA >= max( 1, k )
   // B  for 'N': [ldB * n] matrix,   ldB >= max( 1, k )
   //    for 'T': [ldB * k] matrix,   ldB >= max( 1, n )
   // C         : [ldC * n] matrix,   ldC >= max( 1, m )
   static void gemm( const char *transA,  const char *transB, const int *m, const int *n, const int *k, const double *alpha, const double *A, const int *ldA, const double *B, int *ldB, const double *beta, double *C, const int *ldC );

   // symmetric rank k update. result will be in upper or lower storage format
   // C[:,:] := alpha * ( trans(A[:,:]) * trans(A[:,:]).T ) + beta * C[:,:]
   // uplo        : 'U' : upper storage , 'L' : lower storage
   // trans       : 'N' normal, 'T' transposed
   // A    for 'N': [ldA * k] matrix, ldA >= max( 1, n )
   //      for 'T': [ldA * n] matrix, ldA >= max( 1, k )
   static void syrk( const char *uplo, const char *trans, const int *n, const int *k, const double *alpha, double *A, const int *ldA, const double *beta, double *C, const int *ldC );

   /*
    * Lapack:
    */

   //singular value decomposition (SVD) of a [m * n] matrix A
   static void gesdd( const char *jobz, const int *m, const int *n, double *A, const int *ldA, double *S, double *U, const int *ldU, double *Vt, const int * ldVt, double *work, const int *lwork, int *iwork, int *info );
   
   //symmetric eigenvectors/eigenvalues
   static void syevr( const char *jobz, const char *range, const char *uplo, const int *n, double *A, const int *ldA, const double *vl, const double *vu, const int *il, const int *iu, const double *abstol, int *m, double *W, double *Z, const int *ldZ, int *isuppZ, double *work, const int *lwork, int *iwork, const int *liwork, int *info );
   
   //general matrix least squares solver
   static void gelsy( const int *m, const int *n, const int *nrhs, double *A, const int *ldA, double *B, const int *ldB, int *jpvt, const double *rcond, int *rank, double *work, const int *lwork, int *info );

   /*
    * VML:
    */
   
   // res[:] := Ln(x[:])
   static void vLn( const int n, const double *x, double *res );

   // res[:] := Sqrt(x[:])
   static void vSqrt( const int n, const double *x, double *res );
   };

//global variables:
const int      int_one=1;
const int      int_minusone=-1;
const int      int_zero=0;
const double   double_one=1.0;
const double   double_zero=0.0;
const double   double_minusone=-1.0;
const char     char_V='V';          //output the eigenvalues plus eigenvectors(jobz)
const char     char_A='A';          //compute all eigenvectors (range)
const char     char_L='L';          //the lower part of symmetrix matrix A is being read(uplo)
const char     char_N='N';          //work with the untransposed matrix(trans)
const char     char_T='T';          // Work with the transpose of the matrix(trans)
const char     char_S='S';          //comute the economy size SVD.(jobz)

/*
 * BLAS Level 1:
 */

// vector scaling by a constant
inline void BlasLapack::scal( const int *n, const double *alpha, double *x, const int *incx )
   {
#if defined( USE_MKL )
   dscal( n, alpha, x, incx );
#elif defined ( USE_ACML )
   dscal( *n, *alpha, x, *incx );
#endif
   }

// vector copy
inline void BlasLapack::copy( const int *n, const double *X, const int *incX, double *Y, const int *incY )
   {
#if defined( USE_MKL )
   dcopy( n, X, incX, Y ,incY );
#elif defined( USE_ACML )
   dcopy( *n, (double*)X, *incX, Y, *incY );
#endif
   }

/*
 * BLAS Level 2:
 */
   
// general matrix vector product
inline void BlasLapack::gemv( const char *trans, const int *m, const int *n, const double *alpha, const double *A, const int *ldA, const double *x, const int *incx, const double *beta, double *y, const int *incy )
   {
#if defined( USE_MKL )
   dgemv( trans, m, n, alpha, A, ldA, x, incx, beta, y, incy ); 
#elif defined( USE_ACML )
   dgemv( *trans, *m, *n, *alpha, (double*)A, *ldA, (double*)x, *incx, *beta, y, *incy ); 
#endif
   }

/*
 * BLAS Level 3:
 */

// general matrix matrix product
inline void BlasLapack::gemm( const char *transA,  const char *transB, const int *m, const int *n, const int *k, const double *alpha, const double *A, const int *ldA, const double *B, int *ldB, const double *beta, double *C, const int *ldC )
   {
#if defined( USE_MKL )
   dgemm( transA, transB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC );
#elif defined( USE_ACML )
   dgemm( *transA, *transB, *m, *n, *k, *alpha, (double*)A, *ldA, (double*)B, *ldB, *beta, C, *ldC );
#endif
   }

// symmetric rank k update. result will be in upper or lower storage format
inline void BlasLapack::syrk( const char *uplo, const char *trans, const int *n, const int *k, const double *alpha, double *A, const int *ldA, const double *beta, double *C, const int *ldC )
   {
#if defined( USE_MKL )
   dsyrk( uplo, trans, n, k, alpha, A, ldA, beta, C, ldC );
#elif defined( USE_ACML )
   dsyrk( *uplo, *trans, *n, *k, *alpha, A, *ldA, *beta, C, *ldC );
#endif
   }

/*
 * Lapack:
 */

//singular value decomposition (SVD) of a [m * n] matrix A
inline void BlasLapack::gesdd( const char *jobz, const int *m, const int *n, double *A, const int *ldA, double *S, double *U, const int *ldU, double *Vt, const int * ldVt, double *work, const int *lwork, int *iwork, int *info )
   {
#if defined( USE_MKL )
   dgesdd( jobz, m, n, A, ldA, S, U, ldU, Vt, ldVt, work, lwork, iwork, info );
#elif defined( USE_ACML )
   dgesdd_( (char*)jobz, (int*)m, (int*)n, A, (int*)ldA, S, U, (int*)ldU, Vt, (int*)ldVt, work, (int*)lwork, iwork, info, 1 );
#endif
   }
   
//symmetric eigenvectors/eigenvalues
inline void BlasLapack::syevr( const char *jobz, const char *range, const char *uplo, const int *n, double *A, const int *ldA, const double *vl, const double *vu, const int *il, const int *iu, const double *abstol, int *m, double *W, double *Z, const int *ldZ, int *isuppZ, double *work, const int *lwork, int *iwork, const int *liwork, int *info )
   {
#if defined( USE_MKL )
   dsyevr( jobz, range, uplo, n, A, ldA, vl, vu, il, iu, abstol, m, W, Z, ldZ, isuppZ, work, lwork, iwork, liwork, info );
#elif defined( USE_ACML )
   dsyevr_( (char*)jobz, (char*)range, (char*)uplo, (int*)n, A, (int*)ldA, (double*)vl, (double*)vu, (int*)il, (int*)iu, (double*)abstol, m, W, Z, (int*)ldZ, isuppZ ,work, (int*)lwork, iwork, (int*)liwork, info, 1, 1, 1 );
#endif
   }
   
//general matrix least squares solver
inline void BlasLapack::gelsy( const int *m, const int *n, const int *nrhs, double *A, const int *ldA, double *B, const int *ldB, int *jpvt, const double *rcond, int *rank, double *work, const int *lwork, int *info )
   {
#if defined( USE_MKL )
   dgelsy( m, n, nrhs, A, ldA, B, ldB, jpvt, rcond, rank, work, lwork, info );
#elif defined( USE_ACML )
   dgelsy_( (int*)m, (int*)n, (int*)nrhs, A, (int*)ldA, B, (int*)ldB, (int*)jpvt, (double*)rcond, rank, work, (int*)lwork, info );
#endif
   }

/*
 * VML:
 */
   
// res[:] := Ln(x[:])
inline void BlasLapack::vLn( const int n, const double *x, double *res )
   {
#if defined( USE_MKL )
   vdLn( n, x, res );
#elif defined( USE_ACML )
   for( int i=0; i<n; ++i )
      {
      res[i] = log( x[i] );
      }
   //vrda_log( n, (double*)x, res );
#endif
   }

// res[:] := Sqrt(x[:])
inline void BlasLapack::vSqrt( const int n, const double *x, double *res )
   {
#if defined( USE_MKL )
   vdSqrt( n, x, res );
#elif defined( USE_ACML )
   for ( int i=0; i<n; ++i )
      {
      res[i] = sqrt( x[i] );
      }
#endif
   }

#endif      //BlasLapack_h
