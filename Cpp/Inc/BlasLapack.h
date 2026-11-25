#if !defined( BlasLapack_h )
#define BlasLapack_h
/*
 *******************************************************************
 *
 *    Copyright (c) Microsoft. All rights reserved.
 *
 *    THIS CODE IS PROVIDED *AS IS* WITHOUT WARRANTY OF
 *    ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING ANY
 *    IMPLIED WARRANTIES OF FITNESS FOR A PARTICULAR
 *    PURPOSE, MERCHANTABILITY, OR NON-INFRINGEMENT.
 *
 ******************************************************************
 */

/*
 * This file sets up the appropriate math library to be used by FastLmmC
 *  Currently supported libraries are:
 *    Intel's Math Kernel Library (MKL) http://software.intel.com/en-us/articles/intel-mkl/
 *      or http://software.intel.com/sites/products/documentation/hpc/mkl/mklman/
 *    AMD's Core Math Library (ACML) http://www.amd.com/acml
 *
 *  You must select one of the libraries by defining USE_MKL or USE_ACML when 
 *    compiling FastLmmC
 */
#if defined( USE_MKL )
#define  BLAS_INT MKL_INT
#include "mkl.h"
#elif defined( USE_ACML )
#define  BLAS_INT int
#include "acml.h"
#else
#error USE_MKL or USE_ACML must be defined to set the appropriate math library usage
#endif

//global variables for use with BLAS:
const BLAS_INT int_one = 1;
const BLAS_INT int_minusone = -1;
const BLAS_INT int_zero = 0;
const real     real_one = 1.0;
const real     real_zero = 0.0;
const real     real_minusone = -1.0;
const char     char_V = 'V';        //output the eigenvalues plus eigenvectors(jobz)
const char     char_A = 'A';        //compute all eigenvectors (range)
const char     char_L = 'L';        //the lower part of symmetrix matrix A is being read(uplo)
const char     char_N = 'N';        //work with the untransposed matrix(trans)
const char     char_T = 'T';        // Work with the transpose of the matrix(trans)
const char     char_S = 'S';        //comute the economy size SVD.(jobz)

class BlasLapack
   {
public:
   /*
    * BLAS Level 1:
    */

   // vector scaling by a constant
   // x[::incx] := alpha * x[::incx]
   // x : [1+(n-1)*abs(incx)] vector
   // alpha real constant
   static void scal( const BLAS_INT *n, const real *alpha, real *x, const BLAS_INT *incx );

   // vector copy
   // Y[::incY] := X[::incX]
   // X : [1+(n-1)*abs(incX)] vector
   // Y : [1+(n-1)*abs(incY)] vector
   static void copy( const BLAS_INT *n, const real *X, const BLAS_INT *incX, real *Y, const BLAS_INT *incY );

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
   static void gemv( const char *trans, const BLAS_INT *m, const BLAS_INT *n, const real *alpha, const real *A, const BLAS_INT *ldA, const real *x, const BLAS_INT *incx, const real *beta, real *y, const BLAS_INT *incy );

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
   static void gemm( const char *transA,  const char *transB, const BLAS_INT *m, const BLAS_INT *n, const BLAS_INT *k, const real *alpha, const real *A, const BLAS_INT *ldA, const real *B, const BLAS_INT *ldB, const real *beta, real *C, const BLAS_INT *ldC );

   // symmetric rank k update. result will be in upper or lower storage format
   // C[:,:] := alpha * ( trans(A[:,:]) * trans(A[:,:]).T ) + beta * C[:,:]
   // uplo        : 'U' : upper storage , 'L' : lower storage
   // trans       : 'N' normal, 'T' transposed
   // A    for 'N': [ldA * k] matrix, ldA >= max( 1, n )
   //      for 'T': [ldA * n] matrix, ldA >= max( 1, k )
   static void syrk( const char *uplo, const char *trans, const BLAS_INT *n, const BLAS_INT *k, const real *alpha, real *A, const BLAS_INT *ldA, const real *beta, real *C, const BLAS_INT *ldC );

   /*
    * Lapack:
    */

   //singular value decomposition (SVD) of a [m * n] matrix A
   static void gesdd( const char *jobz, const BLAS_INT *m, const BLAS_INT *n, real *A, const BLAS_INT *ldA, real *S, real *U, const BLAS_INT *ldU, real *Vt, const BLAS_INT * ldVt, real *work, const BLAS_INT *lwork, BLAS_INT *iwork, BLAS_INT *info );
   
   //symmetric eigenvectors/eigenvalues
   static void syevr( const char *jobz, const char *range, const char *uplo, const BLAS_INT *n, real *A, const BLAS_INT *ldA, const real *vl, const real *vu, const BLAS_INT *il, const BLAS_INT *iu, const real *abstol, BLAS_INT *m, real *W, real *Z, const BLAS_INT *ldZ, BLAS_INT *isuppZ, real *work, const BLAS_INT *lwork, BLAS_INT *iwork, const BLAS_INT *liwork, BLAS_INT *info );
   
   //general matrix least squares solver
   static void gelsy( const BLAS_INT *m, const BLAS_INT *n, const BLAS_INT *nrhs, real *A, const BLAS_INT *ldA, real *B, const BLAS_INT *ldB, BLAS_INT *jpvt, const real *rcond, BLAS_INT *rank, real *work, const BLAS_INT *lwork, BLAS_INT *info );

   /*
    * VML:
    */
   
   // res[:] := Ln(x[:])
   static void vLn( const BLAS_INT n, const real *x, real *res );

   // res[:] := Sqrt(x[:])
   static void vSqrt( const BLAS_INT n, const real *x, real *res );

   // cholesky decompostion
   static void potrf(const char *uplo, const BLAS_INT *n, real *a, const BLAS_INT *lda, BLAS_INT *info);

   // solve the linear system using cholesky decompostion
   static void potrs(const char *uplo, const BLAS_INT *n, const BLAS_INT *nrhs, real *a, const BLAS_INT *lda, real *b, const BLAS_INT *ldb, BLAS_INT *info);

   };

/*
 * BLAS Level 1:
 */

// vector scaling by a constant
inline void BlasLapack::scal( const BLAS_INT *n, const real *alpha, real *x, const BLAS_INT *incx )
   {
#if defined( USE_MKL )

   #if defined( USE_FLOAT )
      sscal( n, alpha, x, incx );
   #else
      dscal( n, alpha, x, incx );
   #endif
   
#elif defined ( USE_ACML )
   #if defined( USE_FLOAT )
      sscal( *n, *alpha, x, *incx );
   #else
      dscal( *n, *alpha, x, *incx );
   #endif
#endif
   }

// vector copy
inline void BlasLapack::copy( const BLAS_INT *n, const real *X, const BLAS_INT *incX, real *Y, const BLAS_INT *incY )
   {
#if defined( USE_MKL )
   #if defined( USE_FLOAT )
      scopy( n, X, incX, Y ,incY );
   #else
      dcopy( n, X, incX, Y ,incY );
   #endif
#elif defined( USE_ACML )
   #if defined( USE_FLOAT )
      scopy( *n, (real*)X, *incX, Y, *incY );
   #else
      dcopy( *n, (real*)X, *incX, Y, *incY );
   #endif
#endif
   }

/*
 * BLAS Level 2:
 */
   
// general matrix vector product
//  gemv performs the matrix vector operation y := (alpha * A * x) + (beta * y)
//    alpha and beta are scalars (in)
//    x is an n-element vector (in)
//    A is an m-row x n-column matrix (in)
//    y is an n-element vector (in & out)
//
inline void BlasLapack::gemv( const char *trans, const BLAS_INT *m, const BLAS_INT *n, const real *alpha, const real *A, const BLAS_INT *ldA, const real *x, const BLAS_INT *incx, const real *beta, real *y, const BLAS_INT *incy )
   {
#if defined( USE_MKL )
   #if defined( USE_FLOAT )
      sgemv( trans, m, n, alpha, A, ldA, x, incx, beta, y, incy ); 
   #else
      dgemv( trans, m, n, alpha, A, ldA, x, incx, beta, y, incy ); 
   #endif
#elif defined( USE_ACML )
   #if defined( USE_FLOAT )
      sgemv( *trans, *m, *n, *alpha, (real*)A, *ldA, (real*)x, *incx, *beta, y, *incy );  
   #else
      dgemv( *trans, *m, *n, *alpha, (real*)A, *ldA, (real*)x, *incx, *beta, y, *incy ); 
   #endif
#endif
   }

/*
 * BLAS Level 3:
 */

// general matrix matrix product
//   gemm performs the matrix - matrix operation C = (alpha * op(A) * op(B)) + (beta * C)
//   alpha and beta are scalars (in)
//   A is an m row x k column matrix (in)
//   B is an k row x n column matrix (in)
//   C is an m row x n column matrix (in & out)
//
inline void BlasLapack::gemm( const char *transA,  const char *transB, const BLAS_INT *m, const BLAS_INT *n, const BLAS_INT *k, const real *alpha, const real *A, const BLAS_INT *ldA, const real *B, const BLAS_INT *ldB, const real *beta, real *C, const BLAS_INT *ldC )
   {
#if defined( USE_MKL )
   #if defined( USE_FLOAT )
      sgemm( transA, transB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC ); 
   #else
      dgemm( transA, transB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC );
   #endif
#elif defined( USE_ACML )
   #if defined( USE_FLOAT )
      sgemm( *transA, *transB, *m, *n, *k, *alpha, (real*)A, *ldA, (real*)B, *ldB, *beta, C, *ldC );
   #else
      dgemm( *transA, *transB, *m, *n, *k, *alpha, (real*)A, *ldA, (real*)B, *ldB, *beta, C, *ldC );
   #endif
#endif
   }

// symmetric rank k update. result will be in upper or lower storage format
inline void BlasLapack::syrk( const char *uplo, const char *trans, const BLAS_INT *n, const BLAS_INT *k, const real *alpha, real *A, const BLAS_INT *ldA, const real *beta, real *C, const BLAS_INT *ldC )
   {
#if defined( USE_MKL )
   #if defined( USE_FLOAT )
      ssyrk( uplo, trans, n, k, alpha, A, ldA, beta, C, ldC );
   #else
      dsyrk( uplo, trans, n, k, alpha, A, ldA, beta, C, ldC );
   #endif
#elif defined( USE_ACML )
   #if defined( USE_FLOAT )
      ssyrk( *uplo, *trans, *n, *k, *alpha, A, *ldA, *beta, C, *ldC );
   #else
      dsyrk( *uplo, *trans, *n, *k, *alpha, A, *ldA, *beta, C, *ldC );
   #endif
#endif
   }

/*
 * Lapack:
 */




//computes cholesky factorization of a symmetric positive-definite matrix
inline void BlasLapack::potrf(const char *uplo, const BLAS_INT *n, real *a, const BLAS_INT *lda, BLAS_INT *info)
   {
#if defined( USE_MKL )
   #if defined( USE_FLOAT )
      spotrf(uplo, n, a, lda, info);
   #else
      dpotrf(uplo, n, a, lda, info);
   #endif
#elif defined( USE_ACML )
   cout << "no ACML definition" << endl;
#endif
   }



//solve linear system using cholesky decomposition
inline void BlasLapack::potrs(const char *uplo, const BLAS_INT *n, const BLAS_INT *nrhs, real *a, const BLAS_INT *lda, real *b, const BLAS_INT *ldb, BLAS_INT *info)
   {
#if defined( USE_MKL )
   #if defined( USE_FLOAT )
      spotrs(uplo, n, nrhs, a, lda, b, ldb, info);
   #else
      dpotrs(uplo, n, nrhs, a, lda, b, ldb, info);
   #endif
#elif defined( USE_ACML )
   cout << "no ACML definition" << endl;
#endif
   }


//singular value decomposition (SVD) of a [m * n] matrix A
inline void BlasLapack::gesdd( const char *jobz, const BLAS_INT *m, const BLAS_INT *n, real *A, const BLAS_INT *ldA, real *S, real *U, const BLAS_INT *ldU, real *Vt, const BLAS_INT * ldVt, real *work, const BLAS_INT *lwork, BLAS_INT *iwork, BLAS_INT *info )
   {
#if defined( USE_MKL )
   #if defined( USE_FLOAT )
      sgesdd( jobz, m, n, A, ldA, S, U, ldU, Vt, ldVt, work, lwork, iwork, info );
   #else
      dgesdd( jobz, m, n, A, ldA, S, U, ldU, Vt, ldVt, work, lwork, iwork, info );
   #endif
#elif defined( USE_ACML )
   #if defined( USE_FLOAT )
      sgesdd_( (char*)jobz, (BLAS_INT*)m, (BLAS_INT*)n, A, (BLAS_INT*)ldA, S, U, (BLAS_INT*)ldU, Vt, (BLAS_INT*)ldVt, work, (BLAS_INT*)lwork, iwork, info, 1 );
   #else
      dgesdd_( (char*)jobz, (BLAS_INT*)m, (BLAS_INT*)n, A, (BLAS_INT*)ldA, S, U, (BLAS_INT*)ldU, Vt, (BLAS_INT*)ldVt, work, (BLAS_INT*)lwork, iwork, info, 1 );
   #endif
#endif
   }
   
//symmetric eigenvectors/eigenvalues
inline void BlasLapack::syevr( const char *jobz, const char *range, const char *uplo, const BLAS_INT *n, real *A, const BLAS_INT *ldA, const real *vl, const real *vu, const BLAS_INT *il, const BLAS_INT *iu, const real *abstol, BLAS_INT *m, real *W, real *Z, const BLAS_INT *ldZ, BLAS_INT *isuppZ, real *work, const BLAS_INT *lwork, BLAS_INT *iwork, const BLAS_INT *liwork, BLAS_INT *info )
   {
#if defined( USE_MKL )
   #if defined( USE_FLOAT )
      ssyevr( jobz, range, uplo, n, A, ldA, vl, vu, il, iu, abstol, m, W, Z, ldZ, isuppZ, work, lwork, iwork, liwork, info );
   #else
      dsyevr( jobz, range, uplo, n, A, ldA, vl, vu, il, iu, abstol, m, W, Z, ldZ, isuppZ, work, lwork, iwork, liwork, info );
   #endif
#elif defined( USE_ACML )
   #if defined( USE_FLOAT )
      ssyevr_( (char*)jobz, (char*)range, (char*)uplo, (BLAS_INT*)n, A, (BLAS_INT*)ldA, (real*)vl, (real*)vu, (BLAS_INT*)il, (BLAS_INT*)iu, (real*)abstol, m, W, Z, (BLAS_INT*)ldZ, isuppZ ,work, (BLAS_INT*)lwork, iwork, (BLAS_INT*)liwork, info, 1, 1, 1 );
   #else
      dsyevr_( (char*)jobz, (char*)range, (char*)uplo, (BLAS_INT*)n, A, (BLAS_INT*)ldA, (real*)vl, (real*)vu, (BLAS_INT*)il, (BLAS_INT*)iu, (real*)abstol, m, W, Z, (BLAS_INT*)ldZ, isuppZ ,work, (BLAS_INT*)lwork, iwork, (BLAS_INT*)liwork, info, 1, 1, 1 );
   #endif
#endif
   }
   
//general matrix least squares solver
inline void BlasLapack::gelsy( const BLAS_INT *m, const BLAS_INT *n, const BLAS_INT *nrhs, real *A, const BLAS_INT *ldA, real *B, const BLAS_INT *ldB, BLAS_INT *jpvt, const real *rcond, BLAS_INT *rank, real *work, const BLAS_INT *lwork, BLAS_INT *info )
   {
#if defined( USE_MKL )
   #if defined( USE_FLOAT )
      sgelsy( m, n, nrhs, A, ldA, B, ldB, jpvt, rcond, rank, work, lwork, info );
   #else
      dgelsy( m, n, nrhs, A, ldA, B, ldB, jpvt, rcond, rank, work, lwork, info );
   #endif
#elif defined( USE_ACML )
   #if defined( USE_FLOAT )
      sgelsy_( (BLAS_INT*)m, (BLAS_INT*)n, (BLAS_INT*)nrhs, A, (BLAS_INT*)ldA, B, (BLAS_INT*)ldB, (BLAS_INT*)jpvt, (real*)rcond, rank, work, (BLAS_INT*)lwork, info );
   #else
      dgelsy_( (BLAS_INT*)m, (BLAS_INT*)n, (BLAS_INT*)nrhs, A, (BLAS_INT*)ldA, B, (BLAS_INT*)ldB, (BLAS_INT*)jpvt, (real*)rcond, rank, work, (BLAS_INT*)lwork, info );
   #endif
#endif
   }

/*
 * VML:
 */
   
// res[:] := Ln(x[:])
inline void BlasLapack::vLn( const BLAS_INT n, const real *x, real *res )
   {
#if defined( USE_MKL )
   #if defined( USE_FLOAT )
      vsLn( n, x, res );
   #else
      vdLn( n, x, res );
   #endif
#elif defined( USE_ACML )
   for( BLAS_INT i=0; i<n; ++i )
      {
      res[i] = log( x[i] );
      }
   //vrda_log( n, (real*)x, res );
#endif
   }

// res[:] := Sqrt(x[:])
inline void BlasLapack::vSqrt( const BLAS_INT n, const real *x, real *res )
   {
#if defined( USE_MKL )
   #if defined( USE_FLOAT )
      vsSqrt( n, x, res );
   #else
      vdSqrt( n, x, res );
   #endif
#elif defined( USE_ACML )
   for ( BLAS_INT i=0; i<n; ++i )
      {
      res[i] = sqrt( x[i] );
      }
#endif
   }

/*
 *  Set the MKL library to provide more 'reproducible' floating point results
 *    This reduces performance by limiting MKL to one thread and minmizing
 *    the dynamic dispatch on CPU type.
 */
inline void SetReproducibleFP( bool fDoReproducibleFP )
   {
#if defined( USE_MKL )
#if defined( USE_LINUX_V2x_KERNEL )
/*
 *  TODO:
 *  The MKL v10.x we are using on the Linux 2.x kernel does not support
 *    the MKL v11.x reproducible math functions.  
 *    RedHat RHEL 6 is still using a 2.6 kernel and gcc 4.4 toolset.
 *    Enable SetReproducibleFP() on RHEL after RHEL moves to a 3.x kernel
 *    and the MKL is updated to v11.x+
 */
#else
   // set MKL to us 'reproducible' paths and instructions and restrict thread usage
   // see http://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_userguide_win/index.htm#GUID-497CAEBE-9D78-4C3E-90FE-78437D73C038.htm
   if ( fDoReproducibleFP )
      {
      mkl_cbwr_set( MKL_CBWR_COMPATIBLE );
      mkl_set_num_threads( 1 );
      mkl_set_dynamic( false );     // integer 0 == false
      }
   else
      {
      mkl_cbwr_set( MKL_CBWR_AUTO );
      mkl_set_dynamic( true );      // integer true == 1
      }
#endif   // USE_LINUX_V2x_KERNEL
#endif   // USE_MKL
   }

inline string MathLibName()
   {
   [[noreturn]] void Fatal( const char *szFmt, ... );
   void Warn( const char *szFmt, ... );
#if defined( USE_MKL )
   char szMklVersion[ 64 ];
   MKLVersion version;
   
   mkl_get_version( &version );
   sprintf( szMklVersion, "MKL v%d.%02d.%02d - Build: %s", 
            version.MajorVersion, 
            version.MinorVersion, 
            version.UpdateVersion, 
            version.Build );

#if defined( USE_LINUX_V2x_KERNEL )
   if ( (version.MajorVersion != 10)
     || (version.MinorVersion != 3)
     || (version.UpdateVersion != 7)
     || (_strcmpi( version.Build, "20111003" ) != 0)
      )
      {
      Warn( "Unexpected MKL Version Found.  "
          "\n          Expecting: MKL v10.03.07 - Build: 20111003"
          "\n              Found: %s", szMklVersion );
      if (getenv( "FastLmmUseAnyMklLib" ) == nullptr )
         {
         Fatal( "Exit due to MKL version mismatch"
              "\n              To bypass this check, Set FastLmmUseAnyMklLib=1" );
         }
      else
         {
         Warn( "FastLmmUseAnyMklLib set.  Continuing execution.\n\n" );
         }
      }
#else
   // Pinned to the oneMKL release whose results were validated against the
   //   FaSTLMM 2.07c MKL v11.00.04 (Build 20130517) binary.  oneMKL 2026.2.1
   //   (PyPI mkl-devel/mkl-static) reports itself as v2026.00.02.
   if ( (version.MajorVersion != 2026)
     || (version.MinorVersion != 0)
     || (version.UpdateVersion != 2)
     || (_strcmpi( version.Build, "20260915" ) != 0)
      )
      {
      Warn( "Unexpected MKL Version Found.  "
          "\n          Expecting: MKL v2026.00.02 - Build: 20260915"
          "\n              Found: %s", szMklVersion );
      if (getenv( "FastLmmUseAnyMklLib" ) == nullptr )
         {
         Fatal( "Exit due to MKL version mismatch"
              "\n              To bypass this check, Set FastLmmUseAnyMklLib=1" );
         }
      else
         {
         Warn( "FastLmmUseAnyMklLib set.  Continuing execution.\n\n" );
         }
      }
#endif

   return( szMklVersion );

#elif defined( USE_ACML )
   char szAcmlVersion[ 32 ];
   int  acmlMajor;
   int  acmlMinor;
   int  acmlPatch;

   acmlversion( &acmlMajor, &acmlMinor, &acmlPatch );
   sprintf( szAcmlVersion, "ACML v%d.%02d.%02d", acmlMajor, acmlMinor, acmlPatch );
   if ( (acmlMajor != 4)
     || (acmlMinor != 4)
     || (acmlPatch != 0)
      )
      {
      Warn( "Unexpected ACML Version Found.  "
          "\n          Expecting: ACML v4.04.00"
          "\n              Found: %s", szAcmlVersion );
      if (getenv( "FastLmmUseAnyAcmlLib" ) == nullptr )
         {
         Fatal( "Exit due to version mismatch"
              "\n              To bypass this check, Set FastLmmUseAnyAcmlLib=1" );
         }
      else
         {
         Warn( "FastLmmUseAnyAcmlLib set.  Continuing execution.\n\n" );
         }
      }
   return( szAcmlVersion );
#endif
   }

#endif      //BlasLapack_h
