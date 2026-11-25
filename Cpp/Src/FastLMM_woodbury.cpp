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

#include "FastLmmC.h"

//constructor takes an EigenSym object as input
FastLMM_woodbury::FastLMM_woodbury( EigenSym *eigen ) : FastLMM( eigen )
   {
   if ( this->eigen->n <= this->eigen->rankU )
      {
      Fatal( "Eigen dimensions do not agree with FastLMM_woodbury. rankU=%i, n=%i", eigen->rankU, eigen->n );
      }

   this->rankU = eigen->rankU;
   SU_y = nullptr;         // preset all the pointers to avoid garbage
   y    = nullptr;
   SU_x = nullptr;
   X    = nullptr;
   SUcovar = nullptr;
   covar   = nullptr;

   modelID = "FastLMM_woodbury";
   }

FastLMM_woodbury::~FastLMM_woodbury()
   {
   if ( SU_y )    { delete[] SU_y;    SU_y    = nullptr; }      
   if ( SU_x )    { delete[] SU_x;    SU_x    = nullptr; }        
   if ( SUcovar ) { delete[] SUcovar; SUcovar = nullptr; }

   //TODO: REFACTOR
   // covar, X, and y are really cached pointers to other 
   //   people's data.  A shared pointer would be better
   //   here but don't delete them now
   //if ( covar )   { delete[] covar;   covar   = nullptr; }
   //if ( X )       { delete[] X;       X       = nullptr; }
   //if ( y )       { delete[] y;       y       = nullptr; }
   }

//SNPs [n by nSNPs] holds the nSNPs SNPs to be tested for each of n individuals.
void FastLMM_woodbury::set_SNPs( real *SNPs, size_t nSNPs_, vector< SnpInfo >& snpInfo_ )
   {
   this->snpInfo = snpInfo_;
   this->nSNPs = nSNPs_;
   if ( this->nSNPs != this->snpInfo.size() )
      {
      Fatal( "SNP counts do not match.  nSNPs = %i, snpInfo.size() = %i", this->nSNPs, this->snpInfo.size() );
      }

   if ( SU_x != nullptr ) { delete[] SU_x; }
//   if ( X    != nullptr ) { delete[] X; }

   SU_x = new real[ nSNPs * rankU ];       //as number of SNPs might change, always reallocate memory for UX.
   BLAS_INT n_int = (BLAS_INT)n;
   BLAS_INT rankU_int = (BLAS_INT)rankU;
   BLAS_INT nSNPs_int = (BLAS_INT)nSNPs;

   // SU_x = U.T * SNPs
   BlasLapack::gemm( &char_T, &char_N, &rankU_int, &nSNPs_int, &n_int, &real_one, eigen->U, &n_int, SNPs, &n_int, &real_zero, SU_x, &rankU_int );
   this->X = SNPs;
   }

// phenotype is a vector of length n, holding the (single) phenotype.
void FastLMM_woodbury::set_Phenotype( real *phenotype, const string& phenLabel_ )
   {
   if ( SU_y != 0 ) { delete[] SU_y; }
//   if ( y    != 0 ) { delete[] y; }
   phenLabel = phenLabel_;

   SU_y = new real[ rankU ];              // The phenotypes transformed to low rank
   BLAS_INT n_int = (BLAS_INT)n;
   BLAS_INT rankU_int = (BLAS_INT)rankU;

   // Uy = U.T * phenotype
   BlasLapack::gemv( &char_T, &n_int, &rankU_int, &real_one, eigen->U, &n_int, phenotype, &int_one, &real_zero, SU_y, &int_one );
   this->y = phenotype;
   }

// covariates is an [n*ncovar] matrix holding the covariates, including the bias
void FastLMM_woodbury::set_covariates( real *covar, size_t ncovariates )
   {
   ncovar = ncovariates;                  // the number of covariates (including the fixed effects)
   if ( SUcovar != nullptr )     { delete[] SUcovar; }
//   if ( this->covar != nullptr ) { delete[] this->covar; }

   this->SUcovar = new real[ ncovar * rankU ];     // The covariates transformed to the low rank space
   BLAS_INT n_int = (BLAS_INT)n;
   BLAS_INT rankU_int  = (BLAS_INT)rankU;
   BLAS_INT ncovar_int = (BLAS_INT)ncovar;

   // UX = U.T * SNPs
   BlasLapack::gemm( &char_T, &char_N, &rankU_int, &ncovar_int, &n_int, &real_one, eigen->U, &n_int, covar, &n_int, &real_zero, SUcovar, &rankU_int );
   this->covar = covar;
   }

//find ML parameters for the null model
void FastLMM_woodbury::trainNull()
   {
   real SU_xdummy = 0.0;         // As We don't have a SNP, we pass a dummy variable

   // find optimal delta:
   if ( !_isnan( logDeltaPreset ) )
      {
      delta_0[ 0 ] = logDeltaPreset;
      }
   else
      {
      if ( FastLMM::brentExhaustive )
         {
         // Find the optimal delta
//         delta_0[ 0 ] = FastLMM_woodbury::optimize_delta_Brentintervals( 0, (BLAS_INT)ncovar, (BLAS_INT)n, (BLAS_INT)rankU, SU_y, &SU_xdummy, y, &SU_xdummy, eigen->S, SUcovar, covar );
         delta_0[ 0 ] = FastLMM_woodbury::optimize_delta_Brentintervals( 0, (BLAS_INT)ncovar, (BLAS_INT)n, (BLAS_INT)rankU, SU_y, &SU_xdummy, y, &SU_xdummy, eigen->S, SUcovar, covar, numberOptimizationStartsGlobal, minLogVal, maxLogVal, brentTol );
         }
      else
         {
         // Find the optimal delta
         delta_0[ 0 ] = FastLMM_woodbury::optimize_delta_BrentintervalsCheckTriplets( 0, (BLAS_INT)ncovar, (BLAS_INT)n, (BLAS_INT)rankU, SU_y, &SU_xdummy, y, &SU_xdummy, eigen->S, SUcovar, covar );
         }
      }

   // evaluate the LL at delta_0 and get the ML parameters for the optimal delta
   LL_0[ 0 ] = (real)-1.0 * nLLeval( nullptr, delta_0[0], 0, (BLAS_INT)ncovar, (BLAS_INT)n, (BLAS_INT)rankU, SU_y, &SU_xdummy, y, &SU_xdummy, eigen->S, SUcovar, covar, &sigg2_0[0], &sige2_0[0], &beta_0[0], &rank_beta_0[0], true );

   if ( verboseLevel >= verboseDebug )
      {
      printf( "\n------------------------------------------" );
      printf( "\nNULL: LL:%f", LL_0[ 0 ] );
      printf( "\n  sigg2=%f sige2=%f delta=%f", sigg2_0[ 0 ], sige2_0[ 0 ], delta_0[ 0 ] );
      printf(" beta: ");
      for( size_t i=0; i < ncovar; ++i )
         {
         printf( "%f ",beta_0[ i ] );
         }
      printf("\n-------------------------------------------" );
   }

   // copy the "default" null model parameters
   for ( size_t snp = 1; snp < nSNPs; ++snp )
      {
      delta_0[ snp ] = delta_0[ 0 ];
      sigg2_0[ snp ] = sigg2_0[ 0 ];
      sige2_0[ snp ] = sige2_0[ 0 ];
      rank_beta_0[ snp ] = rank_beta_0[ 0 ];
      LL_0[ snp ] = LL_0[ 0 ];
      for ( size_t i_covar = 0; i_covar < ncovar; ++i_covar )
         {
         beta_0[ (snp * ncovar) + i_covar ] = beta_0[ i_covar ];
         }
      }
   }

// find ML parameters for all SNPs
void FastLMM_woodbury::trainAlt()
   {
   // alternative models:
   size_t d = 1;                // assume that always one SNP is being tested
   size_t d_ = d + ncovar;      // the number of weights (dimensionality of a single beta)

   #pragma omp parallel for
   for ( int iSnp = 0; iSnp < (int)nSNPs; ++iSnp )
      {
      // train the alternative model for the SNP     (parallel for loop)
      // set pointers to the current SNP and ML parameters:
      real* currentSU_x  = SU_x + ( iSnp * ( d * rankU ) );
      real* currentX     = X + ( iSnp * ( n * d ) );
      real* currentBeta  = beta + ( iSnp * d_ );
      real* currentSigg2 = sigg2 + iSnp;
      real* currentSige2 = sige2 + iSnp;
      real* currentDelta = delta + iSnp;
      BLAS_INT* currentRank_beta = rank_beta + iSnp;

      // find delta:
      if ( covLearnType == FULL )
         {
         // Learn delta
         if ( FastLMM::brentExhaustive )
            {
            // Find the optimal delta method 1
//            *currentDelta = FastLMM_woodbury::optimize_delta_Brentintervals( (BLAS_INT)d, (BLAS_INT)ncovar, (BLAS_INT)n, (BLAS_INT)rankU, SU_y, currentSU_x, y, currentX, eigen->S, SUcovar, covar );
            *currentDelta = FastLMM_woodbury::optimize_delta_Brentintervals( (BLAS_INT)d, (BLAS_INT)ncovar, (BLAS_INT)n, (BLAS_INT)rankU, SU_y, currentSU_x, y, currentX, eigen->S, SUcovar, covar, numberOptimizationStartsGlobal, minLogVal, maxLogVal, brentTol );
            }
         else
            {
            // Find the optimal delta method 2
            *currentDelta = FastLMM_woodbury::optimize_delta_BrentintervalsCheckTriplets( (BLAS_INT)d, (BLAS_INT)ncovar, (BLAS_INT)n, (BLAS_INT)rankU, SU_y, currentSU_x, y, currentX, eigen->S, SUcovar, covar );
            }
         }
      else if ( covLearnType == ONCE )
         {
         // Or: use the lock-in. So it equals delta_0
         *currentDelta = delta_0[ iSnp ];
         }

      // evaluate the LL at delta_0 and get the ML parameters
      LL[ iSnp ] = (real)-1.0 * nLLeval( f_stat+(iSnp*d_), *currentDelta, (BLAS_INT)d, (BLAS_INT)ncovar, (BLAS_INT)n, (BLAS_INT)rankU, SU_y, currentSU_x, y, currentX, eigen->S, SUcovar, covar, currentSigg2 , currentSige2, currentBeta, currentRank_beta, true );

      // Compute p-value from LR test:
      real df = *currentRank_beta - (real)rank_beta_0[ iSnp ];
      if ( df > 0.0 )
         {
         // check the degrees of freedom for validity
         if ( use_f_test )
            {
            //Fishers F-test
            //
            // PvalUsingFisherFTestCdf focuses on keeping precision 
            //    for very small pvalues by inverting on the reverse 
            //    path of the Beta::BetaInc() branch
            pvals[ iSnp ] = PvalUsingFisherFTestCdf( f_stat[ iSnp * d_ ], df, (real)(n - *currentRank_beta) );
            //pvals[ iSnp ] = 1.0 - FisherF::Cdf( f_stat[ iSnp * d_ ], df , (real)(n - *currentRank_beta) );
            }
         else
            {
            // Likelihood-ratio test
            pvals[ iSnp ] = PvalFromLikelihoodRatioTest( LL[ iSnp ] - LL_0[ iSnp ], (real)0.5 * df );
            }
         }
      else
         {
         // This happens for example if the SNP is a constant. In this case df = 0
         pvals[ iSnp ] = std::numeric_limits< real >::quiet_NaN();
         }
      if ( verboseLevel >= verboseDebug )
         {
         printf( "\n- LL[%i]:%f p=%f",iSnp,LL[ iSnp ], pvals[ iSnp ] );
         printf( "\n  sigg2=%f sige2=%f delta=%f", *currentSigg2, *currentSige2, *currentDelta );
         printf( " beta: ");
         for( size_t i = 0; i < d_; ++i )
            {
            printf( "%f ", currentBeta[ i ] );
            }
         }
      }
   }

// static function:
// finds delta by optimizing over delta intervals
real FastLMM_woodbury::optimize_delta_BrentintervalsCheckTriplets( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU, const real *SU_y, const real *currentSU_x, const real *y, const real *currentX, const real *S, const real *SUcovar, const real *covar )
   {
   real *negLL_interval = new real[ numberIntervals + 1 ];
   real *log_delta_interval = new real[ numberIntervals + 1 ];

   real best_log_delta = std::numeric_limits< real >::quiet_NaN();
   real best_negLL = std::numeric_limits< real >::infinity();

   real MLparamsdummy = 0.0;
   BLAS_INT rankDummy = 0;
   size_t nevals = 0;

   // evaluate the boundaries of all intervals
   #pragma omp parallel for
   for ( int i = 0; i <= (int)this->numberIntervals; ++i )
      {
      log_delta_interval[ i ] = minLogVal + (((real)i * (maxLogVal - minLogVal)) / numberIntervals);
      negLL_interval[ i ] = FastLMM_woodbury::nLLeval( nullptr, log_delta_interval[ i ], d, ncovar, n, rankU, SU_y, currentSU_x, y, currentX, S, SUcovar, covar, &MLparamsdummy, &MLparamsdummy, &MLparamsdummy, &rankDummy, false );
      ++nevals;
      if ( negLL_interval[ i ] <= best_negLL )
         {
         best_negLL = negLL_interval[ i ];
         best_log_delta = log_delta_interval[ i ];
         }
      }

   brentFunctor_FastLMM_woodbury func = brentFunctor_FastLMM_woodbury( d, ncovar, n, rankU, SU_y, currentSU_x, y, currentX, S, SUcovar, covar );
   #pragma omp parallel for //parallelization might be slow here, as most threads finish immediately
   for ( int i = 1; i < (int)this->numberIntervals; ++i )
      {
      if ( (negLL_interval[ i ] < negLL_interval[ i + 1 ])
        && (negLL_interval[ i ] < negLL_interval[ i - 1 ])
         )
         {
         // check whether a local optimum exists in this triplet
         // If there is a local optimum, optimize over the current triplet:
         real current_negLL = std::numeric_limits< real >::infinity();
         size_t nevals_Int = 0;
         real current_log_delta = BrentC::minimize( func, log_delta_interval[ i - 1 ], log_delta_interval[ i + 1 ], brentTol, current_negLL, nevals_Int, brentMaxIter, true );
         nevals += nevals_Int;
         if ( current_negLL <= best_negLL )
            {
            best_negLL = current_negLL;
            best_log_delta = current_log_delta;
            }
         }
      }

   delete [] log_delta_interval;
   delete [] negLL_interval;

#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
      const char *szFmt1 = "\nnevals: %Iu logdelta: %f\n";
#else                      // Linux/g++ uses a %zu specifier for size_t
      const char *szFmt1 = "\nnevals: %zu logdelta: %f\n";
#endif
   Verbose( verboseDebug, szFmt1, nevals, best_log_delta );
   return( best_log_delta );
   }

// static function:
// finds delta by optimizing over log delta intervals
real FastLMM_woodbury::optimize_delta_Brentintervals( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU, const real *SU_y, const real *currentSU_x, const real *y, const real *currentX, const real *S, const real *SUcovar, const real *covar, const int numberIntervals_, const real minLogVal_, const real maxLogVal_, const real brentTol )
   {
   real best_log_delta = std::numeric_limits< real >::quiet_NaN();
   real best_negLL = std::numeric_limits< real >::infinity();

   // create the Functor. This will be passed to the optmizer
   brentFunctor_FastLMM_woodbury func = brentFunctor_FastLMM_woodbury( d, ncovar, n, rankU, SU_y, currentSU_x, y, currentX, S, SUcovar, covar );

   size_t nevals = 0;      // count the total number of LL-function evaluations

   // Perform a Brent optimization on each interval
   #pragma omp parallel for
   for ( int i = 0; i < numberIntervals_; ++i )
      {
      // Setup min and max of current Brent interval
      real current_log_delta_min = minLogVal_ + i * ( maxLogVal_ - minLogVal_ ) / numberIntervals_;
      real current_log_delta_max = minLogVal_ + ( i + 1 ) * ( maxLogVal_ - minLogVal_ ) / numberIntervals_;
      real current_negLL = std::numeric_limits< real >::infinity();  // This will hold the optimal negative Log Likelihood
      size_t nevals_Int = 0;                                         // The number of evaluations on the current interval
      real current_log_delta = BrentC::minimize( func, current_log_delta_min, current_log_delta_max, this->brentTol, current_negLL, nevals_Int, this->brentMaxIter, true );

      nevals += nevals_Int;            // Increment the total number of evaluations by #evaluations on the current one

      if ( current_negLL <= best_negLL )
         { 
         // Check, whether we have found a better optimum on the current interval
         best_negLL = current_negLL;       // If Yes:   Set the optimum to the current one
         best_log_delta = current_log_delta;
         }
      }

#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
      const char *szFmt1 = "\nnevals: %Iu logdelta: %f\n";
#else                      // Linux/g++ uses a %zu specifier for size_t
      const char *szFmt1 = "\nnevals: %zu logdelta: %f\n";
#endif
   Verbose( verboseDebug, szFmt1, nevals, best_log_delta );
   return( best_log_delta );  //The Brent optimization is in log space. We want the actual delta. => exp()
   }


// static function:
// cleaned up version of nLLeval() for the Woodbury model
real FastLMM_woodbury::nLLeval( real *f_statistic, const real logdelta, const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU, const real *SU_y, const real *SU_x, const real *y, const real *X, const real *S, const real *SUcovar, const real *covar, real *sigg2, real *sige2, real *beta, BLAS_INT *rank_beta, const bool returnMLparams )
   {
   real delta;
   if ( useHeritability )
      {
      delta = (1.0 - logdelta) / logdelta;
      }
   else
      {
      delta = exp( logdelta );
      }

   BLAS_INT d_ = d + ncovar;           // the number of parameters to be estimated is number of SNP coefficients + number of covariates
   real *Sd = new real[ rankU ];       // holds the sqrt of the diagonal values of (S+delta*I) (deleted)
   real *ySdi = new real[ rankU ];     // the residual vecttor. Will be initialized to ySdi, later overwritten by the residual of the regression (deleted)
   real *beta_temp = new real[ d_ ];   // the parameters. Will be initialized to ySdi (deleted)
   real *XSdi = new real[ ((size_t)rankU) * ((size_t)d_) ];    // XSdi This one will be used to compute the residual in the end. (deleted)
   size_t nd = ((size_t)n) * ((size_t)d);
   bool delete_f_stat = false;

   if ( f_statistic == nullptr )
      {
      f_statistic = new double[ d_ ];
      delete_f_stat = true;
      }

   size_t nncovar = ((size_t)n) * ((size_t)ncovar);
   real *Xcovar = new real[ nd + nncovar ];        //Xcovar is the concatenation of X and covar

   //concatenate X and covar in Xcovar
   for ( size_t i = 0; i < nd; i += n )
      {
      BlasLapack::copy( &n, X + i, &int_one, Xcovar + i, &int_one );    //copy the SNPs to Xcovar
      }

   real *Xcovar_ind = Xcovar + nd;
   for ( size_t i = 0; i < nncovar; i += n )
      {
      BlasLapack::copy( &n, covar + i, &int_one, Xcovar_ind + i, &int_one );  //copy the covariates to Xcovar
      }

   real yy = 0.0;                // yy is one of the three squared form variables yy, XX, and Xy
   for ( size_t ind = 0; (BLAS_INT)ind < n; ++ind )
      {
      yy += y[ ind ] * y[ ind ];
      }

   //--------------------------------------------------
   // compute the determinant and ySdi
   real det = (n - rankU) * log(delta);     // the log determinant is the sum of the diagonal of (S + delta * I)

   for ( size_t i = 0; (BLAS_INT)i < rankU; ++i )
      {
      // add delta to S
      Sd[ i ] = S[ i ] + delta;
      }

   BlasLapack::vLn( rankU, Sd, ySdi );    // res holds the Ln( S + deltaI )
   BlasLapack::vSqrt( rankU, Sd, Sd );    // square root of each Sdi element

   for ( size_t i = 0; (BLAS_INT)i < rankU; ++i )
      {
      det += ySdi[ i ];                   // compute the determinant
      ySdi[ i ] = SU_y[ i ] / Sd[ i ];    // scale the transformed ySdi
      yy -= ( ySdi[ i ] * ySdi[ i ] );
      }

   size_t XSdi_ind = 0;
   for ( size_t snp = 0; (BLAS_INT)snp < d; ++snp )
      {
      for ( size_t ind = 0; (BLAS_INT)ind < rankU; ++ind )
         {
         // Scale the transformed X
         XSdi[ XSdi_ind ] = SU_x[ XSdi_ind ] / Sd[ ind ];
         ++XSdi_ind;
         }
      }

   const real* SUcovar_ind = SUcovar;
   for ( size_t cov = 0; (BLAS_INT)cov < ncovar; ++cov )
      {
      for ( size_t ind = 0; (BLAS_INT)ind < rankU; ++ind )
         {
         // scale the transformed covariates
         XSdi[ XSdi_ind ] = *SUcovar_ind / Sd[ ind ];
         ++XSdi_ind;
         ++SUcovar_ind;
         }
      }

   delete[] Sd;  // Sdi is no longer needed

   real *XX = new real[ d_ * d_ ];     // compute XX = (Xcovar.T * Xcovar) - (XSdi.T * XSdi)
   if ( d_ > 0 )
      {
      BlasLapack::syrk( &char_L, &char_T, &d_, &n, &real_one, Xcovar, &n, &real_zero, XX, &d_ );
      BlasLapack::syrk( &char_L, &char_T, &d_, &rankU, &real_minusone, XSdi, &rankU, &real_one, XX, &d_ );
      }

   real *Xy = new real[ d_ ];          // compute Xy = Xcovar.T*y - XSdi.T*y
   BlasLapack::gemv( &char_T, &n, &d_, &real_one, Xcovar, &n, y, &int_one, &real_zero, Xy, &int_one );
   BlasLapack::gemv( &char_T, &rankU, &d_, &real_minusone, XSdi, &rankU, ySdi, &int_one, &real_one, Xy, &int_one );

   delete[] XSdi;          // XSdi is no longer needed.
   delete[] Xcovar;        // Xcovar is no longer needed
   delete[] ySdi;          // ySdi is no longer needed

   // This part should now be the same for all branches, except the sige2/sigg2 issue:

   // solve the least squares problem.
   real logDetXX = 0.0;
   FastLMM::solveLstSqSym( d_, XX, Xy, beta_temp, rank_beta, f_statistic, &logDetXX );

   real sige2_temp = yy;   //the ML sigma_e^2
   for ( size_t dim = 0; (BLAS_INT)dim < d_; ++dim )
      {
      sige2_temp -= (Xy[ dim ] * beta_temp[ dim ]);
      }

   real negLL;    // negative Log Likelihood Temp variable
   real sigg2_temp;

   if ( mlType == ML )
      {
      // compute the negative LL using the ML parameters at delta
      sige2_temp /= n;
      sigg2_temp = sige2_temp / delta;
      negLL = (real)0.5 * (( n * (_log2pi + log( sigg2_temp ) + (real)1.0) ) + det);
      }
   else if ( mlType == REML )
      {
      // compute the negative LL using the REML parameters at delta
      Warn( "The -woodbury lowrank formula for REML has to be verified. In the meantime please use -nowoodbury" );
      sige2_temp /= n-d_;
      sigg2_temp = sige2_temp / delta;
      negLL = (real)0.5 * ( (n - d_) * (_log2pi + log( sigg2_temp ) + (real)1.0) + det + logDetXX - (d_ * log(delta)) );
      }
   else
      {
      Fatal( "Only ML or REML are valid ML-types" );
      }

   if ( returnMLparams )
      {
      // only return the ML params if needed. (e.g. no need during Brents search) This avoids danger of Memory errors
      *sige2 = sige2_temp;
      *sigg2 = sigg2_temp;
      BlasLapack::copy( &d_, beta_temp, &int_one, beta, &int_one );
      }

   delete[] beta_temp;  //delete the ML beta_temp
   delete[] XX;
   delete[] Xy;

   if ( delete_f_stat )
      {
      delete[] f_statistic;
      }
   else
      {
      for ( size_t dim = 0; (BLAS_INT)dim < d_; ++dim )
         {
         f_statistic[ dim ] /= sige2_temp;
         }
      }

   //printf("negLL:%f    sigg2:%f    sige2:%f\n",negLL,sigg2_temp,sigg2_temp*delta);
   return( negLL );
   }

void FastLMM_woodbury::init()
   {
   if ( beta_0 != nullptr )
      {
      delete[] beta_0;
      beta_0 = nullptr;
      }
   if ( sigg2_0 != nullptr )
      {
      delete[] sigg2_0;
      sigg2_0 = nullptr;
      }
   if ( sige2_0 != nullptr )
      {
      delete[] sige2_0;
      sige2_0 = nullptr;
      }
   if ( delta_0 != nullptr )
      {
      delete[] delta_0;
      delta_0 = nullptr;
      }
   if ( rank_beta_0 != nullptr )
      {
      delete[] rank_beta_0;
      rank_beta_0 = nullptr;
      }
   if ( LL_0 != nullptr )
      {
      delete[] LL_0;
      }

   if ( nSNPs > 0 )
      {
      LL_0    = new real[ nSNPs ];
      beta_0  = new real[ ncovar * nSNPs ];
      sigg2_0 = new real[ nSNPs ];
      sige2_0 = new real[ nSNPs ];
      delta_0 = new real[ nSNPs ];
      rank_beta_0 = new BLAS_INT[ nSNPs ];
      }
   else
      {
      LL_0    = new real[ 1 ];
      beta_0  = new real[ ncovar ];
      sigg2_0 = new real[ 1 ];
      sige2_0 = new real[ 1 ];
      delta_0 = new real[ 1 ];
      rank_beta_0 = new BLAS_INT[ 1 ];
      }
    
   // alternative models:
   size_t d = 1;                // assume that always one SNP is being tested
   size_t d_ = d + ncovar;      // the number of weights (dimensionality of a single beta)

   // if the parameters have already been set in a previous iteration delete them.
   if ( beta      != nullptr ) { delete[] beta;      beta = nullptr; }   
   if ( delta     != nullptr ) { delete[] delta;     delta = nullptr; }
   if ( sigg2     != nullptr ) { delete[] sigg2;     sigg2 = nullptr; }
   if ( sige2     != nullptr ) { delete[] sige2;     sige2 = nullptr; }
   if ( LL        != nullptr ) { delete[] LL;        LL = nullptr; }
   if ( rank_beta != nullptr ) { delete[] rank_beta; rank_beta = nullptr; }
   if ( pvals     != nullptr ) { delete[] pvals;     pvals = nullptr; }
   if ( f_stat    != nullptr ) { delete[] f_stat;    f_stat = nullptr; }

   // re-allocate storage for all the variables
   f_stat    = new real[ nSNPs * d_ ];
   beta      = new real[ nSNPs * d_ ];
   delta     = new real[ nSNPs ];
   sigg2     = new real[ nSNPs ];
   sige2     = new real[ nSNPs ];
   LL        = new real[ nSNPs ];
   pvals     = new real[ nSNPs ];
   rank_beta = new BLAS_INT[ nSNPs ];
   }

// Functor holdong all the needed members of FastLMM_woodbury to be handed to the Brent's routine
// the constructor gets the members of FastLMM_woodbury
brentFunctor_FastLMM_woodbury::brentFunctor_FastLMM_woodbury( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU, const real *SU_y, const real *SU_x, const real *y, const real *X, const real *S, const real *SUcovar, const real *covar )
   {
   // member parameters from FastLMM_woodbury:
   this->d        = d;
   this->rankU    = rankU;
   this->ncovar   = ncovar;
   this->n        = n;
   this->SU_y     = SU_y;
   this->SU_x     = SU_x;
   this->y        = y;
   this->X        = X;
   this->covar    = covar;
   this->S        = S;
   this->SUcovar  = SUcovar;
   }

// destructor has nothing to do
brentFunctor_FastLMM_woodbury::~brentFunctor_FastLMM_woodbury(){}

// the operator evaluates nLLeval
real brentFunctor_FastLMM_woodbury::operator() ( const real logdelta )
   {
   // the rest are dummies:
   BLAS_INT rank_beta  = 0;
   real dummy_MLparams = 0.0;

   return FastLMM_woodbury::nLLeval( nullptr, logdelta, d, ncovar, n, rankU, SU_y, SU_x, y, X, S, SUcovar, covar, &dummy_MLparams, &dummy_MLparams, &dummy_MLparams, &rank_beta, false );
   }

#if 0
real *MatrixConcatenate( real *pSrc1, size_t cRows1, size_t cCols1, real *pSrc2, size_t cRows2, size_t cCols2 )
   {
   if ( cCols1 != cCols2 )
      {
      Fatal( "Programming error.  MatrixConcatenate cannot join matrices of different widths."
           "\n    cCols1: %d"
           "\n    cCols2: %d", cCols1, cCols2 );
      }
   real *pReal = new real[ (cRows1 + cRows2) * cCols1 ];
   size_t cEntries1 = cRows1 * cCols1;
   memmove( pReal, pSrc1,  cEntries1 * sizeof( real ) );
   real *pRealT = pReal + cEntries1;
   size_t cEntries2 = cRows2 * cCols2;
   memmove( pRealT, pSrc2, cEntries2 * sizeof( real ) );
   return( pReal );
   }
#endif

/*
 *  Send all the SNP pairs in this run to be processed.
 */
void FastLMM_woodbury::TrainNullAndAltForEpistatic()
   {
   size_t cIndividuals = n;         // use more descriptive names for now
   size_t cCovariates = ncovar;
   size_t cSnps = nSNPs;
   size_t cSnpPairsBetweenDisplay = fClusterIO ? 10000 : 1000;
   char   nl = fClusterIO ? '\n' : '\r';

   // Allocate a working array for covariate data.  We need 
   //   enough for cCovariates * cIndividuals + we need 
   //   to add the two more rows for snp2data and snp12data
   //   NOTE:::  if we change to triples, do we need more than two 
   //            additional rows?  ( cCovariates + 1 + 1 )
   real* rgCovariatesTemp = new real[ (cCovariates + 1 + 1) * cIndividuals ];    // allocate 2d matrix
   real* rgSUcovarTemp = new real[ (cCovariates + 1 + 1) * rankU ];              // space for the covariates transformed to the low rank space
   real* rgTemp = new real[ 2 * rankU ];                                         // work space in TrainNullAndAltForEpistaticWithSnpPairs()

   // Preload the initial covariate data.
   size_t cEntries = cCovariates * cIndividuals;
   memmove( rgCovariatesTemp, covar, cEntries * sizeof( real ) );    // copy original covar data

   size_t crgSUcovarTempEntries = cCovariates * rankU;
   memmove( rgSUcovarTemp, SUcovar, crgSUcovarTempEntries * sizeof( real ) ); // copy origianal SUcovar data

   if ( !snpFilterOpts.FFilterByTask() )
      {
      // We got here w/o the task being setup, so set
      //   1 task to run and this is task 0
      //   say how many snp we have in the block
      snpFilterOpts.SetByTask( 1, 0 ); // setup 1 task to run and run task 0
      snpFilterOpts.SetBlockSize( cSnps );
      }

   /*
    *  Tell the user where the output file is
    *    going and open it up... 
    */
   ProgressNL( "                       Write output: [%s]", fullOutputFilename.c_str() );
   FILE *pFile = fopen( fullOutputFilename.c_str(), "w" );
   if ( pFile == nullptr )
      {
      Fatal( "unable to open output file.  Sorry."
           "\n  CRT Error: %s", strerror( errno ) );
      }

   // Print the header and setup to measure our throughput
   // NOTE:  GCC and VC don't share the same format specifiers for 
   //        size_t, so we convert to use int in this case.
   fprintf( pFile,   "#  Tasks Requested: %9d", (int)snpFilterOpts.cTasksAsk );
   fprintf( pFile, "\n#     Tasks Needed: %9d", (int)snpFilterOpts.cTasksMax );
   fprintf( pFile, "\n#      Task Number: %9d", (int)snpFilterOpts.iTaskAsk );
   fprintf( pFile, "\n#       Partitions: %9d", (int)snpFilterOpts.cBlocksMax );
   fprintf( pFile, "\n#   SNPs/Partition: %9.4lf", snpFilterOpts.cSnpsInBlock );
   fprintf( pFile, "\n# pvalue Threshold: %.2e", pvaluePrintThreshold );

   size_t iBlk1 = snpFilterOpts.iBlock1;
   size_t iBlk2 = snpFilterOpts.iBlock2;
   fprintf( pFile, "\n#      Block1 Info: %9d %9d %9d", (int)iBlk1, (int)snpFilterOpts.IdxFromBlockId( iBlk1 ), (int)snpFilterOpts.CountOfSnpsInBlock( iBlk1 ) );
   fprintf( pFile, "\n#      Block2 Info:" ); // keep a constant number of header lines.
   if ( iBlk1 != iBlk2 )
      {
      fprintf( pFile, " %9d %9d %9d", (int)iBlk2, (int)snpFilterOpts.IdxFromBlockId( iBlk2 ), (int)snpFilterOpts.CountOfSnpsInBlock( iBlk2 ) );
      }

   fprintf( pFile, "\nsnpId_1\tsnpId_2\tpvalue\tSnpPair_Wgt\tSnpPair_OddsRatio\n" );

   size_t iSnp1Base = snpFilterOpts.IdxFromBlockId( iBlk1 );
   size_t iSnp2Base = snpFilterOpts.IdxFromBlockId( iBlk2 );
   size_t cSnps1Max = snpFilterOpts.CountOfSnpsInBlock( iBlk1 );
   size_t cSnps2Max = snpFilterOpts.CountOfSnpsInBlock( iBlk2 );
   size_t cSnpPairs = 0;
   size_t cSnpPairsWritten = 0;
   CTimer snpTimer( true );

   if ( fDoSnpId1 )
      {
      size_t iSnp1 = IndexFromList( snpInfo, snpId1Global );
      if ( iSnp1 == static_cast< size_t >( -1 ) )
         {
         Fatal( "Unable to locate SnpId [%s] in list of SNPs.", snpId1Global.c_str() );
         }

      Progress( "\n       Processing SNP pairs for SNP: %s", snpId1Global.c_str() );
      Progress( "\n                SnpPairs to process: %7d", (int)((cSnps1Max - 1 )) );   // take out the snpId1
      Progress( "\n             pvalue print Threshold: %.2e", pvaluePrintThreshold );
      Progress( "\n SnpPair %9d [%7.2lf SnpPairs/sec]", (int)cSnpPairs, (double)0.0 );

      for ( size_t iSnp2 = 0; iSnp2 < cSnps; ++ iSnp2 )
         {
         if ( iSnp2 != iSnp1 )
            {
            cSnpPairsWritten += TrainNullAndAltForEpistaticWithSnpPairs( pFile, iSnp1Base, iSnp1, iSnp1Base, iSnp2, rgCovariatesTemp, rgSUcovarTemp, rgTemp );
            if ( (++cSnpPairs % cSnpPairsBetweenDisplay) == 0 )
               {
               Progress( "%c SnpPair %9d [%7.2lf SnpPairs/sec]", nl, (int)cSnpPairs, (double)cSnpPairs / snpTimer.Elapsed() );
               }
            }
         }
      }
   else if ( iBlk1 != iBlk2 )
      {  //                                           Block Base count
      Progress( "\n           Number of blocks in File: %7d", (int)snpFilterOpts.cBlocksMax );
      Progress( "\n       Processing 2 blocks for Task: %7d of %7d", (int)snpFilterOpts.iTaskAsk, (int)snpFilterOpts.cTasksMax );
      Progress( "\n                       Block 1 Info: %7d  %7d  %7d", (int)iBlk1, (int)iSnp1Base, (int)cSnps1Max );
      Progress( "\n                       Block 2 Info: %7d  %7d  %7d", (int)iBlk2, (int)iSnp2Base, (int)cSnps2Max );
      Progress( "\n                SnpPairs to process: %7d", (int)(cSnps1Max * cSnps2Max ) );
      Progress( "\n             pvalue print Threshold: %.2e", pvaluePrintThreshold );
      Progress( "\n SnpPair %9d [%7.2lf SnpPairs/sec]", cSnpPairs, (double)0.0 );

      size_t iSnp2BaseOffset = iSnp2Base - cSnps1Max;    // account for Snp2 location in the combined array

      for ( size_t iSnp1 = 0; iSnp1 < cSnps1Max; ++iSnp1 )
         {
         for( size_t iSnp2 = cSnps1Max; iSnp2 < cSnps; ++iSnp2 )
            {
            cSnpPairsWritten += TrainNullAndAltForEpistaticWithSnpPairs( pFile, iSnp1Base, iSnp1, iSnp2BaseOffset, iSnp2, rgCovariatesTemp, rgSUcovarTemp, rgTemp );
            if ( (++cSnpPairs % cSnpPairsBetweenDisplay) == 0 )
               {
               Progress( "%c SnpPair %9d [%7.2lf SnpPairs/sec]", nl, (int)cSnpPairs, (double)cSnpPairs / snpTimer.Elapsed() );
               }
            }
         }
      }
   else
      {  //                                           Block Base count
      Progress( "\n           Number of blocks in File: %7d", (int)snpFilterOpts.cBlocksMax );
      Progress( "\n        Processing 1 block for Task: %7d of %7d", (int)snpFilterOpts.iTaskAsk, (int)snpFilterOpts.cTasksMax );
      Progress( "\n                       Block 1 Info: %7d  %7d  %7d", (int)iBlk1, (int)iSnp1Base, (int)cSnps1Max );
      Progress( "\n                       Block 2 Info:" ); // maintain consistent line count for output
      Progress( "\n                SnpPairs to process: %7d", (int)((cSnps1Max * (cSnps1Max-1)) / 2) );
      Progress( "\n             pvalue print Threshold: %.2e", pvaluePrintThreshold );
      Progress( "\n SnpPair %9d [%7.2lf SnpPairs/sec]", cSnpPairs, (double)0.0 );

      for ( size_t iSnp1 = 0; iSnp1 < cSnps-1; ++iSnp1 )
         {
         for ( size_t iSnp2 = iSnp1 + 1; iSnp2 < cSnps; ++ iSnp2 )
            {
            cSnpPairsWritten += TrainNullAndAltForEpistaticWithSnpPairs( pFile, iSnp1Base, iSnp1, iSnp1Base, iSnp2, rgCovariatesTemp, rgSUcovarTemp, rgTemp );
            if ( (++cSnpPairs % cSnpPairsBetweenDisplay) == 0 )
               {
               Progress( "%c SnpPair %9d [%7.2lf SnpPairs/sec]", nl, (int)cSnpPairs, (double)cSnpPairs / snpTimer.Elapsed() );
               }
            }
         }
      }

   // Report the end numbers
   fprintf( pFile,"\n** End **\n%d,%d\n", (int)cSnpPairsWritten, (int)cSnpPairs );
   fclose( pFile );
   snpTimer.Stop();
   Progress( "\n                    SnpPairs Tested: %7d [%7.2lf SnpPairs/sec]", (int)cSnpPairs, (double)cSnpPairs / snpTimer.Elapsed() );
   Progress( "\n SnpPairs passing pval print filter: %7d\n", (int)cSnpPairsWritten );

   delete [] rgTemp;
   delete [] rgCovariatesTemp;
   delete [] rgSUcovarTemp;
   }

/*
 * Given the index for 2 snps in the snp array, compute the pvalue and 
 *   supporting information for the snp pair
 *  Outputs are LL_null, LL_alt, pval
 * The null model uses snp1 + snp2, while the alternative model uses 
 *   snp1 + snp2 + snp1*snp2
 *
 * NB. No need to account for proximal contamination here because the 
 *     null model here *does* include snp1 + snp2
 */
size_t FastLMM_woodbury::TrainNullAndAltForEpistaticWithSnpPairs( FILE *pFile, size_t iSnp1Base, size_t iSnp1, size_t iSnp2Base, size_t iSnp2, real* rgCovariatesTemp, real* rgSUcovarTemp, real* rgTemp )
   {
   // in order to have the least chance of misusing the univariate SNP code, 
   // I will aribitrarily put snp1 in for the "single" snp in the univariate
   // test, and put the second snp, and the interaction effect with the covariates.
   size_t d = 1;
   size_t cIndividuals = n;         // use more descriptive names for now
   size_t cCovariates = ncovar;

   // pointers to data for Snp1 & Snp2 (one cIndividual's worth)
   real* snp1data = X + (iSnp1 * (cIndividuals * d));
   real* snp2data = X + (iSnp2 * (cIndividuals * d));

   // put the first SNP in here (this is arbitrary, could have been the 
   // second, or none), and yes this repeats above slightly one of the 
   // lines just above
   real* currentSU_x = SU_x + (iSnp1 * (d * rankU));
   real* currentX    = X    + (iSnp1 * (cIndividuals * d));

   // we are re-using the same memory locations since we don't want/need
   // to store all the information for now.
   real* currentBeta  = beta + 0;
   real* currentSigg2 = sigg2 + 0;
   real* currentSige2 = sige2 + 0;
   real* currentDelta = delta + 0;
   BLAS_INT* currentRank_beta = rank_beta + 0;

   //find delta:
   if ( covLearnType == FULL )
      {
      Fatal( "Not Implemented.  covLearnType not supported in this code" );
      }
   else if ( covLearnType == ONCE )
      {
      //Or: use the lock-in. So it equals delta_0
//         *currentDelta = delta_0[ iSnp ];        // original code
      *currentDelta = delta_0[ 0 ];           // TODO:  Review
      }

   //will be modified for each of null and alt, below, one building on the other
   //real* covarTmp;
   //real* SUcovarTmp;

   // Evaluate the null LL (includes snp1 + snp2) at delta_0 and 
   //   get the ML parameters
   // NOTE: could be made more efficient by rotating only each newly 
   //       added feature, but this should be fine for now (re-rotating 
   //       covariates plus new feature) especially for the large Azure 
   //       run which will not have covariates other than a bias term
   //----------------------------------------------------------------------
   // modify SUcovar and covar to also include snp2, by first computing 
   //  covar, then applying the rotation to obtain SUcovar
   // cRows = cIndividuals = n
   // cCols = cSnps = ncovar

   // We filled in the 'updated' covariates array with the original covariate
   //  data in the caller.
   // We now need to append snp2data to the end of the Covariates array we use.
   size_t crgCovariatesEntries = cCovariates * cIndividuals;
   size_t crgSUcovarEntries = cCovariates * rankU;       // and make sure the entries are accounted for in rgSUcovarTemp too.

   //  The initial fill of the working matrix is now done once in calling routine.
   //memmove( rgCovariatesTemp, covar, cEntries * sizeof( real ) );    // copy original covar data
   memmove( rgCovariatesTemp + crgCovariatesEntries, snp2data, cIndividuals * sizeof( real ) );    // add snp2data

   size_t cCovariatesTemp = cCovariates + 1;             // ONLY include snp2data data for now.

   // Obtain SUcovar by (1) rotating _just_ snp2data, 
   //   and then appending the rotated SNP to the previously rotated covariates
   //
   real* SUsnp2data = rgTemp;       // use workspace allocated above rather than allocating in loop new real[ rankU * 2 ];
   BLAS_INT n_int      = (BLAS_INT)cIndividuals;
   BLAS_INT rankU_int  = (BLAS_INT)rankU;
   BLAS_INT ncovar_int = (BLAS_INT)1;
   BlasLapack::gemm( &char_T, &char_N, &rankU_int, &ncovar_int, &n_int, &real_one, eigen->U, &n_int, snp2data, &n_int, &real_zero, SUsnp2data, &rankU_int );

   //  The initial fill of the working matrix rgSUcovarTemp is done once in calling routine.
   memmove( rgSUcovarTemp + crgSUcovarEntries, SUsnp2data, rankU * sizeof( real ) );   // add rotated snp2data to the end of the matrix

   //evaluate the LL at delta_0 and get the ML parameters
   real LL_null = (real)-1.0 * nLLeval( nullptr, 
                                       *currentDelta, (BLAS_INT)d, 
                                       (BLAS_INT)cCovariatesTemp, (BLAS_INT)cIndividuals, 
                                       (BLAS_INT)rankU, SU_y,
                                       currentSU_x, y, 
                                       currentX, eigen->S, 
                                       rgSUcovarTemp, rgCovariatesTemp,
                                       currentSigg2, currentSige2, 
                                       currentBeta, currentRank_beta, 
                                       true );

   //evaluate the alternative  LL (includes snp1 + snp2 + snp1*snp2) at delta_0 and get the ML parameters
   //--------------------------------------------------------------------------------------------
   //modify SUcovar and covar to also include the interaction effect, snp1*snp2
   //specifically, re-use the concatenation we had above for the null, but add in one more feature, snp1*snp2

   // create the SNP interaction effect feature vector as the last row of the rgCovariatesTemp array
   //  this should be set equal to the element-wise product of snp1data and snp2data:  
   //  i.e. for i=1:length(snp1data) snp12data=snp1data[i]*snp2data[i] end
   real* snp12data = (rgCovariatesTemp + crgCovariatesEntries) + cIndividuals;      // point to last row of table e.g. [cCovariates + 1]
   for ( size_t iIndividual = 0; iIndividual < cIndividuals; ++ iIndividual )
      {
      snp12data[ iIndividual ] = snp1data[ iIndividual ] * snp2data[ iIndividual ];
      }

   cCovariatesTemp += 1;         // Account for adding in the snp12data 

   // Obtain SUcovar by rotating the new feature and then concatenating it to existing value of rgSUcovarTemp
   real *SUsnp12data = SUsnp2data;     // re-use the memory allocation of SUsnp2data from above
   n_int = (BLAS_INT)cIndividuals;
   rankU_int  = (BLAS_INT)rankU;
   ncovar_int = (BLAS_INT)1;
   BlasLapack::gemm( &char_T, &char_N, &rankU_int, &ncovar_int, &n_int, &real_one, eigen->U, &n_int, snp12data, &n_int, &real_zero, SUsnp12data, &rankU_int );
   memmove( rgSUcovarTemp + crgSUcovarEntries + rankU, SUsnp12data, rankU * sizeof( real ) );

   //evaluate the LL at delta_0 and get the ML parameters
   real LL_alt  = (real)-1.0 * nLLeval( nullptr, 
                                       *currentDelta, (BLAS_INT)d, 
                                       (BLAS_INT)cCovariatesTemp, (BLAS_INT)cIndividuals, 
                                       (BLAS_INT)rankU, SU_y, 
                                       currentSU_x, y, 
                                       currentX, eigen->S, 
                                       rgSUcovarTemp, rgCovariatesTemp, 
                                       currentSigg2, currentSige2, 
                                       currentBeta, currentRank_beta, 
                                       true );

   // Compute p-value from LR test:
   // real df = *currentRank_beta - (real)rank_beta_0[ snp ];
   real df = 1.0; // TODO: Review
   real pval;

   //Likelihood-ratio test
   pval = Gamma::gammaQ( LL_alt - LL_null, (real)0.5 * df );

   // fix NaN output
   if ( _isnan( pval ) ) 
      {
      Warn( "Pval of SnpPair[%s:%s] NaN.  Set pval to 1.00e-000\n", snpInfo[iSnp1].idSnp.c_str(), snpInfo[iSnp2].idSnp.c_str() );
      pval = 1.00e-000;
      }

   if ( pvaluePrintThreshold < pval )
      {
      // pval did not get below the print threshold . . do not print this one
      return( 0 );   // number of SnpPairs printed ( 0 or 1 )
      }

   // Prepare the odds ratios for printing
   // snpWeight is last element of currentBeta[] (which is one larger than cCovariatesTemp)
   // multiply snp1 & snp2 StdDev to get StdDev for SnpPairs
   real snpPairWeight = currentBeta[ cCovariatesTemp ];
   real snpPairStdDeviation = snpInfo[ iSnp1 ].StdDeviation() * snpInfo[ iSnp2 ].StdDeviation();   
   real snpPairOddsRatio = exp( ::abs( snpPairWeight ) / snpPairStdDeviation );

   fprintf( pFile, "%s\t%s\t%e\t%e\t%e\n", 
                        snpInfo[iSnp1].idSnp.c_str(), snpInfo[iSnp2].idSnp.c_str(),
                        pval,
                        snpPairWeight, snpPairOddsRatio
                        );
   return( 1 );   // number of SnpPairs printed ( 0 or 1 )
   }
