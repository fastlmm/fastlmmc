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

// constructor takes an EigenSym object as input
FastLMM_fullrank::FastLMM_fullrank( EigenSym *eigen ) : FastLMM( eigen )
   {
   if ( this->eigen->n != this->eigen->rankU )
      {
      //Fatal( "\nerror, Eigen dimensions do not agree with FastLMM_fullrank. rankU=%i, n=%i", eigen->rankU, eigen->n );
      }

   Uy = nullptr;              // null pointer
   UX = nullptr;              // null pointer
   Ucovar = nullptr;          // null pointer
   UW = nullptr;              // null pointer
   nExclude = nullptr;
   this->excludeSnpStart = nullptr;
   modelID = "FastLMM_fullrank";
   }


FastLMM_fullrank::~FastLMM_fullrank()
   {
   if (Uy){       delete[] Uy; Uy = nullptr; }
   if (UX){       delete[] UX; UX = nullptr; }
   if (Ucovar){   delete[] Ucovar; Ucovar = nullptr; }
   if (UW){       delete[] UW; UW = nullptr; }
   if (nExclude){ delete[] nExclude; nExclude = nullptr; }
   if (excludeSnpStart){ delete[] excludeSnpStart; excludeSnpStart = nullptr; }
   }

// SNPs [n by nSNPs] holds the nSNPs SNPs to be tested for each of n individuals
void FastLMM_fullrank::set_SNPs( real *SNPs, size_t nSNPs_, vector<SnpInfo>& snpInfo_ )
   {
   if( UX != nullptr )
      {
      delete[] UX;
      }

   this->snpInfo = snpInfo_;
   this->nSNPs = nSNPs_;
   if (this->nSNPs != this->snpInfo.size())
      {
      Fatal("SNP counts do not match.  nSNPs = %i, snpInfo.size() = %i", this->nSNPs, this->snpInfo.size());
      }

   UX = new real[ nSNPs * n ];         // as number of SNPs might change, always reallocate memory for UX.
   BLAS_INT n_int = (BLAS_INT)n;
   BLAS_INT nSNPs_int = (BLAS_INT)nSNPs;
   BlasLapack::gemm(&char_T,  &char_N, &n_int, &nSNPs_int, &n_int, &real_one, eigen->U, &n_int, SNPs, &n_int, &real_zero, UX, &n_int );  //UX=U.T*SNPs
   }

// W [n by nSNPs] holds the nSNPsW SNPs to be tested for each of n individuals
void FastLMM_fullrank::set_SNPsCov(real *W, size_t nSNPsW, vector<SnpInfo>& snpInfoCov)
   {
   if( UW != nullptr )
      {
      delete[] UW;
      }

   if(excludeSnpParams.fileType==FileSim)
      {
      this->snpInfoCovariance = snpInfoCov;
      this->nSNPsW = nSNPsW;
      if (this->nSNPsW != this->snpInfoCovariance.size()){
         Fatal("SNP counts do not match.  nSNPs = %i, snpInfo.size() = %i", this->nSNPsW, this->snpInfoCovariance.size());
      }
      UW=new real[nSNPsW*n];   //as number of SNPs might change, alsways reallocate memory for UW.
      BLAS_INT n_int = (BLAS_INT)n;     //Warning: size_t BLAS_INT conversion
      BLAS_INT nSNPsW_int = (BLAS_INT)nSNPsW;  //Warning: size_t BLAS_INT conersion
      BlasLapack::gemm(&char_T,  &char_N, &n_int, &nSNPsW_int, &n_int, &real_one, eigen->U, &n_int, W, &n_int, &real_zero, UW, &n_int );  //UW=U.T*SNPs
      }
   else if(excludeSnpParams.fileType==TrainOnTest)
      {
      Fatal("These SNPs have undergone different preprocessing.");
      this->UW = this->UX;
      this->nSNPsW = this->nSNPs;
      this->snpInfoCovariance = this->snpInfo;
      }
   else if(excludeSnpParams.fileType==NoFile)
      {
      this->UW = nullptr;
      this->nSNPsW = 0;
      }
   }

// phenotype is a vector of length n, holding the (single) phenotype.
void FastLMM_fullrank::set_Phenotype( real *phenotype, string& phenLabel_ )
   {
   // things can be resized so clear any existing memory allocation and get clean memory
   if( Uy != nullptr ) { delete[] Uy; Uy = nullptr; }
   Uy = new real[ n ];              // get space for result

   phenLabel = phenLabel_;          // ensure we have the phenotype name saved locally
   BLAS_INT n_int = (BLAS_INT)n;    // Ensure we have a BLAS_INT for the call

   // Uy = U.T * phenotype
   BlasLapack::gemv( &char_T,       // trans = 'T'  => y := (alpha * A' * x) + (beta * y)
                     &n_int,        // m     = number of rows in matrix A
                     &n_int,        // n     = number of columns in matrix A
                     &real_one,     // alpha = specifies the scalar alpha
                     eigen->U,      // a     = Array[ lda, n ] leading m x n part of array a must contain matrix coefficients
                     &n_int,        // lda   = first dimension of a as declared in calling routine
                     phenotype,     // x     = array x[ m ] must contain valid vector x
                     &int_one,      // incx  = increment for elements of x
                     &real_zero,    // beta  = scalar beta
                     Uy,            // y     = array y[m] if beta != 0.0, y must contain valid vector y
                     &int_one       // incy  = increment for elements of y 
                     );   
   }

//covariates is an [n*ncovar] matrix holding the covariates, including the bias
void FastLMM_fullrank::set_covariates( real *covar, size_t ncovariates )
   {
   ncovar = ncovariates;
   if( Ucovar==nullptr )
      {
      Ucovar=new real[n*ncovar];
      }
   BLAS_INT n_int = (BLAS_INT)n;  //Warning: size_t BLAS_INT conversion
   BLAS_INT ncovar_int = (BLAS_INT)ncovar;  //Warning: size_t BLAS_INT conversion
   BlasLapack::gemm( &char_T, &char_N, &n_int, &ncovar_int, &n_int, &real_one, eigen->U, &n_int, covar, &n_int, &real_zero, Ucovar, &n_int );  //UX=U.T*SNPs
   }

//find ML parameters for the null model
void FastLMM_fullrank::trainNull()
   {

   BLAS_INT ncovar_int = (BLAS_INT)ncovar;
   BLAS_INT n_int = (BLAS_INT)n;
   real UXdummy = 0.0;  // As we do not plug in a SNP, we need a dummy

   // find optimal delta:
   if ( logDeltaPreset == logDeltaPreset )
      {
      delta_0[0]=logDeltaPreset;
      }
   else
      {
      if ( FastLMM::brentExhaustive )
         {
         // Find the optimal delta
         delta_0[0] = FastLMM_fullrank::optimize_delta_Brentintervals( 0, ncovar_int, n_int, Uy, &UXdummy, eigen->S, Ucovar, &UXdummy, 0 );
         }
      else
         {
         //Find the optimal delta
         delta_0[0] = FastLMM_fullrank::optimize_delta_BrentintervalsCheckTriplets( 0, ncovar_int, n_int, Uy, &UXdummy, eigen->S, Ucovar, &UXdummy, 0 );
         }
      }

   //get the ML parameters for the optimal delta:
   LL_0[0] = (real)-1.0 * nLLeval( nullptr, delta_0[0], 0, ncovar_int, n_int, Uy, &UXdummy, eigen->S, Ucovar, &sigg2_0[0], &sige2_0[0], &beta_0[0], &rank_beta_0[0], true,&UXdummy,0);  //evaluate the LL at delta_0 and get the ML parameters
   if ( verboseLevel >= verboseDebug )
      {
      printf( "\n------------------------------------------" );
      printf( "\nNULL: LL:%f", LL_0[0] );
      printf( "\n  sigg2=%f sige2=%f delta=%f", sigg2_0[0], sige2_0[0], delta_0[0] );
      printf( " beta: ");
      for( size_t i=0; i < ncovar; ++i )
         {
         printf( "%f ", beta_0[i] );
         }
      printf( "\n------------------------------------------" );
      printf( "\n" );
      }
   //copy the "default" null model parameters
   for ( size_t snp = 1; snp < nSNPs; ++snp )
      {
      delta_0[ snp ] = delta_0[ 0 ];
      rank_beta_0[ snp ] = rank_beta_0[ 0 ];
      LL_0[ snp ] = LL_0[ 0 ];
      sigg2_0[ snp ] = sigg2_0[ 0 ];
      sige2_0[ snp ] = sige2_0[ 0 ];
      for (size_t i_covar = 0; i_covar < ncovar; ++i_covar )
         {
         beta_0[ (snp * ncovar) + i_covar ] = beta_0[ i_covar ];
         }
      }
   }

//find ML parameters for all SNPs
void FastLMM_fullrank::trainAlt()
   {
   //alternative models:
   size_t d = 1;           //assume that always one SNP is being tested
   size_t d_= d + ncovar;  //the number of weights (dimensionality of a single beta)

   //train the alternative model for the SNP
   #pragma omp parallel for //This pragma will mess up with the SNP exclusion function
   for ( int snp=0; snp<(int)nSNPs; ++snp )
      { // (parallel for loop)  //Warning: size_t int conversion
         {
         //set pointers to the current SNP and ML parameters:
         real* currentUX = UX + (snp * (n*d));
         real* currentBeta = beta + (snp * d_);
         real* currentSigg2 = sigg2+snp;
         real* currentSige2 = sige2+snp;
         real* currentDelta = delta+snp;
         real* UWcurrent = nullptr;
         BLAS_INT* currentRank_beta = rank_beta + snp;
         
         if ( nExclude[ snp ] > 0 )
            {
            // printf( "SNP: %i nExcl: %i\n", (int)snp, (int)nExclude );
            UWcurrent = UW + (n * this->excludeSnpStart[snp]);
            }

         if ( (nExclude[ snp ] > 0)
           && (!use_f_test)
           && ((covLearnType == REFIT_EXCLUDE) || (covLearnType == FULL)) )
            {
            // update null delta
            if ( FastLMM::brentExhaustive )
               {
               // Find the optimal delta
               delta_0[ snp ] = FastLMM_fullrank::optimize_delta_Brentintervals( 0, (BLAS_INT)ncovar, (BLAS_INT)n, Uy, currentUX, eigen->S, Ucovar, UWcurrent, (BLAS_INT)nExclude[ snp ] );
               }
            else
               {
               // Find the optimal delta
               delta_0[ snp ] = FastLMM_fullrank::optimize_delta_BrentintervalsCheckTriplets( 0, (BLAS_INT)ncovar, (BLAS_INT)n, Uy, currentUX, eigen->S, Ucovar, UWcurrent, (BLAS_INT)nExclude[ snp ] );
               }
            }
         if ( nExclude[ snp ] > 0 )
            {
            //update the null likelihood
            LL_0[snp]=(real)-1.0*nLLeval(nullptr, delta_0[snp], 0, (BLAS_INT)ncovar, (BLAS_INT)n,Uy, currentUX, eigen->S, Ucovar, &sigg2_0[snp], &sige2_0[snp], &beta_0[((size_t)snp)*ncovar], &rank_beta_0[snp], true, UWcurrent, (BLAS_INT)this->nExclude[ snp ] );//evaluate the LL at delta_0 and get the ML parameters   //Warning: size_t BLAS_INT conversion
            }
         if ( (covLearnType == FULL)
           || ((nExclude[ snp ] > 0) && (use_f_test) && (covLearnType == REFIT_EXCLUDE)) )
            {
            //update alternative delta
            if ( FastLMM::brentExhaustive )
               {
               *currentDelta = FastLMM_fullrank::optimize_delta_Brentintervals( (BLAS_INT)d, (BLAS_INT)ncovar, (BLAS_INT)n, Uy,currentUX, eigen->S, Ucovar, UWcurrent, (BLAS_INT)this->nExclude[ snp ] ); //Find the optimal delta  //Warning: size_t BLAS_INT conversion
               }
            else
               {
               *currentDelta = FastLMM_fullrank::optimize_delta_BrentintervalsCheckTriplets( (BLAS_INT)d, (BLAS_INT)ncovar, (BLAS_INT)n, Uy,currentUX, eigen->S, Ucovar, UWcurrent, (BLAS_INT)this->nExclude[snp] ); //Find the optimal delta   //Warning: size_t BLAS_INT conversion
               }
            }
         else if ( (covLearnType == ONCE) 
              || ((nExclude[snp] == 0)  && (covLearnType == REFIT_EXCLUDE))
              || ((nExclude[snp] > 0)  && (covLearnType == REFIT_EXCLUDE) && (!use_f_test)) )
            {
            //use existing null-delta
            *currentDelta = delta_0[ snp ];
            }
         else
            {
            Fatal( "Illegal covLearnType" );
            }
         //compute the alternative likelihood
         LL[ snp ] = (real)-1.0 * nLLeval( f_stat + (snp*d_), *currentDelta, (BLAS_INT)d, (BLAS_INT)ncovar, (BLAS_INT)n, Uy, currentUX, eigen->S, Ucovar, currentSigg2, currentSige2, currentBeta, currentRank_beta, true, UWcurrent, (BLAS_INT)this->nExclude[ snp ] );   //evaluate the LL at delta_0 and get the ML parameters   //Warning: size_t BLAS_INT conversion
         }
      //Compute p-value from LR test:
      real df = rank_beta[ snp ] - ((real)1.0 * rank_beta_0[ snp ]) ;
      if ( df > (real)0.0 )
         {
         //check the degrees of freedom for validity
         if ( use_f_test )
            {
            //Fishers F-test
            //
            // PvalUsingFisherFTestCdf focuses on keeping precision 
            //    for very small pvalues by inverting on the reverse 
            //    path of the Beta::BetaInc() branch
            pvals[ snp ] = PvalUsingFisherFTestCdf( f_stat[ snp * d_ ], df, (real)(n - rank_beta[ snp ]) );
            //pvals[ snp ] = 1.0 - FisherF::Cdf( f_stat[snp*d_], df, (real)(n - rank_beta[snp]) );
            }
         else
            {
            //Likelihood-ratio test
            pvals[ snp ] = PvalFromLikelihoodRatioTest( LL[ snp ] - LL_0[ snp ], ((real)0.5 * df) );
            }
         }
      else
         {
         //This happens for example if the SNP is a constant. In this case df=0
         pvals[snp]=std::numeric_limits<real>::quiet_NaN();
         }
   if ( verboseLevel >= verboseDebug )
         {
         printf( "\n- LL[%i]:%f p=%f", snp, LL[snp], pvals[snp] );
         printf( "\n  sigg2=%f sige2=%f delta=%f", sigg2[snp], sige2[snp], delta[snp] );
         printf( " beta: ");
         for( size_t i=0; i<d_; ++i )
            {
            printf( "%f ", beta[ (snp*d_)+i ] );
            }
         }
      }
   }

//static function:
//finds delta by optimizing over delta intervals
real FastLMM_fullrank::optimize_delta_BrentintervalsCheckTriplets( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const real *Uy, const real *currentUX, const real *S, const real *Ucovar, const real *UW, const BLAS_INT nExclude )
   {
   real *negLL_interval = new real[ this->numberIntervals + 1 ];
   real *log_delta_interval = new real[ this->numberIntervals + 1 ];
   
   real best_log_delta = std::numeric_limits< real >::quiet_NaN();
   real best_negLL = std::numeric_limits< real >::infinity();

   real MLparamsdummy = 0.0;
   
   BLAS_INT rankDummy = 0;
   size_t nevals = 0;

   //evaluate the boundaries of all intervals
   #pragma omp parallel for
   for( int i = 0; i <= (int)this->numberIntervals; ++i )
      {
      log_delta_interval[ i ] = this->minLogVal + ( i * (this->maxLogVal-this->minLogVal) / this->numberIntervals );
      negLL_interval[ i ] = FastLMM_fullrank::nLLeval( nullptr, log_delta_interval[ i ], d, ncovar, n, Uy, currentUX, S, Ucovar, &MLparamsdummy, &MLparamsdummy, &MLparamsdummy, &rankDummy, false, UW, nExclude );
      ++nevals;
      if ( negLL_interval[ i ] <= best_negLL )
         {
         best_negLL = negLL_interval[ i ];
         best_log_delta = log_delta_interval[ i ];
         }
      }

   brentFunctor_FastLMM_fullrank func = brentFunctor_FastLMM_fullrank( d, ncovar, n, Uy,currentUX, S, Ucovar, UW, nExclude );

   #pragma omp parallel for  //This might actually cause more overhead, most of the loops are not calling brents
   for ( int i = 1; i < (int)this->numberIntervals; ++i )
      {
      if ( (negLL_interval[ i ] < negLL_interval[ i + 1 ])
        && (negLL_interval[ i ] < negLL_interval[ i - 1 ])
         )
         {
         //check wether a local optimum exists in this triplet
         real current_negLL = std::numeric_limits< real >::infinity();
         size_t nevals_Int = 0;
         //If there is a local optimum, optimize over the current triplet:
         real current_log_delta = BrentC::minimize( func, log_delta_interval[ i - 1 ], log_delta_interval[ i + 1 ], this->brentTol, current_negLL, nevals_Int, this->brentMaxIter, true );
         nevals += nevals_Int;
         if ( current_negLL <= best_negLL )
            {
            best_negLL = current_negLL;
            best_log_delta = (current_log_delta);
            }
         }
      }

   delete []negLL_interval;
   negLL_interval = 0;
   delete []log_delta_interval;
   log_delta_interval = 0;

#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
   const char *szFmt1 = "\nnevals: %Iu logdelta: %f\n";
#else                      // Linux/g++ uses a %zu specifier for size_t
   const char *szFmt1 = "\nnevals: %zu logdelta: %f\n";
#endif
   Verbose( verboseDebug, szFmt1, nevals, best_log_delta );
   return( best_log_delta );
   }

//static function:
//finds delta by optimizing over delta intervals
real FastLMM_fullrank::optimize_delta_Brentintervals( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const real *Uy, const real *currentUX,const real *S, const real *Ucovar, const real *UW, const BLAS_INT nExclude )
   {
   real best_log_delta = std::numeric_limits< real >::quiet_NaN();
   real best_negLL = std::numeric_limits< real >::infinity();
   brentFunctor_FastLMM_fullrank func = brentFunctor_FastLMM_fullrank( d, ncovar, n, Uy,currentUX, S, Ucovar, UW, nExclude );
   size_t nevals = 0;

   #pragma omp parallel for
   for( int i = 0; i < (int)this->numberIntervals; ++i )
      {
      real current_log_delta_min = this->minLogVal + (i * (this->maxLogVal - this->minLogVal) / this->numberIntervals);
      real current_log_delta_max = this->minLogVal + (( i + 1) * (this->maxLogVal - this->minLogVal ) / this->numberIntervals);
      real current_negLL = std::numeric_limits< real >::infinity();
      size_t nevals_Int = 0;
      real current_log_delta = BrentC::minimize( func, current_log_delta_min, current_log_delta_max, this->brentTol, current_negLL, nevals_Int, this->brentMaxIter, true );
      nevals += nevals_Int;
      if ( current_negLL <= best_negLL )
         {
         best_negLL = current_negLL;
         best_log_delta = (current_log_delta);
         }
      }

#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
   const char *szFmt1 = "\nnevals: %Iu logdelta: %f\n";
#else                      // Linux/g++ uses a %zu specifier for size_t
   const char *szFmt1 = "\nnevals: %zu logdelta: %f\n";
#endif
   Verbose( verboseDebug, szFmt1, nevals, best_log_delta );
   return( best_log_delta );
   }

//static function:
//cleaned up version of nLLeval() for the full rank model
real FastLMM_fullrank::nLLeval(real *f_statistic, const real logdelta,const BLAS_INT d, const BLAS_INT ncovar,const BLAS_INT n,const real *Uy, const real *UX, const real *S, const real *Ucovar, real *sigg2, real *sige2, real *beta, BLAS_INT *rank_beta, const bool returnMLparams, const real *UW, BLAS_INT nExclude)
   {  
   real det = 0.0;                  // the log determinant is the sum of the diagonal of (S+delta*I)
   const BLAS_INT d_ = d + ncovar;  // the number of parameters to be estimated is number of SNP coefficients + number of covariates
   const size_t d_n = ((size_t)d_ * (size_t)n);
   real *Sd = new real[n];          // holds the sqrt of the diagonal values of (S+delta*I) (deleted)
   real *beta_temp = new real[d_];  // the parameters.  (deleted)
   real *XSdi = new real[d_n];      // XSdi This one will be used to compute the residual in the end. (deleted)
   bool delete_f_stat = false;

   if ( f_statistic == nullptr )
      {
      // Allocate space for f_statistic if needed
      // and flag for deletion if we do
      f_statistic = new real[d_];
      delete_f_stat = true;
      }
   
   real delta;
   delta = exp(logdelta);

   real yy=0.0;
   //--------------------------------------------------
   //compute the determinant and ySdi
   if(useHeritability)
      {
      for(size_t i=0;(BLAS_INT)i<n;++i)
         {
         Sd[i] = logdelta*S[i]+(1.0-logdelta);
         }
      }
   else{
      for(size_t i=0;(BLAS_INT)i<n;++i)
         {
         Sd[i]=S[i]+delta;
         }
       }
      real *UX_=new real[d_n];
   //compute XSdi. The inner two for loops concatenate UX and Ucovar
   for (size_t ind=0;(BLAS_INT)ind<n;++ind)
      {
      det+=log(Sd[ind]);
      yy+=Uy[ind]*Uy[ind]/Sd[ind];
      }
   
   size_t ind0=0;
   for(size_t snp=0;(BLAS_INT)snp<d;++snp)
      {
      for(size_t ind=0;(BLAS_INT)ind<n;++ind)
         {
         XSdi[ind0]=UX[ind0]/Sd[ind];
         UX_[ind0]=UX[ind0];
         ++ind0;
         }
      }
   size_t ind1=0;
   for(size_t cov=0;(BLAS_INT)cov<ncovar;++cov)
      {
      for(size_t ind=0;(BLAS_INT)ind<n;++ind)
         {
         XSdi[ind0]=Ucovar[ind1]/Sd[ind];
         UX_[ind0]=Ucovar[ind1];
         ++ind0;
         ++ind1;
         }
      }

   real *XX= new real[d_*d_];   //compute XX = XSdi.T*XSdi
   real *Xy= new real[d_];      //compute Xy = XSdi.T*y
   if (d_>0)
      {
      BlasLapack::gemm(&char_T,&char_N,&d_,&d_,&n,&real_one,XSdi,&n,UX_,&n,&real_zero,XX,&d_);
      BlasLapack::gemv(&char_T, &n, &d_, &real_one, XSdi, &n, Uy, &int_one, &real_zero, Xy, &int_one );
      }

   //Exclude SNPs from the RRM in the likelihood evaluation
   if (nExclude>0)   
      {
      real *WSdi = new real[((size_t)nExclude)*((size_t)n)];
      size_t ind0=0;
      for(size_t snp=0;(BLAS_INT)snp<nExclude;++snp)
         {
         for(size_t ind=0;(BLAS_INT)ind<n;++ind)
            {
            WSdi[ind0]=UW[ind0]/Sd[ind];
            ++ind0;
            }
         }

      //compute (I-WW)
      size_t nExclude2 = ((size_t)nExclude)*((size_t)nExclude);
      real *WW = new real[nExclude2];
      //subtract WW
      BlasLapack::gemm(&char_T,&char_N,&nExclude,&nExclude,&n,&real_minusone,WSdi,&n,UW,&n,&real_zero,WW,&nExclude);

      //add the identity
      for(size_t snp=0;(BLAS_INT)snp<nExclude;++snp)
         {
         ++WW[(snp*((size_t)nExclude)) + snp];
         }

      real *WX = new real[(size_t)nExclude * (size_t)d_];
      real *Wy = new real[nExclude];

      if (d_>0)
         {
         BlasLapack::gemm(&char_T,&char_N,&nExclude,&d_,&n,&real_one,WSdi,&n,UX_,&n,&real_zero,WX,&nExclude);
         }
      BlasLapack::gemv(&char_T, &n, &nExclude, &real_one, WSdi, &n, Uy, &int_one, &real_zero, Wy, &int_one );
      
      //compute the inverse and determinant of WW:
      //1) eigen-decomposition of (I-WW)
      real *U_WW= new real[nExclude2];  //Warning: BLAS_INT product
      real *S_WW= new real[nExclude];
      EigenSym::EV_symmetric(nExclude, WW, U_WW, S_WW);
      //2) determinant and solves using eigen-decomposition
      
      // U_W^T * Wy
      real *UWy = new real[nExclude];
      BlasLapack::gemv(&char_T, &nExclude, &nExclude ,&real_one, U_WW, &nExclude, Wy, &int_one, &real_zero, UWy, &int_one );
      // U_W^T * WX 
      real *UWX = new real[(size_t)nExclude * d_];
      BlasLapack::gemm(&char_T,  &char_N, &nExclude, &d_, &nExclude, &real_one, U_WW, &nExclude, WX, &nExclude, &real_zero, UWX, &nExclude );  //UX=U.T*SNPs

      //compute S_WW^{-1} * UWX
      size_t col = 0;
      for (size_t snp = 0; snp<(size_t)d_; ++snp)
         {
         for (size_t i_ex = 0; i_ex<(size_t)nExclude; ++i_ex)
            {
            //if (S_WW[i_ex]>MathFunctions::eps_rank || S_WW[i_ex]<-MathFunctions::eps_rank)
               {
               WX[col + i_ex] = UWX[col + i_ex]/S_WW[i_ex];
               }
            //else
            //   {
            //   WX[col + i_ex] = (real)0.0;
            //   }
            }
            col+=(size_t)nExclude;
         }

      //compute S_WW^{-1} * UWy and the determinant update
      for (size_t i_ex = 0; i_ex<(size_t)nExclude; ++i_ex)
         {
         //if (S_WW[i_ex]>MathFunctions::eps_rank || S_WW[i_ex]<-MathFunctions::eps_rank)
            {
            Wy[i_ex] = UWy[i_ex]/S_WW[i_ex];
            //update the logdet term in the likelihood
            det += log(S_WW[i_ex]);
            //update yy += yWWy
            yy+=Wy[i_ex] * UWy[i_ex];
            //yy_up+=Wy[i_ex] * UWy[i_ex];
            }
         //else
         //   {
         //   Wy[i_ex] = (real)0.0;
         //   }
         }

      if (d_>0)
         {
         //update XX += XWWX
         BlasLapack::gemm(&char_T,&char_N,&d_,&d_,&nExclude,&real_one,UWX,&nExclude,WX,&nExclude,&real_one,XX,&d_);
         //update Xy += XWWy
         BlasLapack::gemv(&char_T, &nExclude, &d_, &real_one, UWX, &nExclude, Wy, &int_one, &real_one, Xy, &int_one );
         }

      delete[] WSdi;
      delete[] WW;
      delete[] WX;
      delete[] Wy;
      delete[] UWy;
      delete[] UWX;
      delete[] U_WW;
      delete[] S_WW;
      }

   delete[] UX_;
   delete[] Sd;           //Sd is no longer needed
   //This part should now be the same for all branches, except the sige2/sigg2 issue:

   //solve the least squares problem.
   real logDetXX=0.0;
   FastLMM::solveLstSqSym(d_,XX,Xy,beta_temp,rank_beta, f_statistic,&logDetXX);

   real sigg2_temp=yy;  //the ML sigma_g^2
   for (size_t dim=0; (BLAS_INT)dim<d_;++dim)
      {
      sigg2_temp-=Xy[dim]*beta_temp[dim];
      }

   real negLL;       // negative Log Likelihood Temp variable

   if(mlType == ML)
      {
      sigg2_temp /= n;
      negLL=(real)0.5*(n*(_log2pi+log(sigg2_temp)+(real)1.)+det);//compute the negative LL using the ML parameters at delta
      }
   else if(mlType == REML)
      {
      sigg2_temp /= n-d_;
      negLL=(real)0.5*( (n-d_)*(_log2pi+log(sigg2_temp)+(real)1.)+det + logDetXX);//compute the negative LL using the ML parameters at delta
      }
   else
      {
      Fatal("\nOnly ML or REML are valid ML-types");
      }

   delete[] XSdi;          //XSdi is no longer needed.

   if (returnMLparams)//only return the ML params if needed. (e.g. no need during Brents search) This avoids danger of Memory errors
      {
      if ( useHeritability )
         {
         *sigg2=logdelta*sigg2_temp;
         *sige2=(1.0-logdelta)*sigg2_temp;
         }
      else
         {
         *sigg2=sigg2_temp;
         *sige2=delta**sigg2;
         }
      BlasLapack::copy(&d_,beta_temp,&int_one,beta,&int_one);
      }

   delete[] beta_temp;
   delete[] XX;
   delete[] Xy;
   
   if (delete_f_stat){delete[] f_statistic;}
   else
      {
      for (size_t dim = 0;(BLAS_INT)dim<d_;++dim){f_statistic[dim]/=sigg2_temp;}
      }
   return negLL;
   }

// Functor holding all the needed members of FastLMM_fullrank 
//   to be handed to the Brent's routine.  The constructor 
//   gets the members of FastLMM_fullrank
brentFunctor_FastLMM_fullrank::brentFunctor_FastLMM_fullrank( const BLAS_INT d_, const BLAS_INT ncovar_, const BLAS_INT n_, const real *Uy_, const real *UX_, const real *S_, const real *Ucovar_, const real *UW, const BLAS_INT nExclude )
   {
   // member parameters from FastLMM_fullrank:
   d = d_;
   this->nExclude = nExclude;
   this->UW = UW;
   ncovar = ncovar_;
   n = n_;
   Uy = Uy_;
   UX = UX_;
   S = S_;
   Ucovar = Ucovar_;
   }

void FastLMM_fullrank::init()
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
   size_t d  = 1;          //assume that always one SNP is being tested
   size_t d_ = d + ncovar; //the number of weights (dimensionality of a single beta)

      // if the parameters have already been set in a previous iteration delete them.
   if ( beta   != nullptr ) { delete[] beta;   beta   = nullptr; }
   if ( delta  != nullptr ) { delete[] delta;  delta  = nullptr; }
   if ( sigg2  != nullptr ) { delete[] sigg2;  sigg2  = nullptr; }
   if ( sige2  != nullptr ) { delete[] sige2;  sige2  = nullptr; }
   if ( LL     != nullptr ) { delete[] LL;     LL     = nullptr; }
   if ( pvals  != nullptr ) { delete[] pvals;  pvals  = nullptr; }
   if ( f_stat != nullptr ) { delete[] f_stat; f_stat = nullptr; }
   if ( rank_beta != nullptr ) { delete[] rank_beta; rank_beta = nullptr; }

   f_stat = new real[ nSNPs * d_ ];
   beta   = new real[ nSNPs * d_ ];
   delta  = new real[ nSNPs ];
   sigg2  = new real[ nSNPs ];
   sige2  = new real[ nSNPs ];
   LL     = new real[ nSNPs ];
   pvals  = new real[ nSNPs ];
   rank_beta = new BLAS_INT[ nSNPs ];
   }

//destructor has nothing to do
brentFunctor_FastLMM_fullrank::~brentFunctor_FastLMM_fullrank(){}

//the operator evaluates nLLeval
real brentFunctor_FastLMM_fullrank::operator()( const real logdelta )
   {
   //the rest are dummies
   BLAS_INT rank_beta_dummy = 0;
   real dummy_ML_params = 0.0;
   return FastLMM_fullrank::nLLeval( nullptr, logdelta, d, ncovar, n, Uy, UX, S, Ucovar, &dummy_ML_params, &dummy_ML_params, &dummy_ML_params, &rank_beta_dummy, false, UW, nExclude );
   }
