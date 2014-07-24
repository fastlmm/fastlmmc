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
FastLMM_lowrank::FastLMM_lowrank( EigenSym *eigen_ ):FastLMM( eigen_ )
   {
   // initialize all member variables in the constructor to prevent garbage operation
   if ( this->eigen->n <= this->eigen->rankU )
      {
      Fatal( "Eigen dimensions do not agree with FastLMM_lowrank. rankU=%i, n=%i", eigen->rankU, eigen->n );
      }

   modelID = "FastLMM_lowrank";     // in FastLMM sub-class

   rankU    = eigen_->rankU;
   Uy       = nullptr;
   y_UUy    = nullptr;
   UX       = nullptr;
   X_UUX    = nullptr;
   Ucovar   = nullptr;
   covar_UUcovar  = nullptr;
   W_UUW    = nullptr;
   UW       = nullptr;
   }

 
FastLMM_lowrank::~FastLMM_lowrank(void)
   {
   if (Uy){       delete[] Uy;     Uy     = nullptr; }      
   if (UX){       delete[] UX;     UX     = nullptr; }        
   if (Ucovar){   delete[] Ucovar; Ucovar = nullptr; }
   if (covar_UUcovar){ delete[] covar_UUcovar; covar_UUcovar = nullptr; }
   if (X_UUX){    delete[] X_UUX;  X_UUX  = nullptr; }
   if (y_UUy){    delete[] y_UUy;  y_UUy  = nullptr; }
   if (UW){       delete[] UW;     UW     = nullptr; }
   if (W_UUW){    delete[] W_UUW;  W_UUW  = nullptr; }
   }
 
// SNPs [n by nSNPs] holds the nSNPs SNPs to be tested for each of n individuals SNPs will be overwritten by X_UUX.
void FastLMM_lowrank::set_SNPs(real *SNPs, size_t nSNPs_, vector<SnpInfo>& snpInfo_ )
   {   
   this->snpInfo = snpInfo_;
   this->nSNPs = nSNPs_;
   if (this->nSNPs != this->snpInfo.size())
      {
      Fatal("SNP counts do not match.  nSNPs = %i, snpInfo.size() = %i", this->nSNPs, this->snpInfo.size());
      }
   if(UX!=0){        delete[] UX; }
   if(X_UUX != 0){   delete[] X_UUX; }
   
   UX = new real[ nSNPs * rankU ];     //as number of SNPs might change, alsways reallocate memory for UX.
   BLAS_INT n_int = (BLAS_INT)n;
   BLAS_INT nSNPs_int = (BLAS_INT)nSNPs;
   BLAS_INT rankU_int = (BLAS_INT)rankU;

   // UX = U.T * SNPs
   BlasLapack::gemm( &char_T, &char_N, &rankU_int, &nSNPs_int, &n_int, &real_one, eigen->U, &n_int, SNPs, &n_int, &real_zero, UX, &rankU_int );

   size_t sizeX = n * nSNPs;  
   X_UUX = new real[ sizeX ];
#if 1
   for(size_t snp=0;snp<sizeX;snp+=n)
      {
      BlasLapack::copy(&n_int, SNPs + snp, &int_one, X_UUX + snp, &int_one);
      }
#else
   BLAS_INT sizeX_int =(BLAS_INT)sizeX;  //Warning: size_t int converision
   BlasLapack::copy(&sizeX_int, SNPs, &int_one, X_UUX, &int_one);            //Initialize X_UUX with the SNPs. This is an [n * nSNPs] matrix //Warning size_t BLAS_INT conversion product
#endif
   // This matrix multiply performs X_UUX = SNPs - U * UX
   BlasLapack::gemm( &char_N, &char_N, &n_int, &nSNPs_int, &rankU_int, &real_minusone, eigen->U, &n_int, UX, &rankU_int, &real_one, X_UUX, &n_int );
   }

// W [n by nSNPs] holds the nSNPsW SNPs to be tested for each of n individuals
void FastLMM_lowrank::set_SNPsCov( real *W, size_t nSNPsW_, vector<SnpInfo>& snpInfoCov_ )
   {
   if( UW != nullptr )
      {
      delete[] UW;
      }
   if( W_UUW != nullptr )
      {
      delete[] W_UUW;
      }

   if ( excludeSnpParams.fileType == FileSim )
      {
      this->snpInfoCovariance = snpInfoCov_;
      this->nSNPsW = nSNPsW_;
      if (this->nSNPsW != this->snpInfoCovariance.size())
         {
         Fatal( "SNP counts do not match.  nSNPs = %i, snpInfo.size() = %i", this->nSNPsW, this->snpInfoCovariance.size() );
         }

      UW = new real[ nSNPsW_ * this->rankU ];  //as number of SNPs might change, alsways reallocate memory for UW.
      W_UUW = new real[ nSNPsW_ * this->n ];   //as number of SNPs might change, alsways reallocate memory for UW.
      
      BLAS_INT n_int = (BLAS_INT)n;
      BLAS_INT rankU_int = (BLAS_INT)rankU;
      BLAS_INT nSNPsW_int = (BLAS_INT)nSNPsW;
      
      BlasLapack::gemm(&char_T,  &char_N, &rankU_int, &nSNPsW_int, &n_int, &real_one, eigen->U, &n_int, W, &n_int, &real_zero, UW, &rankU_int );   //UX=U.T*SNPs
      for ( size_t i = 0; i< n * nSNPsW_; ++i )
         {
         W_UUW[ i ] = W[ i ];
         }

      // This matrix multiply performs W_UUW = W - (U * UW)
      BlasLapack::gemm( &char_N, &char_N, &n_int, &nSNPsW_int, &rankU_int, &real_minusone, eigen->U, &n_int, UW, &rankU_int, &real_one, W_UUW, &n_int );
      }
   else if ( excludeSnpParams.fileType == TrainOnTest )
      {
      Fatal( "These SNPs have undergone different preprocessing." );
      this->UW = this->UX;
      this->W_UUW = this->X_UUX;
      this->nSNPsW = this->nSNPs;
      this->snpInfoCovariance = this->snpInfo;
      }
   else if ( excludeSnpParams.fileType == NoFile )
      {
      this->UW = nullptr;
      this->W_UUW = nullptr;
      this->nSNPsW = 0;
      }
   }

//phenotype is a vector of length n, holding the (single) phenotype.
void FastLMM_lowrank::set_Phenotype( real *phenotype, string& phenLabel_ )
   {
   if ( Uy != 0 )    { delete[] Uy; }
   if ( y_UUy != 0 ) { delete[] y_UUy; }
   phenLabel = phenLabel_;

   Uy = new real[ rankU ];                // The phenotypes transformed to low rank
   BLAS_INT n_int = (BLAS_INT)n;
   BLAS_INT rankU_int = (BLAS_INT)rankU;

   // Uy = U.T * phenotype
   BlasLapack::gemv( &char_T, &n_int, &rankU_int, &real_one, eigen->U, &n_int, phenotype, &int_one, &real_zero, Uy, &int_one );
   
   // Initialize with original phenotype and then 
   //  U_UUy = phenotype - (U * Uy)
   y_UUy = new real[ n ];
   BlasLapack::copy( &n_int, phenotype, &int_one, y_UUy ,&int_one );
   BlasLapack::gemv( &char_N, &n_int, &rankU_int, &real_minusone, eigen->U, &n_int, Uy, &int_one, &real_one, y_UUy, &int_one );
   }

//covariates is an [n*ncovar] matrix holding the covariates, including the bias
void FastLMM_lowrank::set_covariates(real *covar,size_t ncovariates){
   this->ncovar=ncovariates;                 //the number of covariates (including the fixed effects)
   if(Ucovar!=0){ delete[] Ucovar; }
   if(covar_UUcovar!=0){ delete[] covar_UUcovar; }
   
   Ucovar=new real[this->ncovar*this->rankU];    //The covariates transformed to the low rank space
   
   BLAS_INT n_int =(BLAS_INT)this->n;            //Warning: size_t BLAS_INT conversion
   BLAS_INT ncovar_int=(BLAS_INT)this->ncovar;   //Warning: size_t BLAS_INT conversion
   BLAS_INT rankU_int =(BLAS_INT)this->rankU;    //Warning: size_t BLAS_INT conversion
   BlasLapack::gemm(&char_T,  &char_N, &rankU_int, &ncovar_int, &n_int, &real_one, eigen->U, &n_int, covar, &n_int, &real_zero, Ucovar, &rankU_int );//UX=U.T*SNPs
   
   size_t sizeCovar=this->ncovar*this->n; 
   this->covar_UUcovar=new real[sizeCovar];
#if 1
   for(size_t covar_i=0; covar_i<sizeCovar; covar_i+=this->n)
      {
      BlasLapack::copy(&n_int, &covar[covar_i], &int_one, &covar_UUcovar[covar_i], &int_one);
      }
#else
   BLAS_INT sizeCovar_int=(BLAS_INT)(ncovar*n); //Warning: size_t BLAS_INT conversion product
   BlasLapack::copy(&sizeCovar_int, covar, &int_one, covar_UUcovar, &int_one); //Initialize with the original covariates
#endif
   BlasLapack::gemm(&char_N,  &char_N, &n_int, &ncovar_int, &rankU_int, &real_minusone, eigen->U, &n_int, Ucovar, &rankU_int, &real_one, covar_UUcovar, &n_int );//covar_UUcovar=covar-U*Ucovar
}

//find ML parameters for the null model
void FastLMM_lowrank::trainNull()
   {   
   real UXdummy = 0.0;        //As We don't have a SNP, we pass a dummy variable

   //find optimal delta:
   if ( logDeltaPreset == logDeltaPreset )
      {
      delta_0[ 0 ] = logDeltaPreset;
      }
   else
      {
      if ( FastLMM::brentExhaustive )
         {
         //Find the optimal delta
         delta_0[ 0 ] = FastLMM_lowrank::optimize_delta_Brentintervals( 0, (BLAS_INT)ncovar, (BLAS_INT)n, (BLAS_INT)rankU, Uy, &UXdummy, y_UUy, &UXdummy, eigen->S, Ucovar, covar_UUcovar, &UXdummy, &UXdummy, 0 );
         }
      else
         {
          //Find the optimal delta
         delta_0[ 0 ] = FastLMM_lowrank::optimize_delta_BrentintervalsCheckTriplets( 0, (BLAS_INT)ncovar, (BLAS_INT)n, (BLAS_INT)rankU, Uy, &UXdummy, y_UUy, &UXdummy, eigen->S, Ucovar, covar_UUcovar, &UXdummy, &UXdummy, 0 );
         }
      }

   //get the ML parameters for the optimal delta:
   LL_0[0]= (real)-1.0 * nLLeval( nullptr, delta_0[0], 0,(BLAS_INT)ncovar, (BLAS_INT)n, (BLAS_INT)rankU, Uy, &UXdummy, y_UUy, &UXdummy, eigen->S, Ucovar, covar_UUcovar, &sigg2_0[0], &sige2_0[0], &beta_0[0], &rank_beta_0[0], true, &UXdummy, &UXdummy, 0 );//evaluate the LL at delta_0 and get the ML parameters   //Warning: size_t BLAS_INT conversion
   if ( verboseLevel >= verboseDebug )
      {
      printf( "\n------------------------------------------" );
      printf( "\nNULL: LL:%f", LL_0[0] );
      printf( "\n  sigg2=%f sige2=%f delta=%f", sigg2_0[0], sige2_0[0], delta_0[0] );
      printf( " beta: ");
      for( size_t i=0; i < ncovar; ++i )
         {
         printf( "%f ", beta_0[ i ] );
         }
      printf( "\n------------------------------------------" );
      }

   //copy the "default" null model parameters
   for ( size_t iSnp = 1; iSnp < nSNPs; ++iSnp )
      {
      delta_0[ iSnp ]     = delta_0[ 0 ];
      rank_beta_0[ iSnp ] = rank_beta_0[ 0 ];
      LL_0[ iSnp ]        = LL_0[ 0 ];
      sigg2_0[ iSnp ]     = sigg2_0[ 0 ];
      sige2_0[ iSnp ]     = sige2_0[ 0 ];
      for ( size_t iCovariates = 0; iCovariates < ncovar; ++iCovariates )
         {
         beta_0[ (iSnp*ncovar) + iCovariates ] = beta_0[ iCovariates ];
         }
      }
   }

//find ML parameters for all SNPs
void FastLMM_lowrank::trainAlt(){      
   //alternative models:
   size_t d = 1;                //assume that always one SNP is being tested
   size_t d_ = d + ncovar;      //the number of weights (dimensionality of a single beta)
   
   #pragma omp parallel for //This does not work in the exclude setting as it messes up the exclude SNP setting
   for (int snp=0; snp<(int)nSNPs;++snp){//train the alternative model for the SNP  (parallel for loop) //Warning: size_t int conversion
      real* UWcurrent      = nullptr;
      real* W_UUWcurrent   = nullptr;
         
      if(nExclude[snp]>0)
         {
         UWcurrent      = UW     + (this->rankU * this->excludeSnpStart[snp]);
         W_UUWcurrent   = W_UUW  + (this->n     * this->excludeSnpStart[snp]);
         }

      //set pointers to the current SNP and ML parameters:
      real* currentUX=UX+snp*(d*rankU);
      real* currentX_UUX=X_UUX+snp*(n*d);
      real* currentBeta=beta+snp*d_;
      real* currentSigg2=sigg2+snp;
      real* currentSige2=sige2+snp;
      real* currentDelta=delta+snp;
      BLAS_INT* currentRank_beta=rank_beta+snp;
      
      if (nExclude[snp]>0 && !use_f_test && (covLearnType == REFIT_EXCLUDE || covLearnType == FULL))
            {
            //update null delta
            if(FastLMM::brentExhaustive)
              {
              delta_0[snp]=FastLMM_lowrank::optimize_delta_Brentintervals((BLAS_INT)(d-1), (BLAS_INT)ncovar, (BLAS_INT)n,(BLAS_INT)rankU, Uy,currentUX,y_UUy,currentX_UUX, eigen->S, Ucovar,covar_UUcovar,UWcurrent, W_UUWcurrent,(BLAS_INT)nExclude[snp]); //Find the optimal delta   //Warning: size_t BLAS_INT conversion
              }
            else
              {
              delta_0[snp]=FastLMM_lowrank::optimize_delta_BrentintervalsCheckTriplets((BLAS_INT)(d-1), (BLAS_INT)ncovar, (BLAS_INT)n,(BLAS_INT)rankU, Uy,currentUX,y_UUy,currentX_UUX, eigen->S, Ucovar,covar_UUcovar,UWcurrent, W_UUWcurrent,(BLAS_INT)nExclude[snp]); //Find the optimal delta   //Warning: size_t BLAS_INT conversion
              }
            }
         if(nExclude[snp]>0)
            {
            //update the null likelihood
            LL_0[snp]=(real)-1.0*nLLeval(nullptr, delta_0[snp], (BLAS_INT)(d-1),(BLAS_INT)ncovar, (BLAS_INT)n,(BLAS_INT)rankU,Uy, currentUX,y_UUy,currentX_UUX, eigen->S, Ucovar,covar_UUcovar, &sigg2_0[snp], &sige2_0[snp], &beta_0[((size_t)snp)*ncovar], &rank_beta_0[snp], true,UWcurrent,W_UUWcurrent,(BLAS_INT)this->nExclude[snp]);//evaluate the LL at delta_0 and get the ML parameters   //Warning: size_t BLAS_INT conversion
            }
         if (covLearnType == FULL || (nExclude[snp]>0 && use_f_test && covLearnType == REFIT_EXCLUDE))
            {
            //update alternative delta
            if(FastLMM::brentExhaustive)
               {
               *currentDelta=FastLMM_lowrank::optimize_delta_Brentintervals((BLAS_INT)d, (BLAS_INT)ncovar, (BLAS_INT)n,(BLAS_INT)rankU, Uy,currentUX,y_UUy, currentX_UUX, eigen->S, Ucovar, covar_UUcovar, UWcurrent, W_UUWcurrent, (BLAS_INT)this->nExclude[snp] ); //Find the optimal delta  //Warning: size_t BLAS_INT conversion
               }
            else
               {
               *currentDelta=FastLMM_lowrank::optimize_delta_BrentintervalsCheckTriplets((BLAS_INT)d, (BLAS_INT)ncovar, (BLAS_INT)n,(BLAS_INT)rankU, Uy,currentUX,y_UUy, currentX_UUX, eigen->S, Ucovar, covar_UUcovar, UWcurrent, W_UUWcurrent, (BLAS_INT)this->nExclude[snp] ); //Find the optimal delta   //Warning: size_t BLAS_INT conversion
               }
            }
         else if((covLearnType == ONCE) || (nExclude[snp]==0  && covLearnType == REFIT_EXCLUDE) || (nExclude[snp]>0  && covLearnType == REFIT_EXCLUDE && !use_f_test))
            {
            //use existing null-delta
            *currentDelta=delta_0[snp];
            }
         else
            {
            Fatal("Illegal covLearnType");
            }

      //evaluate the LL and find the parameters
      LL[snp]=(real)-1.0*nLLeval(f_stat+(snp*d_), *currentDelta, (BLAS_INT)d,(BLAS_INT)ncovar, (BLAS_INT)n,(BLAS_INT)rankU,Uy, currentUX,y_UUy,currentX_UUX, eigen->S, Ucovar,covar_UUcovar,currentSigg2 ,currentSige2, currentBeta, currentRank_beta, true, UWcurrent, W_UUWcurrent, (BLAS_INT)this->nExclude[snp]);//evaluate the LL at delta_0 and get the ML parameters   //Warning: size_t BLAS_INT conversion
      
      //Compute p-value from LR test:
      real df=*currentRank_beta-(real)rank_beta_0[snp];
      //check the degrees of freedom for validity
      if ( df > 0.0 )
         {    
         if ( use_f_test )
            {
            //Fishers F-test
            //
            // PvalUsingFisherFTestCdf focuses on keeping precision 
            //    for very small pvalues by inverting on the reverse 
            //    path of the Beta::BetaInc() branch
            pvals[ snp ] = PvalUsingFisherFTestCdf( f_stat[ snp * d_ ], df, (real)(n - *currentRank_beta) );
            //pvals[snp] = 1.0 - FisherF::Cdf( f_stat[ snp * d_ ], df, (real)(n - *currentRank_beta) );
            }
         else
            { //Likelihood-ratio test
            pvals[ snp ] = PvalFromLikelihoodRatioTest( LL[snp] - LL_0[ snp ], (real)0.5 * df );
            }
         }
      else
         {
         //This happens for example if the SNP is a constant. In this case df=0
         pvals[snp]=std::numeric_limits<real>::quiet_NaN();
         }
      if ( verboseLevel >= verboseDebug ){
         printf( "\n- LL[%i]:%f p=%f", snp, LL[ snp ], pvals[ snp ] );
         printf( "\n  sigg2=%f sige2=%f delta=%f", *currentSigg2, *currentSige2, *currentDelta );
         printf( " beta:" );
         for( size_t i = 0; i < d_; ++i )
            {
            printf( "%f ", currentBeta[ i ] );
            }
      }
   }
}

void FastLMM_lowrank::init()
   {
   if (beta_0 != nullptr)
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

   size_t d = 1;              // assume that always one SNP is being tested
   size_t d_ = d + ncovar;    // the number of weights (dimensionality of a single beta)
   
   // if the parameters have already been set in a previous iteration delete them.
   if ( beta   != nullptr ) { delete[] beta;   beta   = nullptr; }
   if ( delta  != nullptr ) { delete[] delta;  delta  = nullptr; }
   if ( sigg2  != nullptr ) { delete[] sigg2;  sigg2  = nullptr; }
   if ( sige2  != nullptr ) { delete[] sige2;  sige2  = nullptr; }
   if ( LL     != nullptr ) { delete[] LL;     LL     = nullptr; }
   if ( pvals  != nullptr ) { delete[] pvals;  pvals  = nullptr; }
   if ( f_stat != nullptr ) { delete[] f_stat; f_stat = nullptr; }
   if ( rank_beta != nullptr ) { delete[] rank_beta; rank_beta = nullptr; }
   
   f_stat = new real[nSNPs*d_];
   beta=new real[nSNPs*d_];       //create all the variables
   delta=new real[nSNPs];
   sigg2=new real[nSNPs];
   sige2=new real[nSNPs];
   LL=new real[nSNPs];
   pvals=new real[nSNPs];
   rank_beta=new BLAS_INT[nSNPs];
   }

//static function:
//finds delta by optimizing over delta intervals
real FastLMM_lowrank::optimize_delta_BrentintervalsCheckTriplets( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n,const BLAS_INT rankU, const real *Uy, const real *currentUX, const real *y_UUy, const real *currentX_UUX,const real *S, const real *Ucovar, const real *covar_UUcovar, const real *UW, const real *W_UUW, const BLAS_INT nExclude )
   {
   real *negLL_interval=new real[this->numberIntervals+1];
   real *log_delta_interval=new real[this->numberIntervals+1];
   
   real best_log_delta=std::numeric_limits<real>::quiet_NaN();
   real best_negLL=std::numeric_limits<real>::infinity();

   real MLparamsdummy=0.0;
   BLAS_INT rankDummy=0;
   size_t nevals=0;

   //evaluate the boundaries of all intervals
   #pragma omp parallel for
   for(int i=0;i<=(int)this->numberIntervals;++i){
      log_delta_interval[i]=this->minLogVal+i*(this->maxLogVal-this->minLogVal)/this->numberIntervals;
      negLL_interval[i]=FastLMM_lowrank::nLLeval(nullptr, log_delta_interval[i],d,ncovar,n,rankU,Uy,currentUX,y_UUy,currentX_UUX,S,Ucovar,covar_UUcovar,&MLparamsdummy,&MLparamsdummy,&MLparamsdummy,&rankDummy,false,UW,W_UUW,nExclude);
      ++nevals;
      if (negLL_interval[i]<=best_negLL){
            best_negLL=negLL_interval[i];
            best_log_delta=log_delta_interval[i];
      }
   }
   
   brentFunctor_FastLMM_lowrank func=brentFunctor_FastLMM_lowrank(d, ncovar, n,rankU, Uy,currentUX,y_UUy,currentX_UUX, S, Ucovar,covar_UUcovar, UW, W_UUW, nExclude);
   #pragma omp parallel for  //might be slow if parallelized as most threads finish quickly
   for(int i=1;i<(int)this->numberIntervals;++i){
      if (negLL_interval[i]<negLL_interval[i+1] && negLL_interval[i]<negLL_interval[i-1]){//check wether a local optimum exists in this triplet
         real current_negLL=std::numeric_limits<real>::infinity();
         size_t nevals_Int=0;
         //If there is a local optimum, optimize over the current triplet:
         real current_log_delta=BrentC::minimize(func,log_delta_interval[i-1],log_delta_interval[i+1],this->brentTol,current_negLL,nevals_Int,this->brentMaxIter,true);
         nevals+=nevals_Int;
         if (current_negLL<=best_negLL){
            best_negLL=current_negLL;
            best_log_delta=(current_log_delta);
         }
      }
   }
   delete []log_delta_interval;
   delete []negLL_interval;

#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
   const char *szFmt1 = "\nnevals: %Iu logdelta: %f\n";
#else                      // Linux/g++ uses a %zu specifier for size_t
   const char *szFmt1 = "\nnevals: %zu logdelta: %f\n";
#endif
   Verbose( verboseDebug, szFmt1, nevals, best_log_delta );
   return( best_log_delta );
   }

//static function:
//finds delta by optimizing over log delta intervals
real FastLMM_lowrank::optimize_delta_Brentintervals(const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n,const BLAS_INT rankU, const real *Uy, const real *currentUX, const real *y_UUy, const real *currentX_UUX,const real *S, const real *Ucovar, const real *covar_UUcovar, const real *UW, const real *W_UUW, const BLAS_INT nExclude){
   real best_log_delta = std::numeric_limits<real>::quiet_NaN();
   real best_negLL = std::numeric_limits<real>::infinity();
   
   //create the Functor. This will be passed to the optmizer
   brentFunctor_FastLMM_lowrank func = brentFunctor_FastLMM_lowrank( d, ncovar, n, rankU, Uy, currentUX, y_UUy, currentX_UUX, S, Ucovar, covar_UUcovar, UW, W_UUW, nExclude );
   
   size_t nevals=0;  //count the total number of LL-function evaluations

   //Perform a Brent optimization on each interval
   #pragma omp parallel for
   for( int i=0; i<(int)this->numberIntervals; ++i )
      {
      // Set the mimimum and maximum of the current Brent interval
      real current_log_delta_min = this->minLogVal + i * (this->maxLogVal - this->minLogVal) / this->numberIntervals;
      real current_log_delta_max = this->minLogVal + (i+1) * (this->maxLogVal - this->minLogVal) / this->numberIntervals;
      real current_negLL = std::numeric_limits<real>::infinity();       // Setup for the optimal Negative Log Likelihood
      size_t nevals_Int = 0;        //The number of evaluations on the current interval
      real current_log_delta = BrentC::minimize( func, current_log_delta_min, current_log_delta_max, this->brentTol, current_negLL, nevals_Int, this->brentMaxIter, true );
      
      // Increment the total number of evaluations by #evaluations on the current one
      nevals += nevals_Int;
      
      // Save the 'better' optimum as we iterate on the current interval
      if ( current_negLL <= best_negLL )
         {
         best_negLL = current_negLL;
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


//static function:
//cleaned up version of nLLeval() for the full rank model
real FastLMM_lowrank::nLLeval(real *f_statistic, const real logdelta,const BLAS_INT d, const BLAS_INT ncovar,const BLAS_INT n, const BLAS_INT rankU,const real *Uy, const real *UX,const real *y_UUy, const real *X_UUX, const real *S, const real *Ucovar, const real *covar_UUcovar, real *sigg2, real *sige2, real *beta, BLAS_INT *rank_beta, const bool returnMLparams, const real *UW, const real *W_UUW, BLAS_INT nExclude)
   {
   real delta = exp( logdelta );
   
   BLAS_INT n_ = n + rankU;               //we regress over both, the n residuals and the rankU latent individuals
   BLAS_INT d_ = d + ncovar;              //the number of parameters to be estimated is number of SNP coefficients + number of covariates
   size_t nd = (size_t)n_* (size_t)d_ ;
   real *res = new real[ n_ ];            //the residual vector. Will be initialized to ySdi, later overwritten by the residual of the regression (deleted)
   real *Sd = new real[ rankU ];          //holds the sqrt of the diagonal values of (S+delta*I) (deleted)
   real *beta_temp = new real[ d_ ];      //the parameters. (deleted)
   real *XSdi = new real[ nd ];           //XSdi This one will be used to compute the residual in the end. (deleted)
   bool delete_f_stat = false;
   real yy = 0.0;

   if ( f_statistic == nullptr )
      {
      f_statistic = new real[ d_ ];
      delete_f_stat = true;
      }

   real *Sd_ = Sd;
   real *res_ = res;
   const real *S_ = S;
   const real *Uy_ = Uy;

   //--------------------------------------------------
   //compute the determinant and ySdi
   real delta_sqrt;
   real det;

   if ( useHeritability )
      {
      delta_sqrt = (real)1.0 / sqrt( 1.0 - logdelta );
      det = ( n - rankU ) * log( 1.0 - logdelta );
      for ( size_t i=0; (BLAS_INT)i < rankU; ++i )
         {
         Sd[ i ] = (logdelta * S[i]) + ( 1.0 - logdelta );
         }
      }
   else
      {
      delta_sqrt = (real)1.0 / sqrt( delta );
      det = ( n - rankU ) * log( delta );       // the log determinant is the sum of the diagonal of (S+delta*I)
      for ( size_t i=0; (BLAS_INT)i < rankU; ++i )
         {
         Sd[ i ] = S[ i ] + delta;
         }
      }

   BlasLapack::vLn( rankU, Sd, res );           //res holds the Ln(S+deltaI)
   BlasLapack::vSqrt( rankU, Sd, Sd );

   for ( size_t i=0; (BLAS_INT)i < rankU; ++i )
      {
      det += *res_;
      *res_ = *Uy_ / *Sd_;
      yy += ((*res_) * (*res_));
      ++Sd_;
      ++res_;
      ++Uy_;
      ++S_;
      }

   BlasLapack::copy( &n, y_UUy, &int_one, res_, &int_one );
   BlasLapack::scal( &n, &delta_sqrt, res_, &int_one );
   for ( size_t i=0; (BLAS_INT)i < n; ++i )
      {
      yy += res_[ i ] * res_[ i ];
      }

   real *XSdi_ind = XSdi;
   const real *UX_ind = UX;
   real *Sd_ind = Sd;
   const real *X_UUX_ind = X_UUX;
   for ( size_t snp=0; (BLAS_INT)snp < d; ++snp )
      {
      Sd_ind = Sd;
      for ( size_t ind=0; (BLAS_INT)ind < rankU; ++ind )
         {
         *XSdi_ind = *UX_ind / *Sd_ind;
         ++UX_ind;
         ++XSdi_ind;
         ++Sd_ind;
         }
      
      BlasLapack::copy( &n, X_UUX_ind, &int_one, XSdi_ind, &int_one );
      BlasLapack::scal( &n, &delta_sqrt, XSdi_ind, &int_one );
      X_UUX_ind += n;
      XSdi_ind += n;
      }
   
   const real *covar_UUcovar_ind = covar_UUcovar;
   size_t ind1 = 0;

   for ( size_t cov=0; (BLAS_INT)cov < ncovar; ++cov )
      {
      for ( size_t ind=0; (BLAS_INT)ind < rankU; ++ind )
         {
         *XSdi_ind = Ucovar[ ind1 ] / Sd[ ind ];
         ++XSdi_ind;
         ++ind1;
         }
      
      BlasLapack::copy( &n, covar_UUcovar_ind, &int_one, XSdi_ind, &int_one );
      BlasLapack::scal( &n, &delta_sqrt, XSdi_ind, &int_one );
      covar_UUcovar_ind += n;
      XSdi_ind += n;
      }

   real *XX = new real[ d_ * d_ ];     //compute XX = XSdi.T*XSdi
   if ( d_ > 0 )
      {
      BlasLapack::syrk( &char_L, &char_T, &d_, &n_, &real_one, XSdi, &n_, &real_zero, XX, &d_ );
      }

   real *Xy = new real[ d_ ];          //compute Xy = XSdi.T*y
   BlasLapack::gemv( &char_T, &n_, &d_, &real_one, XSdi, &n_, res, &int_one, &real_zero, Xy, &int_one );

   //Exclude SNPs from the RRM in the likelihood evaluation
   if ( nExclude > 0 )
      {
      // W
      real *WSdi = new real[ n_ * nExclude ];   //deleted
      real *WSdi_ind = WSdi;
      const real *UW_ind = UW;
      const real *W_UUW_ind = W_UUW;
      for ( size_t snpC = 0; snpC < (size_t)nExclude; ++snpC )
         {
         Sd_ind = Sd;
         for ( size_t ind=0; (BLAS_INT)ind < rankU; ++ind )
            {
            *WSdi_ind = *UW_ind / *Sd_ind;
            ++UW_ind;
            ++WSdi_ind;
            ++Sd_ind;
            }

         BlasLapack::copy( &n, W_UUW_ind, &int_one, WSdi_ind, &int_one );
         BlasLapack::scal( &n, &delta_sqrt, WSdi_ind, &int_one );
         W_UUW_ind += n;
         WSdi_ind += n;
         }

      //compute (I-WW)
      size_t nExclude2 = (size_t)nExclude * (size_t)nExclude;
      real *WW = new real[ nExclude2 ];
      //subtract WW
      BlasLapack::gemm( &char_T, &char_N, &nExclude, &nExclude, &n_, &real_minusone, WSdi, &n_, WSdi, &n_, &real_zero, WW, &nExclude );
      //add the identity
      for( size_t snp = 0; (BLAS_INT)snp < nExclude; ++snp )
         {
         ++WW[ (snp * ((size_t)nExclude)) + snp ];
         }

      real *WX = new real[ (size_t)nExclude * (size_t)d_ ];
      real *Wy = new real[ nExclude ];

      if ( d_ > 0 )
         {
         BlasLapack::gemm( &char_T, &char_N, &nExclude, &d_, &n_, &real_one, WSdi, &n_, XSdi, &n_, &real_zero, WX, &nExclude );
         }

      BlasLapack::gemv( &char_T, &n_, &nExclude, &real_one, WSdi, &n_, res, &int_one, &real_zero, Wy, &int_one );

      //compute the inverse and determinant of WW:
      //1) eigen-decomposition of (I-WW)
      real *U_WW = new real[ nExclude2 ];    //Warning: BLAS_INT product
      real *S_WW = new real[ nExclude ];
      EigenSym::EV_symmetric( nExclude, WW, U_WW, S_WW );
      //2) determinant and solves using eigen-decomposition
      
      // U_W^T * Wy
      real *UWy = new real[ nExclude ];
      BlasLapack::gemv( &char_T, &nExclude, &nExclude ,&real_one, U_WW, &nExclude, Wy, &int_one, &real_zero, UWy, &int_one );
      // U_W^T * WX 
      real *UWX = new real[ (size_t)nExclude * d_ ];
      // UX = U.T * SNPs
      BlasLapack::gemm( &char_T,  &char_N, &nExclude, &d_, &nExclude, &real_one, U_WW, &nExclude, WX, &nExclude, &real_zero, UWX, &nExclude );

      //compute S_WW^{-1} * UWX
      size_t col = 0;
      for ( size_t snp = 0; snp < (size_t)d_; ++snp )
         {
         for ( size_t i_ex = 0; i_ex<(size_t)nExclude; ++i_ex )
            {
            //if (S_WW[i_ex]>MathFunctions::eps_rank || S_WW[i_ex]<-MathFunctions::eps_rank)
               {
               WX[ col + i_ex ] = UWX[col + i_ex] / S_WW[ i_ex ];
               }
            //else
            //   {
            //   WX[ col + i_ex ] = (real)0.0;
            //   }
            }
            col += (size_t)nExclude;
         }

      //compute S_WW^{-1} * UWy and the determinant update
      for ( size_t i_ex = 0; i_ex < (size_t)nExclude; ++i_ex )
         {
         //if (S_WW[i_ex]>MathFunctions::eps_rank || S_WW[i_ex]<-MathFunctions::eps_rank)
            {
            Wy[ i_ex ] = UWy[ i_ex ] / S_WW[ i_ex ];
            //update the logdet term in the likelihood
            det += log( S_WW[ i_ex ] );
            //update yy += yWWy
            yy += Wy[ i_ex ] * UWy[ i_ex ];
            //yy_up += Wy[ i_ex ] * UWy[ i_ex ];
            }
         //else
         //   {
         //   Wy[i_ex] = (real)0.0;
         //   }
         }

      if ( d_ > 0 )
         {
         //update XX += XWWX
         BlasLapack::gemm( &char_T, &char_N, &d_, &d_, &nExclude, &real_one, UWX, &nExclude, WX, &nExclude, &real_one, XX, &d_ );
         //update Xy += XWWy
         BlasLapack::gemv( &char_T, &nExclude, &d_, &real_one, UWX, &nExclude, Wy, &int_one, &real_one, Xy, &int_one );
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

   delete[] Sd;      //Sdi is no longer needed
   delete[] XSdi;    //XSdi is no longer needed.
   delete[] res;     //ySdi is no longer needed

   //This part should now be the same for all branches, except the sige2/sigg2 issue:

   //solve the least squares problem.
   real logDetXX = 0.0;
   FastLMM::solveLstSqSym( d_, XX, Xy, beta_temp, rank_beta, f_statistic, &logDetXX );

   real sigg2_temp = yy;      //the ML sigma_g^2
   for (size_t dim=0; (BLAS_INT)dim < d_; ++dim )
      {
      sigg2_temp -= Xy[ dim ] * beta_temp[ dim ];
      }

   real negLL;    // negative Log Likelihood Temp variable

   if ( mlType == ML )
      {
      //compute the negative LL using the ML parameters at delta
      sigg2_temp /= n;
      negLL = (real)0.5 * ((n * (_log2pi + log( sigg2_temp ) + (real)1.0 )) + det );
      }
   else if ( mlType == REML )
      {
      //compute the negative LL using the REML parameters at delta
      sigg2_temp /= n-d_;
      negLL = (real)0.5 * ( (( n - d_ ) * (_log2pi + log( sigg2_temp ) + (real)1.0 )) + det + logDetXX );
      }
   else
      {
      Fatal( "Only ML or REML are valid ML-types" );
      }
   
   if ( returnMLparams )
      {
      //only return the ML params if needed. (e.g. no need during Brents search) This avoids danger of Memory errors
      if ( useHeritability )
         {
         *sigg2 = logdelta * sigg2_temp;
         *sige2 = (1.0-logdelta) * sigg2_temp;
         }
      else
         {
         *sigg2 = sigg2_temp;
         *sige2 = delta**sigg2;
         }

      BlasLapack::copy( &d_, beta_temp, &int_one, beta, &int_one );
      }
   
   delete[] beta_temp;  //delete the ML beta_temp
   
   if ( delete_f_stat )
      {
      delete[] f_statistic;
      }
   else
      {
      for (size_t dim = 0;(BLAS_INT)dim<d_;++dim)
         {
         f_statistic[dim] /= sigg2_temp;
         }
      }

   //printf("negLL:%f    sigg2:%f    sige2:%f\n",negLL,sigg2_temp,sigg2_temp*delta);
   return negLL;
   }

//Functor holdong all the needed members of FastLMM_lowrank to be handed to the Brent's routine
//the constructor gets the members of FastLMM_lowrank
brentFunctor_FastLMM_lowrank::brentFunctor_FastLMM_lowrank(const BLAS_INT d_, const BLAS_INT ncovar_, const BLAS_INT n_,const BLAS_INT rankU_, const real *Uy_,const  real *UX_, const real *y_UUy_,const  real *X_UUX_, const real *S_, const real *Ucovar_,const real *covar_UUcovar_, const real *UW, const real *W_UUW, const BLAS_INT nExclude){
   //member parameters from FastLMM_lowrank:
   d=d_;
   rankU=rankU_;
   ncovar=ncovar_;
   n=n_;
   Uy=Uy_;
   UX=UX_;
   y_UUy=y_UUy_;
   X_UUX=X_UUX_;
   covar_UUcovar=covar_UUcovar_;
   this->nExclude = nExclude;
   this->UW = UW;
   this->W_UUW = W_UUW;
   S=S_;
   Ucovar=Ucovar_;
}

//destructor has nothing to do
brentFunctor_FastLMM_lowrank::~brentFunctor_FastLMM_lowrank(){}

//the operator evaluates nLLeval
real brentFunctor_FastLMM_lowrank::operator()(const real logdelta){
   
   //the rest are dummies:
   BLAS_INT rank_beta = 0;
   real dummy_MLparams=0.0;

   return FastLMM_lowrank::nLLeval(nullptr, logdelta, d, ncovar, n, rankU, Uy, UX, y_UUy, X_UUX, S, Ucovar, covar_UUcovar, &dummy_MLparams, &dummy_MLparams, &dummy_MLparams, &rank_beta, false, UW, W_UUW, nExclude);
}

