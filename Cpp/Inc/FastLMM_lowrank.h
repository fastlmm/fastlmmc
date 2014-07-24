#if !defined(FastLMM_lowrank_h)
#define FastLMM_lowrank_h
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
 * this class does low rank FastLMM.
 */
class FastLMM_lowrank: public FastLMM
   {
   real *Uy;
   real *y_UUy;
   real *UX;
   real *X_UUX;
   real *Ucovar;
   real *covar_UUcovar;
   size_t rankU;
   real *UW;
   real *W_UUW;

   real optimize_delta_Brentintervals( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU, const real *Uy, const real *currentUX, const real *y_UUy, const real *currentX_UUX,const real *S, const real *Ucovar, const real *covar_UUcovar, const real *UW, const real *W_UUW, const BLAS_INT nExclude );
   real optimize_delta_BrentintervalsCheckTriplets( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU, const real *Uy, const real *currentUX, const real *y_UUy, const real *currentX_UUX, const real *S, const real *Ucovar, const real *covar_UUcovar, const real *UW, const real *W_UUW, const BLAS_INT nExclude );

public:
   FastLMM_lowrank( EigenSym *eigen );
   ~FastLMM_lowrank();
   //void set_Kernel(size_t rankK, real *K);//K is an [n by n] symmetric (lower triangle has to be stored) positve semi-definite matrix. This function performs the rankK-SVD  of K. K will be deleted.
   //void set_Kernel_CovFromSNPsample(size_t nSNPsample, real *SNPsample);//SNPsample is an [n by nSNPsample] matrix. This function performs an SVD of SNPsamples to estimate a min(n,nSNPsample)-rank Covariance Kernel.
   //void set_Kernel_CovFromSNPsample2(size_t nSNPsample, real *SNPsample);
   //void set_Kernel( const size_t rankK, real *U_, real *S_ );  //U [n by rankK] and S [rankK by 1] are the Eigenvectors and Eigenvalues of the symmetric rankK-rank Kernel matrix K.
   void set_SNPs( real *SNPs, size_t nSNPs_, vector<SnpInfo>& snpInfo_ );     // SNPs [n by nSNPs] holds the nSNPs SNPs to be tested for each of n individuals
   void set_Phenotype( real *phenotype, string& phenLabel );   //phenotype is a vector of length n, holding the (single) phenotype.
   void set_covariates( real *covar,size_t ncovariates );      //covariates is an [n*ncovar] matrix holding the covariates, including the bias
   void set_SNPsCov( real *W, size_t nSNPsW, vector<SnpInfo>& snpInfoCov );   //W [n by nSNPs] holds the nSNPsW SNPs to be tested for each of n individuals
   void trainNull();     //find ML parameters for the null model
   void trainAlt();      //find ML parameters for all SNPs
   void init();
   static real nLLeval(real *f_statistic, const real logdelta,const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU,const real *Uy, const real *UX,const real *y_UUy, const real *X_UUX, const real *S, const real *Ucovar, const real *covar_UUcovar, real *sigg2, real *sige2, real *beta, BLAS_INT *rank_beta, const bool returnMLparams, const real* UW, const real* W_UUW, const BLAS_INT nExclude);
   };

//Functor holdong all the needed members of FastLMM_lowrank to be handed to the Brent's routine
class brentFunctor_FastLMM_lowrank: public BrentFunctor
   {
   const real *Uy;
   const real *UX;
   const real *X_UUX;
   const real *y_UUy;
   const real *covar_UUcovar;
   BLAS_INT    rankU;
   const real *S;
   const real *Ucovar;
   const real *UW;
   const real *W_UUW;
   BLAS_INT    d;
   BLAS_INT    ncovar;
   BLAS_INT    n;
   BLAS_INT    nExclude;

public:
   brentFunctor_FastLMM_lowrank( const BLAS_INT d_, const BLAS_INT ncovar_, const BLAS_INT n_,const BLAS_INT rankU_, const real *Uy_,const  real *UX_, const real *y_UUy_,const  real *X_UUX_, const real *S_, const real *Ucovar_,const real *covar_UUcovar_, const real* UW, const real* W_UUW, const BLAS_INT nExclude );
   ~brentFunctor_FastLMM_lowrank();
   virtual real operator()( const real logdelta );
   };

#endif   //FastLMM_lowrank
