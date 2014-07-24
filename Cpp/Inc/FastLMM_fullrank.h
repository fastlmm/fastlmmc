#if !defined( FastLMM_fullrank_h )
#define FastLMM_fullrank_h
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
 * this class does standard FastLMM.
 */
class FastLMM_fullrank: public FastLMM
   {
   real *Uy;
   real *UX;
   real *UW;
   real *Ucovar;

   real optimize_delta_Brentintervals( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const real *Uy, const real *currentUX, const real *S, const real *Ucovar, const real *UW, const BLAS_INT nExclude );
   real optimize_delta_BrentintervalsCheckTriplets( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const real *Uy, const real *currentUX, const real *S, const real *Ucovar, const real *UW, const BLAS_INT nExclude );

public:
   FastLMM_fullrank( EigenSym *eigen );
   ~FastLMM_fullrank();

   //void set_Kernel( real *K );               //K is an [n by n] symmetric (lower triangle has to be stored) positve semi-definite matrix. This function performs the SVD  of K. K will be deleted.
   //void set_Kernel( real *U_, real *S_ );  //U [N by n] and S [n by 1] are the Eigenvectors and Eigenvalues of the symmetric Kernel matrix K.
   void set_SNPs( real *SNPs, size_t nSNPs_, vector<SnpInfo>& snpInfo_ );  //SNPs [n by nSNPs] holds the nSNPs SNPs to be tested for each of n individuals
   void set_Phenotype( real *phenotype, string& phenLabel );    //phenotype is a vector of length n, holding the (single) phenotype.
   void set_covariates( real *covar, size_t ncovariates );   // covariates is an [n*ncovar] matrix holding the covariates, including the bias
   void set_SNPsCov( real *W, size_t nSNPsW, vector<SnpInfo>& snpInfoCov );   // W [n by nSNPs] holds the nSNPsW SNPs to be tested for each of n individuals

   void set_UW( real *pW ) { UW=pW; }        // temporary write access to UW

   void init();
   void trainNull();                         //find ML parameters for the null model
   void trainAlt();                          //find ML parameters for all SNPs
   static real nLLeval( real *f_statistic, const real logdelta, const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const real *Uy, const real *UX, const real *S, const real *Ucovar, real *sigg2, real *sige2, real *beta, BLAS_INT *rank_beta, const bool returnMLparams, const real* UW, const BLAS_INT nExclude );
   };

//Functor passing all the needed members of FastLMM_fullrank to be handed to the Brent's routine
class brentFunctor_FastLMM_fullrank: public BrentFunctor
   {
   const real *Uy;
   const real *UX;
   const real *S;
   const real *Ucovar;
   const real *UW;
   BLAS_INT    d;
   BLAS_INT    ncovar;
   BLAS_INT    n;
   BLAS_INT    nExclude;

public:
   brentFunctor_FastLMM_fullrank( const BLAS_INT d_, const BLAS_INT ncovar_, const BLAS_INT n_, const real *Uy_,const  real *UX_, const real *S_, const real *Ucovar_, const real *UW, const BLAS_INT nExclude );
   ~brentFunctor_FastLMM_fullrank();
   virtual real operator()( const real logdelta );
   };

#endif   //FastLMM_fullrank_h
