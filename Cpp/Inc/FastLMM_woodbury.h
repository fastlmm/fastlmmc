#if !defined(FastLMM_woodbury_h)
#define FastLMM_woodbury_h
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
class FastLMM_woodbury: public FastLMM
   {
   real *SU_y;
   real *y;
   real *SU_x;
   real *X;
   real *SUcovar;
   real *covar;
   size_t rankU;

//   real optimize_delta_Brentintervals( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU, const real *SU_y, const real *currentSU_x, const real *y, const real *currentX, const real *S, const real *SUcovar, const real *covar );
   real optimize_delta_Brentintervals( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU, const real *SU_y, const real *currentSU_x, const real *y, const real *currentX, const real *S, const real *SUcovar, const real *covar, const int numberIntervals, const real minLogVal, const real maxLogVal, const real brentTol );
   real optimize_delta_BrentintervalsCheckTriplets( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU, const real *SU_y, const real *currentSU_x, const real *y, const real *currentX, const real *S, const real *SUcovar, const real *covar );

public:
   FastLMM_woodbury( EigenSym *eigen );
   ~FastLMM_woodbury();

   // K is an [n by n] symmetric (lower triangle has to be stored) positve semi-definite matrix. 
   // This function performs the rankK-SVD  of K. K will be deleted.
   //void set_Kernel(size_t rankK, real *K);

   // SNPsample is an [n by nSNPsample] matrix. 
   // This function performs an SVD of SNPsamples to estimate a min( n, nSNPsample )-rank Covariance Kernel.
   //void set_Kernel_CovFromSNPsample( size_t nSNPsample, real *SNPsample );
   //void set_Kernel_CovFromSNPsample2( size_t nSNPsample, real *SNPsample );

   // U [n by rankK] and S [rankK by 1] are the Eigenvectors and Eigenvalues of the symmetric rankK-rank Kernel matrix K.
   //void set_Kernel(const size_t rankK, real *U_, real *S_);

   // SNPs [n by nSNPs] holds the nSNPs SNPs to be tested for each of n individuals
   void set_SNPs( real *SNPs, size_t nSNPs_, vector<SnpInfo>& snpInfo );

   // phenotype is a vector of length n, holding the (single) phenotype.
   void set_Phenotype( real *phenotype, const string& phenLabel );

   // covariates is an [n*ncovar] matrix holding the covariates, including the bias
   void set_covariates( real *covar, size_t ncovariates );

   void trainNull();     //find ML parameters for the null model
   void trainAlt();      //find ML parameters for all SNPs
   void init();
   static real nLLeval( real *f_statistic, const real delta, const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU, const real *Uy, const real *UX, const real *y_UUy, const real *X_UUX, const real *S, const real *Ucovar, const real *covar_UUcovar, real *sigg2, real *sige2, real *beta, BLAS_INT *rank_beta, const bool returnMLparams );

   void TrainNullAndAltForEpistatic();
//   void TrainNullAndAltForEpistaticWithSnpPairs( FILE *pFile, size_t idxSnp1, size_t idxSnp2 );
//   size_t TrainNullAndAltForEpistaticWithSnpPairs( FILE *pFile, size_t iSnp1, size_t iSnp2, real* rgCovariatesTemp, real* rgSUcovarTmp );
   size_t TrainNullAndAltForEpistaticWithSnpPairs( FILE *pFile, size_t iSnp1Base, size_t iSnp1, size_t iSnp2Base, size_t iSnp2, real* rgCovariatesTemp, real* rgSUcovarTmp, real* rgTemp );
   };

//Functor holdong all the needed members of FastLMM_woodbury to be handed to the Brent's routine
class brentFunctor_FastLMM_woodbury: public BrentFunctor
   {
   const real *SU_y;
   const real *SU_x;
   const real *X;
   const real *y;
   const real *covar;
   BLAS_INT    rankU;
   const real *S;
   const real *SUcovar;
   BLAS_INT    d;
   BLAS_INT    ncovar;
   BLAS_INT    n;

public:
   brentFunctor_FastLMM_woodbury( const BLAS_INT d, const BLAS_INT ncovar, const BLAS_INT n, const BLAS_INT rankU, const real *SU_y, const real *SU_x, const real *y, const real *X, const real *S, const real *SUcovar, const real *covar );
   ~brentFunctor_FastLMM_woodbury();
   virtual real operator()(const real logdelta);
   };

#endif   //FastLMM_woodbury
