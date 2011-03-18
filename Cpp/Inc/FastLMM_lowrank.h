#if !defined(FastLMM_lowrank_h)
#define FastLMM_lowrank_h

/*
 * this class does low rank FastLMM.
 */
class FastLMM_lowrank: public FastLMM{
   double *Uy;
   double *y_UUy;
   double *UX;
   double *X_UUX;
   double *Ucovar;
   double *covar_UUcovar;
   int rankU;
   int* rank_beta;
   int rank_beta_0;
   int ncovar,n, nSNPs;
   double optimize_delta_Brentintervals(const int d, const int ncovar, const int n,const int rankU, const double *Uy, const double *currentUX, const double *y_UUy, const double *currentX_UUX,const double *S, const double *Ucovar, const double *covar_UUcovar,const int numberIntervals, const double minLogVal,const double maxLogVal,const double brentTol);
   double optimize_delta_BrentintervalsCheckTriplets(const int d, const int ncovar, const int n,const int rankU, const double *Uy, const double *currentUX, const double *y_UUy, const double *currentX_UUX,const double *S, const double *Ucovar, const double *covar_UUcovar,const int numberIntervals, const double minLogVal,const double maxLogVal,const double brentTol);
public:
   double *sigg2;
   double *sige2;
   double *delta;
   double *beta;
   double *LL;
   double sigg2_0;
   double sige2_0;
   double delta_0;
   double *beta_0;
   double *pvals;
   double LL_0;
   FastLMM_lowrank(EigenSym *eigen);
   ~FastLMM_lowrank(void);
   //void set_Kernel(int rankK, double *K);//K is an [n by n] symmetric (lower triangle has to be stored) positve semi-definite matrix. This function performs the rankK-SVD  of K. K will be deleted.
   //void set_Kernel_CovFromSNPsample(int nSNPsample, double *SNPsample);//SNPsample is an [n by nSNPsample] matrix. This function performs an SVD of SNPsamples to estimate a min(n,nSNPsample)-rank Covariance Kernel.
   //void set_Kernel_CovFromSNPsample2(int nSNPsample, double *SNPsample);
   //void set_Kernel(const int rankK, double *U_, double *S_);//U [n by rankK] and S [rankK by 1] are the Eigenvectors and Eigenvalues of the symmetric rankK-rank Kernel matrix K.
   void set_SNPs(double *SNPs, int nSNPs_, vector<string> snpLabels);//SNPs [n by nSNPs] holds the nSNPs SNPs to be tested for each of n individuals
   void set_Phenotype(double *Phenotype, string phenLabel);//Phenotype is a vector of length n, holding the (single) phenotype.
   void set_covariates(double *covar,int ncovariates);//covariates is an [n*ncovar] matrix holding the covariates, including the bias
   void trainNull();     //find ML parameters for the null model
   void trainAlt();      //find ML parameters for all SNPs
   ResultsLMM getResults();
   static double nLLeval(const double delta,const int d, const int ncovar,int n, const int rankU,const double *Uy, const double *UX,const double *y_UUy, const double *X_UUX, const double *S, const double *Ucovar, const double *covar_UUcovar, double *sigg2, double *sige2, double *beta, int *rank_beta, const bool returnMLparams);
};

//Functor holdong all the needed members of FastLMM_lowrank to be handed to the Brent's routine
class brentFunctor_FastLMM_lowrank: public BrentFunctor{
   const double *Uy;
   const double *UX;
   const double *X_UUX;
   const double *y_UUy;
   const double *covar_UUcovar;
   int rankU;
   const double *S;
   const double *Ucovar;
   int d,ncovar,n;
public:
   brentFunctor_FastLMM_lowrank(const int d_, const int ncovar_, const int n_,const int rankU_, const double *Uy_,const  double *UX_, const double *y_UUy_,const  double *X_UUX_, const double *S_, const double *Ucovar_,const double *covar_UUcovar_);
   ~brentFunctor_FastLMM_lowrank();
   virtual double operator()(const double logdelta);
};

#endif   //FastLMM_lowrank
