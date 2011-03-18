#if !defined( FastLMM_fullrank_h )
#define FastLMM_fullrank_h

/*
 * this class does standard FastLMM.
 */
class FastLMM_fullrank: public FastLMM{
   double *Uy;
   double *UX;
   double *Ucovar;
   int* rank_beta;
   int rank_beta_0;
   int ncovar,n, nSNPs;
   double optimize_delta_Brentintervals(const int d, const int ncovar, const int n, const double *Uy, const double *currentUX,const double *S, const double *Ucovar,const int numberIntervals, const double minLogVal,const double maxLogVal,const double brentTol);
   double optimize_delta_BrentintervalsCheckTriplets(const int d, const int ncovar, const int n, const double *Uy, const double *currentUX,const double *S, const double *Ucovar,const int numberIntervals, const double minLogVal,const double maxLogVal,const double brentTol);
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
   FastLMM_fullrank(EigenSym *eigen);
   ~FastLMM_fullrank(void);
   //void set_Kernel(double *K);               //K is an [n by n] symmetric (lower triangle has to be stored) positve semi-definite matrix. This function performs the SVD  of K. K will be deleted.
   //void set_Kernel(double *U_, double *S_);  //U [N by n] and S [n by 1] are the Eigenvectors and Eigenvalues of the symmetric Kernel matrix K.
   void set_SNPs(double *SNPs, int nSNPs_, vector<string> snpLabels);  //SNPs [n by nSNPs] holds the nSNPs SNPs to be tested for each of n individuals
   void set_Phenotype(double *Phenotype, string phenLabel);    //Phenotype is a vector of length n, holding the (single) phenotype.
   void set_covariates(double *covar,int ncovariates);   //covariates is an [n*ncovar] matrix holding the covariates, including the bias
   void trainNull();                         //find ML parameters for the null model
   void trainAlt();                          //find ML parameters for all SNPs
   static double nLLeval(const double delta,const int d, const int ncovar,int n,const double *Uy, const double *UX, const double *S, const double *Ucovar, double *sigg2, double *sige2, double *beta, int *rank_beta, const bool returnMLparams);
   ResultsLMM getResults();
};

//Functor hp;dong all the needed members of FastLMM_fullrank to be handed to the Brent's routine
class brentFunctor_FastLMM_fullrank: public BrentFunctor{
   const double *Uy;
   const double *UX;
   const double *S;
   const double *Ucovar;
   int d,ncovar,n;
public:
   brentFunctor_FastLMM_fullrank(const int d_, const int ncovar_, const int n_, const double *Uy_,const  double *UX_, const double *S_, const double *Ucovar_);
   ~brentFunctor_FastLMM_fullrank();
   virtual double operator()(const double logdelta);
};

#endif   //FastLMM_fullrank_h
