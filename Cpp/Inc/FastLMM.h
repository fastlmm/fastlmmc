#if !defined( FastLMM_h )
#define FastLMM_h

class ResultsLMM
{
public:
   ResultsLMM();
   ResultsLMM(string modelID, vector<string> snpLabels, string phenLabel, int n_individuals, int ncovariates, int *index, double *LL_alt, double *LL_0, double *pvals, double *snpWeights_alt, double *bias_alt, double *sigg2_alt, double *sige2_alt, double *delta_alt, double *bias_0, double *sige2_0,  double *sigg2_0, double *delta_0, double *qvals);
   vector<string> snpLabels;  //identifier of the SNP
   string phenLabel;          //the name of the phenotype
   int ncovariates;           //the number of covariates in the study
   int n_individuals;         //the number of individuals in the study
   int *index;                //The results will be printed out, sorted by index
   double *LL_alt;            //Log-likelihoods of the alternative models
   double *LL_0;              //Log-likelihoods of the null models
   double *pvals;             //p-values
   double *snpWeights_alt;    //the SNP weights of the alternative models
   double *bias_alt;          //the bias weights of the alternative models
   double *sigg2_alt;         //the genetic variance of the alternative models
   double *sige2_alt;         //the random vriance of the alternative models
   double *delta_alt;         //the fraction of sige2_alt/sigg2_alt
   double *bias_0;            //the bias terms of the null models
   double *sigg2_0;           //the genetic variance of the null models
   double *sige2_0;           //the random variance of the null models
   double *delta_0;           //fraction of sige2_0/sigg2_0
   double *qvals;             //q-values (minimum FDR, under which the current SNP would be significant)
   void WriteResultFile( string &resultFile );
};

class FastLMM
{
public:
   static const bool brentExhaustive;
   string modelID;
   string phenLabel;
   vector<string> snpLabels;
   EigenSym *eigen;
   virtual void trainNull()=0;  //Train the Null Models
   virtual void trainAlt()=0;   //Train all alternative models
   virtual ResultsLMM getResults()=0;
};

#endif   //FastLMM_h
