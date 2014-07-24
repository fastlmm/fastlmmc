#if !defined( FastLMM_h )
#define FastLMM_h
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

class ResultsLMM
   {
public:
   ResultsLMM();
   ~ResultsLMM();             // need destructor to cleanup the memory allocations in ResultsLMM

   vector<SnpInfo> snpInfo;   //all the SNP related info
   string phenLabel;          //the name of the phenotype
   size_t ncovariates;        //the number of covariates in the study
   size_t n_individuals;      //the number of individuals in the study
   size_t *index;             //The results will be printed out, sorted by index
   size_t *nExclude;
   size_t *excludeSnpStart;
   real *LL_alt;              // Log-likelihoods of the alternative models
   real *LL_0;                // Log-likelihoods of the null models
   real *pvals;               // p-values
   real *snpWeights_alt;      // the SNP weights of the alternative models
   real *bias_alt;            // the bias weights of the alternative models
   real *sigg2_alt;           // the genetic variance of the alternative models
   real *sige2_alt;           // the random vriance of the alternative models
   real *delta_alt;           // the fraction of sige2_alt/sigg2_alt
   real *bias_0;              // the bias terms of the null models
   real *sigg2_0;             // the genetic variance of the null models
   real *sige2_0;             // the random variance of the null models
   real *delta_0;             // fraction of sige2_0/sigg2_0
   real *qvals;               // q-values (minimum FDR, under which the current SNP would be significant)
   real *f_stat;              // f-statistics of the SNPs
   BLAS_INT *dof;             // degrees of freedom

   real *beta_0;              // array of Null Covariate Weights
   real *beta;                // array of Alternative Covariate Weights

   //  Add support for iterative results and printing them
   bool   fIterativeResults;  // Are we using an iterative process to compute GWAS?
   size_t iResultBlock;       // report the iteration this block is printing from.
   size_t cPriorResults;      // the number of SNP results printed in previous blocks

   void Clear();              // release all the allocated memory fields
   void PrepareForResults( size_t cSnps, size_t cCovariates ); // Allocate the memory for these results
   void WriteResultFile( const string& resultFile );

   // remove _V1 & _V2 after validation complete
//TODO: DELETE
//   void WriteResultFile_V1( const string& resultFile );
   void WriteResultFile_V2( const string& resultFile, bool fWriteVerboseOutput, real maxPvaluePrintThreshold );

   void OpenResultsFile( const string& filename_, bool fWriteVerboseOutput_, real pvaluePrintThreshold_, size_t cCovariates_ );
   void WriteResultsBody();         // 
   void WriteResultsSummary();      // write the summary information for 'skipped' SNPs (if any)
   void CloseResultsFile();         // 

private:
   string resultsFilename;                   // save name of results file
   FILE   *pFileResults;                     // and the file pointer

   // Add support for formatting results statistics
   bool   fProduceVerboseOutput;             // print the 'extra' columns for each SNP too
   const char *szSeparator;                  // cache pointer to the separator and format strings
   const char *fmtInt;
   const char *fmtString;
   const char *fmtReal;
   bool   fExcludeSnps;                      // Display the region of excluded SNPs 
   size_t firstCovariate;                    // if biasTerm is used, first covariate is bias so skip it
   size_t cCovariates;                       // count of Covariates in array
   size_t cSnpsFailedPvaluePrintThreshold;   // how many SNPs did not print due to thresholding
   real   pvaluePrintThresholdMax;           // do not print pvals greater than pvaluePrintThreshold
   };

class FastLMM
   {
public:
   FastLMM(EigenSym *eigen);
   ~FastLMM();
   static const bool brentExhaustive;
   string   modelID;
   string   phenLabel;
   vector<SnpInfo> snpInfo;            // Information about each SNP in test data
   vector<SnpInfo> snpInfoCovariance;  // Information about each SNP in W
   size_t   ncovar;                    // number of covariates (cCovariates)
   size_t   n;
   size_t   nSNPs;                     // number of SNPs in test data (cSnps)
   size_t   nSNPsW;                    // number of SNPs in covariate data

   size_t   *excludeSnpStart;
   size_t   *nExclude;                 // the indices of the SNPs to exclude
   EigenSym *eigen;
   real     *pvals;                    // p-values
   real     *sigg2;
   real     *sige2;
   real     *delta;
   real     *beta;
   real     *LL;
   real     *sigg2_0;                  // the genetic variance of the null models
   real     *sige2_0;                  // the random variance of the null models
   real     *delta_0;                  // fraction of sige2_0/sigg2_0
   real     *beta_0;                   // array of Null Covariate Weights
   real     *LL_0;                     // Log-likelihoods of the null models
   real     *f_stat;                   // f-statistics of the SNPs
   BLAS_INT *rank_beta;
   BLAS_INT *rank_beta_0;
   size_t   numberIntervals;           // Brent intervals
   real     minLogVal;                 // Brent min Log Val
   real     maxLogVal;                 // Brent max Log val
   real     brentTol;                  // Brent tolerance
   size_t   brentMaxIter;

   virtual void trainNull() = 0;       //Train the Null Models
   virtual void trainAlt() = 0;        //Train all alternative models

   void getResults( ResultsLMM& resultsLMM );
   void getResultsNullOnly( ResultsLMM& resultsLMM );
   void updateExcludeList();

   // solve the least squares problem
   static BLAS_INT solveLstSq( BLAS_INT n, BLAS_INT d_, real *beta_temp, BLAS_INT *rank_beta, real *XSdi, real *res , real *f_stat );
    
   // solve the least squares problem for precomputed XX and Xy 
   //   in the likelihood function but try to use specialized 
   //   functions for small d_.
   static BLAS_INT solveLstSqSym( BLAS_INT d_, real *XX, real *Xy, real *beta_temp, BLAS_INT *rank_beta, real *f_stat, real *logDetXX );
   };

#endif   //FastLMM_h
