#if !defined( LinReg_h )
#define LinReg_h
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

struct SnpAndPvalue;

class LinReg
   {
public:
   LinReg(StudyData *pStudyData);
   ~LinReg();

   real     sig2_0;
   real     LL_0;          // likelihood of null model
   BLAS_INT rank_beta_0;

   string   modelID;
   string   phenLabel;     // phenotype name

   StudyData *data;        // data read to operate on

   real     *pvals;        // vector of p-values - 1 for each SNP (data->cSnps)
   real     *sig2;         // 
   real     *beta;         // matrix of Covariate Weights - cSnps x cCovariates+1
   real     *beta_0;       // vector of Null Covariate Weights - 1 for each covariate (data->cCovariates)
   real     *LL;           // vector of likelihood of alternative model - 1 for each SNP (data->cSnps)
   real     *f_stat;       // matrix of f-statistics of the SNPs ? - cSNPs x cCovariates+1
   real     *f_stat_0;
   BLAS_INT *rank_beta;

   void     trainNull();  //Train the Null Models
   void     trainAlt();   //Train all alternative models
   real     evalLL( BLAS_INT n, BLAS_INT d_, real* y, real* beta_, real*X, real*sig2_, real *f_stat_, BLAS_INT *rank_beta_ );
   real     evalLL( size_t iSNP, bool altModel );
   vector< size_t > GetOrderedSNPIndex();
   vector< SnpAndPvalue > *GetOrderedSnpAndPvalueInfo( size_t limit );

   void     WriteResultFile( const string& filename );

   //TODO: Remove after new output format is verified
   void     WriteResultFile_V1( const string& filename );
   void     WriteResultFile_V2( const string& filename, bool fVerboseOutput_, real maxPvaluePrintThreshold );

   void     trainPairs();
   size_t   TrainNullAndAltForEpistaticWithSnpPairs( FILE *pFile, size_t iSnp1Base, size_t iSnp1, size_t iSnp2Base, size_t iSnp2, real* rgCovariatesTemp );
   };

#endif   //LinReg_h
