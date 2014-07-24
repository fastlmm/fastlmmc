#if !defined( CLogisticRegression_h )
#define CLogisticRegression_h
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


// TODO: REFACTOR
//       CLogisticRegression is a 'clone' of LinReg until 
//        we can do a better refactoring job.
class CLogisticRegression
   {
public:
   CLogisticRegression( StudyData *pStudyData );
   ~CLogisticRegression();

   real EvalLL( size_t iSNP, bool altModel );   // evaluate the Log Likelihood for each SNP
   void TrainNull();           // Train the Null Model
   void TrainAlt();            // Train all alternative models
   void WriteResultFile( const string& filename, bool fVerboseOutput_, real maxPvaluePrintThreshold );

private:
   // Routine provided to optimizer to maximize LL
   double ProbRoutine( vector< double >::iterator itWgts, vector< double >::iterator itGradient );

   // Storage needed
   real     LL_0;          // likelihood of null model
   BLAS_INT rank_beta_0;   // 

   string   modelID;
   string   phenLabel;     // phenotype name

   StudyData *data;        // data read to operate on

   real     *pvals;        // vector of p-values.  There is 1 for each SNP (data->cSnps)
   real     *beta;         // matrix of Covariate Weights - cSnps x cCovariates+1
   real     *beta_0;       // vector of Null Covariate Weights - 1 for each covariate (data->cCovariates)
   real     *LL;           // vector of Log Likelihood of alternative model - 1 for each SNP (data->cSnps)
   BLAS_INT *rank_beta;    //

   // the additional members for logistic regression support
   size_t   cRows;         // count of rows in matrix (cSnps or cCovariates or ?)
   size_t   cColumns;      // count of columns in matrix (cIndividuals)
   size_t   cVars;         // cSnps + cCovariates
   real     sumLL;         // sum of vLL

   // gcc does not compile the unique_ptr< real[] >
   //unique_ptr< real[] > pMatrixData;      // unique_pointer to raw memory for data
   real     *pMatrixData;                 // raw pointer to memory for data

   vector< vector< double > > vMatrix;    // work area for optimizer 
   vector< real > vPhenotype;             // phenotype storage for optimizer to work on
   vector< real > vWeights;               // weight vector produced by optimizer
   vector< real > vZ;                     // MACC[]
   vector< real > vP;                     // vP[] = 1 / (1 + EXP(-vZ[]) )
   vector< real > vLL;                    // vector of LogLikelihoods
   vector< real > vErr;                   // value for optimizer to minimize
   };

#endif   //CLogisticRegression_h
