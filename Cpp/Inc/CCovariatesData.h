#if !defined( CCovariatesData_h )
#define CCovariatesData_h
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

class CCovariatesData
   {
public:
   size_t         n_individuals;
   size_t         n_covariates;
   string         filename;

   vector<string> indIDs;
   vector<string> famIDs;
   vector<vector<real> > covariates;      // temporary storage of Covariates
   vector<bool>   markedAsConstant;       // indicator of constant covariates
   map< string, size_t > keyIndividualToCovariatesDataIndividualIndex;

   CCovariatesData();
   ~CCovariatesData();
   void clear();                          // clear all CCovariatesData collections
   void ReadCovariatesFile( string& covariatesFile );
   void WriteCovariatesFile( const string& covariatesFile );
//   void CreateColumnMajorArrayFromVectors();
   void CreateColumnMajorArrayFromVectorsAndMapping( const vector<string>& columnLabelOrder , real *covariatesArray, size_t n_covariates );
   void ImputeAndZeroMeanCovariates( const vector<string>& individualsSelected );
   void ImputeAndZeroMeanCovariatesFix( const vector<string>& individualsSelected );
   void MarkConstantCovariates();

   bool     FKeyIndividualInDatFile( const string& keyIndividual );
   size_t   IdxFromKeyIndividual( const string& keyIndividual );
   size_t   IdxFromIdFamilyAndIdIndividual( const string& idFamily, const string& idIndividual );

   real     CovariateValue( size_t iCovariate, const string& keyIndividual );
   void     CovariateValues( size_t iCovariate, const vector< string >& rgKeyIndividuals, real *rgOutArray );
   };

#endif   // CCovariatesData_h
