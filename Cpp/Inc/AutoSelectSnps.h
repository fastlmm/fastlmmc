#if !defined( AutoSelectSnps_h )
#define CSnpInfo_h
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
 *  Declare types, structures, and publics
 */
struct SnpAndPvalue           // struct to hold Snp & p-value association
   {
   string idSnp;
   double pval;               
   };

/*
 * Declare the AutoSelectSnps related functions
 */
void DoAutoSelectSnps();
void InitCovariates( StudyData *studydata, CCovariatesData *covariatesData );
void InitializeStudyDataFromPlinkFile( StudyData *studyData_, const CPlinkFile &plinkFileData, vector< size_t > &rgIdxIndividuals );
void ResetSnpData( StudyData *studyData, size_t *sorted, size_t nSelectedSNPs );
void WriteAutoSelectCrossValidationFile( const string& autoSelectOutputFilenameRoot, size_t cSnpsSelected, const vector<double>& averageMeanSquaredError, const vector<double>& averageLogLikelihood );
void WriteAutoSelectSnpsFile( const string& autoSelectOutputFilenameRoot, const vector<SnpInfo>& rgSnpInfo, size_t cSnpsToPrint, const vector<size_t>& rgSortedIdx );
void WriteSelectedSnpsAsPedFile( const string& pedFileBase, const CPlinkFile& plinkFileT );

SnpFilterOptions *SelectTopK();
vector< size_t > RunLinReg( StudyData& data );
vector< size_t > RunLinRegWithThinning( StudyData& data );
vector< SnpAndPvalue >* RunLinRegSnpAndPvalueInfo( StudyData& data, size_t limit );

#endif   // AutoSelectSnps_h
