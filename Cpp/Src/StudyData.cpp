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
 * Include Files
 */
#include "FastLmmC.h"

StudyData::StudyData()
   { 
   cIndividuals = 0;
   cPhenotypes = 0;
   cSnps = 0;
   cSnpsCovariance = 0;
   cCovariates = 0;
   ind_phenotype = 0;
   phenotypeData = nullptr;
   snpData = nullptr;
   snpDataCovariance = nullptr;
   kernelData = nullptr;
   covariatesData = nullptr;
   permutation = nullptr;
   groupSim = nullptr;
   nperm = 0;
   }
  
StudyData::~StudyData()
   {
   if ( snpData != nullptr )
      {
      // if snpData and snpCovarianceData are pointing to the same
      //   memory, only free it once!
      if ( snpDataCovariance == snpData ) { snpDataCovariance = nullptr; }
      delete[] snpData;
      snpData = nullptr;
      }
   if ( snpDataCovariance != nullptr ) { delete[] snpDataCovariance; snpDataCovariance = nullptr; }
   if ( phenotypeData != nullptr )     { delete[] phenotypeData;     phenotypeData = nullptr; }
   if ( kernelData != nullptr )        { delete[] kernelData;        kernelData = nullptr; }
   if ( covariatesData != nullptr )    { delete[] covariatesData;    covariatesData = nullptr; }
   if ( permutation != nullptr )       { delete[] permutation;       permutation = nullptr; }
   if ( groupSim != nullptr )          { delete[] groupSim;          groupSim = nullptr; }
   }

void StudyData::GeneratePermuations( size_t nperm_ )
   {
   if ( permutation && (nperm != nperm_) )
      {
      Fatal( "Permuations have already been generated." );
      }

   nperm = nperm_;
   permutation = new size_t[ nperm * cIndividuals ];
   size_t *permutationMax = permutation + ( nperm * cIndividuals );
   for (size_t *currentperm = permutation; currentperm < permutationMax; currentperm += cIndividuals )
      {
      MathFunctions::randperm( cIndividuals, currentperm );
      }
   }

void StudyData::Dump( FILE *pf, const string& comment ) const
   {
#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
   const char *szFmt1 = "\n  size_t cIndividuals: %Iu";
   const char *szFmt2 = "\n  size_t cPhenotypes: %Iu";
   const char *szFmt3 = "\n  size_t cSnps: %Iu";
   const char *szFmt4 = "\n  size_t cCovariates: %Iu";
   const char *szFmt5 = "\n  size_t ind_phenotype: %Iu";
   const char *szFmt6 = "\n  size_t nperm: %Iu";
#else                      // Linux/g++ uses a %zu specifier for size_t
   const char *szFmt1 = "\n  size_t cIndividuals: %zu";
   const char *szFmt2 = "\n  size_t cPhenotypes: %zu";
   const char *szFmt3 = "\n  size_t cSnps: %zu";
   const char *szFmt4 = "\n  size_t cCovariates: %zu";
   const char *szFmt5 = "\n  size_t ind_phenotype: %zu";
   const char *szFmt6 = "\n  size_t nperm: %zu";
#endif

   fprintf( pf, "\nStudyData: [%s]", comment.c_str() );
   fprintf( pf, szFmt1, cIndividuals );
   fprintf( pf, szFmt2, cPhenotypes );
   fprintf( pf, szFmt3, cSnps );
   fprintf( pf, szFmt4, cCovariates );
   fprintf( pf, szFmt5, ind_phenotype );
   fprintf( pf, szFmt6, nperm );
   DumpStringVector( pf, 2, individualLabels, "individualLabels" );
   DumpStringVector( pf, 2, phenotypeLabels, "phenotypeLabels" );
   DumpStringVector( pf, 2, covariateLabels, "covariateLabels" );
   DumpSnpInfoVector( pf, 2, snpInfo, "snpInfo" );
   DumpSnpInfoVector( pf, 2, snpInfoCovariance, "snpInfoCovariance" );
   DumpRealArray( pf, 2, phenotypeData, "phenotypeData", cIndividuals, cPhenotypes );
   DumpRealArray( pf, 2, snpData, "snpData", cIndividuals, cSnps );
   DumpRealArray( pf, 2, covariatesData, "covariatesData", cIndividuals, cCovariates );
   DumpRealArray( pf, 2, kernelData, "kernelData", cIndividuals, cIndividuals );
   DumpRealArray( pf, 2, snpDataCovariance, "snpDataCovariance", cIndividuals, cSnps );
   DumpGroupSim( pf, 2, groupSim, "groupSim" );
   fprintf( pf, "\n" );
   }
