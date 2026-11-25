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
 * FastLmmC - {Factored Spectrally Transformed Linear Mixed Models in C++}
 *
 *         File Name:   AutoSelectSnps.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements FastLMM's -autoSelect option
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Include Files
 */
#include "FastLmmC.h"

/*
 * My local classes
 */
class CDnaRegion
   {
public:
   CDnaRegion();
   CDnaRegion( const CDnaRegion &rhs );
   CDnaRegion( size_t iSnp_, const SnpInfo* snpInfo_, real offset_=0.0 );
   bool operator < ( const CDnaRegion& rhs ) const;
   bool const FInRegion( const SnpInfo& snpInfoT ) const;
   bool const FInRegion( int iChromosome_, real geneticDistance_ ) const;
   void UpdateRegionExtents( real offset_ );

//private:
   size_t   iSnp;
   const SnpInfo* snpInfo;
   int      iChromosome;   // 1-22, X, Y, XY, MT, or 0 if unplaced
   real     regionCenter;
   real     regionMin;
   real     regionMax;
   };

CDnaRegion::CDnaRegion()
   {
   iSnp        = 0;
   snpInfo     = nullptr;
   iChromosome = 0;
   regionCenter= 0.0;
   regionMin   = 0.0;
   regionMax   = 0.0;
   };

CDnaRegion::CDnaRegion( const CDnaRegion &rhs )
   {
   iSnp        = rhs.iSnp; 
   snpInfo     = rhs.snpInfo; 
   iChromosome = rhs.iChromosome; 
   regionCenter= rhs.regionCenter;
   regionMin   = rhs.regionMin; 
   regionMax   = rhs.regionMax;
   }

CDnaRegion::CDnaRegion( size_t iSnp_, const SnpInfo* snpInfo_, real offset_ )
   {
   iSnp = iSnp_;
   snpInfo = snpInfo_;
   iChromosome = snpInfo->iChromosome;
   regionCenter = snpInfo->geneticDistance;
   UpdateRegionExtents( offset_ );
   }

bool const CDnaRegion::FInRegion( const SnpInfo& snpInfoT ) const
   {
   return( FInRegion( snpInfoT.iChromosome, snpInfoT.geneticDistance ) );
   }

bool const CDnaRegion::FInRegion( int iChromosome_, real geneticDistance_ ) const
   {
   if ( (iChromosome_ != iChromosome)
     || (geneticDistance_ < regionMin)
     || (geneticDistance_ > regionMax)
      )
      {
      return( false );
      }
   return( true );
   }

inline bool CDnaRegion::operator < ( const CDnaRegion &rhs ) const
   {
   if ( iChromosome < rhs.iChromosome )
      {
      return( true );
      }
   if ( iChromosome > rhs.iChromosome )
      {
      return( false );
      }
   return( regionMax < rhs.regionMax );
   }

void CDnaRegion::UpdateRegionExtents( real offset_ )
   {
   regionMin = std::max( 0.0, regionCenter - offset_ );
   regionMax = regionCenter + offset_;
   }

   /*
 * assume that covariate is already in zero mean state.
 */
void InitCovariates( StudyData *studydata, CCovariatesData *covariatesData )
   {
   if ( covariatesData->n_covariates > 0 )
      {
      //covariates have been set
      if ( useBiasTerm )
         {
         studydata->covariatesData = new real[ (covariatesData->n_covariates + 1) * studydata->cIndividuals ];
         studydata->cCovariates = covariatesData->n_covariates + 1;
         covariatesData->CreateColumnMajorArrayFromVectorsAndMapping( studydata->individualLabels, studydata->covariatesData, covariatesData->n_covariates );
         //zeroMeanColumns( studydata->covariatesData, studydata->cIndividuals, studydata->cCovariates-1, true );

         //append an array of ones -> not needed, if all variables are zero-mean
         for ( size_t i = (studydata->cCovariates - 1) * studydata->cIndividuals; i < studydata->cCovariates * studydata->cIndividuals; ++i )
            {
            studydata->covariatesData[ i ] = 1.0;
            }
         }
      else
         {
         //no bias term (no bias term is added to the covariates)
         studydata->covariatesData = new real[ (covariatesData->n_covariates) * studydata->cIndividuals ];
         studydata->cCovariates = covariatesData->n_covariates;
         covariatesData->CreateColumnMajorArrayFromVectorsAndMapping( studydata->individualLabels, studydata->covariatesData, covariatesData->n_covariates );
         //zeroMeanColumns( studydata->covariatesData, studydata->cIndividuals, studydata->cCovariates, true );     //Imputation for the kernel SNPs
         }
      }
   else
      {
      if ( useBiasTerm )
         {
         //generate the covariates as an array of ones and make them the only covariates
         studydata->covariatesData = new real[ studydata->cIndividuals ];
         for ( size_t i=0; i < studydata->cIndividuals; ++i )
            {
            studydata->covariatesData[i] = 1.0;
            }
         studydata->cCovariates = 1;
         }
      else
         {
         //no bias term, and no covariates
         studydata->covariatesData = new real[ 0 ];
         studydata->cCovariates = 0;
         }
      }
   }

/*
 * initialize the StudyData from CPlinkFile data and 
 *   the specific individuals to include
 */
void InitializeStudyDataFromPlinkFile( StudyData *studyData_, const CPlinkFile &plinkFileData, vector< size_t > &rgIdxIndividuals )
   {
   if ( studyData_ == nullptr )
      {
      Fatal( "Passed nullptr to Initialize StudyDataFromPlinkFile().  This is not good." );
      }

   StudyData* psd = studyData_;     // psd = Pointer to Study Data
   psd->cIndividuals = rgIdxIndividuals.size();
   psd->cPhenotypes = plinkFileData.cPhenotypes;
   psd->cSnps = plinkFileData.cSnps;
   if ( plinkFileData.covariatesFileData != nullptr )
      {
      psd->cCovariates = plinkFileData.covariatesFileData->n_covariates;
      }
   
   // cSnpsCovariance does not come from Plink file

   psd->ind_phenotype = plinkFileData.ind_phenotype;
   for ( vector< size_t >::iterator iter = rgIdxIndividuals.begin(); iter!=rgIdxIndividuals.end(); ++iter )
      {
      psd->individualLabels.push_back( plinkFileData.individualLabels[ *iter ] );
      }

   copy( plinkFileData.phenotypeLabels.begin(), plinkFileData.phenotypeLabels.end(), back_inserter( psd->phenotypeLabels ) );
   copy( plinkFileData.rgSnpInfo.begin(), plinkFileData.rgSnpInfo.end(), back_inserter( psd->snpInfo ) );

   // snpInfoCovariance does not come from Plink file
   // covariateLabels does not come from Plink file - not there

   // Copy phenotype data (cIndividual rows x cPhenotype columns)
   psd->phenotypeData = new real[ psd->cIndividuals * psd->cPhenotypes ];
   for( size_t iPhenotype=0; iPhenotype<psd->cPhenotypes; ++iPhenotype )
      {
      for( size_t iIndividual=0; iIndividual<psd->cIndividuals; ++iIndividual )
         {
         size_t idxDest = ColumnMajorIndex( psd->cIndividuals, psd->cPhenotypes, iIndividual, iPhenotype );
         size_t idxSrc = ColumnMajorIndex( plinkFileData.cIndividuals, plinkFileData.cPhenotypes, rgIdxIndividuals[ iIndividual ], iPhenotype );
         psd->phenotypeData[ idxDest ] = plinkFileData.phenotypeData[ idxSrc ];
         }         
      }

   // copy snpData
   psd->snpData = new real[ psd->cIndividuals * psd->cSnps ];
   for( size_t iIndividual=0; iIndividual<psd->cIndividuals; ++iIndividual )
      {
      for( size_t iSnp=0; iSnp<psd->cSnps; ++iSnp )
         {
         size_t idxDest = ColumnMajorIndex( psd->cIndividuals, psd->cSnps, iIndividual, iSnp );
         size_t idxSrc = ColumnMajorIndex( plinkFileData.cIndividuals, plinkFileData.cSnps, rgIdxIndividuals[ iIndividual ], iSnp );
         psd->snpData[ idxDest ] = plinkFileData.snpData[ idxSrc ];
         }         
      }
   }

/*
 * SelectTopK() is an option that is used when the program is running in
 *   a constrained memory system.  This options reduces the total number
 *   of SNPs to consider by selecting only the 'k' best most influential
 *   SNPs to propagate through the system.
 *   SelectTopK() reads the covariance(Sim) data and then runs a linear
 *   regression test to identify the 'k SNPs with the lowest p-values
 */
SnpFilterOptions *SelectTopK()
   {
   CCovariatesData covariatesData;
   if ( FUseCovariatesFile() )
      {
      covariatesData.ReadCovariatesFile( covariatesFile );
      }

   size_t cbMemoryTotal = GetTotalSystemMemory();
   size_t cbMemoryToUse = (size_t)((double)cbMemoryTotal * mFraction);

   string filetype;
   switch( plinkFileTypeCovariance )
      {
   case Natural:
      filetype = ".ped";
      break;
   case Binary:
      filetype = ".bed";
      break;
   case Transposed:
      filetype = ".bed";
      break;
   case Dosage1:
   case Dosage2:
      filetype = ".dat";
      break;
   case ImputedGenotype:
      filetype = ".big";
      break;
   case UnknownFileType:
   default:
      Fatal( "Unknown Plink covariance filetype. Found [%d]", plinkFileTypeCovariance );
      }

   size_t cbFile = 0;
   string filename = plinkFilenameRootCovarianceVector[0] + filetype;
   FILE *pFile = fopen( filename.c_str(), "rb" );
   if ( pFile == nullptr )
      {
      // If no 'filename', check for a compressed file version too
      filename += ".gz";
      pFile = fopen( filename.c_str(), "rb" );
      }
   if ( pFile != nullptr )
      {
      fseek( pFile, 0, SEEK_END );
      cbFile = ftell( pFile );   // was a shadowing local, which left njobs always 1
      fclose( pFile );
      }

   int njobs = static_cast<int>(cbFile / cbMemoryToUse) + 1;      // force at least 1 iteration

   vector< vector< SnpAndPvalue > *> *snpinfos = new vector< vector< SnpAndPvalue > *>();
   ProgressNL( "  ++Selecting top SNPs Using LinReg: " );
   for ( int ijob = 0; ijob < njobs; ijob++ ) 
      { 
      CPlinkFile *plinkFile = new CPlinkFile();
      SnpFilterOptions snpFilterOpts;
      snpFilterOpts.SetByJob(njobs, ijob);
      ProgressNL( "  ++  Start Loading Covariance Data:" );
      plinkFile->ReadPlinkFiles( plinkFilenameRootCovarianceVector[0], plinkFileTypeCovariance, alternatePhenotypeParameters, snpFilterOpts );
      ProgressNL( "  --    End Loading Covariance Data:" );

      // Eigenstrat Imputation for the kernel snps to be
      zeroMeanColumns( plinkFile->phenotypeData, plinkFile->cIndividuals, plinkFile->cPhenotypes, false );

      // perform preprocessing on the (soon to be) covariance SNPs
      if ( simMissingness )
         {
         replaceSNPsByMissingness( plinkFile->snpData, plinkFile->cIndividuals, plinkFile->cSnps );
         }
      else
         {
         // Zero mean each individual and divide by sqrt(N_SNPs)
         // for autoselect we always want to zero mean and unit variance the SNPs
         ImputeAndZeroMeanSNPs( plinkFile->rgSnpInfo, plinkFile->snpData, plinkFile->cIndividuals, plinkFile->cSnps, unitVarianceSimilaritySNPs, zeroMeanSimilaritySNPs, normalizationInfo, usePseudocounts );
         preprocessSNPsForBuildingCovariance( plinkFile->snpData, plinkFile->cIndividuals, plinkFile->cSnps, plinkFile->cSnps );
         }

      // covariate imputation in the plinkFile
      covariatesData.ImputeAndZeroMeanCovariatesFix( plinkFile->individualLabels );

      StudyData studydata;
      studydata.cIndividuals = plinkFile->cIndividuals;
      studydata.cPhenotypes = plinkFile->cPhenotypes;
      studydata.cSnps = plinkFile->cSnps;
      studydata.individualLabels = plinkFile->individualLabels;
      studydata.phenotypeLabels = plinkFile->phenotypeLabels;
      studydata.ind_phenotype = plinkFile->ind_phenotype;
      studydata.snpInfo = plinkFile->rgSnpInfo;
      studydata.phenotypeData = plinkFile->phenotypeData;
      studydata.snpData = plinkFile->snpData;

      InitCovariates( &studydata, &covariatesData );

      // validate topk value is legal on the high side.
      if ( topk > studydata.cSnps )
         {
         Warn( "TopKbyLinReg value larger than number of SNPs.  Changing TopKbyLinReg to: %d", studydata.cSnps );
         topk = studydata.cSnps;
         }

      vector< SnpAndPvalue > *sortedSnpAndPvalues = RunLinRegSnpAndPvalueInfo( studydata, topk );
      snpinfos->push_back( sortedSnpAndPvalues );
      delete plinkFile;
      }

   SnpFilterOptions *snpFilterOpts = new SnpFilterOptions();
   for ( size_t i = 0; i < topk; ++i ) 
      {
      double minp = 1;
      int min_index = 0;
      for ( int j = 0; j < njobs; ++j ) 
         {
         vector< SnpAndPvalue > *tmp = snpinfos->at( j );
         if ( minp > tmp->at(0).pval ) 
            {
            minp = tmp->at(0).pval;
            min_index = j;
            }
         }
      vector< SnpAndPvalue > *tmp = snpinfos->at( min_index );
      string idSnp = tmp->front().idSnp;
      tmp->erase( tmp->begin() );
      snpFilterOpts->snpIdsToExtract.push_back( idSnp );
      }

   snpFilterOpts->filterType = snpFilterOpts->FilterByFilelist;
   snpFilterOpts->extractFile = "internal";
   snpFilterOpts->cSnpsInJob = static_cast<int>(topk);
   ProgressNL( "    ++  Selection Using LinReg Done: \n" );

   for ( int ijob = 0; ijob < njobs; ++ijob )
      {
      vector< SnpAndPvalue > *tmp = snpinfos->at( ijob );
      delete tmp;
      }

   delete snpinfos;
   return( snpFilterOpts );
   }

/*
 *
 */
void ResetSnpData( StudyData *studyData, size_t *sorted, size_t nSelectedSNPs )
   {
   real *newsnpdata = new real[ studyData->cIndividuals * nSelectedSNPs ];
   for ( size_t isnp = 0; isnp < nSelectedSNPs; ++isnp )
      {
      for ( size_t ind = 0; ind < studyData->cIndividuals; ++ind )
         {
         size_t idxDest = ColumnMajorIndex( studyData->cIndividuals, nSelectedSNPs, ind, isnp );
         size_t idxSrc = ColumnMajorIndex( studyData->cIndividuals, studyData->cSnps, ind, sorted[ isnp ] );
         newsnpdata[ idxDest ] = studyData->snpData[ idxSrc ];
         }
      }

   delete[] studyData->snpData;
   studyData->snpData = newsnpdata;
   studyData->cSnps = nSelectedSNPs;
   }

/*
 * Compute GWAS w/ Linear Regression returning a list of SNP indexes
 */
vector< size_t > RunLinReg( StudyData& data )
   {
   //ProgressNL( "  Compute GWAS w/ Linear Regression:" );
   CTimer linRegTimer( true );
   LinReg lin = LinReg( &data );
   lin.trainNull();
   lin.trainAlt();
   linRegTimer.Stop();
   //linRegTimer.Report( -1, "                LinReg elapsed time: %s" );
   vector< size_t > sorted = lin.GetOrderedSNPIndex();
   return( sorted );
   }

/*
 * Run Linear Regression to get a list of SNPs and p-values ordered by p-values
 */
vector< SnpAndPvalue >* RunLinRegSnpAndPvalueInfo( StudyData& data, size_t limit )
   {
   LinReg lin = LinReg( &data );
   lin.trainNull();
   lin.trainAlt();
   vector< SnpAndPvalue > *ordered = lin.GetOrderedSnpAndPvalueInfo( limit );
   return( ordered );      // Caller must delete this pointer
   }

/*
 * Run Linear Regression with Thinning 
 */
vector< size_t > RunLinRegWithThinning( StudyData& data )
   {
   vector< size_t > rgIndexes = RunLinReg( data );
   vector< size_t > rgSelectedSnps;

   set< CDnaRegion > rgDnaRegions;

   for ( size_t cSnp = 0; cSnp<data.cSnps; ++cSnp )
      {
      size_t iSnp = rgIndexes[ cSnp ];
      SnpInfo* snpInfo = &(data.snpInfo[ iSnp ]);

      /*
       *  Create the region object w/ no extents to search
       *    for the actual GenomicDistance.  After the search
       *    update region extents prior to adding it to the
       *    regions collection
       */
      CDnaRegion dnaRegion( iSnp, snpInfo );
      auto iter = rgDnaRegions.upper_bound( dnaRegion );
      dnaRegion.UpdateRegionExtents( autoSelectRegionDistance );

      if ( (iter != rgDnaRegions.end())
        && (iter->FInRegion( snpInfo->iChromosome, snpInfo->geneticDistance ))
         )
         {
         // Already in a region.  Do Not Add to regions
                  //                                          genetic  region region
                  //                                          distance  min    max
         Verbose( verboseMinor, 
                    "          Snp Region overlap found--"
                  "\n                       Region SnpId: [%-10s]  %3d %8.4f  %8.4f  %8.4f"
                  "\n                      Thinned SnpId: [%-10s]  %3d %8.4f"
                  , iter->snpInfo->idSnp.c_str()
                  , iter->snpInfo->iChromosome
                  , iter->snpInfo->geneticDistance
                  , iter->regionMin
                  , iter->regionMax
                  , snpInfo->idSnp.c_str()
                  , snpInfo->iChromosome
                  , snpInfo->geneticDistance
                  );
         }
      else 
         {
         // this SNP's geneticDistance is not found in any region already
         //   created.  Add it to the regions set
         Verbose( verboseMinor, "          >>> Create region for SNP: [%-10s]  %3d %8.4f", snpInfo->idSnp.c_str(), snpInfo->iChromosome, snpInfo->geneticDistance );
         rgDnaRegions.insert( dnaRegion );
         rgSelectedSnps.push_back( iSnp );
         }
      }
   return( rgSelectedSnps );
   }

/*
 * WriteAutoSelectCrossValidationFile() outputs the results and statistics from the autoselect process
 */
void WriteAutoSelectCrossValidationFile( const string& autoSelectOutputFilenameRoot, size_t cSnpsSelected, const vector<double>& averageMeanSquaredError, const vector<double>& averageLogLikelihood )
   {
   /*
    * Write the Cross Validation File output
    */
   void SelectFormatStrings_V2( const char*& fmtInt_, const char*& fmtString_, const char*& fmtReal_ );

   const char *fmtInt    = nullptr;
   const char *fmtString = nullptr;
   const char *fmtReal   = nullptr;

   SelectFormatStrings_V2( fmtInt, fmtString, fmtReal );
   const char *szSeparator = "\t";

   string autoSelectOutputXvalFile = autoSelectOutputFilenameRoot + ".xval.txt";
   ProgressNL( "      Write Cross Validation output: [%s]", autoSelectOutputXvalFile.c_str() );

   FILE *pfileCrossValidation = fopen( autoSelectOutputXvalFile.c_str(), "w" );
   if ( pfileCrossValidation == nullptr )
      {
      Fatal( "Unable to open Cross Validation output: [%s]"
           "\n        CRT Error: %d [%s]", autoSelectOutputXvalFile.c_str(), errno, strerror( errno ) );
      }

   fprintf( pfileCrossValidation, "Bucket%s#SNPs%sMean_Squared_Error%sLog_Likelihood\n", szSeparator, szSeparator, szSeparator );
   for ( size_t i = 0; i < averageLogLikelihood.size(); ++i )
      {
      // fprintf( pfileCrossValidation, "%d%s%lf%s%lf\n", autoSelectSearchValues[ i ], szSeparator, (double)avg[ i ], szSeparator, (double)avg2[ i ] );
      fprintf( pfileCrossValidation, fmtInt, "", autoSelectSearchValues[ i ] );     // no leading separator
      int cSnps = std::min( autoSelectSearchValues[ i ], static_cast<int>(cSnpsSelected) );
      fprintf( pfileCrossValidation, fmtInt, szSeparator, cSnps );
      fprintf( pfileCrossValidation, fmtReal, szSeparator, averageMeanSquaredError[ i ] );
      fprintf( pfileCrossValidation, fmtReal, szSeparator, averageLogLikelihood[ i ] );
      fprintf( pfileCrossValidation, "\n" );
      }
   fclose( pfileCrossValidation );
   }

/*
 * WriteAutoSelectSnpsFile() outputs the specfic SNPs selected from the autoselect process
 */
void WriteAutoSelectSnpsFile( const string& autoSelectOutputFilenameRoot, const vector<SnpInfo>& rgSnpInfo, size_t cSnpsToPrint, const vector<size_t>& rgSortedIdx )
   {
   /*
    * Write the AutoSelect SNPs File output
    */
   string autoSelectOutputSnpsFile = autoSelectOutputFilenameRoot + ".snps.txt";
   ProgressNL( "       Write AutoSelect SNPs output: [%s]", autoSelectOutputSnpsFile.c_str() );

   FILE *pfileSnps = fopen( autoSelectOutputSnpsFile.c_str(), "w" );
   if ( pfileSnps == nullptr )
      {
      Fatal( "Unable to open AutoSelect SNPs output: [%s]"
           "\n        CRT Error: %d [%s]", autoSelectOutputSnpsFile.c_str(), errno, strerror( errno ) );
      }

   char *pchDebugEmitEntireSnpOrder = getenv( "FastLmmDebugEmitEntireSnpOrder" );
   if ( pchDebugEmitEntireSnpOrder != nullptr )
      {
      // This is a local copy of the number of SNPs to print.
      //   Update it to output the full array of SNPs.
      cSnpsToPrint = rgSortedIdx.size();     // print all SNPs in the vector
      }
   if ( cSnpsToPrint > 0 )
      {
      for ( size_t i = 0; i < cSnpsToPrint; ++i )
         {
         fprintf( pfileSnps, "%s\n", rgSnpInfo[ rgSortedIdx[ i ] ].idSnp.c_str());
         }
      }
   else
      {
      fprintf( pfileSnps, "# 0 SNPs selected\n" );
      }
   
   fclose( pfileSnps );
   }

/*
 *  Write a Plink PED/MAP file pair using the information from
 *    plinkFileT into bedFileBase as the file basename
 */
void WriteSelectedSnpsAsPedFile( const string& pedFileBase, const CPlinkFile& plinkFileT )
   {
   // Write the .ped file which hold the IID information:
   //   FamilyId IndividualId FatherId MotherId Sex Phenotype [genotype array]
   string fullped = pedFileBase + ".ped";
   FILE *pFile = fopen( fullped.c_str(), "wt" );
   for( size_t iIndividual = 0; iIndividual < plinkFileT.cIndividuals; ++iIndividual )
      {
      // print the id + hard code 0 0 1 -9 for IndividualInfo
      fprintf( pFile, "%s\t0\t0\t1\t-9\t", plinkFileT.individualLabels[ iIndividual ].c_str() );
      for( size_t iSnp = 0; iSnp < plinkFileT.cSnps; ++iSnp )
         {
         size_t idx = ColumnMajorIndex( plinkFileT.cIndividuals, plinkFileT.cSnps, iIndividual, iSnp );
         real rT = plinkFileT.snpData[ idx ];
         if ( rT == 2.0 )
            {
            fprintf( pFile, "\t2 2" );
            }
         else if ( rT == 1.0 )
            {
            fprintf( pFile, "\t2 1" );
            }
         else if ( rT == 0.0 )
            {
            fprintf( pFile, "\t1 1" );
            }
         else
            {
            fprintf( pFile, "\t0 0" );
            }
         }
      fprintf( pFile, "\n" );
      }
   fclose( pFile );

   // Write the .map file which holds the SnpInfo information
   //    Chr SnpId GeneticDistance GeneticPosition
   string fullmap = pedFileBase + ".map";
   pFile = fopen( fullmap.c_str(), "wt" );
   for ( size_t iSnpInfo = 0; iSnpInfo < plinkFileT.rgSnpInfo.size(); ++iSnpInfo )
      {
      const SnpInfo& snpInfo = plinkFileT.rgSnpInfo[ iSnpInfo ];
      fprintf( pFile, "%s\t%s\t%f\t%d\n",
                        snpInfo.idChromosome.c_str(),
                        snpInfo.idSnp.c_str(),
                        snpInfo.geneticDistance,
                        snpInfo.basepairPosition );
      }
   fclose( pFile );
   }

/*
 *  SetSnpCovariance() updates the in memory SNP covariance array
 *    with the results from AutoSelect
 */
void SetSnpCovariance( StudyData *data, size_t *idxSNPs, size_t nSNPsUsed )
   {
   StudyData *psd = data;
   delete[] psd->snpDataCovariance;    // delete previous covariance data

   psd->cSnpsCovariance = nSNPsUsed;
   psd->snpDataCovariance = new real[ psd->cIndividuals * nSNPsUsed ];

   for( size_t iIndividual = 0; iIndividual < psd->cIndividuals; ++iIndividual )
      {
      for( size_t iSnp=0; iSnp < nSNPsUsed; ++iSnp )
         {
         size_t idxDest = ColumnMajorIndex( psd->cIndividuals, psd->cSnps, iIndividual, iSnp );
         size_t idxSrc = ColumnMajorIndex( psd->cIndividuals, psd->cSnps, iIndividual, idxSNPs[ iSnp ]);
         psd->snpDataCovariance[ idxDest ] = psd->snpData[ idxSrc ];
         }
      }
   }

/*
 * Ensure the values in the autoSelectSeachrValues vector are consistent with
 *   the number of SNPs for processing (cSnps)
 */
void TrimAutoSelectSearchValues( vector< int > &autoSelectSearchValues_, size_t cSnps )
   {
   if ( cSnps < (size_t)autoSelectSearchValues_.back() )
      {
      while( autoSelectSearchValues_.size() && ( cSnps <= (size_t)autoSelectSearchValues_.back() ) )
         {
         autoSelectSearchValues_.pop_back();
         }
      autoSelectSearchValues_.push_back( static_cast< int >(cSnps) );
      }
   }

/*
 * DoAutoSelectSnps implements the -autoSelect option
 */
void DoAutoSelectSnps()
   {
   useWoodbury = false;
   //
   // Load the snp, individual genotypes, and covariate data from the plink file(s)
   //
   CCovariatesData covariatesData;
   if ( FUseCovariatesFile() )
      {
      covariatesData.ReadCovariatesFile( covariatesFile );
      }

   // Create CPlinkFile and initialize covariates information to SelectIndividuals if we have it
   CPlinkFile plinkFile;
   plinkFile.covariatesFileData = covariatesData.n_individuals ? &covariatesData : nullptr;

   ProgressNL( "  ++        Start Loading Test Data:" );
   if ( selectTopK )
      {
      //REVIEW:  What are the interactions between this SelectTopK 
      //         filter and the ones the user supplied in simFilterOpts?
      SnpFilterOptions * snpFilter = SelectTopK();
      plinkFile.ReadPlinkFiles( plinkFilenameRootCovarianceVector[ 0 ], 
                                plinkFileTypeCovariance, 
                                alternatePhenotypeParameters, 
                                *snpFilter );
      }
   else
      {
      plinkFile.ReadPlinkFiles( plinkFilenameRootCovarianceVector[ 0 ],
                                plinkFileTypeCovariance, 
                                alternatePhenotypeParameters, 
                                simFilterOpts );
      }

   //  
   // At this point we have read in the Plink data and stored it in a CPlinkFile
   //   This is the data that we will use to build the covariance / kernel
   //   data, but first the data needs some 'massaging'
   //  This is similar to the covariance processing in LoadDataAltPhen()
   ProgressNL( "  ++  Start Preprocessing Test Data:\n" );

   // Perform Eigenstrat Imputation on phenotype data
   zeroMeanColumns( plinkFile.phenotypeData, plinkFile.cIndividuals, plinkFile.cPhenotypes, false );

   // perform preprocessing on the (soon to be) covariance SNPs
   if ( simMissingness )
      {
      replaceSNPsByMissingness( plinkFile.snpData, plinkFile.cIndividuals, plinkFile.cSnps );
      }
   else
      {
      //Zero mean each individual and divide by sqrt(N_SNPs)
      ImputeAndZeroMeanSNPs( plinkFile.rgSnpInfo, plinkFile.snpData, plinkFile.cIndividuals, plinkFile.cSnps, unitVarianceSimilaritySNPs, zeroMeanSimilaritySNPs, normalizationInfo, usePseudocounts );
      preprocessSNPsForBuildingCovariance( plinkFile.snpData, plinkFile.cIndividuals, plinkFile.cSnps, 1 );
      }

   // covariate imputation in the plinkFile
   covariatesData.ImputeAndZeroMeanCovariatesFix( plinkFile.individualLabels );

   StudyData studydata;
   studydata.cIndividuals = plinkFile.cIndividuals;
   studydata.cPhenotypes = plinkFile.cPhenotypes;
   studydata.cSnps = plinkFile.cSnps;
   studydata.individualLabels = plinkFile.individualLabels;
   studydata.phenotypeLabels = plinkFile.phenotypeLabels;
   studydata.ind_phenotype = plinkFile.ind_phenotype;
   studydata.snpInfo = plinkFile.rgSnpInfo;
   studydata.phenotypeData = plinkFile.phenotypeData;
   studydata.snpData = plinkFile.snpData;

   InitCovariates( &studydata, &covariatesData );
   TrimAutoSelectSearchValues( autoSelectSearchValues, studydata.cSnps );

   vector< int > autoSelectSearchValuesThinned;
   vector <vector <double> *> Msqs;
   vector <vector <double> *> Mlls;

   //
   // Create a list to index all the individuals and then randomize it
   //
   size_t cIndividuals = plinkFile.cIndividuals;
   size_t cFolds = autoSelectFoldCount;         // Count of folds to create
   if ( cFolds > cIndividuals )
      {
      Warn( "-autoSelectFolds[%d] greater than number of individuals[%d].  Setting Fold count to [%d]", cFolds, cIndividuals, cIndividuals );
      cFolds = cIndividuals;
      }

   ProgressNL( "  ++            Start Fold Creation:" );
   Verbose( verboseMajor, "          Number of folds to Create: %9d", static_cast<int>(cFolds) );

// DO NOT USE srand() RANDOMLY 
//   You have to be able to reproduce your results
//   SO.. Allow the user to set the seed value
//   AND  display the seed value to the user so 
//         he knows how to reproduce the results
//  
//   We initialize rand() on program start up.
//
   vector< size_t > idxIndividuals;
   idxIndividuals.reserve( plinkFile.cIndividuals );
   generate_n( back_inserter( idxIndividuals ), cIndividuals, GenSequence( 0 ) );
   random_shuffle( idxIndividuals.begin(), idxIndividuals.end() );

   //
   // split the randomized list into cFolds sub-lists
   //
   vector< vector< size_t > > rgFolds;          // a vector of Folds (which is a vector of size_t indexes)
   rgFolds.resize( cFolds );
   GenRange range( cIndividuals, cFolds );

   //float cIndividualsPerFold = (float)cIndividuals/(float)cFolds;
   float cIndividualsPerFold = floor((float)cIndividuals/(float)cFolds);
   for ( size_t iFold = 0; iFold < cFolds; ++iFold )
      {
      range.SelectRange( iFold );
      copy( idxIndividuals.begin() + range.begin(), idxIndividuals.begin() + range.end(), back_inserter( rgFolds[ iFold ] ) ); 
      }

   // We have the folds, do we emit them for debugging
   char *pchDebugEmitFolds = getenv( "FastLmmDebugEmitFolds" );
   if ( pchDebugEmitFolds != nullptr )
      {
      FILE *ofile = fopen( "DebugEmitFolds.txt", "w" );
      if ( ofile == nullptr )
         {
         Fatal( "Unable to write Folds.txt file.  CRT Error: %d - %s", errno, strerror( errno ) );
         }
      fprintf( ofile, "CID\tFold\n" );
      for ( size_t iFold = 0; iFold < cFolds; ++iFold )
         {
         for ( vector< size_t >::iterator fold = rgFolds[iFold].begin(); fold != rgFolds[iFold].end(); ++fold )
            {
            string cid = plinkFile.individualLabels[ *fold ];
            fprintf( ofile, "%s\t%d\n", cid.c_str() , (int)iFold );              // iFold should never overflow an int...
            }
         }
      fclose( ofile );
      }

   //
   // We have the Folds, combine them 
   // 
   vector< size_t > testSet;
   vector< size_t > trainSet;
   testSet.reserve( (size_t)(cIndividualsPerFold) );
   trainSet.reserve( cIndividuals - (size_t)cIndividualsPerFold );

   ProgressNL( " Select SNPs using cross validation:" );
   CTimer timerCrossValidation( true );

   for ( size_t iFoldTest = 0; iFoldTest < cFolds; ++iFoldTest )
      {
      Verbose( verboseMajor, "              Start Processing Fold: %4d", iFoldTest + 1 );
      testSet.clear();
      trainSet.clear();

      for ( size_t iFold = 0; iFold < cFolds; ++iFold )
         {
         if ( iFoldTest == iFold )
            {
            // copy to testSet
            copy( rgFolds[ iFold ].begin(), rgFolds[ iFold ].end(), back_inserter( testSet ) );
            }
         else
            {
            // copy to trainSet
            copy( rgFolds[ iFold ].begin(), rgFolds[ iFold ].end(), back_inserter( trainSet ) );
            }
         }

      vector<double> *sqs = new vector<double>();
      vector<double> *lls = new vector<double>();

      StudyData *studyDataTest = new StudyData();
      StudyData *studyDataTrain = new StudyData();
      InitializeStudyDataFromPlinkFile( studyDataTrain, plinkFile, trainSet );
      InitializeStudyDataFromPlinkFile( studyDataTest, plinkFile, testSet );

      InitCovariates( studyDataTrain, &covariatesData );
      InitCovariates( studyDataTest, &covariatesData );

      vector< size_t > sorted;

      if ( fAutoSelectApplyThinning )
         {
         sorted = RunLinRegWithThinning( *studyDataTrain );
         Verbose( verboseMajor, "                SNPs After Thinning: %9d", sorted.size() );
         }
      else
         {
         sorted = RunLinReg( *studyDataTrain );
         }

      autoSelectSearchValuesThinned = autoSelectSearchValues;
      TrimAutoSelectSearchValues( autoSelectSearchValuesThinned, sorted.size() );
      size_t cAutoSelectSearchValues = autoSelectSearchValuesThinned.size();

      ResetSnpData( studyDataTrain, &sorted[0], autoSelectSearchValuesThinned.back() );
      ResetSnpData( studyDataTest,  &sorted[0], autoSelectSearchValuesThinned.back() );

      //copySNPDataToSNPcov( studyDataTrain );
      //copySNPDataToSNPcov( studyDataTest );
      studyDataTrain->snpDataCovariance=studyDataTrain->snpData;
      studyDataTest->snpDataCovariance=studyDataTest->snpData;

      studyDataTrain->snpInfo.resize( autoSelectSearchValuesThinned.back() );
      studyDataTest->snpInfo.resize( autoSelectSearchValuesThinned.back() );

      for ( size_t sl = 0; sl < cAutoSelectSearchValues; ++sl )
         {
         size_t nSNPcovTest = autoSelectSearchValuesThinned[ sl ];

         // we now have testSet and trainSet individuals selected.
         //   Populate StudyData to operate on for the next stage of processing. 
         //   Globals may be a problem.  
         //
         if ( nSNPcovTest == 0 )
            {
            const BLAS_INT *nTrainI = (const BLAS_INT *)&studyDataTrain->cIndividuals;
            const BLAS_INT *nTestI = (const BLAS_INT *)&studyDataTest->cIndividuals;
            const BLAS_INT *nCovar = (const BLAS_INT *)&studyDataTrain->cCovariates;
            const BLAS_INT *nCovars = (const BLAS_INT *)&studyDataTest->cCovariates;

            real *TrainPh = studyDataTrain->phenotypeData;
            real *TestPh = studyDataTest->phenotypeData;
            real *Covar = studyDataTrain->covariatesData;
            real *Covars = studyDataTest->covariatesData;

            real *sts = new real[*nCovar * *nCovar];
            BlasLapack::gemm( &char_T, &char_N, nCovar, nCovar, nTrainI, &real_one, Covar, nTrainI, Covar, nTrainI, &real_zero, sts, nCovar);

            real *sty = new real[*nCovar];
            BlasLapack::gemv( &char_T, nTrainI, nCovar, &real_one, Covar, nTrainI, TrainPh, &int_one, &real_zero, sty, &int_one);

            // Cholesky Decomposition
            int *info = new int[1];
            BlasLapack::potrf( &char_L, nCovar, sts, nCovar, info);

            // Solve the linear system (L L^T)^-1 cen_ph
            BlasLapack::potrs( &char_L, nCovar, &int_one, sts, nCovar, sty, nCovar, info );
            real *solve_sts_sty = sty;
            real *beta = sty;

            real *est_mu = new real[ *nTestI ];
            //REVIEW:  Eun,  I found this call to gemv where you are using Covars which is allocated as nCovars sized
            //         elements and nCovar (without the 's').  However nCovars and nCovar appear to always have the 
            //         same value so I wanted to have you review the change befor committing it.
            //         For now, IF the sizes are different, I want to know about it.
            //
            if ( *nCovars != *nCovar )
               {
               Error( "Expected *nCovars == *nCovar.  Found %d != %d", *nCovars, *nCovar );
               }

            BlasLapack::gemv( &char_N, nTestI, nCovar, &real_one, Covars, nTestI, solve_sts_sty, &int_one, &real_zero, est_mu, &int_one ); 

            real msqerr = ComputeMSQERR( TestPh, est_mu, *nTestI );
            //cout << nSNPcovTest << "    " << msqerr << endl;
            sqs->push_back(msqerr);

            // Now estimate the sigmae2  1/(N-D)*(y-sb)^T(y-sb)
            real *sb = new real[*nTrainI];
            //printMatrix( Covar, *nTrainI, *nCovar, "Covar" );
            //printVector( beta, *nCovar, "beta" );
            BlasLapack::gemv( &char_N, nTrainI, nCovar, &real_one, Covar, nTrainI, beta, &int_one, &real_zero, sb, &int_one );

            vectorSubtractBack( TrainPh, sb, *nTrainI); //y-sb
            double sigmae2 = squared_sum(sb, *nTrainI);
            sigmae2 = sigmae2/(*nTrainI - *nCovar);


            // now we compute LL with TestPh, est_mu, sigma
            //term1 
            real term1 = ( -1.0 / 2.0 ) * (*nTestI) *_log2pi;
            real term2 = ( -1.0 / 2.0 ) * (*nTestI) * log( sigmae2 );
            real term3 = 0;
            term3 = ( -1.0 / 2.0 ) * msqerr /sigmae2;
            real LL = term1 + term2 + term3;

            delete[] sts;
            delete[] sty;
            delete[] est_mu;
            delete[] info;
            delete[] sb;

            //cout << "LL(" << nSNPcovTest << ")" << " =  " << LL << endl;
            lls->push_back( LL );
            continue;
            }

         ResultsLMM results;

         studyDataTrain->cSnpsCovariance = nSNPcovTest;
         studyDataTest->cSnpsCovariance = nSNPcovTest;
         EigenSym *eigen = new EigenSym();
         useWoodbury = false;
         computeEigenSym_SnpSample( "eigen.txt", studyDataTrain, eigen );
         runAlternativeModels = false;
        
         if ( nSNPcovTest >= studyDataTrain->cIndividuals )
            {
            FastLMM_fullrank LMM_fullrank = FastLMM_fullrank( eigen );

            //As the Kernel gets bigger with the number of SNPs in the kernel the optimal delta goes up linearly with the number of SNPs in the kernel as well
            LMM_fullrank.minLogVal=minLogValGlobal+log((real)nSNPcovTest);
            LMM_fullrank.maxLogVal=maxLogValGlobal+log((real)nSNPcovTest);

            LMM_fullrank.set_Phenotype( studyDataTrain->phenotypeData, studyDataTrain->phenotypeLabels[studyDataTrain->ind_phenotype] );
            //LMM_fullrank.set_SNPs( studyDataTrain->snpData, n_SNPs, studyDataTrain->snpInfo );
            //LMM_fullrank.set_SNPsCov(studyDataTrain->snpDataCovariance, studyDataTrain->cSnpsCovariance, studyDataTrain->snpInfoCovariance);
            LMM_fullrank.set_covariates( studyDataTrain->covariatesData, studyDataTrain->cCovariates );

            LMM_fullrank.init();
            LMM_fullrank.trainNull();
            //LMM_fullrank.updateExcludeList();

            //get the results
            LMM_fullrank.getResultsNullOnly( results );

            real log_delta = results.delta_0[0];
            real sigmag2 = results.sigg2_0[0];


            real delta = exp(log_delta);
            const int *nSNP = (const int *)&nSNPcovTest;
            const BLAS_INT *nTrainI = (const BLAS_INT *)&studyDataTrain->cIndividuals;
            const BLAS_INT *nTestI = (const BLAS_INT *)&studyDataTest->cIndividuals;
            const BLAS_INT *nCovar = (const BLAS_INT *)&studyDataTrain->cCovariates;
            const BLAS_INT *nCovars = (const BLAS_INT *)&studyDataTest->cCovariates;
            real *beta = results.beta_0;


            real *TrainPh = studyDataTrain->phenotypeData;
            real *TestPh = studyDataTest->phenotypeData;
            real *Covar = studyDataTrain->covariatesData;
            real *Covars = studyDataTest->covariatesData;

            real *X = studyDataTrain->snpDataCovariance;
            real *Xs = studyDataTest->snpDataCovariance;
             
            // real *U and real *S
            real *S = new real[ *nTrainI ];
            real *U = eigen->U;  // #individual by rankU
   

            for ( int i = 0; i < *nTrainI; ++i )
               {
               *(S+i) = 1.0/(*(eigen->S+i) + exp(log_delta));
               //cout << *(S+i) << endl;
               }
            // Now S become (S + deltaI)-1
   
            // First Part
            real *ssb = new real[ *nTestI ];

            BlasLapack::gemv( &char_N, nTestI, nCovars, &real_one, Covars, nTestI, beta, &int_one, &real_zero, ssb, &int_one);

            //printMatrix(Covars, *nTestI, *nCovars, "Covars");
            //printVector(beta, *nCovar, "beta");
            //printVector(ssb, *nTestI, "ssb");

            // Second Part
            real *xsxt = new real[*nTestI * *nTrainI];
            real *us = new real[*nTrainI * *nTrainI];
            real *usut = new real[*nTrainI * *nTrainI];
            real *xsxtMusut = new real[*nTestI * *nTrainI];
            real *sb = new real[*nTrainI];  //y_sb
            real *xsxtMusutMy_sb = new real[*nTestI];

            BlasLapack::gemm( &char_N, &char_T, nTestI, nTrainI, nSNP, &real_one, Xs, nTestI, X, nTrainI, &real_zero, xsxt, nTestI); 
            multi_mat_diag(U, *nTrainI, *nTrainI, S, us); 
            BlasLapack::gemm( &char_N, &char_T, nTrainI, nTrainI, nTrainI, &real_one, us, nTrainI, U, nTrainI, &real_zero, usut, nTrainI);
            BlasLapack::gemm(&char_N, &char_N, nTestI, nTrainI, nTrainI, &real_one, xsxt, nTestI, usut, nTrainI, &real_zero, xsxtMusut, nTestI);
            BlasLapack::gemv( &char_N, nTrainI, nCovar, &real_one, Covar, nTrainI, beta, &int_one, &real_zero, sb, &int_one); // S beta
            real *y_sb = sb;
            vectorSubtractBack(TrainPh, sb, *nTrainI); // y - S beta
            BlasLapack::gemv( &char_N, nTestI, nTrainI, &real_one, xsxtMusut, nTestI, y_sb, &int_one, &real_zero, xsxtMusutMy_sb, &int_one);

            //Sum
            real *est_mu = new real[*nTestI];
            vectorAddition(ssb, xsxtMusutMy_sb, est_mu, *nTestI);
            //printVector(est_mu, *nTestI, "est_mu");

            double msqerr = ComputeMSQERR(TestPh, est_mu, *nTestI);
            //cout << nSNPcovTest << "    " << msqerr << endl;
            sqs->push_back(msqerr);

            delete[] ssb;
            delete[] us;
            delete[] usut;
            delete[] xsxtMusut;
            delete[] sb;
            delete[] xsxtMusutMy_sb;
            //delete[] est_mu;

            // Computing Sigma
            // part 1
            real *xsxst = new real[*nTestI * *nTestI];
            //printMatrix(Xs, *nTestI, *rU, "Xs");
            BlasLapack::gemm( &char_N, &char_T, nTestI, nTestI, nSNP, &real_one, Xs, nTestI, Xs, nTestI, &real_zero, xsxst, nTestI); // X* X*T
            //printMatrix(xsxst, *nTestI, *nTestI, "xsxst");
            
            add_diag_mat(xsxst, *nTestI, delta, xsxst); // X* X*^T + delta I
            //printMatrix(xsxst, *nTestI, *nTestI, "xsxst + deltaI");

            scale_matrix(xsxst, *nTestI, *nTestI, sigmag2);
            real *part1 = xsxst;
            //printMatrix(part1, *nTestI, *nTestI, "part1");

            //part 2
            //xsxt -> previously computed
            real *xsxtu = new real[*nTestI * *nTrainI];
            //printMatrixF(xsxt, *nTestI, *nTrainI, "xsxt.tab");
            //printMatrixF(U, *nTrainI, *nTrainI, "U_cur.tab");
            BlasLapack::gemm(&char_N, &char_N, nTestI, nTrainI, nTrainI, &real_one, xsxt, nTestI, U, nTrainI, &real_zero, xsxtu, nTestI); // (X* X^T) U
            //printMatrixF(xsxtu, *nTestI, *nTrainI, "xsxtu.tab");
            
            real *xsxtus = new real[*nTestI * *nTrainI];
            multi_mat_diag(xsxtu, *nTestI, *nTrainI, S, xsxtus);
            //printMatrixF(xsxtus, *nTestI, *nTrainI, "xsxtus.tab");

            real *xsxtusMxsxtu = new real[*nTestI * *nTestI];
            BlasLapack::gemm( &char_N, &char_T, nTestI, nTestI, nTrainI, &real_one, xsxtus, nTestI, xsxtu, nTestI, &real_zero, xsxtusMxsxtu, nTestI);
            //printMatrixF(xsxtusMxsxtu, *nTestI, *nTestI, "xsxtusMxsxtu");
            scale_matrix(xsxtusMxsxtu, *nTestI, *nTestI, -sigmag2);
            real *part2 = xsxtusMxsxtu;

            add_matrix(part1, part2, *nTestI, *nTestI, part1);
            real *sigma = part1;

            delete[] xsxt;
            delete[] xsxtusMxsxtu;
            delete[] xsxtus;
            delete[] S;

            // now we compute LL with TestPh, est_mu, sigma
            //term1 
            real term1 = -0.5* *nTestI *_log2pi;
            //cout << "term1 " << term1 << endl;


            // subtract mean
            real *cen_ph = new real[*nTestI];
            real *cen_phT = new real[*nTestI];
            vectorSubtract(TestPh, est_mu, cen_ph, *nTestI);

            for (int i = 0; i < *nTestI; ++i) {
                  *(cen_phT+i) = *(cen_ph+i);
               }

            //printMatrix(sigma, *nTestI, *nTestI, "sigma");
            // Cholesky Decomposition
            int *info = new int[1];
            BlasLapack::potrf( &char_L, nTestI, sigma, nTestI, info);  

            //compute log(det)
            real logdet = ComputeLogDet(sigma, *nTestI);
            //printMatrix(sigma, *nTestI, *nTestI, "chol(sigma)");
            real term2 = -0.5*logdet;

            // Solve the linear system (L L^T)^-1 cen_ph
            BlasLapack::potrs( &char_L, nTestI, &int_one, sigma, nTestI, cen_ph, nTestI, info);
            real term3 = 0;
            for (int i = 0; i < *nTestI; ++i) {
                  term3 = term3 + *(cen_ph+i) * *(cen_phT+i);
                  //cout << *(cen_ph+i) << " " << *(cen_phT+i) << endl;
               }
            term3 = -0.5*term3;

            real LL = term1 + term2 + term3;

            delete[] sigma;
            delete[] est_mu;
            delete[] cen_ph;
            delete[] cen_phT;
            delete[] info;

            //cout << "LL(" << nSNPcovTest << ")" << " =  " << LL << endl;
            lls->push_back(LL);
            }
         else
            {
            size_t rankU = eigen->rankU;
   
            FastLMM_lowrank LMM_lowrank = FastLMM_lowrank( eigen );
            LMM_lowrank.minLogVal = minLogValGlobal + log( (real)nSNPcovTest );
            LMM_lowrank.maxLogVal = maxLogValGlobal + log( (real)nSNPcovTest );

            LMM_lowrank.set_Phenotype( studyDataTrain->phenotypeData, studyDataTrain->phenotypeLabels[studyDataTrain->ind_phenotype] );
            //LMM_lowrank.set_SNPs( studyDataTrain->snpData, n_Snps, studyDataTrain->snpInfo );
            LMM_lowrank.set_covariates( studyDataTrain->covariatesData, studyDataTrain->cCovariates );
            //LMM_lowrank.set_SNPsCov( studyDataTrain->snpDataCovariance, studyDataTrain->cSnpsCovariance, studyDataTrain->snpInfoCovariance );

            LMM_lowrank.init();
            LMM_lowrank.trainNull();
            LMM_lowrank.updateExcludeList();
            LMM_lowrank.getResultsNullOnly( results );

            real log_delta = results.delta_0[0];
            real sigmag2 = results.sigg2_0[0];
            //cout << "sigma e^2 : " << sigmae2 << endl;
            //cout << "sigma g^2 : " << sigmag2 << endl;

            real delta = exp(log_delta);
            real inv_delta = 1.0/delta;
            const real *inv_deltap = (const real *)&inv_delta;
            const BLAS_INT *rU = (const BLAS_INT *)&rankU;
            const BLAS_INT *nTrainI = (const BLAS_INT *)&studyDataTrain->cIndividuals;
            const BLAS_INT *nTestI = (const BLAS_INT *)&studyDataTest->cIndividuals;
            const BLAS_INT *nCovar = (const BLAS_INT *)&studyDataTrain->cCovariates;
            const BLAS_INT *nCovars = (const BLAS_INT *)&studyDataTest->cCovariates;
            real *beta = results.beta_0;
            real *TrainPh = studyDataTrain->phenotypeData;
            real *TestPh = studyDataTest->phenotypeData;
            real *Covar = studyDataTrain->covariatesData;
            real *Covars = studyDataTest->covariatesData;
            real *X = studyDataTrain->snpDataCovariance;
            real *Xs = studyDataTest->snpDataCovariance;
             
            // real *U and real *S
            real *S = new real[ rankU ];
            real *U = eigen->U;  // #individual by rankU
   
            for ( size_t i = 0; i < eigen->rankU; ++i )
               {
               // *(S+i) = 1/(*(eigen->S+i) + exp(log_delta));
               S[ i ] = 1.0 / ( eigen->S[ i ] + exp( log_delta ) );
               }

            // Now S become (S + deltaI)-1
            // First Part
            real *ssb = new real[*nTestI];

            BlasLapack::gemv( &char_N, nTestI, nCovars, &real_one, Covars, nTestI, beta, &int_one, &real_zero, ssb, &int_one);

            //printMatrix(Covars, *nTestI, *nCovars, "Covars");
            //printVector(beta, *nCovar, "beta");
            //printVector(ssb, *nTestI, "ssb");

            // Second Part
            real *u1ty = new real[rankU];
            real *u1ts = new real[rankU*studyDataTrain->cCovariates];
            real *u1tsb = new real[rankU];
            real *u1ty_u1tsb = new real[rankU];
            real *xtu1 = new real[rankU*rankU];
            real *xsxtu1 = new real[studyDataTest->cIndividuals*rankU];
            real *xsxtu1s = new real[studyDataTest->cIndividuals*rankU];
            real *xsxtu1sMu1ty_u1tsb = new real[studyDataTest->cIndividuals];

            BlasLapack::gemv( &char_T, nTrainI, rU, &real_one, U, nTrainI, TrainPh, &int_one, &real_zero, u1ty, &int_one);  // U1^T y
            //printMatrix(U, *nTrainI, *rU, "U");
            //printVector(TrainPh, *nTrainI, "TrainPh");
            //printVector(u1ty, *rU, "u1ty");

            BlasLapack::gemm( &char_T, &char_N, rU, nCovar, nTrainI, &real_one, U, nTrainI, Covar, nTrainI, &real_zero, u1ts, rU); // U1^T S
            //printMatrix(u1ts, *rU, *nCovar, "u1ts");
            BlasLapack::gemv( &char_N, rU, nCovar, &real_one, u1ts, rU, beta, &int_one, &real_zero, u1tsb, &int_one); //U1^T S beta 
            //printVector(u1tsb, *rU, "u1tsb");
            vectorSubtract(u1ty, u1tsb, u1ty_u1tsb, rankU);  // U1^T Y - U1^T S beta
            //printVector(u1ty_u1tsb, *rU, "u1ty_u1tsb");
            BlasLapack::gemm( &char_T, &char_N, rU, rU, nTrainI, &real_one, X, nTrainI, U, nTrainI, &real_zero, xtu1, rU); // X^T U1
            //printMatrix(xtu1, *rU, *rU, "xtu1");
            BlasLapack::gemm( &char_N, &char_N, nTestI, rU, rU, &real_one, Xs, nTestI, xtu1, rU, &real_zero, xsxtu1, nTestI); // X* (X^T U1)
            //printMatrix(xsxtu1, *nTestI, *rU, "xsxtu1");
            multi_mat_diag(xsxtu1, *nTestI, *rU, S, xsxtu1s);
            //printMatrix(xsxtu1s, *nTestI, *rU, "xsxtu1s");

            BlasLapack::gemv( &char_N, nTestI, rU, &real_one, xsxtu1s, nTestI, u1ty_u1tsb, &int_one, &real_zero, xsxtu1sMu1ty_u1tsb, &int_one); // [(X* X^T U1) * S1^-1] x [U1^T Y - U1^T S beta]
            //printVector(xsxtu1sMu1ty_u1tsb, *nTestI, "xsxtu1sMu1ty_u1tsb");

            // Third Part
            real *sb = new real[*nTrainI];
            BlasLapack::gemv( &char_N, nTrainI, nCovar, &real_one, Covar, nTrainI, beta, &int_one, &real_zero, sb, &int_one); // S beta
            real *y_sb = sb;
            vectorSubtractBack(TrainPh, sb, *nTrainI); // y - S beta
            real *xty_sb = new real[*rU];
            BlasLapack::gemv( &char_T, nTrainI, rU, &real_one, X, nTrainI, y_sb, &int_one, &real_zero, xty_sb, &int_one); // X^T (y - S beta)
            real *xsxty_sb = new real[*nTestI];
            BlasLapack::gemv( &char_N, nTestI, rU, &real_one, Xs, nTestI, xty_sb, &int_one, &real_zero, xsxty_sb, &int_one); // X* (X^T (y - S beta))

            real *u1ty_sb = new real[*rU];
            BlasLapack::gemv( &char_T, nTrainI, rU, &real_one, U, nTrainI, y_sb, &int_one, &real_zero, u1ty_sb, &int_one); // U1^T (y - S beta)
            real *xsxtu1Mu1ty_sb = new real[*nTestI];
            //printMatrix(xsxtu1, *nTestI, *rU, "xsxtu1");
            //printVector(u1ty_sb, *rU, "u1ty_sb");
            BlasLapack::gemv( &char_N, nTestI, rU, &real_one, xsxtu1, nTestI, u1ty_sb, &int_one, &real_zero, xsxtu1Mu1ty_sb, &int_one); //  X*(X^T U1) * U1^T (y - S beta)
            //printVector(xsxtu1Mu1ty_sb, *nTestI, "xsxtu1Mu1ty_sb");
            real *xsxty_sb__xsxtu1Mu1ty_sb = xsxtu1Mu1ty_sb;
            vectorSubtractBack(xsxty_sb, xsxty_sb__xsxtu1Mu1ty_sb, *nTestI);
            BlasLapack::scal(nTestI, inv_deltap, xsxty_sb__xsxtu1Mu1ty_sb, &int_one);
            //printVector(xsxtu1Mu1ty_sb, *nTestI, "xsxtu1Mu1ty_sb");

            //Sum
            real *est_mu = new real[*nTestI];
            vectorAddition(ssb, xsxtu1sMu1ty_u1tsb, xsxty_sb__xsxtu1Mu1ty_sb, est_mu, *nTestI);
            //printVector(est_mu, *nTestI, "est_mu");
            //printVector(TestPh, *nTestI, "test phenotype");

            double msqerr = ComputeMSQERR(TestPh, est_mu, *nTestI);
            //cout << nSNPcovTest << "    " << msqerr << endl;
            sqs->push_back(msqerr);

            delete[] S;
            delete[] ssb;
            delete[] u1ty;
            delete[] u1ts; 
            delete[] u1tsb;
            delete[] u1ty_u1tsb;
            delete[] xtu1;
            ////delete[] xsxtu1; 
            ////delete[] xsxtu1s;
            delete[] xsxtu1sMu1ty_u1tsb;
            delete[] sb;
            delete[] xty_sb;
            delete[] xsxty_sb;
            delete[] u1ty_sb;
            delete[] xsxtu1Mu1ty_sb;
            //delete[] est_mu;

            //continue;
            // Now Computing Sigma

            // part 1
            real *xsxst = new real[*nTestI * *nTestI];
            //printMatrix(Xs, *nTestI, *rU, "Xs");
            BlasLapack::gemm( &char_N, &char_T, nTestI, nTestI, rU, &real_one, Xs, nTestI, Xs, nTestI, &real_zero, xsxst, nTestI); // X* X*T
            //printMatrix(xsxst, *nTestI, *nTestI, "xsxst");

            add_diag_mat(xsxst, *nTestI, delta, xsxst); // X* X*^T + delta I

            scale_matrix(xsxst, *nTestI, *nTestI, sigmag2);
            real *part1 = xsxst;

            //printMatrix(part1, *nTestI, *nTestI, "part1");

            //part2 
            //xtu1 previously computed 
            //xsxtu1 previously computed   <---- used alot
            //xsxtu1s * (xsxtu1)^T
            real *xsxtu1sMxsxtu1t = new real[*nTestI * *nTestI];
            BlasLapack::gemm( &char_N, &char_T, nTestI, nTestI, rU, &real_one, xsxtu1s, nTestI, xsxtu1, nTestI, &real_zero, xsxtu1sMxsxtu1t, nTestI);
            scale_matrix(xsxtu1sMxsxtu1t, *nTestI, *nTestI, -1*sigmag2);
            real *part2 = xsxtu1sMxsxtu1t;

            //printMatrix(part2, *nTestI, *nTestI, "part2");

            //part3
            real *xtx = new real[*rU * *rU];
            real *xsxtx = new real[*nTestI * *rU];
            real *xsxtxxst = new real[*nTestI * *nTestI];
            BlasLapack::gemm( &char_T, &char_N, rU, rU, nTrainI, &real_one, X, nTrainI, X, nTrainI, &real_zero, xtx, rU);
            BlasLapack::gemm( &char_N, &char_N, nTestI, rU, rU, &real_one, Xs, nTestI, xtx, rU, &real_zero, xsxtx, nTestI);
            BlasLapack::gemm( &char_N, &char_T, nTestI, nTestI, rU, &real_one, xsxtx, nTestI, Xs, nTestI, &real_zero, xsxtxxst, nTestI);
            scale_matrix(xsxtxxst, *nTestI, *nTestI, -sigmag2/delta);
            real *part3 = xsxtxxst;

            //printMatrix(part3, *nTestI, *nTestI, "part3");

            //part4
            real *xsxtu1Mxsxtu1t = new real[*nTestI * *nTestI];
            BlasLapack::gemm( &char_N, &char_T, nTestI, nTestI, rU, &real_one, xsxtu1, nTestI, xsxtu1, nTestI, &real_zero, xsxtu1Mxsxtu1t, nTestI);
            scale_matrix(xsxtu1Mxsxtu1t, *nTestI, *nTestI, sigmag2/delta);
            real *part4 = xsxtu1Mxsxtu1t;

            //printMatrix(part4, *nTestI, *nTestI, "part4");

            add_square_matrix(part1, part2, part3, part4, *nTestI, part1);
            real *sigma = part1;

            delete[] xsxtu1; 
            delete[] xsxtu1s;
            delete[] xsxtu1sMxsxtu1t;
            delete[] xtx;
            delete[] xsxtx;
            delete[] xsxtxxst;


            // now we compute LL with TestPh, est_mu, sigma
            //term1 
            real term1 = -1.0 / 2.0 * *nTestI * _log2pi;

            // subtract mean
            real *cen_ph = new real[ *nTestI ];
            real *cen_phT = new real[ *nTestI ];
            vectorSubtract( TestPh, est_mu, cen_ph, *nTestI );

            for ( int i = 0; i < *nTestI; ++i )
               {
               cen_phT[ i ] = cen_ph[ i ];
               }

            // Cholesky Decomposition
            int *info = new int[1];
            BlasLapack::potrf( &char_L, nTestI, sigma, nTestI, info);  

            //compute log(det)
            real logdet = ComputeLogDet( sigma, *nTestI );
            //printMatrix(sigma, *nTestI, *nTestI, "chol(sigma)");
            real term2 = -0.5 * logdet;

            // Solve the linear system (L L^T)^-1 cen_ph
            BlasLapack::potrs( &char_L, nTestI, &int_one, sigma, nTestI, cen_ph, nTestI, info);
            real term3 = 0;
            for ( int i = 0; i < *nTestI; ++i )
               {
               term3 = term3 + *(cen_ph+i) * *(cen_phT+i);
               //cout << *(cen_ph+i) << " " << *(cen_phT+i) << endl;
               }

            term3 = -0.5 * term3;
            real LL = term1 + term2 + term3;

            delete[] sigma;
            delete[] est_mu;
            delete[] cen_ph;
            delete[] cen_phT;
            delete[] info;

            //cout << "LL(" << nSNPcovTest << ")" << " =  " << LL << endl;
            lls->push_back( LL );
            }
         }
      delete studyDataTest;
      delete studyDataTrain;

#if defined( OUTPUT_MEAN_SQUARED_ERROR )
      cout << "mean sqared error" << endl;
      for ( int i = 0; i < sqs->size(); ++i )
         {
         cout << sqs[ i ]) << " ";
         }
      cout << endl;
#endif

      Msqs.push_back(sqs);

#if defined( OUTPUT_LOG_LIKELIHOOD )
      cout << "Log Likelihood" << endl;
      for ( int i = 0; i < lls->size(); ++i )
         {
         cout << lls[ i ] << " ";
         }
      cout << endl;
#endif
      Mlls.push_back( lls );
      }

   ProgressNL( "             Cross Validataion Done:" );
   timerCrossValidation.Stop();
   timerCrossValidation.Report( -1, "      Cross Validation Elapsed Time: %s" );

   // we expect these vectors to be initialized to 0.0 and
   //   since we are creating the vector< double > x() with 
   //   a size initializer, the constructor has already 
   //   initialized them to 0.0 
   //
   // Average Mean Squared Error and Likelihood
   vector< double > averageMeanSquaredError( Mlls[0]->size() );
   vector< double > averageLogLikelihood( Mlls[0]->size() );
   vector< vector< double > *>::iterator itr;
   for ( itr = Msqs.begin(); itr < Msqs.end(); ++itr )
      {
      vector< double > *vitr = (vector< double > *)*itr;
      for ( size_t i = 0; i < vitr->size(); ++i )
         {
         averageMeanSquaredError[ i ] = averageMeanSquaredError[ i ] + (*vitr)[ i ];
         }
      }

   for ( size_t i = 0; i < averageMeanSquaredError.size(); ++i )
      {
      averageMeanSquaredError[ i ] = averageMeanSquaredError[ i ] / cFolds;
      }

   // Average Likelihood
   for ( itr = Mlls.begin(); itr < Mlls.end(); ++itr )
      {
      vector< double > *vitr = (vector< double > *)*itr;
      for ( size_t i = 0; i < vitr->size(); ++i )
         {
         averageLogLikelihood[ i ] = averageLogLikelihood[ i ] + (*vitr)[ i ]; 
         }
      }

   for ( size_t i = 0; i < averageLogLikelihood.size(); ++i )
      {
      averageLogLikelihood[ i ] = averageLogLikelihood[ i ] / cFolds;
      }

   size_t idxSelected = 0;    // Which 'bucket' holds the right number of SNPs to use
                              //   initialize as the first element in the loop and
                              //   find the right one.
   /*
    * use MSE to determine the number of SNPs to use but 
    *   keep max_ll around for now in case we change the 
    *   behavior again
    */
   if ( autoSelectCriterion == useMSE )
      {
      /*
       *  Find the min MSE which determines the number of SNPs
       *    to use
       */
      double min_avgMse = averageMeanSquaredError[0];       // load initial values for test loop
      for ( size_t i = 1; i < averageMeanSquaredError.size(); ++i )
         {
         if ( min_avgMse > averageMeanSquaredError[ i ] )
            {
            min_avgMse = averageMeanSquaredError[ i ];
            idxSelected = i;
            }
         }
      }
   else if ( autoSelectCriterion == useLL )
      {
      /*
       *  Find the max LL which determines the number of SNPs
       *    to use
       */
      double max_ll = averageLogLikelihood[0];         // load initial values for test loop
      for ( size_t i = 1; i < averageLogLikelihood.size(); ++i )
         {
         if ( max_ll < averageLogLikelihood[ i ] )
            {
            max_ll = averageLogLikelihood[ i ];
            idxSelected = i;
            }
         }
      }
   else
      {
      Fatal( "Unknown value for autoSelectCriterion.  Found [%d]", autoSelectCriterion );
      }

   vector< size_t > selectedSnps;
   if ( fAutoSelectApplyThinning )
      {
      selectedSnps = RunLinRegWithThinning( studydata );
      }
   else
      {
      selectedSnps = RunLinReg( studydata );
      }

   // Use the original autoSelectSearchValues vector to get the 'k' value
   size_t max_nsnp = std::min( static_cast<size_t>(autoSelectSearchValues[ idxSelected ]), selectedSnps.size() );
   ProgressNL( "            Number of SNPs Selected: %7d", max_nsnp );

   // If we need to output a PED file with selected SNPs, do it now
   if ( fWritePedFile )
      {
      if ( max_nsnp > 0 )
         {

         // Create the list of SnpIds
         Verbose( verboseMajor, "\n        Re-read to produce PED file:" );
         vector< SnpId > rgSnpIds;
         rgSnpIds.reserve( max_nsnp );
         for( size_t iSnp=0; iSnp<max_nsnp; ++iSnp )
            {
            rgSnpIds.push_back( studydata.snpInfo[ selectedSnps[ iSnp ] ].idSnp );
            }

         SnpFilterOptions snpFilterOptions;
         snpFilterOptions.SetByListOfSnpIds( string("$AutoSelect$"), rgSnpIds );
         CPlinkFile plinkFileT;
         plinkFileT.ReadPlinkFiles( plinkFilenameRootCovarianceVector[0], plinkFileTypeCovariance, alternatePhenotypeParameters, snpFilterOptions );
         ProgressNL( "                  Write pedOutFiles: [%s]", pedOutputFile.c_str() );
         WriteSelectedSnpsAsPedFile( pedOutputFile, plinkFileT );

#if 0 // the original code
         topk = max_nsnp;           // SelectTopK() uses global variable topk.  Set it from AutoSelect
         SnpFilterOptions *pSnpFilterOptions = SelectTopK();
         CPlinkFile plinkFileT;
         plinkFileT.ReadPlinkFiles( plinkFilenameRootCovarianceVector[0], plinkFileTypeCovariance, alternatePhenotypeParameters, *pSnpFilterOptions );
         ProgressNL( "                  Write pedOutFiles: [%s]", pedOutputFile.c_str() );
         WriteSelectedSnpsAsPedFile( pedOutputFile, plinkFileT );
         delete pSnpFilterOptions;
#endif
         }
      else
         {
         Warn( "No PED file created because 0 SNPs selected." );
         }
      }

   WriteAutoSelectCrossValidationFile( autoSelectOutputFilenameRoot, selectedSnps.size(), averageMeanSquaredError, averageLogLikelihood );
   WriteAutoSelectSnpsFile( autoSelectOutputFilenameRoot, studydata.snpInfo, max_nsnp, selectedSnps );

   if ( !fullOutputFilename.empty() )
      {
      CTimer gwasTimer( true );

      if ( max_nsnp == 0 )
         {
         LinReg lin = LinReg( &studydata );
         lin.trainNull();
         lin.trainAlt();
         gwasTimer.Stop();
         gwasTimer.Report( -1, "                  GWAS elapsed time: %s" );
         ProgressNL( "                       Write output: [%s]", fullOutputFilename.c_str() );
         lin.WriteResultFile( fullOutputFilename );
         }
      else
         {
         ResultsLMM *results = new ResultsLMM();
         EigenSym *eigen = new EigenSym();
         useWoodbury = false;
         runAlternativeModels = true;

         SetSnpCovariance( &studydata, &selectedSnps[ 0 ], max_nsnp );
         computeEigenSym_SnpSample( "eigen.txt", &studydata, eigen );
         LmmGwas( &studydata, eigen, results );
         gwasTimer.Stop();
         gwasTimer.Report( -1, "                  GWAS elapsed time: %s" );
         ProgressNL( "                       Write output: [%s]", fullOutputFilename.c_str() );
         results->WriteResultFile( fullOutputFilename );
         }
      }
   }
