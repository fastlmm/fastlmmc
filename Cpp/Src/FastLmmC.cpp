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
 *         File Name:   FastLmmC.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements the FastLMM tool to C++
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Include Files
 */
#include "FastLmmC.h"

void DoIteratingFastLmm();
void InitializeDebugSupport();
void WriteEmptyTaskResultFile( const string& ofile, const SnpFilterOptions& sfo );

void printMatrixF( real *mat, size_t cRows, size_t cColumns, const string& name );

FILE *OpenResultFile( const string& resultFile );
void UpdateExcludeLists( ExcludeSNPfromSimType exclusionType, 
                         vector< SnpInfo >& rgSnpInfoTest, 
                         vector< SnpInfo >& rgSnpInfoCovariance, 
                         size_t *rgExcludeStart, 
                         size_t *rgNExclude );

/*
 * The code
 */

/*
 *  This is the main entrypoint to the program
 */
int main(int argc, char* argv[])
   {
   extern const char *szSplash;

   CTimer programTimer( true );

   InitializeDebugSupport();
   InitializeProgram();
   ProgressNL( szSplash, szFastLmmCVersion, MathLibName().c_str() );
   ProgressNL( " ++    Start Processing CommandLine:" );
   ProcessCommandLine( argc, argv );
   ProgressNL( " --      End Processing CommandLine:\n" );

   // CommandLine can change log file location so wait
   //  til after commandline processing to start logging.
   LogMemoryState( "Start MemoryState Log", OverWrite );

   if ( (snpFilterOpts.filterType == SnpFilterOptions::FilterByTask)
     && (snpFilterOpts.iTaskAsk >= snpFilterOpts.cTasksMax)
      )
      {
      //
      // This task is beyond the data but we need to produce the file
      //   so post processing of the files will find the number
      //   of files expected and work 'properly'
      //
      Warn( "No SNPs are associated with this task." );
      WriteEmptyTaskResultFile( fullOutputFilename, snpFilterOpts );
      }
   else if ( fAutoSelectSnps )
      {
      //
      // Automatically select the SNPs to build the kernel
      //
      ProgressNL( " ++            Start DoAutoSelect():" );
      DoAutoSelectSnps();
      ProgressNL( " --              End DoAutoSelect():\n" );

      LogMemoryState( "After DoAutoSelectSnps()" );
      }
   else if ( fDoIteratingFastLmm )
      {
      //
      // This option tells us to not load everything into memory 
      //   and then run.  Instead we will load the 'stuff'
      //   that needs a big matrix all at once and iterate
      //   over the stuff that doesn't need to be all in
      //   memory at one time
      //
      DoIteratingFastLmm();
      }
   else
      {
      // create and initialize the 'data' objects:
      StudyData *data = new StudyData();
      EigenSym *eigen = new EigenSym();
      ResultsLMM *results = new ResultsLMM();

      //load the data object
      ProgressNL( " ++      Start Loading FastLmm Data:" );
      LoadDataAltPhen( data );
      ProgressNL( " --        End Loading FastLmm Data:\n" );

      LogMemoryState( "After LoadDataAltPhen()" );

      if ( fDumpSnpFilterOptions )
         {
         string partitionLog = FileNameFromPath( plinkFilenameRoot );
         partitionLog += ".partition.log";
         snpFilterOpts.Dump( partitionLog );
         }

      if( linReg )
         {
         ProgressNL( "  Compute GWAS w/ Linear Regression:" );
         CTimer linRegTimer( true );
         LinReg lin = LinReg( data );

         LogMemoryState( "After LinReg()" );

         if ( !fDoSnpPairs )
            {
            lin.trainNull();
            lin.trainAlt();

            LogMemoryState( "After lin.trainAlt()" );
            }
         else
            {
            lin.trainPairs();

            LogMemoryState( "After lin.trainPairs()" );
            }

         linRegTimer.Stop();
         linRegTimer.Report( -1, "                LinReg elapsed time: %s" );
         ProgressNL( "                  Write output file: [%s]", fullOutputFilename.c_str() );
         if ( !fDoSnpPairs )
            {
            lin.WriteResultFile( fullOutputFilename );
            }
         }
      else if ( fUseLogisticRegression )
         {
#if defined( USE_LINUX_V2x_KERNEL )
         Fatal( "Logistic Regression is not supported on Linux v2.x kernel/toolset" );
#else
         ProgressNL( "Compute GWAS w/ Logistic Regression:" );
         CTimer logisticRegressionTimer( true );
         CLogisticRegression logReg = CLogisticRegression( data );

         LogMemoryState( "After CLogisticRegressions()" );

         logReg.TrainNull();
         logReg.TrainAlt();

         LogMemoryState( "After logReg.TrainAlt()" );

         logisticRegressionTimer.Stop();
         logisticRegressionTimer.Report( -1, "   Logistic Regression elapsed time: %s" );
         ProgressNL( "                  Write output file: [%s]", fullOutputFilename.c_str() );
         logReg.WriteResultFile( fullOutputFilename, fVerboseOutput, pvaluePrintThreshold );
#endif
         }
      else
         {
         // compute/load EigenSym
         // generate an eigenId from all the files that contribute to Eigen
         ProgressNL( "              Compute/Load EigenSym:" );
         string eigenid = alternatePhenotypeParameters.alternatePhenotypeFilename;
         for ( size_t i_file=0; i_file < plinkFilenameRootCovarianceVector.size(); ++i_file )
            {
            eigenid = eigenid.append( "_" );
            eigenid = eigenid.append( plinkFilenameRootCovarianceVector[i_file] );
            }

         fillEigenSym( eigenid, eigen, data );

         LogMemoryState( "After fillEigenSym()" );

         if ( kernelFileOut.length() )
            {
            eigen->saveReconstructionToDisk( data, kernelFileOut );
            }

         //permute the data as specified
         if( permutationMode != NO_PERM )
            {
            permuteStudyData( eigen, data );
            LogMemoryState( "After permuteStudyData()" );
            }

         //run GWAS
         if( runGwasType == RUN )
            {
            ProgressNL( "             Compute GWAS using LMM:" );

            CTimer gwasTimer( true );
            LmmGwas( data, eigen, results );

            LogMemoryState( "After LmmGwas()" );

            gwasTimer.Stop();
            gwasTimer.Report( -1, "                  GWAS elapsed time: %s" );

            if ( !fDoSnpPairs )
               {
               ProgressNL( "                  Write output file: [%s]", fullOutputFilename.c_str() );
               results->WriteResultFile( fullOutputFilename );
               }
            }
         }

      delete results;   results = nullptr;
      delete eigen;     eigen   = nullptr;
      delete data;      data    = nullptr;
      }

   programTimer.Stop();
   programTimer.Report( -1, "                 Total elapsed time: %s" );
   Progress( "\n" );

   LogMemoryState( "At program exit()" );
   return( 0 );
   }

/*
 *  Take all the steps needed to setup the program
 *    to support any internal debug or diagnostic
 *    capabilities
 */
void InitializeDebugSupport()
   {
#if defined( _DEBUG )
   if ( getenv( "FastLmmLogMemoryDiagnostics" ) != nullptr )
      {
      fLogMemoryDiagnostics = true;
      }
   if ( fLogMemoryDiagnostics )
      {
      _CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF      // enable debug allocation tracking
                  //  | _CRTDBG_DELAY_FREE_MEM_DF // keep freed blocks (consumes memory simulates low memory conditions)
                  //  | _CRTDBG_CHECK_ALWAYS_DF   // call _CrtCheckMemory at every allocation/de-allocation (VERY SLOW...but effective)
                  //  | _CRTDBG_CHECK_CRT_DF      // include _CRT_BLOCKs in memory check/leak detection operations
                  //  | _CRTDBG_LEAK_CHECK_DF     // set DEBUG version to validate memory on exit.
                    );
      _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );    // and emit CRT_WARNings to stderr rather than OutputDebugString()
      _CrtSetReportFile( _CRT_WARN, _CRTDBG_FILE_STDERR );
      }
   // To use heap debugging routines effectively, look for documentation on
   //   _CrtMemCheckpoint()
   //   _CrtMemDifference()
   //   _CrtMemDumpStatistics()
   //   _CrtMemDumpAllObjectsSince()
   //   _CrtDumpMemoryLeaks()
   //
#endif
   }

/*
 *  Load the data from files in prep for running FaSTLMM
 */
void LoadDataAltPhen( StudyData *studydata )
   {
   CCovariatesData *covariatesData = nullptr; //regression covariates

   if ( FUseCovariatesFile() )
      {
      covariatesData = new CCovariatesData();
      covariatesData->ReadCovariatesFile( covariatesFile );
      if ( covariatesData->n_covariates == 0 )
         {
         delete covariatesData;
         covariatesData = nullptr;
         }
      }

   if ( (computeEigenType == NOSAVE) || (computeEigenType == SAVETODISK) )
      {      
      if ( studydata->snpDataCovariance != nullptr )
         {
         Fatal( "snpDataCovariance has already been set. ");
         }

      if ( kernelFile.empty() 
        && !identityKernel 
        && !linReg 
        && !fUseLogisticRegression
        && (simGroups.empty() || groupsFile.empty())
        )
         {
         // Read the id, genotype, and phenotype data for each individual 
         //   in the study.  Covariance file does not participate in the 
         //   Snp filtering (-extract options etc.)
         CPlinkFile plinkFileCovariance;
         plinkFileCovariance.covariatesFileData = covariatesData;

         ProgressNL( "  ++  Start Loading Covariance Data:" );
         plinkFileCovariance.ReadPlinkFiles( plinkFilenameRootCovarianceVector[0], plinkFileTypeCovariance, alternatePhenotypeParameters, simFilterOpts );
         ProgressNL( "  --    End Loading Covariance Data:" );

         // create studydata and copy the right information into it.
         studydata->cIndividuals      = plinkFileCovariance.cIndividuals;
         studydata->cPhenotypes       = plinkFileCovariance.cPhenotypes;
         studydata->cSnpsCovariance   = plinkFileCovariance.cSnps;
         studydata->individualLabels  = std::move( plinkFileCovariance.individualLabels );
         studydata->phenotypeLabels   = std::move( plinkFileCovariance.phenotypeLabels );
         studydata->ind_phenotype     = plinkFileCovariance.ind_phenotype;
         studydata->snpInfoCovariance = std::move( plinkFileCovariance.rgSnpInfo );

         // transfer ownership of memory for phenotypeData and snpData from plinkFileCovariance to studydata
         studydata->phenotypeData          = plinkFileCovariance.phenotypeData;
         plinkFileCovariance.phenotypeData = nullptr;
         studydata->snpDataCovariance      = plinkFileCovariance.snpData;
         plinkFileCovariance.snpData       = nullptr;
         }
      else if ( linReg ){/*no kernel needed for linReg*/}
      else if ( fUseLogisticRegression ) { /* No kernel for Logistic Regression either */ }
      else if ( identityKernel )
         {
         Verbose( verboseMajor, "  **        Use the Identity Kernel:" );
         //load the phenotype
         CPlinkAlternatePhenotypeFile *alternatePhenotypeFile = new CPlinkAlternatePhenotypeFile( alternatePhenotypeParameters.alternatePhenotypeFilename );
         alternatePhenotypeFile->Read();
         if ( alternatePhenotypeParameters.selectedPhenotypeId.empty() )
            {
            alternatePhenotypeFile->SetSelectedPhenotype( alternatePhenotypeParameters.selectedPhenotypeIndex );
            }
         else
            {
            alternatePhenotypeFile->SetSelectedPhenotype( alternatePhenotypeParameters.selectedPhenotypeId );
            }

         Verbose( verboseMinor, "            Using Phenotype Index: %7d", 1 + alternatePhenotypeFile->SelectedPhenotypeIndex() );
         Verbose( verboseMinor, "             Using Phenotype Name: [%s]", alternatePhenotypeFile->SelectedPhenotypeName().c_str() );

         //DEBUG:
         if ( studydata->phenotypeLabels.size() != 0 )
            {
            Fatal( "Expected phenotypLabels.size() to be 0, and it is not!" );
            }
         // We only use 1 phenotype at the moment.  Set it up.
         studydata->phenotypeLabels.push_back( alternatePhenotypeFile->SelectedPhenotypeName() );
         studydata->ind_phenotype = 0;
         studydata->cPhenotypes = 1;

         //get the individualsMapping
         vector<string> individualsIndex = alternatePhenotypeFile->KeysForIndividualsWithValidPhenotypes( alternatePhenotypeParameters.selectedPhenotypeIndex );
         
         //load the kernelData
         studydata->cIndividuals = individualsIndex.size();
         studydata->kernelData = new real[ studydata->cIndividuals*studydata->cIndividuals ];
         for ( size_t i=0; i<studydata->cIndividuals*studydata->cIndividuals; ++i )
            {
            studydata->kernelData[ i ] = 0.0;
            }
         for ( size_t i=0; i<studydata->cIndividuals; ++i )
            {
            studydata->kernelData[ i + (studydata->cIndividuals*i) ] = 1.0;
            }
         }
      else if( !(groupsFile.empty() || simGroups.empty()) )
         {
         Verbose( verboseMajor, "                 Load groups file: [%s]", groupsFile.c_str() );
         GroupSim* sim = new GroupSim( alternatePhenotypeParameters, groupsFile );
         sim->loadKernel( simGroups );
         studydata->groupSim = sim;
         Verbose( verboseMajor, "                            Found: %i individuals in %i groups.", sim->cIndividuals, sim->cGroups );
         }
      else if( !kernelFile.empty() )
         {
         //load the phenotype
         Verbose( verboseMajor, "        Load the kernel from file: [%s]", kernelFile.c_str() );
         CPlinkAlternatePhenotypeFile *alternatePhenotypeFile = new CPlinkAlternatePhenotypeFile( alternatePhenotypeParameters.alternatePhenotypeFilename );
         alternatePhenotypeFile->Read();
         if ( alternatePhenotypeParameters.selectedPhenotypeId.empty() )
            {
            alternatePhenotypeFile->SetSelectedPhenotype( alternatePhenotypeParameters.selectedPhenotypeIndex );
            }
         else
            {
            alternatePhenotypeFile->SetSelectedPhenotype( alternatePhenotypeParameters.selectedPhenotypeId );
            }

         Verbose( verboseMajor, "            Using Phenotype Index: %7d", 1 + alternatePhenotypeFile->SelectedPhenotypeIndex() );
         Verbose( verboseMajor, "             Using Phenotype Name: [%s]", alternatePhenotypeFile->SelectedPhenotypeName().c_str() );

         //DEBUG:
         if ( studydata->phenotypeLabels.size() != 0 )
            {
            Fatal( "Expected phenotypLabels.size() to be 0, and it is not!" );
            }
         // We only use 1 phenotype at the moment.  Set it up.
         studydata->phenotypeLabels.push_back( alternatePhenotypeFile->SelectedPhenotypeName() );
         studydata->ind_phenotype = 0;
         studydata->cPhenotypes = 1;

         //get the individualsMapping
         vector<string> individualsIndex = alternatePhenotypeFile->KeysForIndividualsWithValidPhenotypes( alternatePhenotypeParameters.selectedPhenotypeIndex );
         
         //load the kernelData
         CKernelData* kerneldata = new CKernelData();
         kerneldata->ReadKernelFile( kernelFile );
         kerneldata->CreateColumnMajorArrayFromVectorsAndMapping( individualsIndex );
         studydata->kernelData = kerneldata->kernelArray;
         studydata->cIndividuals = individualsIndex.size();
         //studydata->individualLabels = alternatePhenotypeFile->phenotypeIndividualIdIndex;
         }
      }
   if ( runGwasType==RUN )
      {
      if ( studydata->snpData!=nullptr )
         {
         Fatal( "snpData has already been set." );
         }
      // Read the id, genotype, and phenotype data for each individual in the study
      CPlinkFile plinkFile;
      plinkFile.covariatesFileData = covariatesData;          // pass covariates information to SelectIndividuals
      ProgressNL( "  ++        Start Loading Test Data:" );
      plinkFile.ReadPlinkFiles( plinkFilenameRoot, plinkFileType, alternatePhenotypeParameters, snpFilterOpts );
      ProgressNL( "  --          End Loading Test Data:\n" );

      // create studydata and copy the right information into it.
      studydata->cIndividuals     = plinkFile.cIndividuals;
      studydata->cPhenotypes      = plinkFile.cPhenotypes;
      studydata->cSnps            = plinkFile.cSnps;
      studydata->individualLabels = std::move( plinkFile.individualLabels );
      studydata->phenotypeLabels  = std::move( plinkFile.phenotypeLabels );
      studydata->ind_phenotype    = plinkFile.ind_phenotype;
      studydata->snpInfo          = std::move( plinkFile.rgSnpInfo );

      // transfer ownership of memory for phenotypeData and snpData from plinkFile to studydata
      studydata->phenotypeData = plinkFile.phenotypeData; 
      plinkFile.phenotypeData  = nullptr;
      studydata->snpData       = plinkFile.snpData;
      plinkFile.snpData        = nullptr;
      }

   // Perform Eigenstrat Imputation:
   if ( (studydata->phenotypeData != nullptr) && (!useBiasTerm) )
      {  
      // Eigenstrat Imputation for the kernel SNPs
      // phenotypes have been set and no bias term is used => need to zero-mean the phenotype
      zeroMeanColumns( studydata->phenotypeData, studydata->cIndividuals, studydata->cPhenotypes, false );
      }

   if ( studydata->snpDataCovariance != nullptr )
      {  
      // perform the preprocessing on the covariance SNPs
      if ( simMissingness )
         {
         replaceSNPsByMissingness( studydata->snpDataCovariance, studydata->cIndividuals, studydata->cSnpsCovariance );
         }
      else
         {
         // covariance SNPs
         ImputeAndZeroMeanSNPs( studydata->snpInfoCovariance, studydata->snpDataCovariance, studydata->cIndividuals, studydata->cSnpsCovariance, unitVarianceSimilaritySNPs, zeroMeanSimilaritySNPs, normalizationInfo, usePseudocounts );   //Imputation for the kernel SNPs
         if ( kernelFile.empty() && (plinkFilenameRootCovarianceVector.size() > 1) )
            {
            // Zero mean each individual and divide by 1
            preprocessSNPsForBuildingCovariance( studydata->snpDataCovariance, studydata->cIndividuals, studydata->cSnpsCovariance, 1 );
            }
         else
            {
            // Zero mean each individual and divide by sqrt(N_SNPs)
            preprocessSNPsForBuildingCovariance( studydata->snpDataCovariance, studydata->cIndividuals, studydata->cSnpsCovariance, studydata->cSnpsCovariance );
            }
         }
      }

   if ( studydata->snpData != nullptr )
      {  
      // perform preprocessing on the test snps
      //  test SNPs (optionally not zero-meaned and unit-variance)
      ImputeAndZeroMeanSNPs( studydata->snpInfo, studydata->snpData, (BLAS_INT)studydata->cIndividuals, (BLAS_INT)studydata->cSnps, unitVarianceTestSNPs, zeroMeanTestSNPs, normalizationInfo, usePseudocounts );      //Imputation for the test SNPs
      }

   if ( covariatesData != nullptr )
      {
      // covariates have been set
      if ( useBiasTerm )
         {
         studydata->covariatesData = new real[ (covariatesData->n_covariates+1) * studydata->cIndividuals ];
         studydata->cCovariates = covariatesData->n_covariates + 1;
         covariatesData->CreateColumnMajorArrayFromVectorsAndMapping( studydata->individualLabels, studydata->covariatesData, covariatesData->n_covariates );
         zeroMeanColumns( studydata->covariatesData, studydata->cIndividuals, studydata->cCovariates-1, true );

         // append an array of ones -> not needed, if all variables are zero-mean
         for ( size_t i=(studydata->cCovariates-1)*studydata->cIndividuals; i<studydata->cCovariates*studydata->cIndividuals; ++i )
            {
            studydata->covariatesData[i] = 1.0;
            }
         }
      else
         {
         // no bias term (no bias term is added to the covariates)
         studydata->covariatesData = new real[ (covariatesData->n_covariates) * studydata->cIndividuals ];
         studydata->cCovariates = covariatesData->n_covariates;
         covariatesData->CreateColumnMajorArrayFromVectorsAndMapping( studydata->individualLabels, studydata->covariatesData, covariatesData->n_covariates );
         // Imputation for the kernel SNPs
         zeroMeanColumns( studydata->covariatesData, studydata->cIndividuals, studydata->cCovariates, true );
         }
      }
   else
      {
      if ( useBiasTerm )
         {
         // generate the covariates as an array of ones and make them the only covariates
         studydata->covariatesData = new real[ studydata->cIndividuals ];
         for ( size_t i=0; i<studydata->cIndividuals; ++i )
            {
            studydata->covariatesData[ i ] = 1.0;
            }
         studydata->cCovariates = 1;
         }
      else
         {
         // no bias term, and no covariates
         studydata->covariatesData = new real[ 0 ];
         studydata->cCovariates = 0;
         }
      }
   }

void permuteStudyData( EigenSym *eigen, StudyData *data )
   {
   //generate an index for the permutation
   data->GeneratePermuations( 1 );
   
   if ( data->cIndividuals != eigen->n )
      {
      Fatal( "permuteData : data->cIndividuals does not match eigen->n" );
      }

   if ( permutationMode == PERMUTE_KERNEL )     //Permute only the kernel
      {
      ProgressNL( "permuting the kernel" );
      permuteKernel( eigen, data );
      }
   else if ( permutationMode == PERMUTE_PHEN )  //Have to permute everything except for the SNP
      {
      ProgressNL( "permuting the phenotype, covariates, and the kernel" );
      permuteKernel( eigen, data );
      permutePhen( data );
      permuteCovariates( data );
      }
   }

void permuteKernel( EigenSym *eigen, StudyData *data )
   {
   real *permutedU = new real[data->cIndividuals*eigen->rankU];
   for ( size_t dim = 0; dim<(size_t)(eigen->rankU); ++dim )
      {
      // eigen->rankU must be positive!
      for ( size_t ind=0; ind<data->cIndividuals; ++ind )
         {
         permutedU[dim * data->cIndividuals + ind] = eigen->U[dim * data->cIndividuals + data->permutation[ind]];
         }
      }
   delete[] eigen->U;
   eigen->U = permutedU;
   }

void permutePhen( StudyData *data )
   {
   real *permutedPhen = new real[data->cIndividuals * data->cPhenotypes];
   for ( size_t dim = 0; dim<data->cPhenotypes; ++dim )
      {
      for ( size_t ind=0; ind<data->cIndividuals; ++ind )
         {
         permutedPhen[(dim * data->cIndividuals) + ind] = data->phenotypeData[(dim * data->cIndividuals) + data->permutation[ind]];
         }
      }
   delete[] data->phenotypeData;
   data->phenotypeData = permutedPhen;
   }

void permuteCovariates( StudyData *data )
   {
   real *permutedCovariates = new real[data->cIndividuals * data->cCovariates];
   for ( size_t dim = 0; dim<data->cCovariates; ++dim )
      {
      for ( size_t ind=0; ind<data->cIndividuals; ++ind )
         {
         permutedCovariates[(dim * data->cIndividuals) + ind] = data->covariatesData[(dim * data->cIndividuals) + data->permutation[ind]];
         }
      }
   delete[] data->covariatesData;
   data->covariatesData = permutedCovariates;
   }

/*
 *
 */
void setSNPPhenotypeCovariateFromInds( StudyData *studyData, const CPlinkFile &plinkFileData, const vector< size_t > &rgIdxIndividuals, CCovariatesData *covariatesData )
   {
   StudyData *psd = studyData;

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

   InitCovariates( psd, covariatesData );
   }

/*
 *
 */
void resetSNPData2( StudyData *studyData, size_t nSelectedSNPs )
   {
   real *newsnpdata = new real[ studyData->cIndividuals * nSelectedSNPs ];
   for ( size_t isnp = 0; isnp < nSelectedSNPs; ++isnp )
      {
      for ( size_t ind = 0; ind < studyData->cIndividuals; ++ind )
         {
         size_t idxDest = ColumnMajorIndex( studyData->cIndividuals, nSelectedSNPs, ind, isnp );
         size_t idxSrc = ColumnMajorIndex( studyData->cIndividuals, studyData->cSnps, ind, isnp );
         newsnpdata[ idxDest ] = studyData->snpData[ idxSrc ];
         }
      }
   delete[] studyData->snpData;
   studyData->snpData = newsnpdata;
   studyData->cSnps = nSelectedSNPs;
   }

/*
 *
 */
void copySNPDataToSNPcov( StudyData *study )
   {
   size_t limit = study->cIndividuals * study->cSnps;
   study->snpDataCovariance = new real[ limit ];
   for ( size_t i = 0; i < limit; ++i )
      {
      study->snpDataCovariance[ i ] = study->snpData[ i ];
      }
   }

/*
 *
 */
void printVector( real *vec, size_t n, const string& name )
   {
   Progress( "\n%s", name.c_str() );
   for ( size_t i = 0; i < n; ++i )
      {
      Progress( "\n  %lf", vec[ i ] );
      }
   Progress( "\n" );
   }

/*
 *
 */
void printVectorF( real *vec, size_t n, const string& name ) 
   {
   //
   // REVIEW:  no error handling
   //
   FILE *fp = fopen(name.c_str(), "w");
   for ( size_t i = 0; i < n; ++i )
      {
      fprintf( fp, "%f\n", vec[ i ] );
      }
   fclose( fp );
   }

/*
 *
 */
void printMatrixF( real *mat, size_t cRows, size_t cColumns, const string& name )
   {
   //
   // REVIEW:  no error handling
   //
   FILE *fp = fopen( name.c_str(), "w" );
   for ( size_t iRow = 0; iRow < cRows; ++iRow )
      {
      for ( size_t iColumn =0; iColumn < cColumns; ++iColumn )
         {
         size_t idx = ColumnMajorIndex( cRows, cColumns, iRow, iColumn );
         fprintf( fp, "%f ", mat[ idx ] );
         }
      fprintf( fp, "\n" );
      }
   fclose( fp );
   }

/*
 *
 */
void printVectorF2( size_t *mat, size_t cRows, const string& name )
   {
   //
   // REVIEW:  no error handling
   //
   FILE *fp = fopen( name.c_str(), "w" );
   for ( size_t idx = 0; idx < cRows; ++idx )
      {
         fprintf( fp, "%d ", (int)mat[ idx ] );
      }
   fclose( fp );
   }

/*
 *
 */
void printMatrix( real *mat, size_t m, size_t n, const string& name )
   {
   Progress( "\n%s", name.c_str() );
   for ( size_t i = 0; i < m; ++i )
      {
      for ( size_t j = 0; j < n; ++j )
         {
         size_t idx = ColumnMajorIndex( m, n, i, j );
         Progress( "%lf ", mat[ idx ] );
         }
      Progress( "\n" );
      }
   Progress( "\n" );
   }

/*
 *
 */
real meanSqError( real *mu, real *y, size_t n )
   {
   real err = 0.0;
   for ( size_t i = 0; i < n; ++i )
      {
      real er = mu[ i ] - y[ i ];
      err = err + ( er * er );
      }
   return( err );
   }

int run_lowrank_null_fitting( EigenSym *eigen, ResultsLMM *results, StudyData *data )
   {
   ProgressNL( " --          Start lowrank training:" );

   FastLMM_lowrank LMM_lowrank = FastLMM_lowrank( eigen );

   LMM_lowrank.set_Phenotype( data->phenotypeData, data->phenotypeLabels[ data->ind_phenotype ] );
   LMM_lowrank.set_SNPs( data->snpData, data->cSnps, data->snpInfo );
   LMM_lowrank.set_covariates( data->covariatesData, data->cCovariates );
   LMM_lowrank.set_SNPsCov( data->snpDataCovariance, data->cSnpsCovariance, data->snpInfoCovariance );
   LMM_lowrank.init();
   LMM_lowrank.trainNull();
   LMM_lowrank.updateExcludeList();

   //get the results
   LMM_lowrank.getResults( *results );
   ProgressNL( " --           Lowrank training done:" );
   return( 0 );
   }

void create_eye( real *sigma, size_t m )
   {
   for ( size_t i = 0; i < m; ++i )
      {
      size_t idx = ColumnMajorIndex( m, m, i, i );
      sigma[ idx ] = 1.0;
      }
   }

/*
 *  Compute the pvalue using Likelihood-ratio test as long as
 *    deltaLL is positive.  Otherwise return 1.0
 *    Issue a warning if deltaLL is too far from an expected
 *    value 
 */
real PvalFromLikelihoodRatioTest( real deltaLL, real dof_ )
   {
   real retVal = real1;
   if ( deltaLL > real0 )
      {
      retVal = Gamma::gammaQ( deltaLL, dof_ );
      }
   else
      {
      // Optimizer didn't find the 'best' solution for LL_0
      real absDeltaLL = ::abs( deltaLL );
      if ( absDeltaLL > maxDeltaLL )
         {
         Warn( "deltaLL[%.2e] > maxDeltaLL[%2e].  Setting pval to 1.0"
               "\n      This may indicate optimizer parameters need tuning in the sources",
               absDeltaLL, maxDeltaLL );
         }
      }

   return( retVal );
   }

/*
 *  When running SnpPairs, the number of tasks we actually run must
 *    follow the form (n * (n-1)) as we are processing the diagonal
 *    of a matrix.  
 *  The user gives an approximate task number T and we compute the 
 *    number of tasks T' less than T that meets the above form.
 *    When the user provides a task number above T' and less than
 *    T, we emit an 'empty' file that follows the format so post
 *    processing will not wait for the files.
 */
void WriteEmptyTaskResultFile( const string& ofile, const SnpFilterOptions& sfo )
   {
   // This task is beyond the data but we need to produce the file
   //   so post processing will 'work'
   // Copied from FastLMM_woodbury output to produce an empty file
   FILE *pFile = fopen( ofile.c_str(), "w" );
   if ( pFile == nullptr )
      {
      Fatal( "unable to open output file.  Sorry."
            "\n  CRT Error: %s", strerror( errno ) );
      }

   // NOTE:  GCC and VC don't share the same format specifiers for 
   //        size_t, so we convert to use int in this case.
   // Print the header and setup to measure our throughput
   fprintf( pFile,   "#  Tasks Requested: %9d", (int)sfo.cTasksAsk );
   fprintf( pFile, "\n#     Tasks Needed: %9d", (int)sfo.cTasksMax );
   fprintf( pFile, "\n#      Task Number: %9d", (int)sfo.iTaskAsk );
   fprintf( pFile, "\n#       Partitions: %9d", (int)sfo.cBlocksMax );
   fprintf( pFile, "\nsnpId_1\tsnpId_2\tpvalue\tSnpPair_Wgt\tSnpPair_OddsRatio\n" );

   // Report the end numbers
   fprintf( pFile,"** End **\n%d,%d\n", 0, 0 );
   fclose( pFile );
   }

/*
 * ValidateIteringFastLmmPreconditions() asserts
 *   all the preconditions for running are met.
 */
void ValidateIteratingFastLmmPreconditions()
   {
   if ( useWoodbury ) Fatal( "useWoodbury is not supported in this scenario." );
   if ( selectTopK )  Fatal( "Snp filtering is not supported in this scenario." );
   if ( !alternatePhenotypeParameters.fUseAlternatePhenotype ) Fatal( "Alternate Phenotype file is required for this scenario." );
   if ( plinkFileType != Binary ) Fatal( "-Binary PLINK file (-bfile) is required for this scenario." );
   if ( simMissingness != false ) Fatal( "simMissingness is required to be false for this scenario." );
   if ( permutationMode != NO_PERM ) Fatal( "-permute options are not supported in this scenario." );
   }

/*
 *  Load the snp, individual genotypes, and covariates data from the plink file(s)
 */
void LoadIteratingFastLmmEarlyData( CCovariatesData& covariatesData, CPlinkAlternatePhenotypeFile& apFile, CPlinkFamFile& famFile, CPlinkBimFile& bimFile )
   {
   if ( FUseCovariatesFile() )
      {
      covariatesData.ReadCovariatesFile( covariatesFile );
      }

   // read in the AlternatePhenotype Information
   apFile.Read( alternatePhenotypeParameters.alternatePhenotypeFilename );
   if ( alternatePhenotypeParameters.selectedPhenotypeId.empty() )
      {
      apFile.SetSelectedPhenotype( alternatePhenotypeParameters.selectedPhenotypeIndex );
      }
   else
      {
      apFile.SetSelectedPhenotype( alternatePhenotypeParameters.selectedPhenotypeId );
      }

   // Read in the FAM and BIM files as the individuals part of the Plink binary file load (-bfile)
   famFile.Load( plinkFilenameRoot + ".fam" );
   bimFile.Load( plinkFilenameRoot + ".bim" );
   }

/*
 * GetListOfIndividualKeys() intersects Individuals across PLink, 
 *   AlternatePhenotype, and Covariates data.  The order of processing
 *   the individuals is not important so we use AlternatePhenotype as 
 *   the base individual list as it is frequently smaller.
 */
vector< string > GetListOfIndividualKeys( CCovariatesData& covariatesData, CPlinkAlternatePhenotypeFile& apFile, CPlinkFamFile& famFile )
   {
   vector< string > listOfIndividualKeys;
   size_t cIndividualsMax = apFile.CountOfIndividuals();    // number of individuals read from file
   for ( size_t iIndividual=0; iIndividual<cIndividualsMax; ++iIndividual )
      {
      // Get the key from apFile and verify it has a valid phenotype
      string keyForIndividual = apFile.KeyFromIndividualIndex( iIndividual );
      real pheno = apFile.PhenotypeValue( iIndividual );
      if ( _isnan( pheno ) )
         {
         Warn( "Found 'missing' phenotype marker.  Excluding [%s] from the study", keyForIndividual.c_str() );
         continue;
         }

      // verify the individual exists in the fam file
      IndividualInfo* pIi = famFile.PIndividualInfo( keyForIndividual );
      if ( pIi == nullptr )
         {
         Warn( "Individual from alternate phenotype file not found in .FAM file."
             "\nExcluding [%s] from the study", keyForIndividual.c_str() );
         continue;
         }

      if ( FUseCovariatesFile() && !covariatesData.FKeyIndividualInDatFile( keyForIndividual ) )
         {
         Warn( "Individual from alternate phenotype file not found in Covariates file [%s]"
             "\nExcluding [%s] from the study", covariatesData.filename.c_str(), keyForIndividual.c_str() );
         continue;
         }
      listOfIndividualKeys.push_back( keyForIndividual );
      }
   return( listOfIndividualKeys );
   }

/*
 * FGenotypeVectorHasVariation() returns true if any of the 
 *   non-NaN elements are different 
 */
bool FGenotypeVectorHasVariation( const vector< real >& vGenotypes )
   {
   real previousGenotype = -1.0;
   size_t cGenotypeMax   = vGenotypes.size();
   real const *pGenotype = &vGenotypes[0];

   for ( size_t iGenotype = 0; iGenotype < cGenotypeMax; ++iGenotype )
      {
      real genotype = *pGenotype;
      if ( !_isnan( genotype ) )
         {
         if ( previousGenotype == -1.0 )
            {
            previousGenotype = genotype;
            }
         else
            {
            if ( previousGenotype != genotype )
               {
               return( true );
               }
            }
         }
      ++pGenotype;;
      }
   return( false );
   }

struct IteratingFastLmmParameters
   {
   CCovariatesData covariatesData;
   CPlinkAlternatePhenotypeFile apFile;
   CPlinkFamFile  famFile;
   CPlinkBimFile  bimFile;
   CBedFile       bedFile;
   };

// ********** ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
/*
 *  DoIteratingFastLmm() which is used to reduce the memory footprint
 *    by reading a subset of the SNPs under test in a time.
 */
void DoIteratingFastLmm()
   {
   ValidateIteratingFastLmmPreconditions();

   // Get working variables and load the data from the disk
   CCovariatesData covariatesData;
   CPlinkAlternatePhenotypeFile apFile;
   CPlinkFamFile  famFile;
   CPlinkBimFile  bimFile;

   CTimer iteratingFastLmmTimer( true );
   ProgressNL( "               DoIteratingFastLmm():" );

   // Load all the covariates data, alternate phenotype data, 
   //   individual's data, and the SnpInfo and intersect to get
   //   the individuals we will be working with this time
   LoadIteratingFastLmmEarlyData( covariatesData, apFile, famFile, bimFile );
   vector< string > individualLabels = GetListOfIndividualKeys( covariatesData, apFile, famFile );

   // Extract the 'size of the problem' parameters from all 
   //   the input files and finish preparing the data for
   //   processing
   size_t cPhenotypesMax  = apFile.CountOfPhenotypes();
   size_t cPhenotypes  = 1;                                 // we currently only support 1 phenotype at a time.
   size_t cBiasTerms = useBiasTerm ? 1 : 0;
   size_t cCovariates = covariatesData.n_covariates;
   size_t cIndividualsMax = apFile.CountOfIndividuals();    // number of individuals read from file
   size_t cIndividuals = individualLabels.size();
   size_t cSnpsMax = bimFile.CountOfSnpInfo();
   size_t cSnpsPerBlockMax = min( cSnpsMax, cIteratingBlockSizeMax );        // use default block size for iterating if less than cSnpsMax
   size_t cSnpsPerCovarianceBlockMax = min( cSnpsMax, cSnpsPerBlockMax * 2 );    // double the size of the 

   ProgressNL( "                  Generating Kernel:" );
   ProgressNL( "         Number of Individuals Read: %10d", cIndividualsMax );
   ProgressNL( "     Number of Individuals Selected: %10d", cIndividuals );
   ProgressNL( "          Number of Phenotypes Read: %10d", cPhenotypesMax );
   ProgressNL( "      Number of Phenotypes Selected: %10d", cPhenotypes );
   ProgressNL( "          Number of Covariates Read: %10d", covariatesData.n_covariates );
   ProgressNL( "                    Using bias term: %10s", useBiasTerm ? "True" : "False" );
   ProgressNL( "                Number of SNPs Read: %10d", cSnpsMax );
   ProgressNL( "          SNP block processing size: %10d", cSnpsPerBlockMax );

   // Now open .bed file and read block of data to create the kernel 
   //   using blocks to read, and then iterate over all the read 
   //   blocks of the SNP genotype data in from the .BED file to 
   //   create the kernel
   CBedFile bedFile;
   bedFile.Open( plinkFilenameRoot + ".bed", &famFile, &bimFile );
   if ( bedFile.GetLayoutMode() != LayoutGroupGenotypesBySnp )
      {
      Fatal( "This scenario only supports .bed files organized by SNP.  Found [%s] organized by Individual.", bedFile.Filename().c_str() );
      }

   vector< size_t > listOfIndividualIndexes;
   for ( size_t i=0; i<individualLabels.size(); ++i )
      {
      listOfIndividualIndexes.push_back( famFile.IndexFromKeyIndividual( individualLabels[ i ] ) );
      }

   // Create our 'working' variables
   vector< size_t >  rgIdxSnpsSelected;            rgIdxSnpsSelected.reserve( cSnpsMax );
   vector< SnpInfo > rgSnpInfoSelected;            rgSnpInfoSelected.reserve( cSnpsMax );
   vector< SnpInfo > rgCovarianceSnpInfoInBlock;   rgCovarianceSnpInfoInBlock.reserve( cSnpsPerCovarianceBlockMax );
   vector< real >    rgCovarianceSnpValues( cSnpsPerCovarianceBlockMax * cIndividuals );   // a matrix which is SNPs x IndividualGenotype at that SNP
   vector< real >    rgTestSnpValues( cSnpsPerBlockMax * cIndividuals );  // a matrix which is SNPs x IndividualGenotype at that SNP
   vector< real >    rgKernel( cIndividuals * cIndividuals );  // create the kernel matrix!
   vector< real >    vecGenotypes( cIndividuals );

// ********** ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
   ProgressNL( "     Processing SNPs for covariance:\n" );

   size_t iSnpToRead   = 0;
   size_t cCovarianceSnpsInBlock = 0;
   size_t iBlock = 0;
   do
      {
      // read in this 'block' of the .bed file
      cCovarianceSnpsInBlock = 0;
      for ( ; (cCovarianceSnpsInBlock < cSnpsPerBlockMax) && (iSnpToRead < cSnpsMax); ++iSnpToRead )
         {
         bedFile.ReadGenotypes( iSnpToRead, listOfIndividualIndexes, vecGenotypes );
         if ( FGenotypeVectorHasVariation( vecGenotypes ) )
            {
            memcpy( &rgCovarianceSnpValues[ cCovarianceSnpsInBlock * cIndividuals ], &vecGenotypes[ 0 ], cIndividuals * sizeof( real ) );
            ++cCovarianceSnpsInBlock;
            rgIdxSnpsSelected.push_back( iSnpToRead );
            rgCovarianceSnpInfoInBlock.push_back( *bimFile.PSnpInfo( iSnpToRead ) );
            }
         else
            {
            Warn( "SNP[%s] has no variation.  Filtering SNP from dataset.\n", bimFile.PSnpInfo( iSnpToRead )->idSnp.c_str() );
            }
         }

      // these SNPs are used to build the covariance
      ImputeAndZeroMeanSNPs( rgCovarianceSnpInfoInBlock, &rgCovarianceSnpValues[0], cIndividuals, cCovarianceSnpsInBlock, unitVarianceSimilaritySNPs, zeroMeanSimilaritySNPs, normalizationInfo, usePseudocounts );
      preprocessSNPsForBuildingCovariance( &rgCovarianceSnpValues[0], cIndividuals, cCovarianceSnpsInBlock, 1 );

      // process to create this 'chunk' of the kernel and add to the previous chunk
      BLAS_INT nRows_int = (BLAS_INT)cIndividuals;
      BLAS_INT nCols_int = (BLAS_INT)cCovarianceSnpsInBlock;

      // symmetric rank k update. result will be in upper or lower storage format
      BlasLapack::syrk( &char_L,
                        &char_N,
                        &nRows_int,
                        &nCols_int,
                        &real_one,
                        &rgCovarianceSnpValues[ 0 ],
                        &nRows_int,
                        &real_one,        // rgKernel already initialized to 0.0
                        &rgKernel[ 0 ],
                        &nRows_int
                        );

      move( rgCovarianceSnpInfoInBlock.begin(), rgCovarianceSnpInfoInBlock.end(), back_inserter( rgSnpInfoSelected ) );
      rgCovarianceSnpInfoInBlock.clear();
      ++iBlock;
      if ( (iBlock == 1) || ((iBlock % 100) == 0) )
         {
         Progress( "\r                      Complete pass: %6d", ++iBlock );
         }
      } while ( iSnpToRead < cSnpsMax );

   Progress( "\r                      Complete pass: %6d\n", iBlock );

   size_t cSnps = rgIdxSnpsSelected.size();
   size_t cKernelElements = cIndividuals * cIndividuals;
   real rT = (real)cSnps;

   for ( size_t i = 0; i < cKernelElements; ++i )
      {
      rgKernel[ i ] /= rT;
      }

   // Prepare the TEST DATA information for FastLmm
   string eigenId = alternatePhenotypeParameters.alternatePhenotypeFilename + "_" + plinkFilenameRoot;
   EigenSym eigen;

   //TODO: DEBUG
   //printMatrixF( &rgKernel[ 0 ], cIndividuals, cIndividuals, "Out\\Refactored_Kernel_Matrix_2000" );

   eigen.set_fromKernel( eigenId, &rgKernel[ 0 ], cIndividuals, cIndividuals, individualLabels );

// ********** ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
   // Done creating the Kernel!
   //  Start on setting up the data
   ProgressNL( "           Processing SNPs TestData:" );
   ProgressNL( "                  Output results to: [%s]", fullOutputFilename.c_str() );

   // Create the phenotype data array (a vector has the same 
   //   in memory layout as an array[]) and perform Eigenstrat 
   //   Imputation on the phenotype data (false => DO NOT Mean Imput NaNs)
   vector< string > phenotypeLabels;
   phenotypeLabels.push_back( apFile.SelectedPhenotypeName() );
   vector< real > rgPhenotypeData;
   rgPhenotypeData.reserve( cIndividuals );
   for ( auto it = individualLabels.begin(); it != individualLabels.end(); ++it )
      {
      rgPhenotypeData.push_back( apFile.PhenotypeValue( *it ) );
      }

// ********** ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
   // prepare the covariates and bias data and create the 
   //   associated array for later processing
   covariatesData.ImputeAndZeroMeanCovariatesFix( individualLabels );
   vector< real > rgCovariateData( (cCovariates + cBiasTerms) * cIndividuals );
   for ( size_t iCovariate=0; iCovariate<cCovariates; ++iCovariate )
      {
      real *pRealOut = &rgCovariateData[ iCovariate * cIndividuals ];   // we get cIndividuals each pass
      covariatesData.CovariateValues( iCovariate, individualLabels, pRealOut );
      }
   for ( size_t iBiasTerms = 0; iBiasTerms<cBiasTerms;  ++iBiasTerms )
      {
      real *pRealOut = &rgCovariateData[ (cCovariates + iBiasTerms) * cIndividuals ];
      for ( size_t iIndividual=0; iIndividual<cIndividuals; ++iIndividual )
         {
         pRealOut[ iIndividual ] = 1.0;
         }
      }

// ********** ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
   // Prepare for the LMM work
   FastLMM_fullrank LMM_fullrank( &eigen );
   LMM_fullrank.set_Phenotype( &rgPhenotypeData[ 0 ], phenotypeLabels[ 0 ] );

   vector< SnpInfo > rgTestSnpInfoInBlock;   rgTestSnpInfoInBlock.reserve( cSnpsPerBlockMax );
   vector< size_t > rgExcludeStart( cSnps );
   vector< size_t > rgExcludeCount( cSnps );
   UpdateExcludeLists( excludeSnpParams.distType, rgSnpInfoSelected, rgSnpInfoSelected, &rgExcludeStart[0], &rgExcludeCount[0] );

   // We are now ready to iterate reading the test SNPS and covariance SNPs
   //   to produce the GWAS output.  We KNOW that both test and covariance
   //   SNPs come from the same file so we read them into the same buffer and
   //   then break them out for their individual processing.
   //   We open the output file and write the header
   //   then start re-reading the plink file to get orginal 
   //   data and finally close the file after all is read
   ResultsLMM results;
   results.fIterativeResults = true;         // signal we are iterating to produce results
   results.OpenResultsFile( fullOutputFilename, fVerboseOutput, pvaluePrintThreshold, cCovariates + cBiasTerms );

   bool fExcludeSnps = ( excludeSnpParams.fileType == FileSim ) ? true : false;

   size_t iCovarianceMin = 0;
   size_t iCovarianceMax = 0;
   size_t iTestMin = 0;
   size_t iTestMax = 0;

   do
      {
      // Set loop control variables
      iTestMin = iTestMax;
      iTestMax = min( iTestMax + cSnpsPerBlockMax, cSnps );

      if ( fExcludeSnps )
         {
         iCovarianceMin = rgExcludeStart[ iTestMin ];
         iCovarianceMax = rgExcludeStart[ iTestMax-1 ]  + rgExcludeCount[ iTestMax-1 ];
         }
      else
         {
         iCovarianceMin = iTestMin;
         iCovarianceMax = iTestMax;
         }

      cCovarianceSnpsInBlock = 0;
      size_t cTestSnpsInBlock = 0;

      // Fetch preceding Covariance data
      size_t cSnpsInPreFetch = iTestMin - iCovarianceMin;
      for ( size_t iSnp=0; iSnp<cSnpsInPreFetch; ++iSnp )
         {
         iSnpToRead = rgIdxSnpsSelected[ iSnp + iCovarianceMin ];
         bedFile.ReadGenotypes( iSnpToRead, listOfIndividualIndexes, vecGenotypes );
         memcpy( &rgCovarianceSnpValues[ cCovarianceSnpsInBlock * cIndividuals ], &vecGenotypes[ 0 ], cIndividuals * sizeof( real ) );
         rgCovarianceSnpInfoInBlock.push_back( *bimFile.PSnpInfo( iSnpToRead ) );
         ++cCovarianceSnpsInBlock;
         }

      // Fetch the Test data and put in both Covariance vector and Test vector
      size_t cSnpsInMainBlock = iTestMax - iTestMin;
      for ( size_t iSnp=0; iSnp<cSnpsInMainBlock; ++iSnp )
         {
         iSnpToRead = rgIdxSnpsSelected[ iSnp + iTestMin ];
         bedFile.ReadGenotypes( iSnpToRead, listOfIndividualIndexes, vecGenotypes );
         memcpy( &rgCovarianceSnpValues[ cCovarianceSnpsInBlock * cIndividuals ], &vecGenotypes[ 0 ], cIndividuals * sizeof( real ) );
         rgCovarianceSnpInfoInBlock.push_back( *bimFile.PSnpInfo( iSnpToRead ) );
         ++cCovarianceSnpsInBlock;

         memcpy( &rgTestSnpValues[ cTestSnpsInBlock * cIndividuals ], &vecGenotypes[ 0 ], cIndividuals * sizeof( real ) );
         rgTestSnpInfoInBlock.push_back( *bimFile.PSnpInfo( iSnpToRead ) );
         ComputeSnpInfoStatistics( rgTestSnpInfoInBlock[ cTestSnpsInBlock ], cIndividuals, &rgTestSnpValues[ cTestSnpsInBlock * cIndividuals ] );
         ++cTestSnpsInBlock;
         }

      // Fetch the suffix Covariance data
      size_t cSnpsInPostFetch = iCovarianceMax - iTestMax;
      for ( size_t iSnp=0; iSnp<cSnpsInPostFetch; ++iSnp )
         {
         iSnpToRead = rgIdxSnpsSelected[ iSnp + iTestMax ];
         bedFile.ReadGenotypes( iSnpToRead, listOfIndividualIndexes, vecGenotypes );
         memcpy( &rgCovarianceSnpValues[ cCovarianceSnpsInBlock * cIndividuals ], &vecGenotypes[ 0 ], cIndividuals * sizeof( real ) );
         rgCovarianceSnpInfoInBlock.push_back( *bimFile.PSnpInfo( iSnpToRead ) );
         ++cCovarianceSnpsInBlock;
         }

      // Test SNPs (optionally not zero mean and not unit variance)
      ImputeAndZeroMeanSNPs( rgTestSnpInfoInBlock, &rgTestSnpValues[0], cIndividuals, cTestSnpsInBlock, unitVarianceTestSNPs, zeroMeanTestSNPs, normalizationInfo, usePseudocounts );

      // covariance SNPs (always zero mean and unit variance)
      ImputeAndZeroMeanSNPs( rgCovarianceSnpInfoInBlock, &rgCovarianceSnpValues[0], cIndividuals, cCovarianceSnpsInBlock, unitVarianceSimilaritySNPs, zeroMeanSimilaritySNPs, normalizationInfo, usePseudocounts );
      preprocessSNPsForBuildingCovariance( &rgCovarianceSnpValues[0], cIndividuals, cCovarianceSnpsInBlock, cSnps ); // divide by total number of covariance SNPs

      LMM_fullrank.set_SNPs( &rgTestSnpValues[0], cTestSnpsInBlock, rgTestSnpInfoInBlock );
      LMM_fullrank.set_SNPsCov( &rgCovarianceSnpValues[0], cCovarianceSnpsInBlock, rgCovarianceSnpInfoInBlock );
//      if ( excludeSnpParams.fileType == FileSim )
//         {
//         LMM_fullrank.snpInfoCovariance = rgSnpInfoSelected;
//         LMM_fullrank.set_UW( &rgCovarianceValues[0] );
//         LMM_fullrank.nSNPsW = rgCovarianceValues.size();
//         // we are not creating the SNP genotype values that the normal path does.  
//         //   This could be a problem...
//         }
//      LMM_fullrank.set_SNPsCov( data->snpDataCovariance, cSnpsInBlock, rgSnpInfoInBlock );
//        LMM_fullrank.UW = nullptr;
//        LMM_fullrank.nSNPsW = 0;

      if ( results.iResultBlock == 0 )
         {
         // ONLY run set_covariates, init, and trainNull on the first pass
         LMM_fullrank.set_covariates( &rgCovariateData[0], cCovariates + cBiasTerms );
         LMM_fullrank.init();
         LMM_fullrank.trainNull();
         }

      LMM_fullrank.updateExcludeList();
      if ( runAlternativeModels )
         {
         LMM_fullrank.trainAlt();
         }

      LMM_fullrank.getResults( results );
      results.WriteResultsBody();
      ++results.iResultBlock;
      results.cPriorResults += cTestSnpsInBlock;
      rgTestSnpInfoInBlock.clear();
      rgCovarianceSnpInfoInBlock.clear();
      } while ( iTestMax < cSnps );

   results.WriteResultsSummary();
   results.CloseResultsFile();

   iteratingFastLmmTimer.Stop();
   iteratingFastLmmTimer.Report( 0x02 | -1, "  DoIteratingFastLmm() elapsed time: %s" );
   }

/*
 * SNPs can be excluded from the kernel based on the relative 
 *   positions in the genome to the SNP under test.  
 *   UpdateExcludeLists() fills a pair of arrays associated with 
 *   the SNPs selected to produce the kernel (Covariance) and 
 *   the SNPs under test.  The array pair indicates which SNPs
 *   in the kernel need to be removed to process the specific SNP.
 *   rgExcludeStart holds the index into an array specifying the 
 *   SNPs that make up the kernel and the Covariance array
 */
void UpdateExcludeLists( ExcludeSNPfromSimType exclusionType, 
                         vector< SnpInfo >& rgSnpInfoTest, 
                         vector< SnpInfo >& rgSnpInfoCovariance, 
                         size_t *rgExcludeStart, 
                         size_t *rgNExclude )
   {
   // PRE-CONDITION:
   //   The arrays rgNExclude and rgExcludeStart _must_ be at 
   //   least rgSnpInfoTest.size() or larger

   size_t cSnpsTest = rgSnpInfoTest.size();
   size_t cSnpsCovariance = rgSnpInfoCovariance.size();

   switch( exclusionType )
      {
   case NoDist:
      {
      size_t cb = rgSnpInfoTest.size() * sizeof( size_t );
      memset( rgNExclude, 0, cb );
      memset( rgExcludeStart, 0, cb );
      }
      break;
   case ByID:
      Fatal( "Exlusion by ID is not currnetly supported." );
      break;
   case ByPosition:
      {
      size_t iExclusionStart = 0;
      size_t cExclusions = 0;
      for( size_t iSnp=0; iSnp<cSnpsTest; ++iSnp )
         {
         // Set the lower limit of the exclusion region (a SNP should always excludes at least itself)
         while( 1 )
            {
            // past the end of the covariance data and still below this SNP
            if ( iExclusionStart >= cSnpsCovariance )
               {
               rgExcludeStart[ iSnp ] = cSnpsCovariance;
               break;
               }
            int distance = CPlinkFile::getBpDistance( rgSnpInfoCovariance[ iExclusionStart ], rgSnpInfoTest[ iSnp ] );
            if ( distance >= (0 - excludeSnpParams.bpDist) )
               {
               // Covariance location data is now at or above the lower limit of the
               //   exlusion range.  Set start location 
               rgExcludeStart[ iSnp ] = iExclusionStart;
               break;
               }

            // bump the index to exclude and reduce the number to exclude
            ++iExclusionStart;
            if ( cExclusions )
               {
               // reduce the number of Covariance SNS to exclude
               --cExclusions;
               }
            }

         // Set the upper limit of the exclusion region (a SNP always excludes at least itself)
         while( 1 )
            {
            if ( iExclusionStart + cExclusions >= cSnpsCovariance )
               {
               rgNExclude[ iSnp ] = cSnpsCovariance - iExclusionStart;
               break;
               }

            int distance = CPlinkFile::getBpDistance( rgSnpInfoCovariance[ iExclusionStart + cExclusions ], rgSnpInfoTest[ iSnp ] );
            if ( distance > excludeSnpParams.bpDist )
               {
               // Covariance location data is above the upper limit of the
               //   exlusion range.  Set number to exclude
               rgNExclude[ iSnp ] = cExclusions;
               break;
               }
            ++cExclusions;
            }
         }
      }
      break;
   case ByGeneticDistance:
      {
      size_t iExclusionStart = 0;
      size_t cExclusions = 0;
      for( size_t iSnp=0; iSnp<cSnpsTest; ++iSnp )
         {
         // Set the lower limit of the exclusion region (a SNP should always excludes at least itself)
         while( 1 )
            {
            // past the end of the covariance data and still below this SNP
            if ( iExclusionStart >= cSnpsCovariance )
               {
               rgExcludeStart[ iSnp ] = cSnpsCovariance;
               break;
               }
            real distance = CPlinkFile::getGeneticDistance( rgSnpInfoCovariance[ iExclusionStart ], rgSnpInfoTest[ iSnp ] );
            if ( distance >= (0.0 - excludeSnpParams.geneticDist) )
               {
               // Covariance location data is now at or above the lower limit of the
               //   exlusion range.  Set start location
               rgExcludeStart[ iSnp ] = iExclusionStart;
               break;
               }

            // bump the index to exclude and reduce the number to exclude
            ++iExclusionStart;
            if ( cExclusions )
               {
               // reduce the number of Covariance SNS to exclude
               --cExclusions;
               }
            }

         // Set the upper limit of the exclusion region (a SNP always excludes at least itself)
         while( 1 )
            {
            if ( iExclusionStart + cExclusions >= cSnpsCovariance )
               {
               rgNExclude[ iSnp ] = cSnpsCovariance - iExclusionStart;
               break;
               }

            real distance = CPlinkFile::getGeneticDistance( rgSnpInfoCovariance[ iExclusionStart + cExclusions ], rgSnpInfoTest[ iSnp ] );
            if ( distance > excludeSnpParams.geneticDist )
               {
               // Covariance location data is above the upper limit of the
               //   exlusion range.  Set number to exclude
               rgNExclude[ iSnp ] = cExclusions;
               break;
               }
            ++cExclusions;
            }
         }
      }
      break;
   default:
      Fatal( "UpdateExcludeLists() was passed a bad exclusion type value [%d]", exclusionType );
      break;
      }
   }
