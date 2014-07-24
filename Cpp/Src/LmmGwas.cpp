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
 * LmmGwas - {Do the GWAS work}
 *
 *         File Name:   LmmGwas.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *    Module Purpose:   This file implements the GWAS functionality 
 *                      in FastLmmC
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Include Files
 */
#include "FastLmmC.h"

/*
 * Externs we need
 */

/*
 * Globals we define
 */
/*
 * The code
 */
void LmmGwas( StudyData *data, EigenSym *eigen, ResultsLMM *results )
   {
   if ( fDoSnpPairs 
     && ( (!(eigen->n > eigen->rankU)) || !useWoodbury )
      )
      {
      Fatal( "SnpPairs only implemented in runGWAS_woodbury() which requires eigen->n[%d] > eigen->rankU[%d].", eigen->n, eigen->rankU );
      }

   //run the Gwas given the data
   if ( eigen->n > eigen->rankU )
      {
      if ( useWoodbury )
         {
         runGWAS_woodbury( eigen, results, data );
         }
      else
         {
         runGWAS_lowrank( eigen, results, data );
         }
      }
   else
      {
      runGWAS_fullrank( eigen, results, data );
      }
   }

int computeEigenSym_Kernel( const string& id, StudyData *data, EigenSym *eigen )
   {
   // Computer the full rank EigenSym from the input Kernel matrix
   Verbose( verboseMajor, "  Compute full rank EigenSym Kernel:" );
   eigen->set_fromKernel( id , data->kernelData, data->cIndividuals, data->cIndividuals, data->individualLabels );
   return 0;
   }
      
int computeEigenSym_SnpSample( const string& id, StudyData *data, EigenSym *eigen )
   {
   if( (data->cSnpsCovariance >= data->cIndividuals) || forceFullRank )
      {
      //do the full rank model
      BLAS_INT rankU_int = (BLAS_INT)data->cIndividuals;
      BLAS_INT n_individuals_int=(BLAS_INT)data->cIndividuals;

      // Compute the covariance matrix as the inner product matrix of the Covariance SNPs
      real *Kernel = new real[ data->cIndividuals * data->cIndividuals ];
      ZeroRealArray( Kernel, data->cIndividuals * data->cIndividuals );
      BLAS_INT nCols_int = (BLAS_INT)data->cSnpsCovariance;

      BlasLapack::syrk( &char_L, &char_N, &n_individuals_int, &nCols_int, &real_one, data->snpDataCovariance, &n_individuals_int, &real_zero, Kernel, &n_individuals_int );

      //optionally add kernels from other SNP samples.
      real nSNPsKernelTotal = 1.0 * data->cSnpsCovariance;
      for (size_t i_file = 1 ; kernelFile.empty() && i_file<plinkFilenameRootCovarianceVector.size(); ++i_file)
         {
         // Read the id, genotype, and phenotype data for each individual in the study
         //   The covariance file does not participate in the SNP filtering (-extract)
         Verbose( verboseMajor, "                    Processing File: [%s]", plinkFilenameRootCovarianceVector[i_file].c_str() );
         CPlinkFile plinkFileCovariance;
         plinkFileCovariance.ReadPlinkFiles( plinkFilenameRootCovarianceVector[i_file], plinkFileTypeCovariance, alternatePhenotypeParameters );
         if( data->cIndividuals != plinkFileCovariance.cIndividuals )
            {
            Fatal("Number of individuals does not match up in current covariancefile %s! should be %i is %i", plinkFilenameRootCovarianceVector[i_file].c_str(), data->cIndividuals, plinkFileCovariance.cIndividuals );
            }
         if( simMissingness)
            {
            Fatal( "-simMissingness is only supported for a single SNP-file for computing similarities" );
            replaceSNPsByMissingness( plinkFileCovariance.snpData, plinkFileCovariance.cIndividuals, plinkFileCovariance.cSnps );
            }
         else
            {
            ImputeAndZeroMeanSNPs( plinkFileCovariance.rgSnpInfo, plinkFileCovariance.snpData, plinkFileCovariance.cIndividuals, plinkFileCovariance.cSnps, unitVarianceSimilaritySNPs, zeroMeanSimilaritySNPs, normalizationInfo, usePseudocounts );     //Eigenstrat Imputation for the kernel SNPs
            preprocessSNPsForBuildingCovariance( plinkFileCovariance.snpData, plinkFileCovariance.cIndividuals,plinkFileCovariance.cSnps , 1 );                               //Zero mean each individual and divide by sqrt(N_SNPs)
            }

         BLAS_INT nsnpscurrentkernel_int = (BLAS_INT)plinkFileCovariance.cSnps;   //Warning: size_t BLAS_INT conversion
         BlasLapack::syrk( &char_L, &char_N, &n_individuals_int, &nsnpscurrentkernel_int, &real_one, plinkFileCovariance.snpData, &n_individuals_int, &real_one,  Kernel, &n_individuals_int );
         nSNPsKernelTotal += plinkFileCovariance.cSnps;
         }

      //real sum=0.0;
      if( !simMissingness )
         {
         if ( kernelFile.empty() && (plinkFilenameRootCovarianceVector.size()>1) )
            {
            for(size_t row=0; row<data->cIndividuals;++row)
               {
               for(size_t col=row; col<data->cIndividuals;++col)
                  {
                  Kernel[row*data->cIndividuals + col] /= nSNPsKernelTotal;
                  }
               }
            }
         }

      // compute the Eigenvectors and set EigenSym from the kernel
      eigen->set_fromKernel( id, Kernel, data->cIndividuals, data->cIndividuals, data->individualLabels );
      delete[] Kernel;  // Kernel is no longer needed
      }
   else
      {
      // Compute the Eigenvectors and set EigenSym directly from the covarianceSNPs
      eigen->set_covFromX( id, data->snpDataCovariance, data->cIndividuals, data->cSnpsCovariance, data->individualLabels );
      }

   return 0;
   }
      
int runGWAS_lowrank( EigenSym *eigen, ResultsLMM *results, StudyData *data )
   {
   Verbose( verboseMajor, " --             Do lowrank training:" );

   size_t n_individuals = eigen->n;
   size_t n_Snps = data->cSnps;
   size_t rankU = eigen->rankU;
   
   FastLMM_lowrank LMM_lowrank = FastLMM_lowrank( eigen );
   LMM_lowrank.set_Phenotype( data->phenotypeData, data->phenotypeLabels[data->ind_phenotype] );
   LMM_lowrank.set_SNPs( data->snpData, n_Snps, data->snpInfo );
   LMM_lowrank.set_covariates( data->covariatesData, data->cCovariates );
   LMM_lowrank.set_SNPsCov( data->snpDataCovariance, data->cSnpsCovariance, data->snpInfoCovariance );

   LMM_lowrank.init();

   LMM_lowrank.trainNull();

   LMM_lowrank.updateExcludeList();
   if ( runAlternativeModels )
      {
      LMM_lowrank.trainAlt();
      }

   //get the results
   LMM_lowrank.getResults( *results );
   return 0;
   }

int runGWAS_woodbury( EigenSym *eigen, ResultsLMM *results, StudyData *data )
   {
   Verbose( verboseMajor, " --            Do Woodbury training:" );

   size_t n_individuals = eigen->n;
   size_t n_Snps = data->cSnps;
   size_t rankU = eigen->rankU;

   FastLMM_woodbury LMM_woodbury = FastLMM_woodbury(eigen);
   if ( (data->phenotypeLabels.empty()) 
     || (data->phenotypeLabels.size()==1 && data->phenotypeLabels[0].length()==0) 
     || (data->phenotypeLabels[data->ind_phenotype].length()==0))
      {
      string szNoName = "NoName";
      LMM_woodbury.set_Phenotype( data->phenotypeData, szNoName );
      }
   else
      {
      LMM_woodbury.set_Phenotype( data->phenotypeData, data->phenotypeLabels[data->ind_phenotype] );
      }

   LMM_woodbury.set_SNPs( data->snpData, n_Snps, data->snpInfo );
   LMM_woodbury.set_covariates( data->covariatesData, data->cCovariates );
   LMM_woodbury.init();
   LMM_woodbury.trainNull();

   //LMM_woodbury.updateExcludeList();

   if (runAlternativeModels)
      {
      if ( fDoSnpPairs )
         {
         LMM_woodbury.TrainNullAndAltForEpistatic();
         }
      else
         {
         LMM_woodbury.trainAlt();
         }
      }

   //get the results
   LMM_woodbury.getResults( *results );
   Verbose( verboseMinor, " --          Woodbury training done:" );
   return 0;
   }

int runGWAS_fullrank( EigenSym *eigen, ResultsLMM *results, StudyData *data )
   {
   Verbose( verboseMajor, "               Do fullrank training:" );
   size_t n_individuals = eigen->n;
   size_t n_SNPs = data->cSnps;

   FastLMM_fullrank LMM_fullrank = FastLMM_fullrank( eigen );

   LMM_fullrank.set_Phenotype( data->phenotypeData, data->phenotypeLabels[data->ind_phenotype] );
   LMM_fullrank.set_SNPs( data->snpData, n_SNPs, data->snpInfo );
   LMM_fullrank.set_SNPsCov(data->snpDataCovariance, data->cSnpsCovariance, data->snpInfoCovariance);
   LMM_fullrank.set_covariates( data->covariatesData, data->cCovariates );
   
   LMM_fullrank.init();
   
   LMM_fullrank.trainNull();

   LMM_fullrank.updateExcludeList();
   if (runAlternativeModels)
      {
      LMM_fullrank.trainAlt();
      }

   //get the results
   LMM_fullrank.getResults( *results );
   Verbose( verboseMinor, "             fullrank training done:" );
   return 0;
   }


int fillEigenSym( const string& id, EigenSym *eigen, StudyData *data )
   {
   string fullPath = FullPath( eigenDir );

   if ( (computeEigenType==NOSAVE) || (computeEigenType==SAVETODISK) )
      {
      CTimer timer( true );         // Create with timer running

      if ( data->kernelData )
         {
         //kernelData has been set
         computeEigenSym_Kernel( id, data, eigen );      //create EigenSym from the kernel stored in kernelFile
         }
      else if( data->groupSim )
         {
         eigen->U = new real[ data->groupSim->cIndividuals * data->groupSim->cGroups ];
         eigen->S = new real[ data->groupSim->cGroups ];
         eigen->rankU = data->groupSim->cGroups; 
         eigen->n = data->groupSim->cIndividuals;
         eigen->id = id;
         eigen->columnLabelOrder = data->groupSim->individualsIndex;
         EigenSym::EVsymmetric_ZKZ( data->groupSim->kernelGroups->kernelArray, data->groupSim->Z, data->groupSim->cIndividuals, data->groupSim->cGroups, eigen->U, eigen->S );
         }
      else
         {   //a snpDataCovFile has been specified
         computeEigenSym_SnpSample( id, data, eigen );   //create EigenSym of a covariance kernel from the SNP sample stored in snpDataCovFile
         }

      timer.Stop();
      timer.Report( 0x02, "           Time to compute EigenSym: %s" );

      if ( computeEigenType == SAVETODISK )
         {
         //option to save EigenSym
         Verbose( verboseMajor, "              Save eigenSym to file: [%s]", fullPath.c_str() );
         eigen->saveToDiskBINARY( eigenDir );      //Save the EigenSym object for later use in the directory eigenDir
         }
      }
   else
      {
      if( computeEigenType != LOADFROMDISK )
         {
         Fatal( "computeEigenType at non-legal value." );
         }

      Verbose( verboseMajor, "            Load eigenSym from file: [%s]", fullPath.c_str() );
      eigen->readFromDiskBINARY( eigenDir );   //Load an existing EigenSym object from the directory eigenDir
      }
   return 0;
   }
