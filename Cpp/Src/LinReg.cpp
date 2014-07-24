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

#include "FastLmmC.h"

LinReg::~LinReg()
   {
   if ( pvals )      { delete[] pvals; }
   if ( sig2 )       { delete[] sig2; }
   if ( beta )       { delete[] beta; }
   if ( beta_0 )     { delete[] beta_0; }
   if ( LL )         { delete[] LL; }
   if ( f_stat )     { delete[] f_stat; }
   if ( f_stat_0 )   { delete[] f_stat_0; }
   if ( rank_beta )  { delete[] rank_beta; }
   }

LinReg::LinReg( StudyData *pStudyData )
   {
   sig2_0    = 0.0;
   LL_0      = 0.0;
   rank_beta_0 = (BLAS_INT)0;

   data      = pStudyData;
   phenLabel = pStudyData->phenotypeLabels[0];

   pvals     = nullptr;
   sig2      = nullptr;
   beta      = nullptr;
   beta_0    = nullptr;
   LL        = nullptr;
   f_stat    = nullptr;
   f_stat_0  = nullptr;
   rank_beta = nullptr;
   }

void LinReg::trainPairs()
   {
   size_t cIndividuals = this->data->cIndividuals;         // use more descriptive names for now
   size_t cCovariates = this->data->cCovariates;
   size_t cSnps = this->data->cSnps;
   size_t cSnpPairsBetweenDisplay = fClusterIO ? 10000 : 1000;
   char   nl = fClusterIO ? '\n' : '\r';

   // Allocate a working array for covariate data.  We need 
   //   enough for cCovariates * cIndividuals + we need 
   //   to add the two more rows for snp2data and snp12data
   //   NOTE:::  if we change to triples, do we need more than two 
   //            additional rows?  ( cCovariates + 1 + 1 + 1 )
   real* rgCovariatesTemp = new real[ (cCovariates + 3) * cIndividuals ];    // allocate 2d matrix
   
   // Preload the initial covariate data.
   size_t cEntries = cCovariates * cIndividuals;
   memmove( rgCovariatesTemp, this->data->covariatesData, cEntries * sizeof( real ) );    // copy original covar data

   FILE *pFile = fopen( fullOutputFilename.c_str(), "w" );
   if ( pFile == nullptr )
      {
      Fatal( "unable to open output file.  Sorry."
           "\n  CRT Error: %s", strerror( errno ) );
      }

   // NOTE:  GCC and VC don't share the same format specifiers for 
   //        size_t, so we convert to use int in this case.
   // Print the header and setup to measure our throughput
   fprintf( pFile,   "#  Tasks Requested: %9d", (int)snpFilterOpts.cTasksAsk );
   fprintf( pFile, "\n#     Tasks Needed: %9d", (int)snpFilterOpts.cTasksMax );
   fprintf( pFile, "\n#      Task Number: %9d", (int)snpFilterOpts.iTaskAsk );
   fprintf( pFile, "\n#       Partitions: %9d", (int)snpFilterOpts.cBlocksMax );
   fprintf( pFile, "\n#   SNPs/Partition: %9.4lf", snpFilterOpts.cSnpsInBlock );

   size_t iBlk1 = snpFilterOpts.iBlock1;
   size_t iBlk2 = snpFilterOpts.iBlock2;
   fprintf( pFile, "\n#      Block1 Info: %9d %9d %9d", (int)iBlk1, (int)snpFilterOpts.IdxFromBlockId( iBlk1 ), (int)snpFilterOpts.CountOfSnpsInBlock( iBlk1 ) );
   fprintf( pFile, "\n#      Block2 Info:" ); // keep a constant number of header lines.
   if ( iBlk1 != iBlk2 )
      {
      fprintf( pFile, " %9d %9d %9d", (int)iBlk2, (int)snpFilterOpts.IdxFromBlockId( iBlk2 ), (int)snpFilterOpts.CountOfSnpsInBlock( iBlk2 ) );
      }

   fprintf( pFile, "\nsnpId_1\tsnpId_2\tpvalue\tSnpPair_Wgt\tSnpPair_OddsRatio\n" );

   size_t iSnp1Base = snpFilterOpts.IdxFromBlockId( iBlk1 );
   size_t iSnp2Base = snpFilterOpts.IdxFromBlockId( iBlk2 );
   size_t cSnps1Max = snpFilterOpts.CountOfSnpsInBlock( iBlk1 );
   size_t cSnps2Max = snpFilterOpts.CountOfSnpsInBlock( iBlk2 );
   size_t cSnpPairs = 0;
   size_t cSnpPairsWritten = 0;
   CTimer snpTimer( true );

   if ( fDoSnpId1 )
      {
      size_t iSnp1 = IndexFromList( data->snpInfo, snpId1Global );
      if ( iSnp1 == -1 )
         {
         Fatal( "Unable to locate SnpId [%s] in list of SNPs.", snpId1Global.c_str() );
         }

      Progress( "\n       Processing SNP pairs for SNP: %s", snpId1Global.c_str() );
      Progress( "\n                SnpPairs to process: %7d", (int)((cSnps1Max - 1 )) );   // take out the snpId1
      Progress( "\n SnpPair %9d [%7.2lf SnpPairs/sec]", (int)cSnpPairs, (double)0.0 );

      for ( size_t iSnp2 = 0; iSnp2 < cSnps; ++ iSnp2 )
         {
         if ( iSnp2 != iSnp1 )
            {
            cSnpPairsWritten += TrainNullAndAltForEpistaticWithSnpPairs( pFile, iSnp1Base, iSnp1, iSnp1Base, iSnp2, rgCovariatesTemp );
            if ( (++cSnpPairs % cSnpPairsBetweenDisplay) == 0 )
               {
               Progress( "%c SnpPair %9d [%7.2lf SnpPairs/sec]", nl, (int)cSnpPairs, (double)cSnpPairs / snpTimer.Elapsed() );
               }
            }
         }
      }
   else if ( iBlk1 != iBlk2 )
      {  //                                           Block Base count
      Progress( "\n           Number of blocks in File: %7d", (int)snpFilterOpts.cBlocksMax );
      Progress( "\n       Processing 2 blocks for Task: %7d of %7d", (int)snpFilterOpts.iTaskAsk, (int)snpFilterOpts.cTasksMax );
      Progress( "\n                       Block 1 Info: %7d  %7d  %7d", (int)iBlk1, (int)iSnp1Base, (int)cSnps1Max );
      Progress( "\n                       Block 2 Info: %7d  %7d  %7d", (int)iBlk2, (int)iSnp2Base, (int)cSnps2Max );
      Progress( "\n                SnpPairs to process: %7d", (int)(cSnps1Max * cSnps2Max ) );
      Progress( "\n SnpPair %9d [%7.2lf SnpPairs/sec]", cSnpPairs, (double)0.0 );

      size_t iSnp2BaseOffset = iSnp2Base - cSnps1Max;    // account for Snp2 location in the combined array

      for ( size_t iSnp1 = 0; iSnp1 < cSnps1Max; ++iSnp1 )
         {
         for( size_t iSnp2 = cSnps1Max; iSnp2 < cSnps; ++iSnp2 )
            {
            cSnpPairsWritten += TrainNullAndAltForEpistaticWithSnpPairs( pFile, iSnp1Base, iSnp1, iSnp2BaseOffset, iSnp2, rgCovariatesTemp );
            if ( (++cSnpPairs % cSnpPairsBetweenDisplay) == 0 )
               {
               Progress( "%c SnpPair %9d [%7.2lf SnpPairs/sec]", nl, (int)cSnpPairs, (double)cSnpPairs / snpTimer.Elapsed() );
               }
            }
         }
      }
   else
      {  //                                           Block Base count
      Progress( "\n           Number of blocks in File: %7d", (int)snpFilterOpts.cBlocksMax );
      Progress( "\n        Processing 1 block for Task: %7d of %7d", (int)snpFilterOpts.iTaskAsk, (int)snpFilterOpts.cTasksMax );
      Progress( "\n                       Block 1 Info: %7d  %7d  %7d", (int)iBlk1, (int)iSnp1Base, (int)cSnps1Max );
      Progress( "\n                       Block 2 Info:" ); // maintain consistent line count for output
      Progress( "\n                SnpPairs to process: %7d", (int)((cSnps1Max * (cSnps1Max-1)) / 2) );
      Progress( "\n SnpPair %9d [%7.2lf SnpPairs/sec]", (int)cSnpPairs, (double)0.0 );

      for ( size_t iSnp1 = 0; iSnp1 < cSnps-1; ++iSnp1 )
         {
         for ( size_t iSnp2 = iSnp1 + 1; iSnp2 < cSnps; ++ iSnp2 )
            {
            cSnpPairsWritten += TrainNullAndAltForEpistaticWithSnpPairs( pFile, iSnp1Base, iSnp1, iSnp1Base, iSnp2, rgCovariatesTemp );
            if ( (++cSnpPairs % cSnpPairsBetweenDisplay) == 0 )
               {
               Progress( "%c SnpPair %9d [%7.2lf SnpPairs/sec]", nl, (int)cSnpPairs, (double)cSnpPairs / snpTimer.Elapsed() );
               }
            }
         }
      }

   // Report the end numbers
   fprintf( pFile,"\n** End **\n%d,%d\n", (int)cSnpPairsWritten, (int)cSnpPairs );
   fclose( pFile );
   snpTimer.Stop();
   Progress( "%c SnpPair %9d [%7.2lf SnpPairs/sec]\n", nl, (int)cSnpPairs, (double)cSnpPairs / snpTimer.Elapsed() );

   delete [] rgCovariatesTemp;
   
}

size_t LinReg::TrainNullAndAltForEpistaticWithSnpPairs( FILE *pFile, size_t iSnp1Base, size_t iSnp1, size_t iSnp2Base, size_t iSnp2, real* rgCovariatesTemp ){
   // in order to have the least chance of misusing the univariate SNP code, 
   // I will aribitrarily put snp1 in for the "single" snp in the univariate
   // test, and put the second snp, and the interaction effect with the covariates.
   size_t d = 1;
   size_t cIndividuals = data->cIndividuals;         // use more descriptive names for now
   size_t cCovariates = data->cCovariates;
   real* X = data->snpData;

   // pointers to data for Snp1 & Snp2 (one cIndividual's worth)
   real* snp1data = X + (iSnp1 * (cIndividuals));
   real* snp2data = X + (iSnp2 * (cIndividuals));

   // we are re-using the same memory locations since we don't want/need
   // to store all the information for now.
   real* currentBeta_0  =  new real[ (data->cCovariates + 2) ];
   real* currentBeta  =  new real[ (data->cCovariates + 3) ];
   real currentSige2 = 0.0;
   real * currentF_stat_0 =  new real[ (data->cCovariates + 2) ];
   real * currentF_stat =  new real[ (data->cCovariates + 3) ];
   BLAS_INT currentRank_beta = 0;
   
   // Evaluate the null LL (includes snp1 + snp2) at delta_0 and 
   //   get the ML parameters
   //----------------------------------------------------------------------
   // cRows = cIndividuals = n
   // cCols = cSnps = ncovar

   // We filled in the 'updated' covariates array with the original covariate
   //  data in the caller.
   // We now need to append snp2data to the end of the Covariates array we use.
   size_t crgCovariatesEntries = cCovariates * cIndividuals;
   
   //  The initial fill of the working matrix is now done once in calling routine.
   //memmove( rgCovariatesTemp, covar, cEntries * sizeof( real ) );    // copy original covar data
   memmove( rgCovariatesTemp + crgCovariatesEntries                , snp1data, cIndividuals * sizeof( real ) );    // add snp1data
   memmove( rgCovariatesTemp + crgCovariatesEntries + cIndividuals , snp2data, cIndividuals * sizeof( real ) );    // add snp2data

   size_t cCovariatesTemp = cCovariates + 2;             // include snp2data data for now.

   //evaluate the LL and get the ML parameters
   real LL_null = LinReg::evalLL( (BLAS_INT)cIndividuals, (BLAS_INT)cCovariatesTemp, data->phenotypeData, currentBeta_0, rgCovariatesTemp, &currentSige2, currentF_stat_0, &currentRank_beta );

   //evaluate the alternative  LL (includes snp1 + snp2 + snp1*snp2) and get the ML parameters
   //--------------------------------------------------------------------------------------------
   //modify covar to also include the interaction effect, snp1*snp2
   //specifically, re-use the concatenation we had above for the null, but add in one more feature, snp1*snp2

   // create the SNP interaction effect feature vector as the last row of the rgCovariatesTemp array
   //  this should be set equal to the element-wise product of snp1data and snp2data:  
   //  i.e. for i=1:length(snp1data) snp12data=snp1data[i]*snp2data[i] end
   real* snp12data = (rgCovariatesTemp + crgCovariatesEntries) + cIndividuals + cIndividuals;      // point to last row of table e.g. [cCovariates + 2]
   for ( size_t iIndividual = 0; iIndividual < cIndividuals; ++ iIndividual )
      {
      snp12data[ iIndividual ] = snp1data[ iIndividual ] * snp2data[ iIndividual ];
      }

   cCovariatesTemp += 1;         // Account for adding in the snp12data 

   //evaluate the alternate LL and get the ML parameters
   real LL_alt = LinReg::evalLL( (BLAS_INT)cIndividuals, (BLAS_INT)cCovariatesTemp, data->phenotypeData, currentBeta, rgCovariatesTemp, &currentSige2, currentF_stat, &currentRank_beta );

   // Compute p-value from LR test:
   // real df = *currentRank_beta - (real)rank_beta_0[ snp ];
   real df = 1.0; // TODO: Review
   real pval;

   //Likelihood-ratio test
   pval = Gamma::gammaQ( LL_alt - LL_null, (real)0.5 * df );

   // fix NaN output
   if ( pval != pval )
      {
      Warn( "SnpPair[%s - %s] pval is a NaN.  Setting pval to 1.00e-000", data->snpInfo[iSnp1].idSnp.c_str(), data->snpInfo[iSnp2].idSnp.c_str() );
      pval = 1.00e-000;
      }

   // Prepare the odds ratios for printing
   //   snpWeight is last element of currentBeta[] (which is one larger than cCovariatesTemp)
   //   multiply snp1 & snp2 StdDev to get StdDev for SnpPairs
   real snpPairWeight = currentBeta[ cCovariatesTemp ];  
   real snpPairStdDeviation = data->snpInfo[ iSnp1 ].StdDeviation() * data->snpInfo[ iSnp1 ].StdDeviation();   
   real snpPairOddsRatio = exp( ::abs( snpPairWeight ) / snpPairStdDeviation );

   fprintf( pFile, "%s\t%s\t%e\t%e\t%e\n", 
                        data->snpInfo[iSnp1].idSnp.c_str(), data->snpInfo[iSnp2].idSnp.c_str(),
                        pval,
                        snpPairWeight, snpPairOddsRatio
                        );

   delete[] currentBeta;
   delete[] currentF_stat;
   delete[] currentBeta_0;
   delete[] currentF_stat_0;

   return( 1 );   // number of SnpPairs printed ( 0 or 1 )
   }

void LinReg::trainNull()
   {
   if ( beta_0 )   { delete[] beta_0; }
   if ( f_stat_0 ) { delete[] f_stat_0; }

   beta_0   = new real[ this->data->cCovariates ];
   f_stat_0 = new real[ this->data->cCovariates ];
   LL_0     = this->evalLL( std::numeric_limits<size_t>::quiet_NaN(), false );

   if ( verboseLevel >= verboseDebug )
      {
      printf( "\n------------------------------------------"
              "\nNULL: LL:%f\nsig2=%f",
              LL_0, sig2_0 );
      printf( " beta:" );
      for( size_t i = 0; i < data->cCovariates; ++i )
         {
         printf( "%f ", beta_0[ i ] );
         }
      printf( "\n------------------------------------------"
              "\n" );
      }
   }

void LinReg::trainAlt()
   {
   if ( LL )         { delete[] LL; }
   if ( rank_beta )  { delete[] rank_beta; }
   if ( pvals )      { delete[] pvals; }
   if ( beta )       { delete[] beta; }
   if ( sig2 )       { delete[] sig2; }
   if ( f_stat )     { delete[] f_stat; }
   rank_beta = new BLAS_INT[ this->data->cSnps ];
   f_stat    = new real[ (this->data->cCovariates+1) * this->data->cSnps ];
   beta      = new real[ (this->data->cCovariates+1) * this->data->cSnps ];
   sig2      = new real[ this->data->cSnps];
   pvals     = new real[ this->data->cSnps ];
   LL        = new real[ this->data->cSnps ];

   for ( size_t snp=0; snp < this->data->cSnps; ++snp )
      {
      LL[ snp ] = this->evalLL( snp, true );
      
      //Compute p-value from LR test:
      real df = rank_beta[ snp ] - (real)1.0 * rank_beta_0;
      if ( df > (real)0.0 )
         {
         //check the degrees of freedom for validity
         if ( use_f_test )
            {
            //Fishers F-test
            //
            // PvalUsingFisherFTestCdf focuses on keeping precision 
            //    for very small pvalues by inverting on the reverse 
            //    path of the Beta::BetaInc() branch
            pvals[ snp ] = PvalUsingFisherFTestCdf( f_stat[ snp * (data->cCovariates+1) ], df, (real)(data->cIndividuals - rank_beta[ snp ]) );
            //pvals[ snp ] = 1.0 - FisherF::Cdf( f_stat[ snp * (data->cCovariates+1) ], df, (real)(data->cIndividuals - rank_beta[ snp ]) );
            }
         else
            {
            //Likelihood-ratio test
            pvals[ snp ] = PvalFromLikelihoodRatioTest( LL[ snp ] - LL_0, (real)0.5 * df );
            }
         }
      else
         {
         //This happens for example if the SNP is a constant. In this case df=0
         pvals[ snp ] = std::numeric_limits<real>::quiet_NaN();
         }

      if ( verboseLevel >= verboseDebug )
         {
#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
         const char *szFmt1 = "\n- LL[%Iu]:%f p=%f\nsig2=%f ";
#else                      // Linux/g++ uses a %zu specifier for size_t
         const char *szFmt1 = "\n- LL[%zu]:%f p=%f\nsig2=%f ";
#endif
         printf( szFmt1, snp, LL[snp], pvals[snp], sig2[snp] );
         printf( "beta:\n " );
         for( size_t i=0; i<(data->cCovariates+1); ++i )
            {
            printf( "%f ", beta[ snp*(data->cCovariates+1)+i ] );
            }
         }
      }
   }

real LinReg::evalLL( size_t ind_SNP, bool altModel )
   {
   BLAS_INT n  = (BLAS_INT)data->cIndividuals;
   BLAS_INT d_ = (BLAS_INT)data->cCovariates;
   real *y     = data->phenotypeData;  // point to phenotypeData in LinReg.data
   real *beta_ = nullptr;     // beta hides LinReg::beta.  Use local beta_
   real *X     = nullptr;
   real *sig2_ = nullptr;     // sig2 hides LinReg::sig2.  Use local sig2_
   real *f_stat_ = nullptr;
   BLAS_INT *rank_beta_ = nullptr;

   if ( altModel )
      {
      //alternative model includes a SNP
      ++d_;
      X = new real[ n * d_ ];
      BlasLapack::copy( &n, &data->snpData[ ind_SNP * n ], &int_one, X, &int_one );

      for(size_t i_covar=0; i_covar < (data->cCovariates * data->cIndividuals); i_covar += data->cIndividuals )
         {
         BlasLapack::copy( &n, (data->covariatesData + i_covar), &int_one, (X + data->cIndividuals + i_covar), &int_one );
         }

      beta_ = &this->beta[ (data->cCovariates+1) * ind_SNP ];
      f_stat_ = &this->f_stat[ (data->cCovariates+1) * ind_SNP ];
      sig2_ = &this->sig2[ ind_SNP ];
      rank_beta_ = &this->rank_beta[ ind_SNP ];
      }
   else
      {
      //null model
      X = data->covariatesData;
      beta_ = this->beta_0;
      sig2_ = &this->sig2_0;
      rank_beta_ = &this->rank_beta_0;
      f_stat_ = this->f_stat_0;
      }

   real result =  LinReg::evalLL( n, d_, y, beta_, X, sig2_, f_stat_, rank_beta_ );

   if( altModel )
      {
      // alternative model includes a SNP
      //   and we had to allocate this memory
      delete[] X;
      }

   return result;
   }

real LinReg::evalLL( BLAS_INT n, BLAS_INT d_, real* y, real* beta_, real*X, real*sig2_, real *f_stat_, BLAS_INT *rank_beta_ )
   {
   real *XX = new real[ d_ * d_ ];     //Warning: size_t BLAS_INT conversion product
   memset( XX, 0, (d_ * d_ * sizeof( real )) );
   real *Xy = new real[ d_ ];
   memset( Xy, 0, (d_ * sizeof( real )) );

   if ( d_ > 0 )
      {
      BlasLapack::syrk( &char_L, &char_T, &d_, &n, &real_one, X, &n, &real_zero, XX, &d_ );
      BlasLapack::gemv( &char_T, &n, &d_, &real_one, X, &n, y, &int_one, &real_zero, Xy, &int_one );
      }

   real logDetXX = 0.0;
   FastLMM::solveLstSqSym( d_, XX, Xy, beta_, rank_beta_, f_stat_, &logDetXX );

   real res = 0.0;
   for ( BLAS_INT i = 0; i < n; ++i )
      {
      res += y[ i ] * y[ i ];
      }
   for ( BLAS_INT i = 0; i < d_; ++i )
      {
      res -= Xy[ i ] * beta_[ i ];
      }

   *sig2_ = res / n;
   for ( BLAS_INT dim = 0; dim < d_; ++dim )
      {
      f_stat_[ dim ] /= *sig2_;
      }
   real resultLL = -(real)0.5 * (n * (_log2pi + log( *sig2_ ) + (real)1.0));   //compute the negative LL using the ML parameters at delta
   delete[] XX;
   delete[] Xy;

   return( resultLL );
   }

void LinReg::WriteResultFile( const string& resultFile )
   {
   // CLEANUP: 
   WriteResultFile_V2( resultFile, fVerboseOutput, pvaluePrintThreshold );
   }

vector< string > LinRegColumnNames( bool fVerboseOutput_, bool fUseBiasTerm, size_t cCovariates )
   {
   vector< string > columnNames;
   columnNames.reserve( 20 );          // pre-allocate enough for our currently known cases

   columnNames.push_back( "SNP" );
   columnNames.push_back( "Chromosome" );
   columnNames.push_back( "GeneticDistance" );
   columnNames.push_back( "Position" );

   if ( fVerboseOutput_ )
      {
      columnNames.push_back( "Phenotype" );
      }

   columnNames.push_back( "Pvalue" );
   columnNames.push_back( "Qvalue" );
   columnNames.push_back( "N" );
   columnNames.push_back( "NullLogLike" );
   columnNames.push_back( "AltLogLike" );
   columnNames.push_back( "SnpWeight" );
   columnNames.push_back( "SnpWeightSE" );      // Standard Error
   columnNames.push_back( "OddsRatio" );        // Odds Ratio = exp( abs(SnpWeight) / SnpWeightSE )
   columnNames.push_back( "WaldStat" );

   if ( fUseBiasTerm )
      {
      columnNames.push_back( "NullBias" );
      }

   // if fUseBiasTerm is used, it changes cCovariates.  We need to adjust for it
   size_t adjustForBiasTerm = fUseBiasTerm ? 1 : 0;

   for ( size_t iCovariate = 0; iCovariate < cCovariates - adjustForBiasTerm; ++iCovariate )
      {
      // iterate to generate 'NullCov%02dWeight' column names
      char szNullCovWeight[32];
      sprintf( szNullCovWeight, "NullCov%02dWeight", (int)(iCovariate) );
      columnNames.push_back( szNullCovWeight );
      }

   if ( fVerboseOutput_ )
      {
      for ( size_t iCovariate = 0; iCovariate < cCovariates - adjustForBiasTerm; ++iCovariate )
         {
         // iterate to generate 'Cov%02dWeight' column names
         char szCovWeight[ 32 ];
         sprintf( szCovWeight, "Cov%02dWeight", (int)(iCovariate) );
         columnNames.push_back( szCovWeight );
         }
      }

   columnNames.push_back( "NullVar" );
   return( columnNames );
   }

// Get the index of SNPs ordered by significance
vector< size_t > LinReg::GetOrderedSNPIndex()
   {
   // REVIEW:  
   //   Potential memory leak for qvals and index?
   //     new[] w/o corresponding delete[]
   //     returning new[] index.  Who is calling delete[] on index?
   //
   vector< size_t >index( data->cSnps );
   real *qvals = new real[ data->cSnps ];
   Qval::computeQvals( this->pvals, qvals, &index[0], data->cSnps );
   return( index );
   }
 
vector< SnpAndPvalue > *LinReg::GetOrderedSnpAndPvalueInfo( size_t limit )
   {
   // REVIEW:
   //   Potential memory leak for index, qvals, and ordered?
   //     new[] w/o corresponding delete[]
   //     returning new vector<> ordered.  Who is calling delete ordered?
   //
   vector< SnpAndPvalue > *ordered = new vector< SnpAndPvalue >();
   size_t *index = new size_t[ data->cSnps ];
   real   *qvals = new real[ data->cSnps ];

   Qval::computeQvals( pvals, qvals, index, data->cSnps );

   SnpAndPvalue sinfo;     // declare here to avoid construction/destruction in loop

   // Copy the values to SnpAndPvalue
   for( size_t i=0; i<limit; ++i )
      {
      size_t ind = index[ i ];
      sinfo.idSnp = data->snpInfo[ ind ].idSnp;
      sinfo.pval = pvals[ ind ];
      ordered->push_back( sinfo );
      }

   return( ordered );
   }

void LinReg::WriteResultFile_V2( const string& resultFile, bool fVerboseOutput_, real maxPvaluePrintThreshold )
   {
   /*
    *  Output a tab separated or comma separated file with the results of the FastLmm compuatation
    *  The first row contains column headers and each subsequent row contains
    *  the information for a specific SNP.
    *
    */
   void SelectFormatStrings_V2( const char*& fmtInt_, const char*& fmtString_, const char*& fmtReal_ );  // borrow from FastLmm.cpp
   FILE *OpenResultFile( const string& resultFile );
   void OutputColumnNames( FILE *pFile, const vector< string >& names, const char* separator );

   const char *szSeparator = FIsCsvFile( resultFile ) ? "," : "\t";
   const char *fmtInt    = nullptr;
   const char *fmtString = nullptr;
   const char *fmtReal   = nullptr;
   SelectFormatStrings_V2( fmtInt, fmtString, fmtReal );

   // Compute, store, and index all the final values to be output
   size_t cSnps = data->cSnps;
   size_t *indexSnps = new size_t[ cSnps ];
   real   *qvals = new real[ cSnps ];
   Qval::computeQvals( this->pvals, qvals, indexSnps, cSnps );

   size_t cCovariates = data->cCovariates;
   vector< string > columnNames = LinRegColumnNames( fVerboseOutput_, useBiasTerm, cCovariates );
   FILE *pFile = OpenResultFile( resultFile );
   OutputColumnNames( pFile, columnNames, szSeparator );

   vector< SnpInfo >& rgSnpInfo = data->snpInfo;
   size_t cSnpsFailedPvaluePrintThreshold = 0;

   // write the rows
   for( size_t i=0; i < rgSnpInfo.size(); ++i )
      {
      size_t iSnp = fEmitSorted ? indexSnps[ i ] : i;
      real pvalT = pvals[ iSnp ];
      if ( pvalT > maxPvaluePrintThreshold )
         {
         ++cSnpsFailedPvaluePrintThreshold;
         continue;      // skip this SNP
         }

      SnpInfo& snp = rgSnpInfo[ iSnp ];
      fprintf( pFile, "%s",                 snp.idSnp.c_str() );        // SNP - the SNP label (no leading separator
      fprintf( pFile, fmtInt,  szSeparator, snp.iChromosome );          // Chromosome - the chromosome number
      fprintf( pFile, "%s%f",  szSeparator, snp.geneticDistance );      // GeneticDistance - (morgans)
      fprintf( pFile, fmtInt,  szSeparator, snp.basepairPosition );     // Position - 

      if ( fVerboseOutput_ )
         {
         fprintf( pFile, fmtString, szSeparator, phenLabel.c_str() );
         }

      if ( _isnan( pvalT ) )
         {
         // this is a NaN!
         Warn( "Found a NaN pvalue for SNP [%s] at row %d", snp.idSnp.c_str(), i );
         fprintf( pFile, fmtString, szSeparator, "NA" );    // NaN ?
         }
      else
         {
         // valid pvalue - output the value
         fprintf( pFile, fmtReal, szSeparator, pvalT );
         }

      if ( _isnan( qvals[ iSnp ] ) )
         {
         // this is a NaN!
         Warn( "Found a NaN qvalue for SNP [%s] at row %d", snp.idSnp.c_str(), i );
         fprintf( pFile, fmtString, szSeparator, "NA" );    // NaN ?
         }
      else
         {
         // Qvalue - output proper qvalues
         fprintf( pFile, fmtReal, szSeparator, qvals[ iSnp ] );
         }
      
      fprintf( pFile, fmtInt,  szSeparator, data->cIndividuals );       // N - sample size used in analysis (cIndividuals)
      fprintf( pFile, fmtReal, szSeparator, LL_0 );                     // NullLogLike - likelihood of null model
      fprintf( pFile, fmtReal, szSeparator, LL[ iSnp ] );               // AltLogLike - likelihood of alternative model

      // Compute index into the beta[cCovariates+1][cSnps] matrix to get the snpWeight
      //   it is stored in the 0th element of the iSnp vector
      size_t index = ColumnMajorIndex( cCovariates + 1, cSnps, 0, iSnp );
      real snpWeight = beta[ index ];                                         // SNPWeight -> right from beta vector
      fprintf( pFile, fmtReal, szSeparator, snpWeight );                      // SNPWeight - fixed-effect weight of SNP

      real snpWeightSE = sqrt( (snpWeight * snpWeight) / (f_stat[ index ]) ); // SNPWeightSE -> SNP weight standard error
      fprintf( pFile, fmtReal, szSeparator, snpWeightSE );                    // SNPWeightSE -> SNP weight standard error

      real snpOddsRatio = exp( ::abs( snpWeight ) /  snp.StdDeviation() );
      fprintf( pFile, fmtReal, szSeparator, snpOddsRatio );                   // SNPOddsRatio -> Odds ratio (exp( |snpWeight| / snpStdDeviation )

      real waldStat = (snpWeight * snpWeight) / (snpWeightSE * snpWeightSE);  // WaldStat = (SnpWeight / stdErr)^2
      fprintf( pFile, fmtReal, szSeparator, waldStat );                       // WaldStat = (SnpWeight / stdErr)^2

      if ( useBiasTerm )   
         {
         // if biasTerm is used, last element of beta_0[] is the NullBias so get it
         fprintf( pFile, fmtReal, szSeparator, beta_0[ cCovariates - 1 ] );
         }

      // if fUseBiasTerm is used, it changes cCovariates.  We need to adjust for it
      size_t adjustForBiasTerm = useBiasTerm ? 1 : 0;    

      for ( size_t iCovariate = 0; iCovariate < cCovariates - adjustForBiasTerm; ++iCovariate )
         {
         fprintf( pFile, fmtReal, szSeparator, beta_0[ iCovariate ] );     // NullCov<xx>Weight
         }

      if ( fVerboseOutput_ )
         {
         // iCovariate starts at 1 to account for beta[0] being the SNPWeight 
         // and take useBiasTerm into account at the end where the last element is the biasTerm
         for ( size_t iCovariate = 0; iCovariate < cCovariates - adjustForBiasTerm; ++iCovariate )
            {
            size_t index = ColumnMajorIndex( cCovariates + 1, cSnps, iCovariate + 1, iSnp );    // Cov<xx>Weight
            fprintf( pFile, fmtReal, szSeparator, beta[ index ] );
            }
         }

      fprintf( pFile, fmtReal, szSeparator, sig2_0 );                   // NullVar - 
      fprintf( pFile, "\n" );
      }

   if ( (maxPvaluePrintThreshold != 1.0) && (cSnpsFailedPvaluePrintThreshold > 0) )
      {
#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
      const char *szFmt = "\n"
                          "\nCount of SNPs that passed print threshold: %Iu" 
                          "\nCount of SNPs that failed print threshold: %Iu"
                          "\n";
#else                      // Linux/g++ uses a %zu specifier for size_t
      const char *szFmt = "\n"
                          "\nCount of SNPs that passed print threshold: %zu" 
                          "\nCount of SNPs that failed print threshold: %zu"
                          "\n";
#endif
      fprintf( pFile, szFmt, cSnps - cSnpsFailedPvaluePrintThreshold, cSnpsFailedPvaluePrintThreshold );
      }

   fclose( pFile );
   delete[] qvals;
   delete[] indexSnps;
   }

