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

#if defined( USE_LINUX_V2x_KERNEL )
/*
 *  TODO:
 *  The vec.h code used by Logistic Regression is not supported by gcc until 4.6.+
 *    RedHat RHEL is still using a 2.6 kernel and gcc 4.4 toolset.
 *    Re-enable logistic regression on Linux after RHEL moves to 3.x kernel
 */
#else
#include <numeric>
#include <functional>
#include "minimize.h"

/* ************************************************************************* */
// Tuning variables for the optimizer
static const double lbfgsTolerance = (double)1.0e-12;
static const unsigned int cMaxIter = 200;
static const unsigned int cLbfgsHistory = 20;

/* ************************************************************************* */

double CLogisticRegression::ProbRoutine( vector< double >::iterator itWgts, vector< double >::iterator itGradient )
      {
      // vZ = MACC( vMatrix[], itWgt[] )     // multiply and accumlate
      for ( size_t iIndividual=0; iIndividual<cColumns; ++iIndividual )
         {
         vZ[ iIndividual ] = real0;
         auto itT = itWgts;
         for ( size_t iVar=0; iVar<cVars; ++iVar )
            {
            vZ[ iIndividual ] += vMatrix[ iVar ][ iIndividual ] * itT[ iVar ];
            }
         }

      //  vP[]  = 1 + exp( -vZ[] )
      //  vLL[] = log(  ? vP[] : 1-vP[]
      //  vErr[] = vP[] - vPhenotype[]
      for ( size_t iIndividual=0; iIndividual<cColumns; ++iIndividual )
         {
         vP[ iIndividual ] = real1 / ( real1 + exp( -vZ[ iIndividual ] ) );
         vLL[ iIndividual ] = log( (vPhenotype[ iIndividual ] == real1) ? vP[ iIndividual ] : 1-vP[ iIndividual ] );
         vErr[ iIndividual ] = vP[ iIndividual ] - vPhenotype[ iIndividual ];
         }

      auto itGradientT = itGradient;
      for ( int iWeight = 0; iWeight<cVars; ++iWeight )
         {
         *itGradientT = real0;
         for ( size_t iIndividual=0; iIndividual<cColumns; ++iIndividual )
            {
            *itGradientT += ( vErr[ iIndividual ] * vMatrix[ iWeight ][ iIndividual ] );
            }
         ++itGradientT;
         }

      sumLL = accumulate( vLL.begin(), vLL.end(), real0 );
      return( -sumLL );
      }

/* ************************************************************************* */

CLogisticRegression::~CLogisticRegression()
   {
   if ( pvals )      { delete[] pvals; }
   if ( beta )       { delete[] beta; }
   if ( beta_0 )     { delete[] beta_0; }
   if ( LL )         { delete[] LL; }
   if ( rank_beta )  { delete[] rank_beta; }

   // g++ does not compile the unique_ptr<>
   if ( pMatrixData ) { delete[] pMatrixData; }
   }

bool FHasOnlyAorB( real a, real b, real *p, size_t cReals, size_t *idx = nullptr )
   {
   for ( size_t iReal = 0; iReal<cReals; ++iReal )
      {
      real rT = p[ iReal ];
      if ( (rT != a) && (rT != b)
         )
         {
         if ( idx ) *idx = iReal;
         return( false );
         }
      }
   if ( idx ) *idx = cReals;
   return( true );
   }

CLogisticRegression::CLogisticRegression( StudyData *pStudyData )
   {
   // Clear the local data on initialization
//   sig2_0    = real0;
   LL_0      = real0;
   rank_beta_0 = (BLAS_INT)0;

   data      = pStudyData;
   phenLabel = pStudyData->phenotypeLabels[0];

   pvals     = nullptr;
//   sig2      = nullptr;
   beta      = nullptr;
   beta_0    = nullptr;
   LL        = nullptr;
//   f_stat    = nullptr;
//   f_stat_0  = nullptr;
   rank_beta = nullptr;

   // add the LogisticRegression support
   // Allocate vMatrix, and other vectors
   cRows = data->cCovariates + 1;      // save space for SNP on the end
   cColumns = data->cIndividuals;
   
   //pMatrixData = unique_ptr< real[] >(new real[ cRows * cColumns ]);
   pMatrixData = new real[ cRows * cColumns ];

   vMatrix.resize( cRows );
   for( size_t iRow=0; iRow<cRows; ++iRow )
      {
      vMatrix[ iRow ].resize( cColumns, real0 );
      }

   vPhenotype.resize( cColumns, real0 );
   vWeights.resize( cRows, real0 );
   vZ.resize( cColumns, real0 );
   vP.resize( cColumns, real0 );
   vLL.resize( cColumns, real0 );
   vErr.resize( cColumns, real0 );

   /*
    *  Validate the phenotype data for use in logistic regression
    *  The data going into logistic regression must be 0 or 1
    *  However phenotype data can be 0/1 or 1/2
    *  If it is 1/2 we have to 'adjust' it.
    *  PLINK values can have 1 or 2 for phenotypes meaning 
    *       1 = nonaffected
    *       2 = affected
    *  Fatal if this is NOT the case
    *
    *  Then copy 'adjusted' phenotype value to phenotype vector (for now)
    */

   size_t iIndividual = 0;
   real rDiff = real0;
   if ( !FHasOnlyAorB( real0, real1, data->phenotypeData, data->cIndividuals, &iIndividual ) )
      {
      rDiff = real1;
      if ( !FHasOnlyAorB( real1, real2, data->phenotypeData, data->cIndividuals, &iIndividual ) )
         {
         Fatal( "Malformed phenotype for [%s].  Expected 1 or 2 but found [%.2e]", 
                  data->individualLabels[ iIndividual ].c_str(),
                  data->phenotypeData[ iIndividual ] );
         }
      }

   for( iIndividual=0; iIndividual<data->cIndividuals; ++iIndividual )
      {
      vPhenotype[ iIndividual ] = data->phenotypeData[ iIndividual ] - rDiff;
      }

   // copy covariate data for processing
   vMatrix.resize( data->cCovariates + 1 );    // include room for one row of SNP data
   vector< real >vT;
   vT.resize( data->cIndividuals );
   for ( size_t iCovariate=0; iCovariate<data->cCovariates; ++iCovariate )
      {
      for ( size_t iIndividual=0; iIndividual<data->cIndividuals; ++iIndividual )
         {
         size_t idx = ColumnMajorIndex( data->cIndividuals, data->cCovariates, iIndividual, iCovariate );
         vT[ iIndividual ] = data->covariatesData[ idx ];
         }
      vMatrix[ iCovariate ] = vT;
      }
   }

void CLogisticRegression::TrainNull()
   {
   if ( beta_0 )   { delete[] beta_0; }
//   if ( f_stat_0 ) { delete[] f_stat_0; }

   beta_0   = new real[ data->cCovariates ];
//   f_stat_0 = new real[ data->cCovariates ];
   size_t cbFill = sizeof( real ) * data->cCovariates;
   memset( beta_0, 0, cbFill );
//   memset( f_stat_0, 0, cbFill );

   LL_0 = EvalLL( 0, false );          // EvalLL on the Null Model path

   if ( verboseLevel >= verboseDebug )
      {
      printf( "\n------------------------------------------"
              "\nNULL: LL:%f",
              LL_0 );
      printf( " beta:" );
      for( size_t i = 0; i < data->cCovariates; ++i )
         {
         printf( "%f ", beta_0[ i ] );
         }
      printf( "\n------------------------------------------"
              "\n" );
      }
   }

void CLogisticRegression::TrainAlt()
   {
   // TODO:  WHY would these have a value needing to be deleted???
   if ( LL )         { delete[] LL; }
   if ( rank_beta )  { delete[] rank_beta; }
   if ( pvals )      { delete[] pvals; }
   if ( beta )       { delete[] beta; }

   rank_beta = new BLAS_INT[ data->cSnps ];
   beta      = new real[ (data->cCovariates+1) * data->cSnps ];
   pvals     = new real[ data->cSnps ];
   LL        = new real[ data->cSnps ];

   // initialize them all to 0 or 0.0
   memset( rank_beta, 0, data->cSnps * sizeof( BLAS_INT ) );
   memset( beta, 0, (data->cCovariates+1) * data->cSnps * sizeof( real ) );
   memset( pvals, 0, data->cSnps * sizeof( real ) );
   memset( LL, 0, data->cSnps * sizeof( real ) );

   for ( size_t iSnp=0; iSnp < data->cSnps; ++iSnp )
      {
      LL[ iSnp ] = EvalLL( iSnp, true );
      
      // check the degrees of freedom for validity and use Likelihood-ratio test
      real df = (real)(rank_beta[ iSnp ]) - (real)rank_beta_0;
      if ( df > real0 )
         {
         pvals[ iSnp ] = PvalFromLikelihoodRatioTest( LL[ iSnp ] - LL_0, (real)0.5 * df );
         }
      else
         {
         //This happens for example if the SNP is a constant. In this case df=0
         pvals[ iSnp ] = std::numeric_limits<real>::quiet_NaN();
         }

      if ( verboseLevel >= verboseDebug )
         {
#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
         const char *szFmt1 = "\n- LL[%Iu]:%f p=%f\n ";
#else                      // Linux/g++ uses a %zu specifier for size_t
         const char *szFmt1 = "\n- LL[%zu]:%f p=%f\n ";
#endif
         printf( szFmt1, iSnp, LL[iSnp], pvals[iSnp] );
         printf( "beta:\n " );
         for( size_t i=0; i<(data->cCovariates+1); ++i )
            {
            size_t idxBeta = ColumnMajorIndex( data->cCovariates + 1, data->cSnps, i, iSnp );
            assert( idxBeta == (iSnp*(data->cCovariates+1)+i) );
            printf( "%f ", beta[ idxBeta ] );
            }
         }
      }
   }

real CLogisticRegression::EvalLL( size_t ind_SNP, bool altModel )
   {
   BLAS_INT n  = (BLAS_INT)data->cIndividuals;
   BLAS_INT d_ = (BLAS_INT)data->cCovariates;
   real *y     = data->phenotypeData;  // point to phenotypeData in CLogisticRegression.data
   real *X     = nullptr;
   BLAS_INT *rank_beta_ = &rank_beta_0;   // rank_beta_0 is lines of covariates in null model
   real *beta_ = beta_0;                  // beta_ points to the location to hold the weight value

   // Set cVars for optimization routine to know how many to work on and point to workspace
   cVars = (int)data->cCovariates;
   //X = pMatrixData.get();
   X = pMatrixData;
   if ( vMatrix.size() != cVars )
   if ( altModel )
      {
      //alternative model includes one SNP
      ++d_;    // add 1 SNP to the covariates
      ++cVars;

      // Copy this SNPs data Create the matrix to operate on
      vector< real >vT;
      vT.resize( data->cIndividuals );
      for ( size_t iIndividual=0; iIndividual<data->cIndividuals; ++iIndividual )
         {
         size_t idxSnpData = ColumnMajorIndex( data->cIndividuals, data->cSnps, iIndividual, ind_SNP );
         vT[ iIndividual ] = data->snpData[ idxSnpData ];

         size_t idxMatrixData = ColumnMajorIndex( data->cIndividuals, data->cCovariates, iIndividual, data->cCovariates );
         X[ idxMatrixData ] = data->snpData[ idxSnpData ];
         }
      vMatrix[ data->cCovariates ] = vT;

      rank_beta_ = &rank_beta[ ind_SNP ];
      beta_ = &beta[ (data->cCovariates+1) * ind_SNP ];
      }

   if ( vWeights.size() != cVars )
      {
      vWeights.resize( cVars );
      }

   std::function< double (vector< double >::iterator itWgts, vector< double >::iterator itGrad ) > prob;
   prob = std::bind( &CLogisticRegression::ProbRoutine, this, std::placeholders::_1, std::placeholders::_2 );
   vector< double > vResults = ts::minimize::LBFGS( prob, static_cast< unsigned int >(cVars), vWeights, lbfgsTolerance, cMaxIter, cLbfgsHistory );

   *rank_beta_ = static_cast< int >( cVars );
   for ( size_t iBeta = 0; iBeta < cVars; ++iBeta )
      {
      beta_[ iBeta ] = vWeights[ iBeta ];
      }

   return( -vResults.back() );   // returns negative of the last element of the results array (minimum 
   }


vector< string > CLogisticRegressionColumnNames( bool fVerboseOutput_, bool fUseBiasTerm, size_t cCovariates )
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
//   columnNames.push_back( "SnpWeightSE" );      // Standard Error
//   columnNames.push_back( "OddsRatio" );        // Odds Ratio = exp( abs(snpWeight) * snpStdDeviation )
//   columnNames.push_back( "WaldStat" );

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

   return( columnNames );
   }

void CLogisticRegression::WriteResultFile( const string& resultFile, bool fVerboseOutput_, real maxPvaluePrintThreshold )
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
   vector< string > columnNames = CLogisticRegressionColumnNames( fVerboseOutput_, useBiasTerm, cCovariates );
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
      real snpWeight = beta[ index ];                                   // SNPWeight -> right from beta vector
      fprintf( pFile, fmtReal, szSeparator, snpWeight );                // SNPWeight - fixed-effect weight of SNP

//      real snpWeightSE = sqrt( (snpWeight * snpWeight) / (f_stat[ index ]) ); // SNPWeightSE -> SNP weight standard error
//      fprintf( pFile, fmtReal, szSeparator, snpWeightSE );              // SNPWeightSE -> SNP weight standard error
//
//      real snpOddsRatio = exp( ::abs( snpWeight ) / snp.StdDeviation() );
//      fprintf( pFile, fmtReal, szSeparator, snpOddsRatio );
//
//      real waldStat = (snpWeight * snpWeight) / (snpWeightSE * snpWeightSE);   // WaldStat = (SnpWeight / stdErr)^2
//      fprintf( pFile, fmtReal, szSeparator, waldStat );

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

#endif   // defined( USE_LINUX_V2x_KERNEL )
