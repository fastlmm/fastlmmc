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

const bool FastLMM::brentExhaustive = false;

FastLMM::~FastLMM()
   {
   // eigen is a 'cached' pointer.  DO NOT DELETE IT HERE
   // if ( eigen != nullptr ) { delete eigen; eigen = nullptr; }
   if (sigg2){    delete[] sigg2;  sigg2 = nullptr; }
   if (sige2){    delete[] sige2;  sige2 = nullptr; }
   if (delta){    delete[] delta;  delta = nullptr; }
   if (beta){     delete[] beta;   beta  = nullptr; }
   if (LL){       delete[] LL;     LL    = nullptr; }
   if (LL_0){     delete[] LL_0;   LL_0  = nullptr; }
   if (rank_beta){delete[] rank_beta; rank_beta = nullptr; }
   if (pvals){    delete[] pvals;   pvals   = nullptr; }
   if (sigg2_0){  delete[] sigg2_0; sigg2_0 = nullptr; }
   if (sige2_0){  delete[] sige2_0; sige2_0 = nullptr; }
   if (delta_0){  delete[] delta_0; delta_0 = nullptr; }
   if (beta_0){   delete[] beta_0;  beta_0  = nullptr; }
   if (f_stat){   delete[] f_stat;  f_stat  = nullptr; }
   if (rank_beta_0){ delete[] rank_beta_0; rank_beta_0 = nullptr; }
   }

FastLMM::FastLMM( EigenSym *eigen_ )
   {
   // initialize all member variables in the constructor to prevent garbage operation
   phenLabel = "NA";
   ncovar   = std::numeric_limits<BLAS_INT>::quiet_NaN();   // number of covariates (cCovariates)
   n        = eigen_->n;                                    // ?
   nSNPs    = std::numeric_limits<BLAS_INT>::quiet_NaN();   // number of SNPs in test data
   nSNPsW   = std::numeric_limits<BLAS_INT>::quiet_NaN();   // number of SNPs in covariate data

   eigen             = eigen_;      // cache pointer to EigenSym (shared pointer)

   excludeSnpStart   = nullptr;     // ?
   nExclude          = nullptr;     // indices of SNPs to exclude
   pvals             = nullptr;     // array of pvalues
   sigg2             = nullptr;
   sige2             = nullptr;
   delta             = nullptr;
   beta              = nullptr;
   LL                = nullptr;     // Log Likelihood
   sigg2_0           = nullptr;
   sige2_0           = nullptr;
   delta_0           = nullptr;
   beta_0            = nullptr;
   LL_0              = nullptr;
   f_stat            = nullptr;
   rank_beta         = nullptr;
   rank_beta_0       = nullptr;

   numberIntervals   = numberOptimizationStartsGlobal;
   minLogVal         = minLogValGlobal;
   maxLogVal         = maxLogValGlobal;
   brentTol          = brentTolGlobal;
   brentMaxIter      = brentMaxIterGlobal;
   }
 
void FastLMM::updateExcludeList()
   {
   if ( nExclude )
      {
      delete[] nExclude;
      nExclude = nullptr;
      }
   if ( excludeSnpStart )
      {
      delete[] excludeSnpStart;
      excludeSnpStart = nullptr;
      }
   nExclude = new size_t[nSNPs];
   excludeSnpStart = new size_t[nSNPs];
   nExclude[0] = 0;
   excludeSnpStart[0] = 0;
   for ( size_t snp = 0 ; snp < this->nSNPs; ++snp )
      {
      if ( snp > 0 )
         {
         nExclude[snp] = nExclude[snp-1];
         excludeSnpStart[snp] = excludeSnpStart[snp-1];
         }
      if ( excludeSnpParams.fileType != NoFile )
         {
         if ( excludeSnpParams.distType == ByID )
            {
            Fatal( "TODO build a map holding all the IDs of the kernel SNPs" );//if(this->snpInfo[snp].idSnp)
            }
         else if ( excludeSnpParams.distType == ByGeneticDistance )
            {
            if ( (snp > 0) 
              && (CPlinkFile::getGeneticDistance( this->snpInfo[snp], this->snpInfo[snp-1] ) < 0.0)
               )
               {
               Fatal( "test SNP file is not ordered ascending in genetic distance" );
               }
            while( (this->excludeSnpStart[snp] < this->nSNPsW)
                && (CPlinkFile::getGeneticDistance( this->snpInfoCovariance[this->excludeSnpStart[snp]], this->snpInfo[snp]) < -excludeSnpParams.geneticDist)
                 )
               {
               ++(this->excludeSnpStart[ snp ]);
               if ( nExclude[snp] > 0 )
                  {
                  --nExclude[ snp ];
                  }
               }
            while( (this->excludeSnpStart[snp] + nExclude[snp] < this->nSNPsW)
                && (CPlinkFile::getGeneticDistance(this->snpInfoCovariance[this->excludeSnpStart[snp] + nExclude[snp]], this->snpInfo[snp]) <= excludeSnpParams.geneticDist)
                 )
               {
               ++nExclude[ snp ];
               }
            }
         else if ( excludeSnpParams.distType == ByPosition )
            {
            if ( (snp>0) && (CPlinkFile::getBpDistance(this->snpInfo[snp], this->snpInfo[snp-1]) < 0) )
               {
               Fatal( "test SNP file is not ordered ascending in BP distance" );
               }
            while( (this->excludeSnpStart[snp] < this->nSNPsW)
                && (CPlinkFile::getBpDistance(this->snpInfoCovariance[this->excludeSnpStart[snp]], this->snpInfo[snp]) < -excludeSnpParams.bpDist)
                 )
               {
               ++(this->excludeSnpStart[ snp ]);
               if ( nExclude[ snp ] > 0 )
                  {
                  --nExclude[ snp ];
                  }
               }
            while( (this->excludeSnpStart[snp] + nExclude[snp] < this->nSNPsW)
                && (CPlinkFile::getBpDistance(this->snpInfoCovariance[this->excludeSnpStart[snp] + nExclude[snp]], this->snpInfo[snp]) <= excludeSnpParams.bpDist)
                 )
               {
               ++nExclude[ snp ];
               }
            }
         }
      }
   }


//solve the least squares problem for precomputed XX and Xy in the likelihood function but try to use specialized functions for small d_.
BLAS_INT FastLMM::solveLstSqSym( BLAS_INT d_, real *XX, real *Xy, real *beta_temp, BLAS_INT *rank_beta, real *f_stat, real *logDetXX )
   {
   if ( d_ > 3 )
      {//we trust the hard-coded least square solver functions up to d_=3.
      real *U = new real[ d_ * d_ ];  //Warning: BLAS_INT product
      real *S = new real[ d_ ];
      EigenSym::EV_symmetric( d_, XX, U, S );
      BlasLapack::gemv( &char_T, &d_, &d_, &real_one, U, &d_, Xy, &int_one, &real_zero, XX, &int_one );
      *rank_beta = 0;

      for ( size_t i = 0; (BLAS_INT)i < d_; ++i )
         {
         f_stat[ i ] = 0.0;
         }
      
      *logDetXX = 0.0;
      for (size_t i = 0; (BLAS_INT)i < d_; ++i )
         {
         if( S[i] > MathFunctions::eps_rank )
            {
            *logDetXX += log( S[i] );
            ++*rank_beta;
            XX[ i ] /= S[ i ];
            for( size_t j = 0; (BLAS_INT)j < d_; ++j )
               {
               f_stat[ j ] += U[ (i * d_) + j ] * U[ (i * d_) + j ] / S[i];   //Warning: size_t BLAS_INT conversion product
               }
            }
         else
            {
            //printf("\n lowrank S[%i] = %f", i, S[i]);
            XX[ i ] = 0.0; 
            *logDetXX = -std::numeric_limits<real>::infinity();
            }
         }
      
      BlasLapack::gemv( &char_N, &d_, &d_, &real_one, U, &d_, XX, &int_one, &real_zero, beta_temp, &int_one );
      delete[] U;
      delete[] S;
      for (size_t i = 0; (BLAS_INT)i < d_; ++i )
         {
         f_stat[i] = beta_temp[i] * beta_temp[i] / f_stat[i];
         }
      }  //end d>3
   else if ( d_ == 0 )
      {
      *rank_beta = 0;
      *logDetXX = -std::numeric_limits<real>::infinity();
      }  //end d_==0
   else if ( d_ == 1 )
      {
      if ( XX[0] < MathFunctions::eps_rank )
         {//X is constant => after zero-mean it is 0.0
         *rank_beta = 0;   //res stays constant.
         beta_temp[0] = 0.0;
         f_stat[0] = 0.0;
         *logDetXX = -std::numeric_limits<real>::infinity();
         }
      else
         {
         *rank_beta = 1;
         beta_temp[0] = Xy[0] / XX[0];
         f_stat[0] = Xy[0] * Xy[0] / XX[0];
         *logDetXX =log( XX[ 0 ] );
         }
      }  //end d_==1
   else if ( d_ == 2 )
      {
      //BlasLapack::gemm(&char_T, &char_N, &d_, &d_, &n_, &real_one, XSdi, &n_, XSdi, &n_, &real_zero, XX, &d_);
      *rank_beta=MathFunctions::leastSquares2by2SymLower( XX, Xy, beta_temp, f_stat, logDetXX );
      }
   else if ( d_ == 3 )
      {
      *rank_beta = MathFunctions::solve3by3symmL( XX, Xy, beta_temp, f_stat );
      if ( *rank_beta == 0 )
         {//in the rank deficient case, call the slow standard function. TODO: code up 3-by-3 low rank solver.
         real *eigenvec= new real[9];
         real *eigenval= new real[d_];
         EigenSym::EV_symmetric( d_, XX, eigenvec, eigenval );
         BlasLapack::gemv( &char_T, &d_ ,&d_ ,&real_one, eigenvec, &d_, Xy, &int_one, &real_zero, XX, &int_one );
         *logDetXX = 0.0;
         if( eigenval[0] > MathFunctions::eps_rank)
            {
            ++*rank_beta;
            XX[0] /= eigenval[0];
            f_stat[0] = eigenvec[0] * eigenvec[0] / eigenval[0];
            f_stat[1] = eigenvec[1] * eigenvec[1] / eigenval[0];
            f_stat[2] = eigenvec[2] * eigenvec[2] / eigenval[0];
            *logDetXX += log( eigenval[0] );
            }
         else
            {
            XX[ 0 ] = 0.0;
            f_stat[ 0 ] = 0.0;
            f_stat[ 1 ] = 0.0;
            f_stat[ 2 ] = 0.0;
            *logDetXX = -std::numeric_limits<real>::infinity();
            }
         if( eigenval[1] > MathFunctions::eps_rank)
            {
            ++*rank_beta;
            XX[ 1 ] /= eigenval[ 1 ];
            f_stat[ 0 ] += eigenvec[ 3 ] * eigenvec[ 3 ] / eigenval[ 1 ];
            f_stat[ 1 ] += eigenvec[ 4 ] * eigenvec[ 4 ] / eigenval[ 1 ];
            f_stat[ 2 ] += eigenvec[ 5 ] * eigenvec[ 5 ] / eigenval[ 1 ];
            *logDetXX += log( eigenval[ 0 ] );
            }
         else
            {
            XX[ 1 ] = 0.0;
            *logDetXX = -std::numeric_limits<real>::infinity();
            }
         if( eigenval[ 2 ] > MathFunctions::eps_rank )
            {
            ++*rank_beta;
            XX[ 2 ] /= eigenval[ 2 ];
            f_stat[ 0 ] += eigenvec[ 6 ] * eigenvec[ 6 ] / eigenval[ 2 ];
            f_stat[ 1 ] += eigenvec[ 7 ] * eigenvec[ 7 ] / eigenval[ 2 ];
            f_stat[ 2 ] += eigenvec[ 8 ] * eigenvec[ 8 ] / eigenval[ 2 ];
            *logDetXX += log( eigenval[ 0 ] );
            }
         else
            {
            XX[ 2 ] = 0.0;
            *logDetXX = -std::numeric_limits<real>::infinity();
            }
         BlasLapack::gemv( &char_N, &d_, &d_, &real_one, eigenvec, &d_, XX, &int_one, &real_zero, beta_temp, &int_one );
         delete[] eigenvec;
         delete[] eigenval;
         f_stat[ 0 ] = beta_temp[ 0 ] * beta_temp[ 0 ] / f_stat[ 0 ];
         f_stat[ 1 ] = beta_temp[ 1 ] * beta_temp[ 1 ] / f_stat[ 1 ];
         f_stat[ 2 ] = beta_temp[ 2 ] * beta_temp[ 2 ] / f_stat[ 2 ];
         //FastLMM::dgelsy( n, d_, XSdi_leastsquares, beta_temp, rank_beta, XSdi, res );
         }
      }
   return( 0 );
   }

/*
 * FastLMM::getResults() copies the results from the FastLMM
 */
void FastLMM::getResults( ResultsLMM& resultsLMM )
   {
   resultsLMM.PrepareForResults( snpInfo.size(), ncovar );  // ensure proper memory allocations

   //resultsLMM.modelID       = modelID;    // ResultsLMM does not have a modelID member
                                          // TODO: WORKITEM change to use a pointer avoid the copy?
   resultsLMM.snpInfo         = snpInfo;  // TODO: This COPIES the vector<SnpInfo> to a new vector!
   resultsLMM.phenLabel       = phenLabel;
   resultsLMM.ncovariates     = ncovar;
   resultsLMM.n_individuals   = n;

   Qval::computeQvals( pvals, resultsLMM.qvals, resultsLMM.index, nSNPs );

   for ( size_t iSnp = 0; iSnp < nSNPs; ++iSnp )
      {
      resultsLMM.LL_alt[ iSnp ]    = LL[ iSnp ];
      resultsLMM.LL_0[ iSnp ]      = LL_0[ iSnp ];
      resultsLMM.pvals[ iSnp ]     = pvals[ iSnp ];
      resultsLMM.snpWeights_alt[ iSnp ] = beta[ (1 + ncovar) * iSnp ];
      resultsLMM.bias_alt[ iSnp ]  = beta[ ((1 + ncovar) * iSnp) + ncovar ];
      resultsLMM.sigg2_alt[ iSnp ] = sigg2[ iSnp ];
      resultsLMM.sige2_alt[ iSnp ] = sige2[ iSnp ];
      resultsLMM.delta_alt[ iSnp ] = delta[ iSnp ];
      resultsLMM.bias_0[ iSnp ]    = beta_0[ (ncovar * iSnp) + (ncovar - 1) ];
      resultsLMM.sige2_0[ iSnp ]   = sige2_0[ iSnp ];
      resultsLMM.sigg2_0[ iSnp ]   = sigg2_0[ iSnp ];
      resultsLMM.delta_0[ iSnp ]   = delta_0[ iSnp ];
      resultsLMM.f_stat[ iSnp ]    = f_stat[ (1 + ncovar) * iSnp ];
      resultsLMM.dof[ iSnp ]       = rank_beta[ iSnp ] - rank_beta_0[ iSnp ];

      if ( excludeSnpParams.fileType != NoFile )
         {
         resultsLMM.nExclude[ iSnp ] = nExclude[ iSnp ];
         resultsLMM.excludeSnpStart[ iSnp ] = excludeSnpStart[ iSnp ];
         }
      }

   // beta_res has one additional rank compared to beta_0 
   //  due to the SNP weight so initialize both arrays
   //  with the proper data
   size_t cbT1 = (nSNPs * ncovar * sizeof( real ));
   memmove( resultsLMM.beta_0, beta_0, cbT1 );

   size_t cbT2 = (nSNPs * (1+ncovar) * sizeof( real ));
   memmove( resultsLMM.beta, beta, cbT2 );
   }

void FastLMM::getResultsNullOnly( ResultsLMM& resultsLMM )
   {
   resultsLMM.Clear();
   //resultsLMM.modelID       = modelID;    // ResultsLMM does not have a modelID member
                                          // TODO: WORKITEM change to use a pointer avoid the copy?
   resultsLMM.snpInfo         = snpInfo;  // TODO: This COPIES the vector<SnpInfo> to a new vector!
   resultsLMM.phenLabel       = phenLabel;
   resultsLMM.ncovariates     = ncovar;
   resultsLMM.n_individuals   = n;

   resultsLMM.LL_0      = new real[ 1 ];
   resultsLMM.bias_0    = new real[ 1 ];
   resultsLMM.sige2_0   = new real[ 1 ];
   resultsLMM.sigg2_0   = new real[ 1 ];
   resultsLMM.delta_0   = new real[ 1 ];
   resultsLMM.beta_0    = new real[ ncovar ];
   
   resultsLMM.LL_0[ 0 ]    = LL_0[ 0 ];
   resultsLMM.bias_0[ 0 ]  = beta_0[ (ncovar * 0) + (ncovar - 1) ];
   resultsLMM.sige2_0[ 0 ] = sige2_0[ 0 ];
   resultsLMM.sigg2_0[ 0 ] = sigg2_0[ 0 ];
   resultsLMM.delta_0[ 0 ] = delta_0[ 0 ];
   
   size_t cb = (ncovar * sizeof( real ));
   memmove( resultsLMM.beta_0, beta_0, cb );
   }

/* **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** */
// ResultsLMM Helpers here
/* **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** */
FILE *OpenResultFile( const string& resultFile )
   {
   if ( resultFile.size() == 0 )
      {
      Fatal( "No output filename to open " );
      }

   FILE *pFile = fopen( resultFile.c_str(), "wt" );          // write ascii
   if ( !pFile )
      {
      string fullPath = FullPath( resultFile );
      Fatal( "Cannot open output file [%s].  \n  CRT Error %d: %s", fullPath.c_str(), errno, strerror( errno ) );
      }

   return( pFile );
   }

void SelectFormatStrings_V2( const char*& fmtInt_, const char*& fmtString_, const char*& fmtReal_ )
   {
   static char fmtFastLmmRealsPrintPrecision[16];  // static storage for real format string

   fmtInt_    = "%s%i";       // use the %s so we can create CSV or TAB separated output
   fmtString_ = "%s%s";

   // command line arg -SetOutputPrecision <int> takes 
   //   precedence over environment setting.
   if ( !fSetFastLmmRealsPrintPrecision ) //check for command line argument
      {
      // if no command line argument.  Check for environment variable set to
      // adjust the real output precision
      char *precisionEnvValue = getenv( "FastLmmRealsPrintPrecision" );
      if ( precisionEnvValue != nullptr )
         {
         // found one.  Setup for the new precision
         fastLmmRealsPrintPrecision = atoi( precisionEnvValue );

         if ( fastLmmRealsPrintPrecision > 18 ) fastLmmRealsPrintPrecision = 18;
         if ( fastLmmRealsPrintPrecision <  3 ) fastLmmRealsPrintPrecision =  3;
         }
      }
   sprintf( fmtFastLmmRealsPrintPrecision, "%%s%%.%dE", fastLmmRealsPrintPrecision );
   fmtReal_ = fmtFastLmmRealsPrintPrecision;
   }

vector< string > ResultColumnNames( bool fVerboseOutput_, bool fExcludeSnps_, size_t firstCovariate, size_t cCovariates )
   {
   vector< string > columnNames;
   columnNames.reserve( 32 );                      // pre-allocate enough for 'most' cases.

   // Build vector of column names
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

   if ( fExcludeSnps_ )
      {
      columnNames.push_back( "NumSNPsExcluded" );        // excludebygeneticdistance
      columnNames.push_back( "IndexExclusionStart" );    // excludebygeneticdistance
      }

   if ( fVerboseOutput_ )
      {
      columnNames.push_back( "DOF" );
      }

   columnNames.push_back( "NullLogLike" );
   columnNames.push_back( "AltLogLike" );
   columnNames.push_back( "SNPWeight" );
   columnNames.push_back( "SNPWeightSE" );
   columnNames.push_back( "OddsRatio" );
   columnNames.push_back( "WaldStat" );
   if ( useHeritability )
      {
      columnNames.push_back( "heritability" );
      }
   else
      {
      columnNames.push_back( "NullLogDelta" );
      }
   columnNames.push_back( "NullGeneticVar" );
   columnNames.push_back( "NullResidualVar" );
   columnNames.push_back( "NullBias" );

   // TODO:  What names really go here?
   for( size_t iCovariate = 0; (firstCovariate + iCovariate) < cCovariates; ++iCovariate )
      {
      // iterate to generate 'NullCov%02dWeight' column names
      char szNullCovWeight[32];
      sprintf( szNullCovWeight, "NullCov%02dWeight", (int)(firstCovariate + iCovariate) );
      columnNames.push_back( szNullCovWeight );
      }

   if ( fVerboseOutput_ )
      {
      if ( useHeritability )
         {
         columnNames.push_back( "heritability" );        // verbose
         }
      else
         {
         columnNames.push_back( "LogDelta" );            // verbose
         }
      columnNames.push_back( "GeneticVar" );             // verbose
      columnNames.push_back( "ResidualVar" );            // verbose
      columnNames.push_back( "Bias" );                   // verbose

      // TODO:  What names really go here?
      for( size_t iCovariate = 0; (firstCovariate + iCovariate) < cCovariates; ++iCovariate )
         {
         // iterate to generate 'Cov%02dWeight' column names
         char szCovWeight[32];
         sprintf( szCovWeight, "Cov%02dWeight", (int)(firstCovariate + iCovariate) );
         columnNames.push_back( szCovWeight );
         }

      columnNames.push_back( "SNPIndex" );               // verbose
      columnNames.push_back( "SNPCount" );               // verbose
      }

   return( columnNames );
   }

void OutputColumnNames( FILE *pFile, const vector< string >& names, const char* separator )
   {
   fprintf( pFile, "%s", names[0].c_str() );
   for ( size_t iName = 1; iName < names.size(); ++iName )
      {
      fprintf( pFile, "%s%s", separator, names[ iName ].c_str() );
      }

   fprintf( pFile, "\n" );
   }

/* **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** */
// ResultsLMM method / implementation routines here
/* **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** */
ResultsLMM::ResultsLMM()
   {
   ncovariates    = 0;        // The number of covariates in the study
   n_individuals  = 0;        // The number of individuals in the study
   index          = nullptr;  // The results will be printed out, sorted by index
   nExclude       = nullptr;  // 
   excludeSnpStart= nullptr;  // 
   LL_alt         = nullptr;  // Log-likelihoods of the alternative models
   LL_0           = nullptr;  // Log-likelihoods of the null models
   pvals          = nullptr;  // p-values
   snpWeights_alt = nullptr;  // the SNP weights of the alternative models
   bias_alt       = nullptr;  // the bias weights of the alternative models
   sigg2_alt      = nullptr;  // the genetic variance of the alternative models
   sige2_alt      = nullptr;  // the random vriance of the alternative models
   delta_alt      = nullptr;  // the fraction of sige2_alt/sigg2_alt
   bias_0         = nullptr;  // the bias terms of the null models
   sigg2_0        = nullptr;  // the genetic variance of the null models
   sige2_0        = nullptr;  // the random variance of the null models
   delta_0        = nullptr;  // fraction of sige2_0/sigg2_0
   qvals          = nullptr;  // q-values (minimum FDR, under which the current SNP would be significant)
   f_stat         = nullptr;  // f-statistics of the SNPs
   dof            = nullptr;  // degrees of freedom
   beta_0         = nullptr;  // array of Null Covariate Weights
   beta           = nullptr;  // array of Alternative Covariate Weights

   // Initialize the Iterative Result fields
   fIterativeResults     = false; // signify NOT iterative by default
   iResultBlock          = 0;       // default to zero'th block
   cPriorResults         = 0;       // default to zero prior results

   // Initialize the support for formatting results & statistics
   pFileResults          = nullptr;
   fProduceVerboseOutput = false;   // print the 'extra' columns for each SNP too
   szSeparator           = nullptr; // cache pointers used
   fmtInt                = nullptr; 
   fmtString             = nullptr; 
   fmtReal               = nullptr; 
   fExcludeSnps          = false;   // 
   firstCovariate        = 0;       // if biasTerm is used, first covariate is bias so skip it
   cSnpsFailedPvaluePrintThreshold = 0;
   pvaluePrintThresholdMax = 1.0;   // do not print pvals greater than pvaluePrintThreshold
   }

/*
 * ResultsLMM::~ResultsLMM() - destructor must release raw allocated 
 *                               memory and resources
 */
ResultsLMM::~ResultsLMM()
   {
   Clear();       // Release all the memory held by raw memory pointers in ResultsLMM
   if ( pFileResults )
      {
      CloseResultsFile();
      }
   }

void ResultsLMM::WriteResultFile( const string& resultFile )
   {
   WriteResultFile_V2( resultFile, fVerboseOutput, pvaluePrintThreshold );
   }

/*
 * ResultsLMM:Clear() Releases all of the internally allocated raw memory fields
 */
void ResultsLMM::Clear()
   {
   if ( index           != nullptr ) { delete[] index           ; index           = nullptr; }
   if ( nExclude        != nullptr ) { delete[] nExclude        ; nExclude        = nullptr; }
   if ( excludeSnpStart != nullptr ) { delete[] excludeSnpStart ; excludeSnpStart = nullptr; }
   if ( LL_alt          != nullptr ) { delete[] LL_alt          ; LL_alt          = nullptr; }
   if ( LL_0            != nullptr ) { delete[] LL_0            ; LL_0            = nullptr; }
   if ( pvals           != nullptr ) { delete[] pvals           ; pvals           = nullptr; }
   if ( snpWeights_alt  != nullptr ) { delete[] snpWeights_alt  ; snpWeights_alt  = nullptr; }
   if ( bias_alt        != nullptr ) { delete[] bias_alt        ; bias_alt        = nullptr; }
   if ( sigg2_alt       != nullptr ) { delete[] sigg2_alt       ; sigg2_alt       = nullptr; }
   if ( sige2_alt       != nullptr ) { delete[] sige2_alt       ; sige2_alt       = nullptr; }
   if ( delta_alt       != nullptr ) { delete[] delta_alt       ; delta_alt       = nullptr; }
   if ( bias_0          != nullptr ) { delete[] bias_0          ; bias_0          = nullptr; }
   if ( sigg2_0         != nullptr ) { delete[] sigg2_0         ; sigg2_0         = nullptr; }
   if ( sige2_0         != nullptr ) { delete[] sige2_0         ; sige2_0         = nullptr; }
   if ( delta_0         != nullptr ) { delete[] delta_0         ; delta_0         = nullptr; }
   if ( qvals           != nullptr ) { delete[] qvals           ; qvals           = nullptr; }
   if ( f_stat          != nullptr ) { delete[] f_stat          ; f_stat          = nullptr; }
   if ( dof             != nullptr ) { delete[] dof             ; dof             = nullptr; }
   if ( beta_0          != nullptr ) { delete[] beta_0          ; beta_0          = nullptr; }
   if ( beta            != nullptr ) { delete[] beta            ; beta            = nullptr; }
   }

// Allocate the memory for these results
void ResultsLMM::PrepareForResults( size_t cSnps, size_t cCovariates )
   {
   //TODO: WORKITEM  We might be able to reuse memory if we have
   //                  the same cSnp & cCovariate values.  
   //                  Confirm later and fix.
   Clear();          // deallocates all previously allocated memory

   snpInfo.clear();
   phenLabel.clear();
   ncovariates     = 0;
   n_individuals   = 0;
   index           = new size_t[ cSnps ];  // The results will be printed out, sorted by index
   nExclude        = new size_t[ cSnps ];  // 
   excludeSnpStart = new size_t[ cSnps ];  // 
   LL_alt          = new real[ cSnps ];    // Log-likelihoods of the alternative models
   LL_0            = new real[ cSnps ];    // Log-likelihoods of the null models
   pvals           = new real[ cSnps ];    //
   snpWeights_alt  = new real[ cSnps ];    //
   bias_alt        = new real[ cSnps ];    //
   sigg2_alt       = new real[ cSnps ];    //
   sige2_alt       = new real[ cSnps ];    //
   delta_alt       = new real[ cSnps ];    //
   bias_0          = new real[ cSnps ];    //
   sigg2_0         = new real[ cSnps ];    //
   sige2_0         = new real[ cSnps ];    //
   delta_0         = new real[ cSnps ];    //
   qvals           = new real[ cSnps ];    //
   f_stat          = new real[ cSnps ];    // 
   dof             = new BLAS_INT[ cSnps ];   // degrees of freedom
   beta_0          = new real[ cSnps * cCovariates ];    // Null Covariate Weights
   beta            = new real[ cSnps * (cCovariates+1) ];   // + room for SNP weight
   }

/*
 * ResultsLMM::OpenResultsFile()
 */
void ResultsLMM::OpenResultsFile( const string& resultFile_, bool fWriteVerboseOutput_, real pvaluePrintThreshold, size_t cCovariates_ )
   {
   if ( resultFile_.size() == 0 )
      {
      Fatal( "No output filename to open " );
      }

   resultsFilename = resultFile_;
   pFileResults = fopen( resultsFilename.c_str(), "wt" );   // we create and write an ascii file
   if ( !pFileResults )
      {
      string fullPath = FullPath( resultsFilename );
      Fatal( "Cannot open output file [%s].  \n  CRT Error %d: %s", fullPath.c_str(), errno, strerror( errno ) );
      }

   // cache our results formatting information
   szSeparator = FIsCsvFile( resultsFilename ) ? "," : "\t";   // if .CSV file, separator is <comma>, else <tab>
   SelectFormatStrings_V2( fmtInt, fmtString, fmtReal );
   fProduceVerboseOutput = fWriteVerboseOutput_;
   pvaluePrintThresholdMax = pvaluePrintThreshold;
   fExcludeSnps = FExcludeSnps();                  // going to global variables.
   firstCovariate = useBiasTerm ? 1 : 0;           // if biasTerm is used, first covariate is bias so skip it
   cCovariates = cCovariates_;

   vector< string > columnNames = ResultColumnNames( fWriteVerboseOutput_, fExcludeSnps, firstCovariate, cCovariates );
   OutputColumnNames( pFileResults, columnNames, szSeparator );
   }

/*
 * ResultsLMM::CloseResultsFile() - release the file handle
 */
void ResultsLMM::CloseResultsFile()
   {
   if ( pFileResults )
      {
      fclose( pFileResults );
      pFileResults = nullptr;
      }
   }

void ResultsLMM::WriteResultsBody()
   {
   size_t cSnps = snpInfo.size();

   for( size_t i=0; i<cSnps; ++i )
      {
      size_t iSnp = ( fEmitSorted && runAlternativeModels ) ? index[ i ] : i;
      real pvalT = pvals[ iSnp ];
      if ( pvalT > pvaluePrintThresholdMax )
         {
         ++cSnpsFailedPvaluePrintThreshold;
         continue;      // skip this SNP
         }

      SnpInfo& snp = snpInfo[ iSnp ];
      fprintf( pFileResults, "%s",                 snp.idSnp.c_str() );        // SNP - the SNP label (no leading separator)
      fprintf( pFileResults, fmtInt,  szSeparator, snp.iChromosome );          // Chromosome - the chromosome number
      fprintf( pFileResults, "%s%f",  szSeparator, snp.geneticDistance );      // GeneticDistance - (morgans)
      fprintf( pFileResults, fmtInt,  szSeparator, snp.basepairPosition );     // Position - 
      if ( fProduceVerboseOutput )
         {
         fprintf( pFileResults, fmtString, szSeparator, phenLabel.c_str() );   // Phenotype - Name
         }

      // REVIEW: in old WriteResulteFile() runAlternativeModels was used
      if ( runAlternativeModels )
         {
         if ( _isnan( pvalT ) )
            {
            // this is a NaN!
            Warn( "Found a NaN pvalue for SNP [%s] at input row %d", snp.idSnp.c_str(), iSnp );
            fprintf( pFileResults, fmtString, szSeparator, "NA" );             // NaN ?
            }
         else
            {
            // valid pvalue - output the value
            fprintf( pFileResults, fmtReal, szSeparator, pvalT );
            }
         if ( fIterativeResults )
            {
            fprintf( pFileResults, fmtString, szSeparator, "" );     // blank field (for now)
            }
         else if ( _isnan( qvals[ iSnp ] ) )
            {
            // this is a NaN!
            Warn( "Found a NaN qvalue for SNP [%s] at row %d", snp.idSnp.c_str(), i );
            fprintf( pFileResults, fmtString, szSeparator, "NA" );    // NaN ?
            }
         else
            {
            // Qvalue - output proper qvalues
            fprintf( pFileResults, fmtReal, szSeparator, qvals[ iSnp ] );
            }
         }

      fprintf( pFileResults, fmtInt, szSeparator, n_individuals );                // N - samples size used in analysis (cIndividuals)
      if ( fExcludeSnps )
         {
         fprintf( pFileResults, fmtInt, szSeparator, nExclude[ iSnp ] );          // NumSNPsExcluded - 
         fprintf( pFileResults, fmtInt, szSeparator, excludeSnpStart[ iSnp ] );   // IndexExclusionStart - 
         }

      // REVIEW: in old WriteResulteFile() runAlternativeModels was used as condition
      if ( fProduceVerboseOutput )
         {
         fprintf( pFileResults, fmtInt, szSeparator, dof[ iSnp ] );               // DOF - Degrees of Freedom
         }

      if ( runAlternativeModels )
         {
         fprintf( pFileResults, fmtReal, szSeparator, LL_0[ iSnp ] );             // NullLogLike - likelihood of null model
         fprintf( pFileResults, fmtReal, szSeparator, LL_alt[ iSnp ] );           // AltLogLike - likelihood of alternative model
         fprintf( pFileResults, fmtReal, szSeparator, snpWeights_alt[ iSnp ] );   // SNPWeight - fixed-effect weight of SNP

         real snpWeightSE = sqrt( (snpWeights_alt[ iSnp ] * snpWeights_alt[ iSnp ]) / (f_stat[ iSnp ]) );
         fprintf( pFileResults, fmtReal, szSeparator, snpWeightSE );              // SNPWeightSE -> SNP weight standard error

         //real snpOddsRatio = exp( ::abs( snpWeights_alt[ iSnp ] ) / snp.StdDeviation() );//CL: This line seems to be buggy
         real snpOddsRatio;
         if ( unitVarianceTestSNPs )
            {
            snpOddsRatio = exp( snpWeights_alt[ iSnp ]  / snp.StdDeviation() );
            }
         else
            {
            snpOddsRatio = exp( snpWeights_alt[ iSnp ] );
            }
         fprintf( pFileResults, fmtReal, szSeparator, snpOddsRatio );             // snpOddsRatio -> Odds Ratio = exp( abs(snpWeight) * snpStdDeviation )

         real waldStat = (snpWeights_alt[ iSnp ]*snpWeights_alt[ iSnp ] / (snpWeightSE*snpWeightSE)); // WaldStat = (SnpWeight / stderr)^2
         fprintf( pFileResults, fmtReal, szSeparator, waldStat );

         fprintf( pFileResults, fmtReal, szSeparator, delta_0[ iSnp ] );          // NullLogDelta - / or heritability
         fprintf( pFileResults, fmtReal, szSeparator, sigg2_0[ iSnp ] );          // NullGeneticVar - 
         fprintf( pFileResults, fmtReal, szSeparator, sige2_0[ iSnp ] );          // NullResidualVar - 
         fprintf( pFileResults, fmtReal, szSeparator, bias_0[ iSnp ] );           // NullBias - 

         for( size_t iCovariate = 0; (firstCovariate + iCovariate) < cCovariates; ++iCovariate )
            {
            size_t idx = ColumnMajorIndex( cCovariates, cSnps, (firstCovariate + iCovariate), iSnp );
            fprintf( pFileResults, fmtReal, szSeparator, beta_0[ idx ] );         // NullCov<xx>Weight
            }

         if ( fProduceVerboseOutput )
            {
            fprintf( pFileResults, fmtReal, szSeparator, delta_alt[ iSnp ] );     // LogDelta - / or heritability
            fprintf( pFileResults, fmtReal, szSeparator, sigg2_alt[ iSnp ] );     // GeneticVar - 
            fprintf( pFileResults, fmtReal, szSeparator, sige2_alt[ iSnp ] );     // ResidualVar - 
            fprintf( pFileResults, fmtReal, szSeparator, bias_alt[ iSnp ] );      // Bias - 
            
            for( size_t iCovariate = 0; (firstCovariate + iCovariate) < cCovariates; ++iCovariate )
               {
               size_t idx = ColumnMajorIndex( cCovariates + 1, cSnps, (firstCovariate + iCovariate), iSnp );
               fprintf( pFileResults, fmtReal, szSeparator, beta[ idx ] );        //  Cov<xx>Weight - 
               }
            
            fprintf( pFileResults, fmtInt, szSeparator, (cPriorResults + iSnp + 1) );  // SNPIndex - display as 1 based
            fprintf( pFileResults, fmtInt, szSeparator, (cPriorResults + cSnps) );     // SNPCount - 
            }
         }
      
      fprintf( pFileResults, "\n" );
      }
   }


/*
 * ResultsLMM::WriteResultsSummary() - If the SNP threshold is set, print the summary
 *                                       saying how many SNPs where NOT printed too.
 */
void ResultsLMM::WriteResultsSummary()
   {
   if ( (pvaluePrintThresholdMax != 1.0) && (cSnpsFailedPvaluePrintThreshold > 0) )
      {
      // VC and GCC don't agree on printf() size_t type specifier
      //   convert to unsigned long so both will print
      const char *szFmt = "\n"
                          "\nCount of SNPs that passed print threshold: %lu" 
                          "\nCount of SNPs that failed print threshold: %lu"
                          "\n";
      fprintf( pFileResults, szFmt, (unsigned long)(snpInfo.size() - cSnpsFailedPvaluePrintThreshold), (unsigned long)cSnpsFailedPvaluePrintThreshold );
      }
   }

/*
 *  ResultsLMM::WriteResultFile_V2 - 
 *    Output a tab or comma separated file with the results of the FastLmm 
 *    compuatation. The first row contains column headers and each subsequent 
 *    row contains the information for a specific SNP.
 */
void ResultsLMM::WriteResultFile_V2( const string& resultFile, bool fWriteVerboseOutput_, real maxPvaluePrintThreshold )
   {
   OpenResultsFile( resultFile, fWriteVerboseOutput_, maxPvaluePrintThreshold, ncovariates );
   WriteResultsBody();
   WriteResultsSummary();
   CloseResultsFile();
   }

