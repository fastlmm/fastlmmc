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

/*
 * Parameters: 
 * SNPs [nIndividuals by nSNPs]:
 *                       Matrix stored in column-major order. 
 *                       This will hold the result.
 *                       NaNs will be set to 0.0 in the result.
 * unitVariance: 
 *                true:  Apply unit-variance to the SNPs, 
 *                false: do not unit-variance the SNPs (they are between 0 and 2 - mean)
 * zeroMean:
 *                true:  zero mean the SNPs
 *                false: do not zero Mean the SNPs (leaves them always positive)
 * normalizationInfo
 *                normalizationType:HardyWeinberg: 
 *                       assume SNPs are encoded as 0,1,2 and unitVaraince 
 *                       under Hardy-Weinberg model, 
 *                normalizationType:Beta:
 *                       normalize using SNP /= Beta(MAF, a, b )
 *                normalizationType:Normal:
 *                       use empirical estimator of standard deviation
 * useUnbiasedEstimatorOrBetaBinomial: 
 *                true:  under HardyWeinberg, this add both one positive and 
 *                       one negative pseudo individual in the Beta-Binomial 
 *                       model for the SNP frequency. for empirical std this 
 *                       will use the unbiased estimator.
 *                false: No pseudocounts/biased estimator for std
 */
void ImputeAndZeroMeanSNPs( 
   const vector<SnpInfo>& snpInfo, 
   real *SNPs,                // Array of nIndividuals per SNP
   const size_t nIndividuals, // 
   const size_t nSNPs, 
   const bool fPerformUnitVariance,
   const bool fPerformZeroMean,
   const NormalizationInfo& normalizationInfo_,
   const bool fUseUnbiasedEstimatorOrBetaBinomial
   )
   {
   //pseudo-counts for the Beta prior on the SNP frequency
   const real pseudoCountDenominator = 1.0;
   const real pseudoCountNumerator   = 1.0;
   
   for ( size_t iSnp = 0; iSnp < nSNPs; ++iSnp, SNPs += nIndividuals )  // advance SNPs array pointer too
      {
      const string& idSnp = snpInfo[ iSnp ].idSnp;
      real n_observed = 0.0;
      real sum_s  = 0.0;      //the sum of a SNP over all observed individuals
      real sum2_s = 0.0;      //the sum of the squares of the SNP over all observed individuals
      for( size_t iIndividual = 0; iIndividual < nIndividuals; ++iIndividual )
         {
         if ( !_isnan( SNPs[ iIndividual ] ) )
            {
            // We have 'real' number, not a NaN, continue checking and accumulate sums
            if ( (normalizationInfo.normalizationType == HardyWeinberg)
             &&  (SNPs[iIndividual] != 0.0)
             &&  (SNPs[iIndividual] != 1.0)
             &&  (SNPs[iIndividual] != 2.0)
               )
               {
               Fatal( "Under the Hardy-Weinberg model, SNPs are not allowed to have values other than 0,1,2. SNPs[%i][%s] = %.12f", iIndividual, idSnp.c_str(), SNPs[iIndividual] );
               }

            sum_s += SNPs[ iIndividual ];
            sum2_s+= SNPs[ iIndividual ] * SNPs[ iIndividual ];
            ++n_observed;
            }
         }

      if ( n_observed < 1.0 )
         {
         Fatal( "No individual observed for the SNP[%s].", idSnp.c_str() );
         }

      real mean_s  = sum_s  / n_observed;    //compute the mean over observed individuals for the current SNP
      real mean2_s = sum2_s / n_observed;    //compute the mean of the squared SNP

      real std  = 0.0;                       // the standard deviation
      real freq = 0.0;                       // The SNP frequency

      if ( fUseUnbiasedEstimatorOrBetaBinomial )
         {
         freq = (sum_s + pseudoCountDenominator) / ((n_observed + pseudoCountNumerator) * (real)2.0);
         }
      else
         {
         freq = (sum_s) / (n_observed * (real)2.0);   //compute snp-freq as in the Visscher Height paper (Nat Gen, Yang et al 2010).
         }

      if ( _isnan( freq ) || (freq >= (real)1.0) || (freq <= (real)0.0) )
         {
         Warn( "Illegal SNP frequency: %.2f for SNPs[:][%s]", freq, idSnp.c_str() );
         }

      if ( (normalizationInfo_.normalizationType == Normal) || (normalizationInfo_.normalizationType == Beta) )
         {
         real variance = mean2_s-mean_s * mean_s;        //By the Cauchy Shwartz inequality this should always be positive
         if ( fUseUnbiasedEstimatorOrBetaBinomial && (n_observed > 1.0) )
            {
            if ( n_observed < 2.0 )
               {
               Fatal( "only %i values observed for SNPs[:][%s]. Can not use unbiased variance estimator.", n_observed, idSnp.c_str() );
               }

            variance *= ( n_observed / (n_observed-(real)1.0) );  //replace biased variance estimtor by the unbiased one.
            }

         std = sqrt( variance );
         }
      else if ( normalizationInfo_.normalizationType == HardyWeinberg )
         {
         std = sqrt( (real)2.0 * freq * ((real)1.0-freq));  //compute std under Hardy-Weinberg, given by sqrt(p_i*(1-p_i)) 
         }
      else
         {
         Fatal( "Unknown normalization type.  This is a program error.  Found %d", normalizationInfo_.normalizationType );
         }

      if ( _isnan( std ) || (std <= (real)0.0) )
         {
         // a std == 0.0 means all SNPs have the same value (no variation or Single Nucleotide Constant (SNC))
         //   however ALL SNCs should have been removed in previous filtering steps
         //   This test now prevents a divide by zero error below
         Fatal( "std=.%2f has illegal value for SNPs[:][%s]", std, idSnp.c_str() );
         }

      for ( size_t iIndividual=0; iIndividual < nIndividuals; ++iIndividual )
         {
         //check for NaN
         if ( _isnan( SNPs[ iIndividual ] ) )
            {
            // In the NaN case
            if ( (normalizationInfo_.normalizationType == Normal) || (normalizationInfo_.normalizationType == Beta) )
               {
               // normalizationType == Normal
               if ( fPerformZeroMean )
                  {
                  SNPs[ iIndividual ] = 0.0;
                  }
               else
                  {
                  SNPs[ iIndividual ] = mean_s;
                  }
               }
            else if ( normalizationInfo_.normalizationType == HardyWeinberg )
               {
               // normalizationType == HardyWeinberg
               if ( fPerformZeroMean )
                  {
                  SNPs[ iIndividual ] = mean_s - (2.0 * freq);
                  }
               else
                  {
                  SNPs[ iIndividual ] = mean_s;
                  }
               }
            }
         else
            {
            // Valid SNP value observed case
            if ( fPerformZeroMean )
               {
               if ( normalizationInfo_.normalizationType == Normal )
                  {
                  SNPs[ iIndividual ] -= mean_s;     //subtract the mean from the data
                  }
               else if ( normalizationInfo_.normalizationType == HardyWeinberg )
                  {
                  SNPs[ iIndividual ] -= (2.0 * freq);
                  }
               else if ( normalizationInfo_.normalizationType == Beta )
                  {
                  // ensure freq is Minor Allele Frequencey (MAF)
                  if ( freq > 0.5 )
                     {
                     freq = 1.0 - freq;
                     }
                  real rT = Beta::Pdf( freq, normalizationInfo_.betaA, normalizationInfo_.betaB );
                  SNPs[ iIndividual ] *= rT;
                  }
               else
                  {
                  Fatal( "Program error.  Illegal normalizationType value %d", normalizationInfo_.normalizationType );
                  }
               }
            if ( fPerformUnitVariance  )
               {
               SNPs[ iIndividual ] /= std;        //unit variance as well
               }
            }
         }
      }
   }

void replaceSNPsByMissingness( real *matrix, const size_t nRows, const size_t nColumns )
   {
   for (size_t column=0; column<nColumns; ++column )
      {
      for ( size_t row=0; row<nRows; ++row )
         {
         if ( matrix[ (column*nRows) + row ] == matrix[ (column*nRows) + row ] )
            {
            matrix[ (column*nRows) + row ]= (real)0.0;
            }
         else
            {
            matrix[ (column*nRows) + row ]=(real)2.0;
            }
         }
      }
   }

//Parameters: SNPs [nIndividuals by nSNPs] Matrix stored in column-major order. This will hold the result.
void preprocessSNPsForBuildingCovariance( real *SNPs, const size_t nIndividuals, const size_t nSNPs, const size_t nSNPsdivide )
   {
   real one_over_sqrt_nSNPs = sqrt( (real)nSNPsdivide );
   size_t snpMax = nIndividuals * nSNPs;

   for ( size_t indiv = 0; indiv < nIndividuals; ++indiv )
      {
      real mean_i = (real)0.0;      //the mean of one individual over all SNPs won't be touched if RRM is used.

      if ( kernelType == COVARIANCE )
         {//covariance needs Zero-meaned SNPs
         size_t count=0;
         for ( size_t snp=indiv; snp < snpMax; snp += nIndividuals )
            {
            if ( SNPs[ snp ] != SNPs[ snp ] )
               {
               Fatal( "\nThere is a NaN at pos %i in the covarianceData", snp );
               }

            mean_i += SNPs[ snp ];
            ++count;
            }

         mean_i /= nSNPsdivide;     //compute the mean of the current individual over all SNPs
         }

      for( size_t snp=indiv; snp < snpMax; snp += nIndividuals )
         {
         SNPs[ snp ] -= mean_i;
         SNPs[ snp ] /= one_over_sqrt_nSNPs;
         }
      }
   }

//zero-mean the columns of a column-major matrix [nRows * nColumns].
void zeroMeanColumns( real *matrix, const size_t nRows, const size_t nColumns, bool meanImputeNaN )
   {
   size_t nEntries = nColumns * nRows;

   for ( size_t column = 0; column < nEntries; column += nRows )
      {
      size_t colMax = column + nRows;
      size_t n_obs = 0;
      real colSum = 0.0;
      for (size_t entry = column; entry < colMax; ++entry )
         {
         if( matrix[entry] == matrix[entry] )
            {
            colSum += matrix[ entry ];
            ++n_obs;
            }
         }

      colSum /= n_obs;
      for ( size_t entry = column; entry < colMax; ++entry )
         {
         if( matrix[ entry ] == matrix[ entry ] )
            {
            matrix[ entry ] -= colSum;
            }
         else
            {
            if ( meanImputeNaN )
               {
               matrix[ entry ] = 0.0;
               }
            }
         }
      }
   }
