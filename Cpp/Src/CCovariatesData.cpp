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

CCovariatesData::CCovariatesData()
   {
   n_individuals = 0;
   n_covariates = 0;
   }

CCovariatesData::~CCovariatesData()
   {
   }

/*
 * Clear all the collections from CovariatesData
 */
void CCovariatesData::clear()
   {
   n_individuals = 0;
   n_covariates = 0;
   filename.clear();

   indIDs.clear();
   famIDs.clear();
   covariates.clear();
   markedAsConstant.clear();       // indicator of constant covariates
   keyIndividualToCovariatesDataIndividualIndex.clear();
   }

/*
 *  Read a PLink formated Covariate file into this object
 */
void CCovariatesData::ReadCovariatesFile( string& covariatesFile )
   {
   filename = covariatesFile;
   n_individuals = 0;
   CToken tok;
   vector<real> covariatesCurrent;     // keep constructor out of loop

   CPlinkLexer lex( covariatesFile );
   lex.NextToken( tok );

   while ( tok.type != tokEOF )        // for each line until EOF
      {
      if ( tok.type == tokEOL )        // skip blank line...or comment line or ???
         {
         lex.NextToken( tok );
         continue;
         }

      string idFamily;
      string idIndividual;
      real covarCurrent;
      PhenotypeType covarCurrentType = PHENOTYPE_UNINITIALIZED;

      lex.ExpectId( tok, idFamily, "FamilyID" );
      famIDs.push_back( idFamily );
      lex.ExpectId( tok, idIndividual, "IndividualID" );
      indIDs.push_back( idIndividual );
      
      covariatesCurrent.clear();       // start this iteration w/ a 'clean' vector

      while( (tok.type != tokEOL) && (tok.type != tokEOF) )
         {
         lex.ExpectPhenotype( tok, covarCurrent, covarCurrentType );
         covariatesCurrent.push_back( covarCurrent );
         }

      if ( (n_individuals > 0) && (covariatesCurrent.size() != n_covariates) )
         {
         Fatal( "\nExpected  %i covariates, found %i", n_covariates, covariatesCurrent.size() );
         }
      else if ( n_individuals == 0 )
         {
         n_covariates = covariatesCurrent.size();
         }

      // Validate we have a unique family:individual name before we add it to the index.
      string key = KeyFromIdFamilyAndIdIndividual( idFamily, idIndividual );
      if ( keyIndividualToCovariatesDataIndividualIndex.count( key ) )
         {
         Fatal( "Duplicate Family:Individual id [%s] found in covariates file near lines %d and %d", key.c_str(), keyIndividualToCovariatesDataIndividualIndex[ key ]+1, tok.line );
         }

      keyIndividualToCovariatesDataIndividualIndex[ key ] = covariates.size();
      covariates.push_back( covariatesCurrent );
      ++n_individuals;
      }

   /*
    * If no covariates found, clear ourselves and the covariates 
    *   filename to prevent trying to reload the file.
    */
   if ( n_covariates == 0 )
      {
      Warn( "No covariates found in file %s.  Ignoring file.", covariatesFile.c_str() );
      clear();
      covariatesFile.clear();
      }
   }

// looks for constants among the covariates and removes them from the analysis.
void CCovariatesData::MarkConstantCovariates()
   {
   for ( size_t cov = 0; cov < n_covariates; ++cov )
      {
      bool onlyNaN = true;    // so far we have observed only NaN for the feature
      bool constant = true;   // In the beginning we expect the feature to be constant
      real val;               // this holds the first value found.

      for ( size_t ind = 0; ind < n_individuals; ++ind )
         {
         if( onlyNaN )
            {
            // This is the case, if all individuals up to here are NaN
            if ( !_isnan( covariates[ind][cov] ) )
               {
               // not NaN!  We have the first 'observed value'
               val = covariates[ind][cov];
               onlyNaN = false;
               }
            }
         else
            {
            // There has already been an observed value.
            if ( (!_isnan( covariates[ind][cov] ))
              && (val != covariates[ind][cov])
               )
               {
               // not NaN and the covariate is not constant.
               //   use in the further analysis.
               constant = false;
               break;
               }
            }
         }//end for loop over individuals

      if ( constant )
         {
         //TODO: somehow remove the covariate
         if ( onlyNaN )
            {
            Warn( "Covariate number %i consists only of NaN values", cov );
            }
         else
            {
            Warn( "Covariate number %i is constant", cov );
            }
         }
      }
   }


void CCovariatesData::CreateColumnMajorArrayFromVectorsAndMapping( const vector<string>& individualsLabelOrder, real *covariatesArray, size_t n_covariates_ )
   {
//   if ( individualsLabelOrder.size() != n_individuals )
//      {
//      Warn( "Covariates File contains %i individuals, Only %i observed phenotypes are specified.", n_individuals, individualsLabelOrder.size() );
//      }
   if ( n_covariates_ != n_covariates )
      {
      Fatal( "Covariates file contains %i covariates, size of array is specified as %i covariates", n_covariates_, n_covariates );
      }
   size_t cIndividualsOut = individualsLabelOrder.size();
   string keyIndividual;
   for ( size_t iIndividualOut = 0; iIndividualOut < cIndividualsOut; ++iIndividualOut )
      {
      keyIndividual = individualsLabelOrder[ iIndividualOut ];
      if ( !FKeyIndividualInDatFile( keyIndividual ) )
         {
         Fatal( "Cannot create column major array from Vector and Mapping.  individualIDs of covariates do not contain mapping key [%s]", keyIndividual.c_str() );
         }
      
      size_t idxIndividual = keyIndividualToCovariatesDataIndividualIndex[ keyIndividual ];
      for ( size_t iCovariate = 0; iCovariate < n_covariates; ++iCovariate )
         {
         size_t idxOut = ColumnMajorIndex( cIndividualsOut, n_covariates_, iIndividualOut, iCovariate );
         covariatesArray[ idxOut ] = covariates[idxIndividual][iCovariate];
         }
      }
   }




void CCovariatesData::ImputeAndZeroMeanCovariatesFix( const vector< string >& selectedIndividuals )
   {
   // for each covariate, enumerate the selected people, get the average/mean of that
   //  covariate, then re-state the value
   size_t cIndividuals = selectedIndividuals.size();
   vector< real > rgAverage( n_covariates, 0.0 );
   for ( size_t iCovariate = 0; iCovariate<n_covariates; ++iCovariate )
      {
      size_t cValidReals = 0;
      real  sumOfCovariates = 0.0;
      // get the mean of covariates
      for ( size_t iIndividual=0; iIndividual < cIndividuals; ++iIndividual )
         {
         if ( !FKeyIndividualInDatFile( selectedIndividuals[ iIndividual ] ) )
            {
            Fatal( "Cannot ImputeAndZeroMeanCovariates().  [%s] not found.", selectedIndividuals[ iIndividual ].c_str() );
            }
         size_t idxIndividual = keyIndividualToCovariatesDataIndividualIndex[ selectedIndividuals[ iIndividual ] ];
         if ( !_isnan( covariates[ idxIndividual ][ iCovariate ] ) )
            {
            sumOfCovariates += covariates[ idxIndividual ][ iCovariate ];
            //cout << covariates[ idxIndividual ][ iCovariate ] << endl;
            ++cValidReals;
            }
         }

      real meanOfCovariates = sumOfCovariates / (real)cValidReals;

      // fill out the missing covariates
      for ( size_t iIndividual=0; iIndividual < cIndividuals; ++iIndividual )
         {
         if ( !FKeyIndividualInDatFile( selectedIndividuals[ iIndividual ] ) )
            {
            Fatal( "Cannot ImputeAndZeroMeanCovariates().  [%s] not found.", selectedIndividuals[ iIndividual ].c_str() );
            }
         size_t idxIndividual = keyIndividualToCovariatesDataIndividualIndex[ selectedIndividuals[ iIndividual ] ];
         if ( !_isnan( covariates[ idxIndividual ][ iCovariate ] ) )
            {
            covariates[ idxIndividual ][ iCovariate ] -= meanOfCovariates;
            }
         else
            {
            covariates[ idxIndividual ][ iCovariate ] = 0.0;
            }
         }
      }
   }





void CCovariatesData::ImputeAndZeroMeanCovariates( const vector< string >& selectedIndividuals )
   {
   // for each covariate, enumerate the selected people, get the average/mean of that
   //  covariate, then re-state the value
   size_t cIndividuals = selectedIndividuals.size();
   vector< real > rgAverage( n_covariates, 0.0 );
   for ( size_t iIndividual=0; iIndividual < cIndividuals; ++iIndividual )
      {
      if ( !FKeyIndividualInDatFile( selectedIndividuals[ iIndividual ] ) )
         {
         Fatal( "Cannot ImputeAndZeroMeanCovariates().  [%s] not found.", selectedIndividuals[ iIndividual ].c_str() );
         }
      size_t idxIndividual = keyIndividualToCovariatesDataIndividualIndex[ selectedIndividuals[ iIndividual ] ];
      size_t cValidReals = 0;
      real  sumOfCovariates = 0.0;
      for ( size_t iCovariate = 0; iCovariate<n_covariates; ++iCovariate )
         {
         if ( !_isnan( covariates[ idxIndividual ][ iCovariate ] ) )
            {
            sumOfCovariates += covariates[ idxIndividual ][ iCovariate ];
            ++cValidReals;
            }
         }

      real meanOfCovariates = sumOfCovariates / (real)cValidReals;

      for ( size_t iCovariate = 0; iCovariate<n_covariates; ++iCovariate )
         {
         if ( !_isnan( covariates[ idxIndividual ][ iCovariate ] ) )
            {
            covariates[ idxIndividual ][ iCovariate ] -= meanOfCovariates;
            }
         else
            {
            covariates[ idxIndividual ][ iCovariate ] = 0.0;
            }
         }
      }
   }

size_t CCovariatesData::IdxFromKeyIndividual( const string& keyIndividual )
   {
   if ( keyIndividualToCovariatesDataIndividualIndex.count( keyIndividual ) == 0 )
      {
      Fatal( "Unable to find individual [%s] in Covariates data.", keyIndividual.c_str() );
      }
   size_t idx = keyIndividualToCovariatesDataIndividualIndex[ keyIndividual ];
   return( idx );
   }

size_t CCovariatesData::IdxFromIdFamilyAndIdIndividual( const string& idFamily, const string& idIndividual )
   {
   string key = KeyFromIdFamilyAndIdIndividual( idFamily, idIndividual );
   size_t idx = IdxFromKeyIndividual( key );
   return( idx );
   }

bool CCovariatesData::FKeyIndividualInDatFile( const string& keyIndividual )
   {
   return( keyIndividualToCovariatesDataIndividualIndex.count( keyIndividual ) != 0 );
   }

// Retrieve the covariate value associated with tuple [ keyIndividual ][ iCovariate ]
real CCovariatesData::CovariateValue( size_t iCovariate, const string& keyIndividual )
   {
   if ( iCovariate >= n_covariates )
      {
      Fatal( "Covariate lookup out of range.  Expected less than [%d] and found [%d]", n_covariates, iCovariate );
      }
   size_t idxIndividual = IdxFromKeyIndividual( keyIndividual );
   return( covariates[ idxIndividual ][ iCovariate ] );
   }

// Fill a an array (vector) with the covariate values associated with the 
//   values from the keyIndividual vector and iCovariate
void CCovariatesData::CovariateValues( size_t iCovariate, const vector< string >& rgKeyIndividuals, real *rgOutArray )
   {
   // Possibly hoist iCovariates range check out of loop if needed for performance?
   for ( auto it = rgKeyIndividuals.begin(); it != rgKeyIndividuals.end(); ++it )
      {
      *rgOutArray = CovariateValue( iCovariate, *it );
      ++rgOutArray;
      }
   }
