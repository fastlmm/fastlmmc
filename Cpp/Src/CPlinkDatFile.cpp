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
 * CPlinkDatFile - {PLINK Dat File Access Class}
 *
 *         File Name:   CPlinkDatFile.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:    4 Dec 2011
 *     Revision Date:    4 Dec 2011
 *
 *    Module Purpose:   This file implements the CPlinkDatFile class 
 *                         for FastLmmC
 *
 *    Change History:   
 *
 * Test Files: 
 */

#include "FastLmmC.h"

CPlinkDatFile::CPlinkDatFile( const string& filename_, PlinkFileType fileType_ )
   {
   filename = filename_;
   pFile = nullptr;
   fileType = fileType_;
   }

//CPlinkDatFile::CPlinkDatFile( const string& filename_, size_t cIndividuals_, size_t cSnps_ ) {}

CPlinkDatFile::~CPlinkDatFile()
   {
   if ( pFile != nullptr )
      {
      fclose( pFile );
      pFile = nullptr;
      }
   }

void CPlinkDatFile::Load()
   {
   Load( nullptr );
   }

#if defined( _MSC_VER )      // Visual C++-only pragmas
#pragma optimize( "", off )
#pragma warning( push )
#pragma warning( disable : 4748 )
#endif
void CPlinkDatFile::Load( vector<FamRecord>* prgFam_ )
   {
   SnpDosageData         snpDD;
   vector<SnpDosageData> rgSnpDD;
   DatRecord             dr;
   IndividualId          individualId;
   size_t                cIndividualsExpected = prgFam_ ? prgFam_->size() : 0;
   string allele1Header;
   string allele2Header;

   Verbose( verboseMajor, "                  Loading .DAT file: [%s]", filename.c_str() );
   CTimer timer( true );
   CToken tok;
   CPlinkLexer lex( filename );

   do { // skip leading blank lines and comments.
      lex.NextToken( tok );
      } while ( tok.type == tokEOL );
   if ( (tok.type == tokSymbol) && (tok.text == "SNP") )
      {
      // Header row in file

      lex.NextToken( tok );      // advance past "SNP" and throw away
      lex.ExpectId( tok, allele1Header, "Allele1_Header" );    // advance past "A1" and throw away
      lex.ExpectId( tok, allele2Header, "Allele2_Header" );    // advance past "A2" and throw away
      while( (tok.type != tokEOF) && (tok.type != tokEOL) )
         {
         // Get the Family/Individual ID pairs
         lex.ExpectId( tok, individualId.idFamily, "FamilyID" );
         lex.ExpectId( tok, individualId.idIndividual, "IndividualID" );

         // Validate we have a unique key for this individual
         string key = KeyFromIdFamilyAndIdIndividual( individualId.idFamily, individualId.idIndividual );
         if ( FKeyIndividualInDatFile( key ) )
            {
            Fatal( "Duplicate FamilyId:IndividualId [%s:%s] found in header elements %d and %d", 
                  individualId.idFamily.c_str(),
                  individualId.idIndividual.c_str(),
                  keyIndividualToDatRecordIndividualIndex[ key ]+1, 
                  keyIndividualToDatRecordIndividualIndex.size()+1 );
            }

         keyIndividualToDatRecordIndividualIndex[ key ] = rgIndividualIds.size();
         rgIndividualIds.push_back( individualId );
         }

      cIndividualsExpected = rgIndividualIds.size();
      }

   // No header row or bogus header row
   if ( cIndividualsExpected == 0 )
      {
      Fatal( "Expected a header row starting with \"SNP\" in file [%s]", filename.c_str() );
      }

   // We are either at the first token of the file without a header or the first token of a new line
   while ( tok.type != tokEOF )        // for each line until EOF
      {
      if ( tok.type == tokEOL )        // blank line...or comment line or ???
         {
         lex.NextToken( tok );
         continue;
         }

      rgSnpDD.clear();
      lex.ExpectId( tok, dr.idSnp, "SnpID" );
      lex.ExpectSnpAlleles( tok, dr.majorAllele, dr.minorAllele );
      while( (tok.type != tokEOF) && (tok.type != tokEOL) )
         {
         if ( fileType == Dosage1 )
            {
            lex.ExpectSnpDosage( tok, snpDD.dosage );
            }
         else
            {
            lex.ExpectSnpProbabilities( tok, snpDD );
            }
         rgSnpDD.push_back( snpDD );
         }

      if ( rgSnpDD.size() != cIndividualsExpected )
         {
#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
         const char *szFmt1 = "Expected SNP probability pairs for %Iu individuals on line %Iu.  Found %Iu";
#else                      // Linux/g++ uses a %zu specifier for size_t
         const char *szFmt1 = "Expected SNP probability pairs for %zu individuals on line %zu.  Found %zu";
#endif
         Fatal( szFmt1, cIndividualsExpected, tok.line, rgSnpDD.size() );
         }

      if ( idSnpToDatRecordIndex.count( dr.idSnp ) )
         {
#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
         const char *szFmt2 = "Duplicate SNP Id %s found in SNPs %Iu and %Iu";
#else                      // Linux/g++ uses a %zu specifier for size_t
         const char *szFmt2 = "Duplicate SNP Id %s found in SNPs %zu and %Iu";
#endif
         Fatal( szFmt2, dr.idSnp.c_str(), idSnpToDatRecordIndex[ dr.idSnp ]+1, rgDat.size()+1 );
         }

      idSnpToDatRecordIndex[ dr.idSnp ] = rgDat.size();
      dr.rgSnpDosageData = rgSnpDD;
      rgDat.push_back( dr );
      }

   if ( rgIndividualIds.size() == 0 )
      {
      /*
       * We got here because we have no header row but we do have data.  
       *   Copy the IndividualId from the rgFam information and construct
       *   the needed indexing support.
       */
      for( size_t iIndividual = 0; iIndividual < prgFam_->size(); ++iIndividual )
         {
         FamRecord* pfr = &(prgFam_->at( iIndividual ));
         individualId.idFamily = pfr->idFamily;
         individualId.idIndividual = pfr->idIndividual;
         string key = KeyFromIdFamilyAndIdIndividual( individualId.idFamily, individualId.idIndividual );
         keyIndividualToDatRecordIndividualIndex[ key ] = rgIndividualIds.size();
         rgIndividualIds.push_back( individualId );
         }
      }
   timer.Stop();
   timer.Report( 0x02, "     Loading .DAT file elapsed time: %s" );

   // see if we need to write out a Dosage1 formated file
   char *pchFastLmmEmitSnpDosageData = getenv( "FastLmmEmitSnpDosageData" );
   if ( pchFastLmmEmitSnpDosageData != nullptr )
      {
      if ( *pchFastLmmEmitSnpDosageData == '1' )
         {
         // emit SnpDosageData in Dosage1
         string ofilename = filename + "1";
         FILE* pfile = fopen( ofilename.c_str(), "w" );
         if ( pfile != nullptr )
            {
            Progress( "\n     Writing new SnpDosageData file: [%s]", ofilename.c_str() );
            // if we have a header row, write it to the output file.
            if ( allele1Header.size() )
               {
               fprintf( pfile, "SNP %s %s", allele1Header.c_str(), allele2Header.c_str() );
               for ( vector< IndividualId >::iterator it = rgIndividualIds.begin(); it != rgIndividualIds.end(); ++it )
                  {
                  fprintf( pfile, "  %s %s", it->idFamily.c_str(), it->idIndividual.c_str() ); 
                  }
               fprintf( pfile, "\n" );
               }
            for ( vector< DatRecord >::iterator it = rgDat.begin(); it != rgDat.end(); ++it )
               {
               fprintf( pfile, "%s  %c %c", it->idSnp.c_str(), it->minorAllele, it->majorAllele );
               for ( vector< SnpDosageData >::iterator itDD = it->rgSnpDosageData.begin(); itDD != it->rgSnpDosageData.end(); ++itDD )
                  {
                  real r = SnpDosageData::RealFromSnpDosageData( *itDD, fileType );
                  if ( _isnan( r ) ) { r = missingGenotypeValue; }
                  fprintf( pfile, "  %g", r );
                  }
               fprintf( pfile, "\n" );
               }
            fclose( pfile );
            }
         else
            {
            Warn( "Unable to create DosageData file [%s]", ofilename.c_str() );
            }
         }
      }
   }
#if defined( _MSC_VER )
#pragma warning( pop )
#pragma optimize( "", on )
#endif

bool CPlinkDatFile::FGetSnpInfo( const string& idSnp, SnpInfo& snpInfo_ )
   {
   snpInfo_.Clear();

   if ( idSnpToDatRecordIndex.count( idSnp ) == 0 )
      {
      return( false );
      }

   size_t idx = idSnpToDatRecordIndex[ idSnp ];
   snpInfo_.idSnp = rgDat[ idx ].idSnp;
   snpInfo_.majorAllele = rgDat[ idx ].majorAllele;
   snpInfo_.minorAllele = rgDat[ idx ].minorAllele;
   return( true );
   }

bool CPlinkDatFile::FGetSnpInfo( size_t idx, SnpInfo& snpInfo_ )
   {
   snpInfo_.Clear();

   if ( idx >= rgDat.size() )
      {
      return( false );
      }

   snpInfo_.idSnp = rgDat[ idx ].idSnp;
   snpInfo_.majorAllele = rgDat[ idx ].majorAllele;
   snpInfo_.minorAllele = rgDat[ idx ].minorAllele;
   return( true );
   }

bool CPlinkDatFile::FSnpHasVariation( const string& idSnp )
   {
   size_t idxSnp = IdxFromIdSnp( idSnp );
   return( FSnpHasVariation( idxSnp ) );
   }

bool CPlinkDatFile::FSnpHasVariation( size_t idxSnp )
   {
   DatRecord* dr = DatRecordPointer( idxSnp );
   SnpDosageData dosageData(-1.0, -1.0);          // keep track of values

   for ( size_t iDosageData=0; iDosageData < dr->rgSnpDosageData.size(); ++ iDosageData )
      {
      const SnpDosageData& dosageDataT = dr->rgSnpDosageData[ iDosageData ];
      if ( dosageDataT.probabilityOfHomozygousMinor != missingGenotypeValue )
         {
         // Data is NOT missing, evaluate for variation in the SNP values.
         if ( dosageData.probabilityOfHomozygousMinor == -1.0 )
            {
            dosageData = dosageDataT;
            }
         else if ( dosageData != dosageDataT )
            {
            return( true );
            }
         }
      }
   return( false );
   }

bool CPlinkDatFile::FSnpHasVariation( const string& idSnp, const vector< size_t >& selected )
   {
   size_t idxSnp = IdxFromIdSnp( idSnp );
   return( FSnpHasVariation( idxSnp, selected ) );
   }

bool CPlinkDatFile::FSnpHasVariation( size_t idxSnp, const vector< size_t >& selected )
   {
   DatRecord* dr = DatRecordPointer( idxSnp );
   SnpDosageData dosageData(-1.0, -1.0);          // keep track of values

   for ( size_t iSelected=0; iSelected<selected.size(); ++iSelected )
      {
      size_t iDosageData = selected[ iSelected ];
      const SnpDosageData& dosageDataT = dr->rgSnpDosageData[ iDosageData ];
      if ( dosageDataT.probabilityOfHomozygousMinor != missingGenotypeValue )
         {
         // Data is NOT missing, evaluate for variation in the SNP values.
         if ( dosageData.probabilityOfHomozygousMinor == -1.0 )
            {
            dosageData = dosageDataT;
            }
         else if ( dosageData != dosageDataT )
            {
            return( true );
            }
         }
      }
   return( false );
   }

SnpDosageData CPlinkDatFile::GetSnpDosageData( const string& idSnp, const string& idIndividual, const string& idFamily )
   {
   size_t idxSnp = IdxFromIdSnp( idSnp );
   size_t idxIndividual = IdxFromIdFamilyAndIdIndividual( idFamily, idIndividual );
   return( rgDat[ idxSnp ].rgSnpDosageData[ idxIndividual ] );
   }

SnpDosageData CPlinkDatFile::GetSnpDosageData( const string& idSnp, const string& keyIndividual )
   {
   size_t idxSnp = IdxFromIdSnp( idSnp );
   size_t idxIndividual = IdxFromKeyIndividual( keyIndividual );
   return( rgDat[ idxSnp ].rgSnpDosageData[ idxIndividual ] );
   }

SnpDosageData CPlinkDatFile::GetSnpDosageData( size_t idxSnp, size_t idxIndividual )
   {
   return( rgDat[ idxSnp ].rgSnpDosageData[ idxIndividual ] );
   }

size_t CPlinkDatFile::IdxFromIdSnp( const string& idSnp )
   {
   if ( idSnpToDatRecordIndex.count( idSnp ) == 0 )
      {
      Fatal( "Unable to find SnpId [%s] in .DAT SNP information.", idSnp.c_str() );
      }
   size_t idx = idSnpToDatRecordIndex[ idSnp ];
   return( idx );
   }

size_t CPlinkDatFile::IdxFromKeyIndividual( const string& keyIndividual )
   {
   if ( keyIndividualToDatRecordIndividualIndex.count( keyIndividual ) == 0 )
      {
      Fatal( "Unable to find individual [%s] in .DAT genotype.", keyIndividual.c_str() );
      }
   size_t idx = keyIndividualToDatRecordIndividualIndex[ keyIndividual ];
   return( idx );
   }

size_t CPlinkDatFile::IdxFromIdFamilyAndIdIndividual( const string& idFamily, const string& idIndividual )
   {
   string key = KeyFromIdFamilyAndIdIndividual( idFamily, idIndividual );
   size_t idx = IdxFromKeyIndividual( key );
   return( idx );
   }

bool CPlinkDatFile::FKeyIndividualInDatFile( const string& keyIndividual )
   {
   return( keyIndividualToDatRecordIndividualIndex.count( keyIndividual ) != 0 );
   }

bool CPlinkDatFile::FIdSnpInDatFile( const string& idSnp )
   {
   return( idSnpToDatRecordIndex.count( idSnp ) != 0 );
   }

DatRecord* CPlinkDatFile::DatRecordPointer( const string& idSnp )
   {
   if ( idSnpToDatRecordIndex.count( idSnp ) == 0 )
      {
      return( nullptr );
      }

   size_t idx = idSnpToDatRecordIndex[ idSnp ];
   DatRecord* dr = &rgDat[ idx ];
   return( dr );
   }

DatRecord* CPlinkDatFile::DatRecordPointer( size_t idxSnp )
   {
   if ( idxSnp > rgDat.size() )
      {
      return( nullptr );
      }

   DatRecord* dr = &rgDat[ idxSnp ];
   return( dr );
   }

