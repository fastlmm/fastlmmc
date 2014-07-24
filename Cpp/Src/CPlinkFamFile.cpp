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
 * CPlinkFamFile - {PLINK Fam File Access Class}
 *
 *         File Name:   CPlinkFamFile.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:    4 Dec 2011
 *     Revision Date:    4 Dec 2011
 *
 *    Module Purpose:   This file implements the CPlinkFamFile class 
 *                         for FastLmmC
 *
 *                      A .MAP file contains the following fields related
 *                         to SNPs.  
 *                            chromosome (1-22, X, Y or 0 if unplaced)
 *                            rs# or snp identifier
 *                            Genetic distance (morgans)
 *                            Base-pair position (bp units)
 *                      This is all represented by a SnpInfo struct
 *
 *    Change History:   
 *
 * Test Files: 
 */

#include "FastLmmC.h"

CPlinkFamFile::CPlinkFamFile()
   {
   phenotypeType = PHENOTYPE_UNINITIALIZED;
   }

CPlinkFamFile::CPlinkFamFile( const string& filename_ )
   {
   filename = filename_;
   }

CPlinkFamFile::~CPlinkFamFile()
   {
   }

/*
 *  Load the .FAM file saved in the member filename
 *    into memory and make the elements available for 
 *    later access
 */
void CPlinkFamFile::Load()
   {
   Verbose( verboseMinor, "                  Loading .FAM file: [%s]", filename.c_str() );
   CTimer timer( true );
   IndividualInfo ii;                  // the working storage
   CToken tok;
   CPlinkLexer lex( filename );

   lex.NextToken( tok );
   while ( tok.type != tokEOF )        // for each line until EOF
      {
      if ( tok.type == tokEOL )        // blank line...or comment line or ???
         {
         lex.NextToken( tok );
         continue;
         }
      lex.ExpectId( tok, ii.idFamily, "FamilyID" );
      lex.ExpectId( tok, ii.idIndividual, "IndividualID" );
      lex.ExpectId( tok, ii.idPaternal, "PaternalID" );
      lex.ExpectId( tok, ii.idMaternal, "MaternalID" );
      lex.ExpectSex( tok, ii.sex );
      lex.ExpectPhenotype( tok, ii.phenotype, phenotypeType );
      
      // Validate we have a unique key for this individual and add it to the map
      string key = KeyFromIdFamilyAndIdIndividual( ii.idFamily, ii.idIndividual );
      if ( idKeyToIndividualInfoIndex.count( key ) )
         {
         Fatal( "Duplicate Individual Id [%s:%s] found in lines %d and %d", 
                     ii.idFamily.c_str(), 
                     ii.idIndividual.c_str(), 
                     idKeyToIndividualInfoIndex[ key ] + 1, 
                     idKeyToIndividualInfoIndex.size() + 1
                     );
         }

      idKeyToIndividualInfoIndex[ key ] = rgIndividualInfo.size();
      rgIndividualInfo.push_back( ii );
      }

   timer.Stop();
   timer.Report( 0x02, "     Loading .MAP file elapsed time: %s" );
   }

/*
 *  Load the .FAM file passed to this routine into 
 *    memory and make the elements available for 
 *    later access
 */
void CPlinkFamFile::Load( const string& filename_ )
   {
   filename = filename_;
   Load();
   }

/*
 *  Get Index associated with idKey
 */
size_t CPlinkFamFile::IndexFromKeyIndividual( const string& idKey )
   {
   size_t idxIndividual = -1;

   if ( idKeyToIndividualInfoIndex.count( idKey ) != 0 )
      {
      // idKey is found in the map
      idxIndividual = idKeyToIndividualInfoIndex[ idKey ];
      }
   return( idxIndividual );
   }

/*
 * Get IndividualInfo associated with index
 */
IndividualInfo* CPlinkFamFile::PIndividualInfo( size_t idxIndividualInfo )
   {
   IndividualInfo *pIndividualInfo = nullptr;

   if ( idxIndividualInfo < rgIndividualInfo.size() )
      {
      pIndividualInfo = &(rgIndividualInfo[ idxIndividualInfo ]);
      }
   return( pIndividualInfo );
   }

/*
 *  Get IndividualInfo associated with idKey
 */
IndividualInfo* CPlinkFamFile::PIndividualInfo( const string& idKey )
   {
   return( PIndividualInfo( IndexFromKeyIndividual( idKey ) ) );
   }
