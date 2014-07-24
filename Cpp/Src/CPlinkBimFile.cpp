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
 * CPlinkBimFile - {PLINK Bim File Access Class}
 *
 *         File Name:   CPlinkBimFile.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:    4 Dec 2011
 *     Revision Date:    4 Dec 2011
 *
 *    Module Purpose:   This file implements the CPlinkBimFile class 
 *                         for FastLmmC
 *
 *                      A .MAP file contains the following fields related
 *                         to SNPs.
 *                            chromosome (1-22, X, Y or 0 if unplaced)
 *                            rs# or snp identifier         string w/o whitespace
 *                            Genetic distance (morgans)    floating point number
 *                            Base-pair position (bp units) integral number
 *
 *                      The .BIM file shares these 4 fields and adds two more
 *                            Minor Allele                  one character
 *                            Major Allele                  one character
 *
 *                      This is all represented by a SnpInfo struct
 *
 *    Change History:   
 *
 * Test Files: 
 */

#include "FastLmmC.h"

CPlinkBimFile::CPlinkBimFile()
   {
   }

CPlinkBimFile::CPlinkBimFile( const string& filename_ )
   {
   filename = filename_;
   }

CPlinkBimFile::~CPlinkBimFile()
   {
   }

/*
 *  Load the .FAM file saved in the member filename
 *    into memory and make the elements available for 
 *    later access
 */
void CPlinkBimFile::Load()
   {
   Verbose( verboseMinor, "                  Loading .BIM file: [%s]", filename.c_str() );
   CTimer timer( true );
   SnpInfo si;                         // the working storage
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

      lex.ExpectChromosome( tok, si.idChromosome, si.iChromosome );
      lex.ExpectId( tok, si.idSnp, "SNP Id" );
      lex.ExpectReal( tok, si.geneticDistance, "Genetic Distance" );
      lex.ExpectInt( tok, si.basepairPosition, "Basepair Position" );
      lex.ExpectSnpAlleles( tok, si.majorAllele, si.minorAllele );

      // Validate we have a unique SNP id and add it to the map
      if ( idSnpToSnpInfoIndex.count( si.idSnp ) )
         {
         Fatal( "Duplicate SNP Id [%s] found in lines %d and %d", 
                        si.idSnp.c_str(), 
                        idSnpToSnpInfoIndex[ si.idSnp ] + 1,
                        rgSnpInfo.size() + 1
                        );
         }

      idSnpToSnpInfoIndex[ si.idSnp ] = rgSnpInfo.size();
      rgSnpInfo.push_back( si );
      }

   timer.Stop();
   timer.Report( 0x02, "     Loading .MAP file elapsed time: %s" );
   }

/*
 *  Load the .BIM file passed to this routine into 
 *    memory and make the elements available for 
 *    later access
 */
void CPlinkBimFile::Load( const string& filename_ )
   {
   filename = filename_;
   Load();
   }


/*
 *  Get SnpInfo associated with idSnp 
 */
size_t CPlinkBimFile::IndexFromIdSnp( const string& idSnp )
   {
   size_t idxSnp = -1;

   if ( idSnpToSnpInfoIndex.count( idSnp ) != 0 )
      {
      // idSnp is found in the map
      idxSnp = idSnpToSnpInfoIndex[ idSnp ];
      }
   return( idxSnp );
   }

/*
 * Get SnpInfo associated with index
 */
SnpInfo* CPlinkBimFile::PSnpInfo( size_t idxSnpInfo )
   {
   SnpInfo *pSnpInfo = nullptr;

   if ( idxSnpInfo < rgSnpInfo.size() )
      {
      pSnpInfo = &(rgSnpInfo[ idxSnpInfo ]);
      }
   return( pSnpInfo );
   }

/*
 *  Get SnpInfo associated with idSnp 
 */
SnpInfo* CPlinkBimFile::PSnpInfo( const string& idSnp )
   {
   return( PSnpInfo( IndexFromIdSnp( idSnp ) ) );
   }
