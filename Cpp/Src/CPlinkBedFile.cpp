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
 * CPlinkBedFile - {PLINK BED File Access Class}
 *
 *         File Name:   CPlinkBedFile.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements the CPlinkBedFile 
 *                      class for FastLmmC
 *
 *                      A .BED file contains compressed binary genotype 
 *                         values for for individuals by SNPs.  
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Include Files
 */
#include "FastLmmC.h"

real unknownOrMissing = std::numeric_limits<real>::quiet_NaN();  // now used by SnpInfo
real homozygousPrimaryAllele = 0.0;                // Major Allele
real heterozygousAllele = 1.0;                     
real homozygousSecondaryAllele = 2.0;              // Minor Allele

real mapBedGenotypeToRealAllele[4] = { 
                  homozygousSecondaryAllele,       // look-up 0
                  unknownOrMissing,                // look-up 1
                  heterozygousAllele,              // look-up 2
                  homozygousPrimaryAllele,         // look-up 3
                  };

CBedFile::CBedFile()
   {
   layout = LayoutUnknown;    // layout describes the matrix layout on disk
                              // 0=RowMajor(all snps per individual together);
                              // 1=ColumnMajor(all individuals per SNP together in memory)
   cIndividuals = 0;
   cSnps        = 0;
   cbStride     = 0;

   pFamFile = nullptr;
   pBimFile = nullptr;
   pFile    = nullptr;
   }

CBedFile::~CBedFile()
   {
   if ( pFile )
      {
      fclose( pFile );
      pFile = nullptr;
      }
   }

void CBedFile::Open( const string& filename_, CPlinkFamFile *pFamFile_, CPlinkBimFile *pBimFile_ )
   {
   // cache the associated files used to define the shape of the .BED file.
   pFamFile = pFamFile_;
   pBimFile = pBimFile_;
   Open( filename_, pFamFile->CountOfIndividuals(), pBimFile->CountOfSnpInfo() );
   }

void CBedFile::Open( const string& filename_, size_t cIndividuals_, size_t cSnps_ )
   {
   if ( filename_.empty() )
      {
      Fatal( "Could not create BedFile Reader.  Parameter 'filename' is zero length string" );
      }

   filename = FullPath( filename_ );         // save class local copy of expanded filename
   cIndividuals = cIndividuals_;
   cSnps = cSnps_;

   pFile = fopen( filename.c_str(), "rb" );  // read in binary to ensure ftell works right
   if ( !pFile )
      {
      Fatal( "Cannot open input file [%s].\n  CRT Error %d: %s", filename.c_str(), errno, strerror( errno ) );
      }

   //  Verify 'magic' number
   unsigned char rd1 = NextChar();
   unsigned char rd2 = NextChar();
   if ( (bedFileMagic1 != rd1) || (bedFileMagic2 != rd2))
      {
      Fatal( "Ill-formed BED file [%s]."
           "\n  BED file header is incorrect."
           "\n  Expected magic number of 0x%02x 0x%02x, found 0x%02x 0x%02x", 
               filename.c_str(), bedFileMagic1, bedFileMagic2, rd1, rd2 );
      }

   // Verify 'mode' is valid
   unsigned char rd3 = NextChar();
   switch( rd3 )
      {
   case 0:  // mode = 'IndividualMajor' or RowMajor
      layout = LayoutGroupGenotypesByIndividual;   // all SNPs per individual are sequential in memory
      cbStride = (cSnps + 3)/4;                    // 4 genotypes per byte so round up
      break;
   case 1:  // mode = 'SnpMajor' or ColumnMajor
      layout = LayoutGroupGenotypesBySnp;          // all individuals per SNP are sequential in memory
      cbStride = (cIndividuals + 3)/4;             // 4 genotypes per byte so round up
      break;
   default:
      Fatal( "Ill-formed BED file [%s].  BED file header is incorrect.  Expected mode to be 0 or 1, found %d", filename.c_str(), rd3 );
      break;
      }

   // allocate the read buffer for a SNP
   rgBytes.resize( cbStride );
   rgBedGenotypes.resize( cIndividuals, bedMissingGenotype );
   }

LayoutMode  CBedFile::GetLayoutMode()
   {
   return( layout );
   }

int CBedFile::NextChar()
   {
   int value = fgetc( pFile );
   if ( value == EOF )
      {
      Fatal( "Ill-formed BED file [%s].  Encountered EOF before exepected.", filename.c_str() );
      }
   return( (unsigned char)value );
   }

size_t CBedFile::Read( BYTE *pb, size_t cbToRead )
   {
   size_t cbRead = fread( pb, 1, cbToRead, pFile );
   if ( cbRead != cbToRead )
      {
      if ( feof( pFile ) )
         {
         Fatal( "Encountered EOF before exepected in BED file. Ill-formed BED file [%s]", filename.c_str() );
         }
      int err = ferror( pFile );
      if ( err )
         {
         Fatal( "Encountered a file error %d in BED file [%s]", err, filename.c_str() );
         }
      }
   return( cbRead );
   }

size_t CBedFile::ReadLine( BYTE *pb, size_t idx )
   {
   long long fpos = cbHeader + (idx*cbStride);
   long long fposCur = _ftelli64( pFile );
   if ( fpos != fposCur )
      {
      _fseeki64( pFile, fpos, SEEK_SET );
      }
   
   size_t cbRead = Read( pb, cbStride );
   return( cbRead );
   }

/*
 * Read the genotype for all the individuals in iidList at the SNP specified by iSNP
 *   and store the results in pvOut
 */
void CBedFile::ReadGenotypes( size_t iSnp, const vector< size_t >& idxIndividualList, vector< real >& pvOut )
   {
   ReadLine( &rgBytes[0], iSnp );

   // 'decompress' the genotype information
   size_t iIndividual = 0;
   for ( size_t ib = 0; ib < cbStride; ++ib )
      {
      BYTE genotypeByte = rgBytes[ ib ];

      // manually unrolled loop to decompress this byte
      if ( iIndividual < cIndividuals ) rgBedGenotypes[ iIndividual++ ] = (BedGenotype)( genotypeByte       & 0x03);
      if ( iIndividual < cIndividuals ) rgBedGenotypes[ iIndividual++ ] = (BedGenotype)((genotypeByte >> 2) & 0x03);
      if ( iIndividual < cIndividuals ) rgBedGenotypes[ iIndividual++ ] = (BedGenotype)((genotypeByte >> 4) & 0x03);
      if ( iIndividual < cIndividuals ) rgBedGenotypes[ iIndividual++ ] = (BedGenotype)((genotypeByte >> 6) & 0x03);
      }

   for ( size_t i=0; i<idxIndividualList.size(); ++i )
      {
      size_t idx = idxIndividualList[ i ];
      pvOut[ i ] = mapBedGenotypeToRealAllele[ rgBedGenotypes[ idx ] ];
      }
   }
