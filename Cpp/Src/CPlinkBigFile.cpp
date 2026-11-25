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
 * CPlinkBigFile - {PLINK Big File Access Class}
 *
 *         File Name:   CPlinkBigFile.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements the CPlinkBigFile 
 *                      class for FastLmmC
 *
 *                      A .BIG file contains Binary Imputed Genotype 
 *                         values for individuals by SNPs.  
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Include Files
 */
#include "FastLmmC.h"

#if 0
int ConvertBigEndianInt( BYTE *pb )
   {
   int wT = 0;
   wT  = pb[0];
   wT <<= 8;
   wT += pb[1];
   wT <<= 8;
   wT += pb[2];
   wT <<= 8;
   wT += pb[3];
   return( wT );
   }
#endif

CPlinkBigFile::CPlinkBigFile( const string& filename_, size_t cIndividuals_, size_t cSnps_ )
   {
   filename = FullPath( filename_ );         // save class local copy of expanded filename
   cIndividuals = cIndividuals_;
   cSnps = cSnps_;
   cbStride = 0;

   if ( filename.empty() )
      {
      Fatal( "Could not create BigFile Reader.  Parameter 'filename' is zero length string" );
      }

   pFile = fopen( filename.c_str(), "rb" );  // read in binary to ensure ftell works right
   if ( !pFile )
      {
      Fatal( "Cannot open input file [%s].\n  CRT Error %d: %s", filename.c_str(), errno, strerror( errno ) );
      }

   BYTE rgBuf[32];
   Read( rgBuf, sizeof( rgBuf ) );

   //  Verify 'magic' number
   BYTE magic1 = rgBuf[0];
   BYTE magic2 = rgBuf[1];
   if ( (bigFileMagic1 != magic1) || (bigFileMagic2 != magic2))
      {
      Fatal( "Ill-formed .BIG file [%s]."
           "\n  BIG file header is incorrect."
           "\n  Expected magic number of 0x%02x 0x%02x, found 0x%02x 0x%02x", 
               filename.c_str(), bigFileMagic1, bigFileMagic2, magic1, magic2 );
      }

   // Capture & Verify the version info (4 byte ascii)  expecting 0100
   char rdVer10 = rgBuf[2];
   char rdVer01 = rgBuf[3];
   char rdRel10 = rgBuf[4];
   char rdRel01 = rgBuf[5];

   if ( !isdigit( rdVer10 ) || !isdigit( rdVer01) || !isdigit( rdRel10 ) || !isdigit( rdRel01 ) )
      {
      // poorly formed header as version must be all decimal digits
      Fatal( "Ill-formed .BIG file [%s]."
           "\n  BIG file header is incorrect."
           "\n  Expected version to be 4 decimal digits, found '%c%c%c%c' [0x%02x 0x%02x 0x%02x 0x%02x ]", 
               filename.c_str(), rdVer10, rdVer01, rdRel10, rdRel01, rdVer10, rdVer01, rdRel10, rdRel01 );
      }

   if ( (rdVer10 != '0') || (rdVer01 != '1') )
      {
      Fatal( "Version mismatch in .BIG file [%s]."
           "\n  Expected version '01', found version '%c%c'", filename.c_str(), rdVer10, rdVer01 );
      }

   // Verify we have the Dos EOF marker
   char t = rgBuf[6];
   if ( t != 0x1a )
      {
      Fatal( "Ill-formed .BIG file [%s]."
           "\n  BIG file header is incorrect."
           "\n  Expected to find end of header marker (26) after version, found [%d]", 
           filename.c_str(), t );
      }

   // Verify 'mode' is valid
   BYTE rd3 = rgBuf[7];
   switch( rd3 )
      {
   case 0:  // mode = 'IndividualMajor' or RowMajor
      layout = LayoutGroupGenotypesByIndividual;      // all SNPs per individual are sequential in memory
      cbStride = (cSnps * sizeof( float ) );          // 4 genotypes per byte so round up
      break;
   case 1:  // mode = 'SnpMajor' or ColumnMajor
      layout = LayoutGroupGenotypesBySnp;             // all individuals per SNP are sequential in memory
      cbStride = (cIndividuals * sizeof( float ) );   // 4 genotypes per byte so round up
      break;
   default:
      Fatal( "Ill-formed .BIG file [%s].  BIG file header is incorrect.  Expected mode to be 0 or 1, found %d", filename.c_str(), rd3 );
      break;
      }

   int rdCSnps = *(int*)(rgBuf+8);
   int rdCIndividuals = *(int*)(rgBuf+12);

   if ( static_cast< size_t >( rdCSnps ) != cSnps_ )
      {
      Fatal( "BIG file [%s] formated differently than expected."
         "\n  Expected %d SNPs, found %d SNPs",
         filename.c_str(), cSnps_, rdCSnps );
      }
   if ( static_cast< size_t >( rdCIndividuals ) != cIndividuals_ )
      {
      Fatal( "BIG file [%s] formated differently than expected."
         "\n  Expected %d individuals, found %d individuals",
         filename.c_str(), cIndividuals_, rdCIndividuals );
      }
   }

CPlinkBigFile::~CPlinkBigFile()
   {
   if ( pFile )
      {
      fclose( pFile );
      pFile = nullptr;
      }
   }

int CPlinkBigFile::NextChar()
   {
   int value = fgetc( pFile );
   if ( value == EOF )
      {
      Fatal( "Ill-formed BIG file [%s].  Encountered EOF before exepected.", filename.c_str() );
      }
   return( (unsigned char)value );
   }

float CPlinkBigFile::ReadFloat()
   {
   float value = std::numeric_limits<float>::quiet_NaN();
   Read( (BYTE*)&value, sizeof( value ) );
   return( value );
   }

size_t CPlinkBigFile::Read( BYTE *pb, size_t cbToRead )
   {
   size_t cbRead = fread( pb, 1, cbToRead, pFile );
   if ( cbRead != cbToRead )
      {
      if ( feof( pFile ) )
         {
         Fatal( "Encountered EOF before exepected in BIG file. Ill-formed BIG file [%s]", filename.c_str() );
         }
      int err = ferror( pFile );
      if ( err )
         {
         Fatal( "Encountered a file error %d in BIG file [%s]", err, filename.c_str() );
         }
      }
   return( cbRead );
   }

size_t CPlinkBigFile::ReadLine( BYTE *pb, size_t idx )
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

bool CPlinkBigFile::WriteFile( string& filename, size_t cIndividuals, size_t cSnps, LayoutMode layoutMode, real *pdata )
   {
   FILE *pFile;
   char rgBuf[32] =
            { 
            bigFileMagic1, bigFileMagic2,
            '0', '1', '0', '0',
            0x1a,
            01,
            00,00,00,00, 
            00,00,00,00, 
            00,00,00,00, 
            00,00,00,00, 
            00,00,00,00, 
            00,00,00,00
            };
   if ( sizeof( rgBuf ) != 32 )
      {
      Fatal( "Compile error.  rgBuf must be 32 bytes in size." );
      }

   pFile = fopen( filename.c_str(), "wb" );
   if ( pFile == nullptr )
      {
      return( false );
      }

   // prepare and write the header
   rgBuf[7] = layoutMode;
   *(int*)(&rgBuf[8]) = (int)cSnps;
   *(int*)(&rgBuf[12]) = (int)cIndividuals;
   fwrite( rgBuf, sizeof( rgBuf ), (size_t)1, pFile );

   if ( layoutMode == LayoutGroupGenotypesBySnp )
      {
      float* rgf = new float[ cIndividuals ];
      for ( size_t iSnp=0; iSnp<cSnps; ++iSnp )
         {
         for ( size_t iIndividual=0; iIndividual<cIndividuals; ++iIndividual )
            {
            size_t idx = ColumnMajorIndex( cIndividuals, cSnps, iIndividual, iSnp );
            rgf[ iIndividual ] = (float)pdata[ idx ];
            }
         fwrite( rgf, sizeof( float ), cIndividuals, pFile );
         }
      fclose( pFile );
      delete [] rgf;
      }
   else
      {
      float* rgf = new float[ cSnps ];
      for ( size_t iIndividual=0; iIndividual<cIndividuals; ++iIndividual )
         {
         for ( size_t iSnp=0; iSnp<cSnps; ++iSnp )
            {
            size_t idx = ColumnMajorIndex( cIndividuals, cSnps, iIndividual, iSnp );
            rgf[ iIndividual ] = (float)pdata[ idx ];
            }
         fwrite( rgf, sizeof( float ), cIndividuals, pFile );
         }
      fclose( pFile );
      delete [] rgf;
      }
   return( true );
   }
