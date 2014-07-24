#if !defined( CPlinkBigFile_h )
#define CPlinkBigFile_h
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
 *         File Name:   CPlinkBigFile.h
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   24 Apr 2012
 *     Revision Date:   24 Apr 2012
 *
 *    Module Purpose:   This file defines the CPlinkBigFile class 
 *                         for FastLmmC
 *
 *                      A .BIG file contains binary genotype values for 
 *                         for individuals by SNPs. In this case a single
 *                         4 byte floating point value (float) in the 
 *                         range 0.0 to 2.0 is stored
 *
 *    Change History:   
 *
 * Test Files: 
 */

const BYTE bigFileMagic1 = 'b';        // 0b01100010 or 'b' (lowercase 'L')
const BYTE bigFileMagic2 = 'd';        // 0b01100100 or 'd' <esc>

#if 0
enum LayoutMode
   {
   RowMajor=0                    // all elements of a row are sequential in memory
  ,ColumnMajor=1                 // all elements of a colomn are sequential in memory
  ,GroupGenotypesByIndividual=0  // all SNP genotypes for a specific individual are seqential in memory
  ,GroupGenotypesBySnp=1         // all Individual's genotypes for a specific SNP are sequential in memory
   };
#endif

class CPlinkBigFile
   {
public:
   //CPlinkBigFile( const string& filename_ );
   CPlinkBigFile( const string& filename_, size_t cIndividuals_, size_t cSnps_ );
   ~CPlinkBigFile();

   int      NextChar();
   size_t   Read( BYTE *pb, size_t cbToRead );
   size_t   ReadLine( BYTE *pb, size_t idx );
   float    ReadFloat();
   static bool WriteFile( string& filename, size_t cIndividuals, size_t cSnps, LayoutMode mode, real *pdata );

//private:
   static const size_t   cbHeader = 32;      // 
   string   filename;
   FILE     *pFile;
   
   int      version;
   int      release;
   LayoutMode  layout;        // 0=RowMajor(all snps per individual together);
                              // 1=ColumnMajor(all individuals per SNP together in memory)
   size_t   cIndividuals;
   size_t   cSnps;
   size_t   cbStride;

   void     Init( const string& filename_, size_t cIndividuals_, size_t cSnps_ );
   };

#endif      // CPlinkBigFile_h
