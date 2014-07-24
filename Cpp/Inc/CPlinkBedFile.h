#if !defined( CPlinkBedFile_h )
#define CPlinkBedFile_h
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
 * CPlinkBedFile - {PLINK Bed File Access Class}
 *
 *         File Name:   CPlinkBedFile.h
 *           Version:   1.00
 *            Author:   
 *     Creation Date:    4 Dec 2011
 *     Revision Date:    4 Dec 2011
 *
 *    Module Purpose:   This file defines the CPlinkBedFile class 
 *                         for FastLmmC
 *
 *                      A .BED file contains compressed binary genotype values for 
 *                         for individuals by SNPs.  For FastLmm, we prefer and may
 *                         require the file LayoutMode be LayoutGroupGenotypesBySnp
 *
 *                      The .bed header is three bytes followed immediately by data.
 *                      
 *                         bedFileMagic1 | bedFileMagic2 | LayoutMode
 *                         [... data ... ]
 *
 *    Change History:   
 *
 * Test Files: 
 */

const BYTE bedFileMagic1 = 0x6C;       // 0b01101100 or 'l' (lowercase 'L')
const BYTE bedFileMagic2 = 0x1B;       // 0b00011011 or <esc>

enum LayoutMode
   {
   LayoutUnknown     = -1
  ,LayoutRowMajor    = 0                  // all elements of a row are sequential in memory
  ,LayoutColumnMajor = 1                  // all elements of a colomn are sequential in memory
  ,LayoutGroupGenotypesByIndividual = 0   // all SNP genotypes for a specific individual are seqential in memory
  ,LayoutGroupGenotypesBySnp = 1          // all Individual's genotypes for a specific SNP are sequential in memory
   };

enum BedGenotype     // integer representation of genotype values in Plink's binary .BED file
   {
   bedHomozygousMinor = 0
  ,bedMissingGenotype = 1
  ,bedHeterozygous    = 2
  ,bedHomozygousMajor = 3
   };

extern real unknownOrMissing;  // now used by SnpInfo
extern real homozygousPrimaryAllele;
extern real heterozygousAllele;
extern real homozygousSecondaryAllele;
extern real mapBedGenotypeToRealAllele[4];

class CBedFile
   {
public:
   CBedFile();
   ~CBedFile();

   void Open( const string& filename_, size_t cIndividuals_, size_t cSnps_ );                // validate the file matches the extents we expect
   void Open( const string& filename_, CPlinkFamFile *pFamFile_, CPlinkBimFile *pBimFile_ ); // these three files all work together to form one unit

   // return the layout mode of the file.  
   //   We work better with files that are LayoutGroupGenotypesBySnp
   LayoutMode  GetLayoutMode();

   // return the compressed length of one line of SNP data (related to cIndividuals)
   // TODO:  MAKE THIS PRIVATE!  Consumers should be using cIndividuals or cSnps
   size_t      CbStride() { return( cbStride ); }

   // return the filename associated with this CBedFile
   const string& Filename() { return( filename ); }

   // read the data for one SNP (idxSnp) into the BYTE buffer pb
   size_t   ReadLine( BYTE *pb, size_t idxSnp );

   // read the genotype for all the individuals in 'list' at the SNP specified by iSNP
   void     ReadGenotypes( size_t iSnp, const vector< string >& iidList, vector< real >& pvOut );
   void     ReadGenotypes( size_t iSnp, const vector< size_t >& iIndividualList, vector< real >& pvOut );

private:
   int      NextChar();
   size_t   Read( BYTE *pb, size_t cbToRead );

   static const size_t   cbHeader = 3;         // 
   string   filename;
   FILE     *pFile;
   vector< BYTE > rgBytes;
   vector< BedGenotype > rgBedGenotypes;
   
   LayoutMode  layout;        // 0=RowMajor(all snps per individual together);
                              // 1=ColumnMajor(all individuals per SNP together in memory)
   size_t   cIndividuals;
   size_t   cSnps;
   size_t   cbStride;

//private:
   CPlinkFamFile *pFamFile;
   CPlinkBimFile *pBimFile;
   };

#endif      // CPlinkBedFile_h
