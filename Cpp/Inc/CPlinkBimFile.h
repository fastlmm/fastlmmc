#if !defined( CPlinkBimFile_h )
#define CPlinkBimFile_h
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
 *         File Name:   CPlinkBimFile.h
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

class CPlinkBimFile
   {
public:
   CPlinkBimFile();
   CPlinkBimFile( const string& filename_ );
   ~CPlinkBimFile();

   void     Load();
   void     Load( const string& filename_ );

   size_t   CountOfSnpInfo() { return( rgSnpInfo.size() ); }
   size_t   IndexFromIdSnp( const string& idSnp );
   SnpInfo* PSnpInfo( const string& idSnp );
   SnpInfo* PSnpInfo( size_t idxSnpInfo );

//private:
   string   filename;
   vector< SnpInfo >     rgSnpInfo;
   map< string, size_t > idSnpToSnpInfoIndex;
   };

#endif      // CPlinkBimFile_h
