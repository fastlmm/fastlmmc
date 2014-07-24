#if !defined( CPlinkFamFile_h )
#define CPlinkFamFile_h
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
 *         File Name:   CPlinkFamFile.h
 *           Version:   1.00
 *            Author:   
 *     Creation Date:    4 Dec 2011
 *     Revision Date:    4 Dec 2011
 *
 *    Module Purpose:   This file defines the CPlinkFamFile class 
 *                         for FastLmmC
 *
 *                      A .FAM file contains records for individuals
 *                         covering family and phenotype values.
 *
 *    Change History:   
 *
 * Test Files: 
 */

class CPlinkFamFile
   {
public:
   CPlinkFamFile();
   CPlinkFamFile( const string& filename_ );
   ~CPlinkFamFile();

   void     Load();
   void     Load( const string& filename_ );

   size_t   CountOfIndividuals() { size_t rc = rgIndividualInfo.size(); return( rc ); };
   size_t   IndexFromKeyIndividual( const string& idKey );
   IndividualInfo* PIndividualInfo( const string& idKey );
   IndividualInfo* PIndividualInfo( size_t idxIndividualInfo );

private:
   PhenotypeType phenotypeType;
   string filename;
   vector< IndividualInfo > rgIndividualInfo;
   map< string, size_t > idKeyToIndividualInfoIndex;
   };
#endif      // CPlinkFamFile_h
