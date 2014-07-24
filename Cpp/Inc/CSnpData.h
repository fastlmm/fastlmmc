#if !defined( CSnpData_h )
#define CSnpData_h
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
 * CSnpData - {SNP Data Class}
 *
 *         File Name:   CSnpData.h
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file declares the CSnpData class for FastLmmC
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Get our includes
 */

/*
 * 'Publish' our defines
 */

/*
 * 'Publish' our class declarations / function prototypes
 */
class CSnpData
   {
public:
   vector<string> rowLabels;
   vector<string> columnLabels;
   vector<string> snpVectors;      // temporary storage of SNP reads

   size_t         cRows;
   size_t         cColumns;
   map<string, size_t> columnMap;      // map column to column index
   map<string, size_t> rowMap;         // map column to column index

   real         *snpArray;

   CSnpData();
   ~CSnpData();
   void ReadSnpFile( string& snpFile );
   void WriteSnpFile( string& snpFile );
   void CreateColumnMajorArrayFromVectors();
   void CreateColumnMajorArrayFromVectorsAndMapping( vector<string> columnLabelOrder );

private:
   real realFromChar( char ch );
   };


/*
 * 'Publish' the globals we define
 */

#endif   // CSnpData_h
