#if !defined( CKernelData_h )
#define CKernelData_h
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
 * CKernelData - {Kernel Data Class}
 *
 *         File Name:   CKernelData.h
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file declares the CKernelData class for FastLmmC
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
class CKernelData
   {
public:
   vector<string> columnLabels;
   vector<string> rowLabels;
   vector<vector<real> > kernelVectors;        // temporary storage of Kernel reads
   
   size_t         cColumns;
   size_t         cRows;
   map<string, size_t> columnMap;      // map column to column index
   map<string, size_t> rowMap;         // map column to column index

   real         *kernelArray;

   CKernelData();
   ~CKernelData();
   void ReadKernelFile( const string& kernelFile );
   void WriteKernelFile( const string& kernelFile );
   void CreateColumnMajorArrayFromVectors();
   void CreateColumnMajorArrayFromVectorsAndMapping( const vector<string>& columnLabelOrder );
   };

/*
 * 'Publish' the globals we define
 */

#endif   // CKernelData_h
