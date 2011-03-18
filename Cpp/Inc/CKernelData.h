#if !defined( CKernelData_h )
#define CKernelData_h
/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * CKernelData - {Kernel Data Class}
 *
 *         File Name:   CKernelData.h
 *           Version:   1.00
 *            Author:   BobD
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
   vector<vector<double> > kernelVectors;        // temporary storage of Kernel reads
   
   size_t         cColumns;
   size_t         cRows;
   map<string, size_t> columnMap;      // map column to column index
   map<string, size_t> rowMap;         // map column to column index

   double         *kernelArray;

   CKernelData();
   ~CKernelData();
   void ReadKernelFile( string& kernelFile );
   void WriteKernelFile( string& kernelFile );
   void CreateColumnMajorArrayFromVectors();
   void CreateColumnMajorArrayFromVectorsAndMapping( vector<string> columnLabelOrder );
   };

/*
 * 'Publish' the globals we define
 */

#endif   // CKernelData_h
