#if !defined( CPhenData_h )
#define CPhenData_h
/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * CPhenData - {Phenotype Data Class}
 *
 *         File Name:   CPhenData.h
 *           Version:   1.00
 *            Author:   BobD
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file declares the CPhenData class for FastLmmC
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
class CPhenData
   {
public:
   vector<string> rowLabels;
   vector<string> columnLabels;
   vector<vector<double> > phenVectors;       // temporary storage of Phenotype reads

   size_t         cRows;
   size_t         cColumns;
   map<string, size_t> columnMap;      // map column to column index
   map<string, size_t> rowMap;         // map column to column index

   double         *phenArray;

   CPhenData();
   ~CPhenData();
   void ReadPhenFile( string& phenFile );
   void WritePhenFile( string& phenFile );
   void CreateColumnMajorArrayFromVectors();
   void CreateColumnMajorArrayFromVectorsAndMapping( vector<string> columnLabelOrder );
   };


/*
 * 'Publish' the globals we define
 */

#endif   // CPhenData_h
