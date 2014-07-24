#if !defined( CPhenData_h )
#define CPhenData_h
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
 * CPhenData - {Phenotype Data Class}
 *
 *         File Name:   CPhenData.h
 *           Version:   1.00
 *            Author:   
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
typedef real   Phenotype;        // Phenotypes are 'real' values
                                 // -9=missing, 0=missing, 1=unaffected, 2=affected, or real (anything other than 0,1,2,-9)
enum PhenotypeType               // The type of a phenotype indicates how to treat the phenotype 'value'
   {
   PHENOTYPE_UNINITIALIZED = 0
  ,PHENOTYPE_AFFECTION     = 1          // A phenotype where 1.0 = unaffected and 2.0 = affected
  ,PHENOTYPE_QUANTITATIVE  = 2          // A phenotype with a real value representing 'weighting'
   };

/*
 * 'Publish' our class declarations / function prototypes
 */
class CPhenData
   {
public:
   vector<string> rowLabels;
   vector<string> columnLabels;
   vector<vector<real> > phenVectors;       // temporary storage of Phenotype reads

   size_t         cRows;
   size_t         cColumns;
   map<string, size_t> columnMap;      // map column to column index
   map<string, size_t> rowMap;         // map column to column index

   real         *phenArray;

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
