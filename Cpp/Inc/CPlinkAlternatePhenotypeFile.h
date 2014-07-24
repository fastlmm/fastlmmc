#pragma once
#if !defined( CPlinkAlternatePhenotypFile_h )
#define CPlinkAlternatePhenotypFile_h
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
 * CPlinkAlternatePhenotypeFile - {PLINK Alternate Phenotype File Class}
 *
 *         File Name:   CPlinkAlternatePhenotypeFile.h
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file declare the CPlinkAlternatePhenotypeFile
 *                      class for FastLmmC
 *
 *    Change History:   
 *
 */


/*
 * 'Publish' our defines
 */

/*
 * 'Publish' our class declarations / function prototypes
 */

struct AlternatePhenotypeParameters
   {
   int    fUseAlternatePhenotype;            // Using PLINK alternative phenotype file?
   size_t selectedPhenotypeIndex;
   string selectedPhenotypeId;
   string alternatePhenotypeFilename;

   AlternatePhenotypeParameters() {fUseAlternatePhenotype=0; selectedPhenotypeIndex=0;}
   };

struct AlternatePhenotypeRecord
   {
   string idFamily;
   string idIndividual;
   vector<real> phenotypes;
   };

class CPlinkAlternatePhenotypeFile
   {
public:
   CPlinkAlternatePhenotypeFile( const string& alternatePhenotypeFilename );
   CPlinkAlternatePhenotypeFile();
   ~CPlinkAlternatePhenotypeFile();

   real    PhenotypeValue( const string& idFamily, const string& idIndividual, const string& phenotypeLabel );
   real    PhenotypeValue( const string& idFamily, const string& idIndividual, size_t phenotypeIndex );
   real    PhenotypeValue( const string& idFamily, const string& idIndividual );
   real    PhenotypeValue( const string& individualKey, size_t phenotypeIndex );
   real    PhenotypeValue( size_t individualIndex, size_t phenotypeIndex );
   real    PhenotypeValue( const string& individualKey );
   real    PhenotypeValue( size_t individualIndex );

   size_t  IndexFromPhenotypeLabel( const string& phenotypeLabel );
   size_t  IndexFromIdFamilyAndIdIndividual( const string& idFamily, const string& idIndividual );
   size_t  IndexFromIndividualKey( const string& key );
   void    KeyFromIndividualIndex( size_t idx, string& keyOut );
   string  KeyFromIndividualIndex( size_t idx );

   size_t  CountOfPhenotypes()  {return( cPhenotypes ); }
   size_t  CountOfIndividuals() {return( cIndividuals ); }

   void    Read();
   void    Read( const string& alternatePhenotypeFilename );
   void    SetSelectedPhenotype( string& idPhenotype );
   void    SetSelectedPhenotype( size_t idx );
   size_t  SelectedPhenotypeIndex() { return( selectedPhenotype ); }
   string& SelectedPhenotypeName();
   string& PhenotypeName( size_t idx );

   vector< string > KeysForIndividualsWithValidPhenotypes( size_t idx );
   vector< string > KeysForIndividualsWithValidPhenotypes( const string& idPhenotype );

private:
   string  filename;
   map< string, size_t > phenotypeIndividualIdIndex;
   map< string, size_t > phenotypeLabelIndex;

   vector< string > phenotypeNames;
   vector< AlternatePhenotypeRecord > records;
   size_t  cIndividuals;         // count of individuals in alternate phenotype file
   size_t  cPhenotypes;          // count of Phenotypes per individual
   size_t  selectedPhenotype;    // the 'selected' phenotype is the one we are studying this run

   void Init() {cIndividuals=cPhenotypes=selectedPhenotype=0;}
   };

#endif
