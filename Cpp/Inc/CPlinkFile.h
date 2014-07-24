#if !defined( CPlinkFile_h )
#define CPlinkFile_h
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
 * CPlinkFile - {PLINK File Access Class}
 *
 *         File Name:   CPlinkFile.h
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file declare the CPlinkFile class for FastLmmC
 *
 *    Change History:   
 *
 */

/*
 * 'Publish' our defines
 */
extern real unknownOrMissing;

/*
 * 'Publish' our class declarations / function prototypes
 */
enum Genotype        // integer representation of genotype
   {
   missingGenotype = -1
  ,homozygousMajor =  0
  ,heterozygous    =  1
  ,homozygousMinor =  2
   };

typedef float  BigGenotype;   // BIG file genotype is float in the range of 0.00 - 2.00
                              // std::numeric_limits<float>::quiet_NaN() represents missing
struct PedRecord
   {
#if 0    //bd     // TODO: refactor to IndividualInfo (through phenotype)
   IndividualInfo ii;
#else
   string idFamily;           // can be alpha-numeric / or numeric... :-(
   string idIndividual;
   string idPaternal;
   string idMaternal;
   int    sex;                // 1=male, 2=female, other=unknown
   real   phenotype;          // -9=missing, 0=missing, 1=unaffected, 2=affected, or real (anything other than 0,1,2,-9)
#endif

   vector<SnpNucleotides> rgSnps;        //

   PedRecord() { sex = 0; phenotype = 0.0; };
#if 0
   ~PedRecord();
#endif
   };

struct TPedRecord
   {  // TPed is the same as MapRecord for first 4 elements with the snp data for all individuals appended to each SNP
   SnpInfo snpInfo;           // the TPedRecord conains the SNP description information plus all the SNP data
   vector<SnpNucleotides> rgSnps;
   };

struct DosageRecord           // information from dosage file (.DAT)
   {  // .DAT file contains two pieces of snpInfo then the dosage inforamtion for each individual
   SnpInfo  snpInfo;
   vector<SnpDosageData> rgSnpDosageData;
   };

// TODO:  CPlinkLexer.h Requires the SnpNucleotides definition at the top of this file...
//          fix this later
#include "CPlinkLexer.h"      

string KeyFromIdFamilyAndIdIndividual( const string& idFamily, const string& idIndividual );

class CPlinkFile
   {
public:
   CPlinkFile();
   CPlinkFile( const string& baseFilename_ );
   CPlinkFile( const string& baseFilename_, const AlternatePhenotypeParameters& alternatePhenoptypeParameters_);
   CPlinkFile( const string& basefilename_, const AlternatePhenotypeParameters& alternatePhenotypeParameters_, SnpFilterOptions& snpOptions_ );
   ~CPlinkFile();

   void ReadPlinkFiles( const string& basefilename_, PlinkFileType filetype_, const AlternatePhenotypeParameters& alternatePhenotypeParameters_, SnpFilterOptions& snpOptions_ );
   void ReadPlinkFiles( const string& basefilename_, PlinkFileType filetype_, const AlternatePhenotypeParameters& alternatePhenotypeParameters_ );
   void ReadPlinkFiles( const string& basefilename_, PlinkFileType filetype_ );
   void ReadPlinkFiles( PlinkFileType filetype );

   size_t   cIndividuals;
   size_t   cPhenotypes;
   size_t   cSnps;
   vector<string>  individualLabels;   // labels for individuals selected
   vector<string>  phenotypeLabels;    // labels for phenotypes selected
   vector<SnpInfo> rgSnpInfo;          // Meta-data Infomation for snps selected (rg means rangeof or vectorof)
   size_t   ind_phenotype;

   real     *snpData;                  // snpData[cIndividuals, cSnps] (column major)
   real     *phenotypeData;            // phenotypeData[cIndividuals, cPhenotypes] (column major)

   CCovariatesData *covariatesFileData;

   static inline int compareSnpInfoByPosition(const void *x, const void *y);   //if position (chromosome + basepairposition) of (x>y): 1, (x<y):-1, (x==y):0
   static inline int getBpDistance(SnpInfo &x, SnpInfo &y);       //compute BP distance. returns a large number if on diffrent chromosomes
   static inline real getGeneticDistance(SnpInfo &x, SnpInfo &y); //compute genetic distance. returns infinity if on different chromosomes

private:
   PlinkFileType fileType;             // Natural=use PED/MAP, Binary=use BED/FAM/BIM, Transposed=use TPED/TFAM
   string   baseFilename;
   AlternatePhenotypeParameters  alternatePhenotypeParameters;
   CPlinkAlternatePhenotypeFile  alternatePhenotypeFile;

   SnpFilterOptions snpFilter;         // some PLink files are filtered on SNPs too.

   PhenotypeType  phenotypeType;       // 0=uninitialized, 1=affection phenotype, 2=quantitative phenotype
   size_t   cIndividualsRead;
   size_t   cSnpsRead;                 // number of SNPs read from data files

   vector< size_t > selectedIndividuals;  // index of individuals selected for inclusion in study
   vector< size_t > selectedSnps;         // index of SNPs selected for inclusion in this study

   real   *phenArray;
   real   *snpArray;
//   real   *genotypeArray;                 //

   vector<PedRecord> rgPed;
   vector<MapRecord> rgMap;               // The Map file contains _only_ SnpInfo records
   vector<FamRecord> rgFam;
   //vector<BedRecord> rgBed;             // The Bed file is a clone of rgMap
   //vector<BimRecord> rgBim;
   vector<TPedRecord> rgTPed;
   //vector<TFamRecord> rgTFam;           // The TFam file is clone of rgFam
   vector<DosageRecord> rgDosage;         // paired with rgFam for individual info & alternate phenotype for phenotype info

   vector< vector<real> > rgBinaryGenotype;

   map< string, size_t > indexOfFamilyIdAndIndividualId;
   size_t   IndexFromIdFamilyAndIdIndividual( const string& idFamily, const string& idIndividual );
   void     AddToFamilyIdAndIndividualIdIndex( const string& idFamily, const string& idIndividual, size_t idx );

   map< string, size_t > indexOfSnpIds;
   size_t   IndexFromSnpId( const string& snpId );
   void     AddSnpIdToIndex( const string& snpId, size_t idx );

   void     ReadAltPhenotype( CPlinkAlternatePhenotypeFile& apf, AlternatePhenotypeParameters& altParameters );

   void     ReadNaturalFiles();
   void     ReadPedFile();
   void     WritePedFile( const string& pedFile );
   void     ReadMapFile();
   void     WriteMapFile( const string& mapFile );
   void     ProduceFilteredDataFromNaturalFiles();

#if 0
   void     ReadBinaryFiles();
#endif
//   void     ReadBedFile();

#if 0
   void     ReadBinaryFiles2();
   void     ReadBedFile2();

   void     ReadBinaryFiles3();
   void     ReadBedFile3( vector< size_t >& selectedSnps );
#endif

   void     ReadBinaryFiles4();
   void     ReadBedFile4( vector< size_t >& preSelectedSnps, vector< size_t >& finalSelectedSnps );

//   void     WriteBedFile( const string& bedFile, LayoutMode layout_ );

   void     ReadFamFile();
   void     ReadFamFile( const string& famFile );
   void     WriteFamFile( const string& famFile );
   void     ReadBimFile();
   void     WriteBinFile();
#if 0
   void     ProduceFilteredDataFromBinaryFiles();
#endif
//   void     ProduceFilteredDataFromBinaryFiles3();

   void     ReadTransposedFiles();
   void     ReadTFamFile();
   //void     WriteTFamFile( const string& tfamFile );
   void     ReadTPedFile();
   void     WriteTPedFile( const string& tmapFile );
   void     ProduceFilteredDataFromTransposedFiles();

   void     ReadDosageFiles();
   void     ProduceFilteredDataFromDosageFiles( CPlinkDatFile& datFile, CPlinkMapFile& mapFile );

   void     ReadImputedGenotypeFiles();
   void     ReadBigFile( vector< size_t >& preSelectedSnps, vector< size_t >& finalSelectedSnps );

   void     ComputeNaturalSnpAlleleChars();
   void     CreateNaturalColumnMajorSnpArray();
   void     CreateNaturalColumnMajorPhenArray();

//   void     CreateBinaryColumnMajorSnpArray();
//   void     CreateBinaryColumnMajorSnpArray3();
   void     CreateBinaryColumnMajorPhenArray();

   void     CreateTransposedColumnMajorSnpArray();
   void     CreateTransposedColumnMajorPhenArray();

   void     CreateDosageColumnMajorSnpArray( CPlinkDatFile& datFile );
   void     CreateDosageColumnMajorPhenArray( CPlinkDatFile& datFile );

   void     SelectIndividualsFromRgPed();
   void     SelectIndividualsFromRgFam();
   void     SelectIndividualsFromDosageFiles( CPlinkAlternatePhenotypeFile& altPhenotypeFile, CPlinkDatFile& datFile );

   void     SelectSnps( size_t cSnpsRead, SnpFilterOptions& snpFilter );
   void     SelectSnps( size_t cSnpsRead, SnpFilterOptions& snpFilter, vector< size_t >& selectSnps );
   void     SelectSnpsFromNaturalFiles();
//   void     SelectSnpsFromBinaryFiles();
//   void     SelectSnpsFromBinaryFiles3();
//   void     SelectSnpsFromBinaryFiles3( size_t cSnpsRead, SnpFilterOptions& snpFilter, vector< size_t >& selectedSnps );
   void     SelectSnpsFromTransposedFiles();
   void     SelectSnpsFromDosageFiles( CPlinkDatFile& datFile );

   void     ExtractSnpsFromRgMap();
   void     ExtractSnpsFromRgMap2();      // ExtractSnpsForNaturalFiles()
   void     ExtractSnpsFromRgTPed();
   void     ExtractSnpsFromDosage( CPlinkDatFile& datFile, CPlinkMapFile& mapFile );

   void     DoInit( const string* basefilename, const AlternatePhenotypeParameters *alternatePhenoptypeParameters, SnpFilterOptions *snpFilerOptions );
   real     RealFromSnpNuculeotides( char& majorAllele_, SnpNucleotides& snp );
   real     RealFromSnpDosageData( SnpDosageData& snpDD );
   void     FreePrivateMemory();

   void     WriteDatFile( const string& datFilename );

   bool     FSnpHasVariation( size_t iSnp );
//   bool     FSnpHasVariationBinaryFile( size_t iSnp );
//   bool     FSnpHasVariation3( SnpInfo& snpInfo, vector< BedGenotype >& rgBedGenotypes );
   bool     FSnpHasVariation4( SnpInfo& snpInfo, vector< BedGenotype >& rgBedGenotypes, vector< size_t >& selected );
   bool     FSnpHasVariation( SnpInfo& snpInfo, vector< BigGenotype >& rgBigGenotypes, vector< size_t >& selected );
   bool     FSnpHasVariation4b( SnpInfo& snpInfo, real* genotypes );
   bool     FSnpHasVariationDosageFile( size_t iSnp );
   bool     FSnpHasVariationTransposedFile( size_t iSnp );

   void     ValidateSnpGeneticDistanceInformation( SnpInfo& snpInfo_ );
   void     ValidateSnpPositionInformation( SnpInfo& snpInfo_ );
   };

// compute BP distance. returns a large number if on diffrent chromosomes
inline int CPlinkFile::getBpDistance( SnpInfo &x, SnpInfo &y )
   {
   if ( x.iChromosome == y.iChromosome )
      {
      return( x.basepairPosition - y.basepairPosition );
      }
   else if ( x.iChromosome < y.iChromosome )
      {
      return( -2000000000 );
      }
   else
      {
      return( 2000000000 );
      }
   }

// compute genetic distance. returns infinity if on different chromosomes
inline real CPlinkFile::getGeneticDistance( SnpInfo &x, SnpInfo &y )
   {
   if ( x.iChromosome == y.iChromosome )
      {
      return( x.geneticDistance - y.geneticDistance );
      }
   else if ( x.iChromosome < y.iChromosome )
      {
      return( -std::numeric_limits<real>::infinity() );
      }
   else
      {
      return( std::numeric_limits<real>::infinity() );
      }
   }

inline int CPlinkFile::compareSnpInfoByPosition( const void *x, const void *y )   //if position (chromosome + basepairposition) of (x>y): 1, (x<y):-1, (x==y):0
   {
   //check the chromosome
   SnpInfo* px = (SnpInfo*)x;
   SnpInfo* py = (SnpInfo*)y;

   if ( px->iChromosome > py->iChromosome )
      {
      return( 1 );
      }
   else if( px->iChromosome < py->iChromosome)
      {
      return( -1 );
      }
   else
      {//if the chromosomes match up, compare the position
      if ( px->basepairPosition > py->basepairPosition )
         {
         return( 1 );
         }
      else if( px->basepairPosition < py->basepairPosition )
         {
         return( -1 );
         }
      else
         {
         return( 0 );
         }
      }
   }

/* 
 * ComputeSnpInfoStatistics() fills in SnpInfo SNP statistics for 
 *   this SnpInfo using cSnpValues values pointed to bypSnpValues
 */
void ComputeSnpInfoStatistics( SnpInfo& snpInfo, size_t cSnpValues, real *pSnpValues );

#endif   // CPlinkFile_h
