#if !defined( CSnpInfo_h )
#define CSnpInfo_h
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

struct SnpNucleotides         // what two nucleotides occur here
   {
   char alleles[2];

   SnpNucleotides() {Clear();}
   void Clear() {alleles[0] = alleles[1] = '\0';}
   };

/*
 * Dosage Data has 3 input formats with 1, 2, or 3 real values associated 
 *   We currently support format 1 (Dosage) and format 2 (Probabilities)
 *   As they are the same type of inputs, we use a union to overlay two 
 *   options them in the same location so we can use meaninful names.
 *   Since it is a union, they cannot both be active at the same time.
 */
extern real unknownOrMissing;

struct SnpDosageData
   {
   union
      {
      real dosage;
      real probabilityOfHomozygousMinor;
      };
   real probabilityOfHeterozygous;

   SnpDosageData() {Clear();}
   SnpDosageData( real dosage_ )
      { dosage = dosage_; probabilityOfHeterozygous = 0.0; }
   SnpDosageData( real homozygousMinor, real heterozygous )
      { probabilityOfHomozygousMinor = homozygousMinor; probabilityOfHeterozygous = heterozygous; }
   void Clear()
      { probabilityOfHomozygousMinor = probabilityOfHeterozygous = 0.0; }
   static real RealFromSnpDosageData( SnpDosageData& snpDD, PlinkFileType fileType_ )
      {
      // Return the real representing the dosage
      //  If we are using Dosage1 format, return the single value
      //  If we are using Dosage2 format, either missing or (2 * probability of homozygousMminor) + (probability of heterozygous)
      //
      if ( fileType_ == Dosage1 )
         {
         return( snpDD.dosage );
         }

      // requires a bit of processing to get value when in Dosage2 format
      real rc;
      if ( (snpDD.probabilityOfHomozygousMinor == missingGenotypeValue) || (snpDD.probabilityOfHeterozygous == missingGenotypeValue) )
         {
         rc = unknownOrMissing;
         }
      else
         {
         rc = (snpDD.probabilityOfHomozygousMinor * 2) + snpDD.probabilityOfHeterozygous;
         }
      return( rc );
      }
   };

inline bool operator==( const SnpDosageData &lhs, const SnpDosageData& rhs )
   {
   if ( (lhs.probabilityOfHomozygousMinor == rhs.probabilityOfHomozygousMinor) 
     && (lhs.probabilityOfHeterozygous == rhs.probabilityOfHeterozygous) )
      {
      return( true );
      }
   return( false );
   }

inline bool operator!=( const SnpDosageData &lhs, const SnpDosageData& rhs )
   {
   return( ! (lhs == rhs) );
   }

struct SnpInfo                // information that describes a SNP
   {
   int    iChromosome;        // 1-22, X, Y, XY, MT, or 0 if unplaced
   string idChromosome;
   string idSnp;              // rs# or snp identifier
   real   geneticDistance;    // morgans
   int    basepairPosition;   // bp units
   char   majorAllele;        // A2 when reading/writing
   char   minorAllele;        // A1 when reading/writing

   size_t cSamples;           // number of valid SNP values from individuals
   real   meanValue;          // mean or average SNP value
   real   stdDev;             // sample standard deviation for SNP values

   SnpInfo()    {iChromosome=basepairPosition=0;cSamples=0;geneticDistance=meanValue=stdDev=0.0;majorAllele=minorAllele='\0';}
   void Clear() {iChromosome=basepairPosition=0;cSamples=0;geneticDistance=meanValue=stdDev=0.0;majorAllele=minorAllele='\0';idChromosome.clear();idSnp.clear();}
   real StdDeviation() { return( stdDev ); }
   };

// Helper routines
inline size_t IndexFromList( const vector< SnpInfo >& rgSnpInfo, const string& snpId )
   {
   for ( size_t i=0; i<rgSnpInfo.size(); ++i )
      {
      if ( rgSnpInfo[ i ].idSnp == snpId )
         {
         return( i );
         }
      }
   return( -1 );
   }


#endif   // CSnpInfo_h
