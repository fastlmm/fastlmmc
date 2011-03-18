#if !defined( CPlinkFile_h )
#define CPlinkFile_h
/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * CPlinkFile - {PLINK File Access Class}
 *
 *         File Name:   CPlinkFile.h
 *           Version:   1.00
 *            Author:   BobD
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

/*
 * 'Publish' our class declarations / function prototypes
 */
struct GENOTYPE
   {
   char allele1[2];
   };

struct Ped_Record
   {
   string idFamily;           // alpha-numeric
   string idIndividual;
   string idPaternal;
   string idMaternal;
   int    sex;                // 1=male, 2=female, other=unknown
   int    cGenotype;           // 
   double phenotypeTrait;     // -9=missing, 0=missing, 1=unaffected, 2=affected, or double (anything other than 0,1,2,-9)
   GENOTYPE* rgGenotype;

   Ped_Record();
   ~Ped_Record();
   };

struct Map_Record
   {
   string idChromosome;
   string idSnp;
   unsigned int geneticDistance;
   unsigned int basepairPosition;
   };

class CPlinkFile
   {
public:
   CPlinkFile();
   ~CPlinkFile();
   void ReadFile( string& filename );

private:
   string plinkFilename;
   vector<Ped_Record> rgPed;
   vector<Map_Record> rgMap;

   void ReadPedFile();
   void ReadMapFile();
   };

#endif   // CPlinkFile_h
