/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * CPlinkFile - {PLINK File Access Class}
 *
 *         File Name:   CPlinkFile.cpp
 *           Version:   1.00
 *            Author:   BobD
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements the CPlinkFile class for FastLmmC
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Include Files
 */
#include "FastLmmC.h"
#include "CPlinkFile.h"

Ped_Record::Ped_Record()
   {
   sex = 0;
   cGenotype = 0;
   rgGenotype = nullptr;
   }

Ped_Record::~Ped_Record()
   {
   if ( rgGenotype )
      {
      delete[] rgGenotype;
      }

   cGenotype = 0;
   }

void CPlinkFile::ReadFile( string& readfile )
   {
   plinkFilename = readfile;
   ReadPedFile();
   ReadMapFile();
   }

void CPlinkFile::ReadPedFile()
   {
   string filename = plinkFilename + ".PED";

   }

void CPlinkFile::ReadMapFile()
   {
   string filename = plinkFilename + ".MAP";
   }

