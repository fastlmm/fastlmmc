/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * FastLmmC - {Fast Linear Mixed Model in C++}
 *
 *         File Name:   FastLmmC.cpp
 *           Version:   1.00
 *            Author:   BobD
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements the FastLMM tool to C++
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Include Files
 */
#include "FastLmmC.h"

/*
 * The code
 */
int main(int argc, char* argv[])
   {
   extern char *szSplash;

   InitializeProgram();
   fprintf( stderr, szSplash, szFastLmmCVersion );
   ProcessCommandLine( argc, argv );

   bool dotiming = false;

   if ( dotiming )
      {
      LmmGwas_oldfortiming();
      }
   else
      {
      LmmGwas();
      }

   fprintf( stderr, "\n" );
   return( 0 );
   }

