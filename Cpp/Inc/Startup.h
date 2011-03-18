#if !defined( Startup_h )
#define Startup_h
/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * Utils - {Utility Routines}
 *
 *         File Name:   Startup.h 
 *           Version:   1.00
 *            Author:   BobD
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file declares a number of utility routines for console applications
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * 'Publish' our data structures
 */
enum COV_LEARN_TYPE
   {
   FULL = 0
   ,ONCE = 1
   ,ALT_PARTIAL = 2
   };

enum PARAM_LEARNING_TYPE
   {
   ML = 0
   //,REML = 1
   };

enum COMPUTE_EIGEN_TYPE
   {
   NOSAVE=0
   ,LOADFROMDISK=1
   ,SAVETODISK=2
   };

enum RUN_GWAS_TYPE
   {
   RUN=0
   ,NORUN=1
   };

/*
 * 'Publish' our function prototypes
 */
void InitializeProgram();
void ProcessCommandLine( int argc, char* argv[] );
void Usage();

/*
 * 'Publish' the globals we define
 */
extern char szFastLmmCVersion[];
extern int  wFastLmmCVersion;
extern long long i64FastLmmCFullVersion;

extern string snpDataFile;                // The program options that we use.
extern string phenDataFile;
extern string outputFile;
extern string kernelFile;
extern string snpDataCovFile;
extern string eigenDir;

extern int    numberOptimizationStarts;   // brent grid search options/parameters
extern int    brentMaxIter;
extern double minLogVal;
extern double maxLogVal;
extern double brentTol;

extern COV_LEARN_TYPE covLearnType;
extern COMPUTE_EIGEN_TYPE computeEigenType;
extern RUN_GWAS_TYPE runGwasType; 
extern PARAM_LEARNING_TYPE mlType;

// debug and conditional reporting flags
extern int fCreateReadValidationFiles;

#endif   // Startup_h
