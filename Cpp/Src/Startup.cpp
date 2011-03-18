/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * FastLmmC - {Fast Linear Mixed Model in C++}
 *
 *         File Name:   Startup.cpp
 *           Version:   1.00
 *            Author:   BobD
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements the FastLMM startup code
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
 * Externs we need from Splash.cpp
 */
extern int wVerMajor;
extern int wVerMinor;
extern int wVerBuild;
extern int wVerPatch;
extern char *szSplash;

/*
 * Globals we define
 */
char szComputerName[ MAX_COMPUTERNAME_LENGTH + 1 ] = "";
char szFastLmmCVersion[ 32 ] = "";
int  wFastLmmCVersion = 0;
long long i64FastLmmCFullVersion = (long long)0;

string snpDataFile;                 // Program options and parameters
string phenDataFile;
string outputFile;
string kernelFile;
string snpDataCovFile;
string eigenDir;

int    numberOptimizationStarts;    // brent grid search options/parameters   
int    brentMaxIter;                //    default values are set in ProcessingCommandLine.
double minLogVal;
double maxLogVal;
double brentTol;

COV_LEARN_TYPE covLearnType;
PARAM_LEARNING_TYPE mlType;
COMPUTE_EIGEN_TYPE computeEigenType;
RUN_GWAS_TYPE runGwasType; 

// define our debug flags too
int fCreateReadValidationFiles = 0;

/*
 * The code
 */
void Usage()
   {
   fprintf( stderr, "\n"
                    "\nDescription:"
                    "\n  FastLmmC computes Linear Mixed Model GWAS"
                    "\n"
                    "\nUsage:"
                    "\n  FastLmmC -option <parameter>"
                    "\n    -snpDataFile filename        : file containing the SNP information"
                    "\n    -snpDataCovFile filename     : file containing the SNP information for building the Covariance matrix"
                    "\n    -phenDataFile filename       : file containing the phenotype information"
                    "\n    -kernelFile filename         : file containing the Kernel information"
                    "\n    -OutputFile filename         : file to write output into [OutputFile.txt]"
                    //"\n    -rank [int k]                : [full]"
                    "\n    -CovLearnType.[Full/Once]    : ?? [Full]"
                    "\n    -ParamLearningType.[ML/REML] : ?? [ML] REML not yet implemented"
                    "\n    -Brent.Starts                : numberOptimizationStarts [100]"
                    "\n    -Brent.MaxIter               : brentMaxIter [100000]"
                    "\n    -Brent.MinLogVal             : minLogVal [-5.0]"
                    "\n    -Brent.MaxLogVal             : maxLogVal [ 5.0]"
                    "\n    -Brent.Tol                   : brentTol [1e-6]"
                    "\n    -runGwasType.[Run/NoRun]     : Run GWAS or exit after computing EigenSym [RUN]"
                    "\n    -computeEigenType.[Save/Load] directoryname"
                    "\n                                 : Save/Load the EigenSym object to/from the named directory"
                    "\n  Debug Options:"
                    "\n    -CreateValidationFiles       : Create duplicate files to validate reading data correctly"
                    "\n" );
   exit(1);
   }

void ProcessCommandLine( int argc, char* argv[] )
   {
	char *pch;
   int iArg;

   // set our 'default' values
   covLearnType = ONCE;
   mlType = ML;
   numberOptimizationStarts = 100;
   brentMaxIter = 100000;
   minLogVal = -5.0;
   maxLogVal =  5.0;
   brentTol = 1e-6;
   computeEigenType = NOSAVE;
   runGwasType = RUN; 
   
   // REVIEW:  set debug options to default values 
   //          set fVerbose during development
   fVerbose = 1;
   fCreateReadValidationFiles = 0;

	for ( iArg=1; iArg<argc; ++iArg )
		{
		pch = argv[ iArg ];
		if ( (*pch != '-') && (*pch != '/') )
			{
			break;
			}

		++pch;
		if ( (_strcmpi( pch, "?" ) == 0 ) || (_strcmpi( pch, "help" ) == 0) )
			{
			Usage();
			}
      else if ( (_strcmpi( pch, "verbose" ) == 0 ))
         {
         fVerbose = 1;
         }
      else if ( (_strcmpi( pch, "createvalidationfiles" ) == 0 ))
         {
         fCreateReadValidationFiles = 1;
         }
		else if ( (_strcmpi( pch, "snpdatafile" ) == 0 ) )
			{
         // get the snpDataFile name
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            snpDataFile = argv[ iArg ];
				}
			}
		else if ( (_strcmpi( pch, "phendatafile" ) == 0 ) )
			{
         // get the snpDataFile name
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            phenDataFile = argv[ iArg ];
				}
			}
      else if ( (_strcmpi( pch, "snpdatacovfile" ) == 0 ) )
			{
         // get the snpDataCovFile name
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            snpDataCovFile = argv[ iArg ];
				}
			}
		else if ( (_strcmpi( pch, "outputfile" ) == 0 ) )
			{
         // get the outputFile name
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            outputFile = argv[ iArg ];
				}
			}
		else if ( (_strcmpi( pch, "kernelfile" ) == 0 ) )
			{
         // get the kernelFile name
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            kernelFile = argv[ iArg ];
				}
			}
		else if ( (_strcmpi( pch, "covlearntype.full" ) == 0 ) )
			{
         // get the COV_LEARN_TYPE (and we only accept FULL or ONCE for now)
         covLearnType = FULL;
			}
      else if ( (_strcmpi( pch, "covlearntype.once" ) == 0 ) )
			{
         // get the COV_LEARN_TYPE (and we only accept FULL or ONCE for now)
         covLearnType = ONCE;
			}
		else if ( (_strcmpi( pch, "paramlearningtype.ml" ) == 0 ) )
			{
         // get the PARAM_LEARNING_TYPE (and we only accept ML for now)
         mlType = ML;
         }
      //runGwasType: default = RUN
      else if ( (_strcmpi( pch, "rungwastype.run" ) == 0 ) )
			{
         runGwasType= RUN;
         }
      else if ( (_strcmpi( pch, "rungwastype.norun" ) == 0 ) )
			{
         runGwasType= NORUN;
         }
      //computeEigenType: default = NOSAVE
      else if ( (_strcmpi( pch, "computeeigentype.nosave" ) == 0 ) )
			{
         computeEigenType= NOSAVE;
         }
      else if ( (_strcmpi( pch, "computeeigentype.save" ) == 0 )
             || (_strcmpi( pch, "computeeigentype.savetodisk" ) == 0 )
              )
			{
         computeEigenType= SAVETODISK;
         // get the eigenDir name
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            eigenDir = argv[ iArg ];
				}
			}
      else if ( (_strcmpi( pch, "computeeigentype.load" ) == 0 ) 
             || (_strcmpi( pch, "computeeigentype.loadfromdisk" ) == 0 )
              )
			{
         computeEigenType= LOADFROMDISK;
		   
         // get the eigenDir name
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            eigenDir = argv[ iArg ];
				}	
         }

      // Look for Brent parameters
		else if ( (_strcmpi( pch, "brent.starts" ) == 0 ) )
			{
         // get the numberOptimizationStarts
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            numberOptimizationStarts = atoi( argv[ iArg ] );
            }
         }
		else if ( (_strcmpi( pch, "brent.minlogval" ) == 0 ) )
			{
         // get the minLogVal
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            minLogVal = atof( argv[ iArg ] );
            }
         }
		else if ( (_strcmpi( pch, "brent.maxlogval" ) == 0 ) )
			{
         // get the maxLogVal
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            maxLogVal = atof( argv[ iArg ] );
            }
         }
		else if ( (_strcmpi( pch, "brent.tol" ) == 0 ) )
			{
         // get the brentTol
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            brentTol = atof( argv[ iArg ] );
            }
         }
		else if ( (_strcmpi( pch, "brent.maxiter" ) == 0 ) )
			{
         // get the brentMaxIter
			if ( (++iArg >= argc) )
				{
				Usage();
				}
			else
				{
            brentMaxIter = atoi( argv[ iArg ] );
            }
         }
      // Oops, don't know this option... Report the error
		else
			{
			printf( "\nError:  Unknown option. [%s]", argv[ iArg ] );
			Usage();
			}
      }

   // now verify that all the proper arguments and options have been read and are consistent.
      if ( snpDataFile.length() == 0 && runGwasType==RUN)
         {
         Error( "missing SnpDataFile information." );
         Usage();
         }
      if ( phenDataFile.length() == 0  && runGwasType==RUN)
         {
         Error( "missing PhenDataFile information." );
         Usage();
         }
      if ( outputFile.length() == 0  && runGwasType==RUN)
         {
         Error( "missing OutputFile information." );
         Usage();
         }
      if ( kernelFile.length() == 0 && snpDataCovFile.length() == 0  && computeEigenType!=LOADFROMDISK  )
         {
         Error( "missing KernelFile information. Cannot compute covariance matrix" );
         Usage();
         }
      if ( eigenDir.length() == 0 && computeEigenType!=NOSAVE  )
         {
         Error( "missing computeEigenType.[SAVETODISK/LOADFROMDISK] directory information. Cannot save/load Eigenvectors" );
         Usage();
         }
   }

void InitializeProgram()
   {
   // Build up the 'correct' version string
	if ( wVerBuild )
		{
      sprintf( szFastLmmCVersion, "v%d.%02d.%08d", wVerMajor, wVerMinor, wVerBuild );
		}
	else
		{
   	sprintf( szFastLmmCVersion, "v%d.%02d", wVerMajor, wVerMinor );
		}

   // create the integer representations for the current version
   wFastLmmCVersion = (wVerMajor * 100) + wVerMinor;
	i64FastLmmCFullVersion = (wFastLmmCVersion * (long long)100000000) + wVerBuild;

   // get the host machine name
	unsigned long cchSzComputerName = sizeof( szComputerName );
	if ( !GetComputerName( szComputerName, &cchSzComputerName ) )
		{
		Error( "Unable to retrieve this computer's name...and I want it!" );
		NTLastError();
		}
   }
