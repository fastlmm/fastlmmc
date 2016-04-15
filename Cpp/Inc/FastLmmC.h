#if !defined( FastLmmC_h )
#define FastLmmC_h
/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * Utils - {Utility Routines}
 *
 *         File Name:   Utils.h
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
 * Get our system includes
 */
#if defined( _MSC_VER )          // using Visual C compiler
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#if defined( _MSC_VER )          // using Visual C compiler
#include <io.h>                  // for access()
#else
#include <unistd.h>              // for access() on POSIX/Linux
#endif

#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().

#include <string>
#include <map>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>

using namespace std;
using std::map;
using std::string;
using std::vector;
using std::list;
//using std::min;

#include "Utils.h"
#include "Startup.h"
#include "BrentC.h"
#include "CLexer.h"
#include "TsvLexer.h"
#include "CSnpData.h"
#include "CPhenData.h"
#include "CKernelData.h"
#include "EigenSym.h"
#include "MatrixIO.h"
#include "FastLMM.h"             // must be after EigenSym.h
#include "FastLMM_fullrank.h"    // must follow FastLMM.h
#include "FastLMM_lowrank.h"     // must follow FastLMM.h
#include "Gamma.h"
#include "Qval.h"
#include "LmmGwas.h"
#include "MathFunctions.h"

#include "BlasLapack.h"          // BlasLapack.h isolates FastLmmC's usage of advanced math 
                                 //  routinesfrom Intel's MKL and AMD's ACML libraries.  If 
                                 //  you needto use an alternative BLAS/LAPACK library, 
                                 //  this is where to look.

// Declare our function prototypes

// Parameters: SNPs [nIndividuals by nSNPs] Matrix stored in column-major order. This will hold the result.
// NaNs will be set to 0.0 in the result.
// maxValue: maximum value a SNP can have (usually value for minor-minor)
// unitVariance: true: apply the imputation as per Price et al. 2007, false:only zeroMean the SNPs
void ImputeAndZeroMeanSNPs(double *SNPs, const int nIndividuals, const int nSNPs, const double maxValue, const bool unitVariance );

// Parameters: SNPs [nIndividuals by nSNPs] Matrix stored in column-major order. This will hold the result.
void preprocessSNPsForBuildingCovariance( double *SNPs, const int nIndividuals, const int nSNPs );

void WriteResultFile( string& resultFile, vector<string> snpLabels,string phenLabel, int ncovariates, int *index, double *LL_alt, double LL_0, double *pvals,double *beta, double *sigg2, double *sige2, double *delta, double *beta_0, double sigg2_0, double sige2_0, double delta_0, double *qvals);

//global variable
const double _PI = 2.0*acos(0.0); 
const double _log2pi=log(2.0*_PI);

#endif   // FastLmmC_h
