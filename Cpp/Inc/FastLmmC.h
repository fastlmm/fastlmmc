#if !defined( FastLmmC_h )
#define FastLmmC_h
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
 * FastLmmC - {FastLmmC Primary Header}
 *
 *         File Name:   FastLmmC.h
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *    Module Purpose:   This file defines or includes everything for 
 *                      components of FastLmmC
 *    Change History:   
 *
 * Test Files: 
 */
#if defined( _MSC_VER )
#define  COMPILE_FOR_WINDOWS     // using Visual C compiler
#elif defined( __GNUC__ )
#define  COMPILE_FOR_LINUX       // using gcc/g++ to compiler for Linux
#else
#error We only support VisualC++ on Windows and GCC on Linux.  Unrecognized compiler environment.  Stop Compilation
#endif

/*
 * Get our system includes
 */
#if defined( COMPILE_FOR_WINDOWS )
#define  NOMINMAX
#include <windows.h>

#define  _CRTDBG_MAP_ALLOC_NEW   // Define for heap debugging of new (broken, but...the best we can do for now)
#define  _CRTDBG_MAP_ALLOC       // Define for heap debugging of malloc
#include <stdlib.h>              // MUST BE IN THIS ORDER FOR crtdbg.h to work properly
#include <crtdbg.h>              // To use heap debugging routines effectively, look for documentation on
                                 //   _CrtMemCheckpoint()
                                 //   _CrtMemDifference()
                                 //   _CrtMemDumpStatistics()
                                 //   _CrtMemDumpAllObjectsSince()
                                 //   _CrtDumpMemoryLeaks()
#elif defined( COMPILE_FOR_LINUX )
#define  nullptr        NULL           // gcc does not yet support nullptr
typedef  unsigned char  BYTE;          // common Windows types not found in gcc
typedef  long long      __int64;
typedef  long long      LARGE_INTEGER;
#include <stdlib.h>
#else
#error Unrecognized compiler environment.  Stop Compilation
#endif

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#if defined( COMPILE_FOR_WINDOWS )
#include <io.h>                  // for access()
#elif defined( COMPILE_FOR_LINUX )
#include <unistd.h>              // for access() on POSIX/Linux
#else 
#error Unrecognized compiler environment.  Stop Compilation
#endif

#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#if defined( USE_RANDOM )
#include <random>
#endif

#if defined( USE_FLOAT )
   typedef float real;
#else
   typedef double real;
#endif

using namespace std;
using std::list;
using std::map;
using std::set;
using std::string;
using std::vector;
using std::min;

#include "CPhenData.h"        // PhenotypeValue typedef for IndividualInfo

struct IndividualInfo         // TODO:  Replace FamRecord
   {
   string idFamily;           // can be alpha-numeric / or numeric... :-(
   string idIndividual;       // 
   string idPaternal;
   string idMaternal;
   int    sex;                // 1=male, 2=female, other=unknown
   Phenotype phenotype;       // -9=missing, 0=missing, 1=unaffected, 2=affected, or real (anything other than 0,1,2,-9)
   IndividualInfo() { sex = 0; phenotype = 0.0; };
   };

struct FamRecord              // or struct IndividualInfo... ?
   {
   // used by Fam and TFam
   string idFamily;           // can be alpha-numeric / or numeric... but no spaces
   string idIndividual;       // idIndiviual+idFamily must be unique!
   string idPaternal;         // must match another idIndividual ? (what is Family relationship?)
   string idMaternal;         // must match another idIndividual ? (what is Family relationship?)
   int    sex;                // 1=male, 2=female, other=unknown
   real   phenotype;          // -9=missing, 0=missing, 1=unaffected, 2=affected, or real (anything other than 0,1,2,-9)

   FamRecord() {sex=0;phenotype=0.0;}
   };

//global variables and constants we use
const real _PI = (real)2.0*acos((real)0.0); 
const real _log2pi = log((real)2.0*_PI);

const real real0 = static_cast< real >(0.0);
const real real1 = static_cast< real >(1.0);
const real real2 = static_cast< real >(2.0);
const real maxDeltaLL = static_cast< real >(1.0e-10);   // warn when absDeltaLL > maxDeltaLL 

const size_t NoMatch = (size_t)-1;

#include "Utils.h"               // need some definitions process before BlasLapack.h
#include "BlasLapack.h"          // BlasLapack.h isolates FastLmmC's usage of advanced math 
                                 //  routines from Intel's MKL and AMD's ACML libraries.  If 
                                 //  you need to use an alternative BLAS/LAPACK library, 
                                 //  this is where to look.

#include "Startup.h"
#include "BrentC.h"
#include "CLexer.h"
#include "TsvLexer.h"

#include "CSnpData.h"
#include "CPhenData.h"
#include "CKernelData.h"
#include "CCovariatesData.h"

#include "CSnpFilterOptions.h"
#include "CSnpInfo.h"
#include "CPlinkBimFile.h"
#include "CPlinkDatFile.h"       // must be before CPlinkFile.h
#include "CPlinkFamFile.h"
#include "CPlinkMapFile.h"       // must be before CPlinkFile.h
#include "CPlinkBedFile.h"       // define .BED file processing
#include "CPlinkBigFile.h"       // define .BIG file processing
#include "CPlinkAlternatePhenotypeFile.h" // must be before CPlinkFile.h
#include "CPlinkFile.h"

#include "GroupSim.h"
#include "StudyData.h"

#include "EigenSym.h"
#include "MatrixIO.h"
#include "FastLMM.h"             // must be after EigenSym.h
#include "FastLMM_fullrank.h"    // must follow FastLMM.h
#include "FastLMM_lowrank.h"     // must follow FastLMM.h
#include "FastLMM_woodbury.h"    // must follow FastLMM.h
#include "Gamma.h"
#include "Beta.h"
#include "FisherF.h"
#include "Qval.h"
#include "LmmGwas.h"
#include "MathFunctions.h"
#include "LinReg.h"
#include "CLogisticRegression.h"
#include "AutoSelectSnps.h"

#if defined( COMPILE_FOR_WINDOWS )          // using Visual C compiler
/*
 * TODO: gcc v4.1 does not support OMP.
 *       Re-enable OpenMP on Linux after we upgrade the 
 *       Linux gcc compiler to v4.2 or greater.
 */
#include <omp.h>                 
#endif

// Declare our function prototypes

/*
 * Parameters: 
 *    real *SNPs - pointer to output matrix of reals [nIndividuals x nSNPs]
 *                 The matrix stored in column-major order. 
 *                 NaNs will be set to 0.0 in the result.
 *                 This matrix will receive the result of the computations.
 *    const size_t nIndividuals - number of individual samples along the row axis
 *    const size_t nSNPs - number of SNPs along the column axis
 *    const bool unitVariance - apply imputation as per Price et al. 2007 or not
 *                 true:   apply the imputation as per Price et al. 2007
 *                 false:  only zeroMean the SNPs
 *    zeroMean:    true:   zero mean the SNPs
 *                 false:  do not zero Mean the SNPs (leaves them always positive)
 *    const bool assumeHardyWeinberg - SNPs are encoded under Hardy-Weinber model
 *                 true:   SNPs are encoded as 0,1,2 and unitVariance under Hardy-Weinber
 *                 false:  use empirical estimator of standard deviation
 *    const bool useUnbiasedEstimatorOrBetaBinomial - 
 *                 true:   under HardyWeinberg, this add both one positive and one negative pseudo individual in the Beta-Binomial model for the SNP frequency. 
 *                         for empirical std, this will use the unbiased estimator., 
 *                 false: No pseudocounts/biased estimator for std
 */
void ImputeAndZeroMeanSNPs( const vector<SnpInfo>& snpInfo, real *SNPs, const size_t nIndividuals, const size_t nSNPs, const bool unitVariance, const bool zeroMean, const NormalizationInfo& normalizationInfo_, const bool useUnbiasedEstimatorOrBetaBinomial );

/*
 * Parameters:
 *    real *SNPs - pointer to output matrix of reals [nIndividuals x nSNPs]
 *                 The matrix stored in column-major order. 
 *                 This matrix will receive the result of the computations.
 *    const size_t nIndividuals - number of individual samples along the row axis
 *    const size_t nSNPs - number of SNPs along the column axis
 *    const size_t nSNPsdivide - ???
 */
void preprocessSNPsForBuildingCovariance( real *SNPs, const size_t nIndividuals, const size_t nSNPs, const size_t nSNPsdivide );
 
/*
 *  Compute the pvalue using Likelihood-ratio test as long as
 *    deltaLL is positive.  Otherwise return 1.0
 *    Issue a warning if deltaLL is too far from an expected
 *    value 
 */
real PvalFromLikelihoodRatioTest( real deltaLL, real dof_ );

/*
 *  more function prototypes w/out documentation
 */
void LoadDataAltPhen( StudyData *studydata );
void permuteStudyData( EigenSym *eigen, StudyData *data );
void permuteKernel( EigenSym *eigen, StudyData *data );
void permuteCovariates( StudyData *data );
void permutePhen( StudyData *data );
void replaceSNPsByMissingness( real *matrix, const size_t nRows, const size_t nColumns );
void zeroMeanColumns( real *matrix, const size_t nRows, const size_t nColumns, bool meanImputeNaN );

/*
 *  Debug Dump enums and routines
 */
enum LogFileWriteMode
   {
   OverWrite = 0,
   Append = 1
   };

void DumpLmmGwasInputs( const StudyData *pStudyData, const EigenSym *pEigenSym, const ResultsLMM *pResultsLmm, LogFileWriteMode fOpenAppend=Append );
void DumpStudyData( const StudyData *pStudyData, const string& comment, LogFileWriteMode fOpenAppend=Append );
void DumpEigenSymState( const EigenSym *pEigenSym, const string& comment, LogFileWriteMode fOpenAppend=Append );
void DumpKernelArray( const real *pReal, size_t cRows, size_t cColumns, const string& comment, LogFileWriteMode fOpenAppend=Append );

void DumpStringVector( FILE* pf, int indent, const vector<string>& v, const string& name );
void DumpSnpInfoVector( FILE *pf, int indent, const vector<SnpInfo>& v, const string& name );
void DumpRealArray( FILE *pf, int indent, const real *p, const string& name, size_t rows, size_t columns );
void Dumpreals( FILE *pf, const real *pdbl, size_t count );
void DumpGroupSim( FILE *pf, int indent, const GroupSim* p, const string& name );

void LogMemoryState( const string& comment, LogFileWriteMode fOpenAppend=Append );

#endif   // FastLmmC_h
