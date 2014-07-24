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
 * Startup - {Start the program by processing the command line options}
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

// TODO: obsolete?  remove in future since we only support INDIVIDUAL_ONLY
INDIVIDUAL_ID_TYPE individualIdType = INDIVIDUALID_ONLY;

char szComputerName[ MAX_COMPUTERNAME_LENGTH + 1 ] = "";
char szFastLmmCVersion[ 32 ] = "";
int  wFastLmmCVersion = 0;
long long i64FastLmmCFullVersion = 0ll;

enum  PlinkFileType plinkFileType;
enum  PlinkFileType plinkFileTypeCovariance;

AlternatePhenotypeParameters alternatePhenotypeParameters;

string missingPhenotypeString = "-9";  // default is -9
real   missingPhenotypeValue  = -9.0;
string missingGenotypeString  = "-9";  // default is -9
real   missingGenotypeValue   = -9.0;
HANDLE_NON_PSD_TYPE handle_non_psd = SET_TO_ZERO;

// Program options and parameters
string covariatesFile;              // name of file containing Covariate information

string plinkFilenameRoot;           // Use PLINK Ped/Map/Bed/Fam/Bim files for Test data
vector< string > plinkFilenameRootCovarianceVector;   // Use PLINK Ped/Map files for building the covariance matrix
int    indPheno = 0;                // index of selected phenotype
int    setRank = -1;
PERMUTATION_TYPE permutationMode = NO_PERM;
bool   linReg = false;                    // use linear regression to solve
bool   fUseLogisticRegression = false;    // use logistic regression to solve

string outputFilename;              // store the output filename as given by user or default from input file
string fullOutputFilename;          // keep the full path to the output file around
bool   fVerboseOutput = false;      // Use to output more details in final report

string kernelFile;
string kernelFileOut;               //save kernel to text file
string snpDataCovFile;              //TODO: NOT USED!!!  This is not set anywhere an is always "".  Remove it?
string eigenDir;

/*
 * Define our SnpPair options
 */
bool fDoSnpPairs = false;           // Signal SnpPairs option in effect - compute epistatic effects
bool fDoSnpId1 = false;             // Signal SnpId1 option in effect - restrict epistatic study to 1 snp x N SNPs
string snpId1Global;                // the snpId1
bool fEmitSnpIndexOnly = false;
bool fClusterIO = false;

/*
 * Define our low memory options
 */
bool fDoIteratingFastLmm = false;      // triggered by -bfilesim-
size_t cIteratingBlockSizeMax = 20000; // default block size for iterating

/*
 *  Define our AutoSelect options
 */
bool   fAutoSelectSnps = false;              // preprocess the files to determine the 'best' snps to use
int    autoSelectFoldCount = 10;             // default AutoSelect fold count value = 10
bool   fAutoSelectApplyThinning = false;     // apply thinning/filtering to SearchValues
int    autoSelectMaxSizeForThinning = 2000;  // Restrict thinning to apply only to search values below 2000
real   autoSelectRegionDistance = 0.50;      // width of 'filter' region on each of SNP (+/- 0.5 centiMorgans) (geneticDistance)
AutoSelectCriterion autoSelectCriterion = useLL;   // default autoSelectCriterion = useLL or useMSE
string autoSelectOutputFilenameRoot;         // root of names for AutoSelect output files.
vector< int > autoSelectSearchValues;        // k values to search to yield the highest LL
int    defaultAutoSelectSearchValues[] =
                                    {
                                      0,  1,  2,  3,  4,  5,  6,  7,  8, 9,
                                         10, 20, 30, 40, 50, 60, 70, 80, 90,
                                    100,125,160,200,250,320,400,500,630,800,1000,
                                    -1 // sentinel value
                                    };
bool   selectTopK = false;             // Use memory partitioning during AutoSelect
real   mFraction = 0.20;               // percentMemoryToUse in AutoSelect
size_t topk;                           // 
bool   fWritePedFile = false;          // Signal we are to produce a PED file from AutoSelect
string pedOutputFile = "AutoSelected"; // Root name of plink .PED file produced by AutoSelect
bool   zeroMeanTestSNPs = true;        // subtract mean from test SNPs? (default true)
bool   unitVarianceTestSNPs = true;    // divide SNPs by their standard deviation? (default true)
bool   zeroMeanSimilaritySNPs = true;     // 
bool   unitVarianceSimilaritySNPs = true; //   

ExcludeSnpParams excludeSnpParams;

bool   forceFullRank = false;
bool   identityKernel = false;
int    numberOptimizationStartsGlobal; // brent grid search options/parameters   
int    brentMaxIterGlobal;             //    default values are set in ProcessingCommandLine.
NormalizationInfo normalizationInfo;   // default to 'Normal' normalization
bool   usePseudocounts = false;        // use pseudocounts (Beta-Binomial model)  for computing SNP frequency. default is false.
bool   useBiasTerm;              //use a bias term in the regression instead of zero-meaning the phenotype
bool   useWoodbury;
bool   use_f_test;               //compute an F-test instead of LR test?
bool   simMissingness = false;   //build a Similarity matrix from missingness pattern in SNPs default is false
bool   useHeritability = false;  //by default do optimization in logdelta space

string groupsFile;
string simGroups;

bool   runAlternativeModels = true;

real   minLogValGlobal = -5.0;      // ensure 'reasonable' initialization values
real   maxLogValGlobal =  5.0;
real   brentTolGlobal  = (real)1e-6;
real   logDeltaPreset = std::numeric_limits< real >::quiet_NaN();    //if set, this will be used instead of optimization

COV_LEARN_TYPE covLearnType;
PARAM_LEARNING_TYPE mlType;
COMPUTE_EIGEN_TYPE computeEigenType;
RUN_GWAS_TYPE runGwasType; 
TREAT_MISSING missingType;
KERNEL_TYPE kernelType;

// Process / store the SnpFilterOptions for the data and genetic similarity matrices
SnpFilterOptions snpFilterOpts;
SnpFilterOptions simFilterOpts;

bool fMaximizeWorkingSet = true;             // Allow consumption of most of the system ram
bool fSetFastLmmRealsPrintPrecision = false; // use custom value for output precision? (default: 16)
int  fastLmmRealsPrintPrecision = 16;        // storage to hold the current reals print precision
bool fNoDosageRangeCheck = false;            // turn off validation of dosage values
real pvaluePrintThreshold = 1.0;             // p-values must be less than pvaluePrintThreshold to be output
bool fEmitSorted = true;                     // write output sorted by pvalue? (default: true)

// Setup for Random Number Generation (RNG)
bool           fRngSeedSuppliedByUser = false;  //
unsigned long  randomNumberSeed = 1;            // initialize to our 'default' seed value

#if defined( USE_RANDOM )     // use the C++0x <random> template library
default_random_engine               rng;
uniform_real_distribution<double>   drand( 0.0, 1.0 );
#endif

// define our debug flags too
bool   fCreateReadValidationFiles = false;
bool   fWriteLogFile = false;
string logDir = ".\\Logs";
bool   fDumpSnpFilterOptions = false;
bool   fLogMemoryDiagnostics = false;
bool   fUseReproducibleFP = false;           // Set the math library to minimize FP variation (MKL only for now)

/*
 * The code
 */

/*
 *  SetMlType and associated variables together
 */
void SetMlType( PARAM_LEARNING_TYPE ml_ )
   {
   if ( ml_ == ML )
      {
      mlType = ML;
      use_f_test  = false;
      useBiasTerm = true;
      //TODO:  All instances set useWoodbury = false, and it does
      //         not have all the current features integrated into
      //         it.  Thus the woodbury path should be removed to 
      //         cleanup and simplify the source code.
      //TODO: The SnpPair path requires useWoodbury to be true as 
      //         it is only implemented along the woodbury path.
      //         Before we can remove woodbury, we need to add
      //         snpPair support to the non-woodbury path.
      //         See next function SetDoSnpPairs()
      useWoodbury = false;
      }
   else if ( ml_ == REML )
      {
      // REML is used with an f-test and with a bias term, otherwise the extra terms are affected.
      mlType = REML;
      use_f_test  = true;
      useBiasTerm = true;
      useWoodbury = false;
      }
   else
      {
      Fatal( "Bad parameter to SetMlType().  Expected [%d] or [%d].  Found [%d]", ML, REML, ml_ );
      }
   }

void SetDoSnpPairs()
   {
   fDoSnpPairs = true;
   SetMlType( ML );

   // useWoodbury is true in ML only for snpPairs
   //
   //TODO:  SnpPair support should be updated to not require
   //       useWoodbury as it does not have all the features
   //       implemented along the woodbury path
   useWoodbury = true;
   
   covLearnType = ONCE;
   }

void Usage()
   {
   fprintf( stderr, "\n"
                    "\nDescription:"
                    "\n  FastLmmC computes Linear Mixed Model GWAS"
                    "\n"
                    "\nUsage:"
                    "\n  FastLmmC -option <parameter>"
                    "\n"
                    "\n  -file basefilename"
                    "\n       basename for PLINK's .map and .ped files"
                    "\n"
                    "\n  -bfile basefilename"
                    "\n       basename for PLINK's binary .fam, .bim, and .bed files"
                    "\n"
                    "\n  -tfile basefilename"
                    "\n       basename for PLINK's transposed .tfam and .tped files"
                    "\n"
                    "\n  -dosage1 basefilename"
                    "\n  -dosage2 basefilename"
                    "\n       Dosage supports PLINK dosage file formats 1 and 2.  Use this"
                    "\n       option to define the appropriate format to read and the basename"
                    "\n       for PLINK's dosage .dat, .fam, and optional .map files"
                    "\n"
                    "\n  -afile basefilename [** In Development **]"
                    "\n       basename for imputed genotype files .fam, .bim, and .big"
                    "\n       These files are similar format to PLINK binary files above."
                    "\n"
                    "\n  -pheno filename"
                    "\n       name of phenotype file"
                    "\n"
                    "\n  -mpheno <int>"
                    "\n       Specifies the index for phenotype in -pheno file to process,"
                    "\n       starting at 1 for the first phenotype column. Cannot be used"
                    "\n       together with -pheno-name."
                    "\n       Default: 1."
                    "\n"
                    "\n  -pheno-name phenotypeName"
                    "\n       phenotype name for phenotype in -pheno file to process.  If"
                    "\n       this option is used, the phenotype name must be specified in"
                    "\n       the header row.  Cannot be used together with -mpheno."
                    "\n"
                    "\n  -NoZeroMeanTestSNPs"
                    "\n       Do not subtract the mean from SNPs tested during pre-processing"
                    "\n       Default: false"
                    "\n"
                    "\n  -NoUnitVarTestSNPs"
                    "\n       Do not divide SNPs tested by the standard deviation during"
                    "\n       pre-processing"
                    "\n       Default: false"
                    "\n"
                    "\n  -numjobs <int>"
                    "\n       Partion the SNPS into <int> groups and run FaSTLMM on the"
                    "\n       partition specified by -thisjob"
                    "\n       Requires thisjob to be set too."
                    "\n"
                    "\n  -thisjob <int>"
                    "\n       Specifies which partition of SNPS created by -numjobs to"
                    "\n       process with FaSTLMM"
                    "\n       Requires numjobs to be set too."
                    "\n"
                    "\n  -extract filename"
                    "\n       This is a SNP filter option.  FaSTLMM will only analyze the"
                    "\n       SNPs explicitly listed in the 'filename'"
                    "\n"
                    "\n  -extractSim filename"
                    "\n       This is a genetic similarity SNP filter option.  FaSTLMM will"
                    "\n       only use SNPs explicitly listed in the 'filename' for computing"
                    "\n       genetic similarity."
                    "\n"
                    "\n  -extractSimTopK filename <int>"
                    "\n       Similar to -extractSim, this is a genetic similarity SNP filter"
                    "\n       option.  FaSTLMM will only use the first <int> SNPs explicitly"
                    "\n       listed in the 'filename' for computing genetic similarity."
                    "\n"
                    "\n  -fileSim basefilename"
                    "\n       basename for PLINK's .map and .ped files for computing"
                    "\n       genetic similarity"
                    "\n"
                    "\n  -bfileSim basefilename"
                    "\n       basename for PLINK's binary .fam, .bim, and .bed files for"
                    "\n       building genetic similarity"
                    "\n"
                    "\n  -tfileSim basefilename"
                    "\n       basename for PLINK's transposed .tfam and .tped files for"
                    "\n       building genetic similarity"
                    "\n"
                    "\n  -dosage1Sim basefilename"
                    "\n  -dosage2Sim basefilename"
                    "\n       Dosage supports PLINK dosage file formats 1 and 2.  Use this"
                    "\n       option to define the appropriate format to read and the basename"
                    "\n       for PLINK's dosage .dat, .fam, and optional .map files for"
                    "\n       building genetic similarity"
                    "\n"
                    "\n  -afilesim basefilename"
                    "\n       basename for imputed genotype files .fam, .bim, and .big used"
                    "\n       to build genetic similarity."
                    "\n       These files are similar to PLINK binary files above."
                    "\n"
                    "\n  -sim filename"
                    "\n       file containing the genetic similarity matrix in tab-delimited"
                    "\n       ascii double values (overrides -fileSim)"
                    "\n"
                    "\n  -simOut filename"
                    "\n       write out genetic similarity matrix to filename in tab-delimited"
                    "\n       ascii double values"
                    "\n"
                    "\n  -NoZeroMeanSimSNPs"
                    "\n       Do not subtract the mean from similarity SNPs during pre-processing"
                    "\n       Default: false"
                    "\n"
                    "\n  -NoUnitVarSimSNPs"
                    "\n       Do not divide similarity SNPs by the standard deviation during"
                    "\n       pre-processing"
                    "\n       Default: false"
                    "\n"
                    "\n  -autoSelect filename"
                    "\n       Preprocess the SNP data to determine the best SNPs to include in"
                    "\n       the analysis. Results are written to two files using names derived"
                    "\n       from 'filename' in the following manner."
                    "\n       Using 'Experiment1.txt' as the 'filename' will produce"
                    "\n          'Experiment1.xval.txt'   contains the details of the cross"
                    "\n                                   validation selection process"
                    "\n          'Experiment1.snps.txt'   contains the list of SNPs selected which"
                    "\n                                   can be passed to the -extractsim option."
                    "\n"
                    "\n  -autoSelectFolds <int>"
                    "\n       During AutoSelect divide the individuals into <int> number of"
                    "\n       groups during the cross validation process."
                    "\n       Valid range is: 2 <= <int> <= count_of_individuals"
                    "\n       Default: 10"
                    "\n"
                    "\n  -autoSelectSearchValues <int_list>"
                    "\n       When running AutoSelect, specify which values to use to identify"
                    "\n       the region to search.  The <int_list> is either a filename or a"
                    "\n       list of integers separated by a whitespace or comma.  The list"
                    "\n       should be enclosed in double quotes."
                    "\n       Default: 0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,"
                    "\n                100,125,160,200,250,320,400,500,630,800,1000"
                    "\n"
                    "\n  -autoSelectCriterionLL"
                    "\n       Directs AutoSelect to use out-of-sample log likelihood for the"
                    "\n       selection criteria. [default]"
                    "\n"
                    "\n  -autoSelectCriterionMSE"
                    "\n       Directs AutoSelect to use out-of-sample mean-squared error for the"
                    "\n       selection criteria."
                    "\n"
                    "\n  -autoSelectApplyThinning"
                    "\n       Signal the autoSelect process to apply thinning to the search"
                    "\n       values used in AutoSelect."
                    "\n"
                    "\n  -autoSelectMaxSizeForThin <int>"
                    "\n       [NYI]Specify the largest autoSelectSearchValue where thinning will"
                    "\n       apply."
                    "\n       Default: 2000"
                    "\n"
                    "\n  -autoSelectRegionDistance <dbl>"
                    "\n       Filter SNPs from the same 'region' during the autoselection process."
                    "\n       The <dbl> defines a 'radius' in cMorgans around the SNP to 'thin'"
                    "\n       other SNPs during the autoSelect computation."
                    "\n       Default: 0.50"
                    "\n"
                    "\n  -covar filename"
                    "\n       optional file containing the covariates"
                    "\n"
                    "\n  -missingPhenotype <dbl>"
                    "\n       identifier for missing values.  If the phenotype for an"
                    "\n       individual is missing, then the individual is ignored."
                    "\n       If a covariate value for an individual is missing, then"
                    "\n       it is mean imputed."
                    "\n"
                    "\n  -out filename"
                    "\n       the name of the output file."
                    "\n       Default: [basefilename].out.txt"
                    "\n"
                    "\n  -verboseOutput"
                    "\n       Write more detailed information to the output file."
                    "\n       This includes things like logDelta, GeneticVar, etc."
                    "\n"
//                  "\n    -setrank [int k]                : [full]"
//                  "\n"
                    "\n  -ML"
                    "\n       use maximum likelihood parameter learning (default REML)"
                    "\n"
                    "\n  -REML"
                    "\n       use restricted Maximum likelihood parameter learning (default REML)"
                    "\n"
                    "\n  -simLearnType [FULL/ONCE]"
                    "\n       if set to Once, then delta, the ratio of residual"
                    "\n       to genetic covariance, is optimized only for the null model and"
                    "\n       used for each alternate model. If set to Full (the default), then the ratio"
                    "\n        is re-estimated for each alternative model."
                    "\n"
                    "\n  -simType [RRM/COVARIANCE]"
                    "\n       if set to RRM (the default), then the RRM is used for genetic"
                    "\n       similarity.  If set to COVARIANCE, then the empirical SNP"
                    "\n       covariance matrix is used."
                    "\n"
//hidden option     "\n    -UsePseudoCounts"
//hidden option     "\n       Use pseudocounts when estimating the SNP-standard deviation"
//hidden option     "\n       (as described in Eigenstrat)"
//hidden option     "\n       Default: false"
//                  "\n
//                  "\n    -ParamLearningType.[ML/REML] : ?? [ML] REML not yet implemented"
//                  "\n"
                    "\n  -useHeritability"
                    "\n       perform optimization of variance components in heritability h^2"
                    "\n       instead of delta  (h^2 = sigg2/(sigg2+sige2))"
                    "\n"
                    "\n  -brentStarts <int>"
                    "\n       number of interval boundary points for optimization of delta"
                    "\n       (see Section 2.1 of the Supplemental Information)."
                    "\n       Default: 10."
                    "\n"
                    "\n  -brentMaxIter <int>"
                    "\n       maximum  number of iterations per interval for the"
                    "\n       optimization of delta.  Default: 1e5."
                    "\n"
                    "\n  -brentMinLogVal <dbl>"
                    "\n       lower interval threshold for (log) delta optimization."
                    "\n       Default: -5."
                    "\n"
                    "\n  -brentMaxLogVal <dbl>"
                    "\n       upper interval treshold for (log) delta optimization."
                    "\n       Default: 5."
                    "\n"
                    "\n  -brentTol <dbl>"
                    "\n       convergence tolerance of Brent’s method used to optimize delta."
                    "\n       Default: 1e-6."
                    "\n"
//hidden option     //"\n  -UseBiasTerm"
//hidden option     //"\n       Use a bias term in the regression and output the weight. If"
//hidden option     //"\n       false the phenotype is zero-meaned instead."
//hidden option     //"\n       Default: false"
//hidden option     //"\n"
//hidden option     //"\n  -NoBiasTerm"
//hidden option     //"\n"
                    //"\n  -noHWE"
                    //"\n"
                    "\n  -HardyWeinberg"
                    "\n       use the Hardy-Weinberg estimate for the SNP variance instead"
                    "\n       of empirical estimate. This option is recommended for sample"
                    "\n       sizes that are too small for accurate variance estimation."
                    "\n"
                    "\n  -Beta <dbl_A> <dbl_B>"
                    "\n       Add more detail on Beta here"
                    "\n       Both dbl_A and dbl_B are optional and default to 1.0.  However"
                    "\n       if either dbl_A or dbl_B are specified, both must be entered."
                    "\n"
                    "\n  -runGwasType [RUN/NORUN]"
                    "\n       run GWAS or exit after computing the spectral decomposition of"
                    "\n       the genetic similarity matrix. Use NORUN, to cache the spectral"
                    "\n       decomposition.  This option, in combination with the next, is"
                    "\n       useful for parallelizing the tests of many SNPs. Default: RUN."
                    "\n"
                    "\n  -eigen directoryname"
                    "\n       load the spectral decomposition object from the directoryname."
                    "\n       When specified, the computations leading to the spectral"
                    "\n       decomposition of the genetic similarity matrix are skipped"
                    "\n       (note that that SNP file specifying the genetic similarities"
                    "\n       must still be given)."
                    "\n"
                    "\n  -eigenOut directoryname"
                    "\n       save the spectral decomposition object to the directoryname."
                    "\n       May be used with -runGwasType option."
                    "\n"
                    "\n  -logDelta <dbl>"
                    "\n       use a pre-specified log(delta) for both null and alternative models"
                    "\n"
                    "\n  -linReg"
                    "\n       use a linear regression model instead of a linear mixed model"
                    "\n"
                    "\n  -logReg"
                    "\n       use a logistic regression model instead of a linear mixed model"
                    "\n"
                    "\n  -groups [filename]"
                    "\n       specifies the filename"
                    "\n"
                    //"\n"
                    //"\n  -fileSimGroups [basefilename]"
                    //"\n       basefilename containing group-wise SNPs for building the RRM"
                    //"\n       in PLINK format"
                    //"\n"
                    //"\n  -tfileSimGroups [basefilename]"
                    //"\n       basefilename containing group-wise SNPs for building the RRM"
                    //"\n       in transposed PLINK format"
                    //"\n"
                    //"\n  -bfileSimGroups [basefilename]"
                    //"\n       basefilename containing group-wise SNPs for building the RRM"
                    //"\n       in binary PLINK format"
                    //"\n"
                    "\n  -simGroups [filename]"
                    "\n       file containing the group-by-group similarity matrix"
                    "\n"
                    //"\n  -UsePseudoCounts"
                    //"\n       use pseudocounts (i.e., Bayesian smoothing) when estimating"
                    //"\n       the SNP-standard deviation (as described in Eigenstrat)."
                    //"\n       Default: false."
                    //"\n"
                    "\n  -SnpPairs"
                    "\n         Compute GWAS using SNP pairs rather than individual SNPs."
                    "\n         SnpPairs sets the -ML option.  The default behaior computes"
                    "\n         all pairs of SNPs in the input."
                    "\n         NOTE:  Very compute intensive so use with care."
                    "\n"
                    "\n  -SnpId1 <snpId>"
                    "\n         This option modifies the behavior of SnpPairs to compute"
                    "\n         GWAS using SNP pairs but with Snp1 fixed as 'snpId1'."
                    "\n         -SnpId1 implies -SnpPairs but is not compatible with"
                    "\n         -BlockSNPs or -Tasks."
                    "\n"
                    "\n  -BlockSNPs <NumberOfBlocks> <Block1> <Block2>"
                    "\n         BlockSNPs is used in conjunction with -SnpPairs flag to partition"
                    "\n         the SNPs into NumberOfBlocks groups and then perform the SNP pair"
                    "\n         analysis by pairing each SNP in Block1 with each SNP in Block2"
                    "\n         Block1 and Block2 are 0 based and 0<=Block1<=Block2<NumberOfBlocks"
                    "\n         must hold.  To cover the entire space, FastLmmC should be invoked"
                    "\n         once with each unique block pair e.g."
                    "\n           for ( BlockI=0; BlockI<NumberOfBlocks; ++BlockI ) {"
                    "\n             for ( BlockJ=BlockI; BlockJ<NumberOfBlocks; ++BlockJ ) {"
                    "\n               fastlmmc -BlockSNPs NumberOfBlocks BlockI BlockJ"
                    "\n             }"
                    "\n           }"
                    "\n"
                    "\n  -Tasks <NumberOfTasks> <ThisTask>"
                    "\n         Tasks reshapes the BlockSNPs option to run more cleanly as a"
                    "\n         parametric sweep.  Enter the number of tasks you want break"
                    "\n         the job up into, the which job to run in this invocation and"
                    "\n         FastLmmC will figure out how to partition it.  There will likely"
                    "\n         be 'empty' files produced by the last few tasks at the number"
                    "\n         of tasks is defined by the forumlay (n*(n+1))/2 where n is the" 
                    "\n         number of partitions"
                    "\n"
                    "\n  Configuration Options:"
                    "\n"
                    "\n  -maxChromosomeValue <int>"
                    "\n       PLINK format specified human chromosome numbers and are thus"
                    "\n       limited to 0-26 (0=missing, 1-23, x=23, y=24, xy=25, mt=26)."
                    "\n       Use -MaxChromosomeValue to specify a different limit for the"
                    "\n       maximum chromosome number FastLmm accepts."
                    "\n"
                    "\n  -maxThreads <int>"
                    "\n       Suggests the number of threads for the math library to use"
                    "\n"
                    "\n  -maximizeWorkingSet"
                    "\n       Set program options to use most of available physical ram"
                    "\n       This will penalize other programs running on the system."
                    "\n"
                    "\n  -emitSnpIndexOnly"
                    "\n       Some large outputs can reduce file size by not writing the"
                    "\n       SNP id to the file, but only writing the snp index" 
                    "\n"
                    "\n  -noDosageRangeCheck "
                    "\n       Turn off range checking for dosage data and allow any real"
                    "\n       as the dosage value."
                    "\n"
                    "\n  -setOutputPrecision <int>"
                    "\n       Set the report out floating point precision to be <int>"
                    "\n       digits to the right of the decimal point.  3 <= <int> <= 18"
                    "\n       Default: 16"
                    "\n"
                    "\n  -pvaluePrintThreshold <dbl>"
                    "\n       Filter output to only report SNPs with p-values less than"
                    "\n       or equal to <dbl>"
                    "\n       Default: 1.0 (no filtering))"
                    "\n"
                    "\n  -randomSeed <int>"
                    "\n       Set the seed value for the random number generator."
                    "\n       Default: 1"
                    "\n"
                    "\n  Debug Options:"
                    "\n"
                    "\n  -verbose"
                    "\n       optionally print out log-likelihoods of null and alternative"
                    "\n       models to the standard output"
                    "\n"
                    "\n  -log"
                    "\n       Log tells the program to write more progress and debug data"
                    "\n       to a series of *.log files."
                    "\n"
                    "\n  -logDir [directory name]"
                    "\n       LogDir implies -log and sets the directory to write the output"
                    "\n       Default: .\\Logs"
                    "\n"
                    "\n  -Cluster"
                    "\n       Set when running on a cluster to reduce output and flush i/o"
                    "\n       more frequently to give visibility to headless cluster machines."
                    "\n"
//hidden option     "\n  -CreateValidationFiles"
//hidden option     "\n       create duplicate input files to validate that the data was"
//hidden option     "\n       correctly read"
//hidden option     "\n"
                    "\n-------------------------------------------------------------------"
                    "\nA sample command would be:"
                    "\n  fastlmmc -tfile geno_test -tfileSim geno_cov -pheno pheno.txt -covar covariate.txt"
                    "\n"
                    "\n  This requires the following input files:"
                    "\n"
                    "\n    geno_test.tfam, geno_test.tped"
                    "\n      These use PLINKs transposed file format for the test SNPs"
                    "\n"
                    "\n    geno_cov.tfam, geno_cov.tped"
                    "\n      Same format for the SNPs used to compute the RRM"
                    "\n"
                    "\n    pheno.txt"
                    "\n      PLINKs alternate phenotype file"
                    "\n"
                    "\n    covariate.txt"
                    "\n      The covariate file follows PLINKs covariate file format."
                    "\n      No header line is allowed."
                    "\n"
                    "\n  This produces geno_test.output.txt in the same location"
                    "\n  as geno_test files."
                    "\n" );
   exit(1);
   }

// small helper routines
bool IsLeadinChar( char c )
   {
   if ( (c == '-') 
#if defined( COMPILE_FOR_WINDOWS )
     || (c == '/')            // only detect on Windows.  Valid file leadin on *nix
#endif
     || (c == (char)0x96)
      )
      {
      return( true );
      }
   return( false );
   }

void CheckPlinkFileType( PlinkFileType& filetype, PlinkFileType type )
   {
   if ( filetype != UnknownFileType )
      {
      if ( type != filetype )
         {
         Fatal( "\nIncompatable options.  Cannot have multiple of -file, -bfile, -tfile, -dfile, and -afile options on the command line." );
         }
      }
   else
      {
      filetype = type;
      }
   }

void VerifyAlternativePhenotypeUsage()
   {
   // Verify all the proper alternative phenotype options are consistent
   if ( !alternatePhenotypeParameters.fUseAlternatePhenotype )
      {
      Fatal( "Alternate phenotype file required, but not specified." );
      }
   if ( alternatePhenotypeParameters.alternatePhenotypeFilename.empty() )
      {
      Fatal( "Specified alternate phenotype file but no alternate phenotype file name specified with -pheno option" );
      }
   if ( (alternatePhenotypeParameters.selectedPhenotypeIndex == 0) && alternatePhenotypeParameters.selectedPhenotypeId.empty() )
      {
      //Warn( "Alternate phenotype file was specified without an index or an index of 0."
      //      "\n          Defaulting to index 1"
      //      "\n          Use the -mpheno or -pheno-name options to change the default." );
      alternatePhenotypeParameters.selectedPhenotypeIndex = 1;
      }
   if ( (alternatePhenotypeParameters.selectedPhenotypeIndex != 0) && !alternatePhenotypeParameters.selectedPhenotypeId.empty() )
      {
      Fatal( "\nConflicting arguments specified alternate phenotype file.  Use only one of -mpheno or -pheno-name options" );
      }

   alternatePhenotypeParameters.selectedPhenotypeIndex -= 1;      // make zero based from now on!
   }

void VerifyAutoSelectUsage()
   {
   // Make sure all AutoSelect settings are correct and consistent w/
   //   other parameter settings
   if ( fAutoSelectSnps )
      {
      // -AutoSelect _MUST_ have the ?fileSim
      int cCovarianceFiles = (int)(plinkFilenameRootCovarianceVector.size());
      if ( cCovarianceFiles == 0 )
         {
         Fatal( "AutoSelect must have an input covariance file but we found none.  Use -?fileSim <filename>" );
         }
      else if ( cCovarianceFiles > 1 )
         {
         Fatal( "AutoSelect requires only one input covariance file but we found [%d]", cCovarianceFiles );
         }

      // -AutoSelect is not compatible with any of these options
      if ( !kernelFile.empty()                              // Cannot use a kernel file
//        || !snpDataCovFile.empty()                          //TODO:  OBSOLETE and NOT USED 
        || (computeEigenType == LOADFROMDISK)               // Cannot load similarity matrix from disk
        || identityKernel                                   // Cannot be using an identityKernel
        || linReg                                           // Cannot be running Linear Regression
        || fUseLogisticRegression                           // Cannot be running Logistic Regression
        || !simGroups.empty()                               // 
        || !groupsFile.empty()
         )
         {
         Fatal( "AutoSelect is not compatible with other covariance options" );
         }

      /*
       *  If there is no data file defined nor an output file defined
       *    do not run GWAS.
       */
      if ( !plinkFilenameRoot.empty() )
         {
         Warn( "-AutoSelect option selected.  Will not use test file [%s]", plinkFilenameRoot.c_str() );
         plinkFilenameRoot.clear();
         }
      if ( !outputFilename.empty() )
         {
         Warn( "-AutoSelect option selected.  Will not produce GWAS output file [%s]", outputFilename.c_str() );
         outputFilename.clear();
         }
      if ( outputFilename.empty() && plinkFilenameRoot.empty() )
         {
         runGwasType = NORUN;
         }

      /*
       * If the user did not define the AutoSelectSearchValues,
       *   initialize the list with default values
       *   We use a sentinal of -1 in the defaultAutoSelectSearchValues
       *   array to signal end of initialization data.
       */
      if ( autoSelectSearchValues.empty() )
         {
         for ( int i=0; defaultAutoSelectSearchValues[ i ] != -1; ++i )
            {
            autoSelectSearchValues.push_back( defaultAutoSelectSearchValues[ i ] );
            }
         }
      }
   }

void VerifySnpFilterUsage( int cJobsToRun_, int iJobToRun_ )
   {
   // all SnpFilter information _should_ be set correctly or validated during 
   //   snpFilterOpts.Setxxx _Except_ if we only got partial job information
   if ( !snpFilterOpts.FUseSnpFilter() && !((iJobToRun_ == -1) && (cJobsToRun_ == -1)) )
      {
      // rats.  WE say no filtering, but we got one but not the other Jobs option
      if ( iJobToRun_ == -1 )
         {
         Fatal( "-numjobs option found, but not -thisjob. Usage of -numjobs option requires -thisjob to be specified too." );
         }
      else if ( cJobsToRun_ == -1 )
         {
         Fatal( "-thisjob option found, but not -numjobs. Usage of -thisjob option requires -numjobs to be specified too." );
         }
      else
         {
         Fatal( "Program malfunction processing job info." );
         }
      }

   if ( snpFilterOpts.FFilterByFile() )      // Add filtering to Data file?
      {
      snpFilterOpts.LoadSnpIdsToExtract();
      }

   if ( simFilterOpts.FFilterByFile() )      // Add filtering to Similarity file?
      {
      simFilterOpts.LoadSnpIdsToExtract();
      }
   }

void VerifyLogRegUsage()
   {
   if ( fUseLogisticRegression )
      {
      // Cannot do both Linear and Logistic Regression
      if ( linReg )
         {
         Fatal( "Incompatible options: -linReg is not compatible with -logReg option" );
         }

      // Logistic regression does not use a covariance file
      if ( !plinkFilenameRootCovarianceVector.empty() )
         {
         Fatal( "Incompatible options: -logReg is not compatible with -xFileSim options" );
         }
      }
   }

void VerifyLinRegUsage()
   {
   if ( linReg )
      {
      // Cannot do both Linear and Logistic Regression
      if ( fUseLogisticRegression )
         {
         Fatal( "Incompatible options: -linReg is not compatible with -logReg option" );
         }

      // Linear regression does not use a covariance file
      if ( !plinkFilenameRootCovarianceVector.empty() )
         {
         Fatal( "Incompatible options: -linReg is not compatible with -xFileSim options" );
         }
      }
   }

void VerifySnpPairsUsage()
   {
   if ( fDoSnpPairs )
      {
      if ( mlType != ML )
         {
         Fatal( "Incompatible options: -SnpPair currently requires ML option" );
         }
      if ( covLearnType != ONCE )
         {
         Fatal( "Incompatible options: -SnpPair currently requires covLearnType set to ONCE" );
         }
      if ( fDoSnpId1 )  // 
         {
         if ( snpFilterOpts.FFilterByTask() || snpFilterOpts.FFilterByBlocks() )
            {
            Fatal( "Incompatible options: -SnpId1 is not compatible with -Task or -BlockSNPs" );
            }
         }
      }
   }

void VerifyBrentUsage()
   {
   if ( minLogValGlobal >= maxLogValGlobal )
      {
      Fatal( "Incompatible options: -brentMinLogVal[%.3e] must be less than -brentMaxLogVal[%.3e]", minLogValGlobal, maxLogValGlobal );
      }
   }

/*
 * Ensure we have all the 'right stuff' for running
 *  the iterative FastLmm path
 */
void VerifyOptimizeForLowMemoryUsage()
   {
   if ( fDoIteratingFastLmm )
      {
      //  CANNOT have a *sim file
      if ( plinkFileTypeCovariance != UnknownFileType )
         {
         Fatal( "Incompatible options: -LowMem cannot be used with any other xfilesim option." );
         }
      // CANNOT have any SNP filter options set
      if ( snpFilterOpts.filterType != SnpFilterOptions::FilterByNone )
         {
         Fatal( "Incompatible options: -LowMem cannot be used with SNP filtering." );
         }
      if ( simFilterOpts.filterType != SnpFilterOptions::FilterByNone )
         {
         Fatal( "Incompatible options: -LowMem cannot be used with sim filtering" );
         }
      // CANNOT have -autoSelect
      if ( fAutoSelectSnps )
         {
         Fatal( "Incompatible options: -LowMem cannot be used with AutoSelect" );
         }
      if ( plinkFileType != Binary )
         {
         Fatal( "Incompatible options: -LowMem requires a -bfile" );
         }
      // MUST have fullrank
      if ( !forceFullRank )
         {
         Warn( "-LowMem requires -fullrank.  Setting -fullrank" );
         forceFullRank = true;
         }
      }
   }

void VerifyNormalizationUsage()
   {
   if ( simMissingness && (normalizationInfo.normalizationType == HardyWeinberg) )
      {
      Warn( "simMissingness is not compatible with HardyWeinberg option.  Removing HardyWeinberg." );
      normalizationInfo.normalizationType = Normal;
      }
   }

// return true if pch points to a valid integer list string
bool IsIntListString( const char* pch )
   {
   const char *pchT = pch;
   while( *pchT && ( isspace( *pchT ) || isdigit( *pchT ) || (*pchT == ',') ) )
      {
      ++pchT;
      }
   return( (*pchT == '\0') && (pchT != pch) );  // true if at end-of-string and not zero length
   }

// return true if str is a valid integer list string
bool IsIntListString( const string& str )
   {
   return( IsIntListString( str.c_str() ) );
   }

// populate the vector with the integers found in the integer list 
void PopulateIntListFromString( vector< int >& vec, string& str )
   {
   char *pchDst;
   char rgch[ 32 ];
   const char *pchSrc = str.c_str();

   while( *pchSrc )
      {
      while( isspace( *pchSrc ) || (*pchSrc == ',') ) { ++pchSrc; }     // skip whitespace or ','
      pchDst = rgch;
      while( *pchSrc && isdigit( *pchSrc ) ) { *pchDst++ = *pchSrc++; } // accumulate digits
      if ( pchDst != rgch )
         {
         *pchDst = '\0';
         int iT;
         sscanf( rgch, "%i", &iT );
         vec.push_back( iT );
         }
      else if ( *pchSrc )
         {
         Fatal( "Unrecognized input in integer list [%s]", pchSrc );
         }
      }
   }

void PopulateIntListFromFile( vector< int >& vec, const string& str )
   {
   FILE *pfile = fopen( str.c_str(), "r" );
   if ( pfile == nullptr )
      {
      Fatal( "Expected integer list in file [%s]"
           "\n        CRT error: %d  [%s]", str.c_str(), errno, strerror( errno ) );
      }

   string integerList;
   const size_t bytesToRead = 1024 * 32;
   char rgch[ bytesToRead + 1 ];

   while( !feof( pfile ) )
      {
      size_t cchRead = fread( rgch, sizeof( char ), bytesToRead, pfile );
      if ( ferror( pfile ) )
         {
         Fatal( "Error reading file [%s]"
              "\n        CRT error: %d  [%s]", str.c_str(), errno, strerror( errno ) );
         }
      rgch[ cchRead ] = '\0';
      integerList += rgch;
      }

   fclose( pfile );
   PopulateIntListFromString( vec, integerList );   
   }

void ProcessCommandLine( int argc, char* argv[] )
   {
   char *pch;
   int iArg;

   if ( argc == 1 )
      {
      Usage();
      }

   // set our 'default' values
   plinkFileType = UnknownFileType;
   plinkFileTypeCovariance = UnknownFileType;   

   excludeSnpParams.distType = NoDist;
   excludeSnpParams.fileType = NoFile;
   excludeSnpParams.geneticDist = std::numeric_limits<double>::quiet_NaN();
   excludeSnpParams.bpDist = -1;

// obsolete?  remove in future since we only support INDIVIDUAL_ONLY
   individualIdType = INDIVIDUALID_ONLY;

   missingPhenotypeString = "-9";
   missingPhenotypeValue = (real)atof( missingPhenotypeString.c_str() );
   covLearnType = FULL;

   // mlType is initialized to REML, then we must set use_f_test, useBiasTerm, and clear useWoodbury
   SetMlType( REML );

   missingType = MEAN;              //by default impute by mean value
   numberOptimizationStartsGlobal = 10;
   brentMaxIterGlobal = 100000;
   minLogValGlobal = -5.0;
   maxLogValGlobal =  5.0;
   brentTolGlobal = (real)1e-6;
   computeEigenType = NOSAVE;
   runGwasType = RUN; 
   usePseudocounts = false;
   kernelType = RRM;
   indPheno = 0; //index of the phenotype in phenDataFile to run
   fDoSnpPairs = false;

   // REVIEW:  set debug options to default values 
   //          set verboseLevel during development
   verboseLevel = 0;
   fCreateReadValidationFiles = false;

   // a couple of locals ...
   int iJobToRun = -1;
   int cJobsToRun = -1;

   for ( iArg=1; iArg<argc; ++iArg )
      {
      pch = argv[ iArg ];
      if ( !__isascii(*pch) )
         {
         Error( "Encountered non-ASCII character '%c' [0x%02X] in command line argument %d [%s]", *pch, (unsigned char)*pch, iArg, pch );
         }

      if ( !IsLeadinChar( *pch ) )
         {
         break;
         }

      ++pch;
      if ( (_strcmpi( pch, "?" ) == 0 ) || (_strcmpi( pch, "help" ) == 0) || (_strcmpi( pch, "h" ) == 0) )
         {
         Usage();
         }
      else if ( (_strcmpi( pch, "Verbose" ) == 0 )
             || (_strcmpi( pch, "v" ) == 0 )
              )
         {
         if ( (iArg+1 >= argc)
           || (!isdigit( *argv[iArg+1] )) 
            )
            {
            // no following arguments.  Set verboseLevel to 2 to show Minor Progress
            verboseLevel = verboseMinor;
            }
         else
            {
            ++iArg;
            sscanf( argv[ iArg ], "%i", &verboseLevel );
            }
         }
      else if ( (_strcmpi( pch, "useheritability" ) == 0) )
         {
         useHeritability = true;
         minLogValGlobal = max(minLogValGlobal, (1E-5) );
         maxLogValGlobal = min(maxLogValGlobal, 1.0-(1E-5) );
         }
      else if ( (_strcmpi( pch, "skipalternative" ) == 0) )
         {
         runAlternativeModels = false;
         }
      else if ( (_strcmpi( pch, "SnpPairs" ) == 0) )
         {
         SetDoSnpPairs();
         }
      else if ( (_strcmpi( pch, "SnpId1" ) == 0) )
         {
         if ( iArg+1 >= argc )
            {
            Fatal( "-SnpId1 needs to be followed by a SnpId." );
            }

         fDoSnpId1 = true;
         snpId1Global = argv[ ++iArg ];
         SetDoSnpPairs();
         }
      else if ( (_strcmpi( pch, "BlockSNPs" ) == 0)
             || (_strcmpi( pch, "BS" ) == 0)
              )
         {
         if ( iArg+3 >= argc )
            {
            Fatal( "-BlockSNPs needs to be followed by 3 integers specifying number of blocks and the block pair to analyze." );
            }
         char *pchBlockCount = argv[ iArg+1 ];
         char *pchBlock1 = argv[ iArg+2 ];
         char *pchBlock2 = argv[ iArg+3 ];
         int blockCountT;
         if ( IsLeadinChar( *pchBlockCount ) 
           || (sscanf( pchBlockCount, "%i", &blockCountT ) != 1)
           || (blockCountT <= 0)
            )
            {
            Fatal( "Expeced a positive BlockCount to follow -BlockSNPs.  Found [%s]", pchBlockCount );
            }
//         blockRunParams.cBlocksMax = t;

         int iBlockT1;
         if ( IsLeadinChar( *pchBlock1 ) 
           || (sscanf( pchBlock1, "%i", &iBlockT1 ) != 1)
           || (iBlockT1 < 0)
           || (iBlockT1 >= blockCountT)
            )
            {
            Fatal( "-BlockSNPs needs 0 <= blockId 1[%d] < Block Count[%d].  Found [%s]", blockCountT, pchBlock1 );
            }
//         blockRunParams.iBlock1 = t;

         int iBlockT2;
         if ( IsLeadinChar( *pchBlock2 ) 
           || (sscanf( pchBlock2, "%i", &iBlockT2 ) != 1) 
           || (iBlockT2 < iBlockT1)
           || (iBlockT2 >= blockCountT)
            )
            {
            Fatal( "-BlockSNPs needs blockId 1[%d] <= blockId 2[%d] < Block Count[%d].  Found [%s]", iBlockT1, iBlockT2, blockCountT, pchBlock2 );
            }

         snpFilterOpts.SetByBlocks( blockCountT, iBlockT1, iBlockT2 );
         SetDoSnpPairs();
         iArg += 3;
         }
      else if ( (_strcmpi( pch, "Tasks" ) == 0)
             //|| (_strcmpi( pch, "Ta" ) == 0)
              )
         {
         if ( iArg+2 >= argc )
            {
            Fatal( "-Tasks needs to be followed by 2 integers specifying number of tasks and the specific task to run." );
            }
         char *pchTasks = argv[ iArg+1 ];
         char *pchIdxTask = argv[ iArg+2 ];
         int cTasksT;
         if ( IsLeadinChar( *pchTasks ) 
           || (sscanf( pchTasks, "%i", &cTasksT ) != 1)
           || (cTasksT <= 0)
            )
            {
            Fatal( "Expected a positive number of Tasks to follow -Tasks.  Found [%s]", pchTasks );
            }

         int iTaskT;
         if ( IsLeadinChar( *pchIdxTask ) 
           || (sscanf( pchIdxTask, "%i", &iTaskT ) != 1)
           || (iTaskT < 0)
           || (iTaskT >= cTasksT)
            )
            {
            Fatal( "-Tasks needs 0 <= iTask[%d] < Tasks[%d].  Found [%s]", iTaskT, cTasksT );
            }

         snpFilterOpts.SetByTask( cTasksT, iTaskT );
         SetDoSnpPairs();
         iArg += 2;
         }
      else if ( (_strcmpi( pch, "Timers" ) == 0)
             || (_strcmpi( pch, "t" ) == 0)
              )
         {
         if ( (iArg+1 >= argc)
           || (!isdigit( *argv[iArg+1] )) 
            )
            {
            // no following arguments.  Set fReport times to show everything
            fReportTimes = -1;
            }
         else
            {
            ++iArg;
            sscanf( argv[ iArg ], "%i", &fReportTimes );
            }
         }
      else if ( (_strcmpi( pch, "Log" ) == 0 )
//             || (_strcmpi( pch, "l" ) == 0 )
              )
         {
         fWriteLogFile = true;
         }
      else if ( (_strcmpi( pch, "LogDir" ) == 0 )
             || (_strcmpi( pch, "ld" ) == 0 )
         )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "-LogDir needs to be followed by a path specifying the location to write log files." );
            }
         fWriteLogFile = true;
         logDir = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "MaxThreads" ) == 0)
             || (_strcmpi( pch, "mt" ) == 0)
         )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "-MaxThreads needs to be followed by a number specifying the number of threads." );
            }
         int cThreads;
         sscanf( argv[ iArg ], "%i", &cThreads );
         if ( cThreads <= 0 )
            {
            Fatal( "-MaxThreads argument needs to be positive.  Found %d.", cThreads );
            }
#if defined( USE_MKL )
         //mkl_domain_set_num_threads( cThreads, MKL_DOMAIN_ALL )
         mkl_set_num_threads( cThreads );
#endif
         }
      else if ( (_strcmpi( pch, "SetOutputPrecision" ) == 0)
             || (_strcmpi( pch, "sop" ) == 0)
              )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "-SetOutputPrecision needs to be followed by a number specifying the output precision." );
            }
         int outputPrecisionT;
         sscanf( argv[ iArg ], "%i", &outputPrecisionT );
         if ( outputPrecisionT < 3 )
            {
            Warn( "-SetOutputPrecision must >= 3.  Found %d.  Using 3", outputPrecisionT );
            outputPrecisionT = 3;
            }
         else if ( outputPrecisionT > 18 )
            {
            Warn( "-SetOutputPrecision must <= 18.  Found %d.  Using 18", outputPrecisionT );
            outputPrecisionT = 18;
            }
         fSetFastLmmRealsPrintPrecision = true;
         fastLmmRealsPrintPrecision = outputPrecisionT;
         }
      else if ( (_strcmpi( pch, "NoDosageRangeCheck" ) == 0)
             || (_strcmpi( pch, "ndrc" ) == 0)
              )
         {
         fNoDosageRangeCheck = true;
         }
      else if ( (_strcmpi( pch, "PvaluePrintThreshold" ) == 0)
             || (_strcmpi( pch, "ppt" ) == 0)
              )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "-PvaluePrintThreshold needs to be followed by a number specifying pvalue limit." );
            }

         // use float input since a double can successfully hold any float the user provides
         float pvaluePrintThresholdT = std::numeric_limits<float>::quiet_NaN();
         int rc = sscanf( argv[ iArg ], "%f", &pvaluePrintThresholdT );
         if ( (rc == EOF) || (rc == 0) )
            {
            Fatal( "-PvaluePrintThreshold needs to be followed by a number specifying pvalue limit.  Found [%s]", argv[ iArg ] );
            }
         if ( _isnan( pvaluePrintThreshold ) )
            {
            Fatal( "-PvaluePrintThreshold could not represent your input [%s].  Please simplify.", argv[ iArg ] );
            }
         if ( (pvaluePrintThresholdT <= 0.0) || (pvaluePrintThresholdT > 1.0)  )
            {
            Fatal( "-PvaluePrintThreshold must be in the range 0.0 < PvalueThreadhold <= 1.0.  Found [%e]", pvaluePrintThresholdT );
            }
         pvaluePrintThreshold = pvaluePrintThresholdT;
         }
      else if ( (_strcmpi( pch, "Sort" ) == 0) )
         {
         fEmitSorted = true;
         }
      else if ( (_strcmpi( pch, "Sort-" ) == 0)
             || (_strcmpi( pch, "NoSort" ) == 0)
              )
         {
         fEmitSorted = false;
         }
      else if ( (_strcmpi( pch, "ReproducibleFP" ) == 0)
             || (_strcmpi( pch, "rfp" ) == 0)
              )
         {
         fUseReproducibleFP = true;
         SetReproducibleFP( fUseReproducibleFP );
         }
      else if ( (_strcmpi( pch, "ReproducibleFP-" ) == 0)
             || (_strcmpi( pch, "rfp-" ) == 0)
              )
         {
         fUseReproducibleFP = false;
         SetReproducibleFP( fUseReproducibleFP );
         }
      else if ( (_strcmpi( pch, "RandomSeed" ) == 0)
             || (_strcmpi( pch, "rs" ) == 0)
              )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "-RandomSeed needs to be followed by a number specifying seed value for the random number generator." );
            }

         unsigned long seedT;
         char *pchT = argv[ iArg ];
         if ( FIsDigitString( pchT ) )
            {
            sscanf( pchT, "%lu", &seedT );
            // No range validation is really needed
            }
         else
            {
            seedT = (unsigned long)time( NULL );      // possible overflow as time() return 64bit value now
            }
         randomNumberSeed = seedT;
         fRngSeedSuppliedByUser = true;
         }
      else if ( (_strcmpi( pch, "createvalidationfiles" ) == 0)
             || (_strcmpi( pch, "cvf" ) == 0)
              )
         {
         fCreateReadValidationFiles = true;
         }
      else if ( (_strcmpi( pch, "numjobs" ) == 0)
             || (_strcmpi( pch, "TaskCount" ) == 0)
              )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "%s needs to be followed by an integer number of partions to work on.", pch );
            }
         int i;
         sscanf( argv[ iArg ], "%i", &i );
         if ( i < 1 )
            {
            Error( "%s parameter must be positive.  Found: %d", pch, i );
            }
         else
            {
            cJobsToRun = i;
            if ( iJobToRun != -1 )
               {
               snpFilterOpts.SetByJob( cJobsToRun, iJobToRun );
               }
            }
         }
      else if ( (_strcmpi( pch, "thisjob" ) == 0 )
//             || (_strcmpi( pch, "Tasks" ) == 0)
              )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "%s needs to be followed by an integer number indicating which partion to work on.", pch );
            }
         int i;
         sscanf( argv[ iArg ], "%i", &i );
         if ( i < 0 )
            {
            Error( "%s parameter must be non-negative.  Found: %d", pch, i );
            }
         else
            {
            iJobToRun = i;
            if ( cJobsToRun != -1 )
               {
               snpFilterOpts.SetByJob( cJobsToRun, iJobToRun );
               }
            }
         }
      else if ( (_strcmpi( pch, "pedOutFile" ) == 0)
             || (_strcmpi( pch, "POF" ) == 0)
              )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "-pedOutFile needs to be followed by a filename specifying ped filename." );
            }
         fWritePedFile = true;
         pedOutputFile = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "extract" ) == 0) )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "-extract needs to be followed by a filename specifying the SNPs to work on." );
            }
         if ( snpFilterOpts.FUseSnpFilter() )
            {
            Fatal( "Only one SNP filter can be applied.  [-extract %s ] is an error.", argv[iArg] );
            }
         string s = argv[ iArg ];
         snpFilterOpts.SetByFileList( s );
         }
      else if ( (_strcmpi( pch, "extractSim" ) == 0) )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "-extractSim needs to be followed by a filename specifying the SNPs to use for the genetic similarity matrix." );
            }
         string s = argv[ iArg ];
         simFilterOpts.SetByFileList( s );
         }
      else if ( (_strcmpi( pch, "extractSimTopK" ) == 0)
             || (_strcmpi( pch, "esk" ) == 0)
              )
         {
         if ( iArg+2 >= argc )
            {
            Fatal( "-extractSimTopK needs two arguments, a filename and an integer." );
            }
         string fname = argv[ ++iArg ];
         int i;
         sscanf( argv[ ++iArg ], "%i", &i );
         if ( i <= 0 )
            {
            Fatal( "-extractSimTokK needs a positive integer argument after the filename.  Found [%d]", i );
            }
         simFilterOpts.SetByTopN( fname , i );
         }
      else if ( (_strcmpi( pch, "fullrank" ) == 0) )
         {
         forceFullRank = true;
         }
      else if ( (_strcmpi( pch, "noUnitVarTestSnps" ) == 0)
             || (_strcmpi( pch, "noUnitVarianceTestSnps" ) == 0)
             || (_strcmpi( pch, "nuvt" ) == 0)
              )
         {
         unitVarianceTestSNPs = false;
         }
      else if ( (_strcmpi( pch, "noZeroMeanTestSnps" ) == 0)
             || (_strcmpi( pch, "nzmt" ) == 0)
              )
         {
         zeroMeanTestSNPs = false;
         }
      else if ( (_strcmpi( pch, "noUnitVarSimSnps" ) == 0)
             || (_strcmpi( pch, "noUnitVarianceSimilaritySnps" ) == 0)
             || (_strcmpi( pch, "nuvs" ) == 0)
              )
         {
         unitVarianceSimilaritySNPs = false;
         }
      else if ( (_strcmpi( pch, "noZeroMeanSimSnps" ) == 0)
             || (_strcmpi( pch, "noUnitVarianceSimilaritySnps" ) == 0)
             || (_strcmpi( pch, "nzms" ) == 0)
              )
         {
         zeroMeanSimilaritySNPs = false;
         }
      else if ( (_strcmpi( pch, "identity" ) == 0) )
         {
         identityKernel = true;
         logDeltaPreset = (real)0.0;
         covLearnType = ONCE;
         }
      else if ( (_strcmpi( pch, "linreg" ) == 0) )
         {
         linReg = true;
         }
      else if ( (_strcmpi( pch, "logreg" ) == 0) )
         {
         fUseLogisticRegression = true;
         }
      else if ( (_strcmpi( pch, "excludeById" ) == 0) )
        {
        Fatal( "excludeByID is not implemented yet." );
        excludeSnpParams.distType = ByID;
        //if no filetype has been selected, default to SimFile
        if ( excludeSnpParams.fileType == NoFile )
           {
           excludeSnpParams.fileType = FileSim;
           }
        }
      else if ( (_strcmpi( pch, "excludeByPosition" ) == 0)
             || (_strcmpi( pch, "ebp" ) == 0)
              )
         {
         //forceFullRank = true;
         useWoodbury = false;
         excludeSnpParams.distType = ByPosition;
         //excludeSnpParams.geneticDist = std::numeric_limits<double>::quiet_NaN();
         //excludeSnpParams.bpDist = -1;

         //if no filetype has been selected, default to SimFile
         if ( excludeSnpParams.fileType == NoFile )
            {
            excludeSnpParams.fileType = FileSim;
            }
         // get the distance
         if ( (++iArg >= argc) )
            {
            Fatal( "-ExcludeByPosition must be followed by an integer argument" );
            }
         
         int distance;
         sscanf( argv[ iArg ], "%i", &distance );     // use sscanf() to support other integer input formats (e.g. 0x400)
         if ( distance < 0 )
            {
            Fatal( "-ExcludeByPosition argument must be greater or equal to 0.  Found %d", distance );
            }
         excludeSnpParams.bpDist = distance;
         }
      else if ( (_strcmpi( pch, "excludeByGeneticDistance" ) == 0)
             || (_strcmpi( pch, "ebgd" ) == 0)
              )
         {
         //forceFullRank = true;
         useWoodbury = false;
         excludeSnpParams.distType = ByGeneticDistance;
         //if no filetype has been selected, default to SimFile
         if ( excludeSnpParams.fileType == NoFile )
            {
            excludeSnpParams.fileType = FileSim;
            }
         // get the distance
         if ( (++iArg >= argc) )
            {
            Fatal( "-ExcludeByGeneticDistance must be followed by a numeric argument" );
            }
         if ( !FIsFpString( argv[ iArg ] ) )
            {
            Fatal( "ExcludeByGeneticDistance must be followed by a numeric argument.  Found [%s]", argv[ iArg ] );
            }
         real distance = atof( argv[ iArg ] );
         if ( distance < 0.0 )
            {
            Fatal( "-ExcludeByGeneticDistance argument must be greater or equal to 0.0.  Found %f", distance );
            }
         excludeSnpParams.geneticDist = distance;
         }
      else if ( (_strcmpi( pch, "trainOnTest" ) == 0) )
         {
         Fatal( "trainOnTest is not implemented yet." );
         //note: the normalization of SNPs has to be done to match the kernel SNPs =>divide by #SNPs in Kernel
         excludeSnpParams.fileType = TrainOnTest;
         //if no distance type has been set, default to the SNP-ID
         if ( excludeSnpParams.distType == NoDist )
            {
            excludeSnpParams.distType = ByID;
            }
         }
      else if ( (_strcmpi( pch, "HardyWeinberg" ) == 0)
             || (_strcmpi( pch, "hwe" ) == 0)
              )
         {
         normalizationInfo.normalizationType = HardyWeinberg;
         }
      else if ( (_strcmpi( pch, "Beta" ) == 0) )
         {
         // Set the normalization to Beta
         normalizationInfo.normalizationType = Beta;
         if ( (iArg + 2 < argc) )
            {
            // Have the Beta argument and enough options following to cover both A and B
            if ( FIsFpString( argv[ iArg + 1 ] ) && FIsFpString( argv[ iArg + 2 ] ) )
               {
               real rT = atof( argv[ iArg + 1 ] );
               if ( rT <= 0.0 )
                  {
                  Fatal( "Beta argument 'A' must be a non-zero, positive number.  Found [%.2f]", rT );
                  }
               normalizationInfo.betaA = rT;
               rT = atof( argv[ iArg + 2 ] );
               if ( rT <= 0.0 )
                  {
                  Fatal( "Beta argument 'B' must be a non-zero, positive number.  Found [%.2f]", rT );
                  }
               normalizationInfo.betaB = rT;
               iArg += 2;
               }
            }
         }
      else if ( (_strcmpi( pch, "simMissingness" ) == 0) )
         {
         simMissingness = true;
         }
      else if ( (_strcmpi( pch, "logDelta" ) == 0) )
         {
         // get the delta
         if ( (++iArg >= argc) )
            {
            Fatal( "logDelta needs to be followed by a double argument." );
            }
         if ( !FIsFpString( argv[ iArg ] ) )
            {
            Fatal( "logDelta must be followed by a numeric argument.  Found [%s]", argv[ iArg ] );
            }
         logDeltaPreset = (real)atof( argv[ iArg ] ); //will be used as delta_0
         covLearnType = ONCE;                         //all alternative models use the delta_0
         }
      else if ( (_strcmpi( pch, "setRank" ) == 0) )
         {
         forceFullRank = true;
         // get the rank
         if ( (++iArg >= argc) )
            {
            Fatal( "setRank needs to be followed by an integer argument." );
            }
         sscanf( argv[ iArg ], "%i", &setRank );
         } 
      else if ( (_strcmpi( pch, "ftest" ) == 0) )
         {
         use_f_test = true;
         }
      else if ( (_strcmpi( pch, "ftest-" ) == 0) )
         {
         use_f_test = false;
         }
      else if ( (_strcmpi( pch, "permute.kernel" ) == 0) )
         {
         if(permutationMode!=NO_PERM) { Fatal( "only one permuation can be specified" ); }
         permutationMode = PERMUTE_KERNEL;
         }
      else if ( (_strcmpi( pch, "nonpsd.subtractmin" ) == 0) )
         {
         handle_non_psd = SUBTRACT_MIN;
         }
      else if ( (_strcmpi( pch, "nonpsd.setzero" ) == 0) )
         {
         handle_non_psd = SET_TO_ZERO;
         }
      else if ( (_strcmpi( pch, "permute.phen" ) == 0) )
         {
         if(permutationMode!=NO_PERM) { Fatal( "only one permuation can be specified" ); }
         permutationMode = PERMUTE_PHEN;
         }
      else if ( (_strcmpi( pch, "woodbury" ) == 0) )
         {
         useWoodbury = true;
         }
      else if ( (_strcmpi( pch, "noWoodbury" ) == 0) )
         {
         useWoodbury = false;
         }
      else if ( (_strcmpi( pch, "usePseudoCounts" ) == 0) )
         {
         usePseudocounts = true;
         }
      else if ( (_strcmpi( pch, "covar" ) == 0) )
         {
         // get the covariatesFile name
         if ( (++iArg >= argc) )
            {
            Fatal( "covar must be followed by a filename." );
            }
         covariatesFile = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "out" ) == 0) )
         {
         // get the outputFilename
         if ( (++iArg >= argc) )
            {
            Fatal("out needs to be followed by a filename.");
            }
         outputFilename = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "VerboseOutput" ) == 0)
             || (_strcmpi( pch, "verboseOut" ) == 0)        // obsolete option name
             || (_strcmpi( pch, "vo" ) == 0)
         )
         {
         fVerboseOutput = true;
         }
      else if ( (_strcmpi( pch, "sim" ) == 0) )
         {
         // get the kernelFile name
         if ( (++iArg >= argc) )
            {
            Fatal( "sim needs to be followed by a filename" );
            }
         kernelFile = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "simOut" ) == 0) )
         {
         // get the kernelFile name
         if ( (++iArg >= argc) )
            {
            Fatal( "simOut needs to be followed by a filename" );
            }
         kernelFileOut = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "MemoryFraction" ) == 0) 
             || (_strcmpi( pch, "MF" ) == 0)
             )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "MemoryFraction needs to be followed by a floating point number from 0.1 to 0.9" );
            }

         real rT = atof( argv[ iArg ] );
         if ( rT < 0.10 )
            {
            Warn( "MemoryFraction option is too small.  Setting to 0.1" );
            rT = (real)0.10;
            }
         if ( rT > 0.90 )
            {
            Warn( "MemoryFraction option is too large.  Setting to 0.9" );
            rT = (real)0.90;
            }
         mFraction = rT;
         }
      else if ( (_strcmpi( pch, "TopKbyLinReg" ) == 0) 
             || (_strcmpi( pch, "TK" ) == 0)
             )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "TopKbyLinReg needs to be followed by a number." );
            }
         int iT;
         sscanf( argv[ iArg ], "%i", &iT );        // use sscanf() to support other integer input formats (e.g. 0x400)
         if ( iT < 1 )
            {
            // TopKbyLinReg must be a positive number upto the number of SNPs
            Fatal( "-TopKbyLinReg must be a positive number.  Found [%d]", iT );
            }
         topk = iT;
         selectTopK = true;
         }
      else if ( (_strcmpi( pch, "AutoSelect" ) == 0) 
             || (_strcmpi( pch, "AS" ) == 0)
             )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "AutoSelect needs to be followed by a filename." );
            }
         autoSelectOutputFilenameRoot = argv[ iArg ];
         fAutoSelectSnps = true;
         }
      else if ( (_strcmpi( pch, "AutoSelectSearchValues" ) == 0)
             || (_strcmpi( pch, "ASSV" ) == 0)
         )
         {
         // Get the search values to use or the filename
         //   command line search values are a start with 
         //   a double quote " and are a list of whitespace
         //   separated list of non-negative integers
         if ( ++iArg >= argc )
            {
            Fatal( "AutoSelectSearchValues needs to be followed by a file name or \"list of integers\"." );
            }

         string strT = argv[ iArg ];
         if ( IsIntListString( strT ) )
            {
            PopulateIntListFromString( autoSelectSearchValues, strT );
            }
         else
            {
            PopulateIntListFromFile( autoSelectSearchValues, strT );
            }
         if ( autoSelectOutputFilenameRoot.empty() )
            {
            autoSelectOutputFilenameRoot = "fastlmm.autoselect.results.txt";
            fAutoSelectSnps = true;                      // force AutoSelect too.
            }
         }
      else if ( (_strcmpi( pch, "AutoSelectCriterionLL" ) == 0)
             || (_strcmpi( pch, "ASCL" ) == 0)
              )
         {
         autoSelectCriterion = useLL;
         }
      else if ( (_strcmpi( pch, "AutoSelectCriterionMSE" ) == 0)
             || (_strcmpi( pch, "ASCM" ) == 0)
              )
         {
         autoSelectCriterion = useMSE;
         }
      else if ( (_strcmpi( pch, "AutoSelectFolds" ) == 0)
             || (_strcmpi( pch, "ASF" ) == 0)
             || (_strcmpi( pch, "AutoSelectFoldCount" ) == 0)
             || (_strcmpi( pch, "ASFC" ) == 0)
              )
         {
         // Get the number of folds to use during AutoSelect processing
         if ( ++iArg >= argc )
            {
            Fatal( "-autoSelectFolds must be followed by an integer argument." );
            }
         int i;
         sscanf( argv[ iArg ], "%i", &i );
         if ( i < 2 )
            {
            // Fold Count must be 2 or more, upto the number of individuals
            Fatal( "-autoSelectFolds must be 2 or more.  Found [%d]", i );
            }
         autoSelectFoldCount = i;
         }
      else if ( (_strcmpi( pch, "AutoSelectApplyThinning" ) == 0)
             || (_strcmpi( pch, "ASAT" ) == 0)
              )
         {
         fAutoSelectApplyThinning = true;
         }
      else if ( (_strcmpi( pch, "AutoSelectApplyThinning-" ) == 0)
             || (_strcmpi( pch, "ASAT-" ) == 0)
              )
         {
         fAutoSelectApplyThinning = false;
         }
      else if ( (_strcmpi( pch, "AutoSelectMaxSizeForThinning" ) == 0)
             || (_strcmpi( pch, "ASMSFT" ) == 0)
              )
         {
         // Get the max binsize for Thinning
         if ( ++iArg >= argc )
            {
            Fatal( "-autoSelectMaxSizeForThinning must be followed by an integer argument." );
            }
         int i;
         sscanf( argv[ iArg ], "%i", &i );
         if ( i < 100 )
            {
            // Fold Count must be 2 or more, upto the number of individuals
            Fatal( "-autoSelectMaxSizeForThinning must be 100 or more.  Found [%d]", i );
            }
         autoSelectMaxSizeForThinning = i;
         fAutoSelectApplyThinning = true;
         }
      else if ( (_strcmpi( pch, "AutoSelectRegionDistance" ) == 0)
             || (_strcmpi( pch, "ASRD" ) == 0)
              )
         {
         // get the distance
         if ( (++iArg >= argc) )
            {
            Fatal( "-autoSelectRegionDistance must be followed by a numeric argument" );
            }
         if ( !FIsFpString( argv[ iArg ] ) )
            {
            Fatal( "ExcludeByGeneticDistance must be followed by a numeric argument.  Found [%s]", argv[ iArg ] );
            }
         real distance = atof( argv[ iArg ] );
         if ( distance <= 0.0 )
            {
            Fatal( "-autoSelectRegionDistance argument must be greater than 0.0.  Found %f", distance );
            }
         autoSelectRegionDistance = distance;
         fAutoSelectApplyThinning = true;
         }
      else if ( (_strcmpi( pch, "file" ) == 0 )
             || (_strcmpi( pch, "plink" ) == 0 )
              )
         {
         // Get the basename for the PLINK files
         if ( ++iArg >= argc )
            {
            Fatal( "file needs to be followed by a basefilename." );
            }
         CheckPlinkFileType( plinkFileType, Natural );       // fatal if incompat options
         plinkFilenameRoot = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "bfile" ) == 0)
             || (_strcmpi( pch, "blink" ) == 0)
              )
         {
         // Get the basename for the BLINK files
         if ( ++iArg >= argc )
            {
            Fatal( "bfile needs to be followed by a basefilename." );
            }
         
         CheckPlinkFileType( plinkFileType, Binary );       // fatal if incompat options
         plinkFilenameRoot = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "tfile" ) == 0)
             || (_strcmpi( pch, "tlink" ) == 0)
              )
         {
         // Get the basename for the PLINK transposed files
         if ( ++iArg >= argc )
            {
            Fatal( "tfile needs to be followed by a basefilename." );
            }
         CheckPlinkFileType( plinkFileType, Transposed );       // fatal if incompat options
         plinkFilenameRoot = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "dosage1" ) == 0)
             || (_strcmpi( pch, "dosage2" ) == 0)
             || (_strcmpi( pch, "dfile1" ) == 0)
             || (_strcmpi( pch, "dfile2" ) == 0)
              )
         {
         // Get the basename for the PLINK dosage files
         if ( ++iArg >= argc )
            {
            Fatal( "-%s needs to be followed by a basefilename.", pch );   // pch still points to option
            }
         char *pchT = strchr( pch, '2' );                      // is there a '2' in the string?
         if ( pchT == nullptr )
            {
            CheckPlinkFileType( plinkFileType, Dosage1 );      // fatal if incompat options
            }
         else
            {
            CheckPlinkFileType( plinkFileType, Dosage2 );      // fatal if incompat options
            }
         plinkFilenameRoot = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "afile" ) == 0)
//             || (_strcmpi( pch, "ImputedGenotype" ) == 0)
              )
         {
         // Get the basename for the imputed genotype files
         if ( ++iArg >= argc )
            {
            Fatal( "afile needs to be followed by a basefilename." );
            }
         CheckPlinkFileType( plinkFileType, ImputedGenotype );          // fatal if incompat options
         plinkFilenameRoot = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "pheno" ) == 0) )
         {
         // Get the basename for the PLINK alternate phenotype files
         if ( ++iArg >= argc )
            {
            Fatal( "-pheno needs to be followed by a filename." );
            }
         alternatePhenotypeParameters.fUseAlternatePhenotype = 1;
         alternatePhenotypeParameters.alternatePhenotypeFilename = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "mpheno" ) == 0 ) )
         {
         // Get the index of which phenotype to operate on
         if ( ++iArg >= argc )
            {
            Fatal( "-mpheno must be followed by an integer argument." );
            }
         int i;
         sscanf( argv[ iArg ], "%i", &i );
         alternatePhenotypeParameters.selectedPhenotypeIndex = (size_t)i;
         alternatePhenotypeParameters.fUseAlternatePhenotype = 1;
         }
      else if ( (_strcmpi( pch, "pheno-name" ) == 0) )
         {
         // Get the index of which phenotype to operate on
         if ( ++iArg >= argc )
            {
            Fatal( "-pheno-name needs to be followed by a phenotype name." );
            }
         alternatePhenotypeParameters.fUseAlternatePhenotype = 1;
         alternatePhenotypeParameters.selectedPhenotypeId = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "filesim" ) == 0)
             || (_strcmpi( pch, "plinksim" ) == 0)
              )
         {
         CheckPlinkFileType( plinkFileTypeCovariance, Natural );       // fatal if incompat options
         for (int i=1; (iArg+i < argc) && !IsLeadinChar( *argv[ iArg+i ] ); ++i )
            {
            // Get ALL the basenames for multiple PLINK files
            plinkFilenameRootCovarianceVector.push_back( argv[ iArg+i ] );
            if ( plinkFilenameRootCovarianceVector.size() > 1 )
               {
               Verbose( verboseMajor, "     Additional covariance filename: [%s]", argv[ iArg+i ] );
               }
            }

         iArg += (int)plinkFilenameRootCovarianceVector.size();
         if ( (iArg >= argc) || (plinkFilenameRootCovarianceVector.size()==0) )
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            Fatal( "-fileSim needs to be followed by at least one basefilename." );
            }
         
         if ( plinkFilenameRootCovarianceVector.size()>1 )
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            ProgressNL( "Running fullrank code." );
            forceFullRank = true;            //allow multiple files only with fullrank mode
            }
         }
      else if ( (_strcmpi( pch, "bfilesim-" ) == 0)
             || (_strcmpi( pch, "bfs-" ) == 0)
              )
         {
         fDoIteratingFastLmm = true;
         // allways use full machine memory even in iterating scenario
         //fMaximizeWorkingSet = false;
         forceFullRank = true;
         }
      else if ( (_strcmpi( pch, "LowMem" ) == 0)
             || (_strcmpi( pch, "lm" ) == 0)
              )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "-LowMem needs to be followed by an integer greater than 50." );
            }
         int i;
         sscanf( argv[ iArg ], "%i", &i );
         if ( i < 50 )
            {
            Fatal( "-LowMem needs to be followed by an integer greater than 50. Found [%s]", argv[ iArg ] );
            }
         cIteratingBlockSizeMax = i;
         fDoIteratingFastLmm = true;
         // allways use full machine memory even in iterating scenario
         //fMaximizeWorkingSet = false;
         forceFullRank = true;
         fEmitSorted = false;             // disable sorted output for IteratingFastLmm
         }
      else if ( (_strcmpi( pch, "bfilesim" ) == 0)
             || (_strcmpi( pch, "blinksim" ) == 0)
              )
         {
         CheckPlinkFileType( plinkFileTypeCovariance, Binary );       // fatal if incompat options
         for ( int i=1; (iArg+i < argc) && !IsLeadinChar( *argv[ iArg+i ] ); ++i )
            {
            plinkFilenameRootCovarianceVector.push_back( argv[iArg+i] );         // Get the basenames for the PLINK files
            if ( i > 1 )
               {
               // One file is assumed.  Report any files that are actually appended
            Verbose( verboseMajor, "Append %s to covariance SNP-basefilenames", argv[iArg+i] );
            }
            }

         iArg += (int)plinkFilenameRootCovarianceVector.size();

         if ( (iArg >= argc) || (plinkFilenameRootCovarianceVector.size()==0) )
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            Fatal( "-bfileSim needs to be followed by at least one basefilename." );
            }
         
         if ( plinkFilenameRootCovarianceVector.size() > 1 )
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            ProgressNL( "Running fullrank code." );
            forceFullRank = true;         //allow multiple files only with fullrank mode
            }
         }
      else if ( (_strcmpi( pch, "tfilesim" ) == 0)
             || (_strcmpi( pch, "tlinksim" ) == 0)
              )
         {
         CheckPlinkFileType( plinkFileTypeCovariance, Transposed );       // fatal if incompat options
         for ( int i=1; (iArg+i < argc) && !IsLeadinChar( *argv[ iArg+i ] ); ++i )
            {
            plinkFilenameRootCovarianceVector.push_back( argv[iArg+i] );   // Get the basenames for the PLINK files
            if ( i > 1 )
               {
               // One file is assumed.  Report any files that are actually appended
            Verbose( verboseMajor, "Append %s to covariance SNP-basefilenames", argv[iArg+i] );
            }
            }

         iArg += (int)plinkFilenameRootCovarianceVector.size();
         if ( (iArg >= argc) || (plinkFilenameRootCovarianceVector.size()==0) )
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            Fatal( "-tfileSim needs to be followed by at least one basefilename." );
            }
         if ( plinkFilenameRootCovarianceVector.size() > 1)
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            ProgressNL( "Running fullrank code." );
            forceFullRank=true;     //allow multiple files only with fullrank mode
            }
         }
      else if ( (_strcmpi( pch, "dosage1Sim" ) == 0)
             || (_strcmpi( pch, "dosage2Sim" ) == 0)
             || (_strcmpi( pch, "dfile1Sim" ) == 0)
             || (_strcmpi( pch, "dfile2Sim" ) == 0)
              )
         {
         char *pchT = strchr( pch, '2' );                         // is there a '2' in the string?
         if ( pchT == nullptr )
            {
            CheckPlinkFileType( plinkFileTypeCovariance, Dosage1 );   // fatal if incompat options
            }
         else
            {
            CheckPlinkFileType( plinkFileTypeCovariance, Dosage2 );   // fatal if incompat options
            }

         for ( int i=1; (iArg+i < argc) && !IsLeadinChar( *argv[ iArg+i ] ); ++i )
            {
            plinkFilenameRootCovarianceVector.push_back( argv[iArg+i] );   // Get the basenames for the PLINK files
            if ( i > 1 )
               {
               // One file is assumed.  Report any files that are actually appended
            Verbose( verboseMajor, "Append %s to covariance SNP-basefilenames", argv[iArg+i] );
            }
            }

         iArg += (int)plinkFilenameRootCovarianceVector.size();
         if ( (iArg >= argc) || (plinkFilenameRootCovarianceVector.size()==0) )
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            Fatal( "-%s needs to be followed by at least one basefilename.", pch );    // pch still points to the option
            }
         if ( plinkFilenameRootCovarianceVector.size() > 1)
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            ProgressNL( "Running fullrank code." );
            forceFullRank=true;     //allow multiple files only with fullrank mode
            }
         }
#if 0
      else if ( (_strcmpi( pch, "dfilesim" ) == 0)
             || (_strcmpi( pch, "dosagesim" ) == 0)
              )
         {
         CheckPlinkFileType( plinkFileTypeCovariance, Dosage );       // fatal if incompat options
         for ( int i=1; (iArg+i < argc) && !IsLeadinChar( *argv[ iArg+i ] ); ++i )
            {
            plinkFilenameRootCovarianceVector.push_back( argv[iArg+i] );   // Get the basenames for the PLINK files
            if ( i > 1 )
               {
               // One file is assumed.  Report any files that are actually appended
            Verbose( verboseMajor, "Append %s to covariance SNP-basefilenames", argv[iArg+i] );
            }
            }

         iArg += (int)plinkFilenameRootCovarianceVector.size();
         if ( (iArg >= argc) || (plinkFilenameRootCovarianceVector.size()==0) )
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            Fatal( "-dfileSim needs to be followed by at least one basefilename." );
            }
         if ( plinkFilenameRootCovarianceVector.size() > 1)
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            ProgressNL( "Running fullrank code." );
            forceFullRank=true;     //allow multiple files only with fullrank mode
            }
         }
#endif
      else if ( (_strcmpi( pch, "afilesim" ) == 0)
//             || (_strcmpi( pch, "ImputedGenotypeSim" ) == 0)
              )
         {
         CheckPlinkFileType( plinkFileTypeCovariance, ImputedGenotype );       // fatal if incompat options
         for ( int i=1; (iArg+i < argc) && !IsLeadinChar( *argv[ iArg+i ] ); ++i )
            {
            plinkFilenameRootCovarianceVector.push_back( argv[iArg+i] );   // Get the basenames for the PLINK files
            if ( i > 1 )
               {
               // One file is assumed.  Report any files that are actually appended
            Verbose( verboseMajor, "Append %s to covariance SNP-basefilenames", argv[iArg+i] );
            }
            }

         iArg += (int)plinkFilenameRootCovarianceVector.size();
         if ( (iArg >= argc) || (plinkFilenameRootCovarianceVector.size()==0) )
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            Fatal( "-afileSim needs to be followed by at least one basefilename." );
            }
         if ( plinkFilenameRootCovarianceVector.size() > 1)
            {
            ProgressNL( "Number of Covariance SNP files: %i", plinkFilenameRootCovarianceVector.size() );
            ProgressNL( "Running fullrank code." );
            forceFullRank=true;     //allow multiple files only with fullrank mode
            }
         }
      else if ( (_strcmpi( pch, "missingphenotype" ) == 0) )
         {
         // Get the user defined missingPhenotypeString and 
         //    missingPhenotypeValue for the PLINK files
         if ( ++iArg >= argc )
            {
            Fatal( "missingPhenotype needs to be followed by a numeric argument." );
            }
         missingPhenotypeString = argv[ iArg ];
         missingPhenotypeValue = (real)atof( missingPhenotypeString.c_str() );
         }
      else if ( (_strcmpi( pch, "simlearntype" ) == 0) )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "simLearnType needs to be followed by FULL or ONCE." );
            }
         // get the COV_LEARN_TYPE (and we only accept FULL or ONCE for now)
         if ( (_strcmpi( argv[iArg], "full" ) == 0) )
            {
            covLearnType=FULL;
            }
         else if( (_strcmpi( argv[iArg], "once" ) == 0) )
            {
            covLearnType=ONCE;
            }
         else
            {
            Fatal( "unknown simLearnType :%s\nOnly FULL or ONCE are supported.\n", argv[iArg] );
            }
         }
      else if ( (_strcmpi( pch, "usebiasterm" ) == 0)
             || (_strcmpi( pch, "UBT" ) == 0)
              )
         {
         // run with a bias term and do not zero-mean the data instead?
         useBiasTerm=true;
         }
      else if ( (_strcmpi( pch, "nobiasterm" ) == 0)
             || (_strcmpi( pch, "NBT" ) == 0)
              )
         {
         // run without a bias term? is a little bit faster and produces same result as long as the kernel is zero-mean.
         useBiasTerm=false;
         }
      else if ( (_strcmpi( pch, "simtype" ) == 0) )   //-KernelType.[RRM/COVARIANCE]
         {
         if ( ++iArg >= argc )
            {
            Fatal( "simType needs to be followed by RRM or COVARIANCE." );
            }
            
         // get the COV_LEARN_TYPE (and we only accept FULL or ONCE for now)
         if ( (_strcmpi( argv[iArg], "rrm" ) == 0) )
            {
            kernelType=RRM;
            }
         else if( (_strcmpi( argv[iArg], "covariance" ) == 0) )
            {
            kernelType=COVARIANCE;
            }
         else
            {
            Fatal( "unknown simType :%s\nOnly RRM or COVARIANCE are supported.\n", argv[iArg] );
            }
         }
      else if ( (_strcmpi( pch, "ml" ) == 0) )
         {
         // get the PARAM_LEARNING_TYPE and set associated variables
         SetMlType( ML );
         }
      else if ( (_strcmpi( pch, "reml" ) == 0) )
         {
         // get the PARAM_LEARNING_TYPE and set associated variables
         //   REML is used with an f-test and with a bias term, otherwise the extra terms are affected.
         SetMlType( REML );
         }
      //runGwasType: default = RUN
      else if ( (_strcmpi( pch, "rungwastype" ) == 0) )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "runGwasType needs to be followed by RUN or NORUN." );
            }
         // get the COV_LEARN_TYPE (and we only accept FULL or ONCE for now)
         if ( (_strcmpi( argv[iArg], "run" ) == 0) )
            {
            runGwasType= RUN;
            }
         else if( (_strcmpi( argv[iArg], "norun" ) == 0) )
            {
            runGwasType= NORUN;
            }
         else
            {
            Fatal( "unknown runGwasType :%s\nOnly RUN or NORUN are supported.\n", argv[iArg] );
            }
         }
      //computeEigenType: default = NOSAVE
      else if ( (_strcmpi( pch, "eigen" ) == 0) )
         {
         computeEigenType= LOADFROMDISK;
         // get the eigenDir name
         if ( ++iArg >= argc )
            {
            Fatal( "eigen must be followed by a directoryname" );
            }
         eigenDir = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "eigenout" ) == 0) )
         {
         computeEigenType= SAVETODISK;

         // get the eigenDir name
         if ( ++iArg >= argc )
            {
            Fatal( "eigenOut must be followed by a directoryname" );
            }
         eigenDir = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "groups" ) == 0) )
         {
         // get the groupsFile name
         if ( ++iArg >= argc )
            {
            Fatal( "groups must be followed by a filename" );
            }
         groupsFile = argv[ iArg ];
         }
      else if ( (_strcmpi( pch, "simgroups" ) == 0 ) )
         {
         // get the eigenDir name
         if ( ++iArg >= argc )
            {
            Fatal( "simGroups must be followed by a filename" );
            }
         simGroups = argv[ iArg ];
         }
      // Look for Brent parameters
      else if ( (_strcmpi( pch, "brentStarts" ) == 0 )
             || (_strcmpi( pch, "BSt" ) == 0)
              )
         {
         // get the numberOptimizationStarts
         if ( ++iArg >= argc )
            {
            Fatal( "brentStarts must be followed by an integer argument" );
            }
         sscanf( argv[ iArg ], "%i", &numberOptimizationStartsGlobal );
         }
      else if ( (_strcmpi( pch, "brentMinLogVal" ) == 0)
             || (_strcmpi( pch, "BMnLV" ) == 0)
              )
         {
         // get the minLogVal
         if ( ++iArg >= argc )
            {
            Fatal( "brentMinLogVal must be followed by a numeric argument" );
            }
         if ( !FIsFpString( argv[ iArg ] ) )
            {
            Fatal( "brentMinLogVal must be followed by a numeric argument.  Found [%s]", argv[ iArg ] );
            }
         minLogValGlobal = (real)atof( argv[ iArg ] );
         }
      else if ( (_strcmpi( pch, "brentMaxLogVal" ) == 0)
             || (_strcmpi( pch, "BMxLV" ) == 0)
              )
         {
         // get the maxLogVal
         if ( ++iArg >= argc )
            {
            Fatal( "brentMaxLogVal must be followed by a numeric argument" );
            }
         if ( !FIsFpString( argv[ iArg ] ) )
            {
            Fatal( "brentmaxlogval must be followed by a numeric argument.  Found [%s]", argv[ iArg ] );
            }
         maxLogValGlobal = (real)atof( argv[ iArg ] );
         }
      else if ( (_strcmpi( pch, "brentTol" ) == 0) )
         {
         // get the brentTol
         if ( ++iArg >= argc )
            {
            Fatal( "brentTol must be followed by a numeric argument." );
            }
         if ( !FIsFpString( argv[ iArg ] ) )
            {
            Fatal( "brentTol must be followed by a numeric argument.  Found [%s]", argv[ iArg ] );
            }
         brentTolGlobal = (real)atof( argv[ iArg ] );
         }
      else if ( (_strcmpi( pch, "brentMaxIter" ) == 0) )
         {
         // get the brentMaxIter
         if ( ++iArg >= argc )
            {
            Fatal( "brentMaxIter must be followed by an integer argument." );
            }
         sscanf( argv[ iArg ], "%i", &brentMaxIterGlobal );
         }
      else if ( (_strcmpi( pch, "MaxChromosomeValue" ) == 0)
             || (_strcmpi( pch, "MaxC" ) == 0)
              )
         {
         if ( ++iArg >= argc )
            {
            Fatal( "MaxChromosomeValue must be followed by an integer argument." );
            }
         sscanf( argv[ iArg ], "%i", &maxChromosomeValue );
         }
      else if ( (_strcmpi( pch, "MaximizeWorkingSet" ) == 0)
             || (_strcmpi( pch, "mws" ) == 0)
              )
         {
         fMaximizeWorkingSet = true;
         }
      else if ( (_strcmpi( pch, "MaximizeWorkingSet-" ) == 0)
             || (_strcmpi( pch, "mws-" ) == 0)
              )
         {
         fMaximizeWorkingSet = false;
         }
      else if ( (_strcmpi( pch, "Cluster" ) == 0)
             || (_strcmpi( pch, "cl" ) == 0)
              )
         {
         fClusterIO = true;
         }
      // TODO:  OBSOLETE this option (accept the option but don't set the variable to change it)
      else if ( (_strcmpi( pch, "EmitSnpIndexOnly" ) == 0)
             || (_strcmpi( pch, "esio" ) == 0)
              )
         {
         //fEmitSnpIndexOnly = true;
         }
      // Oops, don't know this option... Report the error
      else
         {
         Fatal( "Unknown option. [%s]\nUse the option -help to view the list of arguments.", argv[ iArg ] );
         }
      }

   /*
    * Now verify that all the proper arguments and options have been read and are consistent.
    */
   if ( iArg < argc )
      {
      // We only have options and parameters, no standalone filenames go back to main.
      Fatal( "The command line is ill formed.  Parameters and options do not match beginning with argument %d[%s].", iArg, argv[ iArg ] );
      }

   VerifyAlternativePhenotypeUsage();
   VerifyAutoSelectUsage();
   VerifySnpFilterUsage( cJobsToRun, iJobToRun );
   VerifyLinRegUsage();
   VerifySnpPairsUsage();
   VerifyBrentUsage();
   VerifyOptimizeForLowMemoryUsage();
   VerifyNormalizationUsage();

   if ( !fAutoSelectSnps && (runGwasType == RUN) )
      {
      if ( plinkFilenameRoot.empty() )
         {
         Fatal( "runGwasType = RUN and missing Plink file information." );
         }
      if ( outputFilename.empty() )
         {
         if ( !plinkFilenameRoot.empty() )
            {
            // Create output name from PLINK root name
            string filename = FileNameFromPath( plinkFilenameRoot );
            outputFilename = filename + ".out.txt";
            }
         }
      }
   if ( (!plinkFilenameRootCovarianceVector.empty()) 
     && ((!snpDataCovFile.empty()) )
      )
      {
      Fatal( "Incompatible options:  Cannot process both -P(T)LinkSim and -snpDataCovFile options" );
      }
   if ( kernelFile.empty() 
//     && snpDataCovFile.empty() 
     && plinkFilenameRootCovarianceVector.empty() 
     && (computeEigenType != LOADFROMDISK) 
     && (!identityKernel) 
     && (!linReg) 
     && (!fUseLogisticRegression)
     && (simGroups.empty() || groupsFile.empty())
     && (!fAutoSelectSnps)
     && (!fDoIteratingFastLmm)
      )
      {
      Fatal( "missing Sim/P(T)LinkSim information. Cannot compute covariance matrix" );
      }
   if ( eigenDir.empty() && computeEigenType!=NOSAVE  )
      {
      Fatal( "missing eigen/eigenOut directory information. Cannot save/load Eigenvectors" );
      }
   if ( fWriteLogFile || fLogMemoryDiagnostics )
         {
      logDir = FullPath( logDir );
      MakeDirectory( logDir );
      }
   if ( !outputFilename.empty() )
      {
      fullOutputFilename = FullPath( outputFilename );
      if ( !FIsFilenameWritable( fullOutputFilename ) )
      {
         Fatal( "Unable to write to output file [%s]", fullOutputFilename.c_str() );
         }
      }

   // abend if we had any problems with parameters
   FatalIfPreviousErrors();

   // and do the last few preparation steps
   PrepareRng();
   if ( fMaximizeWorkingSet )
      {
      MaximizeWorkingSetLimits();
      }
   }

/*
 *  Initialize program by getting all internal variables initialized
 *    and formated with information available now that was not
 *    possible (or easy) to do at compile time.
 */
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
   i64FastLmmCFullVersion = (wFastLmmCVersion * 100000000ll) + wVerBuild;

   // get the host machine name
   unsigned long cchSzComputerName = sizeof( szComputerName );
   if ( !GetComputerName( szComputerName, &cchSzComputerName ) )
      {
      Error( "Unable to retrieve this computer's name...and I want it!" );
      NTLastError();
      }
   }

void PrepareRng()
   {
   // We want to have a consistent, reproducible sequence from
   //   the RNG so we can debug issues, but we want users
   //   to be able to seed the RNG too
   if ( fRngSeedSuppliedByUser )
      {
      Progress( "\n                  Using rand() seed: %lu", randomNumberSeed );
      }
   srand( randomNumberSeed );

#if defined( USE_RANDOM )        // use the C++0x <random> template library
   //
   // Spin up the C++0x RNGs
   //   However, for reproducibility across runs, 
   //   DO NOT generate a random seed
   // 
   //random_device rnd_device;
   //randomNumberSeed = rnd_device();
   if ( fRngSeedSuppliedByUser )
      {
      Progress( "\n                   Using rng<> seed: %lu", randomNumberSeed );
      }

   rng.seed( randomNumberSeed );
#endif
   }
