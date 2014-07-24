#if !defined( Startup_h )
#define Startup_h
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
 * Startup - {Utility Routines}
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
   FULL           = 0
  ,ONCE           = 1
  ,ALT_PARTIAL    = 2
  ,REFIT_EXCLUDE  = 3
   };

enum PERMUTATION_TYPE
   {
   NO_PERM        = 0
  ,PERMUTE_KERNEL = 1
  ,PERMUTE_PHEN   = 2
   };

enum HANDLE_NON_PSD_TYPE
   {
   SUBTRACT_MIN   = 0
  ,SET_TO_ZERO    = 1
   };

enum KERNEL_TYPE
   {
   RRM         = 0
  ,COVARIANCE  = 1
   };

enum TREAT_MISSING
   {
   MEAN     = 0
  ,REMOVE   = 1
   };

enum PARAM_LEARNING_TYPE
   {
   ML    = 0
  ,REML  = 1
   };

enum COMPUTE_EIGEN_TYPE
   {
   NOSAVE         = 0
  ,LOADFROMDISK   = 1
  ,SAVETODISK     = 2
   };

enum RUN_GWAS_TYPE
   {
   RUN   = 0
  ,NORUN = 1
   };

enum PHEN_FORMAT
   {
   PLINK = 0
  ,OLD   = 1
   };

// obsolete?  remove in future since we only support INDIVIDUAL_ONLY
enum INDIVIDUAL_ID_TYPE
  {
  INDIVIDUALID_ONLY       = 0
 ,FAMILYID_INDIVIDUALID   = 1
  };
extern INDIVIDUAL_ID_TYPE individualIdType;

enum PlinkFileType            // what type of PLINK file is used?
   {
   UnknownFileType   = 0
  ,Natural           = 1
  ,Binary            = 2
  ,Transposed        = 3
  ,Dosage1           = 4
  ,Dosage2           = 5
  ,ImputedGenotype   = 6
   };

enum ExcludeSNPfromSimType
   {
   NoDist            = 0
  ,ByID              = 1
  ,ByPosition        = 2
  ,ByGeneticDistance = 3
   };

enum ExcludeFileType
   {
   NoFile      = 0
  ,TrainOnTest = 1
  ,FileSim     = 2
   };

struct ExcludeSnpParams
   {
   ExcludeSNPfromSimType distType;
   ExcludeFileType fileType;
   int bpDist;          //base pair distance (uses chromosome and position)
   double geneticDist;  //genetic distance (uses chromosome and genetic distance)
   };

extern ExcludeSnpParams excludeSnpParams;
inline bool FExcludeSnps() { return( excludeSnpParams.fileType != NoFile ); }
inline bool FExcludeByGeneticDistance() { return( excludeSnpParams.distType == ByGeneticDistance ); }
inline bool FExcludeByPosition() { return( excludeSnpParams.distType == ByPosition ); }

enum NormalizationType     // define method of normalizing the SNPs under test
   {
   Normal = 0              // 
  ,HardyWeinberg = 1       // assume Hardy-Weinberg equilibrium when building RRM/Cov
  ,Beta = 2                // normalize: SNP = rawSNPvalue / Beta::Pdf( MAF, a, b )
   };

struct NormalizationInfo
   {
   NormalizationType  normalizationType;
   real betaA;
   real betaB;
   NormalizationInfo() : normalizationType(Normal), betaA(1.0), betaB(25.0) {}
   };

extern NormalizationInfo normalizationInfo;

/*
 * 'Publish' our function prototypes
 */
void InitializeProgram();
void PrepareRng();
void ProcessCommandLine( int argc, char* argv[] );
void Usage();

/*
 * 'Publish' the globals we define
 */
extern char    szFastLmmCVersion[];
extern int     wFastLmmCVersion;
extern long long i64FastLmmCFullVersion;

extern bool zeroMeanTestSNPs;             // Should we zero-mean the Test or Similarity SNPs? (Has effects on SNP_weight and bias)
extern bool zeroMeanSimilaritySNPs;       //   subtract mean from SNPs? (default true)
extern bool unitVarianceTestSNPs;         // Should we unitVariance the Test or Similarity SNPs? (Has effects on SNP_weight, not on bias)
extern bool unitVarianceSimilaritySNPs;   //   divide SNPs by their standard deviation? (default true)

extern bool   identityKernel;             // run the analysis using an identity kernel
extern PERMUTATION_TYPE permutationMode;
extern HANDLE_NON_PSD_TYPE handle_non_psd;
extern int    setRank;                    // fix rank to a certain value
extern PlinkFileType plinkFileType;
extern PlinkFileType plinkFileTypeCovariance;
extern string missingPhenotypeString;     // default is "-9"
extern real   missingPhenotypeValue;
extern string missingGenotypeString;      // default is "-9"
extern real   missingGenotypeValue;
extern bool   simMissingness;             // build a Similarity matrix from missingness pattern in SNPs default is false

extern string covariatesFile;             // name of file containing Covariate information
inline bool   FUseCovariatesFile() { return( !covariatesFile.empty() ); }

extern string plinkFilenameRoot;          // Use PLINK Ped/Map files
extern vector<string> plinkFilenameRootCovarianceVector;// Use PLINK Ped/Map files for building the covariance matrix

extern struct AlternatePhenotypeParameters alternatePhenotypeParameters;

extern bool runAlternativeModels;

extern string outputFilename;             // store the output filename as given by user or default from input file
extern string fullOutputFilename;         // keep the full path to the output file around
extern bool   fVerboseOutput;
extern string kernelFile;
extern string snpDataCovFile;             //TODO: NOT USED!!!  This is not set anywhere an is always "".  Remove it?
extern string eigenDir;
extern string kernelFileOut;              // save kernel to text file
extern string groupsFile;                 // file that maps individuals to groups
extern string simGroups;                  // file containing the group-wise similarities

/*
 *  Define our AutoSelect options
 */
enum AutoSelectCriterion { useLL = 0, useMSE = 1 };

extern bool   fAutoSelectSnps;               // preprocess the files to determine the 'best' snps to use
extern int    autoSelectFoldCount;           // number of Folds or Groups of individuals used during AutoSelect
extern bool   fAutoSelectApplyThinning;      // apply thinning/filtering to SearchValues
extern int    autoSelectMaxSizeForThinning;  // Restrict thinning to apply only to search values below 2000
extern real   autoSelectRegionDistance;      // width of 'filter' region on each of SNP (1 centiMorgan) (geneticDistance)
extern AutoSelectCriterion autoSelectCriterion;    // default autoSelectCriterion = useLL
extern string autoSelectOutputFilenameRoot;  // root of names for AutoSelect output files.
extern vector<int> autoSelectSearchValues;   // k values to search to yield the highest LL
extern real   mFraction;
extern size_t topk;
extern bool   selectTopK;
extern bool   fWritePedFile;                 // Signal we are to produce a PED file from AutoSelect
extern string pedOutputFile;                 // Root name of plink .PED file produced by AutoSelect
 
extern bool   linReg;
extern bool   fUseLogisticRegression;           // Compute GWAS using Logistic Regression
extern int    numberOptimizationStartsGlobal;   // brent grid search options/parameters
extern int    brentMaxIterGlobal;
extern real   minLogValGlobal;
extern real   maxLogValGlobal;
extern real   brentTolGlobal;
extern KERNEL_TYPE kernelType;
extern bool   usePseudocounts;            // use pseudocounts (Beta-Binomial model)  for computing SNP frequency. default is false.
extern bool   useBiasTerm;                // use a bias term? if false, then all variables are zero-meaned
extern bool   forceFullRank;
extern bool   useWoodbury;
extern bool   use_f_test;                 // compute an F-test instead of LR test?
extern bool   fDoSnpPairs;                // Signal SnpPairs option in effect - compute epistatic effects
extern bool   fDoSnpId1;                  // Signal SnpId1 option in effect - restrict epistatic study to 1 SNP x n SNPs
extern string snpId1Global;               // the snpId1
extern bool   fEmitSnpIndexOnly;          // **OBSOLETE** minimize output.  Display only Snp Index, not the SnpId

extern COV_LEARN_TYPE covLearnType;
extern COMPUTE_EIGEN_TYPE computeEigenType;
extern RUN_GWAS_TYPE runGwasType; 
extern PARAM_LEARNING_TYPE mlType;
extern TREAT_MISSING missingType;
extern real   logDeltaPreset;             // if set, this will be used instead of optimization
extern bool   useHeritability;            // optionally optimize over h^2 instead of logdelta. 
                                          //   by default do optimization in logdelta space

class SnpFilterOptions;                   // Declare the class

// Process / store the SnpFilterOptions for the data and genetic similarity matrices
extern SnpFilterOptions snpFilterOpts;
extern SnpFilterOptions simFilterOpts;

extern bool fMaximizeWorkingSet;             // Allow consumption of most of the system ram
extern bool fSetFastLmmRealsPrintPrecision;  // use custom value for output precision? (16)
extern int  fastLmmRealsPrintPrecision;      // storage to hold the current reals print precision
extern bool fNoDosageRangeCheck;             // turn off validation of dosage values
extern real pvaluePrintThreshold;            // p-values must be less than pvaluePrintThreshold to be output
extern bool fEmitSorted;                     // write output sorted by pvalue? (default: true)

// Setup Random Number Generator <RNG>
extern bool                                fRngSeedSuppliedByUser;
extern unsigned long                       randomNumberSeed;

#if defined( USE_RANDOM )        // use the C++0x <random> template library
extern default_random_engine               rng;
extern uniform_real_distribution<double>   drand;
#endif

// debug and conditional reporting flags
extern bool   fCreateReadValidationFiles;
extern bool   fWriteLogFile;
extern string logDir;
extern bool   fClusterIO;                 // minimum 'noise' output and flush frequently

extern bool   fDoIteratingFastLmm;        // triggered by -bfilesim-
extern size_t cIteratingBlockSizeMax;     // block size for iterating

extern bool   fDumpSnpFilterOptions;      //
extern bool   fLogMemoryDiagnostics;
extern bool   fUseReproducibleFP;         // Set the math library to minimize FP variation

#endif   // Startup_h
