#if !defined( Utils_h )
#define Utils_h
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
 * Utils - {Utility Routines}
 *
 *         File Name:   Utils.h
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *    Module Purpose:   This file declares a number of utility 
 *                      routines for console applications
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * 'Publish' our class declarations / function prototypes
 */
// Console I/O Utility Routines
void Error( const char *szFmt, ... );
void Fatal( const char *szFmt, ... );
void FatalIfPreviousErrors();
void NTLastError( void );
void Progress( const char *szFmt, ... );
void ProgressNL( const char *szFmt, ... );
void Verbose( const char *szFmt, ... );
void Verbose( int lvl, const char *szFmt, ... );
void Warn( const char *szFmt, ... );

// System & Utility Routines
bool   FIsDigitString( const char *pch );
bool   FIsDigitString( const string& strT );
bool   FIsFpString( const char *pch );
bool   FIsFpString( const string& s );
size_t GetTotalSystemMemory();
void   MaximizeWorkingSetLimits();
string SecondsToString( long long sec );


// File / File system Utility Routines
int    FFileExists( const string& name );
int    FDirExists( const string& name );
bool   FIsFilenameWritable( const string& filename );
string FullPath( const string& filename );
void   MakeDirectory( const string& dirname );
string MakePath( const string& drive, const string& dir, const string& fname, const string& ext );
void   SplitPath( const string& filename, string& drive, string& dir, string& fname, string& ext );
string FileNameFromPath( const string& szPathIn );
string FilenameInsertBeforeFileExt( const string& filenameIn, const string& str, const string& defaultExt = string("") );
bool   FIsCsvFile( const string& fname );    // does the filename end in '.csv'

// Some math routines
/*
 * PvalUsingFisherFTestCDF computes the pvalue with better precision
 *   for small pvalues.  It is code cloned from FisherF.cpp and Beta.cpp
 *   but modified to not subtract really small numbers from 1.0
 *
 * xF  Value at which to compute the Cummulative Distribution Function (CDF)
 * d1  Degree of freedom parameter
 * d2  Degree of freedom parameter
 */
real PvalUsingFisherFTestCdf( real xF, real d1, real d2 );
real ComputeLogDet( real *mat, size_t m );
real ComputeMSQERR( real *TestPh, real *est_mu, size_t n );

/*
 * 'Publish' the globals we define
 */
extern int fReportTimes;
extern int verboseLevel;
const  int verboseQuiet  = 0;    // 0 = no Verbose output
const  int verboseMajor  = 1;    // 1 = Major progress output
const  int verboseMinor  = 2;    // 2 = Minor progress output
const  int verboseDebug  = 4;    // 4 = debugging output

const  int Tab = '\t';
const  int NewLine = '\n';

/*
 *  Define the stubs, translation routines, and 'stuff'
 *   that Windows normally provides so we can support 
 *   a bit of cross platform operation
 */
#if defined( COMPILE_FOR_LINUX )
#include <linux/limits.h>
#include <time.h>
#include <libgen.h>

#if defined( USE_LINUX_V2x_KERNEL )
/*
 *  TODO:
 *  The std::isnan() is not supported by gcc until 4.6.+
 *    RedHat RHEL 6 is still using a 2.6 kernel and gcc 4.4 toolset.
 *    Re-enable stdd::isnan() on RHEL after RHEL moves to a 3.x kernel
 */
inline int  _isnan( real x ) { return( x != x ); }
#else
inline int  _isnan( real x ) { return( std::isnan( x ) ); }
#endif

inline real abs( real x ) { return( fabs(x) ); }
inline int min( int m1, int m2 ) { return( ( m1 < m2 ) ? m1 : m2 ); }
inline real min( real m1, real m2 ) { return( (m1 < m2 ) ? m1 : m2 ); }
inline int max( int m1, int m2 ) { return( ( m1 > m2 ) ? m1 : m2 ); }
inline int _strcmpi( const char *s1, const char *s2 ) { return( strcasecmp( s1, s2 ) ); }
inline int _access( const char *path, int mode ) { return( access( path, mode ) ); }
inline long long _ftelli64( FILE *pfile ) { return( (long long)(ftell( pfile )) ); }
inline int _fseeki64( FILE *pfile, long long offset, int origin ) { return( fseek( pfile, offset, origin ) ); }
inline int QueryPerformanceFrequency( LARGE_INTEGER* lpFrequency )
   {
   //
   //  *nix defines clock_gettime() in terms of nano-seconds so report
   //    nano-seconds even though the 'increment' will not be 1
   //
   timespec t;
   int rc = clock_getres( CLOCK_MONOTONIC, &t );
   if ( rc == -1 )
      {
      Fatal( "Unable to Configure() clock.  Error: [%d] [%s]", errno, strerror( errno ) );
      }
   *lpFrequency = (long long)1 * 1000LL * 1000LL * 1000LL;        // use nano second frequency / timerResolution.Nanoseoconds;
   return( 0 );
   }
inline int QueryPerformanceCounter( LARGE_INTEGER* lpPerformanceCount )
   {
   // move all the nano-seconds into one LARGE_INTEGER
   timespec t;
   clock_gettime( CLOCK_MONOTONIC, &t );
   *lpPerformanceCount = ((LARGE_INTEGER)(t.tv_sec) * (1llu * 1000llu * 1000llu * 1000llu)) + t.tv_nsec;
   return( 0 );
   }

inline int _mkdir( char *szDir ) { return( mkdir( szDir, 0777 ) ); }

int  GetComputerName( char *pchName, unsigned long* pcchMax );
#define  MAX_COMPUTERNAME_LENGTH 15

//char *_fullpath( char *absPath, const char *relPath, size_t maxLength );
#define  _MAX_PATH   260

unsigned int _set_output_format( unsigned int format );
#define  _TWO_DIGIT_EXPONENT  1

#endif   // COMPILE_FOR_LINUX

/*
 *  ColumnMajor means individual column elements are located together in linear
 *    memory.  For a 10x10 matrix:
 *       [row,col] 
 *           [0,0] =>  0  (keep all column 0 together in memory)
 *           [1,0] =>  1
 *           [2,0] =>  2
 *           [0,1] => 10
 */
inline size_t ColumnMajorIndex( size_t cRows, size_t cColumns, size_t iRow, size_t iColumn )
   {
   size_t index = (cRows*iColumn) + iRow;
   return( index );
   }

/*
 *  RowMajor means individual row elements are located together in linear
 *    memory.  For a 10x10 matrix:
 *       [row,col] 
 *           [0,0] =>  0  (keep all row 0 together in memory)
 *           [0,1] =>  1
 *           [0,2] =>  2
 *           [1,0] => 10
 */
inline size_t RowMajorIndex( size_t cRows, size_t cColumns, size_t iRow, size_t iColumn )
   {
   size_t index = (cColumns*iRow) + iColumn;
   return( index );
   }

inline void ZeroRealArray( real *pReal, size_t cReals )
   {
   for ( size_t i=0; i<cReals; ++i )
      {
      *pReal = 0.0;
      ++pReal;
      }
   }

class CTimer
   {
public:
   CTimer() { Init( false ); }
   CTimer( bool fStart ) { Init( fStart ); }
   CTimer( const CTimer& src );
   void Start( bool fRestart = false );
   void Stop();
   bool IsRunning();    // 
   double Elapsed();    // elapsed time in seconds
   double Elapsedms();  // elapsed time in milliseconds
   double Elapsedus();  // elapsed time in microseconds
   double Resolution() { return(1.0/(double)m_frequency); } // timer resolution in seconds.
   string ToString();                                    // returns a formatted string with the elapsed time
   string ToString( const string& formatString );        // returns formated string
   void   Report( const string& formatString );          // print the time interval
   void   Report( int lvl, const string& formatString ); // conditionally prints the time interval

   const CTimer& operator=(const CTimer& src);  // Assignment
   const CTimer& operator-(const CTimer& src);  // subtraction

protected:
   void Init( bool fStart );
   void Copy( const CTimer& src );
   void Configure();

private:
   static long long  m_frequency;   // ticks per second
   long long m_start;               // timer start time.  positive if counting, negative if stopped
   };


inline void CTimer::Init( bool fStart )
   {
   if ( m_frequency == 0LL )
      {
      Configure();
      }

   m_start = 0LL;
   if ( fStart )
      {
      Start();
      }
   }

inline CTimer::CTimer( const CTimer& src )
   {
   Copy( src );
   }

inline void CTimer::Copy( const CTimer& src )
   {
   if ( this == &src )
      {
      return;                       // do not copy onto self
      }
   m_start = src.m_start;
   }

inline void CTimer::Start( bool fReset )
   {
   long long t;
   QueryPerformanceCounter((LARGE_INTEGER*)&t );
   if ( !fReset && (m_start < 0 ) )    // counter stopped w/ time on it
      {
      m_start += t;
      }
   else
      {
      m_start = t;
      }
   }

inline void CTimer::Stop()
   {
   if ( m_start <= 0 )     // already stopped
      {
      return;
      }
   long long t;
   QueryPerformanceCounter((LARGE_INTEGER*)&t );
   m_start -= t;
   }

inline double CTimer::Elapsed()
   {
   CTimer t(*this);
   t.Stop();
   double retVal = ((double)(-t.m_start))/(double)m_frequency; 
   return( retVal );
   }

inline double CTimer::Elapsedms()
   {
   CTimer t(*this);
   t.Stop();
   double retVal = ((-t.m_start)*1000.0)/(double)m_frequency; 
   return( retVal );
   }

inline double CTimer::Elapsedus()
   {
   CTimer t(*this);
   t.Stop();
   double retVal = ((-t.m_start)*1000000.0)/(double)m_frequency; 
   return( retVal );
   }

/*
 * GenSequence is a class that produces a sequence of 
 *   numbers increasing by offset using operator()
 */
class GenSequence
   {
private:
   size_t value;
   size_t offset;

public:
   GenSequence( ) : value(0), offset(1) {}
   GenSequence( size_t initialValue, size_t offsetInitializer = 1 ) : value( initialValue ), offset( offsetInitializer ) {}
   size_t operator() () { size_t t = value; value += offset; return( t ); }     // save initial value, increment, return initial value
   };

/*
 * GenerateRange is a class that determines the lower and upper 
 *   limits of each block in a subdivided set
 */
class GenRange
   {
private:
   float  cElementsPerRange;     // elements in a range (not equal numbers but keeps all ranges within 1 element)
   size_t cElements;             // count of elements in all ranges
   size_t cRanges;               // count of ranges to create
   size_t iRange;                // index of current range
   size_t min;                   // start of current range
   size_t max;                   // non-inclusive end of current range

public:
   inline GenRange( size_t cElements_, size_t cRanges_ )
      {
      cElements = cElements_;
      cRanges = cRanges_;
      cElementsPerRange = (float)cElements_ / (float)cRanges;
      SelectRange( 0 );
      }
   inline void SelectRange( size_t iRange_ )
      {
      iRange = iRange_;
      min = (size_t)(cElementsPerRange * (float)iRange);
      max   = (size_t)(cElementsPerRange * (float)(iRange+1));
      }
   inline size_t begin() { return( min ); }
   inline size_t end() { return( max ); }
   };

/*
 * Some simple unoptimized vector and matrix math routines
 */
void add_diag_mat( real *mat, size_t m, real alpha, real *dest );
void add_square_matrix( real *m1, real *m2, real *m3, real *m4, size_t m, real *dest );
void add_matrix( real *m1, real *m2, size_t m, size_t n, real *dest );
void multi_mat_diag( real *mat, size_t m, size_t n, real *S, real *dest );
//void multi_diag_mat( real *S, real *mat, size_t m, size_t n, real *dest );
void scale_matrix( real *mat, size_t m, size_t n, real alpha );
double squared_sum( real *vec, size_t n );

void vectorSubtract( real *a, real *b, real *dest, size_t n );
void vectorSubtractBack( real *a, real *b, size_t n );
void vectorAddition( real *a, real *c, real *dest, size_t n );
void vectorAddition( real *a, real *b, real *c, real *dest, size_t n );

#endif   // Utils_h
