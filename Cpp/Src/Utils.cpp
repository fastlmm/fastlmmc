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
 *         File Name:   Utils.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements a number of utility routines for console applications
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
 * Globals we define
 */
int verboseLevel = verboseMajor; // verboseLevel default is to emitMajor messages
int fReportTimes = 0;            // global flag holds bit values to mask which times will be reported.

static int cErrorsMax = 25;
static int cErrors = 0;

long long CTimer::m_frequency = 0ll;

/*
 * The code
 */
// Issue Fatal error message to stderr and exit the program with error(1)
void Fatal( const char *szFmt, ... )
   {
   va_list va;

   fprintf( stderr, "\nFatal Error : " );

   va_start( va, szFmt );
   vfprintf( stderr, szFmt, va );
   va_end( va );

   fprintf( stderr, "\n" );
   fflush( stderr );
   exit( 1 );
   }

// Issue error message to stderr and Fatal() if error count reaches limit
void Error( const char *szFmt, ... )
   {
   va_list va;

   ++cErrors;
   fprintf( stderr, "\nError : " );

   va_start( va, szFmt );
   vfprintf( stderr, szFmt, va );
   va_end( va );
   if ( cErrors >= cErrorsMax )
      {
      Fatal( "Maximum error count exceeded. %d errors", cErrors );
      }
   }


void FatalIfPreviousErrors()
   {
   if ( cErrors )
      {
      Fatal( "Exit with previous errors." );
      }
   }
   

void Warn( const char *szFmt, ... )
   {
   va_list va;

   fprintf( stderr, "\nWarning : " );

   va_start( va, szFmt );
   vfprintf( stderr, szFmt, va );
   va_end( va );
   fflush( stderr );
   }


void Verbose( const char *szFmt, ... )
   {
   if ( verboseLevel < verboseMinor )     // default message level is 'Minor' 
      {
      return;
      }

   va_list va;

   fprintf( stderr, "\n" );

   va_start( va, szFmt );
   vfprintf( stderr, szFmt, va );
   va_end( va );
   }

void Verbose( int lvl, const char *szFmt, ... )
   {
   if ( verboseLevel < lvl )
      {
      return;
      }

   va_list va;

   fprintf( stderr, "\n" );

   va_start( va, szFmt );
   vfprintf( stderr, szFmt, va );
   va_end( va );
   }

void Progress( const char *szFmt, ... )
   {
   va_list va;

   va_start( va, szFmt );
   vfprintf( stderr, szFmt, va );
   va_end( va );
   fflush( stderr );
   }

void ProgressNL( const char *szFmt, ... )
   {
   va_list va;

   fprintf( stderr, "\n" );

   va_start( va, szFmt );
   vfprintf( stderr, szFmt, va );
   va_end( va );
   fflush( stderr );
   }

#if defined ( COMPILE_FOR_WINDOWS )       // using Visual C compiler

void NTLastError( void )
   {
   LPVOID lpMsgBuf;
   DWORD err = GetLastError();
   FormatMessage( FORMAT_MESSAGE_ALLOCATE_BUFFER | 
                     FORMAT_MESSAGE_FROM_SYSTEM | 
                     FORMAT_MESSAGE_IGNORE_INSERTS,
                  NULL,
                  err,
                  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
                  (LPTSTR) &lpMsgBuf,
                  0,
                  NULL );
   fprintf( stderr, "\nNTError : 0x%08X : %s", err, (const char *)lpMsgBuf );
   fflush( stderr );
   LocalFree( lpMsgBuf );
   exit( 1 );
   }

#elif defined( COMPILE_FOR_LINUX )

void NTLastError( void )
   {
   Fatal( "NTLastError stub for fatal error" );
   }

unsigned int _set_output_format( unsigned int format )
   {
   // stub a _set_output_format() routine for systems that don't have it
   static unsigned int old_format = 0;
   unsigned int rc = old_format;
   old_format = format;
   return( rc );
   }

int GetComputerName( char *pchName, unsigned long* pcchMax )
   {
   static char szName[] = "ThisComputer";
   unsigned long l = *pcchMax;
   *pcchMax = sizeof( szName );
   if ( l < sizeof( szName ) )
      {
      return( 0 );
      }
   strcpy( pchName, szName );
   return( 1 );
   }
#else
#error We only support VisualC++ on Windows and GCC on Linux.  Unrecognized compiler environment.  Stop Compilation
#endif
// Validate pch points to a string that is all digits
bool FIsDigitString( const char *pch )
   {
   // false if nullptr or zero length string
   if ( (pch == nullptr)  || (*pch == '\0') )
      {
      return( false );
      }

   // scan past numeric digits
   while ( *pch && ( isdigit( *pch ) ) )
      {
      ++pch;
      }

   // return true if digits for the entire string
   return( *pch == '\0' );
   }

bool FIsDigitString( const string& strT )
   {
   return( FIsDigitString( strT.c_str() ) );
   }


//Validate string is a floating point number
bool FIsFpString( const char *pch )
   {
   bool fMatch = false;
   int state = 0;

   while( *pch )
      {
      switch( state )
         {
      case 0:  // looking for leadin +/- if there
         if ( ( *pch == '-' ) || ( *pch == '+' ) )
            {
            ++pch;
            }
         state = 1;
         break;
      case 1:  // past sign. Get digit or decimal point
         if ( isdigit( *pch ) )
            {
            ++pch;
            state = 2;  // accumulate digits before decimal or exponent
            fMatch = true;
            }
         else if ( *pch == '.' )
            {
            ++pch;
            state = 3;  // Require digit after decimal point
            }
         else
            {
            return( false );
            }
         break;
      case 2:  // accumulate digits before the decimal point.  Get digit, '.', or 'E'
         if ( isdigit( *pch ) )
            {
            ++pch;      // accumulate digits
            }
         else if ( *pch == '.' )
            {
            ++pch;
            state = 4;  // accumulate digits for fractional part after decimal point
            }
         else if ( (*pch == 'e') || (*pch == 'E') )
            {
            ++pch;
            state = 5;  // scan for exponent
            fMatch = false;
            }
         else
            {
            return( false );
            }
         break;
      case 3:  // we have at no digit before decimal point, get at least one digit
         if ( isdigit( *pch ) )
            {
            ++pch;
            state = 4;  // accumulate digits for fractional part after decimal point
            fMatch = true;
            }
         else
            {
            return( false );
            }
         break;
      case 4:  // accumulate digits for fractional part after decimal point
         if ( isdigit( *pch ) )
            {
            ++pch;
            }
         else if ( (*pch == 'e') || (*pch == 'E') )
            {
            ++pch;
            state = 5;  // scan for exponent
            fMatch = false;
            }
         else
            {
            return( false );
            }
         break;
      case 5:  // scan for exponent.  Look for lead-in +/-
         if ( ( *pch == '-' ) || ( *pch == '+' ) )
            {
            ++pch;
            }
         state = 6;
         break;
      case 6:  // get first exponent digit
         if ( isdigit( *pch ) )
            {
            ++pch;
            state = 7;
            fMatch = true;
            }
         else
            {
            return( false );
            }
         break;
      case 7:  // accumulate exponent digits
         if ( isdigit( *pch ) )
            {
            ++pch;
            }
         else
            {
            return( false );
            }
         break;
      default:
         Fatal( "Programming Error giving invalid state[%d] in FValidateFpString().  Expected 0-7", state );
         break;
         }
      }
   return( fMatch );
   }

bool FIsFpString( const string& str )
   {
   return( FIsFpString( str.c_str() ) );
   }

/*
 *  File Handling helper routines
 */
#if defined ( COMPILE_FOR_WINDOWS )    // using Visual C compiler
#include <direct.h>
#endif

/*
 *  Check for existence of a file.  
 *  Return 1 (true) if file exists
 *  return 0 (false) if no file or it is a directory
 */
static int FCheckFile( const string& name, int mask )
   {
   struct stat status;
   int rc = stat( name.c_str(), &status );
   if ( rc == 0 )
      {
      // status filled out successfully
      if ( (status.st_mode & mask) == static_cast< decltype( status.st_mode ) >( mask ) )
         {
         // file exists and mode matches mask requirements
         return( true );
         }
      }
   else if ( rc == EINVAL )
      {
      Fatal( "Failure in FCheckFile( %s, 0x%04X )"
           "\n  CRT Error: %d : %s", name.c_str(), mask, errno, strerror( errno ) );
      }
   return( false );
   }

bool FIsFilenameWritable( const string& filename )
   {
   static const char szConsole[] = "/dev/stdout";

   bool result = false;
   if ( filename.size() == 0 )
      {
      Fatal( "FIsFileWritable( \"\" ) - No filename to be validated as Writable." );
      }

   if ( _strcmpi( filename.c_str(), szConsole ) == 0 )
      {
      return( true );
      }

   struct stat status;
   int rc = stat( filename.c_str(), &status );
   if ( rc == EINVAL )
      {
      Fatal( "Failure in FIsFileWritable( %s )"
           "\n  CRT Error: %d : %s", filename.c_str(), errno, strerror( errno ) );
      }
   else if ( rc == 0 )
      {
      // the name in the file namespace is used!
      //  File or directory exists.  Is it a writeable file?
      if ( (status.st_mode & (S_IFREG | S_IWRITE)) == (S_IFREG | S_IWRITE) )
         {
         result = true;
         }
      }
   else 
      {
      // stat says filename does not exist in the file namespace.
      //  Can we create it write to it?
      FILE *pFile = fopen( filename.c_str(), "wt" );          // append ascii
      if ( pFile )
         {
         fclose( pFile );
         remove( filename.c_str() );         
         result = true;
         }
      else
         {
         // Cannot create the file
         Fatal( "Failure creating file FIsFileWritable( %s )"
              "\n  CRT Error: %d : %s", filename.c_str(), errno, strerror( errno ) );
         }
      }
   return( result );
   }

#if 0
int FWriteableFile( const string& name )
   {
   return( FCheckFile( name, S_IFREG | S_IWRITE ) );
   }

int FWriteableDir( const string& name )
   {
   return( FCheckFile( name, S_IFDIR | S_IWRITE ) );
   }
#endif

int FFileExists( const string& name )
   {
   return( FCheckFile( name, S_IFREG ) );
   }

int FDirExists( const string& name )
   {
   return( FCheckFile( name, S_IFDIR ) );
   }

void MakeDirectory( const string& dirname )
   {
   if ( _access( dirname.c_str(), 0 ) == 0 )
      {
      // The name exists in the filesystem
      struct stat status;
      stat( dirname.c_str(), &status );
      if ( !(status.st_mode & S_IFDIR) )
         {
         Fatal( "File already exists.  Cannot create directory: [%s]", dirname.c_str() );
         }    
      }    
   else
      {
      int rc = _mkdir( (char *)dirname.c_str() );
      if ( rc == -1 )
         {
         Fatal( "Unable to create directory: [%s]"
              "\n  CRT Error: %d : %s", dirname.c_str(), errno, strerror( errno ) );
         }
      }
   }

string FullPath( const string& filename )
   {
#if defined ( COMPILE_FOR_WINDOWS )    // using Visual C compiler
   // Use Windows _fullpath() to expand the name
   char szFullPath[ _MAX_PATH ];
   _fullpath( szFullPath, filename.c_str(), sizeof( szFullPath ) );
   return( szFullPath );
#elif defined( COMPILE_FOR_LINUX )
   // Use Linux realpath() to expand the name
   char szRealPath[ PATH_MAX ];
   char *szT = realpath( filename.c_str(), szRealPath );
   if ( szT != nullptr )
      {
      return( string( szRealPath ) );
      }
   else if ( errno == ENOENT )  // was it a file not found error?
      {
      /*
       * Yes.  The file was not found, so see if the parent directory
       *    exists by splitting the path and filename parts to verify
       *    we can create the file later
       */
      char *szT1 = strdup( filename.c_str() );
      char *szDirName = dirname( szT1 );
      szT = realpath( szDirName, szRealPath );
      if ( szT != nullptr )
         {
         // valid directory, reconstruct the fullpath
         strcpy( szT1, filename.c_str() );
         char *szBasename = basename( szT1 );
         strcat( szRealPath, "/" );
         strcat( szRealPath, szBasename );
         free( szT1 );
         return( string( szRealPath ) );
         }
      free( szT1 );
      }

   // Rats!  If we get here we have some error other 
   //   than the initial file not found
   Fatal( "Unable to resolve realpath() for filename: [%s]"
        "\n       errno: %6d"
        "\n  strerror(): %s",
        filename.c_str(),
        errno,
        strerror( errno ) );
#endif
   }

string MakePath( const string& drive, const string& dir, const string& fname, const string& ext )
   {
#if defined ( COMPILE_FOR_WINDOWS )       // using Visual C compiler
   char szFullPath[ _MAX_PATH ];
   _makepath( szFullPath, drive.c_str(), dir.c_str(), fname.c_str(), ext.c_str() );
   return( szFullPath );
#else
   string t;
   if ( drive.size() > 0 )
      {
      if ( drive[0] != '/' )
         {
         t = "/";
         }
      t += drive;
      }
   if ( dir.size() > 0 )
      {
      if ( (t.size() > 0) && (dir[0] != '/') )
         {
         t += "/";
         }
      t += dir;
      }
   if ( fname.size() > 0 )
      {
      if ( (t.size() > 0) && (fname[0] != '/') )
         {
         t += "/";
         }
      t += fname;
      }
   if ( ext.size() > 0 )
      {
      t += "." + ext;
      }
   return( t );
#endif
   }

void SplitPath( const string& filename, string& drive, string& dir, string& fname, string& ext )
   {
#if defined ( COMPILE_FOR_WINDOWS )    // using Visual C compiler
   char szDrive[ _MAX_DRIVE ];
   char szDir[ _MAX_DIR ];
   char szFname[ _MAX_FNAME ];
   char szExt[ _MAX_EXT ];

   _splitpath( filename.c_str(), szDrive, szDir, szFname, szExt );
   drive = szDrive;
   dir = szDir;
   fname = szFname;
   ext = szExt;
#else
   // nothing like this on Linux. :-(
   //char *szDrive;
   char *szDir;
   char *szFname;
   char *szExt;

   char *szT = strdup( filename.c_str() );
   char *pchSlash = strrchr( szT, '/' );
   if ( pchSlash == nullptr )
      {
      // it is all the filename
      szFname = szT;
      }
   else
      {
      // we have a path too
      *pchSlash = '\0';
      szDir = szT;
      szFname = pchSlash + 1;
      }
   char *pchDot = strrchr( szFname, '.' );
   if ( pchDot != nullptr )
      {
      *pchDot = '\0';
      szExt = pchDot + 1;
      }

   drive = "";
   dir = szDir;
   fname = szFname;
   ext = szExt;
   free( szT );
#endif
   }

string FileNameFromPath( const string& szPathIn )
   {
#if defined ( COMPILE_FOR_WINDOWS )    // using Visual C compiler
   char szFullPath[ _MAX_PATH ];
   char szDrive[ _MAX_DRIVE ];
   char szDir[ _MAX_DIR ];
   char szFName[ _MAX_FNAME ];
   char szExt[ _MAX_EXT ];
   //char szPath[ _MAX_PATH ];
   char szModule[ _MAX_PATH ];

   _fullpath( szFullPath, szPathIn.c_str(), sizeof( szFullPath ) );
   _splitpath( szFullPath, szDrive, szDir, szFName, szExt );
   _makepath( szModule, nullptr, nullptr, szFName, szExt );
   string filename( szModule );
   return( filename );
#else
   /*
    * The Linux function basename() _CAN_ modify the string passed
    *   into the function.  Create temp copies for basename to work 
    *   with first.
    */
   char *szT = strdup( szPathIn.c_str() );
   char *szBasename = basename( szT );
   string filename( szBasename );
   free( szT );
   return( filename );
#endif
   }

string FilenameInsertBeforeFileExt( const string& filenameIn, const string& str, const string& defaultExt )
   {
   string strNil;
   string strDrive;
   string strDir;
   string strFname;
   string strExt;

   SplitPath( filenameIn, strDrive, strDir, strFname, strExt );
   string strT = MakePath( strDrive, strDir, strFname, str );
   string strT2 = MakePath( strNil, strNil, strT, strExt.empty() ? defaultExt : strExt ); 
   return( strT2 );
   }

bool FIsCsvFile( const string& fname )
   {
   static char szCsvFileType[] = ".csv";
   size_t cchFileType = strlen( szCsvFileType );

   size_t cch = fname.size();
   if ( cch > cchFileType )
      {
      const char *pch = fname.c_str() + fname.size() - cchFileType;
      if ( _strcmpi( pch, szCsvFileType ) == 0 )
         {
         return( true );
         }
      }
   return( false );
   }


/*
 *  Higher resolution timing routines
 */
void CTimer::Configure()
   {
   QueryPerformanceFrequency((LARGE_INTEGER*)(&m_frequency));
   if ( !m_frequency )
      {
      m_frequency = 1;
      }
   // This is where you'd remove overhead of start/stop routines, but
   //  it is too jittery
   }

string CTimer::ToString( const string& formatString )
   {
   char rgch[256];

   sprintf( rgch, formatString.c_str(), ToString().c_str() );
   string val( rgch );
   return( val );
   }

string CTimer::ToString()
   {
   char rgch[ 64 ];

   CTimer timer(*this);
   timer.Stop();
   long long usecT = (-timer.m_start * 1000000LL) / m_frequency;
   int  usec = (int)(usecT % 1000LL);
   int  msec = (int)((usecT / 1000LL) % 1000LL);
   int  seconds = (int)((usecT / (1000LL*1000LL)) % 60LL);
   int  minutes = (int)((usecT / (1000LL*1000LL*60LL)) % 60LL);
   int  hours = (int)((usecT / (1000LL*1000LL*60LL*60LL)));
   if ( hours )
      {
      sprintf( rgch, "%2d:%02d hrs", hours, minutes );
      }
   else if ( minutes )
      {
      sprintf( rgch, "%2d:%02d mins", minutes, seconds );
      }
   else if ( seconds )
      {
      sprintf( rgch, "%2d.%03d sec", seconds, msec );
      }
   else if ( msec )
      {
      sprintf( rgch, "%3d.%03d ms", msec, usec );
      }
   else
      {
      sprintf( rgch, "%3d us", usec );
      }
   string s(rgch);
   return( s );
   }

void CTimer::Report( const string& formatString )
   {
   ProgressNL( formatString.c_str(), ToString().c_str() );
   }

void CTimer::Report( int rptLevel, const string& formatString )
   {
   if ( (fReportTimes & rptLevel) || (rptLevel == -1) )
      {
      ProgressNL( formatString.c_str(), ToString().c_str() );
      }
   }

string SecondsToString( long long sec )
   {
   char rgch[ 64 ];

   int  seconds = (int)(sec % 60LL);
   int  minutes = (int)((sec / 60LL) % 60ll);
   int  hours = (int)((sec / (60LL*60LL)));
   if ( hours )
      {
      sprintf( rgch, "%2d:%02d hrs", hours, minutes );
      }
   else if ( minutes )
      {
      sprintf( rgch, "%2d:%02d mins", minutes, seconds );
      }
   else
      {
      sprintf( rgch, "%2d.%03d sec", seconds, 0 );
      }
   string s(rgch);
   return( s );
   }


/*
 * Other misc helper routines
 */


/*
 *  Debug Dump routines
 */
void BinaryWrite( string fname, void* pv, size_t cb )
   {
   FILE *ofile = fopen( fname.c_str(), "wb");
   if ( ofile == nullptr )
      {
      Fatal( "Unable to open output file [%s]", fname.c_str() );
      }
   if ( pv )
      {
      fwrite( pv, 1 , cb, ofile );
      }
   fclose( ofile );
   }

static const char *openModeWrite = "w";
static const char *openModeWriteAppend = "a";

void DumpLmmGwasInputs( const StudyData *pStudyData, const EigenSym *pEigenSym, const ResultsLMM *pResultsLmm, LogFileWriteMode fOpenAppend )
   {
   if ( !fWriteLogFile )
      {
      return;
      }

   string logFile = MakePath( "", logDir, "FastLmm.GwasInputs", "Log" );
   const char *openMode = (fOpenAppend == Append) ? openModeWriteAppend : openModeWrite;
   FILE *pf = fopen( logFile.c_str(), openMode );
   if ( pf == nullptr )
      {
      Fatal( "Unable to open log file: [%s]"
         "\n  CRT Error [%d]: %s", logFile.c_str(), errno, strerror( errno ) );
      }

   string commentString = "(DumpGwasInputs)";
   pStudyData->Dump( pf, commentString );
   pEigenSym->Dump( pf, commentString );
// pResultsLmm->Dump( pf, commentString );
   fclose( pf );
   }

void DumpStudyData( const StudyData *pStudyData, const string& comment, LogFileWriteMode fOpenAppend )
   {
   if ( !fWriteLogFile )
      {
      return;
      }

   string logFile = MakePath( "", logDir, "FastLmm.StudyData", "Log" );
   const char *openMode = (fOpenAppend == Append) ? openModeWriteAppend : openModeWrite;
   FILE *pf = fopen( logFile.c_str(), openMode );
   if ( pf == nullptr )
      {
      Fatal( "Unable to open log file: [%s]"
         "\n  CRT Error [%d]: %s", logFile.c_str(), errno, strerror( errno ) );
      }

   pStudyData->Dump( pf, comment );
   fclose( pf );
   }

void DumpEigenSymState( const EigenSym *pEigenSym, const string& comment, LogFileWriteMode fOpenAppend )
   {
   if ( !fWriteLogFile )
      {
      return;
      }

   string logFile = MakePath( "", logDir, "FastLmm.EigenState", "Log" );
   const char *openMode = (fOpenAppend == Append) ? openModeWriteAppend : openModeWrite;
   FILE *pf = fopen( logFile.c_str(), openMode );
   if ( pf == nullptr )
      {
      Fatal( "Unable to open log file: [%s]"
         "\n  CRT Error [%d]: %s", logFile.c_str(), errno, strerror( errno ) );
      }

   pEigenSym->Dump( pf, comment );
   fclose( pf );
   }

void DumpKernelArray( const real *pReal, size_t cRows, size_t cColumns, const string& comment, LogFileWriteMode fOpenAppend )
   {
   if ( !fWriteLogFile )
      {
      return;
      }

   string logFile = MakePath( "", logDir, "FastLmm.KernelArray", "Log" );
   const char *openMode = (fOpenAppend == Append) ? openModeWriteAppend : openModeWrite;
   FILE *pf = fopen( logFile.c_str(), openMode );
   if ( pf == nullptr )
      {
      Fatal( "Unable to open log file: [%s]"
         "\n  CRT Error [%d]: %s", logFile.c_str(), errno, strerror( errno ) );
      }

   fprintf( pf, "\nKernelArray: [%s]", comment.c_str() );
   DumpRealArray( pf, 2, pReal, "Kernel", cRows, cColumns );
   fclose( pf );
   }

void DumpStringVector( FILE* pf, int indent, const vector<string>& v, const string& name )
   {
#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
   const char *szFmt1 = "\n%*svector<string> %s: [%Iu]";
   const char *szFmt2 = "\n%*s[%Iu]: %s";
#else                      // Linux/g++ uses a %zu specifier for size_t
   const char *szFmt1 = "\n%*svector<string> %s: [%zu]";
   const char *szFmt2 = "\n%*s[%zu]: %s";
#endif

   size_t cStrings = v.size();
   fprintf( pf, szFmt1, indent, "", name.c_str(), cStrings );
   for( size_t iString=0; iString<cStrings; ++iString )
      {
      fprintf( pf, szFmt2, indent+2, "", iString, v[iString].c_str() );
      }
   }

void DumpSnpInfoVector( FILE *pf, int indent, const vector<SnpInfo>& v, const string& name )
   {
#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
   const char *szFmt1 = "\n%*svector<SnpInfo> %s: [%Iu]";
   const char *szFmt2 = "\n%*s[%Iu]: idSnp: %s";
#else                      // Linux/g++ uses a %zu specifier for size_t
   const char *szFmt1 = "\n%*svector<SnpInfo> %s: [%zu]";
   const char *szFmt2 = "\n%*s[%zu]: idSnp: %s";
#endif

   size_t c = v.size();
   fprintf( pf, szFmt1, indent, "", name.c_str(), c );
   for( size_t iSnpInfo=0; iSnpInfo<c; ++iSnpInfo )
      {
      fprintf( pf, szFmt2, indent+2, "", iSnpInfo, v[iSnpInfo].idSnp.c_str() );
      fprintf( pf, "\n%*sChromosome: %d [%s]", indent+4, "", v[iSnpInfo].iChromosome, v[iSnpInfo].idChromosome.c_str() );
      fprintf( pf, "\n%*sMorgans: %e", indent+4, "", v[iSnpInfo].geneticDistance );
      fprintf( pf, "\n%*sPosition: %d", indent+4, "", v[iSnpInfo].basepairPosition );
      fprintf( pf, "\n%*sAlleles: %c %c", indent+4, "", v[iSnpInfo].minorAllele, v[iSnpInfo].majorAllele );
      }
   }

void DumpRealArray( FILE *pf, int indent, const real *p, const string& name, size_t cRows, size_t cColumns )
   {
#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
   const char *szFmt1 = "\n%*sreal* %s: [cRows: %Iu cColumns: %Iu]";
   const char *szFmt2 = "\n%*s[%Iu]:";
#else                      // Linux/g++ uses a %zu specifier for size_t
   const char *szFmt1 = "\n%*sreal* %s: [cRows: %zu cColumns: %zu]";
   const char *szFmt2 = "\n%*s[%zu]:";
#endif

   fprintf( pf, szFmt1, indent, "", name.c_str(), cRows, cColumns );
   if ( p == nullptr )
      {
      fprintf( pf, "  ** nullptr **" );
      }
   else
      {
      for ( size_t iRow=0; iRow<cRows; ++iRow )
         {
         fprintf( pf, szFmt2, indent+2, "", iRow );
         for ( size_t iColumn=0; iColumn<cColumns; ++iColumn )
            {
            fprintf( pf, " %.17E", p[ColumnMajorIndex( cRows, cColumns, iRow, iColumn )] );
            }
         }
      }
   }

void Dumpreals( FILE *pf, real *pdbl, size_t count )
   {
#if defined( COMPILE_FOR_WINDOWS )     // Windows will zero fill w/ 0x%016p
   const char *szFmt1 = "\n 0x%016p: ";
#else                      // gcc/Linux doesn't like 0x%016p.  Use 0x%16p
   const char *szFmt1 = "\n 0x%16p: ";
#endif
   int i = 0;
   while ( count-- )
      {
      if ( i==0 )
         {
         fprintf( pf, szFmt1, pdbl );
         i = 8;
         }
      fprintf( pf, " %8.3e", *pdbl );
      ++pdbl;
      --i;
      }
   fprintf( pf, "\n" );
   }

void DumpGroupSim( FILE *pf, int indent, const GroupSim* p, const string& name )
   {
   fprintf( pf, "\n%*sGroupSim* %s: ", indent, "", name.c_str() );
   if ( p == nullptr )
      {
      fprintf( pf, "  ** nullptr **" );
      }
   else
      {
      // Dump the rest of GroupSim Class members
      fprintf( pf, "\n%*s%s", indent+2, "", " *** NYI DumpGroupSim ***" );
      }
   }

/*
 *  Write out the CRT heap state to the log file with the 'comment'
 *    prepended to beginning of the output.
 */
void LogMemoryState( const string& comment, LogFileWriteMode fOpenAppend )
   {
#if defined( _DEBUG )
   if ( !fLogMemoryDiagnostics )
      {
      return;
      }

   _CrtMemState crmT;
   _CrtMemCheckpoint( &crmT );

   string logFile = MakePath( "", logDir, "FastLmm.MemoryUsage", "Log" );
   const char *openMode = (fOpenAppend == Append) ? openModeWriteAppend : openModeWrite;
   FILE *pf = fopen( logFile.c_str(), openMode );
   if ( pf == nullptr )
      {
      Fatal( "Unable to open memory usage log file: [%s]"
         "\n  CRT Error [%d]: %s", logFile.c_str(), errno, strerror( errno ) );
      }

   // TODO:  Resolve differences with CRT about meaning of 
   //          lHighWaterCount
   size_t cbNow = 0;
   size_t cBlocksNow = 0;
   for ( int i=0; i<_MAX_BLOCKS; ++i )
      {
      cbNow += crmT.lSizes[ i ];
      cBlocksNow += crmT.lCounts[ i ];
      }
   double dSize = ((double)cbNow) / (1024.0 * 1024.0);

   if ( fOpenAppend == OverWrite )
      {
      fprintf( pf, "%-32.32s"
                   "\t%-12.12s"
                   "\t%-12.12s"
                   "\t%-12.12s"
                   "\t%-12.12s"
                   , "Location"
                   , "High Mem"
                   , "Total Mem"
                   , "cbNow"
                   , "cBlocksNow"
                   );
      }
   fprintf( pf, "\n%-32.32s"
                "\t%9.3lf MB"
                "\t%9.3lf MB"
                "\t%9.3lf MB"
                "\t%12ld"
                , comment.c_str()
                , ((double)(crmT.lHighWaterCount))/(1024.0 * 1024.0)
                , ((double)(crmT.lTotalCount))/(1024.0 * 1024.0)
                , dSize
                , cBlocksNow
                );
   fclose( pf );
#endif
   }

/*
 * GetTotalSystemMemory() returns the amount of physical RAM in bytes
 */
size_t GetTotalSystemMemory()
   { 
#if defined( COMPILE_FOR_WINDOWS )     // using Visual C compiler on Windows
   MEMORYSTATUSEX status;
   status.dwLength = sizeof(status);
   GlobalMemoryStatusEx(&status);
   return( (size_t)status.ullTotalPhys );
#else                            // using GCC and Linux
   long pages = sysconf(_SC_PHYS_PAGES); 
   long page_size = sysconf(_SC_PAGE_SIZE); 
   return pages * page_size; 
#endif
   } 

/*
 *  Set the process to be able to consume 'most' of the physical memory
 *    and not restrict our usage and page against ourselves
 */
void MaximizeWorkingSetLimits()
   {
#if defined( COMPILE_FOR_WINDOWS )
   // conversion from bytes to giga-bytes
   static const double bytesPerGb = (1024.0 * 1024.0 * 1024.0);

   // Remove hard working set limits for FastLmm and setup to use
   //  _most_ of the available RAM.
   //  How much memory / free memory is in the system?
   MEMORYSTATUSEX memoryStatusEx;
   memoryStatusEx.dwLength = sizeof( memoryStatusEx );
   GlobalMemoryStatusEx( &memoryStatusEx );
   //SIZE_T cbPhysicalRam = memoryStatusEx.ullTotalPhys;
   //SIZE_T cbFreeRam = memoryStatusEx.ullAvailPhys;

   Verbose( verboseMajor, "                      Memory Status:" );
   Verbose( verboseMajor, "                 Physical Ram Found: %8.3lf GB", (double)memoryStatusEx.ullTotalPhys / bytesPerGb );
   Verbose( verboseMajor, "                  Physical Ram Free: %8.3lf GB", (double)memoryStatusEx.ullAvailPhys / bytesPerGb );

   HANDLE hProcess;
   SIZE_T cbMinimumWorkingSetSize;
   SIZE_T cbMaximumWorkingSetSize;
   DWORD  workingSet_Flags;

   hProcess = GetCurrentProcess();
   GetProcessWorkingSetSizeEx( hProcess, &cbMinimumWorkingSetSize, &cbMaximumWorkingSetSize, &workingSet_Flags );
   Verbose( verboseMajor, "          Initial MinimumWorkingSet: %8.3lf GB", (double)cbMinimumWorkingSetSize / bytesPerGb );
   Verbose( verboseMajor, "          Initial MaximumWorkingSet: %8.3lf GB", (double)cbMaximumWorkingSetSize / bytesPerGb );

   workingSet_Flags &= ~QUOTA_LIMITS_HARDWS_MAX_ENABLE;
   workingSet_Flags |=  QUOTA_LIMITS_HARDWS_MAX_DISABLE;

   SIZE_T cbNewMaximumWorkingSetSize = (SIZE_T)((double)memoryStatusEx.ullAvailPhys * 0.95);

   if ( (cbNewMaximumWorkingSetSize > cbMaximumWorkingSetSize) 
     && (cbNewMaximumWorkingSetSize < memoryStatusEx.ullAvailPhys )
      )
      {
      // Set max working set to 95% of free memory
      SetProcessWorkingSetSizeEx( hProcess, cbMinimumWorkingSetSize, cbNewMaximumWorkingSetSize, workingSet_Flags );
      GetProcessWorkingSetSizeEx( hProcess, &cbMinimumWorkingSetSize, &cbMaximumWorkingSetSize, &workingSet_Flags );
      }
   
   Verbose( verboseMajor, "              Set MaximumWorkingSet: %8.3lf GB", (double)cbMaximumWorkingSetSize / bytesPerGb );
#else
   // MaximizeWorkingSetLimits() is a NOP in Linux for now
#endif
   }

/*
 * PvalUsingFisherFTestCDF computes the pvalue with better precision
 *   for small pvalues.  It is code cloned from FisherF.cpp and Beta.cpp
 *   but modified to not subtract really small numbers from 1.0
 *
 * xF  Value at which to compute the Cummulative Distribution Function (CDF)
 * d1  Degree of freedom parameter
 * d2  Degree of freedom parameter
 */
real PvalUsingFisherFTestCdf( real xF, real d1, real d2 )
   {
   //  PValue is 1 - CDF.
   //  If the value (xF) for CDF evaluation is 0.0 or 1.0, 
   //    we know the CDF is 0.0 or 1.0 and the 
   //    pvalue is 1.0 - xF, so short circuit
   if ( (xF == 0.0) || (xF == 1.0) )
      {
      return( 1.0 - xF );
      }

   if ( xF < 0.0 )
      {
      Fatal( "Fisher_F::Cdf parameter xF must be > 0.  Found xF=%g", xF );
      }
   if ( (d1 <= 0) || (d2 <= 0) )
      {
      Fatal( "Fisher_F::Cdf parameters d1 and d2 must be > 0.  Found d1=%g d2=%g", d1, d2 );
      }

   real d1x = d1 * xF;
   real xB = d1x / (d1x + d2);
   
   //return Beta::BetaInc( xB, (d1 / 2), (d2 / 2) );
   //real Beta::BetaInc(real xB, real a , real b)
   // So we compute the regularized incomplete beta function now
   real a = d1 / 2;
   real b = d2 / 2;
   real c = exp( MathFunctions::logGamma( a + b ) 
               - MathFunctions::logGamma( a ) 
               - MathFunctions::logGamma( b ) 
               + (a * log( xB ))
               + (b * log( 1 - xB ))
               );

#if 0    // debugging aids (values to inspect in a debugger); unused otherwise
   real rT01 = ( a + 1.0 ) / ( a + b + 2.0 );
   real rMin = std::numeric_limits< real >::min();
   real rEpsilon = std::numeric_limits< real >::epsilon();
   real rEpsilonLD = std::numeric_limits< long double >::epsilon();
#endif
   real rT0a = 0.0;
   real rT0b = 0.0;
   real pVal = 0.0;

   if ( xB < (( a + 1.0 ) / ( a + b + 2.0 )) )
      {
      real cf = Beta::BetaCF( xB, a, b );
      rT0a = cf / a;
      rT0b = c * rT0a;
      pVal = 1 - rT0b;
      //pVal = 1 - (c * cf / a);      // inverse of BetaInc()
      }
   else
      {
      // Use symmetry relation
      real cf = Beta::BetaCF( (1 - xB), b, a );
      rT0a = cf / b;
      pVal = c * rT0a;
      //pVal = c * (cf / b);          // inverse of BetaInc()
      }

   if (pVal >= 1.0)
      {
      pVal = 1.0;
      }
   else if (pVal < 0.0 )
      {
      pVal = 0.0;
      }
   return pVal;
   }

// function assumes that mat is L or L^T from cholesky decomposition
real ComputeLogDet( real *mat, size_t m )
   {
   real det = 0.0;
   for ( size_t i = 0; i < m; ++i )
      {
      size_t idx = ColumnMajorIndex( m, m, i, i );
      det = det + (2 * log( mat[ idx ] ));
      }
   return( det );
   }

real ComputeMSQERR( real *TestPh, real *est_mu, size_t n )
   {
   real result = 0.0;
   for ( size_t i = 0; i < n; ++i )
      {
      real diff = TestPh[ i ] - est_mu[ i ];
      result = result + (diff * diff);
      }
   return( result );
   }


/*
 * Some simple unoptimized vector and matrix math routines
 */
// matrix x diag
void multi_mat_diag( real *mat, size_t m, size_t n, real *S, real *dest )
   {
   for ( size_t j = 0; j < n; ++j )
      {
      real scale = S[ j ];
      for ( size_t i = 0; i < m; ++i )
         {
         size_t idx = ColumnMajorIndex( m, n, i, j );
         dest[ idx ] = mat[ idx ] * scale;
         }
      }
   }

// diag x matrix
void multi_diag_mat( real *S, real *mat, size_t m, size_t n, real *dest )
   {
   for ( size_t i = 0; i < m; ++i )
      {
      real scale = S[ i ];
      for ( size_t j = 0; j < n; ++j )
         {
         size_t idx = ColumnMajorIndex( m, n, i, j );
         dest[ idx ] = mat[ idx ] * scale;
         }
      }
   }

// matrix * alpha
void scale_matrix( real *mat, size_t m, size_t n, real alpha )
   {
   size_t limit = m * n;
   for ( size_t i = 0; i < limit; ++i )
      {
      mat[ i ] *= alpha;
      }
   }

// matrix + diag(a)
// assume that mat is square 
void add_diag_mat( real *mat, size_t m, real alpha, real *dest )
   {
   for ( size_t i = 0; i < m; ++i )
      {
      size_t idx = ColumnMajorIndex( m, m, i, i );
      dest[ idx ] = mat[ idx ] + alpha;
      }
   }

void add_square_matrix( real *m1, real *m2, real *m3, real *m4, size_t m, real *dest )
   {
   size_t limit = m * m;
   for ( size_t i = 0; i < limit; ++i )
      {
      dest[ i ] = m1[ i ] + m2[ i ] + m3[ i ] + m4[ i ];
      }
   }

void add_matrix( real *m1, real *m2, size_t m, size_t n, real *dest )
   {
   size_t limit = m * n;
   for ( size_t i = 0; i < limit; ++i )
      {
      dest[ i ] = m1[ i ] + m2[ i ];
      }
   }

void vectorSubtract( real *a, real *b, real *dest, size_t n )
   {
   for ( size_t i = 0; i < n; ++i )
      {
      dest[ i ] = a[ i ] - b[ i ];
      }
   }

void vectorSubtractBack( real *a, real *b, size_t n )
   {
   for ( size_t i = 0; i < n; ++i )
      {
      b[ i ] = a[ i ] - b[ i ];
      }
   }

void vectorAddition( real *a, real *c, real *dest, size_t n )
   {
   for ( size_t i = 0; i < n; ++i )
      {
      dest[ i ] = a[ i ] + c[ i ];
      }
   }

void vectorAddition( real *a, real *b, real *c, real *dest, size_t n )
   {
   for ( size_t i = 0; i < n; ++i )
      {
      dest[ i ] = a[ i ] + b[ i ] + c[ i ];
      }
   }

// result = sum of all elements ^2
double squared_sum( real *vec, size_t n )
   {
   double result = 0.0;
   for ( size_t i = 0; i < n; ++i )
      {
      result = result + (vec[ i ] * vec[ i ]);
      }
   return( result );
   }

