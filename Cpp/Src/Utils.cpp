/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * Utils - {Utility Routines}
 *
 *         File Name:   Utils.cpp
 *           Version:   1.00
 *            Author:   BobD
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
int fVerbose = 0;             // global flag for verbose output

static int cErrorsMax = 25;
static int cErrors = 0;
static int cWarns = 0;

/*
 * The code
 */
void Fatal( char *szFmt, ... )
   {
   va_list va;

   fprintf( stderr, "\nFatal Error : " );

   va_start( va, szFmt );
   vfprintf( stderr, szFmt, va );
   va_end( va );
   exit( 1 );
   }


void Error( char *szFmt, ... )
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
   

void Warn( char *szFmt, ... )
   {
   va_list va;

   fprintf( stderr, "\nWarning : " );

   va_start( va, szFmt );
   vfprintf( stderr, szFmt, va );
   va_end( va );
   }


void Verbose( char *szFmt, ... )
   {
   if (!fVerbose )
      {
      return;
      }

   va_list va;

   fprintf( stderr, "\n" );

   va_start( va, szFmt );
   vfprintf( stderr, szFmt, va );
   va_end( va );
   }

void Progress( char *szFmt, ... )
   {
   va_list va;

   va_start( va, szFmt );
   vfprintf( stderr, szFmt, va );
   va_end( va );
   }

#if defined ( _MSC_VER )       // using Visual C compiler
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
   fprintf( stderr, "\nNTError : 0x%08X : %s", err, lpMsgBuf );
   LocalFree( lpMsgBuf );
   exit( 1 );
   }

#else    //#if !defined( _MSC_VER )

void NTLastError( void )
   {
   Fatal( "NTLastError stub for fatal error" );
   }

char *_fullpath( char *absPath, const char *relPath, size_t maxLength )
   {
   // stub a _fulllpath() routine for systems that don't have it
   strcpy( absPath, relPath );
   return( absPath );
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

#endif
