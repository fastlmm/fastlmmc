#if !defined( Utils_h )
#define Utils_h
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
 * 'Publish' our class declarations / function prototypes
 */
void Error( char *szFmt, ... );
void Fatal( char *szFmt, ... );
void FatalIfPreviousErrors();
void NTLastError( void );
void Progress( char *szFmt, ... );
void Verbose( char *szFmt, ... );
void Warn( char *szFmt, ... );

/*
 * 'Publish' the globals we define
 */
extern int fVerbose;

/*
 *  Define the stubs, translation routines, and 'stuff'
 *   that Windows normally provides so we can support 
 *   a bit of cross platform operation
 */
#if !defined( _MSC_VER )
inline double abs( double x ) { return( fabs(x) ); }
inline int min( int m1, int m2 ) { return( ( m1 < m2 ) ? m1 : m2 ); }
inline double min( double m1, double m2 ) { return( (m1 < m2 ) ? m1 : m2 ); }
inline int max( int m1, int m2 ) { return( ( m1 > m2 ) ? m1 : m2 ); }
inline int _strcmpi( const char *s1, const char *s2 ) { return( strcasecmp( s1, s2 ) ); }
inline int _access( const char *path, int mode ) { return( access( path, mode ) ); }

int  GetComputerName( char *pchName, unsigned long* pcchMax );
#define  MAX_COMPUTERNAME_LENGTH 15

char *_fullpath( char *absPath, const char *relPath, size_t maxLength );
#define  _MAX_PATH   260

unsigned int _set_output_format( unsigned int format );
#define  _TWO_DIGIT_EXPONENT  1
#endif   // ! _MSC_VER

#endif   // Utils_h
