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
 * CLexer - {Simple Lexer functionality}
 *
 *         File Name:   CLexer.cpp [Simple Lexer Class]
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements the basic lexer/scanner for FastLmmC
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Include Files
 */
#include "FastLmmC.h"
#if !defined( NO_GZIP_SUPPORT )
#include "zlib.h"

#if defined ( COMPILE_FOR_WINDOWS )    // If we are using Visual C compiler
#if defined( USE_ACML )
#pragma comment( lib, "zlibwapi" )     // tell linker to load the appropriate library (dll)
#else
#pragma comment( lib, "zlibstat" )     // tell linker to load the appropriate library (static)
#endif
#endif
#endif

/*
 * Externs we need
 */

/*
 * Globals we define
 */

/*
 * The code
 */
CLexer::CLexer( const string& filename )
   {
   m_fileName = filename;

   if ( m_fileName.size() == 0 )
      {
      Fatal( "Could not create CLexer.  Parameter 'filename' is zero length string" );
      }

   m_pFile = CLEX_OPEN( m_fileName.c_str(), "rb" );         // read in binary to ensure ftell works right, and gzip requires it
   if ( !m_pFile )
      {
      string gzFilename = filename + ".gz";
      m_pFile = CLEX_OPEN( gzFilename.c_str(), "rb" );      // look for gzipped file
      if ( !m_pFile )
         {
         string fullPath = FullPath( m_fileName );
         Fatal( "Cannot open input file [%s].\n  CRT Error %d: %s", fullPath.c_str(), errno, strerror( errno ) );
         }
      Verbose( verboseMinor, "            Loading compressed file: [%s]", gzFilename.c_str() );
      }

   m_chLook = 0;
   m_line = 1ul;
   m_column = 1ul;
   AdvanceCh();            // this will read the first character into m_chLook
                           //   and advance cchInLine 1 with a '\0' in m_chLook
   m_cchInLine = 0;        //   so reset m_cchInLine properly
   m_iFile = 0;            // and m_iFile to reflect the look ahead char is not 'read' yet
   }

CLexer::~CLexer()
   {
   if ( m_pFile )
      {
      int rc = CLEX_CLOSE( m_pFile );
      if ( rc != 0 )
         {
         string fullPath = FullPath( m_fileName );
         Fatal( "fclose( %s ) failed.\n  CRT Error %d: %s", fullPath.c_str(), errno, strerror( errno ) );
         }
      m_pFile = nullptr;
      }

   m_fileName.clear();
   }

void CLexer::AdvanceCh()
   {
   /*
    * Advancing the character pointer so keep track of the 'file' location in row column too
    */
   if ( m_chLook == '\n')
      {
      ++m_iFile;
      m_column = 1;
      m_cchInLine = 0;
      ++m_line;
      }
   else if ( m_chLook == '\t')
      {
      ++m_iFile;
      m_column = ((m_column+8) & (~(8-1)))+1;
      ++m_cchInLine;
      }
   else if ( (m_chLook >= ' ') && (m_chLook <= '~') )         // printable?
      {
      ++m_iFile;
      ++m_cchInLine;
      ++m_column;
      }
   else 
      {
      ++m_iFile;
      ++m_cchInLine;
      // don't know what to do with column (eg. \a or 'bell')
      }

SkipCh:
   m_chLook = CLEX_GETC( m_pFile );                         // advance to next charater?
   if (m_chLook == '\r')
      {
      ++m_iFile;
      m_column = 1;
      m_cchInLine = 0;
      goto SkipCh;                                          // swallow <cr> and let them become <nl> or '\n'
      }

   }

void CLexer::SkipToEol()
   {
   do
      {
      AdvanceCh();
      } while( (LookAheadChar() != EOF) 
            && (LookAheadChar() != '\n') );
   }

void CLexer::SkipWhiteSpace()
   {
   do
      {
      AdvanceCh();
      } while( (LookAheadChar() == ' ') 
            || (LookAheadChar() == '\t') );
   }
