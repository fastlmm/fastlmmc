/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * FastLmmC - {Fast Linear Mixed Model in C++}
 *
 *         File Name:   CLexer.cpp [Tab Separated Values Lexer]
 *           Version:   1.00
 *            Author:   BobD
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

/*
 * Externs we need
 */

/*
 * Globals we define
 */

/*
 * The code
 */
CLexer::CLexer( string filename )
   {
   m_fileName = filename;

   if ( m_fileName.size() == 0 )
      {
      Fatal( "Could not create CLexer.  Parameter 'filename' is zero length string" );
      }

   m_pFile = fopen( m_fileName.c_str(), "rb" );          // read in binary to ensure ftell works right

   if ( !m_pFile )
      {
      char szFullPath[ _MAX_PATH ];
         
      _fullpath( szFullPath, m_fileName.c_str(), sizeof( szFullPath ) );
      Fatal( "Cannot open input file [%s].\n  CRT Error %d: %s", szFullPath, errno, strerror( errno ) );
      }

   m_chLook = 0;
   m_line = 1ul;
   m_column = 1ul;
   AdvanceCh();            // this will read the first character into m_chLook
                           //   and advance cchInLine 1 with a '\0' in m_chLook
   m_cchInLine = 0;        //   so reset m_cchInLine properly

   m_text.clear();

   }

CLexer::~CLexer()
   {
   if ( m_pFile )
      {
      if ( fclose( m_pFile ) != 0 )
         {
         char szFullPath[ _MAX_PATH ];

         _fullpath( szFullPath, m_fileName.c_str(), sizeof( szFullPath ) );
         Fatal( "fclose( %s ) failed.\n  CRT Error %d: %s", szFullPath, errno, strerror( errno ) );
         }
      m_pFile = NULL;
      }

   m_fileName.clear();
   }

#if 0
void CLexer::Init()
   {
   if ( m_fileName.size() == 0 )
      {
      Fatal( "No file to open" );
      }

   m_pFile = fopen( m_fileName.c_str(), "rb" );          // read in binary to ensure ftell works right

   if ( !m_pFile )
      {
      char szFullPath[ _MAX_PATH ];
         
      _fullpath( szFullPath, m_fileName.c_str(), sizeof( szFullPath ) );
      Fatal( "Cannot open input file [%s].\n  CRT Error %d: %s", szFullPath, errno, strerror( errno ) );
      }

   m_chLook = 0;
   m_line = 1;
   m_column = 1;
   AdvanceCh();            // this will read the first character into m_chLook
                           //   and advance cchInLine 1 with a '\0' in m_chLook
   m_cchInLine = 0;        //   so reset m_cchInLine properly
   }
#endif

void CLexer::AdvanceCh()
   {
   /*
    * Advancing the character pointer so keep track of the 'file' location in row column too
    */
   if ( m_chLook == '\n')
      {
      m_column = 1;
      m_cchInLine = 0;
      ++m_line;
      }
   else if ( m_chLook == '\t')
      {
      m_column = ((m_column+8) & (~(8-1)))+1;
      ++m_cchInLine;
      }
   else if ( (m_chLook > ' ') && (m_chLook < '~') )         // printable?
      {
      ++m_cchInLine;
      ++m_column;
      }
   else 
      {
      ++m_cchInLine;
      // don't know what to do with column (eg. \a or 'bell')
      }

SkipCh:
   m_chLook = fgetc( m_pFile );                             // advance to next charater?
   if (m_chLook == '\r')
      {
      m_column = 1;
      m_cchInLine = 0;
      goto SkipCh;                                          // swallow <cr> and let them become <nl> or '\n'
      }
   }

