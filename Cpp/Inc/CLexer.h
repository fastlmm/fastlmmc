#if !defined( CLexer_h )
#define CLexer_h
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
 * CLexer - {Simple Scanner/Lexer Functionality}
 *
 *         File Name:   CLexer.h [Simple Lexer Class]
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file declares a routines for a basic lexer/scanner
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Get our includes
 */
#if defined( NO_GZIP_SUPPORT )
#define  CLEX_FILE   FILE*
#define  CLEX_OPEN   fopen
#define  CLEX_CLOSE  fclose
#define  CLEX_GETC   fgetc
#else
#if defined( COMPILE_FOR_WINDOWS )
typedef struct gzFile_s *gzFile;    /* zlib v1.2.7 on Windows     */
#else
typedef struct gzFile_s *gzFile;    /* zlib v1.2.? on Ubuntu    */
#endif
#define  CLEX_FILE   gzFile
#define  CLEX_OPEN   gzopen
#define  CLEX_CLOSE  gzclose
#define  CLEX_GETC   gzgetc
#endif

/*
 * 'Publish' our defines
 */
typedef enum _TokType
   {
   tokError   = 0
   ,tokTab    = '\t'
   ,tokEOF    = 256
   ,tokEOL
   ,tokSymbol
   ,tokSeparator
   } TokType;

class CToken
   {
public:
   CToken() { type=tokError;offset=0;line=column=0; };

   TokType     type;       // stores Token Type
   long long   offset;     // track file position of token
   size_t      line;
   size_t      column;
   string      text;
   };

class CLexer
   {
public:
   CLexer( const string& filename );
   ~CLexer();

   void        AdvanceCh();
   void        SkipToEol();
   void        SkipWhiteSpace();
   int         LookAheadChar()   { return( m_chLook ); }
   size_t      CurColumn()       { return( m_column ); }
   size_t      CurLine()         { return( m_line ); }
   string&     Filename()        { return( m_fileName ); }
   long long   LexerFileOffset() {  return( m_iFile ); }

private:
   string   m_fileName;
   CLEX_FILE   m_pFile;

   // track file position
   long long     m_iFile;        // current file offset pointer
   unsigned long m_line;
   unsigned long m_column;       // expanded for tabs
   int           m_cchInLine;    // count of characters in the line (

   int           m_chLook;       // look ahead
   };

#endif      // CLexer
