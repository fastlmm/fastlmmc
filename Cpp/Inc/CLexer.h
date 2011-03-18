#if !defined( CLexer_h )
#define CLexer_h
/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * TsvLexer - {Tab Separated Values Scanner/Lexer Routines}
 *
 *         File Name:   TsvLexer.h
 *           Version:   1.00
 *            Author:   BobD
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file declares the TsvScanner routines for FastLmmC
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Get our includes
 */

/*
 * 'Publish' our defines
 */
#define  tokEOF      (256)
#define  tokEOL      (257)
#define  tokSymbol   (258)

#define  tokError    (-1)
#define  tokTab      ('\t')

class CToken
   {
public:
   CToken() { type=0;offset=line=column=0; };

   unsigned type;       // stores tokType
   size_t   offset;     // track file position of token
   size_t   line;
   size_t   column;
   string   text;
   };

class CLexer
   {
public:
   CLexer( string filename );
   ~CLexer();
//   unsigned NextToken( CToken& tok );

   void        AdvanceCh();
   int         LookAheadChar() { return( m_chLook ); }
   size_t      CurColumn() { return( m_column ); }
   size_t      CurLine()   { return( m_line ); }
   string&     Filename()  { return( m_fileName ); }
   long long   LexerFileOffset() {  return( _ftelli64( m_pFile )-1 ); }   // already read the first character

private:
   string   m_fileName;
   FILE     *m_pFile;

   // track file position
   unsigned long m_line;
   unsigned long m_column;       // expanded for tabs
   int      m_cchInLine;         // count of characters in the line (

   int      m_chLook;            // look ahead
   string   m_text;
   };

#endif      // CLexer
