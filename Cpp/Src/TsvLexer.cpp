/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * FastLmmC - {Fast Linear Mixed Model in C++}
 *
 *         File Name:   TsvLexer.cpp [Tab Separated Values Lexer]
 *           Version:   1.00
 *            Author:   BobD
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements the TSV scanner for FastLmmC
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
CTsvLexer::CTsvLexer( string filename ) : CLexer(filename)
   {
   }

CTsvLexer::~CTsvLexer()
   {
   // CLexer does all the cleanup right now
   }


unsigned CTsvLexer::NextToken( CToken& tok )
   {
   tok.type = tokError;
   tok.column = CurColumn();     // get the current position so we start collecting the token --  m_column;
   tok.line = CurLine(); 
   tok.offset = LexerFileOffset();
   tok.text.clear();

   switch( LookAheadChar() )
      {
   case EOF:      // return EOF token
      tok.type = tokEOF;
      break;
   case '\t':     // return a TAB token
      tok.type = tokTab;
      AdvanceCh();
      break;
   case '\n':     // return a EOL token
      tok.type = tokEOL;
      AdvanceCh();
      break;
   default:       // see if it is a Symbol (ascii printable)
      if ( (LookAheadChar() < ' ') || (LookAheadChar() > '~') )
         {
         Fatal( "Found invalid character in file [%s] at Line: %d Column: %d", Filename().c_str(), tok.line, tok.column );
         break;
         }
      
      // accumulate ascii printable characters as a symbol
      while ( (LookAheadChar() >= ' ') && (LookAheadChar() <= '~' ) )
         {
         tok.text += (char)LookAheadChar();
         AdvanceCh();
         }
      tok.type = tokSymbol;
      break;
      }
   return( tok.type );
   }

