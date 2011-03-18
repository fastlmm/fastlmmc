#if !defined( TsvLexer_h )
#define TsvLexer_h
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

/*
 * 'Publish' our class declarations / function prototypes
 */
class CTsvLexer : CLexer
   {
public:
   CTsvLexer( string filename );
   ~CTsvLexer();
   unsigned NextToken( CToken& tok );
   };

/*
 * 'Publish' the globals we define
 */

#endif   // TsvLexer_h
