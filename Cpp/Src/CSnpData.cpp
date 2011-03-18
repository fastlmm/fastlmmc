/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * CSnpData - {SNP Data Class}
 *
 *         File Name:   CSnpData.cpp
 *           Version:   1.00
 *            Author:   BobD
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements the CSnpData class for FastLmmC
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
union _mynan {
   unsigned char chNaN[8];
   double NaN;
   } myNaN = { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f };
/*
 * The code
 */
CSnpData::CSnpData()
   {
   snpArray = NULL;
   }

CSnpData::~CSnpData()
   {
   if (snpArray )
      {
      delete [] snpArray;
      }
   }

void CSnpData::ReadSnpFile( string& snpFile )
   {
   /*
    * The snp file is a tab separated file with the form
    *   var NameXxx NameXxy NameXxz ... <eol>
    *   NameYyy x x ...<eol>
    *   NameYyz x x ...<eol>
    *   ... <eol>
    *   <eof>
    *
    * The first row labels the columns
    * each subsequent row contains the data as a Name followed by a single character value for each column
    */
   cRows = 0;
   cColumns = 0;

   CTsvLexer lex( snpFile );
   CToken tok;

   if ( (lex.NextToken( tok ) != tokSymbol) || (tok.text != "var") )
      {
      Fatal( "Expected the file [%s] to start with 'var' but found '%s'", snpFile.c_str(), tok.text.c_str() );
      }

   lex.NextToken( tok );         // consume the 'var' token

   // Get all the column headers
   while( (tok.type != tokEOF) && (tok.type != tokEOL) )
      {
      if ( tok.type != tokTab )
         {
         Fatal( "Expected <tab> to precede column header %d at Line: %d Column: %d", columnLabels.size(), tok.line, tok.column );
         break;
         }
      if ( lex.NextToken(tok) != tokSymbol )
         {
         Fatal( "Expected column header at Line: %d, Column: %d\n  Found %s", tok.line, tok.column, tok.text.c_str() );
         break;
         }
      else
         {
         if (columnMap.count( tok.text ) == 1)
            {
            Fatal( "Duplicate column header [%s] found at Line: %d, Column: %d", tok.text.c_str(), tok.line, tok.column );
            }

         columnLabels.push_back( tok.text );
         columnMap[ tok.text ] = cColumns;
         ++cColumns;
         lex.NextToken( tok );                  // consume the column header
         }
      }

   // Validate the column headers are OK as best we can
   if ( columnLabels.size() != cColumns )
      {
      Fatal( "Program Error.  columnLabels.size()[%d] != cColumns[%d]", columnLabels.size(), cColumns );
      }
   if ( columnLabels.size() < 1 )
      {
      Fatal( "SNP Data must have 1 or more columns of data." );
      }

   Verbose( "Found %d columns of SNP data", columnLabels.size() );

   lex.NextToken( tok );         // consume the <eol> and get the 'row' data
   while( tok.type != tokEOF )   // for all the rows
      {
      if ( tok.type == tokEOL )  // a blank line is acceptable at the end of the file
         {                       //   consume <eol> until <eof>
         while( lex.NextToken( tok ) != tokEOF )
            {
            if ( tok.type != tokEOL )
               {
               Fatal( "Expecting <eof> and found [%s] at Line: %d, Column: %d", tok.text.c_str(), tok.line, tok.column );
               }
            }
         break;                  // was clean <eol>'s to end of file.  Ignore them and continue
         }
      if ( tok.type != tokSymbol )     // expecting a symbol for RowID
         {
         Fatal( "Expected RowID on Line: %d", tok.line );
         }
      if ( rowMap.count(tok.text ) != 0 )
         {
         Fatal( "Duplicate row ID [%s] found at Line: %d, Column: %d", tok.text.c_str(), tok.line, tok.column );
         }
      rowLabels.push_back( tok.text );
      rowMap[tok.text] = cRows;
      ++cRows;
      lex.NextToken( tok );         // consume the ID token

      if ( tok.type != tokTab )     // accumulate the 
         {
         Fatal( "Expected <tab> after row label on RowID: %s, Line: %d", rowLabels[rowLabels.size()-1].c_str(), tok.line );
         }
      lex.NextToken( tok );         // consume the <tab> token

      if ( tok.type != tokSymbol )
         {
         Fatal( "Expected SNP data for RowID: %s at Line: %d, Column: %d.  Found '%s'.", rowLabels[rowLabels.size()-1].c_str(), tok.line, tok.column, tok.text.c_str() );
         }
      size_t cchTok = tok.text.size();
      if ( cchTok != columnLabels.size() )
         {
         Fatal( "Mismatch in SNP data size on RowID: %s, Line: %d\n  Expected %d elements and found %d elements.", rowLabels[rowLabels.size()-1].c_str(), tok.line, columnLabels.size(), cchTok );
         }
      for ( size_t i=0; i<cchTok; ++i )
         {
         size_t ch = tok.text[i];
         if ( (ch != '0') && (ch != '1') && (ch != '2') && (ch != '?') )
            {
            Fatal( "Invalid SNP data at RowID: %s, Line: Column: %d.\n  Found unrecognized character: [%c]", rowLabels[rowLabels.size()-1].c_str(), tok.line, tok.column, ch );
            }
         }
      snpVectors.push_back( tok.text );
      lex.NextToken( tok );         // consume the SNP token
      if ( (tok.type != tokEOL) && (tok.type != tokEOF) )
         {
         Fatal( "Expected an EOL or EOF at RowID: %s, Line: %d, Column: %d.", rowLabels[rowLabels.size()-1].c_str(), tok.line, tok.column );
         }
      lex.NextToken( tok );   // consume <eol>
      }
      
   if ( rowLabels.size() != cRows )
      {
      Fatal( "Program Error.  columnLabels.size()[%d] != cColumns[%d]", columnLabels.size(), cColumns );
      }
   if ( rowLabels.size() < 1 )
      {
      Fatal( "SNP Data must have 1 or more rows of data." );
      }

   Verbose( "Found %d rows of SNP data", rowLabels.size() );
   FatalIfPreviousErrors();

   if ( fCreateReadValidationFiles )         // test global to write validation file 
      {
      string s = snpFile + ".validate.txt";
      WriteSnpFile( s );
      }
   }


void CSnpData::WriteSnpFile( string& snpFile )
   {
   /*
    * The snp file we write is a tab separated file with the form
    *   var NameXxx NameXxy NameXxz ... <eol>
    *   NameYyy x x ...<eol>
    *   NameYyz x x ...<eol>
    *   ... <eol>
    *   <eof>
    *
    * The first row labels the columns
    * each subsequent row contains the data as an Attribute Name followed by a double value separated by <tab> for each column
    */
   FILE *pFile;

   if ( snpFile.size() == 0 )
      {
      Fatal( "No output filename to open " );
      }

   pFile = fopen( snpFile.c_str(), "wt" );          // write ascii

   if ( !pFile )
      {
      char szFullPath[ _MAX_PATH ];
         
      _fullpath( szFullPath, snpFile.c_str(), sizeof( szFullPath ) );
      Fatal( "Cannot open output file [%s].  \n  CRT Error %d: %s", szFullPath, errno, strerror( errno ) );
      }

   // Write the first row (var\t<columnLabels><eol>)
   fprintf( pFile, "var" );
   for( size_t i=0; i<cColumns; ++i )
      {
      fprintf( pFile, "\t%s", columnLabels[i].c_str() );
      }
   fprintf( pFile, "\n" );

   // write the rest of the rows
   for( size_t i=0; i<cRows; ++i )
      {
      fprintf( pFile, "%s\t%s\n", rowLabels[i].c_str(), snpVectors[i].c_str() );
      }
   fclose( pFile );
   }

void CSnpData::CreateColumnMajorArrayFromVectors()
   {
   // create the appropriate array / memory structure for the snp data.
   //   this creates the array in FORTRAN's and MatLab's column major order
   size_t cArrayElements = cColumns * cRows;
   snpArray = new double[cArrayElements];
   for( size_t i=0; i<cRows; ++i )
      {
      for ( size_t j=0; j<cColumns; ++j )
         {
         snpArray[ (i*cColumns)+j ] = DoubleFromChar( snpVectors[i][j] );
         }
      }
   }

void CSnpData::CreateColumnMajorArrayFromVectorsAndMapping( vector<string> columnLabelOrder )
   {
   // Create a 'Fortran' style column major array but use the columnLabelOrder to order the columns
   if ( cColumns != columnLabelOrder.size() )
      {
      Fatal( "Cannot create SnpData column major array from Vector and Mapping.  Vectors are different lengths.  "
             "SnpColumns[%d] and mapping[%d]", cColumns, columnLabelOrder.size() );
      }
   vector<size_t> rgMap( cColumns );
   for( size_t iCol=0; iCol<cColumns; ++iCol )
      {
      if ( columnMap.count( columnLabelOrder[iCol] ) != 1 )
         {
         Fatal( "Cannot create column major array from Vector and Mapping.  SnpColumns does not contain mapping key [%s]", columnLabelOrder[iCol].c_str() );
         }
      rgMap[iCol] = columnMap[ columnLabelOrder[iCol] ];
      }

   size_t csnp = cColumns * cRows;
   snpArray = new double[ csnp ];
   for ( size_t iRow=0; iRow<cRows; ++iRow )
      {
      for ( size_t iCol=0; iCol<cColumns; ++iCol )
         {
         if ( rgMap[iCol] != iCol )
            {
            Verbose( "\n -- assign snparray[%d]= snpVector[%d][%d]", (iRow*cColumns)+iCol, iRow, rgMap[iCol] );
            }
         char ch = snpVectors[iRow][rgMap[iCol]];
         snpArray[ (iRow*cColumns) + iCol ] = DoubleFromChar( ch );
         }
      }
   }

double CSnpData::DoubleFromChar( char ch )
   {
   double dbl;
   if (ch == '0' )
      {
      dbl = 0.0;
      }
   else if (ch == '1' )
      {
      dbl = 1.0;
      }
   else if (ch == '2' )
      {
      dbl = 2.0;
      }
   else // ( ch == '?' )
      {
      dbl = myNaN.NaN;
      }
   return( dbl );
   }

