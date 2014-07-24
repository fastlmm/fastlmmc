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

#include "FastLmmC.h"

//const string MatrixIO::OS_directoryDelim="\\";     //Windows   //TODO: Make this selection a macro
const string MatrixIO::OS_directoryDelim="/";    //UNIX

int MatrixIO::writeMatrix2dBIN(const string& filename, const string& matrixID, size_t n, size_t m, real *A){
   
   string binfilename=filename;
   binfilename.append(".bin");
   string headerfilename=filename;
   headerfilename.append(".txt");
   //write the binary data:
   std::ofstream binfile;
   binfile.open(binfilename.c_str(), std::ios_base::binary);
   binfile.write(reinterpret_cast<char*>(A), n*m*sizeof(real));
   binfile.close();

   //write the header:
   std::ofstream headerfile;
   headerfile.open(headerfilename.c_str());
   headerfile << "ID:\t";
   headerfile << matrixID;
   headerfile << "\nn:\t";
   headerfile << n;
   headerfile << "\nm:\t";
   headerfile << m;
   headerfile.close();
   return 0;
}

int MatrixIO::querySizeMatrixBIN(const string& filename, string *matrixID, size_t *n, size_t *m){
   string headerfilename=filename;
   headerfilename.append(".txt");
   //ifstream headerfile;
   //headerfile.open(headerfilename);

   CTsvLexer lex( headerfilename );
   CToken tok;

   if ( (lex.NextToken( tok ) != tokSymbol) || (tok.text != "ID:") )
      {
      Fatal( "Expected the file [%s] to start with 'ID:' but found '%s'", headerfilename.c_str(), tok.text.c_str() );
      }
       
   lex.NextToken( tok );                  // consume the 'ID:' token
   if ( lex.NextToken( tok ) != tokSymbol )
      {
         Fatal( "Expected a symbol" );
      }
   *matrixID=tok.text;
   lex.NextToken( tok );                  // consume the matrixID
   if ( (lex.NextToken( tok ) != tokSymbol) || (tok.text != "n:") )
      {
         Fatal( "Expected the second row in file [%s] to start with 'n:' but found '%s'", headerfilename.c_str(), tok.text.c_str() );
      }
   lex.NextToken( tok );                  // consume 'n:'
   if ( lex.NextToken( tok ) != tokSymbol )
      {
         Fatal( "Expected a symbol" );
      }
   *n=atoi(tok.text.c_str()); 
   lex.NextToken( tok );                  // consume n
   if ( (lex.NextToken( tok ) != tokSymbol) || (tok.text != "m:") )
      {
         Fatal( "Expected the second row in file [%s] to start with 'm:' but found '%s'", headerfilename.c_str(), tok.text.c_str() );
      }
   lex.NextToken( tok );                  // consume 'm:'
   if ( lex.NextToken( tok ) != tokSymbol )
      {
         Fatal( "Expected a symbol" );
      }
   *m=atoi(tok.text.c_str());
   return 0;
}

int MatrixIO::readBinary2dBIN(const string& filename, const string& matrixID, size_t n, size_t m, real *A){
   string headerfilename=filename;
   headerfilename.append(".txt");
   std::ifstream headerfile;
   //headerfile.open(headerfilename);

   CTsvLexer lex( headerfilename );
   CToken tok;

   if ( (lex.NextToken( tok ) != tokSymbol) || (tok.text != "ID:") )
      {
      Fatal( "Expected the file [%s] to start with 'ID:' but found '%s'", headerfilename.c_str(), tok.text.c_str() );
      }
       
   lex.NextToken( tok );                  // consume the 'ID:' token
   if ( (lex.NextToken( tok ) != tokSymbol) || (tok.text != matrixID) )
      {
         Fatal( "Expected the matrixID in file [%s] to be '%s' but found '%s'", headerfilename.c_str(),matrixID.c_str(), tok.text.c_str() );
      }

   //printf("%s\n", tok.text.c_str());
   lex.NextToken( tok );                  // consume the matrixID
   if ( (lex.NextToken( tok ) != tokSymbol) || (tok.text != "n:") )
      {
         Fatal( "Expected the second row in file [%s] to start with 'n:' but found '%s'", headerfilename.c_str(), tok.text.c_str() );
      }
   lex.NextToken( tok );                  // consume 'n:'
   if ( (lex.NextToken( tok ) != tokSymbol) || (atoi(tok.text.c_str()) != n) )
      {
         Fatal( "Expected n in file [%s] to be '%i' but found '%s'", headerfilename.c_str(), n, tok.text.c_str() );
      }
   lex.NextToken( tok );                  // consume n
   if ( (lex.NextToken( tok ) != tokSymbol) || (tok.text != "m:") )
      {
         Fatal( "Expected the second row in file [%s] to start with 'm:' but found '%s'", headerfilename.c_str(), tok.text.c_str() );
      }
   lex.NextToken( tok );                  // consume 'm:'
   if ( (lex.NextToken( tok ) != tokSymbol) || (atoi(tok.text.c_str()) != m) )
      {
         Fatal( "Expected m in file [%s] to be '%i' but found '%s'", headerfilename.c_str(), m, tok.text.c_str() );
      }
   
   string binfilename=filename;
   binfilename.append(".bin");
   std::ifstream binfile;
   binfile.open(binfilename.c_str(),std::ios::binary);
   //for(int i=0;i<(n) * (m);++i){printf("%f ", A[i]);}printf("\n");
   binfile.read(reinterpret_cast<char*>(A),n * m * sizeof(real));
   //for(int i=0;i<(n) * (m);++i){printf("%f ", A[i]);}  printf("\n");
   binfile.close();
   return 1;
}

vector<string> MatrixIO::readHeaderFile(const string& filename){
   vector<string> header(0);
   std::ifstream headerfile;
   //headerfile.open(filename);
   
   CTsvLexer lex( filename );
   CToken tok;
   //int i=0;
   lex.NextToken(tok);
   while(tok.type != tokEOF){
      
      //if (tok.type!=tokEOF){
         header.push_back(tok.text);
         lex.NextToken(tok);
         lex.NextToken(tok);
      //}
      
      //printf(header[i].c_str());printf("\n");
      //++i;
      
      //printf(tok.text.c_str());
      //lex.NextToken(tok);
      //printf(tok.text.c_str());
      //if (lex.NextToken( tok ) != tokSymbol) 
      //   {
      //   Fatal( "Expected the file [%s] to have a newline in the current position  but found '%s'", tok.text.c_str() );
      //   }
   }
   return header;
}
int MatrixIO::writeHeaderFile(const string& filename, vector<string>& header){
   //write the header:
   std::ofstream headerfile;
   headerfile.open(filename.c_str());
   for (size_t i=0;i<header.size();++i){
      headerfile << header[i];
      headerfile << "\t";
   }
   headerfile.close();
   return 0;
}

int MatrixIO::createDirectory(const string& dirname)
   {
   if ( _access( dirname.c_str(), 0 ) == 0 )
      {
      struct stat status;
      stat( dirname.c_str(), &status );
      if ( status.st_mode & S_IFDIR )
         {
         //cout << "The directory exists." << endl;
         }        
      else
         {
         Fatal( "Cannot create directory [%s].  A file by that name already exists", dirname.c_str() );
         }    
      }    
   else
      {
      string mkdirectory = "mkdir ";
      mkdirectory.append( dirname );
      int i = system( mkdirectory.c_str() );
      if ( i == -1 )
         {
         Fatal( "system( %s ) returned failure.", mkdirectory.c_str() );
         }
      }

   return( 0 );
   }

int MatrixIO::writeMatrix2dTXT(const string& filename, const string& matrixID, const string& separator, vector<string>& rowhdr, vector<string>& colhdr, size_t rows, size_t cols, real *A)
   {
   /*
    * The kernel file we write is a tab separated file with the form
    *   var NameXxx NameXxy NameXxz ... <eol>
    *   NameYyy x x ...<eol>
    *   NameYyz x x ...<eol>
    *   ... <eol>
    *   <eof>
    *
    * The first row labels the columns
    * each subsequent row contains the data as an Attribute Name followed by a real value separated by <tab> for each column
    */
   FILE *pFile;

   if ( filename.size() == 0 )
      {
      Fatal( "No output filename to open " );
      }

   pFile = fopen( filename.c_str(), "wt" );          // write ascii

   if ( !pFile )
      {
      string fullPath = FullPath( filename );
      Fatal( "Cannot open output file [%s].  \n  CRT Error %d: %s", fullPath.c_str(), errno, strerror( errno ) );
      }

   if( rowhdr.size()!=rows || colhdr.size() != cols) Fatal("\nDimensions do not match up for matrix to be written. rowhdr.size() = %i, rows = %i, colhdr.size() = %i, cols = %i\n");
   
   unsigned int old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);   // to be compatible, use 2 digit exponents

   // Write the first row (var\t<columnLabels><eol>)
   fprintf( pFile, "%s", matrixID.c_str() );
   for( size_t i=0; i<cols; ++i )
      {
      fprintf( pFile, "%s%s", separator.c_str(),colhdr[i].c_str() );
      }
   fprintf( pFile, "\n" );

   // write the rest of the rows
   for( size_t i=0; i<rows; ++i )
      {
      fprintf( pFile, "%s", rowhdr[i].c_str() );
      for ( size_t j=0; j<cols; ++j )
         {
         fprintf( pFile, "%s%.15lG", separator.c_str(), A[j * rows + i] );    // get all the precision using g format
         }
      fprintf( pFile, "\n" );
      }

   _set_output_format( old_exponent_format );                  // restore the default format behavior
   fclose( pFile );
   return( 0 );
   }
