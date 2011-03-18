#include "FastLmmC.h"

const string MatrixIO::OS_directoryDelim="\\";     //Windows   //TODO: Make this selection a macro
//const string MatrixIO::OS_directoryDelim="/";    //UNIX

int MatrixIO::writeMatrix2dBIN(string filename, string matrixID, int n, int m, double *A){
   
   string binfilename=filename;
   binfilename.append(".bin");
   string headerfilename=filename;
   headerfilename.append(".txt");
   //write the binary data:
   ofstream binfile;
   binfile.open(binfilename.c_str(), ios_base::binary);
   binfile.write(reinterpret_cast<char*>(A), n*m*sizeof(double));
   binfile.close();

   //write the header:
   ofstream headerfile;
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

int MatrixIO::querySizeMatrixBIN(string filename, string *matrixID, int *n, int *m){
   string headerfilename=filename;
   headerfilename.append(".txt");
   ifstream headerfile;
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
   return 1;
}

int MatrixIO::readBinary2dBIN(string filename, string matrixID, int n, int m, double *A){
   string headerfilename=filename;
   headerfilename.append(".txt");
   ifstream headerfile;
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
   ifstream binfile;
   binfile.open(binfilename.c_str(),ios::binary);
   //for(int i=0;i<(n) * (m);++i){printf("%f ", A[i]);}printf("\n");
   binfile.read(reinterpret_cast<char*>(A),n * m * sizeof(double));
   //for(int i=0;i<(n) * (m);++i){printf("%f ", A[i]);}  printf("\n");
   binfile.close();
   return 1;
}

vector<string> MatrixIO::readHeaderFile(string filename){
   vector<string> header(0);
   ifstream headerfile;
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
int MatrixIO::writeHeaderFile(string filename, vector<string> header){
   //write the header:
   ofstream headerfile;
   headerfile.open(filename.c_str());
   for (size_t i=0;i<header.size();++i){
      headerfile << header[i];
      headerfile << "\t";
   }
   headerfile.close();
   return 0;
}

int MatrixIO::createDirectory(string dirname){
   if ( _access( dirname.c_str(), 0 ) == 0 )    {
      struct stat status;
      stat( dirname.c_str(), &status );
      if ( status.st_mode & S_IFDIR )        {
         //cout << "The directory exists." << endl;        
      }        
      else        {
         printf("The path you entered is a file.\n");
         throw(1);
      }    
   }    
   else    {
      string mkdirectory="mkdir ";
      mkdirectory.append(dirname);
      system (mkdirectory.c_str());
   }
   return 0;
}
