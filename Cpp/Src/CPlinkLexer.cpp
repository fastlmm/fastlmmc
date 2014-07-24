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
 * CPlinkLexer - {PLINK File Lexer Class}
 *
 *         File Name:   CPlinkLexer.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements the CPlinkLexer class for FastLmmC
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Include Files
 */
#include "FastLmmC.h"

#if 0
#if !defined( __GNUC__ )         // gcc doesn't yet support <regex>
#include <regex>

static std::tr1::regex rxInt( "[-+]?[0-9]+" );
static std::tr1::regex rxFp( "[-+]?([0-9]*\\.)?[0-9]+([eE][-+]?[0-9]+)?" );
static std::tr1::regex rxChromosome( "[0-9]+" );
#endif
#endif

int maxChromosomeValue = 26;

/*
 * Validate string is a floating point number
 */
bool FValidateFpString( const string& s )
   {
   bool fMatch = false;
   const char *pch = s.c_str();
   int state = 0;

   while( *pch )
      {
      switch( state )
         {
      case 0:  // looking for leadin +/- if there
         if ( ( *pch == '-' ) || ( *pch == '+' ) )
            {
            ++pch;
            }
         state = 1;
         break;
      case 1:  // past sign. Get digit or decimal point
         if ( isdigit( *pch ) )
            {
            ++pch;
            state = 2;  // accumulate digits before decimal or exponent
            fMatch = true;
            }
         else if ( *pch == '.' )
            {
            ++pch;
            state = 3;  // Require digit after decimal point
            }
         else
            {
            return( false );
            }
         break;
      case 2:  // accumulate digits before the decimal point.  Get digit, '.', or 'E'
         if ( isdigit( *pch ) )
            {
            ++pch;      // accumulate digits
            }
         else if ( *pch == '.' )
            {
            ++pch;
            state = 4;  // accumulate digits for fractional part after decimal point
            }
         else if ( (*pch == 'e') || (*pch == 'E') )
            {
            ++pch;
            state = 5;  // scan for exponent
            fMatch = false;
            }
         else
            {
            return( false );
            }
         break;
      case 3:  // we have at no digit before decimal point, get at least one digit
         if ( isdigit( *pch ) )
            {
            ++pch;
            state = 4;  // accumulate digits for fractional part after decimal point
            fMatch = true;
            }
         else
            {
            return( false );
            }
         break;
      case 4:  // accumulate digits for fractional part after decimal point
         if ( isdigit( *pch ) )
            {
            ++pch;
            }
         else if ( (*pch == 'e') || (*pch == 'E') )
            {
            ++pch;
            state = 5;  // scan for exponent
            fMatch = false;
            }
         else
            {
            return( false );
            }
         break;
      case 5:  // scan for exponent.  Look for lead-in +/-
         if ( ( *pch == '-' ) || ( *pch == '+' ) )
            {
            ++pch;
            }
         state = 6;
         break;
      case 6:  // get first exponent digit
         if ( isdigit( *pch ) )
            {
            ++pch;
            state = 7;
            fMatch = true;
            }
         else
            {
            return( false );
            }
         break;
      case 7:  // accumulate exponent digits
         if ( isdigit( *pch ) )
            {
            ++pch;
            }
         else
            {
            return( false );
            }
         break;
      default:
         Fatal( "Programming Error giving invalid state[%d] in FValidateFpString().  Expected 0-7", state );
         break;
         }
      }
   return( fMatch );
   }

/*
 * Validate string is an integer number
 */
bool FValidateIntString( const string& s )
   {
   bool fMatch = false;
   const char *pch = s.c_str();
   int state = 0;

   while( *pch )
      {
      switch( state )
         {
      case 0:  // looking for leadin +/- if there
         if ( ( *pch == '-' ) || ( *pch == '+' ) )
            {
            ++pch;
            }
         state = 1;
         break;
      case 1:  // past sign. Get at least one digit
         if ( isdigit( *pch ) )
            {
            ++pch;
            state = 2;  // accumulate digits
            fMatch = true;
            }
         else
            {
            return( false );
            }
         break;
      case 2:  // accumulate remainder of the digits.
         if ( isdigit( *pch ) )
            {
            ++pch;      // accumulate digits
            }
         else
            {
            return( false );
            }
         break;
      default:
         Fatal( "Programming Error giving invalid state[%d] in FValidateFpString().  Expected 0-2", state );
         break;
         }
      }
   return( fMatch );
   }

/*
 * Validate string is an integer number
 */
bool FValidateChromosomeString( const string& s )
   {
   bool fMatch = false;
   const char *pch = s.c_str();
   int state = 0;

   while( *pch )
      {
      switch( state )
         {
      case 0:  // Get at least one digit
         if ( isdigit( *pch ) )
            {
            ++pch;
            state = 1;  // accumulate digits
            fMatch = true;
            }
         else
            {
            return( false );
            }
         break;
      case 1:  // accumulate remainder of the digits.
         if ( isdigit( *pch ) )
            {
            ++pch;      // accumulate digits
            }
         else
            {
            return( false );
            }
         break;
      default:
         Fatal( "Programming Error giving invalid state[%d] in FValidateFpString().  Expected 0-1", state );
         break;
         }
      }
   return( fMatch );
   }

CPlinkLexer::CPlinkLexer( const string& file ) : CLexer( file )
   {
   // CLexer does all the initialization
   }

CPlinkLexer::~CPlinkLexer()
   {
   // CLexer does all the cleanup right now
   }

unsigned CPlinkLexer::NextToken( CToken& tok )
   {
RetryAfterSkip:

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
   case '\t':     // White space - skip it all and return separator
   case ' ':
      SkipWhiteSpace();
      tok.type = tokSeparator;
      goto RetryAfterSkip;       // skip whitespace
      break;
   case '#':      // comment - skip to end of line
      SkipToEol();
      goto RetryAfterSkip;      // could be EOF
      break;
   case '\n':     // return a EOL token
      tok.type = tokEOL;
      AdvanceCh();
      break;
   default:       // see if it is a Symbol (ascii printable but NOT ' ')
      if ( (LookAheadChar() < ' ') || (LookAheadChar() > '~') )
         {
         Fatal( "Found invalid character in file [%s] near line: %d:%d", Filename().c_str(), tok.line, tok.column );
         break;
         }
      
      // accumulate ascii printable characters as a symbol
      while ( (LookAheadChar() > ' ') && (LookAheadChar() <= '~' ) )
         {
         tok.text += (char)LookAheadChar();
         AdvanceCh();
         }
      tok.type = tokSymbol;
      break;
      }
   return( tok.type );
   }

void CPlinkLexer::ExpectId( CToken& tok, string& output, const char *description )
   {
   if ( tok.type != tokSymbol )     // get ID
      {
      Fatal( "Expected %s in file [%s] near line %d:%d.  Found [%s]", description, Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
      }

   output = tok.text;
   NextToken( tok );                // consumed the token so advance
   }

void CPlinkLexer::ExpectInt( CToken& tok, int& output, const char *description )
   {
   if ( tok.type != tokSymbol )     // get string
      {
      Fatal( "Expected %s in file [%s] near line %d:%d.  Found [%s]", description, Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
      }
//#if !defined( __GNUC__ )         // gcc doesn't yet support <regex>
//   //std::tr1::regex rxInt( "[-+]?[0-9]+" );       // move outside loop for perf
//   if ( !std::tr1::regex_match( tok.text, rxInt ) )
   if ( !FValidateIntString( tok.text ) )
      {
      Fatal( "Expected an integer %s in file [%s] near line %d:%d.  Found [%s]", description, Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
      }
//#endif
   output = atoi( tok.text.c_str() );
   if ( output == std::numeric_limits<int>::max() )
      {
      Fatal( "Integer overflow converting %s in file [%s] near line %d:d. Found [%s]", description, Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
      }
   NextToken( tok );
   }

void CPlinkLexer::ExpectReal( CToken& tok, real& output, const char *description )
   {
   if ( tok.type != tokSymbol )     // get string representation
      {
      Fatal( "Expected %s in file [%s] near line %d:%d.  Found [%s]", description, Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
      }
//#if !defined( __GNUC__ )         // gcc doesn't yet support <regex>
//   //std::tr1::regex rxFp( "[-+]?([0-9]*\\.)?[0-9]+([eE][-+]?[0-9]+)?" );     // move outside loop for perf
//   if ( !std::tr1::regex_match( tok.text, rxFp ) )
   if ( !FIsFpString( tok.text ) )
      {
      Fatal( "Expected a floating point number in file [%s] near line: %d:%d.  Found [%s]", Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
      }
//#endif
   output = (real)atof( tok.text.c_str() );
   NextToken( tok );
   }

void CPlinkLexer::ExpectSex( CToken& tok, int& output )
   {
   if ( tok.type != tokSymbol )     // get the fifth column - sex
      {
      Fatal( "Expected a Sex indicator in file [%s] near line %d:%d.  Found [%s]", Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
      }
   if ( tok.text.length() != 1 )
      {
      Fatal( "Expected Sex indicator in file [%s] near line %d:%d to be 1 character long.  Found %d characters [%s]",
         Filename().c_str(), tok.line, tok.column, tok.text.length(), tok.text.c_str() );
      }
   char ch = tok.text[0];
   if ( ch == '1' ) output = 1;
   else if ( ch == '2' ) output = 2;
   else output = 0;
   NextToken( tok );
   }

void CPlinkLexer::ExpectPhenotype( CToken& tok, real& output, PhenotypeType& outputType )
   {
   if ( tok.type != tokSymbol )     //   this is trickier as it can be 0, 1, 2, -9, or an FP number
      {
      Fatal( "Expected phenotype indicator in file [%s] near line %d:%d.  Found [%s]", Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
      }

   outputType = PHENOTYPE_AFFECTION;            
   if ( _strcmpi( tok.text.c_str(), "0" ) == 0 ) output = 0.0;
   else if ( _strcmpi( tok.text.c_str(), "1" ) == 0 ) output = 1.0;
   else if ( _strcmpi( tok.text.c_str(), missingPhenotypeString.c_str() ) == 0 ) output = std::numeric_limits<real>::quiet_NaN();
   else 
      {
      // must be an FP number
//#if !defined( __GNUC__ )         // gcc doesn't yet support <regex>
//      //std::tr1::regex rxFp( "[-+]?([0-9]*\\.)?[0-9]+([eE][-+]?[0-9]+)?" );     // move outside loop for perf
//      if ( !std::tr1::regex_match( tok.text, rxFp ) )
      if ( !FIsFpString( tok.text ) )
         {
         Fatal( "Expected phenotype indicator in file [%s] near line %d:%d to be 0, 1, -9, or a floating point number.  Found [%s]",
            Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
         }
//#endif
      output = (real)atof( tok.text.c_str() );
      if ( output == missingPhenotypeValue )
         {
         output = std::numeric_limits<real>::quiet_NaN();
         }
      outputType = PHENOTYPE_QUANTITATIVE;
      }
   NextToken( tok );
   }

void CPlinkLexer::ExpectAllele( CToken& tok, char& output, const char *description )
   {
   if ( (tok.type != tokSymbol) || (tok.text.length() != 1) )
      {
      Fatal( "Expected %s allele of genotype pair in file [%s] near line %d:%d.  Found %s", description, Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
      }

   output = tok.text[0];
   NextToken( tok );                // consumed the token so advance
   }

void CPlinkLexer::ExpectSnpNucleotides( CToken& tok, char& firstNucleotide, char& secondNucleotide )
   {
   size_t tokenLine = tok.line;
   size_t tokenColumn = tok.column;
   ExpectAllele( tok, firstNucleotide, "First" );
   ExpectAllele( tok, secondNucleotide, "Second" );
   if ( ((firstNucleotide == '0') && (secondNucleotide != '0'))
      ||((firstNucleotide != '0') && (secondNucleotide == '0'))
      )
      {
      // bad genotype data.  Mixed a missing and a non-missing character
      // TODO:  PLink documentation calls out that mixing missing and
      //         not missing characters for the genotype is illegal
      //         and the user should resolve the data to provide a 
      //         homozygous genotype or set both alleles to missing.
      //
      //         I tend to prefer to not parse 'illegal' input as it 
      //         causes 'undefined' behavior to occur.  
      //         In most cases, the user should fix the input.
      //
      //   Fatal
      Warn( "Inconsistent missing genotype information found in file [%s] near line %d:%d."
            "\n  found [%c] [%c]", Filename().c_str(), tokenLine, tokenColumn, firstNucleotide, secondNucleotide );
      Warn( "    Genotype information treated as missing" );
      firstNucleotide = '0';
      secondNucleotide = '0';
      }
   return;
   }

void CPlinkLexer::ExpectSnpAlleles( CToken& tok, char& majorAllele, char& minorAllele )
   {
   size_t tokenLine = tok.line;
   size_t tokenColumn = tok.column;
   ExpectAllele( tok, minorAllele, "minor" );
   ExpectAllele( tok, majorAllele, "major" );
   if ( (majorAllele == '0') || (minorAllele == '0') )
      {
      // bad allele data for SNP.  Must have allele information
      Fatal( "Found missing SNP minor or major allele information in file [%s] near line %d:%d."
            "\n  found [%c] [%c]", Filename().c_str(), tokenLine, tokenColumn, minorAllele, majorAllele );
      }
   return;
   }

void CPlinkLexer::ExpectSnpDosage( CToken& tok, real& output )
   {
   size_t tokenLine = tok.line;
   size_t tokenColumn = tok.column;
   ExpectReal( tok, output, "dosage" );
   if ( (!fNoDosageRangeCheck) && (output != missingGenotypeValue) && ((output < 0.0) || (output > 2.0)) )
      {
      Fatal( "SNP dosage is out of range in file [%s] near line %d:%d."
           "\n  Expected a number between 0.0 and 2.0, but found %.13e",
           Filename().c_str(), tokenLine, tokenColumn, output );
      }
   }

void CPlinkLexer::ExpectProbability( CToken& tok, real& output, const char *description )
   {
   size_t tokenLine = tok.line;
   size_t tokenColumn = tok.column;
   ExpectReal( tok, output, description );
   if ( (output != missingGenotypeValue) && ((output < 0.0) || (output > 1.0)) )
      {
      Fatal( "%s SNP probability out of range in file [%s] near line %d:%d."
           "\n  Expected a number between 0.0 and 1.0, but found %.13e",
           description, Filename().c_str(), tokenLine, tokenColumn, output );
      }
   }

void CPlinkLexer::ExpectSnpProbabilities( CToken& tok, SnpDosageData& output )
   {
   size_t tokenLine = tok.line;
   size_t tokenColumn = tok.column;
   ExpectProbability( tok, output.probabilityOfHomozygousMinor, "Probability Of Homozygous Minor Genotype" );
   ExpectProbability( tok, output.probabilityOfHeterozygous, "Probability Of Heterozygous Genotype" );

   if ( ((output.probabilityOfHomozygousMinor == missingGenotypeValue) && (output.probabilityOfHeterozygous != missingGenotypeValue))
      ||((output.probabilityOfHomozygousMinor != missingGenotypeValue) && (output.probabilityOfHeterozygous == missingGenotypeValue))
      )
      {
      // bad genotype data.  Mixed a missing and a non-missing character
      //
      Warn( "Inconsistent missing genotype information found in file [%s] near line %d:%d."
          "\n  found [%f] [%f]", Filename().c_str(), tokenLine, tokenColumn, output.probabilityOfHomozygousMinor, output.probabilityOfHeterozygous );
      Warn( "    Genotype information treated as missing" );
      output.probabilityOfHomozygousMinor = missingGenotypeValue;
      output.probabilityOfHeterozygous = missingGenotypeValue;
      }
   else
      {
      real sum = output.probabilityOfHomozygousMinor + output.probabilityOfHeterozygous;
      if ( sum > 1.0 )
         {
         Fatal( "Expected sum of SNP probabilities to be no less than 0.0 and no greater than 1.0."
              "\n  Found %f and %f in file [%s] near line %d:%d", 
              output.probabilityOfHomozygousMinor, 
              output.probabilityOfHeterozygous, 
              Filename().c_str(),
              tokenLine, 
              tokenColumn );
         }
      }
   }

void CPlinkLexer::ExpectChromosome( CToken& tok, string& output, int& intOutput )
   {
   if ( tok.type != tokSymbol )     //   this is trickier as it can be and integer 0-26 or 'X', 'Y', 'XY', or 'MT'
      {
      Fatal( "Expected Chromosome in file [%s] near line %d:%d.  Found [%s]", Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
      }

   if ( _strcmpi( tok.text.c_str(), "x" ) == 0 ) intOutput = 23;
   else if ( _strcmpi( tok.text.c_str(), "y" ) == 0 ) intOutput = 24;
   else if ( _strcmpi( tok.text.c_str(), "xy" ) == 0 ) intOutput = 25;
   else if ( _strcmpi( tok.text.c_str(), "mt" ) == 0 ) intOutput = 26;
   else 
      {
//#if !defined( __GNUC__ )         // gcc doesn't yet support <regex>
//      //std::tr1::regex rxChromosome( "[0-9]+" );        // move outside loop for perf
//      if ( !std::tr1::regex_match( tok.text, rxChromosome ) )
      if ( !FValidateChromosomeString( tok.text ) )
         {
         char szFmt[] = "Expected a Chromosome in the form of an integer or a string X, Y, XY, or MT in file [%s] near line %d:%d.  Found [%s]";
         Fatal( szFmt, Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
         }
//#endif
      int i = atoi( tok.text.c_str() );
      if ( i > maxChromosomeValue )
         {
         char szFmt[] = "Expected a Chromosome in the form of an integer from 0-%d in file [%s] near line %d:%d.  Found [%s]";
         Fatal( szFmt, maxChromosomeValue, Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
         }
      intOutput = i;
      }

   output = tok.text;
   NextToken( tok );
   }
