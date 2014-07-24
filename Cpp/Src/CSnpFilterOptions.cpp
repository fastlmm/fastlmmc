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
 * CSnpFilterOptions - {SNP Filtering Options Class}
 *
 *         File Name:   CSnpFilterOptions.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   15 Mar 2012
 *     Revision Date:   15 Mar 2012
 *
 *    Module Purpose:   This file implements the CSnpFilterOptions
 *                         class for FastLmmC
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Include Files
 */
#include "FastLmmC.h"

SnpFilterOptions::SnpFilterOptions()
   {
   /*
    * Initialize the member variables.
    */
   {filterType=FilterByNone; cSnpsInJob=-1; cJobsInFile=-1; thisJobIndex=-1; numberOfJobs=-1; 
                       cSnpsInFile=-1; cBlocksMax=-1; cSnpsInBlock=-1.0; iBlock1=-1; iBlock2=-1;
                       cTasksAsk=0; iTaskAsk=0; cTasksMax=0;}
   filterType   = FilterByNone;  // old init values
   cSnpsInJob   = 0;             // -1
   cJobsInFile  = 0;             // -1
   thisJobIndex = 0;             // -1
   numberOfJobs = 0;             // -1
   cSnpsInFile  = 0;             // -1
   cBlocksMax   = 0;             // -1
   cSnpsInBlock = 0.0;           // -1.0
   iBlock1      = 0;             // -1
   iBlock2      = 0;             // -1
   cTasksAsk    = 0;             // 0
   iTaskAsk     = 0;             // 0
   cTasksMax    = 0;             // 0
   }

// TODO:  fix for job handling too.
size_t SnpFilterOptions::IndexToFirstSnpInJob( size_t cSnpsToPartition ) 
   {
   cSnpsInJob = ((int)cSnpsToPartition + numberOfJobs - 1) / numberOfJobs;
   return( cSnpsInJob * thisJobIndex );
   }

size_t SnpFilterOptions::IndexToLastSnpInJob( size_t cSnpsToPartition )
   {
   cSnpsInJob = ((int)cSnpsToPartition + numberOfJobs - 1) / numberOfJobs;
   size_t iSnpLast = ( cSnpsInJob * (thisJobIndex+1) );
   iSnpLast = min( iSnpLast, cSnpsToPartition );
   return( iSnpLast );
   }

void SnpFilterOptions::SetBlockSize( size_t cSnpsRead )
   {
   cSnpsInFile = cSnpsRead;
   cSnpsInBlock = (double)cSnpsRead / (double) cBlocksMax;
   if ( cSnpsInBlock < 10.0 )
      {
//      Fatal( "Cannot block a file to have less than 10 SNPs in a block.  Found %6.3lf", cSnpsInBlock );
      }
   }

void SnpFilterOptions::BlockRange( size_t idx, size_t& iStart, size_t& iEnd )
   {
   iStart = IdxFromBlockId( idx );
   iEnd   = IdxFromBlockId( idx+1 );
   }

// return the index value in the original file for the beginning of this block
size_t SnpFilterOptions::IdxFromBlockId( size_t blockId )
   {
   if ( (blockId < 0) || (blockId > cBlocksMax) )
      {
      Fatal( "Block out of range.  Expected less than %d, received %d", cBlocksMax, blockId );
      }
   size_t retValue = (size_t)((cSnpsInBlock * (double)blockId) + 0.5);
   retValue = min( retValue, cSnpsInFile );
   return( retValue );
   }

// Block are now a variable size (rounding error... :-()  
// Compute the size of a specific block
size_t SnpFilterOptions::CountOfSnpsInBlock( size_t blockId )
   {
   size_t cSnps = IdxFromBlockId( blockId+1 ) - IdxFromBlockId( blockId );
   return( cSnps );
   }

void SnpFilterOptions::SetByCount( int cSnps )
   {
   if ( cSnps <= 0 )
      {
      Fatal( "FilterOption cSnpsInJob must be a positive number.  Found [%d]", cSnps );
      }
   filterType = FilterByCount;
   cSnpsInJob = cSnps;
   }

void SnpFilterOptions::SetByJob( int cJobs, int iJob )
   {
   if ( cJobs <= 0 )
      {
      Fatal( "FilterOption numberOfJobs must be a positive number. Found [%d]", cJobs );
      }
   filterType = FilterByJob;
   numberOfJobs = cJobs;
   if ( iJob >= numberOfJobs )
      {
      Fatal( "FilterOption thisJob must be less than the numberOfJobs argument [%d]. Found [%d]", numberOfJobs, iJob );
      }
   thisJobIndex = iJob;
   }

void SnpFilterOptions::SetByBlock( size_t cBlocks, size_t iBlock )
   {
   if ( cBlocks <= 0 )
      {
      Fatal( "Number of blocks to operate on must be a positive number. Found [%d]", cBlocks );
      }

   filterType = FilterByBlock;
   cBlocksMax = cBlocks;

   if ( (iBlock < 0) || (iBlock >= cBlocksMax) )
      {
      Fatal( "The block index must be less than the number of blocks to operate on[%d]. Found [%d]", cBlocksMax, iBlock );
      }
   iBlock1 = iBlock;
   }

void SnpFilterOptions::SetByBlocks( size_t cBlocks, size_t iBlk1, size_t iBlk2 )
   {
   if ( cBlocks <= 0 )
      {
      Fatal( "Number of blocks to operate on must be a positive number. Found [%d]", cBlocks );
      }

   filterType = FilterByBlocks;
   cBlocksMax = cBlocks;

   if ( (iBlk1 < 0) || (iBlk1 >= cBlocksMax) )
      {
      Fatal( "The block1 index must be less than the number of blocks to operate on[%d]. Found [%d]", cBlocksMax, iBlk1 );
      }
   iBlock1 = iBlk1;

   if ( (iBlk2 < iBlock1) || (iBlk1 >= cBlocksMax) )
      {
      Fatal( "The block2 index must be greater or equal to block1[%s] and less than the number of blocks[%d]. Found [%d]", iBlock1, cBlocksMax, iBlk2 );
      }
   iBlock2 = iBlk2;
   }

void SnpFilterOptions::SetByFileList( const string& filename )
   {
   filterType = FilterByFilelist;
   extractFile = filename;
   cSnpsInJob = std::numeric_limits<int>::max();      // 
   }

void SnpFilterOptions::SetByTopN( const string& filename, int cSnps )
   {
   if ( filename.size() == 0 )
      {
      Fatal( "FilterOption 'SetByTopN' requires a valid filename.  None found." );
      }
   if ( cSnps <= 0 )
      {
      Fatal( "FilterOption cSnpsInJob must be a positive number.  Found [%d]", cSnps );
      }
   filterType = FilterByTopN;
   cSnpsInJob = cSnps;
   extractFile = filename;
   }

void SnpFilterOptions::SetByTask( int cTasks, int iTask )
   {
   if ( cTasks <= 0 )
      {
      Fatal( "FilterOption numberOfTasks must be a positive number. Found [%d]", cTasks );
      }
   if ( iTask >= cTasks )
      {
      Fatal( "FilterOption iTask must be less than numberOfTasks[%d].  Found [%d]", cTasks, iTask );
      }

   cTasksAsk = cTasks;
   iTaskAsk = iTask;
   cBlocksMax = (int)(sqrt( (real)(cTasksAsk*2) ));    // 
   if ( ((cBlocksMax * ( cBlocksMax + 1 )) / 2 ) > (size_t)cTasks )
      {
      --cBlocksMax;
      }

   cTasksMax = (cBlocksMax * ( cBlocksMax + 1 )) / 2;            // compute how many tasks really needed to complete job
   if ( (size_t)iTask < cTasksMax )
      {
#if 0
      pair< size_t, size_t > taskBlocks = BlockIndexesFromITask( iTask );
      iBlock1 = taskBlocks.first;
      iBlock2 = taskBlocks.second;
#else
      // this is a valid index and needs the block id's set
      iBlock1 = 0;
      iBlock2 = 0;
      size_t iJob = 0; //iTask;
      size_t cBlocksInRow = cBlocksMax;

      for( ; iBlock1<cBlocksMax; ++iBlock1 )
         {
         if ( (iJob + (cBlocksMax - iBlock1)) > (size_t)iTask )
            {
            iBlock2 = iBlock1 + (iTask - iJob);
            break;
            }
         iJob += (cBlocksMax - iBlock1);
         }
#endif
      }
   else
      {
      iBlock1 = iBlock2 = cBlocksMax;    // initialize to end of data
      }
   filterType = FilterByTask;
   }


void SnpFilterOptions::SetByListOfSnpIds( const string& sourceId, const vector< SnpId >& rgSnpIds )
   {
   extractFile = sourceId;
   snpIdsToExtract = rgSnpIds;
   cSnpsInJob = std::numeric_limits<int>::max();      // do not truncate Snps early
   filterType = FilterByListOfSnpIds;
   }

pair< size_t, size_t> SnpFilterOptions::BlockIndexesFromITask( size_t iTask )
   {
   pair< size_t, size_t > taskBlocks( 0, 0 );

   for ( taskBlocks.first = 0; taskBlocks.first < cBlocksMax; ++taskBlocks.first )
      {
      for ( taskBlocks.second = taskBlocks.first; taskBlocks.second < cBlocksMax; ++taskBlocks.second )
         {
         if ( iTask-- == 0 ) return( taskBlocks );
         }
      }
   return( taskBlocks );
   }

void SnpFilterOptions::LoadSnpIdsToExtract()
   {
   if ( !FFilterByFile() )
      {
      return;
      }

   SnpId idSnp;
   map< SnpId, size_t > detectDuplicatSnpIds;

   if ( filterType == FilterByFilelist )
      {
      ProgressNL( "  --               Process -extract: [%s]", extractFile.c_str() );
      }
   else
      {
      ProgressNL( " --            Process -extractTopK: [%s] %d", extractFile.c_str(), cSnpsInJob );
      }
   CPlinkLexer lex( extractFile );
   CToken tok;

   lex.NextToken( tok );               // prime the lexer
   while ( tok.type != tokEOF )        // for each line until EOF
      {
      if ( tok.type == tokEOL )        // blank line...or comment line or ???
         {
         lex.NextToken( tok );
         continue;
         }

      lex.ExpectId( tok, idSnp, "Snp ID" );
      if ( tok.type != tokEOL )
         {
         Fatal( "Expecting <EOL> after Snp ID in file [%s] near line %d:%d.  Found [%s]", lex.Filename().c_str(), tok.line, tok.column, tok.text.c_str() );
         }

      lex.NextToken( tok );      // consume the <EOL>

      if ( detectDuplicatSnpIds.count( idSnp ) )
         {
         Fatal( "Found duplicate Snp ID [%s] in file [%s] near lines %d and %d", idSnp.c_str(), lex.Filename().c_str(), detectDuplicatSnpIds[idSnp]+1, snpIdsToExtract.size()+1 );
         }

      // Valid idSnp  Add it and move on.
      detectDuplicatSnpIds[ idSnp ] = snpIdsToExtract.size();     // add to the duplicate detection table
      snpIdsToExtract.push_back( idSnp );

      if ( snpIdsToExtract.size() >= (size_t)cSnpsInJob )
         {
         break;
         }
      }

   ProgressNL( "  --                   SnpIds found: %d", snpIdsToExtract.size() );
   }

void SnpFilterOptions::Dump( FILE *pFile )
   {
   // GCC and VC don't share the same format specifiers for size_t 
   //   so we convert to use int in this case.
   fprintf( pFile,   "#  FileBaseName: %s", plinkFilenameRoot.c_str() );
   fprintf( pFile, "\n# cSnps In File: %d", (int)cSnpsInFile );
   fprintf( pFile, "\n#    cBlocksMax: %d", (int)cBlocksMax );
   fprintf( pFile, "\n#  cSnpsInBlock: %10.4lf", cSnpsInBlock );
   fprintf( pFile, "\n# PartitionList:" );
   fprintf( pFile, "\n# Partition\tStart\tCount\tMax" );
   for( size_t iPartition=0; iPartition < cBlocksMax; ++iPartition )
      {
      size_t idxPartitionStart = IdxFromBlockId( iPartition );
      size_t csnpPartition = CountOfSnpsInBlock( iPartition );
      size_t idxPartitionMax = IdxFromBlockId( iPartition+1 );
      fprintf( pFile, "\n%d\t%d\t%d\t%d", (int)iPartition, (int)idxPartitionStart, (int)csnpPartition, (int)idxPartitionMax );
      }
   fprintf( pFile, "\n" );
   }

void SnpFilterOptions::Dump( const string& filename )
   {
   FILE *pFile;

   pFile = fopen( filename.c_str(), "w" );
   if ( pFile == nullptr )
      {
      Warn( "Unable to open file [%s]."
          "\n   CRT error: %s", filename.c_str(), strerror( errno ) );
      return;
      }

   Dump( pFile );
   fclose( pFile );
   }
