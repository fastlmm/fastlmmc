#if !defined( CSnpFilterOptions_h )
#define Startup_h
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
 *         File Name:   CSnpFilterOptions.h
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   15 Mar 2012
 *     Revision Date:   15 Mar 2012
 *
 *    Module Purpose:   This file defines the CSnpFilterOptions
 *                         class interface for FastLmmC
 *
 *    Change History:   
 *
 * Test Files: 
 */
typedef  string   SnpId;

class SnpFilterOptions              // hold all the SnpFilter related parameters
   {
public:
   enum FilterType 
      {
      FilterByNone = 0
     ,FilterByCount                 // use 'n' SNPs
     ,FilterByJob
     ,FilterByFilelist
     ,FilterByTopN
     ,FilterByBlock
     ,FilterByBlocks                // use '2' blocks
     ,FilterByTask                  // snpPair but for job/task scheduling
     ,FilterByListOfSnpIds          // use vector< SnpId > snpIdsToExtract 
      }   filterType;

   int    cSnpsInJob;               // how many SNPs/Genotypes to process in one pass
   int    cJobsInFile;              // the number of iterations to complete all jobs
   int    thisJobIndex;             // the iteration we are processing this invocation
   int    numberOfJobs;

   size_t cSnpsInFile;              // Cache how many SNPs we read from the file
   size_t cBlocksMax;               // how many blocks to divide the SNPs in the file into
   double cSnpsInBlock;             // number of SNPs in a full block
   size_t iBlock1;                  // First block to process
   size_t iBlock2;                  // Second block to process

   size_t cTasksAsk;                // the number of tasks asked for 
   size_t iTaskAsk;                 // the task we asked for
   size_t cTasksMax;                // Computed max task.  MUST be less than cTaskAsk
                                    // we set cBlocksMax, iBlock1, iBlock2

   string extractFile;              // file with SNPid's to operate on.
   vector< SnpId > snpIdsToExtract; // list of SNPid's to extract

   SnpFilterOptions();
   void   SetBlockSize( size_t cSnpsRead );
   void   BlockRange( size_t idx, size_t& iStart, size_t& iEnd );
   size_t IdxFromBlockId( size_t blockId );
   size_t CountOfSnpsInBlock( size_t blockId );

   size_t IndexToFirstSnpInJob( size_t cSnpsToPartition );
   size_t IndexToLastSnpInJob( size_t cSnpsToPartition );

   bool FUseSnpFilter() { return( filterType != FilterByNone ); }
   bool FFilterByFile() { return( (filterType == FilterByFilelist) || (filterType == FilterByTopN) ); }
   bool FFilterByBlocks() { return( filterType == FilterByBlocks ); }
   bool FFilterByTask() { return( filterType == FilterByTask ); }
   bool FFilterBySnpIds() { return( (filterType == FilterByListOfSnpIds) || (filterType == FilterByFilelist) || (filterType == FilterByTopN) ); }

   void SetByCount( int cSnps );
   void SetByJob( int cJobs, int iJob );
   void SetByFileList( const string& filename );
   void SetByTopN( const string& filename, int cSnps );
   void SetByBlock( size_t cBlocks, size_t iBlock );
   void SetByBlocks( size_t cBlocks, size_t iBlk1, size_t iBlk2 );
   void SetByTask( int cTasks, int iTask );
   void SetByListOfSnpIds( const string& sourceId, const vector< SnpId >& rgSnpIds );

   pair< size_t, size_t> BlockIndexesFromITask( size_t iTask );      // get the two iBlocks that makeup 'iTask'

   void LoadSnpIdsToExtract();

   void Dump( FILE *pFile );
   void Dump( const string& filename );
   };

#endif   // CSnpFilterOptions_h
