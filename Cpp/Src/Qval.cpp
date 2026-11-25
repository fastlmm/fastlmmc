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

class val_ind
   {
public:
   size_t index;     // index in original array
   real val;         //value
   static void Display( const val_ind* p1, const val_ind* p2 )
      {
      // memcpy (not a pointer cast) to view the bits without breaking strict aliasing
      long long bits1 = 0, bits2 = 0;
      memcpy( &bits1, &p1->val, sizeof( p1->val ) );
      memcpy( &bits2, &p2->val, sizeof( p2->val ) );
      Progress( "\ncompare: p1->index [%d] p1->val [%.16E] 0x%16llx", p1->index, p1->val, bits1 );
      Progress( "\n         p2->index [%d] p2->val [%.16E] 0x%16llx", p2->index, p2->val, bits2 );
      }

   static int compare(const void *x, const void *y)
      {
      val_ind *p1 = (val_ind*)x;
      val_ind *p2 = (val_ind*)y;

      if ( p1->val > p2->val )
         {
         return( 1 );
         }
      else if( p1->val < p2->val )
         {
         return( -1 );
         }

      // index must be unique.  Use it to make this a deterministic sort
      else if ( p1->index > p2->index )
         {
//         Display( p1, p2 );
         return( 1 );
         }
      else if ( p1->index < p2->index )
         {
//         Display( p1, p2 );
         return( -1 );
         }
      else
         {
//         Display( p1, p2 );
         return( 0 );
         }
      }
   };

// static function to compute the qvalues from the pvalues
// pvals: real array of size n
// index: int array of size n: holds the index of the sorted pvalues
void Qval::computeQvals( real *pvals, real *qvals, size_t *index, size_t npvals )
   {
   vector<val_ind> pvals_ind( npvals );
   vector<size_t> pvals_notobs_ind( npvals );

   size_t n_noNaN = 0;     // count the number of SNPs that are not NaN
   size_t n_NaN = 0;       // count the number of SNPs that are NaN
   for( size_t i=0; i<npvals; ++i )
      {
      if ( pvals[i] == pvals[i] )
         {
         //there might be NaNs among the p-values
         pvals_ind[ n_noNaN ].val = pvals[ i ];
         pvals_ind[ n_noNaN ].index = i;        //mapping back to original SNP-index
         ++n_noNaN;
         }
      else
         {
         qvals[ i ] = std::numeric_limits<real>::quiet_NaN();
         pvals_notobs_ind[ n_NaN ] = i;
         ++n_NaN;
         }
      }
   if ( n_noNaN > 0 )
      {
      qsort( &pvals_ind[0], n_noNaN, sizeof(val_ind), val_ind::compare );

      for( size_t i=0; i<n_noNaN; ++i )
         {
         pvals_ind[i].val = min( (real)1.0, ((pvals_ind[i].val * ((real)1.0*n_noNaN)) / ((real)1.0+i)) );
         }

      real currentMin = pvals_ind[ n_noNaN-1 ].val;
      for( int i=(int)n_noNaN-1; i >= 0; --i )
         {
         currentMin = min( pvals_ind[ i ].val, currentMin );
         qvals[pvals_ind[ i ].index] = currentMin;
         }

      //for the argsort put the observed ones to the front in order of pvalue magnitude
      for( size_t i=0; i<n_noNaN; ++i )
         {
         index[ i ] = pvals_ind[ i ].index;
         }
      }

   //append the
   for( size_t i=0; i<n_NaN; ++i )
      {
      index[ i+n_noNaN ] = pvals_notobs_ind[ i ];
      }
   } 
