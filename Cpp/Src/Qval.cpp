#include "FastLmmC.h"

class val_ind
   {
public:
   int index;     // index in original array
   double val;    //value
   static int compare(const void *x, const void *y)
      {
      if ((*(val_ind*)x).val>(*(val_ind*)y).val)
         {
         return( 1 );
         }
      else if((*(val_ind*)x).val<(*(val_ind*)y).val)
         {
         return( -1 );
         }
      else
         {
         return( 0 );
         }
      }
   };

//static function to compute the qvalues from the pvalues
//pvals: double array of size n
//index: int array of size n: holds the index of the sorted pvalues
void Qval::computeQvals( double *pvals, double *qvals, int *index, int npvals,double piNotm, int numTests )
   {
   vector<val_ind> pvals_ind(npvals);
     
   for( int i=0; i<npvals; ++i )
      {
      pvals_ind[i].val = pvals[i];
      pvals_ind[i].index = i;
      }

   qsort( &pvals_ind[0], npvals, sizeof(val_ind), val_ind::compare );

   for(int i=0; i<npvals; ++i )
      {
      pvals_ind[i].val = min( 1.0, ((pvals_ind[i].val * (1.0*numTests))/(1.0+i)) );
      }

   double currentMin = pvals_ind[npvals-1].val;
   for( int i=npvals-1; i>=0; --i)
      {
      currentMin = min( pvals_ind[i].val, currentMin );
      qvals[pvals_ind[i].index] = currentMin;
      }

   for( int i=0; i<npvals; ++i )
      {
      index[i] = pvals_ind[i].index;
      }
   }
