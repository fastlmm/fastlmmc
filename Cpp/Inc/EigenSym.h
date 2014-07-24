#if !defined( EigenSym_h )
#define EigenSym_h
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
 
class EigenSym
   {
   static const string headerfilename;
   static const string Ufilename;
   static const string Sfilename;
public:
   size_t n;      // Dimensionality of the original symmetric matrix ( K \in \R^[n*n] )
   size_t rankU;  // the rank of the reconstruction U*diag(S)*U^T
   real   *U;     // [n*rankU] matrix holding the Eigenvectors
   real   *S;     // [rankU]-vector holding the Eigenvalues   K=U*diag(S)*U^T
   string id;
   vector<string> columnLabelOrder;

   EigenSym();
   EigenSym( const string& id_, real *U_, real *S_, size_t n_, size_t rankU_, const vector<string>& columnLabelOrder_ );
   ~EigenSym();

   int set_fromKernel( const string& id_, real *K, size_t n_, size_t rankU_, const vector<string>& columnLabelOrder_ );
   int set_covFromX( const string& id_, real *X, size_t n_, size_t m, const vector<string>& columnLabelOrder_ );    //X -[n*m] matrix, with K=X*X^T
   int saveToDiskTAB( const string& dirname );           // saves an EigenSym object to disk in TAB delimited format
   int saveToDiskBINARY( const string& dirname );        // saves an EigenSym object to disk in BINARY delimited format
   int readFromDiskTAB( const string& dirname );         // reads an EigenSym Object from disk in TAB delim format
   int readFromDiskBINARY( const string& dirname_ );     // reads an EigenSym object from disk in BINARY format

   int saveReconstructionToDisk( StudyData *data, const string& filename );   //saves the reconstruction to disk

   static int SVD_economy( size_t n_, size_t d, real *A, real *U_, real *S_, real *VT );
   // computes the Eigenvectors of an [n by n] symmetric matrix A=U*S*U.T. Note that A will be overwritten.
   // eigenvalues are returned in ascending order.
   static int EV_symmetric( size_t n_, real *A, real *U_, real *S_ );
   
   // X: g x s
   // U: n x min(g,s)
   // S: min(g,s)
   // Z: n (holds the group indicator for individual i, entries in the range 0..g-1)
   static int EVsymmetric_ZXXZ( real* X, vector<size_t> Z, size_t n, size_t s, size_t g, real *U, real *S );
   
   // K: g x g
   // U: n x g
   // S: g
   // Z: n (holds the group indicator for individual i, entries in the range 0..g-1)
   static int EVsymmetric_ZKZ( real* K, vector<size_t> Z, size_t n, size_t g, real *U, real *S );

   static void treatNonPSD( real *S, size_t rankU, bool orderedEvAscending );
   void Dump( FILE *pf, const string& comment ) const;      // DEBUG support to Dump EigenSym to file
   };

#endif      // EigenSym_h
