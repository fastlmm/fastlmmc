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

const string EigenSym::headerfilename="IndividualID_header.txt";
const string EigenSym::Ufilename="U";
const string EigenSym::Sfilename="S";

void DumpBlasLapackSyevrArgs( const char *pJobz,               // 1 jobz
                              const char *pRange,              // 2 range
                              const char *pUpLo,               // 3 uplo
                              int  *pN,                        // 4 n
                              real *pA,                        // 5 a
                              int  *pLeadingDimensionA,        // 6 lda
                              real *pValLower,                 // 7 vl
                              real *pValUpper,                 // 8 vu
                              int  *pIndexLower,               // 9 il
                              int  *pIndexUpper,               // 10 iu
                              real *pAbsoluteErrorTolerance,   // 11 abstol
                              int  *pEigenValuesFound,         // 12 m
                              real *pSelectedEigenValues,      // 13 w
                              real *pSelectedEigenVectors,     // 14 z
                              int  *pLeadingDimensionZ,        // 15 ldz
                              int  *pDimensionWork,            // 16 isuppz
                              real *pWorkArray,                // 17 work
                        const int  *pLWork,                    // 18 lwork
                              int  *pIWork,                    // 19 iwork
                        const int  *pLIWork,                   // 20 liwork
                              int  *pInfo,                     // 21 info

                              const string& comment
                            //  , bool fOpenAppend
                              );

EigenSym::EigenSym()
   {
   // force initial values to be consistent.
   n = 0;
   rankU = 0;
   U = nullptr;
   S = nullptr;
   id.clear();
   columnLabelOrder.clear();
   }

EigenSym::~EigenSym()
   {
   if ( U != nullptr )
      {
      delete[] U;
      U = nullptr;
      }
   if ( S != nullptr )
      {
      delete[] S;
      S = nullptr;
      }
   rankU = 0;
   n = 0;
   }

//saves an EigenSym object to disk in TAB delimited format
int EigenSym::saveToDiskTAB( const string& Filename )
   {
   Fatal( "EigenSym::saveToDiskTAB() --- Not implemented." );
   return 0;
   }                        

//saves an EigenSym object to disk in BINARY delimited format
int EigenSym::saveToDiskBINARY( const string& dirname_ )
   {
   string dirname = (dirname_.length() != 0 ) ? dirname_ : ".";

   MatrixIO::createDirectory( dirname );

   string filenameHeader = dirname;
   filenameHeader.append( MatrixIO::OS_directoryDelim );
   filenameHeader.append( EigenSym::headerfilename );
   MatrixIO::writeHeaderFile( filenameHeader, this->columnLabelOrder );

   string filenameU = dirname;
   filenameU.append( MatrixIO::OS_directoryDelim );
   filenameU.append( EigenSym::Ufilename );
   
   string filenameS = dirname;
   filenameS.append( MatrixIO::OS_directoryDelim );
   filenameS.append( EigenSym::Sfilename );
   
   MatrixIO::writeMatrix2dBIN( filenameU, id, n, rankU, U );
   MatrixIO::writeMatrix2dBIN( filenameS, id, rankU, 1, S );
   return 0;
   }                     

//reads an EigenSym Object from disk in TAB delim format
int EigenSym::readFromDiskTAB( const string& dirname )
   {
   Fatal( "EigenSym::readFromDiskTAB() --- Not implemented." );
   return 0;
   }           

//reads an EigenSym object from disk in BINARY format
int EigenSym::readFromDiskBINARY( const string& dirname_ )
   {
   string dirname = (dirname_.length() != 0 ) ? dirname_ : ".";

   string filenameHeader = dirname;
   filenameHeader.append( MatrixIO::OS_directoryDelim );
   filenameHeader.append( EigenSym::headerfilename );
   this->columnLabelOrder = MatrixIO::readHeaderFile( filenameHeader );

   //read U
   string filenameU = dirname;
   filenameU.append( MatrixIO::OS_directoryDelim );
   filenameU.append( EigenSym::Ufilename );
   
   MatrixIO::querySizeMatrixBIN( filenameU, &id, &n, &rankU );
   
   if ( this->U != 0)
      {
      delete[] U;
      }
   this->U = new real[ n*rankU ];
   
   MatrixIO::readBinary2dBIN( filenameU, id, n, rankU, U );

   //read S
   string filenameS = dirname;
   filenameS.append( MatrixIO::OS_directoryDelim );
   filenameS.append( EigenSym::Sfilename );
   
   string idS;
   size_t nS = 0;
   size_t mS = 0;
   MatrixIO::querySizeMatrixBIN( filenameS, &idS, &nS, &mS );

   if ( id != idS )
      {
      Fatal( "The ids of U and S do not match up. idU:%s idS:%s", id.c_str(), idS.c_str() );
      }
   if ( nS != rankU )
      {
      Fatal( "The first dimension of matrix S has to be rankU=%i. Found %i", rankU, nS );
      }
   if ( mS != 1)
      {
      Fatal( "The second dimension of matrix S has to be 1. Found %i", mS );
      }
   if ( this->columnLabelOrder.size() != this->n )
      {
      Fatal( "The number of columnLabels=%i has to be n=%i.", this->columnLabelOrder.size(), n );
      }
   if ( this->S != 0 )
      {
      delete[] S;
      }
   this->S = new real[ rankU ];
   
   MatrixIO::readBinary2dBIN( filenameS, idS, rankU, mS, S );
   return 0;
   }        

//K: g x g  symmetric matrix in lower storage
//U: n x g
//S: g
//Z: n (holds the group indicator for individual i, entries in the range 0..g-1)
int EigenSym::EVsymmetric_ZKZ( real* K, vector<size_t> Z, size_t n, size_t g, real *U, real *S )
   {
   real *ZZ= new real[g];     //This matrix holds the group member counts.
   
   for ( size_t group=0; group<g; ++group )
      {
      ZZ[ group ] = 0.0;      // zero init the new array
      }

   for ( size_t ind=0; ind<n; ++ind )
      {
      if ( Z[ind] >= g )
         {
         Fatal( "group index %i of individual %i is larger than number of groups %i", Z[ind], ind, g );
         }
      ++ZZ[Z[ind]];
      }

   for ( size_t group=0; group<g; ++group )
      {
      ZZ[group] = sqrt( ZZ[group] );
      }

   for ( size_t column = 0; column<g; ++column )
      {
      for (size_t row=column; row<g; ++row )
         {
         K[row+(column*g)] *= ZZ[row]*ZZ[column];  //adjust the lower triangular part of K.
         }
      }
   real *ZZU= new real[g*g];
   int ret = EigenSym::EV_symmetric( g, K, ZZU, S );

   EigenSym::treatNonPSD( S, g, true );

   for ( size_t group=0; group<g; ++group )
      {
      for ( size_t ind=0; ind<n; ++ind )
         {
         U[ind+(group*n)] =ZZU[Z[ind]+(group*g)] / ZZ[Z[ind]];

         if ( useWoodbury )
            {
            U[ind+group*n] *= sqrt( S[group] );
            }
         }
      }
   delete[] ZZU;
   delete[] ZZ;
   return ret;
   }

//X: g x s
//U: n x min(g,s)
//S: min(g,s)
//Z: n (holds the group indicator for individual i, entries in the range 0..g-1)
int EigenSym::EVsymmetric_ZXXZ( real* X, vector<size_t> Z, size_t n, size_t s, size_t g, real *U, real *S )
   {
   if ( s >= g )
      {
      real *K= new real[ g*g ];
      BLAS_INT g_int = (BLAS_INT)g; //Warning: size_t BLAS_INT conversion
      BLAS_INT s_int = (BLAS_INT)g; //Warning: size_t BLAS_INT conversion
      
      BlasLapack::syrk( &char_L, &char_N, &g_int, &s_int, &real_one, X, &g_int, &real_zero, K, &g_int );
      int ret = EigenSym::EVsymmetric_ZKZ( K, Z, n, g, U, S );
      delete[] K;
      return ret;
      }

   real *ZZ= new real[g];           //This matrix holds the group member counts.
   real *ZZi= new real[g];          //This matrix holds the inverse group member counts.
   
   for ( size_t group=0; group<g; ++group )
      {
      ZZ[ group ] = 0.0;
      }

   for ( size_t ind=0; ind<n; ++ind )
      {
      if ( Z[ind] >= (size_t)g )
         {
         Fatal( "Group index %i of individual %i is larger than number of groups %i", Z[ind], ind, g );
         }
      ++ZZ[Z[ind]];
      }

   for ( size_t group=0; group<g; ++group )
      {
      Fatal( "\nTODO." );
#if 0
      if (ZZ)
      ZZ[group]=sqrt(ZZ[group]);
      ZZi[group]=(real)1.0/ZZ[group];
#endif
      }

   for ( size_t column=0; column<s; ++column )
      {
      for ( size_t row =0; row<g; ++row )
         {
         X[row+(column*g)] *= ZZ[row];    //adjust the lower triangular part of K.
         }
      }

   real *ZZU= new real[g*s];
   real *ZZVT= new real[s*s];
   int ret = EigenSym::SVD_economy( g, s, X, ZZU, S, ZZVT );
   for ( size_t group=0; group<g; ++group )
      {
      for ( size_t ind=0; ind<n; ++ind )
         {
         U[ind+(group*n)] = ZZU[Z[ind]+(group*g)] * ZZi[Z[ind]];
         }
      }

   delete[] ZZU;
   delete[] ZZ;
   delete[] ZZi;
   delete[] ZZVT;
   return ret;
   }

//X -[n*m] matrix, with K=X*X^T  --->K[n*n];  rankU=min(n,m)
int EigenSym::set_covFromX( const string& id_, real* X, size_t n_, size_t m, const vector<string>& columnLabelOrder_ )
   {
   if ( columnLabelOrder_.size() != n_ )
      {
      Fatal( "The number of columnLabels=%i does not match the number of individuals n_=%i.", columnLabelOrder_.size(), n_ );
      }
   n = n_;
   rankU = min( n, m );    //The number of singular vectors (rankU) is the minimum of the two dimensions of the SNP sample matrix.
   if ( U != nullptr )
      {
      delete[] U;
      U = nullptr;
      }
   if ( S != nullptr )
      {
      delete[] S;
      S = nullptr;
      }
   if ( rankU == n )
      {  //full rank case
      ProgressNL( "        Create covariance from SNPs: ..." );
      real *K = new real[ n * n ];
      BLAS_INT n_int = (BLAS_INT)n;
      BLAS_INT m_int = (BLAS_INT)m;
      BlasLapack::syrk( &char_L, &char_N, &n_int, &m_int, &real_one, X, &n_int, &real_zero, K, &n_int );
      
      Progress( " done." );
      set_fromKernel( id_, K, n, rankU, columnLabelOrder_ );
      delete[] K;
      }
   else
      {   
      // low rank case
      if ( useWoodbury )
         {
         ////ProgressNL( "          Compute Woodbury low rank: ..." );
         columnLabelOrder = columnLabelOrder_;
         id = id_;
         U = new real[ n * m ];           //the left singular vectors of the SNP sample matrix are the Eigenvectors of the Covariance Kernel matrix.
         S = new real[ m ];

         // K is stored in U.
         BLAS_INT n_int = (BLAS_INT)n;
         BLAS_INT m_int = (BLAS_INT)m;
         BlasLapack::syrk( &char_L, &char_T, &m_int, &n_int, &real_one, X, &n_int, &real_zero, U, &m_int );
         real *V = new real[ m * m ];     //the right singular vectors of the SNP sample matrix

         //Note, that 'U' here is really 'K'
         EigenSym::EV_symmetric( m, U, V, S );
         BlasLapack::gemm( &char_N, &char_N, &n_int, &m_int, &m_int, &real_one, X, &n_int, V, &m_int, &real_zero, U, &n_int );
         delete[] V;    //We only need the left singular vectors, so this part can go.
         
         if ( useWoodbury == false )
            {
            for ( size_t col=0; col<m; ++col )
               {
               real sqrtS = (real)0.0;
               if ( S[col] > (real)MathFunctions::eps_rank )
                  {
                  sqrtS = sqrt( S[col] );
                  }
               for ( size_t row=0; row<n; ++row )
                  {
                  U[ row + (col * n) ] /= sqrtS;
                  }
               }
            }
         Progress( " done." );
         }
      else
         {
         //ProgressNL( "                   Compute low rank: ..." );
         columnLabelOrder = columnLabelOrder_;
         id = id_;
         U = new real[ n * rankU ];             // the left singular vectors of the SNP sample matrix are the Eigenvectors of the Covariance Kernel matrix.
         S = new real[ rankU ];
         real *VT = new real[ rankU * m ];      // the right singular vectors of the SNP sample matrix
         size_t size_x = (size_t)n * (size_t)m;

         real *X_ = new real[ size_x ];
         for ( size_t i=0 ; i<size_x; ++i )
            {
            X_[ i ] = X[ i ];
            }
         EigenSym::SVD_economy( n, m, X_, U, S, VT ); // The economy SVD returns the first rankU Singular vectors and Singular values
          
         delete[] VT;                                 // We only need the left singular vectors, so this part can go.
         delete[] X_;
         for ( size_t i=0; i<rankU; ++i )
            {
            S[ i ] = S[ i ] * S[ i ];                 // the Eigenvalues of the Kernel are the squared Singular values of the SNP-matrix
            }
         }
      }
   return 0;
   }  

EigenSym::EigenSym( const string& id_, real *U_, real *S_, size_t n_, size_t rankU_, const vector<string>& columnLabelOrder_ )
   {
   if ( columnLabelOrder_.size() != n_ )
      {
      Fatal( "The number of columnLabels=%i does not match the number of individuals n=%i.",columnLabelOrder_.size(), n_ );
      }
   columnLabelOrder = columnLabelOrder_;
   id = id_;
   if ( rankU_ > n_ )
      {
      Fatal( "rankU=%i>n%i. This might be a bug!", rankU_, n_ );
      }
   rankU = rankU_;
   n = n_;
   U = U_;
   S = S_;
   }

int EigenSym::set_fromKernel( const string& id_, real *K, size_t n_, size_t rankU_, const vector<string>& columnLabelOrder_ )
   {
   if ( columnLabelOrder_.size() != n_ )
      {
      Fatal( "The number of columnLabels=%i does not match the number of individuals n=%i.", columnLabelOrder_.size(), n_ );
      }
   id = id_;
   columnLabelOrder = columnLabelOrder_;
   //ProgressNL( "   Compute Eigenvectors from Kernel: ..." );
   if ( U != nullptr )
      {
      delete[] this->U;
      }
   if ( S != nullptr )
      {
      delete[] this->S;
      }
   n  =  n_;
   rankU =   min( n, rankU_ );
   if ( n == rankU )
      {
      U = new real[ n * n ];
      S = new real[ n ];
      ZeroRealArray( U, n * n );
      ZeroRealArray( S, n );

      //Verbose( "Running eigen" );
      EigenSym::EV_symmetric( n, K, U, S );
      }
   else
      {
      Fatal( "Low rank Eigenvectors of Kernel are not implemented yet. Use Lanczos algorithm here." );
      }
   
   //optionally set the rank of the kernel to a fixed value
   if ( setRank > 0 )
      {
      for ( size_t i=0; i<n-setRank; ++i)
         {
         //The eigenvalues have inverse ordering
         S[ i ] = 0.0;
         }
      }

   treatNonPSD( S, rankU, true );
//   ProgressNL( "           End compute EigenVectors:" );
   return 0;
   }

//computes the economy sized SVD of an [n by d] matrix A=U*S*VT. Note that A will be overwritten.
//eigenvalues are returned in descending order.
int EigenSym::SVD_economy( size_t n_, size_t d, real *A, real *U_, real *S_, real *VT )
   {
   BLAS_INT n_int = (BLAS_INT)n_;
   BLAS_INT d_int = (BLAS_INT)d;
   BLAS_INT dmin  = (BLAS_INT)min( n_, d );
   real work_svd_query;
   BLAS_INT lwork_svd  = -1;
   BLAS_INT *iwork_svd =new BLAS_INT[ 8 * dmin ];     // to be deleted

   BLAS_INT info_svd;
   BlasLapack::gesdd( &char_S, &n_int, &d_int, A, &n_int, S_, U_, &n_int, VT, &dmin, &work_svd_query,  &lwork_svd, iwork_svd, &info_svd ); //workspace query
   lwork_svd = (BLAS_INT)work_svd_query;
   real *work_svd = new real[ lwork_svd ];            // to be deleted
   BlasLapack::gesdd( &char_S, &n_int, &d_int, A, &n_int, S_, U_, &n_int, VT, &dmin, work_svd,  &lwork_svd, iwork_svd, &info_svd );        //execute the SVD

   delete[] work_svd;
   delete[] iwork_svd;

   return( info_svd );
   }

// Computes the Eigenvectors of an [n by n] symmetric matrix A=U*S*U.T. 
// Note that A will be overwritten.
// eigenvalues are returned in ascending order.
int EigenSym::EV_symmetric( size_t n_, real *A, real *U_, real *S_ )
   {
   CTimer eigenTimer( true );
   BLAS_INT int_n = (BLAS_INT)n_; //Warning: size_t BLAS_INT conversion
   real vl = 0.0;                // the range of eigenvalues (not used in our scenario)
   real vu = 0.0;
   BLAS_INT il = 0;              // the indices of eigenvalues (not used in our scenario)
   BLAS_INT iu = 0;
   real abstol = (real)3.0E-8;            //tolerance of the EV algorithm
 //TODO:  bobd  - confirm what real precision should be
 //  real abstol = (real)3.0E-13;           //tolerance of the EV algorithm

   BLAS_INT m = 0;                        // out: the number of eigenvalues found, should be n_.
   BLAS_INT* isuppz = new BLAS_INT[2*n_];  // ---: second integer workspace (deleted)
   real workdummy = 0.0;                  // out: holds the optimal size of the workspace
   BLAS_INT iworkdummy = 0;               // out: holds the optimal size of the BLAS_INT workspace
   BLAS_INT info = 0;                     // out: holds status/error indicator syevr()
   memset( isuppz, 0, (n_ * 2 * sizeof( BLAS_INT )) );

 DumpBlasLapackSyevrArgs( &char_V, &char_A, &char_L, &int_n, A, &int_n, &vl, &vu, &il, &iu, &abstol, &m, S_, U_, &int_n, isuppz, &workdummy, &int_minusone,  &iworkdummy, &int_minusone, &info, "before optimal workspace call"  );
   BlasLapack::syevr( &char_V, &char_A, &char_L, &int_n, A, &int_n, &vl, &vu, &il, &iu, &abstol, &m, S_, U_, &int_n, isuppz, &workdummy, &int_minusone, &iworkdummy, &int_minusone, &info );//optimal workspace request
 DumpBlasLapackSyevrArgs( &char_V, &char_A, &char_L, &int_n, A, &int_n, &vl, &vu, &il, &iu, &abstol, &m, S_, U_, &int_n, isuppz, &workdummy, &int_minusone,  &iworkdummy, &int_minusone, &info, "After optimal workspace call"  );

   BLAS_INT lwork=(BLAS_INT)workdummy;    //optimal size of the workspace
   BLAS_INT liwork=iworkdummy;            //optimal size of the BLAS_INTeger workspace
   real* work=new real[lwork];            //workspace of optimal size (deleted)
   BLAS_INT* iwork= new BLAS_INT[liwork]; //integer workspace of optimal size (deleted)

   ZeroRealArray( work, lwork );
   memset( iwork, 0, (liwork*(sizeof(BLAS_INT))) );

 //DumpEigenSymState( this, "Before syevr", false );
 DumpBlasLapackSyevrArgs( &char_V, &char_A, &char_L, &int_n, A, &int_n, &vl, &vu, &il, &iu, &abstol, &m, S_, U_, &int_n, isuppz, work, &lwork,  iwork, &liwork, &info, "Before Eigenvectors computed"  );
   BlasLapack::syevr( &char_V, &char_A, &char_L, &int_n, A, &int_n, &vl, &vu, &il, &iu, &abstol, &m, S_, U_, &int_n, isuppz, work, &lwork,  iwork, &liwork, &info );//Eigenvectors and Eigenvalues are computed
 DumpBlasLapackSyevrArgs( &char_V, &char_A, &char_L, &int_n, A, &int_n, &vl, &vu, &il, &iu, &abstol, &m, S_, U_, &int_n, isuppz, work, &lwork,  iwork, &liwork, &info, "After Eigenvectors computed"  );
 //DumpEigenSymState( this, "After syevr", true );

   delete[] work;       //delete workspace
   delete[] iwork;      //delete integer workspace
   delete[] isuppz;     //delete second integer workspace

   eigenTimer.Stop();
   eigenTimer.Report( 0x10, "     Compute symmetric Eigenproblem: %s" );
   return( info );      //return EV routine error message (info=0 means everything fine)
   }

//saves the reconstruction to disk
int EigenSym::saveReconstructionToDisk( StudyData *data, const string& filename )
   {
   ProgressNL( "      Save reconstruction of Kernel: [%s]", FullPath( filename ).c_str() );
   if ( !data->kernelData )
      {
      data->kernelData = new real[data->cIndividuals * data->cIndividuals];
      }
   if ( (this->rankU == this->n) || (useWoodbury == false) )
      {
      real *U2= new real[data->cIndividuals * data->cIndividuals];
      for ( size_t i=0; i<(size_t)this->rankU; ++i )
         {
         for( size_t j=0; j<data->cIndividuals; ++j )
            {
            U2[(i*data->cIndividuals)+j] = this->U[(i*data->cIndividuals)+j] * this->S[i];
            }
         }
      BLAS_INT n_individuals = (BLAS_INT)data->cIndividuals;   //Warning: size_t BLAS_INT conversion
      BLAS_INT n_cols = (BLAS_INT)this->rankU;                 //Warning: size_t BLAS_INT conversion
      BlasLapack::gemm( &char_N, &char_T, &n_individuals, &n_individuals, &n_cols, &real_one, U2, &n_individuals, this->U, &n_individuals, &real_zero, data->kernelData, &n_individuals );
      delete[] U2;
      }
   else if (this->rankU > this->n )
      {
      Fatal( "rankU=%i is larger than n=%i", this->rankU, this->n );
      }
   else
      {
      BLAS_INT n_individuals_int = (BLAS_INT)data->cIndividuals;  //Warning: size_t BLAS_INT conversion
      BLAS_INT n_cols_int = (BLAS_INT)this->rankU;                //Warning: size_t BLAS_INT conversion
      BlasLapack::gemm( &char_N, &char_T, &n_individuals_int, &n_individuals_int, &n_cols_int, &real_one, this->U, &n_individuals_int, this->U, &n_individuals_int, &real_zero,  data->kernelData, &n_individuals_int );
      }
   MatrixIO::writeMatrix2dTXT( filename, "var", "\t", data->individualLabels, data->individualLabels, data->cIndividuals, data->cIndividuals, data->kernelData );
   return 0;
   }

//check, if the kernel is low rank
void EigenSym::treatNonPSD( real *S_, size_t rankU_, bool orderedEvAscending )
   {
   if ( !orderedEvAscending )
      {
      Fatal( "non-ascending ordering is not implemented for EigenSym::treatNonPSD()" );
      }

   if ( S_[0] < 0.0 )
      {
      //The esigenvalues have inverse ordering
      real evMin = S_[0];
      if ( evMin < -MathFunctions::eps_rank )
         {
         Warn( "Kernel is not positive semi definite!" );
         Warn( "The smallest Eigenvalue is %.4e, the largest Eigenvalue is %.4e.", evMin, S_[rankU_-1] );
         }
      else
         {
         Warn(          "The kernel has a few Eigenvalues that are a tiny bit smaller"
             "\n           than zero and are considered to be numerically zero."
             );
         }

      if ( handle_non_psd == SUBTRACT_MIN )
         {
         Warn( "Make the kernel positive semi definite." );
         Warn( "Subtracting the smallest negative eigenvalue (%.4e) from all eigenvalues.", evMin );
         for ( size_t i = 0; i<rankU_; ++i )
            {
            S_[i] -= evMin;
            }
         }
      else if ( handle_non_psd == SET_TO_ZERO )
         {
         Warn( "Make the kernel positive semi definite." );
         Warn( "Setting negative eigenvalues to zero." );
         for ( size_t i = 0; i<rankU_; ++i )
            {
            if( S_[i] < 0.0 )
               {
               S_[i] = 0.0;
               }
            }
         }
      }
   }

void DumpBlasLapackSyevrArgs( const char *pJobz,               // 1 jobz
                              const char *pRange,              // 2 range
                              const char *pUpLo,               // 3 uplo
                              int  *pN,                        // 4 n
                              real *pA,                        // 5 a
                              int  *pLeadingDimensionA,        // 6 lda
                              real *pValLower,                 // 7 vl
                              real *pValUpper,                 // 8 vu
                              int  *pIndexLower,               // 9 il
                              int  *pIndexUpper,               // 10 iu
                              real *pAbsoluteErrorTolerance,   // 11 abstol
                              int  *pEigenValuesFound,         // 12 m
                              real *pSelectedEigenValues,      // 13 w
                              real *pSelectedEigenVectors,     // 14 z
                              int  *pLeadingDimensionZ,        // 15 ldz
                              int  *pDimensionWork,            // 16 isuppz
                              real *pWorkArray,                // 17 work
                        const int  *pLWork,                    // 18 lwork
                              int  *pIWork,                    // 19 iwork
                        const int  *pLIWork,                   // 20 liwork
                              int  *pInfo,                     // 21 info

                              const string& comment
                            //  , bool fOpenAppend
                              )
   {
   static bool fPreviouslyOpened = false;
   if ( !fWriteLogFile )
      {
      return;
      }

   string logFile = MakePath( "", logDir, "FastLmm.SyevrArgs", "Log" );
   string openMode = fPreviouslyOpened ? "a" : "w";
   FILE *pf = fopen( logFile.c_str(), openMode.c_str() );
   if ( pf == nullptr )
      {
      Fatal( "Unable to open log file: [%s]"
           "\n  CRT Error [%d]: %s", logFile.c_str(), errno, strerror( errno ) );
      }
   fPreviouslyOpened = true;

   fprintf( pf, "\ndsyevr() Arguments: [%s]", comment.c_str() );
   fprintf( pf, "\n  char *Jobz: %c", *pJobz  );
   fprintf( pf, "\n  char *Range: %c", *pRange );
   fprintf( pf, "\n  char *UpLo: %c", *pUpLo );
   fprintf( pf, "\n  int  *N: %d", *pN );
   fprintf( pf, "\n  real *A:" );   // array 2d [lda, n]
     DumpRealArray( pf, 4, pA, "A", *pLeadingDimensionA, *pN );
   fprintf( pf, "\n  int  *pLeadingDimensionA: %d:", *pLeadingDimensionA );
   fprintf( pf, "\n  real *pValLower: %.17e", *pValLower );
   fprintf( pf, "\n  real *pValUpper: %.17e", *pValUpper );
   fprintf( pf, "\n  int  *pIndexlower: %d", *pIndexLower );
   fprintf( pf, "\n  int  *pIndexUpper: %d", *pIndexUpper );
   fprintf( pf, "\n  real *pAbsoluteErrorTolerance: %.17E", *pAbsoluteErrorTolerance );
   fprintf( pf, "\n  int  *pEigenValuesFound: %d", *pEigenValuesFound );
   fprintf( pf, "\n  real *pSelectedEigenValues:" );
     DumpRealArray( pf, 4, pSelectedEigenValues, "w", 1, *pN );
   fprintf( pf, "\n  real *pSelectedEigenVectors:" );
     DumpRealArray( pf, 4, pSelectedEigenVectors, "z", *pLeadingDimensionZ, *pEigenValuesFound );
   fprintf( pf, "\n  int  *pLeadingDimensionZ: %d", *pLeadingDimensionZ );
   fprintf( pf, "\n  int  *pDimensionWork: %d", *pDimensionWork );
   fprintf( pf, "\n  real *pWorkArray: %.17E", *pWorkArray );
   fprintf( pf, "\n  int  *pLWork: %d", *pLWork );
   fprintf( pf, "\n  int  *pIWork: %d", *pIWork );
   fprintf( pf, "\n  int  *pLIWork: %d", *pLIWork );
   fprintf( pf, "\n  int  *pInfo: %d", *pInfo );
   fclose( pf );
   }

void EigenSym::Dump( FILE* pf, const string& comment ) const
   {
#if defined( COMPILE_FOR_WINDOWS )     // Windows/VC uses a %Iu specifier for size_t
   const char *szFmt1 = "\n  size_t n: %Iu";
   const char *szFmt2 = "\n  size_t rankU: %Iu";
#else                      // Linux/g++ uses a %zu specifier for size_t
   const char *szFmt1 = "\n  size_t n: %zu";
   const char *szFmt2 = "\n  size_t rankU: %zu";
#endif

   fprintf( pf, "\nEigenSym: [%s]", comment.c_str() );
   fprintf( pf, "\n  string headerfilename: [%s]", headerfilename.c_str() );
   fprintf( pf, "\n  string Ufilename: [%s]", Ufilename.c_str() );
   fprintf( pf, "\n  string Sfilename: [%s]", Sfilename.c_str() );
   fprintf( pf, "\n  string id: [%s]", id.c_str() );
   fprintf( pf, szFmt1, n );
   fprintf( pf, szFmt2, rankU );
   DumpRealArray( pf, 2, U, "U", n, rankU );
   DumpRealArray( pf, 2, S, "S", rankU, 1 );
   fprintf( pf, "\n" );
   }
