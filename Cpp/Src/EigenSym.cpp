#include "FastLmmC.h"

const string EigenSym::headerfilename="IndividualID_header.txt";
const string EigenSym::Ufilename="U";
const string EigenSym::Sfilename="S";

EigenSym::EigenSym(void)
{
   this->U=0;
   this->S=0;
}

EigenSym::~EigenSym(void)
{
   delete[] U;
   delete[] S;
}

//saves an EigenSym object to disk in TAB delimited format
int EigenSym::saveToDiskTAB(string Filename){
   printf("Not implemented.");
   throw(1);
   return 0;
}                        

//saves an EigenSym object to disk in BINARY delimited format
int EigenSym::saveToDiskBINARY(string dirname){
   if(dirname.length()==0){dirname=".";};

   MatrixIO::createDirectory(dirname);
   
   string filenameHeader=dirname;
   filenameHeader.append(MatrixIO::OS_directoryDelim);
   filenameHeader.append(EigenSym::headerfilename);
   MatrixIO::writeHeaderFile(filenameHeader, this->columnLabelOrder);

   string filenameU=dirname;
   filenameU.append(MatrixIO::OS_directoryDelim);
   filenameU.append(EigenSym::Ufilename);
   
   string filenameS=dirname;
   filenameS.append(MatrixIO::OS_directoryDelim);
   filenameS.append(EigenSym::Sfilename);
   
   MatrixIO::writeMatrix2dBIN(filenameU,id,n,rankU,U);
   MatrixIO::writeMatrix2dBIN(filenameS,id,rankU,1,S);
   return 0;
}                     

//reads an EigenSym Object from disk in TAB delim format
int EigenSym::readFromDiskTAB(string dirname){
   printf("Not implemented.");
   throw(1);
   return 0;
}           

//reads an EigenSym object from disk in BINARY format
int EigenSym::readFromDiskBINARY(string dirname){
   if(dirname.length()==0){dirname=".";};

   string filenameHeader=dirname;
   filenameHeader.append(MatrixIO::OS_directoryDelim);
   filenameHeader.append(EigenSym::headerfilename);
   this->columnLabelOrder = MatrixIO::readHeaderFile(filenameHeader);

   //read U
   string filenameU=dirname;
   filenameU.append(MatrixIO::OS_directoryDelim);
   filenameU.append(EigenSym::Ufilename);
   
   MatrixIO::querySizeMatrixBIN( filenameU, &id, &n, &rankU );
   
   if (this->U!=0){delete[] U;}
   this->U = new double[n*rankU];
   
   MatrixIO::readBinary2dBIN( filenameU, id, n, rankU, U );

   //read S
   string filenameS=dirname;
   filenameS.append(MatrixIO::OS_directoryDelim);
   filenameS.append(EigenSym::Sfilename);
   
   string idS;
   int nS=0;
   int mS=0;
   MatrixIO::querySizeMatrixBIN( filenameS, &idS, &nS, &mS );

   if(id!=idS){
      Fatal("The ids of U and S do not match up. idU:%s idS:%s",id.c_str(),idS.c_str());
   }
   if(nS!=rankU){
      Fatal("The first dimension of matrix S has to be rankU=%i. Found %i",rankU,nS);
   }
   if(mS!=1){
      Fatal("The second dimension of matrix S has to be 1. Found %i",mS);
   }
   if (this->columnLabelOrder.size()!=this->n){
      Fatal("The number of columLabels=%i has to be n=%i.",this->columnLabelOrder.size(),n);
   }

   if (this->S!=0){delete[] S;}
   this->S = new double[rankU];
   
   MatrixIO::readBinary2dBIN( filenameS, idS, rankU, mS, S );
   
   return 0;
}        

//X -[n*m] matrix, with K=X*X^T  --->K[n*n];  rankU=min(n,m)
int EigenSym::set_covFromX(string id, double* X, int n, int m, vector<string> columnLabelOrder){
   if (columnLabelOrder.size()!=n){
      printf("the number of columnLabels=%i does not math the number of individuals n=%i.",columnLabelOrder.size(),n);
      throw(1);
   }
   this->n=n;
   this->rankU=min(n,m);     //The number of singular vectors (rankU) is the minimum of the two dimensions of the SNP sample matrix.
      if (this->U!=0){
         delete[] U;
         this->U=0;
      }
      if (this->S!=0){
         delete[] this->S;
         this->S=0;
      }
   if(rankU==n){  //full rank case
      printf("Creating covariance Matrix from SNP sample...");
      double *K=new double[n*n];
      BlasLapack::syrk(&char_L,&char_N,&n,&m,&double_one, X, &n, &double_zero,  K, &n);
      printf("done.\n");
      this->set_fromKernel(id, K,n,rankU,columnLabelOrder);
      delete[] K;
   }else{   //low rank case
      printf("Computing the low rank Eigenvectors of the Kernelmatrix...");
      this->columnLabelOrder=columnLabelOrder;
      this->id=id;
      this->U=new double[n*rankU];            //the left singular vectors of the SNP sample matrix are the Eigenvectors of the Covariance Kernel matrix.
      this->S=new double[rankU];
      double *VT=new double[rankU*m];           //the right singular vectors of the SNP sample matrix
      EigenSym::SVD_economy(n, m, X, U,S, VT);  //The economy SVD returns the first rankU Singular vectors and Singular values
      delete[] VT;                                 //We only need the left singular vectors, so this part can go.

      for (int i=0;i<rankU;++i){
         S[i]=S[i]*S[i];      //the Eigenvalues of the Kernel are the squared Singular values of the SNP-matrix
      }
      this->id=id;
      printf("done.\n");
   }
   return 0;
}  

EigenSym::EigenSym(string id, double *U, double *S, int n, int rankU, vector<string> columnLabelOrder){
   if (columnLabelOrder.size()!=n){
      printf("the number of columnLabels=%i does not math the number of individuals n=%i.",columnLabelOrder.size(),n);
      throw(1);
   }
   this->columnLabelOrder=columnLabelOrder;
   this->id=id;
   if(rankU>n){
      printf("rankU=%i>n%i. This is a bug!",rankU,n);
      throw(1);
   }
   this->rankU=rankU;
   this->n=n;
   this->U=U;
   this->S=S;
}


int EigenSym::set_fromKernel(string id, double *K, int n, int rankU, vector<string> columnLabelOrder){
   if (columnLabelOrder.size()!=n){
      printf("the number of columnLabels=%i does not math the number of individuals n=%i.",columnLabelOrder.size(),n);
      throw(1);
   }
   this->id=id;
   this->columnLabelOrder=columnLabelOrder;
   printf("Computing the Eigenvectors of the Kernelmatrix...");
   if (this->U!=0){ delete[] this->U; }
   if (this->S!=0){ delete[] this->S; }
   this->n=n;
   this->rankU=min(n,rankU);
   if (this->n==this->rankU){
      this->U=new double[n*n];
      this->S=new double[n];

	   EigenSym::EV_symmetric(this->n,K,this->U,this->S);
   }
   else{
      printf("not implemented.");
      throw(1);
   }
   printf("done.\n");
   return 0;
}

//computes the economy sized SVD of an [n by d] matrix A=U*S*VT. Note that A will be overwritten.
//eigenvalues are returned in descending order.
int EigenSym::SVD_economy(int n, int d, double *A, double *U,double *S, double *VT){
	int nd=n*d;
	int dmin=min(n,d);
	double work_svd_query;
	int lwork_svd =-1;
	int *iwork_svd=new int[8*dmin];           //deleted
	int info_svd;
	BlasLapack::gesdd(&char_S, &n, &d, A, &n, S, U, &n, VT, &dmin, &work_svd_query,  &lwork_svd, iwork_svd, &info_svd ); //workspace query
	lwork_svd=(int)work_svd_query;
	double *work_svd=new double[lwork_svd];   //deleted
	BlasLapack::gesdd(&char_S, &n, &d, A, &n, S, U, &n, VT, &dmin, work_svd,  &lwork_svd, iwork_svd, &info_svd );        //execute the SVD
	delete[] work_svd;
	delete[] iwork_svd;
	return info_svd;
}

//computes the Eigenvectors of an [n by n] symmetric matrix A=U*S*U.T. Note that A will be overwritten.
//eigenvalues are returned in ascending order.
int EigenSym::EV_symmetric(int n,double *A, double *U,double *S){
	double vl,vu;        //the range of eigenvalues (not referenced)
	int il,iu;           //the indices of eigenvalues (not referenced)
	double abstol=3E-8;  //tolerance of the EV algorithm
	int m;               //out: the number of eigenvalues found, should be n.
	int* isuppz=new int[2*n];  //second integer workspace (deleted)
	double workdummy;    //will hold the optimal size of the workspace
	int iworkdummy;      //will hold the optimal size of the int workspace
	int info;            //error message from Eigenvector routine
	BlasLapack::syevr( &char_V, &char_A, &char_L, &n, A, &n, &vl, &vu, &il, &iu, &abstol, &m, S, U, &n, isuppz, &workdummy, &int_minusone,  &iworkdummy, &int_minusone, &info );//optimal workspace request
	int lwork=(int)workdummy;        //optimal size of the workspace
	int liwork=iworkdummy;           //optimal size of the integer workspace
	double* work=new double[lwork];  //workspace of optimal size (deleted)
	int* iwork= new int[liwork];     //integer workspace of optimal size (deleted)
	BlasLapack::syevr( &char_V, &char_A, &char_L, &n, A, &n, &vl, &vu, &il, &iu, &abstol, &m, S, U, &n, isuppz, work, &lwork,  iwork, &liwork, &info );//Eigenvectors and Eigenvalues are computed
	delete[] work;       //delete workspace
	delete[] iwork;      //delete integer workspace
   delete[] isuppz;     //delete second integer workspace
	return info;         //return EV routine error message (info=0 means evrth fine)
}
