#if !defined( EigenSym_h )
#define EigenSym_h

class EigenSym
{
   static const string headerfilename;
   static const string Ufilename;
   static const string Sfilename;
public:
   int n;      //Dimensionality of the original symmetric matrix ( K \in \R^[n*n] )
   int rankU;  //the rank of the reconstruction U*diag(S)*U^T
   double *U;  // [n*rankU] matrix holding the Eigenvectors
   double *S;  // [rankU]-vector holding the Eigenvalues   K=U*diag(S)*U^T
   string id;
   vector<string> columnLabelOrder;

   EigenSym(void);
   EigenSym(string id, double *U, double *S, int n, int rankU, vector<string> columnLabelOrder);
   int set_fromKernel(string id, double *K, int n, int rankU, vector<string> columnLabelOrder);
   int set_covFromX(string id, double* X, int n, int m, vector<string> columnLabelOrder);      //X -[n*m] matrix, with K=X*X^T
   int saveToDiskTAB(string dirname);             //saves an EigenSym object to disk in TAB delimited format
   int saveToDiskBINARY(string dirname);          //saves an EigenSym object to disk in BINARY delimited format
   int readFromDiskTAB(string dirname);           //reads an EigenSym Object from disk in TAB delim format
   int readFromDiskBINARY(string dirname);        //reads an EigenSym object from disk in BINARY format
   

   static int SVD_economy(int n, int d, double *A, double *U,double *S, double *VT);
   //computes the Eigenvectors of an [n by n] symmetric matrix A=U*S*U.T. Note that A will be overwritten.
   //eigenvalues are returned in ascending order.
   static int EV_symmetric(int n,double *A, double *U,double *S);
   ~EigenSym(void);
};

#endif      //EigenSym_h
