#if !defined(MatrixIO_h)
#define MatrixIO_h

class MatrixIO
   {
public:
   const static string OS_directoryDelim;
   static int writeMatrix2dBIN(string filename, string matrixID, int n, int m, double *A);
   static int querySizeMatrixBIN(string filename, string *matrixID, int *n, int *m);
   static int readBinary2dBIN(string filename, string matrixID, int n, int m, double *A);
   static vector<string> readHeaderFile(string filename);
   static int writeHeaderFile(string filename, vector<string> header);
   static int createDirectory(string dirname);
   };

#endif      //MatriIO_h
