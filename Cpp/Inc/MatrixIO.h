#if !defined(MatrixIO_h)
#define MatrixIO_h
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

class MatrixIO
   {
public:
   const static string OS_directoryDelim;
   static int writeMatrix2dBIN(const string& filename, const string& matrixID, size_t n, size_t m, real *A);
   static int querySizeMatrixBIN(const string& filename, string *matrixID, size_t *n, size_t *m);
   static int readBinary2dBIN(const string& filename, const string& matrixID, size_t n, size_t m, real *A);
   static int writeMatrix2dTXT(const string& filename, const string& matrixID, const string& separator, vector<string>& rowhdr, vector<string>& colhdr, size_t rows, size_t cols, real *A);
   static vector<string> readHeaderFile(const string& filename);
   static int writeHeaderFile(const string& filename, vector<string>& header);
   static int createDirectory(const string& dirname);
   };

#endif      //MatriIO_h
