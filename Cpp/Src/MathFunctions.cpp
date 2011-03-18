#include "FastLmmC.h"

//const double _PI = 2.0*acos(0.0);          // now defined in FastLmmC.h
const double _halflog2pi=0.5*log(2.0*_PI);
const double coeffsForLogGamma[] = { 12.0, -360.0, 1260.0, -1680.0, 1188.0 };

// Gamma and LogGamma

// Use GalenA's implementation of LogGamma - it's faster!
/// <summary>Returns the log of the gamma function</summary>
/// <param name="x">Argument of function</param>
/// <returns>Log Gamma(x)</returns>
/// <remarks>Accurate to eight digits for all x.</remarks>
double MathFunctions::logGamma(double x)
{
   if (x <= 0.0){
      printf("LogGamma arg=%f must be > 0.",x);
      throw(1);
   }

   double res = 0.0;
   if (x < 6.0)
   {
         int toAdd = (int)floor(7 - x);
         double v2 = 1.0;
         for (int i = 0; i < toAdd; i++)
         {
            v2 *= (x + i);
         }
         res = -log(v2);
         x += toAdd;
   }
   x -= 1.0;

   res += _halflog2pi + (x + 0.5) * log(x) - x;

   // correction terms
   double xSquared = x * x;
   double pow = x;
   for (int i=0; i<5; ++i)   //the length of the coefficient array is 5.
   {
         double newRes = res + 1.0 / (coeffsForLogGamma[i] * pow);
         if (newRes == res)
         {
            return res;
         }
         res = newRes;
         pow *= xSquared;
   }

   return res;
}

//solve the quadratic least squares system of a positive semi-definite symmetric matrix X in lower storage. On exit X holds the inverse of X and y holds the solution.
//returns the rank of X
int MathFunctions::leastSquares2by2SymLower(double *X, double *y){
   double determinant=X[0]*X[3]-X[1]*X[1];
   if ( determinant > 0.0 ){//full rank case. just solve
      //inverse of X
      determinant=1.0/determinant;
      X[1]*=-1.0*determinant;
      X[2]=X[0];
      X[0]=determinant*X[3];
      X[3]=X[2]*determinant;
      X[2]=y[0];
      //solve y
      y[0]=X[0]*y[0]+X[1]*y[1];
      y[1]=X[1]*y[0]+X[3]*y[1];
      X[2]=X[1];
      return 2;//rank2
   }else{// see http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html for derivations of Eigenvalues and eigenvectors
      double trace=X[0]+X[3]; //In this case the trace equals the nonzero eigenvalue L_1 in the derivation
      if (trace == 0.0){
         y[0]= 0.0;
         y[1]= 0.0;
         return 0;
      }
      else{
         trace=1.0/trace;  // the inverse of the trace
         if (X[1] == 0.0){
            y[0]*=trace;
            y[1]=0.0;
         }else{
            trace*=X[0];
            X[2]=y[0]; //cache y[0]
            y[0]=(X[0]*X[2]+X[1]*y[1])*trace;
            y[1]=(X[1]*X[2]+X[0]*y[1])*trace;
         }
         return 1;//rank1
      }
   }
}
