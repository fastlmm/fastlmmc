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

//const real _PI = 2.0*acos(0.0);          // now defined in FastLmmC.h
const real _halflog2pi=(real)0.5*log((real)2.0*_PI);
const real coeffsForLogGamma[] = { 12.0, -360.0, 1260.0, -1680.0, 1188.0 };

const real MathFunctions::eps_rank=(real)3E-8;

// Gamma and LogGamma

// Use GalenA's implementation of LogGamma - it's faster!
/// <summary>Returns the log of the gamma function</summary>
/// <param name="x">Argument of function</param>
/// <returns>Log Gamma(x)</returns>
/// <remarks>Accurate to eight digits for all x.</remarks>
real MathFunctions::logGamma(real x)
{
   if (x <= (real)0.0){
      printf("LogGamma arg=%f must be > 0.",x);
      throw(1);
   }

   real res = (real)0.0;
   if (x < (real)6.0)
   {
         int toAdd = (int)floor(7 - x);
         real v2 = (real)1.0;
         for (int i = 0; i < toAdd; i++)
         {
            v2 *= (x + i);
         }
         res = -log(v2);
         x += toAdd;
   }
   x -= (real)1.0;

   res += _halflog2pi + (x + (real)0.5) * log(x) - x;

   // correction terms
   real xSquared = x * x;
   real pow = x;
   for (int i=0; i<5; ++i)   //the length of the coefficient array is 5.
   {
         real newRes = res + (real)1.0 / (coeffsForLogGamma[i] * pow);
         if (newRes == res)
         {
            return res;
         }
         res = newRes;
         pow *= xSquared;
   }

   return res;
}

// Beta and LogBeta
/// <summary>Computes the log beta function</summary>
double MathFunctions::LogBeta(real x, real y)
{
   if (x <= 0.0 || y <= 0.0){
      printf("LogBeta args must be > 0.");
      throw(1);
   }
   return MathFunctions::logGamma(x) + MathFunctions::logGamma(y) - MathFunctions::logGamma(x + y);
}

int MathFunctions::randperm(const size_t n, const real *in, real *out){
   if (n>0){
      out[0] = in[0];
      size_t j;
      for (size_t i=1; i<n; ++i){
         j = rand() % (i+1);
         out[i] = out[j];
         out[j] = in[i];
      }
   }
   return 0;
}

int MathFunctions::randperm(const size_t n, size_t *randperm){
   if (n>0){
      randperm[0] = 0;
      size_t j;
      for (size_t i=1; i<n; ++i){
         j = rand() % (i+1);
         randperm[i] = randperm[j];
         randperm[j] = i;
      }
   }
   return 0;
}
