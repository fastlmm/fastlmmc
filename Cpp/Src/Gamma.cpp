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

// This is gammainc(x,a,'upper') in Matlab
/// <summary>q-gamma function</summary>
/// <remarks>This gammainc(x,a,'upper')in Matlab</remarks>
real Gamma::gammaQ(real x, real a)
{
   if (x < a + (real)1.0)
   {
         return (real)1.0 - Gamma::gammaIncLower(x, a);
   }
   else
   {
         return Gamma::gammaIncUpper(x, a);
   }
}

// P(a,x) using series expansion
// See 6.5.4 and 6.5.29:
// y(a,x) = x^a * exp(-x) * sum (Gamma(a) * z^n / Gamma(a+n+1))
// where sum goes from 0 to infinity
// We use 6.1.15 & 6.1.16 to compute Gamma(a)/Gamma(a+n)
// Gamma(a) = (a-1)!
// Gamma(a+n) = (n-1+a)*(n-2+a)...(1+a)*a! = (n-1+a)*(n-2+a)...(1+a)*a*(a-1)! 
// Gamma(a)/Gamma(a+n) = (n-1+a)*(n-2+a)...(1+a)*a


// This converges quickly for x < a + 1
/// <summary>Lower incomplete gamma function</summary>
real Gamma::gammaIncLower(real x, real a)
{
   real sum = (real)1.0 / a;
   real frac = sum;
   real an = a;

   for (int i = 1; i <= MAXITER; ++i)
   {
         an++;
         frac = x * frac / an;
         sum += frac;

         // Stopping criterion: sum won't change
         if (::abs(frac) < ::abs(sum * EPS))
            break;
   }

   return sum * exp(-x - MathFunctions::logGamma(a) + a * log(x));
}

// Q(a,x) using continued fractions
// This converges quickly for x > a + 1
// See 6.5.31 in ref [2]
// Using the modified Lentz's method, described in section 5.2 of 
// Numerical Recipes
/// <summary>Upper incomplete gamma function</summary>
real Gamma::gammaIncUpper(real x, real a)
{
   real f = TINY;
   real c = f;
   real d = 0;
   real aj = 1;
   real bj = x + 1 - a;
            
   for (int i = 1; i <= MAXITER; ++i)
   {
         d = bj + aj * d;
         if (::abs(d) < TINY) d = TINY;
         c = bj + aj / c;
         if (::abs(c) < TINY) c = TINY;
         d = 1 / d;
         real delta = c * d;
         f *= delta;

         if (::abs(delta - 1) < EPS)
            break;

         bj += 2;
         aj = -i * (i - a);
   }

   return f * exp(-x - MathFunctions::logGamma(a) + a * log(x));
}

/// <summary>Cumulative distribution function</summary>
/// <param name="x">Value at which to compute the cdf</param>
/// <param name="k">Shape parameter</param>
/// <param name="theta">Scale parameter</param>
real Gamma::cdf(real x, real k, real theta)
{
   if (x < 0){
      printf("Gamma::cdf argument x=%f must be >= 0",x);
      throw(1);
   }
   if (k <= 0 || theta <= 0){
      printf("Gamma::cdf arguments k=%f and theta=%f must be > 0",k,theta);
      throw(1);
   }
   return Gamma::gammaP(x / theta, k);
}

// Incomplete gamma functions

/// <summary>p-gamma function</summary>
/// <remarks>This gammainc(x,a,'lower') or gammainc(x,a) in Matlab</remarks>
real Gamma::gammaP(real x, real a)
{
   if (x < a + 1.0)
   {
         return Gamma::gammaIncLower(x, a);
   }
   else
   {
         return (real)1.0 - Gamma::gammaIncUpper(x, a);
   }
}

/// <summary>Inverse cumulative distribution function</summary>
/// <param name="p">Probability at which to compute the inverse cdf</param>
/// <param name="k">Shape parameter</param>
/// <param name="theta">Scale parameter</param>
real Gamma::inv(real p, real k, real theta)
{
   // We use an iterative search to find x - this is Excel's algorithm (I think)

   if (p == 0) return 0.0;
   if (p == 1) return std::numeric_limits<real>::infinity( );

   if (p < 0 || p > 1){
         printf("Gamma::inv parameter p must be between 0 and 1.");
         throw(1);
   }
   if (k <= 0 || theta <= 0){
         printf("Gamma::inv parameters a and b must be greater than 0.");
         throw(1);
   }

   real guesslo = 0.0;
   real ab = k * theta;
   // What's a good initial guess when a*b < 1 ???
   real guesshi = ab > (real)1.0 ? (real)10.0 * ab : (real)10.0;

   real guess = 0.0;
   real pguess;
   for (int i = 0; i < 400; ++i)
   {
         guess = (guesslo + guesshi) * (real)0.5;
         pguess = Gamma::cdf(guess, k, theta);

         if (::abs(p - pguess) < (real)1e-16)
            break;
                
         if (pguess > p)
            guesshi = guess;
         else
            guesslo = guess;
   }

   return guess;
}

        /// <summary>Computes the psi (digamma) function</summary>
double Gamma::Psi(real k)
{
   real dig = 0;

   if (k < 8)
   {
         dig = Gamma::Psi(k + 1);
         dig -= 1.0 / k;
   }
   else
   {
         real ksq = k * k;
         real k3 = ksq * k;
         real k4 = ksq * ksq;
         real k5 = k4 * k;
         real k6 = ksq * k4;
         real k7 = k6 * k;

         dig = log(k) - 1.0 / (2.0 * k) - 1.0 / (12.0 * ksq) + 1.0 / (120.0 * k4) - 1.0 / (252.0 * k6);
   }
   return dig;
}
