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

/// <summary>Brent's method for minimizing a 1d function</summary>
#include "FastLmmC.h"

// Machine eps
real BrentC::MACHEPS = (real)2.2204460492503131e-016;
real BrentC::MACHEPS_SQRT = sqrt(MACHEPS);
real BrentC::c = ((real)3.0 - sqrt((real)5.0)) / (real)2.0;

   
/// <summary>Minimize the function over the interval [a, b]</summary>
/// <param name="f">Function to minimize</param>
/// <param name="a">Left side of the bracket</param>
/// <param name="b">Right side of the bracket</param>
/// <param name="eps">Stopping tolerance</param>
/// <param name="funcx">Function evaluated at the minimum</param>
/// <param name="numiter">number of function evaluations</param>
/// <param name="maxIter">maximum number of function evaluations allowed</param>
/// <param name="quiet">print out function evaluations?</param>
/// <returns>Point that minimizes the function</returns>
/// <remarks>This implements the algorithm from Brent's book, "Algorithms for Minimization without Derivatives"</remarks>
real BrentC::minimize( BrentFunctor &f, real a, real b, real eps, real &funcx, size_t &numiter, size_t maxIter, bool quiet )
   {
   if ( a >= b )
      {
      printf( "Exception: a must be < b" );
      throw( 1 );
      }

   real x = a + c * (b - a);
   real v = x;
   real w = x;
   real e = 0;

   real fx = f(x);
   real fv = fx;
   real fw = fx;

   numiter = 0;

   while ( true )
      {
      real m = (real)0.5 * (a + b);
      real tol = MACHEPS_SQRT * ::abs(x) + eps;
      real tol2 = (real)2 * tol;
                
      // Check the stopping criterion
      if ( ::abs(x - m) <= (tol2 - (0.5 * (b - a))) )
         {
         break;
         }

      // Stop if we've exceeded the maximum number of iterations
      if ( ++numiter > maxIter )
         {
         printf( "Exception: Exceeded maximum number of iterations." );
         throw( 2 );
         }

      real p = 0.0, q = 0.0, r = 0.0;
      real d = 0.0;
      real u = 0.0;

      if ( ::abs( e ) > tol )
         {
         // Fit parabola
         r = (x - w) * (fx - fv);
         q = (x - v) * (fx - fw);
         p = (x - v) * q - (x - w)*r;
         q = (real)2.0 * (q - r);
         if ( q > (real)0.0 )
            {
            p = -p;
            }
         else
            {
            q = -q;
            }
         r = e;
         e = d;
         }

      if ( ( ::abs(p) < ::abs( (real)0.5 * q * r ) )
        && (p < (q * ( a - x )) )
        && (p < (q * ( b - x )) )
         )
         {
         // Parabolic interpolation step
         d = p / q;
         u = x + d;
         // f must not be evaluated too close to a or b
         if ( ((u - a) < tol2) || ((b - u) < tol2) )
            {
            d = (x < m) ? tol : -tol;
            }
         }
      else
         {
         // Golden section step
         e = (x < m) ? (b - x) : (a - x);
         d = c * e;
         }

      // f must not be evaluated too close to x
      if ( ::abs( d ) >= tol )
         {
         u = x + d;
         }
      else if ( d > 0.0 )
         {
         u = x + tol;
         }
      else
         {
         u = x - tol;
         }

      real funcu = f( u );
                
      // Update
      if ( funcu <= fx )
         {
         if ( u < x )
            {
            b = x;
            }
         else
            {
            a = x;
            }

         v = w; fv = fw;
         w = x; fw = fx;
         x = u; fx = funcu;
         }
      else
         {
         if ( u < x )
            {
            a = u;
            }
         else
            {
            b = u;
            }

         if ( (funcu <= fw) || (w == x) )
            {
            v = w;
            fv = fw;
            w = u;
            fw = funcu;
            }
         else if ( (funcu <= fv) || (v == x) || (v == w) )
            {
            v = u;
            fv = funcu;
            }
         }

      if ( !quiet )
         {
#if defined( COMPILE_FOR_WINDOWS )    // Windows/VC uses a %Iu specifier for size_t
         const char *szFmt1 = "Iteration %Iu, min_x = %f, f(min_x) = %f, ";
#else                      // Linux/g++ uses a %zu specifier for size_t
         const char *szFmt1 = "Iteration %zu, min_x = %f, f(min_x) = %f, ";
#endif
         printf( szFmt1, numiter, x, fx );
         }
      }

   funcx = f( x );
   return( x );
   }
