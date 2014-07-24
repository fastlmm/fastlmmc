#if !defined( BrentC_h )
#define BrentC_h
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

//abstract base class for Brent functors:
class BrentFunctor{
public:
   //the operator (function evaluation) takes a real as argumant and returns a real
   virtual real operator()(const real x)=0;
};

/// <summary>Brent's method for minimizing a 1d function</summary>
class BrentC{
   // Machine eps
   static real MACHEPS;
   static real MACHEPS_SQRT;
   static real c;

public:
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
   static real minimize(BrentFunctor &f, real a, real b, real eps, real &funcx, size_t &numiter, size_t maxIter, bool quiet);
};

#endif      //BrentC_h
