#if !defined( BrentC_h )
#define BrentC_h

//abstract base class for Brent functors:
class BrentFunctor{
public:
   //the operator (function evaluation) takes a double as argumant and returns a double
   virtual double operator()(const double x)=0;
};

/// <summary>Brent's method for minimizing a 1d function</summary>
class BrentC{
   // Machine eps
   static double MACHEPS;
   static double MACHEPS_SQRT;
   static double c;

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
   static double minimize(BrentFunctor &f, double a, double b, double eps, double &funcx, int &numiter, int maxIter, bool quiet);
};

#endif      //BrentC_h
