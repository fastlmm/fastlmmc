
#include "FastLmmC.h"

// Beta distribution
// See: http://en.wikipedia.org/wiki/Beta_distribution
    /// <summary>Beta distribution</summary>

        /// <summary>Cumulative distribution function</summary>
        /// <param name="x">Value at which to compute the cdf</param>
        /// <param name="a">Shape parameter (alpha)</param>
        /// <param name="b">Shape parameter (beta)</param>
real Beta::Cdf(real x, real a, real b){
   if (a <= 0 || b <= 0){
      printf("Beta.Cdf parameters, a and b, must be > 0");
      throw(1);
   }
   if (x < 0 || x > 1){
      printf("Beta.Cdf parameter x must be between 0 and 1");
      throw(1);
   }
   return Beta::BetaInc(x, a, b);
}
        
/// <summary>Probability distribution function</summary>
/// <param name="x">Value at which to compute the pdf</param>
/// <param name="a">Shape parameter (alpha)</param>
/// <param name="b">Shape parameter (beta)</param>
real Beta::Pdf(real x, real a, real b){
   if (a <= 0 || b <= 0){
      printf("Beta.Pdf parameters, a and b, must be > 0");
      throw(1);
   }

   if (x > 1) return 0;
   if (x < 0) return 0;

   real lnb = MathFunctions::LogBeta(a, b);
   return exp((a - 1) * log(x) + (b - 1) * log(1 - x) - lnb);
}



/// <summary>Mean and variance</summary>
/// <param name="a">Shape parameter (alpha)</param>
/// <param name="b">Shape parameter (beta)</param>
/// <param name="var">Output: variance</param>
/// <returns>Mean</returns>
real Beta::Stats(real a, real b, real *var)
{
   if (a <= 0 || b <= 0){
         printf("Beta.Stats parameters, a and b, must be > 0");
         throw(1);
   }
   real mean = a / (a + b);
   real ab = a + b;
   *var = (a * b) / (ab * ab * (ab + 1));
   return mean;
}

/// <summary>Inverse beta cumulative probability distribution</summary>
/// <param name="p">Probability at which to compute the inverse cdf</param>
/// <param name="a">Shape parameter (alpha)</param>
/// <param name="b">Shape parameter (beta)</param>
real Beta::Inv(real p, real a, real b)
{
   // We use an iterative search to find x
   // I'm sure there are better ways of doing this

   if (p < 0 || p > 1){
      printf("Beta.Inv parameter p must be between 0 and 1.");
      throw(1);
   }
   if (a <= 0 || b <= 0){
         printf("Beta.Inv parameters a and b must be greater than 0.");
         throw(1);
   }
   if (p == 0) return 0;
   if (p == 1) return 1;

   // We know that x is between 0 and 1, guess in log-odds space
   real guesslo = -709.0;
   real guesshi = 36.0;

   real guess = 0.0;
   real x = 0.0;
   for (int i = 0; i < MAXITER; ++i)
   {
         // Make a guess half way between the lo and hi
         guess = (guesslo + guesshi) / 2;
         x = 1.0 / (1.0 + exp(-guess));
         if (BetaInc(x, a, b) > p)
            guesshi = guess;    // we guessed too high
         else
            guesslo = guess;    // we guessed too low

         // Convergence!
         if (guesshi - guesslo < 1e-10)
            break;
   }

   return x;
}

real Beta::BetaInc(real x, real a , real b)
{
   if (x == 0 || x == 1)
   {
         return x;
   }
   else
   {
         real c = exp(MathFunctions::logGamma(a + b) - MathFunctions::logGamma(a) - MathFunctions::logGamma(b) +
                  a * log(x) + b * log(1 - x));
         real p;

         if (x < (a+1)/(a+b+2))
         {
            real cf = BetaCF(x, a, b);
            p= c * cf/a;
         }
         else
         {
            // Use symmetry relation
            real cf = BetaCF(1 - x, b, a);
            p = 1 - c* cf / b;
         }

         if (p >= 1.0) p = 1.0;
         return p;
   }
}


// Using the modified Lentz's method, described in section 5.2 of 
// Numerical Recipes
real Beta::BetaCF(real x, real a, real b)
{
   real ap1 = a + 1;
   real am1 = a - 1;
   real m2 = 2;
   real ab = a + b;

   real d = 1 - ab * x / ap1;
   if (::abs(d) < TINY) d = TINY;
   real c = 1;
   d = 1.0 / d;
   real cf = d;

   for (int m = 1; m <= MAXITER; ++m)
   {
         real aj = m * (b - m) * x / ((am1 + m2) * (a + m2));
         d = 1.0 + aj * d;
         if (::abs(d) < TINY) d = TINY;
         c = 1 + aj / c;
         if (::abs(c) < TINY) c = TINY;
         d = 1 / d;
         real delta = c * d;
         cf *= delta;

         aj = -(a + m) * (ab + m) * x/((ap1 + m2) * (a + m2));
         d = 1.0 + aj * d;
         if (::abs(d) < TINY) d = TINY;
         c = 1 + aj / c;
         if (::abs(c) < TINY) c = TINY;
         d = 1 / d;
         delta = c * d;
         cf *= delta;

         if (::abs(delta - 1) < EPS)
            break;

         m2 += 2;
   }

   return cf;
}

real Beta::Eval(real y, real lo, real hi, real s)
{
   real a = lo;
   real b = hi;
   real fa = y - Gamma::Psi(a) + Gamma::Psi(a + s);
   if (fa == 0) return a;
   real fb = y - Gamma::Psi(b) + Gamma::Psi(b + s);
   if (fb == 0) return b;

   real c = a;
   real fc = fa;

   real tol = 1e-10;

   for (int i = 0; i < MAXITER; ++i)
   {
         real diff = b - a;
         if (::abs(fc) < ::abs(fb))
         {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
         }

         real tolact = 2.0 * EPS * ::abs(b) + tol / 2.0;
         real nextstep = (c - b) / 2;

         if (::abs(nextstep) <= tolact || fb == 0.0)
         {
            return b;
         }

         if (::abs(diff) >= tolact && ::abs(fa) > ::abs(fb))
         {
            real p, q;
            real cb = c-b;
            if (a == c)
            {
               real t1 = fa/fb;
               p = cb*t1;
               q = 1.0 - t1;
            }
            else
            {
               q = fa/fc; 
               real t1 = fb/fc;
               real t2 = fb/fa;
               p = t2*(cb*q*(q-t1) - (b-a)*(t1-1.0));
               q = (q-1.0)*(t1-1.0)*(t2-1.0);
            }

            if (p > 0)
               q = -q;
            else
               p = -p;

            if (p < (0.75*cb*q-::abs(tolact*q)/2) && p < ::abs(diff*q/2))
               nextstep = p/q;

         }

         if (::abs(nextstep) < tolact)
         {
            if (nextstep > 0)
               nextstep = tolact;
            else
               nextstep = -tolact;
         }
                
         // Previous approximation
         a = b;
         fa = fb;
         b += nextstep;
         fb = y - Gamma::Psi(b) + Gamma::Psi(b + s);
         if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0))
         {
            c = a;
            fc = fa;
         }

   }

   return b;
}

/*
// Test case: [.9,.1,.3,.23]
/// <summary>Maximum likely estimate</summary>
/// <param name="dataIn">Data to fit</param>
/// <param name="bhat">Shape estimate (beta)</param>
/// <returns>Shape estimate (alpha)</returns>
real Beta::MLE(real *dataIn, real *bhat, int data_size)
{
   real *data = dataIn;
   if (data_size < 2){
         printf("Beta.MLE parameter, data, must be of length > 1");
         throw(1);
   }
   real meanlnp = 0;
   real meanlp1 = 0;

   for (int i = 0; i < data_size; ++i)
   {
         if (data[i] <= 0 || data[i] >= 1){
            printf("Data for Beta.MLE must be > 0 and < 1.");
            throw(1);
         }
         meanlnp += log(data[i]);
         meanlp1 += log(1 - data[i]);
   }

   real cnt = (real)data_size;
   meanlnp /= cnt;
   meanlp1 /= cnt;

   real p = data.Mean();
   real p1 = 1 - p;
   real var = data.Var();
   real pp = p * p1 / var - 1;

   // Starting guesses
   real ahat = pp * p;
   *bhat = pp * p1;

   real toler = 1.0e-6;
            
   real alo = 0.001;
   real ahi = ahat * 2;
   real blo = 0.001;
   real bhi = *bhat * 2;
            
   for (int i = 0; i < MAXITER; ++i)
   {
         real ahatnew = Eval(meanlnp, alo, ahi, *bhat);
         real bhatnew = Eval(meanlp1, blo, bhi, ahat);

         real adiff = ahatnew - ahat;
         real bdiff = bhatnew - bhat;
         real dist = sqrt(adiff * adiff + bdiff * bdiff);
         ahat = ahatnew;
         *bhat = bhatnew;
         if (dist <= toler)
         {
            bool fTryAgain = false;
            // If we bump against the value of ahi or bhi, we increase it and
            // try again. It's possible for the change in a or bhat to get very small but
            // not have converged.
            if (abs(ahat - ahi) <= toler)
            {
               ahi *= 2;
               fTryAgain = true;
            }
            if (abs(*bhat - bhi) <= toler)
            {
               bhi *= 2;
               fTryAgain = true;
            }

            if (!fTryAgain)
               break;
         }
   }

   return ahat;
}
*/
