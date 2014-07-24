#if !defined(Beta_h)
#define Beta_h

class Beta{
public:
   static real BetaCF(real x, real a, real b);
   static real BetaInc(real x, real a , real b);
   static real Cdf(real x, real a, real b);
   static real Pdf(real x, real a, real b);
   static real Stats(real a, real b, real* var);
   static real Inv(real p, real a, real b);
   static real Eval(real y, real lo, real hi, real s);
   //static real MLE(real *dataIn, real *bhat);
};


const int MAXITER = 400;
const real TINY = (real)1.0e-30;
const real EPS = (real)2.2204e-016;

#endif //Beta_h
