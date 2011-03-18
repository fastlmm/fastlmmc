#if !defined(Gamma_h)
#define Gamma_h

class Gamma
{
public:
   static double gammaQ(double x, double a);
   static double gammaIncLower(double x, double a);
   static double gammaIncUpper(double x, double a);
   static double gammaP(double x, double a);
   static double cdf(double x, double k, double theta);
   static double inv(double p, double k, double theta);
};

#endif  //Gamma_h
