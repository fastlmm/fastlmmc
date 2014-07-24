#if !defined(FisherF_h)
#define FisherF_h
#include "FastLmmC.h"
class FisherF{
   public:
      static real Cdf(real x, real v1, real v2);
      static real Pdf(real x, real v1, real v2);
      static real Inv(real p, real v1, real v2);
      static real Stats(real v1, real v2, real *var);
};

#endif //Fisher_F_h
