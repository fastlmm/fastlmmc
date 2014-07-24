#if !defined(Gamma_h)
#define Gamma_h
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


class Gamma
{
public:
   static real gammaQ(real x, real a);
   static real gammaIncLower(real x, real a);
   static real gammaIncUpper(real x, real a);
   static real gammaP(real x, real a);
   static real cdf(real x, real k, real theta);
   static real inv(real p, real k, real theta);
   static real Psi(double k);
};

#endif  //Gamma_h
