#include "FastLmmC.h"

//Parameters: SNPs [nIndividuals by nSNPs] Matrix stored in column-major order. This will hold the result.
//NaNs will be set to 0.0 in the result.
//maxValue: maximum value a SNP can have (usually value for minor-minor)
//unitVariance: true: apply the imputation as per Price et al. 2007, false:only zeroMean the SNPs
void ImputeAndZeroMeanSNPs(double *SNPs, const int nIndividuals, const int nSNPs, const double maxValue, const bool unitVariance){
   for (int snp=0;snp++<nSNPs;){
      double n_observed=0.0;
      double mean_s=0.0;//the mean of a SNP over all observed individuals
      for(int ind=0;ind<nIndividuals;++ind){
         if (SNPs[ind]==SNPs[ind]){//check for not NaN
            mean_s+=SNPs[ind];
            ++n_observed;
         }
      }
      double oneover_sqrt_pi=(1.0+mean_s)/(2.0+maxValue*n_observed);//compute p_i as in the Eigenstrat paper (Price et al 2007) Note that mean_s still holds the sum.
      oneover_sqrt_pi=1.0/sqrt(oneover_sqrt_pi*(1-oneover_sqrt_pi));//compute inverse of sqrt(p_i*(1-p_i)) (multiplication is faster than division...)
      mean_s/=n_observed;//compute the mean over observed individuals for the current SNP
      for(int ind=0;ind<nIndividuals;++ind){
         if (SNPs[ind]!=SNPs[ind]){//check for NaN
            SNPs[ind]=0.0;
         }
         else{
            SNPs[ind]-=mean_s;//subtract the mean from the data
            if(unitVariance){
               SNPs[ind]*=oneover_sqrt_pi;//unit variance as well
            }
         }
      }
      SNPs+=nIndividuals;
   }
}

//Parameters: SNPs [nIndividuals by nSNPs] Matrix stored in column-major order. This will hold the result.
void preprocessSNPsForBuildingCovariance(double *SNPs, const int nIndividuals, const int nSNPs){
   double one_over_sqrt_nSNPs=1.0/sqrt(1.0*nSNPs);
   //printf("%f\n",one_over_sqrt_nSNPs);
   int snpMax=nIndividuals*nSNPs;
   for (int indiv=0;indiv<nIndividuals;++indiv){
      double mean_i=0.0;//the mean of one individual over all SNPs
      //int snpMax=1+indiv+nSNPs*(indiv-1);
      int count=0;
      for(int snp=indiv;snp<snpMax;snp+=nIndividuals){
         mean_i+=SNPs[snp];++count;
      }
      //printf("%i ",count);
      mean_i/=nSNPs;//compute the mean of the current individual over all SNPs
      //printf("%f\n",mean_i);
      for(int snp=indiv;snp<snpMax;snp+=nIndividuals){
         SNPs[snp]-=mean_i;
         SNPs[snp]*=one_over_sqrt_nSNPs;
      }
   }
}
