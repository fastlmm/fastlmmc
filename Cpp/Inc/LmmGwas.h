#if !defined(LmmGwas_h)
#define LmmGwas_h

void LmmGwas();
void LmmGwas_oldfortiming();
int computeEigenSym_Kernel(string kernelFile, EigenSym *eigen);
int computeEigenSym_SnpSample(string snpDataCovFile, EigenSym *eigen);
int runGWAS_lowrank(EigenSym *eigen, ResultsLMM *results, CSnpData *snpData, CPhenData *phenData, double *covariates, int n_covariates);
int runGWAS_fullrank(EigenSym *eigen, ResultsLMM *results, CSnpData *snpData, CPhenData *phenData, double *covariates, int n_covariates);

#endif   //LmmGwas
