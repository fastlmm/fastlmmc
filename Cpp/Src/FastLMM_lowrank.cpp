#include "FastLmmC.h"

//constructor takes an EigenSym object as input
FastLMM_lowrank::FastLMM_lowrank(EigenSym *eigen){
   if (eigen->n<=eigen->rankU){
      printf("error, Eigen dimensions do not agree with FastLMM_lowrank. rankU=%i, n=%i",eigen->rankU,eigen->n);
      throw(1);
   }
   this->eigen=eigen;
   this->n=eigen->n;
   this->rankU=eigen->rankU;
   Uy=0;             //null pointer
   y_UUy=0;          //null pointer
   UX=0;             //null pointer
   X_UUX=0;          //null pointer
   Ucovar=0;         //null pointer
   covar_UUcovar=0;  //null pointer
   sigg2=0;          //null pointer
   sige2=0;          //null pointer
   delta=0;          //null pointer
   beta=0;           //null pointer
   LL=0;             //null pointer
   rank_beta=0;      //null pointer
   pvals=0;          //null pointer
   beta_0=0;         //null pointer
   rank_beta_0=0;
   sigg2_0=-1.0;
   sige2_0=-1.0;
   delta_0=-1.0;
   LL_0=-1.0e100;
   rank_beta_0=0;
   ncovar=0;
   nSNPs=0;
   modelID="FastLMM_lowrank";
}


FastLMM_lowrank::~FastLMM_lowrank(void)
{
   delete[] Uy;       
   delete[] UX;       
   delete eigen;     
   delete[] Ucovar;
   delete[] covar_UUcovar;
   delete[] X_UUX;
   delete[] y_UUy;
   delete[] sigg2;   
   delete[] sige2;   
   delete[] delta;  
   delete[] beta;    
   delete[] LL;     
   delete[] rank_beta;
   delete[] pvals;     
   delete[] beta_0;
}

//SNPs [n by nSNPs] holds the nSNPs SNPs to be tested for each of n individuals SNPs will be overwritten by X_UUX.
void FastLMM_lowrank::set_SNPs(double *SNPs, int nSNPs_, vector<string> snpLabels){
   nSNPs=nSNPs_;
   this->snpLabels=snpLabels;
   
   if(UX!=0){ delete[] UX; }
   if(X_UUX != 0){ delete[] X_UUX; }
   
   UX=new double[nSNPs*rankU];      //as number of SNPs might change, alsways reallocate memory for UX.
   BlasLapack::gemm(&char_T,  &char_N, &rankU, &nSNPs, &n, &double_one, eigen->U, &n, SNPs, &n, &double_zero, UX, &rankU );   //UX=U.T*SNPs

   int sizeX=n*nSNPs;
   X_UUX=new double[sizeX];
   BlasLapack::copy(&sizeX, SNPs, &int_one, X_UUX, &int_one);            //Initialize X_UUX with the SNPs. This is an [n * nSNPs] matrix
   BlasLapack::gemm(&char_N,  &char_N, &n, &nSNPs, &rankU, &double_minusone, eigen->U, &n, UX, &rankU, &double_one, X_UUX, &n ); // This matrix multiply performs X_UUX= SNPs-U*UX
}

//Phenotype is a vector of length n, holding the (single) phenotype.
void FastLMM_lowrank::set_Phenotype(double *Phenotype, string phenLabel){
   if(Uy!=0){ delete[] Uy; }
   if(y_UUy!=0){ delete[] y_UUy; }
   this->phenLabel=phenLabel;

   Uy=new double[rankU];         //The phenotypes transformed to low rank
   BlasLapack::gemv(&char_T,&n,&rankU,&double_one, eigen->U, &n,Phenotype,&int_one,&double_zero,Uy,&int_one); // Uy = U.T*Phenotype
   
   y_UUy=new double[n];
   BlasLapack::copy(&n, Phenotype, &int_one, y_UUy ,&int_one);        //Initialize with original Phenotype
   BlasLapack::gemv(&char_N, &n, &rankU, &double_minusone, eigen->U, &n, Uy, &int_one, &double_one, y_UUy, &int_one); // U_UUy = Phenotype-U*Uy
}

//covariates is an [n*ncovar] matrix holding the covariates, including the bias
void FastLMM_lowrank::set_covariates(double *covar,int ncovariates){
   ncovar=ncovariates;                 //the number of covariates (including the fixed effects)
   if(Ucovar!=0){ delete[] Ucovar; }
   if(covar_UUcovar!=0){ delete[] covar_UUcovar; }
   
   Ucovar=new double[ncovar*rankU];    //The covariates transformed to the low rank space
   BlasLapack::gemm(&char_T,  &char_N, &rankU, &ncovar, &n, &double_one, eigen->U, &n, covar, &n, &double_zero, Ucovar, &rankU );//UX=U.T*SNPs
   
   int sizeCovar=ncovar*n;
   covar_UUcovar=new double[sizeCovar];
   BlasLapack::copy(&sizeCovar, covar, &int_one, covar_UUcovar, &int_one); //Initialize with the original covariates
   BlasLapack::gemm(&char_N,  &char_N, &n, &ncovar, &rankU, &double_minusone, eigen->U, &n, Ucovar, &rankU, &double_one, covar_UUcovar, &n );//covar_UUcovar=covar-U*Ucovar
}

//find ML parameters for the null model
void FastLMM_lowrank::trainNull(){
   if(beta_0!=0){ delete[] beta_0; }
   beta_0=new double[ncovar];

   double UXdummy=0.0;     //As We don't have a SNP, we pass a dummy variable
   
   //find optimal delta:
   if(FastLMM::brentExhaustive){
      delta_0=FastLMM_lowrank::optimize_delta_Brentintervals(0, ncovar, n,rankU, Uy,&UXdummy,y_UUy,&UXdummy, eigen->S, Ucovar,covar_UUcovar,numberOptimizationStarts,minLogVal,maxLogVal,brentTol); //Find the optimal delta
   }
   else{
      delta_0=FastLMM_lowrank::optimize_delta_BrentintervalsCheckTriplets(0, ncovar, n,rankU, Uy,&UXdummy,y_UUy,&UXdummy, eigen->S, Ucovar,covar_UUcovar,numberOptimizationStarts,minLogVal,maxLogVal,brentTol); //Find the optimal delta
   }
   //get the ML parameters for the optimal delta:
   LL_0=-1.0*nLLeval(delta_0, 0,ncovar, n,rankU,Uy, &UXdummy,y_UUy,&UXdummy, eigen->S, Ucovar,covar_UUcovar,&sigg2_0 ,&sige2_0, beta_0, &rank_beta_0, true);//evaluate the LL at delta_0 and get the ML parameters
   printf("------------------------------------------\nNULL: LL:%f\nsigg2=%f sige2=%f delta=%f",LL_0,sigg2_0,sige2_0,delta_0);printf(" beta: ");for(int i=0;i<ncovar;++i)printf("%f ",beta_0[i]);printf("\n---------------------------------------------\n");
}

//find ML parameters for all SNPs
void FastLMM_lowrank::trainAlt(){      
   //alternative models:
   int d = 1;                //assume that always one SNP is being tested
   int d_ = d + ncovar;      //the number of weights (dimensionality of a single beta)
   
   if(beta!=0){   delete[] beta; }   //if the parameters have already been set in a previous iteration delete them.
   if(delta!=0){  delete[] delta; }
   if(sigg2!=0){  delete[] sigg2; }
   if(sige2!=0){  delete[] sige2; }
   if(LL!=0){     delete[] LL; }
   if(rank_beta!=0){ delete[] rank_beta; }
   if(pvals!=0){  delete[] pvals; }

   beta=new double[nSNPs*d_];       //create all the variables
   delta=new double[nSNPs];
   sigg2=new double[nSNPs];
   sige2=new double[nSNPs];
   LL=new double[nSNPs];
   pvals=new double[nSNPs];
   rank_beta=new int[nSNPs];
   
   //set pointers to the current SNP and ML parameters:
   double* currentUX=UX;
   double* currentX_UUX=X_UUX;
   double* currentBeta=beta;
   double* currentSigg2=sigg2;
   double* currentSige2=sige2;
   double* currentLL=LL;
   double* currentDelta=delta;
   double *currentPval=pvals;
   int* currentRank_beta=rank_beta;
   for (int snp=0; snp<nSNPs;++snp){//train the alternative model for the SNP   
      //find delta:
      if (covLearnType == FULL){       //Learn delta
         if(FastLMM::brentExhaustive){
            *currentDelta=FastLMM_lowrank::optimize_delta_Brentintervals(d, ncovar, n,rankU, Uy,currentUX,y_UUy,currentX_UUX, eigen->S, Ucovar,covar_UUcovar,numberOptimizationStarts,minLogVal,maxLogVal,brentTol); //Find the optimal delta
         }
         else{
            *currentDelta=FastLMM_lowrank::optimize_delta_BrentintervalsCheckTriplets(d, ncovar, n,rankU, Uy,currentUX,y_UUy,currentX_UUX, eigen->S, Ucovar,covar_UUcovar,numberOptimizationStarts,minLogVal,maxLogVal,brentTol); //Find the optimal delta
         }
      }else if(covLearnType == ONCE){  //Or: use the lock-in. So it equals delta_0
         *currentDelta=delta_0;
      }

      //evaluate the LL and find the parameters
      *currentLL=-1.0*nLLeval(*currentDelta, d,ncovar, n,rankU,Uy, currentUX,y_UUy,currentX_UUX, eigen->S, Ucovar,covar_UUcovar,currentSigg2 ,currentSige2, currentBeta, currentRank_beta, true);//evaluate the LL at delta_0 and get the ML parameters
      
      //Compute p-value from LR test:
      double df=1.0;
      *currentPval=Gamma::gammaQ(*currentLL-LL_0,0.5*df);
      printf("LL[%i]:%f p=%f\nsigg2=%f sige2=%f delta=%f",snp,*currentLL,*currentPval,*currentSigg2,*currentSige2,*currentDelta);printf("beta:\n ");for(int i=0;i<d_;++i)printf("%f ",currentBeta[i]);printf("\n");

      //move pointers to the next SNP and ML parameters:
      currentUX+=rankU*d;
      currentX_UUX+=n*d;
      currentBeta+=d_;
      ++currentSigg2;
      ++currentSige2;
      ++currentLL;
      ++currentDelta;
      ++currentRank_beta;
      ++currentPval;
   }
}

//static function:
//finds delta by optimizing over delta intervals
double FastLMM_lowrank::optimize_delta_BrentintervalsCheckTriplets(const int d, const int ncovar, const int n,const int rankU, const double *Uy, const double *currentUX, const double *y_UUy, const double *currentX_UUX,const double *S, const double *Ucovar, const double *covar_UUcovar,const int numberIntervals, const double minLogVal,const double maxLogVal,const double brentTol){
   double *nLL_interval=new double[numberIntervals+1];
   double *log_delta_interval=new double[numberIntervals+1];
   
   double best_log_delta=-1.0;
   double best_nLL=1.0E100;

   double MLparamsdummy=0.0;
   int rankDummy=0;
   int nevals=0;

   //evaluate the boundaries of all intervals
   for(int i=0;i<=numberIntervals;++i){
      log_delta_interval[i]=minLogVal+i*(maxLogVal-minLogVal)/numberIntervals;
      nLL_interval[i]=FastLMM_lowrank::nLLeval(exp(log_delta_interval[i]),d,ncovar,n,rankU,Uy,currentUX,y_UUy,currentX_UUX,S,Ucovar,covar_UUcovar,&MLparamsdummy,&MLparamsdummy,&MLparamsdummy,&rankDummy,false);
      ++nevals;
      if (nLL_interval[i]<=best_nLL){
            best_nLL=nLL_interval[i];
            best_log_delta=log_delta_interval[i];
      }
   }
   
   brentFunctor_FastLMM_lowrank func=brentFunctor_FastLMM_lowrank(d, ncovar, n,rankU, Uy,currentUX,y_UUy,currentX_UUX, S, Ucovar,covar_UUcovar);
   for(int i=1;i<numberIntervals;++i){
      if (nLL_interval[i]<nLL_interval[i+1] && nLL_interval[i]<nLL_interval[i-1]){//check wether a local optimum exists in this triplet
         double current_nLL=1.0E+100;
         int nevals_Int=0;
         //If there is a local optimum, optimize over the current triplet:
         double current_log_delta=BrentC::minimize(func,log_delta_interval[i-1],log_delta_interval[i+1],brentTol,current_nLL,nevals_Int,brentMaxIter,true);
         nevals+=nevals_Int;
         if (current_nLL<=best_nLL){
            best_nLL=current_nLL;
            best_log_delta=(current_log_delta);
         }
      }
   }
   delete []log_delta_interval;
   delete []nLL_interval;
   printf("\nnevals: %i logdelta: %f\n",nevals,best_log_delta);
   return exp(best_log_delta);
}

//static function:
//finds delta by optimizing over log delta intervals
double FastLMM_lowrank::optimize_delta_Brentintervals(const int d, const int ncovar, const int n,const int rankU, const double *Uy, const double *currentUX, const double *y_UUy, const double *currentX_UUX,const double *S, const double *Ucovar, const double *covar_UUcovar,const int numberIntervals, const double minLogVal,const double maxLogVal,const double brentTol){
   double best_log_delta=-1.0;
   double best_nLL=1.0E100;
   
   //create the Functor. This will be passed to the optmizer
   brentFunctor_FastLMM_lowrank func=brentFunctor_FastLMM_lowrank(d, ncovar, n,rankU, Uy,currentUX,y_UUy,currentX_UUX, S, Ucovar,covar_UUcovar);
   
   int nevals=0;  //count the total number of LL-function evaluations

   //Perform a Brent optimization on each interval
   for(int i=0;i<numberIntervals;++i){
      double current_log_delta_min=minLogVal+i*(maxLogVal-minLogVal)/numberIntervals;        //The mimimum of the current Brent interval
      double current_log_delta_max=minLogVal+(i+1)*(maxLogVal-minLogVal)/numberIntervals;    //The maximum of the current Brent interval
      double current_nLL=1.0E+100;     //This will hold the optimal nLL
      int nevals_Int=0;                //The number of evaluations on the current interval
      double current_log_delta=BrentC::minimize(func,current_log_delta_min,current_log_delta_max,brentTol,current_nLL,nevals_Int,brentMaxIter,true);
      
      nevals+=nevals_Int;              //Increment the total number of evaluations by #evaluations on the current one
      
      if (current_nLL<=best_nLL){      //Check, whether we have found a better optimum on the current interval
         best_nLL=current_nLL;         //If Yes:   Set the optimum to the current one
         best_log_delta=(current_log_delta);
      }
   }
   printf("\nnevals: %i logdelta: %f\n",nevals,best_log_delta);
   return exp(best_log_delta);         //The Brent optimization is in log space. We want the actual delta. => exp()
}


//static function:
//cleaned up version of the nLL evaluation for the full rank model
double FastLMM_lowrank::nLLeval(const double delta,const int d, const int ncovar,int n, const int rankU,const double *Uy, const double *UX,const double *y_UUy, const double *X_UUX, const double *S, const double *Ucovar, const double *covar_UUcovar, double *sigg2, double *sige2, double *beta, int *rank_beta, const bool returnMLparams){
   int n_=n+rankU;                              //we regress over both, the n residuals and the rankU latent individuals
	int d_=d+ncovar;                             //the number of parameters to be estimated is number of SNP coefficients + number of covariates
   int nd=n_*d_;
	double *Sdi= new double[rankU];              //holds the sqrt of the diagonal values of (S+delta*I)^-1 (deleted)
	double *res=new double[n_];                  //the residual vecttor. Will be initialized to ySdi, later overwritten by the residual of the regression (deleted)
	double *beta_temp=new double[max(n_,d_)];    //the parameters. Will be initialized to ySdi (deleted)
   double *XSdi=new double[nd];                 //XSdi This one will be used to compute the residual in the end. (deleted)
	double *XSdi_leastsquares=new double[nd];    //another XSdi. This one will be overwritten in the ML learning. (deleted)
   //--------------------------------------------------
   //compute the determinant and ySdi
   double delta_sqrt=1.0/sqrt(delta);
   double det=(n-rankU)*log(delta);             // the log determinant is the sum of the diagonal of (S+delta*I)
	double *Sdi_=Sdi;
	double *res_=res;
	double *ySdi_=beta_temp;
	const double *S_=S;
	const double *Uy_=Uy;

   for(int i=0;i<rankU;++i){
      Sdi[i]=S[i]+delta;
   }
   BlasLapack::vLn(rankU,Sdi,res);              //res holds the Ln(S+deltaI)
  
   for(int i=0;i<rankU;++i){
      Sdi[i]=1.0/Sdi[i];                        //invert Sdi. Sdi is alway > than zero, as S[i]>=0.0 && delta>0.0
   }
   
   BlasLapack::vSqrt( rankU, Sdi, Sdi );

	double *beta_=beta;
	for (int i=0;i++<rankU;){
		det+=*res_;
		*res_=*Uy_* *Sdi_;
		++Sdi_;
		++res_;
		++Uy_;
		++S_;
	}
   BlasLapack::copy(&n,y_UUy,&int_one,res_,&int_one);
	BlasLapack::scal(&n,&delta_sqrt,res_,&int_one);
	BlasLapack::copy(&n_,res,&int_one,beta_temp,&int_one);
	
   double *XSdi_ind = XSdi;
	const double *UX_ind=UX;
	double *Sdi_ind=Sdi;
	const double *X_UUX_ind=X_UUX;
	for(int snp=0;snp++<d;){
		Sdi_ind=Sdi;
		for(int ind=0;ind++<rankU;){
			*XSdi_ind=*UX_ind * *Sdi_ind;
			++UX_ind;
			++XSdi_ind;
			++Sdi_ind;
		}
		BlasLapack::copy(&n,X_UUX_ind,&int_one,XSdi_ind,&int_one);
		BlasLapack::scal(&n,&delta_sqrt,XSdi_ind,&int_one);
		X_UUX_ind+=n;
		XSdi_ind+=n;
	}
	const double *covar_UUcovar_ind=covar_UUcovar;
	int ind1=0;
	for(int cov=0;cov<ncovar;++cov){
		for(int ind=0;ind<rankU;++ind){
			*XSdi_ind=Ucovar[ind1]*Sdi[ind];
			++XSdi_ind;
			++ind1;
		}
		BlasLapack::copy(&n,covar_UUcovar_ind,&int_one,XSdi_ind,&int_one);
		BlasLapack::scal(&n,&delta_sqrt,XSdi_ind,&int_one);
		covar_UUcovar_ind+=n;
		XSdi_ind+=n;
	}
   BlasLapack::copy(&nd,XSdi,&int_one,XSdi_leastsquares,&int_one);

   //double *XX= new double[d_*d_];
   //BlasLapack::gemm(&char_T, &char_N, &d_, &d_, &n_, &double_one, XSdi, &n_, XSdi, &n_, &double_zero, XX, &d_);
   
	delete[] Sdi;  //Sdi is no longer needed
	
   //-----------------------------------------------
   //ML learning of beta: This could be speeded up, by assigning a fixed workspace to avoid reallocation at each call (see commented out part)
   int INFO=0;             //error messages from the solver go here
	double RCOND=3E-8;      //the rank cutoff on (XSdi.T*XSdi) of the least squares solver
	int*JPVT=new int[n];    //some integer workspace of the least squares (deleted)
	for (int i=0;i<n;++i){
		JPVT[i]=0;
	}
	int LDB=max(d_,n_);

   double work_beta_dummy=0.0;                  //Storage for the optimal workspace size
	BlasLapack::gelsy(&n_,&d_,&int_one,XSdi_leastsquares, &n_,beta_temp, &LDB,JPVT,&RCOND, rank_beta,&work_beta_dummy, &int_minusone, &INFO);   //workspace request
	int lwork_beta= (int)work_beta_dummy;        //the optimal size of the workspace
	double *work_betanew=new double[lwork_beta]; //create the workspace
	BlasLapack::gelsy(&n_,&d_,&int_one,XSdi_leastsquares, &n_,beta_temp, &LDB,JPVT,&RCOND, rank_beta,work_betanew, &lwork_beta, &INFO);   //call Least squares solver
	delete[] work_betanew;                       //delete the workspace
   
   delete[]JPVT;                    //this was only needed for the least squares
	delete[] XSdi_leastsquares;      //this has been overwritten in the regression
   //----------------------------------------------------------------------------

   //compute the residual vector res: (ySdi-XSdi*beta_temp):
	BlasLapack::gemv(&char_N, &n_,&d_,&double_minusone,XSdi,&n_, beta_temp, &int_one,&double_one,res,&int_one );
	delete[] XSdi;          //XSdi is no longer needed.

	double sigg2_temp=0.0;  //the ML sigma_g^2
	for (int i=0;i<n_;++i){
		sigg2_temp+=res[i]*res[i];
	}
	delete[] res;        //residual is no longer needed.
	sigg2_temp/=n;       //divide by the original n, not n_
   if (returnMLparams){ //only return the ML params if needed. (e.g. no need during Brents search) This avoids danger of Memory errors
      *sigg2=sigg2_temp;
      *sige2=delta**sigg2;
      BlasLapack::copy(&d_,beta_temp,&int_one,beta,&int_one);
   }
   
   delete[] beta_temp;  //delete the ML beta_temp
   double nLL = 0.5 * ( n * ( _log2pi + log(sigg2_temp) + 1. ) + det);  //compute the negative LL usning the ML parameters at delta
   
   //printf("nLL:%f    sigg2:%f    sige2:%f\n",nLL,sigg2_temp,sigg2_temp*delta);
	return nLL;
}

//Functor holdong all the needed members of FastLMM_lowrank to be handed to the Brent's routine
//the constructor gets the members of FastLMM_lowrank
brentFunctor_FastLMM_lowrank::brentFunctor_FastLMM_lowrank(const int d_, const int ncovar_, const int n_,const int rankU_, const double *Uy_,const  double *UX_, const double *y_UUy_,const  double *X_UUX_, const double *S_, const double *Ucovar_,const double *covar_UUcovar_){
   //member parameters from FastLMM_lowrank:
   d=d_;
   rankU=rankU_;
   ncovar=ncovar_;
   n=n_;
   Uy=Uy_;
   UX=UX_;
   y_UUy=y_UUy_;
   X_UUX=X_UUX_;
   covar_UUcovar=covar_UUcovar_;
   S=S_;
   Ucovar=Ucovar_;
}

//destructor has nothing to do
brentFunctor_FastLMM_lowrank::~brentFunctor_FastLMM_lowrank(){}

//the operator evaluates nLLeval
double brentFunctor_FastLMM_lowrank::operator()(const double logdelta){
   double delta = exp(logdelta);
   
   //the rest are dummies:
   int rank_beta = 0;
   double dummy_MLparams=0.0;
   
   return FastLMM_lowrank::nLLeval(delta, d, ncovar, n, rankU, Uy, UX, y_UUy, X_UUX, S, Ucovar, covar_UUcovar, &dummy_MLparams, &dummy_MLparams, &dummy_MLparams, &rank_beta, false);
}

ResultsLMM FastLMM_lowrank::getResults(){
   int *index_res=new int[nSNPs];
   double *qvals_res = new double[nSNPs];
   Qval::computeQvals(this->pvals,qvals_res,index_res,nSNPs,1.0,nSNPs);
   string modelID_res=this->modelID;
   string phenLabel_res=this->phenLabel;
   vector<string> snpLabels_res(nSNPs);
   double *LL_0_res = new double[nSNPs];
   double *LL_alt_res=new double[nSNPs];
   double *pvals_res = new double[nSNPs];
   double *snpWeights_alt_res = new double[nSNPs];
   double *bias_alt_res = new double[nSNPs];
   double *sigg2_alt_res = new double[nSNPs];
   double *sige2_alt_res = new double[nSNPs];
   double *delta_alt_res = new double[nSNPs];
   double *bias_0_res = new double[nSNPs];
   double *sige2_0_res = new double[nSNPs];
   double *sigg2_0_res = new double[nSNPs];
   double *delta_0_res = new double[nSNPs];
   for (int i=0;i<nSNPs;++i){
      snpLabels_res[i]=this->snpLabels[i];
      LL_alt_res[i]=this->LL[i];
      LL_0_res[i] =  this->LL_0;
      pvals_res[i] = this->pvals[i];
      snpWeights_alt_res[i] = this->beta[(1+ncovar)*i];
      bias_alt_res[i] = this->beta[(1+ncovar)*i+ncovar];
      sigg2_alt_res[i] = this->sigg2[i];
      sige2_alt_res[i] = this->sige2[i];
      delta_alt_res[i] = this->delta[i];
      bias_0_res[i] = this->beta_0[(ncovar-1)];
      sige2_0_res[i] = this->sige2_0;
      sigg2_0_res[i] = this->sigg2_0;
      delta_0_res[i] = this->delta_0;
   }
   //ResultsLMM(string modelID, vector<string> snpLabels, string phenLabel, int n_individuals, int ncovariates, int *index, double *LL_alt, double *LL_0, double *pvals, double *snpWeights_alt, double *bias_alt, double *sigg2_alt, double *sige2_alt, double *delta_alt, double *bias_0, double *sige2_0, double *delta_0, double *qvals);
   return ResultsLMM(modelID_res, snpLabels_res, phenLabel_res, n, ncovar, index_res, LL_alt_res, LL_0_res, pvals_res, snpWeights_alt_res, bias_alt_res, sigg2_alt_res, sige2_alt_res, delta_alt_res, bias_0_res, sige2_0_res, sigg2_0_res, delta_0_res, qvals_res);
}

