#include "FastLmmC.h"

//constructor takes an EigenSym object as input
FastLMM_fullrank::FastLMM_fullrank( EigenSym *eigen )
   {
   if ( eigen->n != eigen->rankU )
      {
      Fatal( "error, Eigen dimensions do not agree with FastLMM_fullrank. rankU=%i, n=%i", eigen->rankU, eigen->n );
      }

   this->eigen = eigen;
   this->n = eigen->n;
   Uy = 0;              // null pointer
   UX = 0;              // null pointer
   Ucovar = 0;          // null pointer
   sigg2 = 0;           // null pointer
   sige2 = 0;           // null pointer
   delta = 0;           // null pointer
   beta = 0;            // null pointer
   LL = 0;              // null pointer
   rank_beta = 0;       // null pointer
   pvals = 0;           // null pointer
   beta_0 = 0;          // null pointer
   rank_beta_0 = 0;
   sigg2_0 = -1.0;
   sige2_0 = -1.0;
   delta_0 = -1.0;
   LL_0 = -1.0e100;
   rank_beta_0 = 0;
   ncovar = 0;
   nSNPs = 0;
   modelID = "FastLMM_fullrank";
   }


FastLMM_fullrank::~FastLMM_fullrank()
   {
   delete eigen;
   delete[] Uy;       
   delete[] UX;            
   delete[] Ucovar;  
   delete[] sigg2;   
   delete[] sige2;   
   delete[] delta;  
   delete[] beta;    
   delete[] LL;     
   delete[] rank_beta;
   delete[] pvals;  
   delete[] beta_0;
   }

//SNPs [n by nSNPs] holds the nSNPs SNPs to be tested for each of n individuals
void FastLMM_fullrank::set_SNPs( double *SNPs, int nSNPs_, vector<string> snpLabels )
   {
   if( UX != 0 )
      {
      delete[] UX;
      }

   this->snpLabels = snpLabels;
   nSNPs = nSNPs_;

   UX=new double[nSNPs*n];    //as number of SNPs might change, alsways reallocate memory for UX.
   BlasLapack::gemm(&char_T,  &char_N, &n, &nSNPs, &n, &double_one, eigen->U, &n, SNPs, &n, &double_zero, UX, &n );  //UX=U.T*SNPs
   }

//Phenotype is a vector of length n, holding the (single) phenotype.
void FastLMM_fullrank::set_Phenotype( double *Phenotype, string phenLabel )
   {
   if( Uy != 0 )
      {
      delete[] Uy;
      }

   this->phenLabel = phenLabel;

   Uy = new double[n];
   BlasLapack::gemv(&char_T,&n,&n,&double_one,eigen->U,&n,Phenotype,&int_one,&double_zero,Uy,&int_one);  //Uy=U.T*Phenotype
   }

//covariates is an [n*ncovar] matrix holding the covariates, including the bias
void FastLMM_fullrank::set_covariates( double *covar, int ncovariates )
   {
   ncovar = ncovariates;
   if( Ucovar==0 )
      {
      Ucovar=new double[n*ncovar];
      }
   BlasLapack::gemm( &char_T, &char_N, &n, &ncovar, &n, &double_one, eigen->U, &n, covar, &n, &double_zero, Ucovar, &n );  //UX=U.T*SNPs
   }

//find ML parameters for the null model
void FastLMM_fullrank::trainNull()
   {
   if ( beta_0 != 0 )
      {
      delete[] beta_0;
      }

   beta_0=new double[ncovar];

   double UXdummy=0.0;  //As we do not plug in a SNP, we need a dummy
   //find optimal delta:
   if(FastLMM::brentExhaustive)
      {
      delta_0=FastLMM_fullrank::optimize_delta_Brentintervals(0, ncovar, n, Uy,&UXdummy, eigen->S, Ucovar,numberOptimizationStarts,minLogVal,maxLogVal,brentTol); //Find the optimal delta
      }
   else
      {
      delta_0 = FastLMM_fullrank::optimize_delta_BrentintervalsCheckTriplets(0, ncovar, n, Uy, &UXdummy, eigen->S, Ucovar, numberOptimizationStarts, minLogVal, maxLogVal, brentTol ); //Find the optimal delta
      }
   //get the ML parameters for the optimal delta:
   LL_0 = -1.0 * nLLeval( delta_0, 0, ncovar, n, Uy, &UXdummy, eigen->S, Ucovar, &sigg2_0, &sige2_0, beta_0, &rank_beta_0, true );  //evaluate the LL at delta_0 and get the ML parameters
   printf( "\n------------------------------------------"
           "\nNULL: LL:%f\nsigg2=%f sige2=%f delta=%f",
           LL_0, sigg2_0, sige2_0, delta_0 );
   printf( " beta:");
   for( int i=0; i<ncovar; ++i )
      {
      printf( "%f ", beta_0[i] );
      }
   printf( "\n------------------------------------------"
           "\n" );
   }

//find ML parameters for all SNPs
void FastLMM_fullrank::trainAlt(){      
   //alternative models:
   int d=1;          //assume that always one SNP is being tested
   int d_=d+ncovar;  //the number of weights (dimensionality of a single beta)
   if (beta!=0){     //if the parameters have already been set in a previous iteration delete them.
      delete[] beta;
      delete[] delta;
      delete[] sigg2;
      delete[] sige2;
      delete[] LL;
      delete[] rank_beta;
      delete[] pvals;
   }
   beta=new double[nSNPs*d_];
   delta=new double[nSNPs];
   sigg2=new double[nSNPs];
   sige2=new double[nSNPs];
   LL=new double[nSNPs];
   pvals=new double[nSNPs];
   rank_beta=new int[nSNPs];
   
   //set pointers to the current SNP and ML parameters:
   double* currentUX=UX;
   double* currentBeta=beta;
   double* currentSigg2=sigg2;
   double* currentSige2=sige2;
   double* currentLL=LL;
   double* currentDelta=delta;
   double *currentPval=pvals;
   int* currentRank_beta=rank_beta;

   //train the alternative model for the SNP
   for (int snp=0; snp<nSNPs;++snp){
      if (covLearnType == FULL){
         if(FastLMM::brentExhaustive){
            *currentDelta=FastLMM_fullrank::optimize_delta_Brentintervals(d, ncovar, n, Uy,currentUX, eigen->S, Ucovar,numberOptimizationStarts,minLogVal,maxLogVal,brentTol); //Find the optimal delta
         }
         else{
            *currentDelta=FastLMM_fullrank::optimize_delta_BrentintervalsCheckTriplets(d, ncovar, n, Uy,currentUX, eigen->S, Ucovar,numberOptimizationStarts,minLogVal,maxLogVal,brentTol); //Find the optimal delta
         }
      }else if(covLearnType == ONCE){
         *currentDelta=delta_0;
      }
      *currentLL=-1.0*nLLeval(*currentDelta, d,ncovar, n,Uy, currentUX, eigen->S, Ucovar,currentSigg2 ,currentSige2, currentBeta, currentRank_beta, true);//evaluate the LL at delta_0 and get the ML parameters
      
      //Compute p-value from LR test:
      double df=1.0;
      *currentPval=Gamma::gammaQ(*currentLL-LL_0,0.5*df);
      printf("LL[%i]:%f p=%f\nsigg2=%f sige2=%f delta=%f",snp,*currentLL,*currentPval,*currentSigg2,*currentSige2,*currentDelta);printf("beta:\n ");for(int i=0;i<d_;++i)printf("%f ",currentBeta[i]);printf("\n");

      //move pointers to the next SNP and ML parameters:
      currentUX+=n*d;
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
double FastLMM_fullrank::optimize_delta_BrentintervalsCheckTriplets(const int d, const int ncovar, const int n, const double *Uy, const double *currentUX,const double *S, const double *Ucovar,const int numberIntervals, const double minLogVal,const double maxLogVal,const double brentTol){
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
      nLL_interval[i]=FastLMM_fullrank::nLLeval(exp(log_delta_interval[i]),d,ncovar,n,Uy,currentUX,S,Ucovar,&MLparamsdummy,&MLparamsdummy,&MLparamsdummy,&rankDummy,false);
      ++nevals;
      if (nLL_interval[i]<=best_nLL){
            best_nLL=nLL_interval[i];
            best_log_delta=log_delta_interval[i];
      }
   }
   
   brentFunctor_FastLMM_fullrank func=brentFunctor_FastLMM_fullrank(d, ncovar, n, Uy,currentUX, S, Ucovar);
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
   delete []nLL_interval;
   delete []log_delta_interval;
   printf("\nnevals: %i logdelta: %f\n",nevals,best_log_delta);
   return exp(best_log_delta);
}

//static function:
//finds delta by optimizing over delta intervals
double FastLMM_fullrank::optimize_delta_Brentintervals(const int d, const int ncovar, const int n, const double *Uy, const double *currentUX,const double *S, const double *Ucovar,const int numberIntervals, const double minLogVal,const double maxLogVal,const double brentTol){
   double best_log_delta=-1.0;
   double best_nLL=1.0E100;
   brentFunctor_FastLMM_fullrank func=brentFunctor_FastLMM_fullrank(d, ncovar, n, Uy,currentUX, S, Ucovar);
   int nevals=0;
   for(int i=0;i<numberIntervals;++i){
      double current_log_delta_min=minLogVal+i*(maxLogVal-minLogVal)/numberIntervals;
      double current_log_delta_max=minLogVal+(i+1)*(maxLogVal-minLogVal)/numberIntervals;
      double current_nLL=1.0E+100;
      int nevals_Int=0;
      double current_log_delta=BrentC::minimize(func,current_log_delta_min,current_log_delta_max,brentTol,current_nLL,nevals_Int,brentMaxIter,true);
      nevals+=nevals_Int;
      if (current_nLL<=best_nLL){
         best_nLL=current_nLL;
         best_log_delta=(current_log_delta);
      }
   }
   printf("\nnevals: %i logdelta: %f\n",nevals,best_log_delta);
   return exp(best_log_delta);
}

//static function:
//cleaned up version of the nLL evaluation for the full rank model
double FastLMM_fullrank::nLLeval(const double delta,const int d, const int ncovar,int n,const double *Uy, const double *UX, const double *S, const double *Ucovar, double *sigg2, double *sige2, double *beta, int *rank_beta, const bool returnMLparams){
   //printf("delta: %f\n",delta);
   double det=0.0;                  // the log determinant is the sum of the diagonal of (S+delta*I)
	const int d_=d+ncovar;           //the number of parameters to be estimated is number of SNP coefficients + number of covariates
	const int d_n=d_*n;
   double *Sdi= new double[n];      //holds the sqrt of the diagonal values of (S+delta*I)^-1 (deleted)
	double *res=new double[n];       //the residual vecttor. Will be initialized to ySdi, later overwritten by the residual of the regression (deleted)
	double *beta_temp=new double[max(n,d_)];  //the parameters. Will be initialized to ySdi (deleted)
   double *XSdi=new double[d_n];   //XSdi This one will be used to compute the residual in the end. (deleted)
	double *XSdi_leastsquares=new double[d_n];  //another XSdi. This one will be overwritten in the ML learning. (deleted)

   //--------------------------------------------------
   //compute the determinant and ySdi
	
   for(int i=0;i<n;++i){
      Sdi[i]=S[i]+delta;
   }
   BlasLapack::vLn(n,Sdi,res);        //res holds the Ln(S+deltaI)
   
   for(int i=0;i<n;++i){
      Sdi[i]=1.0/Sdi[i];
   }
   
   BlasLapack::vSqrt( n, Sdi, Sdi );

   double *Sdi_=Sdi;
	double *res_=res;
	double *ySdi_=beta_temp;
	const double *S_=S;
	const double *Uy_=Uy;
	for (int i=0;i<n;++i){
		det+=*res_;          //res holds the Ln(S+deltaI)
		*res_=*Uy_* *Sdi_;
		*ySdi_=*res_;        //res is overwritten with ySdi
		++Sdi_;
		++res_;
		++ySdi_;
		++Uy_;
		++S_;
	}

   //compute XSdi. The inner two for loops concatenate UX and Ucovar
	int ind0=0;
	for(int snp=0;snp<d;++snp){
		for(int ind=0;ind<n;++ind){
			XSdi[ind0]=UX[ind0]*Sdi[ind];
			XSdi_leastsquares[ind0]=XSdi[ind0];
			++ind0;
		}
	}
	int ind1=0;
	for(int cov=0;cov<ncovar;++cov){
		for(int ind=0;ind<n;++ind){
			XSdi[ind0]=Ucovar[ind1]*Sdi[ind];
			XSdi_leastsquares[ind0]=XSdi[ind0];
			++ind0;
			++ind1;
		}
	}


	delete[] Sdi;           //Sdi is no longer needed

   //-----------------------------------------------
   //ML learning of beta: This could be speeded up, by assigning a fixed workspace to avoid reallocation at each call (see commented out part)
   int INFO=0;//error messages from the solver go here
	double RCOND=3E-8;//the rank cutoff on (XSdi.T*XSdi) of the least squares solver
	int*JPVT=new int[n];//some integer workspace of the least squares (deleted)
	for (int i=0;i<n;++i){
		JPVT[i]=0;
	}
	int LDB=max(d_,n);
   double work_beta_request=0.0;
   BlasLapack::gelsy(&n,&d_,&int_one,XSdi_leastsquares, &n,beta_temp, &LDB,JPVT,&RCOND, rank_beta, &work_beta_request, &int_minusone, &INFO);//workspace request
	int lwork_beta= (int)work_beta_request;//the optimal size of the workspace
	double *work_betanew=new double[lwork_beta];//create the workspace
	BlasLapack::gelsy(&n,&d_,&int_one,XSdi_leastsquares, &n,beta_temp, &LDB,JPVT,&RCOND, rank_beta,work_betanew, &lwork_beta, &INFO);//call Least squares solver
	delete[] work_betanew;//delete the workspace
	
   delete[]JPVT;           //this was only needed for the least squares
	delete[] XSdi_leastsquares;//this has been overwritten in the regression
   //----------------------------------------------------------------------------
	
   //compute the residual vector res: (ySdi-XSdi*beta_temp)
	BlasLapack::gemv(&char_N, &n,&d_,&double_minusone,XSdi,&n, beta_temp, &int_one,&double_one,res,&int_one );
	delete[] XSdi;          //XSdi is no longer needed.
	
   double sigg2_temp=0.0;  //the ML sigma_g^2
	for (int i=0;i<n;++i){
		sigg2_temp+=res[i]*res[i];
	}
	delete[] res;           //residual is no longer needed.
	sigg2_temp/=n;
   if (returnMLparams){//only return the ML params if needed. (e.g. no need during Brents search) This avoids danger of Memory errors
      *sigg2=sigg2_temp;
      *sige2=delta**sigg2;
      BlasLapack::copy(&d_,beta_temp,&int_one,beta,&int_one);
   }

   delete[] beta_temp;
   double nLL=0.5*(n*(_log2pi+log(sigg2_temp)+1.)+det);//compute the negative LL usning the ML parameters at delta

	return nLL;
}


ResultsLMM FastLMM_fullrank::getResults(){
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
   return ResultsLMM(modelID_res, snpLabels_res, phenLabel_res, n, ncovar, index_res, LL_alt_res, LL_0_res, pvals_res, snpWeights_alt_res, bias_alt_res, sigg2_alt_res, sige2_alt_res, delta_alt_res, bias_0_res, sige2_0_res, sigg2_0_res, delta_0_res, qvals_res);
}


//Functor holdong all the needed members of FastLMM_fullrank to be handed to the Brent's routine
//the constructor gets the members of FastLMM_fullrank
brentFunctor_FastLMM_fullrank::brentFunctor_FastLMM_fullrank(const int d_, const int ncovar_, const int n_, const double *Uy_,const  double *UX_, const double *S_, const double *Ucovar_){
   //member parameters from FastLMM_fullrank:
   d=d_;
   ncovar=ncovar_;
   n=n_;
   Uy=Uy_;
   UX=UX_;
   S=S_;
   Ucovar=Ucovar_;
}

//destructor has nothing to do
brentFunctor_FastLMM_fullrank::~brentFunctor_FastLMM_fullrank(){}

//the operator evaluates nLLeval
double brentFunctor_FastLMM_fullrank::operator()(const double logdelta)
   {
   double delta=exp( logdelta );
   
   //the rest are dummies
   int rank_beta_dummy=0;
   double dummy_ML_params=0.0;
   
   return FastLMM_fullrank::nLLeval( delta, d, ncovar, n, Uy, UX, S, Ucovar, &dummy_ML_params, &dummy_ML_params, &dummy_ML_params, &rank_beta_dummy, false );
   }
