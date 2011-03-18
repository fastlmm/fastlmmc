/*
 * Copyright (c) Microsoft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * FastLmmC - {Fast Linear Mixed Model in C++}
 *
 *         File Name:   LmmGwas.cpp
 *           Version:   1.00
 *            Author:   BobD
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Module Purpose:   This file implements the GWAS functionality in FastLmmC
 *
 *    Change History:   
 *
 * Test Files: 
 */

/*
 * Include Files
 */
#include "FastLmmC.h"

/*
 * Externs we need
 */
double nLLeval_lowrank(double delta, int d, int ncovar, int n, int rank_U, double *Uy, double *y_UUy, double *UX, double *X_UUX, double *S, double *Ucovar, double *covar_UUcovar, double *sigg2, double *sige2, double *beta, int*rank_beta, int*lwork_beta, double*work_beta);

void test_EV_symm();
void test_SVD_econ();
void test_LL_evals();
void test_SVD_full();
void test_LL_evals2();

/*
 * Globals we define
 */
bool runAlternativeModels = true;
bool forceFullRank = false;

/*
 * The code
 */
void LmmGwas()
   {
   EigenSym *eigen = new EigenSym();   //create an empty EigenSym object. This will be filled in the following.

   // Read all the input files
   if (computeEigenType==NOSAVE || computeEigenType==SAVETODISK){
      printf("computing eigensym\n");
      if ( kernelFile.size()>0 ){   //a kernelFile has been specified
         computeEigenSym_Kernel(kernelFile, eigen);         //create EigenSym from the kernel stored in kernelFile
      }else{   //a snpDataCovFile has been specified
         computeEigenSym_SnpSample(snpDataCovFile, eigen);  //create EigenSym of a covariance kernel from the SNP sample stored in snpDataCovFile
      }

      if (computeEigenType==SAVETODISK){     //option to save EigenSym
         printf("saving eigensym to directory %s\n", eigenDir.c_str());
         eigen->saveToDiskBINARY(eigenDir);  //Save the EigenSym object for later use in the directory eigenDir
      }else{
         printf("proceeding without saving eigensym\n");
      }
   }
   else{
      //computeEigenType==LOADFROMDISK 
      printf("loading eigensym from directory %s\n", eigenDir.c_str());
      eigen->readFromDiskBINARY(eigenDir);   //Load an existing EigenSym object from the directory eigenDir
   }
   

   //run the GWAS SNP testing
   if (runGwasType==RUN){
      
      int n_individuals= eigen->n;
      ResultsLMM *results= new ResultsLMM();

      //generate the covariates as an array of ones
      int n_covariates=1;
      double *covariates = new double[n_individuals*n_covariates];
      for (int i=0; i<n_individuals; ++i){
         covariates[i]=1.0;
      }

      //load the phenotypes
      printf("loading the test data\n");
      CSnpData *snpData=new CSnpData();
      snpData->ReadSnpFile(snpDataFile);
      int n_SNPs = (int)snpData->cRows;
      printf("preprocessing the test data\n");
      snpData->CreateColumnMajorArrayFromVectorsAndMapping(eigen->columnLabelOrder);
      
      //perform the preprocessing on the SNPs to be tested.
      ImputeAndZeroMeanSNPs(snpData->snpArray, n_individuals, n_SNPs, 2.0, true);//Eigenstrat Imputation

      //load  and preprocess the SNPs
      CPhenData *phenData = new CPhenData();
      phenData->ReadPhenFile(phenDataFile);
      phenData->CreateColumnMajorArrayFromVectorsAndMapping(eigen->columnLabelOrder);

      //run he Gwas given the data
      if (eigen->n>eigen->rankU){
         runGWAS_lowrank(eigen, results,snpData,phenData,covariates,n_covariates);
      }else{
         runGWAS_fullrank(eigen, results,snpData,phenData,covariates,n_covariates);
      }

      results->WriteResultFile(outputFile);
   }
   Warn( "LmmGwas() done.");   
}

void LmmGwas_oldfortiming()
   {
   clock_t startLoading,startPreprocessing,startSVD,saveSVD,loadSVD,startTransform,startNull,startAlt,finished;
   startLoading=clock();

   // Read all the input files
   CSnpData snpData;
   snpData.ReadSnpFile( snpDataFile );
   snpData.CreateColumnMajorArrayFromVectors();

   CPhenData phenData;
   phenData.ReadPhenFile( phenDataFile );
   phenData.CreateColumnMajorArrayFromVectorsAndMapping( snpData.columnLabels );

   CKernelData kernelData;
   if ( kernelFile.size()>0 )
      {
      kernelData.ReadKernelFile( kernelFile );
      kernelData.CreateColumnMajorArrayFromVectorsAndMapping( snpData.columnLabels );
      }

   // We need this second copy to build the covariance matrix for the low rank model:
   CSnpData covarianceSNPs;
   
   //preprocess the data:
   int n_individuals=(int)phenData.cColumns;
   int n_SNPs=(int)snpData.cRows;
   int n_phenotypes=(int)phenData.cRows;
   int rankU = n_individuals; //=(int)covarianceSNPs.cRows;
   int ncovariates=1;         //assume that there is only the bias term
   double *Kernel=0;          //null pointer
   string kernelID;
   if (kernelFile.size()>0){
      kernelID="Kernel";
      Kernel=kernelData.kernelArray;
      printf("\nPreprocessing the data.\n");
      startPreprocessing=clock();
   }else{
      kernelID="covarianceKernel";
      // We need these SNPs to build a covariance matrix for the low rank model:
      covarianceSNPs.ReadSnpFile( snpDataFile );
      covarianceSNPs.CreateColumnMajorArrayFromVectors();
      printf("\nPreprocessing the data.\n");
      rankU = min(rankU,(int)covarianceSNPs.cRows);
      //perform the preprocessing on the covariance SNPs
      startPreprocessing=clock();
      ImputeAndZeroMeanSNPs(covarianceSNPs.snpArray, n_individuals, (int)covarianceSNPs.cRows, 2.0, true);     //Eigenstrat Imputation for the kernel SNPs
      preprocessSNPsForBuildingCovariance(covarianceSNPs.snpArray, n_individuals, (int)covarianceSNPs.cRows);  //Zero mean each individual and divide by sqrt(N_SNPs)
      if((int)covarianceSNPs.cRows*2>=n_individuals ||forceFullRank){//do the full rank model
         //compute the covariance matrix as the inner product matrix of the Covariance SNPs
         Kernel= new double[n_individuals*n_individuals];
         int nCols=(int)covarianceSNPs.cRows;
         BlasLapack::syrk(&char_L,&char_N,&n_individuals,&nCols,&double_one, covarianceSNPs.snpArray, &n_individuals, &double_zero,  Kernel, &n_individuals);
      }
   }

   //perform the preprocessing on the SNPs to be tested.
   ImputeAndZeroMeanSNPs(snpData.snpArray, n_individuals, n_SNPs, 2.0, true);//Eigenstrat Imputation
   
   //Create the Bias Vector:
   int n_covariates=1;
   double *covariates=new double[n_individuals];
   for(int ind=0;ind<n_individuals;++ind){
      covariates[ind]=1.0;
   }

   // Do the work
   if (n_phenotypes!=1){
      Warn( "LmmGwas() for multiple phenotypes is Not Yet Implemented. Will only use the first Phenotype in the following.\a" );            // send the message with a 'beep' (the \a)
   }

   ResultsLMM result;   //the result variables

   if(Kernel!=0){//run the full rank model
      printf("Running the full rank code.\n");
      printf("Computing the SVD.\n");
      //Do the association mapping:
      startSVD=clock();
      EigenSym *eigen= new EigenSym();
      eigen->set_fromKernel(kernelID, Kernel,n_individuals,n_individuals, snpData.columnLabels);

      saveSVD=clock();
      printf("Saving the SVD to disk.\n");
      eigen->saveToDiskBINARY("svd_test");
      delete eigen;

      loadSVD=clock();
      printf("Loading the SVD from disk.\n");
      EigenSym *eigenReload=new EigenSym();
      eigenReload->readFromDiskBINARY("svd_test");

      FastLMM_fullrank LMM_fullrank = FastLMM_fullrank(eigenReload);
      //LMM_fullrank.set_Kernel(kernelData.kernelArray);//use the kernel provided in the file
      
      //LMM_fullrank.set_Kernel(Kernel);//use the covariance matrix
      printf("Transforming the data.\n");
      startTransform=clock();
      LMM_fullrank.set_Phenotype(phenData.phenArray, phenData.rowLabels[0]);
      LMM_fullrank.set_SNPs(snpData.snpArray,n_SNPs,snpData.rowLabels);
      LMM_fullrank.set_covariates(covariates,n_covariates);
      printf("Starting the training.\n");
      startNull=clock();
      LMM_fullrank.trainNull();
      startAlt=clock();
      if (runAlternativeModels){
         LMM_fullrank.trainAlt();//This takes a while.
         //printf("LL_alt: ");for (int i=0;i<n_SNPs;++i)printf("%f ",LMM_fullrank.LL[i]);printf("\n");
      }
      result=LMM_fullrank.getResults();
   }else{//run the low rank model
      //Do the association mapping:
      printf("Fewer SNPs than individuals in Covariance Matrix. Running the low rank code.\n");
      startSVD=clock();
      EigenSym *eigen=new EigenSym();
      eigen->set_covFromX(kernelID, covarianceSNPs.snpArray, n_individuals, rankU, snpData.columnLabels);

      saveSVD=clock();
      printf("Saving the SVD to disk.\n");
      eigen->saveToDiskBINARY("svd_test");
      delete eigen;

      loadSVD=clock();
      printf("Loading the SVD from disk.\n");
      EigenSym *eigenReload=new EigenSym();
      eigenReload->readFromDiskBINARY("svd_test");
      
      FastLMM_lowrank LMM_lowrank = FastLMM_lowrank(eigenReload);
      startTransform=clock();
      printf("Transforming the data.\n");
      LMM_lowrank.set_Phenotype(phenData.phenArray, phenData.rowLabels[0]);
      LMM_lowrank.set_SNPs(snpData.snpArray,n_SNPs,snpData.rowLabels);
      LMM_lowrank.set_covariates(covariates,n_covariates);
      printf("Starting the training.\n");
      startNull=clock();
      LMM_lowrank.trainNull();
      startAlt=clock();
      if (runAlternativeModels){
         LMM_lowrank.trainAlt();//This takes a while.
         //printf("LL_alt: ");for (int i=0;i<n_SNPs;++i)printf("%f ",LMM_fullrank.LL[i]);printf("\n");
      }
      //printf("LL_null: %f\n",LMM_fullrank.LL_0);
      result=LMM_lowrank.getResults();
      //printf("qvals: ");for(int i=0;i<n_SNPs;++i){  printf("%f ",qvals[i]); }printf("\n");
   }
   
   finished=clock();
   printf("Time to load the data:                  %.2fs\n", (1.0*(startPreprocessing-startLoading)) / (1.0*CLOCKS_PER_SEC));
   printf("Time to preprocess the data:            %.2fs\n", (1.0*(startSVD-startPreprocessing)) / (1.0*CLOCKS_PER_SEC));
   printf("Time to compute the SVD:                %.2fs\n", (1.0*(saveSVD-startSVD)) / (1.0*CLOCKS_PER_SEC));
   printf("Time to save the SVD to disk:           %.2fs\n", (1.0*(loadSVD-saveSVD)) / (1.0*CLOCKS_PER_SEC));
   printf("Time to load the SVD from disk:         %.2fs\n", (1.0*(startTransform-loadSVD)) / (1.0*CLOCKS_PER_SEC));
   printf("Time to transform the data:             %.2fs\n", (1.0*(startNull-startTransform)) / (1.0*CLOCKS_PER_SEC));
   printf("Time to compute null model:             %.2fs\n", (1.0*(startAlt-startNull)) / (1.0*CLOCKS_PER_SEC));
   printf("Time to compute alternative models:     %.2fs\n", (1.0*(finished-startAlt)) / (1.0*CLOCKS_PER_SEC));
   printf("Total Time:                             %.2fs\n", (1.0*(finished-startLoading)) / (1.0*CLOCKS_PER_SEC));

   //Save the results:
   result.WriteResultFile(outputFile);


   Warn( "LmmGwas() results have been saved to %s.\a\n",outputFile.c_str() );   
};

int computeEigenSym_Kernel(string kernelFile, EigenSym *eigen){
   printf("computing full rank EigenSym from input Kernel matrix\n");
   CKernelData kernelData;
   if ( kernelFile.size()>0 )
      {
      kernelData.ReadKernelFile( kernelFile );
      kernelData.CreateColumnMajorArrayFromVectorsAndMapping( kernelData.columnLabels );
      }
   eigen->set_fromKernel(kernelFile, kernelData.kernelArray, (int)kernelData.cColumns, (int)kernelData.cColumns, kernelData.columnLabels);
   return 0;
}
      
int computeEigenSym_SnpSample(string snpDataCovFile, EigenSym *eigen){
   printf("Reading the data for Cov building\n");
   CSnpData covarianceSNPs;
   // We need these SNPs to build a covariance matrix for the low rank model:
   covarianceSNPs.ReadSnpFile( snpDataFile );
   covarianceSNPs.CreateColumnMajorArrayFromVectors();
   printf("Preprocessing the data for Cov building.\n");
   
   //perform the preprocessing on the covariance SNPs
   ImputeAndZeroMeanSNPs(covarianceSNPs.snpArray, (int)covarianceSNPs.cColumns, (int)covarianceSNPs.cRows, 2.0, true);     //Eigenstrat Imputation for the kernel SNPs
   preprocessSNPsForBuildingCovariance(covarianceSNPs.snpArray, (int)covarianceSNPs.cColumns, (int)covarianceSNPs.cRows);  //Zero mean each individual and divide by sqrt(N_SNPs)
   
   if((int)covarianceSNPs.cRows>=(int)covarianceSNPs.cColumns || forceFullRank){//do the full rank model
      printf("computing full rank EigenSym\n");
      int rankU = (int)covarianceSNPs.cColumns;
      int n_individuals=(int)covarianceSNPs.cColumns;
      
      //compute the covariance matrix as the inner product matrix of the Covariance SNPs
      double *Kernel= new double[(int)covarianceSNPs.cColumns*(int)covarianceSNPs.cColumns];
      int nCols=(int)covarianceSNPs.cRows;
      BlasLapack::syrk(&char_L,&char_N,&n_individuals,&nCols,&double_one, covarianceSNPs.snpArray, &n_individuals, &double_zero,  Kernel, &n_individuals);

      //compute the Eigenvectors and set EigenSym from the kernel
      eigen->set_fromKernel(snpDataCovFile, Kernel,n_individuals,n_individuals, covarianceSNPs.columnLabels);
      delete[] Kernel;  //this is no longer needed
   }else{
      printf("computing low rank EigenSym\n");
      int rankU = (int)covarianceSNPs.cRows;
      int n_individuals=(int)covarianceSNPs.cColumns;

      //compute the Eigenvectors and set EigenSym directly from the covarianceSNPs
      eigen->set_covFromX(snpDataCovFile, covarianceSNPs.snpArray,n_individuals,rankU,covarianceSNPs.columnLabels);
   }
   
   return 0;
}
      
int runGWAS_lowrank(EigenSym *eigen, ResultsLMM *results, CSnpData *snpData, CPhenData *phenData, double *covariates, int n_covariates){
   printf("starting lowrank training\n");
   int n_individuals =eigen->n;
   int n_Snps =(int)snpData->cRows;
   int rankU=eigen->rankU;

   FastLMM_lowrank LMM_lowrank = FastLMM_lowrank(eigen);
   LMM_lowrank.set_Phenotype(phenData->phenArray, phenData->rowLabels[0]);
   LMM_lowrank.set_SNPs(snpData->snpArray, n_Snps, snpData->columnLabels);
   LMM_lowrank.set_covariates(covariates, n_covariates);

   LMM_lowrank.trainNull();

   if (runAlternativeModels){
      LMM_lowrank.trainAlt();
   }

   //get the results
   *results=LMM_lowrank.getResults();
   printf("done lowrank training\n");
   return 0;
}

int runGWAS_fullrank(EigenSym *eigen, ResultsLMM *results, CSnpData *snpData, CPhenData *phenData, double *covariates, int n_covariates){
   printf("starting fullrank training\n");
   int n_individuals=eigen->n;
   int n_SNPs= (int)snpData->cRows;

   FastLMM_fullrank LMM_fullrank = FastLMM_fullrank(eigen);
 
   LMM_fullrank.set_Phenotype(phenData->phenArray, phenData->rowLabels[0]);
   LMM_fullrank.set_SNPs(snpData->snpArray, n_SNPs, snpData->rowLabels);
   LMM_fullrank.set_covariates(covariates,n_covariates);
   
   LMM_fullrank.trainNull();
   if (runAlternativeModels){
      LMM_fullrank.trainAlt();
   }

   //get the results
   *results=LMM_fullrank.getResults();
   printf("done fullrank training\n");
   return 0;
}
