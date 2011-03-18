#include "FastLmmC.h"
//#include "FastLMM.h"
//#include "Utils.h"

const bool FastLMM::brentExhaustive=false;

ResultsLMM::ResultsLMM(){}

ResultsLMM::ResultsLMM(string modelID, vector<string> snpLabels, string phenLabel, int n_individuals, int ncovariates, int *index, double *LL_alt, double *LL_0, double *pvals, double *snpWeights_alt, double *bias_alt, double *sigg2_alt, double *sige2_alt, double *delta_alt, double *bias_0, double *sige2_0, double *sigg2_0, double *delta_0, double *qvals){
   this->snpLabels = snpLabels;           //identifier of the SNP
   this->phenLabel = phenLabel;           //the name of the phenotype
   this->ncovariates = ncovariates;       //the number of covariates in the study
   this->n_individuals = n_individuals;   //the number of individuals in the study
   this->index = index;                   //The results will be printed out, sorted by index
   this->LL_alt = LL_alt;                 //Log-likelihoods of the alternative models
   this->LL_0 = LL_0;                     //Log-likelihoods of the null models
   this->pvals = pvals;                   //p-values
   this->snpWeights_alt = snpWeights_alt; //the SNP weights of the alternative models
   this->bias_alt = bias_alt;             //the bias weights of the alternative models
   this->sigg2_alt = sigg2_alt;           //the genetic variance of the alternative models
   this->sige2_alt = sige2_alt;           //the random vriance of the alternative models
   this->delta_alt = delta_alt;           //the fraction of sige2_alt/sigg2_alt
   this->bias_0 = bias_0;                 //the bias terms of the null models
   this->sigg2_0 = sigg2_0;               //the genetic variance of the null models
   this->sige2_0 = sige2_0;               //the random variance of the null models
   this->delta_0 = delta_0;               //fraction of sige2_0/sigg2_0
   this->qvals = qvals;                   //q-values (minimum FDR, under which the current SNP would be significant)
}

void ResultsLMM::WriteResultFile( string& resultFile)
{
   /*
    * The output file we write is a tab separated file with the form
    *   rowIndex rowCount nullIndex SNPId ... <eol>
    *   NameYyy x x ...<eol>
    *   NameYyz x x ...<eol>
    *   ... <eol>
    *   <eof>
    *
    * The first row labels the columns
    * each subsequent row contains the data as an Attribute Name followed by a value separated by <tab> for each column
    */
   string columnLabels="rowIndex\trowCount\tnullIndex\tSNPId\tphenId\tlogLikeAlt\tlogLikeNull\tpValue\tSnpWeight\tbias\tgeneticVarSq\tresidualVarSq\tdelta\tnullBias\tnullGeneticVar\tnullResidualVar\tnullDelta\tnumCidsRemoved\tqValue\n";
   FILE *pFile;

   if ( resultFile.size() == 0 )
      {
      Fatal( "No output filename to open " );
      }

   pFile = fopen( resultFile.c_str(), "wt" );          // write ascii

   if ( !pFile )
      {
      char szFullPath[ _MAX_PATH ];
         
      _fullpath( szFullPath, resultFile.c_str(), sizeof( szFullPath ) );
      Fatal( "Cannot open output file [%s].  \n  CRT Error %d: %s", szFullPath, errno, strerror( errno ) );
      }

   size_t cRows = snpLabels.size();

   unsigned int old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);   // to be compatible, use 2 digit exponents

   // Write the first row (var\t<columnLabels><eol>)
   fprintf( pFile, columnLabels.c_str() );
   
   // write the rest of the rows
   for( size_t i=0; i<cRows; ++i )
   {
      fprintf(pFile,"%i\t",index[i]);                    //row index
      fprintf(pFile,"%i\t",cRows);                       //number of rows
      fprintf(pFile,"-1\t");                             //null  index is always -1
      fprintf(pFile,"%s\t",snpLabels[index[i]].c_str()); //the snp label
      fprintf(pFile,"%s\t",phenLabel.c_str());           //the phenotype name
      fprintf(pFile,"%.13E\t",LL_alt[index[i]]);         //The likelihood of the alternative model
      fprintf(pFile,"%.13E\t",LL_0[index[i]]);           //the likelihood of the null model
      fprintf(pFile,"%.13E\t",pvals[index[i]]);          //pvalue
      fprintf(pFile,"%.13E\t",snpWeights_alt[index[i]]); //SNP weight
      fprintf(pFile,"%.13E\t",bias_alt[index[i]]);       //bias
      fprintf(pFile,"%.13E\t",sigg2_alt[index[i]]);      //genetic variance squared
      fprintf(pFile,"%.13E\t",sige2_alt[index[i]]);      //residual variance squared
      fprintf(pFile,"%.13E\t",delta_alt[index[i]]);      //delta
      fprintf(pFile,"%.13E\t",bias_0[index[i]]);         //nullBias
      fprintf(pFile,"%.13E\t",sigg2_0[index[i]]);        //nullGeneticVar
      fprintf(pFile,"%.13E\t",sige2_0[index[i]]);        //nullResidualVar
      fprintf(pFile,"%.13E\t",delta_0[index[i]]);        //nullDelta
      fprintf(pFile,"0\t");                              //numCidsRemoved
      fprintf(pFile,"%.13E\n",qvals[index[i]]);          //qValue
   }

   _set_output_format( old_exponent_format );            // restore the default format behavior
   fclose( pFile );
}
