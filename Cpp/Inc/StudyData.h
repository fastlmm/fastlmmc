#if !defined( StudyData_h )
#define StudyData_h
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
 
class StudyData
   {
public:
   size_t cIndividuals;          // count of individuals
   size_t cPhenotypes;           // count of phenotypes (1 in PLink)
   size_t cSnps;                 // count of SNPs
   size_t cCovariates;           // count of covariate (minimum 1 for bias term)
   size_t cSnpsCovariance;       // count of SNPs used to estimate covariance
   size_t ind_phenotype;         // index of the phenotype to analyze
   vector<string> individualLabels;    // Name/Id for each individual
   vector<string> phenotypeLabels;     // Name/Id for each phenotype
   vector<SnpInfo> snpInfo;            // Associated SNP info for each SNP
   vector<SnpInfo> snpInfoCovariance;  // Associated SNP info for each SNP used to estimate covariance
   vector<string> covariateLabels;     // Name/Id for each covariate
   real *phenotypeData;          // size cIndividuals * cPhenotypes
   real *snpData;                // size cIndividuals * cSnps
   real *covariatesData;         // size cIndividuals * cCovariates
   real *kernelData;             // size cIndividuals * cIndividuals
   real *snpDataCovariance;      // size cIndividuals * cSnps
   GroupSim *groupSim;           // similarity matrix for groups
   size_t *permutation;          // size nperm * cIndividuals
   size_t nperm;                 // number of permutations

   StudyData();                  // constructor/initializer
   ~StudyData();                 // destructor must clean up allocations

   void GeneratePermuations( size_t nperm );             // generate random permutations
   void Dump( FILE *pf, const string& comment ) const;   // DEBUG support to Dump StudyData to file
   };

#endif //StudyData_h
