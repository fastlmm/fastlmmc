#if !defined( GroupSim_h )
#define GroupSim_h
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

class GroupSim
   {
   void loadIndividualsIndexFromPhen(AlternatePhenotypeParameters& alternatePhenotypeParameters);
   void loadGroupIndex(const string& groupFile);
public:
   CKernelData *kernelGroups;    //group-by-group kernel
   map<string,size_t> indMap;    //The group of each individual (individualID to groupIndex)
   map<string,size_t> groupIndexMap;   //the groupIndex of each group
   vector<string> groupIds;      //size: cGroups; sorted by groupIndex
   vector<size_t> groupSizes;
   size_t cGroups;         //number of groups
   size_t cIndividuals;    //number of individuals in all groups
   vector<size_t> Z;       //the Z-matrix, in integer form (size: cIndividuals)
   vector<string> individualsIndex; //contains the unique identifiers for each individual
   GroupSim(AlternatePhenotypeParameters& alternatePhenotypeParameters, const string& groupFile);
   ~GroupSim();
   
   void loadKernel(const string& simGroups);
   };


#endif     //GroupSim_h
