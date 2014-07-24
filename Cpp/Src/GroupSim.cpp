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

#include "FastLmmC.h"

GroupSim::GroupSim(AlternatePhenotypeParameters& alternatePhenotypeParameters, const string& groupFile)
   {
   kernelGroups = nullptr;
   this->loadIndividualsIndexFromPhen(alternatePhenotypeParameters);
   this->loadGroupIndex(groupFile);
   }

GroupSim::~GroupSim()
   {
   delete []kernelGroups;
   }

void GroupSim::loadKernel(const string& simGroups)
   {
   this->kernelGroups = new CKernelData();
   this->kernelGroups->ReadKernelFile(simGroups);
   this->kernelGroups->CreateColumnMajorArrayFromVectorsAndMapping(this->groupIds);
   }

void GroupSim::loadIndividualsIndexFromPhen(AlternatePhenotypeParameters& alternatePhenotypeParameters)
   {
   CPlinkAlternatePhenotypeFile phen = CPlinkAlternatePhenotypeFile(alternatePhenotypeParameters.alternatePhenotypeFilename);   
   phen.Read();

   if ( !alternatePhenotypeParameters.selectedPhenotypeId.empty() )
      {
      // lookup the phenotype name and set the index
      alternatePhenotypeParameters.selectedPhenotypeIndex = phen.IndexFromPhenotypeLabel( alternatePhenotypeParameters.selectedPhenotypeId );
      Verbose( verboseMinor, "             Using Phenotype Name: [%s]", alternatePhenotypeParameters.selectedPhenotypeId.c_str() );
      }
   if ( !(alternatePhenotypeParameters.selectedPhenotypeIndex < phen.CountOfPhenotypes()) )
      {
      Fatal( "alternatePhenotypeIndex of %d is out of range.  Max is %d.", 1+alternatePhenotypeParameters.selectedPhenotypeIndex, phen.CountOfPhenotypes() );
      }
   Verbose( verboseMinor, "            Using Phenotype Index: %7d", 1+alternatePhenotypeParameters.selectedPhenotypeIndex );
   //get the individualsMapping
   this->individualsIndex = phen.KeysForIndividualsWithValidPhenotypes(alternatePhenotypeParameters.selectedPhenotypeIndex);
   this->cIndividuals=this->individualsIndex.size();
   }

void GroupSim::loadGroupIndex(const string& groupFile)
   {
   //this->cIndividuals = 0;
   this->cGroups=0;
   CPlinkLexer lex( groupFile );
   CToken tok;
   lex.NextToken( tok );
   while ( tok.type != tokEOF )        // for each line until EOF
      {
      if ( tok.type == tokEOL )        // blank line...or comment line or ???
         {
         lex.NextToken( tok );
         continue;
         }
      string idFamily;
      string idIndividual;
      string idGroup;
      size_t indexGroup;
      
      lex.ExpectId( tok, idFamily, "FamilyID" );
      //famIDs.push_back(idFamily);
      lex.ExpectId( tok, idIndividual, "IndividualID" );
      //indIDs.push_back(idIndividual);
      
      string currentIndId;

      if ( individualIdType==INDIVIDUALID_ONLY )
         {
         currentIndId=idIndividual;
         }
      else if ( individualIdType==FAMILYID_INDIVIDUALID )
         {
         currentIndId = KeyFromIdFamilyAndIdIndividual( idFamily, idIndividual );
         }
      else
         {
         Fatal( "Illegal individualIdType specified :[%i]",individualIdType );
         }

      

      if(indMap.count(idIndividual)){Fatal("\nIndividual %s appears twice in GroupFile",idIndividual.c_str());}
      
      size_t cColumns_=0;
      
      while( tok.type != tokEOL && tok.type!=tokEOF)
         {
         ++cColumns_;
         lex.ExpectId( tok, idGroup , "IndividualID");
         }
      if (cColumns_>1)
         {
         Fatal("\ngroupFile has too many columns. (%i)", cColumns_);
         }

      if (!groupIndexMap.count(idGroup))
         {
         groupIndexMap[idGroup]=cGroups;
         groupIds.push_back(idGroup);
         groupSizes.push_back(0);
         ++cGroups;
         }
      indexGroup=groupIndexMap[idGroup];
      indMap[idIndividual]=indexGroup;
      }

   //filter out the used individuals
   for(size_t i=0; i<this->individualsIndex.size();++i)
      {
      if(indMap.count(individualsIndex[i])!=1)
         {
         Fatal( "Cannot create group mapping. kernelFile contains [%i] copies of mapping key [%s] (epected 1 copy)", indMap.count(individualsIndex[i]), individualsIndex[i].c_str() );
         }
      size_t iGroup = indMap[individualsIndex[i]];
      Z.push_back(iGroup);
      ++groupSizes[iGroup];
      }
   }
