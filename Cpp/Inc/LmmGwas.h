#if !defined(LmmGwas_h)
#define LmmGwas_h
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

void LmmGwas( StudyData *data, EigenSym *eigen, ResultsLMM *results );

int computeEigenSym_SnpSample( const string& id, StudyData *data, EigenSym *eigen );
int runGWAS_lowrank( EigenSym *eigen, ResultsLMM *results, StudyData *data );
int runGWAS_woodbury( EigenSym *eigen, ResultsLMM *results, StudyData *data );
int runGWAS_fullrank( EigenSym *eigen, ResultsLMM *results, StudyData *data );
int fillEigenSym( const string& id, EigenSym *eigen, StudyData *data );
#endif   //LmmGwas
