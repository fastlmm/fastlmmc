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

/*
 * Splash - {Splash screen information}
 *
 *         File Name:   Splash.cpp
 *           Version:   1.00
 *            Author:   
 *     Creation Date:   18 Nov 2010
 *     Revision Date:   18 Nov 2010
 *
 *    Modlue Purpose:   This file defines the version string and uses 
 *                      the__DATE__ and __TIME__ macros.  We MUST 
 *                      compile this file for every build for proper 
 *                      operation.  This is accomplished by having a 
 *                      'pre-link step' in the build that runs the 
 *                      following command:
 *
 *     ShowTime // Generated file **DO NOT EDIT**\nint wVerBuild = $y$J; \n#define szBuildHost \"$H\"\n >BuildVer.h
 */

/*
 * Requires the following to be visible
extern int wVerMajor;
extern int wVerMinor;
extern int wVerBuild;
extern int wVerPatch;
extern char *szSplash;
 */

#include "BuildVer.h"      // get the updated build / host information (generated during build in VisualStudio)

int wVerMajor = 2;         // 2.06 -> transition to VS2012 & MKL 11.0
int wVerMinor =  8;        // odd number is 'beta' / even number is 'release'
                           // 2.08 -> 2026 revival: oneMKL 2026.2.1, g++ 11 / VS2022 (2.07 was the 2014 release)
//int wVerBuild = 0;       // Now set in BuildVer.h (above)
int wVerPatch =  0;

#if defined ( _MSC_VER )   // COMPILE_FOR_WINDOWS using Visual C compiler
#define  myRuntimeHost  "Windows"
#elif defined( __GNUC__ )  // COMPILE_FOR_LINUX using g++
#if defined( USE_LINUX_V2x_KERNEL )
#define  myRuntimeHost  "Linux v2.x kernel"
#else
#define  myRuntimeHost  "Linux v3.x kernel"
#endif
#else
#error We only support VisualC++ on Windows and GCC on Linux.  Unrecognized compiler environment.  Stop Compilation
#endif

const char *szSplash = "\nFastLmmC %s - Factored Spectrally Transformed Linear Mixed Models"
#if defined( _DEBUG )
      " [Debug]"
#else
      " [Release]"
#endif
      "\n  Copyright Microsoft Corporation"
#if !defined( COMMERCIAL_USE )
      " -- Licensed Only for Non-Commercial use."
#endif
      "\n  Compiled " __DATE__ " at " __TIME__ " by " szBuildHost " for " myRuntimeHost
      "\n  using %s"
      "\n";

