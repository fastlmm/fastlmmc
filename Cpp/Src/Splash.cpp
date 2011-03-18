/*
 * Copyright (c) 2010 Microsft Corporation -- All rights reserved
 *
 * Microsoft makes no representations about the fitness nor
 * suitability of this software for any purpose.  It is provided
 * "as is" with neither expressed nor implied warranty.
 *
 * File: Splash.cpp
 * Date: 05/01/2010
 * Owner: BobD
 * Description: This file defines the version string and uses the
 *		  __DATE__ and __TIME__ macros.  We MUST compile this
 *		  file for every build for proper operation.  This is 
 *                accomplished by having a 'pre-link step' in the 
 *                build that runs the following command:
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

int wVerMajor = 1;
int wVerMinor = 0;
//int wVerBuild = 0;       // Now set in BuildVer.h (above)
int wVerPatch = 0;

char *szSplash = "\nFastLmmC %s - Fast Linear Mixed Model"
      "\n  Copyright Microsoft Corporation"
      "\n  Compiled " __DATE__ " at " __TIME__ " by " szBuildHost 
#if defined( USE_MKL )
      " using MKL"
#elif defined( USE_ACML )
      " using ACML"
#endif
      "\n";

