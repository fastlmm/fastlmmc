FaST-LMM, (Factored Spectrally Transformed Linear Mixed Models) is 
a program for performing genome-wide association studies (GWAS) on 
large data sets. It runs on both Windows and Linux systems, and has 
been tested on data sets with over 120,000 individuals.

This file is located in the root directory for the FaSTLMM source
tree.  To get started quickly with a precompiled version, see the 
user's manual in the Doc directory.

All the sources for the C++ versions are included in this tree.
However, there are additional external dependencies that you will 
need to download and install to build the code.

This document describes the directories within the tree and their 
usage/function.  The directories and their functions are:

  Cpp -
      The C++ project is contained within the Cpp tree.  All the 
      files for editing and starting a build with Visual Studio 
      are present.

  Data -
      The Data directory holds files for running simple tests
      to confirm proper installation of FaSTLMM.  The data may
      also be used as a sanity test after re-building from 
      sources or downloading copy of the program

  Doc - 
      All documentation files for the project are located under the 
      Doc directory

  Externals -
      All external dependencies used under their own license
      are located under the Externals directory
	  
  Test
      Unit Tests that drive the different FaSTLMM executables 
      and validate their behavior are located under this directory.
      
	  
---------- Building the code ----------
In addition to the code in this tree, the following external dependencies
that must be installed to build FaSTLMM.

Building for Windows.
  The C++ versions require Visual Studio 2010 (VS)- 
      Any version of Visual Studio (Express through Universal) is 
      capable of building FaSTLMM.  If you do not already have a copy
      of Visual Studio, the Visual Studio 2010 Express edition can be 
      freely downloaded from 
        http://www.microsoft.com/express/downloads/

  For the C++ version 
    FastLmmC uses a 3rd party math library for advanced math funtions
    and performance.  FastLmmC can use either Intel's MKL or AMD's 
    ACML math libraries.  Once you have installed the appropriate 
    library, use the Visual Studio IDE to open the solution file
    Cpp\FastLmmC.sln and select the configuration matching your
    library from the solution and build.  ACML requires an 
    additional step to tell Visual Studio where it is located.  
    You must set the environment variable ACML_ROOT to point to 
    your install location or libraries - for example
	  C>set ACML_ROOT=C:\AMD\acml4.4.0
	You can find more about the math libraries at their respective
	web sites:
    http://software.intel.com/en-us/articles/intel-mkl/
    http://developer.amd.com/libraries/acml/pages/default.aspx

  Building the C++ version for Linux
    FastLmmC uses a 3rd party math library for advanced math functions 
    and performance.  The program has been run on Ubuntu Linux and can 
    use either Intel's MKL or AMD's ACML math libraries for Linux.  Once 
    you have selected and installed the appropriate library, you can 
    then build using the appropriate script file located in the Cpp 
    directory.  Review of the two files, DoMKL_linux and DoAcml_linux, 
    will show very simple scripts to compile the program using g++ and 
    then link the .o files with the appropriate math library.  The *.o 
    files are written to version specific directories, so it is necessary 
    to create the appropriate directory prior to running the script.  
    For more details, see the script.

----------
KNOWN ISSUES
----------
