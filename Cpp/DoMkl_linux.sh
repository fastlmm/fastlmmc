#!/bin/bash

if [ "$1" = "v2.x" ]; then
  LINUX_KERNEL=-DUSE_LINUX_V2x_KERNEL
  echo "Compiling for RedHat Enterprise Linux with v2.x kernel (LINUX_KERNEL=$LINUX_KERNEL)"
else
  LINUX_KERNEL=
  echo "Compiling for Linux with V3.x kernel (LINUX_KERNEL=$LINUX_KERNEL)"
fi

echo "Generating Inc/BuildVer.h"
echo "// Generated File  ** DO NOT EDIT **" >Inc/BuildVer.h
echo "int wVerBuild = `date +%Y%m%d`;"     >>Inc/BuildVer.h
echo "#define szBuildHost \"`uname -n`\""  >>Inc/BuildVer.h

echo "Compiling sources with MKL Library"
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/AutoSelectSnps.cpp       -o ./obj-mkl/AutoSelectSnps.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/Beta.cpp                 -o ./obj-mkl/Beta.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/BrentC.cpp               -o ./obj-mkl/BrentC.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CCovariatesData.cpp      -o ./obj-mkl/CCovariatesData.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CKernelData.cpp          -o ./obj-mkl/CKernelData.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CLexer.cpp               -o ./obj-mkl/CLexer.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CLogisticRegression.cpp  -o ./obj-mkl/CLogisticRegression.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CPhenData.cpp            -o ./obj-mkl/CPhenData.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CPlinkAlternatePhenotypeFile.cpp -o ./obj-mkl/CPlinkAlternatePhenotypeFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CPlinkBedFile.cpp        -o ./obj-mkl/CPlinkBedFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CPlinkBigFile.cpp        -o ./obj-mkl/CPlinkBigFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CPlinkBimFile.cpp        -o ./obj-mkl/CPlinkBimFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CPlinkDatFile.cpp        -o ./obj-mkl/CPlinkDatFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CPlinkFamFile.cpp        -o ./obj-mkl/CPlinkFamFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CPlinkFile.cpp           -o ./obj-mkl/CPlinkFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CPlinkLexer.cpp          -o ./obj-mkl/CPlinkLexer.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CPlinkMapFile.cpp        -o ./obj-mkl/CPlinkMapFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CSnpData.cpp             -o ./obj-mkl/CSnpData.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/CSnpFilterOptions.cpp    -o ./obj-mkl/CSnpFilterOptions.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/EigenStratImputation.cpp -o ./obj-mkl/EigenStratImputation.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/EigenSym.cpp             -o ./obj-mkl/EigenSym.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/FastLMM.cpp              -o ./obj-mkl/FastLMM.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/FastLMM_fullrank.cpp     -o ./obj-mkl/FastLMM_fullrank.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/FastLMM_lowrank.cpp      -o ./obj-mkl/FastLMM_lowrank.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/FastLMM_woodbury.cpp     -o ./obj-mkl/FastLMM_woodbury.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/FastLmmC.cpp             -o ./obj-mkl/FastLmmC.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/FisherF.cpp              -o ./obj-mkl/FisherF.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/Gamma.cpp                -o ./obj-mkl/Gamma.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/GroupSim.cpp             -o ./obj-mkl/GroupSim.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/LinReg.cpp               -o ./obj-mkl/LinReg.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/LmmGwas.cpp              -o ./obj-mkl/LmmGwas.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/MathFunctions.cpp        -o ./obj-mkl/MathFunctions.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/MatrixIO.cpp             -o ./obj-mkl/MatrixIO.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/Qval.cpp                 -o ./obj-mkl/Qval.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/Splash.cpp               -o ./obj-mkl/Splash.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/Startup.cpp              -o ./obj-mkl/Startup.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/StudyData.cpp            -o ./obj-mkl/StudyData.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/TsvLexer.cpp             -o ./obj-mkl/TsvLexer.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64 Src/Utils.cpp                -o ./obj-mkl/Utils.o

echo "Linking fastlmmc_mkl"
g++ -o ./fastlmmc_mkl		\
	obj-mkl/AutoSelectSnps.o	\
	obj-mkl/Beta.o			\
	obj-mkl/BrentC.o		\
	obj-mkl/CCovariatesData.o	\
	obj-mkl/CKernelData.o	\
	obj-mkl/CLexer.o		\
	obj-mkl/CLogisticRegression.o   \
	obj-mkl/CPhenData.o		\
	obj-mkl/CPlinkAlternatePhenotypeFile.o \
	obj-mkl/CPlinkBedFile.o	\
	obj-mkl/CPlinkBigFile.o	\
	obj-mkl/CPlinkBimFile.o \
	obj-mkl/CPlinkDatFile.o	\
	obj-mkl/CPlinkFamFile.o \
	obj-mkl/CPlinkFile.o	\
	obj-mkl/CPlinkLexer.o	\
	obj-mkl/CPlinkMapFile.o	\
	obj-mkl/CSnpData.o		\
	obj-mkl/CSnpFilterOptions.o \
	obj-mkl/EigenStratImputation.o 	\
	obj-mkl/EigenSym.o		\
	obj-mkl/FastLMM.o		\
	obj-mkl/FastLmmC.o		\
	obj-mkl/FastLMM_fullrank.o 	\
	obj-mkl/FastLMM_lowrank.o 	\
	obj-mkl/FastLMM_woodbury.o 	\
	obj-mkl/FisherF.o		\
	obj-mkl/Gamma.o 		\
	obj-mkl/GroupSim.o      \
	obj-mkl/LinReg.o        \
	obj-mkl/LmmGwas.o	 	\
	obj-mkl/MathFunctions.o \
	obj-mkl/MatrixIO.o		\
	obj-mkl/Qval.o			\
	obj-mkl/Splash.o		\
	obj-mkl/Startup.o		\
	obj-mkl/StudyData.o     \
	obj-mkl/TsvLexer.o		\
	obj-mkl/Utils.o			\
	-Wl,--start-group		\
	/opt/intel/mkl/lib/intel64/libmkl_blas95_lp64.a		\
	/opt/intel/mkl/lib/intel64/libmkl_lapack95_lp64.a	\
	/opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a		\
	/opt/intel/mkl/lib/intel64/libmkl_intel_thread.a	\
	/opt/intel/mkl/lib/intel64/libmkl_core.a		\
	/opt/intel/lib/intel64/libiomp5.a	\
	-Wl,--end-group	\
	-lpthread	\
	-lz		\
	-ldl	\
	-lrt	\
	-static

echo "Copy to fastlmmc and zip the executable file"
cp fastlmmc_mkl fastlmmc
bzip2 -k -f fastlmmc 
