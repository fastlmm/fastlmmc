#!/bin/bash

if [ "$1" = "v2.x" ]; then
  LINUX_KERNEL=-DUSE_LINUX_V2x_KERNEL
  echo "Compiling for RedHat Enterprise Linux with v2.x kernel (LINUX_KERNEL=$LINUX_KERNEL)"
else
  LINUX_KERNEL=
  echo "Compiling for Linux with V3.x kernel (LINUX_KERNEL=$LINUX_KERNEL)"
fi

# oneMKL location. Defaults to a standard oneAPI install; override MKLROOT to use
# another copy, e.g. one installed from PyPI with:
#   uv pip install --target ~/.cache/fastlmmcpp/mkl-2026.2.1 mkl-devel==2026.2.1 mkl-static==2026.2.1
#   MKLROOT=~/.cache/fastlmmcpp/mkl-2026.2.1 ./DoMkl_linux.sh
MKLROOT=${MKLROOT:-/opt/intel/oneapi/mkl/latest}
MKL_INCLUDE=${MKL_INCLUDE:-$MKLROOT/include}
MKL_LIB=${MKL_LIB:-$MKLROOT/lib}
# Intel OpenMP runtime: in oneAPI it ships with the compiler; the PyPI layout puts it beside MKL.
if [ -z "$IOMP5_LIB" ]; then
  if [ -f "$MKL_LIB/libiomp5.a" ]; then
    IOMP5_LIB=$MKL_LIB/libiomp5.a
  else
    IOMP5_LIB=$MKLROOT/../../compiler/latest/lib/libiomp5.a
  fi
fi
for f in "$MKL_INCLUDE/mkl.h" "$MKL_LIB/libmkl_core.a" "$IOMP5_LIB"; do
  if [ ! -f "$f" ]; then echo "Missing $f (set MKLROOT or IOMP5_LIB)"; exit 1; fi
done
echo "Using MKLROOT=$MKLROOT"
set -e
mkdir -p obj-mkl

echo "Generating Inc/BuildVer.h"
echo "// Generated File  ** DO NOT EDIT **" >Inc/BuildVer.h
echo "int wVerBuild = `date +%Y%m%d`;"     >>Inc/BuildVer.h
echo "#define szBuildHost \"`uname -n`\""  >>Inc/BuildVer.h

echo "Compiling sources with MKL Library"
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/AutoSelectSnps.cpp       -o ./obj-mkl/AutoSelectSnps.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/Beta.cpp                 -o ./obj-mkl/Beta.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/BrentC.cpp               -o ./obj-mkl/BrentC.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CCovariatesData.cpp      -o ./obj-mkl/CCovariatesData.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CKernelData.cpp          -o ./obj-mkl/CKernelData.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CLexer.cpp               -o ./obj-mkl/CLexer.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CLogisticRegression.cpp  -o ./obj-mkl/CLogisticRegression.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CPhenData.cpp            -o ./obj-mkl/CPhenData.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CPlinkAlternatePhenotypeFile.cpp -o ./obj-mkl/CPlinkAlternatePhenotypeFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CPlinkBedFile.cpp        -o ./obj-mkl/CPlinkBedFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CPlinkBigFile.cpp        -o ./obj-mkl/CPlinkBigFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CPlinkBimFile.cpp        -o ./obj-mkl/CPlinkBimFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CPlinkDatFile.cpp        -o ./obj-mkl/CPlinkDatFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CPlinkFamFile.cpp        -o ./obj-mkl/CPlinkFamFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CPlinkFile.cpp           -o ./obj-mkl/CPlinkFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CPlinkLexer.cpp          -o ./obj-mkl/CPlinkLexer.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CPlinkMapFile.cpp        -o ./obj-mkl/CPlinkMapFile.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CSnpData.cpp             -o ./obj-mkl/CSnpData.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/CSnpFilterOptions.cpp    -o ./obj-mkl/CSnpFilterOptions.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/EigenStratImputation.cpp -o ./obj-mkl/EigenStratImputation.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/EigenSym.cpp             -o ./obj-mkl/EigenSym.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/FastLMM.cpp              -o ./obj-mkl/FastLMM.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/FastLMM_fullrank.cpp     -o ./obj-mkl/FastLMM_fullrank.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/FastLMM_lowrank.cpp      -o ./obj-mkl/FastLMM_lowrank.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/FastLMM_woodbury.cpp     -o ./obj-mkl/FastLMM_woodbury.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/FastLmmC.cpp             -o ./obj-mkl/FastLmmC.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/FisherF.cpp              -o ./obj-mkl/FisherF.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/Gamma.cpp                -o ./obj-mkl/Gamma.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/GroupSim.cpp             -o ./obj-mkl/GroupSim.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/LinReg.cpp               -o ./obj-mkl/LinReg.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/LmmGwas.cpp              -o ./obj-mkl/LmmGwas.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/MathFunctions.cpp        -o ./obj-mkl/MathFunctions.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/MatrixIO.cpp             -o ./obj-mkl/MatrixIO.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/Qval.cpp                 -o ./obj-mkl/Qval.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/Splash.cpp               -o ./obj-mkl/Splash.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/Startup.cpp              -o ./obj-mkl/Startup.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/StudyData.cpp            -o ./obj-mkl/StudyData.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/TsvLexer.cpp             -o ./obj-mkl/TsvLexer.o
g++ -c -O3 -std=c++0x -DUSE_MKL $LINUX_KERNEL -I./Inc -I"$MKL_INCLUDE" Src/Utils.cpp                -o ./obj-mkl/Utils.o

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
	"$MKL_LIB/libmkl_blas95_lp64.a"		\
	"$MKL_LIB/libmkl_lapack95_lp64.a"	\
	"$MKL_LIB/libmkl_intel_lp64.a"		\
	"$MKL_LIB/libmkl_intel_thread.a"	\
	"$MKL_LIB/libmkl_core.a"		\
	"$IOMP5_LIB"	\
	-Wl,--end-group	\
	-lpthread	\
	-lz		\
	-ldl	\
	-lrt	\
	-static

echo "Copy to fastlmmc and zip the executable file"
cp fastlmmc_mkl fastlmmc
bzip2 -k -f fastlmmc 
