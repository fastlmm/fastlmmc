#!/bin/bash

echo "// Generated File  ** DO NOT EDIT **" >Inc/BuildVer.h
echo "int wVerBuild = `date +%Y%m%d`;"     >>Inc/BuildVer.h
echo "#define szBuildHost \"`uname -n`\""  >>Inc/BuildVer.h

echo "Compiling sources with ACML Library"
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/AutoSelectSnps.cpp        -o ./obj-acml/AutoSelectSnps.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/Beta.cpp                 -o ./obj-acml/Beta.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/BrentC.cpp               -o ./obj-acml/BrentC.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CCovariatesData.cpp      -o ./obj-acml/CCovariatesData.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CKernelData.cpp          -o ./obj-acml/CKernelData.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CLexer.cpp               -o ./obj-acml/CLexer.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CLogisticRegression.cpp  -o ./obj-acml/CLogisticRegression.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CPhenData.cpp            -o ./obj-acml/CPhenData.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CPlinkAlternatePhenotypeFile.cpp -o ./obj-acml/CPlinkAlternatePhenotypeFile.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CPlinkBedFile.cpp        -o ./obj-acml/CPlinkBedFile.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CPlinkBigFile.cpp        -o ./obj-acml/CPlinkBigFile.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CPlinkBimFile.cpp        -o ./obj-acml/CPlinkBimFile.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CPlinkDatFile.cpp        -o ./obj-acml/CPlinkDatFile.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CPlinkFamFile.cpp        -o ./obj-acml/CPlinkFamFile.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CPlinkFile.cpp           -o ./obj-acml/CPlinkFile.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CPlinkLexer.cpp          -o ./obj-acml/CPlinkLexer.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CPlinkMapFile.cpp        -o ./obj-acml/CPlinkMapFile.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CSnpData.cpp             -o ./obj-acml/CSnpData.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/CSnpFilterOptions.cpp    -o ./obj-acml/CSnpFilterOptions.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/EigenStratImputation.cpp -o ./obj-acml/EigenStratImputation.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/EigenSym.cpp             -o ./obj-acml/EigenSym.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/FastLMM.cpp              -o ./obj-acml/FastLMM.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/FastLMM_fullrank.cpp     -o ./obj-acml/FastLMM_fullrank.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/FastLMM_lowrank.cpp      -o ./obj-acml/FastLMM_lowrank.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/FastLMM_woodbury.cpp     -o ./obj-acml/FastLMM_woodbury.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/FastLmmC.cpp             -o ./obj-acml/FastLmmC.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/FisherF.cpp              -o ./obj-acml/FisherF.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/Gamma.cpp                -o ./obj-acml/Gamma.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/GroupSim.cpp             -o ./obj-acml/GroupSim.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/LinReg.cpp               -o ./obj-acml/LinReg.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/LmmGwas.cpp              -o ./obj-acml/LmmGwas.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/MathFunctions.cpp        -o ./obj-acml/MathFunctions.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/MatrixIO.cpp             -o ./obj-acml/MatrixIO.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/Qval.cpp                 -o ./obj-acml/Qval.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/Splash.cpp               -o ./obj-acml/Splash.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/Startup.cpp              -o ./obj-acml/Startup.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/StudyData.cpp            -o ./obj-acml/StudyData.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/TsvLexer.cpp             -o ./obj-acml/TsvLexer.o
g++ -c -O3 -DUSE_ACML -I./Inc -I/opt/acml4.4.0/gfortran64/include Src/Utils.cpp                -o ./obj-acml/Utils.o

echo "Linking fastlmmc_acml"
g++ -o ./fastlmmc_acml		\
	obj-acml/AutoSelectSnps.o	\
	obj-acml/Beta.o         \
	obj-acml/BrentC.o		\
	obj-acml/CCovariatesData.o	\
	obj-acml/CKernelData.o	\
	obj-acml/CLexer.o		\
	obj-acml/CLogisticRegression.o   \
	obj-acml/CPhenData.o	\
	obj-acml/CPlinkAlternatePhenotypeFile.o	\
	obj-acml/CPlinkBedFile.o	\
	obj-acml/CPlinkBigFile.o	\
	obj-acml/CPlinkBimFile.o	\
	obj-acml/CPlinkDatFile.o	\
	obj-acml/CPlinkFamFile.o	\
	obj-acml/CPlinkFile.o	\
	obj-acml/CPlinkLexer.o	\
	obj-acml/CPlinkMapFile.o	\
	obj-acml/CSnpData.o		\
	obj-acml/CSnpFilterOptions.o \
	obj-acml/EigenStratImputation.o	\
	obj-acml/EigenSym.o		\
	obj-acml/FastLMM.o		\
	obj-acml/FastLmmC.o		\
	obj-acml/FastLMM_fullrank.o	\
	obj-acml/FastLMM_lowrank.o	\
	obj-acml/FastLMM_woodbury.o	\
	obj-acml/FisherF.o		\
	obj-acml/Gamma.o		\
	obj-acml/GroupSim.o		\
	obj-acml/LinReg.o		\
	obj-acml/LmmGwas.o		\
	obj-acml/MathFunctions.o 	\
	obj-acml/MatrixIO.o		\
	obj-acml/Qval.o			\
	obj-acml/Splash.o		\
	obj-acml/Startup.o		\
	obj-acml/StudyData.o	\
	obj-acml/TsvLexer.o		\
	obj-acml/Utils.o		\
	/opt/acml4.4.0/gfortran64/lib/libacml.a	\
	/usr/lib/gcc/x86_64-linux-gnu/4.4/libgfortran.a	\
	-lz						\
	-static -lrt

echo "Copy to fastlmmc and zip the executable file"
cp fastlmmc_acml fastlmmc
bzip2 -k -f fastlmmc 
