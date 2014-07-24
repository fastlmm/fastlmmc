@Echo Off
@if DEFINED ECHO  Echo %ECHO%
@SetLocal
@if NOT DEFINED PscpPswd   Set PscpPswd=%1
@if NOT DEFINED PscpPswd   goto :Usage
@if NOT DEFINED PscpLinux  Set PscpLinux=Erg00Lx
@if NOT DEFINED PscpUser   Set PscpUser=%USERNAME%

Set SrcDir=%CD%
Set DstDir=%PscpLinux%:Projects/fastlmm

@Echo.
@Echo Copy License and Script files
@Echo    from Windows: %SrcDir%
@Echo        to Linux: %DstDir%
@Echo.
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\..\LICENSE.TXT  %DstDir%/LICENSE.TXT
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\..\README.TXT   %DstDir%/README.TXT
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\..\NOTICE.TXT   %DstDir%/NOTICE.TXT
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\DoAcml_linux.sh %DstDir%/DoAcml_linux.sh
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\DoMkl_linux.sh  %DstDir%/DoMkl_linux.sh
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\DoTest_linux.sh %DstDir%/DoTest_linux.sh

Echo.
Echo Copy Sources
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\AutoSelectSnps.cpp      %DstDir%/Src/AutoSelectSnps.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\Beta.cpp                %DstDir%/Src/Beta.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\BrentC.cpp              %DstDir%/Src/BrentC.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CCovariatesData.cpp     %DstDir%/Src/CCovariatesData.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CKernelData.cpp         %DstDir%/Src/CKernelData.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CLexer.cpp              %DstDir%/Src/CLexer.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CLogisticRegression.cpp %DstDir%/Src/CLogisticRegression.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CPhenData.cpp           %DstDir%/Src/CPhenData.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CPlinkAlternatePhenotypeFile.cpp %DstDir%/Src/CPlinkAlternatePhenotypeFile.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CPlinkBedFile.cpp                %DstDir%/Src/CPlinkBedFile.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CPlinkBigFile.cpp                %DstDir%/Src/CPlinkBigFile.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CPlinkBimFile.cpp                %DstDir%/Src/CPlinkBimFile.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CPlinkDatFile.cpp                %DstDir%/Src/CPlinkDatFile.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CPlinkFamFile.cpp                %DstDir%/Src/CPlinkFamFile.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CPlinkFile.cpp                   %DstDir%/Src/CPlinkFile.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CPlinkLexer.cpp                  %DstDir%/Src/CPlinkLexer.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CPlinkMapFile.cpp                %DstDir%/Src/CPlinkMapFile.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CSnpData.cpp                     %DstDir%/Src/CSnpData.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\CSnpFilterOptions.cpp            %DstDir%/Src/CSnpFilterOptions.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\EigenStratImputation.cpp         %DstDir%/Src/EigenStratImputation.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\EigenSym.cpp         %DstDir%/Src/EigenSym.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\FastLMM.cpp          %DstDir%/Src/FastLMM.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\FastLmmC.cpp         %DstDir%/Src/FastLmmC.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\FastLMM_fullrank.cpp %DstDir%/Src/FastLMM_fullrank.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\FastLMM_lowrank.cpp  %DstDir%/Src/FastLMM_lowrank.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\FastLMM_woodbury.cpp %DstDir%/Src/FastLMM_woodbury.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\FisherF.cpp          %DstDir%/Src/FisherF.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\Gamma.cpp            %DstDir%/Src/Gamma.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\GroupSim.cpp         %DstDir%/Src/GroupSim.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\LinReg.cpp           %DstDir%/Src/LinReg.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\LmmGwas.cpp          %DstDir%/Src/LmmGwas.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\MathFunctions.cpp    %DstDir%/Src/MathFunctions.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\MatrixIO.cpp         %DstDir%/Src/MatrixIO.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\Qval.cpp             %DstDir%/Src/Qval.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\Splash.cpp           %DstDir%/Src/Splash.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\Startup.cpp          %DstDir%/Src/Startup.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\StudyData.cpp        %DstDir%/Src/StudyData.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\TsvLexer.cpp         %DstDir%/Src/TsvLexer.cpp
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Src\Utils.cpp            %DstDir%/Src/Utils.cpp

Echo.
Echo Copy Headers
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\AutoSelectSnps.h      %DstDir%/Inc/AutoSelectSnps.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\Beta.h                %DstDir%/Inc/Beta.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\BlasLapack.h          %DstDir%/Inc/BlasLapack.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\BrentC.h              %DstDir%/Inc/BrentC.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\BuildVer.h            %DstDir%/Inc/BuildVer.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CCovariatesData.h     %DstDir%/Inc/CCovariatesData.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CKernelData.h         %DstDir%/Inc/CKernelData.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CLexer.h              %DstDir%/Inc/CLexer.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CLogisticRegression.h %DstDir%/Inc/CLogisticRegression.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CPhenData.h           %DstDir%/Inc/CPhenData.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CPLinkAlternatePhenotypeFile.h %DstDir%/Inc/CPlinkAlternatePhenotypeFile.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CPlinkBedFile.h      %DstDir%/Inc/CPlinkBedFile.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CPlinkBigFile.h      %DstDir%/Inc/CPlinkBigFile.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CPlinkBimFile.h      %DstDir%/Inc/CPlinkBimFile.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CPlinkDatFile.h      %DstDir%/Inc/CPlinkDatFile.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CPlinkFamFile.h      %DstDir%/Inc/CPlinkFamFile.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CPlinkFile.h         %DstDir%/Inc/CPlinkFile.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CPlinkLexer.h        %DstDir%/Inc/CPlinkLexer.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CPlinkMapFile.h      %DstDir%/Inc/CPlinkMapFile.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CSnpData.h           %DstDir%/Inc/CSnpData.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CSnpFilterOptions.h  %DstDir%/Inc/CSnpFilterOptions.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\CSnpInfo.h           %DstDir%/Inc/CSnpInfo.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\EigenSym.h           %DstDir%/Inc/EigenSym.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\FastLMM.h            %DstDir%/Inc/FastLMM.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\FastLmmC.h           %DstDir%/Inc/FastLmmC.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\FastLMM_fullrank.h   %DstDir%/Inc/FastLMM_fullrank.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\FastLMM_lowrank.h    %DstDir%/Inc/FastLMM_lowrank.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\FastLMM_woodbury.h   %DstDir%/Inc/FastLMM_woodbury.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\FisherF.h            %DstDir%/Inc/FisherF.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\Gamma.h              %DstDir%/Inc/Gamma.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\GroupSim.h           %DstDir%/Inc/GroupSim.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\LinReg.h             %DstDir%/Inc/LinReg.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\LmmGwas.h            %DstDir%/Inc/LmmGwas.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\MathFunctions.h      %DstDir%/Inc/MathFunctions.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\MatrixIO.h           %DstDir%/Inc/MatrixIO.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\minimize.h           %DstDir%/Inc/minimize.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\Qval.h               %DstDir%/Inc/Qval.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\Startup.h            %DstDir%/Inc/Startup.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\StudyData.h          %DstDir%/Inc/StudyData.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\TsvLexer.h           %DstDir%/Inc/TsvLexer.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\Utils.h              %DstDir%/Inc/Utils.h
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\Inc\vec.h                %DstDir%/Inc/vec.h

goto :EOF

:Usage
:HelpMsg
Echo.
Echo CopySrcFiles2Linux copies the files necessary to build FastLmmC 
Echo    to a Linux system.  The files will be copied to the users 
Echo    ~/Project/fastlmm directory and it is assumed that the
Echo    default directory tree structure is already in place.  
Echo    This includes giving the directories and files the
Echo    appropriate atributes and ownership flags too.  
Echo    The sub-directories expected in the ~Project/fastlmm 
Echo    directory include:
Echo      Inc
Echo      Src
Echo      obj-acml
Echo      obj-mkl
Echo      Test
Echo      Data/FastLmmC/PLink.100
Echo    The following environment variable may also be used to
Echo    modify the default behavior:
Echo      PscpLinux - Linux system name [Erg00Lx]
Echo      PscpUser  - User name on Linux system [%USERNAME%]
Echo      PscpPswd  - password associated with PscpUser []
Echo.
Echo   Example:
Echo     CopySrcFiles2Linux pswd
Echo.
Echo   After the copy completes, you can build on the Linux box using
Echo     ./DoMkl_linux
Echo   and
Echo     ./DoAcml_linux
Echo.
goto :EOF
