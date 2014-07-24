@Echo Off
@if DEFINED ECHO  Echo %ECHO%
@SetLocal
@if NOT DEFINED PscpLinux  Set PscpLinux=Erg00Lx
@if NOT DEFINED PscpUser   Set PscpUser=%USERNAME%
@if NOT DEFINED PscpPswd   Set PscpPswd=%1
@if NOT DEFINED PscpPswd   goto :Usage

Set SrcDir=%CD%
Set DstDir=%PscpLinux%:Projects/fastlmm

Echo Copy the files required to run tests on Linux
Echo    from Windows: %SrcDir%
Echo        to Linux: %DstDir%
Echo.
Echo Copy the Test Script file
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\DoTest_linux.sh %DstDir%/DoTest_linux.sh

Set SrcDir=%CD%\..\Data\FastLmmC\PLink.100
Set DstDir=%PscpLinux%:Projects/fastlmm/Data/FastLmmC/PLink.100

Echo.
Echo Copy all the PLink.100.* files
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\AutoSelectSearchValues.txt %DstDir%/AutoSelectSearchValues.txt
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\PLink.100.bed  %DstDir%/PLink.100.bed
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\PLink.100.big  %DstDir%/PLink.100.big
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\PLink.100.bim  %DstDir%/PLink.100.bim
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\PLink.100.covariates.txt   %DstDir%/PLink.100.covariates.txt
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\PLink.100.dat  %DstDir%/PLink.100.dat
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\PLink.100.fam  %DstDir%/PLink.100.fam
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\PLink.100.map  %DstDir%/PLink.100.map
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\PLink.100.ped  %DstDir%/PLink.100.ped
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\PLink.100.pheno    %DstDir%/PLink.100.pheno
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\plink.100.tfam %DstDir%/plink.100.tfam
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\plink.100.tped %DstDir%/plink.100.tped
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\SnpSelection.00.04.txt %DstDir%/SnpSelection.00.04.txt
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\SnpSelection.00.99.txt %DstDir%/SnpSelection.00.99.txt
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\SnpSelection.55.59.txt %DstDir%/SnpSelection.55.59.txt
pscp -pw %PscpPswd% -l %PscpUser% %SrcDir%\SnpSelection.90.99.txt %DstDir%/SnpSelection.90.99.txt
goto :eof

:Usage
:HelpMsg
Echo.
Echo DoPscpData copies the files necessary to run FastLmmC tests on
Echo   a Linux system.  The files will all be copied to 
Echo   Projects/fastlmm/Data/FastLmmC/PLink.100 and are used by
Echo   the test script on the Linux system
Echo.
Echo   Environment variables used
Echo     PscpLinux - Linux machine to communicate with [%PscpLinux%]
Echo     PscpUser  - Linux username [%PscpUser%]
Echo     PscpPswd  - your linux password []
Echo.
Echo   Example:
Echo      CopyTestFiles2Linux [password]
Echo.
Echo    After the copy is complete you can run the FastLmmC 
Echo      test on the Linux box using
Echo      ./DoTest_linux
Echo.
