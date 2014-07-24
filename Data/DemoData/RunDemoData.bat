@Echo off
if DEFINED ECHO Echo %ECHO%
SetLocal
REM This file runs FaST-LMM through the commands in the user manual
REM Please call runsampledata.bat from the DemoData directory.
REM
REM Add fastlmmc.exe to the path
REM
if     "%1" == ""  Set PATH=..\..\Bin\Windows\Cpp_MKL;%PATH%
if NOT "%1" == ""  Set PATH=%1;%PATH%
@Echo on
fastlmmc -verboseOutput -bfile snps -bfilesim snps -pheno pheno.txt -out snps.dummy.txt
@if errorlevel 1 goto :Fatal

fastlmmc -autoselect ASout -autoSelectSearchValues "0 1 2 4 8 16 32 64 128 256 512 1024 2048 10100" -randomSeed 1 -autoselectFolds 10 -bfilesim snps -pheno pheno.txt -covar covariate.txt
@if errorlevel 1 goto :Fatal

fastlmmc -bfile snps -bfilesim snps -extractSim ASout.snps.txt -pheno pheno.txt -covar covariate.txt -out snps.insample.txt -excludebyPosition 1
@if errorlevel 1 goto :Fatal


@Rem This command covers every pair for the 10100 SNPs in snps.bim. (50,999,950 pair)
@Rem   so it takes a longer than a 'normal' demo run would take.
@Rem
@if NOT "%2%" == "" fastlmmc -snpPairs -bfile snps -bfileSim snps -extractSim ASout.snps.txt -pheno pheno.txt -covar covariate.txt -out snps.pairs.txt
@if errorlevel 1 goto :Fatal

@Rem
@Rem You can focus the analysis on a smaller set of snps by extracting them from the test set
@Rem  using the -extract option.  
fastlmmc -snpPairs -bfile snps -extract topsnps.txt -bfileSim snps -extractSim ASout.snps.txt -pheno pheno.txt -covar covariate.txt -out topsnps.pairs.txt
@if errorlevel 1 goto :Fatal

fastlmmc -snpPairs -bfile snps -bfileSim snps -extractSim ASout.snps.txt -pheno pheno.txt -covar covariate.txt -out snps.pairs.250.000.txt -tasks 250 0
@if errorlevel 1 goto :Fatal

goto :EOF

:Fatal
Echo.
Echo fastlmmc returned an error.  Stopping script.
Echo.

