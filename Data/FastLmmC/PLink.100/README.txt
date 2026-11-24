PLink.100 test inputs for Cpp/DoTest_linux.sh

Copied from the 2012 Escience/RegressionTest/inputs folder (carlk3BioVirusCount snapshot):
PLink.100.{bed,bim,fam,ped,map,dat,big,pheno,1.pheno}, plink.100.{tped,tfam},
SnpSelection.{00.04,00.99,50.54,90.99}.txt

Reconstructed (the originals were not found):
  SnpSelection.55.59.txt - SnpSelection.50.54.txt with the window moved to entries 55-59.

Stand-ins (synthetic; the originals were not found, so results that use them
cannot be compared with any historical output, only between builds):
  AutoSelectSearchValues.txt - FastLmmC's built-in default search values, capped at 100 SNPs.
  PLink.100.covariates.txt   - two covariates per individual (Gaussian, 0/1), Python random.Random(20261124).
