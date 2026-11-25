#!/bin/bash
# Regression cases added in 2026 beyond Cpp/DoTest_linux.sh: they exercise the low-rank,
# Woodbury and forced full-rank paths, REML, SNP pairs, and the .tped, .dat (dosage2) and .big readers.
# Same conventions as DoTest_linux.sh: run from the work dir that holds ./fastlmmc and Data/.
D=../Data/FastLmmC/PLink.100
P="-pheno $D/PLink.100.pheno -mpheno 1"
pushd Test
# 5-SNP kernel on 99 individuals: rank-deficient, so the low-rank path (ML default) ...
../fastlmmc -vo -sort- -file $D/PLink.100 $P -filesim $D/PLink.100 -extractSim $D/SnpSelection.00.04.txt              -out ./Out/x.lowrank.out.txt
# ... the Woodbury path ...
../fastlmmc -vo -sort- -file $D/PLink.100 $P -filesim $D/PLink.100 -extractSim $D/SnpSelection.00.04.txt -woodbury    -out ./Out/x.woodbury.out.txt
# ... and the same kernel forced through the full-rank path.
../fastlmmc -vo -sort- -file $D/PLink.100 $P -filesim $D/PLink.100 -extractSim $D/SnpSelection.00.04.txt -fullrank    -out ./Out/x.fullrank.out.txt
../fastlmmc -vo -sort- -file $D/PLink.100 $P -filesim $D/PLink.100 -REML                                               -out ./Out/x.reml.out.txt
../fastlmmc -vo -sort- -file $D/PLink.100 $P -filesim $D/PLink.100 -extractSim $D/SnpSelection.00.04.txt -snpPairs    -out ./Out/x.snppairs.out.txt
../fastlmmc -vo -sort- -tfile $D/plink.100 $P -tfilesim $D/plink.100                                                   -out ./Out/x.tfile.out.txt
../fastlmmc -vo -sort- -bfile $D/PLink.100 $P -dosage2sim $D/PLink.100                                                 -out ./Out/x.dosage2sim.out.txt
../fastlmmc -vo -sort- -afile $D/PLink.100 $P -bfilesim $D/PLink.100                                                   -out ./Out/x.afile.out.txt
popd
