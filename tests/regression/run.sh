#!/bin/bash
# Regression suite for fastlmmc: runs Cpp/DoTest_linux.sh, extra_tests.sh and the demo commands
# from Data/DemoData/RunDemoData.bat, then compares the outputs with expected/ (see README.md).
#
# Usage: tests/regression/run.sh <fastlmmc binary> [--update]
#   --update  replace expected/ with this run's outputs instead of comparing.
# Work files go under $FASTLMM_TEST_WORK (default ~/.cache/fastlmmcpp/regression).
set -euo pipefail

if [ $# -lt 1 ] || [ $# -gt 2 ] || { [ $# -eq 2 ] && [ "$2" != "--update" ]; }; then
  echo "Usage: $0 <fastlmmc binary> [--update]" >&2
  exit 2
fi
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
bin=$(realpath "$1")
work=${FASTLMM_TEST_WORK:-$HOME/.cache/fastlmmcpp/regression}/$(date +%Y%m%d-%H%M%S)-$$

mkdir -p "$work/Test/Out" "$work/Demo"
cp "$bin" "$work/fastlmmc"
# A Windows FastLmmC.exe (run through WSL interop) needs the Intel OpenMP DLL beside it.
if [ -f "$(dirname "$bin")/libiomp5md.dll" ]; then cp "$(dirname "$bin")/libiomp5md.dll" "$work/"; fi
# A copy, not a symlink, so a Windows fastlmmc.exe run through WSL interop can read it too.
cp -r "$repo/Data" "$work/Data"
cd "$work"
echo "Binary:   $bin"
echo "Work dir: $work"

echo "Running Cpp/DoTest_linux.sh ..."
bash "$repo/Cpp/DoTest_linux.sh" > dotest.log 2>&1
echo "Running extra_tests.sh ..."
bash "$here/extra_tests.sh" > extra.log 2>&1

echo "Running demo commands ..."
cp "$repo"/Data/DemoData/* Demo/
(
  cd Demo
  ../fastlmmc -verboseOutput -bfile snps -bfilesim snps -pheno pheno.txt -out snps.dummy.txt
  ../fastlmmc -autoselect ASout -autoSelectSearchValues "0 1 2 4 8 16 32 64 128 256 512 1024 2048 10100" \
      -randomSeed 1 -autoselectFolds 10 -bfilesim snps -pheno pheno.txt -covar covariate.txt
  ../fastlmmc -bfile snps -bfilesim snps -extractSim ASout.snps.txt -pheno pheno.txt -covar covariate.txt \
      -out snps.insample.txt -excludebyPosition 1
  ../fastlmmc -snpPairs -bfile snps -extract topsnps.txt -bfileSim snps -extractSim ASout.snps.txt \
      -pheno pheno.txt -covar covariate.txt -out topsnps.pairs.txt
) > demo.log 2>&1 || echo "A demo command failed; see $work/demo.log"

mkdir -p out/demo
cp -r Test/Out out/dotest
cp Demo/snps.dummy.txt Demo/ASout.snps.txt Demo/ASout.xval.txt Demo/snps.insample.txt Demo/topsnps.pairs.txt out/demo/ \
  || echo "Some demo outputs are missing"

if grep -l 'Fatal Error' dotest.log extra.log demo.log; then
  echo "FAILED: fastlmmc reported a fatal error in the log(s) above (in $work)."
  exit 1
fi

if [ "${2:-}" = "--update" ]; then
  rm -rf "$here/expected"
  cp -r out "$here/expected"
  echo "Updated $here/expected from this run."
else
  ignore=()
  case "$bin" in
    *.exe|*.EXE)
      # AutoSelect's cross-validation folds come from std::random_shuffle, i.e. the C library's
      # rand(), whose sequence differs between the Microsoft and GNU runtimes. expected/ is
      # from Linux, so a Windows build's fold-by-fold results (*.xval.txt) cannot match.
      echo "Windows binary: skipping AutoSelect cross-validation tables (*.xval.txt)."
      ignore=(--ignore '*.xval.txt') ;;
  esac
  python3 -I "$here/compare.py" "$here/expected" out "${ignore[@]}"
fi
