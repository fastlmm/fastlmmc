# Regression tests for fastlmmc

`run.sh` runs a `fastlmmc` binary on three groups of cases and compares every output with
`expected/`:

- `Cpp/DoTest_linux.sh`: the 26 original FaSTLMM test commands, on `Data/FastLmmC/PLink.100`.
- `extra_tests.sh`: cases added in 2026 for the low-rank, Woodbury and forced full-rank
  paths, REML, SNP pairs, and the `.tped`, `.dat` (dosage2) and `.big` readers.
- The commands from `Data/DemoData/RunDemoData.bat` (3,000 individuals, 10,100 SNPs), except
  the all-pairs runs, which are too large for routine testing.

```bash
tests/regression/run.sh Cpp/fastlmmc          # or: make -C Cpp test MKLROOT=...
```

A full run takes a few minutes. Work files go under `~/.cache/fastlmmcpp/regression`
(override with `FASTLMM_TEST_WORK`). A Windows `FastLmmC.exe` can be tested from WSL the same
way; put the work folder on a Windows drive, for example
`FASTLMM_TEST_WORK=/mnt/c/deldir/fastlmmcpp/regression-work`.

For a Windows `.exe`, `run.sh` skips the AutoSelect cross-validation tables (`*.xval.txt`).
AutoSelect assigns folds with `std::random_shuffle`, which uses the C library's `rand()`, and
the Microsoft and GNU `rand()` sequences differ for the same `-randomSeed`. The SNPs that
AutoSelect finally chooses still match.

## Tolerances

`compare.py` accepts numbers that agree within a relative 1e-4 (absolute 1e-8 near zero), and
lets rows swap order when near-tied p-values sort differently. Eigenvectors are compared up to
sign, skipping repeated eigenvalues, whose eigenvectors are not unique. Different MKL versions
and CPU code paths typically agree to about 1e-14 in log-likelihoods and 1e-6 in p-values.

## Where `expected/` came from

`expected/` holds the outputs of the FaSTLMM 2.07c release binary (`Bin/Linux_MKL/fastlmmc`,
built in 2014 with MKL 11.0.4), regenerated with:

```bash
MKL_CBWR=COMPATIBLE tests/regression/run.sh <2.07c fastlmmc> --update
```

`MKL_CBWR=COMPATIBLE` is required: without it, that MKL crashes with "Illegal instruction"
on current CPUs (seen on an AMD Ryzen 9 7950X).

Some inputs are stand-ins; see `Data/FastLmmC/PLink.100/README.txt`. Results that use them
check consistency between builds, not agreement with any historical run.
