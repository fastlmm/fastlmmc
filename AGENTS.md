# fastlmmcpp

The FaST-LMM C++ program (`fastlmmc`): history imported from 2011 and 2016 snapshots and the
2014 FaSTLMM 2.07c release, then updated in 2026 to build with current compilers and oneMKL.
The code is under the Microsoft Research License Agreement (`LICENSE.TXT`), not an open-source
license.

## Build

Linux (g++, oneMKL, static binary):

```bash
make -C Cpp -j MKLROOT=~/.cache/fastlmmcpp/mkl-2026.2.1
```

Without root, oneMKL comes from PyPI:
`uv pip install --target ~/.cache/fastlmmcpp/mkl-2026.2.1 mkl-devel==2026.2.1 mkl-static==2026.2.1`.
`make OPENMP=1` builds an experimental `fastlmmc-omp` with the source's OpenMP loops turned on.

Windows: `Cpp/FastLmmC.sln` has two configurations, `MKL-Release|x64` and `MKL-Debug|x64`
(Visual Studio 2022, toolset v143). Set `MKLROOT` (a folder with `include\` and `lib\`; for the
PyPI `win_amd64` packages, their `Library\` folder) and `ZLIB_ROOT` (`include\`, plus a static
zlib as `lib\zlibstat.lib` built with `/MT` and `lib\Debug\zlibstat.lib` built with `/MTd`).
The 32-bit and AMD ACML configurations were removed in 2.08; the `USE_ACML` code paths remain.

## Test

```bash
tests/regression/run.sh Cpp/fastlmmc
```

See `tests/regression/README.md`. Run it after any change that could affect numerical results.

## Rules

- The version is `wVerMajor`/`wVerMinor` in `Cpp/Src/Splash.cpp` plus the VERSIONINFO block in
  `Cpp/FastLmmC.rc` (UTF-16); change both together. 2.07 is the 2014 release; 2.08 is the 2026
  update. Write the minor number without a leading zero (`08` is invalid octal in C++).

- Sources use CRLF line endings, and `FastLmmC.rc` and `resource.h` are UTF-16; preserve both.
- `Inc/BlasLapack.h` refuses to run with an MKL version other than the validated one unless
  `FastLmmUseAnyMklLib=1` is set. When moving to a new MKL, run the regression suite, then
  update the pinned version there.
