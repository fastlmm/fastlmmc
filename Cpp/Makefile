# FaST-LMM C++ (fastlmmc) Linux build against Intel oneMKL, statically linked.
#
#   make -j MKLROOT=~/.cache/fastlmmcpp/mkl-2026.2.1   # build ./fastlmmc
#   make test MKLROOT=...                               # build, then run the regression suite
#   make OPENMP=1 MKLROOT=...                           # experimental: ./fastlmmc-omp (see below)
#   make clean
#
# MKLROOT defaults to a standard oneAPI install. Without root, oneMKL can come from PyPI:
#   uv pip install --target ~/.cache/fastlmmcpp/mkl-2026.2.1 mkl-devel==2026.2.1 mkl-static==2026.2.1
#
# OPENMP=1 compiles the '#pragma omp parallel for' loops in FastLMM_{fullrank,lowrank,woodbury}.cpp
# (the 2.07c Windows release build had OpenMP on; the 2.07c Linux build did not) and links MKL's
# GNU-threading layer so only one OpenMP runtime (libgomp) is linked.

MKLROOT     ?= /opt/intel/oneapi/mkl/latest
MKL_INCLUDE ?= $(MKLROOT)/include
MKL_LIB     ?= $(MKLROOT)/lib
# Intel OpenMP runtime: beside MKL in the PyPI layout; with the compiler in a oneAPI install.
IOMP5_LIB   ?= $(firstword $(wildcard $(MKL_LIB)/libiomp5.a) $(MKLROOT)/../../compiler/latest/lib/libiomp5.a)

OPT    ?= -O3
OPENMP ?= 0

ifeq ($(OPENMP),1)
  BIN     ?= fastlmmc-omp
  OBJ_DIR := obj-mkl-omp
  OMPFLAG := -fopenmp
  MKL_THREAD_LIBS := $(MKL_LIB)/libmkl_gnu_thread.a
else
  BIN     ?= fastlmmc
  OBJ_DIR := obj-mkl
  OMPFLAG :=
  MKL_THREAD_LIBS := $(MKL_LIB)/libmkl_intel_thread.a $(IOMP5_LIB)
endif

# $(OBJ_DIR) comes first so the generated BuildVer.h wins over any stale Inc/BuildVer.h.
CPPFLAGS := -DUSE_MKL -I$(OBJ_DIR) -IInc -I$(MKL_INCLUDE)
CXXFLAGS := $(OPT) -std=c++11 -Wall $(OMPFLAG) -MMD -MP
LDFLAGS  := -static $(OMPFLAG)
LDLIBS   := -Wl,--start-group $(MKL_LIB)/libmkl_intel_lp64.a $(MKL_THREAD_LIBS) $(MKL_LIB)/libmkl_core.a -Wl,--end-group \
            -lpthread -lz -ldl -lrt -lm

SRCS := $(wildcard Src/*.cpp)
OBJS := $(patsubst Src/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

.PHONY: all test clean FORCE
all: $(BIN)

ifeq ($(filter clean,$(MAKECMDGOALS)),)
  ifeq ($(wildcard $(MKL_INCLUDE)/mkl.h),)
    $(error oneMKL not found ($(MKL_INCLUDE)/mkl.h). Set MKLROOT; see the top of this Makefile)
  endif
endif

$(BIN): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(OBJ_DIR)/%.o: Src/%.cpp | $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Build date and host shown in the splash text. Rewritten only when its content changes,
# so Splash.cpp recompiles at most once a day.
$(OBJ_DIR)/BuildVer.h: FORCE | $(OBJ_DIR)
	@printf '// Generated file **DO NOT EDIT**\nint wVerBuild = %s;\n#define szBuildHost "%s"\n' \
	    "$$(date +%Y%m%d)" "$$(uname -n)" > $@.tmp
	@if cmp -s $@.tmp $@; then rm $@.tmp; else mv $@.tmp $@; fi
$(OBJ_DIR)/Splash.o: $(OBJ_DIR)/BuildVer.h

$(OBJ_DIR):
	mkdir -p $@

test: $(BIN)
	../tests/regression/run.sh ./$(BIN)

clean:
	rm -rf obj-mkl obj-mkl-omp fastlmmc fastlmmc-omp

-include $(OBJS:.o=.d)
