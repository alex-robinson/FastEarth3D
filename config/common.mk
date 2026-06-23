# Shared dependency wiring for FastEarth3D.
#
# Loaded by config/Makefile *after* the compiler fragment, so it may reference
# FFLAGS_BASE, MODFLAGS, CPPFLAGS_PP, INC_NC and LIB_NC defined there.
#
# All numerical dependencies (FFTW, SHTns) and the fesm-utils helper library are
# provided by the fesm-utils package, expected at the repo root as a symlink:
#
#     ln -s ../fesm-utils fesm-utils
#
# SHTns (spherical-harmonic transforms) is built into fesm-utils with:
#     cd fesm-utils && ./build.py -m <machine> -c <compiler> --component shtns
# It links FFTW, which fesm-utils also builds. Serial variants are used by
# default; `make openmp=1` swaps in the OpenMP variants below.

# --- fesm-utils helper library (ncio, nml, mapping_scrip, ...) ---------------
FESMUTILSROOT = fesm-utils/utils
INC_FESMUTILS = -I$(FESMUTILSROOT)/include-serial
LIB_FESMUTILS = -L$(FESMUTILSROOT)/include-serial -lfesmutils

# --- FFTW --------------------------------------------------------------------
FFTWROOT = fesm-utils/fftw-serial
INC_FFTW = -I$(FFTWROOT)/include
LIB_FFTW = -L$(FFTWROOT)/lib -lfftw3 -lm

# --- SHTns (provides shtns.f03 Fortran 2003 interface + libshtns.a) ----------
SHTNSROOT = fesm-utils/shtns-serial
INC_SHTNS = -I$(SHTNSROOT)/include
LIB_SHTNS = -L$(SHTNSROOT)/lib -lshtns

# --- LIS (Library of Iterative Solvers; provides lisf.h + liblis.a) ----------
# Backs the per-degree banded saddle-point solve in fe_radial_fe. The Fortran
# interface is the preprocessor header lisf.h (needs -cpp, set as CPPFLAGS_PP);
# INC_LIS is what lets `#include "lisf.h"` resolve. Real-scalar build.
LISROOT = fesm-utils/lis-serial
INC_LIS = -I$(LISROOT)/include
LIB_LIS = -L$(LISROOT)/lib -llis

# --- OpenMP ------------------------------------------------------------------
# `make openmp=1` adds -fopenmp (via config/Makefile, appended to FFLAGS) to thread
# the per-degree loop in fe_response (begin_step / commit_step). The dependency
# libraries deliberately stay SERIAL: we thread over independent per-degree systems
# ourselves, each solved by the re-entrant banded LU (fe_band); LIS is used only for
# the single j=1 system (run by one thread) and SHTns/FFTW are called outside the
# parallel regions. The omp library variants are intentionally NOT linked — lis-omp
# would nest threads inside our parallel loop, and concurrent LIS solves are not
# re-entrant anyway (which is why j>=2 moved to the banded LU in the first place).

# --- Final flag sets ---------------------------------------------------------
# MODFLAGS (-I/-J objdir) and FFLAGS_BASE come from the compiler fragment.
# INC_SHTNS is what lets `include 'shtns.f03'` in src/fe_sht.f90 be found.
CPPFLAGS_FE = $(CPPFLAGS_PP)
FFLAGS_FE   = $(FFLAGS_BASE) $(MODFLAGS) $(INC_NC) $(INC_FESMUTILS) $(INC_FFTW) $(INC_SHTNS) $(INC_LIS)

# Static archives resolve left-to-right, so a library must precede the libraries
# it depends on: SHTns before FFTW (SHTns calls FFTW), fesm-utils before netCDF.
# LIS is self-contained (only needs libm, already pulled in by FFTW/system).
LFLAGS_FE   = $(LIB_FESMUTILS) $(LIB_SHTNS) $(LIB_FFTW) $(LIB_LIS) $(LIB_NC) $(LFLAGS_EXTRA)
