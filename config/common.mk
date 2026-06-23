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

# --- (LIS removed) -----------------------------------------------------------
# The per-degree solve is now a dependency-free pivoted banded LU (fe_band); LIS
# is no longer linked. Keeping INC_LIS / LIB_LIS empty so the flag lists below
# (and any external references) stay valid.
INC_LIS =
LIB_LIS =

# --- OpenMP ------------------------------------------------------------------
# `make openmp=1` adds -fopenmp (via config/Makefile, appended to FFLAGS) to thread
# the per-degree loop in fe_response (begin_step / commit_step). The dependency
# libraries deliberately stay SERIAL: we thread over independent per-degree systems
# ourselves, each solved by the re-entrant banded LU (fe_band); SHTns/FFTW are
# called only outside the parallel regions. There is no LIS to reconcile — the
# iterative solver was removed in favour of the direct banded LU, which is also why
# this build no longer depends on a serial-vs-OpenMP LIS variant at all.

# --- Final flag sets ---------------------------------------------------------
# MODFLAGS (-I/-J objdir) and FFLAGS_BASE come from the compiler fragment.
# INC_SHTNS is what lets `include 'shtns.f03'` in src/fe_sht.f90 be found.
CPPFLAGS_FE = $(CPPFLAGS_PP)
FFLAGS_FE   = $(FFLAGS_BASE) $(MODFLAGS) $(INC_NC) $(INC_FESMUTILS) $(INC_FFTW) $(INC_SHTNS) $(INC_LIS)

# Static archives resolve left-to-right, so a library must precede the libraries
# it depends on: SHTns before FFTW (SHTns calls FFTW), fesm-utils before netCDF.
LFLAGS_FE   = $(LIB_FESMUTILS) $(LIB_SHTNS) $(LIB_FFTW) $(LIB_LIS) $(LIB_NC) $(LFLAGS_EXTRA)
