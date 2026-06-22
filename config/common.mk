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

# --- OpenMP variants ---------------------------------------------------------
# `make openmp=1`: swap the serial dependency builds for their OpenMP variants.
# Note the differing library names (libshtns.a vs libshtns_omp.a, matching FFTW).
ifeq ($(openmp),1)
	INC_FESMUTILS = -I$(FESMUTILSROOT)/include-omp
	LIB_FESMUTILS = -L$(FESMUTILSROOT)/include-omp -lfesmutils

	FFTWROOT = fesm-utils/fftw-omp
	INC_FFTW = -I$(FFTWROOT)/include
	LIB_FFTW = -L$(FFTWROOT)/lib -lfftw3_omp -lfftw3 -lm

	SHTNSROOT = fesm-utils/shtns-omp
	INC_SHTNS = -I$(SHTNSROOT)/include
	LIB_SHTNS = -L$(SHTNSROOT)/lib -lshtns_omp
endif

# --- Final flag sets ---------------------------------------------------------
# MODFLAGS (-I/-J objdir) and FFLAGS_BASE come from the compiler fragment.
# INC_SHTNS is what lets `include 'shtns.f03'` in src/fe_sht.f90 be found.
CPPFLAGS_FE = $(CPPFLAGS_PP)
FFLAGS_FE   = $(FFLAGS_BASE) $(MODFLAGS) $(INC_NC) $(INC_FESMUTILS) $(INC_FFTW) $(INC_SHTNS)

# Static archives resolve left-to-right, so a library must precede the libraries
# it depends on: SHTns before FFTW (SHTns calls FFTW), fesm-utils before netCDF.
LFLAGS_FE   = $(LIB_FESMUTILS) $(LIB_SHTNS) $(LIB_FFTW) $(LIB_NC) $(LFLAGS_EXTRA)
