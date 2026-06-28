#!/usr/bin/env bash
#
# run_modal_vs_ve.sh — stage/submit the experiment set that quantifies how well
# the reduced MODAL response approximates the full viscoelastic (VE) solver, in
# both accuracy and cost, using runme (https://github.com/.../runme).
#
# Two sets, both forced by the Tarasov deglaciation and measured against the VE
# run in the *same* set (the ground truth):
#
#   1. radial   — 1-D (radially-symmetric) viscosity. The clean limit where modal
#                 with n_modes=all converges to VE exactly; isolates the accuracy
#                 of the mode-count / ranking dial (no lateral approximation).
#   2. deglac3d — full run with laterally-varying (3-D) viscosity. The LGM spin-up
#                 is done with the cheap 1-D solver in every case (spinup_1d=true),
#                 then the transient runs on the 3-D path. This is where modal is a
#                 genuine approximation to VE (design §4).
#
# Each set sweeps the modal dial: earth_response=modal x n_modes x mode_rank, plus
# n_modes=all, plus the single VE reference. Compare each modal out.nc (rsl/bsl)
# against the set's ve/out.nc, and the wall time from out.out ([PROFILE] lines) or
# the SLURM accounting.
#
# Run it from anywhere; it cd's to the repo root so runme finds .runme/.
#
#   ./scripts/run_modal_vs_ve.sh
#
# Override any setting from the environment, e.g.
#   FORCING=/work/ice.nc VISC3D=/work/bagge.nc RUNME_FLAGS="-s -r" ./scripts/run_modal_vs_ve.sh
#
set -euo pipefail
cd "$(dirname "$0")/.."                      # repo root (where .runme/ lives)

# ============================================================================
# MACHINE-SPECIFIC paths — EDIT THESE for the target cluster (absolute paths).
# ============================================================================
FORCING=${FORCING:-/Users/alrobi001/models/climber-x/input/geo_ice_tarasov_deglac.nc}   # ice_thickness(lon,lat,time)
VISC3D=${VISC3D:-/Users/alrobi001/models/isostasy_data/earth_structure/viscosity/bagge2021.nc}  # log10(eta)(lon,lat,r)

# ============================================================================
# Experiment knobs.
# ============================================================================
LMAX=${LMAX:-64}                             # spherical-harmonic degree (ref file must match)
T0=${T0:--26000.0}                           # transient start [yr] (LGM)
T1=${T1:-0.0}                                # transient end   [yr] (present)
DT_COUPLE=${DT_COUPLE:-100.0}                # coupling interval [yr] (forcing cadence)
DT_EQUIL=${DT_EQUIL:-10000.0}                # LGM-memory spin-up [yr/pass]
OMP=${OMP:-8}                                # OpenMP threads per run
EXP=${EXP:-runs/modal_vs_ve}                 # experiment root (under gitignored runs/)

# Modal dial swept in each set (comma lists => runme ensemble dimensions).
N_MODES=${N_MODES:-1,2,4,8}                  # modes kept per degree (truncated)
RANKS=${RANKS:-isostatic,rate,residue}       # mode_rank metric

# How runme launches each (group) of runs:
#   "-s -r"  prepare SLURM scripts AND submit   (HPC, default)
#   "-s"     prepare SLURM scripts, do NOT submit (inspect, submit by hand)
#   "-r"     run locally in the background
#   ""       stage the run dirs only (no run)
RUNME_FLAGS=${RUNME_FLAGS--s -r}     # note: `-` not `:-`, so RUNME_FLAGS="" means stage-only

REF=data/reference/rtopo_gauss_l${LMAX}.nc   # present-day reference (i_eq=1); resolved via the rundir 'data' symlink

# Turning on 3-D viscosity + 1-D spin-up. runme writes booleans quoted ('true'),
# but the model's nml reader parses 'true'/'false' as logicals, so -p is fine.
VISC3D_ON=(fe3d.l_visc_3d=true fe3d.spinup_1d=true fe3d.visc_3d_file="$VISC3D")

# Parameters common to every run (machine paths + the shared deglaciation setup).
COMMON=(
  fe3d.lmax="$LMAX"
  fe3d.file_forcing="$FORCING"
  fe3d.name_ice=ice_thickness
  fe3d.i_eq=1
  fe3d.z_bed_ref_file="$REF"
  fe3d.h_ice_ref_file="$REF"
  fe3d.dt_couple="$DT_COUPLE"
  fe3d.dt_equil="$DT_EQUIL"
  fe3d.time_init="$T0"
  fe3d.time_end="$T1"
  fe3d.file_out=out.nc
)

# launch <outdir> <ensemble:0|1> <extra -p args...>
launch() {
  local out=$1 ens=$2; shift 2
  local aflag=()
  [ "$ens" = "1" ] && aflag=(-a)                     # name ensemble member dirs from their params
  echo ">>> $out   ($*)"
  runme -o "$out" -e main --omp "$OMP" $RUNME_FLAGS \
        ${aflag[@]+"${aflag[@]}"} \
        -p "${COMMON[@]}" "$@"
}

echo "lmax=$LMAX  window=[$T0,$T1]  dt_couple=$DT_COUPLE  dt_equil=$DT_EQUIL  omp=$OMP"
echo "exp root: $EXP    runme flags: '$RUNME_FLAGS'"

# ---------------------------------------------------------------------------
# Set 1 — idealized: 1-D (radial) viscosity. modal(n_modes=all) -> VE exactly.
# ---------------------------------------------------------------------------
launch "$EXP/radial/ve"        0 fe3d.earth_response=ve    fe3d.scheme=fe          fe3d.l_visc_3d=false
launch "$EXP/radial/modal"     1 fe3d.earth_response=modal fe3d.n_modes="$N_MODES" fe3d.mode_rank="$RANKS" fe3d.l_visc_3d=false
launch "$EXP/radial/modal_all" 0 fe3d.earth_response=modal fe3d.n_modes=-1         fe3d.l_visc_3d=false

# ---------------------------------------------------------------------------
# Set 2 — full deglaciation: 3-D viscosity, 1-D spin-up.
# ---------------------------------------------------------------------------
launch "$EXP/deglac3d/ve"        0 fe3d.earth_response=ve    fe3d.scheme=fe          "${VISC3D_ON[@]}"
launch "$EXP/deglac3d/modal"     1 fe3d.earth_response=modal fe3d.n_modes="$N_MODES" fe3d.mode_rank="$RANKS" "${VISC3D_ON[@]}"
launch "$EXP/deglac3d/modal_all" 0 fe3d.earth_response=modal fe3d.n_modes=-1         "${VISC3D_ON[@]}"

echo "done."
