# FastEarth3D — performance note

> ⚠️ **Superseded numbers (see [`performance-assessment.md`](performance-assessment.md)
> §Measured results).** The 162.5 s / 216 ms-step anchor below and the
> "effectively serial" diagnosis do **not** reproduce on a controlled rebuild:
> measured `-O2` E2 is **58.8 s / 78 ms-step**, and the run **parallelizes ~4×**
> (`user/real ≈ 4`), not serial. The lmax-scaling table here is therefore optimistic
> on absolute times and should be rebased before use. Kept for the scaling *shape*
> and the caveats, which still hold.

End-to-end timing of the coupled solid-earth + sea-level solver, and an
extrapolation to production-scale global domains. For the per-degree solver
micro-optimisations (band LU, degree-grouped memory, skip-negligible, OpenMP over
the degree loop) see the "Performance" section of `doc/design.md`; this note is
about the *whole-step* cost and what it means for real runs.

## Measured anchor

Martinec-2018 case **E2** (the full migrating-coastline SLE: subgrid coast, basin
topography, a growing off-pole ice cap), run as `test_benchmark_sle.x E2`:

| quantity | value |
|---|---|
| resolution | lmax = 128 (≈156 km, grid 256×512) |
| steps | 750 (T1 history, dt = 20 yr) |
| wall time | **162.5 s** → **~216 ms / step** |
| threads | `OMP_NUM_THREADS=8` |
| machine | 10-core Apple-Silicon laptop (4 performance + 6 efficiency) |

**Key observation:** `user` time (159.6 s) ≈ `real` time (162.5 s). The run is
**effectively serial**. The band-LU work made the per-degree solves cheap
(`begin_step` ≈ 58 ms at lmax 128), so the bottleneck has moved to the **SLE
fixed-point's spherical-harmonic transforms** — roughly 3 SHTs (one analysis, two
syntheses) per inner iteration, up to `n_outer × n_inner` iterations per step — and
SHTns is linked in its **serial** variant (deliberately, to avoid OpenMP nesting
inside a host model). So at lmax 128 the model is SHT-bound and single-core-bound;
the OpenMP over the degree loop barely moves the wall clock at this resolution.

## Extrapolation to a global domain

The dominant cost (the SHTs) scales as ≈ O(lmax³) on a Gauss grid; the step count is
(timespan / dt). Taking the measured 216 ms/step at lmax 128 and scaling the SHT
term by (lmax/128)³:

| resolution | ≈ ms/step | glacial cycle (120 kyr, dt 50 yr → 2400 steps) | deglaciation (20 kyr, dt 20 yr → 1000 steps) |
|---|---|---|---|
| **lmax 128** (~156 km) | 216 (measured) | **~9 min** | ~3.6 min |
| lmax 256 (~78 km)  | ~1.7 s (×8)  | ~1.2 h  | ~30 min |
| lmax 512 (~39 km)  | ~14 s (×64)  | ~9 h    | ~4 h |

(At dt = 20 yr the glacial-cycle numbers are ~2.5× larger.)

## Verdict

- **lmax 128 (~150 km) is comfortably production-fast** — a full glacial cycle in
  ~10–20 min wall on a laptop. For coupling into CLIMBER-X (a fast EMIC that runs
  glacial cycles in hours), the solid-earth overhead is negligible, and ~150 km is a
  reasonable GIA resolution.
- **lmax 256 (~78 km): fine for standalone / offline (~1 h/cycle), acceptable but
  noticeable when coupled.**
- **lmax 512 (~39 km): heavy (~hours–day); would need optimisation before routine
  use.**

## Caveats and levers (in priority order)

1. **SLE iteration count on real coastlines — the biggest uncertainty.** Per-step
   cost is roughly linear in `n_outer × n_inner`. The benchmark converges the
   migrating coastline in `n_outer = 3`; the pure-eustatic subgrid test
   (`test_sle_subgrid`, null response, strongly-migrating shallow basin, no VE
   stabilisation) needed up to ~60. A real, complex coastline likely sits between,
   so the table above could move 2–3×. Measure this on a real ICE-6G-style run
   before trusting the extrapolation.
2. **It is serial-bound, so a bigger machine does not help as-is.** The clean lever
   is **parallelising the SHTs**: SHTns has an OpenMP variant we currently avoid for
   host-nesting reasons, but for a *standalone* driver it could be switched on for a
   several-× speedup (pushing lmax 256 into comfortable range). Cutting the SLE
   iteration count (tighter convergence, the skip-negligible trick already in
   `ve_response`) is the complementary lever.
3. **Time step / stability.** Larger dt cuts the step count proportionally, but the
   explicit forward-Euler Maxwell scheme has a stability limit (the ETD0 /
   exponential-memory alternative was tried and abandoned — see the project notes),
   so the usable dt must be checked against the fastest relaxation mode at the
   target resolution.

## Reproducing

```sh
rm -rf obj && python config.py config/macbook_gfortran && make openmp=1 test_benchmark_sle
OMP_NUM_THREADS=8 /usr/bin/time -p bin/test_benchmark_sle.x E2
```

A firmer scaling curve would come from recompiling at lmax 256 (a compile-time
parameter in `tests/test_benchmark_sle.f90`) and timing a short segment, and/or
instrumenting `begin_step` vs the SLE-loop SHT time to confirm the serial-bound
diagnosis directly.
