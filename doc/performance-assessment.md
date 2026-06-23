# FastEarth3D — performance assessment & roadmap

A whole-model performance review aimed at **real transient runs** (a full glacial
cycle: thousands of time steps with a migrating coastline and a laterally varying
load), *not* the one-shot validation benchmarks. It records where the time goes,
the prioritized options, what has been implemented, and what is deferred.

For the per-step measured anchor and the lmax-scaling extrapolation see
[`performance.md`](performance.md); for the per-degree solver micro-optimisations
(band LU, degree-grouped memory, skip-negligible, OpenMP over the degree loop) see
§Performance of [`design.md`](design.md). This note sits on top of both.

## Cost model

For a transient run the wall time is

```
total  =  (timespan / dt)  ×  per_step
per_step  ≈  begin_step  +  SLE_fixed_point
SLE_fixed_point  ≈  n_outer × n_inner × (~3 SHTs per inner iteration)
```

- `begin_step` (`fe_response`) does ~2 real banded solves per active `(l,m)`; the
  band LU + OpenMP-over-degree work already made this cheap (~58 ms at lmax 128).
- The **SLE fixed point's spherical-harmonic transforms now dominate.** At lmax 128
  the run is *effectively serial* and **SHT-bound** (~216 ms/step, `user ≈ real`),
  because SHTns is linked serial and the per-step cost is roughly linear in
  `n_outer × n_inner`. SHT cost scales ≈ O(lmax³).

So the two highest-value levers for transient runs are **(a) the step count** (the
`dt` lever) and **(b) the SLE iteration count × per-SHT cost**. The benchmarks hide
both: they don't pay the long step count, and their clean coastlines converge the
SLE in few iterations.

## Prioritized options (3-agent review)

Ranked by return on effort for transient runs. ✅ = implemented now; ⏭ = deferred.

| # | Option | Payoff | Effort/risk | Status |
|---|--------|--------|-------------|--------|
| 1 | **Optimization compiler flags** (`-O3 -mcpu=native -funroll-loops -ffast-math`) | 1.5–3× (free) | low | ✅ |
| 2 | **Warm-start the SLE fixed point** from the previous step's `rsl` | ~2× (cuts iteration count) | low | ✅ |
| 3 | **ETD0 / exponential memory update** → larger stable `dt` | 2–5× (step count) | med–high | ⏭ |
| 4 | **OpenMP SHTns as the default** (offline *and* coupled) | several× at lmax ≥ 256 | low | ⏭ |
| 5 | **Fuse the two syntheses** to one synthesis of `N_lm − u_lm` | ~33% of inner-loop SHTs | low | ⏭ |
| 6 | **Batched multi-RHS band solve** over all `m`/re-im at fixed `l`; kill the dense degree-1 LU via nullspace projection | high at production lmax | med | ⏭ |
| 7 | **Multi-rate cadence** — re-solve the SLE less often than the cheap memory advance | ~2× | med | ⏭ |
| — | Memory footprint ~O(lmax²): ~0.5 GB @128, ~2 GB @256, ~8 GB @512 | feasibility gate | — | note |
| — | Restart I/O: write the full `tau_*` memory state only at checkpoints, not per step | avoid a multi-GB/step footgun | low | note |

Items 1–3 alone plausibly compound to ~5–15× on a real glacial cycle.

## Implemented now

### 1. Optimization compiler flags — `config/macbook_gfortran`

The only optimization flag was `-O2`; no architecture targeting, unrolling, or
vectorization. The hot kernels (`dissipative_rhs`, `advance_memory`, the band
solve) are tight double-precision loops that benefit directly. Now:

```make
DFLAGS_NODEBUG = -O3 -mcpu=native -funroll-loops -ffast-math
```

- `-mcpu=native` is the aarch64 / Apple-Silicon idiom (the x86 `-march=native`
  equivalent); it enables NEON vectorization. On other machine fragments use the
  platform-appropriate flag (`-march=native` on x86-64).
- **`-ffast-math` is ON by design**, for parity with the CLIMBER-X production
  build (the host enables it, and the coupled build cannot reliably be configured
  differently). It relaxes IEEE semantics, so this is the one change that **must
  be validated on the target machine**: re-run `make check` and confirm the SLE
  mass-conservation tests (~1e-16) still pass. (This container has neither
  gfortran nor the `fesm-utils` deps, so the build could not be exercised here.)

### 2. Warm-start the SLE fixed point — `fe_sle.f90`, `fe_coupling.f90`

`sle_solve` previously did `rsl = 0.0` at the top of *every* call, discarding the
previous step's converged solution. For a transient run the coastline and RSL
barely move between adjacent steps, so that solution is a near-converged seed.

- New `sle_solver%warm_start` flag (default `.false.` → cold start, so **every
  existing benchmark/test is byte-for-byte unchanged**).
- `rsl` changed from `intent(out)` to `intent(inout)`; the cold-start reset is now
  guarded by `.not. warm_start`. With the flag off the reset still runs, so the
  semantics for existing callers are identical.
- The coupling driver (`solid_earth_init`) turns it **on**: `self%rsl` persists
  across sub-steps and coupling intervals (seeded to 0 at init), so each solve
  starts from the last converged field.

The converged fixed point is unique (the SLE is a Fredholm equation of the second
kind → a contraction) and mass is rebalanced every iteration via Δφ, so the
**answer is unchanged to tolerance** — only the inner iteration count drops
(typically to 1–2, and it tames the pathological strongly-migrating-coastline case
that can otherwise need tens of iterations).

## Deferred — roadmap (next sessions)

**3. ETD0 / exponential memory update (highest upside).** `dt` is currently capped
by the *stability* limit of the forward-Euler Maxwell recurrence
(`Δt ≲ 2η_min/μ`, the reason for the viscosity floor), not by accuracy. The fix is
a one-line, fully local substitution in `advance_memory`
(`fe_viscoelastic.f90:216-235`): `(1−M) → exp(−M)` and `2μM → 2μ(1−exp(−M))`. This
is unconditionally stable and exact for piecewise-constant strain; the operator,
banded LU, and the `begin_step/apply/commit` structure are untouched, and the
elastic (M→0) / fluid (μ=0) limits are preserved. Plausibly lifts `dt` from
2.5–20 yr toward 25–100 yr.

> ⚠️ **The claim that ETD0 "was tried and abandoned" has no trace in this repo** —
> no commit (`git log -S ETD`), branch, or note, only two lines of prose in
> `performance.md`/`design.md`. Treat it as unsubstantiated and revisit. The real
> risk is *accuracy* at large `dt` (time truncation through the multi-element
> coupled solve), not stability — so implement behind a flag and run a
> `dt`-convergence sweep against `test_relax` and a Martinec-2018 case before
> trusting it. This is the single highest-value experiment available.

**4. OpenMP SHTns as the default.** SHTns is currently linked serial
(`config/common.mk`) — a deliberate choice to avoid OpenMP nesting inside
CLIMBER-X. **Decision (this review): make the OpenMP SHTns variant the default for
both the standalone/offline driver and the coupled build.** The transforms are the
serial bottleneck and the per-SHT cost is O(lmax³), so threading them is the clean
lever that moves lmax 256 into comfortable range. The `fe_sht` wrapper needs no
source change — only the build wiring (link the `shtns-omp` variant; it must be
built in `fesm-utils`). The model's own degree-loop OpenMP and the SHT calls run in
*different phases* of the step (never nested concurrently), so they can share one
thread pool safely. (Requires building the OpenMP SHTns variant in `fesm-utils`,
hence deferred to a build/dependency pass.)

**5. Fuse the two syntheses.** The inner loop only consumes `Sraw = N − u`;
synthesize the combined spectrum `N_lm − u_lm` once instead of `u` and `N`
separately, cutting inner-loop SHTs from 3 to 2 (~33%). Synthesize `u`/`N`
individually only on the converged iterate, for diagnostics. Pure reordering.
(`fe_sle.f90:177-178`.)

**6. Batched multi-RHS band solve + degree-1 nullspace projection.** `begin_step`
calls the band solver one RHS at a time, twice per coefficient (re/im), but for a
fixed degree `l` the factored LU is identical across all `m` and both re/im.
Batching them into one multi-RHS triangular solve (`fe_band.f90:114`) turns
memory-bound work into cache-friendly BLAS-3-style reuse — biggest at lmax ≥ 256.
Separately, degree 1's dense KKT border makes its solve O(nr³) and serializes one
thread; replace with a Sherman–Morrison / nullspace projection on the narrow band
(`fe_radial_fe.f90`, bordered path). Validate against the degree-1 benchmark.

**7. Multi-rate cadence.** The viscous memory evolves on millennial timescales; the
expensive migrating-coastline SLE need not be re-converged at the same fine cadence
as the cheap Maxwell memory advance. The affine `begin_step/apply` split already
separates frozen drift from the current load, so the SLE can run on a coarser
cadence (e.g. per coupling interval) with the load frozen between. Complementary to
warm-starting (#2).

### Notes / guardrails (not bottlenecks today)

- **Memory footprint scales O(lmax²)** (six `(NLAM, ne, nk)` memory-stress arrays):
  ~0.5 GB at lmax 128, ~2 GB at 256, **~8 GB at 512** — a genuine feasibility wall
  at 512 on a laptop, and a real cost when this state is checkpointed. lmax 128 is
  the right coupled default; lmax > 256 would want a sparse/compressed memory layout
  (the skip-negligible machinery already avoids *solving* negligible coefficients
  but still *stores* them all).
- **Restart I/O.** The coupling `update` does no per-step I/O (good). But a full
  snapshot writes the entire `tau_*` memory state (multi-GB at high resolution).
  Enforce that heavy `tau_*` restart writes happen only at checkpoints; per-step
  diagnostics should use the small lon×lat field subset (`fe_io` already supports a
  variable subset via its `nms` argument).

## Validation status of this change set

- **Warm-start:** default-off preserves all existing benchmark/test behavior
  exactly; only the coupling path opts in. The `intent(out) → intent(inout)` change
  is benign for existing callers (the cold-start reset still defines `rsl`).
- **Compiler flags:** could not be exercised here (no toolchain/deps in this
  container). **Run `make check` on the target machine after pulling**, paying
  particular attention to the SLE mass-conservation tolerances under `-ffast-math`.
