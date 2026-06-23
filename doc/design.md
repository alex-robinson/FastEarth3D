# FastEarth3D — design

This document records the design decisions, the method comparison that led to
them, the validation plan, and the implementation pitfalls the GIA literature
flags. It is the reference for filling in the module stubs.

## 1. Goal and constraints

Replace **VILMA** (the closed-source solid-Earth model coupled to CLIMBER-X)
with an open-source Fortran model that is:

- **state of the art** — full sea-level equation, rotational feedback, 3D
  (laterally varying) viscosity;
- **simple and fast** — no normal-mode root finding, no heavy 3D-FEM
  infrastructure;
- **3D-ready from the start**, validated first in 1D.

We reproduce **VILMA's physics** (confirmed from Albrecht, Bagge & Klemann 2024,
*The Cryosphere* 18:4233, and Bagge et al. 2021): an **incompressible,
self-gravitating, Maxwell** viscoelastic sphere, solved by the
**spectral–finite-element, time-domain** method of **Martinec (2000)**, *GJI*
142:117, with rotational feedback (Martinec & Hagedoorn 2014) and a full
migrating-coastline sea-level equation. This is a clean-room reimplementation
from the published literature — we do not read VILMA source.

## 2. Why this method (landscape)

| Family | Examples | 3D viscosity | Speed | Fit |
|---|---|---|---|---|
| Normal-mode + Love numbers + SLE | SELEN, TABOO, ALMA3 | no (1D only) | fast (1D) | dead end for a 3D model |
| **Time-domain spectral-FE** | **VILMA** (Martinec 2000) | **yes** | **fast** | **chosen** |
| Full 3D FE/FV | CitcomSVE, ASPECT, Elmer | yes | slow, heavy | too complex |

The decisive structural fact: with the **explicit** Maxwell time scheme, the
stiffness solve has **no angular coupling**, so each spherical-harmonic degree
`l` decouples into an independent small banded radial system — embarrassingly
parallel and cheap. Lateral (3D) viscosity does **not** couple the harmonics in
the solve; it enters only through the explicitly-known **memory stress**,
evaluated pointwise on the Gauss-Legendre grid and transformed back to spectral
each step (pseudo-spectral). Therefore the **1D and 3D code paths are the same
solver** — 3D viscosity is just a spatially varying field instead of a constant.
That is what makes "3D-ready from day one, validated in 1D" cheap rather than a
rewrite.

## 3. Architecture (modules)

| Module | Role | Status |
|---|---|---|
| `fe_precision` | working precision (= C double) | done |
| `fe_constants` | physical / reference constants (benchmark conventions) | done |
| `fe_sht` | SHTns wrapper — the transform kernel | done + tested |
| `fe_earth_structure` | radial layers + optional 3D viscosity field | types |
| `fe_radial_integrals` | Appendix C P1/P0 element integrals | done + tested |
| `fe_band` | pivoted banded LU (dependency-free, re-entrant) | done + tested |
| `fe_radial_fe` | per-degree saddle-point operator + banded-LU solve | done + tested |
| `fe_viscoelastic` | Maxwell memory-stress explicit time stepping (1-D) + shared kernel | done + tested |
| `fe_response` | surface-load → (uplift, geoid) operator: elastic + viscoelastic field driver | done + tested |
| `fe_gravity` | self-gravitation / Poisson coupling | stub |
| `fe_sle` | sea-level equation (ocean function, migration) | done + tested (elastic + VE) |
| `fe_rotation` | rotational feedback / TPW | stub |
| `fe_coupling` | CLIMBER-X-compatible init/update/finalize API | stub |
| `fastearth` | umbrella re-export | done |

Convention in `fe_sht`: fully-normalized real harmonics, no Condon-Shortley
phase, Gauss grid, phi-contiguous layout; spectral arrays hold `m >= 0`.

## 4. Coupling contract (from CLIMBER-X `src/geo/vilma.F90`)

- **in:** ice thickness `h_ice` [m]; **out:** relative sea level `rsl` [m] and
  bedrock `z_bed = z_bed_eq - rsl` [m].
- Host (CLIMBER-X) maps between its grid and the model's lat-lon Gauss grid with
  conservative (in) / bilinear (out) SCRIP weights; **all SH transforms stay
  inside this model**.
- Reference VILMA resolution: SH degree 170; SLE grid 1024×2048; radial FE 5 km
  (→420 km) / 10 km (→670 km) / 40–60 km (→CMB); explicit Δt 2.5 yr; coupling
  cadence `n_year_geo` (default 10 yr); viscosity floor ~1e19 Pa·s.

## 5. Validation ladder (each rung = a published benchmark)

1. **SH transform** round trip — `tests/test_sht.f90` (done).
2. **Love numbers** (elastic + viscoelastic) vs **Spada et al. (2011)**, *GJI*
   185:106, tests 2/1, 3/1, 5/1, 7/1 on Earth model **M3–L70–V01**; cross-check
   **ALMA3**.
3. **Deformation / geoid** for a disc load vs Spada (2011) tests 1/2–2/2.
4. **Sea-level equation** vs **Martinec et al. (2018)**, *GJI* 215:389, cases
   A→E (the migrating-coastline benchmark VILMA itself passed).
5. **Rotation** vs Spada (2011) test 3/2 (polar motion).
6. **3D** vs ASPECT / TABOO cross-checks (no closed published 3D benchmark yet;
   Klemann et al. 2022 effort ongoing).

## 6. Implementation pitfalls to design in (not patch later)

- **Degree-1 is frame-dependent** — do not hard-code `h1=l1=1, k1=0`. Work in the
  **CM frame** (Blewitt 2003). Geocenter motion *is* the degree-1 signal.
- **Rotation is degree-2** and uses **tidal** Love numbers; fix the fluid Love
  number `k_f` from the observed flattening (Mitrovica et al. 2005) to avoid the
  lithosphere-thickness paradox.
- **Lithosphere = η→∞ exactly**, never large-but-finite (a finite value injects a
  spurious slow relaxation mode that contaminates multi-cycle runs).
- **Pseudo-spectral SLE**: convolution in spectral space, ocean-function product
  on the spatial grid — this is what kills Gibbs ringing at coastlines.
- **Conserve mass, not volume**: apply `rho_ice/rho_water` exactly once each in
  the flotation criterion, the eustatic conversion, and the SLE; fix one density
  set (benchmark: ρ_ice = 931, ρ_water = 1000).
- The SLE is a **Fredholm equation of the second kind** → nested fixed-point
  iteration (~3 outer × 3 inner; Kendall, Mitrovica & Milne 2005).
- **Explicit-scheme stability**: enforce a viscosity floor; the time step is
  bounded by the shortest Maxwell time.

## 7. Build system

`configme`-style: `config.py` inserts a machine fragment
(`config/<machine>_<compiler>`) into the template `config/Makefile`; shared
dependency wiring is in `config/common.mk`, source/rule lists in
`config/Makefile_fastearth.mk`. Dependencies come from a `fesm-utils` symlink
(FFTW, SHTns, `fesmutils`) plus system netCDF (`nf-config`). SHTns was added to
fesm-utils' `build.py` as a first-class component (serial + OpenMP variants).

## 8. Key references

- Martinec (2000), *GJI* 142:117 — the spectral-FE method.
- Albrecht, Bagge & Klemann (2024), *The Cryosphere* 18:4233 — PISM–VILMA,
  confirms VILMA physics.
- Martinec & Hagedoorn (2014), *GJI* 199:1823 — time-domain rotational feedback.
- Kendall, Mitrovica & Milne (2005), *GJI* 161:679 — SLE with moving shorelines.
- Spada et al. (2011), *GJI* 185:106 — GIA benchmark.
- Martinec et al. (2018), *GJI* 215:389 — SLE benchmark.
- Mitrovica et al. (2005), *GJI* 161:491 — revised rotation theory.
- Blewitt (2003), *JGR* 108:2103 — degree-1 / reference frames.

## 9. Radial solver formulation (rung 2) — working notes

Verified from Martinec's *Continuum Mechanics* lecture notes (Ch. 9), the Hanyk
PhD thesis, and the VEGA benchmark description (Martinec et al. 2018). Items
marked **(unverified)** await the paywalled Martinec (2000) PDF.

**Governing equations** (incompressible, self-gravitating, quasi-static):
- Incremental momentum balance `Div t^L + ρ₀ f^L + ρ^L f₀ + (pre-stress) = 0`,
  body force `f^L = −∇φ₁ + ∇u·∇φ₀` (Martinec Eq. 9.32, 9.152).
- Poisson via advected density: `∇²φ₁ = 4πG ρ₁`, `ρ₁ = −∇·(ρ₀u)` (Eq. 9.98).
- Incompressibility `Div u = 0`; isotropic stress Π is the Lagrange multiplier
  (K→∞), deviatoric part `τ^dev = 2μ ε^dev` (Maxwell in time, §D below).

**Spectral/radial:** spheroidal-only for 1D loading. Radial scalar set per degree
`{U, V, F}` nodal + pressure Π; tractions are natural BCs of the weak form.
**P1 ("tent") radial basis** (VEGA: 165 elements). Mesh = done (§3, `fe_radial_fe`).

**Incompressibility discretization (B4) — RESOLVED.** Confirmed from the
Martinec (2000) PDF: **P1 (U,V,F nodal) / P0 (Π per element)**, `Div u = 0`
enforced weakly via the pressure block (eq 82). Assembled and validated — see
[formulation.md](formulation.md) and `fe_radial_fe%build_dense_operator`.

**Boundary conditions — settled (whole-sphere mesh).** The fluid core IS meshed
(μ=0 region), so there is **no explicit CMB boundary condition** — free-slip
emerges. No explicit centre BC either: Martinec meshes through r=0 and the r²
weighting handles regularity (the singular `I⁷` term vanishes via `R₁=0`).
Surface: the (j+1)F(a) exterior match on the F-F diagonal + the load RHS.
Density-jump interfaces are natural conditions of the weak form. (The earlier
Wu & Peltier CMB-BC plan assumed an un-meshed core and is superseded.)

**Time scheme (rung 3) — DONE (1-D).** Explicit ω=1 Maxwell scheme (Martinec
2000 eqs 23-25), implemented in `fe_viscoelastic%ve_degree`. Memory stress
`τ^{V,i} = (1−M)τ^{V,i-1} − 2μM ε^i`, `M = μΔt/η`; it enters the SAME elastic
operator as the dissipative RHS forcing `−∫τ^V:δε dV` (radial Gauss-2 + the
spectral double-dot over the four spheroidal tensor components, eqs 94/110).
Elastic layers (η→∞) freeze, fluid layers (μ=0) carry no memory. Stability
`Δt ≲ 2η_min/μ` ⇒ viscosity floor; VEGA Δt = 20 yr. The fixed operator is built
and LU-factored once (`fe_band`, banded) and reused every step (~20 µs/solve).
Validated (`test_relax`): held load relaxes elastic→fluid, `t_relax ∝ η`.

**Love numbers (rung 2) — DONE.** `h_n = g₀ U(a)/φ^L`, `l_n = g₀ V(a)/φ^L`,
`k_n = −F(a)/φ^L − 1`, with `φ^L = 4πG a σ/(2n+1)` (`fe_radial_fe%loading_love`).
The fluid (`h_n→−(2n+1)/3`, `k_n→−1`) and rigid (`→0`) limits are reproduced to
~1e-5 (`test_love`). `l_n` sign/normalization and the quantitative Spada (2011)
Test 2/1 match are the remaining calibration items.

## 10. Sea-level equation (rung 4) — working notes

**Response operator (`fe_response`).** The SLE is built on an abstraction:
`response_operator%apply(σ_lm) → (u_lm, N_lm)` maps a spectral surface mass load
to surface uplift `u` and geoid height `N`, so the elastic and viscoelastic Earth
responses are swappable. **Geoid mapping: `N(a) = −F(a)/g`** (Bruns' formula; the
geopotential perturbation is `−F` since Martinec's `φ₁ = F → −φ^L` rigid). This
uses only `U` and `F`, *not* the unresolved horizontal Love number `l`, so the
SLE is not blocked by that open item. Three implementations:
- `elastic_response` — per-degree gains `gu(l)=U(a)`, `gn(l)=−F(a)/g` precomputed
  once from a unit-load solve;
- `null_response` — `u≡N≡0` (rigid, non-self-gravitating) → eustatic baseline;
- `ve_response` — viscoelastic field driver (below).

**SLE fixed point (`fe_sle`) — DONE.** Pseudo-spectral (KMM 2005): the load →
response convolution is spectral, the ocean-function product `C·S` is pointwise on
the Gauss grid (kills coastline Gibbs ringing). For a fixed coastline,
`S = C·(N − u + Δφ)`, `N,u` = response to `L = ρ_i ΔI + ρ_w C·S`, and the spatial
constant `Δφ` is fixed *each iteration* by ocean-mass conservation
`ρ_w ∫C·S dA = −ρ_i ∫ΔI dA` ⇒ mass conserved to machine precision by construction.
Inner loop iterates `S`; outer loop migrates the coastline from `topo0 − S`.
Validated (`test_sle`): eustatic limit → uniform barystatic rise (mass resid
~4e-16); self-gravitating elastic M3 → mass resid 0, fixed point converges, real
spatial structure.

**Viscoelastic field driver (`ve_response`) — DONE.** The SLE needs a *field*
load with an **independent memory history per (l,m)** (each coefficient has its
own load time-series). Key insight: the response at a fixed time is **affine** in
the current load, `u_lm = gu(l)·σ_lm + drift_lm`, where `drift_lm` comes from the
*frozen* past-relaxation memory — so the SLE fixed point can call `apply()`
repeatedly without corrupting state. The step is bracketed: `begin_step` freezes
the drift (one memory-forcing solve per (l,m)), `commit_step` advances the Maxwell
memory with the converged load. Each complex (l,m) history is two real histories
(re/im) since the operator and `M = μΔt/η` are real. The per-element Maxwell
kernel (`strain_coeffs`, `ve_strain_constants`, `dissipative_rhs`,
`advance_memory`) is shared with the 1-D stepper — no duplicated algorithm.
Validated (`test_ve_response`): first-step gains == `elastic_response` exactly; a
held degree-2 load reproduces the 1-D `ve_degree` history to ~1e-13 relative.
End-to-end (`test_sle_ve`): a held 2 km ice cap → ocean mass conserved ~1e-16 per
step, eustatic mean held, ~16 m of viscoelastic relaxation over 500 yr.

**Degree-1 (geocenter) is now carried in the field driver — DONE.** Previously
skipped (`gu(1)=gn(1)=0`, `l < 2` guards) as a workaround for the *dense* j=1
operator hanging the stepper. The sparse KKT degree-1 operator (merged from
`degree1-sparse`) fixed that at the source — j=1 solves in ~2 GMRES iters in the
CM frame (wᵀd = 0 rigid-translation removal, Blewitt 2003). The `l < 2` guards
in `fe_response` are now `l < 1`: degree 1 joins the unit-load loop in
`ve_response_init` (assembles `ops(1)`, real gains + nodal fields like any
degree), and carries memory drift through begin/apply/commit; only degree 0
(monopole geoid, no deformation, no operator) stays special. Validated
(`test_ve_response`): degree-1 first-step gains == `elastic_response` exactly
(0.0 diff), and the held (1,0) field response reproduces the 1-D `ve_degree` at
j=1 to ~1e-11 relative; `test_sle_ve` mass conservation unchanged (~6.6e-16).
The frame is inherited from the operator (same CM frame as `elastic_response`
and `ve_degree`), so the SLE's `N=−F/g` at degree 1 is self-consistent with the
displacement `u`.

**Open / next** (priority order for a fresh session):
0. ~~Fix the elastic low-degree Love-number bug.~~ **DONE.** Benchmark data
   in-repo (`data/benchmarks/love_M3-L70-V01/`; `test_benchmark_love`). Root cause:
   a transposed index in the self-gravity potential-gradient force in
   `build_dense_operator` (the U-F coupling used `I²_αβ` instead of `I²_βα`),
   which broke the U↔F symmetry the energy functional requires. Found by
   re-deriving the continuous gravity form (Martinec eq 65) term-by-term and
   cross-checking the shear block against the `fe_viscoelastic` strain
   representation (eqs 85–88). Fix: `i2(ia,ib) → i2(ib,ia)`. Now elastic AND fluid
   M3-L70-V01 Love numbers match the benchmark to ~0.1% at every degree; the
   operator is exactly symmetric (`test_assembly`). Closes the disc offset (rungs
   2/3) at the source — a direct disc re-run to confirm <1% is a quick follow-up.
   (The earlier "lift the degree-1 skip in `ve_response`" item is now DONE on main
   — the geocenter degree-1 response is carried in the field driver.)
1. **`fe_coupling` wiring** — the CLIMBER-X contract (reference state: z_bed_eq +
   reference ice/topo; host-grid mapping). A deliberate interface decision; swap
   the `visco(:)` member for a `ve_response`, drive `sle%solve` per Δt across the
   coupling interval, return `z_bed = z_bed_eq − rsl`.
2. **Grounded-ice flotation** in the ocean function (currently `topo < 0` only):
   a cell is ocean only where it is below sea level AND ice does not ground
   (`ρ_i I < −ρ_w·topo`).
3. **Martinec et al. (2018) cases A–E** quantitative match — the REFERENCE
   output curves are now in-repo (`data/benchmarks/sle_martinec2018/`, cases
   A/C2/D3/E2/F1, figs 10–13: u, v_θ, v_φ, F, sea-surface, SLE). The load/topo
   SPEC is analytic (giapy `tests/sle_test.py`: ice L1–L3 at (θ₀,φ₀,h₀), topo
   B0–B3 exponential basins, time T1–T3) — build these inputs and compare. The
   elastic bug (item 0) that fed the SLE response is now fixed, so this is unblocked.
4. **Performance — DONE.** `begin_step` does 2 real solves per (l,m) per step
   (~O(nlm) solves), the cost driver at VILMA resolution. Four changes, all exact
   (results unchanged) except the threshold-controlled skip:
   - **Banded LU** (`fe_band`) replaces the iterative solver (LIS GMRES+ILU) on the
     per-degree solve, and **LIS is removed entirely**. The operator is banded
     (half-bandwidth ~6) and the equilibrated system is effectively direct (1 GMRES
     iter), so a pivoted band LU — factor once at assemble, band solve per RHS — is
     far faster (~20 µs vs ~700 µs/solve) and cache-light (no ILU fill to evict).
     Pivoting is required (zero pressure (Π,Π) block). j=1 carries the dense KKT
     border (rigid-mode removal), so that one degree has ~full bandwidth and factors
     as a dense LU — still `fe_band`, just wide. Dependency-free, **re-entrant** (so
     no serial-vs-OpenMP LIS variant to reconcile when linked into a larger host).
   - **Degree-grouped storage** of the per-(l,m) Maxwell memory/drift (slot `k`,
     `lm↔k` map) so the loop is contiguous and per-degree.
   - **Skip-negligible**: coefficients whose memory is < `skip_tol`×max are not
     solved (drift ≈ 0); ~2× for a localized cap.
   - **OpenMP** over the degree loop (`make openmp=1`, serial deps + `-fopenmp`);
     safe via the re-entrant band LU. ~5.4× at 8 threads.
   Net: `begin_step` at lmax 128 ≈ 7.9 s (orig LIS, lmax 64 extrapolated ~30 s) →
   ~58 ms (8 threads) — well over 100× single-thread from the band LU alone.
