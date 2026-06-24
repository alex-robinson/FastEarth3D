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
| `fe_params` | `fe_param_class` + `fe_par_load` (one `&fe3d` nml group) | done + tested |
| `fe_sht` | SHTns wrapper — the transform kernel | done + tested |
| `fe_earth_structure` | radial layers (`build_earth`: named / custom) + optional 3D viscosity | done + tested |
| `fe_radial_integrals` | Appendix C P1/P0 element integrals | done + tested |
| `fe_band` | pivoted banded LU (dependency-free, re-entrant) | done + tested |
| `fe_radial_fe` | per-degree saddle-point operator + banded-LU solve | done + tested |
| `fe_viscoelastic` | Maxwell memory-stress time stepping (1-D, scheme-pluggable) + shared kernel | done + tested |
| `fe_response` | surface-load → (uplift, geoid) operator: elastic + viscoelastic field driver | done + tested |
| `fe_gravity` | self-gravitation / Poisson coupling | stub |
| `fe_sle` | sea-level equation (ocean function, migration) | done + tested (elastic + VE) |
| `fe_timestep` | adaptive-Δt controller (step-doubling on the Maxwell memory) | done + tested |
| `fe_rotation` | rotational feedback / TPW (degree-2 Liouville + tidal VE channels, SLE-coupled) | done + tested (5a/5b/5c) |
| `fe_coupling` | CLIMBER-X-compatible init/update/finalize API (adaptive-coupled) | done + tested |
| `fe_io` | netCDF restart + step output (yelmo variable-table convention) | done + tested |
| `fe_drive` | standalone forced-run loop (`program fastearth`) | done + tested |
| `fastearth3d` | umbrella re-export (incl. `fe_param_class`, `fastearth_run`) | done |

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
5. **Rotation** vs Spada (2011) test 3/2 (polar motion) — **DONE** (`test_rotation`,
   `test_rotation_sle`): degree-2 Liouville polar motion `|m(t)|` matches Table 14
   (Cw=0) to <1% for cap + disc at t=0–20 kyr; the SLE-coupled feedback (5c) matches
   Adhikari et al. (2016) eq. 8 to 1e-14 and conserves ocean mass. See §11.
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
5. **Adaptive time stepping (§3c) — DONE.** The 2nd-order lever is the
   **trapezoidal** memory rule (Crank–Nicolson, order 2.00 vs FE order 1; the
   coupling iteration is its implicit solver) with the start-of-step load `σ_n`
   tracked. `fe_timestep`'s `adaptive_stepper` crosses a coupling interval with the
   ice load linearly interpolated, choosing Δt by step-doubling on the Maxwell
   memory. Δt enters only as `Mk=(μ/η)Δt` (cheap rescale, no operator re-factor).
   Pays off on dynamic-range (ice-age) loads; ~1.6× wall there. See
   `doc/performance-assessment.md` §3c.
6. **Parameter type + nml + standalone driver — DONE.** One `fe_param_class`
   loaded from a single `&fe3d` namelist (`fe_params`, yelmo `defaults_file`
   overlay; `fastearth.nml` is the complete defaults; time fields in years→s).
   `build_earth(p)` (named built-in / custom layers). `solid_earth%init(p, …)` /
   `update(h_ice, dt)` distribute the knobs and run the adaptive controller per
   interval (fixed substeps removed). Restart persists the **full** integrator
   state (`dt_try` + `σ_n`) → bit-for-bit continuation. Umbrella module renamed
   `fastearth`→`fastearth3d`; standalone `program fastearth` (`fe_drive`) runs a
   forced simulation from `&fe3d` + an ice forcing already on the Gauss grid.
7. **Real ice forcing + lon-lat → Gauss remapping — DEFERRED (after the physics).**
   The standalone driver assumes inputs already on the model Gauss grid. Real
   forcing (e.g. CLIMBER-X `geo_ice_tarasov_deglac.nc`, 1×1° deglaciation ice+bed)
   needs a conservative regridding layer (fesm-utils `mapping_scrip`) — reusable
   for any lon-lat input, but built later. **Physics first.**

**Rung 6 (3D laterally-varying viscosity) — IN PROGRESS (6a DONE; §12).** Rung 5
(rotation) is DONE — 5a/5b/5c (§11).

## 11. Rotational feedback / TPW (rung 5) — working notes

**Formulation (Spada et al. 2011 §2.1.1; time-domain à la Martinec & Hagedoorn 2014,
i.e. VILMA).** Equatorial polar motion `m = m₁ + i m₂` from the GIA (quasi-static)
Liouville equation with the Chandler wobble neglected (eq. 7, justified since the
Chandler period ≪ GIA timescales):

    [1 − k^T(t)/k_s] ∗ m(t) = Ψ_L(t),   Ψ_L = I/(C−A),   I = [δ+k^L] ∗ I_rigid,

with `k^T,k^L` the degree-2 tidal/loading Love numbers and `k_s ≡ k^T_f` the secular
(fluid) tidal Love number (eq. 11). No explicit Ω — it is absorbed into `m = ω/Ω`
and `k_s`. Validated against the Cw-excluded column of Table 14 (the regime this
quasi-static form models; the Cw≠0 Chandler transient is deliberately out of scope).

**5a — tidal forcing path (`fe_radial_fe`).** An external degree-j potential reuses
the SAME per-degree operator as a surface load; only the natural surface term differs.
In Martinec's φ₁ convention the external potential couples to F(a) with the SAME sign
as the load's own potential, `−(a/4πG)(2j+1)φ_t`, but with NO U-traction (a load
subsides; a tide-raising potential uplifts). `tidal_rhs` + `tidal_love`
(`k^T=−F/φ_t−1`). Validated (`test_tidal`) against the homogeneous-sphere Kelvin
tidal Love numbers: fluid `k^T_f=3/(2(n−1))`, `h^T_f=(2n+1)/(2(n−1))`; rigid →0;
degree-2 elastic `(3/2)/(1+μ̃)`, `μ̃=19μ/(2ρga)`.

**5b — Liouville solve (`fe_rotation`), self-contained, degree-2 only.** Two compact
degree-2 *complex* viscoelastic channels (reusing the `fe_viscoelastic` Maxwell
kernel; no normal modes, no convolution quadrature) carry the convolutions as memory:
a LOADING channel returns `(1+k^L)∗I_rigid = I(t)`; a TIDAL channel (`tidal_rhs`,
forced by the centrifugal potential ∝ `m`) returns `k^T∗m`. The feedback makes the
step ALGEBRAIC in `m` (the affine begin/apply/commit structure of the field driver):

    m_n = [ Ψ_L,n − dF_tidal/k_s ] / [ 1 − k^T_e/k_s ],

then both channels' memory is advanced. The rigid inertia `I₁₃+iI₂₃ =
−a⁴∫σ sinθcosθ e^{iφ}dΩ` is a DIRECT Gauss-grid quadrature of the load (3-D-ready —
no spherical-harmonic normalization assumption; verified by reproducing the paper's
published `G_cap/G_disc` to <0.5%). `k_s = k^T_f` from fluidizing the Maxwell mantle.

**Pitfall (designed-in): the lithosphere-thickness paradox (Mitrovica et al. 2005).**
The secular polar-motion slope is *pathologically* sensitive to `k_s` — a 0.5% change
in `k_s` moved the t=20 kyr `|m|` by ~6%. So `k_s` must be the consistent fluid limit:
fluidize ONLY `RHEOL_MAXWELL` layers (the viscous mantle); keeping the lithosphere
elastic is essential (fluidizing it by mistake inflated `k^T_f` 0.967→0.975 and
under-drove the late-time motion ~11%). For deep-time runs `k_s` is a parameter,
overridable to the observed-flattening value (§5c). Validated (`test_rotation`): cap + disc
`|m(t)|` match Table 14 (Cw=0) to <1% at t=0–20 kyr; `k^T_e=0.303`, `k_s=0.967`.

**5c — feedback into the SLE.** The centrifugal potential of `m` perturbs the sea
surface and deforms the solid, adding a degree-2 contribution to relative sea level
`s_rot = N_rot − u_rot` with (Adhikari et al. 2016, eq. 8) `N_rot = (1+k^T)Λ/g`,
`u_rot = h^T Λ/g`, `Λ = Ω²a² sinθcosθ(m₁cosφ+m₂sinφ)`. The VE `(1+k^T)`,`h^T` reuse
the 5b tidal channel (which now exposes the uplift readout); the rotational fields
reuse the `m`-forced channel exactly. `s_rot` enters the SLE geometry (`Sraw`) but NOT
the surface mass load — the rotational potential forces the Earth through the tidal
channel, not as a load — and `Δφ` is recomputed so ocean mass stays conserved. The
rotation ↔ SLE coupling is a fixed point (`m` responds to the ice+ocean load);
`fe_rotation` is split into begin_step / solve_m / s_rot / commit (affine, no memory
advance until commit) so it can be iterated. In the coupling driver it is applied at
the interval level (a predictor: `s_rot` held across the interval, `m` refreshed from
the end load; the explicit-FE channels are sub-stepped to the Maxwell stability
ceiling `dt_fe_max`). `k_s` exposes two values: the model fluid limit `k_s_fluid`
(default, reproduces Spada) and the observed-flattening closed form
`k_s_flat = 3G(C−A)/(a⁵Ω²) = 0.943` (Adhikari/Mitrovica — the recommended deep-time
value). Validated (`test_rotation_sle`, elastic): hook-off is bit-for-bit the
no-rotation SLE; `s_rot` matches Adhikari eq. 8 pointwise to 1.9e-14; ocean-mass
residual 3e-16; the fixed point converges (ocean feedback on `m` ~0.7%);
`|s_rot| ≈ 1.9 m`. End-to-end (`test_coupling`): an off-axis cap drives `|m| = 0.30°`,
mass conserved, the bed shifts up to 6.7 m vs rotation-off.

## 12. 3D laterally-varying viscosity (rung 6) — working notes

**The one structural change.** Lateral viscosity makes the Maxwell factor
`M = μΔt/η` a *field* `M(θ,φ)`, so the memory update `τ⁺ = (1−M)τ − 2μM·ε` has a
pointwise lateral product that *couples* harmonics. Everything else is untouched:
the per-degree elastic solve and the dissipative RHS `−∫τ^V:δε` are linear in the
memory (no lateral product), so they stay the exact 1-D code path. Only **η** varies
laterally — **μ and ρ stay radial** — which is why the operator/LU factorisation and
the whole speed argument survive. This is the "3D-ready from day one" payoff (§2).

**6a — pseudo-spectral memory advance — DONE (`test_response_3d`).** A new path
`advance_memory_3d` (in `fe_response`, which owns the SHT grid + slot↔lm map;
`fe_viscoelastic` stays grid-unaware) replaces the per-slot scalar advance when
`lat_visc` is set. Per radial element `e` and tensor component `λ`, each radial
shape-coefficient (A,B,C) of the memory `τ` and of the current strain `ε` (built
from the nodal displacement `σ·xUn + drift` per slot, exactly as the FE path) is
**synthesised** to the Gauss grid, advanced pointwise `τ⁺ = (1−M)τ − 2μM·ε` with the
lateral field `Mk3(:,:,e)`, and **analysed** back to spectral. Cost ≈ 36·ne SHTs per
step (NLAM=4 × {A,B,C} × {2 synth + 1 analysis}); serial for now (OpenMP over the
element loop is the perf follow-up). Degree 0 carries no memory slot and is dropped
(no deformation channel to act on). The lateral field is set via
`ve%enable_lateral_visc(sht, pert_elem)` — a per-element `(nphi,nlat,ne)` log10
viscosity perturbation, `η_eff = η·10^pert` ⇒ `MkPerDt3 = (μ/η)·10^(−pert)`; elastic
/fluid elements keep `MkPerDt=0` so the **lithosphere stays exactly elastic**. `set_dt`
rescales `Mk3 = MkPerDt3·Δt` like the 1-D `Mk`. Only the explicit FE scheme is wired
in 3D; the implicit trapezoidal 3D path is a guarded `error stop` (later sub-step).
Validated: a laterally-UNIFORM field reduces the pseudo-spectral advance to the 1-D
advance (SHT round-trip is exact on band-limited fields) — zero perturbation matches
the 1-D `ve_response` memory + uplift trajectory to ~5e-13 rel; a uniform `p` matches
a 1-D run with Maxwell η scaled by `10^p` to ~1e-12 rel.

**Open / next.**
- **6b** — lateral viscosity field builder (synthetic low-viscosity zone) + the
  FastIsostasy test-3b cross-check: a cylindrical disc load (R=100 km, H=100 m,
  ramped 100 yr) over an LVZ (η=1e19 in a 100-km column vs 1e21, 70–170 km depth),
  comparing central-uplift vs a FastIsostasy.jl test3b run (and ASPECT/Abaqus as a
  secondary anchor; geometry caveat: FI/ASPECT are flat-Cartesian, curvature
  negligible at 100 km). **De-aliasing matters here** (sharp lateral η step) —
  use `nlat ~ 3·lmax/2`.
- **6c** — load a real 3D viscosity field from netCDF into `earth%visc_3d`
  (node-based); bridge node→element by log10-mean of the bracketing nodes.
- **TRAP-3D** — extend the pseudo-spectral advance to the implicit trapezoidal
  rule (needed before the adaptive controller runs with lateral viscosity).
