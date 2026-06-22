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
| `fe_radial_fe` | per-degree banded radial FE operators | stub |
| `fe_viscoelastic` | Maxwell memory-stress explicit time integration | stub |
| `fe_gravity` | self-gravitation / Poisson coupling | stub |
| `fe_sle` | sea-level equation (ocean function, migration) | stub |
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

**Incompressibility discretization (B4) — OPEN.** The inf-sup-critical choice.
Candidate (unverified): P1 displacement/potential + P0 pressure, `Div u = 0`
enforced weakly via the pressure block. Must be confirmed against Martinec (2000)
— could be P1–P0, a stabilized pair, Taylor–Hood-like, or pressure-elimination.
**Do not implement the FE assembly until this is pinned.** Need the PDF
(`doc/refs/martinec2000.pdf`, AWI OUP access).

**Boundary conditions:** fluid core NOT meshed → CMB boundary condition
(Wu & Peltier 1982): continuity of normal displacement, normal traction,
vanishing tangential traction, potential + flux. Surface: traction balances the
load, potential-flux jump carries the `4πGσ` source. Density-jump interfaces are
natural conditions of the weak form.

**Time scheme (rung 3):** explicit forward-Euler memory-stress update; previous
memory stress relaxed by `(1 − Δt/τ_M)`, `τ_M = η/μ` (Hanyk Eq. 4.36). Stability
`Δt ≲ 2 η_min/μ` ⇒ impose a viscosity floor. VEGA uses Δt = 20 yr.

**Love numbers (rung 2 target):** from surface values for a unit degree-n load,
`h_n = g₀ N_n U(a)/Φ_n`, `l_n = g₀ N_n V(a)/Φ_n`, `k_n = −N_n F(a)/Φ_n`
(Hanyk Eq. 4.40, Farrell normalization). Fluid limits `h_n → −(2n+1)/3`,
`k_n → −1`. Validate vs Spada (2011) Test 2/1.
