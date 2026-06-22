# Spectral–finite-element formulation (Martinec 2000)

Implementation reference for the radial solver, extracted from Martinec (2000),
*GJI* 142:117 (full PDF held locally, gitignored). Equation numbers below are the
paper's. This is the spec for `fe_radial_fe` assembly and Love-number extraction.

## Problem

Self-gravitating, incompressible, Maxwell-viscoelastic sphere `B` of radius `a`
under a surface mass load `σ`. Governing equations (incremental, quasi-static):

- Momentum (1): `div τ − ρ₀∇φ₁ + div(ρ₀u)∇φ₀ − ∇(ρ₀u·∇φ₀) = 0`
- Poisson (2): `∇²φ₁ + 4πG div(ρ₀u) = 0`
- Constitutive (3,4): `τ = ΠI + 2με`, elastic; Maxwell adds the memory term (rung 3)
- Incompressibility (5): `div u = 0`
- φ₁ = perturbed gravitational potential (interior + the load's exterior potential)

## Weak form (energy functional, §4)

Find `(u, φ₁, Π) ∈ V` with `δE = δF` for all test functions (eq 47).
`E = E_press + E_shear + E_grav + E_uniq` (eqs 30–33); forcing `F = F_surf` (elastic;
the dissipative `F_diss` is rung 3). Variations: eqs 40–46.

- `E_press = ∫ Π div u dV`  (30) — couples pressure ↔ incompressibility
- `E_shear = ∫ μ (ε:ε) dV`  (31)
- `E_grav`  (32) — self-gravity: body force + `(1/8πG)∫|∇φ₁|²`
- `E_uniq`  (33) — removes rigid translation+rotation; **degree j=1 only**
- `F_surf = ∫_∂B (b₀·u + b₁ φ₁) dS` (36), `b₀=−g₀(a)σ eᵣ` (37), `b₁` (38)

## Spectral reduction (§5)

Vector/scalar SH expansion (55–57), `J ≡ j(j+1)`:
- `u = Σ_{j≥1} [U_jm(r) S⁽⁻¹⁾ + V_jm(r) S⁽¹⁾ + W_jm(r) S⁽⁰⁾]` — **no j=0** (incompressibility)
- `φ₁ = Σ_{j≥0} F_jm(r) Y_jm`,  `Π = Σ_{j≥0} Π_jm(r) Y_jm`
- **Spheroidal-only for 1D loading: W=0** (the W block decouples, eq 110).
- `div u = Σ (U' + 2U/r − J V/r) Y_jm`  (58)

Each degree `j` decouples (radially symmetric μ); solve a 1D radial problem per j.

## Discretization (§6)

- P1 "tent" basis ψ_k(r) (71); `U,V,F = Σ_k [·]_k ψ_k` (72) — **nodal**, P+1 nodes.
- P0 basis ξ_k(r) (74); `Π = Σ_k Π_k ξ_k` (73) — **per-element**, P values.
- Per element, μ_k, ρ_k constant (75); gravity `g₀(r) = (4πG/3)(ρ_k r + R_k/r²)`
  (76), with `R_k = Σ_{i≤k}(ρ_{i−1}−ρ_i)r_i³` (77), `R₁=0`.
- Mesh = whole sphere ⟨0,a⟩; fluid core has μ=0. Done in `fe_radial_fe` (218 nodes).

### Assembled bilinear forms (the element matrices to build)

Use the Appendix C element integrals (implemented + tested in
`fe_radial_integrals`): `I1..I7`, `K1..K3`.

- `δE_press` (82): pressure ↔ `(U' + 2U/r − JV/r)` coupling. Uses K1,K2 (per the
  P0×P1 products `∫ξ ψ' r²`, `∫ξ ψ r`, `∫ξ ψ r`-type → K-integrals).
- `δE_shear` (80): U,V (and W, dropped) stiffness. Uses I1,I3,I6 with factor μ_k.
- `δE_grav` (81): self-gravity coupling among U,V,F. Uses I2,I4,I5,I7 (the I7/1/r
  term carries R_k; **skip when R_k=0**, i.e. innermost element) and `1/(4πG)` for
  the `∫|∇φ₁|²` block (I1,I6).
- `δE_uniq` (83): degree-1 rigid-mode removal. Uses K2 (∫ψ r) on degree-1 dofs.

### Surface forcing (eq 84)

`δF_surf = −(a/4πG) Σ (j+1) F^{P+1}_jm δF^{P+1}_jm − a² Σ σ_jm [g₀(a) δU^{P+1}_jm + δF^{P+1}_jm]`

- The `(j+1)F(a)` term = exterior-potential matching (φ₁ ∝ r^{−(j+1)} outside),
  added to the **F–F diagonal at the surface node**.
- The `σ` terms are the **RHS**: forcing on U(a) (`−a² σ g₀(a)`) and on F(a) (`−a² σ`).

## System structure & solve (§10)

Per degree j, DOFs `(U_k, V_k, F_k)` nodal (k=1..P+1) + `Π_k` per element (k=1..P).
Saddle-point (mixed) system, symmetric, **band-diagonal** (P1/P0 local supports):

```
[ A   Bᵀ ] [d ]   [f]          A = shear+grav stiffness (U,V,F)
[ B   0  ] [Π ] = [0]          B = pressure/incompressibility coupling (eq 82)
```

Solve by banded LU (Martinec uses Numerical Recipes BANMUL/BANBKS; we can use a
banded LAPACK `dgbsv`/`zgbsv` or LIS). The matrix is **degree-dependent only
through J**; factor once per j, reuse for all loads/time steps.

### Boundary / regularity conditions
- **Centre r=0:** regularity; for j≥1 displacement → 0 at the centre. The r²
  weighting suppresses node-1 contributions; check whether U₁=V₁=0 must be imposed.
- **Surface:** the (j+1)F(a) exterior term + load RHS above.
- **Uniqueness (j=1):** E_uniq removes the rigid translation null space.
- **Fluid core:** μ=0 region; free-slip emerges (no shear stress). No explicit CMB
  BC (Martinec meshes through the centre).

## Love numbers (§11; normalization to calibrate)

For a unit point-mass load, `σ_jm = (1/a²)√((2j+1)/4π) δ_{m0}` (eq 115). From the
surface coefficients `U_jm(a), V_jm(a), F_jm(a)`, the loading Love numbers follow
(Farrell 1972 normalization):
`h_j ∝ g₀(a) U(a)`, `l_j ∝ g₀(a) V(a)`, `k_j ∝ −F(a)`, scaled by the load's direct
potential `∝ 4πGa/(2j+1)·σ_j`. **Exact constants to be calibrated against the
Spada (2011) Test 2/1 published h,l,k** (and the fluid limits `h_j→−(2j+1)/3`,
`k_j→−1`).

## Validation targets
1. Fluid (t→∞) limits of U,V (paper Fig 6/8 numbers; e.g. for model E,
   `U₂(a)=−278.99e-20 m`). Independent check of the self-gravity/incompressibility.
2. Elastic loading Love numbers h,l,k vs **Spada (2011) Test 2/1**, model M3-L70-V01,
   degrees 2–256.
3. Internal: symmetric matrix; degree-1 uniqueness; mesh refinement convergence.
