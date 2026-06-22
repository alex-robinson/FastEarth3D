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

Per degree j, DOFs `(U_k, V_k, F_k)` nodal (k=1..P+1) + `Π_k` per element (k=1..P),
laid out **node-interleaved** `[U V F | Π]` so the operator stays band-diagonal:

```
[ A   Bᵀ ] [d ]   [f]          A = shear+grav stiffness (U,V,F)
[ B   0  ] [Π ] = [0]          B = pressure/incompressibility coupling (eq 82)
```

The shear block (eq 80) is **symmetric** (it is `∫μ ε:ε`); the **B/Bᵀ** pressure
block is symmetric by construction; but the **self-gravity `I²` U↔F coupling
(eq 81) is not symmetric**, so the assembled operator is non-symmetric overall —
which is why Martinec solves it with a *general* banded LU (BANMUL/BANBKS), not a
symmetric one. Implemented in `fe_radial_fe%build_dense_operator`, transcribed
term-by-term from the PDF and verified against the analytic limits below.

**Solve (no LAPACK): LIS** (`fe_lis`). The physical entries span ~20 orders of
magnitude (`μr²/h` vs the pressure couplings vs `1/4πG`), so the operator is
geometric-mean **row/column equilibrated** before the solve; restarted GMRES with
an ILU(1) preconditioner then converges effectively as a direct solve (1 GMRES
iteration, residual ~1e-13) on this band-diagonal system. The matrix is
**degree-dependent only through J**; assemble + equilibrate once per j, reuse for
all orders m / loads / time steps (precon-reuse optimisation still TODO).

### Boundary / regularity conditions
- **Centre r=0:** **no explicit BC** — Martinec meshes through the centre and the
  r² weighting handles regularity; the singular `I⁷` term is killed by `R₁=0`
  (skipped to avoid 0·∞). Verified: no empty rows, solver well-posed for j≥2.
- **Surface:** the (j+1)F(a) exterior term + load RHS above.
- **Uniqueness (j=1):** E_uniq removes the rigid translation null space.
- **Fluid core:** μ=0 region; free-slip emerges (no shear stress). No explicit CMB
  BC (Martinec meshes through the centre).

## Love numbers (§11; conventions verified)

For a degree-j surface load of coefficient `σ`, the load's own potential at the
surface is `φ^L = 4πG a σ/(2j+1)`. From the surface coefficients (Farrell 1972
normalization), implemented in `fe_radial_fe%loading_love`:

`h_j = g₀(a) U(a)/φ^L`,  `l_j = g₀(a) V(a)/φ^L`,  `k_j = −F(a)/φ^L − 1`.

The `k` form was **pinned empirically by two analytic limits**: Martinec's `φ₁`
(=`F`) is the *total* perturbation potential and carries the load's direct
potential with the **opposite sign** to `φ^L` (`F→−φ^L` for a rigid sphere), so
the induced potential is `−F−φ^L`. `σ` cancels in every ratio (use σ=1).
**`l` still needs its sign / S⁽¹⁾-normalization factor calibrated** against the
published Spada `l` (h and k are fully pinned).

## Validation targets
1. **Fluid limit** (μ→0, homogeneous sphere): `h_j → −(2j+1)/3` and `k_j → −1`.
   ✅ reproduced to ~1e-5 (`test_love`). Independent check of self-gravity +
   incompressibility + Poisson + the load forcing.
2. **Rigid limit** (μ→∞): `h_j, l_j, k_j → 0`. ✅ to ~1e-5 (`test_love`).
   Checks the shear block and the F sign convention.
3. Elastic loading Love numbers h,l,k vs **Spada (2011) Test 2/1**, model
   M3-L70-V01, degrees 2–256. Current output is physical (k₂≈−0.37, decaying with
   j); a quantitative match awaits the published table (not in-repo). 🔶
4. Internal: operator finite (centre I⁷ guard); B=Bᵀ; degree-1 uniqueness;
   gravity R_k reconstruction (`test_assembly`).
