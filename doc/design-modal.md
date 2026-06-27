# FastEarth3D — modal viscoelastic response (design)

This document specifies a new viscoelastic response representation, `RESP_MODAL`,
that reduces the per-degree relaxation to a chosen number of normal modes `K`.
It is the design reference for the implementation on branch `modal`.

It complements `design.md`; the full Martinec time-domain FE method is the
*reference* this approximation is validated against and, with `K = all`, the
limit it reproduces.

## 1. Motivation

FastEarth3D today has two endpoints with nothing between them:

- **Full VE (`RESP_VE`, Martinec FE)** — accurate, but pays a banded radial
  solve every step, carries a heavy radial memory-stress state
  (`Are/Aim/Bre/Bim/Cre/Cim(NLAM, ne, nk)`), is bounded by a Maxwell `Δt`
  stability ceiling (hence sub-stepping), and routes lateral viscosity through
  the expensive 6-component tensor-SH path (`fe_tensor_sh`).
- A hypothetical **single-mode** model (LV-ELVA, as in FastIsostasy) — cheap and
  unconditionally stable, but keeps only one relaxation mode: lossy, with no
  knob to recover accuracy.

The real Earth's loading response, per harmonic degree, is a **sum of a handful
of exponential relaxation modes**. The natural design is therefore neither
endpoint but a **dial**: choose how many modes `K` to keep per degree.

- `K = 1`   → the dominant mode ≈ LV-ELVA;
- `K = all` → reproduces the full FE model (exactly, for radial viscosity);
- `1 < K < all` → tune fidelity to the modelling context — fast continental
  deglaciation at low `K`, careful regional GIA at high `K` — in one code path,
  one validation harness, one coupling interface.

This is the **LV-ELVA representation translated into our global
spherical-harmonic framework**. We do *not* import FastIsostasy's regional
Cartesian / FFT machinery; only its physical idea (a depth-collapsed,
mode-limited relaxation with stacked lateral viscosity).

### Relation to the "no normal-mode root finding" constraint

`design.md §1` rules out normal-mode root finding — meaning the analytic
secular-determinant approach of SELEN/TABOO/ALMA. We do **not** do that. We
**numerically diagonalize the already-assembled FE radial operator** per degree.
This keeps the existing FE infrastructure as the single source of truth: the
modal basis is a property of the same discretized operator the full model uses,
so `K = all` matches `RESP_VE` by construction. The eigensolve is a one-time
init cost, not a per-step root find.

## 2. Physical basis — how many modes

For an **incompressible**, self-gravitating, layered Maxwell Earth the
relaxation spectrum per degree is **finite and discrete** (this is exactly why
incompressibility is assumed — it keeps the count finite; compressibility would
make it a continuum). Roughly one mode per interface:

- **M0** — the fundamental mantle mode (dominant for surface loads);
- one **buoyancy (M-type)** mode per internal density discontinuity;
- one **transient (T)** mode per viscosity contrast / Maxwell layer;
- **C0** — the core–mantle-boundary mode.

For an M3-L70-V01-class structure (lithosphere + 3 mantle Maxwell layers + fluid
core) that is **~10 physical modes, of which typically 1–3 dominate** at any
given degree (which ones is degree-dependent). The numerical eigensolve of the
FE operator yields many more eigenvalues, but the surplus are numerical with
negligible residue. Hence:

> **`n_modes = K` keeps the `K` modes of largest residue (modal strength).**

Physical content saturates well before `K` reaches double digits.

## 3. Core idea

Per degree `l`, the FE-discretized viscoelastic system is linear time-invariant.
Its load→surface response, after diagonalization, is

```
u_lm(t) = u_el(l)·σ_lm  +  Σ_{k=1..K} r^u_k(l) · ξ_k,lm(t)
N_lm(t) = N_el(l)·σ_lm  +  Σ_{k=1..K} r^N_k(l) · ξ_k,lm(t)
V_lm(t) = V_el(l)·σ_lm  +  Σ_{k=1..K} r^V_k(l) · ξ_k,lm(t)

ξ̇_k,lm = −ξ_k,lm / τ_k(l)  +  σ_lm        (one scalar ODE per mode per coeff)
```

- `u_el, N_el, V_el` — instantaneous **elastic** surface gains (the `t→0⁺`
  limit; same quantities `RESP_ELASTIC`/`RESP_VE` already compute).
- `{τ_k(l), r^{u,N,V}_k(l)}` — relaxation times and **modal residues** from the
  per-degree eigensolve (§5).
- The **isostatic (fluid) limit** is reproduced automatically:
  `u_∞ = u_el + Σ_k r^u_k τ_k` etc. — a check, not a separately-imposed datum.

**State shrinks** from the full radial memory tensor to **`K` complex scalars per
(l,m)** — `ξ_k,lm`. Each obeys a scalar linear ODE.

**Exact time integration.** Over a step with constant load `σ`, each mode
advances by its analytic solution:

```
ξ_k ← e^{−Δt/τ_k} · ξ_k  +  τ_k (1 − e^{−Δt/τ_k}) · σ
```

so the scheme is **unconditionally stable — no Δt ceiling, no sub-stepping**.
For the SLE-coupled / fast-load case the same closed form with the trapezoidal
(linear-in-time σ) variant is used, giving 2nd order in the load without a
stability constraint. (This is *not* the rejected ETD-on-FE-memory route: that
failed because `exp(−M)` per element in the coupled FE basis does not decouple.
Diagonalizing **first** is precisely the decoupling that makes the exponential
exact and cheap.)

## 4. Lateral viscosity (the "LV")

The model is built to run with laterally-varying viscosity. Lateral η breaks the
clean per-degree diagonalization (it couples harmonics), so it is handled the way
FastIsostasy handles its single-mode case, generalized to `K` modes:

1. **Modal basis from a reference profile.** Diagonalize once using the
   **laterally-averaged** radial viscosity profile → `{τ_k(l), residues}`. This
   is the one-time eigensolve; lateral variation does **not** trigger
   re-diagonalization.
2. **Lateral variation as pointwise rate modulation (split-operator).** Each
   mode's relaxation rate is modulated by the local viscosity on the Gauss grid:
   multiply by the local rate factor in **real space**, apply the degree/mode
   kernel in **spectral space**, transform back with `fe_sht`. This reuses the
   existing 3D-path pattern but acts on `K` **scalar** fields, not a 6-component
   memory tensor — much cheaper.
3. **Depth-weighted per-mode modulation.** Each mode `k` is modulated by the
   lateral viscosity in the **depth band that mode samples**, weighted by the
   mode's radial strain-energy profile (available from its eigenvector). M0 is
   driven by deep-mantle structure, transient modes by their own layers — a
   fidelity single-mode LV-ELVA cannot represent, and a direct benefit of the
   multi-mode form.

**Exactness ladder (be honest about it):**

- *Radial* η → `K = all` equals the full FE model exactly.
- *Lateral* η → modal is an **approximation at every `K`** (true lateral coupling
  is richer than rate modulation). The full tensor-SH path (`RESP_VE`) remains
  the ground-truth reference for validation; modal-LV is the fast, tunable
  approximation to it, cheaper than tensor-SH at all `K`.

## 5. Eigensolve (the one-time cost)

**Decision: numerically diagonalize the FE operator** (not analytic
root-finding). Per degree `l`:

- The semi-discrete viscoelastic system is a generalized linear relaxation
  problem in the radial DOFs. Form the relaxation pencil from the same assembled
  per-degree operator the full model uses (`radial_operator`, `fe_radial_fe`) and
  the Maxwell rate structure (`μ/η` per element).
- Solve the generalized eigenproblem (LAPACK `dggev` / `dgeev` as appropriate)
  → eigenvalues `s_k = −1/τ_k` and eigenvectors. Incompressibility/pressure makes
  it a descriptor (DAE) system: the algebraic constraints appear as
  infinite/spurious eigenvalues, filtered out by residue.
- For each retained mode compute its surface residues `r^{u,N,V}_k(l)` (from the
  eigenvector's surface DOFs and the load projection) and its radial
  strain-energy weight (for §4.3 depth-weighting).
- Sort by `|residue|`, keep `n_modes` (or all). Store `τ_k(l)`, residues, weights.

**Cost.** Per degree, a dense generalized eigensolve is `O(N_r³)` with `N_r` ≈ a
few hundred radial DOFs; embarrassingly parallel over the `lmax+1` degrees (the
degree loop is already OpenMP-parallel). Estimate: **~1–10 s wall at lmax128**
for a few-hundred-DOF mesh, up to **~tens of seconds** for a fine radial mesh —
**once, at init**, comparable to the existing `se_init` auto-tune. If only `K`
modes are wanted, a partial/Krylov eigensolver drops this to `O(K·N_r²)`.

**Runtime per step** is then negligible: `K` exact-exponential scalar updates per
(l,m), versus the current per-step banded radial solve. Lateral adds `K` scalar
forward/inverse SHTs per step — cheaper than the 6-component tensor-SH advance.
Net: **faster runtime than full VE, modest one-time init.**

## 6. Architecture / integration

A new response kind `RESP_MODAL` in `fe_response.f90`, beside `RESP_VE`.
Everything downstream is untouched because it goes through the abstract
`response` interface (`response_apply`, `response_horizontal`,
`response_begin/commit_step`, the endpoint brackets, `response_set_dt`,
save/restore, restart).

- **Init** (`response_init_modal`): assemble per-degree operators (existing),
  eigensolve (§5), select modes, store `{τ_k(l), r^{u,N,V}_k(l), weight}`.
  Reference profile = laterally-averaged η.
- **State**: `ξ_k,lm` — `K` complex scalars per deforming (l,m), degree-grouped
  like the existing `k` layout. Replaces the memory-stress arrays.
- **apply / horizontal**: affine — elastic gain · σ + Σ_k residue · ξ_k.
- **commit_step / endpoint brackets**: advance `ξ_k` by the exact exponential
  (constant-σ) or its trapezoidal (linear-σ) form; SLE co-convergence uses the
  same closed form re-evaluated against the converged σ.
- **set_dt**: only rescales `e^{−Δt/τ_k}` — no operator re-factor (as today).
- **Lateral** (`response_enable_lateral_visc*`): scalar split-operator via
  `fe_sht` with depth-weighting (§4); reuses the `visc3d_tol`-style activity
  bookkeeping where useful.
- **Unchanged**: SLE (`fe_sle`), geoid (the per-degree modal sum feeds `gn`),
  rotation/TPW (`fe_rotation`), coupling API (`fe_coupling`), I/O/restart
  (`fe_io`) — restart persists the `ξ_k` fields + the modal basis metadata.
- **Namelist** (`fe_params`, `&fe3d`): `n_modes` (integer; `-1`/`all` = full),
  `mode_residue_tol` (drop modes below this fraction of the largest residue),
  and a selector for the modal scheme.

## 7. End goal

`RESP_MODAL` with `K = all` is intended to **subsume `RESP_VE`**: the full model
becomes the all-modes case of the modal model. `RESP_VE` is retained as the
validation reference until the equivalence is demonstrated bit-close on the
benchmarks, then it can be retired (or kept as the lateral-η ground truth).

## 8. Validation

1. **`K = all`, radial η** reproduces the existing FE results bit-close
   (M3-L70-V01 Love numbers; 1-D deglaciation; the §3c benchmarks).
2. **Mode-count sweep** of `K` against full tensor-SH lateral runs (M3-L70-V01,
   deglaciation) → an accuracy/cost curve quantifying the approximation per `K`.
   This curve is the deliverable that justifies "hone to context".
3. **Limit checks**: elastic (`t→0⁺`) and isostatic (`t→∞`) limits from the modal
   sum match the directly-computed elastic and fluid Love numbers.

## 9. Implementation plan (incremental, small commits)

1. **Eigensolve + radial-only modal response.** `RESP_MODAL`, `response_init_modal`,
   the per-degree generalized eigensolve, mode selection, modal state `ξ_k`,
   exact-exponential `apply`/`commit`. **Validate `K = all` = `RESP_VE`** on the
   1-D benchmarks before anything else. Expose `n_modes`.
2. **Lateral η via depth-weighted split-operator rate modulation.**
3. **SLE / endpoint co-convergence + restart** wiring for the modal scheme.
4. **Accuracy/cost sweep + docs**; begin subsuming `RESP_VE`.
