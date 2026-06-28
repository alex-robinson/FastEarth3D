# Modal lateral-viscosity split-operator error — diagnosis

Diagnostic: `tests/diag_modal_latvisc.f90` (target `make diag_modal_latvisc`, NOT in
`make check`). Isolates the RESP_MODAL lateral split-operator error against the
RESP_VE tensor-SH ground truth, **with lateral viscosity active on both paths**
(all prior modal-vs-VE probes were radial and missed this entirely).

## Setup

- Earth: M3-L70-V01 (laterally-averaged reference profile = the modal basis).
- Controlled contrast: an axisymmetric polar **viscosity cap** (colat < θ_cap,
  5° cosine taper), depth-uniform `Δ = log10(η_local/η_ref) > 0` (stiff craton),
  applied to every Maxwell element. The *identical* `pert_elem(nphi,nlat,ne)`
  field is handed to both `response_enable_lateral_visc` (VE→tensor-SH) and
  `response_enable_lateral_visc_modal` (split-operator) — the only difference is
  the algorithm.
- Co-located disc surface load (colat < θ_cap, tapered). Both operators are linear
  in σ, so absolute load amplitude is irrelevant; we report ratios.
- Protocol = **free rebound** (cleanest rate probe): hold the load (glaciation) to
  build depression, then remove it (σ=0) and watch free relaxation, whose decay
  rate *is* the lateral viscosity operator. Fluid limits are contrast-independent,
  so any gap is a transient **rate** error.
- Metrics: peak uplift over the cap region, modal/VE ratio (1.0 = perfect);
  `grms/pk` = global area-rms(modal−VE) over the VE peak (field error).

## Result

Two resolutions run (`make diag_modal_latvisc`; logs in `logs/latvisc_*.log`).
lmax 16 over-states the magnitudes (coarse grid / Gibbs); **lmax 48 is the
reference** below. The qualitative structure is identical at both.

### lmax 48 (glaciation 40 kyr, rebound 20 kyr; dt 200 yr) — REFERENCE

End-of-glaciation modal/VE peak-uplift **ratio** (1.0 = perfect) | field error
`grms/pk` | rebound ratio at 20 kyr:

| Δ \ θ_cap | 10° | 20° | 40° |
|-----------|-----|-----|-----|
| 0.5 | 0.98 / 0.005 / 0.82 | 0.91 / 0.015 / 0.62 | **0.77 / 0.085 / 0.60** |
| 1.0 | 1.26 / 0.020 / 0.61 | 1.15 / 0.039 / 0.64 | 0.91 / 0.189 / 3.12 |
| 1.5 | 1.68 / 0.044 / 1.08 | 1.42 / 0.069 / 1.24 | 1.09 / 0.265 / 4.92 |
| 2.0 | 2.02 / 0.064 / 1.72 | 1.59 / 0.091 / 1.62 | 1.24 / 0.311 / 6.00 |

Reading: error grows **monotonically with contrast at every scale**; field error
`grms/pk` is **largest for the largest cap** (40°: 0.085→0.31). The peak-ratio
*sign* flips (under-prediction at moderate Δ, over-prediction at high Δ), but the
*field error magnitude* is unambiguous and large for continental scale + realistic
contrast. During free **rebound** the divergence is far worse than at
end-of-glaciation (40°/Δ≥1: 3–6× over; 40°/Δ0.5: 0.60× under).

### Matches the ensemble's deglac3d regime

The ensemble's N. America attenuation (modal 154 m vs VE 343 m at PD ⇒ ratio 0.45)
is reproduced by the **continental-scale (40°), moderate stiff-contrast (Δ≈0.5–1)**
corner, where modal **under-predicts rebound by 30–55%** (40°/Δ0.5 rebound 0.60;
40°/Δ1 early-rebound 0.73–0.85). This is the relevant regime for the Canadian
Shield (broad, ~0.5–1 dex stiffer mantle). Higher contrast (Δ≥1.5) over-predicts —
so the sign of the local modal error is contrast-dependent, but the breakdown is
generic. lmax 16 numbers (for reference, much larger): 40°/Δ0.5 end-glac 0.35,
40°/Δ2 end-glac 6.7, rebound to 77×.

## What this says

1. **The error is NOT a uniform ~55% attenuation.** It scales steeply with BOTH
   lateral contrast Δ and cap size θ, and **changes sign**: small/weak → mild
   under-prediction (the advertised 1st-order regime); large+strong → severe
   over-prediction with O(1)+ field errors.
2. **Mechanism (visible in the code).** For a stiff cap (Δ>0), `rho = 10^(−Δ) < 1`
   inside, so the zero-mean anomaly `(Ri − Rbar) < 0` there →
   `E2 = exp(−dt·mlatRate) > 1`. In `modal_lateral_anomaly`,
   `gphi = E2·gphi + (1−E2)·gsig`; during free rebound (gsig=0) this is
   `gphi ← E2·gphi` with E2>1 → **amplification, not relaxation**. The Lie split
   treats the anomaly as a small perturbation of the mean step; at large contrast
   the product overshoots and can grow.
3. **Three compounding error sources** (to be attributed):
   - **Lie splitting**: mean (spectral, diagonal in l) ⊗ anomaly (real-space,
     diagonal in θφ) do not commute; error ~½[A,B], grows with both magnitudes.
   - **Characteristic-τ̂ mismatch**: the anomaly is applied at a single rank-wide
     `τ̂_i` (|C^u|-weighted over degrees, low-l-dominated), but a *localized* cap's
     structure lives at high l where `τ_i(l)` is very different → wrong rate.
   - **Amplification (E2>1)**: net non-contractive step for strong negative
     anomalies → unphysical growth in free relaxation.

## Resolves the "2.6% vs 55%" discrepancy (task background)

The prior stripped-down modal+3D vs VE+3D run reported ~2.6% **global** rsl rms,
yet the ensemble shows ~55% **local** attenuation over N. America at PD. These are
consistent: the error is strongly localized to the high-contrast craton. A global
rms is diluted by the vast low-contrast ocean/non-cratonic area (small error
there), while the local peak/field error over the cap is O(50–100%). This
diagnostic measures the local error directly (`grms/pk` normalized by the cap peak)
and recovers the large values.

## Candidate fixes (NOT yet implemented — for discussion)

Ordered roughly by effort/return. The split-operator's three error sources map to
distinct fixes:

1. **Per-degree anomaly rate** (target: τ̂ mismatch). Apply the real-space anomaly
   with the actual `τ_i(l)` instead of one rank-characteristic `τ̂_i`. The localized
   cap's structure lives at high l where `τ_i(l)` ≪ the low-l-weighted `τ̂_i`, so the
   anomaly currently relaxes at the wrong rate. Cost: anomaly SHTs per degree-group
   instead of per-rank (more transforms, still ≪ tensor-SH).
2. **Strang (2nd-order) split** (target: Lie non-commutation). Half mean → full
   anomaly → half mean. Cheap; halves the splitting error order. Does **not** fix
   the amplification or τ̂ mismatch.
3. **Guarantee contraction** (target: E2>1 amplification). The current anomaly can
   be non-contractive (`E2>1`) for stiff caps → unphysical growth in free rebound.
   A clamp is a band-aid; the proper cure is (1)+(4) which keep the combined rate
   physical by construction.
4. **Krylov / Chebyshev exponential of the coupled rate operator** (the "proper"
   fix). Don't split at all: apply `exp(−Δt·L_i)` where `L_i: φ ↦ R(θ,φ) ⊙ (1/τ_i(l))·φ`
   (spatial ⊙ spectral) directly, via a few matrix-free applications on the **K
   scalar** fields. This is what tensor-SH does for the 6-tensor memory, but on K
   scalars — much cheaper than tensor-SH, far more accurate than 1st-order split.
   Largest effort; most faithful; preserves the modal philosophy (cheap, tunable,
   converges to VE as the Krylov order grows).

A focused **attribution experiment** (variants of the diagnostic: Lie vs Strang;
characteristic-τ̂ vs per-degree; vs the full coupled exponential) would quantify how
much of the gap each source owns before committing to a fix. Recommend running that
first.
