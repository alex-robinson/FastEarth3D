# 3-D viscosity performance — profile + design (next-session pickup)

Status after the VILMA-parity rework (commits `506f372..1b0463d`):

- **1-D viscosity: at/beyond VILMA parity.** Default `scheme=fe` (explicit, 1 SLE solve
  per coupling step) + a-priori stability sub-stepping. Full LGM→present Tarasov
  deglaciation: lmax 64 ≈ 41 s (with 3-pass spin-up), lmax 128 ≈ 53 s.
- **3-D viscosity: still slow for a genuinely-3-D field.** The 1-D/3-D layer split
  (`feat(3d)`) only collapses *laterally uniform* layers; a field that varies at all
  mantle depths (Pan et al. 2022) keeps nearly every element on the pseudo-spectral path.

## Profile (Pan 2022, lmax 64, openmp)

Pan et al. 2022 lateral log10(η), M3-L70-V01 radial reference, lmax 64, openmp:

- `ne3d / ne` = **116 / 217** elements classified laterally 3-D. The split collapsed the
  other 101 (elastic lithosphere, fluid core, and laterally-uniform layers) to the 1-D
  path — so the split helps (~46 % of elements skip the transform), but a Pan-type field
  varies across most of the mantle, so 116 still pay.
- `se%update` per step: **1-D = 62 ms**, **3-D = 4196 ms** → ~**68× slower**.
- The whole 3-D-vs-1-D gap is the memory advance (`advance_memory_3d`): the SLE solves and
  band back-subs are unchanged from the 1-D path, so per-step cost ≈ 1-D cost + the
  tensor-SH advance over the 116 3-D elements (~4.1 s of the 4.2 s).
- Compute-bound on the transforms (CPU-saturated across cores). NOTE: a full 3-D run also
  exposes a SEPARATE bottleneck — per-step output writes the large memory-stress arrays —
  which inflated an earlier wall estimate; track that independently of this transform work.

## Root cause — transform strategy, not the layer split

`advance_memory_3d` (`fe_response.f90`) advances the Maxwell memory pointwise on the Gauss
grid, but it issues the spherical-harmonic transforms **per element, per radial shape-
coefficient (A,B,C), per dyadic component**:

```
do e in e3d:                         ! ~ne3d elements
  gather_tensor_coeffs(...)          ! 3 shape-coeffs × 4 strain components
  advance_shape_tensor × 3           ! A, B, C
     -> tensor_sh%synth (memory)     ! ~ a handful of scalar/vector SHTns calls
     -> tensor_sh%synth (strain)     !
     -> pointwise (1-M)τ - 2μM·ε on the grid (6 dyadic components)
     -> tensor_sh%analysis           !
```

So per 3-D advance ≈ `ne3d × 3 × (2 synth + 1 analysis) × ~6 scalar SHTns` ≈ **O(54·ne3d)
scalar transforms**, each O(lmax³). The transforms dominate.

VILMA (`mod_tevolauxsub.f90: taudel2m / harsye2 / harane2m`) instead:

- transforms **once per 3-D layer** for the *whole* stress tensor, with **hand-fused
  multi-field FFTs** (`cdfft13` does 26 arrays in one bit-reversal pass);
- is OpenMP-parallel **over radial layers** `k2`;
- carries the lateral-η update as a single pointwise grid operation per layer.

That fused, per-layer, batched-FFT structure — not the k1p/k2p split alone — is what keeps
VILMA's 3-D at "a few seconds/year at lmax 128".

## Proposed design (to implement next session)

Replace the per-element / per-shape-coeff transform loop with a **per-3-D-layer, batched**
transform:

1. **Batch the components into one multi-field transform.** Synthesize all dyadic
   components (6) × shape-coeffs (3) for a layer in a single batched `SH_to_spat`
   call (SHTns supports many-field transforms), instead of ~18 separate calls. Same for
   analysis. Reuse one plan/config per thread.
2. **Loop and parallelize over 3-D layers** (not elements×coeffs), mirroring VILMA's `k2`
   loop — one synth + pointwise η-update + analysis per layer.
3. **Pointwise η-update stays as is** (`(1-M)τ - 2μM·ε` with the lateral M-field) — it is
   already the cheap part.
4. **Keep the 1-D/3-D split** in front (already merged): only the 3-D layers enter this path.

Expected: collapse ~`54·ne3d` transforms to ~`(synth+analysis)·n_3d_layers` batched
transforms — the same asymptotic work VILMA does, with the constant factor cut by batching
and plan reuse.

### Verification plan
- Reuse `test_response_3d` (uniform field must still reduce to 1-D) + `test_visc_load`
  (Pan 2022 finite/non-trivial) for correctness.
- Re-profile: target 3-D `se%update` within a small factor of 1-D at lmax 64/128.
- Run the full 3-D deglaciation (`runs/deglac_fe64_3d`, `l_visc_3d=.true.`) end-to-end.

### Open questions
- Does the current `fe_tensor_sh` API expose batched/many-field SHTns transforms, or does
  it need extending? (Likely needs a batched synth/analysis entry point.)
- Memory-layout: store the per-layer memory `(component, shape-coeff, lm)` contiguously so a
  batched transform sees a single strided buffer.
