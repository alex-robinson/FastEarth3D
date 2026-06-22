# Benchmark reference data — provenance

Reference datasets for validating FastEarth3D against the published GIA
community benchmark. None of this is FastEarth3D output; it is external
reference data, vendored here so the validation tests are self-contained.

## `love_M3-L70-V01/mod_M3-L70-V01`

Normal-mode loading Love numbers for the **M3-L70-V01** earth model
(incompressible, self-gravitating, layered Maxwell; 70 km elastic
lithosphere, three mantle layers, inviscid core), degrees 1–256.

- **Source:** the `giapy` package (S. B. Kachuck), file
  `giapy/data/earth/mod_M3-L70-V01`. https://github.com/skachuck/giapy
- **License:** MIT (giapy) — redistribution with attribution is permitted.
- **Underlying benchmark:** the Charles University GIA Benchmark
  (https://geofjv.troja.mff.cuni.cz/GIABenchmark), i.e. the model defined in
  Spada et al. (2011), *A benchmark study for glacial isostatic adjustment
  codes*, Geophys. J. Int. 185, 106–132.
- **How it was generated (not by us):** these are *normal-mode*
  (analytic / semi-analytic) viscoelastic Love numbers — for an incompressible
  layered Maxwell sphere the response is a finite sum of decaying exponentials
  (here 9 modes per degree). They are the output of a normal-mode code such as
  TABOO or ALMA (Spada's codes), not a time-domain FE code. FastEarth3D
  computes the *same physics* by a different route (per-degree FE solve for the
  elastic/fluid limits, time-domain memory-stress integration for the transient),
  so this table is an independent cross-check, and TABOO
  (https://github.com/danielemelini/TABOO, GPLv3) can regenerate it externally
  if needed.

### File format (per the giapy loader)

- Lines 1–5: the earth model — `layer  r[m]  rho[kg/m^3]  mu[Pa]  eta[Pa·s]`
  (η = 1e44 ≈ ∞ marks the elastic lithosphere; μ = η = 0 marks the inviscid core).
- Then, for each degree `n`:
  - header `n  nmodes  k_e  h_e  l_e`  — **elastic** (instantaneous) Love numbers,
  - `nmodes` lines of normal-mode residues + relaxation,
  - one line `... h_f l_f k_f` — the **fluid** (t→∞) Love numbers.
  - (column order in the elastic/fluid lines is k, h, l after the leading index.)

## `sle_martinec2018/*_SBK.dat`

Reference spatial response curves for the sea-level-equation benchmark,
cases A, C2, D3, E2, F1 (figs 10–13), at truncation j256.

- **Source:** `giapy`, `tests/sbk_benchmark_data.tar.gz` (MIT).
- **Underlying benchmark:** the same Charles University GIA Benchmark; the SLE
  intercomparison published as Martinec et al. (2018), *A benchmark study of
  numerical implementations of the sea level equation in GIA modelling*,
  Geophys. J. Int. 215, 389–414.
- **Columns:** `col1` colatitude/longitude [deg]; `col2` vertical displacement
  [m]; `col3`,`col4` θ- and φ-horizontal displacement [m]; `col5` gravitational
  potential increment [SI]; `col6` sea-surface variation w.r.t. h_UF [m];
  `col7` sea-level-equation result [m].
- **Benchmark spec** (from giapy `tests/sle_test.py`): ice cases L1/L2/L3,
  topographies B0–B3 (exponential basins), time histories T1/T2/T3.
