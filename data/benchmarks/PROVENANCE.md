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

## `disc_spada2011/{u,n,dudt}_disc.txt`

Spatial response of the **disc load** test of the Spada et al. (2011) benchmark:
vertical displacement `u`, geoid `n`, and uplift rate `dudt`, for the
**M3-L70-V01** earth model. Load: a 10°-radius disc of 1000 m ice
(ρ_ice = 931 kg m⁻³) on elevated bedrock (no ocean load).

- **Source:** the `isostasy_data` repository (J. Jereczek),
  `model_outputs/Spada-2011/`. https://github.com/JanJereczek/isostasy_data
- **Underlying benchmark:** Spada et al. (2011), *A benchmark study for glacial
  isostatic adjustment codes*, Geophys. J. Int. 185, 106–132 (the Charles
  University GIA Benchmark). Same earth model as the Love table above.
- **Grid:** 201 rows, colatitude θ = 0:0.1:20°; 6 columns, times
  t = [0, 1, 2, 5, 10, 100] kyr. **Column 1 (t=0) is the elastic response.**
- **GOTCHA — the `n_*` (geoid) files are stored in REVERSED θ order.** The
  upstream loader (`isostasy_data`/FastIsostasy.jl `dataloaders.jl`) applies
  `reverse!(X, dims=1)` to the `n_` cases only; `u_*` and `dudt_*` are NOT
  reversed. So `n_disc.txt` row 1 is θ=20° and row 201 is θ=0°; reverse the rows
  before pairing with the θ grid. (See `test_benchmark_disc`, which reverses `n`.)
- **Degree-1 / geoid frame:** the displacement is in the CE-like gauge (h₁≈0,
  geocenter), the geoid in the CM frame (N₁=0). FastEarth3D reproduces both: u to
  ~1% near-field, n to ~1% once the degree-1 geoid is referenced to CM (N₁=0; see
  fe_response). The far-field (θ≳12°) forebulge in `u` is small-amplitude and
  shows larger relative differences (low-degree-truncation sensitive).

## `sle_martinec2018/*_SBK.dat`

Reference spatial response curves for the sea-level-equation benchmark,
cases A, C2, D3, E2, F1 (figs 10–13), at truncation j256.

- **Source:** `giapy`, `tests/sbk_benchmark_data.tar.gz` (MIT).
- **Underlying benchmark:** the same Charles University GIA Benchmark; the SLE
  intercomparison published as Martinec et al. (2018), *A benchmark study of
  numerical implementations of the sea level equation in GIA modelling*,
  Geophys. J. Int. 215, 389–414.
- **Two file formats** (the case-A loading file differs from the SLE-case files):
  - **`A_fig10_SBK.dat`** (case A only): 4 columns, header
    `# Colat(deg) Uplift(m) Horizontal(m) Geoid(m)`; rows are colatitude
    180°→0° at 0.25° spacing (721 rows). This case has NO ocean (B0=0): a pure
    viscoelastic cap-loading response, with degrees n=0,1 dropped (giapy jmin=2).
  - **`{C2,D3,E2,F1}_fig{10..13}_SBK.dat`** (the SLE cases): 7 data columns —
    `col1` longitude/colatitude [deg]; `col2` vertical displacement [m];
    `col3`,`col4` θ- and φ-horizontal displacement [m]; `col5` gravitational
    potential increment [SI, = −geoid·9.815]; `col6` sea-surface variation w.r.t.
    h_UF [m]; `col7` sea-level-equation result [m]. Profiles are along circles of
    constant lon or lat (figs 10–13 use lon=75, lat=25, lon=±/25, lat=35/100).
- **Benchmark spec** (from giapy `tests/sle_test.py`): ice cases L1/L2/L3
  (spherical caps h(δ)=h0·√[(cosδ−cosα)/(1−cosα)], α=10°, at given centre),
  topographies B0–B3 (B = bmax − b0·exp(−δ²/2σ²) exponential basins, σ=26°),
  time histories T1 (Heaviside at 10 kyr) / T2 (10-kyr linear growth) / T3
  (linear decay). giapy ice/water/sea densities 931/1000/1000 kg m⁻³, g=9.815.
  Case A = L1+T1+B0; the SLE cases combine L2/L3 + T2 + B1/B2/B3.
