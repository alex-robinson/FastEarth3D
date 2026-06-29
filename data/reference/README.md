# Reference state (present-day topography + ice) on the Gauss grid

Prebaked present-day reference fields for the driver's `i_eq=1` (and `=3`) modes.
With `i_eq=1` the relaxed reference `z_bed_eq` (= the SLE `topo0`) and `h_ice_ref`
are this present-day state, so relative sea level is measured against today and
`rsl ≈ 0` at the present day.

`read_ref2d` (fe_drive) auto-detects that these files are already on the Gauss
grid (their `lon`/`lat` dims equal the run's `nphi`/`nlat`) and reads them
**directly**, skipping the expensive online conservative-map build (~160 s for the
0.125° RTopo, done twice per run). Point a run at the file matching its `lmax`:

| file                   | lmax | Gauss (nphi×nlat) |
|------------------------|------|-------------------|
| `rtopo_gauss_l32.nc`   | 32   | 128×66            |
| `rtopo_gauss_l64.nc`   | 64   | 256×130           |
| `rtopo_gauss_l96.nc`   | 96   | 384×194           |
| `rtopo_gauss_l128.nc`  | 128  | 512×258           |

## Provenance

Source: `RTopo-2.0.1_0.125deg_DRThydrocorr.nc` (CLIMBER-X input; `bedrock_topography`,
`ice_thickness`). Pipeline (offline, one-time):

1. `rtopo_0.5deg.nc` — the 0.125° RTopo coarsened to 0.5° by an exact 4×4 block
   mean (conservative on a regular grid; the 0.5° `topo_05x05.nc` in CLIMBER-X
   carries bed only, so the full RTopo is coarsened here to keep ice too).
2. `rtopo_gauss_l{32,64,96,128}.nc` — `bin/fastearth_mkref.x prebake_l<L>.nml fastearth.nml`
   conservatively remaps the 0.5° bed (as-is) and ice (mass-conserving) onto each
   Gauss grid with the same `fe_remap` engine the online path uses.

Any integer lmax works (the grid is nphi=4·lmax, nlat=2·lmax+2) — each level is
remapped independently from the 0.5° source, not aggregated from a finer Gauss grid,
so the resolutions need not be powers of two. Regenerate any level with
`make fastearth_mkref` then running its `prebake_l<L>.nml`. A finer source is
unnecessary up to ~lmax-180: the 0.5° topo (~l180) and the lateral viscosity
(bagge2021, 512×256 ≈ l128) both run out of real structure before then.
