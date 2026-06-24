# fastearth restart / output variables

Variable-io table (yelmo convention) for FastEarth3D netCDF I/O, used by
`fe_io` for both restart files and diagnostic `write_step` output. Each row gives
the netCDF variable name, its dimensions, units, and a long_name. The time axis
(unlimited) lets several snapshots live in one file.

The Maxwell memory-stress fields `tau_*` and the adaptive controller's next-step
Δt seed `dt_try` are the prognostic state restored on restart; the reference
fields `z_bed_eq`/`h_ice_ref` are static (written once) and checked on read; the
rest are diagnostic.

| id | variable     | dimensions    | units  | long_name                                        |
|----|--------------|---------------|--------|--------------------------------------------------|
|  1 | tau_a_re     | nlam, ne, nk  | Pa     | Maxwell memory stress, component A (real part)   |
|  2 | tau_a_im     | nlam, ne, nk  | Pa     | Maxwell memory stress, component A (imag part)   |
|  3 | tau_b_re     | nlam, ne, nk  | Pa     | Maxwell memory stress, component B (real part)   |
|  4 | tau_b_im     | nlam, ne, nk  | Pa     | Maxwell memory stress, component B (imag part)   |
|  5 | tau_c_re     | nlam, ne, nk  | Pa     | Maxwell memory stress, component C (real part)   |
|  6 | tau_c_im     | nlam, ne, nk  | Pa     | Maxwell memory stress, component C (imag part)   |
|  7 | z_bed_eq     | lon, lat      | m      | Reference (equilibrium) bedrock elevation        |
|  8 | h_ice_ref    | lon, lat      | m      | Reference grounded-ice thickness                 |
|  9 | h_ice        | lon, lat      | m      | Grounded-ice thickness                           |
| 10 | rsl          | lon, lat      | m      | Relative sea-level change (full field)           |
| 11 | z_bed        | lon, lat      | m      | Bedrock elevation (z_bed_eq - rsl)               |
| 12 | C_ocean      | lon, lat      | 1      | Ocean function (1 ocean / 0 land)                |
| 13 | dt_try       | time          | s      | Adaptive time-stepping next-step Δt suggestion   |
| 14 | sigma_n_re   | nlm, time     | kg m-2 | Trapezoidal start-of-step load σ_n (real part)   |
| 15 | sigma_n_im   | nlm, time     | kg m-2 | Trapezoidal start-of-step load σ_n (imag part)   |
| 16 | sigma_primed | time          | 1      | Flag: σ_n is tracked (1) or not yet primed (0)   |
