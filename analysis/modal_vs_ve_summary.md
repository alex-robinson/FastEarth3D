# Modal vs VE — accuracy & cost summary

Each modal run is compared against the VE run in the same set. `rsl RMSE` is area-weighted (cos lat) over all space-time; `rsl PD` at present day; `cost` is the per-coupling-step solver time; `speedup` = VE/modal.

### radial

VE reference: solver cost = 45.6 ms/step, n_solve = 1.6 sub-steps/step

VE scales: |rsl|max(PD) = 146.0 m, BSL min = -118.9 m, BSL(PD) = 1.51 m

| modal | rsl RMSE | rsl PD | rsl max | bsl RMSE | bsl PD | cost | speedup | nsolve |
|---|---|---|---|---|---|---|---|---|
| all | 1.183 m | 0.885 m | 9.653 m | 0.010 m | 0.004 m | 8.0 ms | 5.7x | 1.0 |
| n1/isostatic | 17.420 m | 6.026 m | 196.720 m | 0.101 m | 0.000 m | 7.4 ms | 6.2x | 1.0 |
| n1/rate | 19.362 m | 6.191 m | 214.367 m | 0.154 m | 0.002 m | 7.3 ms | 6.2x | 1.0 |
| n1/residue | 23.944 m | 7.198 m | 273.533 m | 0.126 m | 0.004 m | 7.3 ms | 6.2x | 1.0 |
| n2/isostatic | 5.776 m | 3.786 m | 59.044 m | 0.038 m | 0.002 m | 7.4 ms | 6.2x | 1.0 |
| n2/rate | 5.455 m | 4.912 m | 50.402 m | 0.037 m | 0.001 m | 7.5 ms | 6.1x | 1.0 |
| n2/residue | 7.559 m | 4.928 m | 63.008 m | 0.074 m | 0.001 m | 7.4 ms | 6.2x | 1.0 |
| n4/isostatic | 1.133 m | 0.838 m | 9.095 m | 0.009 m | 0.000 m | 7.6 ms | 6.0x | 1.0 |
| n4/rate | 1.393 m | 1.147 m | 9.968 m | 0.009 m | 0.004 m | 7.7 ms | 5.9x | 1.0 |
| n4/residue | 1.117 m | 0.815 m | 9.126 m | 0.009 m | 0.004 m | 7.7 ms | 5.9x | 1.0 |
| n8/isostatic | 1.178 m | 0.882 m | 9.588 m | 0.010 m | 0.004 m | 8.0 ms | 5.7x | 1.0 |
| n8/rate | 1.173 m | 0.870 m | 9.526 m | 0.010 m | 0.000 m | 7.9 ms | 5.8x | 1.0 |
| n8/residue | 1.174 m | 0.881 m | 9.564 m | 0.010 m | 0.004 m | 8.0 ms | 5.7x | 1.0 |

### deglac3d

VE reference: solver cost = 1754.2 ms/step, n_solve = 8.6 sub-steps/step

VE scales: |rsl|max(PD) = 342.6 m, BSL min = -118.9 m, BSL(PD) = 1.51 m

| modal | rsl RMSE | rsl PD | rsl max | bsl RMSE | bsl PD | cost | speedup | nsolve |
|---|---|---|---|---|---|---|---|---|
| all | 20.423 m | 22.960 m | 393.184 m | 0.183 m | 0.000 m | 79.4 ms | 22.1x | 1.0 |
| n1/isostatic | 32.395 m | 28.298 m | 507.385 m | 0.238 m | 0.002 m | 14.1 ms | 124.4x | 1.0 |
| n1/rate | 33.339 m | 28.293 m | 508.758 m | 0.268 m | 0.002 m | 13.6 ms | 129.0x | 1.0 |
| n1/residue | 36.395 m | 28.397 m | 520.796 m | 0.256 m | 0.006 m | 13.7 ms | 128.0x | 1.0 |
| n2/isostatic | 24.419 m | 26.460 m | 450.154 m | 0.202 m | 0.001 m | 19.8 ms | 88.6x | 1.0 |
| n2/rate | 25.524 m | 28.165 m | 468.984 m | 0.204 m | 0.005 m | 19.9 ms | 88.2x | 1.0 |
| n2/residue | 26.946 m | 28.310 m | 472.937 m | 0.225 m | 0.007 m | 19.8 ms | 88.6x | 1.0 |
| n4/isostatic | 20.511 m | 22.997 m | 394.466 m | 0.184 m | 0.000 m | 31.9 ms | 55.0x | 1.0 |
| n4/rate | 21.127 m | 23.668 m | 408.528 m | 0.185 m | 0.000 m | 32.1 ms | 54.6x | 1.0 |
| n4/residue | 21.744 m | 23.591 m | 409.909 m | 0.187 m | 0.003 m | 32.0 ms | 54.8x | 1.0 |
| n8/isostatic | 20.431 m | 22.968 m | 393.284 m | 0.183 m | 0.000 m | 55.5 ms | 31.6x | 1.0 |
| n8/rate | 20.955 m | 23.407 m | 406.327 m | 0.181 m | 0.001 m | 55.6 ms | 31.6x | 1.0 |
| n8/residue | 21.625 m | 23.524 m | 408.684 m | 0.185 m | 0.003 m | 55.6 ms | 31.6x | 1.0 |

