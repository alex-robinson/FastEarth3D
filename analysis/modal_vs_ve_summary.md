# Modal vs VE — accuracy & cost summary

Each modal run is compared against the VE run in the same set. `rsl RMSE` is area-weighted (cos lat) over all space-time; `rsl PD` at present day; `cost` is the per-coupling-step solver time; `speedup` = VE/modal.

### radial

VE reference: solver cost = 42.9 ms/step, n_solve = 1.6 sub-steps/step

VE scales: |rsl|max(PD) = 146.0 m, BSL min = -118.9 m, BSL(PD) = 1.51 m

| modal | rsl RMSE | rsl PD | rsl max | bsl RMSE | bsl PD | cost | speedup | nsolve |
|---|---|---|---|---|---|---|---|---|
| all | 0.812 m | 0.342 m | 7.790 m | 0.009 m | 0.000 m | 54.3 ms | 0.8x | 10.3 |
| n1/isostatic | 17.453 m | 6.049 m | 195.197 m | 0.102 m | 0.002 m | 27.4 ms | 1.6x | 5.2 |
| n1/rate | 19.354 m | 6.196 m | 213.002 m | 0.153 m | 0.002 m | 44.3 ms | 1.0x | 8.7 |
| n1/residue | 23.962 m | 7.217 m | 271.688 m | 0.127 m | 0.004 m | 24.2 ms | 1.8x | 4.5 |
| n2/isostatic | 6.081 m | 3.810 m | 61.558 m | 0.041 m | 0.002 m | 39.1 ms | 1.1x | 7.5 |
| n2/rate | 5.555 m | 4.927 m | 51.273 m | 0.039 m | 0.001 m | 38.2 ms | 1.1x | 7.3 |
| n2/residue | 7.657 m | 4.957 m | 62.168 m | 0.074 m | 0.001 m | 41.5 ms | 1.0x | 8.1 |
| n4/isostatic | 0.919 m | 0.360 m | 7.872 m | 0.011 m | 0.000 m | 39.7 ms | 1.1x | 7.9 |
| n4/rate | 1.043 m | 0.723 m | 8.287 m | 0.012 m | 0.005 m | 47.4 ms | 0.9x | 9.9 |
| n4/residue | 0.933 m | 0.387 m | 7.865 m | 0.011 m | 0.000 m | 48.0 ms | 0.9x | 9.6 |
| n8/isostatic | 0.813 m | 0.353 m | 7.764 m | 0.009 m | 0.000 m | 43.4 ms | 1.0x | 9.0 |
| n8/rate | 0.812 m | 0.341 m | 7.808 m | 0.009 m | 0.000 m | 51.5 ms | 0.8x | 10.8 |
| n8/residue | 0.813 m | 0.342 m | 7.793 m | 0.009 m | 0.000 m | 50.2 ms | 0.9x | 10.6 |

### deglac3d

VE reference: solver cost = 1769.4 ms/step, n_solve = 8.6 sub-steps/step

VE scales: |rsl|max(PD) = 342.6 m, BSL min = -118.9 m, BSL(PD) = 1.51 m

| modal | rsl RMSE | rsl PD | rsl max | bsl RMSE | bsl PD | cost | speedup | nsolve |
|---|---|---|---|---|---|---|---|---|
| all | 18.532 m | 23.739 m | 323.657 m | 0.143 m | 0.002 m | 185.3 ms | 9.5x | 22.0 |
| n1/isostatic | 30.352 m | 29.234 m | 421.542 m | 0.200 m | 0.003 m | 41.2 ms | 42.9x | 6.8 |
| n1/rate | 31.683 m | 28.172 m | 448.833 m | 0.264 m | 0.002 m | 71.3 ms | 24.8x | 12.3 |
| n1/residue | 32.726 m | 24.396 m | 403.089 m | 0.207 m | 0.007 m | 46.0 ms | 38.5x | 7.4 |
| n2/isostatic | 22.252 m | 26.689 m | 373.164 m | 0.162 m | 0.006 m | 59.0 ms | 30.0x | 9.4 |
| n2/rate | 20.589 m | 23.681 m | 362.593 m | 0.177 m | 0.007 m | 77.0 ms | 23.0x | 12.3 |
| n2/residue | 20.996 m | 22.814 m | 339.363 m | 0.177 m | 0.001 m | 59.9 ms | 29.5x | 10.4 |
| n4/isostatic | 18.636 m | 23.760 m | 325.696 m | 0.145 m | 0.002 m | 61.3 ms | 28.9x | 9.8 |
| n4/rate | 16.772 m | 19.788 m | 317.705 m | 0.158 m | 0.007 m | 159.1 ms | 11.1x | 26.3 |
| n4/residue | 16.840 m | 19.332 m | 293.908 m | 0.140 m | 0.001 m | 111.3 ms | 15.9x | 18.0 |
| n8/isostatic | 18.533 m | 23.740 m | 323.689 m | 0.143 m | 0.002 m | 162.9 ms | 10.9x | 22.2 |
| n8/rate | 16.657 m | 19.601 m | 315.618 m | 0.155 m | 0.007 m | 182.7 ms | 9.7x | 25.3 |
| n8/residue | 16.709 m | 19.269 m | 293.148 m | 0.138 m | 0.001 m | 228.2 ms | 7.8x | 31.8 |

