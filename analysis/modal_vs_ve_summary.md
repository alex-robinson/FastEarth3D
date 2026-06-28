# Modal vs VE — accuracy & cost summary

Each modal run is compared against the VE run in the same set. `rsl RMSE` is area-weighted (cos lat) over all space-time; `rsl PD` at present day; `cost` is the per-coupling-step solver time; `speedup` = VE/modal.

### radial

VE reference: solver cost = 43.9 ms/step, n_solve = 1.6 sub-steps/step

VE scales: |rsl|max(PD) = 146.0 m, BSL min = -118.9 m, BSL(PD) = 1.51 m

| modal | rsl RMSE | rsl PD | rsl max | bsl RMSE | bsl PD | cost | speedup | nsolve |
|---|---|---|---|---|---|---|---|---|
| all | 1.183 m | 0.885 m | 9.653 m | 0.010 m | 0.004 m | 8.0 ms | 5.5x | 1.0 |
| n1/isostatic | 17.420 m | 6.026 m | 196.720 m | 0.101 m | 0.000 m | 7.8 ms | 5.6x | 1.0 |
| n1/rate | 19.362 m | 6.191 m | 214.367 m | 0.154 m | 0.002 m | 7.3 ms | 6.0x | 1.0 |
| n1/residue | 23.944 m | 7.198 m | 273.533 m | 0.126 m | 0.004 m | 7.4 ms | 5.9x | 1.0 |
| n2/isostatic | 5.776 m | 3.786 m | 59.044 m | 0.038 m | 0.002 m | 7.4 ms | 5.9x | 1.0 |
| n2/rate | 5.455 m | 4.912 m | 50.402 m | 0.037 m | 0.001 m | 7.5 ms | 5.9x | 1.0 |
| n2/residue | 7.559 m | 4.928 m | 63.008 m | 0.074 m | 0.001 m | 7.4 ms | 5.9x | 1.0 |
| n4/isostatic | 1.133 m | 0.838 m | 9.095 m | 0.009 m | 0.000 m | 7.3 ms | 6.0x | 1.0 |
| n4/rate | 1.393 m | 1.147 m | 9.968 m | 0.009 m | 0.004 m | 7.3 ms | 6.0x | 1.0 |
| n4/residue | 1.117 m | 0.815 m | 9.126 m | 0.009 m | 0.004 m | 7.6 ms | 5.8x | 1.0 |
| n8/isostatic | 1.178 m | 0.882 m | 9.588 m | 0.010 m | 0.004 m | 7.8 ms | 5.6x | 1.0 |
| n8/rate | 1.173 m | 0.870 m | 9.526 m | 0.010 m | 0.000 m | 7.8 ms | 5.6x | 1.0 |
| n8/residue | 1.174 m | 0.881 m | 9.564 m | 0.010 m | 0.004 m | 7.8 ms | 5.6x | 1.0 |

### deglac3d

VE reference: solver cost = 2030.6 ms/step, n_solve = 8.6 sub-steps/step

VE scales: |rsl|max(PD) = 342.6 m, BSL min = -118.9 m, BSL(PD) = 1.51 m

| modal | rsl RMSE | rsl PD | rsl max | bsl RMSE | bsl PD | cost | speedup | nsolve |
|---|---|---|---|---|---|---|---|---|
| all | 17.923 m | 23.192 m | 313.439 m | 0.137 m | 0.001 m | 11.7 ms | 173.6x | 1.0 |
| n1/isostatic | 30.360 m | 29.398 m | 421.241 m | 0.196 m | 0.002 m | 7.9 ms | 257.0x | 1.0 |
| n1/rate | 31.639 m | 28.219 m | 447.402 m | 0.261 m | 0.004 m | 7.7 ms | 263.7x | 1.0 |
| n1/residue | 32.738 m | 24.559 m | 401.706 m | 0.204 m | 0.007 m | 7.9 ms | 257.0x | 1.0 |
| n2/isostatic | 21.915 m | 26.563 m | 367.667 m | 0.158 m | 0.006 m | 8.4 ms | 241.7x | 1.0 |
| n2/rate | 20.446 m | 23.664 m | 359.963 m | 0.174 m | 0.007 m | 8.2 ms | 247.6x | 1.0 |
| n2/residue | 20.974 m | 22.961 m | 337.045 m | 0.172 m | 0.002 m | 8.3 ms | 244.7x | 1.0 |
| n4/isostatic | 18.033 m | 23.236 m | 315.221 m | 0.138 m | 0.002 m | 9.1 ms | 223.1x | 1.0 |
| n4/rate | 16.074 m | 19.113 m | 307.395 m | 0.152 m | 0.009 m | 8.9 ms | 228.2x | 1.0 |
| n4/residue | 16.244 m | 18.823 m | 283.494 m | 0.134 m | 0.003 m | 9.0 ms | 225.6x | 1.0 |
| n8/isostatic | 17.924 m | 23.193 m | 313.480 m | 0.137 m | 0.001 m | 13.8 ms | 147.1x | 1.0 |
| n8/rate | 15.935 m | 18.893 m | 305.512 m | 0.148 m | 0.008 m | 13.3 ms | 152.7x | 1.0 |
| n8/residue | 16.106 m | 18.761 m | 282.025 m | 0.132 m | 0.003 m | 13.7 ms | 148.2x | 1.0 |

