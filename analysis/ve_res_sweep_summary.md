# VE resolution & sub-step sweep — accuracy & cost

Reference: lmax=128, cfl=1. Every run is regridded onto the reference grid; `rsl RMSE` is area-weighted (cos lat) over all space-time, `rsl PD` at present day.

Reference cost: solid_earth_update = 14239.3 ms/step (drift 271.9 + memory 12802.4), n_solve = 8.6

Reference scales: |rsl|max(PD) = 340.5 m, BSL min = -120.0 m, BSL(PD) = 1.44 m

| run | lmax | cfl | rsl RMSE | rsl PD | rsl max | bsl PD | cost | drift | mem | nsolve | speedup |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lmax.32 | 32 | 1 | 10.543 m | 4.203 m | 169.175 m | 0.159 m | 392.4 ms | 50.6 | 227.8 | 8.4 | 36.29x |
| lmax.64 | 64 | 1 | 2.223 m | 1.736 m | 61.654 m | 0.078 m | 2481.9 ms | 166.5 | 1910.4 | 8.6 | 5.74x |
| lmax.128 | 128 | 1 | 0.000 m | 0.000 m | 0.000 m | 0.000 m | 14239.3 ms | 271.9 | 12802.4 | 8.6 | 1.00x |
| lmax.128.cfl0.5 | 128 | 0.5 | 0.064 m | 0.023 m | 0.816 m | 0.000 m | 25133.5 ms | 561.5 | 22132.5 | 16.1 | 0.57x |
| lmax.128.cfl1.5 | 128 | 1.5 | 0.075 m | 0.011 m | 0.735 m | 0.000 m | 6860.7 ms | 163.2 | 6003.0 | 5.4 | 2.08x |

