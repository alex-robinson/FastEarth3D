# FastEarth3D

A state-of-the-art but simple and fast **3D solid-Earth model** — a
visco-elastic deformation model coupled with the sea-level equation — intended
as an **open-source replacement for VILMA** within the CLIMBER-X climate model.

The method is a clean-room reimplementation of the **spectral–finite-element,
time-domain** approach of Martinec (2000): spherical harmonics horizontally,
finite elements radially, an incompressible Maxwell rheology integrated
explicitly in time, a self-consistent sea-level equation with migrating
coastlines, and rotational feedback. It is built **3D-ready from the start**
(laterally varying viscosity) but validated first in 1D against the published
GIA benchmarks.

See [doc/design.md](doc/design.md) for the design rationale, the method
comparison, the validation ladder, and the implementation details that the
literature flags as easy to get wrong.

## Status

Scaffold. The build system, module architecture, and the spherical-harmonic
transform kernel (SHTns) are in place and tested; the physics modules are
interface stubs being filled in along the validation ladder.

## Dependencies

All provided by [fesm-utils](https://github.com/fesm-org/fesm-utils) (FFTW,
SHTns, the `fesmutils` helper library) plus a system netCDF. SHTns is built
through fesm-utils' `build.py`:

```bash
# in the fesm-utils checkout
./build.py -m macbook -c gfortran --component shtns --variant serial
```

## Build

```bash
# point a symlink at your fesm-utils checkout
ln -s ../fesm-utils fesm-utils

# generate the Makefile for your machine/compiler, then build + test
python config.py config/macbook_gfortran
make check
```

`make` switches: `debug=0|1|2`, `openmp=0|1`.

## Configure & run

All runtime parameters live in a single namelist group `&fe3d`, loaded into the
`fe_param_class` record by `fe_par_load`. [`fastearth.nml`](fastearth.nml) is the
complete, documented defaults set; a run can pass a sparse file overlaid on it
(yelmo `defaults_file` convention), overriding only what it needs. Time fields
(`dt_*`, `time_*`) are given in **years** and converted to SI seconds on load.

The earth structure is selected by `earth` — a named built-in (e.g.
`"M3-L70-V01"`) or `"custom"` to assemble from the surface-first layer arrays.
The memory integrator (`scheme = "trap"`) is advanced with the adaptive
time-stepper (`fe_timestep`), whose tolerances are also in `&fe3d`.

Build and run the standalone forced-run driver:

```bash
make fastearth                       # -> bin/fastearth.x
./bin/fastearth.x [config.nml]       # default config: fastearth.nml
```

It reads a reference equilibrium state (`file_ref`) and an ice-thickness forcing
(`file_forcing`, `h_ice(lon,lat,time)` on the model Gauss grid), marches the
model across the forcing, and writes the diagnostic surface fields to `file_out`.

Embedding the model in a host (the CLIMBER-X coupling path) uses the same API
behind a single `use fastearth3d`:

```fortran
use fastearth3d
type(fe_param_class) :: par
type(solid_earth)    :: se
call fe_par_load(par, "fastearth.nml")
call se%init(par, sht, z_bed_eq, h_ice_ref)
call se%update(h_ice, dt)            ! advance time -> time+dt; reads se%rsl, se%z_bed
```
