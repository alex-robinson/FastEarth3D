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
