module fastearth
   !! FastEarth3D umbrella module — single `use fastearth` entry point that
   !! re-exports the public API of every component. Host models and tests should
   !! depend on this rather than the individual fe_* modules.
   use fe_precision,       only: wp, sp, dp
   use fe_constants
   use fe_sht,             only: sht_grid
   use fe_earth_structure, only: earth_model
   use fe_radial_fe,       only: radial_operator
   use fe_viscoelastic,    only: viscoelastic_state
   use fe_gravity,         only: potential_perturbation
   use fe_sle,             only: sle_solver
   use fe_rotation,        only: rotation_state
   use fe_coupling,        only: solid_earth
   implicit none
   public

#ifndef VERSION
#define VERSION "unknown"
#endif
   character(len=*), parameter :: fastearth_version = VERSION

end module fastearth
