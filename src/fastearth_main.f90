program fastearth
   !! Standalone FastEarth3D driver. Reads a &fe3d configuration (default
   !! "fastearth.nml", or the path given as the first argument) and runs a forced
   !! simulation over the ice-thickness forcing named in it.
   !!
   !!   ./bin/fastearth.x [config.nml] [defaults.nml]
   !!
   !! With one argument the config must be complete (it is its own defaults set).
   !! With two, the first is a SPARSE config overlaid on the complete defaults in
   !! the second (e.g. a run config over fastearth.nml) -- the convenient form for
   !! experiments that override only a few knobs.
   use fe_drive, only: fastearth_run
   implicit none
   character(len=512) :: cfg, defs

   if (command_argument_count() >= 1) then
      call get_command_argument(1, cfg)
   else
      cfg = "fastearth.nml"
   end if

   if (command_argument_count() >= 2) then
      call get_command_argument(2, defs)
      call fastearth_run(trim(cfg), defaults_file=trim(defs))
   else
      call fastearth_run(trim(cfg))
   end if
end program fastearth
