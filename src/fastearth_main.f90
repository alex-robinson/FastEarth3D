program fastearth
   !! Standalone FastEarth3D driver. Reads a &fe3d configuration (default
   !! "fastearth.nml", or the path given as the first argument) and runs a forced
   !! simulation over the ice-thickness forcing named in it. The config file is
   !! used as its own complete defaults set.
   !!
   !!   ./bin/fastearth.x [config.nml]
   use fe_drive, only: fastearth_run
   implicit none
   character(len=512) :: cfg

   if (command_argument_count() >= 1) then
      call get_command_argument(1, cfg)
   else
      cfg = "fastearth.nml"
   end if

   call fastearth_run(trim(cfg))
end program fastearth
