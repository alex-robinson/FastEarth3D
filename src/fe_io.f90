module fe_io
   !! netCDF input/output for the coupling state, following the yelmo convention:
   !! a general, table-driven write_step writes a chosen subset of variables to a
   !! file with an unlimited time axis, so several snapshots accumulate in one
   !! file (handy for debugging). A restart file is just the full variable set
   !! written this way; fe_restart_read restores the prognostic state.
   !!
   !! The variable registry — names, dimensions, units, long_names — lives in a
   !! markdown table (input/fastearth-variables.md) parsed by variable_io. The
   !! table is both the I/O metadata source and the human documentation.
   !!
   !! Prognostic (restored on restart): the Maxwell memory-stress fields tau_*
   !! (nlam,ne,nk), the model time, and the adaptive controller's next-step Δt
   !! suggestion dt_try (so a restarted run sub-steps the same path → bit-for-bit
   !! continuation). Static (written once, checked on read):
   !! the reference state z_bed_eq, h_ice_ref. Diagnostic: h_ice, rsl, z_bed,
   !! C_ocean (lon,lat) — written for inspection and restored if present.
   use fe_precision,    only: wp
   use fe_constants,    only: rad2deg, sec_per_year
   use fe_viscoelastic, only: NLAM
   use fe_response,     only: ve_response
   use fe_coupling,     only: solid_earth
   use ncio
   use variable_io
   implicit none
   private

   public :: fe_restart_write, fe_restart_read, fe_write_step, fe_io_set_table

   ! The full time-varying variable set written each step (a restart writes all
   ! of these). Static reference fields are written once in the init branch.
   character(len=12), parameter :: ALL_VARS(15) = [character(len=12) :: &
        "tau_a_re",   "tau_a_im",   "tau_b_re", "tau_b_im", "tau_c_re", "tau_c_im", &
        "h_ice",      "rsl",        "z_bed",    "C_ocean",  "dt_try",   &
        "sigma_n_re", "sigma_n_im", "sigma_primed", "bsl"]

   character(len=256), save              :: table_file = "input/fastearth-variables.md"
   type(var_io_type), allocatable, save  :: vtable(:)

contains

   subroutine fe_io_set_table(filename)
      !! Override the variable-io table path (default input/fastearth-variables.md)
      !! and force a reload on next use.
      character(len=*), intent(in) :: filename
      table_file = filename
      if (allocated(vtable)) deallocate(vtable)
   end subroutine fe_io_set_table

   subroutine ensure_table()
      if (.not. allocated(vtable)) call load_var_io_table(vtable, trim(table_file))
   end subroutine ensure_table

   ! --- writing ---------------------------------------------------------------

   subroutine fe_restart_write(self, filename, time, init)
      !! Write the full coupling state as a restart snapshot at `time`. With
      !! init=.true. (default) the file is created; with init=.false. the
      !! snapshot is appended at a new time index.
      class(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename
      real(wp),           intent(in) :: time
      logical, optional,  intent(in) :: init
      call fe_write_step(self, filename, time, init=init)
   end subroutine fe_restart_write

   subroutine fe_write_step(self, filename, time, nms, init)
      !! General step writer. Writes the variables named in nms (default: all
      !! time-varying variables) at the time index for `time`, appending along
      !! the unlimited time axis. The model time should be passed as `time` for
      !! a restart to round-trip the clock.
      class(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename
      real(wp),           intent(in) :: time
      character(len=*), optional, intent(in) :: nms(:)
      logical,          optional, intent(in) :: init
      logical :: do_init
      integer :: ncid, n, q
      real(wp) :: tyr

      call ensure_table()
      do_init = .true.;  if (present(init)) do_init = init

      ! SI is internal; the output time axis is in YEARS (item §14a).
      tyr = time/sec_per_year

      if (do_init) call write_init(self, filename, time)

      call nc_open(filename, ncid, writable=.true.)
      if (do_init) then
         n = 1                                    ! the file was created with slice 1
      else
         n = nc_time_index(filename, "time", tyr, ncid)
         call nc_write(filename, "time", tyr, dim1="time", start=[n], count=[1], ncid=ncid)
      end if

      if (present(nms)) then
         do q = 1, size(nms);     call write_one(self, filename, trim(nms(q)), n, ncid);     end do
      else
         do q = 1, size(ALL_VARS); call write_one(self, filename, trim(ALL_VARS(q)), n, ncid); end do
      end if

      call nc_close(ncid)
   end subroutine fe_write_step

   subroutine write_init(self, filename, time)
      !! Create the file and its dimensions (lon, lat, nlam, ne, nk, time) and
      !! write the static reference fields.
      class(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename
      real(wp),           intent(in) :: time
      real(wp), allocatable :: lon_deg(:), lat_deg(:)

      lon_deg = self%sht%lon * rad2deg
      lat_deg = 90.0_wp - self%sht%colat * rad2deg

      call nc_create(filename, overwrite=.true.)
      call nc_write_dim(filename, "lon",  x=lon_deg, units="degrees_east")
      call nc_write_dim(filename, "lat",  x=lat_deg, units="degrees_north")
      call nc_write_dim(filename, "nlam", x=1, dx=1, nx=NLAM,          units="1")
      call nc_write_dim(filename, "ne",   x=1, dx=1, nx=self%resp%ne,  units="1")
      ! nk = # deforming (l>=1) coefficients, in the ve_response degree-grouped order
      call nc_write_dim(filename, "nk",   x=1, dx=1, nx=self%resp%nk,  units="1")
      ! nlm = # spherical-harmonic (l,m) coefficients (carries the σ_n load vector)
      call nc_write_dim(filename, "nlm",  x=1, dx=1, nx=self%resp%nlm, units="1")
      call nc_write_dim(filename, "time", x=time/sec_per_year, dx=1.0_wp, nx=1, &
           units="years", unlimited=.true.)

      call put2d(self, filename, "z_bed_eq",  self%z_bed_eq,  static=.true.)
      call put2d(self, filename, "h_ice_ref", self%h_ice_ref, static=.true.)
   end subroutine write_init

   subroutine write_one(self, filename, name, n, ncid)
      !! Dispatch a single variable name to its array + dimensions.
      class(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename, name
      integer,            intent(in) :: n, ncid
      select case (name)
      case ("tau_a_re"); call put3d(self, filename, name, self%resp%Are, n, ncid)
      case ("tau_a_im"); call put3d(self, filename, name, self%resp%Aim, n, ncid)
      case ("tau_b_re"); call put3d(self, filename, name, self%resp%Bre, n, ncid)
      case ("tau_b_im"); call put3d(self, filename, name, self%resp%Bim, n, ncid)
      case ("tau_c_re"); call put3d(self, filename, name, self%resp%Cre, n, ncid)
      case ("tau_c_im"); call put3d(self, filename, name, self%resp%Cim, n, ncid)
      case ("h_ice");    call put2d(self, filename, name, self%h_ice, n=n, ncid=ncid)
      case ("rsl");      call put2d(self, filename, name, self%rsl,   n=n, ncid=ncid)
      case ("z_bed");    call put2d(self, filename, name, self%z_bed, n=n, ncid=ncid)
      case ("C_ocean");  call put2d(self, filename, name, self%C,     n=n, ncid=ncid)
      case ("bsl");      call put_scalar(filename, name, self%bsl, n, ncid)
      case ("dt_try");   call put_scalar(filename, name, self%stepper%dt_try/sec_per_year, n, ncid)
      case ("sigma_n_re"); call put_sigma(self, filename, name, want_re=.true.,  n=n, ncid=ncid)
      case ("sigma_n_im"); call put_sigma(self, filename, name, want_re=.false., n=n, ncid=ncid)
      case ("sigma_primed")
         call put_scalar(filename, name, &
              merge(1.0_wp, 0.0_wp, allocated(self%resp%sigma_n) .and. self%resp%sigma_primed), &
              n, ncid)
      case default
         error stop 'fe_io: unknown variable "'//trim(name)//'"'
      end select
   end subroutine write_one

   subroutine put3d(self, filename, name, dat, n, ncid)
      !! Write a (nlam,ne,nk) Maxwell-memory field at time slice n.
      class(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename, name
      real(wp),           intent(in) :: dat(:,:,:)
      integer,            intent(in) :: n, ncid
      type(var_io_type) :: v
      call find_var_io_in_table(v, name, vtable, with_error=.true.)
      call nc_write(filename, name, dat, ncid=ncid, &
           dim1="nlam", dim2="ne", dim3="nk", dim4="time", &
           start=[1,1,1,n], count=[NLAM, self%resp%ne, self%resp%nk, 1], &
           units=trim(v%units), long_name=trim(v%long_name))
   end subroutine put3d

   subroutine put_sigma(self, filename, name, want_re, n, ncid)
      !! Write the real or imaginary part of the trapezoidal start-of-step load σ_n
      !! (a spectral (nlm) vector) at time slice n; zeros if σ_n is not yet tracked.
      class(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename, name
      logical,            intent(in) :: want_re
      integer,            intent(in) :: n, ncid
      type(var_io_type)     :: v
      real(wp), allocatable :: dat(:)
      integer :: nlm
      nlm = self%resp%nlm
      allocate(dat(nlm))
      if (allocated(self%resp%sigma_n)) then
         if (want_re) then;  dat = real(self%resp%sigma_n, wp)
         else;               dat = aimag(self%resp%sigma_n)
         end if
      else
         dat = 0.0_wp
      end if
      call find_var_io_in_table(v, name, vtable, with_error=.true.)
      call nc_write(filename, name, dat, ncid=ncid, dim1="nlm", dim2="time", &
           start=[1,n], count=[nlm,1], units=trim(v%units), long_name=trim(v%long_name))
   end subroutine put_sigma

   subroutine put_scalar(filename, name, val, n, ncid)
      !! Write a single time-varying scalar (controller state) at time slice n.
      character(len=*), intent(in) :: filename, name
      real(wp),         intent(in) :: val
      integer,          intent(in) :: n, ncid
      type(var_io_type) :: v
      call find_var_io_in_table(v, name, vtable, with_error=.true.)
      call nc_write(filename, name, val, ncid=ncid, dim1="time", &
           start=[n], count=[1], units=trim(v%units), long_name=trim(v%long_name))
   end subroutine put_scalar

   subroutine put2d(self, filename, name, dat, n, ncid, static)
      !! Write a (lon,lat) field — static (no time) or at time slice n.
      class(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename, name
      real(wp),           intent(in) :: dat(:,:)
      integer, optional,  intent(in) :: n, ncid
      logical, optional,  intent(in) :: static
      type(var_io_type) :: v
      logical :: is_static
      is_static = .false.;  if (present(static)) is_static = static
      call find_var_io_in_table(v, name, vtable, with_error=.true.)
      if (is_static) then
         call nc_write(filename, name, dat, dim1="lon", dim2="lat", &
              units=trim(v%units), long_name=trim(v%long_name))
      else
         call nc_write(filename, name, dat, ncid=ncid, &
              dim1="lon", dim2="lat", dim3="time", &
              start=[1,1,n], count=[self%sht%nphi, self%sht%nlat, 1], &
              units=trim(v%units), long_name=trim(v%long_name))
      end if
   end subroutine put2d

   ! --- reading ---------------------------------------------------------------

   subroutine fe_restart_read(self, filename, time)
      !! Restore the prognostic state into an already-initialised solid_earth
      !! (init() must have rebuilt the operators and grid). Reads the Maxwell
      !! memory + model time, and the diagnostic state if present, at the time
      !! slice matching `time` (default: the last slice). Validates that the
      !! file dimensions and reference fields match the initialised model.
      class(solid_earth), intent(inout) :: self
      character(len=*),   intent(in)    :: filename
      real(wp), optional, intent(in)    :: time
      real(wp), allocatable :: tvals(:), ref(:,:)
      integer :: nt, n, np, nl, ne, nk
      real(wp) :: tol

      np = self%sht%nphi;  nl = self%sht%nlat
      ne = self%resp%ne;   nk = self%resp%nk

      ! dimension validation
      if (nc_size(filename, "nlam") /= NLAM .or. nc_size(filename, "ne")  /= ne  .or. &
          nc_size(filename, "nk")   /= nk   .or. nc_size(filename, "lon") /= np  .or. &
          nc_size(filename, "lat")  /= nl   .or. nc_size(filename, "nlm") /= self%resp%nlm) &
         error stop 'fe_restart_read: file dimensions do not match the initialised model'

      ! select the time slice. The file stores YEARS (item §14a); convert to SI
      ! seconds immediately so the clock/comparison stay internally consistent.
      nt = nc_size(filename, "time")
      allocate(tvals(nt));  call nc_read(filename, "time", tvals)
      tvals = tvals*sec_per_year
      if (present(time)) then
         n = minloc(abs(tvals - time), dim=1)
         if (abs(tvals(n) - time) > 1.0e-6_wp*max(abs(time), 1.0_wp)) &
            error stop 'fe_restart_read: requested time not present in file'
      else
         n = nt
      end if

      ! reference-state validation (memory must not be paired with a different ref)
      allocate(ref(np,nl))
      tol = 1.0e-3_wp
      call nc_read(filename, "z_bed_eq", ref)
      if (maxval(abs(ref - self%z_bed_eq)) > tol) &
         error stop 'fe_restart_read: z_bed_eq does not match the initialised model'
      call nc_read(filename, "h_ice_ref", ref)
      if (maxval(abs(ref - self%h_ice_ref)) > tol) &
         error stop 'fe_restart_read: h_ice_ref does not match the initialised model'

      ! prognostic Maxwell memory at slice n
      call get3d(filename, "tau_a_re", self%resp%Are, ne, nk, n)
      call get3d(filename, "tau_a_im", self%resp%Aim, ne, nk, n)
      call get3d(filename, "tau_b_re", self%resp%Bre, ne, nk, n)
      call get3d(filename, "tau_b_im", self%resp%Bim, ne, nk, n)
      call get3d(filename, "tau_c_re", self%resp%Cre, ne, nk, n)
      call get3d(filename, "tau_c_im", self%resp%Cim, ne, nk, n)

      ! diagnostic current state (so the restarted object reports correctly)
      call get2d(filename, "h_ice",   self%h_ice, np, nl, n)
      call get2d(filename, "rsl",     self%rsl,   np, nl, n)
      call get2d(filename, "z_bed",   self%z_bed, np, nl, n)
      call get2d(filename, "C_ocean", self%C,     np, nl, n)

      ! restore the clock and the adaptive controller's step-size seed (so the
      ! restarted run sub-steps the same path as the uninterrupted one)
      self%resp%time = tvals(n)
      self%time      = tvals(n)
      call nc_read(filename, "dt_try", self%stepper%dt_try, start=[n], count=[1])
      self%stepper%dt_try = self%stepper%dt_try*sec_per_year   ! file stores years

      ! restore the trapezoidal start-of-step load σ_n (the other prognostic piece
      ! of the implicit scheme); without it the first step would re-derive σ_n to
      ! only the SLE solver tolerance, breaking bit-for-bit continuation.
      call restore_sigma(self, filename, n)
   end subroutine fe_restart_read

   subroutine restore_sigma(self, filename, n)
      !! Restore σ_n at time slice n into the response, marking it primed, if the
      !! snapshot recorded a tracked σ_n (sigma_primed flag set).
      class(solid_earth), intent(inout) :: self
      character(len=*),   intent(in)    :: filename
      integer,            intent(in)    :: n
      real(wp), allocatable :: sre(:), sim(:)
      real(wp) :: primed
      integer  :: nlm
      call nc_read(filename, "sigma_primed", primed, start=[n], count=[1])
      if (primed <= 0.5_wp) return
      nlm = self%resp%nlm
      allocate(sre(nlm), sim(nlm))
      call nc_read(filename, "sigma_n_re", sre, start=[1,n], count=[nlm,1])
      call nc_read(filename, "sigma_n_im", sim, start=[1,n], count=[nlm,1])
      call self%resp%prime_sigma(cmplx(sre, sim, wp))
   end subroutine restore_sigma

   subroutine get3d(filename, name, dat, ne, nk, n)
      character(len=*), intent(in)  :: filename, name
      real(wp),         intent(out) :: dat(:,:,:)
      integer,          intent(in)  :: ne, nk, n
      call nc_read(filename, name, dat, start=[1,1,1,n], count=[NLAM, ne, nk, 1])
   end subroutine get3d

   subroutine get2d(filename, name, dat, np, nl, n)
      character(len=*), intent(in)  :: filename, name
      real(wp),         intent(out) :: dat(:,:)
      integer,          intent(in)  :: np, nl, n
      call nc_read(filename, name, dat, start=[1,1,n], count=[np, nl, 1])
   end subroutine get2d

end module fe_io
