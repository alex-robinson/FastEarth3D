program test_restart
   !! netCDF restart round-trip for the coupling state (yelmo-convention I/O via
   !! fe_io: a time axis lets several snapshots share one file). Checks:
   !!
   !!   (1) direct state restore — read the last snapshot into a fresh model and
   !!       its z_bed/rsl match what was written;
   !!   (2) bit-for-bit continuation — restore the EARLIER snapshot's Maxwell
   !!       memory and step forward; the trajectory reproduces the uninterrupted
   !!       run exactly (the prognostic memory is all that is needed);
   !!   (3) multi-snapshot file — two snapshots at different times coexist, and a
   !!       specific earlier time can be selected on read.
   use fe_precision,       only: wp
   use fe_constants,       only: kyr, pi
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_sht,             only: sht_grid
   use fe_coupling,        only: solid_earth
   use fe_io,              only: fe_restart_write, fe_restart_read
   use ncio,               only: nc_size
   implicit none

   integer, parameter :: LMAX = 12, K1 = 3, K2 = 3
   character(len=*), parameter :: FILE = "obj/test_restart.nc"
   type(sht_grid), target :: sht
   type(earth_model)      :: e
   type(solid_earth)      :: a, b, c
   real(wp), allocatable  :: z_bed_eq(:,:), h_ice_ref(:,:), h_ice(:,:)
   real(wp), allocatable  :: a6_zbed(:,:), a6_rsl(:,:)
   real(wp) :: dt_couple, dt_step, t1, t2, d_restore, d_continue
   integer  :: step, nt
   logical  :: ok

   ok = .true.
   call sht%init(LMAX, nlat=2*LMAX, nphi=4*LMAX)
   e = build_M3L70V01()
   allocate(z_bed_eq(sht%nphi,sht%nlat), h_ice_ref(sht%nphi,sht%nlat), &
            h_ice(sht%nphi,sht%nlat))
   call make_fields(z_bed_eq, h_ice_ref, h_ice)

   dt_couple = 1.0_wp*kyr;  dt_step = 0.5_wp*kyr      ! n_sub = 2

   ! === reference run A =======================================================
   call a%init(e, sht, z_bed_eq, h_ice_ref, dt_couple, dt_step)
   do step = 1, K1
      call a%update(h_ice)
   end do
   t1 = a%time
   call fe_restart_write(a, FILE, t1, init=.true.)       ! snapshot 1 (memory @ K1)

   do step = 1, K2
      call a%update(h_ice)
   end do
   t2 = a%time
   call fe_restart_write(a, FILE, t2, init=.false.)      ! snapshot 2 (state @ K1+K2)
   a6_zbed = a%z_bed;  a6_rsl = a%rsl

   write(*,'(a,i0,a,f6.2,a,f6.2)') ' restart: lmax=', LMAX, &
        '   t1=', t1/kyr, ' kyr   t2=', t2/kyr

   ! === (3) multi-snapshot file ==============================================
   nt = nc_size(FILE, "time")
   write(*,'(a,i0)') '   snapshots in file: ', nt
   if (nt /= 2) then
      write(*,'(a)') '   FAIL: expected two time snapshots in one file'; ok = .false.
   end if

   ! === (1) direct state restore (default = last snapshot, t2) ================
   call b%init(e, sht, z_bed_eq, h_ice_ref, dt_couple, dt_step)
   call fe_restart_read(b, FILE)
   d_restore = max(maxval(abs(b%z_bed - a6_zbed)), maxval(abs(b%rsl - a6_rsl)))
   write(*,'(a,es11.2)') '   (1) state restore  max|B - A|     =', d_restore
   if (d_restore > 1.0e-9_wp) then
      write(*,'(a)') '   FAIL: restored state does not match the written state'; ok = .false.
   end if

   ! === (2) bit-for-bit continuation from the earlier snapshot (t1) ===========
   call c%init(e, sht, z_bed_eq, h_ice_ref, dt_couple, dt_step)
   call fe_restart_read(c, FILE, time=t1)                ! restore memory @ K1
   do step = 1, K2
      call c%update(h_ice)                               ! same load, K2 steps
   end do
   d_continue = max(maxval(abs(c%z_bed - a6_zbed)), maxval(abs(c%rsl - a6_rsl)))
   write(*,'(a,es11.2)') '   (2) continuation   max|C - A|     =', d_continue
   if (d_continue > 1.0e-9_wp) then
      write(*,'(a)') '   FAIL: restarted trajectory diverges from the uninterrupted run'
      ok = .false.
   end if

   call a%finalize();  call b%finalize();  call c%finalize()
   call sht%destroy()

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: netCDF restart restores state, continues bit-for-bit,'
      write(*,'(a)') '       and stores several time snapshots in one file'
   else
      write(*,'(a)') ' FAIL: restart round-trip did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine make_fields(z_bed_eq, h_ice_ref, h_ice)
      !! Polar land cap (colat<50°, +500 m) over deep ocean (−4000 m); no
      !! reference ice; a 2 km grounded ice load on colat<30°.
      real(wp), intent(out) :: z_bed_eq(:,:), h_ice_ref(:,:), h_ice(:,:)
      integer  :: i, j
      real(wp) :: thd
      do j = 1, sht%nlat
         thd = sht%colat(j)*180.0_wp/pi
         do i = 1, sht%nphi
            z_bed_eq(i,j) = merge(500.0_wp, -4000.0_wp, thd < 50.0_wp)
            h_ice(i,j)    = merge(2000.0_wp, 0.0_wp,     thd < 30.0_wp)
         end do
      end do
      h_ice_ref = 0.0_wp
   end subroutine make_fields

end program test_restart
