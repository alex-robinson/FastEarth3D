program diag_modal_sle
   !! DIAGNOSTIC (not pass/fail): does the modal-vs-VE gap appear in the
   !! SLE-COUPLED rsl, the quantity the ensemble actually scores? The per-degree
   !! probes (diag_modal_pblock, diag_modal_ramp) showed the displacement response
   !! u matches VE to ~1e-6 m, so the ensemble's ~1 m radial gap must live in the
   !! parts those probes skip: the sea-level equation (geoid N, ocean-load fixed
   !! point) and the LGM initial state. Note rsl = N − u + Δφ — and the earlier
   !! probes never checked N at all.
   !!
   !! This drives the FULL SLE (fe_sle) with a deglaciation history — a grounded
   !! cap unloaded to zero — through RESP_MODAL and RESP_VE in lockstep on the same
   !! earth (M3-L70-V01) and step cadence, and reports the area-weighted rsl gap
   !! (the ensemble metric). It runs the deglaciation TWICE to bisect the
   !! equilibration suspect:
   !!   * COLD  — memory at rest at LGM (no spin-up);
   !!   * SPUN  — the LGM cap held for n_spin steps first (relaxed initial state).
   !! If the gap is large COLD but small SPUN (or vice-versa), the LGM equilibration
   !! is implicated; if both are large, it is the SLE/geoid coupling itself; if both
   !! are ~1e-6 m, the gap is temporal (VE sub-steps, modal does not) or rotational.
   !!
   !! Rotation is OFF here (s_rot absent) — that is the next bisection if this one
   !! comes back clean.
   !!
   !!   usage:  diag_modal_sle.x [lmax] [n_spin] [n_degla] [dt_yr]
   !!   default:                  16    200      260       100
   use fe_precision,       only: wp
   use fe_constants,       only: kyr, pi, rho_ice, rho_water
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response, response_init_ve, response_init_modal, &
                                 response_set_dt, response_destroy
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, &
                                 sht_grid_surface_integral
   use fe_sle,             only: sle_solve, sle_solver, sle_result
   implicit none

   integer            :: lmax, n_spin, n_degla, k
   real(wp)           :: dt
   character(len=64)  :: arg
   type(sht_grid)     :: sht
   type(earth_model)  :: e
   real(wp), allocatable :: topo0(:,:), ice_lgm(:,:)
   real(wp) :: gap_pd_cold, gap_mx_cold, rlx_cold
   real(wp) :: gap_pd_spun, gap_mx_spun, rlx_spun

   lmax = 16;  n_spin = 200;  n_degla = 260;  dt = 100.0_wp
   if (command_argument_count() >= 1) then; call get_command_argument(1,arg); read(arg,*) lmax;    end if
   if (command_argument_count() >= 2) then; call get_command_argument(2,arg); read(arg,*) n_spin;  end if
   if (command_argument_count() >= 3) then; call get_command_argument(3,arg); read(arg,*) n_degla; end if
   if (command_argument_count() >= 4) then; call get_command_argument(4,arg); read(arg,*) dt;      end if
   dt = dt*1.0e-3_wp*kyr

   call sht_grid_init(sht, lmax, nlat=2*lmax, nphi=4*lmax)
   e = build_M3L70V01()
   allocate(topo0(sht%nphi,sht%nlat), ice_lgm(sht%nphi,sht%nlat))
   call make_fields(topo0, ice_lgm)

   write(*,'(a)') ' === modal vs VE through the SLE — M3-L70-V01, deglaciation rsl ==='
   write(*,'(a,i0,a,i0,a,i0,a,f6.1,a)') '   lmax=', lmax, '   n_spin=', n_spin, &
        '   n_degla=', n_degla, '   dt=', dt/kyr*1.0e3_wp, ' yr   (rotation off)'
   write(*,'(a)') ''

   call run_case(.false., gap_pd_cold, gap_mx_cold, rlx_cold)
   call run_case(.true.,  gap_pd_spun, gap_mx_spun, rlx_spun)

   write(*,'(a)') ' Area-weighted rsl gap  |modal - VE|  [m]  (ensemble metric):'
   write(*,'(a)') '   case            VE relax span   gap rms @PD    gap max @PD'
   write(*,'(a,es14.3,es15.3,es15.3)') '   COLD (no spin)', rlx_cold, gap_pd_cold, gap_mx_cold
   write(*,'(a,es14.3,es15.3,es15.3)') '   SPUN (relaxed)', rlx_spun, gap_pd_spun, gap_mx_spun
   write(*,'(a)') ''
   write(*,'(a)') ' Read: "VE relax span" is how far the VE rsl moves over the run (the signal'
   write(*,'(a)') ' scale). If a gap rms is a sizeable fraction of it, the SLE-coupled rsl is'
   write(*,'(a)') ' where modal departs from VE — and COLD vs SPUN says whether LGM'
   write(*,'(a)') ' equilibration drives it. ~1e-6 m in both => not SLE/equilibration (temporal/rotation).'

   call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   subroutine run_case(do_spin, gap_pd_rms, gap_pd_max, vrelax)
      !! Fresh modal + VE responses; optional LGM spin-up; deglaciate to ice-free in
      !! lockstep. Returns the area-weighted rsl gap at present day and the VE relax span.
      logical,  intent(in)  :: do_spin
      real(wp), intent(out) :: gap_pd_rms, gap_pd_max, vrelax
      type(response)   :: md, ve
      type(sle_solver) :: sle_md, sle_ve
      type(sle_result) :: res
      real(wp), allocatable :: d_ice(:,:), ice(:,:), Smd(:,:), Sve(:,:), Cmd(:,:), Cve(:,:)
      real(wp), allocatable :: Sve0(:,:)
      real(wp) :: f
      integer  :: j

      call response_init_ve(ve, e, sht, dt);  call response_set_dt(ve, dt)
      call response_init_modal(md, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
      call response_set_dt(md, dt)

      allocate(d_ice(sht%nphi,sht%nlat), ice(sht%nphi,sht%nlat))
      allocate(Smd(sht%nphi,sht%nlat), Sve(sht%nphi,sht%nlat))
      allocate(Cmd(sht%nphi,sht%nlat), Cve(sht%nphi,sht%nlat), Sve0(sht%nphi,sht%nlat))

      ! LGM spin-up: hold the full cap so the memory relaxes toward equilibrium.
      if (do_spin) then
         d_ice = ice_lgm;  ice = ice_lgm
         do j = 1, n_spin
            call sle_solve(sle_ve, sht, ve, d_ice, ice, topo0, Sve, Cve, res)
            call sle_solve(sle_md, sht, md, d_ice, ice, topo0, Smd, Cmd, res)
         end do
      end if
      Sve0 = Sve                       ! VE rsl entering the deglaciation

      ! Deglaciation: unload the cap linearly to zero.
      gap_pd_rms = 0.0_wp;  gap_pd_max = 0.0_wp
      do k = 1, n_degla
         f = 1.0_wp - real(k,wp)/real(n_degla,wp)
         d_ice = f*ice_lgm;  ice = d_ice
         call sle_solve(sle_ve, sht, ve, d_ice, ice, topo0, Sve, Cve, res)
         call sle_solve(sle_md, sht, md, d_ice, ice, topo0, Smd, Cmd, res)
      end do

      gap_pd_rms = warea_rms(Smd - Sve)
      gap_pd_max = maxval(abs(Smd - Sve))
      vrelax     = warea_rms(Sve - Sve0)

      call response_destroy(md);  call response_destroy(ve)
   end subroutine run_case

   real(wp) function warea_rms(d) result(r)
      !! Area-weighted (cos lat, via the SHT surface integral) RMS of a grid field.
      real(wp), intent(in) :: d(:,:)
      r = sqrt(sht_grid_surface_integral(sht, d*d) / (4.0_wp*pi))
   end function warea_rms

   subroutine make_fields(topo0, ice_lgm)
      !! Land cap (colat<60°, +500 m) over ocean (−4000 m); a 2 km grounded LGM ice
      !! cap (colat<40°), unloaded to zero over the deglaciation.
      real(wp), intent(out) :: topo0(:,:), ice_lgm(:,:)
      integer  :: i, j
      real(wp) :: th
      do j = 1, sht%nlat
         th = sht%colat(j)
         do i = 1, sht%nphi
            topo0(i,j)   = merge(500.0_wp, -4000.0_wp, th < 60.0_wp*pi/180.0_wp)
            ice_lgm(i,j) = merge(2000.0_wp, 0.0_wp,    th < 40.0_wp*pi/180.0_wp)
         end do
      end do
   end subroutine make_fields

end program diag_modal_sle
