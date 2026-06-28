program diag_modal_latvisc
   !! DIAGNOSTIC (not pass/fail): isolate the modal LATERAL-VISCOSITY split-operator
   !! error against the VE tensor-SH ground truth, with lateral viscosity ACTIVE on
   !! BOTH paths. Earlier modal-vs-VE probes were all RADIAL and so missed this
   !! entirely; here a controlled-contrast synthetic earth drives modal+LV vs VE+LV
   !! through the identical pert_elem field (the only difference is the algorithm:
   !! depth-weighted split-operator vs full 6-component tensor-SH).
   !!
   !! Geometry: an axisymmetric polar viscosity cap (colat < theta_cap) of contrast
   !! Δ = log10(η/η_ref) applied to every Maxwell element (depth-uniform), with a
   !! co-located disc surface load (colat < theta_cap, cosine-tapered edge).
   !!
   !! Protocol (free-rebound = cleanest rate probe):
   !!   A. GLACIATION: hold the disc load for n_glac steps → builds depression.
   !!   B. REBOUND:    remove the load (σ = 0) for n_rebound steps → free relaxation;
   !!      the decay rate IS the lateral viscosity operator. A stiff cap (Δ>0) rebounds
   !!      slowly; the split-operator must reproduce VE's spatially-varying rate.
   !!
   !! Both fluid limits are contrast-independent, so any gap is a TRANSIENT RATE error.
   !! Reports, per (Δ, theta_cap): cap-region peak uplift modal vs VE and the ratio
   !! (the ensemble's "154 vs 343 m" analogue), plus the global area-rms gap, at the
   !! end of glaciation and at rebound checkpoints. mrbar range + active anomaly ranks
   !! are printed for insight.
   !!
   !! Rotation OFF. Earth = M3-L70-V01 (laterally-averaged reference profile).
   !!
   !!   usage:  diag_modal_latvisc.x [lmax] [dt_yr] [n_glac] [n_rebound]
   !!   default:                      48    200     500      200
   use fe_precision,       only: wp
   use fe_constants,       only: kyr, pi
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response, response_init_ve, response_init_modal, &
                                 response_set_dt, response_destroy, &
                                 response_begin_step, response_apply, response_commit_step, &
                                 response_enable_lateral_visc, response_enable_lateral_visc_modal
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, &
                                 sht_grid_synthesis, sht_grid_analysis, &
                                 sht_grid_surface_integral
   implicit none

   integer            :: lmax, n_glac, n_rebound, ic, it, ne
   real(wp)           :: dt
   character(len=64)  :: arg
   type(sht_grid)     :: sht
   type(earth_model)  :: e
   complex(wp), allocatable :: load_lm(:)
   real(wp),    allocatable :: load_g(:,:), pert(:,:,:)
   real(wp) :: theta_load
   ! sweep grids
   integer,  parameter :: NC = 4, NT = 3
   real(wp), parameter :: contrast(NC) = [0.5_wp, 1.0_wp, 1.5_wp, 2.0_wp]
   real(wp), parameter :: thcap_deg(NT) = [10.0_wp, 20.0_wp, 40.0_wp]
   ! rebound checkpoints (kyr after unload) to print
   integer,  parameter :: NCHK = 3
   real(wp) :: chk_kyr(NCHK) = [2.0_wp, 5.0_wp, 20.0_wp]

   lmax = 48;  dt = 200.0_wp;  n_glac = 500;  n_rebound = 200
   if (command_argument_count() >= 1) then; call get_command_argument(1,arg); read(arg,*) lmax;      end if
   if (command_argument_count() >= 2) then; call get_command_argument(2,arg); read(arg,*) dt;        end if
   if (command_argument_count() >= 3) then; call get_command_argument(3,arg); read(arg,*) n_glac;    end if
   if (command_argument_count() >= 4) then; call get_command_argument(4,arg); read(arg,*) n_rebound; end if
   dt = dt*1.0e-3_wp*kyr

   call sht_grid_init(sht, lmax, nlat=2*lmax, nphi=4*lmax)
   e = build_M3L70V01()
   allocate(load_lm(sht%nlm), load_g(sht%nphi,sht%nlat))

   write(*,'(a)')  ' === modal+LV vs VE+LV: lateral split-operator vs tensor-SH (M3-L70-V01) ==='
   write(*,'(a,i0,a,f6.1,a,i0,a,f6.1,a,i0,a,f6.1,a)') &
        '   lmax=', lmax, '  dt=', dt/kyr*1.0e3_wp, ' yr   glaciation=', n_glac, &
        ' steps (', n_glac*dt/kyr, ' kyr)   rebound=', n_rebound, ' steps (', n_rebound*dt/kyr, ' kyr)'
   write(*,'(a)')  '   stiff cap Δ=log10(η/η_ref)>0 over colat<theta_cap; co-located disc load.'
   write(*,'(a)')  '   peak = max uplift over cap region; ratio = modal/VE (1.0 = perfect).'
   write(*,'(a)')  ''

   do ic = 1, NT
      theta_load = thcap_deg(ic)*pi/180.0_wp
      call make_disc_load(sht, theta_load, load_g)
      call sht_grid_analysis(sht, load_g, load_lm)
      write(*,'(a,f5.1,a)') ' --- theta_cap = theta_load = ', thcap_deg(ic), ' deg ---------------------------------'
      write(*,'(a)') '   Δdex   end-glaciation                rebound peak ratio (modal/VE)       mrbar     ranks'
      write(*,'(a)') '          VE_pk     ratio   grms/pk |  2kyr    5kyr   20kyr  |  mean-mod          active'
      do it = 1, NC
         call run_one(contrast(it), thcap_deg(ic), load_lm)
      end do
      write(*,'(a)') ''
   end do

   call sht_grid_destroy(sht);  call radial_fe_finalize()
   deallocate(load_lm, load_g)

contains

   subroutine run_one(delta, thcap_deg_, load_lm)
      !! Build the cap-contrast field, enable LV on both responses, run glaciation +
      !! free rebound, and print the modal-vs-VE comparison row.
      real(wp),    intent(in) :: delta, thcap_deg_
      complex(wp), intent(in) :: load_lm(:)
      type(response)   :: md, ve
      complex(wp), allocatable :: zero_lm(:)
      real(wp) :: theta_cap
      real(wp) :: ve_pk_gl, md_pk_gl, gl_grms
      real(wp) :: rb_ratio(NCHK)
      real(wp) :: mrbar_lo, mrbar_hi
      integer  :: k, ichk, nchk_step(NCHK)
      real(wp) :: ve_pk, md_pk

      theta_cap = thcap_deg_*pi/180.0_wp
      call response_init_modal(md, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
      call response_init_ve(ve, e, sht, dt)
      call response_set_dt(md, dt);  call response_set_dt(ve, dt)
      ne = md%ne
      if (.not. allocated(pert)) allocate(pert(sht%nphi, sht%nlat, ne))
      call make_cap_pert(sht, ne, theta_cap, delta, pert)
      call response_enable_lateral_visc_modal(md, sht, pert)
      call response_enable_lateral_visc(ve, sht, pert)
      mrbar_lo = minval(md%mrbar);  mrbar_hi = maxval(md%mrbar)

      allocate(zero_lm(sht%nlm));  zero_lm = (0.0_wp, 0.0_wp)

      ! ---- A. glaciation: hold the load ------------------------------------
      do k = 1, n_glac
         call step(md, load_lm)
         call step(ve, load_lm)
      end do
      call peak_over_cap(md, theta_cap, md_pk_gl)
      call peak_over_cap(ve, theta_cap, ve_pk_gl)
      gl_grms = grms_gap(md, ve)

      ! ---- B. free rebound: remove the load --------------------------------
      do ichk = 1, NCHK
         nchk_step(ichk) = nint(chk_kyr(ichk)*kyr/dt)
      end do
      rb_ratio = -1.0_wp   ! sentinel: checkpoint not reached
      do k = 1, n_rebound
         call step(md, zero_lm)
         call step(ve, zero_lm)
         do ichk = 1, NCHK
            if (k == nchk_step(ichk)) then
               call peak_over_cap(md, theta_cap, md_pk)
               call peak_over_cap(ve, theta_cap, ve_pk)
               rb_ratio(ichk) = md_pk/max(abs(ve_pk), tiny(1.0_wp))
            end if
         end do
      end do

      write(*,'(f6.2,2x,es9.2,1x,f7.3,1x,f8.4,a,3(1x,f7.3),a,2x,f6.3,a,f6.3,3x,i0)') &
           delta, ve_pk_gl, md_pk_gl/max(abs(ve_pk_gl),tiny(1.0_wp)), &
           gl_grms/max(abs(ve_pk_gl),tiny(1.0_wp)), ' |', &
           rb_ratio(1), rb_ratio(2), rb_ratio(3), ' |', mrbar_lo, '-', mrbar_hi, md%nrank3d

      call response_destroy(md);  call response_destroy(ve)
      deallocate(zero_lm)
   end subroutine run_one

   subroutine step(resp, sigma_lm)
      !! One held-load step: apply (read state) then commit (advance φ / memory).
      type(response), intent(inout) :: resp
      complex(wp),    intent(in)    :: sigma_lm(:)
      complex(wp), allocatable :: u_lm(:), n_lm(:)
      allocate(u_lm(sht%nlm), n_lm(sht%nlm))
      call response_begin_step(resp, sht)
      call response_apply(resp, sht, sigma_lm, u_lm, n_lm)
      call response_commit_step(resp, sht, sigma_lm)
      deallocate(u_lm, n_lm)
   end subroutine step

   subroutine peak_over_cap(resp, theta_cap, pk)
      !! Max |uplift| over the cap region (colat < theta_cap), in load-σ units.
      type(response), intent(inout) :: resp
      real(wp),       intent(in)    :: theta_cap
      real(wp),       intent(out)   :: pk
      complex(wp), allocatable :: u_lm(:), n_lm(:)
      real(wp),    allocatable :: ug(:,:)
      integer :: i, j
      allocate(u_lm(sht%nlm), n_lm(sht%nlm), ug(sht%nphi,sht%nlat))
      call response_begin_step(resp, sht)
      call response_apply(resp, sht, mkzero(), u_lm, n_lm)   ! read current state (load irrelevant to apply)
      call sht_grid_synthesis(sht, u_lm, ug)
      pk = 0.0_wp
      do j = 1, sht%nlat
         if (sht%colat(j) >= theta_cap) cycle
         do i = 1, sht%nphi
            pk = max(pk, abs(ug(i,j)))
         end do
      end do
      deallocate(u_lm, n_lm, ug)
   end subroutine peak_over_cap

   function mkzero() result(z)
      complex(wp), allocatable :: z(:)
      allocate(z(sht%nlm));  z = (0.0_wp, 0.0_wp)
   end function mkzero

   real(wp) function grms_gap(md, ve) result(r)
      !! Global area-weighted rms of the modal-minus-VE uplift field [σ units].
      type(response), intent(inout) :: md, ve
      complex(wp), allocatable :: um(:), nm(:), uv(:), nv(:)
      real(wp),    allocatable :: gm(:,:), gv(:,:)
      allocate(um(sht%nlm), nm(sht%nlm), uv(sht%nlm), nv(sht%nlm))
      allocate(gm(sht%nphi,sht%nlat), gv(sht%nphi,sht%nlat))
      call response_begin_step(md, sht);  call response_apply(md, sht, mkzero(), um, nm)
      call response_begin_step(ve, sht);  call response_apply(ve, sht, mkzero(), uv, nv)
      call sht_grid_synthesis(sht, um, gm);  call sht_grid_synthesis(sht, uv, gv)
      r = sqrt(sht_grid_surface_integral(sht, (gm-gv)**2) / (4.0_wp*pi))
      deallocate(um, nm, uv, nv, gm, gv)
   end function grms_gap

   subroutine make_disc_load(sht, theta_load, g)
      !! Axisymmetric disc load: amplitude 1 inside colat<theta_load, cosine-tapered
      !! to 0 over a 5° edge. (Amplitude arbitrary — the operators are linear in σ.)
      type(sht_grid), intent(in)  :: sht
      real(wp),       intent(in)  :: theta_load
      real(wp),       intent(out) :: g(:,:)
      real(wp) :: taper, th, x
      integer  :: i, j
      taper = 5.0_wp*pi/180.0_wp
      do j = 1, sht%nlat
         th = sht%colat(j)
         if (th < theta_load) then
            x = 1.0_wp
         else if (th < theta_load + taper) then
            x = 0.5_wp*(1.0_wp + cos(pi*(th-theta_load)/taper))
         else
            x = 0.0_wp
         end if
         do i = 1, sht%nphi
            g(i,j) = x
         end do
      end do
   end subroutine make_disc_load

   subroutine make_cap_pert(sht, ne, theta_cap, delta, pert)
      !! Depth-uniform polar viscosity cap: pert = Δ inside colat<theta_cap (cosine-
      !! tapered over 5°), 0 outside, identical across all elements (non-Maxwell
      !! elements are ignored by both enablers). Δ = log10(η_local/η_ref).
      type(sht_grid), intent(in)  :: sht
      integer,        intent(in)  :: ne
      real(wp),       intent(in)  :: theta_cap, delta
      real(wp),       intent(out) :: pert(:,:,:)
      real(wp) :: taper, th, x
      integer  :: i, j
      taper = 5.0_wp*pi/180.0_wp
      do j = 1, sht%nlat
         th = sht%colat(j)
         if (th < theta_cap) then
            x = 1.0_wp
         else if (th < theta_cap + taper) then
            x = 0.5_wp*(1.0_wp + cos(pi*(th-theta_cap)/taper))
         else
            x = 0.0_wp
         end if
         do i = 1, sht%nphi
            pert(i,j,:) = delta*x
         end do
      end do
   end subroutine make_cap_pert

end program diag_modal_latvisc
