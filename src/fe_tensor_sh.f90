module fe_tensor_sh
   !! Axisymmetric tensor spherical-harmonic dyadic transforms (Martinec 2000,
   !! Appendix B), the machinery rung 6 (laterally-varying viscosity) needs.
   !!
   !! The Maxwell memory stress τ and the strain ε are second-order symmetric
   !! tensors on the sphere, expanded in the spheroidal tensor spherical harmonics
   !! Z^λ, λ∈{1,2,5,6} (the toroidal λ=3,4 drop). For lateral viscosity the update
   !! τ⁺ = (1−M)τ − 2μM·ε is pointwise in PHYSICAL space, so the tensor must be
   !! reconstructed on the grid via its dyadic components (eqs 90/91, B10/B11) and
   !! projected back — scalar-synthesising the Z^λ coefficients is wrong (they are
   !! tensor, not scalar, harmonics).
   !!
   !! For an AXISYMMETRIC field (m=0) the φ-derivative harmonics vanish (F≡H≡0), so
   !! only four dyadic components are nonzero and each is a θ-only profile:
   !!   rr  =  Σ_j T¹_j Y_j        (Z¹: e_rr,  scalar Y)
   !!   rθ  =  Σ_j T²_j E_j        (Z²: e_rθ,  E=∂_θY)
   !!   tr  =  θθ+φφ = Σ_j (−2j(j+1)) T⁵_j Y_j   (Z⁵ trace)
   !!   df  =  θθ−φφ = Σ_j 2 T⁶_j G_j            (Z⁶, G=(∂_θθ−cotθ∂_θ)Y for m=0)
   !! with θθ = (tr+df)/2, φφ = (tr−df)/2 and rφ = θφ = 0. Each channel uses an
   !! orthogonal degree basis, so synthesis is a matrix·coeffs and analysis the
   !! Gauss-quadrature adjoint (coeff = ⟨field, basis⟩_w / ‖basis‖²_w) — exact,
   !! and a round trip is the identity (test_tensor_sh).
   !!
   !! The bases are bootstrapped from SHTns itself (Y via scalar synth, E via the
   !! spheroidal vector synth, G = −j(j+1)Y − 2cotθ·E from ∇₁²Y = −j(j+1)Y), so the
   !! normalisation is automatically consistent with the rest of FastEarth3D.
   use fe_precision, only: wp
   use fe_sht,       only: sht_grid
   implicit none
   private

   public :: tensor_sh
   ! Local component order 1..4 maps to Martinec's λ = 1,2,5,6 (as in strain_coeffs).
   integer, parameter, public :: TLAM = 4

   type :: tensor_sh
      !! Axisymmetric (m=0) tensor-SH dyadic transformer tied to one sht_grid.
      integer :: lmax = 0, nlat = 0
      ! Per-channel degree bases on the Gauss latitudes, (nlat, lmax):
      real(wp), allocatable :: Brr(:,:)   !! rr  channel: Y_j
      real(wp), allocatable :: Brt(:,:)   !! rθ  channel: E_j = ∂_θ Y_j
      real(wp), allocatable :: Btr(:,:)   !! tr  channel: −2j(j+1) Y_j
      real(wp), allocatable :: Bdf(:,:)   !! df  channel: 2 G_j
      ! Analysis norms ‖basis‖²_w = Σ_i w_i B(i,j)², (lmax). Zero ⇒ that degree
      ! has no such harmonic (Z⁶ at j=1) and its coefficient is identically zero.
      real(wp), allocatable :: nrr(:), nrt(:), ntr(:), ndf(:)
      real(wp), allocatable :: w(:)       !! Gauss weights (nlat), Σ = 2
   contains
      procedure :: init     => tensor_sh_init
      procedure :: synth    => tensor_sh_synth      !! coeffs (TLAM,lmax) -> dyad fields (nlat,4)
      procedure :: analysis => tensor_sh_analysis   !! dyad fields -> coeffs
      procedure :: destroy  => tensor_sh_destroy
   end type tensor_sh

   ! Dyadic-field column order returned by synth / consumed by analysis.
   integer, parameter, public :: DY_RR = 1, DY_RT = 2, DY_TR = 3, DY_DF = 4

contains

   subroutine tensor_sh_init(self, sht)
      !! Build the four degree bases from SHTns, axisymmetric (m=0).
      class(tensor_sh), intent(out) :: self
      type(sht_grid),   intent(in)  :: sht
      complex(wp), allocatable :: slm(:)
      real(wp),    allocatable :: y(:,:), vth(:,:), vph(:,:)
      real(wp) :: jj, cott
      integer  :: j, i, lm

      self%lmax = sht%lmax;  self%nlat = sht%nlat
      allocate(self%Brr(self%nlat,self%lmax), self%Brt(self%nlat,self%lmax))
      allocate(self%Btr(self%nlat,self%lmax), self%Bdf(self%nlat,self%lmax))
      allocate(self%nrr(self%lmax), self%nrt(self%lmax))
      allocate(self%ntr(self%lmax), self%ndf(self%lmax))
      allocate(self%w(self%nlat));  self%w = sht%gauss_w
      allocate(slm(sht%nlm), y(sht%nphi,sht%nlat), vth(sht%nphi,sht%nlat), &
               vph(sht%nphi,sht%nlat))

      do j = 1, self%lmax
         jj = real(j,wp)*real(j+1,wp)
         lm = sht%lmidx(j, 0)
         slm = (0.0_wp, 0.0_wp);  slm(lm) = (1.0_wp, 0.0_wp)
         call sht%synthesis(slm, y)             ! Y_j(θ)         (column 1 = axisym profile)
         call sht%sph_synthesis(slm, vth, vph)  ! E_j = ∂_θ Y_j  (vth)
         do i = 1, self%nlat
            cott = cos(sht%colat(i))/sin(sht%colat(i))
            self%Brr(i,j) = y(1,i)
            self%Brt(i,j) = vth(1,i)
            self%Btr(i,j) = -2.0_wp*jj*y(1,i)
            ! G_j = (∂_θθ − cotθ ∂_θ)Y_j = −j(j+1)Y_j − 2cotθ ∂_θ Y_j  (m=0)
            self%Bdf(i,j) = 2.0_wp*( -jj*y(1,i) - 2.0_wp*cott*vth(1,i) )
         end do
      end do

      ! Per-channel norms; degrees with a vanishing harmonic (Z⁶ at j=1, where the
      ! analytic norm is exactly 0 but roundoff leaves ~1e-30) are zeroed relative to
      ! the channel scale so analysis treats them as absent (coefficient ≡ 0) instead
      ! of dividing by a numerical-noise norm.
      self%nrr = clean(anorm(self%Brr, self%w))
      self%nrt = clean(anorm(self%Brt, self%w))
      self%ntr = clean(anorm(self%Btr, self%w))
      self%ndf = clean(anorm(self%Bdf, self%w))

      deallocate(slm, y, vth, vph)
   end subroutine tensor_sh_init

   pure function clean(d) result(dc)
      !! Zero entries negligible relative to the channel's largest norm.
      real(wp), intent(in) :: d(:)
      real(wp) :: dc(size(d))
      dc = merge(d, 0.0_wp, d > 1.0e-12_wp*maxval(d))
   end function clean

   pure function anorm(B, w) result(d)
      !! Per-degree quadrature norm Σ_i w_i B(i,j)².
      real(wp), intent(in) :: B(:,:), w(:)
      real(wp) :: d(size(B,2))
      integer  :: j
      do j = 1, size(B,2)
         d(j) = sum(w*B(:,j)*B(:,j))
      end do
   end function anorm

   subroutine tensor_sh_synth(self, c, dyad)
      !! Tensor-harmonic coefficients c(λ=1..4, degree) -> dyadic grid fields
      !! dyad(nlat, {rr,rθ,tr,df}).
      class(tensor_sh), intent(in)  :: self
      real(wp),         intent(in)  :: c(:,:)      !! (TLAM, lmax): rows λ=1,2,5,6
      real(wp),         intent(out) :: dyad(:,:)   !! (nlat, 4)
      dyad(:,DY_RR) = matmul(self%Brr, c(1,:))
      dyad(:,DY_RT) = matmul(self%Brt, c(2,:))
      dyad(:,DY_TR) = matmul(self%Btr, c(3,:))
      dyad(:,DY_DF) = matmul(self%Bdf, c(4,:))
   end subroutine tensor_sh_synth

   subroutine tensor_sh_analysis(self, dyad, c)
      !! Dyadic grid fields -> tensor-harmonic coefficients (orthogonal projection,
      !! Gauss quadrature). Degrees with a zero basis norm (Z⁶ at j=1) give 0.
      class(tensor_sh), intent(in)  :: self
      real(wp),         intent(in)  :: dyad(:,:)   !! (nlat, 4)
      real(wp),         intent(out) :: c(:,:)      !! (TLAM, lmax)
      call project(dyad(:,DY_RR), self%Brr, self%w, self%nrr, c(1,:))
      call project(dyad(:,DY_RT), self%Brt, self%w, self%nrt, c(2,:))
      call project(dyad(:,DY_TR), self%Btr, self%w, self%ntr, c(3,:))
      call project(dyad(:,DY_DF), self%Bdf, self%w, self%ndf, c(4,:))
   end subroutine tensor_sh_analysis

   pure subroutine project(field, B, w, nrm, c)
      !! c_j = ⟨field, B(:,j)⟩_w / ‖B(:,j)‖²_w, or 0 where the norm vanishes.
      real(wp), intent(in)  :: field(:), B(:,:), w(:), nrm(:)
      real(wp), intent(out) :: c(:)
      integer :: j
      do j = 1, size(B,2)
         if (nrm(j) > 0.0_wp) then
            c(j) = sum(w*field*B(:,j)) / nrm(j)
         else
            c(j) = 0.0_wp
         end if
      end do
   end subroutine project

   subroutine tensor_sh_destroy(self)
      class(tensor_sh), intent(inout) :: self
      if (allocated(self%Brr)) deallocate(self%Brr, self%Brt, self%Btr, self%Bdf)
      if (allocated(self%nrr)) deallocate(self%nrr, self%nrt, self%ntr, self%ndf)
      if (allocated(self%w))   deallocate(self%w)
      self%lmax = 0;  self%nlat = 0
   end subroutine tensor_sh_destroy

end module fe_tensor_sh
