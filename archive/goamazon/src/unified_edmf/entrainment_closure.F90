! entrainment_closure.F90
!
! STAGE 3 of the "unified MYNN-EDMF" exploration: the actual dynamical
! and turbulent entrainment/detrainment closure of Cohen et al. (2020,
! JAMES), Section 3, replacing the constant/imposed eps,delta used in
! Stages 1-2. Restricted here to the two-subdomain case (one updraft,
! i, and its environment, j=0), which is all this project's plume
! model uses -- the general Eij/Delta_ij notation of the paper always
! specializes to i0/0i below.
!
! Equations implemented, with their numbers in Cohen et al. (2020):
!   - inverse time scale lambda_i0                      Eq. 30
!   - dry-mixing function D_i0 (sigmoid of buoyancy      Eq. 33-37
!     mixing fraction), and D_0i = 1 - D_i0 (from the
!     stated antisymmetry mu_ij = -mu_ji)
!   - evaporative enhancement function M_i0, M_0i        Eq. 38
!   - fractional entrainment/detrainment rates           Eq. 39-40
!   - turbulent entrainment rate                          Eq. 42-46
!
! Parameter values are the ones calibrated and reported in Table 2 of
! Cohen et al. (2020) -- not re-tuned here.
!
! New code of ours -- not extracted from MPAS. Deliberately still
! decoupled from plume_thermo.F90/plume_velocity.F90 (Stage 2): this
! module only computes eps/delta/eps_turb from buoyancy, velocity,
! area, TKE, relative humidity and a plume-height scale, however the
! caller obtains those. Wiring this into the Stage 1-2 modules (so
! that eps/delta stop being imposed constants and start being computed
! self-consistently from the evolving plume state) is the natural next
! step, not yet done here.

module entrainment_closure

  use mpas_kind_types, only : RKIND

  implicit none
  private

  public :: compute_entrainment_detrainment

  ! ---- Table 2 parameters (Cohen et al. 2020), used as reported ----
  real(kind=RKIND), parameter :: c_eps   = 0.13_RKIND   ! entrainment scaling constant
  real(kind=RKIND), parameter :: c_delta = 0.52_RKIND   ! detrainment scaling constant
  real(kind=RKIND), parameter :: c_lambda= 0.3_RKIND    ! weight of the TKE term in lambda_i0
  real(kind=RKIND), parameter :: beta_rh = 2.0_RKIND    ! detrainment relative-humidity power law
  real(kind=RKIND), parameter :: mu0     = 4.0e-4_RKIND ! s-1, sigmoid smoothness for D_i0
  real(kind=RKIND), parameter :: chi_i   = 0.25_RKIND   ! updraft-air fraction in the dry mixture
  real(kind=RKIND), parameter :: c_gamma = 0.075_RKIND  ! turbulent entrainment scaling constant

contains

  !===================================================================
  subroutine compute_entrainment_detrainment(b_i, b_0, w_i, w_0, a_i, &
       &     e0, RH_i, RH_0, H_i, eps_i0, delta_i0, eps_turb_i0)
  !===================================================================
    ! b_i, b_0    : updraft, environment buoyancy, m/s2
    ! w_i, w_0    : updraft, environment vertical velocity, m/s
    ! a_i         : updraft area fraction, 0-1 (environment is 1-a_i)
    ! e0          : environmental TKE, m2/s2
    ! RH_i, RH_0  : updraft, environment relative humidity, 0-1
    !               (RH_i==1 signals a saturated/cloudy updraft, the
    !               trigger for the evaporative M_i0/M_0i enhancement)
    ! H_i         : updraft height scale for turbulent entrainment, m
    !               (the real closure uses the updraft's own top
    !               height, tracked by the caller across time; this
    !               module takes it as given rather than tracking it
    !               itself, consistent with the rest of this project's
    !               "given profiles in, rates out" module design)
    real(kind=RKIND), intent(in) :: b_i, b_0, w_i, w_0, a_i, e0, RH_i, RH_0, H_i
    real(kind=RKIND), intent(out) :: eps_i0, delta_i0, eps_turb_i0

    real(kind=RKIND) :: lambda_i0, mu_i0, D_i0, D_0i, M_i0, M_0i
    real(kind=RKIND) :: cand1, cand2

    ! --- inverse time scale, Eq. 30: lambda_i0 = smin(|b|/|w|, c_lambda*|b|/sqrt(e0)) ---
    cand1 = abs(b_i - b_0) / max(abs(w_i - w_0), 1.0e-6_RKIND)
    cand2 = c_lambda * abs(b_i - b_0) / sqrt(max(e0, 1.0e-6_RKIND))
    lambda_i0 = smooth_min2(cand1, cand2)

    ! --- dry-mixing function, Eq. 33-37 ---
    ! mu_i0 = (b_mix - b_ij)/(w_i-w_0), with b_mix-b_ij = (b_i-b_0)*(chi_i - a_i)
    ! (using a_i+a_0=1, the two-subdomain specialization of Eq. 36)
    mu_i0 = (b_i - b_0) * (chi_i - a_i) / max(w_i - w_0, 1.0e-6_RKIND)
    D_i0 = 1.0_RKIND / (1.0_RKIND + exp(-mu_i0/mu0))
    ! D_0i = 1 - D_i0, from the stated antisymmetry mu_ij = -mu_ji (a
    ! sigmoid of the negated argument is exactly one minus the original)
    D_0i = 1.0_RKIND - D_i0

    ! --- evaporative enhancement, Eq. 38 ---
    ! M_i0 (used in entrainment E_i0): conditioned on the ENVIRONMENT (j=0)
    ! being saturated -- essentially never true, kept for completeness
    if (RH_0 >= 1.0_RKIND) then
      M_i0 = (max(RH_0**beta_rh - RH_i**beta_rh, 0.0_RKIND))**(1.0_RKIND/beta_rh)
    else
      M_i0 = 0.0_RKIND
    end if
    ! M_0i (used in detrainment Delta_i0): conditioned on the UPDRAFT (i)
    ! being saturated -- the physically important case (evaporative
    ! detrainment of a cloudy updraft mixing into a drier environment)
    if (RH_i >= 1.0_RKIND) then
      M_0i = (max(RH_i**beta_rh - RH_0**beta_rh, 0.0_RKIND))**(1.0_RKIND/beta_rh)
    else
      M_0i = 0.0_RKIND
    end if

    ! --- fractional entrainment/detrainment rates, Eq. 39-40 ---
    eps_i0   = (lambda_i0/max(w_i,1.0e-6_RKIND)) * (c_eps*D_i0 + c_delta*M_i0)
    delta_i0 = (lambda_i0/max(w_i,1.0e-6_RKIND)) * (c_eps*D_0i + c_delta*M_0i)

    ! --- turbulent entrainment, Eq. 42-46 ---
    eps_turb_i0 = 2.0_RKIND * c_gamma * sqrt(max(e0,0.0_RKIND)) / max(H_i, 100.0_RKIND)

  end subroutine compute_entrainment_detrainment

  !-----------------------------------------------------------------
  function smooth_min2(x1, x2) result(smin)
  !-----------------------------------------------------------------
    ! Smooth minimum of two positive candidates, following Lopez-Gomez
    ! et al. (2020) Appendix A, Eqs. A1-A4 (N=2 specialization), with
    ! epsilon=0.1 as used throughout that paper. The length-specific
    ! "at least 1 m" floor from Eq. A5 does NOT apply here (lambda_i0
    ! has units of 1/time, not length; that floor is particular to the
    ! mixing-length use of this same smooth-min function) -- a small
    ! floor is used instead purely to avoid division by zero.
    real(kind=RKIND), intent(in) :: x1, x2
    real(kind=RKIND) :: smin
    real(kind=RKIND) :: xmin, xt1, xt2, s1, s2, lambda_reg, w_lambert
    real(kind=RKIND), parameter :: eps_smin = 0.1_RKIND

    xmin = min(x1, x2)
    if (xmin <= 1.0e-12_RKIND) then
      smin = max(xmin, 0.0_RKIND)
      return
    end if

    ! Lambert W(1/e) for N=2 (n-1=1, so argument is 1/e exactly):
    ! W(1/e) = 0.27846454276...  (a fixed numerical constant here,
    ! since N=2 always -- no need for a general Lambert-W solver)
    w_lambert = 0.27846454276_RKIND
    lambda_reg = max(eps_smin/w_lambert * xmin, 1.0e-9_RKIND)

    xt1 = x1 - xmin
    xt2 = x2 - xmin
    s1 = exp(-xt1/lambda_reg)
    s2 = exp(-xt2/lambda_reg)
    smin = (x1*s1 + x2*s2) / (s1 + s2)

  end function smooth_min2

end module entrainment_closure
