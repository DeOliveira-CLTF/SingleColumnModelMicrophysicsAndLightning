! plume_condensation.F90
!
! STAGE 4: condensation inside the plume, via saturation adjustment on
! the plume's diluted conserved variables (theta_l, liquid-ice
! potential temperature, and q_t, total water -- the same pair
! plume_thermo.F90 already dilutes, and the ones Cohen et al. 2020 use
! throughout). This is the piece flagged as missing at every earlier
! stage: without it, the plume's buoyancy stays "dry" and can never
! sustain deep, latent-heat-driven ascent, and the entrainment
! closure's RH_i trigger (Eq. 38) never has a genuine saturated state
! to act on.
!
! Numerics: classical fixed-point saturation adjustment (the same
! family as Bryan & Fritsch, 2002, cited by Cohen et al. 2020 for
! their own cutoff microphysics) --
!
!   T_(n+1) = theta_l * (p/p0)^(Rd/cp) * exp(+Lv*ql_(n)/(cp*T_(n)))
!   ql_(n+1) = max(qt - qvs(T_(n+1), p), 0)
!
! iterated from the unsaturated first guess (ql=0) to convergence.
! Saturation vapor pressure uses the Bolton (1980) formula, a standard
! Clausius-Clapeyron-based approximation.
!
! Simplifications flagged here, not yet relaxed: liquid water only (no
! ice -- so theta_l above is really theta_l for warm-rain physics, not
! the full liquid-ice potential temperature of Eq. 7 in Cohen et al.,
! which folds in q_i too); constant L_v and c_p (the papers note L_v
! is more properly a linear function of T, Romps 2008, easily added
! later).
!
! New code of ours -- not extracted from MPAS.

module plume_condensation

  use mpas_kind_types, only : RKIND

  implicit none
  private

  public :: saturation_adjustment, plume_buoyancy_moist, qsat_bolton

  real(kind=RKIND), parameter :: Rd_local = 287.04_RKIND    ! J/kg/K
  real(kind=RKIND), parameter :: cp_local = 1004.64_RKIND   ! J/kg/K
  real(kind=RKIND), parameter :: Lv_local = 2.5e6_RKIND     ! J/kg (constant -- see header note)
  real(kind=RKIND), parameter :: p0_local = 1.0e5_RKIND     ! Pa, reference pressure

contains

  !===================================================================
  subroutine saturation_adjustment(theta_l, qt, p, T, qv, ql, theta, niter_out)
  !===================================================================
    ! theta_l, qt, p : inputs -- liquid potential temperature (K),
    !                  total water specific humidity (kg/kg), pressure (Pa)
    ! T              : output, temperature, K
    ! qv, ql         : output, vapor and liquid specific humidity, kg/kg
    ! theta          : output, ordinary (dry) potential temperature, K
    ! niter_out      : optional, iterations used (for diagnostics/testing)
    !
    ! Numerics: Newton-Raphson on F(T) = T*exner*exp(-Lv*ql(T)/(cp*T)) -
    ! theta_l, with ql(T) = max(qt-qsat_bolton(T,p), 0) and a centered
    ! finite-difference derivative. An earlier version of this routine
    ! used a direct fixed-point recursion T_(n+1) = theta_l*exner*
    ! exp(+Lv*ql_(n)/(cp*T_(n))); that recursion turned out to be
    ! numerically UNSTABLE here (its local derivative exceeds 1 in
    ! magnitude near the root for these T/Lv/cp scales), oscillating
    ! between an under-saturated and an over-saturated branch without
    ! ever settling -- caught by this module's own test
    ! (test_plume_condensation.F90) checking the residual of the
    ! defining equation, not just "did the loop exit". Newton-Raphson
    ! does not share that instability (quadratic local convergence).
    real(kind=RKIND), intent(in) :: theta_l, qt, p
    real(kind=RKIND), intent(out) :: T, qv, ql, theta
    integer, intent(out), optional :: niter_out

    real(kind=RKIND) :: exner, F_here, Fp, Fm, dFdT, T_step
    real(kind=RKIND) :: qvs_here
    integer :: iter
    integer, parameter :: max_iter = 50
    real(kind=RKIND), parameter :: tol_F = 1.0e-8_RKIND   ! K, on the residual F(T)
    real(kind=RKIND), parameter :: dT_fd = 0.01_RKIND     ! K, finite-difference step

    exner = (p/p0_local)**(Rd_local/cp_local)

    ! first guess: fully unsaturated
    T = theta_l * exner
    qvs_here = qsat_bolton(T, p)

    if (qt <= qvs_here) then
      ! genuinely unsaturated: done, no iteration needed
      ql = 0.0_RKIND
      qv = qt
      theta = T / exner
      if (present(niter_out)) niter_out = 0
      return
    end if

    ! supersaturated in the dry guess: Newton-Raphson on F(T)=0
    do iter = 1, max_iter
      F_here = residual_F(T)
      if (abs(F_here) < tol_F) exit
      Fp = residual_F(T + dT_fd)
      Fm = residual_F(T - dT_fd)
      dFdT = (Fp - Fm) / (2.0_RKIND*dT_fd)
      T_step = F_here / dFdT
      T = T - T_step
    end do

    qvs_here = qsat_bolton(T, p)
    ql = max(qt - qvs_here, 0.0_RKIND)
    qv = qt - ql
    theta = T / exner
    if (present(niter_out)) niter_out = iter

  contains

    function residual_F(T_in) result(Fval)
      real(kind=RKIND), intent(in) :: T_in
      real(kind=RKIND) :: Fval, qvs_in, ql_in
      qvs_in = qsat_bolton(T_in, p)
      ql_in = max(qt - qvs_in, 0.0_RKIND)
      Fval = (T_in/exner)*exp(-Lv_local*ql_in/(cp_local*T_in)) - theta_l
    end function residual_F

  end subroutine saturation_adjustment

  !===================================================================
  function qsat_bolton(T, p) result(qvs)
  !===================================================================
    ! Saturation specific humidity via Bolton (1980)'s saturation
    ! vapor pressure formula (T in K, es in Pa), then the standard
    ! mixing-ratio approximation qvs = 0.622*es/(p-0.378*es).
    real(kind=RKIND), intent(in) :: T, p
    real(kind=RKIND) :: qvs, es, Tc
    Tc = T - 273.15_RKIND
    es = 611.2_RKIND * exp(17.67_RKIND*Tc/(Tc + 243.5_RKIND))
    qvs = 0.622_RKIND*es / max(p - 0.378_RKIND*es, 1.0_RKIND)
  end function qsat_bolton

  !===================================================================
  function plume_buoyancy_moist(theta_u, qv_u, ql_u, theta_env, qv_env, ql_env, g) &
       &     result(B)
  !===================================================================
    ! Buoyancy from virtual potential temperature WITH condensate
    ! loading: theta_v = theta*(1 + 0.608*qv - ql) -- water vapor
    ! makes air less dense (positive term), condensate loading makes
    ! it denser (negative term, ql subtracted). This replaces
    ! plume_thermo.F90's plume_buoyancy(), which is dry (ql=0 always,
    ! implicitly) -- that function is left in place for callers that
    ! deliberately want the dry-only Stage 2 behavior; this is the
    ! moist superset used from Stage 4 onward.
    real(kind=RKIND), intent(in) :: theta_u, qv_u, ql_u, theta_env, qv_env, ql_env, g
    real(kind=RKIND) :: B, thetav_u, thetav_env
    thetav_u   = theta_u   * (1.0_RKIND + 0.608_RKIND*qv_u   - ql_u)
    thetav_env = theta_env * (1.0_RKIND + 0.608_RKIND*qv_env - ql_env)
    B = g * (thetav_u - thetav_env) / thetav_env
  end function plume_buoyancy_moist

end module plume_condensation
