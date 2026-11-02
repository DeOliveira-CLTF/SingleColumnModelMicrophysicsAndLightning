! plume_thermo.F90
!
! STAGE 2 (part 1) of the "unified MYNN-EDMF" exploration: the plume
! scalar dilution equation
!
!   d(phi_u)/dz = -eps(z) * (phi_u - phi_env(z))
!
! -- entrainment mixes environmental air into the plume, pulling its
! properties (here: liquid-water potential temperature theta_l and
! total water q_t) toward the environment at a fractional rate eps(z)
! per unit height. This is the standard EDMF plume dilution equation
! (Tan et al., 2018; and going back to the classical entraining-plume
! literature). New code of ours -- not extracted from MPAS.
!
! DRY, NO CONDENSATION YET: buoyancy here is computed straight from the
! (diluted) virtual potential temperature, theta_v = theta*(1+0.608*qv)
! -- the same formula used throughout this project. Phase change
! inside the plume (the "4th change" flagged in the earlier gap
! analysis, needed before this can represent real convective clouds)
! is deliberately NOT part of this stage.
!
! Numerics: exact per-layer integration. Within one layer, eps and
! phi_env are held at their layer-center value (piecewise-constant
! approximation of whatever continuous profiles the caller supplies).
! Under that assumption the ODE above is LINEAR with constant
! coefficients and admits a closed-form solution (integrating factor),
! which this module evaluates exactly rather than by a generic
! finite-difference stepper -- this both sidesteps any stability
! concern and gives, for the validation test (test_plume_thermo.F90,
! where eps and the environmental lapse rate really are piecewise
! constant), agreement with the analytical solution to machine
! precision rather than merely first-order convergence.

module plume_thermo

  use mpas_kind_types, only : RKIND

  implicit none
  private

  public :: dilute_plume_scalar, plume_buoyancy

contains

  !===================================================================
  subroutine dilute_plume_scalar(nz, dz, eps, phi_env, phi_u)
  !===================================================================
    ! Integrates d(phi_u)/dz = -eps*(phi_u - phi_env) upward from
    ! phi_u(1) (already set by the caller as the surface/launch value)
    ! through the column, overwriting phi_u(2:nz).
    !
    ! nz            : number of layers
    ! dz(nz)        : layer thickness, m
    ! eps(nz)       : fractional entrainment rate, m-1 (per unit HEIGHT
    !                 here, not per unit time -- this module works in
    !                 z, not t; a caller with a rate in s-1 divides by
    !                 w first)
    ! phi_env(nz)   : environmental value of the scalar, layer center
    ! phi_u(nz)     : plume value; phi_u(1) is the launch/surface value
    !                 on entry, phi_u(2:nz) are overwritten on exit
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: dz(nz), eps(nz), phi_env(nz)
    real(kind=RKIND), intent(inout) :: phi_u(nz)

    integer :: k
    real(kind=RKIND) :: gamma_local, decay

    do k = 2, nz
      ! local environmental lapse rate over this layer (piecewise-
      ! constant approximation, using the two neighbouring layer
      ! centers -- one-sided at k=2 since there is no phi_env(0))
      if (k == 2) then
        gamma_local = (phi_env(k) - phi_env(k-1)) / dz(k-1)
      else
        gamma_local = (phi_env(k) - phi_env(k-2)) / (dz(k-1)+dz(k-2))
      end if

      if (eps(k) > 1.0e-12_RKIND) then
        decay = exp(-eps(k)*dz(k-1))
        ! exact solution of dpsi/dz + eps*psi = -gamma_local, psi=phi_u-phi_env,
        ! evaluated over one layer thickness dz(k-1) (from level k-1 to level k)
        phi_u(k) = phi_env(k) - gamma_local/eps(k) + &
             &     (phi_u(k-1) - phi_env(k-1) + gamma_local/eps(k)) * decay
      else
        ! no entrainment: plume value simply carries its previous value
        ! forward unchanged (eps->0 limit of the exact solution above)
        phi_u(k) = phi_u(k-1)
      end if
    end do

  end subroutine dilute_plume_scalar

  !===================================================================
  function plume_buoyancy(theta_u, qv_u, theta_env, qv_env, g) result(B)
  !===================================================================
    ! B = g*(thetav_u - thetav_env)/thetav_env -- dry virtual-
    ! temperature buoyancy, no condensate loading or ice yet (see
    ! header note).
    real(kind=RKIND), intent(in) :: theta_u, qv_u, theta_env, qv_env, g
    real(kind=RKIND) :: B, thetav_u, thetav_env
    thetav_u   = theta_u   * (1.0_RKIND + 0.608_RKIND*qv_u)
    thetav_env = theta_env * (1.0_RKIND + 0.608_RKIND*qv_env)
    B = g * (thetav_u - thetav_env) / thetav_env
  end function plume_buoyancy

end module plume_thermo
