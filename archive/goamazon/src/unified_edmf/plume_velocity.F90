! plume_velocity.F90
!
! STAGE 2 (part 2): the plume vertical velocity equation, in the
! standard steady (quasi-equilibrium) entraining-plume form used
! throughout the EDMF literature (Siebesma et al., 2007; Tan et al.,
! 2018, before their own extension to a fully time-prognostic momentum
! budget, which additionally needs a perturbation-pressure closure --
! deliberately left for a later stage, see header note at the end):
!
!   w * dw/dz = a1*B - a2*eps*w**2
!
! -- buoyancy B accelerates the updraft; entrainment of slower
! environmental air acts as a drag term proportional to eps*w**2.
! a1 (~virtual mass effect, order 1) and a2 (drag coefficient, order 1)
! are the same nondimensional constants used across the EDMF papers
! cited to the user; both are left as caller-supplied arguments here
! rather than hardcoded, since calibrating them is exactly the kind of
! thing later stages (validation against LES) would revisit.
!
! Substituting u=w**2 makes this a LINEAR first-order ODE in u:
!
!   du/dz + 2*a2*eps*u = 2*a1*B
!
! which (for eps, B piecewise-constant over a layer) has the exact
! closed-form solution this module evaluates layer by layer -- same
! exact-per-layer-integration approach as plume_thermo.F90, for the
! same reasons (no stability concern, machine-precision agreement with
! the analytical test rather than merely first-order convergence).
!
! New code of ours -- not extracted from MPAS.

module plume_velocity

  use mpas_kind_types, only : RKIND

  implicit none
  private

  public :: integrate_velocity

contains

  !===================================================================
  subroutine integrate_velocity(nz, dz, B, eps, a1, a2, w)
  !===================================================================
    ! Integrates w*dw/dz = a1*B - a2*eps*w**2 upward from w(1) (already
    ! set by the caller as the launch value) through the column,
    ! overwriting w(2:nz). Stops (sets w=0 for all levels above) at the
    ! first level where w**2 would go negative -- the plume has lost
    ! all its vertical momentum and stops rising there, same
    ! termination convention as the real dmp_mf.
    !
    ! nz      : number of layers
    ! dz(nz)  : layer thickness, m
    ! B(nz)   : plume buoyancy, layer center, m/s2 (from plume_thermo's
    !           plume_buoyancy(), or any other source)
    ! eps(nz) : fractional entrainment rate, m-1 (same convention as
    !           plume_thermo -- per unit height)
    ! a1, a2  : nondimensional virtual-mass and drag coefficients
    ! w(nz)   : plume vertical velocity, layer center; w(1) is the
    !           launch value on entry, w(2:nz) overwritten on exit
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: dz(nz), B(nz), eps(nz), a1, a2
    real(kind=RKIND), intent(inout) :: w(nz)

    integer :: k
    real(kind=RKIND) :: u(nz), u_eq, decay_rate, stopped_above

    u(1) = w(1)**2
    stopped_above = 0.0_RKIND   ! 0 = still rising, 1 = plume has stopped

    do k = 2, nz
      if (stopped_above > 0.5_RKIND) then
        u(k) = 0.0_RKIND
        cycle
      end if

      if (eps(k) > 1.0e-12_RKIND) then
        decay_rate = 2.0_RKIND*a2*eps(k)
        u_eq = a1*B(k) / (a2*eps(k))   ! equilibrium (terminal) w**2:
                                        ! buoyancy production balances
                                        ! entrainment drag
        u(k) = u_eq + (u(k-1) - u_eq) * exp(-decay_rate*dz(k-1))
      else
        u(k) = u(k-1) + 2.0_RKIND*a1*B(k)*dz(k-1)
      end if

      if (u(k) <= 0.0_RKIND) then
        u(k) = 0.0_RKIND
        stopped_above = 1.0_RKIND
      end if
    end do

    w = sqrt(max(u, 0.0_RKIND))

  end subroutine integrate_velocity

end module plume_velocity

! SCOPE LEFT OUT of this stage (see the gap analysis discussed with
! the user): (1) the fully time-PROGNOSTIC momentum budget
! d(a*w)/dt + ... = ... , which needs its own perturbation-pressure
! closure (a separate piece of the extended-EDMF literature beyond
! Tan et al. 2018 itself) -- this module only solves the steady/
! quasi-equilibrium reduction; (2) condensation inside the plume (the
! buoyancy B passed in is still dry, from plume_thermo's current
! scope) -- until that lands, this w(z) cannot by itself sustain deep,
! phase-change-driven ascent, exactly the limitation flagged for the
! original MYNN dmp_mf, just reproduced here in the new formulation
! rather than solved; (3) the real Cohen et al. (2020) entrainment/
! detrainment closure -- eps here is still a caller-supplied profile,
! not yet computed from B/w**2.
