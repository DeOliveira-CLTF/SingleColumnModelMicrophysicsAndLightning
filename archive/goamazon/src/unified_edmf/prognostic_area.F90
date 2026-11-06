! prognostic_area.F90
!
! STAGE 1 of the "unified MYNN-EDMF" exploration discussed with the
! user: implements ONLY the prognostic plume-area continuity equation
! from Tan et al. (2018, JAMES) --
!
!   d(rho*a)/dt + d(rho*a*w)/dz = rho*a*(E - D)
!
! -- as a standalone, general-purpose finite-volume solver, decoupled
! from any thermodynamics, buoyancy, or microphysics. This is
! deliberately NOT wired into mynn_offline_coupling.F90 or any
! extracted MPAS file: it is new code of ours, meant to be validated
! on its own against a known analytical solution (see
! test_prognostic_area.F90) before any coupling is attempted. The
! real, future entrainment/detrainment closure (Cohen et al., 2020)
! and the buoyancy-generating vertical velocity equation are NOT part
! of this stage -- E, D, w, rho are all just INPUT profiles here,
! however they end up being computed later.
!
! Numerics: first-order upwind finite-volume in space (donor-cell,
! valid for w>=0 -- updrafts only, no downdraft case yet), explicit
! forward Euler in time for the flux-divergence term, combined with an
! EXACT exponential update for the (E-D) reaction term (exact for any
! dt when E-D is held fixed over the step, which removes one avoidable
! source of numerical error and isolates the advection scheme's own
! error in the validation test).
!
! Grid convention (matches bl_mynn_subroutines.F90's own zw/dz1
! convention, see mynn_offline_coupling.F90): a(k) is a layer-CENTER
! value, k=1..nz, bottom to top. w_edge(k), k=0..nz, is the vertical
! velocity at the BOTTOM edge of layer k (w_edge(0) = surface,
! w_edge(nz) = domain top); rho_edge(k) likewise. This lets a future
! caller reuse exactly the same z/zw/dz arrays already built for the
! MYNN coupling (build_vertical_grid in mp_offline_driver.F90).

module prognostic_area

  use mpas_kind_types, only : RKIND

  implicit none
  private

  public :: advance_area

contains

  !===================================================================
  subroutine advance_area(nz, dz, dt, rho_edge, w_edge, E, D, flux_in, a)
  !===================================================================
    ! nz            : number of layers
    ! dz(nz)        : layer thickness, m
    ! dt            : timestep, s
    ! rho_edge(0:nz): air density at layer edges, kg/m3
    ! w_edge(0:nz)  : plume vertical velocity at layer edges, m/s
    !                 (w_edge(k)>=0 assumed throughout -- updraft only)
    ! E(nz), D(nz)  : fractional entrainment/detrainment rate, s-1
    !                 (already turned into a true RATE -- i.e. already
    !                 multiplied by whatever w or 1/length scale the
    !                 eventual closure uses; this module takes them at
    !                 face value and does not care where they came from)
    ! flux_in       : rho*a*w flux of NEW plume area injected at the
    !                 surface (k=0 edge) this step, kg/m2/s -- this is
    !                 how a plume is "launched"; zero means no new
    !                 plume area enters from below this step
    ! a(nz)         : plume area fraction, layer-center, updated in place
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: dz(nz), dt
    real(kind=RKIND), intent(in) :: rho_edge(0:nz), w_edge(0:nz)
    real(kind=RKIND), intent(in) :: E(nz), D(nz)
    real(kind=RKIND), intent(in) :: flux_in
    real(kind=RKIND), intent(inout) :: a(nz)

    real(kind=RKIND) :: flux(0:nz)   ! rho*a*w at each edge, kg/m2/s
    real(kind=RKIND) :: rho_center(nz), a_new(nz)
    integer :: k

    ! layer-center density, simple average of the two bounding edges
    ! (consistent with a layer-center a(k) and edge-defined rho/w)
    do k = 1, nz
      rho_center(k) = 0.5_RKIND*(rho_edge(k-1) + rho_edge(k))
    end do

    ! upwind (donor-cell) fluxes: with w_edge>=0 everywhere, the flux
    ! through edge k is carried by the layer BELOW it (k, using 1-based
    ! layer indexing, sits below edge k) -- except edge 0 (the surface),
    ! which has no layer below it inside the domain and is instead the
    ! prescribed injection flux_in (how a new plume is launched)
    flux(0) = flux_in
    do k = 1, nz-1
      flux(k) = rho_edge(k)*w_edge(k)*a(k)
    end do
    flux(nz) = rho_edge(nz)*w_edge(nz)*a(nz)   ! outflow at the domain top

    ! finite-volume update of the flux-divergence (advection) term only
    do k = 1, nz
      a_new(k) = a(k) - (dt/(rho_center(k)*dz(k))) * (flux(k) - flux(k-1))
    end do

    ! exact exponential update of the (E-D) reaction term, held fixed
    ! over this step -- see header note on why this is done separately
    ! from (and after) the advection sub-step
    do k = 1, nz
      a_new(k) = a_new(k) * exp((E(k)-D(k))*dt)
    end do

    a = max(a_new, 0.0_RKIND)   ! area fraction cannot be negative

  end subroutine advance_area

end module prognostic_area
