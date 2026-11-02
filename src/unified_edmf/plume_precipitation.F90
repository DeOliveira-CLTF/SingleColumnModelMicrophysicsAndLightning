! plume_precipitation.F90
!
! STAGE 5 (first piece): a simple Kessler (1969)-type autoconversion,
! removing liquid water from the plume once it exceeds a critical
! threshold -- the piece flagged as missing at the end of the coupled
! test (test_coupled_plume.F90 let ql_u grow to ~8 g/kg with nothing
! ever falling out, an unrealistically heavy parcel).
!
! Physically: precipitation fallout removes WATER MASS from the
! parcel, but does not un-release the latent heat that already warmed
! it -- the temperature T is unaffected at the instant of removal.
! Since theta_l is DEFINED from T, p and ql (Eq. 7 in Cohen et al.,
! 2020; the same formula plume_condensation.F90 solves), removing ql
! at fixed T changes what theta_l must be to represent that state
! correctly. This module removes ql (and the corresponding mass from
! qt), holds T fixed, and recomputes theta_l directly from the
! defining equation -- deliberately NOT calling saturation_adjustment
! again, since T is already known here and does not need re-solving.
!
! New code of ours -- not extracted from MPAS. This is a first,
! simple, spatially-formulated autoconversion (rate per unit HEIGHT
! rather than per unit time, dividing the usual Kessler rate by the
! plume's own vertical velocity w -- a faster-rising parcel spends
! less time per unit height and so converts less condensate per unit
! height, at the same rate constant). It is not tied to any of the
! Thompson/WSM6 autoconversion formulations elsewhere in this project;
! reconciling the two (or feeding this plume's precipitation into the
! grid-scale microphysics as a genuine mass source) is not attempted
! here.

module plume_precipitation

  use mpas_kind_types, only : RKIND

  implicit none
  private

  public :: apply_precipitation

  real(kind=RKIND), parameter :: Rd_local = 287.04_RKIND
  real(kind=RKIND), parameter :: cp_local = 1004.64_RKIND
  real(kind=RKIND), parameter :: Lv_local = 2.5e6_RKIND
  real(kind=RKIND), parameter :: p0_local = 1.0e5_RKIND

  ! ---- Kessler (1969) autoconversion parameters, standard values ----
  real(kind=RKIND), parameter :: k1_kessler    = 1.0e-3_RKIND   ! s-1
  real(kind=RKIND), parameter :: ql_crit_kessler = 1.0e-3_RKIND ! kg/kg (1 g/kg)

contains

  !===================================================================
  subroutine apply_precipitation(T, p, dz, w, qt, ql, theta_l, precip_removed)
  !===================================================================
    ! T, p     : temperature (K), pressure (Pa) -- T is INPUT ONLY,
    !            unchanged on exit (see header note)
    ! dz, w    : layer thickness (m), plume vertical velocity (m/s) --
    !            used only to convert the time-rate Kessler constant
    !            into a spatial (per-height) removal amount
    ! qt, ql   : total water, liquid water (kg/kg) -- updated in place
    ! theta_l  : liquid potential temperature (K) -- RECOMPUTED in
    !            place from the new ql at the same T (see header note)
    ! precip_removed : the ql mass removed this layer (kg/kg), a
    !            diagnostic for the caller (e.g. to accumulate a rain
    !            rate) -- not otherwise used inside this module
    real(kind=RKIND), intent(in) :: T, p, dz, w
    real(kind=RKIND), intent(inout) :: qt, ql, theta_l
    real(kind=RKIND), intent(out) :: precip_removed

    real(kind=RKIND) :: exner, dt_equiv, removal_rate_time

    if (ql <= ql_crit_kessler) then
      precip_removed = 0.0_RKIND
      return
    end if

    dt_equiv = dz / max(w, 1.0e-3_RKIND)   ! time-equivalent of crossing this layer
    removal_rate_time = k1_kessler * (ql - ql_crit_kessler)   ! kg/kg/s, Kessler's own rate
    precip_removed = min(removal_rate_time * dt_equiv, ql - ql_crit_kessler)
    ! (capped so ql never drops below ql_crit_kessler in one step, and
    ! never goes negative regardless of how large dt_equiv gets)

    ql = ql - precip_removed
    qt = qt - precip_removed   ! the water leaves the parcel entirely

    exner = (p/p0_local)**(Rd_local/cp_local)
    theta_l = (T/exner) * exp(-Lv_local*ql/(cp_local*T))

  end subroutine apply_precipitation

end module plume_precipitation
