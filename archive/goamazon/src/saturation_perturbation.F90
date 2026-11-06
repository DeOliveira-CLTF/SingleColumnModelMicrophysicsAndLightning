! saturation_perturbation.F90
!
! Given a real background column (temperature, pressure, water vapor
! -- typically from era5_bridge.F90) and a chosen layer range, forces
! that layer to saturation (or a chosen supersaturation fraction),
! converting any excess vapor beyond saturation directly into cloud
! water (qc) or cloud ice (qi), depending on temperature. This is the
! "artificially saturated column" building block for the ERA5 x
! saturation-level ensemble: it turns one real ERA5 background profile
! into an idealized initial condition for testing Thompson's cloud
! formation response, independent of whether the real large-scale
! forcing at that point/time would ever actually produce saturation
! there.
!
! Uses RSLF/RSIF (module mpas_atmphys_functions.F, extracted unmodified
! from MPAS-Model) -- the EXACT saturation-mixing-ratio formulas
! Thompson itself uses internally -- so "saturated" here means exactly
! what Thompson itself considers saturated, not a separately-chosen
! formula that could disagree with it at the margin.
!
! Phase choice at the perturbed layer: liquid (rslf) for T > 253.15 K
! (-20C), ice (rsif) at or below that -- a simple, defensible cutover
! matching the temperature range where Thompson treats cloud water as
! the dominant condensate phase versus cloud ice; this is a modeling
! choice for BUILDING the initial condition, independent of Thompson's
! own internal mixed-phase treatment during the run that follows.

module saturation_perturbation

  use mpas_kind_types, only : RKIND
  use mpas_atmphys_functions, only : rslf, rsif

  implicit none
  private

  public :: saturate_layer

contains

  !===================================================================
  subroutine saturate_layer(nz, p, t, k_lo, k_hi, rh_target, qv, qc, qi)
  !===================================================================
    ! p(nz)  : pressure, Pa
    ! t(nz)  : temperature, K
    ! k_lo,k_hi : inclusive layer-index range to saturate (k_lo<=k_hi)
    ! rh_target : fraction of saturation to impose (1.0 = exactly
    !             saturated, 1.05 = 5% supersaturated -- Thompson will
    !             condense the excess on its very first timestep either
    !             way, so >1 mostly just guarantees condensate forms
    !             immediately rather than waiting on evaporation/mixing
    !             to nibble away a marginal perturbation first)
    ! qv,qc,qi(nz) : INOUT -- modified only for k_lo:k_hi; excess vapor
    !             above the target is moved into qc (T>-20C) or qi
    !             (T<=-20C), everything else (qr,qs,qg,ni,nr, and qv/qc/qi
    !             outside the layer) is left untouched by this routine
    integer, intent(in) :: nz, k_lo, k_hi
    real(kind=RKIND), intent(in) :: p(nz), t(nz), rh_target
    real(kind=RKIND), intent(inout) :: qv(nz), qc(nz), qi(nz)

    real(kind=RKIND) :: qsat, excess
    integer :: k

    do k = k_lo, k_hi
      if (t(k) > 253.15_RKIND) then
        qsat = rslf(p(k), t(k))
      else
        qsat = rsif(p(k), t(k))
      end if
      excess = max(qv(k) - rh_target*qsat, 0.0_RKIND)
      qv(k)  = qv(k) - excess
      if (t(k) > 253.15_RKIND) then
        qc(k) = qc(k) + excess
      else
        qi(k) = qi(k) + excess
      end if
      ! if the background was already sub-saturated at rh_target, top
      ! it up to exactly rh_target*qsat rather than leaving it dry --
      ! this is what makes the layer become saturated, not just "not
      ! supersaturated"
      qv(k) = rh_target*qsat
    end do
  end subroutine saturate_layer

end module saturation_perturbation
