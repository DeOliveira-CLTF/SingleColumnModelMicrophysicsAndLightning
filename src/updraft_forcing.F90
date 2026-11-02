! updraft_forcing.F90
!
! Combines ERA5's own (diagnostic, synoptic-scale) vertical velocity
! with a synthetic convective-scale updraft pulse, for the saturation-
! perturbation ensemble. See saturation_perturbation.F90's header for
! why a pulse is needed at all: ERA5's real w at these scales is
! typically far too weak (~cm/s) to sustain a convective updraft long
! enough for Thompson to develop the full warm-rain -> ice -> graupel
! chain within a short test run -- exactly the same issue found with
! GOAMAZON's own prescribed wa.
!
! The pulse is a separable time x height shape: a half-cosine bump in
! time (rises from 0 to w_pulse_max and back to 0 over duration_s,
! smooth at both ends so it never shocks the column with a
! discontinuity) times a Gaussian bump in height centered on the
! saturated layer's own midpoint (so the pulse naturally reinforces
! ascent through the layer that was just saturated, tapering off
! above and below it).
!
! Combination with ERA5's own w: simple sum, not a max or blend -- the
! synoptic-scale w is small enough in practice that this choice rarely
! matters, but addition is the simplest, most defensible choice (real
! large-scale ascent/subsidence still contributes exactly as much as
! it would on its own; the pulse adds convective-scale lift on top of
! it, not instead of it).

module updraft_forcing

  use mpas_kind_types, only : RKIND

  implicit none
  private

  public :: combined_w

contains

  !===================================================================
  function combined_w(nz, z, w_era5, t_sec, z_center, sigma_z, &
       &               duration_s, w_pulse_max) result(w_total)
  !===================================================================
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: z(nz), w_era5(nz)
    real(kind=RKIND), intent(in) :: t_sec           ! seconds since run start
    real(kind=RKIND), intent(in) :: z_center, sigma_z ! pulse height center/width, m
    real(kind=RKIND), intent(in) :: duration_s      ! total pulse duration, s
    real(kind=RKIND), intent(in) :: w_pulse_max     ! peak amplitude, m/s
    real(kind=RKIND) :: w_total(nz)

    real(kind=RKIND), parameter :: pi = 3.14159265358979_RKIND
    real(kind=RKIND) :: time_shape
    integer :: k

    if (t_sec < 0.0_RKIND .or. t_sec > duration_s) then
      time_shape = 0.0_RKIND
    else
      ! half-cosine bump: 0 at t=0 and t=duration_s, 1 at the midpoint
      time_shape = 0.5_RKIND*(1.0_RKIND - cos(2.0_RKIND*pi*t_sec/duration_s))
    end if

    do k = 1, nz
      w_total(k) = w_era5(k) + w_pulse_max*time_shape* &
           &       exp(-0.5_RKIND*((z(k)-z_center)/sigma_z)**2)
    end do

  end function combined_w

end module updraft_forcing
