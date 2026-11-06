! Historical cosine-in-time x Gaussian-in-height perturbation only.
! No saturation adjustment and no transport/integrator changes.
module updraft_pulse
  use mpas_kind_types, only: RKIND
  implicit none
  private
  public :: pulse_profile
contains
  function pulse_profile(z,t,start,duration,center,sigma,amplitude) result(w)
    real(RKIND),intent(in) :: z(:),t,start,duration,center,sigma,amplitude
    real(RKIND) :: w(size(z)),tau,shape
    real(RKIND),parameter :: pi=3.14159265358979_RKIND
    w=0.0_RKIND
    if (amplitude==0.0_RKIND) return
    if (duration<=0.0_RKIND.or.sigma<=0.0_RKIND.or.amplitude<0.0_RKIND) stop 'Invalid pulse parameters'
    tau=t-start
    if (tau<0.0_RKIND.or.tau>duration) return
    shape=0.5_RKIND*(1.0_RKIND-cos(2.0_RKIND*pi*tau/duration))
    w=amplitude*shape*exp(-0.5_RKIND*((z-center)/sigma)**2)
  end function
end module
