! Stage accounting with fixed INITIAL dry-air layer masses. Read-only.
! Exact accounting != a conservative atmospheric water budget: the existing
! advective integrator and prescribed pressure do not solve mass continuity.
module sensitivity_budget
  use mpas_kind_types, only: RKIND
  implicit none
  private
  public :: sb_open,sb_begin,sb_transport,sb_before_mp,sb_finish,sb_close
  real(RKIND),allocatable :: mass(:)
  integer :: unit
  real(RKIND) :: before,after_adv,before_mp,horizontal,large_scale,pulse,clip
contains
  subroutine sb_open(path,nz,dz,p,theta,qv)
    character(*),intent(in) :: path
    integer,intent(in) :: nz
    real(RKIND),intent(in) :: dz(nz),p(nz),theta(nz),qv(nz)
    allocate(mass(nz))
    mass=p/(287.0_RKIND*theta*(p/1.e5_RKIND)**(2.0_RKIND/7.0_RKIND)*(1+qv/(287.0_RKIND/461.6_RKIND)))*dz
    open(newunit=unit,file=path,status='replace',action='write')
    write(unit,'(A)') 'time_s,water_start_ref,water_change_ref,horizontal_ref,LS_vertical_ref,pulse_vertical_ref,adv_clip_ref,adv_actual_ref,MYNN_stage_ref,MP_stage_ref,precip_kg_m2,MP_plus_precip_ref,accounting_residual_ref,adv_reconstruction_error_ref'
  end subroutine
  subroutine sb_begin(nz,state,dt,dz,rvadv,wls,wpulse)
    integer,intent(in) :: nz
    real(RKIND),intent(in) :: state(6,nz),dt,dz(nz),rvadv(nz),wls(nz),wpulse(nz)
    real(RKIND) :: gradient(nz),trial(nz)
    integer :: k
    before=sum(sum(state,dim=1)*mass)
    do k=2,nz-1
      gradient(k)=(state(1,k+1)-state(1,k-1))/(dz(k)+dz(k-1))
    enddo
    gradient(1)=(state(1,2)-state(1,1))/dz(1)
    gradient(nz)=(state(1,nz)-state(1,nz-1))/dz(nz-1)
    horizontal=sum(mass*dt*rvadv)
    large_scale=-sum(mass*dt*wls*gradient)
    pulse=-sum(mass*dt*wpulse*gradient)
    trial=state(1,:)+dt*(rvadv-(wls+wpulse)*gradient)
    clip=sum(mass*max(-trial,0.0_RKIND))
  end subroutine
  subroutine sb_transport(state)
    real(RKIND),intent(in) :: state(:,:)
    after_adv=sum(sum(state,dim=1)*mass)
  end subroutine
  subroutine sb_before_mp(state)
    real(RKIND),intent(in) :: state(:,:)
    before_mp=sum(sum(state,dim=1)*mass)
  end subroutine
  subroutine sb_finish(t,state,precip)
    real(RKIND),intent(in) :: t,state(:,:),precip
    real(RKIND) :: after,adv,mynn,mp,residual,error
    after=sum(sum(state,dim=1)*mass)
    adv=after_adv-before;mynn=before_mp-after_adv;mp=after-before_mp
    residual=after-before-adv-mynn-mp
    error=adv-horizontal-large_scale-pulse-clip
    write(unit,'(ES24.16,13(",",ES24.16))') t,before,after-before,horizontal,large_scale,pulse,clip,adv,mynn,mp,precip,mp+precip,residual,error
  end subroutine
  subroutine sb_close()
    close(unit)
    deallocate(mass)
  end subroutine
end module
