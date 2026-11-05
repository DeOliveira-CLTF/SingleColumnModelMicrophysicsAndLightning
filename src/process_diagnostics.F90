! Observational hooks for the serial offline column; disabled by default.
! Species order: vapor, cloud liquid, rain, cloud ice, snow, graupel.
module process_diagnostics
  use mpas_kind_types, only: RKIND
  implicit none
  logical :: pd_enabled=.false.
  integer, parameter :: npmax=96
  integer :: pd_unit, nprocess=0
  character(len=48) :: names(npmax)
  real(RKIND), allocatable :: increments(:,:,:), marker(:,:), marker_sum(:,:), initial(:,:)
contains
  subroutine pd_open(path,nz)
    character(*),intent(in) :: path
    integer,intent(in) :: nz
    if(allocated(increments)) deallocate(increments,marker,marker_sum,initial)
    allocate(increments(6,nz,npmax),marker(6,nz),marker_sum(6,nz),initial(6,nz))
    increments=0; marker=0; initial=0; nprocess=0; names=''
    open(newunit=pd_unit,file=path,status='replace',action='write')
    write(pd_unit,'(A)') 'time_s,k,process,dqv,dqc,dqr,dqi,dqs,dqg'
    pd_enabled=.true.
  end subroutine
  subroutine pd_begin(state)
    real(RKIND),intent(in) :: state(:,:)
    initial=state; increments=0; marker=0
  end subroutine
  subroutine pd_add(k,name,delta)
    integer,intent(in) :: k
    character(*),intent(in) :: name
    real(RKIND),intent(in) :: delta(6)
    integer :: j
    if(.not.pd_enabled) return
    do j=1,nprocess
      if(names(j)==name) exit
    end do
    if(j>nprocess) then
      if(j>npmax) error stop 'Too many process diagnostics'
      nprocess=j; names(j)=name
    end if
    increments(:,k,j)=increments(:,k,j)+delta
  end subroutine
  subroutine pd_transfer(k,name,donor,receiver,amount)
    integer,intent(in) :: k,donor,receiver
    character(*),intent(in) :: name
    real(RKIND),intent(in) :: amount
    real(RKIND) :: delta(6)
    if(.not.pd_enabled) return
    delta=0; delta(donor)=-amount; delta(receiver)=amount
    call pd_add(k,name,delta)
  end subroutine
  subroutine pd_mark(k,state)
    integer,intent(in) :: k
    real(RKIND),intent(in) :: state(6)
    if(pd_enabled) then
      marker(:,k)=state
      marker_sum(:,k)=sum(increments(:,k,:),dim=2)
    end if
  end subroutine
  subroutine pd_adjust(k,name,state)
    integer,intent(in) :: k
    character(*),intent(in) :: name
    real(RKIND),intent(in) :: state(6)
    if(pd_enabled) call pd_add(k,name,state-marker(:,k)-(sum(increments(:,k,:),dim=2)-marker_sum(:,k)))
  end subroutine
  subroutine pd_difference(k,name,state)
    integer,intent(in) :: k
    character(*),intent(in) :: name
    real(RKIND),intent(in) :: state(6)
    if(pd_enabled) call pd_add(k,name,state-marker(:,k))
  end subroutine
  subroutine pd_reconcile_tendency(k,state)
    integer,intent(in) :: k
    real(RKIND),intent(in) :: state(6)
    if(pd_enabled) call pd_add(k,'tendency_adjustments',state-sum(increments(:,k,:),dim=2))
  end subroutine
  subroutine pd_finish(time,state)
    real(RKIND),intent(in) :: time,state(:,:)
    integer :: k,j
    ! A diagnostic residual, never silently attributed to a physical process.
    do k=1,size(state,2)
      call pd_add(k,'unattributed_residual',state(:,k)-initial(:,k)-sum(increments(:,k,:),dim=2))
    end do
    do j=1,nprocess
      do k=1,size(state,2)
        if(any(increments(:,k,j)/=0._RKIND)) then
          write(pd_unit,'(ES25.16E3,",",I0,",",A,6(",",ES25.16E3))') time,k,trim(names(j)),increments(:,k,j)
        end if
      end do
    end do
  end subroutine
  subroutine pd_close()
    close(pd_unit); pd_enabled=.false.
  end subroutine
end module
