! test_plume_velocity.F90
!
! Validates integrate_velocity() against the exact analytical solution
! for constant B and eps: with u=w**2,
!
!   u(z) = u_eq + (u0 - u_eq)*exp(-2*a2*eps*z),  u_eq = a1*B/(a2*eps)
!
! and, for the eps=0 limit, u(z) = u0 + 2*a1*B*z (pure buoyant
! acceleration, no drag). Both are checked; the eps>0 case also
! confirms w approaches the terminal velocity sqrt(u_eq) far from the
! launch level, a physically meaningful sanity check in its own right
! (a buoyant, entraining plume in a uniform environment should level
! off at a terminal rise speed, not accelerate forever).

program test_plume_velocity

  use mpas_kind_types, only : RKIND
  use plume_velocity, only : integrate_velocity

  implicit none

  real(kind=RKIND), parameter :: a1 = 0.667_RKIND      ! virtual-mass coefficient (~1/(1+0.5))
  real(kind=RKIND), parameter :: a2 = 1.0_RKIND         ! drag coefficient
  real(kind=RKIND), parameter :: B0 = 0.02_RKIND        ! m/s2, constant buoyancy
  real(kind=RKIND), parameter :: eps0 = 0.0015_RKIND    ! m-1
  real(kind=RKIND), parameter :: w0 = 1.0_RKIND         ! m/s, launch vertical velocity
  real(kind=RKIND), parameter :: ztop = 4000.0_RKIND

  integer :: nz, k
  real(kind=RKIND) :: dz_uniform, err_max, u_eq
  real(kind=RKIND), allocatable :: z(:), dz(:), B(:), eps(:), w(:), w_exact(:)

  ! --- case 1: eps0 > 0 ---
  do nz = 80, 320, 80   ! three resolutions

    allocate(z(nz), dz(nz), B(nz), eps(nz), w(nz), w_exact(nz))
    dz_uniform = ztop/real(nz,kind=RKIND)
    do k = 1, nz
      dz(k) = dz_uniform
      z(k)  = (real(k,kind=RKIND)-0.5_RKIND)*dz_uniform
    end do
    B = B0
    eps = eps0
    w(1) = w0

    call integrate_velocity(nz, dz, B, eps, a1, a2, w)

    u_eq = a1*B0/(a2*eps0)
    do k = 1, nz
      w_exact(k) = sqrt(u_eq + (w0**2 - u_eq)*exp(-2.0_RKIND*a2*eps0*(z(k)-z(1))))
    end do
    err_max = maxval(abs(w(2:nz) - w_exact(2:nz)))
    write(*,'(A,I4,A,ES12.4,A,ES12.4,A,ES12.4)') 'CASO eps>0  nz=', nz, &
         '  dz=', dz_uniform, '  err_max=', err_max, '  w_terminal=', sqrt(u_eq)

    deallocate(z, dz, B, eps, w, w_exact)
  end do

  ! --- case 2: eps=0, pure buoyant acceleration, w(z)=sqrt(w0^2+2*a1*B0*z) ---
  nz = 100
  allocate(z(nz), dz(nz), B(nz), eps(nz), w(nz), w_exact(nz))
  dz_uniform = ztop/real(nz,kind=RKIND)
  do k = 1, nz
    dz(k) = dz_uniform
    z(k)  = (real(k,kind=RKIND)-0.5_RKIND)*dz_uniform
  end do
  B = B0
  eps = 0.0_RKIND
  w(1) = w0
  call integrate_velocity(nz, dz, B, eps, a1, a2, w)
  do k = 1, nz
    w_exact(k) = sqrt(w0**2 + 2.0_RKIND*a1*B0*(z(k)-z(1)))
  end do
  err_max = maxval(abs(w(2:nz) - w_exact(2:nz)))
  write(*,'(A,ES12.4)') 'CASO eps=0  err_max=', err_max

end program test_plume_velocity
