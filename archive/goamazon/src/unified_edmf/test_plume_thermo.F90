! test_plume_thermo.F90
!
! Validates dilute_plume_scalar() against the exact analytical
! solution for constant entrainment rate eps0 and a linearly
! stratified environment phi_env(z) = phi_env0 + gamma*z:
!
!   phi_u(z) = phi_env(z) - gamma/eps0
!              + [phi_u(0) - phi_env0 + gamma/eps0] * exp(-eps0*z)
!
! (derived in the plume_thermo.F90 header). Also checks the eps0->0
! limit (no entrainment: plume value stays exactly at its launch
! value, phi_u(z)=phi_u(0) for all z).

program test_plume_thermo

  use mpas_kind_types, only : RKIND
  use plume_thermo, only : dilute_plume_scalar

  implicit none

  real(kind=RKIND), parameter :: eps0 = 0.001_RKIND    ! m-1
  real(kind=RKIND), parameter :: phi_env0 = 300.0_RKIND ! K (theta_env at z=0)
  real(kind=RKIND), parameter :: gamma = 0.003_RKIND    ! K/m (environmental lapse rate)
  real(kind=RKIND), parameter :: phi_u0 = 302.0_RKIND   ! K (plume launch value, 2K warmer)
  real(kind=RKIND), parameter :: ztop = 3000.0_RKIND

  integer :: nz, k
  real(kind=RKIND) :: dz_uniform, err_max
  real(kind=RKIND), allocatable :: z(:), dz(:), phi_env(:), phi_u(:), phi_exact(:), eps(:)

  ! --- case 1: eps0 > 0, check against the full exact formula ---
  do nz = 60, 240, 60   ! three resolutions

    allocate(z(nz), dz(nz), phi_env(nz), phi_u(nz), phi_exact(nz), eps(nz))
    dz_uniform = ztop/real(nz,kind=RKIND)
    do k = 1, nz
      dz(k) = dz_uniform
      z(k)  = (real(k,kind=RKIND)-0.5_RKIND)*dz_uniform
    end do
    phi_env = phi_env0 + gamma*z
    eps = eps0
    phi_u(1) = phi_u0

    call dilute_plume_scalar(nz, dz, eps, phi_env, phi_u)

    do k = 1, nz
      phi_exact(k) = phi_env(k) - gamma/eps0 + &
           &         (phi_u0 - phi_env(1) + gamma/eps0) * exp(-eps0*(z(k)-z(1)))
    end do
    ! layer 1 keeps the caller-set launch value in both, so compare from k=2
    err_max = maxval(abs(phi_u(2:nz) - phi_exact(2:nz)))
    write(*,'(A,I4,A,ES12.4,A,ES12.4)') 'CASO eps>0  nz=', nz, &
         '  dz=', dz_uniform, '  err_max=', err_max

    deallocate(z, dz, phi_env, phi_u, phi_exact, eps)
  end do

  ! --- case 2: eps=0, plume value must stay exactly at the launch value ---
  nz = 100
  allocate(z(nz), dz(nz), phi_env(nz), phi_u(nz), eps(nz))
  dz_uniform = ztop/real(nz,kind=RKIND)
  do k = 1, nz
    dz(k) = dz_uniform
    z(k)  = (real(k,kind=RKIND)-0.5_RKIND)*dz_uniform
  end do
  phi_env = phi_env0 + gamma*z
  eps = 0.0_RKIND
  phi_u(1) = phi_u0
  call dilute_plume_scalar(nz, dz, eps, phi_env, phi_u)
  err_max = maxval(abs(phi_u - phi_u0))
  write(*,'(A,ES12.4)') 'CASO eps=0  err_max (deve ser ~0)=', err_max
  deallocate(z, dz, phi_env, phi_u, eps)

end program test_plume_thermo
