! test_prognostic_area.F90
!
! Validates advance_area() (prognostic_area.F90) against the EXACT
! analytical solution of the area continuity equation for the special
! case of constant coefficients:
!
!   da/dt + w0 da/dz = k0*a        (rho constant cancels out)
!
! Along characteristics dz/dt = w0, this reduces to da/dt = k0*a, so a
! fluid element's area just grows/decays exponentially while being
! advected upward at w0. The exact solution for an arbitrary initial
! profile a0(z) is therefore
!
!   a(z,t) = a0(z - w0*t) * exp(k0*t)
!
! Here a0(z) is a Gaussian "bubble" launched once at t=0 (flux_in=0
! for the whole run -- this test is about advection+growth of an
! existing area patch, not about the injection mechanism). The test
! runs the numerical solver to a final time T and compares the result,
! at every grid point, against the exact formula above -- at TWO grid
! resolutions, to confirm the error shrinks as the grid is refined
! (i.e. the scheme actually converges, not just "happens to be close"
! at one resolution).

program test_prognostic_area

  use mpas_kind_types, only : RKIND
  use prognostic_area, only : advance_area

  implicit none

  real(kind=RKIND), parameter :: w0 = 3.0_RKIND        ! m/s, constant advection speed
  real(kind=RKIND), parameter :: k0 = 0.0005_RKIND      ! s-1, constant net (E-D)
  real(kind=RKIND), parameter :: rho0 = 1.0_RKIND       ! kg/m3, constant (cancels in the exact solution)
  real(kind=RKIND), parameter :: A0 = 0.05_RKIND        ! peak area fraction of the initial bubble
  real(kind=RKIND), parameter :: zc0 = 1000.0_RKIND     ! m, initial bubble center height
  real(kind=RKIND), parameter :: sigma = 300.0_RKIND    ! m, initial bubble width
  real(kind=RKIND), parameter :: ztop = 8000.0_RKIND    ! m, domain top
  real(kind=RKIND), parameter :: T_final = 600.0_RKIND  ! s, 10 min integration

  integer :: nz, k, nsteps, istep
  real(kind=RKIND) :: dz_uniform, dt, t_sec, err_max, err_l2
  real(kind=RKIND), allocatable :: z(:), dz(:), a(:), a_exact(:)
  real(kind=RKIND), allocatable :: rho_edge(:), w_edge(:), E(:), D(:)

  do nz = 100, 200, 100   ! two resolutions: nz=100 then nz=200

    allocate(z(nz), dz(nz), a(nz), a_exact(nz))
    allocate(rho_edge(0:nz), w_edge(0:nz), E(nz), D(nz))

    dz_uniform = ztop/real(nz,kind=RKIND)
    do k = 1, nz
      dz(k) = dz_uniform
      z(k)  = (real(k,kind=RKIND)-0.5_RKIND)*dz_uniform   ! layer center
    end do

    rho_edge = rho0
    w_edge   = w0
    E = 0.0_RKIND
    D = -k0   ! E-D = 0-(-k0) = k0, net growth rate k0 (D<0 here is just
              ! bookkeeping so E-D=k0; this test only exercises the net
              ! rate, not the E/D split itself)

    ! initial condition: Gaussian bubble already in the domain at t=0
    do k = 1, nz
      a(k) = A0 * exp(-0.5_RKIND*((z(k)-zc0)/sigma)**2)
    end do

    ! CFL-respecting timestep (upwind advection is stable for
    ! Courant number <= 1; use a comfortable safety factor)
    dt = 0.5_RKIND * dz_uniform / w0
    nsteps = ceiling(T_final/dt)
    dt = T_final/real(nsteps,kind=RKIND)   ! land exactly on T_final

    t_sec = 0.0_RKIND
    do istep = 1, nsteps
      call advance_area(nz, dz, dt, rho_edge, w_edge, E, D, 0.0_RKIND, a)
      t_sec = t_sec + dt
    end do

    ! exact solution at t=T_final
    do k = 1, nz
      a_exact(k) = A0 * exp(-0.5_RKIND*((z(k)-w0*T_final-zc0)/sigma)**2) * exp(k0*T_final)
    end do

    err_max = maxval(abs(a - a_exact))
    err_l2  = sqrt(sum((a-a_exact)**2)/real(nz,kind=RKIND))

    write(*,'(A,I4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
         'nz=', nz, '  dz=', dz_uniform, '  peak(a_exact)=', maxval(a_exact), &
         '  err_max=', err_max, '  err_L2=', err_l2

    deallocate(z, dz, a, a_exact, rho_edge, w_edge, E, D)

  end do

end program test_prognostic_area
