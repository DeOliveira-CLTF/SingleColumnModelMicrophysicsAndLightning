! test_plume_condensation.F90
!
! Validates saturation_adjustment() by self-consistency (no closed-form
! solution exists for this nonlinear problem, so instead of comparing
! to an exact answer, this checks that the converged state actually
! SATISFIES the two equations it is supposed to solve):
!
!   (a) qv + ql = qt                              (water conservation --
!       exact by construction, checked anyway as a safety net)
!   (b) qv == qsat_bolton(T,p)                     (on the saturated
!       branch: the parcel sits exactly on the saturation curve)
!   (c) T*exner*exp(-Lv*ql/(cp*T)) == theta_l      (the defining
!       equation for theta_l itself, evaluated at the returned T,ql)
!
! Also checks the trivial unsaturated branch (qt small -> ql=0, T
! exactly theta_l*exner, zero iterations used), and that
! plume_buoyancy_moist reduces to the dry buoyancy formula already
! validated in plume_thermo.F90 when ql=0 for both parcel and
! environment.

program test_plume_condensation

  use mpas_kind_types, only : RKIND
  use plume_condensation, only : saturation_adjustment, plume_buoyancy_moist, qsat_bolton
  use plume_thermo, only : plume_buoyancy

  implicit none

  real(kind=RKIND), parameter :: Rd_local = 287.04_RKIND
  real(kind=RKIND), parameter :: cp_local = 1004.64_RKIND
  real(kind=RKIND), parameter :: Lv_local = 2.5e6_RKIND
  real(kind=RKIND), parameter :: p0_local = 1.0e5_RKIND
  real(kind=RKIND), parameter :: g_local = 9.81_RKIND

  real(kind=RKIND) :: T, qv, ql, theta, exner, resid_thetal, resid_qv
  real(kind=RKIND) :: theta_l, qt, p
  integer :: niter
  logical :: all_ok
  real(kind=RKIND) :: B_moist, B_dry
  all_ok = .true.

  ! --- case 1: trivial unsaturated branch ---
  theta_l = 300.0_RKIND
  qt = 0.001_RKIND   ! 1 g/kg, well below saturation at these T,p
  p = 90000.0_RKIND
  call saturation_adjustment(theta_l, qt, p, T, qv, ql, theta, niter)
  exner = (p/p0_local)**(Rd_local/cp_local)
  write(*,'(A,ES14.6,A,I3)') 'CASO 1 (nao saturado): |T - theta_l*exner| = ', &
       abs(T - theta_l*exner), '  niter=', niter
  if (abs(T - theta_l*exner) > 1.0e-9_RKIND .or. ql /= 0.0_RKIND .or. niter /= 0) then
    write(*,*) 'FALHA: ramo nao saturado nao deveria precisar de iteracao'
    all_ok = .false.
  end if

  ! --- case 2: clearly supersaturated in the dry guess ---
  theta_l = 300.0_RKIND
  qt = 0.018_RKIND   ! 18 g/kg -- supersaturated in the dry first guess at this p
  p = 90000.0_RKIND
  call saturation_adjustment(theta_l, qt, p, T, qv, ql, theta, niter)

  resid_qv = abs(qv - qsat_bolton(T,p))
  resid_thetal = abs( (T/exner)*exp(-Lv_local*ql/(cp_local*T)) - theta_l )

  write(*,'(A,I3,A,ES14.6,A,ES14.6)') 'CASO 2 (saturado): niter=', niter, &
       '  |qv-qsat|=', resid_qv, '  |residuo theta_l|=', resid_thetal
  write(*,'(A,F10.4,A,F10.6,A,F10.6)') '  T=', T, ' K   qv=', qv*1000.0_RKIND, &
       ' g/kg   ql=', ql*1000.0_RKIND, ' g/kg'

  if (resid_qv > 1.0e-8_RKIND) then
    write(*,*) 'FALHA: qv nao esta na curva de saturacao'
    all_ok = .false.
  end if
  if (resid_thetal > 1.0e-6_RKIND) then
    write(*,*) 'FALHA: theta_l nao satisfeito pelo T,ql retornados'
    all_ok = .false.
  end if
  if (abs((qv+ql) - qt) > 1.0e-12_RKIND) then
    write(*,*) 'FALHA: agua nao conservada (qv+ql != qt)'
    all_ok = .false.
  end if

  ! --- case 3: plume_buoyancy_moist reduces to the dry formula when ql=0 ---
  B_moist = plume_buoyancy_moist(302.0_RKIND, 0.012_RKIND, 0.0_RKIND, &
       &                          300.0_RKIND, 0.010_RKIND, 0.0_RKIND, g_local)
  B_dry   = plume_buoyancy(302.0_RKIND, 0.012_RKIND, 300.0_RKIND, 0.010_RKIND, g_local)
  write(*,'(A,ES14.6,A,ES14.6)') 'CASO 3: B_moist(ql=0)=', B_moist, '  B_dry=', B_dry
  if (abs(B_moist - B_dry) > 1.0e-12_RKIND) then
    write(*,*) 'FALHA: buoyancy umida com ql=0 deveria bater exatamente com a seca'
    all_ok = .false.
  end if

  if (all_ok) then
    write(*,*) 'TODOS OS CHECKS PASSARAM'
  else
    write(*,*) 'ALGUM CHECK FALHOU -- ver acima'
  end if

end program test_plume_condensation
