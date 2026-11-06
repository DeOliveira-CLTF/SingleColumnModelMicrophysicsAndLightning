! test_plume_precipitation.F90
!
! Checks: (1) below the critical threshold, nothing changes at all;
! (2) above threshold, removal never exceeds what would push
! ql below the threshold in one step; (3) qt decreases by EXACTLY
! the same amount as ql (water leaves the parcel, not just changes
! phase); (4) the recomputed theta_l, evaluated at the same T against
! the NEW ql, satisfies the defining equation to machine precision
! (this is the same self-consistency check style used in
! test_plume_condensation.F90).

program test_plume_precipitation

  use mpas_kind_types, only : RKIND
  use plume_precipitation, only : apply_precipitation

  implicit none

  real(kind=RKIND), parameter :: Rd_local = 287.04_RKIND
  real(kind=RKIND), parameter :: cp_local = 1004.64_RKIND
  real(kind=RKIND), parameter :: Lv_local = 2.5e6_RKIND
  real(kind=RKIND), parameter :: p0_local = 1.0e5_RKIND

  real(kind=RKIND) :: T, p, dz, w, qt, ql, theta_l, precip_removed
  real(kind=RKIND) :: qt0, ql0, exner, resid
  logical :: all_ok
  all_ok = .true.

  ! --- case 1: below threshold -- nothing should happen ---
  T = 275.0_RKIND; p = 60000.0_RKIND; dz = 300.0_RKIND; w = 10.0_RKIND
  qt = 0.008_RKIND; ql = 0.0005_RKIND   ! 0.5 g/kg, below the 1 g/kg critical value
  theta_l = 310.0_RKIND   ! arbitrary -- must come back unchanged
  call apply_precipitation(T, p, dz, w, qt, ql, theta_l, precip_removed)
  write(*,'(A,ES12.4,A,F8.5,A,F10.5)') 'CASO 1 (abaixo do limiar): precip_removed=', &
       precip_removed, '  ql=', ql, '  theta_l=', theta_l
  if (precip_removed /= 0.0_RKIND .or. ql /= 0.0005_RKIND .or. theta_l /= 310.0_RKIND &
       &  .or. qt /= 0.008_RKIND) then
    write(*,*) 'FALHA: nada deveria mudar abaixo do limiar'
    all_ok = .false.
  end if

  ! --- case 2: well above threshold ---
  T = 275.0_RKIND; p = 60000.0_RKIND; dz = 300.0_RKIND; w = 10.0_RKIND
  qt0 = 0.012_RKIND; ql0 = 0.004_RKIND   ! 4 g/kg, well above 1 g/kg critical
  qt = qt0; ql = ql0
  ! compute theta_l consistent with this (T,p,ql) BEFORE precipitation,
  ! so the "before" state is itself physically self-consistent
  exner = (p/p0_local)**(Rd_local/cp_local)
  theta_l = (T/exner) * exp(-Lv_local*ql/(cp_local*T))

  call apply_precipitation(T, p, dz, w, qt, ql, theta_l, precip_removed)

  resid = abs( (T/exner)*exp(-Lv_local*ql/(cp_local*T)) - theta_l )
  write(*,'(A,ES12.4,A,F8.5,A,ES14.6)') 'CASO 2 (acima do limiar): precip_removed=', &
       precip_removed, ' g/kg  ql_final=', ql*1000.0_RKIND, '  residuo theta_l=', resid

  if (abs((qt0-qt) - precip_removed) > 1.0e-12_RKIND) then
    write(*,*) 'FALHA: qt nao caiu exatamente pela agua removida'
    all_ok = .false.
  end if
  if (abs((ql0-ql) - precip_removed) > 1.0e-12_RKIND) then
    write(*,*) 'FALHA: ql nao caiu exatamente pela agua removida'
    all_ok = .false.
  end if
  if (resid > 1.0e-9_RKIND) then
    write(*,*) 'FALHA: theta_l recalculado nao satisfaz a equacao definidora'
    all_ok = .false.
  end if
  if (ql < 1.0e-3_RKIND - 1.0e-12_RKIND) then
    write(*,*) 'FALHA: ql nao deveria cair abaixo do limiar num unico passo'
    all_ok = .false.
  end if

  if (all_ok) then
    write(*,*) 'TODOS OS CHECKS PASSARAM'
  else
    write(*,*) 'ALGUM CHECK FALHOU -- ver acima'
  end if

end program test_plume_precipitation
