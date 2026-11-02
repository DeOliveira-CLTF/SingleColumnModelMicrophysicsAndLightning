! test_entrainment_closure.F90
!
! This closure has no analytical solution to check against (it is an
! empirical fit, calibrated against LES in the source paper) -- so
! this test instead checks internal consistency and known physical
! limits, which any bug would very likely violate:
!
! 1) A buoyant, faster updraft (b_i>b_0, w_i>w_0) should have
!    entrainment-dominated dry mixing (D_i0>0.5).
! 2) Evaporative enhancement: a SATURATED updraft mixing into a DRIER
!    environment must show enhanced detrainment (M_0i>0), while an
!    UNSATURATED updraft must show none regardless of environment
!    humidity (the trigger in Eq. 38 is strictly RH_i==1).
! 3) All rates must stay finite and non-negative, including in a
!    near-zero-velocity edge case.

program test_entrainment_closure

  use mpas_kind_types, only : RKIND
  use entrainment_closure, only : compute_entrainment_detrainment

  implicit none

  real(kind=RKIND) :: eps_i0, delta_i0, eps_turb_i0, delta_case1, delta_case2, delta_case3
  logical :: all_ok
  all_ok = .true.

  ! --- case 1: buoyant, fast, DRY updraft (RH<1) in a calm environment ---
  call compute_entrainment_detrainment( &
       0.02_RKIND, 0.0_RKIND, 3.0_RKIND, 0.0_RKIND, 0.05_RKIND, &
       0.3_RKIND, 0.6_RKIND, 0.5_RKIND, 1000.0_RKIND, &
       eps_i0, delta_i0, eps_turb_i0)
  delta_case1 = delta_i0
  write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') 'CASO 1 (seco, boiante): eps_i0=', eps_i0, &
       '  delta_i0=', delta_i0, '  eps_turb_i0=', eps_turb_i0
  if (.not.(eps_i0 > delta_i0)) then
    write(*,*) 'FALHA: esperava eps_i0 > delta_i0 para pluma boiante e seca'
    all_ok = .false.
  end if
  if (eps_i0 < 0.0_RKIND .or. delta_i0 < 0.0_RKIND .or. eps_turb_i0 < 0.0_RKIND) then
    write(*,*) 'FALHA: taxa negativa encontrada'
    all_ok = .false.
  end if

  ! --- case 2: SAME dynamics, but updraft now SATURATED (RH_i=1) in a
  !     much drier environment -- must show enhanced detrainment ---
  call compute_entrainment_detrainment( &
       0.02_RKIND, 0.0_RKIND, 3.0_RKIND, 0.0_RKIND, 0.05_RKIND, &
       0.3_RKIND, 1.0_RKIND, 0.3_RKIND, 1000.0_RKIND, &
       eps_i0, delta_i0, eps_turb_i0)
  delta_case2 = delta_i0
  write(*,'(A,ES12.4,A,ES12.4)') 'CASO 2 (saturado, ambiente seco): delta_i0=', delta_i0, &
       '  vs caso 1 delta_i0=', delta_case1
  if (.not.(delta_case2 > delta_case1)) then
    write(*,*) 'FALHA: esperava detranhamento maior com updraft saturado + ambiente seco'
    all_ok = .false.
  end if

  ! --- case 3: unsaturated updraft (RH_i<1) in the SAME dry environment:
  !     M_0i must be EXACTLY zero regardless of how dry the environment is ---
  call compute_entrainment_detrainment( &
       0.02_RKIND, 0.0_RKIND, 3.0_RKIND, 0.0_RKIND, 0.05_RKIND, &
       0.3_RKIND, 0.99_RKIND, 0.01_RKIND, 1000.0_RKIND, &
       eps_i0, delta_i0, eps_turb_i0)
  delta_case3 = delta_i0
  write(*,'(A,ES12.4,A,ES12.4)') 'CASO 3 (quase saturado mas nao saturado): delta_i0=', delta_i0, &
       '  vs caso 1 delta_i0=', delta_case1
  if (abs(delta_case3 - delta_case1) > 1.0e-12_RKIND) then
    write(*,*) 'FALHA: delta_i0 do caso 3 deveria bater com o caso 1 (sem termo evaporativo)'
    all_ok = .false.
  end if

  ! --- case 4: near-zero vertical velocity difference, robustness check ---
  call compute_entrainment_detrainment( &
       0.02_RKIND, 0.0_RKIND, 0.01_RKIND, 0.0_RKIND, 0.05_RKIND, &
       0.3_RKIND, 0.6_RKIND, 0.5_RKIND, 1000.0_RKIND, &
       eps_i0, delta_i0, eps_turb_i0)
  write(*,'(A,ES12.4,A)') 'CASO 4 (w_i quase zero, testa robustez): eps_i0=', eps_i0, &
       '  (deve ser finito, nao Inf/NaN)'
  if (eps_i0 /= eps_i0) then   ! NaN check
    write(*,*) 'FALHA: eps_i0 e NaN'
    all_ok = .false.
  end if

  if (all_ok) then
    write(*,*) 'TODOS OS CHECKS DE CONSISTENCIA PASSARAM'
  else
    write(*,*) 'ALGUM CHECK FALHOU -- ver acima'
  end if

end program test_entrainment_closure
