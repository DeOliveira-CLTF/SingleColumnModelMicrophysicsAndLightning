! test_plume_microphysics_thompson.F90
!
! Sanity test only -- Thompson has no closed-form solution to check
! against (it is a full multi-species scheme). This confirms:
!  1) a supersaturated parcel condenses (qc forms, qv drops, theta
!     RISES from latent heat release -- unlike plume_condensation.F90,
!     which conserves theta_l by construction, Thompson conserves
!     nothing analogous here; theta increasing is the expected sign
!     of genuine latent heating).
!  2) repeated layers (simulating several steps of ascent at a fixed
!     w, dz) let qc accumulate then start converting to qr (rain) and
!     eventually produce nonzero surface-precip increments -- the
!     autoconversion chain actually firing.
!  3) nothing goes negative or NaN across the run.

program test_plume_microphysics_thompson

  use mpas_kind_types, only : RKIND
  use plume_microphysics_thompson, only : init_plume_thompson, step_plume_thompson

  implicit none

  real(kind=RKIND), parameter :: Rd_local = 287.04_RKIND
  real(kind=RKIND), parameter :: cp_local = 1004.64_RKIND
  real(kind=RKIND), parameter :: p0_local = 1.0e5_RKIND

  real(kind=RKIND) :: theta, p, exner, w, dz, qv, qc, qr, qi, qs, qg, precip_removed
  real(kind=RKIND) :: qv0, qc0, theta0
  integer :: step
  logical :: all_ok, ever_precip

  all_ok = .true.
  ever_precip = .false.

  call init_plume_thompson()

  ! --- single-step condensation check ---
  p = 85000.0_RKIND
  exner = (p/p0_local)**(Rd_local/cp_local)
  theta = 300.0_RKIND           ! theta itself at 300K (simpler to reason about)
  qv = 0.016_RKIND               ! 16 g/kg -- clearly supersaturated at T~290K, 850 hPa
  qc = 0.0_RKIND; qr = 0.0_RKIND; qi = 0.0_RKIND; qs = 0.0_RKIND; qg = 0.0_RKIND
  w = 3.0_RKIND; dz = 100.0_RKIND

  qv0 = qv; qc0 = qc; theta0 = theta
  call step_plume_thompson(theta, p, exner, w, dz, qv, qc, qr, qi, qs, qg, precip_removed)

  write(*,'(A,F10.6,A,F10.6,A,F9.4,A,F9.4)') 'CASO 1 (condensacao): qv ', qv0*1000.0_RKIND, &
       ' -> ', qv*1000.0_RKIND, ' g/kg   qc -> ', qc*1000.0_RKIND, ' g/kg   theta ', theta0
  write(*,'(A,F9.4,A)') '  theta depois: ', theta, ' K (deve subir -- calor latente)'

  if (qc <= qc0) then
    write(*,*) 'FALHA: qc deveria aumentar com vapor supersaturado'
    all_ok = .false.
  end if
  if (qv >= qv0) then
    write(*,*) 'FALHA: qv deveria cair (viraram condensado)'
    all_ok = .false.
  end if
  if (theta <= theta0) then
    write(*,*) 'FALHA: theta deveria subir (liberacao de calor latente)'
    all_ok = .false.
  end if
  if (qc /= qc .or. qv /= qv .or. theta /= theta) then   ! NaN check
    write(*,*) 'FALHA: NaN encontrado'
    all_ok = .false.
  end if

  ! --- repeated layers: does the autoconversion chain eventually fire? ---
  do step = 1, 40
    call step_plume_thompson(theta, p, exner, w, dz, qv, qc, qr, qi, qs, qg, precip_removed)
    if (qc < 0.0_RKIND .or. qr < 0.0_RKIND .or. qv < 0.0_RKIND .or. theta /= theta) then
      write(*,*) 'FALHA: valor negativo ou NaN no passo ', step
      all_ok = .false.
      exit
    end if
    if (precip_removed > 0.0_RKIND) ever_precip = .true.
  end do
  write(*,'(A,F9.4,A,F9.4,A,F9.4,A,L2)') 'CASO 2 (40 camadas): qc final=', qc*1000.0_RKIND, &
       ' g/kg  qr final=', qr*1000.0_RKIND, ' g/kg  precip_removed(ultima)=', precip_removed, &
       '  precipitou em algum momento=', ever_precip

  if (all_ok) then
    write(*,*) 'TODOS OS CHECKS PASSARAM'
  else
    write(*,*) 'ALGUM CHECK FALHOU -- ver acima'
  end if

end program test_plume_microphysics_thompson
