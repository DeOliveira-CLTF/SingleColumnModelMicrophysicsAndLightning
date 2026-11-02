! test_coupled_plume_thompson.F90
!
! Rewrite of test_coupled_plume.F90, replacing plume_condensation.F90 +
! plume_precipitation.F90 (the from-scratch, warm-rain-only physics)
! with the REAL Thompson microphysics (plume_microphysics_thompson.F90
! -> module_mp_thompson.F, extracted unmodified from MPAS-Model), one
! layer at a time. See plume_microphysics_thompson.F90's own header
! for the Eulerian-vs-Lagrangian architectural note this rests on.
!
! STATE CHANGE from the Kessler version: since Thompson does its own
! phase-change physics (condensation, freezing, riming, autoconversion,
! sedimentation, all at once), the plume no longer carries theta_l/qt
! as its conserved variables. Instead, it carries theta (ordinary dry
! potential temperature) and qv, qc, qr, qi, qs, qg SEPARATELY, each
! diluted by entrainment independently (mixing toward the
! environment's own value -- zero for every hydrometeor, since the
! environment is assumed clear). Buoyancy now includes loading from
! ALL FIVE hydrometeor species, not just liquid ql -- and, unlike the
! Kessler version, genuine ICE-PHASE physics (freezing, riming, graupel
! production) is now possible, since Thompson carries qi/qs/qg as real
! prognostic species rather than this project's own qc/ql substitute.
!
! Same steady-state bulk-plume scope as test_coupled_plume.F90 (see
! that file's header) -- not yet a genuinely time-dependent, multi-
! timestep column coupled to MYNN/the rest of this project's SCM.

program test_coupled_plume_thompson

  use mpas_kind_types, only : RKIND
  use plume_thermo, only : dilute_plume_scalar
  use plume_velocity, only : integrate_velocity
  use entrainment_closure, only : compute_entrainment_detrainment
  use plume_microphysics_thompson, only : init_plume_thompson, step_plume_thompson
  use module_mp_thompson, only : RSLF

  implicit none

  real(kind=RKIND), parameter :: Rd_local = 287.04_RKIND
  real(kind=RKIND), parameter :: cp_local = 1004.64_RKIND
  real(kind=RKIND), parameter :: p0_local = 1.0e5_RKIND
  real(kind=RKIND), parameter :: g_local  = 9.81_RKIND
  real(kind=RKIND), parameter :: a1_vel = 0.667_RKIND, a2_vel = 1.0_RKIND
  real(kind=RKIND), parameter :: e0_const = 0.3_RKIND

  integer, parameter :: nz = 400
  real(kind=RKIND), parameter :: ztop = 16000.0_RKIND

  real(kind=RKIND) :: z(nz), p(nz), theta_env(nz), qv_env(nz), T_env(nz), exner(nz)
  real(kind=RKIND) :: theta_u(nz), qv_u(nz), qc_u(nz), qr_u(nz), qi_u(nz), qs_u(nz), qg_u(nz)
  real(kind=RKIND) :: w_u(nz), a_u(nz), b_u(nz), RH_u(nz), RH_env(nz), M_flux(nz)
  real(kind=RKIND) :: eps_dyn, delta_i0, eps_turb, eps_total, H_i, dz_here, dz_uniform
  real(kind=RKIND) :: dz2(1), eps2(2), phi_env2(2), phi_u2(2), B2(2)
  real(kind=RKIND) :: precip_removed, precip_accum(nz), qcond_total
  real(kind=RKIND) :: thetav_u, thetav_env
  integer :: k, first_cloudy_k, stop_k
  real(kind=RKIND) :: p_surf

  dz_uniform = ztop/real(nz,kind=RKIND)
  call init_plume_thompson()

  ! ---- environment (same sounding as test_coupled_plume.F90: weak
  !      tropospheric lapse rate + stratospheric inversion above 11 km,
  !      exponential moisture decay, saturation-adjusted for self-
  !      consistency -- see that file's header for why each piece is
  !      there) ----
  p_surf = 95000.0_RKIND
  do k = 1, nz
    z(k) = (real(k,kind=RKIND)-0.5_RKIND)*dz_uniform
    if (z(k) < 11000.0_RKIND) then
      theta_env(k) = 300.0_RKIND + 0.003_RKIND*z(k)
    else
      theta_env(k) = (300.0_RKIND + 0.003_RKIND*11000.0_RKIND) + 0.010_RKIND*(z(k)-11000.0_RKIND)
    end if
    qv_env(k) = max(0.014_RKIND * exp(-z(k)/2000.0_RKIND), 0.0005_RKIND)
    p(k) = p_surf * exp(-z(k)/8400.0_RKIND)
    exner(k) = (p(k)/p0_local)**(Rd_local/cp_local)
    T_env(k) = theta_env(k) * exner(k)
    ! cap environmental humidity at RSLF (Thompson's own saturation
    ! function, for consistency with everything else in this driver)
    ! -- the environment itself is assumed clear (no separate
    ! condensation physics), this only prevents the unphysical
    ! RH_env>100% seen in an earlier version of this test
    qv_env(k) = min(qv_env(k), RSLF(p(k), T_env(k)))
    RH_env(k) = qv_env(k) / RSLF(p(k), T_env(k))
  end do

  ! ---- launch conditions ----
  theta_u(1) = theta_env(1) + 2.0_RKIND
  qv_u(1)    = qv_env(1) + 0.001_RKIND
  qc_u(1) = 0.0_RKIND; qr_u(1) = 0.0_RKIND; qi_u(1) = 0.0_RKIND
  qs_u(1) = 0.0_RKIND; qg_u(1) = 0.0_RKIND
  w_u(1) = 1.0_RKIND
  a_u(1) = 0.05_RKIND

  thetav_u   = theta_u(1)   * (1.0_RKIND + 0.608_RKIND*qv_u(1))
  thetav_env = theta_env(1) * (1.0_RKIND + 0.608_RKIND*qv_env(1))
  b_u(1) = g_local * (thetav_u-thetav_env) / thetav_env
  RH_u(1) = qv_u(1) / RSLF(p(1), theta_u(1)*exner(1))
  M_flux(1) = a_u(1) * w_u(1)
  precip_accum(1) = 0.0_RKIND

  first_cloudy_k = -1
  stop_k = -1

  do k = 2, nz
    if (stop_k > 0) then
      w_u(k) = 0.0_RKIND; a_u(k) = 0.0_RKIND
      theta_u(k) = theta_env(k); qv_u(k) = qv_env(k)
      qc_u(k) = 0.0_RKIND; qr_u(k) = 0.0_RKIND; qi_u(k) = 0.0_RKIND
      qs_u(k) = 0.0_RKIND; qg_u(k) = 0.0_RKIND
      b_u(k) = 0.0_RKIND; RH_u(k) = RH_env(k); M_flux(k) = 0.0_RKIND
      precip_accum(k) = precip_accum(k-1)
      cycle
    end if

    dz_here = z(k) - z(k-1)
    H_i = max(z(k-1), 100.0_RKIND)

    call compute_entrainment_detrainment( &
         b_u(k-1), 0.0_RKIND, w_u(k-1), 0.0_RKIND, a_u(k-1), &
         e0_const, RH_u(k-1), RH_env(k-1), H_i, &
         eps_dyn, delta_i0, eps_turb)
    eps_total = eps_dyn + eps_turb

    dz2(1) = dz_here
    eps2 = (/ 0.0_RKIND, eps_total /)

    ! --- dilute theta, qv, and every hydrometeor species separately
    !     (hydrometeors always relax toward the clear environment's
    !     zero value) ---
    phi_env2 = (/ theta_env(k-1), theta_env(k) /)
    phi_u2(1) = theta_u(k-1)
    call dilute_plume_scalar(2, dz2, eps2, phi_env2, phi_u2); theta_u(k) = phi_u2(2)

    phi_env2 = (/ qv_env(k-1), qv_env(k) /)
    phi_u2(1) = qv_u(k-1)
    call dilute_plume_scalar(2, dz2, eps2, phi_env2, phi_u2); qv_u(k) = phi_u2(2)

    phi_env2 = (/ 0.0_RKIND, 0.0_RKIND /)
    phi_u2(1) = qc_u(k-1)
    call dilute_plume_scalar(2, dz2, eps2, phi_env2, phi_u2); qc_u(k) = phi_u2(2)
    phi_u2(1) = qr_u(k-1)
    call dilute_plume_scalar(2, dz2, eps2, phi_env2, phi_u2); qr_u(k) = phi_u2(2)
    phi_u2(1) = qi_u(k-1)
    call dilute_plume_scalar(2, dz2, eps2, phi_env2, phi_u2); qi_u(k) = phi_u2(2)
    phi_u2(1) = qs_u(k-1)
    call dilute_plume_scalar(2, dz2, eps2, phi_env2, phi_u2); qs_u(k) = phi_u2(2)
    phi_u2(1) = qg_u(k-1)
    call dilute_plume_scalar(2, dz2, eps2, phi_env2, phi_u2); qg_u(k) = phi_u2(2)

    ! --- Thompson microphysics: condensation, freezing, riming,
    !     autoconversion, sedimentation, all at once, over this layer ---
    call step_plume_thompson(theta_u(k), p(k), exner(k), w_u(k-1), dz_here, &
         &                   qv_u(k), qc_u(k), qr_u(k), qi_u(k), qs_u(k), qg_u(k), &
         &                   precip_removed)
    precip_accum(k) = precip_accum(k-1) + precip_removed

    qcond_total = qc_u(k) + qr_u(k) + qi_u(k) + qs_u(k) + qg_u(k)
    if (qcond_total > 1.0e-8_RKIND .and. first_cloudy_k < 0) first_cloudy_k = k
    RH_u(k) = qv_u(k) / RSLF(p(k), theta_u(k)*exner(k))

    ! --- moist buoyancy, loading from ALL hydrometeor species ---
    thetav_u   = theta_u(k)   * (1.0_RKIND + 0.608_RKIND*qv_u(k) - qcond_total)
    thetav_env = theta_env(k) * (1.0_RKIND + 0.608_RKIND*qv_env(k))
    b_u(k) = g_local * (thetav_u-thetav_env) / thetav_env

    B2 = (/ b_u(k-1), b_u(k) /)
    eps2 = (/ 0.0_RKIND, eps_total /)
    phi_u2(1) = w_u(k-1)
    call integrate_velocity(2, dz2, B2, eps2, a1_vel, a2_vel, phi_u2)
    w_u(k) = phi_u2(2)

    if (w_u(k) <= 0.0_RKIND) then
      w_u(k) = 0.0_RKIND; stop_k = k; a_u(k) = 0.0_RKIND; M_flux(k) = 0.0_RKIND
      cycle
    end if

    M_flux(k) = M_flux(k-1) * exp((eps_total-delta_i0)*dz_here)
    a_u(k) = M_flux(k) / w_u(k)

  end do

  ! ---- report ----
  write(*,'(A)') 'z(m)  a_u  w_u(m/s)  b_u(m/s2)  qc(g/kg)  qr  qi  qs  qg  RH_u  chuva(mm)'
  do k = 1, nz, 10
    write(*,'(F7.1,2X,ES9.2,2X,F7.3,2X,ES9.2,2X,F6.3,2X,F6.3,2X,F6.3,2X,F6.3,2X,F6.3,2X,F5.3,2X,F8.3)') &
         z(k), a_u(k), w_u(k), b_u(k), qc_u(k)*1000.0_RKIND, qr_u(k)*1000.0_RKIND, &
         qi_u(k)*1000.0_RKIND, qs_u(k)*1000.0_RKIND, qg_u(k)*1000.0_RKIND, RH_u(k), precip_accum(k)
  end do

  write(*,*)
  if (first_cloudy_k > 0) then
    write(*,'(A,F8.1,A)') 'Primeiro nivel com condensado: z = ', z(first_cloudy_k), ' m'
  else
    write(*,'(A)') 'Nenhuma condensacao ocorreu'
  end if
  if (stop_k > 0) then
    write(*,'(A,F8.1,A)') 'Pluma parou (w chegou a zero) em z = ', z(stop_k), ' m'
  else
    write(*,'(A,F8.1,A)') 'Pluma ainda subindo no topo do dominio, z = ', ztop, ' m'
  end if

end program test_coupled_plume_thompson
