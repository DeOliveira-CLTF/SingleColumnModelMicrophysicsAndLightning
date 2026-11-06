! test_coupled_plume.F90
!
! STAGE "coupling": the first genuinely self-consistent integration of
! all four modules built so far (area/mass-flux, dilution, velocity,
! the real entrainment/detrainment closure, and condensation),
! replacing the imposed/synthetic eps, delta, B used in each stage's
! own isolated validation test.
!
! SCOPE, stated plainly: this is a STEADY-STATE bulk-plume rise
! calculation (the classic entraining-plume family: prescribe an
! environmental sounding, launch a plume at the surface, march it
! upward one level at a time, each level's eps/delta/B computed from
! the plume's OWN state at the level below -- not yet a genuinely
! time-dependent, multi-timestep column coupled to MYNN/Thompson).
! Consistent with that scope, the area/mass-flux update below uses the
! STEADY, height-based specialization of the Stage-1 area equation
! (d(a*w)/dz = a*w*(eps-delta), integrated exactly per layer as an
! exponential update of the mass flux M=a*w) rather than calling
! prognostic_area.F90's own time-dependent finite-volume solver, which
! is built for a different problem (a full column, evolving over many
! timesteps) not yet exercised here. dilute_plume_scalar() and
! integrate_velocity() ARE called directly, each on a 2-level slice
! per step -- both were designed for exactly this kind of per-layer
! stepping, so no adaptation was needed there.
!
! Environment: a simple, weakly stable, moistening-with-height sounding
! (not from any real case) -- this run exists to check that the
! coupled system BEHAVES sensibly (condensation switching on where
! expected, the entrainment closure responding to it, the plume
! stopping somewhere reasonable), not to reproduce any specific
! observed profile.

program test_coupled_plume

  use mpas_kind_types, only : RKIND
  use plume_thermo, only : dilute_plume_scalar
  use plume_velocity, only : integrate_velocity
  use plume_condensation, only : saturation_adjustment, plume_buoyancy_moist, qsat_bolton
  use entrainment_closure, only : compute_entrainment_detrainment
  use plume_precipitation, only : apply_precipitation

  implicit none

  real(kind=RKIND), parameter :: Rd_local = 287.04_RKIND
  real(kind=RKIND), parameter :: cp_local = 1004.64_RKIND
  real(kind=RKIND), parameter :: p0_local = 1.0e5_RKIND
  real(kind=RKIND), parameter :: g_local  = 9.81_RKIND
  real(kind=RKIND), parameter :: a1_vel = 0.667_RKIND, a2_vel = 1.0_RKIND  ! Stage 2 velocity coefficients
  real(kind=RKIND), parameter :: e0_const = 0.3_RKIND   ! m2/s2, prescribed TKE (see header note)

  integer, parameter :: nz = 400
  real(kind=RKIND), parameter :: ztop = 16000.0_RKIND

  real(kind=RKIND) :: z(nz), p(nz), theta_env(nz), qv_env(nz), T_env(nz), exner(nz)
  real(kind=RKIND) :: theta_l_u(nz), qt_u(nz), w_u(nz), a_u(nz), b_u(nz)
  real(kind=RKIND) :: T_u(nz), qv_u(nz), ql_u(nz), theta_u(nz)
  real(kind=RKIND) :: RH_u(nz), RH_env(nz), M_flux(nz)
  real(kind=RKIND) :: eps_dyn, delta_i0, eps_turb, eps_total, H_i, dz_here, dz_uniform
  real(kind=RKIND) :: dz2(1), eps2(2), phi_env2(2), phi_u2(2), B2(2)
  real(kind=RKIND) :: precip_removed, precip_accum(nz)
  integer :: k, first_cloudy_k, stop_k
  real(kind=RKIND) :: p_surf

  dz_uniform = ztop/real(nz,kind=RKIND)

  ! ---- build the environment (simple, not a real sounding -- see header) ----
  ! IMPORTANT FIX: the environment must also pass through saturation
  ! adjustment. An earlier version of this test computed RH_env directly
  ! from a prescribed qv_env(z) with no saturation check, and that
  ! synthetic sounding turned out to reach RH_env > 130% above ~1500 m --
  ! a physically impossible "environment", since nothing stopped it from
  ! going supersaturated. Left uncorrected, that inflates the plume's
  ! buoyancy once it entrains such air, and could have made the "deep
  ! convection" result below spurious rather than genuine. Reusing
  ! saturation_adjustment() here (treating the prescribed qv_env(z) as
  ! total water) keeps the environment itself physically consistent
  ! (RH_env<=1 by construction) at no extra modeling cost.
  p_surf = 95000.0_RKIND
  do k = 1, nz
    z(k) = (real(k,kind=RKIND)-0.5_RKIND)*dz_uniform
    ! troposphere (weak lapse rate) capped by a stratospheric inversion
    ! above 11 km -- the earlier version of this test had NO tropopause
    ! at all (theta_env just kept increasing at the same weak tropospheric
    ! rate all the way up), so there was no physical mechanism left in
    ! the environment that could ever stop a sufficiently buoyant plume
    ! within a plausible domain height; this is the fix.
    if (z(k) < 11000.0_RKIND) then
      theta_env(k) = 300.0_RKIND + 0.003_RKIND*z(k)
    else
      theta_env(k) = (300.0_RKIND + 0.003_RKIND*11000.0_RKIND) + 0.010_RKIND*(z(k)-11000.0_RKIND)
    end if
    ! exponential humidity decay (scale height 2000 m) -- the earlier
    ! linear profile decayed too slowly with height, leaving RH_env
    ! pinned at exactly 100% (a permanently overcast "environment")
    ! above ~1800 m once the saturation-adjustment fix above was
    ! applied; that left no genuine dry/moist contrast for the
    ! evaporative detrainment term (M_0i, Eq. 38) to ever respond to,
    ! since it needs a REALLY drier environment than the (saturated)
    ! plume to activate. A more realistic free-tropospheric moisture
    ! scale height restores that contrast.
    qv_env(k)    = max(0.014_RKIND * exp(-z(k)/2000.0_RKIND), 0.0005_RKIND)
    ! hydrostatic-ish pressure profile (dry, constant scale height ~8.4km)
    p(k) = p_surf * exp(-z(k)/8400.0_RKIND)
    exner(k) = (p(k)/p0_local)**(Rd_local/cp_local)
    block
      real(kind=RKIND) :: ql_env_dummy, theta_env_dummy
      call saturation_adjustment(theta_env(k), qv_env(k), p(k), &
           &                     T_env(k), qv_env(k), ql_env_dummy, theta_env_dummy)
      ! qv_env(k) is overwritten here with its saturation-adjusted value
      ! (equal to the original if unsaturated, capped at qsat otherwise);
      ! any ql_env_dummy>0 (environmental cloud/fog) is discarded for
      ! this test -- a further simplification flagged, not solved, here
    end block
    RH_env(k) = qv_env(k) / qsat_bolton(T_env(k), p(k))
  end do

  ! ---- launch conditions at the surface ----
  theta_l_u(1) = theta_env(1) + 2.0_RKIND   ! 2K warmer than environment: a moderate thermal
  qt_u(1)      = qv_env(1) + 0.001_RKIND    ! 1 g/kg moister than environment
  w_u(1)       = 1.0_RKIND
  a_u(1)       = 0.05_RKIND

  call saturation_adjustment(theta_l_u(1), qt_u(1), p(1), T_u(1), qv_u(1), ql_u(1), theta_u(1))
  b_u(1) = plume_buoyancy_moist(theta_u(1), qv_u(1), ql_u(1), theta_env(1), qv_env(1), 0.0_RKIND, g_local)
  RH_u(1) = qv_u(1) / qsat_bolton(T_u(1), p(1))
  M_flux(1) = a_u(1) * w_u(1)
  precip_accum(1) = 0.0_RKIND

  first_cloudy_k = -1
  stop_k = -1

  ! ---- march the plume upward, one level at a time ----
  do k = 2, nz
    if (stop_k > 0) then
      ! plume already stopped below: everything above is quiescent
      w_u(k) = 0.0_RKIND; a_u(k) = 0.0_RKIND; theta_l_u(k) = theta_env(k); qt_u(k) = qv_env(k)
      T_u(k) = T_env(k); qv_u(k) = qv_env(k); ql_u(k) = 0.0_RKIND; theta_u(k) = theta_env(k)
      b_u(k) = 0.0_RKIND; RH_u(k) = RH_env(k); M_flux(k) = 0.0_RKIND
      precip_accum(k) = precip_accum(k-1)
      cycle
    end if

    dz_here = z(k) - z(k-1)
    H_i = max(z(k-1), 100.0_RKIND)   ! updraft height reached so far (Table 2 convention)

    ! entrainment/detrainment from the plume's OWN state at level k-1
    call compute_entrainment_detrainment( &
         b_u(k-1), 0.0_RKIND, w_u(k-1), 0.0_RKIND, a_u(k-1), &
         e0_const, RH_u(k-1), RH_env(k-1), H_i, &
         eps_dyn, delta_i0, eps_turb)
    eps_total = eps_dyn + eps_turb   ! see header note: turbulent entrainment adds to
                                     ! total mixing, not treated as separate detrainment

    ! --- dilute theta_l, qt over this one layer (2-level slice call) ---
    dz2(1) = dz_here
    eps2 = (/ 0.0_RKIND, eps_total /)   ! eps2(1) unused by dilute_plume_scalar (see its own header)

    phi_env2 = (/ theta_env(k-1), theta_env(k) /)
    phi_u2(1) = theta_l_u(k-1)
    call dilute_plume_scalar(2, dz2, eps2, phi_env2, phi_u2)
    theta_l_u(k) = phi_u2(2)

    phi_env2 = (/ qv_env(k-1), qv_env(k) /)
    phi_u2(1) = qt_u(k-1)
    call dilute_plume_scalar(2, dz2, eps2, phi_env2, phi_u2)
    qt_u(k) = phi_u2(2)

    ! --- condensation: where does this diluted state actually sit? ---
    call saturation_adjustment(theta_l_u(k), qt_u(k), p(k), T_u(k), qv_u(k), ql_u(k), theta_u(k))
    if (ql_u(k) > 1.0e-8_RKIND .and. first_cloudy_k < 0) first_cloudy_k = k

    ! --- precipitation: remove condensate above the Kessler threshold,
    !     using w_u(k-1) as the layer's velocity scale (w_u(k) itself
    !     is not known yet at this point in the loop -- same
    !     "state at k-1 sets this layer's physics" convention already
    !     used for eps/delta above). This can change theta_l_u(k) (see
    !     plume_precipitation.F90's header), so it must run BEFORE the
    !     buoyancy below, which depends on ql_u(k) directly, and before
    !     theta_l_u(k) is used by the NEXT layer's dilution step. ---
    call apply_precipitation(T_u(k), p(k), dz_here, w_u(k-1), &
         &                   qt_u(k), ql_u(k), theta_l_u(k), precip_removed)
    precip_accum(k) = precip_accum(k-1) + precip_removed
    ! qv_u(k)/theta_u(k) stay consistent with the (T_u(k), post-precip
    ! ql_u(k)) state directly, without a second saturation_adjustment call
    qv_u(k) = qt_u(k) - ql_u(k)

    RH_u(k) = qv_u(k) / qsat_bolton(T_u(k), p(k))

    ! --- moist buoyancy, now WITH condensate loading ---
    b_u(k) = plume_buoyancy_moist(theta_u(k), qv_u(k), ql_u(k), theta_env(k), qv_env(k), &
         &                        0.0_RKIND, g_local)

    ! --- velocity (2-level slice call) ---
    B2 = (/ b_u(k-1), b_u(k) /)
    eps2 = (/ 0.0_RKIND, eps_total /)
    phi_u2(1) = w_u(k-1)
    call integrate_velocity(2, dz2, B2, eps2, a1_vel, a2_vel, phi_u2)
    w_u(k) = phi_u2(2)

    if (w_u(k) <= 0.0_RKIND) then
      w_u(k) = 0.0_RKIND
      stop_k = k
      a_u(k) = 0.0_RKIND
      M_flux(k) = 0.0_RKIND
      cycle
    end if

    ! --- area, via the steady mass-flux specialization of Eq. 18
    !     (see header note): M=a*w obeys dM/dz = M*(eps_total-delta),
    !     integrated exactly per layer as an exponential update ---
    M_flux(k) = M_flux(k-1) * exp((eps_total-delta_i0)*dz_here)
    a_u(k) = M_flux(k) / w_u(k)

  end do

  ! ---- report ----
  write(*,'(A)') 'z(m)  a_u  w_u(m/s)  b_u(m/s2)  ql_u(g/kg)  RH_u  RH_env  chuva_acum(g/kg)'
  do k = 1, nz, 10
    write(*,'(F7.1,2X,ES9.2,2X,F7.3,2X,ES9.2,2X,F7.3,2X,F6.3,2X,F6.3,2X,F8.3)') &
         z(k), a_u(k), w_u(k), b_u(k), ql_u(k)*1000.0_RKIND, RH_u(k), RH_env(k), &
         precip_accum(k)*1000.0_RKIND
  end do

  write(*,*)
  if (first_cloudy_k > 0) then
    write(*,'(A,F8.1,A)') 'Primeiro nivel com condensado (LCL efetivo): z = ', z(first_cloudy_k), ' m'
  else
    write(*,'(A)') 'Nenhuma condensacao ocorreu -- pluma ficou seca o tempo todo'
  end if
  if (stop_k > 0) then
    write(*,'(A,F8.1,A)') 'Pluma parou (w chegou a zero) em z = ', z(stop_k), ' m'
  else
    write(*,'(A,F8.1,A)') 'Pluma ainda subindo no topo do dominio, z = ', ztop, ' m'
  end if

end program test_coupled_plume
