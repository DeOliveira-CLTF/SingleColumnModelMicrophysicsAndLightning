! plume_microphysics_thompson.F90
!
! Wraps mp_gt_driver (module_mp_thompson.F, extracted unmodified from
! MPAS-Model -- same file mp_offline_driver.F90 already uses and
! validates) for use as the PLUME's own microphysics, one layer at a
! time, replacing plume_condensation.F90 + plume_precipitation.F90
! together in the coupled plume test.
!
! Architectural note (why this needed its own module rather than a
! one-line swap -- see the conversation this was built in): mp_gt_driver
! is written for a fixed EULERIAN column, with kts:kte spanning many
! levels and its own internal sedimentation solver moving mass from
! level k to k-1 WITHIN one call. The plume model here is LAGRANGIAN
! (a single parcel's trajectory, growing one new level per step, not a
! fixed column). The chosen compromise: call mp_gt_driver with a
! ONE-LEVEL "column" (kts=kte=kms=kme=1) at every plume step, with
! dt_in = dz/w (the time-equivalent of the parcel crossing this one
! layer) -- conceptually "given this parcel's state entering the
! layer, and no air above it yet, what does one microphysics timestep
! do to it, including whatever falls out". This mirrors exactly how
! plume_precipitation.F90's own dt_equiv=dz/w device already worked;
! the difference is that Thompson's own condensation/freezing/riming/
! autoconversion/sedimentation now do this physics, in place of this
! project's simplified saturation-adjustment + Kessler.
!
! Consequence: because Thompson handles the phase change itself
! (unlike plume_condensation.F90's theta_l trick), entrainment dilution
! ahead of this call must mix ordinary theta and qv, qc, qr, qi, qs, qg
! SEPARATELY (each simply relaxing toward the environment's own value,
! typically zero for the hydrometeors) -- not theta_l/qt. That change
! belongs in the coupled driver, not here; this module only wraps the
! microphysics call itself.
!
! New code of ours -- not extracted from MPAS. The extracted physics
! it calls (module_mp_thompson.F) is unmodified.

! BUG FOUND AND FIXED (important): mp_gt_driver loops over columns with
! i_end = MIN(ite, ide-1) and j_end = MIN(jte, jde-1). An earlier
! version of this wrapper (and mp_offline_driver.F90 itself) passed
! ide=jde=1, making i_end=j_end=0 -- the column loop NEVER EXECUTED and
! Thompson silently did nothing. ide/jde must be at least ite+1.
!
module plume_microphysics_thompson

  use mpas_kind_types, only : RKIND

  implicit none
  private

  public :: init_plume_thompson, step_plume_thompson, plume_thompson_assume_initialized
  integer(8), public, save :: n_plume_mp_calls = 0
  ! diagnostic split of the cost: calls where the parcel already holds condensate or is supersaturated
  ! (Thompson does real work) versus calls that hit the no_micro early exit
  integer(8), public, save :: n_mp_active = 0, n_mp_idle = 0
  real(8),    public, save :: t_mp_active = 0.d0, t_mp_idle = 0.d0

  logical, save :: initialized = .false.
  integer, save :: itimestep_saved = 0

  integer, parameter :: ids=1,ide=2, jds=1,jde=2, kds=1,kde=2
  integer, parameter :: ims=1,ime=1, jms=1,jme=1, kms=1,kme=1
  integer, parameter :: its=1,ite=1, jts=1,jte=1, kts=1,kte=1

contains

  !===================================================================
  subroutine init_plume_thompson()
  !===================================================================
    use module_mp_thompson, only : thompson_init
    if (.not. initialized) then
      call thompson_init(.false.)
      initialized = .true.
      itimestep_saved = 0
    end if
  end subroutine init_plume_thompson

  subroutine plume_thompson_assume_initialized()
    ! the host has already called thompson_init (grid microphysics = Thompson): do not re-initialize the tables
    initialized = .true.
  end subroutine plume_thompson_assume_initialized

  !===================================================================
  subroutine step_plume_thompson(theta, p, exner, w, dz, qv, qc, qr, qi, qs, qg, &
       &                         precip_removed, fall3)
  !===================================================================
    ! theta, p, exner, w, dz : in/out and in as noted below -- ordinary
    !            (dry) potential temperature (K), pressure (Pa), Exner
    !            function, plume vertical velocity (m/s, IN only), and
    !            this layer's thickness (m, IN only)
    ! qv,qc,qr,qi,qs,qg : hydrometeor mixing ratios (kg/kg), updated
    !            in place by one Thompson microphysics step
    ! precip_removed : accumulated surface-precip increment this step
    !            (mm liquid-water equivalent, RAINNCV's own units) --
    !            a diagnostic for the caller, analogous to
    !            plume_precipitation.F90's own output of the same name
    !            (see note at the end of this subroutine on units)
    use module_mp_thompson, only : mp_gt_driver, RSLF

    real(kind=RKIND), intent(inout) :: theta, qv, qc, qr, qi, qs, qg
    real(kind=RKIND), intent(in) :: p, exner, w, dz
    real(kind=RKIND), intent(out) :: precip_removed
    ! fall3 (optional): precipitation that left this layer, mm = kg/m2, by phase: (1) liquid rain, (2) snow + cloud ice,
    ! (3) graupel. Thompson's RAINNCV is the TOTAL (rain+snow+graupel+ice), SNOWNCV = snow+ice, GRAUPELNCV = graupel.
    real(kind=RKIND), intent(out), optional :: fall3(3)

    real(kind=RKIND) :: qv3(ims:ime,kms:kme,jms:jme), qc3(ims:ime,kms:kme,jms:jme)
    real(kind=RKIND) :: qr3(ims:ime,kms:kme,jms:jme), qi3(ims:ime,kms:kme,jms:jme)
    real(kind=RKIND) :: qs3(ims:ime,kms:kme,jms:jme), qg3(ims:ime,kms:kme,jms:jme)
    real(kind=RKIND) :: ni3(ims:ime,kms:kme,jms:jme), nr3(ims:ime,kms:kme,jms:jme)
    real(kind=RKIND) :: th3(ims:ime,kms:kme,jms:jme), pii3(ims:ime,kms:kme,jms:jme)
    real(kind=RKIND) :: p3(ims:ime,kms:kme,jms:jme), w3(ims:ime,kms:kme,jms:jme)
    real(kind=RKIND) :: dz3(ims:ime,kms:kme,jms:jme)
    real(kind=RKIND) :: re_cloud(ims:ime,kms:kme,jms:jme), re_ice(ims:ime,kms:kme,jms:jme)
    real(kind=RKIND) :: re_snow(ims:ime,kms:kme,jms:jme)
    real(kind=RKIND) :: rainnc(ims:ime,jms:jme), rainncv(ims:ime,jms:jme)
    real(kind=RKIND) :: snownc(ims:ime,jms:jme), snowncv(ims:ime,jms:jme)
    real(kind=RKIND) :: graupelnc(ims:ime,jms:jme), graupelncv(ims:ime,jms:jme)
    real(kind=RKIND) :: sr(ims:ime,jms:jme)
    real(kind=RKIND) :: rainprod(ims:ime,kms:kme,jms:jme), evapprod(ims:ime,kms:kme,jms:jme)
    real(kind=RKIND) :: ntc(ims:ime,jms:jme), muc(ims:ime,jms:jme)
    integer :: has_reqc, has_reqi, has_reqs
    real(kind=RKIND) :: dt_in
    logical :: active_call
    integer(8) :: tc0, tc1, tc_rate

    if (.not. initialized) call init_plume_thompson()

    ! No seeding needed: above 0 C Thompson's own ssati uses liquid
    ! saturation (qvsi=qvs), so its no_micro shortcut clears itself for a
    ! supersaturated warm parcel. (An earlier version of this wrapper
    ! seeded qc on the mistaken belief that no_micro blocked this case;
    ! the real problem was the ide/jde bug described in the header.)

    qv3(1,1,1) = qv; qc3(1,1,1) = qc; qr3(1,1,1) = qr
    qi3(1,1,1) = qi; qs3(1,1,1) = qs; qg3(1,1,1) = qg
    ! number concentrations start at zero, exactly as in
    ! mp_offline_driver.F90 -- Thompson diagnoses its own the first
    ! time a species with ni/nr=0 is found alongside qi/qr>0
    ni3 = 0.0_RKIND; nr3 = 0.0_RKIND
    th3(1,1,1) = theta
    pii3(1,1,1) = exner
    p3(1,1,1) = p
    w3(1,1,1) = w
    dz3(1,1,1) = dz
    ntc = 300.0e6_RKIND   ! continental cloud droplet number, same value
    muc = min(15.0_RKIND, 1000.0e6_RKIND/ntc(1,1) + 2.0_RKIND)  ! as mp_offline_driver.F90
    has_reqc = 0; has_reqi = 0; has_reqs = 0
    rainnc = 0.0_RKIND; rainncv = 0.0_RKIND
    snownc = 0.0_RKIND; snowncv = 0.0_RKIND
    graupelnc = 0.0_RKIND; graupelncv = 0.0_RKIND
    sr = 0.0_RKIND

    ! same dt_equiv convention as plume_precipitation.F90; capped at 300 s so a plume that is nearly stalled
    ! (w -> 0) does not ask Thompson for an hours-long single step (no effect for w > dz/300)
    dt_in = min(dz / max(w, 1.0e-3_RKIND), 300.0_RKIND)
    n_plume_mp_calls = n_plume_mp_calls + 1
    itimestep_saved = itimestep_saved + 1

    active_call = (qc+qr+qi+qs+qg > 0.0_RKIND) .or. (qv > RSLF(p, theta*exner))
    call system_clock(tc0,tc_rate)
    call mp_gt_driver( &
         qv=qv3, qc=qc3, qr=qr3, qi=qi3, qs=qs3, qg=qg3, ni=ni3, nr=nr3, &
         th=th3, pii=pii3, p=p3, w=w3, dz=dz3, &
         dt_in=dt_in, itimestep=itimestep_saved, &
         RAINNC=rainnc, RAINNCV=rainncv, SNOWNC=snownc, SNOWNCV=snowncv, &
         GRAUPELNC=graupelnc, GRAUPELNCV=graupelncv, SR=sr, &
         rainprod=rainprod, evapprod=evapprod, &
         re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
         has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs, &
         ntc=ntc, muc=muc, &
         ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde, &
         ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme, &
         its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte )

    call system_clock(tc1)
    if (active_call) then
      n_mp_active = n_mp_active + 1; t_mp_active = t_mp_active + real(tc1-tc0,8)/real(tc_rate,8)
    else
      n_mp_idle = n_mp_idle + 1;     t_mp_idle = t_mp_idle + real(tc1-tc0,8)/real(tc_rate,8)
    end if
    qv = qv3(1,1,1); qc = qc3(1,1,1); qr = qr3(1,1,1)
    qi = qi3(1,1,1); qs = qs3(1,1,1); qg = qg3(1,1,1)
    theta = th3(1,1,1)
    precip_removed = rainncv(1,1)
    if (present(fall3)) then
      fall3(2) = max(snowncv(1,1), 0.0_RKIND)
      fall3(3) = max(graupelncv(1,1), 0.0_RKIND)
      fall3(1) = max(rainncv(1,1) - snowncv(1,1) - graupelncv(1,1), 0.0_RKIND)
    end if
    ! precip_removed is RAINNCV as Thompson returns it (mm liquid-water
    ! equivalent over this dt_in) -- NOT the same units as
    ! plume_precipitation.F90's own kg/kg output. Converting it to a
    ! kg/kg-equivalent loss would need this layer's air density and
    ! depth; left in mm here since the caller only needs it for an
    ! accumulated-rain diagnostic, not to balance the qt budget --
    ! qv+qc+qr+qi+qs+qg already reflects the true water content
    ! directly here, unlike the Kessler module where qt needed manual
    ! correction after removing ql.

  end subroutine step_plume_thompson

end module plume_microphysics_thompson
