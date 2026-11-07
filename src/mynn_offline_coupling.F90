! mynn_offline_coupling.F90
!
! Wraps the two real, unmodified MPAS/MMM-physics schemes --
! sf_mynn_run (surface layer, sf_mynn.F90) and bl_mynn_run (PBL,
! bl_mynn.F90) -- behind a single mynn_coupling_step() call that
! mp_offline_driver.F90 invokes once per timestep, BEFORE
! mp_gt_driver. This module owns all of MYNN's persistent
! (call-to-call) state (TKE, mixing length, PBL height, similarity
! parameters) as SAVEd arrays, so the main driver only has to supply
! the current thermodynamic column and surface forcing and apply back
! the two tendencies it actually needs (rthblten, rqvblten).
!
! FORCING CONVENTION (this is the "Option 1" decision from README.md
! "Status"): sf_mynn_run does not accept a prescribed surface flux --
! it always DIAGNOSES hfx/qfx itself from the bulk-aerodynamic
! temperature/humidity difference between the surface (tsk) and the
! first model level, exactly as MPAS's own mpas_atmphys_driver_sfclayer.F
! calls it. This module therefore drives it with GOAMAZON's ts_forc
! (skin temperature) and z0 (roughness length), NOT with hfss/hfls --
! the resulting diagnosed flux will not exactly match GOAMAZON's own
! reported hfss/hfls (those come from independent
! observations/LES, this from a local bulk formula), by design. This
! mirrors a genuine, confirmed limitation of MPAS's own MYNN
! integration: MPAS's surface-layer driver has an SCM
! prescribed-flux mode ("scm_force_flux"), but it is wired up only for
! the classic Monin-Obukhov surface-layer scheme (sfclay/sfclayrev),
! never passed to module_sf_mynn.F -- so even a real MPAS run could not
! force MYNN's surface layer by prescribed flux either.
!
! WIND: GOAMAZON's ua_nud/va_nud are used directly as the column wind
! at every level (interpolated in time/height by goamazon_bridge), NOT
! integrated as a real prognostic momentum equation -- this driver has
! none. rublten/rvblten (MYNN's own wind tendencies from vertical
! mixing) are computed internally (mandatory dummy arguments) but
! discarded, not applied -- see README.md "Status".
!
! SCOPE LEFT OUT (see README.md "Status" for the full list): the EDMF
! mass-flux (shallow-cumulus) sub-scheme is left ON at MPAS's own
! Registry.xml default (bl_mynn_edmf=1), since its arrays are
! mandatory dummy arguments anyway; its own real, unmodified code runs,
! and its heat/moisture/cloud-water/cloud-ice tendencies (rthblten,
! rqvblten, rqcblten, rqiblten) ARE returned and applied by the caller
! -- see rqcblten_out/rqiblten_out below, the actual EDMF-to-
! microphysics coupling. Momentum tendencies (rublten/rvblten) are
! still not applied (no momentum equation in this driver -- wind is
! prescribed, not prognostic). Radiative heating (rthraten) is set to
! zero (no radiation scheme in this microphysics+PBL-only driver) --
! MYNN's optional radiation-aware entrainment enhancement therefore
! sees no radiative cooling at cloud top, a real simplification worth
! keeping in mind when interpreting any stratocumulus-like behavior.

module mynn_offline_coupling

  use mpas_kind_types, only : RKIND
  use bl_mynn, only : bl_mynn_run, bl_mynn_init
  use sf_mynn, only : sf_mynn_run
  use mynn_audit, only: ma_enabled,ma_main,ma_dt,ma_skip_shadow,ma_write,ma_boundary
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan

  implicit none
  private

  public :: mynn_coupling_init
  public :: mynn_coupling_step
  logical, save, public :: boundary_fix = .false. ! versioned sensitivity; retain historical reference by default
  logical, save, public :: surface_flux_memory=.false.
  integer, save, public :: surface_flux_mode=0 ! 0 bulk; 1 prescribed scalar-flux sensitivity, bulk momentum transfer

  ! ---- physical constants (match MPAS's own mpas_atmphys_constants.F
  !      and the standard WRF svp1/2/3/svpt0/ep1/ep2 Bolton constants) ----
  real(kind=RKIND), parameter :: cp_c     = 1004.5_RKIND
  real(kind=RKIND), parameter :: g_c      = 9.80616_RKIND
  real(kind=RKIND), parameter :: Rd_c     = 287.0_RKIND
  real(kind=RKIND), parameter :: Rv_c     = 461.6_RKIND
  real(kind=RKIND), parameter :: xlv_c    = 2.5e6_RKIND
  real(kind=RKIND), parameter :: karman_c = 0.4_RKIND
  real(kind=RKIND), parameter :: svp1_c   = 0.6112_RKIND
  real(kind=RKIND), parameter :: svp2_c   = 17.67_RKIND
  real(kind=RKIND), parameter :: svp3_c   = 29.65_RKIND
  real(kind=RKIND), parameter :: svpt0_c  = 273.15_RKIND
  real(kind=RKIND), parameter :: ep1_c    = Rv_c/Rd_c - 1.0_RKIND
  real(kind=RKIND), parameter :: ep2_c    = Rd_c/Rv_c

  ! MPAS Registry.xml defaults for the config_mynn_* namelist options
  ! (confirmed against src/core_atmosphere/Registry.xml):
  integer,  parameter :: bl_mynn_cloudpdf_c  = 2
  integer,  parameter :: bl_mynn_mixlength_c = 2
  integer,  parameter :: bl_mynn_stfunc_c    = 1
  integer,  parameter :: bl_mynn_topdown_c   = 0
  integer,  parameter :: bl_mynn_scaleaware_c= 1
  integer,  parameter :: bl_mynn_dheatopt_c  = 1
  integer,  save, public :: bl_mynn_edmf_c = 1  ! MPAS's own default: ON; 0 = EDMF off (control experiments only, set by the driver)
  integer,  parameter :: bl_mynn_edmf_dd_c   = 0
  integer,  parameter :: bl_mynn_edmf_mom_c  = 0
  integer,  parameter :: bl_mynn_edmf_tke_c  = 0
  integer,  parameter :: bl_mynn_mixscalars_c= 1
  integer,  parameter :: bl_mynn_cloudmix_c  = 1
  integer,  parameter :: bl_mynn_mixqt_c     = 0
  integer,  parameter :: bl_mynn_tkebudget_c = 0
  logical,  parameter :: bl_mynn_tkeadvect_c = .false.
  real(kind=RKIND), parameter :: bl_mynn_closure_c = 2.5_RKIND
  integer,  parameter :: icloud_bl_c = 1
  ! sf_mynn's own driver-side constants (mpas_atmphys_driver_sfclayer.F):
  integer,  parameter :: isfflx_c   = 1  ! surface heat/moisture flux ON
  integer,  parameter :: isftcflx_c = 0  ! Carlson-Boland/Charnock
  integer,  parameter :: iz0tlnd_c  = 0

  logical :: initialized = .false.
  integer :: nz_saved = 0

  ! ---- persistent (call-to-call) state, all shape (nz) for the one column ----
  real(kind=RKIND), allocatable, save :: qke(:), qke_adv(:), tsq(:), qsq(:), cov(:)
  real(kind=RKIND), allocatable, save :: el_pbl(:), sh(:), sm(:)
  real(kind=RKIND), allocatable, save :: cldfra_bl(:), qc_bl(:), qi_bl(:)
  real(kind=RKIND), save :: pblh_s = 0.0_RKIND
  integer,          save :: kpbl_s = 1
  real(kind=RKIND), save :: ust_s = 0.1_RKIND, rmol_s = 0.0_RKIND, mol_s = 0.0_RKIND
  real(kind=RKIND), save :: qsfc_s = 0.0_RKIND
  real(kind=RKIND), save :: hfx_s=0._RKIND,qfx_s=0._RKIND

contains

  !===================================================================
  subroutine mynn_coupling_init(nz)
  !===================================================================
    integer, intent(in) :: nz
    character(len=256) :: errmsg
    integer :: errflg
    real(kind=RKIND), parameter :: xlf_c = 3.34e5_RKIND   ! latent heat of fusion, J/kg
    real(kind=RKIND), parameter :: xls_c = xlv_c + xlf_c  ! latent heat of sublimation, J/kg
    real(kind=RKIND), parameter :: cpv_c = 1885.0_RKIND   ! specific heat of water vapor, J/kg/K
    real(kind=RKIND), parameter :: cliq_c = 4190.0_RKIND  ! specific heat of liquid water, J/kg/K
    real(kind=RKIND), parameter :: cice_c = 2106.0_RKIND  ! specific heat of ice, J/kg/K

    ! one-time initialization of bl_mynn_common's physical constants
    ! (karman, cp, xlv, ...) -- WITHOUT this call they stay at their
    ! Fortran default of 0.0, which silently turned every length-scale
    ! and TKE formula in mym_initialize/mym_length into a division by
    ! zero (this was the exact cause of the Infinity/NaN in qke found
    ! previously -- not a physics bug, a missing init call on my part,
    ! exactly analogous to Thompson's own thompson_init).
    call bl_mynn_init(con_cp=cp_c, con_cpv=cpv_c, con_cice=cice_c, con_cliq=cliq_c, &
         &            con_ep1=ep1_c, con_ep2=ep2_c, con_grav=g_c, con_karman=karman_c, &
         &            con_p0=100000.0_RKIND, con_rd=Rd_c, con_rv=Rv_c, &
         &            con_svp1=svp1_c, con_svp2=svp2_c, con_svp3=svp3_c, con_svpt0=svpt0_c, &
         &            con_xlf=xlf_c, con_xls=xls_c, con_xlv=xlv_c, errmsg=errmsg, errflg=errflg)
    if (errflg /= 0) then
      write(0,'(A)') 'mynn_coupling_init: bl_mynn_init failed: '//trim(errmsg)
      stop 1
    end if

    nz_saved = nz
    allocate(qke(nz), qke_adv(nz), tsq(nz), qsq(nz), cov(nz))
    allocate(el_pbl(nz), sh(nz), sm(nz))
    allocate(cldfra_bl(nz), qc_bl(nz), qi_bl(nz))
    qke = 0.0_RKIND; qke_adv = 0.0_RKIND; tsq = 0.0_RKIND; qsq = 0.0_RKIND; cov = 0.0_RKIND
    el_pbl = 0.0_RKIND; sh = 0.0_RKIND; sm = 0.0_RKIND
    cldfra_bl = 0.0_RKIND; qc_bl = 0.0_RKIND; qi_bl = 0.0_RKIND
    pblh_s = 100.0_RKIND  ! first-guess PBL depth, m -- MYNN re-diagnoses this every call
    kpbl_s = 2
    ust_s = 0.1_RKIND; rmol_s = 0.0_RKIND; mol_s = 0.0_RKIND
    hfx_s=0._RKIND;qfx_s=0._RKIND;qsfc_s=0._RKIND
    initialized = .true.
  end subroutine mynn_coupling_init

  !===================================================================
  subroutine mynn_coupling_step(nz, dt, z, dz, u_col, v_col, th, tt, p, exner, rho, &
       &                        qv, qc, qi, qs, ts, z0_in, ps_in, itimestep, &
       &                        rthblten, rqvblten, rqcblten_out, rqiblten_out, &
       &                        rqcblten_loc_out, rqiblten_loc_out, &
       &                        tke_out, mixing_length_out, exch_h_out, exch_m_out, &
       &                        cldfra_bl_out, qc_bl_out, qi_bl_out, &
       &                        mf_a_out, mf_w_out, pblh_out,hfss_in,hfls_in)
  !===================================================================
    integer, intent(in) :: nz, itimestep
    real(kind=RKIND), intent(in) :: dt, ts, z0_in, ps_in
    real(kind=RKIND), intent(in) :: z(nz), dz(nz), u_col(nz), v_col(nz)
    real(kind=RKIND), intent(in) :: th(nz), tt(nz), p(nz), exner(nz), rho(nz)
    real(kind=RKIND), intent(in) :: qv(nz), qc(nz), qi(nz), qs(nz)
    real(kind=RKIND), intent(in), optional :: hfss_in,hfls_in ! W/m2 upward, complete time forcing
    real(kind=RKIND), intent(out) :: rthblten(nz), rqvblten(nz)
    real(kind=RKIND), intent(out), optional :: rqcblten_out(nz), rqiblten_out(nz)

    ! ---- additional PBL/turbulence diagnostics, all from bl_mynn_run's
    !      own output arrays (qke2, el_pbl2, exch_h2/m2, cldfra_bl2,
    !      qc_bl2/qi_bl2, edmf_a2/edmf_w2) -- computed internally on
    !      every call already, just not returned to the caller before
    !      this extension.
    real(kind=RKIND), intent(out), optional :: tke_out(nz)           ! m2/s2 (=qke/2)
    real(kind=RKIND), intent(out), optional :: mixing_length_out(nz) ! m (el_pbl, master length scale)
    real(kind=RKIND), intent(out), optional :: exch_h_out(nz)        ! m2/s, eddy diffusivity, heat
    real(kind=RKIND), intent(out), optional :: exch_m_out(nz)        ! m2/s, eddy diffusivity, momentum
    real(kind=RKIND), intent(out), optional :: cldfra_bl_out(nz)     ! 0-1, PBL-diagnosed subgrid cloud fraction
    real(kind=RKIND), intent(out), optional :: qc_bl_out(nz)         ! kg/kg, PBL-diagnosed subgrid cloud water
    real(kind=RKIND), intent(out), optional :: qi_bl_out(nz)         ! kg/kg, PBL-diagnosed subgrid cloud ice
    real(kind=RKIND), intent(out), optional :: mf_a_out(nz)          ! 0-1, EDMF updraft area fraction
    real(kind=RKIND), intent(out), optional :: mf_w_out(nz)          ! kg m-2 s-1, rhoz * area-weighted plume velocity
    real(kind=RKIND), intent(out), optional :: pblh_out              ! m, boundary-layer height (scalar)
    ! LOCAL-ONLY condensate tendencies: what bl_mynn_run returns for qc/qi when its
    ! EDMF mass-flux component is switched OFF, from the SAME input state. Computed by
    ! a shadow call (see the ipass loop below). Full - local is the response of the
    ! complete MYNN call to EDMF, including cloud/mixing feedbacks.
    real(kind=RKIND), intent(out), optional :: rqcblten_loc_out(nz), rqiblten_loc_out(nz)
    integer :: ipass, npass
    logical :: edmf_flag
    real(kind=RKIND) :: rqcloc(nz), rqiloc(nz)

    ! ---- surface-layer (sf_mynn_run) scalars, its:ite = 1:1 ----
    real(kind=RKIND) :: mavail(1), pblh1(1), xland(1), tsk(1), psfcpa(1), qcg(1), snowh(1), dx1(1)
    real(kind=RKIND) :: u1d(1), v1d(1), u1d2(1), v1d2(1), qv1d(1), p1d(1), t1d(1), dz8w1d(1), dz2w1d(1), rho1d(1)
    real(kind=RKIND) :: rstoch1d(1)
    real(kind=RKIND) :: regime(1), hfx(1), qfx(1), lh(1), mol1(1), rmol1(1), qgh(1), qsfc1(1)
    real(kind=RKIND) :: znt(1), zol(1), ust1(1), cpm(1), chs2(1), cqs2(1), chs(1), ch(1)
    real(kind=RKIND) :: flhc(1), flqc(1), gz1oz0(1), wspd(1), br(1), psim(1), psih(1)
    real(kind=RKIND) :: u10(1), v10(1), th2(1), t2(1), q2(1), wstar(1), qstar(1)
    character(len=256) :: errmsg
    integer :: errflg

    ! ---- PBL (bl_mynn_run) column arrays, its:ite=1:1 -> shape (1,nz) ----
    real(kind=RKIND) :: u2(1,nz), v2(1,nz), w2(1,nz), th2d(1,nz), tt2d(1,nz), p2(1,nz), &
         &              exner2(1,nz), rho2(1,nz), rthraten2(1,nz)
    real(kind=RKIND) :: sqv2(1,nz), sqc2(1,nz), sqi2(1,nz), sqs2(1,nz)
    real(kind=RKIND) :: qnc2(1,nz), qni2(1,nz), qnwfa2(1,nz), qnifa2(1,nz), qnbca2(1,nz), qozone2(1,nz)
    real(kind=RKIND) :: pattern_spp2(1,nz)
    real(kind=RKIND) :: cldfra_bl2(1,nz), qc_bl2(1,nz), qi_bl2(1,nz)
    real(kind=RKIND) :: el_pbl2(1,nz), qke2(1,nz), qke_adv2(1,nz), cov2(1,nz), qsq2(1,nz), tsq2(1,nz), &
         &              sh2(1,nz), sm2(1,nz)
    real(kind=RKIND) :: rublten2(1,nz), rvblten2(1,nz), rthblten2(1,nz), rqvblten2(1,nz), &
         &              rqcblten2(1,nz), rqiblten2(1,nz), rqsblten2(1,nz), rqncblten2(1,nz), &
         &              rqniblten2(1,nz), rqnifablten2(1,nz), rqnwfablten2(1,nz), rqnbcablten2(1,nz), &
         &              rqozblten2(1,nz)
    real(kind=RKIND) :: edmf_a2(1,nz), edmf_w2(1,nz), edmf_qt2(1,nz), edmf_thl2(1,nz), edmf_ent2(1,nz), &
         &              edmf_qc2(1,nz), sub_thl2(1,nz), sub_sqv2(1,nz), det_thl2(1,nz), det_sqv2(1,nz)
    real(kind=RKIND) :: exch_h2(1,nz), exch_m2(1,nz)
    real(kind=RKIND) :: maxwidth1(1), maxmf1(1), ztop_plume1(1)
    integer :: kpbl1(1), ktop_plume1(1)
    real(kind=RKIND) :: pblh_bl1(1)

    ! dz on the (1,nz) shape bl_mynn_run expects
    real(kind=RKIND) :: dz_arr(1,nz)

    integer :: k
    real(RKIND) :: audit_profiles(nz,20),audit_surface(20),input_hfx,input_qfx,bulk_hfx,bulk_qfx,forcing_hfx,forcing_lh

    if (.not. initialized) call mynn_coupling_init(nz)

    !=================================================================
    ! 1) surface layer (sf_mynn_run): diagnoses ust/rmol/ch/hfx/qfx from
    !    the temperature/humidity contrast between tsk (=ts, GOAMAZON's
    !    ts_forc) and the first two model levels -- see this module's
    !    header for why flux is diagnosed here, not prescribed.
    !=================================================================
    mavail  = 1.0_RKIND     ! inherited wet-limit assumption; not independently validated for Argentina
    pblh1   = pblh_s
    xland   = 1.0_RKIND     ! land
    tsk     = ts
    psfcpa  = ps_in
    qcg     = 0.0_RKIND
    snowh   = 0.0_RKIND
    dx1     = 5000.0_RKIND  ! nominal single-column horizontal scale, m (only used by scale-aware options)

    u1d  = u_col(1); v1d  = v_col(1)
    u1d2 = u_col(min(2,nz)); v1d2 = v_col(min(2,nz))
    qv1d = qv(1); p1d = p(1); t1d = tt(1)
    dz8w1d = dz(1)
    dz2w1d = z(min(2,nz)) - z(1)
    rho1d  = rho(1)
    rstoch1d = 0.0_RKIND

    ! first-guess/previous-step values for the iterative similarity solve:
    ust1 = ust_s; rmol1 = rmol_s; mol1 = mol_s; znt = z0_in; qsfc1 = qsfc_s
    ! No land model updates qsfc in this offline driver. In the existing wet-limit
    ! boundary, request saturation at CURRENT Tskin/ps rather than latch first value.
    if(boundary_fix) qsfc1=0._RKIND
    regime = 0.0_RKIND; hfx = 0.0_RKIND; qfx = 0.0_RKIND; qgh = 0.0_RKIND
    ! The SF scheme reads incoming fluxes for Beljaars wstar BEFORE recomputing
    ! bulk fluxes. Retain the last actually applied flux only in the opt-in mode.
    if(surface_flux_memory) then
      hfx = hfx_s; qfx = qfx_s
    endif
    forcing_hfx=ieee_value(0._RKIND,ieee_quiet_nan);forcing_lh=forcing_hfx
    if(present(hfss_in)) forcing_hfx=hfss_in
    if(present(hfls_in)) forcing_lh=hfls_in
    if(surface_flux_mode==1) then
      if(.not.present(hfss_in).or..not.present(hfls_in)) error stop 'Prescribed surface flux needs both forcing fields'
      hfx=forcing_hfx;qfx=forcing_lh/xlv_c
    else if(surface_flux_mode/=0) then
      error stop 'Unsupported surface flux mode'
    endif
    input_hfx=hfx(1);input_qfx=qfx(1)
    cpm = cp_c; chs2 = 0.0_RKIND; cqs2 = 0.0_RKIND; chs = 0.0_RKIND; ch = 0.0_RKIND
    flhc = 0.0_RKIND; flqc = 0.0_RKIND; gz1oz0 = 0.0_RKIND; wspd = 0.0_RKIND; br = 0.0_RKIND
    psim = 0.0_RKIND; psih = 0.0_RKIND; zol = 0.0_RKIND

    call sf_mynn_run( &
         u1d=u1d, v1d=v1d, t1d=t1d, qv1d=qv1d, p1d=p1d, dz8w1d=dz8w1d, rho1d=rho1d, &
         u1d2=u1d2, v1d2=v1d2, dz2w1d=dz2w1d, cp=cp_c, g=g_c, rovcp=Rd_c/cp_c, r=Rd_c, xlv=xlv_c, &
         psfcpa=psfcpa, chs=chs, chs2=chs2, cqs2=cqs2, cpm=cpm, pblh=pblh1, rmol=rmol1, &
         znt=znt, ust=ust1, mavail=mavail, zol=zol, mol=mol1, regime=regime, psim=psim, &
         psih=psih, xland=xland, hfx=hfx, qfx=qfx, tsk=tsk, u10=u10, v10=v10, th2=th2, &
         t2=t2, q2=q2, flhc=flhc, flqc=flqc, snowh=snowh, qgh=qgh, qsfc=qsfc1, lh=lh, &
         gz1oz0=gz1oz0, wspd=wspd, br=br, isfflx=isfflx_c, dx=dx1, svp1=svp1_c, svp2=svp2_c, &
         svp3=svp3_c, svpt0=svpt0_c, ep1=ep1_c, ep2=ep2_c, karman=karman_c, ch=ch, qcg=qcg, &
         itimestep=itimestep, wstar=wstar, qstar=qstar, &
         spp_pbl=.false., rstoch1d=rstoch1d, isftcflx=isftcflx_c, iz0tlnd=iz0tlnd_c, &
         its=1, ite=1, errmsg=errmsg, errflg=errflg )

    bulk_hfx=hfx(1);bulk_qfx=qfx(1)
    if(surface_flux_mode==1) then
      ! Diagnostic boundary experiment only. Keep SF bulk momentum/stability;
      ! prescribe heat/vapor fluxes to PBL, which recomputes its own 1/L.
      ! Do not claim a fully flux-driven Monin-Obukhov surface-layer solution.
      hfx=forcing_hfx;qfx=forcing_lh/xlv_c;lh=forcing_lh
    endif
    hfx_s=hfx(1);qfx_s=qfx(1)
    if(ma_enabled) ma_boundary=[input_hfx,input_qfx,bulk_hfx,bulk_qfx,hfx(1),qfx(1),lh(1), &
      forcing_hfx,forcing_lh,wstar(1),ust1(1),rmol1(1),rho(1),pblh1(1),ts,qsfc1(1),mavail(1), &
      real(surface_flux_mode,RKIND),merge(1._RKIND,0._RKIND,surface_flux_memory),dt]
    ust_s = ust1(1); rmol_s = rmol1(1); mol_s = mol1(1); qsfc_s = qsfc1(1); pblh_s = pblh1(1)

    !=================================================================
    ! 2) PBL (bl_mynn_run): whole-column TKE-based vertical mixing,
    !    using the ust/rmol/ch just diagnosed as its surface boundary
    !    condition.
    !=================================================================
    ! Two passes when the local-only tendencies are requested: pass 1 = EDMF OFF (shadow,
    ! results kept only as rqcloc/rqiloc), pass 2 = the real call. bl_mynn_run's in/out
    ! arrays are all local copies re-initialised from the persistent module state below,
    ! and that persistent state is only updated AFTER the loop, so both passes see
    ! identical inputs.
    npass = 1
    if (present(rqcblten_loc_out) .or. present(rqiblten_loc_out)) npass = 2
    if(ma_skip_shadow) npass=1
    rqcloc=0.;rqiloc=0.
    do ipass = 3-npass, 2
    ma_main=(ipass==2)
    edmf_flag = (bl_mynn_edmf_c == 1)
    if (npass == 2 .and. ipass == 1) edmf_flag = .false.
    u2(1,:) = u_col; v2(1,:) = v_col; w2(1,:) = 0.0_RKIND
    dz_arr(1,:) = dz
    th2d(1,:) = th; tt2d(1,:) = tt; p2(1,:) = p; exner2(1,:) = exner; rho2(1,:) = rho
    rthraten2 = 0.0_RKIND   ! no radiation scheme in this driver -- see header note
    sqv2(1,:) = qv; sqc2(1,:) = qc; sqi2(1,:) = qi; sqs2(1,:) = qs
    if(boundary_fix) then
      ! MYNN expects specific vapor and condensate per moist-gas mass; SCM uses
      ! dry-air mixing ratios. Match MYNN's inverse q/(1-sqv) exactly.
      sqv2(1,:) = qv/(1._RKIND+qv)
      sqc2(1,:) = qc/(1._RKIND+qv)
      sqi2(1,:) = qi/(1._RKIND+qv)
      sqs2(1,:) = qs/(1._RKIND+qv)
    endif
    qnc2 = 0.0_RKIND; qni2 = 0.0_RKIND; qnwfa2 = 0.0_RKIND; qnifa2 = 0.0_RKIND; qnbca2 = 0.0_RKIND
    qozone2 = 0.0_RKIND; pattern_spp2 = 0.0_RKIND

    cldfra_bl2(1,:) = cldfra_bl; qc_bl2(1,:) = qc_bl; qi_bl2(1,:) = qi_bl
    el_pbl2(1,:) = el_pbl; qke2(1,:) = qke; qke_adv2(1,:) = qke_adv
    cov2(1,:) = cov; qsq2(1,:) = qsq; tsq2(1,:) = tsq; sh2(1,:) = sh; sm2(1,:) = sm

    rublten2 = 0.0_RKIND; rvblten2 = 0.0_RKIND; rthblten2 = 0.0_RKIND; rqvblten2 = 0.0_RKIND
    rqcblten2 = 0.0_RKIND; rqiblten2 = 0.0_RKIND; rqsblten2 = 0.0_RKIND
    rqncblten2 = 0.0_RKIND; rqniblten2 = 0.0_RKIND; rqnifablten2 = 0.0_RKIND
    rqnwfablten2 = 0.0_RKIND; rqnbcablten2 = 0.0_RKIND; rqozblten2 = 0.0_RKIND
    edmf_a2 = 0.0_RKIND; edmf_w2 = 0.0_RKIND; edmf_qt2 = 0.0_RKIND; edmf_thl2 = 0.0_RKIND
    edmf_ent2 = 0.0_RKIND; edmf_qc2 = 0.0_RKIND; sub_thl2 = 0.0_RKIND; sub_sqv2 = 0.0_RKIND
    det_thl2 = 0.0_RKIND; det_sqv2 = 0.0_RKIND
    exch_h2 = 0.0_RKIND; exch_m2 = 0.0_RKIND
    kpbl1 = kpbl_s; ktop_plume1 = 0; maxwidth1 = 0.0_RKIND; maxmf1 = 0.0_RKIND; ztop_plume1 = 0.0_RKIND
    pblh_bl1 = pblh_s

    call bl_mynn_run( &
         initflag=merge(1,0,itimestep<=1), restart=.false., cycling=.false., &
         delt=dt, dz=dz_arr, dx=dx1, znt=znt, u=u2, v=v2, w=w2, th=th2d, &
         sqv=sqv2, sqc=sqc2, sqi=sqi2, sqs=sqs2, qnc=qnc2, qni=qni2, qnwfa=qnwfa2, &
         qnifa=qnifa2, qnbca=qnbca2, qozone=qozone2, p=p2, exner=exner2, rho=rho2, &
         tt=tt2d, xland=xland, ts=tsk, qsfc=qsfc1, ps=psfcpa, ust=ust1, ch=ch, hfx=hfx, &
         qfx=qfx, rmol=rmol1, wspd=wspd, uoce=(/0.0_RKIND/), voce=(/0.0_RKIND/), &
         qke=qke2, qke_adv=qke_adv2, tsq=tsq2, qsq=qsq2, cov=cov2, sh=sh2, sm=sm2, &
         rublten=rublten2, rvblten=rvblten2, rthblten=rthblten2, rqvblten=rqvblten2, &
         rqcblten=rqcblten2, rqiblten=rqiblten2, rqsblten=rqsblten2, rqncblten=rqncblten2, &
         rqniblten=rqniblten2, rqnwfablten=rqnwfablten2, rqnifablten=rqnifablten2, &
         rqnbcablten=rqnbcablten2, rqozblten=rqozblten2, exch_h=exch_h2, exch_m=exch_m2, &
         pblh=pblh_bl1, kpbl=kpbl1, el_pbl=el_pbl2, &
         qc_bl=qc_bl2, qi_bl=qi_bl2, cldfra_bl=cldfra_bl2, icloud_bl=icloud_bl_c, &
         bl_mynn_tkeadvect=bl_mynn_tkeadvect_c, bl_mynn_tkebudget=merge(.true.,.false.,bl_mynn_tkebudget_c==1), &
         bl_mynn_cloudpdf=bl_mynn_cloudpdf_c, bl_mynn_mixlength=bl_mynn_mixlength_c, &
         bl_mynn_closure=bl_mynn_closure_c, bl_mynn_stfunc=bl_mynn_stfunc_c, &
         bl_mynn_topdown=merge(.true.,.false.,bl_mynn_topdown_c==1), &
         bl_mynn_edmf=edmf_flag, &
         bl_mynn_edmf_dd=merge(.true.,.false.,bl_mynn_edmf_dd_c==1), &
         bl_mynn_edmf_mom=merge(.true.,.false.,bl_mynn_edmf_mom_c==1), &
         bl_mynn_edmf_tke=merge(.true.,.false.,bl_mynn_edmf_tke_c==1), &
         bl_mynn_mixscalars=merge(.true.,.false.,bl_mynn_mixscalars_c==1), &
         bl_mynn_output=.true., &
         bl_mynn_cloudmix=merge(.true.,.false.,bl_mynn_cloudmix_c==1), &
         bl_mynn_mixqt=merge(.true.,.false.,bl_mynn_mixqt_c==1), &
         bl_mynn_scaleaware=merge(.true.,.false.,bl_mynn_scaleaware_c==1), &
         bl_mynn_dheatopt=merge(.true.,.false.,bl_mynn_dheatopt_c==1), &
         edmf_a=edmf_a2, edmf_w=edmf_w2, edmf_qt=edmf_qt2, edmf_thl=edmf_thl2, edmf_ent=edmf_ent2, &
         edmf_qc=edmf_qc2, sub_thl=sub_thl2, sub_sqv=sub_sqv2, det_thl=det_thl2, det_sqv=det_sqv2, &
         maxwidth=maxwidth1, maxmf=maxmf1, ztop_plume=ztop_plume1, ktop_plume=ktop_plume1, &
         spp_pbl=0, pattern_spp_pbl=pattern_spp2, rthraten=rthraten2, &
         flag_qc=.true., flag_qi=.true., flag_qs=.true., flag_qnc=.false., flag_qni=.false., &
         flag_qnwfa=.false., flag_qnifa=.false., flag_qnbca=.false., flag_qoz=.false., &
         its=1, ite=1, kts=1, kte=nz, kme=nz, errmsg=errmsg, errflg=errflg )
    if(boundary_fix) then
      ! Convert complete finite-step end states, not just multiply tendencies
      ! by a constant: vapor changes the common denominator during the step.
      if(any(1._RKIND-sqv2(1,:)-dt*rqvblten2(1,:)<=0._RKIND)) error stop 'Invalid MYNN end humidity'
      rqcblten2(1,:)=((sqc2(1,:)+dt*rqcblten2(1,:))/(1._RKIND-sqv2(1,:)-dt*rqvblten2(1,:))-qc)/dt
      rqiblten2(1,:)=((sqi2(1,:)+dt*rqiblten2(1,:))/(1._RKIND-sqv2(1,:)-dt*rqvblten2(1,:))-qi)/dt
      rqvblten2(1,:)=((sqv2(1,:)+dt*rqvblten2(1,:))/(1._RKIND-sqv2(1,:)-dt*rqvblten2(1,:))-qv)/dt
    endif
    if (npass == 2 .and. ipass == 1) then
      rqcloc = rqcblten2(1,:); rqiloc = rqiblten2(1,:)
    end if
    end do   ! ipass
    if (present(rqcblten_loc_out)) rqcblten_loc_out = rqcloc
    if (present(rqiblten_loc_out)) rqiblten_loc_out = rqiloc

    ! carry TKE/mixing-length/PBL-height state forward to the next call
    qke = qke2(1,:); qke_adv = qke_adv2(1,:); tsq = tsq2(1,:); qsq = qsq2(1,:); cov = cov2(1,:)
    el_pbl = el_pbl2(1,:); sh = sh2(1,:); sm = sm2(1,:)
    cldfra_bl = cldfra_bl2(1,:); qc_bl = qc_bl2(1,:); qi_bl = qi_bl2(1,:)
    pblh_s = pblh_bl1(1); kpbl_s = kpbl1(1)

    rthblten = rthblten2(1,:)
    rqvblten = rqvblten2(1,:)
    ! rublten2/rvblten2 (wind mixing tendencies) are still not returned
    ! (no momentum equation in this driver -- wind is prescribed from
    ! GOAMAZON's ua_nud/va_nud, so a wind tendency has nowhere to act).
    !
    ! rqcblten2/rqiblten2 are TOTAL local mixing plus EDMF response, not
    ! isolated detrainment. They may be nonzero when EDMF is inactive.
    ! Full minus shadow includes feedbacks of the complete MYNN call.
    ! rqcblten2/rqiblten2 (cloud-water/ice tendencies)
    ! ARE returned to the microphysics coupling. Local mixing remains available
    ! even when no plume is active; these fields are not plume-only sources.
    if (present(rqcblten_out)) rqcblten_out = rqcblten2(1,:)
    if (present(rqiblten_out)) rqiblten_out = rqiblten2(1,:)

    ! turbulence/mass-flux/cloud-fraction diagnostics requested by the
    ! caller (all already computed above; this just exposes them)
    if (present(tke_out))           tke_out           = 0.5_RKIND*qke2(1,:)
    if (present(mixing_length_out)) mixing_length_out = el_pbl2(1,:)
    if (present(exch_h_out))        exch_h_out        = exch_h2(1,:)
    if (present(exch_m_out))        exch_m_out        = exch_m2(1,:)
    if (present(cldfra_bl_out))     cldfra_bl_out     = cldfra_bl2(1,:)
    if (present(qc_bl_out))         qc_bl_out         = qc_bl2(1,:)
    if (present(qi_bl_out))         qi_bl_out         = qi_bl2(1,:)
    ! mass flux itself (kg/m2/s) is mf_a*mf_w; mf_w ALREADY contains rhoz.
    ! Characteristic plume speed is mf_w/rhoz on diagnosed interfaces.
    if (present(mf_a_out))          mf_a_out          = edmf_a2(1,:)
    if (present(mf_w_out))          mf_w_out          = edmf_w2(1,:)
    if (present(pblh_out))          pblh_out          = pblh_bl1(1)
    if(ma_enabled) then
      audit_profiles(:,1)=.5_RKIND*qke2(1,:)
      audit_profiles(:,2)=rqcloc;audit_profiles(:,3)=rqcblten2(1,:)
      audit_profiles(:,4)=rqiloc;audit_profiles(:,5)=rqiblten2(1,:)
      if(ma_skip_shadow) then
        audit_profiles(:,2)=ieee_value(0._RKIND,ieee_quiet_nan)
        audit_profiles(:,4)=ieee_value(0._RKIND,ieee_quiet_nan)
      endif
      audit_profiles(:,6)=tt;audit_profiles(:,7)=p;audit_profiles(:,8)=rho
      audit_profiles(:,9)=qv;audit_profiles(:,10)=qc;audit_profiles(:,11)=qi
      audit_profiles(:,12)=rthblten2(1,:);audit_profiles(:,13)=rqvblten2(1,:)
      audit_profiles(:,14)=el_pbl2(1,:);audit_profiles(:,15)=exch_h2(1,:);audit_profiles(:,16)=exch_m2(1,:)
      audit_profiles(:,17)=edmf_a2(1,:);audit_profiles(:,18)=edmf_w2(1,:)
      audit_profiles(:,19)=cldfra_bl2(1,:);audit_profiles(:,20)=edmf_ent2(1,:)
      audit_surface=[ts,qsfc1(1),hfx(1),lh(1),qfx(1),ust1(1),rmol1(1),pblh_bl1(1),znt(1), &
        mavail(1),dx1(1),pblh1(1),maxwidth1(1),maxmf1(1),ztop_plume1(1),wstar(1),ps_in,ts/exner(1),qv(1),dt]
      ma_dt=dt
      call ma_write(audit_surface,audit_profiles)
    endif

  contains

    ! dz on the (1,nz) shape bl_mynn_run expects, built from the 1-D
    ! grid array without a separate persistent copy
    function dz2() result(r)
      real(kind=RKIND) :: r(1,nz)
      r(1,:) = dz
    end function dz2

  end subroutine mynn_coupling_step

end module mynn_offline_coupling
