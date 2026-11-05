! mp_offline_driver.F90
!
! Standalone offline driver for the MPAS Thompson cloud microphysics
! scheme (module_mp_thompson.F, extracted unmodified from MPAS-Model),
! forced by the DEPHY-SCM GOAMAZON/SINGLE case (goamazon_bridge.F90).
! Mirrors the architecture of the companion rad_offline_driver.F90,
! adapted for a PROGNOSTIC (time-integrating) scheme instead of a
! diagnostic one -- see README.md for why that changes the driver
! shape (forcing must be added to the column each timestep, not just
! read and handed to the scheme once per hour).
!
! Single column, single process: ids=ide=its=ite=1, jds=jde=jts=jte=1,
! kds=1, kde=kme=kte=nz (the driver's own vertical grid -- see
! build_vertical_grid below), no padding levels needed (unlike RRTMG's
! documented kme=nlev+1 fix: mp_gt_driver has no half-level flux
! arrays that need an extra top interface).
!
! mp_gt_driver's field dummy arguments are explicit-shape
! (ims:ime,kms:kme,jms:jme) = (1,nz,1) 3-D arrays, even for this
! single column -- every thermodynamic/hydrometeor array below is
! declared that way (NOT as a plain (nz) 1-D array) so the actual
! argument shapes match exactly.
!
! FORCING APPLIED EACH MICROPHYSICS TIMESTEP dt_microp, BEFORE calling
! mp_gt_driver (see goamazon_bridge.F90 header for the GCSS/DEPHY
! "revealed forcing" convention this follows):
!   dtheta = [ tnta_adv - w*d(theta)/dz ] * dt_microp
!   drv    = [ tnrv_adv - w*d(rv)/dz    ] * dt_microp
! rt_adv (total water) is read but NOT separately applied by default
! (apply_rt_adv=.false.) -- see README.md "Status": only the vapor
! tendency is applied by default; condensate evolves from microphysics
! plus the vertical-advection term above acting on qc/qi. Vertical
! advection of qc/qr/qi/qs/qg by w is likewise not yet applied. Both
! are flagged as open validation items, not silent omissions.
!
! Pressure is NOT integrated prognostically/hydrostatically: at every
! step it is taken directly from goamazon_bridge's time+height
! interpolation of the file's own pa_forc field (forc_z=1 convention:
! the file's forcing is defined on fixed height levels with pressure
! evolving underneath them, exactly as the input data already
! represents it). pii (Exner function) is recomputed from that
! pressure every step.

program mp_offline_driver

  use mpas_kind_types, only : RKIND
  use mpas_atmphys_constants, only : gravity, cp, cpv, R_d, R_v, svpt0, ep_1, ep_2, epsilon, &
       &                             xls, xlv, xlv0, xlf, rho_a, rho_w, rho_s, cliq, cice, psat, &
       &                             re_qc_bg, re_qi_bg, re_qs_bg, re_qc_max, re_qi_max, re_qs_max
  use module_mp_thompson, only : mp_gt_driver, thompson_init
  use module_mp_wsm6, only : wsm6
  use mp_wsm6, only : mp_wsm6_init
  use goamazon_bridge
  use output_mp_netcdf
  use mynn_offline_coupling
  use bl_mynn_common, only : xlvcp, xlscp
  use plume_microphysics_thompson, only : step_plume_thompson, init_plume_thompson, plume_thompson_assume_initialized, n_plume_mp_calls, &
       &   n_mp_active, n_mp_idle, t_mp_active, t_mp_idle
  use bl_mynn_subroutines, only : edmf_type, uedmf_exc_mult, uedmf_kfac_max, uedmf_orig_gate, uedmf_thermo, uedmf_c0_auto, uedmf_qc_crit, &
       &                             uedmf_e0_floor, uedmf_econv, uedmf_econv_l, uedmf_hmode, uedmf_mp_mode, uedmf_mp_step, uedmf_t_core, uedmf_aw_flux, uedmf_sp_flux, &
       &                             uedmf_src, uedmf_src_th, uedmf_T_flux, uedmf_c_flux, uedmf_thl_flux, uedmf_diag_fo_pos, uedmf_diag_fo_neg, &
       &                             uedmf_c_eps, uedmf_c_delta, uedmf_c_lambda, uedmf_beta, uedmf_mu0, uedmf_chi, uedmf_c_gamma
  use mpas_atmphys_lightning, only : driver_lightning
  use process_diagnostics, only : pd_open, pd_begin, pd_finish, pd_close
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

  implicit none

  ! ---- physical constants: p0/Rd/rcp derived locally; cp/gravity/etc.
  ! now come from the REAL mpas_atmphys_constants (imported above) --
  ! P0 there is also 100000 Pa (same value, avoid the name clash by
  ! keeping our own lowercase p0 as a distinct local derived constant)
  real(kind=RKIND), parameter :: p0  = 100000.0_RKIND   ! reference pressure, Pa
  real(kind=RKIND), parameter :: rcp = R_d/cp

  ! ---- namelist-configurable run parameters (see namelist.mp_offline.sample) ----
  character(len=256) :: path_input  = 'data/GOAMAZON_SINGLE_SCM_driver.nc'
  character(len=256) :: path_output = 'output_mp_offline.nc'
  integer  :: nz         = 60       ! driver's own vertical levels
  real(kind=RKIND) :: ztop = 20000.0_RKIND   ! m, must not exceed the file's own lid
  real(kind=RKIND) :: dz0  = 50.0_RKIND      ! m, spacing of the lowest layer (stretched above)
  real(kind=RKIND) :: dt_microp = 30.0_RKIND ! s, microphysics timestep (sub-cycled between forcing times)
  logical :: apply_rt_adv = .false.          ! see header note above -- off by default in this first version
  logical :: legacy_temperature_tendency = .false. ! historical tnta-as-theta interface
  logical :: legacy_lightning = .false. ! retired EP and fixed IC/CG partition
  logical :: write_process_diagnostics = .false. ! observational hooks only
  logical :: use_mynn = .true.               ! MYNN surface layer + PBL, see mynn_offline_coupling.F90
  logical :: mynn_condensate_local_only = .false. ! apply only the EDMF-OFF (local mixing) part of MYNN's qc/qi tendencies
  integer :: mynn_edmf_type = 0                ! 0 = MPAS DMP_mf plumes (default); 1 = unified plume core inside MYNN (uedmf_core)
  real(kind=RKIND) :: mynn_uedmf_exc_mult = 1.0_RKIND, mynn_uedmf_kfac_max = 3.0_RKIND  ! unified core tunables (uncalibrated)
  real(kind=RKIND) :: mynn_uedmf_e0_floor = 1.0e-4_RKIND, mynn_uedmf_econv_l = 300.0_RKIND   ! min TKE seen by the closure (m2/s2); dissipation length of e_conv (m)
  ! closure parameters of the unified core: defaults = Table 2 of Cohen et al. (2020); hmode 0 = local height, 1 = plume top (paper)
  real(kind=RKIND) :: mynn_uedmf_c_eps=0.13_RKIND, mynn_uedmf_c_delta=0.52_RKIND, mynn_uedmf_c_lambda=0.3_RKIND, mynn_uedmf_beta=2.0_RKIND, &
       &              mynn_uedmf_mu0=4.0e-4_RKIND, mynn_uedmf_chi=0.25_RKIND, mynn_uedmf_c_gamma=0.075_RKIND
  integer :: mynn_uedmf_hmode = 0
  logical :: plume_channel_conserve = .false.   ! mode-2 channels: .false. = naive (clips at zero, CREATES mass: upper bound); .true. = scale the positive inputs so the column conserves mass (lower bound)
  integer :: mynn_uedmf_mp = 0                  ! 1 = explicit Thompson microphysics INSIDE the plume (feasibility option; see bl_mynn_subroutines.F90)
  integer(8) :: clk_start, clk_end, clk_rate
  real(kind=RKIND) :: bal(8) = 0.0_RKIND   ! mode-3 balance diagnostics (see apply_plume_balance)
  real(kind=RKIND) :: ext_net = 0.0_RKIND, ext_abs = 0.0_RKIND, ext_clip = 0.0_RKIND, ext_qg = 0.0_RKIND   ! conservation diagnostics of the plume species channels
  integer :: mynn_uedmf_econv = 0                ! 1 = convective TKE from plume detrainment (two-pass), see bl_mynn_subroutines.F90
  integer :: mynn_uedmf_thermo = 0               ! 0 = MYNN thermo; 1 = consistent theta=T/Exner; 2 = + latent heat of fusion
  real(kind=RKIND) :: mynn_uedmf_c0_auto = 0.0_RKIND, mynn_uedmf_qc_crit = 5.0e-4_RKIND  ! plume condensate sink (1/m; kg/kg); 0 = off
  integer :: mynn_uedmf_orig_gate = 0          ! 1 = unified core also requires MPAS's original activation criteria (control)
  integer :: mynn_edmf = 1                   ! 1 = MYNN EDMF mass-flux ON (MPAS default), 0 = OFF (control)
  logical :: couple_mynn_condensate = .true. ! apply MYNN's qc/qi tendencies (EDMF+mixing) before microphysics; .false. = control run
  character(len=32) :: microp_scheme = 'thompson'  ! 'thompson' or 'wsm6'

  namelist /mp_offline/ path_input, path_output, nz, ztop, dz0, dt_microp, apply_rt_adv, use_mynn, couple_mynn_condensate, mynn_condensate_local_only, mynn_edmf, &
       &                       mynn_edmf_type, mynn_uedmf_exc_mult, mynn_uedmf_kfac_max, mynn_uedmf_orig_gate, &
       &                       mynn_uedmf_thermo, mynn_uedmf_c0_auto, mynn_uedmf_qc_crit, &
       &                       mynn_uedmf_e0_floor, mynn_uedmf_econv, mynn_uedmf_econv_l, mynn_uedmf_hmode, mynn_uedmf_mp, plume_channel_conserve, &
       &                       mynn_uedmf_c_eps, mynn_uedmf_c_delta, mynn_uedmf_c_lambda, mynn_uedmf_beta, mynn_uedmf_mu0, mynn_uedmf_chi, mynn_uedmf_c_gamma, &
       &                microp_scheme, legacy_temperature_tendency, legacy_lightning, write_process_diagnostics

  ! ---- vertical grid (1-D, driver's own bookkeeping) ----
  real(kind=RKIND), allocatable :: z(:), dz(:)

  ! ---- prognostic column state, MPAS/WRF shape (ims:ime,kms:kme,jms:jme) = (1,nz,1) ----
  real(kind=RKIND), allocatable :: th(:,:,:), qv(:,:,:), qc(:,:,:), qr(:,:,:), &
       &                           qi(:,:,:), qs(:,:,:), qg(:,:,:), ni(:,:,:), nr(:,:,:)
  real(kind=RKIND), allocatable :: pres(:,:,:), pii(:,:,:), w3(:,:,:)
  real(kind=RKIND), allocatable :: re_cloud(:,:,:), re_ice(:,:,:), re_snow(:,:,:)

  ! ---- forcing at the current instant, on the driver grid (1-D) ----
  real(kind=RKIND), allocatable :: ta_adv(:), rv_adv(:), rt_adv(:)
  real(kind=RKIND) :: hfss_t, hfls_t, ps_t, ts_t, z0_t
  real(kind=RKIND) :: flash_rate, flash_rate_cg, flash_rate_ic
  real(kind=RKIND), allocatable :: u_col(:), v_col(:)
  real(kind=RKIND), allocatable :: tt(:), rho(:), rthblten(:), rqvblten(:)
  real(kind=RKIND), allocatable :: rqcblten(:), rqiblten(:), rqcblten_loc(:), rqiblten_loc(:)
  real(kind=RKIND), allocatable :: tke_diag(:), mixlen_diag(:), exch_h_diag(:), exch_m_diag(:)
  real(kind=RKIND), allocatable :: cldfra_bl_diag(:), qc_bl_diag(:), qi_bl_diag(:)
  real(kind=RKIND), allocatable :: mf_a_diag(:), mf_w_diag(:)
  real(kind=RKIND) :: pblh_diag

  ! ---- MPAS-style dimension bookkeeping (single column) ----
  integer :: ids,ide, jds,jde, kds,kde, ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte

  ! ---- precip accumulators and cloud-droplet-number arguments, shape (ims:ime,jms:jme) = (1,1) ----
  real(kind=RKIND) :: rainnc(1,1), rainncv(1,1), sr(1,1)
  real(kind=RKIND) :: snownc(1,1), snowncv(1,1), graupelnc(1,1), graupelncv(1,1)
  real(kind=RKIND) :: ntc(1,1), muc(1,1)
  integer :: has_reqc, has_reqi, has_reqs
  real(kind=RKIND), allocatable :: rainprod(:,:,:), evapprod(:,:,:)
  character(len=256) :: wsm6_errmsg
  integer :: wsm6_errflg

  type(goamazon_file_type) :: gfile
  real(kind=RKIND) :: t_sec, t_end
  real(kind=RKIND), allocatable :: ua0(:), va0(:)  ! unused outputs from goamazon_get_initial_profile (wind not needed by Thompson)
  real(kind=RKIND) :: ps0
  integer :: k, itimestep, iu
  logical :: file_exists

  ! ------------------------------------------------------------------
  ! namelist
  ! ------------------------------------------------------------------
  inquire(file='namelist.mp_offline', exist=file_exists)
  if (file_exists) then
    open(newunit=iu, file='namelist.mp_offline', status='old', action='read')
    read(iu, nml=mp_offline)
    close(iu)
  else
    write(*,'(A)') 'mp_offline_driver: namelist.mp_offline not found, using built-in defaults.'
  end if

  ! ------------------------------------------------------------------
  ! vertical grid: stretched, lowest layer dz0, geometric growth to
  ! reach ztop in nz layers. First-pass grid -- see README.md "Status":
  ! not yet checked against MONAN/MPAS's own real vertical spacing.
  ! ------------------------------------------------------------------
  allocate(z(nz), dz(nz))
  call build_vertical_grid(nz, dz0, ztop, z, dz)

  ! ------------------------------------------------------------------
  ! allocate column state at MPAS/WRF shape (1,nz,1)
  ! ------------------------------------------------------------------
  allocate(th(1,nz,1), qv(1,nz,1), qc(1,nz,1), qr(1,nz,1), qi(1,nz,1), &
       &   qs(1,nz,1), qg(1,nz,1), ni(1,nz,1), nr(1,nz,1))
  allocate(pres(1,nz,1), pii(1,nz,1), w3(1,nz,1))
  allocate(re_cloud(1,nz,1), re_ice(1,nz,1), re_snow(1,nz,1))
  allocate(ta_adv(nz), rv_adv(nz), rt_adv(nz))
  allocate(ua0(nz), va0(nz))
  allocate(u_col(nz), v_col(nz), tt(nz), rho(nz), rthblten(nz), rqvblten(nz))
  allocate(rqcblten(nz), rqiblten(nz), rqcblten_loc(nz), rqiblten_loc(nz))
  rqcblten_loc = 0.0_RKIND; rqiblten_loc = 0.0_RKIND
  allocate(tke_diag(nz), mixlen_diag(nz), exch_h_diag(nz), exch_m_diag(nz))
  allocate(cldfra_bl_diag(nz), qc_bl_diag(nz), qi_bl_diag(nz))
  allocate(mf_a_diag(nz), mf_w_diag(nz))
  pblh_diag = 0.0_RKIND
  tke_diag = 0.0_RKIND; mixlen_diag = 0.0_RKIND; exch_h_diag = 0.0_RKIND; exch_m_diag = 0.0_RKIND
  cldfra_bl_diag = 0.0_RKIND; qc_bl_diag = 0.0_RKIND; qi_bl_diag = 0.0_RKIND
  mf_a_diag = 0.0_RKIND; mf_w_diag = 0.0_RKIND
  rqcblten = 0.0_RKIND; rqiblten = 0.0_RKIND

  ! ------------------------------------------------------------------
  ! read the GOAMAZON forcing file and interpolate the initial profile
  ! onto our grid. Sections like pres(1,:,1) of a (1,nz,1)-shaped array
  ! are contiguous (the size-1 leading dimension costs nothing), so
  ! passing them to goamazon_bridge's plain (nz) dummy arguments is a
  ! straight aliased copy, no transpose, no temporary.
  ! ------------------------------------------------------------------
  call goamazon_open(trim(path_input), gfile, legacy_temperature_tendency)
  call goamazon_get_initial_profile(gfile, nz, z, pres(1,:,1), th(1,:,1), qv(1,:,1), &
       &                            qc(1,:,1), qi(1,:,1), ua0, va0, ps0)
  qr = 0.0_RKIND; qs = 0.0_RKIND; qg = 0.0_RKIND; ni = 0.0_RKIND; nr = 0.0_RKIND
  pii = (pres/p0)**rcp
  u_col = ua0; v_col = va0
  tt = th(1,:,1)*pii(1,:,1)
  rho = pres(1,:,1) / (R_d*tt*(1.0_RKIND + 0.61_RKIND*qv(1,:,1)))

  ! ------------------------------------------------------------------
  ! single-column MPAS-style dimensions
  ! ------------------------------------------------------------------
  ids=1; ide=2; jds=1; jde=2; kds=1; kde=nz   ! ide/jde=2: mp_gt_driver loops i=its..MIN(ite,ide-1); ide=1 made it a no-op
  ims=1; ime=1; jms=1; jme=1; kms=1; kme=nz
  its=1; ite=1; jts=1; jte=1; kts=1; kte=nz
  allocate(rainprod(1,nz,1), evapprod(1,nz,1)); rainprod = 0.0_RKIND; evapprod = 0.0_RKIND

  ! ------------------------------------------------------------------
  ! Thompson initialization -- once, before the time loop. l_mp_tables
  ! = .false. means the lookup tables are computed on the fly (no
  ! external MP_THOMPSON_*.DBL files needed -- see README.md).
  ! ntc/muc (per-column, passed as ARGUMENTS to mp_gt_driver, NOT the
  ! module's own internal, private Nt_c/mu_c): GOAMAZON is a
  ! continental (land) case, so following
  ! init_thompson_clouddroplets_forMPAS's own landmask.eq.1 branch,
  ! ntc=300.e6 and muc=min(15, 1000.e6/ntc + 2) -- DEPHY has no CCN
  ! field here to do better than this fixed continental default.
  ! ------------------------------------------------------------------
  if (trim(microp_scheme) == 'thompson') call thompson_init(.false.)
  ntc = 300.0e6_RKIND
  muc = min(15.0_RKIND, 1000.0e6_RKIND/ntc(1,1) + 2.0_RKIND)
  if (trim(microp_scheme) == 'wsm6') then
    ! hail_opt=0 matches MPAS's own hardcoded default (see
    ! mpas_atmphys_driver_microphysics.F's "hail_opt" parameter) --
    ! not exposed as a namelist option there either.
    call mp_wsm6_init(den0=rho_a, denr=rho_w, dens=rho_s, cl=cliq, cpv=cpv, hail_opt=0, &
         &            errmsg=wsm6_errmsg, errflg=wsm6_errflg)
    if (wsm6_errflg /= 0) then
      write(0,'(A)') 'mp_offline_driver: mp_wsm6_init failed: '//trim(wsm6_errmsg)
      stop 1
    end if
  end if
  bl_mynn_edmf_c = mynn_edmf
  edmf_type = mynn_edmf_type; uedmf_exc_mult = mynn_uedmf_exc_mult; uedmf_kfac_max = mynn_uedmf_kfac_max; uedmf_orig_gate = mynn_uedmf_orig_gate
  uedmf_thermo = mynn_uedmf_thermo; uedmf_c0_auto = mynn_uedmf_c0_auto; uedmf_qc_crit = mynn_uedmf_qc_crit
  uedmf_e0_floor = mynn_uedmf_e0_floor; uedmf_econv = mynn_uedmf_econv; uedmf_econv_l = mynn_uedmf_econv_l
  uedmf_hmode = mynn_uedmf_hmode
  uedmf_mp_mode = mynn_uedmf_mp
  if (mynn_uedmf_mp >= 1) then
    uedmf_mp_step => step_plume_thompson
    ! grid microphysics = Thompson already initialised the tables in this run; otherwise initialise them here
    if (trim(microp_scheme) == 'thompson') then
      call plume_thompson_assume_initialized()
    else
      call init_plume_thompson()
    end if
  end if
  call system_clock(clk_start, clk_rate)
  uedmf_c_eps = mynn_uedmf_c_eps; uedmf_c_delta = mynn_uedmf_c_delta; uedmf_c_lambda = mynn_uedmf_c_lambda; uedmf_beta = mynn_uedmf_beta
  uedmf_mu0 = mynn_uedmf_mu0; uedmf_chi = mynn_uedmf_chi; uedmf_c_gamma = mynn_uedmf_c_gamma
  if (use_mynn) call mynn_coupling_init(nz)

  has_reqc = 0; has_reqi = 0; has_reqs = 0   ! effective-radius diagnostics not requested
  rainnc = 0.0_RKIND; rainncv = 0.0_RKIND; sr = 0.0_RKIND
  snownc = 0.0_RKIND; snowncv = 0.0_RKIND; graupelnc = 0.0_RKIND; graupelncv = 0.0_RKIND

  call output_open(trim(path_output), nz, z)
  if (write_process_diagnostics) call pd_open(trim(path_output)//'.processes.csv', nz)

  ! ------------------------------------------------------------------
  ! time loop: sub-cycle at dt_microp between the file's own (coarser,
  ! irregular) forcing times, exactly like MPAS's own n_microp
  ! sub-cycling of mp_gt_driver within one host physics timestep.
  ! ------------------------------------------------------------------
  t_sec = gfile%time(1)
  t_end = gfile%time(gfile%ntime)
  itimestep = 0

  call output_write_timestep(t_sec, nz, th(1,:,1), qv(1,:,1), qc(1,:,1), qr(1,:,1), &
       &                      qi(1,:,1), qs(1,:,1), qg(1,:,1), pres(1,:,1), &
       &                      rainncv(1,1), rainnc(1,1))

  do while (t_sec < t_end)

    call goamazon_get_forcing_at_time(gfile, t_sec, nz, z, pres(1,:,1), ta_adv, rv_adv, rt_adv, &
         &                            w3(1,:,1), hfss_t, hfls_t, ps_t, ts_t, z0_t, u_col, v_col)
    pii = (pres/p0)**rcp
    tt  = th(1,:,1)*pii(1,:,1)
    rho = pres(1,:,1) / (R_d*tt*(1.0_RKIND + 0.61_RKIND*qv(1,:,1)))

    ! --- MYNN surface layer + PBL: thermodynamic mixing tendencies,
    !     added to th/qv alongside the large-scale forcing below (see
    !     mynn_offline_coupling.F90 header for what this call does and
    !     leaves out)
    if (use_mynn) then
      call mynn_coupling_step(nz, dt_microp, z, dz, u_col, v_col, th(1,:,1), tt, pres(1,:,1), &
           &                   pii(1,:,1), rho, qv(1,:,1), qc(1,:,1), qi(1,:,1), qs(1,:,1), &
           &                   ts_t, z0_t, ps_t, itimestep+1, rthblten, rqvblten, &
           &                   rqcblten_out=rqcblten, rqiblten_out=rqiblten, &
           &                   rqcblten_loc_out=rqcblten_loc, rqiblten_loc_out=rqiblten_loc, &
           &                   tke_out=tke_diag, mixing_length_out=mixlen_diag, &
           &                   exch_h_out=exch_h_diag, exch_m_out=exch_m_diag, &
           &                   cldfra_bl_out=cldfra_bl_diag, qc_bl_out=qc_bl_diag, qi_bl_out=qi_bl_diag, &
           &                   mf_a_out=mf_a_diag, mf_w_out=mf_w_diag, pblh_out=pblh_diag)
    else
      rthblten = 0.0_RKIND; rqvblten = 0.0_RKIND
      rqcblten = 0.0_RKIND; rqiblten = 0.0_RKIND
    end if

    ! --- apply forcing: horizontal advection (from the file) plus
    !     vertical advection -w*d(.)/dz using the column's own current
    !     gradient (see header note on the revealed-forcing convention)
    call apply_forcing_and_advance(nz, dz, dt_microp, w3(1,:,1), ta_adv, rv_adv, th(1,:,1), qv(1,:,1))
    th(1,:,1) = th(1,:,1) + dt_microp*rthblten
    qv(1,:,1) = max(qv(1,:,1) + dt_microp*rqvblten, 0.0_RKIND)
    ! --- EXPERIMENTAL, NOT MASS-CONSERVING (upper-bound sensitivity only; see the report): the flux form lacks the
    !     plume's source term, so clipping at zero creates mass (~14% of the activity in the GOAMAZON run).
    ! --- plume species channels (uedmf_mp_mode=2): rain/ice/snow/graupel carried by the plume are applied to the
    !     grid state as the divergence of (plume flux - environment transport), F = sum(rho a w q_u) - sum(rho a w)*q_env,
    !     explicit forward Euler, before the grid microphysics. Interface j is the bottom of layer j (F(1)=0 at the surface).
    if (edmf_type == 1 .and. uedmf_mp_mode == 2 .and. use_mynn) then
      call apply_plume_species(nz, dz, rho, dt_microp, qr(1,:,1), qi(1,:,1), qs(1,:,1), qg(1,:,1), &
           &                   ext_net, ext_abs, ext_clip, ext_qg, plume_channel_conserve)
    end if
    ! EDMF-to-microphysics coupling: cloud water/ice detrained by MYNN's
    ! mass-flux updraft plumes becomes a genuine qc/qi source term here,
    ! feeding directly into whichever scheme (Thompson/WSM6) is about to
    ! run -- see mynn_offline_coupling.F90's header for what this is and
    ! why it is exactly zero whenever the EDMF scheme itself is inactive
    ! (its own internal size/flux criteria unmet), which is the
    ! physically correct outcome, not a placeholder.
    if (couple_mynn_condensate) then
      if (mynn_condensate_local_only) then
        qc(1,:,1) = max(qc(1,:,1) + dt_microp*rqcblten_loc, 0.0_RKIND)
        qi(1,:,1) = max(qi(1,:,1) + dt_microp*rqiblten_loc, 0.0_RKIND)
      else
        qc(1,:,1) = max(qc(1,:,1) + dt_microp*rqcblten, 0.0_RKIND)
        qi(1,:,1) = max(qi(1,:,1) + dt_microp*rqiblten, 0.0_RKIND)
      end if
    end if

    ! --- mode 3: species channels WITH the plume's source terms and precipitation fallout (balanced version of mode 2)
    if (edmf_type == 1 .and. uedmf_mp_mode == 3 .and. use_mynn) then
      call apply_plume_balance(nz, dz, rho, pii(1,:,1), dt_microp, th(1,:,1), qv(1,:,1), qc(1,:,1), &
           &                   qr(1,:,1), qi(1,:,1), qs(1,:,1), qg(1,:,1), bal)
    end if
    if (apply_rt_adv) then
      ! total-water tendency, added to qv only (see README.md open item
      ! -- this is a simplification, not yet partitioned to qc/qi)
      do k = 1, nz
        qv(1,k,1) = qv(1,k,1) + rt_adv(k) * dt_microp
      end do
    end if

    itimestep = itimestep + 1
    if (write_process_diagnostics) call pd_begin(process_state())
    if (trim(microp_scheme) == 'wsm6') then
      call wsm6( &
           th=th, q=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, den=spread(spread(rho,1,1),3,1), &
           pii=pii, p=pres, &
           delz=spread(spread(dz,1,1),3,1), delt=dt_microp, g=gravity, cpd=cp, cpv=cpv, &
           rd=R_d, rv=R_v, t0c=svpt0, ep1=ep_1, ep2=ep_2, qmin=epsilon, &
           xls=xls, xlv0=xlv, xlf0=xlf, den0=rho_a, denr=rho_w, cliq=cliq, cice=cice, psat=psat, &
           rain=rainnc, rainncv=rainncv, snow=snownc, snowncv=snowncv, &
           graupel=graupelnc, graupelncv=graupelncv, sr=sr, &
           has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs, &
           re_qc_bg=re_qc_bg, re_qi_bg=re_qi_bg, re_qs_bg=re_qs_bg, &
           re_qc_max=re_qc_max, re_qi_max=re_qi_max, re_qs_max=re_qs_max, &
           re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
           ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde, &
           ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme, &
           its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte, &
           errmsg=wsm6_errmsg, errflg=wsm6_errflg )
      if (wsm6_errflg /= 0) then
        write(0,'(A)') 'mp_offline_driver: wsm6 error: '//trim(wsm6_errmsg)
        stop 1
      end if
    else
      call mp_gt_driver( &
           qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr,           &
           th=th, pii=pii, p=pres, w=w3, dz=spread(spread(dz,1,1),3,1), &
           dt_in=dt_microp, itimestep=itimestep, &
           RAINNC=rainnc, RAINNCV=rainncv, SNOWNC=snownc, SNOWNCV=snowncv, &
           GRAUPELNC=graupelnc, GRAUPELNCV=graupelncv, SR=sr, &
           rainprod=rainprod, evapprod=evapprod, &
           re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
           has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs, &
           ntc=ntc, muc=muc, &
           ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde, &
           ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme, &
           its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte )
    end if

    ! diagnostic lightning flash rate (mpas_atmphys_lightning.F90,
    ! P. Kubota/INPE-CPTEC prototype, MET-576-4) -- generic column
    ! diagnostic, unmodified, runs identically after Thompson or WSM6
    ! (it only reads qg/qs/qc/w/T/rho, not scheme-internal state).
    ! cu_active=.false. always: this SCM has no cumulus parameterization
    ! (microphysics is resolved explicitly, not sub-grid), so
    ! driver_lightning always takes its EP (electrification potential)
    ! branch, never the Price-Rind/cloud-top-height branch -- dx_km and
    ! cldtop_km below are therefore unused, passed only because the
    ! interface requires them.
    if (write_process_diagnostics) call pd_finish(t_sec + dt_microp, process_state())
    flash_rate = ieee_value(0.0_RKIND, ieee_quiet_nan)
    flash_rate_cg = flash_rate
    flash_rate_ic = flash_rate
    if (legacy_lightning) then
      tt = th(1,:,1)*pii(1,:,1)
      call driver_lightning(nz, z, tt-273.15_RKIND, rho, w3(1,:,1), qg(1,:,1), qs(1,:,1), qc(1,:,1), &
         &                 1.0_RKIND, .false., 0.0_RKIND, &
         &                 flash_rate, flash_rate_cg, flash_rate_ic)
    end if

    t_sec = t_sec + dt_microp

    call output_write_timestep(t_sec, nz, th(1,:,1), qv(1,:,1), qc(1,:,1), qr(1,:,1), &
         &                      qi(1,:,1), qs(1,:,1), qg(1,:,1), pres(1,:,1), &
         &                      rainncv(1,1), rainnc(1,1), flash_rate, flash_rate_cg, flash_rate_ic, &
         &                      snownc(1,1), snowncv(1,1), &
         &                      tke_diag, mixlen_diag, exch_h_diag, exch_m_diag, &
         &                      cldfra_bl_diag, qc_bl_diag, qi_bl_diag, mf_a_diag, mf_w_diag, pblh_diag)

  end do

  call output_close()
  if (write_process_diagnostics) call pd_close()
  call goamazon_close(gfile)

  call system_clock(clk_end)
  if (edmf_type == 1) write(*,'(A,F10.2,A,F10.2,A,I12)') 'TIMING total(s)=', real(clk_end-clk_start,8)/real(clk_rate,8), &
       '  uedmf_core(s)=', uedmf_t_core, '  chamadas microfisica-na-pluma=', n_plume_mp_calls
  if (edmf_type == 1 .and. uedmf_mp_mode == 1) write(*,'(A,I10,A,ES10.3,A,I10,A,ES10.3)') 'CUSTO ativas: n=', n_mp_active, ' t(s)=', t_mp_active, '   ociosas: n=', n_mp_idle, ' t(s)=', t_mp_idle

  if (edmf_type == 1 .and. uedmf_mp_mode == 2) write(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
       'CANAIS: |soma coluna|/soma|.| = ', ext_net/max(ext_abs,1.0e-30_RKIND), '   perda por corte em zero/soma|.| = ', &
       ext_clip/max(ext_abs,1.0e-30_RKIND), '   graupel entregue (kg/m2, acum.) = ', ext_qg, '   soma|.| (kg/m2) = ', ext_abs

  if (edmf_type == 1 .and. uedmf_mp_mode == 3) then
    write(*,'(A)') 'BALANCO (modo 3), residuos relativos a soma dos termos:'
    write(*,'(A,ES10.3)') '  agua, fechamento da coluna (emulando o solver)      = ', bal(1)/max(bal(2),1.0e-30_RKIND)
    write(*,'(A,ES10.3)') '  theta_il, fechamento da coluna (emulando o solver)  = ', bal(3)/max(bal(4),1.0e-30_RKIND)
    write(*,'(A,ES10.3)') '  calor latente do Thompson vs constantes do MYNN     = ', bal(5)/max(bal(6),1.0e-30_RKIND)
    write(*,'(A,ES10.3)') '  massa criada por corte em zero / soma dos termos    = ', bal(7)/max(bal(8),1.0e-30_RKIND)
    write(*,'(A,ES10.3)') '  massa criada pelo Thompson na camada / queda total  = ', uedmf_diag_fo_neg/max(uedmf_diag_fo_pos,1.0e-30_RKIND)
    write(*,'(A,ES10.3,A,ES10.3)') '  massa criada por corte em zero (kg/m2 acum.) = ', bal(7), '   soma |tendencias aplicadas| (kg/m2) = ', bal(8)
  end if

contains

  function process_state() result(state)
    real(kind=RKIND) :: state(6,nz)
    state(1,:) = qv(1,:,1)
    state(2,:) = qc(1,:,1)
    state(3,:) = qr(1,:,1)
    state(4,:) = qi(1,:,1)
    state(5,:) = qs(1,:,1)
    state(6,:) = qg(1,:,1)
  end function process_state

  subroutine apply_plume_balance(nz, dz, rho, pii, dt, th, qv, qc, qr, qi, qs, qg, bal)
    ! Grid-side half of the plume species channels (uedmf_mp_mode=3). The solver has already delivered, through
    ! qt / thl / qc, the plume's total water, ice-liquid theta and liquid cloud water. Here:
    !   external species s in {qr,qi,qs,qg}: tend = -d(F_env)/dz + src(s)   (flux form + the plume's production, which
    !                                        already includes the precipitation that fell out of the plume, deposited here)
    !   vapour: += src(qv) - B_ext, where B_ext = raw flux divergence of the external species (what the solver put in the
    !                                        vapour channel because it carries them inside qt) is taken back out
    !   cloud water: += src(qc)
    !   theta: += src_th + sum_s L_s*B_s/(cp*Pi)  (latent heating in the plume, plus the heat the ice-liquid theta hid
    !                                        in the external species)
    ! The bal(:) accumulators emulate the solver with the same exported fluxes (explicit flux form) to check that the
    ! column closes in total water and in ice-liquid theta.
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: dz(nz), rho(nz), pii(nz), dt
    real(kind=RKIND), intent(inout) :: th(nz), qv(nz), qc(nz), qr(nz), qi(nz), qs(nz), qg(nz), bal(8)
    real(kind=RKIND) :: q(nz,4), qenv(nz+1), Fenv(nz+1), Fraw(nz+1,4), A(nz,4), B(nz,4), tend(nz,4)
    real(kind=RKIND) :: Lsp(4), qbar, qnew, w, tv, tc, tth, thl_bar(nz), Fv(nz+1), Fc(nz+1), Fthl(nz+1)
    real(kind=RKIND) :: Wsum, Wabs, Esum, Eabs, LHabs, LHref, dthil, dth_tot, dtheta_mynn, latmis, cliprm, appabs
    real(kind=RKIND) :: Bext, dqv_c, dqc_c
    integer :: k, s, j

    q(:,1) = qr; q(:,2) = qi; q(:,3) = qs; q(:,4) = qg
    Lsp = (/ xlvcp, xlscp, xlscp, xlscp /)        ! qr is liquid; qi, qs, qg are ice (the same split the plume's thl uses)

    do s = 1, 4
      do j = 1, nz+1
        Fraw(j,s) = uedmf_sp_flux(j,s)
      end do
      do j = 1, nz+1
        if (j == 1) then
          qbar = q(1,s)
        else if (j == nz+1) then
          qbar = q(nz,s)
        else
          qbar = 0.5_RKIND*(q(j-1,s) + q(j,s))
        end if
        Fenv(j) = Fraw(j,s) - uedmf_aw_flux(j)*qbar
      end do
      do k = 1, nz
        A(k,s) = -(Fenv(k+1) - Fenv(k))/(rho(k)*dz(k))
        B(k,s) = -(Fraw(k+1,s) - Fraw(k,s))/(rho(k)*dz(k))
        tend(k,s) = A(k,s) + uedmf_src(k,2+s)
      end do
    end do

    ! ---- budget emulation of the solver part (explicit flux form with the exported hand-off fluxes) ----
    do j = 1, nz+1
      if (j == 1) then
        Fv(j) = 0.0_RKIND; Fc(j) = 0.0_RKIND; Fthl(j) = 0.0_RKIND
      else
        if (j == nz+1) then
          qbar = qv(nz)
        else
          qbar = 0.5_RKIND*(qv(j-1) + qv(j))
        end if
        Fv(j) = (uedmf_T_flux(j) - uedmf_c_flux(j)) - uedmf_aw_flux(j)*qbar
        if (j == nz+1) then
          qbar = qc(nz)
        else
          qbar = 0.5_RKIND*(qc(j-1) + qc(j))
        end if
        Fc(j) = uedmf_c_flux(j) - uedmf_aw_flux(j)*qbar
      end if
    end do
    do k = 1, nz
      thl_bar(k) = th(k) - xlvcp*qc(k)/pii(k) - xlscp*qi(k)/pii(k)
    end do
    do j = 2, nz+1
      if (j == nz+1) then
        qbar = thl_bar(nz)
      else
        qbar = 0.5_RKIND*(thl_bar(j-1) + thl_bar(j))
      end if
      Fthl(j) = uedmf_thl_flux(j) - uedmf_aw_flux(j)*qbar
    end do

    Wsum = 0.0_RKIND; Wabs = 0.0_RKIND; Esum = 0.0_RKIND; Eabs = 0.0_RKIND; LHabs = 0.0_RKIND; LHref = 0.0_RKIND
    do k = 1, nz
      w = rho(k)*dz(k)
      tv  = -(Fv(k+1) - Fv(k))/w
      tc  = -(Fc(k+1) - Fc(k))/w
      tth = -(Fthl(k+1) - Fthl(k))/w
      Bext = B(k,1) + B(k,2) + B(k,3) + B(k,4)
      dqv_c = uedmf_src(k,1) - Bext
      dqc_c = uedmf_src(k,2)
      ! water: solver part + every grid-side term
      Wsum = Wsum + w*(tv + tc + tend(k,1)+tend(k,2)+tend(k,3)+tend(k,4) + dqv_c + dqc_c)
      Wabs = Wabs + w*(abs(tv)+abs(tc)+abs(tend(k,1))+abs(tend(k,2))+abs(tend(k,3))+abs(tend(k,4))+abs(dqv_c)+abs(dqc_c))
      ! ice-liquid theta: theta as the solver would reconstruct it, plus the grid-side heat, minus every species' L*q/(cp*Pi)
      dtheta_mynn = tth + xlvcp*tc/pii(k)
      dth_tot = dtheta_mynn + uedmf_src_th(k) + (xlvcp*B(k,1) + xlscp*(B(k,2)+B(k,3)+B(k,4)))/pii(k)
      dthil = dth_tot - (xlvcp/pii(k))*(tc + uedmf_src(k,2) + tend(k,1)) &
           &          - (xlscp/pii(k))*(tend(k,2)+tend(k,3)+tend(k,4))
      Esum = Esum + w*dthil
      Eabs = Eabs + w*(abs(tth) + abs(uedmf_src_th(k)) + abs(xlvcp*uedmf_src(k,2)/pii(k)) &
           &          + abs(xlvcp*tend(k,1)/pii(k)) + abs(xlscp*(tend(k,2)+tend(k,3)+tend(k,4))/pii(k)))
      ! latent heating of the microphysics call vs L * (condensate produced), with the MYNN's own constants
      latmis = uedmf_src_th(k) - (xlvcp*(uedmf_src(k,2) + uedmf_src(k,3)) + xlscp*(uedmf_src(k,4)+uedmf_src(k,5)+uedmf_src(k,6)))/pii(k)
      LHabs = LHabs + w*abs(latmis)
      LHref = LHref + w*abs(uedmf_src_th(k))
    end do
    bal(1) = bal(1) + abs(Wsum)*dt; bal(2) = bal(2) + Wabs*dt
    bal(3) = bal(3) + abs(Esum)*dt; bal(4) = bal(4) + Eabs*dt
    bal(5) = bal(5) + LHabs*dt;     bal(6) = bal(6) + LHref*dt

    ! ---- apply ----
    cliprm = 0.0_RKIND; appabs = 0.0_RKIND
    do k = 1, nz
      w = rho(k)*dz(k)
      Bext = B(k,1) + B(k,2) + B(k,3) + B(k,4)
      do s = 1, 4
        qnew = q(k,s) + dt*tend(k,s)
        appabs = appabs + w*dt*abs(tend(k,s))
        if (qnew < 0.0_RKIND) cliprm = cliprm + w*(-qnew)
        q(k,s) = max(qnew, 0.0_RKIND)
      end do
      qnew = qv(k) + dt*(uedmf_src(k,1) - Bext)
      appabs = appabs + w*dt*abs(uedmf_src(k,1) - Bext)
      if (qnew < 0.0_RKIND) cliprm = cliprm + w*(-qnew)
      qv(k) = max(qnew, 0.0_RKIND)
      qnew = qc(k) + dt*uedmf_src(k,2)
      appabs = appabs + w*dt*abs(uedmf_src(k,2))
      if (qnew < 0.0_RKIND) cliprm = cliprm + w*(-qnew)
      qc(k) = max(qnew, 0.0_RKIND)
      th(k) = th(k) + dt*(uedmf_src_th(k) + (xlvcp*B(k,1) + xlscp*(B(k,2)+B(k,3)+B(k,4)))/pii(k))
    end do
    bal(7) = bal(7) + cliprm; bal(8) = bal(8) + appabs
    qr = q(:,1); qi = q(:,2); qs = q(:,3); qg = q(:,4)
  end subroutine apply_plume_balance

  subroutine apply_plume_species(nz, dz, rho, dt, qr, qi, qs, qg, net, absum, clip, qg_del, conserve)
    integer, intent(in) :: nz
    logical, intent(in) :: conserve
    real(kind=RKIND) :: ppos, cneg, sfac, qnew2
    real(kind=RKIND), intent(in) :: dz(nz), rho(nz), dt
    real(kind=RKIND), intent(inout) :: qr(nz), qi(nz), qs(nz), qg(nz), net, absum, clip, qg_del
    real(kind=RKIND) :: q(nz,4), Fj(nz+1,4), qenv, tend(nz,4), colnet, colabs, qnew
    integer :: k, is_
    q(:,1) = qr; q(:,2) = qi; q(:,3) = qs; q(:,4) = qg
    do is_ = 1, 4
      Fj(1,is_) = 0.0_RKIND
      do k = 2, nz
        qenv = (q(k-1,is_)*dz(k) + q(k,is_)*dz(k-1))/(dz(k-1)+dz(k))
        Fj(k,is_) = uedmf_sp_flux(k,is_) - uedmf_aw_flux(k)*qenv
      end do
      qenv = q(nz,is_)
      Fj(nz+1,is_) = uedmf_sp_flux(nz+1,is_) - uedmf_aw_flux(nz+1)*qenv
      do k = 1, nz
        tend(k,is_) = -(Fj(k+1,is_) - Fj(k,is_))/(rho(k)*dz(k))
      end do
    end do
    do is_ = 1, 4
      colnet = 0.0_RKIND; colabs = 0.0_RKIND
      ! mass created by clipping at zero (negative results set to 0) and the positive inputs, per column
      ppos = 0.0_RKIND; cneg = 0.0_RKIND
      do k = 1, nz
        qnew = q(k,is_) + dt*tend(k,is_)
        if (qnew < 0.0_RKIND) cneg = cneg + rho(k)*dz(k)*(-qnew)
        if (tend(k,is_) > 0.0_RKIND) ppos = ppos + rho(k)*dz(k)*dt*tend(k,is_)
      end do
      sfac = 1.0_RKIND
      if (conserve .and. ppos > 0.0_RKIND) sfac = max(0.0_RKIND, 1.0_RKIND - cneg/ppos)
      do k = 1, nz
        colnet = colnet + rho(k)*dz(k)*tend(k,is_)
        colabs = colabs + rho(k)*dz(k)*abs(tend(k,is_))
        if (conserve .and. tend(k,is_) > 0.0_RKIND) then
          qnew = q(k,is_) + sfac*dt*tend(k,is_)
        else
          qnew = q(k,is_) + dt*tend(k,is_)
        end if
        if (qnew < 0.0_RKIND) clip = clip + rho(k)*dz(k)*(-qnew)
        q(k,is_) = max(qnew, 0.0_RKIND)
        if (is_ == 4) then
          if (conserve .and. tend(k,is_) > 0.0_RKIND) then
            qg_del = qg_del + sfac*dt*rho(k)*dz(k)*tend(k,is_)
          else
            qg_del = qg_del + dt*rho(k)*dz(k)*max(tend(k,is_),0.0_RKIND)
          end if
        end if
      end do
      net = net + abs(colnet)*dt; absum = absum + colabs*dt
    end do
    qr = q(:,1); qi = q(:,2); qs = q(:,3); qg = q(:,4)
  end subroutine apply_plume_species

  !===================================================================
  subroutine build_vertical_grid(nz, dz0, ztop, z, dz)
  !===================================================================
    ! Geometric stretching: dz(k) = dz0 * ratio**(k-1), solved so that
    ! sum(dz(1:nz)) = ztop. Layer k spans [zw(k), zw(k+1)]; z(k) is the
    ! layer-CENTER height (what th/qv/... live at).
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: dz0, ztop
    real(kind=RKIND), intent(out) :: z(nz), dz(nz)
    real(kind=RKIND) :: ratio, lo, hi, mid, s
    real(kind=RKIND) :: zw(nz+1)
    integer :: k, iter

    ! bisect for the growth ratio that makes the geometric series sum to ztop
    lo = 1.0_RKIND; hi = 1.20_RKIND
    do iter = 1, 60
      mid = 0.5_RKIND*(lo+hi)
      if (abs(mid-1.0_RKIND) < 1.0e-12_RKIND) then
        s = dz0 * real(nz,kind=RKIND)
      else
        s = dz0 * (mid**nz - 1.0_RKIND) / (mid - 1.0_RKIND)
      end if
      if (s < ztop) then
        lo = mid
      else
        hi = mid
      end if
    end do
    ratio = 0.5_RKIND*(lo+hi)

    zw(1) = 0.0_RKIND
    do k = 1, nz
      zw(k+1) = zw(k) + dz0*ratio**(k-1)
    end do
    do k = 1, nz
      dz(k) = zw(k+1) - zw(k)
      z(k)  = 0.5_RKIND*(zw(k)+zw(k+1))
    end do
  end subroutine build_vertical_grid

  !===================================================================
  subroutine apply_forcing_and_advance(nz, dz, dt, w, ta_adv, rv_adv, th, qv)
  !===================================================================
    ! theta/qv update in place: horizontal advection (ta_adv/rv_adv,
    ! already on the driver grid) plus vertical advection -w*d(.)/dz
    ! using a centered difference on the column's OWN current profile
    ! (one-sided at the two end levels).
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: dz(nz), dt, w(nz), ta_adv(nz), rv_adv(nz)
    real(kind=RKIND), intent(inout) :: th(nz), qv(nz)
    real(kind=RKIND) :: dthdz(nz), dqvdz(nz)
    integer :: k

    do k = 2, nz-1
      dthdz(k) = (th(k+1)-th(k-1)) / (dz(k)+dz(k-1)) ! ~ centered, non-uniform grid
      dqvdz(k) = (qv(k+1)-qv(k-1)) / (dz(k)+dz(k-1))
    end do
    dthdz(1)  = (th(2)-th(1))   / dz(1)
    dqvdz(1)  = (qv(2)-qv(1))   / dz(1)
    dthdz(nz) = (th(nz)-th(nz-1)) / dz(nz-1)
    dqvdz(nz) = (qv(nz)-qv(nz-1)) / dz(nz-1)

    do k = 1, nz
      th(k) = th(k) + dt * (ta_adv(k) - w(k)*dthdz(k))
      qv(k) = qv(k) + dt * (rv_adv(k) - w(k)*dqvdz(k))
      qv(k) = max(qv(k), 0.0_RKIND)
    end do
  end subroutine apply_forcing_and_advance

end program mp_offline_driver
