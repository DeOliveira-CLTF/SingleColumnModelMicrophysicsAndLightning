! ensemble_driver.F90
!
! Runs the (unmodified) Thompson microphysics scheme against a full
! factorial ensemble of ERA5-based background columns x saturation
! perturbations, to test its cloud-formation/precipitation response
! across a range of REAL atmospheric backgrounds (rather than one
! hand-picked idealized sounding) without depending on whatever the
! real large-scale forcing at that time/place would or would not have
! actually produced -- see saturation_perturbation.F90 and
! updraft_forcing.F90 for why/how.
!
! No MYNN, no large-scale advective forcing time series here (each
! member is built from a single ERA5 instant, not a time-varying
! forcing file like GOAMAZON): the background profile is held FIXED
! except for vertical advection through its own gradient, driven by
! combined_w() (real ERA5 w, from its own omega, plus a synthetic
! pulse). This is deliberately simpler than mp_offline_driver.F90 --
! the point here is isolating Thompson's response to imposed ascent
! through a real thermodynamic background, not reproducing a
! particular historical event end-to-end.
!
! Member list: read from ensemble.nml (see ensemble.nml.sample) --
! NPOINTS_MAX x NSATCONFIGS_MAX entries, cross-producted (every point
! run against every saturation config).
!
! Output: one record per member appended to a single NetCDF summary
! (path_summary, via output_ensemble_netcdf.F90 -- a "member" dimension,
! point/satconfig labels as character variables, one variable per
! summary quantity) -- peak qc/qr/qi/qs/qg reached anywhere in the
! column, time-to-first-precip, accumulated surface precipitation,
! lightning and snow diagnostics -- plus, optionally, a full profile
! NetCDF per member (same output_mp_netcdf.F90 module the GOAMAZON
! driver uses) if write_profiles=.true., named
! "<path_profile_prefix>_NNN.nc".

program ensemble_driver

  use mpas_kind_types, only : RKIND
  use mpas_atmphys_constants, only : gravity, cp, cpv, R_d, R_v, svpt0, ep_1, ep_2, epsilon, &
       &                             xls, xlv, xlv0, xlf, rho_a, rho_w, rho_s, cliq, cice, psat, &
       &                             re_qc_bg, re_qi_bg, re_qs_bg, re_qc_max, re_qi_max, re_qs_max
  use module_mp_thompson, only : mp_gt_driver, thompson_init
  use module_mp_wsm6, only : wsm6
  use mp_wsm6, only : mp_wsm6_init
  use era5_bridge
  use saturation_perturbation, only : saturate_layer
  use updraft_forcing, only : combined_w
  use output_mp_netcdf
  use output_ensemble_netcdf
  use mpas_atmphys_lightning, only : driver_lightning

  implicit none

  integer, parameter :: NPOINTS_MAX = 10
  integer, parameter :: NSATCONFIGS_MAX = 10

  ! ---- ensemble member definitions (namelist) ----
  integer :: npoints = 0, nsatconfigs = 0
  character(len=256) :: path_pl(NPOINTS_MAX) = ''
  character(len=256) :: path_sl(NPOINTS_MAX) = ''
  character(len=64)  :: point_label(NPOINTS_MAX) = ''
  real(kind=RKIND)   :: target_lat(NPOINTS_MAX) = 0.0_RKIND, target_lon(NPOINTS_MAX) = 0.0_RKIND
  integer :: itime_pick(NPOINTS_MAX) = 16   ! default: index 16 = 15:00 UTC (hour 15, 1-based)

  real(kind=RKIND) :: sat_z_lo(NSATCONFIGS_MAX) = 0.0_RKIND   ! m
  real(kind=RKIND) :: sat_z_hi(NSATCONFIGS_MAX) = 0.0_RKIND   ! m
  real(kind=RKIND) :: sat_rh_target(NSATCONFIGS_MAX) = 1.0_RKIND
  character(len=32) :: satconfig_label(NSATCONFIGS_MAX) = ''

  ! ---- run parameters ----
  integer  :: nz = 60
  real(kind=RKIND) :: ztop = 20000.0_RKIND, dz0 = 50.0_RKIND
  real(kind=RKIND) :: dt_microp = 20.0_RKIND
  real(kind=RKIND) :: run_duration_s = 7200.0_RKIND   ! 2h per member
  real(kind=RKIND) :: pulse_sigma_z = 1000.0_RKIND    ! m
  real(kind=RKIND) :: pulse_duration_s = 1800.0_RKIND ! 30 min
  real(kind=RKIND) :: pulse_w_max = 3.0_RKIND         ! m/s, synthetic component
  character(len=256) :: path_summary = 'ensemble_summary.nc'
  character(len=256) :: path_profile_prefix = 'member'
  logical :: write_profiles = .false.
  character(len=32) :: microp_scheme = 'thompson'  ! 'thompson' or 'wsm6'

  namelist /ensemble/ npoints, nsatconfigs, path_pl, path_sl, point_label, &
       &              target_lat, target_lon, itime_pick, &
       &              sat_z_lo, sat_z_hi, sat_rh_target, satconfig_label, &
       &              nz, ztop, dz0, dt_microp, run_duration_s, &
       &              pulse_sigma_z, pulse_duration_s, pulse_w_max, &
       &              path_summary, path_profile_prefix, write_profiles, microp_scheme

  real(kind=RKIND), allocatable :: z(:), dz(:)
  real(kind=RKIND), allocatable :: th(:,:,:), qv(:,:,:), qc(:,:,:), qr(:,:,:), &
       &                           qi(:,:,:), qs(:,:,:), qg(:,:,:), ni(:,:,:), nr(:,:,:)
  real(kind=RKIND), allocatable :: pres(:,:,:), pii(:,:,:), w3(:,:,:)
  real(kind=RKIND), allocatable :: re_cloud(:,:,:), re_ice(:,:,:), re_snow(:,:,:)
  real(kind=RKIND), allocatable :: w_era5(:)

  real(kind=RKIND), parameter :: p0 = 100000.0_RKIND
  real(kind=RKIND), parameter :: rcp = R_d/cp

  integer :: ids,ide, jds,jde, kds,kde, ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte
  real(kind=RKIND) :: rainnc(1,1), rainncv(1,1), sr(1,1), ntc(1,1), muc(1,1)
  real(kind=RKIND) :: snownc(1,1), snowncv(1,1), graupelnc(1,1), graupelncv(1,1)
  integer :: has_reqc, has_reqi, has_reqs
  real(kind=RKIND), allocatable :: rainprod(:,:,:), evapprod(:,:,:)
  character(len=256) :: wsm6_errmsg
  integer :: wsm6_errflg

  type(era5_file_type) :: efile
  type(era5_column_type) :: col

  integer :: ip, is, iu, k, k_lo, k_hi, itimestep, member_count
  real(kind=RKIND) :: t_sec
  real(kind=RKIND) :: qc_peak, qr_peak, qi_peak, qs_peak, qg_peak, time_to_precip
  real(kind=RKIND) :: flash_rate, flash_rate_cg, flash_rate_ic
  real(kind=RKIND) :: flash_rate_peak, flash_rate_cg_peak, flash_rate_ic_peak
  real(kind=RKIND) :: flash_count_total, flash_count_cg_total, flash_count_ic_total
  logical :: file_exists
  integer :: ilon, ilat
  character(len=512) :: profile_path

  inquire(file='ensemble.nml', exist=file_exists)
  if (.not. file_exists) stop 'ensemble_driver: ensemble.nml not found -- see ensemble.nml.sample'
  open(newunit=iu, file='ensemble.nml', status='old', action='read')
  read(iu, nml=ensemble)
  close(iu)

  if (npoints < 1 .or. npoints > NPOINTS_MAX) stop 'ensemble_driver: npoints out of range'
  if (nsatconfigs < 1 .or. nsatconfigs > NSATCONFIGS_MAX) stop 'ensemble_driver: nsatconfigs out of range'

  allocate(z(nz), dz(nz))
  call build_vertical_grid(nz, dz0, ztop, z, dz)

  allocate(th(1,nz,1), qv(1,nz,1), qc(1,nz,1), qr(1,nz,1), qi(1,nz,1), &
       &   qs(1,nz,1), qg(1,nz,1), ni(1,nz,1), nr(1,nz,1))
  allocate(pres(1,nz,1), pii(1,nz,1), w3(1,nz,1))
  allocate(re_cloud(1,nz,1), re_ice(1,nz,1), re_snow(1,nz,1))
  allocate(w_era5(nz))

  ids=1; ide=2; jds=1; jde=2; kds=1; kde=nz   ! ide/jde=2: mp_gt_driver loops i=its..MIN(ite,ide-1); ide=1 made it a no-op
  ims=1; ime=1; jms=1; jme=1; kms=1; kme=nz
  its=1; ite=1; jts=1; jte=1; kts=1; kte=nz
  allocate(rainprod(1,nz,1), evapprod(1,nz,1)); rainprod = 0.0_RKIND; evapprod = 0.0_RKIND

  call thompson_init(.false.)
  if (trim(microp_scheme) == 'wsm6') then
    call mp_wsm6_init(den0=rho_a, denr=rho_w, dens=rho_s, cl=cliq, cpv=cpv, hail_opt=0, &
         &            errmsg=wsm6_errmsg, errflg=wsm6_errflg)
    if (wsm6_errflg /= 0) then
      write(0,'(A)') 'ensemble_driver: mp_wsm6_init failed: '//trim(wsm6_errmsg)
      stop 1
    end if
  end if

  call ensemble_summary_open(trim(path_summary))

  member_count = 0
  do ip = 1, npoints

    call era5_open(trim(path_pl(ip)), trim(path_sl(ip)), efile)
    call nearest_index(efile%lon, size(efile%lon), target_lon(ip), ilon)
    call nearest_index(efile%lat, size(efile%lat), target_lat(ip), ilat)
    call era5_read_column(efile, ilon, ilat, itime_pick(ip), col)

    do is = 1, nsatconfigs

      member_count = member_count + 1

      ! regrid the ERA5 column (TOA-to-surface order, see era5_bridge.F90)
      ! onto our own driver grid, surface-up
      call regrid_era5_to_driver(col, nz, z, pres(1,:,1), th(1,:,1), qv(1,:,1), &
           &                      qc(1,:,1), qi(1,:,1), w_era5)
      qr = 0.0_RKIND; qs = 0.0_RKIND; qg = 0.0_RKIND; ni = 0.0_RKIND; nr = 0.0_RKIND
      pii = (pres/p0)**rcp

      ! locate the driver-grid layers spanning [sat_z_lo,sat_z_hi] and saturate them
      call bracket_layers(nz, z, sat_z_lo(is), sat_z_hi(is), k_lo, k_hi)
      call saturate_layer(nz, pres(1,:,1), th(1,:,1)*pii(1,:,1), k_lo, k_hi, sat_rh_target(is), &
           &               qv(1,:,1), qc(1,:,1), qi(1,:,1))

      ntc = 300.0e6_RKIND
      muc = min(15.0_RKIND, 1000.0e6_RKIND/ntc(1,1) + 2.0_RKIND)
      has_reqc = 0; has_reqi = 0; has_reqs = 0
      rainnc = 0.0_RKIND; rainncv = 0.0_RKIND; sr = 0.0_RKIND
      snownc = 0.0_RKIND; snowncv = 0.0_RKIND; graupelnc = 0.0_RKIND; graupelncv = 0.0_RKIND
      qc_peak = 0.0_RKIND; qr_peak = 0.0_RKIND; qi_peak = 0.0_RKIND
      qs_peak = 0.0_RKIND; qg_peak = 0.0_RKIND; time_to_precip = -1.0_RKIND
      flash_rate_peak = 0.0_RKIND
      flash_rate_cg_peak = 0.0_RKIND; flash_rate_ic_peak = 0.0_RKIND
      flash_count_total = 0.0_RKIND; flash_count_cg_total = 0.0_RKIND; flash_count_ic_total = 0.0_RKIND

      if (write_profiles) then
        write(profile_path,'(A,A,I0,A)') trim(path_profile_prefix), '_', member_count, '.nc'
        call output_open(trim(profile_path), nz, z)
        call output_write_timestep(0.0_RKIND, nz, th(1,:,1), qv(1,:,1), qc(1,:,1), qr(1,:,1), &
             &                      qi(1,:,1), qs(1,:,1), qg(1,:,1), pres(1,:,1), 0.0_RKIND, 0.0_RKIND)
      end if

      t_sec = 0.0_RKIND
      itimestep = 0
      do while (t_sec < run_duration_s)
        w3(1,:,1) = combined_w(nz, z, w_era5, t_sec, 0.5_RKIND*(sat_z_lo(is)+sat_z_hi(is)), &
             &                 pulse_sigma_z, pulse_duration_s, pulse_w_max)
        pii = (pres/p0)**rcp
        call apply_vertical_advection(nz, dz, dt_microp, w3(1,:,1), th(1,:,1), qv(1,:,1))

        itimestep = itimestep + 1
        if (trim(microp_scheme) == 'wsm6') then
          block
            real(kind=RKIND) :: tt_col(nz), rho_col(nz)
            tt_col  = th(1,:,1)*pii(1,:,1)
            rho_col = pres(1,:,1) / (R_d*tt_col*(1.0_RKIND+0.61_RKIND*qv(1,:,1)))
            call wsm6( &
                 th=th, q=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, &
                 den=spread(spread(rho_col,1,1),3,1), pii=pii, p=pres, &
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
              write(0,'(A)') 'ensemble_driver: wsm6 error: '//trim(wsm6_errmsg)
              stop 1
            end if
          end block
        else
          call mp_gt_driver( &
               qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, ni=ni, nr=nr, &
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

        t_sec = t_sec + dt_microp

        block
          real(kind=RKIND) :: tt_lit(nz), rho_lit(nz)
          tt_lit  = th(1,:,1)*pii(1,:,1)
          rho_lit = pres(1,:,1) / (R_d*tt_lit*(1.0_RKIND+0.61_RKIND*qv(1,:,1)))
          call driver_lightning(nz, z, tt_lit-273.15_RKIND, rho_lit, w3(1,:,1), &
               &                 qg(1,:,1), qs(1,:,1), qc(1,:,1), 1.0_RKIND, .false., 0.0_RKIND, &
               &                 flash_rate, flash_rate_cg, flash_rate_ic)
        end block
        flash_rate_peak = max(flash_rate_peak, flash_rate)
        flash_rate_cg_peak = max(flash_rate_cg_peak, flash_rate_cg)
        flash_rate_ic_peak = max(flash_rate_ic_peak, flash_rate_ic)
        ! rectangle-rule accumulation: flash_rate is in flashes/min, dt_microp
        ! in seconds -- (dt_microp/60) converts the step's contribution to
        ! flashes/min * min = flashes
        flash_count_total    = flash_count_total    + flash_rate    * (dt_microp/60.0_RKIND)
        flash_count_cg_total = flash_count_cg_total + flash_rate_cg * (dt_microp/60.0_RKIND)
        flash_count_ic_total = flash_count_ic_total + flash_rate_ic * (dt_microp/60.0_RKIND)

        qc_peak = max(qc_peak, maxval(qc(1,:,1))); qr_peak = max(qr_peak, maxval(qr(1,:,1)))
        qi_peak = max(qi_peak, maxval(qi(1,:,1))); qs_peak = max(qs_peak, maxval(qs(1,:,1)))
        qg_peak = max(qg_peak, maxval(qg(1,:,1)))
        if (time_to_precip < 0.0_RKIND .and. rainncv(1,1) > 0.0_RKIND) time_to_precip = t_sec

        if (write_profiles) then
          call output_write_timestep(t_sec, nz, th(1,:,1), qv(1,:,1), qc(1,:,1), qr(1,:,1), &
               &                      qi(1,:,1), qs(1,:,1), qg(1,:,1), pres(1,:,1), &
               &                      rainncv(1,1), rainnc(1,1), flash_rate, flash_rate_cg, flash_rate_ic, &
               &                      snownc(1,1), snowncv(1,1))
        end if
      end do

      if (write_profiles) call output_close()

      block
        real(kind=RKIND), parameter :: snow_to_liquid_ratio = 10.0_RKIND  ! 10:1 by depth -- typical
                                                                           ! fresh/dry powder; denser
                                                                           ! wind-packed snow can be
                                                                           ! much lower (5:1 or less) --
                                                                           ! see README.md caveat
        real(kind=RKIND) :: snow_depth_cm
        snow_depth_cm = (snownc(1,1) / 10.0_RKIND) * snow_to_liquid_ratio  ! mm LWE -> cm LWE -> cm snow

        call ensemble_summary_write_member( &
             &  point_label(ip), satconfig_label(is), &
             &  qc_peak, qr_peak, qi_peak, qs_peak, qg_peak, &
             &  time_to_precip, rainnc(1,1), &
             &  flash_rate_peak, flash_rate_cg_peak, flash_rate_ic_peak, &
             &  flash_count_total, flash_count_cg_total, flash_count_ic_total, &
             &  snownc(1,1), snow_depth_cm)
      end block

      write(*,'(A,I0,A,A,A,A)') 'member ', member_count, ' (', trim(point_label(ip)), '/', &
           &                      trim(satconfig_label(is))//') done'

    end do

    call era5_close(efile)

  end do

  call ensemble_summary_close()

contains

  !===================================================================
  subroutine nearest_index(arr, n, target, idx)
  !===================================================================
    integer, intent(in) :: n
    real(kind=RKIND), intent(in) :: arr(n), target
    integer, intent(out) :: idx
    integer :: i
    real(kind=RKIND) :: best
    idx = 1
    best = abs(arr(1)-target)
    do i = 2, n
      if (abs(arr(i)-target) < best) then
        best = abs(arr(i)-target)
        idx = i
      end if
    end do
  end subroutine nearest_index

  !===================================================================
  subroutine bracket_layers(nz, z, z_lo, z_hi, k_lo, k_hi)
  !===================================================================
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: z(nz), z_lo, z_hi
    integer, intent(out) :: k_lo, k_hi
    integer :: k
    k_lo = 1; k_hi = nz
    do k = 1, nz
      if (z(k) >= z_lo) then
        k_lo = k
        exit
      end if
    end do
    do k = nz, 1, -1
      if (z(k) <= z_hi) then
        k_hi = k
        exit
      end if
    end do
    if (k_hi < k_lo) k_hi = k_lo
  end subroutine bracket_layers

  !===================================================================
  subroutine build_vertical_grid(nz, dz0, ztop, z, dz)
  !===================================================================
    ! identical stretching approach to mp_offline_driver.F90 -- see
    ! that file's own copy for the derivation comment
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: dz0, ztop
    real(kind=RKIND), intent(out) :: z(nz), dz(nz)
    real(kind=RKIND) :: ratio, lo, hi, mid, s
    real(kind=RKIND) :: zw(nz+1)
    integer :: k, iter
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
  subroutine regrid_era5_to_driver(col, nz, z, p_out, th_out, qv_out, qc_out, qi_out, w_out)
  !===================================================================
    ! era5_column_type is TOA-to-surface (col%z DEcreasing with index);
    ! flip to surface-up before the same linear-in-height interpolation
    ! goamazon_bridge.F90 uses, then convert T->theta.
    type(era5_column_type), intent(in) :: col
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: z(nz)
    real(kind=RKIND), intent(out) :: p_out(nz), th_out(nz), qv_out(nz), qc_out(nz), qi_out(nz), w_out(nz)
    real(kind=RKIND) :: zcol(col%nlev), pcol(col%nlev), tcol(col%nlev), qvcol(col%nlev), &
         &              qccol(col%nlev), qicol(col%nlev), wcol(col%nlev), thcol(col%nlev)
    integer :: k

    zcol  = col%z(col%nlev:1:-1);  pcol = col%p(col%nlev:1:-1)
    tcol  = col%t(col%nlev:1:-1);  qvcol = col%qv(col%nlev:1:-1)
    qccol = col%qc(col%nlev:1:-1); qicol = col%qi(col%nlev:1:-1)
    wcol  = col%w(col%nlev:1:-1)
    thcol = tcol * (p0/pcol)**rcp

    call interp1(col%nlev, zcol, pcol,  nz, z, p_out)
    call interp1(col%nlev, zcol, thcol, nz, z, th_out)
    call interp1(col%nlev, zcol, qvcol, nz, z, qv_out)
    call interp1(col%nlev, zcol, qccol, nz, z, qc_out)
    call interp1(col%nlev, zcol, qicol, nz, z, qi_out)
    call interp1(col%nlev, zcol, wcol,  nz, z, w_out)
  end subroutine regrid_era5_to_driver

  !===================================================================
  subroutine interp1(n, x, y, nout, xout, yout)
  !===================================================================
    integer, intent(in) :: n, nout
    real(kind=RKIND), intent(in) :: x(n), y(n), xout(nout)
    real(kind=RKIND), intent(out) :: yout(nout)
    integer :: k, j
    real(kind=RKIND) :: wgt
    do k = 1, nout
      if (xout(k) <= x(1)) then
        yout(k) = y(1)
      else if (xout(k) >= x(n)) then
        yout(k) = y(n)
      else
        j = 1
        do while (x(j+1) < xout(k))
          j = j + 1
        end do
        wgt = (xout(k)-x(j)) / (x(j+1)-x(j))
        yout(k) = (1.0_RKIND-wgt)*y(j) + wgt*y(j+1)
      end if
    end do
  end subroutine interp1

  !===================================================================
  subroutine apply_vertical_advection(nz, dz, dt, w, th, qv)
  !===================================================================
    ! same as mp_offline_driver.F90's apply_forcing_and_advance, minus
    ! the horizontal-advective-tendency terms (none here -- see header)
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: dz(nz), dt, w(nz)
    real(kind=RKIND), intent(inout) :: th(nz), qv(nz)
    real(kind=RKIND) :: dthdz(nz), dqvdz(nz)
    integer :: k
    do k = 2, nz-1
      dthdz(k) = (th(k+1)-th(k-1)) / (dz(k)+dz(k-1))
      dqvdz(k) = (qv(k+1)-qv(k-1)) / (dz(k)+dz(k-1))
    end do
    dthdz(1)  = (th(2)-th(1))   / dz(1)
    dqvdz(1)  = (qv(2)-qv(1))   / dz(1)
    dthdz(nz) = (th(nz)-th(nz-1)) / dz(nz-1)
    dqvdz(nz) = (qv(nz)-qv(nz-1)) / dz(nz-1)
    do k = 1, nz
      th(k) = th(k) - dt*w(k)*dthdz(k)
      qv(k) = max(qv(k) - dt*w(k)*dqvdz(k), 0.0_RKIND)
    end do
  end subroutine apply_vertical_advection

end program ensemble_driver
