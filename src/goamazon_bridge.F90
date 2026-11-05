! goamazon_bridge.F90
!
! Reads the DEPHY-SCM common-format forcing file for the GOAMAZON/SINGLE
! case (Amazon deep convection, 2014-10-05, Tian & Zhang 2025 GRL) and
! builds the plain Fortran column arrays mp_gt_driver() expects at one
! instant in time -- i.e. this is the offline, GOAMAZON-driven
! replacement for MPAS's mpas_atmphys_interface.F "*_from_MPAS" packing
! subroutines, for the microphysics-only case (no dynamical core, no
! other physics parameterizations).
!
! Expected input file (GdR-DEPHY/DEPHY-SCM repository, file
! GOAMAZON/SINGLE/GOAMAZON_SINGLE_SCM_driver.nc), DEPHY common format
! v1 (see DEPHY-SCM_CommonFormat_v1.0.pdf in that repository):
!
!   Dimensions: t0(1), time(ntime), lev(nlev)
!
!   lev is NOT a generic index: it IS height above ground in metres,
!   a fixed common grid shared by every DEPHY-SCM case (confirmed for
!   GOAMAZON: lev(1)=0, lev(2)=10, ..., lev(2001)=20000 -- 10 m spacing
!   to 20 km). Both initial-condition variables (dims t0,lev) and
!   forcing variables (dims time,lev) are given on this SAME fixed
!   height grid; pressure itself is one of the time-varying forcing
!   fields (pa_forc), not a coordinate.
!
!   Initial condition (t0,lev): pa, theta, rv, rl, ri, ua, va
!   Forcing (time,lev):         pa_forc, tnta_adv, tnrv_adv, tnrt_adv, wa
!   Forcing (time):             hfss, hfls, ps_forc
!
!   This case has forc_wa=1 (real prescribed large-scale vertical
!   velocity), adv_ta=1/adv_rv=1/adv_rt=1 (real prescribed advective
!   tendencies) -- unlike LBA/REF, which has none of these and is
!   meant to drive an explicit CRM/LES instead of a forced column.
!
! Working variables chosen in MIXING-RATIO basis (rv/rl/ri/rt, not the
! specific-humidity qv/ql/qi/qt also present in the file) because that
! is what WRF/MPAS-convention microphysics (Thompson included) expects
! for qv/qc/qi. rl (liquid) is used to initialize qc; rr/qs/qg have no
! DEPHY equivalent and start at zero, same as ni/nr.
!
! FORCING APPLICATION CONVENTION (this is a modeling choice, not
! something read from the file -- see README.md "Status" section):
! following the standard GCSS/DEPHY "revealed forcing" convention for
! forc_wa=1 cases, tnta_adv/tnrv_adv/tnrt_adv represent HORIZONTAL
! advection only; vertical advection is added separately by the driver
! each timestep as -wa*d(.)/dz using the column's OWN current profile
! gradient (see mp_offline_driver.F90). This module supplies both
! pieces (the horizontal tendency and wa itself, both already
! interpolated to the driver's vertical grid) and leaves that
! combination to the driver.
!
! Vertical interpolation: linear in height, from the file's native
! 10 m grid onto the driver's own (coarser, stretched) grid. Time
! interpolation: linear between the two bracketing forcing times;
! t_sec before the first or after the last forcing time is clamped to
! the nearest end (no extrapolation).

module goamazon_bridge

  use mpas_kind_types, only : RKIND
  use netcdf

  implicit none
  private

  public :: goamazon_file_type
  public :: goamazon_open, goamazon_close
  public :: goamazon_get_initial_profile
  public :: goamazon_get_forcing_at_time

  type :: goamazon_file_type
    integer :: ncid  = -1
    integer :: nlev  = 0
    integer :: ntime = 0
    real(kind=RKIND), allocatable :: lev(:)          ! (nlev) fixed height grid, m
    real(kind=RKIND), allocatable :: time(:)         ! (ntime) forcing times, s since t0

    ! initial condition, t0, on the native height grid:
    real(kind=RKIND), allocatable :: pa0(:), theta0(:), rv0(:), rl0(:), ri0(:), ua0(:), va0(:)
    real(kind=RKIND) :: ps0 = 0.0_RKIND

    ! forcing, native height grid x time:
    real(kind=RKIND), allocatable :: pa_forc(:,:)    ! (nlev,ntime) Pa
    real(kind=RKIND), allocatable :: tnta_adv(:,:)   ! (nlev,ntime) K/s
    real(kind=RKIND), allocatable :: tnrv_adv(:,:)   ! (nlev,ntime) kg/kg/s
    real(kind=RKIND), allocatable :: tnrt_adv(:,:)   ! (nlev,ntime) kg/kg/s
    real(kind=RKIND), allocatable :: wa(:,:)         ! (nlev,ntime) m/s
    real(kind=RKIND), allocatable :: hfss(:), hfls(:), ps_forc(:)  ! (ntime)

    ! surface boundary condition + wind, for the MYNN surface-layer/PBL coupling:
    real(kind=RKIND), allocatable :: ts_forc(:)      ! (ntime) K, skin temperature
    real(kind=RKIND), allocatable :: z0(:)           ! (ntime) m, roughness length
    real(kind=RKIND), allocatable :: ua_nud(:,:)     ! (nlev,ntime) m/s, prescribed wind (see README.md:
    real(kind=RKIND), allocatable :: va_nud(:,:)     ! (nlev,ntime) m/s  used directly, not actually "nudged")
  end type goamazon_file_type

contains

  subroutine nc_check(status, where)
    integer, intent(in) :: status
    character(len=*), intent(in) :: where
    if (status /= nf90_noerr) then
      write(0,'(A)') 'goamazon_bridge: NetCDF error at '//trim(where)//': '// &
           &          trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine nc_check

  !=================================================================
  subroutine goamazon_open(path, f, legacy_temperature_tendency)
  !=================================================================
    character(len=*), intent(in)  :: path
    logical, intent(in), optional :: legacy_temperature_tendency
    type(goamazon_file_type), intent(out) :: f

    integer :: dimid, varid

    call nc_check(nf90_open(trim(path), nf90_nowrite, f%ncid), 'open '//trim(path))

    call nc_check(nf90_inq_dimid(f%ncid, 'lev', dimid), 'inq dim lev')
    call nc_check(nf90_inquire_dimension(f%ncid, dimid, len=f%nlev), 'inquire lev')
    call nc_check(nf90_inq_dimid(f%ncid, 'time', dimid), 'inq dim time')
    call nc_check(nf90_inquire_dimension(f%ncid, dimid, len=f%ntime), 'inquire time')

    allocate(f%lev(f%nlev), f%time(f%ntime))
    allocate(f%pa0(f%nlev), f%theta0(f%nlev), f%rv0(f%nlev), f%rl0(f%nlev), &
             f%ri0(f%nlev), f%ua0(f%nlev), f%va0(f%nlev))
    allocate(f%pa_forc(f%nlev,f%ntime), f%tnta_adv(f%nlev,f%ntime), &
             f%tnrv_adv(f%nlev,f%ntime), f%tnrt_adv(f%nlev,f%ntime), &
             f%wa(f%nlev,f%ntime))
    allocate(f%hfss(f%ntime), f%hfls(f%ntime), f%ps_forc(f%ntime))
    allocate(f%ts_forc(f%ntime), f%z0(f%ntime))
    allocate(f%ua_nud(f%nlev,f%ntime), f%va_nud(f%nlev,f%ntime))

    call nc_check(nf90_inq_varid(f%ncid, 'lev', varid), 'inq lev')
    call nc_check(nf90_get_var(f%ncid, varid, f%lev), 'get lev')
    call nc_check(nf90_inq_varid(f%ncid, 'time', varid), 'inq time')
    call nc_check(nf90_get_var(f%ncid, varid, f%time), 'get time')

    ! initial condition (dims t0,lev in the file -- t0 has length 1, so
    ! a plain (nlev) read works directly)
    call read_1d_t0lev(f%ncid, 'pa',    f%nlev, f%pa0)
    call read_1d_t0lev(f%ncid, 'theta', f%nlev, f%theta0)
    call read_1d_t0lev(f%ncid, 'rv',    f%nlev, f%rv0)
    call read_1d_t0lev(f%ncid, 'rl',    f%nlev, f%rl0)
    call read_1d_t0lev(f%ncid, 'ri',    f%nlev, f%ri0)
    call read_1d_t0lev(f%ncid, 'ua',    f%nlev, f%ua0)
    call read_1d_t0lev(f%ncid, 'va',    f%nlev, f%va0)
    call nc_check(nf90_inq_varid(f%ncid, 'ps', varid), 'inq ps')
    call nc_check(nf90_get_var(f%ncid, varid, f%ps0), 'get ps')

    ! forcing (dims time,lev in the file -> (nlev,ntime) in Fortran's
    ! reversed storage order, read directly with no transpose needed)
    call read_2d_timelev(f%ncid, 'pa_forc',   f%nlev, f%ntime, f%pa_forc)
    ! Integrator advances theta, not T. Legacy mode reproduces old interface.
    if (present(legacy_temperature_tendency)) then
      if (legacy_temperature_tendency) then
        call read_2d_timelev(f%ncid, 'tnta_adv', f%nlev, f%ntime, f%tnta_adv)
      else
        call read_2d_timelev(f%ncid, 'tntheta_adv', f%nlev, f%ntime, f%tnta_adv)
      end if
    else
      call read_2d_timelev(f%ncid, 'tntheta_adv', f%nlev, f%ntime, f%tnta_adv)
    end if
    call read_2d_timelev(f%ncid, 'tnrv_adv',  f%nlev, f%ntime, f%tnrv_adv)
    call read_2d_timelev(f%ncid, 'tnrt_adv',  f%nlev, f%ntime, f%tnrt_adv)
    call read_2d_timelev(f%ncid, 'wa',        f%nlev, f%ntime, f%wa)
    call read_2d_timelev(f%ncid, 'ua_nud',    f%nlev, f%ntime, f%ua_nud)
    call read_2d_timelev(f%ncid, 'va_nud',    f%nlev, f%ntime, f%va_nud)

    call nc_check(nf90_inq_varid(f%ncid, 'ts_forc', varid), 'inq ts_forc')
    call nc_check(nf90_get_var(f%ncid, varid, f%ts_forc), 'get ts_forc')
    call nc_check(nf90_inq_varid(f%ncid, 'z0', varid), 'inq z0')
    call nc_check(nf90_get_var(f%ncid, varid, f%z0), 'get z0')

    call nc_check(nf90_inq_varid(f%ncid, 'hfss', varid), 'inq hfss')
    call nc_check(nf90_get_var(f%ncid, varid, f%hfss), 'get hfss')
    call nc_check(nf90_inq_varid(f%ncid, 'hfls', varid), 'inq hfls')
    call nc_check(nf90_get_var(f%ncid, varid, f%hfls), 'get hfls')
    call nc_check(nf90_inq_varid(f%ncid, 'ps_forc', varid), 'inq ps_forc')
    call nc_check(nf90_get_var(f%ncid, varid, f%ps_forc), 'get ps_forc')

  end subroutine goamazon_open

  !-----------------------------------------------------------------
  subroutine read_1d_t0lev(ncid, name, nlev, arr)
  !-----------------------------------------------------------------
    integer, intent(in) :: ncid, nlev
    character(len=*), intent(in) :: name
    real(kind=RKIND), intent(out) :: arr(nlev)
    integer :: varid
    call nc_check(nf90_inq_varid(ncid, name, varid), 'inq '//trim(name))
    call nc_check(nf90_get_var(ncid, varid, arr, start=(/1,1/), count=(/nlev,1/)), &
         &        'get '//trim(name))
  end subroutine read_1d_t0lev

  !-----------------------------------------------------------------
  subroutine read_2d_timelev(ncid, name, nlev, ntime, arr)
  !-----------------------------------------------------------------
    integer, intent(in) :: ncid, nlev, ntime
    character(len=*), intent(in) :: name
    real(kind=RKIND), intent(out) :: arr(nlev,ntime)
    integer :: varid
    call nc_check(nf90_inq_varid(ncid, name, varid), 'inq '//trim(name))
    ! file storage order is (time,lev); nf90_get_var into a Fortran
    ! array declared (nlev,ntime) with matching start/count gives the
    ! REVERSE-order read netCDF-Fortran always does (fastest-varying
    ! netCDF dim -> first Fortran dim), i.e. exactly what we want with
    ! no transpose.
    call nc_check(nf90_get_var(ncid, varid, arr, start=(/1,1/), count=(/nlev,ntime/)), &
         &        'get '//trim(name))
  end subroutine read_2d_timelev

  !=================================================================
  subroutine goamazon_close(f)
  !=================================================================
    type(goamazon_file_type), intent(inout) :: f
    if (f%ncid >= 0) call nc_check(nf90_close(f%ncid), 'close')
    f%ncid = -1
  end subroutine goamazon_close

  !=================================================================
  ! Linear interpolation of a (nlev)-on-f%lev profile onto an
  ! arbitrary target height grid z_drv(1:nz), clamped (no
  ! extrapolation) at the ends -- both ends are physically safe here
  ! since z_drv always lies inside [0,20000] m by construction
  ! (mp_offline_driver.F90 builds z_drv that way).
  !
  ! lev(1:nlev) only needs to be MONOTONICALLY INCREASING, not
  ! uniformly spaced: the bracketing level j is found by bisection.
  ! An earlier version of this routine assumed uniform 10 m spacing
  ! (computing j directly from lev(2)-lev(1)), which is true of the
  ! original GOAMAZON_SINGLE_SCM_driver.nc (2001 levels, 10 m grid)
  ! but is NOT true of a DEPHY-SCM file built directly on ERA5's own
  ! pressure-level-derived heights (e.g. driver_DEF_era5.py's
  ! GOAMAZON_SINGLE_DEF_driver.nc, ~30-40 irregularly-spaced levels)
  ! -- on such a file the old code silently returned the WRONG
  ! bracket and wrong interpolated values, with no error of any
  ! kind. Bisection is correct for both cases and is not measurably
  ! slower at these array sizes.
  !=================================================================
  subroutine interp_z(nlev, lev, prof, nz, z_drv, out)
    integer, intent(in) :: nlev, nz
    real(kind=RKIND), intent(in) :: lev(nlev), prof(nlev), z_drv(nz)
    real(kind=RKIND), intent(out) :: out(nz)
    integer :: k, j, lo, hi, mid
    real(kind=RKIND) :: wgt

    do k = 1, nz
      if (z_drv(k) <= lev(1)) then
        out(k) = prof(1)
      else if (z_drv(k) >= lev(nlev)) then
        out(k) = prof(nlev)
      else
        ! bisection: find j such that lev(j) <= z_drv(k) < lev(j+1)
        lo = 1; hi = nlev
        do while (hi - lo > 1)
          mid = (lo + hi) / 2
          if (lev(mid) <= z_drv(k)) then
            lo = mid
          else
            hi = mid
          end if
        end do
        j = max(1, min(lo, nlev-1))
        wgt = (z_drv(k) - lev(j)) / (lev(j+1) - lev(j))
        out(k) = (1.0_RKIND - wgt) * prof(j) + wgt * prof(j+1)
      end if
    end do
  end subroutine interp_z

  !=================================================================
  subroutine goamazon_get_initial_profile(f, nz, z_drv, pa, theta, rv, rl, ri, ua, va, ps0)
  !=================================================================
    type(goamazon_file_type), intent(in) :: f
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in)  :: z_drv(nz)
    real(kind=RKIND), intent(out) :: pa(nz), theta(nz), rv(nz), rl(nz), ri(nz), ua(nz), va(nz)
    real(kind=RKIND), intent(out) :: ps0

    call interp_z(f%nlev, f%lev, f%pa0,    nz, z_drv, pa)
    call interp_z(f%nlev, f%lev, f%theta0, nz, z_drv, theta)
    call interp_z(f%nlev, f%lev, f%rv0,    nz, z_drv, rv)
    call interp_z(f%nlev, f%lev, f%rl0,    nz, z_drv, rl)
    call interp_z(f%nlev, f%lev, f%ri0,    nz, z_drv, ri)
    call interp_z(f%nlev, f%lev, f%ua0,    nz, z_drv, ua)
    call interp_z(f%nlev, f%lev, f%va0,    nz, z_drv, va)
    ps0 = f%ps0
  end subroutine goamazon_get_initial_profile

  !=================================================================
  ! Forcing at one instant t_sec (seconds since t0): linear in time
  ! between the two bracketing forcing times (clamped at the ends),
  ! THEN linear in height onto z_drv -- i.e. time-interpolate on the
  ! native height grid first, so the two interpolations don't need to
  ! commute (they do here, since both are linear, but doing it in this
  ! order keeps the height interpolation weights reusable across all
  ! five fields for one call).
  !=================================================================
  subroutine goamazon_get_forcing_at_time(f, t_sec, nz, z_drv, &
       &                                  pa_forc, ta_adv, rv_adv, rt_adv, w, &
       &                                  hfss_t, hfls_t, ps_t, ts_t, z0_t, u_col, v_col)
    type(goamazon_file_type), intent(in) :: f
    real(kind=RKIND), intent(in) :: t_sec
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in)  :: z_drv(nz)
    real(kind=RKIND), intent(out) :: pa_forc(nz), ta_adv(nz), rv_adv(nz), rt_adv(nz), w(nz)
    real(kind=RKIND), intent(out) :: hfss_t, hfls_t, ps_t, ts_t, z0_t
    real(kind=RKIND), intent(out) :: u_col(nz), v_col(nz)

    integer :: it0, it1, it
    real(kind=RKIND) :: wgt
    real(kind=RKIND) :: prof_pa(f%nlev), prof_ta(f%nlev), prof_rv(f%nlev), &
                         prof_rt(f%nlev), prof_w(f%nlev), prof_u(f%nlev), prof_v(f%nlev)

    ! bracket t_sec in f%time (clamped, no extrapolation)
    if (t_sec <= f%time(1)) then
      it0 = 1; it1 = 1; wgt = 0.0_RKIND
    else if (t_sec >= f%time(f%ntime)) then
      it0 = f%ntime; it1 = f%ntime; wgt = 0.0_RKIND
    else
      it0 = 1
      do it = 1, f%ntime-1
        if (f%time(it+1) > t_sec) then
          it0 = it
          exit
        end if
      end do
      it1 = it0 + 1
      wgt = (t_sec - f%time(it0)) / (f%time(it1) - f%time(it0))
    end if

    prof_pa = (1.0_RKIND-wgt)*f%pa_forc(:,it0)  + wgt*f%pa_forc(:,it1)
    prof_ta = (1.0_RKIND-wgt)*f%tnta_adv(:,it0) + wgt*f%tnta_adv(:,it1)
    prof_rv = (1.0_RKIND-wgt)*f%tnrv_adv(:,it0) + wgt*f%tnrv_adv(:,it1)
    prof_rt = (1.0_RKIND-wgt)*f%tnrt_adv(:,it0) + wgt*f%tnrt_adv(:,it1)
    prof_w  = (1.0_RKIND-wgt)*f%wa(:,it0)       + wgt*f%wa(:,it1)
    prof_u  = (1.0_RKIND-wgt)*f%ua_nud(:,it0)   + wgt*f%ua_nud(:,it1)
    prof_v  = (1.0_RKIND-wgt)*f%va_nud(:,it0)   + wgt*f%va_nud(:,it1)

    call interp_z(f%nlev, f%lev, prof_pa, nz, z_drv, pa_forc)
    call interp_z(f%nlev, f%lev, prof_ta, nz, z_drv, ta_adv)
    call interp_z(f%nlev, f%lev, prof_rv, nz, z_drv, rv_adv)
    call interp_z(f%nlev, f%lev, prof_rt, nz, z_drv, rt_adv)
    call interp_z(f%nlev, f%lev, prof_w,  nz, z_drv, w)
    call interp_z(f%nlev, f%lev, prof_u,  nz, z_drv, u_col)
    call interp_z(f%nlev, f%lev, prof_v,  nz, z_drv, v_col)

    hfss_t = (1.0_RKIND-wgt)*f%hfss(it0)    + wgt*f%hfss(it1)
    hfls_t = (1.0_RKIND-wgt)*f%hfls(it0)    + wgt*f%hfls(it1)
    ps_t   = (1.0_RKIND-wgt)*f%ps_forc(it0) + wgt*f%ps_forc(it1)
    ts_t   = (1.0_RKIND-wgt)*f%ts_forc(it0) + wgt*f%ts_forc(it1)
    z0_t   = (1.0_RKIND-wgt)*f%z0(it0)      + wgt*f%z0(it1)

  end subroutine goamazon_get_forcing_at_time

end module goamazon_bridge
