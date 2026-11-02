! output_mp_netcdf.F90
!
! Writes the mp_offline_driver column time series to a single NetCDF
! file: dimensions time(unlimited), lev(nz). Profile fields (theta,
! qv, qc, qr, qi, qs, qg, pressure) plus two surface scalars per step
! (incremental and accumulated surface precipitation, RAINNCV/RAINNC
! straight from mp_gt_driver, mm). Minimal by design -- mirrors
! output_netcdf.F90 from the companion radiation project only in
! spirit (open/write_timestep/close pattern), not in file layout,
! since the two drivers' outputs share no fields.

module output_mp_netcdf

  use mpas_kind_types, only : RKIND
  use netcdf

  implicit none
  private

  public :: output_open, output_write_timestep, output_close

  integer :: ncid, dimid_time, dimid_lev
  integer :: id_time, id_lev, id_theta, id_qv, id_qc, id_qr, id_qi, id_qs, id_qg, &
       &     id_pres, id_rainncv, id_rainnc, id_flashrate, id_flashrate_cg, id_flashrate_ic, &
       &     id_snownc, id_snowncv
  integer :: id_tke, id_mixlen, id_exch_h, id_exch_m, id_cldfra_bl, id_qc_bl, id_qi_bl, &
       &     id_mf_a, id_mf_w, id_pblh
  integer :: itime = 0

contains

  subroutine nc_check(status, where)
    integer, intent(in) :: status
    character(len=*), intent(in) :: where
    if (status /= nf90_noerr) then
      write(0,'(A)') 'output_mp_netcdf: NetCDF error at '//trim(where)//': '// &
           &          trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine nc_check

  !=================================================================
  subroutine output_open(path, nz, z)
  !=================================================================
    character(len=*), intent(in) :: path
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: z(nz)

    call nc_check(nf90_create(trim(path), nf90_clobber, ncid), 'create '//trim(path))

    call nc_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_time), 'def dim time')
    call nc_check(nf90_def_dim(ncid, 'lev',  nz,             dimid_lev),  'def dim lev')

    call nc_check(nf90_def_var(ncid, 'time', nf90_double, (/dimid_time/), id_time), 'def time')
    call nc_check(nf90_put_att(ncid, id_time, 'units', 'seconds since t0'), 'att time')

    call nc_check(nf90_def_var(ncid, 'lev', nf90_double, (/dimid_lev/), id_lev), 'def lev')
    call nc_check(nf90_put_att(ncid, id_lev, 'units', 'm'), 'att lev')
    call nc_check(nf90_put_att(ncid, id_lev, 'long_name', 'driver layer-center height'), 'att lev longname')

    call def_profile_var('theta', 'K',     'potential temperature', id_theta)
    call def_profile_var('qv',    'kg kg-1','water vapor mixing ratio', id_qv)
    call def_profile_var('qc',    'kg kg-1','cloud water mixing ratio', id_qc)
    call def_profile_var('qr',    'kg kg-1','rain mixing ratio', id_qr)
    call def_profile_var('qi',    'kg kg-1','cloud ice mixing ratio', id_qi)
    call def_profile_var('qs',    'kg kg-1','snow mixing ratio', id_qs)
    call def_profile_var('qg',    'kg kg-1','graupel mixing ratio', id_qg)
    call def_profile_var('pressure', 'Pa', 'pressure', id_pres)

    call nc_check(nf90_def_var(ncid, 'rainncv', nf90_double, (/dimid_time/), id_rainncv), 'def rainncv')
    call nc_check(nf90_put_att(ncid, id_rainncv, 'units', 'mm'), 'att rainncv')
    call nc_check(nf90_put_att(ncid, id_rainncv, 'long_name', &
         &        'incremental surface precip this timestep (mp_gt_driver RAINNCV)'), 'att rainncv longname')

    call nc_check(nf90_def_var(ncid, 'rainnc', nf90_double, (/dimid_time/), id_rainnc), 'def rainnc')
    call nc_check(nf90_put_att(ncid, id_rainnc, 'units', 'mm'), 'att rainnc')
    call nc_check(nf90_put_att(ncid, id_rainnc, 'long_name', &
         &        'accumulated surface precip since t0 (mp_gt_driver RAINNC)'), 'att rainnc longname')

    call nc_check(nf90_def_var(ncid, 'flash_rate', nf90_double, (/dimid_time/), id_flashrate), 'def flash_rate')
    call nc_check(nf90_put_att(ncid, id_flashrate, 'units', 'flashes min-1'), 'att flash_rate units')
    call nc_check(nf90_put_att(ncid, id_flashrate, 'long_name', &
         &        'diagnostic total flash rate (mpas_atmphys_lightning, driver_lightning)'), 'att flash_rate longname')
    call nc_check(nf90_def_var(ncid, 'flash_rate_cg', nf90_double, (/dimid_time/), id_flashrate_cg), 'def flash_rate_cg')
    call nc_check(nf90_put_att(ncid, id_flashrate_cg, 'units', 'flashes min-1'), 'att flash_rate_cg units')
    call nc_check(nf90_def_var(ncid, 'flash_rate_ic', nf90_double, (/dimid_time/), id_flashrate_ic), 'def flash_rate_ic')
    call nc_check(nf90_put_att(ncid, id_flashrate_ic, 'units', 'flashes min-1'), 'att flash_rate_ic units')

    call nc_check(nf90_def_var(ncid, 'snownc', nf90_double, (/dimid_time/), id_snownc), 'def snownc')
    call nc_check(nf90_put_att(ncid, id_snownc, 'units', 'mm (liquid-water equivalent)'), 'att snownc units')
    call nc_check(nf90_put_att(ncid, id_snownc, 'long_name', &
         &        'accumulated surface snow+cloud-ice fallout since t0, LIQUID-WATER EQUIVALENT '// &
         &        '(NOT physical snow depth -- see README.md for the depth conversion caveat)'), 'att snownc longname')
    call nc_check(nf90_def_var(ncid, 'snowncv', nf90_double, (/dimid_time/), id_snowncv), 'def snowncv')
    call nc_check(nf90_put_att(ncid, id_snowncv, 'units', 'mm (liquid-water equivalent)'), 'att snowncv units')

    ! ---- MYNN turbulence/mass-flux/cloud-fraction diagnostics (only
    !      meaningful when use_mynn=.true.; zero-filled otherwise, see
    !      output_write_timestep) ----
    call def_profile_var('tke',          'm2 s-2',  'turbulent kinetic energy (MYNN, qke/2)', id_tke)
    call def_profile_var('mixing_length','m',        'MYNN master mixing length (el_pbl)', id_mixlen)
    call def_profile_var('exch_h',       'm2 s-1',  'eddy diffusivity, heat (MYNN)', id_exch_h)
    call def_profile_var('exch_m',       'm2 s-1',  'eddy diffusivity, momentum (MYNN)', id_exch_m)
    call def_profile_var('cldfra_bl',    '0-1',      'MYNN-diagnosed subgrid-scale (PBL) cloud fraction', id_cldfra_bl)
    call def_profile_var('qc_bl',        'kg kg-1', 'MYNN-diagnosed subgrid-scale cloud water', id_qc_bl)
    call def_profile_var('qi_bl',        'kg kg-1', 'MYNN-diagnosed subgrid-scale cloud ice', id_qi_bl)
    call def_profile_var('mf_a',         '0-1',      'EDMF updraft area fraction (MYNN mass-flux scheme)', id_mf_a)
    call def_profile_var('mf_w',         'm s-1',   'EDMF updraft vertical velocity (MYNN mass-flux scheme)', id_mf_w)

    call nc_check(nf90_def_var(ncid, 'pblh', nf90_double, (/dimid_time/), id_pblh), 'def pblh')
    call nc_check(nf90_put_att(ncid, id_pblh, 'units', 'm'), 'att pblh units')
    call nc_check(nf90_put_att(ncid, id_pblh, 'long_name', 'MYNN-diagnosed boundary-layer height'), 'att pblh longname')

    call nc_check(nf90_put_att(ncid, nf90_global, 'source', &
         &        'mp_offline_driver: MPAS Thompson microphysics (module_mp_thompson.F, unmodified) '// &
         &        'forced by DEPHY-SCM GOAMAZON/SINGLE'), 'att global source')

    call nc_check(nf90_enddef(ncid), 'enddef')
    call nc_check(nf90_put_var(ncid, id_lev, z), 'put lev')

    itime = 0
  end subroutine output_open

  !-----------------------------------------------------------------
  subroutine def_profile_var(name, units, long_name, varid)
  !-----------------------------------------------------------------
    character(len=*), intent(in) :: name, units, long_name
    integer, intent(out) :: varid
    call nc_check(nf90_def_var(ncid, name, nf90_double, (/dimid_lev, dimid_time/), varid), &
         &        'def '//trim(name))
    call nc_check(nf90_put_att(ncid, varid, 'units', units), 'att units '//trim(name))
    call nc_check(nf90_put_att(ncid, varid, 'long_name', long_name), 'att longname '//trim(name))
  end subroutine def_profile_var

  !=================================================================
  subroutine output_write_timestep(t_sec, nz, theta, qv, qc, qr, qi, qs, qg, pres, &
       &                           rainncv_t, rainnc_t, flash_rate_t, flash_rate_cg_t, flash_rate_ic_t, &
       &                           snownc_t, snowncv_t, &
       &                           tke_t, mixing_length_t, exch_h_t, exch_m_t, &
       &                           cldfra_bl_t, qc_bl_t, qi_bl_t, mf_a_t, mf_w_t, pblh_t)
  !=================================================================
    real(kind=RKIND), intent(in) :: t_sec
    integer, intent(in) :: nz
    real(kind=RKIND), intent(in) :: theta(nz), qv(nz), qc(nz), qr(nz), qi(nz), qs(nz), qg(nz), pres(nz)
    real(kind=RKIND), intent(in) :: rainncv_t, rainnc_t
    real(kind=RKIND), intent(in), optional :: flash_rate_t, flash_rate_cg_t, flash_rate_ic_t
    real(kind=RKIND), intent(in), optional :: snownc_t, snowncv_t
    real(kind=RKIND), intent(in), optional :: tke_t(nz), mixing_length_t(nz), exch_h_t(nz), exch_m_t(nz)
    real(kind=RKIND), intent(in), optional :: cldfra_bl_t(nz), qc_bl_t(nz), qi_bl_t(nz)
    real(kind=RKIND), intent(in), optional :: mf_a_t(nz), mf_w_t(nz)
    real(kind=RKIND), intent(in), optional :: pblh_t

    itime = itime + 1

    call nc_check(nf90_put_var(ncid, id_time, t_sec, start=(/itime/)), 'put time')
    call nc_check(nf90_put_var(ncid, id_theta, theta, start=(/1,itime/), count=(/nz,1/)), 'put theta')
    call nc_check(nf90_put_var(ncid, id_qv,    qv,    start=(/1,itime/), count=(/nz,1/)), 'put qv')
    call nc_check(nf90_put_var(ncid, id_qc,    qc,    start=(/1,itime/), count=(/nz,1/)), 'put qc')
    call nc_check(nf90_put_var(ncid, id_qr,    qr,    start=(/1,itime/), count=(/nz,1/)), 'put qr')
    call nc_check(nf90_put_var(ncid, id_qi,    qi,    start=(/1,itime/), count=(/nz,1/)), 'put qi')
    call nc_check(nf90_put_var(ncid, id_qs,    qs,    start=(/1,itime/), count=(/nz,1/)), 'put qs')
    call nc_check(nf90_put_var(ncid, id_qg,    qg,    start=(/1,itime/), count=(/nz,1/)), 'put qg')
    call nc_check(nf90_put_var(ncid, id_pres,  pres,  start=(/1,itime/), count=(/nz,1/)), 'put pres')
    call nc_check(nf90_put_var(ncid, id_rainncv, rainncv_t, start=(/itime/)), 'put rainncv')
    call nc_check(nf90_put_var(ncid, id_rainnc,  rainnc_t,  start=(/itime/)), 'put rainnc')

    if (present(flash_rate_t)) then
      call nc_check(nf90_put_var(ncid, id_flashrate, flash_rate_t, start=(/itime/)), 'put flash_rate')
    else
      call nc_check(nf90_put_var(ncid, id_flashrate, 0.0_RKIND, start=(/itime/)), 'put flash_rate')
    end if
    if (present(flash_rate_cg_t)) then
      call nc_check(nf90_put_var(ncid, id_flashrate_cg, flash_rate_cg_t, start=(/itime/)), 'put flash_rate_cg')
    else
      call nc_check(nf90_put_var(ncid, id_flashrate_cg, 0.0_RKIND, start=(/itime/)), 'put flash_rate_cg')
    end if
    if (present(flash_rate_ic_t)) then
      call nc_check(nf90_put_var(ncid, id_flashrate_ic, flash_rate_ic_t, start=(/itime/)), 'put flash_rate_ic')
    else
      call nc_check(nf90_put_var(ncid, id_flashrate_ic, 0.0_RKIND, start=(/itime/)), 'put flash_rate_ic')
    end if

    if (present(snownc_t)) then
      call nc_check(nf90_put_var(ncid, id_snownc, snownc_t, start=(/itime/)), 'put snownc')
    else
      call nc_check(nf90_put_var(ncid, id_snownc, 0.0_RKIND, start=(/itime/)), 'put snownc')
    end if
    if (present(snowncv_t)) then
      call nc_check(nf90_put_var(ncid, id_snowncv, snowncv_t, start=(/itime/)), 'put snowncv')
    else
      call nc_check(nf90_put_var(ncid, id_snowncv, 0.0_RKIND, start=(/itime/)), 'put snowncv')
    end if

    ! MYNN turbulence/mass-flux/cloud-fraction diagnostics -- zero-filled
    ! (not left undefined) when use_mynn=.false., so the variable is
    ! always present and unambiguous to read regardless of run config
    call put_profile_opt(id_tke,         nz, tke_t)
    call put_profile_opt(id_mixlen,      nz, mixing_length_t)
    call put_profile_opt(id_exch_h,      nz, exch_h_t)
    call put_profile_opt(id_exch_m,      nz, exch_m_t)
    call put_profile_opt(id_cldfra_bl,   nz, cldfra_bl_t)
    call put_profile_opt(id_qc_bl,       nz, qc_bl_t)
    call put_profile_opt(id_qi_bl,       nz, qi_bl_t)
    call put_profile_opt(id_mf_a,        nz, mf_a_t)
    call put_profile_opt(id_mf_w,        nz, mf_w_t)

    if (present(pblh_t)) then
      call nc_check(nf90_put_var(ncid, id_pblh, pblh_t, start=(/itime/)), 'put pblh')
    else
      call nc_check(nf90_put_var(ncid, id_pblh, 0.0_RKIND, start=(/itime/)), 'put pblh')
    end if

  end subroutine output_write_timestep

  !-----------------------------------------------------------------
  subroutine put_profile_opt(varid, nz, prof)
  !-----------------------------------------------------------------
    ! writes prof at the current itime if present, else a zero profile
    integer, intent(in) :: varid, nz
    real(kind=RKIND), intent(in), optional :: prof(nz)
    real(kind=RKIND) :: zeros(nz)
    if (present(prof)) then
      call nc_check(nf90_put_var(ncid, varid, prof, start=(/1,itime/), count=(/nz,1/)), 'put profile')
    else
      zeros = 0.0_RKIND
      call nc_check(nf90_put_var(ncid, varid, zeros, start=(/1,itime/), count=(/nz,1/)), 'put profile (zero-filled)')
    end if
  end subroutine put_profile_opt

  !=================================================================
  subroutine output_close()
  !=================================================================
    call nc_check(nf90_close(ncid), 'close output')
  end subroutine output_close

end module output_mp_netcdf
