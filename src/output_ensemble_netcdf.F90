! output_ensemble_netcdf.F90
!
! Writes the ensemble_driver summary (one row per member: point label,
! saturation-config label, peak hydrometeor contents, precipitation,
! lightning and snow diagnostics) to a single NetCDF file, replacing
! the earlier CSV format -- same dimensioned/attributed convention as
! output_mp_netcdf.F90 and the companion rad_offline project's
! output_netcdf.F90 (NF90_CLOBBER classic format, units/long_name on
! every variable, one open/write_member/close call sequence).
!
! Dimensions: member (unlimited, one record per ensemble member),
! name_strlen (fixed-length character dimension for the two label
! variables). Point and saturation-config labels are stored as
! NF90_CHAR arrays (name_strlen, member) -- the portable, classic-
! format way to store per-record strings, rather than NF90_STRING
! (which needs NetCDF-4 and would break the NF90_CLOBBER/classic
! convention this project otherwise sticks to throughout).

module output_ensemble_netcdf

  use mpas_kind_types, only : RKIND
  use netcdf

  implicit none
  private

  public :: ensemble_summary_open, ensemble_summary_write_member, ensemble_summary_close

  integer, parameter :: NAME_STRLEN = 64

  integer :: ncid, dimid_member, dimid_strlen
  integer :: id_point_label, id_satconfig_label
  integer :: id_qc_peak, id_qr_peak, id_qi_peak, id_qs_peak, id_qg_peak
  integer :: id_time_to_precip, id_rainnc_final
  integer :: id_flashrate_peak, id_flashrate_cg_peak, id_flashrate_ic_peak
  integer :: id_flashcount_total, id_flashcount_cg_total, id_flashcount_ic_total
  integer :: id_snownc_final, id_snowdepth_final
  integer :: imember = 0

contains

  subroutine nc_check(status, where)
    integer, intent(in) :: status
    character(len=*), intent(in) :: where
    if (status /= nf90_noerr) then
      write(0,'(A)') 'output_ensemble_netcdf: NetCDF error at '//trim(where)//': '// &
           &          trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine nc_check

  !=================================================================
  subroutine ensemble_summary_open(path)
  !=================================================================
    character(len=*), intent(in) :: path

    call nc_check(nf90_create(trim(path), NF90_CLOBBER, ncid), 'create '//trim(path))

    call nc_check(nf90_def_dim(ncid, 'member', NF90_UNLIMITED, dimid_member), 'def dim member')
    call nc_check(nf90_def_dim(ncid, 'name_strlen', NAME_STRLEN, dimid_strlen), 'def dim name_strlen')

    call nc_check(nf90_def_var(ncid, 'point_label', NF90_CHAR, &
         &        (/dimid_strlen, dimid_member/), id_point_label), 'def point_label')
    call nc_check(nf90_put_att(ncid, id_point_label, 'long_name', &
         &        'ERA5 point/case label'), 'att point_label longname')

    call nc_check(nf90_def_var(ncid, 'satconfig_label', NF90_CHAR, &
         &        (/dimid_strlen, dimid_member/), id_satconfig_label), 'def satconfig_label')
    call nc_check(nf90_put_att(ncid, id_satconfig_label, 'long_name', &
         &        'saturation-perturbation configuration label'), 'att satconfig_label longname')

    call def_member_var('qc_peak', 'kg kg-1', &
         &  'peak cloud water mixing ratio reached anywhere in the column during the run', id_qc_peak)
    call def_member_var('qr_peak', 'kg kg-1', &
         &  'peak rain mixing ratio reached anywhere in the column during the run', id_qr_peak)
    call def_member_var('qi_peak', 'kg kg-1', &
         &  'peak cloud ice mixing ratio reached anywhere in the column during the run', id_qi_peak)
    call def_member_var('qs_peak', 'kg kg-1', &
         &  'peak snow mixing ratio reached anywhere in the column during the run', id_qs_peak)
    call def_member_var('qg_peak', 'kg kg-1', &
         &  'peak graupel mixing ratio reached anywhere in the column during the run', id_qg_peak)
    call def_member_var('time_to_precip_s', 's', &
         &  'time of first nonzero surface precipitation (-1 if none during the run)', id_time_to_precip)
    call def_member_var('rainnc_final_mm', 'mm', &
         &  'accumulated surface precipitation at the end of the run', id_rainnc_final)
    call def_member_var('flash_rate_peak', 'flashes min-1', &
         &  'peak total lightning flash rate during the run', id_flashrate_peak)
    call def_member_var('flash_rate_cg_peak', 'flashes min-1', &
         &  'peak cloud-to-ground flash rate during the run', id_flashrate_cg_peak)
    call def_member_var('flash_rate_ic_peak', 'flashes min-1', &
         &  'peak intra-cloud flash rate during the run', id_flashrate_ic_peak)
    call def_member_var('flash_count_total', 'flashes', &
         &  'time-integrated total flash count over the run', id_flashcount_total)
    call def_member_var('flash_count_cg_total', 'flashes', &
         &  'time-integrated cloud-to-ground flash count over the run', id_flashcount_cg_total)
    call def_member_var('flash_count_ic_total', 'flashes', &
         &  'time-integrated intra-cloud flash count over the run', id_flashcount_ic_total)
    call def_member_var('snownc_final_mm_lwe', 'mm (liquid-water equivalent)', &
         &  'accumulated surface snow+cloud-ice fallout at the end of the run, liquid-water equivalent', &
         &  id_snownc_final)
    call def_member_var('snow_depth_final_cm_10to1', 'cm', &
         &  'snow_depth_final_cm_10to1: snownc_final_mm_lwe converted to physical snow depth '// &
         &  'assuming a fixed 10:1 snow-to-liquid ratio -- see README.md caveat, NOT a direct '// &
         &  'model output', id_snowdepth_final)

    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'source', &
         &        'ensemble_driver: ERA5 x saturation-perturbation ensemble summary '// &
         &        '(replaces the earlier CSV format)'), 'att global source')

    call nc_check(nf90_enddef(ncid), 'enddef')
    imember = 0
  end subroutine ensemble_summary_open

  !-----------------------------------------------------------------
  subroutine def_member_var(name, units, long_name, varid)
  !-----------------------------------------------------------------
    character(len=*), intent(in) :: name, units, long_name
    integer, intent(out) :: varid
    call nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimid_member, varid), 'def '//trim(name))
    call nc_check(nf90_put_att(ncid, varid, 'units', units), 'att units '//trim(name))
    call nc_check(nf90_put_att(ncid, varid, 'long_name', long_name), 'att longname '//trim(name))
  end subroutine def_member_var

  !=================================================================
  subroutine ensemble_summary_write_member(point_label, satconfig_label, &
       &     qc_peak, qr_peak, qi_peak, qs_peak, qg_peak, &
       &     time_to_precip, rainnc_final, &
       &     flash_rate_peak, flash_rate_cg_peak, flash_rate_ic_peak, &
       &     flash_count_total, flash_count_cg_total, flash_count_ic_total, &
       &     snownc_final, snow_depth_final)
  !=================================================================
    character(len=*), intent(in) :: point_label, satconfig_label
    real(kind=RKIND), intent(in) :: qc_peak, qr_peak, qi_peak, qs_peak, qg_peak
    real(kind=RKIND), intent(in) :: time_to_precip, rainnc_final
    real(kind=RKIND), intent(in) :: flash_rate_peak, flash_rate_cg_peak, flash_rate_ic_peak
    real(kind=RKIND), intent(in) :: flash_count_total, flash_count_cg_total, flash_count_ic_total
    real(kind=RKIND), intent(in) :: snownc_final, snow_depth_final

    imember = imember + 1

    call nc_check(nf90_put_var(ncid, id_point_label, trim(point_label), &
         &        start=(/1,imember/), count=(/min(len_trim(point_label),NAME_STRLEN),1/)), 'put point_label')
    call nc_check(nf90_put_var(ncid, id_satconfig_label, trim(satconfig_label), &
         &        start=(/1,imember/), count=(/min(len_trim(satconfig_label),NAME_STRLEN),1/)), &
         &        'put satconfig_label')

    call nc_check(nf90_put_var(ncid, id_qc_peak, qc_peak, start=(/imember/)), 'put qc_peak')
    call nc_check(nf90_put_var(ncid, id_qr_peak, qr_peak, start=(/imember/)), 'put qr_peak')
    call nc_check(nf90_put_var(ncid, id_qi_peak, qi_peak, start=(/imember/)), 'put qi_peak')
    call nc_check(nf90_put_var(ncid, id_qs_peak, qs_peak, start=(/imember/)), 'put qs_peak')
    call nc_check(nf90_put_var(ncid, id_qg_peak, qg_peak, start=(/imember/)), 'put qg_peak')
    call nc_check(nf90_put_var(ncid, id_time_to_precip, time_to_precip, start=(/imember/)), 'put time_to_precip')
    call nc_check(nf90_put_var(ncid, id_rainnc_final, rainnc_final, start=(/imember/)), 'put rainnc_final')
    call nc_check(nf90_put_var(ncid, id_flashrate_peak, flash_rate_peak, start=(/imember/)), 'put flash_rate_peak')
    call nc_check(nf90_put_var(ncid, id_flashrate_cg_peak, flash_rate_cg_peak, start=(/imember/)), &
         &        'put flash_rate_cg_peak')
    call nc_check(nf90_put_var(ncid, id_flashrate_ic_peak, flash_rate_ic_peak, start=(/imember/)), &
         &        'put flash_rate_ic_peak')
    call nc_check(nf90_put_var(ncid, id_flashcount_total, flash_count_total, start=(/imember/)), &
         &        'put flash_count_total')
    call nc_check(nf90_put_var(ncid, id_flashcount_cg_total, flash_count_cg_total, start=(/imember/)), &
         &        'put flash_count_cg_total')
    call nc_check(nf90_put_var(ncid, id_flashcount_ic_total, flash_count_ic_total, start=(/imember/)), &
         &        'put flash_count_ic_total')
    call nc_check(nf90_put_var(ncid, id_snownc_final, snownc_final, start=(/imember/)), 'put snownc_final')
    call nc_check(nf90_put_var(ncid, id_snowdepth_final, snow_depth_final, start=(/imember/)), &
         &        'put snow_depth_final')

  end subroutine ensemble_summary_write_member

  !=================================================================
  subroutine ensemble_summary_close()
  !=================================================================
    call nc_check(nf90_close(ncid), 'close ensemble summary')
  end subroutine ensemble_summary_close

end module output_ensemble_netcdf
