! mpas_stubs_offline.F90
!
! Minimal, serial (no-MPI) replacements for the MPAS infrastructure
! modules that module_mp_thompson.F, mp_radar.F90 (the real one, from
! the external NCAR/MMM-physics repository -- see README.md) and
! mpas_atmphys_functions.F "use". Confirmed by grepping the extracted
! sources for every "use mpas_*"/"USE mpas_*" line; this is the
! complete list. Compiling with -Dmpas selects the branch of
! module_mp_thompson.F's own aerosol-activation-table subroutine
! table_ccnAct (needing module_domain/module_dm) out entirely -- see
! "#if !defined(mpas)" there.
!
! Same pattern as rrtmg_offline/src/mpas_stubs_offline.F90 from the
! companion radiation project: this replaces (a) MPI broadcasts with
! no-ops, since an offline single-process driver has nothing to
! broadcast to, and (b) MPAS's logging/unit-management/error-handling
! helpers with trivial plain-Fortran equivalents. None of this touches
! microphysics physics.

module mpas_kind_types
  implicit none
  integer, parameter :: R4KIND  = selected_real_kind(6)
  integer, parameter :: R8KIND  = selected_real_kind(12)
  integer, parameter :: RKIND   = selected_real_kind(12)
  integer, parameter :: I8KIND  = selected_int_kind(18)
  integer, parameter :: StrKIND = 512
end module mpas_kind_types


module mpas_derived_types
  implicit none
  integer, parameter :: MPAS_LOG_ERR = 3
  type :: dm_info
    integer :: my_proc_id = 0
    integer :: nprocs     = 1
  end type dm_info
end module mpas_derived_types


module mpas_dmpar
  use mpas_derived_types, only : dm_info
  implicit none
  integer, parameter :: IO_NODE = 0
contains
  subroutine mpas_dmpar_bcast_int(dminfo, i)
    type(dm_info), intent(in) :: dminfo
    integer, intent(inout) :: i
  end subroutine mpas_dmpar_bcast_int

  subroutine mpas_dmpar_bcast_real(dminfo, r)
    use mpas_kind_types, only : RKIND
    type(dm_info), intent(in) :: dminfo
    real(kind=RKIND), intent(inout) :: r
  end subroutine mpas_dmpar_bcast_real

  subroutine mpas_dmpar_bcast_reals(dminfo, n, r)
    use mpas_kind_types, only : RKIND
    type(dm_info), intent(in) :: dminfo
    integer, intent(in) :: n
    real(kind=RKIND), dimension(*), intent(inout) :: r
  end subroutine mpas_dmpar_bcast_reals
end module mpas_dmpar


module mpas_atmphys_utilities
  implicit none
contains
  ! Real MPAS's physics_message just forwards to mpas_log_write;
  ! reproduced narrowly for that one purpose.
  subroutine physics_message(str)
    character(len=*), intent(in) :: str
    write(*,'(A)') trim(str)
  end subroutine physics_message

  subroutine physics_error_fatal(mess)
    character(len=*), intent(in) :: mess
    write(0,'(A)') 'FATAL: '//trim(mess)
    stop 1
  end subroutine physics_error_fatal
end module mpas_atmphys_utilities


module mpas_log
  implicit none
contains
  subroutine mpas_log_write(mess, intArgs, realArgs, charArgs, messageType)
    character(len=*), intent(in) :: mess
    integer, dimension(:), intent(in), optional :: intArgs
    real,    dimension(:), intent(in), optional :: realArgs
    character(len=*), dimension(:), intent(in), optional :: charArgs
    integer, intent(in), optional :: messageType

    if (present(intArgs) .and. len_trim(mess) > 0) then
      write(*,'(A,1x,I0)') trim(mess), intArgs(1)
    else if (len_trim(mess) > 0) then
      write(*,'(A)') trim(mess)
    end if
  end subroutine mpas_log_write
end module mpas_log


module mpas_io_units
  implicit none
contains
  subroutine mpas_new_unit(unit_num, unformatted)
    integer, intent(out) :: unit_num
    logical, intent(in), optional :: unformatted
    integer :: i
    logical :: is_open
    do i = 20, 999
      inquire(unit=i, opened=is_open)
      if (.not. is_open) then
        unit_num = i
        return
      end if
    end do
    unit_num = -1
  end subroutine mpas_new_unit

  subroutine mpas_release_unit(unit_num)
    integer, intent(in) :: unit_num
  end subroutine mpas_release_unit
end module mpas_io_units
