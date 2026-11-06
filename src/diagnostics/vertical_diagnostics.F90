! Read-only snapshots at the MYNN call: no state or tendency is returned.
module vertical_diagnostics
  use mpas_kind_types, only: RKIND
  use netcdf
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: vd_open, vd_write, vd_close
  integer :: nc, td, ld, zd, rec=0, ti, zi, li, dzi
  integer :: ids(17)
contains
  subroutine check(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      write(*,*) 'vertical diagnostics: ',trim(nf90_strerror(status))
      stop 1
    endif
  end subroutine

  subroutine define(name,units,dims,id)
    character(*), intent(in) :: name,units
    integer, intent(in) :: dims(:)
    integer, intent(out) :: id
    call check(nf90_def_var(nc,name,nf90_double,dims,id,deflate_level=1))
    call check(nf90_put_att(nc,id,'units',units))
  end subroutine

  subroutine vd_open(path,nz,z,dz,edmf_type,mp_mode)
    character(*), intent(in) :: path
    integer, intent(in) :: nz,edmf_type,mp_mode
    real(RKIND), intent(in) :: z(nz),dz(nz)
    real(RKIND) :: interfaces(nz)
    integer :: k
    call check(nf90_create(path,nf90_netcdf4,nc))
    call check(nf90_def_dim(nc,'time',nf90_unlimited,td))
    call check(nf90_def_dim(nc,'lev',nz,ld))
    ! MYNN edmf_* k=1..nz-1 is at the top of layer k. The surface
    ! and model lid have no diagnosed plume-mean velocity; omit both.
    call check(nf90_def_dim(nc,'interface',nz-1,zd))
    call define('time','s',[td],ti)
    call define('lev','m',[ld],li)
    call define('layer_thickness','m',[ld],dzi)
    call define('z_interface','m',[zd],zi)
    call define('w_LS','m s-1',[ld,td],ids(1))
    call define('w_LS_interface','m s-1',[zd,td],ids(2))
    call define('rho_MYNN','kg m-3',[ld,td],ids(3))
    call define('rhoz','kg m-3',[zd,td],ids(4))
    call define('mf_w','kg m-2 s-1',[zd,td],ids(5))
    call define('mf_a','1',[zd,td],ids(6))
    call define('w_EDMF','m s-1',[zd,td],ids(7))
    call define('plume_active','1',[zd,td],ids(8))
    call define('temperature','K',[ld,td],ids(9))
    call define('pressure','Pa',[ld,td],ids(10))
    call define('qv','kg kg-1',[ld,td],ids(11))
    call define('qc','kg kg-1',[ld,td],ids(12))
    call define('qr','kg kg-1',[ld,td],ids(13))
    call define('qi','kg kg-1',[ld,td],ids(14))
    call define('qs','kg kg-1',[ld,td],ids(15))
    call define('qg','kg kg-1',[ld,td],ids(16))
    call define('graupel_plume_flux_available','1',[td],ids(17))
    call check(nf90_put_att(nc,nf90_global,'sample_phase', &
      'At forcing time, immediately after MYNN call, before applying tendencies, transport or microphysics'))
    call check(nf90_put_att(nc,nf90_global,'velocity_frame', &
      'MYNN receives w2=0; w_EDMF is its area-weighted plume velocity. w_LS prescribed separately; never summed.'))
    call check(nf90_put_att(nc,nf90_global,'rhoz_definition', &
      '(rho(k)*dz(k+1)+rho(k+1)*dz(k))/(dz(k+1)+dz(k)); exact MYNN interpolation'))
    call check(nf90_put_att(nc,nf90_global,'inactive_plume_policy', &
      'mf_a<=0: plume_active=0, w_EDMF=NaN; raw mf_a and mf_w retained'))
    call check(nf90_put_att(nc,nf90_global,'mf_a_definition', &
      'Native area fraction after Psig_w; mf_w divided by unscaled sum(UPA) inside MYNN, so divide mf_w by rhoz only'))
    call check(nf90_put_att(nc,nf90_global,'hydrometeor_support', &
      'Grid-column dry-air mixing ratios at MYNN entry, NOT plume mixing ratios or covariant plume fluxes'))
    call check(nf90_put_att(nc,nf90_global,'species_flux_status', &
      'uedmf_sp_flux channels qr,qi,qs,qg require unified edmf_type=1 and explicit plume MP >=2; unavailable in main type=0 mode=0. Other modes not certified by this exporter.'))
    call check(nf90_put_att(nc,nf90_global,'mynn_edmf_type',edmf_type))
    call check(nf90_put_att(nc,nf90_global,'mynn_uedmf_mp',mp_mode))
    call check(nf90_enddef(nc))
    interfaces(1)=dz(1)
    do k=2,nz
      interfaces(k)=interfaces(k-1)+dz(k)
    enddo
    call check(nf90_put_var(nc,li,z))
    call check(nf90_put_var(nc,dzi,dz))
    call check(nf90_put_var(nc,zi,interfaces(:nz-1)))
    rec=0
  end subroutine

  subroutine vd_write(time,nz,dz,w,rho,area,mfw,t,p,qv,qc,qr,qi,qs,qg)
    integer, intent(in) :: nz
    real(RKIND), intent(in) :: time,dz(nz),w(nz),rho(nz),area(nz),mfw(nz)
    real(RKIND), intent(in) :: t(nz),p(nz),qv(nz),qc(nz),qr(nz),qi(nz),qs(nz),qg(nz)
    real(RKIND) :: rz(nz-1),wi(nz-1),we(nz-1),active(nz-1),zero
    integer :: k
    do k=1,nz-1
      rz(k)=(rho(k)*dz(k+1)+rho(k+1)*dz(k))/(dz(k+1)+dz(k))
      wi(k)=(w(k)*dz(k+1)+w(k+1)*dz(k))/(dz(k+1)+dz(k))
      we(k)=ieee_value(0.0_RKIND,ieee_quiet_nan)
      active(k)=0.0_RKIND
      if (area(k)>0.0_RKIND.and.rz(k)>0.0_RKIND) then
        we(k)=mfw(k)/rz(k)
        active(k)=1.0_RKIND
      endif
    enddo
    rec=rec+1
    call check(nf90_put_var(nc,ti,[time],start=[rec],count=[1]))
    call put(ids(1),w);call put(ids(2),wi);call put(ids(3),rho);call put(ids(4),rz)
    call put(ids(5),mfw(:nz-1));call put(ids(6),area(:nz-1))
    call put(ids(7),we);call put(ids(8),active)
    call put(ids(9),t);call put(ids(10),p);call put(ids(11),qv)
    call put(ids(12),qc);call put(ids(13),qr);call put(ids(14),qi)
    call put(ids(15),qs);call put(ids(16),qg)
    zero=0.0_RKIND
    call check(nf90_put_var(nc,ids(17),[zero],start=[rec],count=[1]))
  end subroutine

  subroutine put(id,x)
    integer,intent(in) :: id
    real(RKIND),intent(in) :: x(:)
    call check(nf90_put_var(nc,id,x,start=[1,rec],count=[size(x),1]))
  end subroutine

  subroutine vd_close()
    call check(nf90_close(nc))
  end subroutine
end module
