! Read-only MYNN audit. Disabled unless explicitly opened by the driver.
module mynn_audit
  use mpas_kind_types, only: RKIND
  use netcdf
  implicit none
  logical :: ma_enabled=.false., ma_main=.true.
  logical :: ma_skip_shadow=.false.
  real(RKIND) :: ma_time=0., ma_dt=0., ma_gate(16)=0., ma_term(10,8)=0.
  real(RKIND) :: ma_boundary(20)=0.
  integer :: boundary_dim,boundary_id
  integer :: nc,rec=0,ld,td,gd,pd,cd,ti,zi,si,gi,pi,ids(20),aid(6)
  real(RKIND),allocatable :: qc_entry(:),qi_entry(:),qc_local(:),qi_local(:),qc_full(:),qi_full(:)
contains
  subroutine ck(status)
    integer,intent(in)::status
    if(status/=nf90_noerr) then
      print *, 'MYNN audit: ',trim(nf90_strerror(status));stop 1
    endif
  end subroutine
  subroutine define(name,unit,dims,id)
    character(*),intent(in)::name,unit
    integer,intent(in)::dims(:)
    integer,intent(out)::id
    call ck(nf90_def_var(nc,name,nf90_double,dims,id,deflate_level=1))
    call ck(nf90_put_att(nc,id,'units',unit))
  end subroutine
  subroutine ma_open(path,nz,z)
    character(*),intent(in)::path
    integer,intent(in)::nz
    real(RKIND),intent(in)::z(nz)
    character(32)::names(20),units(20)
    integer::j
    names=[character(32)::'tke','qc_tendency_local','qc_tendency_full','qi_tendency_local','qi_tendency_full', &
      'temperature','pressure','rho','qv_entry','qc_entry','qi_entry','theta_tendency','qv_tendency', &
      'mixing_length','diffusivity_heat','diffusivity_momentum','edmf_area','edmf_rho_w','cloud_fraction','edmf_entrainment']
    units=[character(32)::'m2 s-2','kg kg-1 s-1','kg kg-1 s-1','kg kg-1 s-1','kg kg-1 s-1', &
      'K','Pa','kg m-3','kg kg-1','kg kg-1','kg kg-1','K s-1','kg kg-1 s-1','m','m2 s-1','m2 s-1','1','kg m-2 s-1','1','m-1']
    allocate(qc_entry(nz),qi_entry(nz),qc_local(nz),qi_local(nz),qc_full(nz),qi_full(nz))
    call ck(nf90_create(path,nf90_netcdf4,nc))
    call ck(nf90_def_dim(nc,'lev',nz,ld));call ck(nf90_def_dim(nc,'time',nf90_unlimited,td))
    call ck(nf90_def_dim(nc,'surface_field',20,si));call ck(nf90_def_dim(nc,'gate_field',16,gd))
    call ck(nf90_def_dim(nc,'terminal_field',10,cd));call ck(nf90_def_dim(nc,'plume',8,pd))
    call ck(nf90_def_dim(nc,'boundary_field',20,boundary_dim))
    call define('surface_boundary','see boundary_fields',[boundary_dim,td],boundary_id)
    call ck(nf90_put_att(nc,nf90_global,'boundary_fields', &
      'incoming_hfx_W_m2 incoming_qfx_kg_m2_s bulk_hfx_W_m2 bulk_qfx_kg_m2_s applied_hfx_W_m2 applied_qfx_kg_m2_s applied_lh_W_m2 forcing_hfx_W_m2 forcing_lh_W_m2 wstar_m_s ust_m_s rmol_bulk_m-1 rho_surface_kg_m3 previous_pblh_m tsk_K qsfc_specific mavail mode memory_enabled dt_s'))
    call define('time','s',[td],ti);call define('lev','m AGL',[ld],zi)
    call define('surface','see surface_fields',[si,td],j);si=j
    call define('gate','see gate_fields',[gd,td],gi)
    call define('termination','see terminal_fields',[cd,pd,td],pi)
    do j=1,20
      call define(trim(names(j)),trim(units(j)),[ld,td],ids(j))
    enddo
    call define('qc_applied_increment','kg kg-1',[ld,td],aid(1))
    call define('qi_applied_increment','kg kg-1',[ld,td],aid(2))
    call define('qc_clipping_increment','kg kg-1',[ld,td],aid(3))
    call define('qi_clipping_increment','kg kg-1',[ld,td],aid(4))
    call define('qc_expected_increment','kg kg-1',[ld,td],aid(5))
    call define('qi_expected_increment','kg kg-1',[ld,td],aid(6))
    call ck(nf90_put_att(nc,nf90_global,'surface_fields', &
      'tsk_K qsfc_specific hfx_W_m2 lh_W_m2 qfx_kg_m2_s ust_m_s rmol_sf_m-1 pblh_m z0_m mavail dx_m previous_pblh_m maxwidth_m maxmf_kg_m2_s plume_top_m wstar_m_s psfc_Pa surface_theta_K lowest_rv_kg_kg dt_s'))
    call ck(nf90_put_att(nc,nf90_global,'gate_fields', &
      'fltv_effective_K_m_s superadiabatic maxwidth_before_gate_m minwidth_m dx_cap_m pbl_cap_m cloud_cap_m flux_cap_m Psig_w Psig_shcu max_resolved_w_asymmetric_m_s nup_after_gate flux_limiter ktop entry_fltv_K_m_s wind_area_factor'))
    call ck(nf90_put_att(nc,nf90_global,'terminal_fields', &
      'reason attempted_interface_height_m last_positive_w_m_s buoyancy_m_s2 entrainment_m-1 buoyancy_velocity_increment_m_s drag_velocity_increment_m_s unclipped_next_w_m_s diameter_m attempted_k'))
    call ck(nf90_put_att(nc,nf90_global,'termination_codes', &
      '0=not_launched;1=failed_first_level;2=momentum_stop_negative_buoyancy_with_drag;3=entrainment_drag_stop_despite_positive_buoyancy;4=indeterminate_stop;5=reached_native_loop_top;6=alternative_core_not_instrumented'))
    call ck(nf90_put_att(nc,nf90_global,'delta_scope', &
      'full minus EDMF-OFF shadow: response of whole MYNN call, not an isolated covariant plume species flux'))
    call ck(nf90_put_att(nc,nf90_global,'local_tendencies_available',merge(0,1,ma_skip_shadow)))
    call ck(nf90_enddef(nc));call ck(nf90_put_var(nc,zi,z))
    ma_enabled=.true.;rec=0
  end subroutine
  subroutine ma_write(surface,profiles)
    real(RKIND),intent(in)::surface(20),profiles(:,:)
    integer::j
    rec=rec+1
    call ck(nf90_put_var(nc,ti,[ma_time],start=[rec],count=[1]))
    call ck(nf90_put_var(nc,si,surface,start=[1,rec],count=[20,1]))
    call ck(nf90_put_var(nc,boundary_id,ma_boundary,start=[1,rec],count=[20,1]))
    call ck(nf90_put_var(nc,gi,ma_gate,start=[1,rec],count=[16,1]))
    call ck(nf90_put_var(nc,pi,ma_term,start=[1,1,rec],count=[10,8,1]))
    do j=1,20
      call put(ids(j),profiles(:,j))
    enddo
    qc_entry=profiles(:,10);qi_entry=profiles(:,11)
    qc_local=profiles(:,2);qc_full=profiles(:,3);qi_local=profiles(:,4);qi_full=profiles(:,5)
  end subroutine
  subroutine ma_applied(qc,qi,mode)
    real(RKIND),intent(in)::qc(:),qi(:)
    integer,intent(in)::mode
    real(RKIND)::dc(size(qc)),di(size(qi))
    dc=0.;di=0.
    if(mode==1) then
      dc=ma_dt*qc_local;di=ma_dt*qi_local
    else if(mode==2) then
      dc=ma_dt*qc_full;di=ma_dt*qi_full
    endif
    call put(aid(1),qc-qc_entry);call put(aid(2),qi-qi_entry)
    call put(aid(3),qc-qc_entry-dc);call put(aid(4),qi-qi_entry-di)
    call put(aid(5),dc);call put(aid(6),di)
  end subroutine
  subroutine put(id,values)
    integer,intent(in)::id
    real(RKIND),intent(in)::values(:)
    call ck(nf90_put_var(nc,id,values,start=[1,rec],count=[size(values),1]))
  end subroutine
  subroutine ma_close()
    call ck(nf90_close(nc));ma_enabled=.false.
  end subroutine
end module
