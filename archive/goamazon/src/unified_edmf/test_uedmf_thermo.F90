! test_uedmf_thermo.F90 -- consistency tests for the plume thermodynamics added inside MYNN
! (uedmf_condensation / uedmf_thil_from_T in bl_mynn_subroutines.F90). No closed-form solution exists,
! so this checks that the converged state satisfies the equations it claims to solve.
program test_uedmf_thermo
  use ccpp_kind_types, only : kind_phys
  use bl_mynn_common
  use bl_mynn, only : bl_mynn_init
  use bl_mynn_subroutines, only : uedmf_condensation, uedmf_thil_from_T
  use mynn_shared, only : qsat_blend
  implicit none
  real(kind_phys) :: THV1,THV2,THV0,Tm1,Tm2,QC1,QC2,QC0,exn,th,P,QT,THL,leff,fi,tsat
  integer :: i, nfail
  character(len=256) :: errmsg; integer :: errflg
  nfail = 0
  call bl_mynn_init(con_cp=1004.5_kind_phys, con_cpv=1885._kind_phys, con_cice=2106._kind_phys, con_cliq=4190._kind_phys, &
       con_ep1=0.608_kind_phys, con_ep2=0.622_kind_phys, con_grav=9.81_kind_phys, con_karman=0.4_kind_phys, &
       con_p0=100000._kind_phys, con_rd=287._kind_phys, con_rv=461.5_kind_phys, con_svp1=0.6112_kind_phys, &
       con_svp2=17.67_kind_phys, con_svp3=29.65_kind_phys, con_svpt0=273.15_kind_phys, con_xlf=3.34e5_kind_phys, &
       con_xls=2.834e6_kind_phys, con_xlv=2.5e6_kind_phys, errmsg=errmsg, errflg=errflg)

  ! levels: (P, THL, QT) chosen supersaturated at several temperatures (warm, mixed, cold)
  call check(90000._kind_phys, 300._kind_phys, 0.0180_kind_phys, 'quente (~290 K)')
  call check(60000._kind_phys, 305._kind_phys, 0.0100_kind_phys, 'misto  (~265 K)')
  call check(30000._kind_phys, 320._kind_phys, 0.0040_kind_phys, 'frio   (~240 K)')
  call check(20000._kind_phys, 335._kind_phys, 0.0010_kind_phys, 'muito frio (~220 K)')
  if (nfail == 0) then
     print *, 'TODOS OS CHECKS PASSARAM'
  else
     print *, 'FALHAS:', nfail
  endif
contains
  subroutine check(P,THL,QT,lab)
    real(kind_phys), intent(in) :: P,THL,QT
    character(len=*), intent(in) :: lab
    real(kind_phys) :: r1,r2,r3,exn
    exn = (P/p1000mb)**rcp
    QC1=0.; QC2=0.; QC0=0.
    call uedmf_condensation(QT,THL,P,2000._kind_phys,THV1,QC1,Tm1,1)
    call uedmf_condensation(QT,THL,P,2000._kind_phys,THV2,QC2,Tm2,2)
    ! (a) defining equation of the ice-liquid theta_l at the returned (T,QC):  THL == thil_from_T
    r1 = abs(uedmf_thil_from_T(Tm1,QC1,exn,1) - THL)
    r2 = abs(uedmf_thil_from_T(Tm2,QC2,exn,2) - THL)
    ! (b) saturation: QT - QC == qsat(T) when condensate exists
    r3 = 0.
    if (QC2 > 0.) r3 = abs((QT-QC2) - qsat_blend(Tm2,t0c,tice,P))
    print '(A,A,4X,A,F8.3,A,F8.3,A,ES9.2,A,ES9.2,A,ES9.2)', trim(lab), ':', 'T1=',Tm1,' T2=',Tm2,'  res_thil(1)=',r1,' res_thil(2)=',r2,' res_sat=',r3
    print '(A,2F9.5,A,2F9.3)', '     QC modo1/modo2 (g/kg):', QC1*1000., QC2*1000., '   THV modo1/modo2:', THV1, THV2
    if (r1 > 1.e-4 .or. r2 > 1.e-4 .or. r3 > 5.e-5) then
       print *, '   FALHA: residuo de autoconsistencia grande'; nfail = nfail + 1
    endif
    ! (c) fusion must never cool the parcel: T2 >= T1 (extra latent heat), and be equal when fi=0 (warm)
    if (Tm2 < Tm1 - 1.e-9) then
       print *, '   FALHA: fusao esfriou a parcela'; nfail = nfail + 1
    endif
    if (Tm1 > t0c+0.5 .and. abs(Tm2-Tm1) > 1.e-9) then
       print *, '   FALHA: acima de 0 C os modos 1 e 2 deveriam coincidir'; nfail = nfail + 1
    endif
  end subroutine check
end program test_uedmf_thermo
