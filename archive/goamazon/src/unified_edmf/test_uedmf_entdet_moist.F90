! test_uedmf_entdet_moist.F90 -- entrainment/detrainment of a SATURATED (cloudy) plume vs environment humidity,
! and its dependence on beta (Eq. 38). Uses the actual Fortran closure inside MYNN.
program test_uedmf_entdet_moist
  use ccpp_kind_types, only : kind_phys
  use bl_mynn, only : bl_mynn_init
  use bl_mynn_subroutines, only : uedmf_entdet, uedmf_beta
  implicit none
  real(kind_phys) :: eps,det,epst,rh(5)
  integer :: i,ib
  real(kind_phys) :: betas(3)
  rh = (/0.3,0.5,0.7,0.9,0.99/)
  betas = (/1.0,2.0,4.0/)
  print '(A)', 'pluma nublada (RH_i=1), b=0.06 m/s2, w=15 m/s, a=1e-4, e0=0.05, H=8 km'
  print '(A)', ' beta  RH_env   eps_dyn(1e-5/m)  delta(1e-5/m)  delta/eps'
  do ib=1,3
    uedmf_beta = betas(ib)
    do i=1,5
      call uedmf_entdet(0.06_kind_phys,0._kind_phys,15._kind_phys,0._kind_phys,1.0e-4_kind_phys,0.05_kind_phys,1.0_kind_phys,rh(i),8000._kind_phys,eps,det,epst)
      print '(F5.1,F8.2,2F16.2,F12.2)', betas(ib), rh(i), 1.0e5*eps, 1.0e5*det, det/eps
    enddo
  enddo
end program
