! test_uedmf_entdet_map.F90 -- evaluates the ACTUAL Fortran closure (uedmf_entdet, inside MYNN) on a grid of
! plume buoyancy b and velocity w, writes a table. Used to characterize how eps_dyn depends on (b,w).
program test_uedmf_entdet_map
  use ccpp_kind_types, only : kind_phys
  use bl_mynn, only : bl_mynn_init
  use bl_mynn_subroutines, only : uedmf_entdet
  implicit none
  real(kind_phys), parameter :: bs(8) = (/0.005,0.01,0.02,0.04,0.06,0.08,0.10,0.15/)
  real(kind_phys), parameter :: ws(8) = (/1.,2.,5.,10.,15.,20.,25.,30./)
  real(kind_phys) :: eps,det,epst
  integer :: i,j,u
  open(newunit=u,file='entdet_map.txt',status='replace')
  do i=1,8
    do j=1,8
      ! dry updraft (RH_i<1): a_i=1e-4, TKE e0=0.05 m2/s2, plume height 8 km
      call uedmf_entdet(bs(i),0.,ws(j),0.,1.0e-4_kind_phys,0.05_kind_phys,0.6_kind_phys,0.5_kind_phys,8000._kind_phys,eps,det,epst)
      write(u,'(2F8.3,3ES14.5)') bs(i),ws(j),eps,det,epst
    enddo
  enddo
  close(u)
end program
