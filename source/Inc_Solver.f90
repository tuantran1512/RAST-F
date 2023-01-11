
      module Inc_Solver
      implicit none

      real(8), allocatable :: rvdelt     (:,:,:)
      real(8), allocatable :: rvdelt_old (:,:,:)
      real(8), allocatable :: betap      (:,:)
      real(8), allocatable :: betap_old  (:,:)

      integer:: nkincomp=1
      integer, allocatable :: kincomp    (:)
      real(8), allocatable :: lambda_d_I (:,:)
      real(8), allocatable :: tvelo      (:,:)
      real(8), allocatable :: tbeta      (:,:)

      real(8), allocatable :: cappa       (:,:,:)
      real(8), allocatable :: omegam      (:,:,:)
      real(8), allocatable :: omega0      (:,:,:)
      real(8), allocatable :: omegap      (:,:,:)
      real(8), allocatable :: u_omega     (:,:,:)
      real(8), allocatable :: flux_extra  (:,:,:)
      real(8), allocatable :: src_tr_fac1 (:,:,:)
      real(8), allocatable :: src_tr_fac2 (:,:,:)

      real(8), allocatable :: iPgenT  (:)
      real(8), allocatable :: rhocoef (:)
      real(8), allocatable :: Pbeta   (:)
      real(8), allocatable :: Plambda (:,:)
      real(8), allocatable :: Pzeta   (:,:)

      real(8):: cetak
      real(8):: cetaf
      real(8):: cetac
      real(8):: nordprec
      real(8):: cetak0
      real(8):: cetakb
      real(8):: cetakr
      real(8):: cetafb
      real(8):: cetacb
      real(8):: cetadelt

      end module Inc_Solver
