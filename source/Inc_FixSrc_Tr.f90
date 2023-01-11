
      MODULE Inc_FixSrc_Tr

      IMPLICIT NONE

      real(8) :: gamma
      real(8) :: rgammap1
      real(8) :: gammam1
      real(8) :: rgamma
      real(8) :: cetakrl
      real(8) :: cappab
      real(8) :: rcdeltl

      REAL(8) :: omegalm
      REAL(8) :: omegal0
      REAL(8) :: omegalp
      REAL(8) :: omegalpd

      REAL(8) :: sd
      REAL(8) :: sdn
      REAL(8) :: sdnt
      REAL(8) :: spnt

      REAL(8), DIMENSION(:), ALLOCATABLE :: rvdtvol
      REAL(8), DIMENSION(:), ALLOCATABLE :: rvdtdvol

      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: rldt
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: rldtgp1
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: cappap1
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: capbrldt
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: capbrldt2

      END MODULE Inc_FixSrc_Tr
