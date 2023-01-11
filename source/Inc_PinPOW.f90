
      MODULE Inc_PinPOW

      IMPLICIT NONE

      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: HFF
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: PinPOW_2D
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: PinQ_2D
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: PinBU_2D
      REAL(8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: PinPOW_3D
      REAL(8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: PinQ_3D
      REAL(8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: PinBU_3D
      REAL(8), DIMENSION(:), ALLOCATABLE :: PinMTU_2D
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: PinMTU_3D
      real(8), dimension(:,:,:,:,:), allocatable :: pinpow_3d_old
      real(8), dimension(:,:,:,:,:), allocatable :: pinq_3d_old

      REAL(8) :: Fq, Fxy, Fdh, Fr

      logical(1) :: flag_savepbu=.false.

      LOGICAL(1) :: Flag_PinPOW, Flag_NPinOdd

      real(8) :: Fq_val
      real(8) :: Fr_val
      real(8) :: Fz_val
      real(8) :: Fz_ave
      real(8) :: FdH_val
      real(8) :: max_Fxy_val
      integer :: Fq_loc(5)
      integer :: Fr_loc(4)
      integer :: Fz_loc(1)
      integer :: FdH_loc(4)
      integer :: max_Fxy_loc(3)

      real(8) :: max_PinBU_2D
      real(8) :: max_PinPOW_2D
      real(8) :: max_PinPOW_3D
      real(8) :: max_PinQ_2D
      real(8) :: max_PinQ_3D
      integer :: max_PinBU_2D_loc(4)
      integer :: max_PinPOW_2D_loc(4)
      integer :: max_PinPOW_3D_loc(5)
      integer :: max_PinQ_2D_loc(4)
      integer :: max_PinQ_3D_loc(5)

      real(8) :: maxPBU_2D
      real(8) :: maxPPD_2D
      real(8) :: maxPPD_3D
      real(8) :: maxPPW_2D
      real(8) :: maxPPW_3D
      integer :: PBU_2D_loc(4)
      integer :: PPD_2D_loc(4)
      integer :: PPD_3D_loc(5)
      integer :: PPW_2D_loc(4)
      integer :: PPW_3D_loc(5)

      ! file name = AnmpM - anm parameters
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: akratio   !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: akappa    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: amu       !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: ar        !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: as        !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: ardet     !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: asnmu     !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: acnmu     !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: asnka     !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: acnka     !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: asnmux    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: acnmux    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: asnkax    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: acnkax    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: asnmuy    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: acnmuy    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: asnkay    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: acnkay    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: asnmuz    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: acnmuz    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: asnkaz    !(nxy,nz)
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: acnkaz    !(nxy,nz)
      LOGICAL(1), DIMENSION(:,:), ALLOCATABLE :: kflag     !(nxy,nz)
      REAL(8):: epsanm = 0.005D0

      logical(1) :: flag_first_pinpow=.false.

      END MODULE Inc_PinPOW
