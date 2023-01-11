
      MODULE Inc_3D

      IMPLICIT NONE

      real(8) :: keff
      real(8) :: keff_Old
      real(8) :: keff_OOld

      real(8), allocatable :: FisSrc          (:,:)
      real(8), allocatable :: FisSrc_Old      (:,:)
      real(8), allocatable :: FisSrc_OOld     (:,:)
      real(8), allocatable :: FisSrc_Iout     (:,:)
      real(8), allocatable :: Power           (:,:)
      real(8), allocatable :: Normal_Power    (:,:)
      real(8), allocatable :: BU              (:,:)
      real(8), allocatable :: BU_Old          (:,:)
      real(8), allocatable :: BU_Predictor    (:,:)
      real(8), allocatable :: T_Fuel          (:,:)
      real(8), allocatable :: T_Mod           (:,:)
      real(8), allocatable :: D_Mod           (:,:)
      real(8), allocatable :: T_Fuel_Old_XSFB (:,:)
      real(8), allocatable :: T_Mod_Old_XSFB  (:,:)
      real(8), allocatable :: D_Mod_Old_XSFB  (:,:)
      real(8), allocatable :: MTU             (:,:)

      real(8), allocatable :: Flux          (:,:,:)
      real(8), allocatable :: Flux_Old      (:,:,:)
      real(8), allocatable :: Flux_OOld     (:,:,:)
      real(8), allocatable :: Flux_Old_XSFB (:,:,:)


      real(8) :: Avg_Power
      real(8), allocatable :: Avg_Flux(:)

      real(8), allocatable :: Flux_det_old (:,:,:)
      real(8), allocatable :: Flux_det     (:,:,:)
      integer :: OPT_Det
      integer :: print_numden=0
      integer :: NN_TIME
      integer :: NN_TIME2
      integer :: OPT_DET_Power
      integer :: I_ITER_JR = 0

      real(8), allocatable :: flux_adj (:,:,:)

      real(8), allocatable :: t_fuel_bk (:,:)
      real(8), allocatable :: t_fuel_ss (:,:)
      real(8), allocatable :: t_mod_bk  (:,:)
      real(8), allocatable :: t_mod_ss  (:,:)

      END MODULE Inc_3D
