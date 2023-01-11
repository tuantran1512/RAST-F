
      MODULE Inc_RST

      IMPLICIT NONE

      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: BU_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Normal_Power_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: T_Fuel_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: T_Mod_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: D_Mod_RST

      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Flux_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: maXS_f_3D_RST

      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_U34_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_U35_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_U36_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_U37_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_U38_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Np37_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Np38_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Np39_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Pu38_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Pu39_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Pu40_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Pu41_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Pu42_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Pu43_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Am41_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_As42_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Am42_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Am43_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Am44_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Cm42_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Cm43_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Cm44_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_I35_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Xe35_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Nd47_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Nd48_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Nd49_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Pm47_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Ps48_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Pm48_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Pm49_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Sm47_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Sm48_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Sm49_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Gd52_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Gd54_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Gd55_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Gd56_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Gd57_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Gd58_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_Gd60_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: N_tmp_RST

      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_U34_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_U35_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_U36_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_U37_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_U38_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Np37_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Np38_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Np39_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Pu38_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Pu39_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Pu40_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Pu41_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Pu42_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Pu43_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Am41_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_As42_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Am42_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Am43_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Am44_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Cm42_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Cm43_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Cm44_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_U34_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_U35_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_U36_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_U37_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_U38_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Np37_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Np38_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Np39_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Pu38_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Pu39_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Pu40_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Pu41_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Pu42_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Pu43_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Am41_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_As42_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Am42_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Am43_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Am44_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Cm42_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Cm43_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_f_Cm44_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_I35_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Xe35_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Nd47_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Nd48_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Nd49_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Pm47_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Ps48_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Pm48_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Pm49_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Sm47_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Sm48_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Sm49_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Gd52_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Gd54_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Gd55_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Gd56_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Gd57_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Gd58_RST
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: miXS_a_Gd60_RST

      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Y_I35_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Y_Xe35_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Y_Nd47_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Y_Nd48_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Y_Nd49_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Y_Pm47_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Y_Ps48_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Y_Pm48_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Y_Pm49_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Y_Sm49_RST
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Y_Pm49_eff_RST

      INTEGER, DIMENSION(:, :), ALLOCATABLE :: I_4N_RST
      REAL(8), DIMENSION(:), ALLOCATABLE :: CR_Bot_RST0
      REAL(8), DIMENSION(:), ALLOCATABLE :: dTemp
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: dCR

      REAL(8) :: TM_In_RST0

      REAL(8), DIMENSION(:, :), ALLOCATABLE :: T_Fuel_RST0
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: T_Mod_RST0
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: D_Mod_RST0

      INTEGER :: N_Branch
      INTEGER :: I_Branch
      INTEGER, DIMENSION(:), ALLOCATABLE :: Type_Branch
#ifdef tuan_fr
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE:: N_FR_RST

#endif
      END MODULE Inc_RST
