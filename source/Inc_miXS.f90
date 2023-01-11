
      module Inc_miXS
      implicit none

      real(8), allocatable :: miXS_tr_U34  (:,:,:)
      real(8), allocatable :: miXS_tr_U35  (:,:,:)
      real(8), allocatable :: miXS_tr_U36  (:,:,:)
      real(8), allocatable :: miXS_tr_U37  (:,:,:)
      real(8), allocatable :: miXS_tr_U38  (:,:,:)
      real(8), allocatable :: miXS_tr_Np37 (:,:,:)
      real(8), allocatable :: miXS_tr_Np38 (:,:,:)
      real(8), allocatable :: miXS_tr_Np39 (:,:,:)
      real(8), allocatable :: miXS_tr_Pu38 (:,:,:)
      real(8), allocatable :: miXS_tr_Pu39 (:,:,:)
      real(8), allocatable :: miXS_tr_Pu40 (:,:,:)
      real(8), allocatable :: miXS_tr_Pu41 (:,:,:)
      real(8), allocatable :: miXS_tr_Pu42 (:,:,:)
      real(8), allocatable :: miXS_tr_Pu43 (:,:,:)
      real(8), allocatable :: miXS_tr_Am41 (:,:,:)
      real(8), allocatable :: miXS_tr_As42 (:,:,:)
      real(8), allocatable :: miXS_tr_Am42 (:,:,:)
      real(8), allocatable :: miXS_tr_Am43 (:,:,:)
      real(8), allocatable :: miXS_tr_Am44 (:,:,:)
      real(8), allocatable :: miXS_tr_Cm42 (:,:,:)
      real(8), allocatable :: miXS_tr_Cm43 (:,:,:)
      real(8), allocatable :: miXS_tr_Cm44 (:,:,:)

      real(8), allocatable :: miXS_a_U34  (:,:,:)
      real(8), allocatable :: miXS_a_U35  (:,:,:)
      real(8), allocatable :: miXS_a_U36  (:,:,:)
      real(8), allocatable :: miXS_a_U37  (:,:,:)
      real(8), allocatable :: miXS_a_U38  (:,:,:)
      real(8), allocatable :: miXS_a_Np37 (:,:,:)
      real(8), allocatable :: miXS_a_Np38 (:,:,:)
      real(8), allocatable :: miXS_a_Np39 (:,:,:)
      real(8), allocatable :: miXS_a_Pu38 (:,:,:)
      real(8), allocatable :: miXS_a_Pu39 (:,:,:)
      real(8), allocatable :: miXS_a_Pu40 (:,:,:)
      real(8), allocatable :: miXS_a_Pu41 (:,:,:)
      real(8), allocatable :: miXS_a_Pu42 (:,:,:)
      real(8), allocatable :: miXS_a_Pu43 (:,:,:)
      real(8), allocatable :: miXS_a_Am41 (:,:,:)
      real(8), allocatable :: miXS_a_As42 (:,:,:)
      real(8), allocatable :: miXS_a_Am42 (:,:,:)
      real(8), allocatable :: miXS_a_Am43 (:,:,:)
      real(8), allocatable :: miXS_a_Am44 (:,:,:)
      real(8), allocatable :: miXS_a_Cm42 (:,:,:)
      real(8), allocatable :: miXS_a_Cm43 (:,:,:)
      real(8), allocatable :: miXS_a_Cm44 (:,:,:)

      real(8), allocatable :: miXS_a_U34_Old  (:,:,:)
      real(8), allocatable :: miXS_a_U35_Old  (:,:,:)
      real(8), allocatable :: miXS_a_U36_Old  (:,:,:)
      real(8), allocatable :: miXS_a_U37_Old  (:,:,:)
      real(8), allocatable :: miXS_a_U38_Old  (:,:,:)
      real(8), allocatable :: miXS_a_Np37_Old (:,:,:)
      real(8), allocatable :: miXS_a_Np38_Old (:,:,:)
      real(8), allocatable :: miXS_a_Np39_Old (:,:,:)
      real(8), allocatable :: miXS_a_Pu38_Old (:,:,:)
      real(8), allocatable :: miXS_a_Pu39_Old (:,:,:)
      real(8), allocatable :: miXS_a_Pu40_Old (:,:,:)
      real(8), allocatable :: miXS_a_Pu41_Old (:,:,:)
      real(8), allocatable :: miXS_a_Pu42_Old (:,:,:)
      real(8), allocatable :: miXS_a_Pu43_Old (:,:,:)
      real(8), allocatable :: miXS_a_Am41_Old (:,:,:)
      real(8), allocatable :: miXS_a_As42_Old (:,:,:)
      real(8), allocatable :: miXS_a_Am42_Old (:,:,:)
      real(8), allocatable :: miXS_a_Am43_Old (:,:,:)
      real(8), allocatable :: miXS_a_Am44_Old (:,:,:)
      real(8), allocatable :: miXS_a_Cm42_Old (:,:,:)
      real(8), allocatable :: miXS_a_Cm43_Old (:,:,:)
      real(8), allocatable :: miXS_a_Cm44_Old (:,:,:)

      real(8), allocatable :: miXS_f_U34  (:,:,:)
      real(8), allocatable :: miXS_f_U35  (:,:,:)
      real(8), allocatable :: miXS_f_U36  (:,:,:)
      real(8), allocatable :: miXS_f_U37  (:,:,:)
      real(8), allocatable :: miXS_f_U38  (:,:,:)
      real(8), allocatable :: miXS_f_Np37 (:,:,:)
      real(8), allocatable :: miXS_f_Np38 (:,:,:)
      real(8), allocatable :: miXS_f_Np39 (:,:,:)
      real(8), allocatable :: miXS_f_Pu38 (:,:,:)
      real(8), allocatable :: miXS_f_Pu39 (:,:,:)
      real(8), allocatable :: miXS_f_Pu40 (:,:,:)
      real(8), allocatable :: miXS_f_Pu41 (:,:,:)
      real(8), allocatable :: miXS_f_Pu42 (:,:,:)
      real(8), allocatable :: miXS_f_Pu43 (:,:,:)
      real(8), allocatable :: miXS_f_Am41 (:,:,:)
      real(8), allocatable :: miXS_f_As42 (:,:,:)
      real(8), allocatable :: miXS_f_Am42 (:,:,:)
      real(8), allocatable :: miXS_f_Am43 (:,:,:)
      real(8), allocatable :: miXS_f_Am44 (:,:,:)
      real(8), allocatable :: miXS_f_Cm42 (:,:,:)
      real(8), allocatable :: miXS_f_Cm43 (:,:,:)
      real(8), allocatable :: miXS_f_Cm44 (:,:,:)

      real(8), allocatable :: miXS_f_U34_Old  (:,:,:)
      real(8), allocatable :: miXS_f_U35_Old  (:,:,:)
      real(8), allocatable :: miXS_f_U36_Old  (:,:,:)
      real(8), allocatable :: miXS_f_U37_Old  (:,:,:)
      real(8), allocatable :: miXS_f_U38_Old  (:,:,:)
      real(8), allocatable :: miXS_f_Np37_Old (:,:,:)
      real(8), allocatable :: miXS_f_Np38_Old (:,:,:)
      real(8), allocatable :: miXS_f_Np39_Old (:,:,:)
      real(8), allocatable :: miXS_f_Pu38_Old (:,:,:)
      real(8), allocatable :: miXS_f_Pu39_Old (:,:,:)
      real(8), allocatable :: miXS_f_Pu40_Old (:,:,:)
      real(8), allocatable :: miXS_f_Pu41_Old (:,:,:)
      real(8), allocatable :: miXS_f_Pu42_Old (:,:,:)
      real(8), allocatable :: miXS_f_Pu43_Old (:,:,:)
      real(8), allocatable :: miXS_f_Am41_Old (:,:,:)
      real(8), allocatable :: miXS_f_As42_Old (:,:,:)
      real(8), allocatable :: miXS_f_Am42_Old (:,:,:)
      real(8), allocatable :: miXS_f_Am43_Old (:,:,:)
      real(8), allocatable :: miXS_f_Am44_Old (:,:,:)
      real(8), allocatable :: miXS_f_Cm42_Old (:,:,:)
      real(8), allocatable :: miXS_f_Cm43_Old (:,:,:)
      real(8), allocatable :: miXS_f_Cm44_Old (:,:,:)

      real(8), allocatable :: nu_miXS_f_U34  (:,:,:)
      real(8), allocatable :: nu_miXS_f_U35  (:,:,:)
      real(8), allocatable :: nu_miXS_f_U36  (:,:,:)
      real(8), allocatable :: nu_miXS_f_U37  (:,:,:)
      real(8), allocatable :: nu_miXS_f_U38  (:,:,:)
      real(8), allocatable :: nu_miXS_f_Np37 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Np38 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Np39 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Pu38 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Pu39 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Pu40 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Pu41 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Pu42 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Pu43 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Am41 (:,:,:)
      real(8), allocatable :: nu_miXS_f_As42 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Am42 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Am43 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Am44 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Cm42 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Cm43 (:,:,:)
      real(8), allocatable :: nu_miXS_f_Cm44 (:,:,:)

      real(8), allocatable :: kap_miXS_f_U34  (:,:,:)
      real(8), allocatable :: kap_miXS_f_U35  (:,:,:)
      real(8), allocatable :: kap_miXS_f_U36  (:,:,:)
      real(8), allocatable :: kap_miXS_f_U37  (:,:,:)
      real(8), allocatable :: kap_miXS_f_U38  (:,:,:)
      real(8), allocatable :: kap_miXS_f_Np37 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Np38 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Np39 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Pu38 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Pu39 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Pu40 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Pu41 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Pu42 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Pu43 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Am41 (:,:,:)
      real(8), allocatable :: kap_miXS_f_As42 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Am42 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Am43 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Am44 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Cm42 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Cm43 (:,:,:)
      real(8), allocatable :: kap_miXS_f_Cm44 (:,:,:)

      real(8), allocatable :: miXS_s_U34  (:,:,:)
      real(8), allocatable :: miXS_s_U35  (:,:,:)
      real(8), allocatable :: miXS_s_U36  (:,:,:)
      real(8), allocatable :: miXS_s_U37  (:,:,:)
      real(8), allocatable :: miXS_s_U38  (:,:,:)
      real(8), allocatable :: miXS_s_Np37 (:,:,:)
      real(8), allocatable :: miXS_s_Np38 (:,:,:)
      real(8), allocatable :: miXS_s_Np39 (:,:,:)
      real(8), allocatable :: miXS_s_Pu38 (:,:,:)
      real(8), allocatable :: miXS_s_Pu39 (:,:,:)
      real(8), allocatable :: miXS_s_Pu40 (:,:,:)
      real(8), allocatable :: miXS_s_Pu41 (:,:,:)
      real(8), allocatable :: miXS_s_Pu42 (:,:,:)
      real(8), allocatable :: miXS_s_Pu43 (:,:,:)
      real(8), allocatable :: miXS_s_Am41 (:,:,:)
      real(8), allocatable :: miXS_s_As42 (:,:,:)
      real(8), allocatable :: miXS_s_Am42 (:,:,:)
      real(8), allocatable :: miXS_s_Am43 (:,:,:)
      real(8), allocatable :: miXS_s_Am44 (:,:,:)
      real(8), allocatable :: miXS_s_Cm42 (:,:,:)
      real(8), allocatable :: miXS_s_Cm43 (:,:,:)
      real(8), allocatable :: miXS_s_Cm44 (:,:,:)

      real(8), allocatable :: miXS_n2n_U35  (:,:,:)
      real(8), allocatable :: miXS_n2n_U38  (:,:,:)
      real(8), allocatable :: miXS_n2n_Np37 (:,:,:)
      real(8), allocatable :: miXS_n2n_Pu39 (:,:,:)

      real(8), allocatable :: miXS_n2n_U35_Old  (:,:,:)
      real(8), allocatable :: miXS_n2n_U38_Old  (:,:,:)
      real(8), allocatable :: miXS_n2n_Np37_Old (:,:,:)
      real(8), allocatable :: miXS_n2n_Pu39_Old (:,:,:)

      real(8), allocatable :: miXS_tr_I35  (:,:,:)
      real(8), allocatable :: miXS_tr_Xe35 (:,:,:)
      real(8), allocatable :: miXS_tr_Nd47 (:,:,:)
      real(8), allocatable :: miXS_tr_Nd48 (:,:,:)
      real(8), allocatable :: miXS_tr_Nd49 (:,:,:)
      real(8), allocatable :: miXS_tr_Pm47 (:,:,:)
      real(8), allocatable :: miXS_tr_Ps48 (:,:,:)
      real(8), allocatable :: miXS_tr_Pm48 (:,:,:)
      real(8), allocatable :: miXS_tr_Pm49 (:,:,:)
      real(8), allocatable :: miXS_tr_Sm47 (:,:,:)
      real(8), allocatable :: miXS_tr_Sm48 (:,:,:)
      real(8), allocatable :: miXS_tr_Sm49 (:,:,:)

      real(8), allocatable :: miXS_a_I35  (:,:,:)
      real(8), allocatable :: miXS_a_Xe35 (:,:,:)
      real(8), allocatable :: miXS_a_Nd47 (:,:,:)
      real(8), allocatable :: miXS_a_Nd48 (:,:,:)
      real(8), allocatable :: miXS_a_Nd49 (:,:,:)
      real(8), allocatable :: miXS_a_Pm47 (:,:,:)
      real(8), allocatable :: miXS_a_Ps48 (:,:,:)
      real(8), allocatable :: miXS_a_Pm48 (:,:,:)
      real(8), allocatable :: miXS_a_Pm49 (:,:,:)
      real(8), allocatable :: miXS_a_Sm47 (:,:,:)
      real(8), allocatable :: miXS_a_Sm48 (:,:,:)
      real(8), allocatable :: miXS_a_Sm49 (:,:,:)

      real(8), allocatable :: miXS_a_I35_Old  (:,:,:)
      real(8), allocatable :: miXS_a_Xe35_Old (:,:,:)
      real(8), allocatable :: miXS_a_Nd47_Old (:,:,:)
      real(8), allocatable :: miXS_a_Nd48_Old (:,:,:)
      real(8), allocatable :: miXS_a_Nd49_Old (:,:,:)
      real(8), allocatable :: miXS_a_Pm47_Old (:,:,:)
      real(8), allocatable :: miXS_a_Ps48_Old (:,:,:)
      real(8), allocatable :: miXS_a_Pm48_Old (:,:,:)
      real(8), allocatable :: miXS_a_Pm49_Old (:,:,:)
      real(8), allocatable :: miXS_a_Sm47_Old (:,:,:)
      real(8), allocatable :: miXS_a_Sm48_Old (:,:,:)
      real(8), allocatable :: miXS_a_Sm49_Old (:,:,:)

      real(8), allocatable :: miXS_s_I35  (:,:,:)
      real(8), allocatable :: miXS_s_Xe35 (:,:,:)
      real(8), allocatable :: miXS_s_Nd47 (:,:,:)
      real(8), allocatable :: miXS_s_Nd48 (:,:,:)
      real(8), allocatable :: miXS_s_Nd49 (:,:,:)
      real(8), allocatable :: miXS_s_Pm47 (:,:,:)
      real(8), allocatable :: miXS_s_Ps48 (:,:,:)
      real(8), allocatable :: miXS_s_Pm48 (:,:,:)
      real(8), allocatable :: miXS_s_Pm49 (:,:,:)
      real(8), allocatable :: miXS_s_Sm47 (:,:,:)
      real(8), allocatable :: miXS_s_Sm48 (:,:,:)
      real(8), allocatable :: miXS_s_Sm49 (:,:,:)

      real(8), allocatable :: miXS_s_I35_Old  (:,:,:)
      real(8), allocatable :: miXS_s_Xe35_Old (:,:,:)
      real(8), allocatable :: miXS_s_Nd47_Old (:,:,:)
      real(8), allocatable :: miXS_s_Nd48_Old (:,:,:)
      real(8), allocatable :: miXS_s_Nd49_Old (:,:,:)
      real(8), allocatable :: miXS_s_Pm47_Old (:,:,:)
      real(8), allocatable :: miXS_s_Ps48_Old (:,:,:)
      real(8), allocatable :: miXS_s_Pm48_Old (:,:,:)
      real(8), allocatable :: miXS_s_Pm49_Old (:,:,:)
      real(8), allocatable :: miXS_s_Sm47_Old (:,:,:)
      real(8), allocatable :: miXS_s_Sm48_Old (:,:,:)
      real(8), allocatable :: miXS_s_Sm49_Old (:,:,:)

      real(8), allocatable :: miXS_tr_Gd52 (:,:,:)
      real(8), allocatable :: miXS_tr_Gd54 (:,:,:)
      real(8), allocatable :: miXS_tr_Gd55 (:,:,:)
      real(8), allocatable :: miXS_tr_Gd56 (:,:,:)
      real(8), allocatable :: miXS_tr_Gd57 (:,:,:)
      real(8), allocatable :: miXS_tr_Gd58 (:,:,:)
      real(8), allocatable :: miXS_tr_Gd60 (:,:,:)

      real(8), allocatable :: miXS_a_Gd52 (:,:,:)
      real(8), allocatable :: miXS_a_Gd54 (:,:,:)
      real(8), allocatable :: miXS_a_Gd55 (:,:,:)
      real(8), allocatable :: miXS_a_Gd56 (:,:,:)
      real(8), allocatable :: miXS_a_Gd57 (:,:,:)
      real(8), allocatable :: miXS_a_Gd58 (:,:,:)
      real(8), allocatable :: miXS_a_Gd60 (:,:,:)

      real(8), allocatable :: miXS_a_Gd52_Old (:,:,:)
      real(8), allocatable :: miXS_a_Gd54_Old (:,:,:)
      real(8), allocatable :: miXS_a_Gd55_Old (:,:,:)
      real(8), allocatable :: miXS_a_Gd56_Old (:,:,:)
      real(8), allocatable :: miXS_a_Gd57_Old (:,:,:)
      real(8), allocatable :: miXS_a_Gd58_Old (:,:,:)
      real(8), allocatable :: miXS_a_Gd60_Old (:,:,:)

      real(8), allocatable :: miXS_s_Gd52 (:,:,:)
      real(8), allocatable :: miXS_s_Gd54 (:,:,:)
      real(8), allocatable :: miXS_s_Gd55 (:,:,:)
      real(8), allocatable :: miXS_s_Gd56 (:,:,:)
      real(8), allocatable :: miXS_s_Gd57 (:,:,:)
      real(8), allocatable :: miXS_s_Gd58 (:,:,:)
      real(8), allocatable :: miXS_s_Gd60 (:,:,:)

      real(8), allocatable :: miXS_tr_H2O    (:,:,:)
      real(8), allocatable :: miXS_tr_B0     (:,:,:)
      real(8), allocatable :: miXS_a_H2O     (:,:,:)
      real(8), allocatable :: miXS_a_B0      (:,:,:)
      real(8), allocatable :: miXS_a_H2O_Old (:,:,:)
      real(8), allocatable :: miXS_a_B0_Old  (:,:,:)
      real(8), allocatable :: miXS_s_H2O     (:,:,:)
      real(8), allocatable :: miXS_s_B0      (:,:,:)

      real(8), allocatable :: Y_I35      (:,:)
      real(8), allocatable :: Y_Xe35     (:,:)
      real(8), allocatable :: Y_Nd47     (:,:)
      real(8), allocatable :: Y_Nd48     (:,:)
      real(8), allocatable :: Y_Nd49     (:,:)
      real(8), allocatable :: Y_Pm47     (:,:)
      real(8), allocatable :: Y_Ps48     (:,:)
      real(8), allocatable :: Y_Pm48     (:,:)
      real(8), allocatable :: Y_Pm49     (:,:)
      real(8), allocatable :: Y_Sm49     (:,:)
      real(8), allocatable :: Y_Xe35_eff (:,:)
      real(8), allocatable :: Y_Pm49_eff (:,:)

      real(8), allocatable :: Y_I35_Old      (:,:)
      real(8), allocatable :: Y_Xe35_Old     (:,:)
      real(8), allocatable :: Y_Nd47_Old     (:,:)
      real(8), allocatable :: Y_Nd48_Old     (:,:)
      real(8), allocatable :: Y_Nd49_Old     (:,:)
      real(8), allocatable :: Y_Pm47_Old     (:,:)
      real(8), allocatable :: Y_Ps48_Old     (:,:)
      real(8), allocatable :: Y_Pm48_Old     (:,:)
      real(8), allocatable :: Y_Pm49_Old     (:,:)
      real(8), allocatable :: Y_Sm49_Old     (:,:)
      real(8), allocatable :: Y_Xe35_eff_Old (:,:)
      real(8), allocatable :: Y_Pm49_eff_Old (:,:)

      real(8), allocatable :: miXS_a_A1      (:,:,:)
      real(8), allocatable :: miXS_a_A2      (:,:,:)
      real(8), allocatable :: miXS_a_A2m     (:,:,:)
      real(8), allocatable :: miXS_a_A3      (:,:,:)
      real(8), allocatable :: miXS_a_B2      (:,:,:)
      real(8), allocatable :: miXS_a_B3      (:,:,:)
      real(8), allocatable :: miXS_a_A4      (:,:,:)
      real(8), allocatable :: miXS_a_V1      (:,:,:)
      real(8), allocatable :: miXS_a_V2      (:,:,:)
      real(8), allocatable :: miXS_a_Cr      (:,:,:)

      real(8), allocatable :: DET_miXS_a_A1  (:,:,:)
      real(8), allocatable :: DET_miXS_a_A2  (:,:,:)
      real(8), allocatable :: DET_miXS_a_A2m (:,:,:)
      real(8), allocatable :: DET_miXS_a_A3  (:,:,:)
      real(8), allocatable :: DET_miXS_a_B2  (:,:,:)
      real(8), allocatable :: DET_miXS_a_B3  (:,:,:)
      real(8), allocatable :: DET_miXS_a_A4  (:,:,:)
      real(8), allocatable :: DET_miXS_a_V1  (:,:,:)
      real(8), allocatable :: DET_miXS_a_V2  (:,:,:)
      real(8), allocatable :: DET_miXS_a_Cr  (:,:,:)

      real(8), allocatable :: miXS_a_B10      (:,:,:)
      real(8), allocatable :: miXS_a_B10_Old  (:,:,:)
      real(8), allocatable :: XS_DEP_B10      (:,:,:,:)
      real(8), allocatable :: XS_DEP_B10_BASE (:,:,:)
      real(8), allocatable :: temp_DEP_B10    (:)

#ifdef tuan_fr

      REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE:: miXS_Hex_ref
      REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE:: miXS_Hex
      integer, DIMENSION(:), ALLOCATABLE:: ZAID
      REAL(8), DIMENSION(:,:), ALLOCATABLE:: Dens_Int
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE:: N_FR
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE:: N_FR_Old
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE:: N_FR_Predictor
      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE:: XS_residual_ref
      REAL(8), DIMENSION(:,:,:),   ALLOCATABLE:: XS_residual
!      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE:: XS_residual


      integer  ::  N_XSfile=0
      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: miXS_n2n
      integer  ::  N_Nuclide_FR=0


#endif

#ifdef tuan_fr_tdep

      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE:: N_FR_T
      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE:: N_FR_T_Old
      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE:: N_FR_T_Predictor
#endif



      end module Inc_miXS
