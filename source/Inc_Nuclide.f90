
      module Inc_Nuclide
      ! *** Variable ***
      ! [Heavy Nuclide]
      ! 2-dimension :: ( Ixy, Iz )
      ! @ = 22 Heavy Nuclide
      !   1)  U-234
      !   2)  U-235
      !   3)  U-236
      !   4)  U-237
      !   5)  U-238
      !   6)  Np-237
      !   7)  Np-238
      !   8)  Np-239
      !   9)  Pu-238
      !   10) Pu-239
      !   11) Pu-240
      !   12) Pu-241
      !   13) Pu-242
      !   14) Pu-243
      !   15) Am-241
      !   16) Am-242  (Stable)
      !   17) Am-242m (Meta Stable)
      !   18) Am-243
      !   19) Am-244
      !   20) Cm-242
      !   21) Cm-243
      !   22) Cm-244
      !
      ! 1) N_Heavy :: Total Number of Heavy Nuclide
      ! 2) N_@     :: Number Density of @ Heavy Nuclide
      ! 3) N_@_Old :: Previous Burnup Step Number Density of @ Heavy Nuclide
      !
      ! [Fission Product]
      ! 2-dimension :: ( Ixy, Iz )
      ! @ = 12 Fission Product
      !   1)  I-135
      !   2)  Xe-135
      !   3)  Nd-147
      !   4)  Nd-148
      !   5)  Nd-149
      !   6)  Pm-147
      !   7)  Pm-148  (Stable)
      !   8)  Pm-148m (Meta Stable)
      !   9)  Pm-149
      !   10) Sm-147
      !   11) Sm-148
      !   12) Sm-149
      !
      ! 1) N_NdSm          :: Total Number of Nd-Sm Chain Nuclide
      ! 2) N_@             :: Number Density of @ Fission Product
      ! 3) N_@_Old         :: Previous Burnup Step Number Density of @ Fission Product
      !
      ! [Burnable Absorber]
      ! 2-dimension :: ( Ixy, Iz )
      ! @ = 7 Burnable Absorber
      !   1) Gd-152
      !   2) Gd-154
      !   3) Gd-155
      !   4) Gd-156
      !   5) Gd-157
      !   6) Gd-158
      !   7) Gd-160
      !
      ! 1) N_Gd     :: Total Number of Gd in Linear Depletion Chain of BP (Gd-154~158; = 5)
      ! 2) N_@      :: Number Density of @ Burnable Absorber
      ! 3) N_@_Old  :: Previous Burnup Step N_@
      ! 4) N_@_OOld :: 2nd Previous Burnup Step N_@
      ! 5) N_@_u    :: Post Corrected N_@_p
      !
      ! *** Constant ***
      ! 1)  M_B0      :: Molar Mass of Natural Boron (B-0)
      ! 2)  M_UO2     :: Molar Mass of UO2   [g/mol]
      ! 3)  D_UO2     :: Initial UO2 Density [g/cc ]
      ! 4)  w_Ps48    :: Branching Ratio [Pm-147  + n ] => [Pm-148 ] (Absorption)
      ! 5)  w_Pm48    :: Branching Ratio [Pm-147  + n ] => [Pm-148m] (Absorption)
      ! 6)  w_Ps48_IT :: Branching Ratio [Pm-148m - E ] => [Pm-148 ] (Isometric Transition)
      ! 7)  w_Sm48    :: Branching Ratio [Pm-148m - e-] => [Sm-148 ] (Beta Decay)
      ! 8)  w_U37     :: Branching Ratio [Pu-241  - He] => [ U-237 ] (Alpha Decay)
      ! 9)  w_Am41    :: Branching Ratio [Pu-241  - e-] => [Am-241 ] (Beta Decay)
      ! 10) w_As42    :: Branching Ratio [Am-241  + n ] => [Am-242 ] (Absorption)
      ! 11) w_Am42    :: Branching Ratio [Am-241  + n ] => [Am-242m] (Absorption)
      ! 12) w_Pu42    :: Branching Ratio [Am-242  + e-] => [Pu-242 ] (Electron Capture)
      ! 13) w_Cm42    :: Branching Ratio [Am-242  - e-] => [Cm-242 ] (Beta Decay)
      ! 14) w_As42_IT :: Branching Ratio [Am-242m - E ] => [Am-242 ] (Isometric Transition)
      ! 15) w_Np38    :: Branching Ratio [Am-242m - He] => [Np-238 ] (Alpha Decay)
      ! 16) beta_@    :: Beta  Decay Constant of @ Nuclide
      ! 17) alpha_@   :: Alpha Decay Constant of @ Nuclide
      ! 18) EC_@      :: Electron  Capture    Decay Constant of @ Nuclide
      ! 19) IT_@      :: Isometric Transition Decay Constant of @ Nuclide
      ! xx) w_Np37    :: Branching Ratio [Am-241  - He] => [Np-237 ] (Alpha Decay)

      implicit none

      real(8), dimension(:,:), allocatable :: N_U34
      real(8), dimension(:,:), allocatable :: N_U35
      real(8), dimension(:,:), allocatable :: N_U36
      real(8), dimension(:,:), allocatable :: N_U37
      real(8), dimension(:,:), allocatable :: N_U38
      real(8), dimension(:,:), allocatable :: N_Np37
      real(8), dimension(:,:), allocatable :: N_Np38
      real(8), dimension(:,:), allocatable :: N_Np39
      real(8), dimension(:,:), allocatable :: N_Pu38
      real(8), dimension(:,:), allocatable :: N_Pu39
      real(8), dimension(:,:), allocatable :: N_Pu40
      real(8), dimension(:,:), allocatable :: N_Pu41
      real(8), dimension(:,:), allocatable :: N_Pu42
      real(8), dimension(:,:), allocatable :: N_Pu43
      real(8), dimension(:,:), allocatable :: N_Am41
      real(8), dimension(:,:), allocatable :: N_As42
      real(8), dimension(:,:), allocatable :: N_Am42
      real(8), dimension(:,:), allocatable :: N_Am43
      real(8), dimension(:,:), allocatable :: N_Am44
      real(8), dimension(:,:), allocatable :: N_Cm42
      real(8), dimension(:,:), allocatable :: N_Cm43
      real(8), dimension(:,:), allocatable :: N_Cm44

      real(8), dimension(:,:), allocatable :: N_U34_Old
      real(8), dimension(:,:), allocatable :: N_U35_Old
      real(8), dimension(:,:), allocatable :: N_U36_Old
      real(8), dimension(:,:), allocatable :: N_U37_Old
      real(8), dimension(:,:), allocatable :: N_U38_Old
      real(8), dimension(:,:), allocatable :: N_Np37_Old
      real(8), dimension(:,:), allocatable :: N_Np38_Old
      real(8), dimension(:,:), allocatable :: N_Np39_Old
      real(8), dimension(:,:), allocatable :: N_Pu38_Old
      real(8), dimension(:,:), allocatable :: N_Pu39_Old
      real(8), dimension(:,:), allocatable :: N_Pu40_Old
      real(8), dimension(:,:), allocatable :: N_Pu41_Old
      real(8), dimension(:,:), allocatable :: N_Pu42_Old
      real(8), dimension(:,:), allocatable :: N_Pu43_Old
      real(8), dimension(:,:), allocatable :: N_Am41_Old
      real(8), dimension(:,:), allocatable :: N_As42_Old
      real(8), dimension(:,:), allocatable :: N_Am42_Old
      real(8), dimension(:,:), allocatable :: N_Am43_Old
      real(8), dimension(:,:), allocatable :: N_Am44_Old
      real(8), dimension(:,:), allocatable :: N_Cm42_Old
      real(8), dimension(:,:), allocatable :: N_Cm43_Old
      real(8), dimension(:,:), allocatable :: N_Cm44_Old

      real(8), dimension(:,:), allocatable :: N_U34_Predictor
      real(8), dimension(:,:), allocatable :: N_U35_Predictor
      real(8), dimension(:,:), allocatable :: N_U36_Predictor
      real(8), dimension(:,:), allocatable :: N_U37_Predictor
      real(8), dimension(:,:), allocatable :: N_U38_Predictor
      real(8), dimension(:,:), allocatable :: N_Np37_Predictor
      real(8), dimension(:,:), allocatable :: N_Np38_Predictor
      real(8), dimension(:,:), allocatable :: N_Np39_Predictor
      real(8), dimension(:,:), allocatable :: N_Pu38_Predictor
      real(8), dimension(:,:), allocatable :: N_Pu39_Predictor
      real(8), dimension(:,:), allocatable :: N_Pu40_Predictor
      real(8), dimension(:,:), allocatable :: N_Pu41_Predictor
      real(8), dimension(:,:), allocatable :: N_Pu42_Predictor
      real(8), dimension(:,:), allocatable :: N_Pu43_Predictor
      real(8), dimension(:,:), allocatable :: N_Am41_Predictor
      real(8), dimension(:,:), allocatable :: N_As42_Predictor
      real(8), dimension(:,:), allocatable :: N_Am42_Predictor
      real(8), dimension(:,:), allocatable :: N_Am43_Predictor
      real(8), dimension(:,:), allocatable :: N_Am44_Predictor
      real(8), dimension(:,:), allocatable :: N_Cm42_Predictor
      real(8), dimension(:,:), allocatable :: N_Cm43_Predictor
      real(8), dimension(:,:), allocatable :: N_Cm44_Predictor

      integer :: N_Heavy

      real(8), dimension(:,:), allocatable :: N_I35
      real(8), dimension(:,:), allocatable :: N_Xe35
      real(8), dimension(:,:), allocatable :: N_Nd47
      real(8), dimension(:,:), allocatable :: N_Nd48
      real(8), dimension(:,:), allocatable :: N_Nd49
      real(8), dimension(:,:), allocatable :: N_Pm47
      real(8), dimension(:,:), allocatable :: N_Ps48
      real(8), dimension(:,:), allocatable :: N_Pm48
      real(8), dimension(:,:), allocatable :: N_Pm49
      real(8), dimension(:,:), allocatable :: N_Sm47
      real(8), dimension(:,:), allocatable :: N_Sm48
      real(8), dimension(:,:), allocatable :: N_Sm49

      real(8), dimension(:,:), allocatable :: N_I35_Old
      real(8), dimension(:,:), allocatable :: N_Xe35_Old
      real(8), dimension(:,:), allocatable :: N_Nd47_Old
      real(8), dimension(:,:), allocatable :: N_Nd48_Old
      real(8), dimension(:,:), allocatable :: N_Nd49_Old
      real(8), dimension(:,:), allocatable :: N_Pm47_Old
      real(8), dimension(:,:), allocatable :: N_Ps48_Old
      real(8), dimension(:,:), allocatable :: N_Pm48_Old
      real(8), dimension(:,:), allocatable :: N_Pm49_Old
      real(8), dimension(:,:), allocatable :: N_Sm47_Old
      real(8), dimension(:,:), allocatable :: N_Sm48_Old
      real(8), dimension(:,:), allocatable :: N_Sm49_Old

      real(8), dimension(:,:), allocatable :: N_I35_Predictor
      real(8), dimension(:,:), allocatable :: N_Xe35_Predictor
      real(8), dimension(:,:), allocatable :: N_Nd47_Predictor
      real(8), dimension(:,:), allocatable :: N_Nd48_Predictor
      real(8), dimension(:,:), allocatable :: N_Nd49_Predictor
      real(8), dimension(:,:), allocatable :: N_Pm47_Predictor
      real(8), dimension(:,:), allocatable :: N_Ps48_Predictor
      real(8), dimension(:,:), allocatable :: N_Pm48_Predictor
      real(8), dimension(:,:), allocatable :: N_Pm49_Predictor
      real(8), dimension(:,:), allocatable :: N_Sm47_Predictor
      real(8), dimension(:,:), allocatable :: N_Sm48_Predictor
      real(8), dimension(:,:), allocatable :: N_Sm49_Predictor

      integer :: N_NdSm

      real(8), dimension(:,:), allocatable :: N_Gd52
      real(8), dimension(:,:), allocatable :: N_Gd54
      real(8), dimension(:,:), allocatable :: N_Gd55
      real(8), dimension(:,:), allocatable :: N_Gd56
      real(8), dimension(:,:), allocatable :: N_Gd57
      real(8), dimension(:,:), allocatable :: N_Gd58
      real(8), dimension(:,:), allocatable :: N_Gd60

      real(8), dimension(:,:), allocatable :: N_Gd52_Old
      real(8), dimension(:,:), allocatable :: N_Gd54_Old
      real(8), dimension(:,:), allocatable :: N_Gd55_Old
      real(8), dimension(:,:), allocatable :: N_Gd56_Old
      real(8), dimension(:,:), allocatable :: N_Gd57_Old
      real(8), dimension(:,:), allocatable :: N_Gd58_Old
      real(8), dimension(:,:), allocatable :: N_Gd60_Old

      real(8), dimension(:,:), allocatable :: N_Gd52_Predictor
      real(8), dimension(:,:), allocatable :: N_Gd54_Predictor
      real(8), dimension(:,:), allocatable :: N_Gd55_Predictor
      real(8), dimension(:,:), allocatable :: N_Gd56_Predictor
      real(8), dimension(:,:), allocatable :: N_Gd57_Predictor
      real(8), dimension(:,:), allocatable :: N_Gd58_Predictor
      real(8), dimension(:,:), allocatable :: N_Gd60_Predictor

      integer :: N_Gd

      real(8), dimension(:,:), allocatable :: N_H2O
      real(8), dimension(:,:), allocatable :: N_B0


      real(8), parameter :: M_B0    = 10.811004D0
      real(8), parameter :: M_UO2   = 267.033984D0
      real(8), parameter :: D_UO2   = 10.3402D0
      real(8), parameter :: M_Gd2O3 = 362.497D0

      real(8) :: w_Sm48    = 1d0-0.049d0
      real(8) :: w_Ps48    = 0.5330d0
      real(8) :: w_Pm48    = 1d0-0.5330d0
      real(8) :: w_U37     = 2.45006D-05
      real(8) :: w_As42_IT = 0.995D0 ! (w_As42_IT+w_Np38=1)
      real(8) :: w_Ps48_IT = 0.049D0
      real(8) :: w_As42    = 0.9190D0
      real(8) :: w_Am42    = 1d0 - 0.9190D0
      real(8) :: w_Am41    = 1d0 - 2.45006D-05
      real(8) :: w_Np38    = 0.005D0
      real(8) :: w_Pu42    = 0.173D0 ! (w_Pu42+w_Cm42=1)
      real(8) :: w_Cm42    = 0.827D0
      real(8) :: w_Np37    = 1d0

      real(8) :: beta_I35
      real(8) :: beta_Xe35
      real(8) :: beta_Nd47
      real(8) :: beta_Nd48
      real(8) :: beta_Nd49
      real(8) :: beta_Pm47
      real(8) :: beta_Ps48
      real(8) :: beta_Pm48
      real(8) :: beta_Pm49
      real(8) :: beta_Sm47
      real(8) :: beta_Sm48
      real(8) :: beta_Sm49
      real(8) :: beta_U37
      real(8) :: beta_Np38
      real(8) :: beta_Np39
      real(8) :: beta_Pu41
      real(8) :: beta_Pu43
      real(8) :: beta_As42
      real(8) :: beta_Am44

      real(8) :: alpha_Pu38
      real(8) :: alpha_Pu41
      real(8) :: alpha_Am42
      real(8) :: alpha_Cm42
      real(8) :: alpha_Cm43
      real(8) :: alpha_Cm44
      real(8) :: alpha_Am41

      real(8) :: EC_As42
      real(8) :: IT_Am42
      real(8) :: IT_Pm48

      real(8), dimension(:), allocatable :: HMMASS

      real(8), dimension(:,:), allocatable :: DET_N_A1
      real(8), dimension(:,:), allocatable :: DET_N_A2
      real(8), dimension(:,:), allocatable :: DET_N_A3
      real(8), dimension(:,:), allocatable :: DET_N_A4
      real(8), dimension(:,:), allocatable :: DET_N_A2m
      real(8), dimension(:,:), allocatable :: DET_N_B2
      real(8), dimension(:,:), allocatable :: DET_N_B3

      real(8), dimension(:,:), allocatable :: DET_N_A1_Old
      real(8), dimension(:,:), allocatable :: DET_N_A2_Old
      real(8), dimension(:,:), allocatable :: DET_N_A3_Old
      real(8), dimension(:,:), allocatable :: DET_N_A4_Old
      real(8), dimension(:,:), allocatable :: DET_N_A2m_Old
      real(8), dimension(:,:), allocatable :: DET_N_B2_Old
      real(8), dimension(:,:), allocatable :: DET_N_B3_Old

      real(8), dimension(:,:), allocatable :: DET_N_V1
      real(8), dimension(:,:), allocatable :: DET_N_V2
      real(8), dimension(:,:), allocatable :: DET_N_Cr52

      real(8), dimension(:,:), allocatable :: DET_N_V1_Old
      real(8), dimension(:,:), allocatable :: DET_N_V2_Old
      real(8), dimension(:,:), allocatable :: DET_N_Cr52_Old

      real(8), dimension(:,:), allocatable :: DET_I_beta_A2
      real(8), dimension(:,:), allocatable :: DET_I_beta_A2m
      real(8), dimension(:,:), allocatable :: DET_I_beta_A3
      real(8), dimension(:,:), allocatable :: DET_I_beta_B2
      real(8), dimension(:,:), allocatable :: DET_I_beta_B3
      real(8), dimension(:,:), allocatable :: DET_I_beta_rh
      real(8), dimension(:,:), allocatable :: DET_I_beta_rh_p
      real(8), dimension(:,:), allocatable :: DET_I_beta_rh_d

      real(8), dimension(:,:), allocatable :: DET_I_beta_V2
      real(8), dimension(:,:), allocatable :: DET_I_beta_si

      real(8), dimension(:,:), allocatable :: N_A1
      real(8), dimension(:,:), allocatable :: N_A2
      real(8), dimension(:,:), allocatable :: N_A3
      real(8), dimension(:,:), allocatable :: N_A4
      real(8), dimension(:,:), allocatable :: N_A2m
      real(8), dimension(:,:), allocatable :: N_B2
      real(8), dimension(:,:), allocatable :: N_B3

      real(8), dimension(:,:), allocatable :: N_A1_Old
      real(8), dimension(:,:), allocatable :: N_A2_Old
      real(8), dimension(:,:), allocatable :: N_A3_Old
      real(8), dimension(:,:), allocatable :: N_A4_Old
      real(8), dimension(:,:), allocatable :: N_A2m_Old
      real(8), dimension(:,:), allocatable :: N_B2_Old
      real(8), dimension(:,:), allocatable :: N_B3_Old

      real(8), dimension(:,:), allocatable :: N_V1
      real(8), dimension(:,:), allocatable :: N_V2
      real(8), dimension(:,:), allocatable :: N_Cr52

      real(8), dimension(:,:), allocatable :: N_V1_Old
      real(8), dimension(:,:), allocatable :: N_V2_Old
      real(8), dimension(:,:), allocatable :: N_Cr52_Old

      real(8), dimension(:,:), allocatable :: I_beta_A2
      real(8), dimension(:,:), allocatable :: I_beta_A2m
      real(8), dimension(:,:), allocatable :: I_beta_A3
      real(8), dimension(:,:), allocatable :: I_beta_B2
      real(8), dimension(:,:), allocatable :: I_beta_B3
      real(8), dimension(:,:), allocatable :: I_beta_rh
      real(8), dimension(:,:), allocatable :: I_beta_rh_p
      real(8), dimension(:,:), allocatable :: I_beta_rh_d

      real(8), dimension(:,:), allocatable :: I_beta_V2
      real(8), dimension(:,:), allocatable :: I_beta_si

      integer :: flag_jumpin_2d = 1
      integer :: flag_jumpin_4n = 0
      real(8), dimension(:,:), allocatable :: inp_jump_bu
      real(8), dimension(:,:), allocatable :: jump_bu
      real(8), dimension(:,:,:), allocatable :: burned_nden


      end module Inc_Nuclide
