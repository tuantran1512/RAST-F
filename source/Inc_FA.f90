
      MODULE Inc_FA

      IMPLICIT NONE

      REAL(8) :: h_FA
      REAL(8) :: N_Pin
      REAL(8) :: N_GT
      REAL(8) :: N_GT_NoCR

      REAL(8) :: R_Fuel
      REAL(8) :: R_Gap
      REAL(8) :: R_Clad
      REAL(8) :: h_Clad
      REAL(8) :: R_GT_CR
      REAL(8) :: R_GT_IClad
      REAL(8) :: R_GT_OClad

      REAL(8) :: VF_H2O_NoCR
      REAL(8) :: VF_H2O_InCR

      logical(1) :: flag_auto_fa = .false.

      END MODULE Inc_FA
