
      MODULE Inc_MATFB

      IMPLICIT NONE

      ! Fuel thermal conductivity as burnup[1]
      real(8), parameter :: FKB_A  = 0.0452_8
      real(8), parameter :: FKB_sA = 1.1599_8
      real(8), parameter :: FKB_B  = 2.46e-4_8
      real(8), parameter :: FKB_E  = 3.5e+9_8
      real(8), parameter :: FKB_F  = 16361._8
      real(8), parameter :: FKB_Q  = 6380._8
      REAL(8) :: FKB_h

      ! Oxide generation
      real(8), parameter :: OX_c1 = 1.9599_8
      real(8), parameter :: OX_c2 = 2.41e-4_8
      real(8), parameter :: OX_c3 = 6.43e-7_8
      real(8), parameter :: OX_c4 = 1.946e-10_8

      ! Gap conductance fitting coefficient
      real(8) :: GK_a1 = -7e-5_8
      real(8) :: GK_b1 = 0.0086_8
      real(8) :: GK_c1 = -0.4174_8
      real(8) :: GK_d1 = 9.9863_8
      real(8) :: GK_e1 = -120.61_8
      real(8) :: GK_f1 = 817.48_8
      real(8) :: GK_g1 = 6715.3_8
      real(8) :: GK_a2 = 4e-5_8
      real(8) :: GK_b2 = -0.013_8
      real(8) :: GK_c2 = 1.7706_8
      real(8) :: GK_d2 = -125.57_8
      real(8) :: GK_e2 = 4918.5_8
      real(8) :: GK_f2 = -101402._8
      real(8) :: GK_g2 = 877660._8
      REAL(8) :: GK_BU
      REAL(8), parameter :: GK_fit_BU_cut = 36.05_8

      END MODULE Inc_MATFB
