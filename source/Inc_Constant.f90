
      MODULE Inc_Constant

      IMPLICIT NONE

      REAL(8), PARAMETER :: D0   =  0.0D0
      REAL(8), PARAMETER :: D1   =  1.0D0
      REAL(8), PARAMETER :: D2   =  2.0D0
      REAL(8), PARAMETER :: D3   =  3.0D0
      REAL(8), PARAMETER :: D4   =  4.0D0
      REAL(8), PARAMETER :: D5   =  5.0D0
      REAL(8), PARAMETER :: D6   =  6.0D0
      REAL(8), PARAMETER :: D7   =  7.0D0
      REAL(8), PARAMETER :: D8   =  8.0D0
      REAL(8), PARAMETER :: D9   =  9.0D0
      REAL(8), PARAMETER :: D10  = 10.0D0
      REAL(8), PARAMETER :: D11  = 11.0D0
      REAL(8), PARAMETER :: D12  = 12.0D0
      REAL(8), PARAMETER :: D13  = 13.0D0
      REAL(8), PARAMETER :: D14  = 14.0D0
      REAL(8), PARAMETER :: D15  = 15.0D0
      REAL(8), PARAMETER :: D16  = 16.0D0
      REAL(8), PARAMETER :: D17  = 17.0D0
      REAL(8), PARAMETER :: D18  = 18.0D0
      REAL(8), PARAMETER :: D19  = 19.0D0
      REAL(8), PARAMETER :: D20  = 20.0D0
      REAL(8), PARAMETER :: D21  = 21.0D0
      REAL(8), PARAMETER :: D22  = 22.0D0
      REAL(8), PARAMETER :: D23  = 23.0D0
      REAL(8), PARAMETER :: D24  = 24.0D0
      REAL(8), PARAMETER :: DP1  = 1.0D+01
      REAL(8), PARAMETER :: DP2  = 1.0D+02
      REAL(8), PARAMETER :: DP3  = 1.0D+03
      REAL(8), PARAMETER :: DP4  = 1.0D+04
      REAL(8), PARAMETER :: DP5  = 1.0D+05
      REAL(8), PARAMETER :: DP6  = 1.0D+06
      REAL(8), PARAMETER :: DP7  = 1.0D+07
      REAL(8), PARAMETER :: DP8  = 1.0D+08
      REAL(8), PARAMETER :: DP9  = 1.0D+09
      REAL(8), PARAMETER :: DP10 = 1.0D+10
      REAL(8), PARAMETER :: DP11 = 1.0D+11
      REAL(8), PARAMETER :: DP12 = 1.0D+12
      REAL(8), PARAMETER :: DP13 = 1.0D+13
      REAL(8), PARAMETER :: DP14 = 1.0D+14
      REAL(8), PARAMETER :: DP15 = 1.0D+15
      REAL(8), PARAMETER :: DP16 = 1.0D+16
      REAL(8), PARAMETER :: DP17 = 1.0D+17
      REAL(8), PARAMETER :: DP18 = 1.0D+18
      REAL(8), PARAMETER :: DP19 = 1.0D+19
      REAL(8), PARAMETER :: DP20 = 1.0D+20
      REAL(8), PARAMETER :: DP21 = 1.0D+21
      REAL(8), PARAMETER :: DP22 = 1.0D+22
      REAL(8), PARAMETER :: DP23 = 1.0D+23
      REAL(8), PARAMETER :: DP24 = 1.0D+24
      REAL(8), PARAMETER :: DM1  = 1.0D-01
      REAL(8), PARAMETER :: DM2  = 1.0D-02
      REAL(8), PARAMETER :: DM3  = 1.0D-03
      REAL(8), PARAMETER :: DM4  = 1.0D-04
      REAL(8), PARAMETER :: DM5  = 1.0D-05
      REAL(8), PARAMETER :: DM6  = 1.0D-06
      REAL(8), PARAMETER :: DM7  = 1.0D-07
      REAL(8), PARAMETER :: DM8  = 1.0D-08
      REAL(8), PARAMETER :: DM9  = 1.0D-09
      REAL(8), PARAMETER :: DM10 = 1.0D-10
      REAL(8), PARAMETER :: DM11 = 1.0D-11
      REAL(8), PARAMETER :: DM12 = 1.0D-12
      REAL(8), PARAMETER :: DM13 = 1.0D-13
      REAL(8), PARAMETER :: DM14 = 1.0D-14
      REAL(8), PARAMETER :: DM15 = 1.0D-15
      REAL(8), PARAMETER :: DM16 = 1.0D-16
      REAL(8), PARAMETER :: DM17 = 1.0D-17
      REAL(8), PARAMETER :: DM18 = 1.0D-18
      REAL(8), PARAMETER :: DM19 = 1.0D-19
      REAL(8), PARAMETER :: DM20 = 1.0D-20
      REAL(8), PARAMETER :: DM21 = 1.0D-21
      REAL(8), PARAMETER :: DM22 = 1.0D-22
      REAL(8), PARAMETER :: DM23 = 1.0D-23
      REAL(8), PARAMETER :: DM24 = 1.0D-24

      REAL(8), PARAMETER :: Half   = 0.5D0
      REAL(8), PARAMETER :: PI     = 3.141592653589793D0
      REAL(8), PARAMETER :: inch   = 2.54D0
      REAL(8), PARAMETER :: DegToK = 273.15D0
      REAL(8), PARAMETER :: barn   = 1.0D-24
      REAL(8), PARAMETER :: NA     = 6.02252D23
      real(8), parameter :: n_mass = 1.00866491588d0 ! amu
      real(8), parameter :: n_mass_mev = 939.565d0 ! MeV/c^2
      real(8), parameter :: v_light = 299792458d0 ! m/s

      INTEGER :: NBD
      INTEGER :: NBI
      INTEGER :: NBS
      INTEGER :: NBF
      INTEGER :: ng2
      INTEGER :: ndecgrp
      LOGICAL:: TRUE
      LOGICAL:: FALSE
      REAL(8) :: CKELVIN
      REAL(8) :: PI_P2R
      REAL(8) :: cinf
      REAL(8) :: zero
      REAL(8) :: one
      REAL(8) :: half_P2R
      REAL(8) :: big
      PARAMETER (NBD=8,NBI=4,NBS=4)
      PARAMETER (NBF=NBD)
      PARAMETER (ng2=2,ndecgrp=6)
      PARAMETER (TRUE=.true.,FALSE=.false.)
      PARAMETER (CKELVIN=273.15,PI_P2R=3.1415927)
      PARAMETER (cinf=1.0e30,zero=0,one=1,half_P2R=0.5,big=cinf)
      INTEGER :: ndegrp

      END MODULE Inc_Constant
