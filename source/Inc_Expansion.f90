#ifdef tuan_fr_TherEx
    MODULE Inc_Expansion
        IMPLICIT NONE

        LOGICAL(1) :: Flag_Card_Thermal_Expansion = .false.
        INTEGER :: N_FuelRegion                     ! Number of fuel regions
!        INTEGER(8), dimension(:), allocatable :: Iz_region
        REAL(8), dimension(:), allocatable :: Iz_region
        REAL(8), dimension(:), allocatable :: Fuel_length
        REAL(8), dimension(:), allocatable :: Fuel_length_Exp
        real(8), allocatable ,dimension(:) :: F_ID(:), Exp_Coef(:), Initial_T_F(:), Expect_T_F(:), ABC(:)
        INTEGER (8) :: N_FuelType
        REAL(8), dimension(:), allocatable :: exp_value
!       radial expansion
!        REAL(8) :: R1
!        REAL(8) :: R2
        INTEGER (8) :: N_CooMAT
        INTEGER (8), dimension(:), allocatable :: CooMAT_ID
        INTEGER (8), dimension(:,:,:), allocatable :: N_FR_Exp
        

!
        LOGICAL(1) :: miXS_exp = .false.
        LOGICAL(1) :: Axial_exp = .false.
        LOGICAL(1) :: Radial_exp = .false.
        LOGICAL(1) :: maXS_exp = .false.
        
        
    END MODULE 

#endif
