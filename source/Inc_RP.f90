
      MODULE Inc_RP
      ! *** Variable ***
      ! 1)  Name_LP      :: Loading Pattern Type Name
      ! 2)  N_LP         :: Total Number of Loading Pattern Type
      ! 3)  AxialComp    :: Composition Number (= Material Number or XS Table)
      !     (I_LP_1N, Iz)   along Axial Node of a Loading Pattern Type
      ! 4)  I_Comp       :: Composition (= Material or XS Table) Index
      !     (Ixy, Iz)
      ! 5)  LP_MTU       :: Loaded Uranium in a Loading Pattern Type (Metric Ton Uranium)  [ton]
      ! 6)  Flag_BP      :: BP FA Flag
      !     (I_Tab)
      ! 7)  I_LP_1N      :: Loading Pattern Type Index (1N/1FA)
      ! 8)  Nxy_FA       :: Total Number of Fuel Assembly Node
      ! 9)  Nxy_RF       :: Total Number of Reflector Node
      ! 10) Nxy_FA_1N    :: Total Number of Fuel Assembly Node (1N/1FA)
      ! 11) Nxy_RF_1N    :: Total Number of Reflector Node (1N/1FA)
      ! 12) I_FA         :: Fuel Assembly Node Index
      !     (Nxy_FA)          (Ixy_FA -> Ixy)
      ! 13) I_RF         :: Reflector Node Index
      !     (Nxy_RF)          (Ixy_RF -> Ixy)
      ! 14) I_FA_1N      :: Fuel Assembly Node Index (1N/1FA)
      !                       (Ixy_FA_1N -> Ixy_1N)
      ! 15) I_RF_1N      :: Reflector Node Index (1N/1FA)
      !                       (Ixy_RF_1N -> Ixy_1N)
      ! 16) I_FARF_1N    :: (Fuel Assembly/Reflector) Node Index (1N/1FA)
      !     (Ixy_1N)          0 = None, 1 = Reflector, 2 >= Fuel Assembly
      ! 17) I_FARF_1N_2D :: (Fuel Assembly/Reflector) Node Index (1N/1FA)
      !     (Ix_1N, Iy_1N)
      !                       0 = None
      !                       1 = Reflector
      !                       2 = Fuel
      ! 18) I_FARF_1N_3D :: I_FARF_1N [3D]
      !     (Ixy_1N, Iz)    If a Node [3D] has Zero nu_maXS_f,
      !                     then this Node is Reflector Node
      ! 19) I_FARF_2D    :: (Fuel Assembly/Reflector) Node Index
      !     (Ix, Iy)

      IMPLICIT NONE

      real(8), dimension(:), allocatable :: LP_MTU
      integer :: N_LP
      integer :: Nxy_FA
      integer :: Nxy_RF
      integer :: Nxy_FA_1N
      integer :: Nxy_RF_1N

      integer, dimension(:), allocatable :: I_LP_1N
      integer, dimension(:), allocatable :: I_FA
      integer, dimension(:), allocatable :: I_RF
      integer, dimension(:), allocatable :: I_FA_1N
      integer, dimension(:), allocatable :: I_RF_1N
      integer, dimension(:), allocatable :: I_FARF_1N

      real(8), dimension(:,:), allocatable :: I_Gd_FA

      integer, dimension(:,:), allocatable :: AxialComp
      integer, dimension(:,:), allocatable :: I_Comp
      integer, dimension(:,:), allocatable :: I_FARF_1N_2D
      integer, dimension(:,:), allocatable :: I_FARF_1N_3D
      integer, dimension(:,:), allocatable :: I_FARF_2D

      integer, dimension(:), allocatable :: index_LP
      character(10), dimension(:), allocatable :: Name_LP
      logical(1), dimension(:), allocatable :: Flag_BP

      integer, dimension(:), allocatable :: ixytoifa
#ifdef tuan_fr
      integer, dimension(:), allocatable :: I_LP_1N_tmp
      integer, dimension(20,100) :: axialcomp_tmp
#endif
      END MODULE Inc_RP
