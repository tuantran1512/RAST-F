
      MODULE Inc_Geometry
      ! *** Variable ***
      ! 1)  h_x             :: Node Size along X-Axis
      !     (Ixy)
      ! 2)  h_y             :: Node Size along Y-Axis
      !     (Ixy)
      ! 3)  h_z             :: Node Size along Z-Axis
      !     (Iz)
      ! 4)  Mesh_x          :: Number of Grid in each Column Grid(= 1N/1FA)
      ! 5)  Mesh_y          :: Number of Grid in each Row    Grid
      ! 6)  Mesh_z          :: Number of Grid in each Axial  Grid
      ! 7)  GridSize_x      :: Grid Size along X-Axis
      !     (Ix_1N)
      ! 8)  GridSize_y      :: Grid Size along Y-Axis
      !     (Iy_1N)
      ! 9)  GridSize_y      :: Grid Size along Z-Axis
      !     (Iz)
      ! 10) MeshSize_x      :: Mesh Size along X-Axis
      !     (Ix)
      ! 11) MeshSize_y      :: Mesh Size along Y-Axis
      !     (Iy)
      ! 12) MeshSize_z      :: Mesh Size along Z-Axis
      !     (Iz)
      ! 13) Nx              :: Total Number of Column
      ! 14) Ny              :: Total Number of Row
      ! 15) Nx_y            :: Total Number of Column in each Row
      ! 16) Nx_1N           :: Total Number of Column (1N/1FA)
      ! 17) Ny_1N           :: Total Number of Row (1N/1FA)
      ! 18) Nx_y_1N         :: Total Number of Column in each Row (1N/1FA)
      ! 19) Nz              :: Total Number of Axial(1-D, Z) Node
      ! 20) Nxy             :: Total Number of Radial(2-D, XY) Node
      ! 21) Nxy_1N          :: Total Number of Radial(2-D, XY) Node (1N/1FA)
      ! 22) Nxyz            :: Total Number of Entire(3-D, XYZ) Node
      ! 23) Ixy             :: Radial Node Index
      ! 24) Ixy_1N          :: Radial Node Index (1N/1FA)
      ! 25) Ixy_Start_y     :: Start Number Ixy of FA or RF in each Row
      ! 26) Ixy_End_y       :: End   Number Ixy of FA or RF in each Row
      ! 27) Ix_Start_y      :: Start Number Ix of FA or RF in each Row (1N/1FA)
      ! 28) Ix_End_y        :: End   Number Ix of FA or RF in each Row (1N/1FA)
      ! 29) Ix_Start_y_1N   :: Start Number Ix of FA or RF in each Row (1N/1FA)
      ! 30) Ix_End_y_1N     :: End   Number Ix of FA or RF in each Row (1N/1FA)
      ! 31) Ix_StartFA_y_1N :: Start Number Ix of FA in each Row (1N/1FA)
      ! 32) Ix_EndFA_y_1N   :: End   Number Ix of FA in each Row (1N/1FA)
      ! 33) Iy_Start_x      :: Start Number Iy of FA or RF in each Column (1N/1FA)
      ! 34) Iy_End_x        :: End   Number Iy of FA or RF in each Column (1N/1FA)
      ! 35) Iy_Start_x_1N   :: Start Number Iy of FA or RF in each Column (1N/1FA)
      ! 36) Iy_End_x_1N     :: End   Number Iy of FA or RF in each Column (1N/1FA)
      ! 37) Iy_StartFA_x_1N :: Start Number Iy of FA in each Column (1N/1FA)
      ! 38) Iy_EndFA_x_1N   :: End   Number Iy of FA in each Column (1N/1FA)
      ! 39) Ix              :: X-Axis Index
      ! 40) Iy              :: Y-Axis Index
      ! 41) Iz              :: Z-Axis Index
      ! 42) IzFuelBot       :: Z-Axis Index of Fuel Bottom Region
      ! 43) IzFuelTop       :: Z-Axis Index of Fuel Top Region
      ! 44) Ix_1N           :: X-Axis Index (1N/1FA)
      ! 45) Iy_1N           :: Y-Axis Index (1N/1FA)
      ! 46) Ix_4Nto1N       :: Ix => Ix_1N
      ! 47) Iy_4Nto1N       :: Iy => Iy_1N
      ! 48) N_ID            :: Node ID; Type of Node Depend on Radial Exposed Surface (4-Direction)
      !                        It Use the Units(1) Digit. Negative(-) Means Opposite Case for Positive(+).
      !
      !                            X   [Radial Expose Type]
      !                            -->
      !                        Y | 1 | 2 | 3    1: Lx, Ly  2: Ly     3: Rx, Ly
      !                          V --|---|--
      !                            4 | 5 | 6    4: Lx      5: No     6: Rx
      !                            --|---|--
      !                            7 | 8 | 9    7: Lx, Ry  8: Ry     9: Rx, Ry
      !
      !                                       - 2: Lx, Rx, Ry  - 4: Rx, Ly, Ry
      !
      !                                       - 6: Lx, Ly, Ry  - 8: Lx, Rx, Ly
      !
      !                                         456: Lx, Rx        258: Ly, Ry
      !
      !                        Axial Expose Type will be Treated Case by Case
      !                        Bottom(Iz = 1): Lz    Top(Iz = Nz): Rz
      ! 49) I_Rx            :: Right (X-Axis) Node Index
      ! 50) I_Lx            :: Left  (X-Axis) Node Index
      ! 51) I_Ry            :: Right (Y-Axis) Node Index
      ! 52) I_Ly            :: Left  (Y-Axis) Node Index
      ! 53) I_Rz            :: Right (Z-Axis) Node Index
      ! 54) I_Lz            :: Left  (Z-Axis) Node Index
      !                            X   [Adjacent Node Index]
      !                            -->
      !                        Y |    | Ly |         ^  | Rz |
      !                          V ---|----|---    Z |  |----|
      !                            Lx | ND | Rx         | ND |
      !                            ---|----|---         |----|
      !                               | Ry |            | Lz |
      ! 55) I_4Nto1N        :: Ixy    => Ixy_1N (4N/1FA or #N/1FA => 1N/1FA)
      ! 56) I_1Nto4N        :: Ixy_1N => Ixy    (1N/1FA => 4N/1FA or #N/1FA)
      ! 57) Ig              :: Neutron Energy Group Index
      ! 58) Core_Height     :: Active Core Hight (Fuel Region)
      ! 59) Tot_Vol         :: Total Volume of (Fuel + Reflector) Region
      ! 60) Tot_FuelVol     :: Total Volume of Fuel Region
      ! 61) NodeVolume      :: Node-wise Volume
      ! 62) BC_Lx           :: Boundary Condition of Left  Surface along X-Axis
      !                          1 = Net Current Zero (Reflected)
      !                          2 = Incoming Partial Current Zero (Vacuum)
      !                          3 = Flux Zero
      ! 63) BC_Rx           :: Boundary Condition of Right Surface along X-Axis
      ! 64) BC_Ly           :: Boundary Condition of Left  Surface along Y-Axis
      ! 65) BC_Ry           :: Boundary Condition of Right Surface along Y-Axis
      ! 66) BC_Lz           :: Boundary Condition of Left  Surface along Z-Axis
      ! 67) BC_Rz           :: Boundary Condition of Right Surface along Z-Axis
      ! 68) OPT_Core        :: Option of Core Geometry
      !                          1 = Full         Geometry
      !                          4 = Quarter(1/4) Geometry
      ! 69) Flag_4N1FA      :: 4N/1FA Flag
      ! 70) Flag_HalfCenter :: Half Center Line Flag in 4N/1FA + Quarter Core
      ! 71) IxIyToIxy       :: (Ix, Iy) -> (Ixy)
      !     (Ix, Iy)           (None Region) -> 0
      ! 72) IxIyToIxy_1N    :: (Ix, Iy) -> (Ixy_1N)
      !     (Ix, Iy)           (None Region) -> 0
      ! 73) IxIy_1NToIxy    :: (Ix_1N, Iy_1N) -> (Ixy)
      !     (Ix_1N, Iy_1N)     (None Region) -> 0
      ! 74) IxIy_1NToIxy_1N :: (Ix_1N, Iy_1N) -> (Ixy_1N)
      !     (Ix_1N, Iy_1N)     (None Region) -> 0
      ! 75) IxToIx_1N       :: (Ix) -> (Ix_1N)
      !     (Ix)
      ! 76) IyToIy_1N       :: (Iy) -> (Iy_1N)
      !     (Iy)
      ! 77) Core_N_FARF     :: Nxy_1N    in Full Core
      ! 78) Core_N_FA       :: Nxy_FA_1N in Full Core

      IMPLICIT NONE

      INTEGER, DIMENSION(:), ALLOCATABLE :: ixytoix
      INTEGER, DIMENSION(:), ALLOCATABLE :: ixytoiy

      REAL(8), DIMENSION(:), ALLOCATABLE :: h_x
      REAL(8), DIMENSION(:), ALLOCATABLE :: h_y
      REAL(8), DIMENSION(:), ALLOCATABLE :: h_z
      REAL(8), DIMENSION(:), ALLOCATABLE :: GridSize_x
      REAL(8), DIMENSION(:), ALLOCATABLE :: GridSize_y
      REAL(8), DIMENSION(:), ALLOCATABLE :: GridSize_z
      REAL(8), DIMENSION(:), ALLOCATABLE :: MeshSize_x
      REAL(8), DIMENSION(:), ALLOCATABLE :: MeshSize_y
      REAL(8), DIMENSION(:), ALLOCATABLE :: MeshSize_z
      INTEGER :: Nx_1N
      INTEGER :: Ny_1N
      INTEGER :: Nz
      INTEGER, DIMENSION(:), ALLOCATABLE :: Nx_y_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Mesh_x
      INTEGER, DIMENSION(:), ALLOCATABLE :: Mesh_y
      INTEGER, DIMENSION(:), ALLOCATABLE :: Mesh_z
      INTEGER :: BC_Rx
      INTEGER :: BC_Lx
      INTEGER :: BC_Ry
      INTEGER :: BC_Ly
      INTEGER :: BC_Rz
      INTEGER :: BC_Lz
      INTEGER :: OPT_Core
      REAL(8) :: Core_Height
      REAL(8) :: Tot_Vol
      REAL(8) :: Tot_FuelVol
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: NodeVolume
      INTEGER :: Nx
      INTEGER :: Ny
      INTEGER :: Nxy
      INTEGER :: Nxyz
      INTEGER :: Nxy_1N
      INTEGER :: Core_N_FARF
      INTEGER :: Core_N_FA
      INTEGER :: IzFuelBot
      INTEGER :: IzFuelTop
      INTEGER, DIMENSION(:), ALLOCATABLE :: Nx_y
      INTEGER, DIMENSION(:), ALLOCATABLE :: N_ID
      INTEGER, DIMENSION(:), ALLOCATABLE :: Ixy_Start_y
      INTEGER, DIMENSION(:), ALLOCATABLE :: Ixy_End_y
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_Lx
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_Rx
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_Ly
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_Ry
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_Rz
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_Lz
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_4Nto1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Ix_4Nto1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Iy_4Nto1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Ix_Start_y_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Ix_End_y_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Ix_StartFA_y_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Ix_EndFA_y_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Iy_Start_x_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Iy_End_x_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Iy_StartFA_x_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Iy_EndFA_x_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Ix_Start_y
      INTEGER, DIMENSION(:), ALLOCATABLE :: Ix_End_y
      INTEGER, DIMENSION(:), ALLOCATABLE :: Iy_Start_x
      INTEGER, DIMENSION(:), ALLOCATABLE :: Iy_End_x
      INTEGER, DIMENSION(:), ALLOCATABLE :: IxToIx_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: IyToIy_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_4Nin1FA

      INTEGER, DIMENSION(:, :), ALLOCATABLE :: I_1Nto4N
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: IxIyToIxy
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: IxIyToIxy_1N
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: IxIy_1NToIxy
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: IxIy_1NToIxy_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: ixy_1ntoix_1n
      INTEGER, DIMENSION(:), ALLOCATABLE :: ixy_1ntoiy_1n

      LOGICAL(1) :: Flag_4N1FA
      LOGICAL(1) :: Flag_HalfCenter

      character, dimension(:), allocatable :: x_fa_index
      integer, dimension(:), allocatable :: y_fa_index
      character(5), dimension(:), allocatable :: xy_fa_index
      character(5), dimension(:), allocatable :: xy_fa_index_org
      character, dimension(:,:), allocatable :: x_fa_findex
      integer, dimension(:,:), allocatable :: y_fa_findex
      character(5), dimension(:,:), allocatable :: xy_fa_findex
      character(5), dimension(:,:), allocatable :: xy_fa_findex_org
      integer :: fa_istart

      INTEGER:: Nxy_FA_P2R
      INTEGER,ALLOCATABLE :: I_Lx_P2R(:)
      INTEGER,ALLOCATABLE :: I_Rx_P2R(:)
      INTEGER,ALLOCATABLE :: I_Ly_P2R(:)
      INTEGER,ALLOCATABLE :: I_Ry_P2R(:)
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: nodelfa
      INTEGER:: idom
      INTEGER:: ndomx
      INTEGER:: ndomy
      INTEGER:: ndomz
      INTEGER:: ndom
      INTEGER:: idomx
      INTEGER:: idomy
      INTEGER:: idomz
      INTEGER:: idomxy
      INTEGER:: npfa
      INTEGER, ALLOCATABLE :: nrnx(:)
      INTEGER, ALLOCATABLE :: ltola(:)
      INTEGER, ALLOCATABLE :: latol(:)
      INTEGER, ALLOCATABLE :: ltox(:)
      INTEGER, ALLOCATABLE :: ltoy(:)
      INTEGER, ALLOCATABLE :: ltolf(:)
      INTEGER, ALLOCATABLE :: ltolfa(:)
      INTEGER, ALLOCATABLE :: lfatol(:)
      INTEGER, ALLOCATABLE :: lfaptr(:)
      INTEGER, ALLOCATABLE :: nodel(:,:)
      INTEGER, ALLOCATABLE :: nodef(:)
      REAL(8), ALLOCATABLE :: znode(:)
      REAL(8), ALLOCATABLE :: zcent(:)
      REAL(8), ALLOCATABLE :: zb(:)
      REAL(8):: hyz
      REAL(8):: hactive
      real(8):: hbot_active ! bot position of active fuel region [cm]
      real(8):: htop_active ! top position of active fuel region [cm]

      logical(1) :: flag_axial_mirror=.true.

      logical(1) :: flag_auto_geom = .false.
      integer :: opt_rx_type = 0
      real(8) :: coef_the_z = 1d0
      integer, dimension( 9, 9) :: auto_farf_1_q
      integer, dimension(10,10) :: auto_farf_2_q
      integer, dimension( 9, 9) :: auto_farf_3_q
      integer, dimension( 8, 8) :: auto_farf_4_q
      integer, dimension(17,17) :: auto_farf_1_f
      integer, dimension(19,19) :: auto_farf_2_f
      integer, dimension(17,17) :: auto_farf_3_f
      integer, dimension(15,15) :: auto_farf_4_f
      data auto_farf_1_q(1:9,1:9)                          &
      / 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        2, 2, 2, 2, 2, 2, 2, 2, 1, &
        2, 2, 2, 2, 2, 2, 2, 2, 1, &
        2, 2, 2, 2, 2, 2, 2, 1, 1, &
        2, 2, 2, 2, 2, 2, 2, 1, 0, &
        2, 2, 2, 2, 2, 2, 1, 1, 0, &
        2, 2, 2, 2, 2, 1, 1, 0, 0, &
        2, 2, 2, 1, 1, 1, 0, 0, 0, &
        1, 1, 1, 1, 0, 0, 0, 0, 0  /
      data auto_farf_2_q(1:10,1:10)                        &
      / 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        2, 2, 2, 2, 2, 2, 2, 2, 1, 1, &
        2, 2, 2, 2, 2, 2, 2, 2, 1, 0, &
        2, 2, 2, 2, 2, 2, 2, 1, 1, 0, &
        2, 2, 2, 2, 2, 2, 1, 1, 0, 0, &
        2, 2, 2, 2, 1, 1, 1, 0, 0, 0, &
        1, 1, 1, 1, 1, 0, 0, 0, 0, 0  /
      data auto_farf_3_q(1:9,1:9)                          &
      / 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        2, 2, 2, 2, 2, 2, 2, 2, 1, &
        2, 2, 2, 2, 2, 2, 2, 1, 1, &
        2, 2, 2, 2, 2, 2, 2, 1, 0, &
        2, 2, 2, 2, 2, 2, 1, 1, 0, &
        2, 2, 2, 2, 2, 1, 1, 0, 0, &
        2, 2, 2, 2, 1, 1, 0, 0, 0, &
        2, 2, 1, 1, 1, 0, 0, 0, 0, &
        1, 1, 1, 0, 0, 0, 0, 0, 0  /
      data auto_farf_4_q(1:8,1:8)                          &
      / 2, 2, 2, 2, 2, 2, 2, 1, &
        2, 2, 2, 2, 2, 2, 2, 1, &
        2, 2, 2, 2, 2, 2, 1, 1, &
        2, 2, 2, 2, 2, 2, 1, 0, &
        2, 2, 2, 2, 2, 1, 1, 0, &
        2, 2, 2, 2, 1, 1, 0, 0, &
        2, 2, 1, 1, 1, 0, 0, 0, &
        1, 1, 1, 0, 0, 0, 0, 0  /
      data auto_farf_1_f(1:17,1:17)                        &
      / 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &
        0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, &
        0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0, &
        0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, &
        0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, &
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, &
        0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, &
        0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, &
        0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0, &
        0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, &
        0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0  /
      data auto_farf_2_f(1:19,1:19)                        &
      / 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &
        0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, &
        0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0, &
        0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, &
        0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, &
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, &
        0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, &
        0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, &
        0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0, &
        0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0, &
        0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0  /
      data auto_farf_3_f(1:17,1:17)                        &
      / 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0, &
        0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0, 0, &
        0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0, &
        0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, &
        0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, &
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, &
        0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, &
        0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, &
        0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0, &
        0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0, 0, &
        0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0  /
      data auto_farf_4_f(1:15,1:15)                        &
      / 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &
        0, 0, 0, 1, 1, 1, 2, 2, 2, 1, 1, 1, 0, 0, 0, &
        0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0, &
        0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, &
        0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, &
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, &
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, &
        0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, &
        0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, &
        0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0, &
        0, 0, 0, 1, 1, 1, 2, 2, 2, 1, 1, 1, 0, 0, 0, &
        0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0  /

      END MODULE Inc_Geometry
