
      MODULE Inc_CR

      IMPLICIT NONE

      REAL(8) :: Step_CR, Speed_CR, OverLap

      real(8), dimension(:), allocatable :: Length_CR
      real(8), dimension(:), allocatable :: Length_CR_Tip
      real(8), dimension(:), allocatable :: CR_Bot
      real(8), dimension(:), allocatable :: CR_Top
      real(8), dimension(:), allocatable :: dCR_MoveCR
      real(8), dimension(:), allocatable :: v_MoveCR
      real(8), dimension(:), allocatable :: PDIL
      real(8), dimension(:), allocatable :: PDWL
      real(8), dimension(:), allocatable :: Reg_VF
      real(8), dimension(:), allocatable :: MoveCR_Bot

      REAL(8), DIMENSION(:, :), ALLOCATABLE :: T_MoveCR

      INTEGER :: N_CR
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_CR_1N
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_CR_4N
      INTEGER, DIMENSION(:), ALLOCATABLE :: Comp_CR
#ifdef tuan_tr_test
      integer(4)                        :: N_MoveCR
      integer(4),         allocatable   :: I_MoveCR                   (:)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: N_MoveCR
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: I_MoveCR
#endif
      CHARACTER(5), DIMENSION(:), ALLOCATABLE :: Name_CR
      LOGICAL(1), DIMENSION(:), ALLOCATABLE :: Flag_Move
      LOGICAL(1) :: Flag_MoveCR

      INTEGER ::  rodmove = 1

      real(8), dimension(:,:), allocatable :: CR_mat_frac

#ifdef tuan_fr_crm
      real(8)   :: CR_Part
      integer(8)   :: Rep_MAT
      real(8) :: NoMove
      real(8), dimension(:), allocatable :: Abs_Bot
      real(8), dimension(:), allocatable :: BotFol_Bot
      real(8), dimension(:), allocatable :: TopFol_Bot
      real(8), dimension(:), allocatable :: Abs_len
      real(8) :: botfol_len
      integer(8), dimension(:), allocatable :: BotFol_Mat
      integer(8), dimension(:), allocatable :: Abs_Mat
      integer(8), dimension(:), allocatable :: TopFol_Mat
      integer(8), dimension(:), allocatable :: Ix_CR
      integer(8), dimension(:), allocatable :: Iy_CR
      integer(8), dimension(:,:), allocatable :: Ix_surFA
      integer(8), dimension(:,:), allocatable :: Iy_surFA
      integer(8), dimension(:), allocatable :: Nxy_CR
      integer(8), dimension(:,:), allocatable :: Nxy_surFA
      real(8), dimension(:,:,:), allocatable :: SPH_F
      real(8), dimension(:,:,:), allocatable :: SPH_CR
     LOGICAL(1) :: flag_decusping = .false.
     LOGICAL(1) :: flag_decusping_update = .false.
     LOGICAL(1) :: flag_decusping_update2 = .false.
     LOGICAL(1) :: flag_SPH_F = .false.
     LOGICAL(1) :: flag_SPH_CR = .false.
     integer(8)   :: N_CtrolType = 0
     integer(8)   :: N_FuelType = 0
      
#endif
#ifdef siarhei_tr_hex
      ! integer, dimension(:,:), allocatable :: AxialComp ! this variable from Inc_RP is what we affect in CR criticality search
      logical                              :: crit_search_CR = .false.
      logical                              :: converged_crit_CR = .false.
      real(8)                              :: keff_cr_search_old ! used for saving keff for the last CR position for further interpolation criticality search
      real(8), dimension(:),   allocatable :: CR_Position        ! copy CR_Bot initially and then update using interpolation for criticality search
      real(8), dimension(:),   allocatable :: CR_Position_old    ! saves previous CR positions for interpolated CR criticality search
      character(10),           allocatable :: Name_CR_FA(:)
      integer, dimension(:,:), allocatable :: CR_Comp_Ins      ! record fully inserted composition state
      integer, dimension(:,:), allocatable :: CR_Comp_Wtr      ! record fully withdrawn composition state
      integer, dimension(:,:), allocatable :: CR_Comp_HEX      ! record current CR state (criticality search)
      integer, dimension(:,:), allocatable :: CR_Comp_HEX_old  ! record previous step CR state
      real(8), dimension(:),   allocatable :: cumulat_height_z ! copies and sums values from GridSize_z for determining CR position
      integer, dimension(:),   allocatable :: Name_LP_to_CR    ! reverse conversion from LP to CR equivalents
      integer                              :: bottom_pos_cr = 0
      real(8), dimension(:), allocatable   :: CR_GridSize_z    ! same as GridSize_z but with 0th element = 0.0
      ! real(8), dimension(:), allocatable   :: position_of_CR   
      real(8), dimension(:,:), allocatable :: history_of_CR_search ! used to save all iteration steps and then print as a table

#endif

      END MODULE Inc_CR
