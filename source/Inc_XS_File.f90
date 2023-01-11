
      MODULE Inc_XS_File

      IMPLICIT NONE

      integer :: N_BU_Ref_Max
      integer :: N_Brch

      REAL(8), DIMENSION(:), ALLOCATABLE :: BU_Ref
      REAL(8), DIMENSION(:), ALLOCATABLE :: BU_Brch
      REAL(8), DIMENSION(:), ALLOCATABLE :: Var_PPM
      REAL(8), DIMENSION(:), ALLOCATABLE :: Var_TF
      REAL(8), DIMENSION(:), ALLOCATABLE :: Var_TM
      REAL(8), DIMENSION(:), ALLOCATABLE :: Var_DM
      REAL(8), DIMENSION(:), ALLOCATABLE :: Var_ref_TF
      REAL(8), DIMENSION(:), ALLOCATABLE :: Var_ref_TM
      REAL(8), DIMENSION(:), ALLOCATABLE :: Var_ref_DM

      REAL(8), DIMENSION(:, :), ALLOCATABLE :: BU_Ref_NoBP
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: BU_Brch_NoBP
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: BU_Brch_wBP
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: BU_Ref_wBP
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Ini_NumDen

      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: W_SP
      INTEGER, DIMENSION(:), ALLOCATABLE :: Type_Tab
      INTEGER, DIMENSION(:), ALLOCATABLE :: N_BU_Ref_wBP
      INTEGER :: I_Tab
      INTEGER :: I_SP
      INTEGER :: I_BU
      INTEGER :: I_Type
      INTEGER :: I_Reg

      logical(1) :: flag_1nhff = .false.
      logical(1) :: flag_1nadf = .false.

      logical(1) :: flag_asymfa = .false.
      integer, dimension(:,:), allocatable :: asym_tab
      integer, dimension(:,:), allocatable :: asym2itab
      integer, dimension(:), allocatable :: asym_rot

      INTEGER :: N_XS_Table=0
      INTEGER :: N_SP
      INTEGER :: XS_File_NbyN
      INTEGER :: N_BU_Ref
      INTEGER :: N_BU_Brch
      INTEGER :: N_BU_Ref_NoBP
      INTEGER :: N_BU_Brch_NoBP
      INTEGER :: N_BU_Brch_wBP
      INTEGER :: N_BU_SDC
      INTEGER :: Buff_N_BU
      INTEGER :: N_SP_FA
      INTEGER :: N_SP_RF

      CHARACTER(5), DIMENSION(:), ALLOCATABLE :: Type_RegBrch
      LOGICAL(1), DIMENSION(:), ALLOCATABLE :: Flag_Reg
      REAL(8), DIMENSION(:), ALLOCATABLE :: XSset
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: XSset_Table
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: XSset_Table_Base
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: XSset_Node
      INTEGER :: N_Region, N_XS_Kind
      INTEGER :: ReflectorDF
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: HFFset
      REAL(8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: HFF_Table
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: HFF_Table_Base
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Ref_maXS
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: dXS_dPPM
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: dXS_dTF
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: dXS_dTM
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: dXS_dDM
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: dXS_dCR
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Ref_HFF
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: dHFF_dPPM
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: dHFF_dTF
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: dHFF_dTM
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: dHFF_dDM
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: dHFF_dCR
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: ADF_Table
      LOGICAL(1) :: Flag_ADF
      LOGICAL(1) :: Flag_RefDF

#ifdef jr_vver
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Ref_maXS_Hex
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: ref_smat_Hex
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: ref_smat
#endif

#ifdef tuan_fr
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: XSset_Table_Hex
      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: XSset_Hex
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: StateID
      INTEGER ::  N_XS_Kind_Hex
      REAL(8) :: FUELTEMP_INT
      REAL(8) :: COOLDENS_INT
      INTEGER  :: N_XS_Table_tmp

#endif


#ifdef siarhei_tr_hex
      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: XSset_Hex_CR ! for recording all rodded cases
#endif

#ifdef tuan_tr_test
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Ref_chid

      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: dxs_hex      ! index, icomp, ig
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: dsmat_hex ! index, icomp, ig, igg 
      !! index 1 = PPM, 2 = TMO, 3 = DM, 4 =TFU 
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: dxs_dCR_hex
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: dsmat_dCR   !! for condense of multi-group XS
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: dsmat_dCR_hex
#endif

      type macro_xs_table_
         real(8), allocatable :: xs(:,:)
         real(8), allocatable :: hff(:,:,:)
      end type macro_xs_table_
      type(macro_xs_table_), allocatable :: macro_xs_table(:,:)
      real(8), allocatable :: MacroXSset(:)

      LOGICAL(1) :: Flag_Leakage
      LOGICAL(1) :: Flag_ReflDF
      LOGICAL(1) :: Flag_LC_mean
      LOGICAL(1) :: Flag_LC_micro
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: Leakage_Table
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: leak_ratio

      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: macxs_b0
      REAL(8), DIMENSION(:), ALLOCATABLE :: W_BU_Ref
      REAL(8), DIMENSION(:), ALLOCATABLE :: W_BU_Brch
      REAL(8), DIMENSION(:), ALLOCATABLE :: W_SP_1
      REAL(8), DIMENSION(:), ALLOCATABLE :: W_PPM
      REAL(8), DIMENSION(:), ALLOCATABLE :: W_TF
      REAL(8), DIMENSION(:), ALLOCATABLE :: W_TM
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_BU_Ref
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_BU_Brch
      REAL(8), DIMENSION(:), ALLOCATABLE :: WC_TF
      REAL(8), DIMENSION(:), ALLOCATABLE :: WC_TM
      REAL(8), DIMENSION(:), ALLOCATABLE :: WC_DM

      INTEGER(4) :: if_th     = 4
      INTEGER(4) :: if_branch = 2
      INTEGER(4) :: if_state  = 1

      END MODULE Inc_XS_File

