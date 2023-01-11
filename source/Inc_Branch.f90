module inc_branch

   logical(1) :: flag_tmfb_back
   logical(1) :: flag_tffb_back
   logical(1) :: flag_xsfb_back
   logical(1) :: flag_thfb_back

   integer(4) :: opt_mode_back
   integer(4) :: opt_xe_back
   integer(4) :: opt_sm_back
   integer(4) :: opt_gd_back
   logical(4) :: OPT_findtavg_back

   real(8)    :: ppm_back
   real(8)    :: ppower_back
   real(8)    :: core_massflow_back
   real(8)    :: tm_in_back
   real(8)    :: keff_back
   real(8)    :: tf_in_back

   real(8), allocatable :: cr_bot_back(:)

   REAL(8), DIMENSION(:, :), ALLOCATABLE :: BU_back

   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_U34_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_U35_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_U36_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_U37_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_U38_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Np37_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Np38_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Np39_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu38_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu39_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu40_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu41_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu42_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu43_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Am41_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_As42_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Am42_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Am43_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Am44_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Cm42_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Cm43_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Cm44_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_I35_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Xe35_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Nd47_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Nd48_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Nd49_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pm47_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Ps48_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pm48_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pm49_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Sm47_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Sm48_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Sm49_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd52_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd54_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd55_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd56_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd57_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd58_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd60_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_B0_back

   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_U34_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_U35_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_U36_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_U37_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_U38_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Np37_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Np38_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Np39_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu38_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu39_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu40_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu41_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu42_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pu43_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Am41_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_As42_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Am42_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Am43_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Am44_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Cm42_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Cm43_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Cm44_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_I35_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Xe35_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Nd47_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Nd48_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Nd49_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pm47_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Ps48_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pm48_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Pm49_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Sm47_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Sm48_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Sm49_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd52_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd54_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd55_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd56_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd57_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd58_old_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: N_Gd60_old_back

   REAL(8), DIMENSION(:, :), ALLOCATABLE :: T_mod_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: D_mod_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: T_fuel_back

   REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: flux_back
   REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: flux_old_back

   REAL(8), DIMENSION(:, :), ALLOCATABLE :: FisSrc_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: Power_back
   REAL(8), DIMENSION(:, :), ALLOCATABLE :: Normal_power_back
   logical(1)                            :: if_mtc_bu = .false.
   logical(1)                            :: if_ndr_logic = .false.
   logical(1)                            :: if_mtc_ppm       =.false.
   logical(1)                            :: if_ftc_pow       =.false.
   logical(1)                            :: if_dpc_pow       =.false.
   logical(1)                            :: if_dpd_pow       =.false.
   logical(1)                            :: if_sdm_pow       =.false.

   type :: level1
      character(30) :: case_name=''
      logical(1),allocatable :: if_step(:)
      real(8),allocatable :: keff(:)
      real(8),allocatable :: ndr_bu(:)
      real(8),allocatable :: TF_Avg(:)
      real(8),allocatable :: TM_Avg(:)
      real(8),allocatable :: DM_Avg(:)
      real(8),allocatable :: ppm_out(:)

      logical(1) :: if_refresh       =.true.

      logical(1) :: if_flag_tmfb     =.false.
      logical(1) :: if_flag_tffb     =.false.
      logical(1) :: if_flag_xsfb     =.false.
      logical(1) :: if_flag_thfb     =.false.
      logical(1) :: if_opt_mode      =.false.
      logical(1) :: if_opt_xe        =.false.
      logical(1) :: if_opt_sm        =.false.
      logical(1) :: if_opt_gd        =.false.
      logical(4) :: if_OPT_findtavg  =.false.
      logical(1) :: if_ppm           =.false.
      logical(1) :: if_ppower        =.false.
      logical(1) :: if_core_massflow =.false.
      logical(1) :: if_tm_in         =.false.
      logical(1) :: if_tf_in         =.false.
      logical(1) :: if_d_TF          =.false.
      logical(1) :: if_d_TM          =.false.
      logical(1) :: if_cr_bot        =.false.
      logical(1) :: if_czp           =.false.
      logical(1) :: if_hzp           =.false.
      integer(4) :: if_ndr           = 0

      integer(4) :: opt_mtc          = 1

      logical(1) :: if_macro_xs      =.false.

      logical(1) :: flag_tmfb
      logical(1) :: flag_tffb
      logical(1) :: flag_xsfb
      logical(1) :: flag_thfb
      integer(4) :: opt_mode = 1
      integer(4) :: opt_xe
      integer(4) :: opt_sm
      integer(4) :: opt_gd
      integer(4) :: n_ppm
      integer(4) :: n_mt
      integer(4) :: n_ppower
      integer(4) :: case_step
      integer(4) :: n_crw
      integer(4) :: dcrm
      integer(4) ,DIMENSION(:),    allocatable:: i_crw
      logical(4) :: OPT_findtavg

      real(8)    :: ppm
      real(8)    :: ppower
      real(8)    :: core_massflow
      real(8)    :: tm_in
      real(8)    :: tf_in
      real(8)    :: d_TF
      real(8)    :: d_TM
      real(8)    :: unc_sdm !shutdown margin uncertainty
      real(8), allocatable :: cr_bot   (:)
      real(8), allocatable:: d_TM1     (:)
      real(8), allocatable:: NDR_ppower(:)
      real(8), allocatable:: NDR_ppm   (:)
      real(8), allocatable:: NDR_mt    (:)
      real(8), allocatable:: NDR_ft    (:)
      real(8), allocatable:: ndr1      (:)
      real(8), allocatable:: ndr2      (:,:)
      real(8), allocatable:: ndr2_1    (:,:)
   end type level1

   integer(4) :: n_case=0
   type(level1),allocatable :: branch_case(:)

   real(8) :: stack_keff(10)
   logical(1) :: flag_stack = .false.
   logical(1) :: flag_coef = .false.
   integer(4) :: crw_cases = 0
   character(3) :: c_coef=''
   integer(4)   :: c_case=0

end module inc_branch
