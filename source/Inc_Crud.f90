      MODULE Inc_Crud
      IMPLICIT NONE
      LOGICAL :: OPT_CRUD=.false.
      LOGICAL :: OPT_READCRUD=.false.
      LOGICAL :: OPT_CRUDBURN=.TRUE.
      LOGICAL :: OPT_CRUDASI=.FALSE.
      LOGICAL :: OPT_CRUDFAASI=.FALSE.
      LOGICAL :: OPT_CRUDCBC=.TRUE.
      real(8) :: EPS_CRUDASI=0.005d0
      REAL(8), ALLOCATABLE :: NB0_Table(:,:) ! N_XS_Table, n_SP
      REAL(8), ALLOCATABLE :: N_crudb0(:,:),n_crudb0_old(:,:),Hist_n_crudb0(:,:,:),hist_n_crudb0_z(:,:)
      character(200) :: crud_file=''
!! BOA
      integer(4) :: n_day_boa=0
      real(8),allocatable :: day_boa(:)
      LOGICAL :: OPT_CRUDBOA=.false.
      character(200) :: crudboa_file=''
      real(8),allocatable :: Boa_b0(:,:,:)  ! (NXY, IZ, BU)
      real(8),allocatable :: Boa_dr(:,:,:)  ! (NXY, IZ, BU)
      integer(4) :: nz_boa
      real(8),allocatable :: z_boa(:)
!! BOA

      real(8),allocatable :: Crud_table(:,:,:)  ! (NXY, IZ, BU)

      character(200) :: cruddr_file=''
      real(8),allocatable :: Cruddr_table(:,:,:)  ! (NXY, IZ, BU)
      character(200) :: FAASI_file=''
      real(8),allocatable :: FAASI_INP(:,:,:)  ! (IX,IY,BU)
      REAL(8),allocatable :: itr_FA_BD (:,:,:) !
      REAL(8),allocatable :: itr_FA_ASI(:,:,:) !
      REAL(8),allocatable :: slope_FA_BD(:,:) !
      real(8) :: crtf=0.0
      real(8) :: burn_factor=1.0d0

      REAL(8), ALLOCATABLE :: cruddr(:,:)
      REAL(8), ALLOCATABLE :: cruddr_old(:,:)
      REAL(8), ALLOCATABLE :: cruddtf(:,:)
      REAL(8) :: LHS=250.0d0
      REAL(8) :: Crud_BNDEN=2.5844190d-02*1d+24 ! (#/cm3)
      REAL(8),allocatable :: TD_crud_bnden(:)
      REAL(8) :: crud_poro=0.75d0 ! Default: 75 %
      REAL(8) :: crud_drf=80d0/3600d0 ! 80 micro-m at 60 burnup
      REAL(8),allocatable :: TD_crud_drf(:)
      REAL(8),allocatable :: TD_crud_ASI(:)
      logical(1),allocatable :: IF_ASI_search(:)
      real(8),allocatable :: B0_RRa(:,:)


      real(8) :: ASI_o=0.0d0,ASI_oo=0.0d0
      real(8) :: NB_o=0.0d0,NB_oo=0.0d0
      integer(4):: n_CRUDASI=0
      real(8),allocatable :: total_Crud(:)

      real(8), allocatable :: Ind_BU(:,:) ! Independent BU, not cumulative
      real(8), allocatable :: Ind_BU_Old(:,:) ! Independent BU, not cumulative
      real(8), allocatable :: Ind_BU_Predictor(:,:) ! Independent BU, not cumulative

      integer :: opt_crud_model = 1           ! 0=old / 1=revised version implemented by KHNP-CRI
      real(8) :: cruddr_thres = 0.d0          ! threshold of crud thickness (micron) for precipitation of LiBO2
      real(8), allocatable :: SNB_mass(:,:)   ! mass of subcooled nucleate boiling for each node (kg)
      real(8) :: max_SNB_mass = 0.d0          ! maximum mass of subcooled nucleate boiling (kg)
      real(8) :: boiling_factor = 1.d-7       ! growth rate of crud thickness by SNB_mass (micron/kg)

#ifdef hj_training
      logical :: opt_crud_trset = .false. !training set
      logical :: opt_crud_depl
      integer :: crud_nfa = 0
      integer, allocatable :: crud_ifa(:)
      real(8), allocatable :: bogus_bu(:)
      real(8), allocatable :: crudmult(:)
      real(8), allocatable :: Hist_crudt(:,:,:), Hist_crudt_Z(:,:)
#endif

      END MODULE Inc_Crud
