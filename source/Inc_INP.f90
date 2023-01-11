      MODULE Inc_INP

      ! 1) I_Rotate_1N        :: Rotation Index (1N/1FA)
      !                          90*Index Degree Clockwise Rotation for 4N/1FA
      !                               X    [1 Index Rotation]
      !                               -->
      !                               -------             -------
      !                          Y | | 1 | 2 |    =>    | 3 | 1 |
      !                            V |---|---|           |---|---|
      !                              | 3 | 4 |           | 4 | 2 |
      !                               -------             -------

      IMPLICIT NONE

      LOGICAL(1) :: Flag_Card_Title
      LOGICAL(1) :: Flag_Card_File
      LOGICAL(1) :: Flag_Card_Option
      LOGICAL(1) :: Flag_Card_Control
      LOGICAL(1) :: Flag_Card_Geometry
      LOGICAL(1) :: Flag_Card_RP
      LOGICAL(1) :: Flag_Card_Shuffling
      LOGICAL(1) :: Flag_Card_Rotation
      LOGICAL(1) :: Flag_Card_CR
      LOGICAL(1) :: Flag_Card_FA_Data
      LOGICAL(1) :: Flag_Card_TH_Data
      LOGICAL(1) :: Flag_Card_ATF_TCD
      LOGICAL(1) :: Flag_Card_Depletion
      LOGICAL(1) :: Flag_Card_Crud
      LOGICAL(1) :: Flag_Card_maXS
      LOGICAL(1) :: Flag_Card_Transient
      LOGICAL(1) :: Flag_Card_Branch


      integer, allocatable :: I_MRI_cycle(:)
      INTEGER, ALLOCATABLE :: I_Shuffle_1N(:)
      INTEGER, ALLOCATABLE :: I_Rotate_1N(:)
      INTEGER, ALLOCATABLE :: I_ATF_TCD_1N(:,:)
      INTEGER, ALLOCATABLE :: I_ATF_TCD(:,:)
      INTEGER, ALLOCATABLE :: I_ATF_TCD_FA_1N(:)
      INTEGER, ALLOCATABLE :: I_ATF_TCD_FA(:)
      INTEGER, ALLOCATABLE :: I_ATF_TCD_FA_RF(:)

      ! FOR TUNING
      logical(1) :: flag_axpow_cor=.false.
      real(8), allocatable :: axpow_inp(:)

      logical(1) :: flag_force_asi=.false.
      logical(1) :: flag_force_ao=.false.
      real(8) :: force_asi
      real(8), allocatable :: sigcor(:)
      real(8), allocatable :: bk_xsset_table(:,:,:,:)
      real(8), allocatable :: bk_xsset_table_base(:,:,:)

      logical(1) :: flag_force_tf=.false.
      real(8) :: force_tf=1d0
      real(8) :: force_dtf=0d0
      logical(1) :: flag_force_tm=.false.
      real(8) :: force_tm=1d0
      real(8) :: force_dtm=0d0
      logical(1) :: flag_force_tdop=.false.
      real(8) :: force_tdop=1d0
      real(8) :: force_dtdop=0d0

      logical(1) :: flag_force_ftc=.false.
      real(8) :: force_ftc=0d0
      logical(1) :: flag_force_mtc=.false.
      real(8) :: force_mtc=0d0

      logical(1) :: flag_force_ftctr=.false.
      real(8) :: force_ftctr=1d0
      real(8) :: force_dftctr=0d0
      logical(1) :: flag_force_mtctr=.false.
      real(8) :: force_mtctr=1d0
      real(8) :: force_dmtctr=0d0

      logical(1) :: flag_crw_cor=.false.
      real(8) :: crw_cor=1d0
      real(8) :: w_cr=1d0
      real(8) :: w_ucr=1d0
      logical(1) :: flag_crw_cor_each=.false.
      real(8), allocatable :: crw_cor_each(:)

      real(8) :: pu239_con     = 1.0d0 !1.017d0
      real(8) :: u235_con      = 1.0d0 !1.017d0
      real(8) :: scatup_con    = 1.0d0 !1.017d0
      real(8) :: effgd_con     = 1.0d0 !1.017d0
      real(8) :: abs2_residual = 1.0d0

      real(8) :: fuelk_con = 1d0
      real(8) :: gaph_con = 1d0

      logical(1) :: flag_force_beff=.false.
      real(8) :: force_beff(1:6)=1d0
      logical(1) :: flag_force_lamb=.false.
      real(8) :: force_lamb(1:6)=1d0
      logical(1) :: flag_force_velo=.false.
      real(8) :: force_velo(1:2)=1d0

      ! MKFCN_INPUT
      logical(1) :: flag_mkfcn=.false.
      integer :: mkfcn_cycle
      integer :: mkfcn_step
      integer :: mkfcn_ncycle
      character(30) :: folder_fcn
      character(30), allocatable :: name_fcninp(:)
      integer, allocatable :: cycle_fcninp(:)

      character(50) :: fcn_baseinp
      integer :: nline_fcnbase
      character(500), allocatable :: fcnbase(:)

      integer, allocatable :: fcn_ixiytoixy(:)
      integer :: fcn_nxya
      integer :: fcn_nxa
      integer :: fcn_nya
      integer :: fcn_npin
      integer :: fcn_nz
      real(8), allocatable :: fcn_height(:)
      real(8) :: fcn_pressure
      real(8) :: fcn_massflow
      real(8) :: fcn_inlettemp
      integer, allocatable :: fcn_burned(:)

      type :: ixya_type
         real(8), allocatable :: qmpy(:,:)
         real(8), allocatable :: qf(:,:,:)
      end type ixya_type
      type :: istep_type
         real(8) :: day
         type(ixya_type), allocatable :: ixya(:)
      end type istep_type
      type :: ifcn_type
         integer :: nstep
         integer, allocatable :: ixa(:)
         integer, allocatable :: iya(:)
         integer, allocatable :: ixya(:)
         integer, allocatable :: lp(:)
         integer, allocatable :: farf(:)
         integer, allocatable :: shuffle(:)
         integer, allocatable :: rotate(:)
         type(istep_type), allocatable :: step(:)
      end type ifcn_type
      type(ifcn_type), allocatable :: ifcn(:)

      type :: ostep_type
         real(8) :: day
         real(8), allocatable :: qmpy(:,:)
         real(8), allocatable :: qf(:,:,:)
      end type ostep_type
      type :: ofcn_type
         integer :: nstep
         type(ostep_type), allocatable :: step(:)
      end type ofcn_type
      type(ofcn_type), allocatable :: ofcn(:)

      ! MKFTN INPUT
      logical(1) :: flag_mkftn=.false.

      character(50) :: ftn_baseinp
      integer :: nline_ftnbase
      character(500), allocatable :: ftnbase(:)

      character(50) :: folder_ftn
      character(50) :: dir_fcn_rst
      integer :: ntrest
      integer, allocatable :: trest_index(:)
      real(8), allocatable :: trest(:)
      integer, allocatable :: trest_map(:)
      real(8), allocatable :: map2trest(:)

      logical(1) :: flag_ftnth_detail=.false.
      logical(1) :: switch_saveth=.false.
      real(8), allocatable :: ctf4ftninp(:,:,:,:,:,:,:)


      END MODULE Inc_INP

