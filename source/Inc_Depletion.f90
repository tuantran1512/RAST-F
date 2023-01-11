
      module Inc_Depletion
      implicit none

      integer :: N_BU
      integer :: n_day
      real(8), allocatable :: Cycle_BU(:)
      real(8), allocatable :: cycle_day(:)
      real(8) :: dBU
      real(8) :: RST_BU
      real(8), allocatable :: MRST_BU(:)
      real(8) :: Tot_MTU

      ! input history
      real(8), allocatable :: inp_dT(:)
      real(8), allocatable :: inp_dBU(:)
      real(8), allocatable :: inp_ppower(:)
      real(8), allocatable :: inp_core_power(:)
      real(8), allocatable :: inp_fa_power(:)

      logical(1) :: if_buinput=.true.
      logical(1) :: Flag_miDepl
      logical(1) :: Flag_stop_neg_ppm =.true.
      logical(1) :: flag_effgad=.true.
      integer :: dep_ver = 0

      ! PC
      real(8), parameter :: PC_w = 0.5d0
      logical(1) :: Flag_PC
      integer :: opt_pc=1 ! 0: Full-PC (Original)
                          ! 1: No PC (Only P)
                          ! 2: Semi-PC
                          ! 3: Full-PC (new) - default

      ! CRAM
      real(8) :: alpha_0
      integer :: Order_CRAM
      complex(8), allocatable :: Mat_I_HN(:,:)
      complex(8), allocatable :: Mat_I_FP(:,:)
      complex(8), allocatable :: alpha(:)
      complex(8), allocatable :: theta(:)

      ! BU adaptation
      logical(1) :: Flag_BU_Adap=.FALSE.
      logical(1) :: Flag_BU_AdapINP=.FALSE.
      logical(1) :: Flag_BU_A_StandAlone = .FALSE.
      logical(1), allocatable :: Flag_Ada_FA(:)
      integer, allocatable :: ID_SubBatch(:,:)
      integer, allocatable :: N_FA_LP(:)
      integer :: I_BU_Ada
      integer :: OPT_ADA
      integer :: N_SubBatch
      real(8) :: MB
      real(8) :: MB_Min
      real(8) :: MB_Max
      real(8) :: dMB
      real(8), allocatable :: MB_SubBatch(:)
      integer, allocatable :: ID_BU_A(:)
      logical(1), allocatable :: Flag_Ada_BU_Step(:)
      real(8), allocatable :: Power_1N_Ref(:,:)

      ! B10 depletion
      logical(1) :: Flag_b10dep=.FALSE.
      real(8) :: B10PCT = 19.8D0 ! STREAM default
      real(8) :: B10DEP = 19.8d0 ! if B10DEP <= 0, B10DEP=B10PCT (at BOC)
      real(8) :: BORINB10 = 19.90D0
      real(8) :: PPMB10 = 0d0
      real(8) :: RCSVOL = 81634.9 ! amount of coolant in RCS for OPR-1000, in gallons
      real(8) :: CORVOL = 4665.4  ! ratio (CORVOL/RCSVOL) is required... for OPR-1000, in gallons
      real(8) :: GFACT = 1.00d0    ! should be modified... see sample input
      real(8) :: AVGRHO = 0.70090d0   ! Core average water density
      real(8) :: BORINFLOW = 0d0      ! Borated water added into RCS over time step
      real(8) :: BORINCONC = 0.7815d0 ! Soluble boron level in borated water
      real(8) :: MAKEFLOW  = 0d0      ! Total water (borated and unborated) added

      ! 10ppm search
      logical(1) :: Flag_10ppm = .false.
      logical(1) :: if_get_interpol = .false.
      logical(1) :: if_conv_10ppm = .false.
      real(8) :: targetppm = 10d0

      ! CR depletion
      logical(1) :: flag_crdep = .false.


      end module Inc_Depletion
