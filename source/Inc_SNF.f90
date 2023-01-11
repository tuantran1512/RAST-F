#ifdef JR_SNF

      MODULE Inc_SNF
      IMPLICIT NONE

      REAL(8)                               :: PCF_96242
      REAL(8)                               :: correct_flux
      LOGICAL(1)                            :: FLAG_SNF_VOL       = .FALSE.
      LOGICAL(1)                            :: FLAG_SNF_COO       = .FALSE.
      LOGICAL(1)                            :: FLAG_NO_PCF       = .FALSE.
      LOGICAL(1)                            :: Flag_SNF          = .FALSE.
      LOGICAL(1)                            :: Flag_SNF_FIRST    = .FALSE.
      LOGICAL(1)                            :: Flag_SNF_FA       = .FALSE.    ! To use print HFF
      LOGICAL(1)                            :: Flag_SNF_PIN      = .FALSE.    ! To use print HFF
      LOGICAL(1)                            :: Flag_SNF_PRINT    = .FALSE.    ! To use print HFF
      LOGICAL(1)                            :: Flag_SNF_ND_PRINT = .TRUE.
#ifdef JR_SNF_INTERPOLATION
      LOGICAL(1)                            :: Flag_calculation  = .TRUE.   ! TRUE: Lagrange, FALSE: Newton-Raphson
#endif
#ifdef JR_SNF_TMO
      REAL(8)                               :: Flag_tmo_order  = 4   ! TRUE: Lagrange, FALSE: Newton-Raphson
#endif
      LOGICAL(1)                            :: flag_hrst         = .FALSE.
      LOGICAL(1)                            :: flag_hri          = .FALSE.
      CHARACTER(200)                        :: NAME_ND
      LOGICAL(1)                            :: Flag_SNF_NFA      = .FALSE.
      INTEGER(4)                            :: SNF_NFA           = 1
      INTEGER(4),DIMENSION(:),ALLOCATABLE   :: SNF_NZ
      INTEGER(4),DIMENSION(:),ALLOCATABLE   :: SNF_Ixy_1N_M                 ! SNF POSITION : Ixy_1N
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: SNF_k_M
      INTEGER(4)                            :: SNF_Ixy_1N                 ! SNF POSITION : Ixy_1N
      INTEGER(4)                            :: SNF_k
      INTEGER(4)                            :: SNF_N_BU       = 0          ! Default, END
      INTEGER(4)                            :: SNF_N_END_CYBU = 0          ! Default, END
      INTEGER(4)                            :: SNF_PIN_TYPE   = 1          ! 1: PIN BY PIN (9GB/1FA), 2: ND distribution (NDD)
      REAL(8), ALLOCATABLE                  :: SNF_BU(:)

      !! decay
      REAL(8)                               :: SNF_DECAY_YEAR = 0d0
      INTEGER(4)                            :: SNF_ENDF  = 1 ! default ENDF70 (0:E68, 1:E70, 2:E71)
      LOGICAL(1)                            :: Flag_SNF_decay = .FALSE.

      INTEGER(4)                            :: SNF_PRINT_INDEX = 0
      INTEGER(4)                            :: SNF_PRINT_NN = 0
      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE   :: phihom_SNF    !(i_step,npin,npin,ng)       !(npin,npin,Nxy_1N,nz,ng)
      REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: phihom_M    !(i_step,npin,npin,ng)       !(npin,npin,Nxy_1N,nz,ng)
      REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: phihom_M_Old    !(i_step,npin,npin,ng)       !(npin,npin,Nxy_1N,nz,ng)
      character(200) :: Name_SNF
      character(200) :: Name_SNF_INTER
      character(200) :: Name_SNF_HRST
      character(200) :: Name_SNF_HRI
      INTEGER        :: W_SNF         = 9001     ! generates source term file
      INTEGER        :: R_ND          = 9002     !
      INTEGER        :: W_INTER       = 9003     ! generates number density file with interpolation method
      INTEGER        :: W_HRST        = 9004     ! generates number density file with interpolation method
      INTEGER        :: R_HRI         = 9005     ! generates number density file with interpolation method

      LOGICAL(1)                  :: FLAG_SNF_BU_INTER    ! Need to interpolation method
      INTEGER(4)                  :: SNF_N_BU_UP
      INTEGER(4)                  :: SNF_N_BU_DOWN
      REAL(8)                     :: SNF_H = 1.D0         ! HEIGHT used in  branch calculation
      REAL(8)                     :: SNF_MIN_VAL = 0d0

      ! INDEX
      REAL(8)                     :: SNF_TMO_INDEX
      REAL(8)                     :: SNF_TFU_INDEX
      REAL(8)                     :: SNF_BOR_INDEX

      ! VOLUME RATIO
      REAL(8), DIMENSION(:,:), ALLOCATABLE      :: SNF_VOL_RATIO
      REAL(8)       :: SNF_VOL

      ! NDENTOGCC
      REAL(8), DIMENSION(:), ALLOCATABLE        :: HNDEN2GCC

      ! PIN MAT
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE     :: SNF_PIN_MAT
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE     :: SNF_PIN_POSITION ! ix,iy
      INTEGER(4)                                :: SNF_MAX_NR
      INTEGER(4)                                :: SNF_MIN_MAT
      INTEGER(4)                                :: SNF_MAX_MAT

      ! CONDITION, TFU/TMO/BOR
      INTEGER(4)                  :: SNF_BRAN_OPT
      REAL(8)                     :: SNF_BRAN_TFU
      REAL(8)                     :: SNF_BRAN_TMO
      REAL(8)                     :: SNF_BRAN_BOR
      INTEGER(4)                  :: SNF_N_ISO
      INTEGER(4)                  :: SNF_N_STEP
      REAL(8)                     :: SNF_BRAN_P
      INTEGER(4)                  :: SNF_N_INDEXING
      REAL(8),DIMENSION(:,:),ALLOCATABLE        :: PCF_RESULT

      ! Power correction factor
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: SNF_PCF_LIST

      ! BRANCH FLUX, BRANCH XS(0 BU)
      REAL(8),DIMENSION(:,:,:),allocatable     :: SNF_BRAN_FLUX
      REAL(8),DIMENSION(:),allocatable       :: SNF_XS   ! ORIGEN LIBRARY

      ! Read .ND FILE
      REAL(8),DIMENSION(:),allocatable       :: SNF_BRAN_BUHIST
      REAL(8),DIMENSION(:),allocatable       :: SNF_BRAN_DAYHIST
      REAL(8),DIMENSION(:),allocatable       :: SNF_LAMBDA_ARRAY  ! (N_ISO)
      INTEGER(4),DIMENSION(:),allocatable       :: SNF_ISO_ARRAY  ! (N_ISO)
      REAL(8),DIMENSION(:,:,:),allocatable     :: SNF_REF_ARRAY  ! (U/D, N_ISO)
      REAL(8),DIMENSION(:,:,:,:),allocatable   :: SNF_TMO_ARRAY  ! (U/D, N_ISO, N_HIST)
      REAL(8),DIMENSION(:,:,:,:),allocatable   :: SNF_TFU_ARRAY  ! (U/D, N_ISO, N_HIST)
      REAL(8),DIMENSION(:,:,:,:),allocatable   :: SNF_BOR_ARRAY  ! (U/D, N_ISO, N_HIST)
      REAL(8),DIMENSION(:,:),allocatable       :: SNF_RESULT_ARRAY  ! (U/D, N_ISO, N_HIST)

      !REAL(8),DIMENSION(:,:,:),allocatable   :: SNF_BRAN_NDD   ! (BU, N_ISO, PIN
      !REAL(8),DIMENSION(:,:,:),allocatable   :: SNF_R2_NDD   ! (BU, N_ISO, PIN

      ! GEOMETRY - NR
      !INTEGER(4),DIMENSION(:,:),allocatable   :: SNF_PIN_POSITION ! (FA #, x_pin, y_pin)
      !INTEGER(4),DIMENSION(:,:),allocatable   :: SNF_PIN_RAD      ! RADIUS INFORMATION ! MAX NR:15


      INTEGER(4),DIMENSION(:),allocatable :: SNF_ISO_LIST
      INTEGER(4),DIMENSION(:),allocatable :: SNF_ISO_TYPE
      INTEGER(4),DIMENSION(:),allocatable :: SNF_PRE_1
      INTEGER(4),DIMENSION(:),allocatable :: SNF_PRE_2

      ! PRINT SUMMARY
      REAL(8),DIMENSION(:,:),allocatable  :: SNF_PRINT_SUM


      real(8), dimension(:,:), allocatable :: SNF_N_U34
      real(8), dimension(:,:), allocatable :: SNF_N_U35
      real(8), dimension(:,:), allocatable :: SNF_N_U36
      real(8), dimension(:,:), allocatable :: SNF_N_U37
      real(8), dimension(:,:), allocatable :: SNF_N_U38
      real(8), dimension(:,:), allocatable :: SNF_N_Np37
      real(8), dimension(:,:), allocatable :: SNF_N_Np38
      real(8), dimension(:,:), allocatable :: SNF_N_Np39
      real(8), dimension(:,:), allocatable :: SNF_N_Pu38
      real(8), dimension(:,:), allocatable :: SNF_N_Pu39
      real(8), dimension(:,:), allocatable :: SNF_N_Pu40
      real(8), dimension(:,:), allocatable :: SNF_N_Pu41
      real(8), dimension(:,:), allocatable :: SNF_N_Pu42
      real(8), dimension(:,:), allocatable :: SNF_N_Pu43
      real(8), dimension(:,:), allocatable :: SNF_N_Am41
      real(8), dimension(:,:), allocatable :: SNF_N_As42
      real(8), dimension(:,:), allocatable :: SNF_N_Am42
      real(8), dimension(:,:), allocatable :: SNF_N_Am43
      real(8), dimension(:,:), allocatable :: SNF_N_Am44
      real(8), dimension(:,:), allocatable :: SNF_N_Cm42
      real(8), dimension(:,:), allocatable :: SNF_N_Cm43
      real(8), dimension(:,:), allocatable :: SNF_N_Cm44

      real(8), dimension(:,:), allocatable :: SNF_N_U34_Old
      real(8), dimension(:,:), allocatable :: SNF_N_U35_Old
      real(8), dimension(:,:), allocatable :: SNF_N_U36_Old
      real(8), dimension(:,:), allocatable :: SNF_N_U37_Old
      real(8), dimension(:,:), allocatable :: SNF_N_U38_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Np37_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Np38_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Np39_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Pu38_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Pu39_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Pu40_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Pu41_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Pu42_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Pu43_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Am41_Old
      real(8), dimension(:,:), allocatable :: SNF_N_As42_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Am42_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Am43_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Am44_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Cm42_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Cm43_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Cm44_Old

      integer :: SNF_N_Heavy

      real(8), dimension(:,:), allocatable :: SNF_N_I35
      real(8), dimension(:,:), allocatable :: SNF_N_Xe35
      real(8), dimension(:,:), allocatable :: SNF_N_Nd47
      real(8), dimension(:,:), allocatable :: SNF_N_Nd48
      real(8), dimension(:,:), allocatable :: SNF_N_Nd49
      real(8), dimension(:,:), allocatable :: SNF_N_Pm47
      real(8), dimension(:,:), allocatable :: SNF_N_Ps48
      real(8), dimension(:,:), allocatable :: SNF_N_Pm48
      real(8), dimension(:,:), allocatable :: SNF_N_Pm49
      real(8), dimension(:,:), allocatable :: SNF_N_Sm47
      real(8), dimension(:,:), allocatable :: SNF_N_Sm48
      real(8), dimension(:,:), allocatable :: SNF_N_Sm49

      real(8), dimension(:,:), allocatable :: SNF_N_I35_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Xe35_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Nd47_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Nd48_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Nd49_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Pm47_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Ps48_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Pm48_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Pm49_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Sm47_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Sm48_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Sm49_Old

      integer :: SNF_N_NdSm

      real(8), dimension(:,:), allocatable :: SNF_N_Gd52
      real(8), dimension(:,:), allocatable :: SNF_N_Gd54
      real(8), dimension(:,:), allocatable :: SNF_N_Gd55
      real(8), dimension(:,:), allocatable :: SNF_N_Gd56
      real(8), dimension(:,:), allocatable :: SNF_N_Gd57
      real(8), dimension(:,:), allocatable :: SNF_N_Gd58
      real(8), dimension(:,:), allocatable :: SNF_N_Gd60

      real(8), dimension(:,:), allocatable :: SNF_N_Gd52_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Gd54_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Gd55_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Gd56_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Gd57_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Gd58_Old
      real(8), dimension(:,:), allocatable :: SNF_N_Gd60_Old

      integer :: SNF_N_Gd

      !real(8), dimension(:,:), allocatable :: SNF_N_H2O
      !real(8), dimension(:,:), allocatable :: SNF_N_B0

      ! Save data
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Xe35
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Sm49
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_I35
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Nd47
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Nd48
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Nd49
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Pm47
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Ps48
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Pm48
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Pm49
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Sm47
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Sm48
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Gd52
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Gd54
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Gd55
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Gd56
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Gd57
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Gd58
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Gd60
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_U34
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_U35
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_U36
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_U37
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_U38
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Np37
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Np38
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Np39
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Pu38
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Pu39
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Pu40
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Pu41
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Pu42
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Pu43
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Am41
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_As42
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Am42
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Am43
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Am44
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Cm42
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Cm43
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SNF_Hist_N_Cm44

      ! HRI
      INTEGER(4)                           :: SNF_HRI_N
      REAL(8)                           :: SNF_HRI_LAST_DAY
      REAL(8), DIMENSION(:), ALLOCATABLE   :: SNF_HRI_BU
      REAL(8), DIMENSION(:), ALLOCATABLE   :: SNF_HRI_DAY
      REAL(8), DIMENSION(:), ALLOCATABLE   :: SNF_HRI_TMO
      REAL(8), DIMENSION(:), ALLOCATABLE   :: SNF_HRI_TFU
      REAL(8), DIMENSION(:), ALLOCATABLE   :: SNF_HRI_BOR
      REAL(8), DIMENSION(:), ALLOCATABLE   :: SNF_HRI_FLUX

      ! COOLING :: DEPLETION
#ifdef JR_SRCTRM
      REAL(8),PARAMETER                             :: BARN=1.0d-24
      CHARACTER(300)                                :: DECAY_BR=''
      CHARACTER(300)                                :: DECAY_MAT=''
      INTEGER(4)                                    :: dep_ver=68
      !! Decay mode
      ! For E6 decay lib
      INTEGER(4),PARAMETER                        :: FR_FBX=1                         ! Negatron beta decay transition to excited state
      INTEGER(4),PARAMETER                        :: FR_FPEC=2                        ! Positron emission or electron capture to ground state
      INTEGER(4),PARAMETER                        :: FR_FPECX=3                       ! Positron emission or electron capture to excited state
      INTEGER(4),PARAMETER                        :: FR_FA=4                          ! Alpha decay to ground state
      INTEGER(4),PARAMETER                        :: FR_FIT=5                         ! Decay from Excited state to ground state
      INTEGER(4),PARAMETER                        :: FR_FSF=6                         ! Spontaneous fission
      INTEGER(4),PARAMETER                        :: FR_FN=7                          ! (Beta+neutron) decay to ground state
      INTEGER(4),PARAMETER                        :: FR_FB=8                          ! Negatron beta decay transition to ground state
      ! For E7 decay lib
      INTEGER(4),PARAMETER                        :: FR_FBG=1                         ! Negatron beta decay transition, ground state daughter
      INTEGER(4),PARAMETER                        :: FR_FBE1=2                        ! Negatron beta decay transition, 1st exicited state daughter
      INTEGER(4),PARAMETER                        :: FR_FBE2=3                        ! Negatron beta decay transition, 2nd exicited state daughter
      INTEGER(4),PARAMETER                        :: FR_FBB=4                         ! Negatron beta decay followed by a beta decay
      INTEGER(4),PARAMETER                        :: FR_FBA=5                         ! Negatron beta decay followed by alpha emission
      INTEGER(4),PARAMETER                        :: FR_FBN=6                         ! Negatron beta decay followed by neutron emission (delayed neutron de
      INTEGER(4),PARAMETER                        :: FR_FBNA=7                        ! Negatron beta decay followed by neutron + alpha emission (delayed ne
      INTEGER(4),PARAMETER                        :: FR_FB2N=8                        ! Negatron beta decay followed by 2 neutron emission (delayed neutron
      INTEGER(4),PARAMETER                        :: FR_FB3N=9                        ! Negatron beta decay followed by 3 neutron emission (delayed neutron
      INTEGER(4),PARAMETER                        :: FR_FEC=10                        ! Positron emission or electron capture ground state daughter
      INTEGER(4),PARAMETER                        :: FR_FEC1=11                       ! Positron emission or electron capture 1st exicited state daughter
      INTEGER(4),PARAMETER                        :: FR_FECA=12                       ! Positron emission or electron capture followed by alpha emission
      INTEGER(4),PARAMETER                        :: FR_FECP=13                       ! Positron emission or electron capture by proton emission
      INTEGER(4),PARAMETER                        :: FR_FECPP=14                      ! Positron emission or electron capture by 2 proton emission
      INTEGER(4),PARAMETER                        :: FR_FIS=15                        ! isomeric transition
      INTEGER(4),PARAMETER                        :: FR_FISE=16                       ! isomeric transition 1st excited state daughter
      INTEGER(4),PARAMETER                        :: FR_FA_7=17                       ! Alpha decay
      INTEGER(4),PARAMETER                        :: FR_FAE=18                        ! Alpha decay 1st excited state daughter
      INTEGER(4),PARAMETER                        :: FR_FN_7=19                       ! Neutron emission
      INTEGER(4),PARAMETER                        :: FR_FNN=20                        ! Neutron emission followed by a neutron emission
      INTEGER(4),PARAMETER                        :: FR_FSF_7=21                      ! Spontaneous fission
      INTEGER(4),PARAMETER                        :: FR_FP=22                         ! proton emission
      INTEGER(4),PARAMETER                        :: FR_FPP=23                        ! proton emission followed by a proton emission

      !For E7.1 Decay lib
      INTEGER(4),PARAMETER                        :: FR_FB_71=1                       ! Negatron beta decay
      INTEGER(4),PARAMETER                        :: FR_FBX_71=2                      ! Negatron beta decay to exited state(1st)
      INTEGER(4),PARAMETER                        :: FR_FB2X_71=3                     ! Negatron beta decay to exited state(2nd)
      INTEGER(4),PARAMETER                        :: FR_F2B_71=4                      ! Negatron beta decay followed by beta decay
      INTEGER(4),PARAMETER                        :: FR_FBA_71=5                      ! Negatron beta decay followed by alpha emission
      INTEGER(4),PARAMETER                        :: FR_FBN_71=6                      ! Negatron beta decay followed by neutron emission
      INTEGER(4),PARAMETER                        :: FR_FB2N_71=7                     ! Negatron beta decay followed by neutron emission(2n)
      INTEGER(4),PARAMETER                        :: FR_FB3N_71=8                     ! Negatron beta decay followed by neutron emission(3n)
      INTEGER(4),PARAMETER                        :: FR_FB4N_71=9                     ! Negatron beta decay followed by neutron emission(4n)
      INTEGER(4),PARAMETER                        :: FR_FE_71=10                      ! Decay by positron emission or electron capture
      INTEGER(4),PARAMETER                        :: FR_FEX_71=11                     ! Decay by positron emission or electron capture to exited state(1st)
      INTEGER(4),PARAMETER                        :: FR_FE2X_71=12                    ! Decay by positron emission or electron capture to exited state(2nd)
      INTEGER(4),PARAMETER                        :: FR_FEA_71=13                     ! Decay by positron emission or electron capture followed by alpha emission
      INTEGER(4),PARAMETER                        :: FR_FESF_71=14                    ! Decay by positron emission or electron capture followed by spontaneous fission
      INTEGER(4),PARAMETER                        :: FR_FEP_71=15                     ! Decay by positron emission or electron capture followed by proton emission
      INTEGER(4),PARAMETER                        :: FR_FE2P_71=16                    ! Decay by positron emission or electron capture followed by proton emission(2p)
      INTEGER(4),PARAMETER                        :: FR_FITG_71=17                    ! Isomeric transition from exited state to ground state
      INTEGER(4),PARAMETER                        :: FR_FITX_71=18                    ! Isomeric transition from exited state to lower exited state(1st)
      INTEGER(4),PARAMETER                        :: FR_FA_71=19                      ! Alpha decay
      INTEGER(4),PARAMETER                        :: FR_FAX_71=20                     ! Alpha decay to exited state(1st)
      INTEGER(4),PARAMETER                        :: FR_FN_71=21                      ! Decay by neutron emission
      INTEGER(4),PARAMETER                        :: FR_F2N_71=22                     ! Decay by neutron emission (2n)
      INTEGER(4),PARAMETER                        :: FR_FSF_71=23                     ! Spontaneous fission
      INTEGER(4),PARAMETER                        :: FR_FP_71=24                      ! Decay by proton emission
      INTEGER(4),PARAMETER                        :: FR_F2P_71=25                     ! Decay by proton emission(2p)

!      INTEGER(4),ALLOCATABLE                      :: m_rex(:)                             !
      INTEGER(4),PARAMETER                        :: n_CRAM_LIB_XS=8
      INTEGER(4)                                  :: optim_val=68                     !                             @@@ DO NOT END
!      CHARACTER(300)                              :: RXfile=''                          !                             @@@ DO NOT END
!      INTEGER(4),ALLOCATABLE                      :: type_rex(:)                        !
      INTEGER(4),PARAMETER                        :: XS_ALPHA=4                       ! (n,alpha) reaction to ground state
      INTEGER(4),PARAMETER                        :: XS_FIS=5                         ! (n,fission)
      INTEGER(4),PARAMETER                        :: XS_N2NX=8                        ! (n,2n) reaction to excited state
      INTEGER(4),PARAMETER                        :: XS_N2N=2                         ! (n,2n) reaction to ground state
      INTEGER(4),PARAMETER                        :: XS_N3N=3                         ! (n,3n) reaction to ground state
      INTEGER(4),PARAMETER                        :: XS_NGX=7                         ! (n,gamma) reaction to excited state
      INTEGER(4),PARAMETER                        :: XS_NG=1                          ! (n,gamma) reaction to ground state
      INTEGER(4),PARAMETER                        :: XS_NP=6                          ! (n,proton) reaction to ground state
#endif

      END MODULE Inc_SNF

#endif
