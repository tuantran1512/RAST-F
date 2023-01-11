
      MODULE Inc_File

      IMPLICIT NONE

      INTEGER :: Len_INP, N_VIPRE_LINE = 0, N_VIPRE_STEP = 1
      character(200) :: Name_XS
      character(200) :: Name_RI
      character(200) :: Name_QS
      character(200) :: Name_OUT
      character(200) :: Name_ANC
      character(200) :: Name_SUM
      character(200) :: Name_RST
      character(200) :: Name_POW
      character(200) :: Name_FLUX
      character(200) :: Name_TH
      character(200) :: Name_maXS
      character(200) :: Name_INP
      character(200) :: INP_Title
      character(200) :: Name_DRFB
      character(200) :: Name_DRFM
      character(200) :: Name_DRFT
      character(200) :: Name_BU_A
      character(200) :: Name_ADA
      character(200) :: Name_OCEAN
      character(200) :: Name_VIPRE_OUT
      character(200) :: Name_VIPRE_IN
      character(200) :: name_inp2
      character(200) :: Name_MAP
      character(200), dimension(:), allocatable ::  Name_MRI
      integer, dimension(:), allocatable ::  Cycle_MRI
      integer :: No_MRI, Cycle_Now = 1

      character(200) :: name_sim
      logical(1) :: flag_sim=.false.
      integer :: no_axial = 0
      integer, allocatable :: axial_cont(:)

      character(200) :: name_adj
      character(200) :: name_pkd
      character(200) :: name_rho

      LOGICAL(1) :: Flag_XS
      LOGICAL(1) :: Flag_RI
      LOGICAL(1) :: Flag_QS
      LOGICAL(1) :: Flag_OUT
      LOGICAL(1) :: Flag_ANC
      LOGICAL(1) :: Flag_SUM
      LOGICAL(1) :: Flag_RST
#ifdef tuan_fr_rst
      LOGICAL(1) :: Flag_RST_HEX
      LOGICAL(1) :: Flag_RI_HEX
#endif
      LOGICAL(1) :: Flag_maXS
      LOGICAL(1) :: Flag_DRFB
      LOGICAL(1) :: Flag_DRFM
      LOGICAL(1) :: Flag_DRFT
      LOGICAL(1) :: Flag_BU_A
      LOGICAL(1) :: Flag_OCEAN
      LOGICAL(1) :: Flag_VIPRE
      logical(1) :: Flag_POW  = .false.
      logical(1) :: Flag_FLUX = .false.
      logical(1) :: Flag_TH   = .false.
      logical(1) :: Flag_MAP = .false.
      logical(1) :: flag_ocean_first = .TRUE.
      LOGICAL(1) :: Flag_MRI=.FALSE.

      logical(1) :: flag_print_core_param = .false.
      logical(1), dimension(:), allocatable :: flag_out_print

      logical(1) :: flag_fcn=.false.
      character(200) :: Name_fcn
      integer :: W_fcn = 2041

      INTEGER :: R_INP   = 1001
      INTEGER :: r_inp2
      INTEGER :: R_XS    = 1002
      INTEGER :: R_RI    = 1003
      INTEGER :: R_QS    = 1004
      INTEGER :: R_VIPRE = 1008

      integer,dimension(:),allocatable :: R_MRI

      INTEGER :: W_OUT   = 2001
      INTEGER :: W_SUM   = 2002
      INTEGER :: W_RST   = 2003
      INTEGER :: W_POW   = 2004
      INTEGER :: W_FLUX  = 2004
      INTEGER :: W_TH    = 2005
      INTEGER :: W_maXS  = 2006
      INTEGER :: W_ADA   = 2007
      INTEGER :: W_OCEAN = 2008
      INTEGER :: W_VIPRE = 2009
      integer :: W_MAP   = 2011
      INTEGER :: W_SIM   = 2010
      INTEGER :: W_ANC   = 2012
#ifdef tuan_fr_rst
      INTEGER :: W_RST_HEX = 2023
#endif


      INTEGER, PARAMETER :: W_ADJ = 2099
      INTEGER, PARAMETER :: W_PKD = 2098
      INTEGER, PARAMETER :: W_RHO = 2097

      LOGICAL(1) :: Flag_Open
#ifdef tuan_fr
      character(200) :: Name_XS_HEX
      character(1), DIMENSION(:), ALLOCATABLE :: Fiss
      integer :: NuNum
      integer :: NGroup
      integer :: NState
      INTEGER(4),ALLOCATABLE  :: N_TF(:)                   
      INTEGER(4),ALLOCATABLE  :: N_DC(:)                   
!      integer :: N_TF
!      integer :: N_DC
      LOGICAL(1) :: Flag_XS_Hex

!      integer :: NState
      Real(8), DIMENSION(:), ALLOCATABLE :: FuelTemp
      Real(8), DIMENSION(:), ALLOCATABLE :: CoolDens

#endif


      END MODULE Inc_File
