
      MODULE Inc_Option

      USE Inc_Constant

      IMPLICIT NONE

      integer :: OPT_Mode
      integer :: OPT_RST
      integer :: OPT_Xe
      integer :: OPT_Sm
      integer :: OPT_SmChain
      integer :: OPT_Nodal
      integer :: OPT_TLA
      integer :: OPT_TFcal
      integer :: N_Group

      integer :: OPT_Gd = 1

      logical(4) :: OPT_Jumpin

      INTEGER, DIMENSION(:) :: OPT_TCD(3) = 0

      INTEGER :: OPT_Tburn = 0
      INTEGER :: OPT_OXcal = 0
      INTEGER :: OPT_BUh = 0

      INTEGER :: OPT_RST_FOLD
      INTEGER :: OPT_RST_FOLD_BC=1

      logical(1) :: flag_print_design=.false.
      logical(1) :: flag_print_kp=.false.

      END MODULE Inc_Option
