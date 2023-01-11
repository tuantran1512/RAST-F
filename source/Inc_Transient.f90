
      MODULE Inc_Transient

      IMPLICIT NONE

      real(8) :: sec_old
      logical(1) :: flag_expand
      logical(1) :: flag_final

      INTEGER :: OPT_Xe_inTr
      INTEGER :: OPT_Sm_inTr
      INTEGER :: OPT_SmChain_inTr

      LOGICAL(1) :: Flag_Transient
      LOGICAL(1) :: Flag_DepBP_inTr

      REAL(8), DIMENSION(:), ALLOCATABLE :: dPPM_MovePPM
      REAL(8), DIMENSION(:), ALLOCATABLE :: MovePPM_Final
      REAL(8), DIMENSION(:), ALLOCATABLE :: v_MovePPM

      REAL(8), DIMENSION(:,:), ALLOCATABLE :: T_MovePPM
      INTEGER, DIMENSION(:), ALLOCATABLE :: N_MovePPM
      LOGICAL(1) :: Flag_MovePPM

      LOGICAL(1) :: Flag_Scram=.FALSE.
      LOGICAL(1) :: trip=.false. , scram=.false.
      REAL(8) :: powtrip, delaydel, scramdelt
      REAL(8) :: tripbeg, delayt
      REAL(8) :: scramstep
      real(8) :: scrambeg = 0d0
      logical(1) :: flag_otrip = .false.
      integer :: opt_trip = 0 ! 0: no trip / 1: power level / 2: power change / 3: 1 & 2
      integer :: trip_param = 0 ! 0: power / 1: excore detector
      real(8) :: dpowtrip,dttrip
      real(8),allocatable :: trip_powhist(:)
      real(8) :: trip_refexd, trip_refpow ! reference excore signal / corresponding power
      integer :: n_crstock = 0
      integer,allocatable :: i_crstock(:)

      logical(1) :: flag_print_kp_anc=.false.

      REAL(8) :: L_omega
      REAL(8), DIMENSION(:), ALLOCATABLE :: L_omega_Surf

      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: Q_Surf
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: Q_Surf_Old
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: C_d_I
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: C_d_I_Old
      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: C_d_I_Surf
      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: C_d_I_Surf_Old

      REAL(8), DIMENSION(:), ALLOCATABLE :: S_0
      real(8), dimension(:), allocatable :: S_1x
      real(8), dimension(:), allocatable :: S_1y
      real(8), dimension(:), allocatable :: S_1z
      real(8), dimension(:), allocatable :: S_2x
      real(8), dimension(:), allocatable :: S_2y
      real(8), dimension(:), allocatable :: S_2z

      REAL(8), DIMENSION(:,:), ALLOCATABLE :: S_Surf
      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: Flux_Surf
      REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: Flux_Surf_Old

      REAL(8), DIMENSION(:), ALLOCATABLE :: chi_p
      REAL(8), DIMENSION(:), ALLOCATABLE :: chi_d

      logical(1) :: flag_mslb=.false.
      integer :: slb_chan, slb_time, slb_chanp
      integer, dimension(:), allocatable :: i_chan
      integer, dimension(:), allocatable :: i_chan_fa
      integer, dimension(:), allocatable :: i_chan_fa_4n
      integer, dimension(:), allocatable :: i_chan_rf
      integer, dimension(:), allocatable :: i_chan_4n
      logical(1) :: flag_scen_tm=.false.
      logical(1) :: flag_scen_p=.false.
      logical(1) :: flag_chan_avg=.false.
      logical(1) :: flag_limit_tm=.false.
      logical(1) :: flag_limit_p=.false.
      character(200) :: name_scen_tm
      character(200) :: name_scen_p
      real(8), dimension(:,:), allocatable :: scen_tm
      real(8), dimension(:,:), allocatable :: scen_p
      real(8) :: limit_tm, limit_p
      integer :: opt_ctf=0

      logical(1) :: flag_tr_ocean = .FALSE.
      logical(1) :: flag_tr_ocean_CR = .FALSE.
      logical(1) :: flag_tr_ocean_PPM = .FALSE.
      logical(1) :: flag_tr_ocean_TM = .FALSE.

      integer :: ntime_troc
      real(8), dimension(:), allocatable :: bin_time_troc
      real(8), dimension(:,:), allocatable :: bin_cr_troc
      real(8), dimension(:), allocatable :: bin_ppm_troc
      real(8), dimension(:), allocatable :: bin_tm_troc

      logical(1) :: opt_tr_noftc=.false.
      logical(1) :: opt_tr_nomtc=.false.

      END MODULE Inc_Transient

