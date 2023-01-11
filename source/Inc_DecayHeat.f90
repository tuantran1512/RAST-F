
      module Inc_DecayHeat

      real(8), allocatable                :: Decay_Power(:, :) ! -/-
      real(8), allocatable                :: DH_precurs(:, :, :, :) ! -/-
      real(8), allocatable                :: Flux_dh_old(:, :, :)
      real(8)                             :: total_DH
      real(8)                             :: total_FP
      real(8)                             :: tn,tn_Old
      real(8), dimension(1:6)             :: C_tn
      real(8), dimension(1:6)             :: C_tn_Old
      real(8), dimension(1:6)             :: dC_tn
      real(8), dimension(1:6)             :: dC_tn_Old
      real(8), dimension(1:6)             :: dP_tn
      real(8), dimension(1:6)             :: dP_tn_Old
      real(8)                             :: P_dec_heat
      real(8)                             :: beta_DH_ratio
      real(8)                             :: beta_parcs = 0.07015986
      real(8)                             :: beta_dunn  = 0.0697094
      real(8), parameter                  :: NONE       = 0.0
      real(8)                             :: curr_time_dh
      real(8)                             :: curr_time_dh_old
      logical                             :: use_decay_heat = .false.
      logical                             :: init_decay_heat = .false.
      logical                             :: print_decay_heat = .false.
      integer                             :: print_dh
      integer                             :: mode_dh
      integer                             :: init_time_dh
      real(8), allocatable                :: dep_total_DH(:, :) ! save value before transient


      ! constants taken from given sources =-=-=-=-=-=-=-=-=-=-=-=-=

      ! user input
      real(8), dimension(1:6)             :: lam6_user = (/0.091825, &
         0.0068277,0.00048149,4.5839E-5,3.3416E-6,1.0734E-11/)
      real(8), dimension(1:6)             :: beta6_user = (/0.025032, &
         0.017649,0.013445,0.0066959,0.0035552,0.0033323/)
      real(8)                             :: beta_user = 0.0

      ! from PARCS manual
      real(8), dimension(1:6)             :: beta6_parcs = (/2.35402E-2, &
         1.89077E-2,1.39236E-2,6.90315E-3,3.56888E-3,3.31633E-3/)

      real(8), dimension(1:6)             :: lam6_parcs = (/1.05345E-1, &
         8.37149E-3,5.20337E-4,4.73479E-5,3.28153E-6,1.17537E-11/)

      ! from F. Dunn paper
      real(8), dimension(1:6)             :: beta6_dunn = (/0.025032, &
         0.017649,0.013445,0.0066959,0.0035552,0.0033323/)

      real(8), dimension(1:6)             :: lam6_dunn = (/0.091825, &
         0.0068277,0.00048149,4.5839E-5,3.3416E-6,1.0734E-11/)
      ! NOTE: constants for 12 group (F.Dunn) are not available because
      ! F.Dunn stated 12 lambdas and only 9 betas in his article
      ! (i.e. not matching data/not enough data) =-=-=-=-=-=-=-=-=-=

      end module Inc_DecayHeat
