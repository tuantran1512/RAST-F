
      MODULE Inc_TH
      ! *** Variable ***
      ! 1)  Core_Power        :: Core Thermal Power                  [W]
      ! 2)  Core_Power_0      :: Initial Core_Power
      ! 3)  Core_Power_100    :: 100% Core Thermal Power             [W]
      ! 4)  PPower            :: (Core_Power/Core_Power_100)         [%]
      ! 5)  PPower_0          :: Initial PPower
      ! 6)  PPower_Old        :: Previous Burnup/Time Step PPower
      ! 7)  Avg_CorePower     :: Average Power for Fuel Region       [W/cm^3]
      ! 8)  TM_In             :: Core Inlet   Moderator Temperature  [C]
      ! 9)  TM_Out            :: Core Outlet  Moderator Temperature  [C]
      ! 10) TF_In             :: Core Initial Fuel      Temperature  [C]
      ! 11) Core_Pressure     :: Core Pressure                       [bar]
      ! 12) Core_MassFlow     :: Core Total Mass Flow                [kg/sec]
      ! 13) FA_Power          :: Fuel Assembly Thermal Power         [W]
      ! 14) FA_Power_0        :: Initial FA_Power
      ! 15) FA_Power_100      :: 100% Fuel Assembly Thermal Power    [W]
      ! 16) FA_MassFlow       :: Fuel Assembly Mass Flow             [kg/sec]
      ! 17) PPM               :: Boron Concentration                 [ppm]
      ! 18) T_Mod_Out         :: Outlet Moderator Temperature        [K] [2D]
      ! 19) Frac_Gamma        :: Fraction of Heat Deposited Directly in Coolant (Direct Heating Fraction)
      ! 20) h_Gap             :: Gap Conductance                     [J/m^2-C]
      ! 21) LevFactor         :: Power/Flux Level Balancing Factor for Design Power
      ! 22) LinPOW            :: Linear Power Density                [W/m]
      !     (Ixy, Iz)
      implicit none

      real(8) :: Core_Power
      real(8) :: Core_Power_0
      real(8) :: Core_Power_100
      real(8) :: PPower
      real(8) :: PPower_0
      real(8) :: PPower_Old
      real(8) :: Avg_CorePower
      real(8) :: TM_In
      real(8) :: TM_Out
      real(8) :: Core_Pressure=155d0 ! bar
      real(8) :: Core_MassFlow
      real(8) :: inp_MassFlow
      real(8) :: FA_Power
      real(8) :: FA_Power_0
      real(8) :: FA_Power_100
      real(8) :: FA_MassFlow
      real(8) :: PPM
      real(8) :: Frac_Gamma
      real(8) :: h_Gap
      real(8) :: LevFactor
      real(8) :: TF_In
      real(8) :: Avg_CorePower0
      real(8), dimension(:), allocatable :: T_Mod_Out
      real(8), dimension(:,:), allocatable :: LinPOW
#ifdef tuan_fr
      real(8) :: TC_In       ! initial clad temperature from input
#endif

      real(8) :: LinPOW_ave
      real(8) :: LinPOW_max

      logical(1) :: opt_fenthl=.false.
      real(8) :: max_fenthl=0d0
      real(8) :: max_fenthl_rise=0d0

      logical(4) :: OPT_findtavg=.false.
      logical(4) :: OPT_findtavg_init=.false.
      real(8) :: target_tavg=585.0d0
      real(8) :: itr_tavg(1:3)
      real(8) :: itr_fa_flowrate(3)
      real(8) :: chanvf_save=0.0d0

      logical(1) :: flag_th_chanwise=.false.
      integer, allocatable :: i_chan_1n(:)
      integer :: max_i_chan
      integer :: min_i_chan
      real(8), allocatable :: chanwise_t_inlet(:)
      real(8), allocatable :: chanwise_g_inlet(:)

      logical(1) :: opt_tfutab=.false.
      real(8) :: table_bu(10)
      real(8) :: table_pow(5)
      real(8) :: table_tfu(5,10)
      real(8) :: segtfu(3)
      data table_bu(1:10) &
      / 0.00d+0, 5.00d+0, 1.00d+1, 1.50d+1, 2.00d+1, &
      & 3.00d+1, 4.00d+1, 5.00d+1, 6.50d+1, 8.00d+1 /
      data table_pow(1:5) &
      / 5.00d-1, 1.00d+0, 1.50d+0, 2.00d+0, 2.50d+0 /
      data table_tfu(1:5,1:10) &
      / 295.6d0, 291.2d0, 288.2d0, 285.6d0, 283.3d0, &
      & 282.0d0, 280.7d0, 280.9d0, 282.2d0, 284.0d0, &
      & 269.6d0, 271.4d0, 273.8d0, 277.2d0, 281.2d0, &
      & 260.6d0, 264.3d0, 268.0d0, 272.4d0, 277.0d0, &
      & 253.9d0, 258.7d0, 262.8d0, 267.3d0, 272.5d0, &
      & 243.9d0, 252.5d0, 263.0d0, 273.9d0, 284.9d0, &
      & 261.2d0, 272.3d0, 282.0d0, 292.3d0, 302.2d0, &
      & 280.9d0, 291.6d0, 300.8d0, 310.6d0, 319.6d0, &
      & 314.2d0, 324.6d0, 333.6d0, 343.2d0, 350.8d0, &
      & 355.7d0, 366.1d0, 375.7d0, 385.4d0, 390.7d0 /
      data segtfu(1:3) &
      / 0.000000d+00, 0.000000d+00, 1.751609d+01 /

      logical(1) :: flag_cal_dnbr=.false.
      real(8), allocatable :: qchf(:,:)
      real(8), allocatable :: qwat(:,:)
      real(8), allocatable :: dnbr(:,:)
      real(8) :: mqchf
      real(8) :: mqwat
      real(8) :: mdnbr

      ! th geom
      integer:: nchan, nzth
      integer:: nr, nrp1, nrp2, nrp3, nrp4, nrp5
      real(8):: acf                                  !coolant flow area
      real(8):: afp                                  !fuel pellet area
      real(8):: xi                                   !wetted preimeter
      real(8):: zeta                                 !heated perimeter
      real(8):: zetap                                !heated perimeter density
      real(8):: deq                                  !equivalent diameter
      real(8):: delr                                 !radial mesh spacing in the pellet region
      real(8):: delrw                                !radial mesh spacing in the cladding region
      real(8):: tworm                                !tw over rm(=rg+0.5*tw)
      real(8):: delr2                                !delr^2
      real(8):: delrw2                               !delrw^2
      integer, dimension(:), allocatable :: ltochan  !neutronic node number to channel number
      integer, dimension(:), allocatable :: lchanptr !ch no. to neut. node no.
      integer, dimension(:), allocatable :: lchantol !chan. no. to neut. node n.
      real(8), dimension(:), allocatable :: r        !radial mesh coordinate

      ! th data
      real(8) :: wfcl
      real(8) :: wfsurf
      real(8) :: akfuel(0:5)
      real(8) :: akclad(0:3)
      real(8) :: arcpfuel(0:3)
      real(8) :: arcpclad(0:3)
      real(8) :: tdoplmax
      real(8), dimension(:,:,:), allocatable :: tfuel
      real(8), dimension(:,:), allocatable :: tdopl
      real(8), dimension(:,:), allocatable :: tcool
      real(8), dimension(:,:), allocatable :: dcool
      real(8), dimension(:,:), allocatable :: hcool

      real(8), dimension(:,:), allocatable :: qflux
      real(8), dimension(:,:), allocatable :: qvol
      real(8), dimension(:,:), allocatable :: qeff
      real(8), dimension(:,:), allocatable :: htcoef
      real(8), dimension(:,:), allocatable :: rhou
      real(8), dimension(:,:), allocatable :: rhohu
      real(8), dimension(:,:), allocatable :: u
      real(8), dimension(:,:), allocatable :: ud

#ifdef js_mpc
      real(8), dimension(:,:,:), allocatable :: tfuel_bk
      real(8), dimension(:,:), allocatable :: tdopl_bk
      real(8), dimension(:,:), allocatable :: tcool_bk
      real(8), dimension(:,:), allocatable :: dcool_bk
      real(8), dimension(:,:), allocatable :: hcool_bk
#endif

      ! th option
      real(8) :: fracdf      !fraction of heat in fuel
      real(8) :: din         !inlet density
      real(8) :: tdopin      !inlet doppler temperature
      real(8) :: toutavg     !average outlet temperature
      real(8) :: tfmax       !maximum fuel centerline temperature
      integer :: fluid_type  ! type of the fluid (water or heavy water =>>> fluid_type = 0 water, =1 heavy water)
#ifdef tuan_fr
      REAL(8) :: t_in             ! fuel temperature from input
      integer :: chole_type       ! type of the chole
      integer :: fuel_type        ! type of the fuel
      integer :: gap_type         ! type of the gap
      integer :: clad_type        ! type of the cladding
      integer :: coolant_type     ! type of the coolant

#endif
#ifdef tuan_tr_test
      real(8),            allocatable   :: tdopl0                     (:,:)
      real(8),            allocatable   :: tcool0                     (:,:)
      real(8),            allocatable   :: dcool0                     (:,:)
#endif


      logical(1) :: flag_th_init_first=.false.

#ifdef hj_timectrl
      type :: tfcal_type
         real(8),allocatable,dimension(:) :: kf,kfb,kfm,kfmb
         real(8),allocatable,dimension(:) :: x,rhocp,ad,al,au,b
      endtype tfcal_type
      type(tfcal_type) :: trtf2
#endif

      END MODULE Inc_TH
