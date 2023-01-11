      MODULE Mod_Manage

      USE Inc_Constant
      USE Inc_Geometry
      USE Inc_RP
      USE Inc_Option

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE Initializing

      USE Inc_Control
      USE Inc_XYZ, only: ASI
      USE Inc_CR
      USE Inc_Depletion
      USE Inc_Detector
      USE Inc_DF
      USE Inc_FA
      USE Inc_File
      USE Inc_Flag
      USE Inc_INP
      USE Inc_Kinetics, ONLY: N_Group_d
      USE Inc_3d, only: avg_power
      USE Inc_maXS, ONLY: nu, kappa
      USE Inc_3d, only: keff, keff_old, keff_oold
      USE Inc_Nuclide, ONLY: N_Heavy, N_NdSm, N_Gd
      USE Inc_PinPOW
      USE Inc_RST
      USE Inc_TH
      USE Inc_Time
      USE Inc_Transient
      USE Inc_Utile
      USE Inc_XS_File
      use inc_option, only: opt_jumpin
      USE Mod_Alloc
      use mod_charedit, only: print_msg
      IMPLICIT NONE
      INTEGER ::  i


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Initializing] in Mod_Manage'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif
#ifdef siarhei_fr
      call print_msg(0,"=============================================================================")
      call print_msg(0,"| oooooooooo.        .o.        .ooooooo.  ooooooooooooo      ooooooooooooo |")
      call print_msg(0,"| `888   `Y88.      .888.      d8P'    `YP 8'   888   `8      `88888888888' |")
      call print_msg(0,"|  888   .d88'     .88888.     Y88bo.           888            888          |")
      call print_msg(0,"|  888ooo88P'     .8' `888.     `Y88888o.       888      oooo  888888888P   |")
      call print_msg(0,"|  888`88b.      .88ooo8888.        `Y888b      888            888      '   |")
      call print_msg(0,"|  888  `88b.   .8'     `888.  db.    .d8P      888            888          |")
      call print_msg(0,"| o888o  o888o o88o     o8888o `88888888P'     o888o          o888o         |")
      call print_msg(0,"=============================================================================")
#else
      call print_msg(0,"=============================================================================")
      call print_msg(0,"| oooooooooo.        .o.        .ooooooo.  ooooooooooooo      ooooo   ooooo |")
      call print_msg(0,"| `888   `Y88.      .888.      d8P'    `YP 8'   888   `8      `888   .d8P'  |")
      call print_msg(0,"|  888   .d88'     .88888.     Y88bo.           888            888o.od8P'   |")
      call print_msg(0,"|  888ooo88P'     .8' `888.     `Y88888o.       888      oooo  888`Yo88     |")
      call print_msg(0,"|  888`88b.      .88ooo8888.        `Y888b      888            888  `Y8b.   |")
      call print_msg(0,"|  888  `88b.   .8'     `888.  db.    .d8P      888            888   `Y8b.  |")
      call print_msg(0,"| o888o  o888o o88o     o8888o `88888888P'     o888o          o888o   o888o |")
      call print_msg(0,"=============================================================================")
#endif

#ifdef siarhei_fr
      call print_msg(0,"                                                                      v.0.9.0")
#else
      call print_msg(0,"                                                                      v.2.2.0")
#endif
      call print_msg(0,"")

      EPS_keff_St   = DM6
      EPS_Linf_St   = DM6
      EPS_Linf_Tr   = DM5
      EPS_PowLev    = DM5!5.D-5
      EPS_Flux      = DM6
      EPS_TF        = DM3
      EPS_TM        = DM3
      EPS_Xe        = DM4
      EPS_Sm        = DM4

#ifdef js_dbg
      EPS_Critical  = 1d-4
      EPS_Global    = 1d-3
      EPS_Local     = 5.D-2
#else
      EPS_Critical  = 1d-5
      EPS_Global    = 1d-5
      EPS_Local     = 5.D-4
#endif
      EPS_Residual  = DM3
      EPS_ERF       = 5.D-3
      EPS_DoppTF    = DM3

      Flag_Conv_Critical = .FALSE.
      Flag_Conv_PowLev   = .FALSE.

      !Iout_Max   = 31
      Iout_Max   = 500
      !Iin_Max    = 1
      Iin_Max    = 10
      I_CriSearch_Max = 20

      I_CRCH_Max = 500

      I_XSFB_Max = 1

      Period_TH   = 3
      Period_NonL = 3
      Period_CMFD = 3

      ini_keff = D1
      ini_Flux = D1
      ini_jOut = D1
      ini_TL   = D0

      ASI = D0

      N_CR = 0

      Speed_CR = D0
      OverLap  = D0

      Flag_MoveCR = .FALSE.

      Flag_Card_CR = .FALSE.

      N_BU = 1

      Flag_miDepl = .FALSE.
      Flag_10ppm  = .FALSE.

      Flag_Card_Depletion = .FALSE.

      dBU     = D0
      RST_BU  = -1d0
      Tot_MTU = D0

      Flag_PC = .FALSE.

      alpha_0 = 2.124853710495237488D-16

      Order_CRAM = 16

      c_Fuzzy        = D1
      nu_per_kap_avg = D1

      h_FA = 21.606D0

      N_Pin     = 264.D0
      N_GT      = 25.D0
      N_GT_NoCR = 1.D0

      R_Fuel     = 0.41195D0        * DM2  ! [cm] => [m]
      R_Gap      = 0.41875D0        * DM2  ! [cm] => [m]
      R_Clad     = 0.47585D0        * DM2  ! [cm] => [m]
      h_Clad     = (R_Clad - R_Gap) * DM2  ! [cm] => [m]
      R_GT_CR    = 0.42175D0        * DM2  ! [cm] => [m]
      R_GT_IClad = 0.57240D0        * DM2  ! [cm] => [m]
      R_GT_OClad = 0.61295D0        * DM2  ! [cm] => [m]

      VF_H2O_NoCR = 0.59D0
      VF_H2O_InCR = 0.56D0

      Len_INP = 1

      Name_INP  = "RAST-K.INP "
      Name_XS   = "RAST-K.XS  "
      Name_RI   = "RAST-K.RI  "
      Name_QS   = "RAST-K.QS  "
      Name_DRFB = "RAST-K.DRFB"
      Name_DRFM = "RAST-K.DRFM"
      Name_DRFT = "RAST-K.DRFT"
      Name_OUT  = "RAST-K.OUT "
      Name_SUM  = "RAST-K.SUM "
      Name_RST  = "RAST-K.RST "
      Name_POW  = "RAST-K.POW "
      Name_FLUX = "RAST-K.FLUX"
      Name_TH   = "RAST-K.TH  "
      Name_maXS = "RAST-K.maXS"
      name_adj  = "RAST-K.adj"
      name_pkd  = "RAST-K.pkd"
      name_rho  = "RAST-K.rho"

      INP_Title = "RAST-K"

      Flag_XS   = .FALSE.
      Flag_RI   = .FALSE.
      Flag_QS   = .FALSE.
      Flag_DRFB = .FALSE.
      Flag_DRFM = .FALSE.
      Flag_DRFT = .FALSE.
      Flag_OUT  = .FALSE.
      Flag_ANC  = .FALSE.
      Flag_SUM  = .FALSE.
      Flag_RST  = .FALSE.
#ifdef tuan_fr_rst
      Flag_RST_HEX  = .FALSE.
#endif

      Flag_POW  = .FALSE.
      Flag_FLUX = .FALSE.
      Flag_TH   = .FALSE.
      Flag_maXS = .FALSE.
      Flag_OCEAN = .FALSE.
      Flag_VIPRE = .FALSE.
      Flag_MAP = .FALSE.

      call alloc(flag_out_print,25)
      flag_out_print(1:18)=.true.
      flag_out_print(19:20)=.false.
      flag_out_print(21:25)=.false.

      Flag_XSFB  = .TRUE.
      !js+Flag_XSFB  = .FALSE.
      !js+Flag_THFB  = .TRUE.
      Flag_THFB  = .FALSE.
      Flag_InDET = .FALSE.
      Flag_ExDET = .FALSE.

      Nx_1N = 1
      Ny_1N = 1
      Nz    = 1

      BC_Rx = 1
      BC_Lx = 1
      BC_Ry = 1
      BC_Ly = 1
      BC_Rz = 1
      BC_Lz = 1

      OPT_Core = 4

      Nx        = 1
      Ny        = 1
      Nxy       = 1
      Nxy_1N    = 1
      Nxyz      = 1

      IzFuelBot = 1
      IzFuelTop = 1

      Core_Height = D1
      Tot_Vol     = D0
      Tot_FuelVol = D0

      Flag_4N1FA      = .FALSE.
      Flag_HalfCenter = .FALSE.
      Flag_Card_Geometry = .FALSE.
      Flag_Card_Title = .FALSE.
      Flag_Card_File = .FALSE.
      Flag_Card_Option = .FALSE.
      Flag_Card_Control = .FALSE.
      Flag_Card_Geometry = .FALSE.
      Flag_Card_RP = .FALSE.
      Flag_Card_Shuffling = .FALSE.
      Flag_Card_Rotation = .FALSE.
      Flag_Card_CR = .FALSE.
      Flag_Card_FA_Data = .FALSE.
      Flag_Card_TH_Data = .FALSE.
      Flag_Card_ATF_TCD = .FALSE.
      Flag_Card_Depletion = .FALSE.
      Flag_Card_maXS = .FALSE.
      Flag_Card_Transient = .FALSE.
      Flag_Card_Branch = .FALSE.
      Flag_RefDF = .FALSE.

      N_Group_d = 6

      N_Heavy = 22
      N_NdSm  = 10
      N_Gd    = 5

      N_LP      = 0
      Nxy_FA    = 0
      Nxy_RF    = 0
      Nxy_FA_1N = 0
      Nxy_RF_1N = 0

      Flag_Card_RP = .FALSE.

      Avg_Power = D0

      nu    = D0
      kappa = D0

      !b    = D0
      !c    = D0
      !s    = D0

      keff      = D1
      keff_Old  = D1
      keff_OOld = D1

      OPT_Mode    = 1
      OPT_RST     = 1
      OPT_Nodal   = 2
      OPT_Xe      = 0
      OPT_Sm      = 0
      OPT_SmChain = 1
      OPT_TLA     = 1
      OPT_TFcal   = 1

      N_Group = 2

      Flag_PinPOW  = .FALSE.
      Flag_NPinOdd = .FALSE.

      OPT_jumpin=.false.

      TM_In_RST0 = D0

      N_Branch = 0
      I_Branch = 0

      Core_Power     = 2775.D6
      Core_Power_0   = Core_Power
      Core_Power_100 = Core_Power
      PPower         = D1
      PPower_0       = PPower
      PPower_Old     = PPower
      Avg_CorePower  = D1
      TM_In          = D0
      TM_Out         = TM_In
      TF_In          = TM_In
      Core_Pressure  = 155.D0
      Core_MassFlow  = 12893.D0
      FA_Power       = 17.67516D0
      FA_Power_0     = FA_Power
      FA_Power_100   = FA_Power
      FA_MassFlow    = 82.12102D0
      PPM            = D0
      Frac_Gamma     = 0.019D0
      h_Gap          = D4
      LevFactor      = D1

      Flag_Card_TH_Data = .FALSE.

      Sec      = D0
      Hour     = D0
      Day      = D0
      dT       = D0
      dT_Tr    = dT
      dT_Old   = dT
      T_Tot    = 5.D0
      T_Switch = D1
      T_Expand = DP1

      I_Time = 1
      N_Time = 1

      OPT_Xe_inTr      = 0
      OPT_Sm_inTr      = 0
      OPT_SmChain_inTr = 1

      Flag_MovePPM    = .FALSE.
      Flag_DepBP_inTr = .FALSE.

      Flag_Transient = .FALSE.

      L_Omega = D0

      N_Dummy = 100

      CALL Alloc( Dum_Real , N_Dummy )
      CALL Alloc( Dum_Int  , N_Dummy )
      CALL Alloc( Dum_Char , N_Dummy )
      CALL Alloc( Dum_Char2 , N_Dummy )
      Dum_Real = - D1
      Dum_Int  = - 1
      Dum_Char = ''

      Buff_Line = ''
      Buff_MSG  = ''
      Buff_MSG2 = ''

      DO i = 1, 1024
         Div_Line  (i:i) = '-'
         Div_Line2 (i:i) = '='
         Div_Line3 (i:i) = '~'
      END DO

      Blank_Line = ''

      I_Tab  = 1
      I_Type = 1
      I_SP   = 1
      I_BU   = 1
      I_Reg  = 1

      N_XS_Table   = 1
      N_SP         = 16
      N_SP_FA      = 16
      N_SP_RF      = 8
      N_BU         = 1
      N_BU_Brch    = 1
      Buff_N_BU    = N_BU
      XS_File_NbyN = 1

      N_BU_Ref_Max   = 0
      N_BU_Ref_NoBP  = 42
      N_BU_Brch_NoBP = 14
      N_BU_Brch_wBP  = 15
      N_BU_SDC       = 64

      N_Brch    = 0
      N_Region  = 1
      N_XS_Kind = 12

      Flag_ADF = .FALSE.

      RETURN
      END SUBROUTINE Initializing


!!!#ifdef siarhei_delete 
      SUBROUTINE Set_Memory
      USE Inc_Solver, ONLY: nkincomp
      USE Inc_3D
      USE Inc_CR
      USE Inc_Depletion
      USE Inc_DF
      USE Inc_File
      USE Inc_Flag
      USE Inc_History
      USE Inc_INP
      USE Inc_Kinetics
      USE Inc_maXS
      USE Inc_miXS
      USE Inc_Nuclide
      USE Inc_Detector
      use Inc_Crud
      USE Inc_RST
      USE Inc_TH
      USE Inc_Time
      USE Inc_Transient
      USE Inc_Utile
      USE Inc_XS_File
      USE Inc_XYZ
      USE Mod_Alloc
      USE Inc_Option
#ifdef siarhei_tr_hex
      use Inc_Control, only: iout_Max ! siarhei_check
#endif

#ifdef js_r2mpi
      use inc_parallel, only: comm
#endif
#ifdef JR_SRCTRM
      use inc_snf
      use inc_pinvar, only: npin
#endif

#ifdef tuan_fr
      use inc_tpen
#endif
      IMPLICIT NONE
      REAL(8) :: Buff_REAL
      INTEGER :: Buff_INT
      INTEGER :: i

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Set_Memory] in Mod_Manage'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif



#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      Buff_REAL = D0
      Buff_INT = 0
      i        = 0

      Allocate( FisSrc      (1: Nxy , 1:Nz ))
      Allocate( FisSrc_Old  (1: Nxy , 1:Nz ))
      Allocate( FisSrc_OOld (1: Nxy , 1:Nz ))
      Allocate( FisSrc_Iout (1: Nxy , 1:Nz ))
      Allocate( Power       (1: Nxy , 1:Nz ))
      Allocate( Normal_Power(1: Nxy , 1:Nz ))

      IF ( .not. Flag_RI ) THEN
         Allocate( BU (1: Nxy , 1:Nz ))
         Allocate( BU_Old(1: Nxy , 1:Nz ))
      END IF
        BU = 0
        BU_Old = 0
        
      Allocate(BU_Predictor    (1:Nxy , 1:Nz ))
      Allocate(T_Fuel          (1:Nxy , 1:Nz ))
      Allocate(T_Fuel_Old_XSFB (1:Nxy , 1:Nz ))
      Allocate(T_Mod           (1:Nxy , 1:Nz ))
      Allocate(T_Mod_Old_XSFB  (1:Nxy , 1:Nz ))
      Allocate(D_Mod           (1:Nxy , 1:Nz ))
      Allocate(D_Mod_Old_XSFB  (1:Nxy , 1:Nz ))

      Allocate( Flux (0: Nxy + 1 , 1: Nz , 1: N_Group ))
      !CALL Alloc( Flux          , Nxy , Nz , N_Group )
      Allocate( Flux_Old      ( 1:Nxy , 1:Nz , 1:N_Group ))
      Allocate( Flux_OOld     ( 1:Nxy , 1:Nz , 1:N_Group ))
      Allocate( Flux_Old_XSFB ( 1:Nxy , 1:Nz , 1:N_Group ))
      Allocate(flux_adj(1:nxy,1:nz,1:n_group))

      Allocate(Reg_VF ( 1:N_Region ))

      IF ( N_CR /= 0 ) THEN
         Allocate( Flag_Move (1: N_CR ))
      END IF

      allocate( Mat_I_HN( N_Heavy, N_Heavy ) )
      Mat_I_HN = (0d0,0d0)
      do i=1,N_Heavy
         Mat_I_HN(i,i)=(1d0,0d0)
      enddo

      allocate( Mat_I_FP( N_NdSm , N_NdSm ) )
      Mat_I_FP = (0d0,0d0)
      do i=1,N_NdSm
         Mat_I_FP(i,i)=(1d0,0d0)
      enddo

      allocate( alpha ( Order_CRAM / 2 ) )
      alpha = ( 0d0 , 0d0 )
      alpha(1) = ( - 5.0901521865224915650D-7 , - 2.4220017652852287970D-5 )
      alpha(2) = ( + 2.1151742182466030907D-4 , + 4.3892969647380673918D-3 )
      alpha(3) = ( + 1.1339775178483930527D+2 , + 1.0194721704215856450D+2 )
      alpha(4) = ( + 1.5059585270023467528D+1 , - 5.7514052776421819979D+0 )
      alpha(5) = ( - 6.4500878025539646595D+1 , - 2.2459440762652096056D+2 )
      alpha(6) = ( - 1.4793007113557999718D+0 , + 1.7686588323782937906D+0 )
      alpha(7) = ( - 6.2518392463207918892D+1 , - 1.1190391094283228480D+1 )
      alpha(8) = ( + 4.1023136835410021273D-2 , - 1.5743466173455468191D-1 )

      allocate( theta ( Order_CRAM / 2 ) )
      theta = ( 0d0 , 0d0 )
      theta(1) = ( - 1.0843917078696988026D+1 , + 1.9277446167181652284D+1 )
      theta(2) = ( - 5.2649713434426468895D+0 , + 1.6220221473167927305D+1 )
      theta(3) = ( + 5.9481522689511774808D+0 , + 3.5874573620183222829D+0 )
      theta(4) = ( + 3.5091036084149180974D+0 , + 8.4361989858843750826D+0 )
      theta(5) = ( + 6.4161776990994341923D+0 , + 1.1941223933701386874D+0 )
      theta(6) = ( + 1.4193758971856659786D+0 , + 1.0925363484496722585D+1 )
      theta(7) = ( + 4.9931747377179963991D+0 , + 5.9968817136039422260D+0 )
      theta(8) = ( - 1.4139284624888862114D+0 , + 1.3497725698892745389D+1 )

      Allocate(ADF_Lx       (1:Nxy, 1: Nz ,1: N_Group ))
      Allocate(ADF_Rx       (1:Nxy, 1: Nz ,1: N_Group ))
      Allocate(ADF_Ly       (1:Nxy, 1: Nz ,1: N_Group ))
      Allocate(ADF_Ry       (1:Nxy, 1: Nz ,1: N_Group ))
      Allocate(ADF_Lz       (1:Nxy, 1: Nz ,1: N_Group ))
      Allocate(ADF_Rz       (1:Nxy, 1: Nz ,1: N_Group ))
      Allocate(ADF_Avg      (1:Nxy, 1: Nz ,1: N_Group ))
      Allocate(Buff_ADF_Lx  (1:Nxy, 1: Nz ,1: N_Group ))
      Allocate(Buff_ADF_Rx  (1:Nxy, 1: Nz ,1: N_Group ))
      Allocate(Buff_ADF_Ly  (1:Nxy, 1: Nz ,1: N_Group ))
      Allocate(Buff_ADF_Ry  (1:Nxy, 1: Nz ,1: N_Group ))
      ADF_Lx      = D1
      ADF_Rx      = D1
      ADF_Ly      = D1
      ADF_Ry      = D1
      ADF_Lz      = D1
      ADF_Rz      = D1
      ADF_Avg     = D1
      Buff_ADF_Lx = D1
      Buff_ADF_Rx = D1
      Buff_ADF_Ly = D1
      Buff_ADF_Ry = D1

      Allocate( CDF_LxLy (1: Nxy , 1:Nz , 1:N_Group ))
      Allocate( CDF_LxRy (1: Nxy , 1:Nz , 1:N_Group ))
      Allocate( CDF_RxLy (1: Nxy , 1:Nz , 1:N_Group ))
      Allocate( CDF_RxRy (1: Nxy , 1:Nz , 1:N_Group ))
      CDF_LxLy = D1
      CDF_LxRy = D1
      CDF_RxLy = D1
      CDF_RxRy = D1

      IF ( OPT_Mode == 1 ) THEN
         Buff_INT = N_Branch + 1
      ELSE IF ( OPT_Mode == 3 ) THEN
         Buff_INT = N_Time

         ! Due to SIZE(Clcle_BU) = N_BU = 1 in Transient Case
          !  dummy_filler = 1 ! @$^ siarhei_plot 
         ! There is Error to Save Hist_Cycle_BU
         Buff_REAL = Cycle_BU(1)

         DEALLOCATE( Cycle_BU )
         Allocate( Cycle_BU  (1:Buff_INT ))
         Cycle_BU = Buff_REAL
      ELSE
         Buff_INT = N_BU + 30
      END IF

      IF ( ALLOCATED( Cycle_BU ) .EQV. .FALSE. ) THEN
         Allocate( Cycle_BU ( 1:Buff_INT ))
      END IF

      Allocate( Hist_Day  (1: Buff_INT ))

      IF ( ALLOCATED( Hist_Sec ) .EQV. .FALSE. ) THEN
         Allocate( Hist_Sec (1: Buff_INT ))
      END IF

      IF ( ALLOCATED( Hist_Hour ) .EQV. .FALSE. ) THEN
         Allocate( Hist_Hour (1: Buff_INT ))
      END IF

      Allocate( Hist_Time ( 1:Buff_INT ))

      Allocate( Hist_Cycle_BU   (1: Buff_INT ))
      Allocate( Hist_Accum_BU   (1: Buff_INT ))
      Allocate( Hist_keff       (1: Buff_INT ))
      Allocate( Hist_Reactivity (1: Buff_INT ))

      if (opt_mode==3) then
         opt_fenthl=.true.
         Allocate( Hist_MaxFenthl(1: Buff_INT ))
         Allocate( Hist_MaxFenthl_Rise(1: Buff_INT ))
      endif

      IF ( ALLOCATED( Hist_PPower ) .EQV. .FALSE. ) THEN
         Allocate( Hist_PPower (1: Buff_INT ))
      END IF

      IF ( ALLOCATED( Hist_PPM ) .EQV. .FALSE. ) THEN
         Allocate( Hist_PPM    (1: Buff_INT ))
         Allocate( Hist_ADJPPM (1: Buff_INT ))
         Allocate( Hist_PPMB10 (1: Buff_INT ))
         Allocate( Hist_B10DEP (1: Buff_INT ))
      END IF

      Allocate( Hist_ASI      (1: Buff_INT ))
      Allocate( Hist_PF_1D    (1: Buff_INT ))
      Allocate( Hist_PP_1D    (1: Buff_INT ))
      Allocate( Hist_PF_2D    (1: Buff_INT ))
      Allocate( Hist_PP_2D    (1: Buff_INT ))
      Allocate( Hist_PF_3D    (1: Buff_INT ))
      Allocate( Hist_PP_3D    (1: Buff_INT ))
      Allocate( Hist_PBU_FA   (1: Buff_INT ))
      Allocate( Hist_PPBU_FA  (1: Buff_INT ))
      Allocate( Hist_PBU_Pin  (1: Buff_INT ))
      Allocate( Hist_PPBU_Pin (1: Buff_INT ))

      if (flag_ocean) then
         Allocate( Hist_Fxy      (1: Buff_INT ))
         Allocate( Hist_Fq       (1: Buff_INT ))
         Allocate( Hist_Fr       (1: Buff_INT ))
      endif

      Allocate( Hist_Fq_val       (1: Buff_INT ))
      Allocate( Hist_Fr_val       (1: Buff_INT ))
      Allocate( Hist_Fz_val       (1: Buff_INT ))
      Allocate( Hist_FdH_val      (1: Buff_INT ))
      Allocate( Hist_max_Fxy_val  (1: Buff_INT ))

      Allocate( HIST_LinPOW_ave   (1: Buff_INT ))
      Allocate( HIST_LinPOW_max   (1: Buff_INT ))
      Allocate( Hist_Fz_ave       (1: Buff_INT ))
      Allocate( Hist_Fxy_val      (1: Buff_INT , 1:Nx_1N, 1:Ny_1N))
      Allocate( Hist_FdH_val_XY   (1: Buff_INT , 1:Nx_1N, 1:Ny_1N))
      Allocate( FdH_val_XY        (1: Nx_1N,1: Ny_1N))
      Allocate( Hist_kinf         (1: Buff_INT ))
      Allocate( HIST_a_1          (1: Buff_INT ))
      Allocate( HIST_a_2          (1: Buff_INT ))
      Allocate( HIST_D_1          (1: Buff_INT ))
      Allocate( HIST_D_2          (1: Buff_INT ))
      Allocate( HIST_nf_1         (1: Buff_INT ))
      Allocate( HIST_nf_2         (1: Buff_INT ))
      Allocate( HIST_r_1          (1: Buff_INT ))
      Allocate( HIST_kf_1         (1: Buff_INT ))
      Allocate( HIST_kf_2         (1: Buff_INT ))
      Allocate( HIST_f_1          (1: Buff_INT ))
      Allocate( HIST_f_2          (1: Buff_INT ))
      Allocate( Hist_FA_kinf      (1: Buff_INT , 1:Nxy_1N))
      Allocate( Hist_FA_kinf_XY   (1: Buff_INT , 1:Nx_1N,1: Ny_1N))
      Allocate( Hist_FA_kinf_4N   (1: Buff_INT , 1:Nxy))
      Allocate( Hist_FA_kinf_XY_4N(1: Buff_INT , 1:Nx, 1:Ny))

      Allocate( Hist_Fq_loc       (1: Buff_INT , 1:5 ))
      Allocate( Hist_Fr_loc       (1: Buff_INT , 1:4 ))
      Allocate( Hist_Fz_loc       (1: Buff_INT , 1:1 ))
      Allocate( Hist_FdH_loc      (1: Buff_INT , 1:4 ))
      Allocate( Hist_max_Fxy_loc  (1: Buff_INT , 1:3 ))

      Allocate( Hist_maxPBU_2D    (1:Buff_INT ))
      Allocate( Hist_maxPPD_2D    (1:Buff_INT ))
      Allocate( Hist_maxPPD_3D    (1:Buff_INT ))
      Allocate( Hist_maxPPW_2D    (1:Buff_INT ))
      Allocate( Hist_maxPPW_3D    (1:Buff_INT ))
      Allocate( Hist_PBU_2D_loc   (1:Buff_INT, 1:4))
      Allocate( Hist_PPD_2D_loc   (1:Buff_INT, 1:4))
      Allocate( Hist_PPD_3D_loc   (1:Buff_INT, 1:5))
      Allocate( Hist_PPW_2D_loc   (1:Buff_INT, 1:4))
      Allocate( Hist_PPW_3D_loc   (1:Buff_INT, 1:5))

      IF ( ALLOCATED( Hist_CR ) .EQV. .FALSE. ) THEN
         Allocate( Hist_CR (1: Buff_INT ,1: N_CR , 1:3 ))
      END IF

      Allocate( Hist_TF_Avg    (1: Buff_INT ))
      Allocate( Hist_TF_Max    (1: Buff_INT ))
      Allocate( Hist_TF_Min    (1: Buff_INT ))
      Allocate( Hist_TFcen_Max (1: Buff_INT ))
      Allocate( Hist_TM_Avg    (1: Buff_INT ))
      Allocate( Hist_TM_Max    (1: Buff_INT ))
      Allocate( Hist_TM_Min    (1: Buff_INT ))

      IF ( ALLOCATED( Hist_TM_In ) .EQV. .FALSE. ) THEN
         Allocate( Hist_TM_In (1: Buff_INT ))
      END IF

      Allocate( Hist_TM_Out    (1: Buff_INT ))
      Allocate( Hist_DM_Avg    (1: Buff_INT ))
      Allocate( Hist_DM_Max    (1: Buff_INT ))
      Allocate( Hist_DM_Min    (1: Buff_INT ))

      flag_cal_dnbr=.true.
      if (flag_cal_dnbr) then
          Allocate(hist_qchf (1: buff_int,1: nxy,1: nz) )
          Allocate(hist_qwat (1: buff_int,1: nxy,1: nz) )
          Allocate(hist_dnbr (1: buff_int,1: nxy,1: nz) )
          Allocate(hist_mqchf(1: buff_int))
          Allocate(hist_mqwat(1: buff_int))
          Allocate(hist_mdnbr(1: buff_int))
      endif

!      Allocate( Hist_Xe35     (1: Buff_INT ))
!      Allocate( Hist_Xe35_ASI (1: Buff_INT ))
!      Allocate( Hist_Sm49     (1: Buff_INT ))
!      Allocate( Hist_Sm49_ASI (1: Buff_INT ))
!      if (flag_out_print(5).or.flag_out_print(6).or.flag_out_print(21).or.flag_out_print(22)) then
!         Allocate( Hist_I35      (1: Buff_INT ))
!         Allocate( Hist_Nd47     (1: Buff_INT ))
!         Allocate( Hist_Nd48     (1: Buff_INT ))
!         Allocate( Hist_Nd49     (1: Buff_INT ))
!         Allocate( Hist_Pm47     (1: Buff_INT ))
!         Allocate( Hist_Ps48     (1: Buff_INT ))
!         Allocate( Hist_Pm48     (1: Buff_INT ))
!         Allocate( Hist_Pm49     (1: Buff_INT ))
!         Allocate( Hist_Sm47     (1: Buff_INT ))
!         Allocate( Hist_Sm48     (1: Buff_INT ))
!
!         Allocate( Hist_Gd52     (1:Buff_INT ))
!         Allocate( Hist_Gd54     (1:Buff_INT ))
!         Allocate( Hist_Gd55     (1:Buff_INT ))
!         Allocate( Hist_Gd56     (1:Buff_INT ))
!         Allocate( Hist_Gd57     (1:Buff_INT ))
!         Allocate( Hist_Gd58     (1:Buff_INT ))
!         Allocate( Hist_Gd60     (1:Buff_INT ))
!
!         Allocate( Hist_U34  (1: Buff_INT ))
!         Allocate( Hist_U35  (1: Buff_INT ))
!         Allocate( Hist_U36  (1: Buff_INT ))
!         Allocate( Hist_U37  (1: Buff_INT ))
!         Allocate( Hist_U38  (1: Buff_INT ))
!         Allocate( Hist_Np37 (1: Buff_INT ))
!         Allocate( Hist_Np38 (1: Buff_INT ))
!         Allocate( Hist_Np39 (1: Buff_INT ))
!         Allocate( Hist_Pu38 (1: Buff_INT ))
!         Allocate( Hist_Pu39 (1: Buff_INT ))
!         Allocate( Hist_Pu40 (1: Buff_INT ))
!         Allocate( Hist_Pu41 (1: Buff_INT ))
!         Allocate( Hist_Pu42 (1: Buff_INT ))
!         Allocate( Hist_Pu43 (1: Buff_INT ))
!         Allocate( Hist_Am41 (1: Buff_INT ))
!         Allocate( Hist_As42 (1: Buff_INT ))
!         Allocate( Hist_Am42 (1: Buff_INT ))
!         Allocate( Hist_Am43 (1: Buff_INT ))
!         Allocate( Hist_Am44 (1: Buff_INT ))
!         Allocate( Hist_Cm42 (1: Buff_INT ))
!         Allocate( Hist_Cm43 (1: Buff_INT ))
!         Allocate( Hist_Cm44 (1: Buff_INT ))
!      endif

      Allocate( Hist_Normal_Power_Z   (1: Buff_INT, 1:Nz ))
      Allocate( Hist_Normal_Power_XY  (1: Buff_INT, 1:Nx_1N, 1:Ny_1N ))
      Allocate( Hist_Normal_Power_XYZ (1: Buff_INT, 1:Nx_1N, 1:Ny_1N , 1:Nz ))

      Allocate( Hist_BU_Z   (1: Buff_INT ,1:Nz ))
      Allocate( Hist_BU_XY  (1: Buff_INT ,1:Nx_1N, 1:Ny_1N ))
      Allocate( Hist_BU_XYZ (1: Buff_INT ,1:Nx_1N, 1:Ny_1N , 1:Nz ))

      Allocate( Hist_T_Fuel_Z   (1: Buff_INT ,1: Nz ))
      Allocate( Hist_T_Fuel_XY  (1: Buff_INT ,1: Nx_1N, 1:Ny_1N ))
      Allocate( Hist_T_Fuel_XYZ (1: Buff_INT ,1: Nx_1N, 1:Ny_1N , 1:Nz ))
      Allocate( Hist_T_Mod_Z    (1: Buff_INT ,1: Nz ))
      Allocate( Hist_T_Mod_XY   (1: Buff_INT ,1: Nx_1N, 1:Ny_1N ))
      Allocate( Hist_T_Mod_XYZ  (1: Buff_INT ,1: Nx_1N, 1:Ny_1N ,1: Nz ))
      Allocate( Hist_D_Mod_Z    (1: Buff_INT ,1: Nz ))
      Allocate( Hist_D_Mod_XY   (1: Buff_INT ,1: Nx_1N, 1:Ny_1N ))
      Allocate( Hist_D_Mod_XYZ  (1: Buff_INT ,1: Nx_1N, 1:Ny_1N , 1:Nz ))

      if (flag_ocean) then
      Allocate( Hist_N_Xe35_XYZ  (1: Buff_INT , 1:Nx_1N, 1:Ny_1N , 1:Nz ))
      Allocate( Hist_N_Sm49_XYZ  (1: Buff_INT , 1:Nx_1N, 1:Ny_1N , 1:Nz ))
      Allocate( Hist_TFlux_XYZ  (1: Buff_INT , 1:Nx_1N, 1:Ny_1N , 1:Nz ))
      Allocate( Hist_FFlux_XYZ  (1: Buff_INT , 1:Nx_1N, 1:Ny_1N , 1:Nz ))
      Allocate( Hist_Normal_Power_XYZ_4N( 1:Buff_INT, 1:Nx, Ny, 1:Nz ))
      Allocate( Hist_BU_XYZ_4N( 1:Buff_INT, 1:Nx, 1:Ny, 1:Nz ))
      Allocate( Hist_T_Fuel_XYZ_4N(1: Buff_INT,1: Nx, 1:Ny, 1:Nz ))
      Allocate( Hist_T_Mod_XYZ_4N (1: Buff_INT,1: Nx, 1:Ny, 1:Nz ))
      Allocate( Hist_N_Xe35_XYZ_4N(1: Buff_INT,1: Nx, 1:Ny, 1:Nz ))
      Allocate( Hist_N_Sm49_XYZ_4N(1: Buff_INT,1: Nx, 1:Ny, 1:Nz ))
      Allocate( Hist_FFlux_XYZ_4N (1: Buff_INT,1: Nx, 1:Ny, 1:Nz ))
      Allocate( Hist_TFlux_XYZ_4N (1: Buff_INT,1: Nx, 1:Ny, 1:Nz ))
      endif

      Allocate( Hist_Normal_Power_XY_4N( 1:Buff_INT, 1:Nx, 1:Ny ))
      Allocate( Hist_BU_XY_4N( 1:Buff_INT, 1:Nx, 1:Ny ))

      Allocate( FA_ASI(1: Nx_1N, 1:Ny_1N ))
      Allocate( Hist_FA_ASI(1: Nx_1N, 1:Ny_1N , 1:Buff_INT)) ! (x,y,step)

      Allocate( CR_mat_frac(1: Nxy,1: Nz ))
      CR_mat_frac = 1d0


      if (flag_print_kp) then
         if (opt_mode==2.or.opt_mode==4) then
            Allocate( hist_kp_beta    (1: Buff_INT,1: N_Group_d ))
            Allocate( hist_kp_lambda  (1: Buff_INT,1: N_Group_d ))
            Allocate( hist_kp_zeta    (1: Buff_INT,1: N_Group_d ))
            Allocate( hist_kp_gentime (1: Buff_INT ))
         else
            write(*,*) 'Print Kinetics Parameter Option only for ST or DP mode ...'
            stop
         endif
      endif

      Allocate( Hist_Power  (1: Buff_INT ,1: Nxy , 1:Nz ))
      Allocate( Hist_N_Xe35 (1: Buff_INT ,1: Nxy , 1:Nz ))
      Allocate( Hist_N_Sm49 (1: Buff_INT ,1: Nxy , 1:Nz ))

      if (flag_out_print(21)) then
         Allocate( Hist_Power_FA  (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Xe35_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Sm49_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_Burnup_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_T_Fuel_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_T_Mod_FA  (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_D_Mod_FA  (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_I35_FA  (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Nd47_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Nd48_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Nd49_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Pm47_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Ps48_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Pm48_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Pm49_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Sm47_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Sm48_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Gd52_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Gd54_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Gd55_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Gd56_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Gd57_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Gd58_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Gd60_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_U34_FA  (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_U35_FA  (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_U36_FA  (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_U37_FA  (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_U38_FA  (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Np37_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Np38_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Np39_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Pu38_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Pu39_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Pu40_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Pu41_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Pu42_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Pu43_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Am41_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_As42_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Am42_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Am43_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Am44_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Cm42_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Cm43_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
         Allocate( Hist_N_Cm44_FA (1: Buff_INT ,1: Nxy_1N , 1:Nz ))
      endif
      if (flag_out_print(22)) then
         Allocate( Hist_Burnup (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_T_Fuel (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_T_Mod  (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_D_Mod  (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_I35  (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Nd47 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Nd48 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Nd49 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Pm47 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Ps48 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Pm48 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Pm49 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Sm47 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Sm48 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Gd52 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Gd54 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Gd55 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Gd56 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Gd57 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Gd58 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Gd60 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_U34  (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_U35  (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_U36  (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_U37  (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_U38  (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Np37 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Np38 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Np39 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Pu38 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Pu39 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Pu40 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Pu41 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Pu42 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Pu43 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Am41 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_As42 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Am42 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Am43 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Am44 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Cm42 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Cm43 (1: Buff_INT ,1: Nxy , 1:Nz ))
         Allocate( Hist_N_Cm44 (1: Buff_INT ,1: Nxy , 1:Nz ))
      endif

      if (flag_out_print(23)) then
         Allocate( Hist_Flux          (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_maXS_tr_3D    (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_D_3D          (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_maXS_a_3D     (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_maXS_f_3D     (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_nu_maXS_f_3D  (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_kap_maXS_f_3D (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_maXS_s_3D     (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_maXS_r_3D     (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_ADF_Lx        (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_ADF_Rx        (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_ADF_Ly        (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_ADF_Ry        (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_CDF_LxLy      (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_CDF_LxRy      (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_CDF_RxLy      (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
         Allocate( Hist_CDF_RxRy      (1: Buff_INT, 1:Nxy, 1:Nz, 1:N_Group ))
      endif


      IF ( Flag_ExDET ) THEN
         Allocate( Hist_DR     ( 1:Buff_INT ,1: 3 ))
         Allocate( Hist_ASI_DR ( 1:Buff_INT     ))
      END IF

      if (flag_indet .and. flag_det_sig) then
         Allocate( Hist_det_sig_rr(1: BUff_INT, 1:DET_N_XY, 1:DET_N_Z ))
         if (opt_prt_detpow) then
            Allocate( Hist_det_detpow(1: BUff_INT, 1:DET_N_XY, 1:DET_N_Z ))
         endif
         if (opt_prt_wprime) then
            Allocate( Hist_det_wprime(1: BUff_INT, 1:DET_N_XY, 1:DET_N_Z ))
         endif
      endif

      if (flag_det_pow) then
         allocate(DET_A(Nxy*Nz*N_Group,Nxy*Nz*N_Group))
         allocate(DET_F(Nxy*Nz*N_Group,Nxy*Nz*N_Group))
         allocate(DET_R(Nxy*Nz*N_Group,Nxy*Nz*N_Group))
         allocate(DET_RINV(Nxy*Nz*N_Group,Nxy*Nz*N_Group))
         allocate(DET_D(DET_NPOW,Nxy*Nz*N_Group))
         allocate(DET_C(Nxy,Nz,N_Group))
         DET_A = 0d0
         DET_F = 0d0
         DET_R = 0d0
         DET_D = 0d0
         DET_C = 0d0
         DET_RINV=0d0
      endif

      Allocate(hist_6factor_eps  (1: buff_int))
      Allocate(hist_6factor_p    (1: buff_int))
      Allocate(hist_6factor_eta  (1: buff_int))
      Allocate(hist_6factor_f    (1: buff_int))
      Allocate(hist_6factor_fnl  (1: buff_int))
      Allocate(hist_6factor_tnl  (1: buff_int))
      Allocate(hist_6factor_kinf (1: buff_int))
      Allocate(hist_core_avg_nu  (1: buff_int))
      Allocate(hist_mfp          (1: buff_int))
      Allocate(hist_mfp_f        (1: buff_int))
      Allocate(hist_mfp_t        (1: buff_int))
      Allocate(hist_n_density    (1: buff_int))
      Allocate(hist_n_density_f  (1: buff_int))
      Allocate(hist_n_density_t  (1: buff_int))
      Allocate(hist_n_flux       (1: buff_int))
      Allocate(hist_n_flux_f     (1: buff_int))
      Allocate(hist_n_flux_t     (1: buff_int))
      Allocate(hist_n_speed      (1: buff_int))
      Allocate(hist_n_speed_f    (1: buff_int))
      Allocate(hist_n_speed_t    (1: buff_int))
      Allocate(hist_spectral_idx (1: buff_int))
      Allocate(hist_dif_length   (1: buff_int))
      Allocate(hist_dif_length_f (1: buff_int))
      Allocate(hist_dif_length_t (1: buff_int))
      Allocate(hist_mig_length   (1: buff_int))
      Allocate(hist_n_lifetime   (1: buff_int))
      Allocate(hist_n_lifetime_f (1: buff_int))
      Allocate(hist_n_lifetime_t (1: buff_int))
      Allocate(hist_n_gentime    (1: buff_int))
      Allocate(hist_n_gentime_f  (1: buff_int))
      Allocate(hist_n_gentime_t  (1: buff_int))
      Allocate(hist_n_energy     (1: buff_int))
      Allocate(hist_n_energy_f   (1: buff_int))
      Allocate(hist_n_energy_t   (1: buff_int))

      Allocate( v_Inv      (1: Nxy , 1:Nz , 1:N_Group   ))
      Allocate( beta_d     (1: Nxy , 1:Nz , 1:N_Group_d ))
      Allocate( beta_d_eff (1: Nxy , 1:Nz , 1:N_Group_d ))
      Allocate( beta_d_Tot (1: Nxy , 1:Nz               ))
      Allocate( lambda_d   (1: Nxy , 1:Nz , 1:N_Group_d ))

      nkincomp = 1

      IF ( .NOT. ALLOCATED(Flag_BP) ) THEN
         Allocate( Flag_BP (1: N_XS_Table ))
      ENDIF

      Allocate(Avg_Flux(1:N_Group))
      Avg_Flux = D0

      !ALLOCATE(D(N_Group))
      !ALLOCATE(maXS_tr(N_Group))
      !ALLOCATE(maXS_a(N_Group))
      !ALLOCATE(maXS_f(N_Group))
      !ALLOCATE(nu_maXS_f(N_Group))
      !ALLOCATE(kap_maXS_f(N_Group))
      !ALLOCATE(maXS_s(N_Group))
      !ALLOCATE(maXS_r(N_Group))
      !D          = D0
      !maXS_tr    = D0
      !maXS_a     = D0
      !maXS_f     = D0
      !nu_maXS_f  = D0
      !kap_maXS_f = D0
      !maXS_s     = D0
      !maXS_r     = D0

      Allocate(maXS_chi_3D(1: Nxy,1: Nz,1: N_Group))
      Allocate(maXS_chid_3D(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_scat_3D(1: N_Group,1:N_Group,1:Nxy,1:Nz))
      maXS_chi_3D(:,:,1)   = D1
      maXS_chi_3D(:,:,2)   = D0
      maXS_chid_3D(:,:,1)  = D1
      maXS_chid_3D(:,:,2)  = D0
      maXS_scat_3D          = D0

       Allocate(D_3D( -Nx-Ny+1:Nxy+Nx+Ny, 0:Nz+1, 1:N_Group))
      Allocate(maXS_tr_3D(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_a_3D(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_f_3D(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(nu_maXS_f_3D(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(kap_maXS_f_3D(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_s_3D(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_r_3D(1: Nxy, 1:Nz, 1:N_Group))
#ifdef tuan_fr
      Allocate(D_3D_MG( -Nx-Ny+1:Nxy+Nx+Ny, 0:Nz+1, 1:N_Group))
      Allocate(maXS_tr_3D_MG(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_a_3D_MG(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_f_3D_MG(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(nu_maXS_f_3D_MG(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(kap_maXS_f_3D_MG(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_s_3D_MG(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_r_3D_MG(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_chi_3D_MG(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_chid_3D_MG(1: Nxy, 1:Nz, 1:N_Group))
      Allocate(maXS_scat_3D_MG(1: N_Group,1:N_Group,1:Nxy,1:Nz))
#endif

#ifdef siarhei_tr_hex
      call alloc0(D_3D_CR, 1,N_CR, 0,Nz+1, 1,N_Group)
      call alloc(maXS_tr_3D_CR, N_CR, Nz, N_Group)
      call alloc(maXS_a_3D_CR, N_CR, Nz, N_Group)
      call alloc(maXS_f_3D_CR, N_CR, Nz, N_Group)
      call alloc(nu_maXS_f_3D_CR, N_CR, Nz, N_Group)
      call alloc(kap_maXS_f_3D_CR, N_CR, Nz, N_Group)
      call alloc(maXS_s_3D_CR, N_CR, Nz, N_Group)
      call alloc(maXS_r_3D_CR, N_CR, Nz, N_Group)
      call alloc(maXS_chi_3D_CR, N_CR, Nz, N_Group)
      call alloc(maXS_chid_3D_CR, N_CR, Nz, N_Group)
      call alloc(maXS_scat_3D_CR, N_Group,N_Group,N_CR,Nz)
      if (.not.allocated(XSset_Hex_CR ))  &
         CALL Alloc(XSset_Hex_CR, N_CR, Nz, N_Group + 6, N_Group ) ! for recording all rodded cases
      if (crit_search_CR) then
         call alloc(history_of_CR_search,Iout_Max,1 + N_CR + 1) ! number of iterations (max), iter - all CR positions - keff
         history_of_CR_search = -444.4_8
      endif


      D_3D_CR          = 0d0
      maXS_tr_3D_CR    = 0d0
      maXS_a_3D_CR     = 0d0
      maXS_f_3D_CR     = 0d0
      nu_maXS_f_3D_CR  = 0d0
      kap_maXS_f_3D_CR = 0d0
      maXS_s_3D_CR     = 0d0
      maXS_r_3D_CR     = 0d0
      maXS_chi_3D_CR   = 0d0
      maXS_chid_3D_CR  = 0d0
      maXS_scat_3D_CR  = 0d0

      XSset_Hex_CR = 0d0
#endif


#ifdef jr_vver
      if(opt_nodal==4) then
         If (.NOT. ALLOCATED(maXS_r_3D)) Allocate(maXS_r_3D(1: Nxy, 1:Nz, 1:N_Group))
         maXS_r_3D = 0d0
      endif
#endif
      D_3D          = D0
      maXS_tr_3D    = D0
      maXS_a_3D     = D0
      maXS_f_3D     = D0
      nu_maXS_f_3D  = D0
      kap_maXS_f_3D = D0
      maXS_s_3D     = D0
      maXS_r_3D     = D0

      Allocate(maXS_f_3D_Old(1: Nxy, 1:Nz, 1:N_Group))
      maXS_f_3D_Old = D0

      if (flag_leakage) then
         Allocate(D_3D_lc( -Nx-Ny+1:Nxy+Nx+Ny, 0:Nz+1, 1:N_Group))

         Allocate(maXS_a_3D_lc(1: Nxy, 1:Nz, 1:N_Group))
         Allocate(maXS_s_3D_lc(1: Nxy, 1:Nz, 1:N_Group))
         Allocate(buck2_lc(1: Nxy, 1:Nz, 1:N_Group))
         Allocate(dbuck2_lc(1: Nxy, 1:Nz, 1:N_Group))
         Allocate(leak_ratio(1: Nxy, 1:Nz, 1:4))
         D_3D_lc = D0
         maXS_a_3D_lc = D0
         maXS_s_3D_lc = D0
         buck2_lc = D0
         dbuck2_lc = D0
         leak_ratio = D0
      endif

      if (flag_det_sig) then
         Allocate(DET_FRSET(1:N_Group))
         Allocate(DET_FR(1:Nxy,1:Nz,1:N_Group))
         Allocate(DET_Flux_XYZ_G(Nxy_1N,1:Nz,1:N_Group))
         Allocate(DET_Flux_XYZ_G_Old(1:Nxy_1N,1:Nz,1:N_Group))
         Allocate(DET_I_FA(1:Nxy_1N))
         DET_FRSET = D0
         DET_FR = D0
         DET_Flux_XYZ_G = D0
         DET_Flux_XYZ_G_Old = D0
         DET_I_FA = 0
      endif

      if (flag_midepl.or.flag_XSFB) then
         Allocate( miXS_tr_U34  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_U35  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_U36  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_U37  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_U38  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Np37 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Np38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Np39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Pu38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Pu39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Pu40 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Pu41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Pu42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Pu43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Am41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_As42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Am42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Am43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Am44 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Cm42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Cm43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Cm44 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_a_U34  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_U35  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_U36  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_U37  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_U38  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Np37 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Np38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Np39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu40 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Am41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_As42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Am42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Am43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Am44 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Cm42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Cm43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Cm44 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_a_U34_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_U35_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_U36_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_U37_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_U38_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Np37_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Np38_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Np39_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu38_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu39_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu40_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu41_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu42_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pu43_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Am41_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_As42_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Am42_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Am43_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Am44_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Cm42_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Cm43_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Cm44_Old (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_f_U34  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_U35  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_U36  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_U37  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_U38  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Np37 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Np38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Np39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu40 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Am41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_As42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Am42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Am43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Am44 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Cm42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Cm43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Cm44 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_f_U34_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_U35_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_U36_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_U37_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_U38_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Np37_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Np38_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Np39_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu38_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu39_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu40_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu41_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu42_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Pu43_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Am41_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_As42_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Am42_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Am43_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Am44_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Cm42_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Cm43_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_f_Cm44_Old (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( nu_miXS_f_U34  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_U35  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_U36  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_U37  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_U38  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Np37 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Np38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Np39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Pu38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Pu39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Pu40 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Pu41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Pu42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Pu43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Am41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_As42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Am42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Am43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Am44 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Cm42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Cm43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( nu_miXS_f_Cm44 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( kap_miXS_f_U34  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_U35  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_U36  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_U37  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_U38  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Np37 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Np38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Np39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Pu38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Pu39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Pu40 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Pu41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Pu42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Pu43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Am41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_As42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Am42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Am43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Am44 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Cm42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Cm43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( kap_miXS_f_Cm44 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_s_U34  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_U35  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_U36  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_U37  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_U38  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Np37 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Np38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Np39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Pu38 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Pu39 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Pu40 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Pu41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Pu42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Pu43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Am41 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_As42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Am42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Am43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Am44 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Cm42 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Cm43 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Cm44 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_tr_I35  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Xe35 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Nd47 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Nd48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Nd49 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Pm47 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Ps48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Pm48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Pm49 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Sm47 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Sm48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Sm49 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_a_I35  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Xe35 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Nd47 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Nd48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Nd49 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pm47 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Ps48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pm48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pm49 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Sm47 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Sm48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Sm49 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_a_I35_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Xe35_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Nd47_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Nd48_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Nd49_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pm47_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Ps48_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pm48_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Pm49_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Sm47_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Sm48_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Sm49_Old (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_s_I35  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Xe35 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Nd47 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Nd48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Nd49 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Pm47 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Ps48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Pm48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Pm49 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Sm47 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Sm48 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Sm49 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_tr_Gd52 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Gd54 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Gd55 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Gd56 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Gd57 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Gd58 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_Gd60 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_a_Gd52 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd54 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd55 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd56 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd57 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd58 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd60 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_a_Gd52_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd54_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd55_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd56_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd57_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd58_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Gd60_Old (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_s_Gd52 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Gd54 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Gd55 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Gd56 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Gd57 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Gd58 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_Gd60 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_tr_H2O (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_tr_B0  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_H2O (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_B0  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_H2O_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_B0_Old  (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_a_B10_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_B10 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_s_H2O (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_s_B0  (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_n2n_U35  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_n2n_U38  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_n2n_Np37 (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_n2n_Pu39 (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_n2n_U35_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_n2n_U38_Old  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_n2n_Np37_Old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_n2n_Pu39_Old (1: Nxy ,1: Nz ,1: N_Group ))
      endif
      !js+if (flag_midepl) then ! need to check later... for transient case
         Allocate( Y_I35      (1: Nxy ,1: Nz ))
         Allocate( Y_Xe35     (1: Nxy ,1: Nz ))
         Allocate( Y_Nd47     (1: Nxy ,1: Nz ))
         Allocate( Y_Nd48     (1: Nxy ,1: Nz ))
         Allocate( Y_Nd49     (1: Nxy ,1: Nz ))
         Allocate( Y_Pm47     (1: Nxy ,1: Nz ))
         Allocate( Y_Ps48     (1: Nxy ,1: Nz ))
         Allocate( Y_Pm48     (1: Nxy ,1: Nz ))
         Allocate( Y_Pm49     (1: Nxy ,1: Nz ))
         Allocate( Y_Sm49     (1: Nxy ,1: Nz ))
         Allocate( Y_Xe35_eff (1: Nxy ,1: Nz ))
         Allocate( Y_Pm49_eff (1: Nxy ,1: Nz ))

         Allocate( Y_I35_Old      (1: Nxy ,1: Nz ))
         Allocate( Y_Xe35_Old     (1: Nxy ,1: Nz ))
         Allocate( Y_Nd47_Old     (1: Nxy ,1: Nz ))
         Allocate( Y_Nd48_Old     (1: Nxy ,1: Nz ))
         Allocate( Y_Nd49_Old     (1: Nxy ,1: Nz ))
         Allocate( Y_Pm47_Old     (1: Nxy ,1: Nz ))
         Allocate( Y_Ps48_Old     (1: Nxy ,1: Nz ))
         Allocate( Y_Pm48_Old     (1: Nxy ,1: Nz ))
         Allocate( Y_Pm49_Old     (1: Nxy ,1: Nz ))
         Allocate( Y_Sm49_Old     (1: Nxy ,1: Nz ))
         Allocate( Y_Xe35_eff_Old (1: Nxy ,1: Nz ))
         Allocate( Y_Pm49_eff_Old (1: Nxy ,1: Nz ))
      !js+endif

      IF ( .not. Flag_RI ) THEN
         Allocate( N_U36  (1: Nxy ,1: Nz ))
         Allocate( N_U37  (1: Nxy ,1: Nz ))
         Allocate( N_Np37 (1: Nxy ,1: Nz ))
         Allocate( N_Np38 (1: Nxy ,1: Nz ))
         Allocate( N_Np39 (1: Nxy ,1: Nz ))
         Allocate( N_Pu38 (1: Nxy ,1: Nz ))
         Allocate( N_Pu39 (1: Nxy ,1: Nz ))
         Allocate( N_Pu40 (1: Nxy ,1: Nz ))
         Allocate( N_Pu41 (1: Nxy ,1: Nz ))
         Allocate( N_Pu42 (1: Nxy ,1: Nz ))
         Allocate( N_Pu43 (1: Nxy ,1: Nz ))
         Allocate( N_Am41 (1: Nxy ,1: Nz ))
         Allocate( N_As42 (1: Nxy ,1: Nz ))
         Allocate( N_Am42 (1: Nxy ,1: Nz ))
         Allocate( N_Am43 (1: Nxy ,1: Nz ))
         Allocate( N_Am44 (1: Nxy ,1: Nz ))
         Allocate( N_Cm42 (1: Nxy ,1: Nz ))
         Allocate( N_Cm43 (1: Nxy ,1: Nz ))
         Allocate( N_Cm44 (1: Nxy ,1: Nz ))

         Allocate( N_U36_Old  (1: Nxy ,1: Nz ))
         Allocate( N_U37_Old  (1: Nxy ,1: Nz ))
         Allocate( N_Np37_Old (1: Nxy ,1: Nz ))
         Allocate( N_Np38_Old (1: Nxy ,1: Nz ))
         Allocate( N_Np39_Old (1: Nxy ,1: Nz ))
         Allocate( N_Pu38_Old (1: Nxy ,1: Nz ))
         Allocate( N_Pu39_Old (1: Nxy ,1: Nz ))
         Allocate( N_Pu40_Old (1: Nxy ,1: Nz ))
         Allocate( N_Pu41_Old (1: Nxy ,1: Nz ))
         Allocate( N_Pu42_Old (1: Nxy ,1: Nz ))
         Allocate( N_Pu43_Old (1: Nxy ,1: Nz ))
         Allocate( N_Am41_Old (1: Nxy ,1: Nz ))
         Allocate( N_As42_Old (1: Nxy ,1: Nz ))
         Allocate( N_Am42_Old (1: Nxy ,1: Nz ))
         Allocate( N_Am43_Old (1: Nxy ,1: Nz ))
         Allocate( N_Am44_Old (1: Nxy ,1: Nz ))
         Allocate( N_Cm42_Old (1: Nxy ,1: Nz ))
         Allocate( N_Cm43_Old (1: Nxy ,1: Nz ))
         Allocate( N_Cm44_Old (1: Nxy ,1: Nz ))

         Allocate( N_I35  (1: Nxy ,1: Nz ))
         Allocate( N_Xe35 (1: Nxy ,1: Nz ))
         Allocate( N_Nd47 (1: Nxy ,1: Nz ))
         Allocate( N_Nd48 (1: Nxy ,1: Nz ))
         Allocate( N_Nd49 (1: Nxy ,1: Nz ))
         Allocate( N_Pm47 (1: Nxy ,1: Nz ))
         Allocate( N_Ps48 (1: Nxy ,1: Nz ))
         Allocate( N_Pm48 (1: Nxy ,1: Nz ))
         Allocate( N_Pm49 (1: Nxy ,1: Nz ))
         Allocate( N_Sm47 (1: Nxy ,1: Nz ))
         Allocate( N_Sm48 (1: Nxy ,1: Nz ))
         Allocate( N_Sm49 (1: Nxy ,1: Nz ))

         Allocate( N_I35_Old  (1: Nxy ,1: Nz ))
         Allocate( N_Xe35_Old (1: Nxy ,1: Nz ))
         Allocate( N_Nd47_Old (1: Nxy ,1: Nz ))
         Allocate( N_Nd48_Old (1: Nxy ,1: Nz ))
         Allocate( N_Nd49_Old (1: Nxy ,1: Nz ))
         Allocate( N_Pm47_Old (1: Nxy ,1: Nz ))
         Allocate( N_Ps48_Old (1: Nxy ,1: Nz ))
         Allocate( N_Pm48_Old (1: Nxy ,1: Nz ))
         Allocate( N_Pm49_Old (1: Nxy ,1: Nz ))
         Allocate( N_Sm47_Old (1: Nxy ,1: Nz ))
         Allocate( N_Sm48_Old (1: Nxy ,1: Nz ))
         Allocate( N_Sm49_Old (1: Nxy ,1: Nz ))

         Allocate( N_H2O (1: Nxy ,1: Nz ))
         Allocate( N_B0  (1: Nxy ,1: Nz ))

      END IF


      IF ( Flag_Card_maXS ) THEN
         Allocate( N_U34  (1: Nxy ,1: Nz ))
         Allocate( N_U35  (1: Nxy ,1: Nz ))
         Allocate( N_U38  (1: Nxy ,1: Nz ))

         Allocate( N_U34_Old  (1: Nxy ,1: Nz ))
         Allocate( N_U35_Old  (1: Nxy ,1: Nz ))
         Allocate( N_U38_Old  (1: Nxy ,1: Nz ))

         Allocate( N_Gd52 (1: Nxy ,1: Nz ))
         Allocate( N_Gd54 (1: Nxy ,1: Nz ))
         Allocate( N_Gd55 (1: Nxy ,1: Nz ))
         Allocate( N_Gd56 (1: Nxy ,1: Nz ))
         Allocate( N_Gd57 (1: Nxy ,1: Nz ))
         Allocate( N_Gd58 (1: Nxy ,1: Nz ))
         Allocate( N_Gd60 (1: Nxy ,1: Nz ))

         Allocate( N_Gd52_Old (1: Nxy ,1: Nz ))
         Allocate( N_Gd54_Old (1: Nxy ,1: Nz ))
         Allocate( N_Gd55_Old (1: Nxy ,1: Nz ))
         Allocate( N_Gd56_Old (1: Nxy ,1: Nz ))
         Allocate( N_Gd57_Old (1: Nxy ,1: Nz ))
         Allocate( N_Gd58_Old (1: Nxy ,1: Nz ))
         Allocate( N_Gd60_Old (1: Nxy ,1: Nz ))
      END IF

      if (flag_sim) then
         if (allocated(macxs_b0)) deallocate(macxs_b0)
         Allocate(macxs_b0(1:2,1:nxy,1:nz))
         macxs_b0(:,:,:)=0d0
      endif

      if (flag_det_sig) then
         Allocate( DET_miXS_a_A1  (1: DET_N_XY ,1: DET_N_Z , 1:N_Group ))
         Allocate( DET_miXS_a_A2  (1: DET_N_XY ,1: DET_N_Z , 1:N_Group ))
         Allocate( DET_miXS_a_A2m (1: DET_N_XY ,1: DET_N_Z , 1:N_Group ))
         Allocate( DET_miXS_a_A3  (1: DET_N_XY ,1: DET_N_Z , 1:N_Group ))
         Allocate( DET_miXS_a_A4  (1: DET_N_XY ,1: DET_N_Z , 1:N_Group ))
         Allocate( DET_miXS_a_B2  (1: DET_N_XY ,1: DET_N_Z , 1:N_Group ))
         Allocate( DET_miXS_a_B3  (1: DET_N_XY ,1: DET_N_Z , 1:N_Group ))

         Allocate( DET_miXS_a_V1  (1: DET_N_XY ,1: DET_N_Z , 1:N_Group ))
         Allocate( DET_miXS_a_V2  (1: DET_N_XY ,1: DET_N_Z , 1:N_Group ))
         Allocate( DET_miXS_a_Cr  (1: DET_N_XY ,1: DET_N_Z , 1:N_Group ))

         Allocate( DET_N_A1  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_A2  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_A2m (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_A3  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_A4  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_B2  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_B3  (1: DET_N_XY , 1:DET_N_Z ))

         Allocate( DET_N_A1_Old  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_A2_Old  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_A2m_Old (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_A3_Old  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_A4_Old  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_B2_Old  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_B3_Old  (1: DET_N_XY , 1:DET_N_Z ))

         Allocate( DET_N_V1_Old   (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_V2_Old   (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_Cr52_Old (1: DET_N_XY , 1:DET_N_Z ))

         Allocate( DET_N_V1   (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_V2   (1: DET_N_XY , 1:DET_N_Z ))
         Allocate( DET_N_Cr52 (1: DET_N_XY , 1:DET_N_Z ))

         Allocate(DET_I_beta_A2    (1: DET_N_XY , 1:DET_N_Z ))
         Allocate(DET_I_beta_A2m   (1: DET_N_XY , 1:DET_N_Z ))
         Allocate(DET_I_beta_A3    (1: DET_N_XY , 1:DET_N_Z ))
         Allocate(DET_I_beta_B2    (1: DET_N_XY , 1:DET_N_Z ))
         Allocate(DET_I_beta_B3    (1: DET_N_XY , 1:DET_N_Z ))
         Allocate(DET_I_beta_rh    (1: DET_N_XY , 1:DET_N_Z ))
         Allocate(DET_I_beta_rh_p  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate(DET_I_beta_rh_d  (1: DET_N_XY , 1:DET_N_Z ))
         Allocate(DET_I_beta_V2    (1: DET_N_XY , 1:DET_N_Z ))
         Allocate(DET_I_beta_si    (1: DET_N_XY , 1:DET_N_Z ))
      endif

      if (Flag_InDET) then
         Allocate( miXS_a_A1  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_A2  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_A2m (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_A3  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_A4  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_B2  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_B3  (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( miXS_a_V1  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_V2  (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( miXS_a_Cr  (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate( N_A1 (1: Nxy ,1: Nz ))
         Allocate( N_A2 (1: Nxy ,1: Nz ))
         Allocate( N_A2m(1: Nxy ,1: Nz ))
         Allocate( N_A3 (1: Nxy ,1: Nz ))
         Allocate( N_A4 (1: Nxy ,1: Nz ))
         Allocate( N_B2 (1: Nxy ,1: Nz ))
         Allocate( N_B3 (1: Nxy ,1: Nz ))

         Allocate( N_A1_Old (1: Nxy ,1: Nz ))
         Allocate( N_A2_Old (1: Nxy ,1: Nz ))
         Allocate( N_A2m_Old(1: Nxy ,1: Nz ))
         Allocate( N_A3_Old (1: Nxy ,1: Nz ))
         Allocate( N_A4_Old (1: Nxy ,1: Nz ))
         Allocate( N_B2_Old (1: Nxy ,1: Nz ))
         Allocate( N_B3_Old (1: Nxy ,1: Nz ))

         Allocate( Flux_det_old (1: Nxy ,1: Nz ,1: N_Group ))
         Allocate( Flux_det     (1: Nxy ,1: Nz ,1: N_Group ))

         Allocate(I_beta_A2  ( 1:Nxy , 1:Nz ))
         Allocate(I_beta_A2m ( 1:Nxy , 1:Nz ))
         Allocate(I_beta_A3  ( 1:Nxy , 1:Nz ))
         Allocate(I_beta_B2  ( 1:Nxy , 1:Nz ))
         Allocate(I_beta_B3  ( 1:Nxy , 1:Nz ))
         Allocate(I_beta_rh  ( 1:Nxy , 1:Nz ))
         Allocate(I_beta_rh_p( 1:Nxy , 1:Nz ))
         Allocate(I_beta_rh_d( 1:Nxy , 1:Nz ))

         Allocate(  N_V1_Old   (1: Nxy ,1: Nz ))
         Allocate(  N_V2_Old   (1: Nxy ,1: Nz ))
         Allocate(  N_Cr52_Old (1: Nxy ,1: Nz ))

         Allocate(  N_V1   (1: Nxy ,1: Nz ))
         Allocate(  N_V2   (1: Nxy ,1: Nz ))
         Allocate(  N_Cr52 (1: Nxy ,1: Nz ))

         Allocate(I_beta_V2(1: Nxy ,1: Nz ))
         Allocate(I_beta_si(1: Nxy ,1: Nz ))
      endif

      Allocate( T_Mod_Out (1: Nxy) )
      Allocate( LinPOW    ( 1:Nxy , 1:Nz ) )

      Allocate( L_omega_Surf   (1: 6          ))
      Allocate( Q_Surf         (1: Nxy, 1:Nz, 1:6 ))
      Allocate( Q_Surf_Old     (1: Nxy, 1:Nz, 1:6 ))
      Allocate( C_d_I          (1: Nxy, 1:Nz, 1:N_Group_d))
      Allocate( C_d_I_Old      (1: Nxy, 1:Nz, 1:N_Group_d))
      Allocate( C_d_I_Surf     (1: Nxy, 1:Nz, 1:N_Group_d, 1:6 ))
      Allocate( C_d_I_Surf_Old (1: Nxy, 1:Nz, 1:N_Group_d, 1:6 ))
      L_omega_Surf = D0
      Q_Surf     = D0
      Q_Surf_Old = D0
      C_d_I     = D0
      C_d_I_Old = D0
      C_d_I_Surf     = D0
      C_d_I_Surf_Old = D0

      Allocate( S_0  (1: N_Group ))
      Allocate( S_1x (1: N_Group ))
      Allocate( S_1y (1: N_Group ))
      Allocate( S_1z (1: N_Group ))
      Allocate( S_2x (1: N_Group ))
      Allocate( S_2y (1: N_Group ))
      Allocate( S_2z (1: N_Group ))

      Allocate( S_Surf        (1: N_Group, 1:6 ))
      Allocate( Flux_Surf     (1: Nxy, 1:Nz, 1:N_Group, 1:6 ))
      Allocate( Flux_Surf_Old (1: Nxy, 1:Nz, 1:N_Group, 1:6 ))
      Allocate( chi_p         (1: N_Group))
      Allocate( chi_d         (1: N_Group))
      S_Surf = D0
      Flux_Surf     = D0
      Flux_Surf_Old = D0
      chi_p(1) = D1
      chi_p(2) = D0
      chi_d(1) = D1
      chi_d(2) = D0

      Allocate( Flag_Reg ( 1:N_Region ))

      Allocate( Normal_Power_Z   (1: Nz                 ))
      Allocate( Normal_Power_XY  (1: Nx_1N , 1:Ny_1N      ))
      Allocate( Normal_Power_XYZ (1: Nx_1N , 1:Ny_1N , 1:Nz ))
      Allocate( BU_Z             (1: Nz                 ))
      Allocate( BU_XY            (1: Nx_1N , 1:Ny_1N      ))
      Allocate( BU_XYZ           (1: Nx_1N , 1:Ny_1N , 1:Nz ))
      Allocate( T_Fuel_Z         (1: Nz                 ))
      Allocate( T_Fuel_XY        (1: Nx_1N , 1:Ny_1N      ))
      Allocate( T_Fuel_XYZ       (1: Nx_1N , 1:Ny_1N , 1:Nz ))
      Allocate( T_Mod_Z          (1: Nz                 ))
      Allocate( T_Mod_XY         (1: Nx_1N , 1:Ny_1N      ))
      Allocate( T_Mod_XYZ        (1: Nx_1N , 1:Ny_1N , 1:Nz ))
      Allocate( D_Mod_Z          (1: Nz                 ))
      Allocate( D_Mod_XY         (1: Nx_1N , 1:Ny_1N      ))
      Allocate( D_Mod_XYZ        (1: Nx_1N , 1:Ny_1N , 1:Nz ))

      Allocate( Fxy_val_XY       (1: Nx_1N , 1:Ny_1N      ))
      Allocate( Hist_BU_XYZ_4N_ANC(1: Nx,    1:Ny, 1:Nz ))

      if (flag_ocean) then
         Allocate( N_Xe35_XYZ      (1: Nx_1N ,1: Ny_1N ,1: Nz ))
         Allocate( N_Sm49_XYZ      (1: Nx_1N ,1: Ny_1N ,1: Nz ))
         Allocate( FFlux_XYZ       (1: Nx_1N ,1: Ny_1N ,1: Nz ))
         Allocate( TFlux_XYZ       (1: Nx_1N ,1: Ny_1N ,1: Nz ))
      endif

      Allocate( fFlux_adj_XYZ   ( 1:nx_1N, 1:ny_1N, 1:nz ))
      Allocate( tFlux_adj_XYZ   ( 1:nx_1N, 1:ny_1N, 1:nz ))
      Allocate( fFlux_adj_XY    ( 1:nx_1N, 1:ny_1N ))
      Allocate( tFlux_adj_XY    ( 1:nx_1N, 1:ny_1N ))

     Allocate( Normal_Power_XY_4N  (1: Nx , 1:Ny      ))
     Allocate( BU_XY_4N  ( 1:Nx , 1:Ny      ))

      If (.NOT. ALLOCATED(FA_ASI)) Allocate( FA_ASI(1: Nx_1N , 1:Ny_1N ))

      RETURN
      END SUBROUTINE Set_Memory
!!!#endif 


      SUBROUTINE Check_File(Chr)
      use mod_charedit, only: print_msg
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: Chr
      LOGICAL(1) :: Exist_File


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Check_File] in Mod_Manage'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      INQUIRE(FILE = Chr, EXIST = Exist_File)

      IF (Exist_File .EQV. .TRUE.) THEN
         RETURN
      ELSE
         call print_msg(3,chr)
      END IF

      RETURN
      END SUBROUTINE Check_File


!!!#ifdef siarhei_delete 
      SUBROUTINE ErrMSG(ErrCase, ErrClue_Real, ErrClue_Char)

      IMPLICIT NONE

      REAL(8), INTENT(IN) :: ErrCase, ErrClue_Real
      CHARACTER(*), INTENT(IN) :: ErrClue_Char

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [ErrMSG] in Mod_Manage'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF (ErrCase == 1.01D0) THEN
         WRITE(*, "(A, F4.2, 2A)") "*** [ ERROR_", ErrCase, " ] ",  &
            "Check File Existence or Name"
         WRITE(*, "(2A)") "*** [  LOCATION  ] ",  &
            "Main/Manage/Check_File"
         WRITE(*, "(2A)") "*** [    CLUE    ] ", ErrClue_Char

      ELSE IF (ErrCase == 1.02D0) THEN

      ELSE IF (ErrCase == 2.00D0) THEN
         WRITE(*, "(A, F4.2, 2A)") "*** [ ERROR_", ErrCase, " ] ",  &
            "Check Card Name'"
         WRITE(*, "(2A)") "*** [  LOCATION  ] ",  &
            "RW_File/Read_I/Read_Inp"
         WRITE(*, "(3A)") "*** [    CLUE    ] ", "Current Card Name is ", ErrClue_Char

      ELSE IF (ErrCase == 2.10D0) THEN
         WRITE(*, "(A, F4.2, 2A)") "*** [ ERROR_", ErrCase, " ] ",  &
            "Check SubCard Name"
         WRITE(*, "(2A)") "*** [  LOCATION  ] ",  &
            "RW_File/Read_I/Read_Inp/Card_maXS"
         WRITE(*, "(3A)") "*** [    CLUE    ] ", "Current SubCard Name is ", ErrClue_Char

      ELSE IF (ErrCase == 2.11D0) THEN
         WRITE(*, "(A, F4.2, 2A)") "*** [ ERROR_", ErrCase, " ] ",  &
            "Check SubCard Name"
         WRITE(*, "(2A)") "*** [  LOCATION  ] ",  &
            "RW_File/Read_I/Read_Inp/Card_Transient"
         WRITE(*, "(3A)") "*** [    CLUE    ] ", "Current SubCard Name is ", ErrClue_Char

      ELSE IF (ErrCase == 2.12D0) THEN
         WRITE(*, "(A, F4.2, 2A)") "*** [ ERROR_", ErrCase, " ] ",  &
            "Check SubCard Name"
         WRITE(*, "(2A)") "*** [  LOCATION  ] ",  &
            "RW_File/Read_I/Read_Inp/Card_Branch"
         WRITE(*, "(3A)") "*** [    CLUE    ] ", "Current SubCard Name is ", ErrClue_Char

      ELSE IF (ErrCase == 2.13D0) THEN
         WRITE(*, "(A, F4.2, 2A)") "*** [ ERROR_", ErrCase, " ] ",  &
            "Check Type_Branch Name. It must be 'FTC' or 'MTC' or 'CRW'"
         WRITE(*, "(2A)") "*** [  LOCATION  ] ",  &
            "RW_File/Read_I/Read_Inp/Card_Branch"
         WRITE(*, "(3A)") "*** [    CLUE    ] ", "Current Type_Branch Name is ", ErrClue_Char

      ELSE IF (ErrCase == 3.01D0) THEN
         WRITE(*, "(A, F4.2, 2A)") "*** [ ERROR_", ErrCase, " ] ",  &
            "Check Flag in 'CALL Get_Avg'. It must be '0 ~ 2'"
         WRITE(*, "(2A)") "*** [  LOCATION  ] ",  &
            "Function/Get_Something/Get_Avg"
         WRITE(*, "(2A, F4.2)") "*** [    CLUE    ] ", "Current Flag is ", ErrClue_Real

      END IF

      STOP

      RETURN
      END SUBROUTINE ErrMSG
!!!#endif 


      subroutine alloc_lscoef
      use inc_lscoef
      use inc_option, only: n_group
      use mod_alloc

#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: iproc2nxy
      use inc_parallel, only: nnodelmpi
#endif
      implicit none
      integer :: nxy2nx, nzp1
      integer :: nxy_bk, nxy2nx_bk


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [alloc_lscoef] in Mod_Manage'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      nxy2nx=nxy+nx+ny
      nzp1=nz+1

      nxy_bk=nxy
      nxy2nx_bk=nxy2nx
#ifdef js_mpi
      if (comm%usempi) then
         nxy_bk=iproc2nxy(iproc)
         nxy2nx_bk=nnodelmpi
      endif
#endif

      ! alloc_lscoef
      call alloc (ccw  , n_group , nxy_bk    , nz   )
      call alloc (cce  , n_group , nxy_bk    , nz   )
      call alloc (ccn  , n_group , nxy_bk    , nz   )
      call alloc (ccs  , n_group , nxy_bk    , nz   )
      call alloc (ccb  , n_group , nxy_bk    , nz   )
      call alloc (cct  , n_group , nxy_bk    , nz   )

      call alloc (dfw  , n_group , nxy2nx_bk , nz   )
      call alloc (dfn  , n_group , nxy2nx_bk , nz   )
      call alloc (dfb  , n_group , nxy_bk    , nzp1 )


      call alloc (dnw  , n_group , nxy2nx_bk , nz   )
      call alloc (dnn  , n_group , nxy2nx_bk , nz   )
      call alloc (dnb  , n_group , nxy_bk    , nzp1 )

      call alloc (am    , n_group ,           nxy_bk , nz )
      call alloc (amcc  , n_group ,           nxy_bk , nz )
      call alloc (af    , n_group ,           nxy_bk , nz )
      call alloc (af2   ,                     nxy_bk , nz )
      call alloc (scat  ,                     nxy_bk , nz )
      call alloc (scatv , n_group , n_group , nxy_bk , nz )

      return
      end subroutine alloc_lscoef


      subroutine alloc_lucoef
      use inc_lscoef
      use mod_alloc
      implicit none

      ! alloc_lu

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [alloc_lucoef] in Mod_Manage'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call alloc (del    , 4 , nx )
      call alloc (au     , 4 , nx )
      call alloc (ainvl  , 4 , nx )
      call alloc (ainvu  , 4 , nx )
      call alloc (ainvd  , 4 , nx )
      call alloc (delinv , 4 , nxy , nz )
      call alloc (al     , 4 , nxy , nz )
      call alloc (deliau , 4 , nxy , nz )

      return
      end subroutine alloc_lucoef


      subroutine alloc_flux
      use inc_fluxvar
      use inc_option, only: n_group

      use mod_alloc
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: iproc2nxy
      use inc_parallel, only: nnodelmpi
#endif
      implicit none
      integer :: nx2, nx2m1, nx4m1, nzp1, nxy2nx
      integer :: nxy_bk, nxy2nx_bk


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [alloc_flux] in Mod_Manage'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      nxy2nx=nxy+nx+ny
      nx2=nx+ny
      nx2m1=nx2-1
      nx4m1=nx2m1+nx2
      nzp1=nz+1

      nxy_bk=nxy
      nxy2nx_bk=nxy2nx
#ifdef js_mpi
      if (comm%usempi) then
         nxy_bk=iproc2nxy(iproc)
         nxy2nx_bk=nnodelmpi
         nx4m1=0
         nx2m1=0
      endif
#endif


      call alloc (flux_add , nxy , nz , n_group )
      call alloc (vr       , nxy , nz , n_group )
      call alloc (src   ,   n_group ,        nxy ,   nz )
      call alloc0(src_tr, 1,n_group , -nx2m1,nxy2nx_bk , 0,nz ) !TODO...mpi

      call alloc (j_net_x_3D , n_group , nxy2nx_bk , nz   )
      call alloc (j_net_y_3D , n_group , nxy2nx_bk , nz   )
      call alloc (j_net_z_3D , n_group , nxy_bk    , nzp1 )

      call alloc0(L_0x_3D , 1,n_group , -nx4m1,nxy2nx_bk , -1,nzp1 )
      call alloc0(L_1x_3D , 1,n_group , -nx2m1,nxy2nx_bk ,  1,nz   )
      call alloc0(L_2x_3D , 1,n_group , -nx2m1,nxy2nx_bk ,  1,nz   )
      call alloc0(L_0y_3D , 1,n_group , -nx4m1,nxy2nx_bk , -1,nzp1 )
      call alloc0(L_1y_3D , 1,n_group , -nx2m1,nxy2nx_bk ,  1,nz   )
      call alloc0(L_2y_3D , 1,n_group , -nx2m1,nxy2nx_bk ,  1,nz   )
      call alloc0(L_0z_3D , 1,n_group , -nx4m1,nxy2nx_bk ,  1,nz   ) !TODO...mpi
      call alloc0(L_1z_3D , 1,n_group ,      1,nxy_bk    ,  0,nz   )
      call alloc0(L_2z_3D , 1,n_group ,      1,nxy_bk    ,  0,nz   )

      return
      end subroutine alloc_flux


      subroutine alloc_bicgs
      use inc_option, only: n_group
      use inc_bicg
      use mod_alloc
      implicit none


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [alloc_bicgs] in Mod_Manage'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call alloc ( vr0,   nxy   ,   nz ,   n_group )
      call alloc ( vp ,   nxy   ,   nz ,   n_group )
      call alloc ( vv ,   nxy   ,   nz ,   n_group )
      call alloc ( vs ,   nxy   ,   nz ,   n_group )
      call alloc ( vt ,   nxy   ,   nz ,   n_group )
      call alloc0( vy , 0,nxy+1 , 1,nz , 1,n_group )
      call alloc0( vz , 0,nxy+1 , 1,nz , 1,n_group )

      call alloc0(s    , 1,n_group , -nx+1 , nxy+nx )
      call alloc (b0   ,   n_group ,         nxy    )
      call alloc (s1dl ,   n_group ,  nx            )
      call alloc (b01d ,   n_group ,  nx            )
      call alloc (y    ,   n_group ,  nx            )

      return
      end subroutine alloc_bicgs

      subroutine alloc_tran
      use inc_xs_file, only: n_xs_table
      use inc_option, only: n_group, opt_mode
      use inc_kinetics, only: n_group_d
      use inc_solver
      use mod_alloc
      use inc_option, only: flag_print_kp
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: iproc2nxy
#endif
#ifdef tuan_tr_test
      use Inc_FixSrc_Tr
#endif
      implicit none
      logical(1), save :: flag_first
      integer :: nxy_bk
      if (.not.flag_first) then

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [alloc_tran] in Mod_Manage'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         flag_first=.true.
      else
         return
      endif

      nxy_bk=nxy
#ifdef js_mpi
      if (comm%usempi) then
         nxy_bk=iproc2nxy(iproc)
      endif
#endif

#ifdef jr_test
write(*,*) 'hh 1'
#endif
      call alloc (rvdelt     , n_group , nxy_bk , nz )
#ifdef jr_test
write(*,*) 'hh 2'
#endif
      call alloc (rvdelt_old , n_group , nxy_bk , nz )
#ifdef jr_test
write(*,*) 'hh 3'
#endif
      call alloc (betap      ,           nxy_bk , nz )
#ifdef jr_test
write(*,*) 'hh 4'
#endif
      call alloc (betap_old  ,           nxy_bk , nz )
#ifdef jr_test
write(*,*) 'hh 5'
#endif

      if (flag_print_kp.or.opt_mode==3) then
#ifdef jr_test
write(*,*) 'hh in'
#endif


         call alloc (kincomp    , n_xs_table )
         call alloc (lambda_d_I , n_group_d , nkincomp )
         call alloc (tvelo      , n_group   , nkincomp )
         call alloc (tbeta      , n_group_d , nkincomp )

         call alloc (cappa       , nxy , nz , n_group_d )
         call alloc (omegam      , n_group_d , nxy , nz )
         call alloc (omega0      , n_group_d , nxy , nz )
         call alloc (omegap      , n_group_d , nxy , nz )
         call alloc (u_omega     , n_group   , nxy , nz )
         call alloc (flux_extra  , n_group   , nxy , nz )
         call alloc (src_tr_fac1 , n_group   , nxy , nz )
         call alloc (src_tr_fac2 , n_group   , nxy , nz )

         call alloc (iPgenT  ,3 )
         call alloc (rhocoef ,6 )
         call alloc (Pbeta   ,n_group_d )
         call alloc (Plambda ,n_group_d , 3 )
         call alloc (Pzeta   ,n_group_d , 3 )
      endif
#ifdef tuan_tr_test

      call alloc ( rldt      , nxy_bk , nz , n_group_d )
      call alloc ( rldtgp1   , nxy_bk , nz , n_group_d )
      call alloc ( capbrldt  , nxy_bk , nz , n_group_d )
      call alloc ( cappap1   , nxy_bk , nz , n_group_d )
      call alloc ( capbrldt2 , nxy_bk , nz , n_group_d )
      call alloc ( rvdtvol   ,               n_group   )
      call alloc ( rvdtdvol  ,               n_group   )
#endif

      return
      end subroutine alloc_tran


      subroutine alloc_th
      use mod_alloc
      use inc_th

#ifdef js_mpi
      use inc_parallel, only: comm
      use inc_parallel, only: iproc
      use inc_parallel, only: iproc2nxy_fa
#endif
      implicit none
      integer :: nchan_bk


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [alloc_th] in Mod_Manage'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call alloc (r , nr+5 )

      nchan_bk=nchan
#ifdef js_mpi
      if (comm%usempi) then
         nchan_bk=iproc2nxy_fa(iproc)
      endif
#endif

      call alloc (htcoef        ,   nzth ,   nchan_bk )
      call alloc (hcool         ,   nzth ,   nchan_bk+1 )
      call alloc (tfuel  , nr+5 ,   nzth ,   nchan_bk )
      call alloc0(tdopl         , 0,nzth , 0,nchan_bk+1 )
      call alloc0(tcool         , 0,nzth , 0,nchan_bk+1 )
      call alloc0(dcool         , 0,nzth , 0,nchan_bk+1 )
      call alloc0(rhou          , 0,nzth , 1,nchan_bk   )
      call alloc0(rhohu         , 0,nzth , 1,nchan_bk   )
      call alloc0(u             , 0,nzth , 1,nchan_bk   )
      call alloc0(ud            , 0,nzth , 1,nchan_bk   )

      call alloc (qflux, nzth , nchan_bk )
      call alloc (qvol , nzth , nchan_bk )
      call alloc (qeff , nzth , nchan_bk )

      return
      end subroutine alloc_th

#ifdef siarhei_delete 
      subroutine alloc_pin
      use inc_option, only: n_group
      use inc_pinvar
      use inc_pinpow
      use mod_alloc
#ifdef JR_SRCTRM
      USE Inc_SNF
      USE Inc_depletion, only: N_BU
      USE Inc_Nuclide , only: N_Gd
#endif
      implicit none

      call alloc (pploc   , nxy_1n    )
      call alloc (nxcs    , ny+1      )
      call alloc (nxce    , ny+1      )
      call alloc (nneighc , ncorn     )
      call alloc (lcnw    , nxy       )
      call alloc (lcsw    , nxy       )
      call alloc (lcne    , nxy       )
      call alloc (lcse    , nxy       )
      call alloc0(lctox   , 0,ncorn   )
      call alloc0(lctoy   , 0,ncorn   )
      call alloc (lcc     , 8 , ncorn )
      call alloc (lcn     , 4 , ncorn )

      call alloc (asnmut   , nxy , nz )
      call alloc (asnmutl  , nxy , nz )
      call alloc (asnmuttl , nxy , nz )
      call alloc (acnmut   , nxy , nz )
      call alloc (acnmutl  , nxy , nz )
      call alloc (acnmuttl , nxy , nz )
      call alloc (asnkat   , nxy , nz )
      call alloc (asnkatl  , nxy , nz )
      call alloc (asnkattl , nxy , nz )
      call alloc (acnkat   , nxy , nz )
      call alloc (acnkatl  , nxy,  nz )
      call alloc (acnkattl , nxy,  nz )
      call alloc (bctka    , nxy , nz )
      call alloc (bctmu    , nxy , nz )
      call alloc0(powvalr  , 1,npin , 1,npin , 0,nz )
      call alloc (powval   , npin, npin, nz, N_Group)
      call alloc (RevNorm  , Nx_1N, Ny_1N)
      call alloc (pinpitch , nxy , nz)
      call alloc (xoffset  , nxy , nz)
      call alloc (yoffset  , nxy , nz)
      call alloc0(nodec, 0, Nx+2, 0, ny+2)
      call alloc0(pppeak, 1, Nxy_1N, 0, nz)
      call alloc0(Peak_X, 1, Nxy_1N, 0, nz)
      call alloc0(Peak_Y, 1, Nxy_1N, 0, nz)

      call alloc0(cc11   ,0,8,1,ncorn,1,nz)
      call alloc0(cc12   ,0,8,1,ncorn,1,nz)
      call alloc0(cc21   ,0,8,1,ncorn,1,nz)
      call alloc0(cc22   ,0,8,1,ncorn,1,nz)
      call alloc0(cornfka,0,2,1,nxy,1,nz)
      call alloc0(cornfmu,0,2,1,nxy,1,nz)
      call alloc (cur11  ,2,nxy,nz)
      call alloc (cur12  ,2,nxy,nz)
      call alloc (cur21  ,2,nxy,nz)
      call alloc (cur22  ,2,nxy,nz)
      call alloc (cpbsrc ,N_Group,ncorn)
      call alloc (phicorn,N_Group,ncorn,nz)
      call alloc (phihom,npin,npin,Nxy_1N,nz,N_Group)

#ifdef JR_SRCTRM
      if(flag_snf_pin) then
         CALL Alloc( phihom_M_old,npin,npin,Nxy_1N,nz,N_Group)
         CALL Alloc( phihom_M,npin,npin,Nxy_1N,nz,N_Group)
         CALL Alloc( phihom_SNF,N_BU,npin,npin,N_Group)
         CALL Alloc(SNF_N_U36  ,int(npin),int(npin) )
         CALL Alloc(SNF_N_U37  ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Np37 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Np38 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Np39 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu38 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu39 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu40 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu41 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu42 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu43 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Am41 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_As42 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Am42 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Am43 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Am44 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Cm42 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Cm43 ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Cm44 ,int(npin),int(npin) )

         CALL Alloc(SNF_N_U36_Old  ,int(npin),int(npin) )
         CALL Alloc(SNF_N_U37_Old  ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Np37_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Np38_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Np39_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu38_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu39_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu40_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu41_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu42_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Pu43_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Am41_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_As42_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Am42_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Am43_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Am44_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Cm42_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Cm43_Old ,int(npin),int(npin) )
         CALL Alloc(SNF_N_Cm44_Old ,int(npin),int(npin) )

         CALL Alloc( SNF_N_I35  ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Xe35 ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Nd47 ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Nd48 ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Nd49 ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Pm47 ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Ps48 ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Pm48 ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Pm49 ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Sm47 ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Sm48 ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Sm49 ,int(npin),int(npin) )

         CALL Alloc( SNF_N_I35_Old  ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Xe35_Old ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Nd47_Old ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Nd48_Old ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Nd49_Old ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Pm47_Old ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Ps48_Old ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Pm48_Old ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Pm49_Old ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Sm47_Old ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Sm48_Old ,int(npin),int(npin) )
         CALL Alloc( SNF_N_Sm49_Old ,int(npin),int(npin) )

      endif
#endif

      call alloc (akratio,nxy,nz)
      call alloc (akappa ,nxy,nz)
      call alloc (amu    ,nxy,nz)
      call alloc (ar     ,nxy,nz)
      call alloc (as     ,nxy,nz)
      call alloc (ardet  ,nxy,nz)
      call alloc (asnmu  ,nxy,nz)
      call alloc (acnmu  ,nxy,nz)
      call alloc (asnka  ,nxy,nz)
      call alloc (acnka  ,nxy,nz)
      call alloc (asnmux ,nxy,nz)
      call alloc (acnmux ,nxy,nz)
      call alloc (asnkax ,nxy,nz)
      call alloc (acnkax ,nxy,nz)
      call alloc (asnmuy ,nxy,nz)
      call alloc (acnmuy ,nxy,nz)
      call alloc (asnkay ,nxy,nz)
      call alloc (acnkay ,nxy,nz)
      call alloc (asnmuz ,nxy,nz)
      call alloc (acnmuz ,nxy,nz)
      call alloc (asnkaz ,nxy,nz)
      call alloc (acnkaz ,nxy,nz)

      call alloc (kflag  ,nxy,nz)

      end subroutine alloc_pin
#endif 


#ifdef tuan_tr_test
      SUBROUTINE Alloc_FBRHO

      USE Inc_Lscoef
      USE Inc_Option, ONLY: N_Group
      USE Mod_Alloc
      USE Inc_TH
      USE Inc_Adjoint
      Use Inc_CR

      IMPLICIT NONE

      call alloc(t_Fuel0, nxy,nz)
      call alloc(t_Mod0,  nxy,nz)
      call alloc(d_Mod0,  nxy,nz)
      call alloc(CR_Bot0, n_cr)
      call alloc(CR_Top0, n_cr)

      call alloc(init_miXS_a_Xe35,nxy,nz,n_group)
      call alloc(init_miXS_a_Sm49,nxy,nz,n_group)
      call alloc(init_N_Xe35,nxy,nz)
      call alloc(init_N_Sm49,nxy,nz)

      CALL alloc(dnw0, N_Group,nxy+nx+ny,nz)   ! initial x-drection d-hat
      CALL alloc(dnn0, N_Group,nxy+nx+ny,nz)   ! initial y-drection d-hat
      CALL alloc(dnb0, N_Group,nxy,nz+1)       ! initial z-drection d-hat
      CALL alloc(xsxea0, N_Group,nxy,nz)       ! initial Xe micro
      CALL alloc(xssma0, N_Group,nxy,nz)       ! initial Sm micro
      CALL alloc(rnxe0, nxy,nz)                ! initial Xe number density
      CALL alloc(rnsm0, nxy,nz)                ! initial Sm number density

      END SUBROUTINE Alloc_FBRHO
      
#elif siarhei_delete 
      SUBROUTINE Alloc_FBRHO

      USE Inc_Lscoef
      USE Inc_Option, ONLY: N_Group
      USE Mod_Alloc
      USE Inc_TH
      USE Inc_Adjoint
      Use Inc_CR

      IMPLICIT NONE

      call alloc(init_tFuel, nxy,nz)
      call alloc(init_tMod,  nxy,nz)
      call alloc(init_dMod,  nxy,nz)
      call alloc(init_CR_Bot, n_cr)
      call alloc(init_CR_Top, n_cr)

      call alloc(init_miXS_a_Xe35,nxy,nz,n_group)
      call alloc(init_miXS_a_Sm49,nxy,nz,n_group)
      call alloc(init_N_Xe35,nxy,nz)
      call alloc(init_N_Sm49,nxy,nz)

      CALL alloc(dnw0, N_Group,nxy+nx+ny,nz)   ! initial x-drection d-hat
      CALL alloc(dnn0, N_Group,nxy+nx+ny,nz)   ! initial y-drection d-hat
      CALL alloc(dnb0, N_Group,nxy,nz+1)       ! initial z-drection d-hat
      CALL alloc(xsxea0, N_Group,nxy,nz)       ! initial Xe micro
      CALL alloc(xssma0, N_Group,nxy,nz)       ! initial Sm micro
      CALL alloc(rnxe0, nxy,nz)                ! initial Xe number density
      CALL alloc(rnsm0, nxy,nz)                ! initial Sm number density

      END SUBROUTINE Alloc_FBRHO
#endif 

      END MODULE Mod_Manage

