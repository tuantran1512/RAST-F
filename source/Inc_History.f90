
      MODULE Inc_History

      IMPLICIT NONE

      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Sec
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Hour
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Day
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Time

      real(8), dimension(:), allocatable :: Hist_Cycle_BU
      real(8), dimension(:), allocatable :: Hist_Accum_BU
      real(8), dimension(:), allocatable :: Hist_PPower
      real(8), dimension(:), allocatable :: Hist_PPM
      real(8), dimension(:), allocatable :: Hist_keff
      real(8), dimension(:), allocatable :: Hist_Reactivity
      real(8), dimension(:), allocatable :: hist_tmm_inp
      real(8), dimension(:), allocatable :: Hist_MaxFenthl_Rise
      real(8), dimension(:), allocatable :: Hist_MaxFenthl

      REAL(8), DIMENSION(:), ALLOCATABLE ::  Hist_ADJPPM
      REAL(8), DIMENSION(:), ALLOCATABLE ::  Hist_PPMB10
      REAL(8), DIMENSION(:), ALLOCATABLE ::  Hist_B10DEP

      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_ASI
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_PF_1D
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_PF_2D
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_PF_3D
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_PBU_FA
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_PBU_Pin

      REAL(8), DIMENSION(:), ALLOCATABLE ::  Hist_Fxy
      REAL(8), DIMENSION(:), ALLOCATABLE ::  Hist_Fq
      REAL(8), DIMENSION(:), ALLOCATABLE ::  Hist_Fr

      real(8), dimension(:), allocatable :: Hist_Fq_val
      real(8), dimension(:), allocatable :: Hist_Fr_val
      real(8), dimension(:), allocatable :: Hist_Fz_val
      real(8), dimension(:), allocatable :: Hist_FdH_val
      real(8), dimension(:), allocatable :: Hist_max_Fxy_val

      real(8), dimension(:),     allocatable :: Hist_LinPOW_ave
      real(8), dimension(:),     allocatable :: Hist_LinPOW_max
      real(8), dimension(:),     allocatable :: Hist_Fz_ave
      real(8), dimension(:,:,:), allocatable :: Hist_Fxy_val
      real(8), dimension(:,:,:), allocatable :: Hist_FdH_val_XY
      real(8), dimension(:),     allocatable :: Hist_kinf
      real(8), dimension(:),     allocatable :: HIST_a_1
      real(8), dimension(:),     allocatable :: HIST_a_2
      real(8), dimension(:),     allocatable :: HIST_D_1
      real(8), dimension(:),     allocatable :: HIST_D_2
      real(8), dimension(:),     allocatable :: HIST_nf_1
      real(8), dimension(:),     allocatable :: HIST_nf_2
      real(8), dimension(:),     allocatable :: HIST_r_1
      real(8), dimension(:),     allocatable :: HIST_kf_1
      real(8), dimension(:),     allocatable :: HIST_kf_2
      real(8), dimension(:),     allocatable :: HIST_f_1
      real(8), dimension(:),     allocatable :: HIST_f_2
      real(8), dimension(:,:),   allocatable :: Hist_FA_kinf
      real(8), dimension(:,:,:), allocatable :: Hist_FA_kinf_XY
      real(8), dimension(:,:),   allocatable :: Hist_FA_kinf_4N
      real(8), dimension(:,:,:), allocatable :: Hist_FA_kinf_XY_4N
      real(8), dimension(:,:,:), allocatable :: HIST_BU_XYZ_4N_ANC

      integer, dimension(:,:), allocatable :: Hist_Fq_loc
      integer, dimension(:,:), allocatable :: Hist_Fr_loc
      integer, dimension(:,:), allocatable :: Hist_Fz_loc
      integer, dimension(:,:), allocatable :: Hist_FdH_loc
      integer, dimension(:,:), allocatable :: Hist_max_Fxy_loc

      real(8), dimension(:), allocatable :: Hist_maxPBU_2D
      real(8), dimension(:), allocatable :: Hist_maxPPD_2D
      real(8), dimension(:), allocatable :: Hist_maxPPD_3D
      real(8), dimension(:), allocatable :: Hist_maxPPW_2D
      real(8), dimension(:), allocatable :: Hist_maxPPW_3D
      integer, dimension(:,:), allocatable :: Hist_PBU_2D_loc
      integer, dimension(:,:), allocatable :: Hist_PPD_2D_loc
      integer, dimension(:,:), allocatable :: Hist_PPD_3D_loc
      integer, dimension(:,:), allocatable :: Hist_PPW_2D_loc
      integer, dimension(:,:), allocatable :: Hist_PPW_3D_loc

      INTEGER, DIMENSION(:), ALLOCATABLE :: Hist_PP_1D
      INTEGER, DIMENSION(:), ALLOCATABLE :: Hist_PP_2D
      INTEGER, DIMENSION(:), ALLOCATABLE :: Hist_PP_3D
      INTEGER, DIMENSION(:), ALLOCATABLE :: Hist_PPBU_FA
      INTEGER, DIMENSION(:), ALLOCATABLE :: Hist_PPBU_Pin

      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_CR

      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_TF_Avg
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_TF_Max
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_TF_Min
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_TM_Avg
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_TM_Max
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_TM_Min
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_TM_In
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_TM_Out
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_DM_Avg
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_DM_Max
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_DM_Min
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_TFcen_Max
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Xe35
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Xe35_ASI
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Sm49
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Sm49_ASI
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_I35
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Nd47
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Nd48
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Nd49
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Pm47
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Ps48
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Pm48
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Pm49
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Sm47
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Sm48
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Gd52
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Gd54
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Gd55
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Gd56
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Gd57
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Gd58
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Gd60

      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_U34
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_U35
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_U36
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_U37
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_U38
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Np37
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Np38
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Np39
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Pu38
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Pu39
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Pu40
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Pu41
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Pu42
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Pu43
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Am41
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_As42
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Am42
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Am43
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Am44
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Cm42
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Cm43
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_Cm44

      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Hist_Normal_Power_Z
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_Normal_Power_XY
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_Normal_Power_XYZ
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Hist_BU_Z
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_BU_XY
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_BU_XYZ
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Hist_T_Fuel_Z
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Hist_T_Mod_Z
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Hist_D_Mod_Z
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_T_Fuel_XY
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_T_Mod_XY
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_D_Mod_XY
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_T_Fuel_XYZ
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_T_Mod_XYZ
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_D_Mod_XYZ

      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_N_Xe35_XYZ
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_N_Sm49_XYZ
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_FFlux_XYZ
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_TFlux_XYZ
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_Normal_Power_XYZ_4N
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_BU_XYZ_4N
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_T_Fuel_XYZ_4N
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_T_Mod_XYZ_4N
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_N_Xe35_XYZ_4N
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_N_Sm49_XYZ_4N
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_FFlux_XYZ_4N
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_TFlux_XYZ_4N

      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_Normal_Power_XY_4N
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_BU_XY_4N

      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Hist_crud_nden_Z
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Hist_cruddr_Z
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Hist_cruddtf_Z
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_crud_nden_XY
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_cruddr_XY
      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_cruddtf_XY
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_crud_nden_XYZ
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_cruddr_XYZ
      REAL(8), DIMENSION(:, :, :, :), ALLOCATABLE :: Hist_cruddtf_XYZ

      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Hist_FA_ASI
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_FTC
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_MTC
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_ITC
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_CRW
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_XEW
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_SMW
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_GDW

      real(8), dimension(:,:,:), allocatable :: Hist_Power
      real(8), dimension(:,:,:), allocatable :: Hist_Burnup
      real(8), dimension(:,:,:), allocatable :: Hist_T_Fuel
      real(8), dimension(:,:,:), allocatable :: Hist_T_Mod
      real(8), dimension(:,:,:), allocatable :: Hist_D_Mod
      real(8), dimension(:,:,:), allocatable :: Hist_N_I35
      real(8), dimension(:,:,:), allocatable :: Hist_N_Xe35
      real(8), dimension(:,:,:), allocatable :: Hist_N_Nd47
      real(8), dimension(:,:,:), allocatable :: Hist_N_Nd48
      real(8), dimension(:,:,:), allocatable :: Hist_N_Nd49
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pm47
      real(8), dimension(:,:,:), allocatable :: Hist_N_Ps48
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pm48
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pm49
      real(8), dimension(:,:,:), allocatable :: Hist_N_Sm47
      real(8), dimension(:,:,:), allocatable :: Hist_N_Sm48
      real(8), dimension(:,:,:), allocatable :: Hist_N_Sm49
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd52
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd54
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd55
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd56
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd57
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd58
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd60
      real(8), dimension(:,:,:), allocatable :: Hist_N_U34
      real(8), dimension(:,:,:), allocatable :: Hist_N_U35
      real(8), dimension(:,:,:), allocatable :: Hist_N_U36
      real(8), dimension(:,:,:), allocatable :: Hist_N_U37
      real(8), dimension(:,:,:), allocatable :: Hist_N_U38
      real(8), dimension(:,:,:), allocatable :: Hist_N_Np37
      real(8), dimension(:,:,:), allocatable :: Hist_N_Np38
      real(8), dimension(:,:,:), allocatable :: Hist_N_Np39
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu38
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu39
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu40
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu41
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu42
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu43
      real(8), dimension(:,:,:), allocatable :: Hist_N_Am41
      real(8), dimension(:,:,:), allocatable :: Hist_N_As42
      real(8), dimension(:,:,:), allocatable :: Hist_N_Am42
      real(8), dimension(:,:,:), allocatable :: Hist_N_Am43
      real(8), dimension(:,:,:), allocatable :: Hist_N_Am44
      real(8), dimension(:,:,:), allocatable :: Hist_N_Cm42
      real(8), dimension(:,:,:), allocatable :: Hist_N_Cm43
      real(8), dimension(:,:,:), allocatable :: Hist_N_Cm44

      real(8), dimension(:,:,:), allocatable :: Hist_Power_FA
      real(8), dimension(:,:,:), allocatable :: Hist_Burnup_FA
      real(8), dimension(:,:,:), allocatable :: Hist_T_Fuel_FA
      real(8), dimension(:,:,:), allocatable :: Hist_T_Mod_FA
      real(8), dimension(:,:,:), allocatable :: Hist_D_Mod_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_I35_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Xe35_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Nd47_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Nd48_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Nd49_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pm47_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Ps48_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pm48_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pm49_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Sm47_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Sm48_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Sm49_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd52_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd54_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd55_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd56_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd57_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd58_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Gd60_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_U34_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_U35_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_U36_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_U37_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_U38_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Np37_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Np38_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Np39_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu38_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu39_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu40_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu41_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu42_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Pu43_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Am41_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_As42_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Am42_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Am43_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Am44_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Cm42_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Cm43_FA
      real(8), dimension(:,:,:), allocatable :: Hist_N_Cm44_FA

      real(8), dimension(:,:,:,:), allocatable :: Hist_Flux
      real(8), dimension(:,:,:,:), allocatable :: Hist_maXS_tr_3D
      real(8), dimension(:,:,:,:), allocatable :: Hist_D_3D
      real(8), dimension(:,:,:,:), allocatable :: Hist_maXS_a_3D
      real(8), dimension(:,:,:,:), allocatable :: Hist_maXS_f_3D
      real(8), dimension(:,:,:,:), allocatable :: Hist_nu_maXS_f_3D
      real(8), dimension(:,:,:,:), allocatable :: Hist_kap_maXS_f_3D
      real(8), dimension(:,:,:,:), allocatable :: Hist_maXS_s_3D
      real(8), dimension(:,:,:,:), allocatable :: Hist_maXS_r_3D
      real(8), dimension(:,:,:,:), allocatable :: Hist_ADF_Lx
      real(8), dimension(:,:,:,:), allocatable :: Hist_ADF_Rx
      real(8), dimension(:,:,:,:), allocatable :: Hist_ADF_Ly
      real(8), dimension(:,:,:,:), allocatable :: Hist_ADF_Ry
      real(8), dimension(:,:,:,:), allocatable :: Hist_CDF_LxLy
      real(8), dimension(:,:,:,:), allocatable :: Hist_CDF_LxRy
      real(8), dimension(:,:,:,:), allocatable :: Hist_CDF_RxLy
      real(8), dimension(:,:,:,:), allocatable :: Hist_CDF_RxRy

      REAL(8), DIMENSION(:, :), ALLOCATABLE :: Hist_DR
      REAL(8), DIMENSION(:), ALLOCATABLE :: Hist_ASI_DR
      real(8), dimension(:,:,:), allocatable :: Hist_det_sig_rr
      real(8), dimension(:,:,:), allocatable :: Hist_det_detpow
      real(8), dimension(:,:,:), allocatable :: Hist_det_wprime

      REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Hist_PinPOW_3D
      REAL(8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Hist_PinBU_3D
      REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Hist_PinPOW_2D
      REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Hist_PinBU_2D
      REAL(8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Hist_PinQ_2D

      real(8), dimension(:), allocatable :: hist_6factor_eps
      real(8), dimension(:), allocatable :: hist_6factor_p
      real(8), dimension(:), allocatable :: hist_6factor_eta
      real(8), dimension(:), allocatable :: hist_6factor_f
      real(8), dimension(:), allocatable :: hist_6factor_fnl
      real(8), dimension(:), allocatable :: hist_6factor_tnl
      real(8), dimension(:), allocatable :: hist_6factor_kinf
      real(8), dimension(:), allocatable :: hist_core_avg_nu
      real(8), dimension(:), allocatable :: hist_mfp
      real(8), dimension(:), allocatable :: hist_mfp_f
      real(8), dimension(:), allocatable :: hist_mfp_t
      real(8), dimension(:), allocatable :: hist_n_density
      real(8), dimension(:), allocatable :: hist_n_density_f
      real(8), dimension(:), allocatable :: hist_n_density_t
      real(8), dimension(:), allocatable :: hist_n_flux
      real(8), dimension(:), allocatable :: hist_n_flux_f
      real(8), dimension(:), allocatable :: hist_n_flux_t
      real(8), dimension(:), allocatable :: hist_n_speed
      real(8), dimension(:), allocatable :: hist_n_speed_f
      real(8), dimension(:), allocatable :: hist_n_speed_t
      real(8), dimension(:), allocatable :: hist_spectral_idx
      real(8), dimension(:), allocatable :: hist_dif_length
      real(8), dimension(:), allocatable :: hist_dif_length_f
      real(8), dimension(:), allocatable :: hist_dif_length_t
      real(8), dimension(:), allocatable :: hist_mig_length
      real(8), dimension(:), allocatable :: hist_n_lifetime
      real(8), dimension(:), allocatable :: hist_n_lifetime_f
      real(8), dimension(:), allocatable :: hist_n_lifetime_t
      real(8), dimension(:), allocatable :: hist_n_gentime
      real(8), dimension(:), allocatable :: hist_n_gentime_f
      real(8), dimension(:), allocatable :: hist_n_gentime_t
      real(8), dimension(:), allocatable :: hist_n_energy
      real(8), dimension(:), allocatable :: hist_n_energy_f
      real(8), dimension(:), allocatable :: hist_n_energy_t

      real(8), allocatable :: hist_qchf(:,:,:)
      real(8), allocatable :: hist_qwat(:,:,:)
      real(8), allocatable :: hist_dnbr(:,:,:)
      real(8), allocatable :: hist_mqchf(:)
      real(8), allocatable :: hist_mqwat(:)
      real(8), allocatable :: hist_mdnbr(:)

      real(8), allocatable :: hist_kp_beta(:,:)
      real(8), allocatable :: hist_kp_lambda(:,:)
      real(8), allocatable :: hist_kp_zeta(:,:)
      real(8), allocatable :: hist_kp_gentime(:)

      END MODULE Inc_History

