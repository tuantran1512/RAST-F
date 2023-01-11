
      MODULE Mod_Save

      USE Inc_Constant
      USE Inc_Geometry
      USE Inc_RP
      USE Inc_Option
      use Mod_GetNode, only: new_asym_itab

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

      CONTAINS

!#ifdef siarhei_delete 
!#ifdef tuan_tr_test
      subroutine Save_Old
      use inc_solver, only: betap, betap_Old, rvdelt, rvdelt_Old
      use inc_3d
      use inc_depletion
      use inc_kinetics, only: N_Group_d
      use inc_maxs, only: maXS_f_3D, maXS_f_3D_Old
      use inc_mixs
      use inc_nuclide
      use inc_th, only: PPower, PPower_Old
      use inc_time
      use inc_transient
      use inc_xs_file, only: I_Tab, I_BU
      use inc_detector
      use inc_flag, only: flag_InDet
      use inc_crud, only: opt_crud, ind_bu, ind_bu_old
#ifdef JR_SRCTRM
      use Inc_Pinvar, only: npin
      use Inc_SNF
#endif
#ifdef js_mpc
      use inc_pinpow, only: pinq_3d, pinq_3d_old
      use inc_multiphysics, only: flag_mpc, flag_mp_oneway
#endif
#ifdef js_mpi
      use inc_parallel, only: comm
      use inc_parallel, only: iproc
      use inc_parallel, only: ixyip2ixy
#endif

#ifdef tuan_fr
      use Inc_miXS
      use Inc_TPEN, ONLY: fluxf,fluxf_old
      use inc_flag, only: FLag_miXS
#endif
      implicit none
      integer :: ixy, iz, ig
      integer :: n, m
      integer :: l0
#ifdef JR_SRCTRM
      INTEGER :: ip, jp
#endif

#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

#ifdef JR_SRCTRM
      IF(Flag_SNF_PIN) THEN
         DO ip=1,int(npin)
            DO jp=1,int(npin)
               phihom_M_old(ip,jp,SNF_Ixy_1N,SNF_k,:) = phihom_M(ip,jp,SNF_Ixy_1N,SNF_k,:)
               SNF_N_U34_Old (ip,jp) = SNF_N_U34 (ip,jp)
               SNF_N_U35_Old (ip,jp) = SNF_N_U35 (ip,jp)
               SNF_N_U36_Old (ip,jp) = SNF_N_U36 (ip,jp)
               SNF_N_U37_Old (ip,jp) = SNF_N_U37 (ip,jp)
               SNF_N_U38_Old (ip,jp) = SNF_N_U38 (ip,jp)
               SNF_N_Np37_Old(ip,jp) = SNF_N_Np37(ip,jp)
               SNF_N_Np38_Old(ip,jp) = SNF_N_Np38(ip,jp)
               SNF_N_Np39_Old(ip,jp) = SNF_N_Np39(ip,jp)
               SNF_N_Pu38_Old(ip,jp) = SNF_N_Pu38(ip,jp)
               SNF_N_Pu39_Old(ip,jp) = SNF_N_Pu39(ip,jp)
               SNF_N_Pu40_Old(ip,jp) = SNF_N_Pu40(ip,jp)
               SNF_N_Pu41_Old(ip,jp) = SNF_N_Pu41(ip,jp)
               SNF_N_Pu42_Old(ip,jp) = SNF_N_Pu42(ip,jp)
               SNF_N_Pu43_Old(ip,jp) = SNF_N_Pu43(ip,jp)
               SNF_N_Am41_Old(ip,jp) = SNF_N_Am41(ip,jp)
               SNF_N_As42_Old(ip,jp) = SNF_N_As42(ip,jp)
               SNF_N_Am42_Old(ip,jp) = SNF_N_Am42(ip,jp)
               SNF_N_Am43_Old(ip,jp) = SNF_N_Am43(ip,jp)
               SNF_N_Am44_Old(ip,jp) = SNF_N_Am44(ip,jp)
               SNF_N_Cm42_Old(ip,jp) = SNF_N_Cm42(ip,jp)
               SNF_N_Cm43_Old(ip,jp) = SNF_N_Cm43(ip,jp)
               SNF_N_Cm44_Old(ip,jp) = SNF_N_Cm44(ip,jp)

               SNF_N_I35_Old (ip,jp) = SNF_N_I35 (ip,jp)
               SNF_N_Xe35_Old(ip,jp) = SNF_N_Xe35(ip,jp)
               SNF_N_Nd47_Old(ip,jp) = SNF_N_Nd47(ip,jp)
               SNF_N_Nd48_Old(ip,jp) = SNF_N_Nd48(ip,jp)
               SNF_N_Nd49_Old(ip,jp) = SNF_N_Nd49(ip,jp)
               SNF_N_Pm47_Old(ip,jp) = SNF_N_Pm47(ip,jp)
               SNF_N_Ps48_Old(ip,jp) = SNF_N_Ps48(ip,jp)
               SNF_N_Pm48_Old(ip,jp) = SNF_N_Pm48(ip,jp)
               SNF_N_Pm49_Old(ip,jp) = SNF_N_Pm49(ip,jp)
               SNF_N_Sm47_Old(ip,jp) = SNF_N_Sm47(ip,jp)
               SNF_N_Sm48_Old(ip,jp) = SNF_N_Sm48(ip,jp)
               SNF_N_Sm49_Old(ip,jp) = SNF_N_Sm49(ip,jp)

               I_Tab = AxialComp( I_LP_1N( SNF_Ixy_1N ), SNF_k)
               I_Tab = new_asym_itab(I_Tab,I_1Nto4N(SNF_Ixy_1N,1))

               SNF_N_Gd52_Old(ip,jp) = SNF_N_Gd52(ip,jp)
               SNF_N_Gd54_Old(ip,jp) = SNF_N_Gd54(ip,jp)
               SNF_N_Gd55_Old(ip,jp) = SNF_N_Gd55(ip,jp)
               SNF_N_Gd56_Old(ip,jp) = SNF_N_Gd56(ip,jp)
               SNF_N_Gd57_Old(ip,jp) = SNF_N_Gd57(ip,jp)
               SNF_N_Gd58_Old(ip,jp) = SNF_N_Gd58(ip,jp)
               SNF_N_Gd60_Old(ip,jp) = SNF_N_Gd60(ip,jp)

            ENDDO
         ENDDO
      ENDIF
#endif

      PPower_old = PPower
      dT_old = dT
      do ixy=1,nxy
         do iz=1,nz
            i_tab=axialcomp(i_lp_1n(i_4nto1n(ixy)),iz)
            i_tab=new_asym_itab(i_tab,ixy)
            do ig=1,n_group
               flux_old(ixy,iz,ig)=flux(ixy,iz,ig)
#ifdef tuan_fr
            if (.not.flag_transient) then 
                  fluxf_old(ixy,iz,ig) = fluxf(ixy,iz,ig)
            endif
#endif
            enddo
            if (opt_mode==3) then
               do ig=1,n_group
                  if (.not.flag_transient.and.I_BU==1) then
                     flux_oold(ixy,iz,ig)=flux(ixy,iz,ig)
                  else
                     flux_oold(ixy,iz,ig)=flux_old(ixy,iz,ig)
                  endif
               enddo

               if (.not.flag_transient.and.I_BU==1) then
                  fissrc_oold(ixy,iz)=fissrc(ixy,iz)
               else
                  fissrc_oold(ixy,iz)=fissrc_old(ixy,iz)
               endif
               fissrc_old(ixy,iz)=fissrc(ixy,iz)
            endif

            BU_old(ixy,iz)=BU(ixy,iz)
 

            l0=ixy
#ifdef js_mpi
            if (comm%usempi) then
               l0=ixyip2ixy(ixy,iproc)
            endif
#endif
            if ((flag_transient).or.(.not.flag_transient.and.I_BU==1)) then
               do ig=1,n_group
                  rvdelt_old(ig,l0,iz)=rvdelt(ig,l0,iz)
               enddo
               betap_old(l0,iz)=betap(l0,iz)
            endif

            if (flag_transient) then
               do ig=1,N_Group_d
                  C_d_I_old(ixy,iz,ig)=C_d_I(ixy,iz,ig)
               enddo
            endif

            if (flag_InDet .and. opt_mode>4) then
               N_A1_old  (ixy,iz)= N_A1  (ixy,iz)
               N_A2_old  (ixy,iz)= N_A2  (ixy,iz)
               N_A2m_old (ixy,iz)= N_A2m (ixy,iz)
               N_A3_old  (ixy,iz)= N_A3  (ixy,iz)
               N_B2_old  (ixy,iz)= N_B2  (ixy,iz)
               N_B3_old  (ixy,iz)= N_B3  (ixy,iz)
               N_V1_old  (ixy,iz)= N_V1  (ixy,iz)
               N_V2_old  (ixy,iz)= N_V2  (ixy,iz)
               N_Cr52_old(ixy,iz)= N_Cr52(ixy,iz)
               N_A4_old  (ixy,iz)= N_A4  (ixy,iz)
            endif
         enddo
      enddo
         
      if (flag_InDet .and. opt_mode>4) then
         do m=1,DET_N_Z
            do n=1,DET_N_XY
               DET_N_A1_old   (n,m) = DET_N_A1   (n,m)
               DET_N_A2_old   (n,m) = DET_N_A2   (n,m)
               DET_N_A2m_old  (n,m) = DET_N_A2m  (n,m)
               DET_N_A3_old   (n,m) = DET_N_A3   (n,m)
               DET_N_B2_old   (n,m) = DET_N_B2   (n,m)
               DET_N_B3_old   (n,m) = DET_N_B3   (n,m)
               DET_N_V1_old   (n,m) = DET_N_V1   (n,m)
               DET_N_V2_old   (n,m) = DET_N_V2   (n,m)
               DET_N_Cr52_old (n,m) = DET_N_Cr52 (n,m)
               DET_N_A4_old   (n,m) = DET_N_A4   (n,m)
            enddo
         enddo
      endif

#ifdef js_mpc
      if (flag_mpc.or.flag_mp_oneway) then
         pinq_3d_old(:,:,:,:,:)=pinq_3d(:,:,:,:,:)
      endif
#endif

      N_U34_old (:,:) = N_U34 (:,:)
      N_U35_old (:,:) = N_U35 (:,:)
      N_U36_old (:,:) = N_U36 (:,:)
      N_U36_old (:,:) = N_U36 (:,:)
      N_U37_old (:,:) = N_U37 (:,:)
      N_U38_old (:,:) = N_U38 (:,:)
      N_Np37_old(:,:) = N_Np37(:,:)
      N_Np38_old(:,:) = N_Np38(:,:)
      N_Np39_old(:,:) = N_Np39(:,:)
      N_Pu38_old(:,:) = N_Pu38(:,:)
      N_Pu39_old(:,:) = N_Pu39(:,:)
      N_Pu40_old(:,:) = N_Pu40(:,:)
      N_Pu41_old(:,:) = N_Pu41(:,:)
      N_Pu42_old(:,:) = N_Pu42(:,:)
      N_Pu43_old(:,:) = N_Pu43(:,:)
      N_Am41_old(:,:) = N_Am41(:,:)
      N_As42_old(:,:) = N_As42(:,:)
      N_Am42_old(:,:) = N_Am42(:,:)
      N_Am43_old(:,:) = N_Am43(:,:)
      N_Am44_old(:,:) = N_Am44(:,:)
      N_Cm42_old(:,:) = N_Cm42(:,:)
      N_Cm43_old(:,:) = N_Cm43(:,:)
      N_Cm44_old(:,:) = N_Cm44(:,:)
      N_I35_old (:,:) = N_I35 (:,:)
      N_Xe35_old(:,:) = N_Xe35(:,:)
      N_Nd47_old(:,:) = N_Nd47(:,:)
      N_Nd48_old(:,:) = N_Nd48(:,:)
      N_Nd49_old(:,:) = N_Nd49(:,:)
      N_Pm47_old(:,:) = N_Pm47(:,:)
      N_Ps48_old(:,:) = N_Ps48(:,:)
      N_Pm48_old(:,:) = N_Pm48(:,:)
      N_Pm49_old(:,:) = N_Pm49(:,:)
      N_Sm47_old(:,:) = N_Sm47(:,:)
      N_Sm48_old(:,:) = N_Sm48(:,:)
      N_Sm49_old(:,:) = N_Sm49(:,:)
      N_Gd52_old(:,:) = N_Gd52(:,:)
      N_Gd54_old(:,:) = N_Gd54(:,:)
      N_Gd55_old(:,:) = N_Gd55(:,:)
      N_Gd56_old(:,:) = N_Gd56(:,:)
      N_Gd57_old(:,:) = N_Gd57(:,:)
      N_Gd58_old(:,:) = N_Gd58(:,:)
      N_Gd60_old(:,:) = N_Gd60(:,:)

      maXS_f_3D_old    (:,:,:) = maXS_f_3D    (:,:,:)

      miXS_a_U34_old   (:,:,:) = miXS_a_U34   (:,:,:)
      miXS_a_U35_old   (:,:,:) = miXS_a_U35   (:,:,:)
      miXS_a_U36_old   (:,:,:) = miXS_a_U36   (:,:,:)
      miXS_a_U37_old   (:,:,:) = miXS_a_U37   (:,:,:)
      miXS_a_U38_old   (:,:,:) = miXS_a_U38   (:,:,:)
      miXS_a_Np37_old  (:,:,:) = miXS_a_Np37  (:,:,:)
      miXS_a_Np38_old  (:,:,:) = miXS_a_Np38  (:,:,:)
      miXS_a_Np39_old  (:,:,:) = miXS_a_Np39  (:,:,:)
      miXS_a_Pu38_old  (:,:,:) = miXS_a_Pu38  (:,:,:)
      miXS_a_Pu39_old  (:,:,:) = miXS_a_Pu39  (:,:,:)
      miXS_a_Pu40_old  (:,:,:) = miXS_a_Pu40  (:,:,:)
      miXS_a_Pu41_old  (:,:,:) = miXS_a_Pu41  (:,:,:)
      miXS_a_Pu42_old  (:,:,:) = miXS_a_Pu42  (:,:,:)
      miXS_a_Pu43_old  (:,:,:) = miXS_a_Pu43  (:,:,:)
      miXS_a_Am41_old  (:,:,:) = miXS_a_Am41  (:,:,:)
      miXS_a_As42_old  (:,:,:) = miXS_a_As42  (:,:,:)
      miXS_a_Am42_old  (:,:,:) = miXS_a_Am42  (:,:,:)
      miXS_a_Am43_old  (:,:,:) = miXS_a_Am43  (:,:,:)
      miXS_a_Am44_old  (:,:,:) = miXS_a_Am44  (:,:,:)
      miXS_a_Cm42_old  (:,:,:) = miXS_a_Cm42  (:,:,:)
      miXS_a_Cm43_old  (:,:,:) = miXS_a_Cm43  (:,:,:)
      miXS_a_Cm44_old  (:,:,:) = miXS_a_Cm44  (:,:,:)

      miXS_f_U34_old   (:,:,:) = miXS_f_U34   (:,:,:)
      miXS_f_U35_old   (:,:,:) = miXS_f_U35   (:,:,:)
      miXS_f_U36_old   (:,:,:) = miXS_f_U36   (:,:,:)
      miXS_f_U37_old   (:,:,:) = miXS_f_U37   (:,:,:)
      miXS_f_U38_old   (:,:,:) = miXS_f_U38   (:,:,:)
      miXS_f_Np37_old  (:,:,:) = miXS_f_Np37  (:,:,:)
      miXS_f_Np38_old  (:,:,:) = miXS_f_Np38  (:,:,:)
      miXS_f_Np39_old  (:,:,:) = miXS_f_Np39  (:,:,:)
      miXS_f_Pu38_old  (:,:,:) = miXS_f_Pu38  (:,:,:)
      miXS_f_Pu39_old  (:,:,:) = miXS_f_Pu39  (:,:,:)
      miXS_f_Pu40_old  (:,:,:) = miXS_f_Pu40  (:,:,:)
      miXS_f_Pu41_old  (:,:,:) = miXS_f_Pu41  (:,:,:)
      miXS_f_Pu42_old  (:,:,:) = miXS_f_Pu42  (:,:,:)
      miXS_f_Pu43_old  (:,:,:) = miXS_f_Pu43  (:,:,:)
      miXS_f_Am41_old  (:,:,:) = miXS_f_Am41  (:,:,:)
      miXS_f_As42_old  (:,:,:) = miXS_f_As42  (:,:,:)
      miXS_f_Am42_old  (:,:,:) = miXS_f_Am42  (:,:,:)
      miXS_f_Am43_old  (:,:,:) = miXS_f_Am43  (:,:,:)
      miXS_f_Am44_old  (:,:,:) = miXS_f_Am44  (:,:,:)
      miXS_f_Cm42_old  (:,:,:) = miXS_f_Cm42  (:,:,:)
      miXS_f_Cm43_old  (:,:,:) = miXS_f_Cm43  (:,:,:)
      miXS_f_Cm44_old  (:,:,:) = miXS_f_Cm44  (:,:,:)

      miXS_a_I35_old   (:,:,:) = miXS_a_I35   (:,:,:)
      miXS_a_Xe35_old  (:,:,:) = miXS_a_Xe35  (:,:,:)
      miXS_a_Nd47_old  (:,:,:) = miXS_a_Nd47  (:,:,:)
      miXS_a_Nd48_old  (:,:,:) = miXS_a_Nd48  (:,:,:)
      miXS_a_Nd49_old  (:,:,:) = miXS_a_Nd49  (:,:,:)
      miXS_a_Pm47_old  (:,:,:) = miXS_a_Pm47  (:,:,:)
      miXS_a_Ps48_old  (:,:,:) = miXS_a_Ps48  (:,:,:)
      miXS_a_Pm48_old  (:,:,:) = miXS_a_Pm48  (:,:,:)
      miXS_a_Pm49_old  (:,:,:) = miXS_a_Pm49  (:,:,:)
      miXS_a_Sm47_old  (:,:,:) = miXS_a_Sm47  (:,:,:)
      miXS_a_Sm48_old  (:,:,:) = miXS_a_Sm48  (:,:,:)
      miXS_a_Sm49_old  (:,:,:) = miXS_a_Sm49  (:,:,:)

      miXS_a_Gd52_old  (:,:,:) = miXS_a_Gd52  (:,:,:)
      miXS_a_Gd54_old  (:,:,:) = miXS_a_Gd54  (:,:,:)
      miXS_a_Gd55_old  (:,:,:) = miXS_a_Gd55  (:,:,:)
      miXS_a_Gd56_old  (:,:,:) = miXS_a_Gd56  (:,:,:)
      miXS_a_Gd57_old  (:,:,:) = miXS_a_Gd57  (:,:,:)
      miXS_a_Gd58_old  (:,:,:) = miXS_a_Gd58  (:,:,:)
      miXS_a_Gd60_old  (:,:,:) = miXS_a_Gd60  (:,:,:)

      miXS_a_H2O_old   (:,:,:) = miXS_a_H2O   (:,:,:)
      miXS_a_B0_old    (:,:,:) = miXS_a_B0    (:,:,:)

      miXS_n2n_U35_old (:,:,:) = miXS_n2n_U35 (:,:,:)
      miXS_n2n_U38_old (:,:,:) = miXS_n2n_U38 (:,:,:)
      miXS_n2n_Np37_old(:,:,:) = miXS_n2n_Np37(:,:,:)
      miXS_n2n_Pu39_old(:,:,:) = miXS_n2n_Pu39(:,:,:)

      miXS_a_B10_old   (:,:,:) = miXS_a_B10   (:,:,:)

      Y_I35_old      (:,:) = Y_I35      (:,:)
      Y_Xe35_old     (:,:) = Y_Xe35     (:,:)
      Y_Nd47_old     (:,:) = Y_Nd47     (:,:)
      Y_Nd48_old     (:,:) = Y_Nd48     (:,:)
      Y_Nd49_old     (:,:) = Y_Nd49     (:,:)
      Y_Pm47_old     (:,:) = Y_Pm47     (:,:)
      Y_Ps48_old     (:,:) = Y_Ps48     (:,:)
      Y_Pm48_old     (:,:) = Y_Pm48     (:,:)
      Y_Pm49_old     (:,:) = Y_Pm49     (:,:)
      Y_Sm49_old     (:,:) = Y_Sm49     (:,:)
      Y_Xe35_eff_old (:,:) = Y_Xe35_eff (:,:)
      Y_Pm49_eff_old (:,:) = Y_Pm49_eff (:,:)

#ifdef tuan_fr
      if (FLag_miXS) then
          N_FR_old (:,:,:) = N_FR(:,:,:)
      endif
      
#endif

#ifdef tuan_fr_tdep
      if (FLag_miXS) then
          N_FR_T_old (:,:,:,:) = N_FR_T(:,:,:,:)
      endif
      
#endif

      RETURN
      END SUBROUTINE Save_Old
!#endif 


      SUBROUTINE Save_History(I_Step)
      USE Inc_3D, ONLY: BU, T_Fuel, T_Mod, D_Mod, Power, Normal_Power, Flux
      USE Inc_CR
      USE Inc_Depletion, ONLY: Cycle_BU, PPMB10, B10PCT, BORINB10, B10DEP, N_BU
      USE Inc_Detector
      USE Inc_Flag
      USE Inc_History
      USE Inc_INP
      USE Inc_3D, ONLY: keff
      USE Inc_Nuclide
      USE Inc_TH
      USE Inc_Time
      USE Inc_XYZ
      USE Mod_GetSome
      USE Inc_PinPOW
      use inc_file, only: flag_ocean
      use inc_file, only: flag_fcn
      use Inc_File, only: flag_out_print
      use inc_transient, only: flag_transient
      use Inc_PinVar, only: npin
      use Inc_maXS
      use Inc_miXS
      use Inc_DF
      use Inc_Branch, only: stack_keff
      use mod_thdriver, only: cal_dnbr
!      use mod_adjoint, only: adjoint_ss, calkp_ss
      use inc_solver, only: Pbeta, Plambda, Pzeta, iPgenT
      use inc_kinetics, only: n_group_d
#ifdef JR_SRCTRM
      use Inc_SNF
      !use Inc_Pinvar, only: phihom
#endif

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I_Step
      INTEGER :: Ixy, Ixy_FA, Iz, Ig, I_CR, Ix_1N, Iy_1N, Iz_Reg
      INTEGER :: Ix, Iy, m, Ixy_1N, Iy_in, Ix_in
      INTEGER :: M1D(1), M2D(2), M3D(3)
      LOGICAL(1) :: Flag_First
      REAL(8) :: tmp1
      real(8), dimension(50) :: Sum_Real
      integer :: FoldCount
      real(8) :: FoldFactor, FoldFactor_1N
      real(8) :: tsum, vsum
#ifdef JR_SRCTRM
      integer(4) :: ip,jp
#endif


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Save_History] in Mod_Save'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      Flag_First = .FALSE.

      Hist_Sec (I_Step) = Sec
      Hist_Hour(I_Step) = Hour
      Hist_Day (I_Step) = Day

      Hist_Cycle_BU  (I_Step) = Cycle_BU(I_Step)
      Hist_Accum_BU  (I_Step) = Cycle_BU(I_Step)
      Hist_PPower    (I_Step) = PPower * DP2
      Hist_PPM       (I_Step) = PPM
      Hist_keff      (I_Step) = keff
      Hist_Reactivity(I_Step) = D1 - ( 1/keff )

      Hist_ADJPPM    (I_Step) = PPM*(B10PCT/BORINB10)
      Hist_PPMB10    (I_Step) = PPMB10
      Hist_B10DEP    (I_Step) = B10DEP

      if (flag_ocean) then
         Hist_Fxy(I_Step) = Fxy
         Hist_Fq (I_Step) = Fq
         Hist_Fr (I_Step) = Fr
      endif
      if (flag_pinpow) then
         Hist_Fq_val (I_Step) = Fq_val
         Hist_Fr_val (I_Step) = Fr_val
         Hist_Fz_val (I_Step) = Fz_val
         Hist_FdH_val(I_Step) = FdH_val
         Hist_max_Fxy_val(I_Step) = max_Fxy_val
         Hist_Fq_loc (I_Step,:) = Fq_loc (:)
         Hist_Fr_loc (I_Step,:) = Fr_loc (:)
         Hist_Fz_loc (I_Step,:) = Fz_loc (:)
         Hist_FdH_loc(I_Step,:) = FdH_loc(:)
         Hist_max_Fxy_loc(I_Step,:) = max_Fxy_loc(:)
         Hist_Fz_ave (I_Step) = Fz_ave
         Hist_Fxy_val(I_Step,:,:)    = Fxy_val_XY(:,:)
         Hist_FdH_val_XY(I_Step,:,:) = FdH_val_XY(:,:)

         tsum=0d0
         do ixy_fa=1,nxy_fa
            ixy=i_fa(ixy_fa)
            do iz=izfuelbot,izfueltop
               tsum=tsum+linpow(ixy,iz)
            enddo
         enddo
         linpow_ave=tsum/nxy_fa/(izfueltop-izfuelbot-1)
         linpow_max=maxval(linpow)
         Hist_linpow_ave(i_step) = linpow_ave
         Hist_linpow_max(i_step) = linpow_max
      endif

      if (flag_pinpow) then
         Hist_maxPBU_2D(I_Step) = max_PinBU_2D
         Hist_maxPPD_2D(I_Step) = max_PinPOW_2D
         Hist_maxPPD_3D(I_Step) = max_PinPOW_3D
         Hist_maxPPW_2D(I_Step) = max_PinQ_2D
         Hist_maxPPW_3D(I_Step) = max_PinQ_3D
         Hist_PBU_2D_loc(I_Step,:) = max_PinBU_2D_loc(:)
         Hist_PPD_2D_loc(I_Step,:) = max_PinPOW_2D_loc(:)
         Hist_PPD_3D_loc(I_Step,:) = max_PinPOW_3D_loc(:)
         Hist_PPW_2D_loc(I_Step,:) = max_PinQ_2D_loc(:)
         Hist_PPW_3D_loc(I_Step,:) = max_PinQ_3D_loc(:)
      endif

#ifdef JR_SRCTRM
      if (.not.(flag_pinpow .and. flag_out_print(19)) .and. flag_snf_pin) then
         call Alloc(Hist_PinPOW_3D, N_BU, Nx_1N, Ny_1N, Nz, npin, npin)
         Hist_PinPow_3D(I_Step,:,:,:,:,:) = PinPOW_3D(:,:,:,:,:)
      endif
#endif


      if (flag_pinpow.and.(flag_fcn.or.flag_mkftn.or.flag_out_print(19))) then
         if (.not. allocated(Hist_PinPow_3D)) then
            if (flag_transient .or. (OPT_Mode==3)) then
               call Alloc(Hist_PinPOW_3D, N_Time, Nx_1N, Ny_1N, Nz, npin, npin)
            else
               call Alloc(Hist_PinPOW_3D, N_BU, Nx_1N, Ny_1N, Nz, npin, npin)
            endif
         endif
         Hist_PinPow_3D(I_Step,:,:,:,:,:) = PinPOW_3D(:,:,:,:,:)
         if (.not. allocated(Hist_PinPow_2D)) then
            if (flag_transient .or. (OPT_Mode==3)) then
               call Alloc(Hist_PinPOW_2D, N_Time, Nx_1N, Ny_1N, npin, npin)
            else
               call Alloc(Hist_PinPOW_2D, N_BU, Nx_1N, Ny_1N, npin, npin)
            endif
         endif
         Hist_PinPow_2D(I_Step,:,:,:,:) = PinPOW_2D(:,:,:,:)
         if (.not. allocated(Hist_PinQ_2D)) then
            if (flag_transient .or. (OPT_Mode==3)) then
               call Alloc(Hist_PinQ_2D, N_Time, Nx_1N, Ny_1N, npin, npin)
            else
               call Alloc(Hist_PinQ_2D, N_BU, Nx_1N, Ny_1N, npin, npin)
            endif
         endif
         Hist_PinQ_2D(I_Step,:,:,:,:) = PinQ_2D(:,:,:,:)
      endif
      if (flag_pinpow .and. flag_out_print(20)) then
         if (.not. allocated(Hist_PinBU_3D)) then
            if (flag_transient .or. (OPT_Mode==3)) then
               call Alloc(Hist_PinBU_3D, N_Time, Nx_1N, Ny_1N, Nz, npin, npin)
            else
               call Alloc(Hist_PinBU_3D, N_BU, Nx_1N, Ny_1N, Nz, npin, npin)
            endif
         endif
         Hist_PinBU_3D(I_Step,:,:,:,:,:) = PinBU_3D(:,:,:,:,:)
         if (.not. allocated(Hist_PinBU_2D)) then
            if (flag_transient .or. (OPT_Mode==3)) then
               call Alloc(Hist_PinBU_2D, N_Time, Nx_1N, Ny_1N, npin, npin)
            else
               call Alloc(Hist_PinBU_2D, N_BU, Nx_1N, Ny_1N, npin, npin)
            endif
         endif
         Hist_PinBU_2D(I_Step,:,:,:,:) = PinBU_2D(:,:,:,:)
      endif

      CALL Get_ASI( ASI, Power )
      Hist_ASI(I_Step) = ASI

      CALL Get_FA_ASI
      Hist_FA_ASI(:,:,i_Step)=FA_ASI(:,:)

      Hist_PF_3D(I_Step) = MAXVAL( Normal_Power_XYZ )
      Hist_PF_2D(I_Step) = MAXVAL( Normal_Power_XY  )
      Hist_PF_1D(I_Step) = MAXVAL( Normal_Power_Z   )

      M3D = MAXLOC(Normal_Power_XYZ)
      M2D = MAXLOC(Normal_Power_XY )
      M1D = MAXLOC(Normal_Power_Z  )

      Hist_PP_3D(I_Step) = 10000*M3D(1) + 100*M3D(2) + M3D(3)
      Hist_PP_2D(I_Step) = 10000*M2D(1) + 100*M2D(2)
      Hist_PP_1D(I_Step) = M1D(1)

      Hist_PBU_FA (I_Step) = MAXVAL( Hist_BU_XY )
      !Hist_PBU_Pin(I_Step)

      M2D = MAXLOC( BU_XY )

      Hist_PPBU_FA = 100000*M2D(1) + 10000*M2D(2)
      !Hist_PPBU_Pin

      IF ( N_CR /= 0 ) THEN
#ifdef tuan_fr_crm
         DO I_CR = 1, N_CR
            Hist_CR(I_Step, I_CR, 1) = Abs_Bot(I_CR)
            Hist_CR(I_Step, I_CR, 2) = D0  ! PDIL
            Hist_CR(I_Step, I_CR, 3) = D0  ! PDWL
         END DO
#else
         DO I_CR = 1, N_CR
            Hist_CR(I_Step, I_CR, 1) = CR_Bot(I_CR)
            Hist_CR(I_Step, I_CR, 2) = D0  ! PDIL
            Hist_CR(I_Step, I_CR, 3) = D0  ! PDWL
         END DO
#endif

      END IF

      !js+CALL Get_Avg( Hist_TF_Avg(I_Step), T_Fuel, 0 )
      tsum=0d0
      vsum=0d0
      do Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)
         do Iz = IzFuelBot, IzFuelTop

            tsum=tsum+t_fuel(ixy,iz)*nodevolume(ixy,iz)
            vsum=vsum+nodevolume(ixy,iz)
         enddo
      enddo
      Hist_TF_Avg(I_Step)=tsum/vsum
      CALL Get_Avg( Hist_TM_Avg(I_Step), T_Mod , 0 )
      CALL Get_Avg( Hist_DM_Avg(I_Step), D_Mod , 0 )

      Hist_TF_Max(I_Step) = MAXVAL(T_Fuel)
      Hist_TF_Min(I_Step) = MAXVAL(T_Fuel)
      DO Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)
         DO Iz = IzFuelBot, IzFuelTop

            IF ( T_Fuel(Ixy, Iz) < Hist_TF_Min(I_Step) ) THEN
               Hist_TF_Min(I_Step) = T_Fuel(Ixy, Iz)
            END IF
         END DO
      END DO

      Hist_TFcen_Max(I_Step) = tfmax + DegToK

      Hist_TM_Max(I_Step) = MAXVAL(T_Mod)
      Hist_TM_Min(I_Step) = MINVAL(T_Mod)
      Hist_TM_In (I_Step) = TM_In + DegToK
      Hist_TM_Out(I_Step) = 0d0
      tmp1 = 0d0
      DO Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)
         tmp1 = tmp1 + NodeVolume(Ixy, Nz)
         Hist_TM_Out(I_Step) = Hist_TM_Out(I_Step) + T_Mod(Ixy, Nz) * NodeVolume(Ixy, Nz)
      END DO
      Hist_TM_Out(I_Step) = Hist_TM_Out(I_Step) / tmp1

      IF ( Hist_TM_Min(I_Step) == D0 ) THEN
         Hist_TM_In(I_Step) = D0
      END IF

      Hist_DM_Max(I_Step) = MAXVAL(D_Mod)
      Hist_DM_Min(I_Step) = MINVAL(D_Mod)

      if (flag_cal_dnbr.and.flag_thfb) then
         call cal_dnbr

         do iz=1,nz
            do ixy=1,nxy
               hist_qchf(i_step,ixy,iz)=qchf(ixy,iz)
               hist_qwat(i_step,ixy,iz)=qwat(ixy,iz)
               hist_dnbr(i_step,ixy,iz)=dnbr(ixy,iz)
            enddo
         enddo
         hist_mqchf(i_step)=mqchf
         hist_mqwat(i_step)=mqwat
         hist_mdnbr(i_step)=mdnbr
      endif

      if (opt_fenthl) then
         Hist_MaxFenthl(I_Step) = max_fenthl
         Hist_MaxFenthl_Rise(I_Step) = max_fenthl_rise
      endif

#ifdef tuan_delete      
     if (flag_out_print(10)) then
      CALL Get_ASI( Hist_Xe35_ASI(I_Step), N_Xe35 )
      CALL Get_ASI( Hist_Sm49_ASI(I_Step), N_Sm49 )
      CALL Get_Avg( Hist_I35 (I_Step), N_I35  , 0 )
      CALL Get_Avg( Hist_Xe35(I_Step), N_Xe35 , 0 )
      CALL Get_Avg( Hist_Nd47(I_Step), N_Nd47 , 0 )
      CALL Get_Avg( Hist_Nd48(I_Step), N_Nd48 , 0 )
      CALL Get_Avg( Hist_Nd49(I_Step), N_Nd49 , 0 )
      CALL Get_Avg( Hist_Pm47(I_Step), N_Pm47 , 0 )
      CALL Get_Avg( Hist_Ps48(I_Step), N_Ps48 , 0 )
      CALL Get_Avg( Hist_Pm48(I_Step), N_Pm48 , 0 )
      CALL Get_Avg( Hist_Pm49(I_Step), N_Pm49 , 0 )
      CALL Get_Avg( Hist_Sm47(I_Step), N_Sm47 , 0 )
      CALL Get_Avg( Hist_Sm48(I_Step), N_Sm48 , 0 )
      CALL Get_Avg( Hist_Sm49(I_Step), N_Sm49 , 0 )
      CALL Get_Avg( Hist_Gd52(I_Step), N_Gd52 , 0 )
      CALL Get_Avg( Hist_Gd54(I_Step), N_Gd54 , 0 )
      CALL Get_Avg( Hist_Gd55(I_Step), N_Gd55 , 0 )
      CALL Get_Avg( Hist_Gd56(I_Step), N_Gd56 , 0 )
      CALL Get_Avg( Hist_Gd57(I_Step), N_Gd57 , 0 )
      CALL Get_Avg( Hist_Gd58(I_Step), N_Gd58 , 0 )
      CALL Get_Avg( Hist_Gd60(I_Step), N_Gd60 , 0 )
     endif

     if (flag_out_print(11)) then
      CALL Get_Avg( Hist_U34 (I_Step), N_U34  , 0 )
      CALL Get_Avg( Hist_U35 (I_Step), N_U35  , 0 )
      CALL Get_Avg( Hist_U36 (I_Step), N_U36  , 0 )
      CALL Get_Avg( Hist_U37 (I_Step), N_U37  , 0 )
      CALL Get_Avg( Hist_U38 (I_Step), N_U38  , 0 )
      CALL Get_Avg( Hist_Np37(I_Step), N_Np37 , 0 )
      CALL Get_Avg( Hist_Np38(I_Step), N_Np38 , 0 )
      CALL Get_Avg( Hist_Np39(I_Step), N_Np39 , 0 )
      CALL Get_Avg( Hist_Pu38(I_Step), N_Pu38 , 0 )
      CALL Get_Avg( Hist_Pu39(I_Step), N_Pu39 , 0 )
      CALL Get_Avg( Hist_Pu40(I_Step), N_Pu40 , 0 )
      CALL Get_Avg( Hist_Pu41(I_Step), N_Pu41 , 0 )
      CALL Get_Avg( Hist_Pu42(I_Step), N_Pu42 , 0 )
      CALL Get_Avg( Hist_Pu43(I_Step), N_Pu43 , 0 )
      CALL Get_Avg( Hist_Am41(I_Step), N_Am41 , 0 )
      CALL Get_Avg( Hist_As42(I_Step), N_As42 , 0 )
      CALL Get_Avg( Hist_Am42(I_Step), N_Am42 , 0 )
      CALL Get_Avg( Hist_Am43(I_Step), N_Am43 , 0 )
      CALL Get_Avg( Hist_Am44(I_Step), N_Am44 , 0 )
      CALL Get_Avg( Hist_Cm42(I_Step), N_Cm42 , 0 )
      CALL Get_Avg( Hist_Cm43(I_Step), N_Cm43 , 0 )
      CALL Get_Avg( Hist_Cm44(I_Step), N_Cm44 , 0 )
     endif
#endif
      DO Iz = 1, Nz
         Hist_Normal_Power_Z(I_Step, Iz) = Normal_Power_Z(Iz)
      END DO

      DO Iy_1N = 1, Ny_1N
         DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
            Hist_Normal_Power_XY(I_Step, Ix_1N, Iy_1N) = Normal_Power_XY(Ix_1N, Iy_1N)
         END DO
      END DO

      DO Iz = 1, Nz
         DO Iy_1N = 1, Ny_1N
            DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
               Hist_Normal_Power_XYZ(I_Step, Ix_1N, Iy_1N, Iz) = Normal_Power_XYZ(Ix_1N, Iy_1N, Iz)
            END DO
         END DO
      END DO

      DO Iz = 1, Nz
         Hist_BU_Z(I_Step, Iz) = BU_Z(Iz)
      END DO

      DO Iy_1N = 1, Ny_1N
         DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
            Hist_BU_XY(I_Step, Ix_1N, Iy_1N) = BU_XY(Ix_1N, Iy_1N)
         END DO
      END DO

      DO Iz = 1, Nz
         DO Iy_1N = 1, Ny_1N
            DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
               Hist_BU_XYZ(I_Step, Ix_1N, Iy_1N, Iz) = BU_XYZ(Ix_1N, Iy_1N, Iz)
            END DO
         END DO
      END DO

      DO Iz = 1, Nz
         Hist_T_Fuel_Z(I_Step, Iz) = T_Fuel_Z(Iz)
         Hist_T_Mod_Z (I_Step, Iz) = T_Mod_Z (Iz)
         Hist_D_Mod_Z (I_Step, Iz) = D_Mod_Z (Iz)
      END DO

      DO Iy_1N = 1, Ny_1N
         DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
            Hist_T_Fuel_XY(I_Step, Ix_1N, Iy_1N) = T_Fuel_XY(Ix_1N, Iy_1N)
            Hist_T_Mod_XY (I_Step, Ix_1N, Iy_1N) = T_Mod_XY (Ix_1N, Iy_1N)
            Hist_D_Mod_XY (I_Step, Ix_1N, Iy_1N) = D_Mod_XY (Ix_1N, Iy_1N)
         END DO
      END DO

      DO Iz = 1, Nz
         DO Iy_1N = 1, Ny_1N
            DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
               Hist_T_Fuel_XYZ(I_Step, Ix_1N, Iy_1N, Iz) = T_Fuel_XYZ(Ix_1N, Iy_1N, Iz)
               Hist_T_Mod_XYZ (I_Step, Ix_1N, Iy_1N, Iz) = T_Mod_XYZ (Ix_1N, Iy_1N, Iz)
               Hist_D_Mod_XYZ (I_Step, Ix_1N, Iy_1N, Iz) = D_Mod_XYZ (Ix_1N, Iy_1N, Iz)
            END DO
         END DO
      END DO

      if ( flag_print_design ) then
         Hist_FTC(I_Step) = (1d0/stack_keff(1) - 1d0/stack_keff( 2))*1d5 / 5d0
         Hist_MTC(I_Step) = (1d0/stack_keff(1) - 1d0/stack_keff( 3))*1d5 / 5d0
         Hist_ITC(I_Step) = (1d0/stack_keff(1) - 1d0/stack_keff( 4))*1d5 / 5d0
         Hist_CRW(I_Step) = (1d0/stack_keff(6) - 1d0/stack_keff( 5))*1d5
         Hist_XEW(I_Step) = (1d0/stack_keff(7) - 1d0/stack_keff( 8))*1d5
         Hist_SMW(I_Step) = (1d0/stack_keff(7) - 1d0/stack_keff( 9))*1d5
         Hist_GDW(I_Step) = (1d0/stack_keff(7) - 1d0/stack_keff(10))*1d5
      endif

#ifdef JR_SRCTRM
     if (flag_SNF) then
      do Ig = 1, N_Group
         do Ixy = 1, Nxy
            do Iz = IzFuelBot, IzFuelTop
               Hist_Flux          (I_Step, Ixy, Iz, Ig) = Flux          (Ixy, Iz, Ig)
            enddo
         enddo
      enddo
      do Ig=1, N_group
         if (flag_snf_pin) then
            DO jp = 1, npin
               DO ip = 1, npin
                  phihom_SNF(I_step,ip,jp,ig)=phihom_M(ip,jp,SNF_Ixy_1N,SNF_k,Ig)
               enddo
            enddo
         endif
      enddo
      if (flag_SNF_pin) then
       DO ip = 1, int(npin)
          DO jp = 1, int(npin)
             SNF_Hist_N_Xe35(I_Step, ip,jp) = SNF_N_Xe35(ip,jp)
             SNF_Hist_N_Sm49(I_Step, ip,jp) = SNF_N_Sm49(ip,jp)
             SNF_Hist_N_I35 (I_Step, ip,jp) = SNF_N_I35 (ip,jp)
             SNF_Hist_N_Nd47(I_Step, ip,jp) = SNF_N_Nd47(ip,jp)
             SNF_Hist_N_Nd48(I_Step, ip,jp) = SNF_N_Nd48(ip,jp)
             SNF_Hist_N_Nd49(I_Step, ip,jp) = SNF_N_Nd49(ip,jp)
             SNF_Hist_N_Pm47(I_Step, ip,jp) = SNF_N_Pm47(ip,jp)
             SNF_Hist_N_Ps48(I_Step, ip,jp) = SNF_N_Ps48(ip,jp)
             SNF_Hist_N_Pm48(I_Step, ip,jp) = SNF_N_Pm48(ip,jp)
             SNF_Hist_N_Pm49(I_Step, ip,jp) = SNF_N_Pm49(ip,jp)
             SNF_Hist_N_Sm47(I_Step, ip,jp) = SNF_N_Sm47(ip,jp)
             SNF_Hist_N_Sm48(I_Step, ip,jp) = SNF_N_Sm48(ip,jp)
             SNF_Hist_N_Gd52(I_Step, ip,jp) = SNF_N_Gd52(ip,jp)
             SNF_Hist_N_Gd54(I_Step, ip,jp) = SNF_N_Gd54(ip,jp)
             SNF_Hist_N_Gd55(I_Step, ip,jp) = SNF_N_Gd55(ip,jp)
             SNF_Hist_N_Gd56(I_Step, ip,jp) = SNF_N_Gd56(ip,jp)
             SNF_Hist_N_Gd57(I_Step, ip,jp) = SNF_N_Gd57(ip,jp)
             SNF_Hist_N_Gd58(I_Step, ip,jp) = SNF_N_Gd58(ip,jp)
             SNF_Hist_N_Gd60(I_Step, ip,jp) = SNF_N_Gd60(ip,jp)
             SNF_Hist_N_U34 (I_Step, ip,jp) = SNF_N_U34 (ip,jp)
             SNF_Hist_N_U35 (I_Step, ip,jp) = SNF_N_U35 (ip,jp)
             SNF_Hist_N_U36 (I_Step, ip,jp) = SNF_N_U36 (ip,jp)
             SNF_Hist_N_U37 (I_Step, ip,jp) = SNF_N_U37 (ip,jp)
             SNF_Hist_N_U38 (I_Step, ip,jp) = SNF_N_U38 (ip,jp)
             SNF_Hist_N_Np37(I_Step, ip,jp) = SNF_N_Np37(ip,jp)
             SNF_Hist_N_Np38(I_Step, ip,jp) = SNF_N_Np38(ip,jp)
             SNF_Hist_N_Np39(I_Step, ip,jp) = SNF_N_Np39(ip,jp)
             SNF_Hist_N_Pu38(I_Step, ip,jp) = SNF_N_Pu38(ip,jp)
             SNF_Hist_N_Pu39(I_Step, ip,jp) = SNF_N_Pu39(ip,jp)
             SNF_Hist_N_Pu40(I_Step, ip,jp) = SNF_N_Pu40(ip,jp)
             SNF_Hist_N_Pu41(I_Step, ip,jp) = SNF_N_Pu41(ip,jp)
             SNF_Hist_N_Pu42(I_Step, ip,jp) = SNF_N_Pu42(ip,jp)
             SNF_Hist_N_Pu43(I_Step, ip,jp) = SNF_N_Pu43(ip,jp)
             SNF_Hist_N_Am41(I_Step, ip,jp) = SNF_N_Am41(ip,jp)
             SNF_Hist_N_As42(I_Step, ip,jp) = SNF_N_As42(ip,jp)
             SNF_Hist_N_Am42(I_Step, ip,jp) = SNF_N_Am42(ip,jp)
             SNF_Hist_N_Am43(I_Step, ip,jp) = SNF_N_Am43(ip,jp)
             SNF_Hist_N_Am44(I_Step, ip,jp) = SNF_N_Am44(ip,jp)
             SNF_Hist_N_Cm42(I_Step, ip,jp) = SNF_N_Cm42(ip,jp)
             SNF_Hist_N_Cm43(I_Step, ip,jp) = SNF_N_Cm43(ip,jp)
             SNF_Hist_N_Cm44(I_Step, ip,jp) = SNF_N_Cm44(ip,jp)
          END DO
       END DO
      endif
     endif
#endif

     if (flag_out_print(23)) then
      do Ig = 1, N_Group
         do Ixy = 1, Nxy
            do Iz = IzFuelBot, IzFuelTop
               Hist_Flux          (I_Step, Ixy, Iz, Ig) = Flux          (Ixy, Iz, Ig)
               Hist_maXS_tr_3D    (I_Step, Ixy, Iz, Ig) = maXS_tr_3D    (Ixy, Iz, Ig)
               Hist_D_3D          (I_Step, Ixy, Iz, Ig) = D_3D          (Ixy, Iz, Ig)
               Hist_maXS_a_3D     (I_Step, Ixy, Iz, Ig) = maXS_a_3D     (Ixy, Iz, Ig)
               Hist_maXS_f_3D     (I_Step, Ixy, Iz, Ig) = maXS_f_3D     (Ixy, Iz, Ig)
               Hist_nu_maXS_f_3D  (I_Step, Ixy, Iz, Ig) = nu_maXS_f_3D  (Ixy, Iz, Ig)
               Hist_kap_maXS_f_3D (I_Step, Ixy, Iz, Ig) = kap_maXS_f_3D (Ixy, Iz, Ig)
               Hist_maXS_s_3D     (I_Step, Ixy, Iz, Ig) = maXS_s_3D     (Ixy, Iz, Ig)
               Hist_maXS_r_3D     (I_Step, Ixy, Iz, Ig) = maXS_r_3D     (Ixy, Iz, Ig)
               Hist_ADF_Lx        (I_Step, Ixy, Iz, Ig) = ADF_Lx        (Ixy, Iz, Ig)
               Hist_ADF_Rx        (I_Step, Ixy, Iz, Ig) = ADF_Rx        (Ixy, Iz, Ig)
               Hist_ADF_Ly        (I_Step, Ixy, Iz, Ig) = ADF_Ly        (Ixy, Iz, Ig)
               Hist_ADF_Ry        (I_Step, Ixy, Iz, Ig) = ADF_Ry        (Ixy, Iz, Ig)
               Hist_CDF_LxLy      (I_Step, Ixy, Iz, Ig) = CDF_LxLy      (Ixy, Iz, Ig)
               Hist_CDF_LxRy      (I_Step, Ixy, Iz, Ig) = CDF_LxRy      (Ixy, Iz, Ig)
               Hist_CDF_RxLy      (I_Step, Ixy, Iz, Ig) = CDF_RxLy      (Ixy, Iz, Ig)
               Hist_CDF_RxRy      (I_Step, Ixy, Iz, Ig) = CDF_RxRy      (Ixy, Iz, Ig)
            enddo
         enddo
      enddo
     endif

      DO Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)
         DO Iz = IzFuelBot, IzFuelTop
            Hist_Power (I_Step, Ixy, Iz) = Power (Ixy, Iz)
            Hist_N_Xe35(I_Step, Ixy, Iz) = N_Xe35(Ixy, Iz)
            Hist_N_Sm49(I_Step, Ixy, Iz) = N_Sm49(Ixy, Iz)
         END DO
      END DO
     if (flag_out_print(22)) then
      DO Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)
         DO Iz = IzFuelBot, IzFuelTop
            Hist_Burnup(I_Step, Ixy, Iz) = BU    (Ixy, Iz)
            Hist_T_Fuel(I_Step, Ixy, Iz) = T_Fuel(Ixy, Iz)
            Hist_T_Mod (I_Step, Ixy, Iz) = T_Mod (Ixy, Iz)
            Hist_D_Mod (I_Step, Ixy, Iz) = D_Mod (Ixy, Iz)
            Hist_N_I35 (I_Step, Ixy, Iz) = N_I35 (Ixy, Iz)
            Hist_N_Nd47(I_Step, Ixy, Iz) = N_Nd47(Ixy, Iz)
            Hist_N_Nd48(I_Step, Ixy, Iz) = N_Nd48(Ixy, Iz)
            Hist_N_Nd49(I_Step, Ixy, Iz) = N_Nd49(Ixy, Iz)
            Hist_N_Pm47(I_Step, Ixy, Iz) = N_Pm47(Ixy, Iz)
            Hist_N_Ps48(I_Step, Ixy, Iz) = N_Ps48(Ixy, Iz)
            Hist_N_Pm48(I_Step, Ixy, Iz) = N_Pm48(Ixy, Iz)
            Hist_N_Pm49(I_Step, Ixy, Iz) = N_Pm49(Ixy, Iz)
            Hist_N_Sm47(I_Step, Ixy, Iz) = N_Sm47(Ixy, Iz)
            Hist_N_Sm48(I_Step, Ixy, Iz) = N_Sm48(Ixy, Iz)
            Hist_N_Gd52(I_Step, Ixy, Iz) = N_Gd52(Ixy, Iz)
            Hist_N_Gd54(I_Step, Ixy, Iz) = N_Gd54(Ixy, Iz)
            Hist_N_Gd55(I_Step, Ixy, Iz) = N_Gd55(Ixy, Iz)
            Hist_N_Gd56(I_Step, Ixy, Iz) = N_Gd56(Ixy, Iz)
            Hist_N_Gd57(I_Step, Ixy, Iz) = N_Gd57(Ixy, Iz)
            Hist_N_Gd58(I_Step, Ixy, Iz) = N_Gd58(Ixy, Iz)
            Hist_N_Gd60(I_Step, Ixy, Iz) = N_Gd60(Ixy, Iz)
            Hist_N_U34 (I_Step, Ixy, Iz) = N_U34 (Ixy, Iz)
            Hist_N_U35 (I_Step, Ixy, Iz) = N_U35 (Ixy, Iz)
            Hist_N_U36 (I_Step, Ixy, Iz) = N_U36 (Ixy, Iz)
            Hist_N_U37 (I_Step, Ixy, Iz) = N_U37 (Ixy, Iz)
            Hist_N_U38 (I_Step, Ixy, Iz) = N_U38 (Ixy, Iz)
            Hist_N_Np37(I_Step, Ixy, Iz) = N_Np37(Ixy, Iz)
            Hist_N_Np38(I_Step, Ixy, Iz) = N_Np38(Ixy, Iz)
            Hist_N_Np39(I_Step, Ixy, Iz) = N_Np39(Ixy, Iz)
            Hist_N_Pu38(I_Step, Ixy, Iz) = N_Pu38(Ixy, Iz)
            Hist_N_Pu39(I_Step, Ixy, Iz) = N_Pu39(Ixy, Iz)
            Hist_N_Pu40(I_Step, Ixy, Iz) = N_Pu40(Ixy, Iz)
            Hist_N_Pu41(I_Step, Ixy, Iz) = N_Pu41(Ixy, Iz)
            Hist_N_Pu42(I_Step, Ixy, Iz) = N_Pu42(Ixy, Iz)
            Hist_N_Pu43(I_Step, Ixy, Iz) = N_Pu43(Ixy, Iz)
            Hist_N_Am41(I_Step, Ixy, Iz) = N_Am41(Ixy, Iz)
            Hist_N_As42(I_Step, Ixy, Iz) = N_As42(Ixy, Iz)
            Hist_N_Am42(I_Step, Ixy, Iz) = N_Am42(Ixy, Iz)
            Hist_N_Am43(I_Step, Ixy, Iz) = N_Am43(Ixy, Iz)
            Hist_N_Am44(I_Step, Ixy, Iz) = N_Am44(Ixy, Iz)
            Hist_N_Cm42(I_Step, Ixy, Iz) = N_Cm42(Ixy, Iz)
            Hist_N_Cm43(I_Step, Ixy, Iz) = N_Cm43(Ixy, Iz)
            Hist_N_Cm44(I_Step, Ixy, Iz) = N_Cm44(Ixy, Iz)
         END DO
      END DO
     endif
     if (flag_out_print(21)) then
      FoldFactor    = D0
      FoldFactor_1N = D0
      FoldCount     = 0
      DO Iy_1N = 1, Ny_1N
         DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
            Ixy_1N = IxIy_1NToIxy_1N(Ix_1N, Iy_1N)
            DO Iz = 1, Nz
               Sum_Real   = D0
               FoldCount  = 0
               DO m = 1, 4
                  IF (I_1Nto4N(Ixy_1N, m) == 0) THEN
                     FoldCount = FoldCount + 1
                     CYCLE
                  ELSE
                     Sum_Real(01) = Sum_Real(01) + Power (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(02) = Sum_Real(02) + N_Xe35(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(03) = Sum_Real(03) + N_Sm49(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(04) = Sum_Real(04) + BU    (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(05) = Sum_Real(05) + T_Fuel(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(06) = Sum_Real(06) + T_Mod (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(07) = Sum_Real(07) + D_Mod (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(08) = Sum_Real(08) + N_I35 (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(09) = Sum_Real(09) + N_Nd47(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(10) = Sum_Real(10) + N_Nd48(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(11) = Sum_Real(11) + N_Nd49(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(12) = Sum_Real(12) + N_Pm47(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(13) = Sum_Real(13) + N_Ps48(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(14) = Sum_Real(14) + N_Pm48(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(15) = Sum_Real(15) + N_Pm49(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(16) = Sum_Real(16) + N_Sm47(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(17) = Sum_Real(17) + N_Sm48(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(18) = Sum_Real(18) + N_Gd52(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(19) = Sum_Real(19) + N_Gd54(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(20) = Sum_Real(20) + N_Gd55(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(21) = Sum_Real(21) + N_Gd56(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(22) = Sum_Real(22) + N_Gd57(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(23) = Sum_Real(23) + N_Gd58(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(24) = Sum_Real(24) + N_Gd60(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(25) = Sum_Real(25) + N_U34 (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(26) = Sum_Real(26) + N_U35 (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(27) = Sum_Real(27) + N_U36 (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(28) = Sum_Real(28) + N_U37 (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(29) = Sum_Real(29) + N_U38 (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(30) = Sum_Real(30) + N_Np37(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(31) = Sum_Real(31) + N_Np38(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(32) = Sum_Real(32) + N_Np39(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(33) = Sum_Real(33) + N_Pu38(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(34) = Sum_Real(34) + N_Pu39(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(35) = Sum_Real(35) + N_Pu40(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(36) = Sum_Real(36) + N_Pu41(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(37) = Sum_Real(37) + N_Pu42(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(38) = Sum_Real(38) + N_Pu43(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(39) = Sum_Real(39) + N_Am41(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(40) = Sum_Real(40) + N_As42(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(41) = Sum_Real(41) + N_Am42(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(42) = Sum_Real(42) + N_Am43(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(43) = Sum_Real(43) + N_Am44(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(44) = Sum_Real(44) + N_Cm42(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(45) = Sum_Real(45) + N_Cm43(I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(46) = Sum_Real(46) + N_Cm44(I_1Nto4N(Ixy_1N, m), Iz)
                  END IF
               END DO
               FoldFactor = D4/(D4 - REAL(FoldCount, 8))
               Hist_Power_FA (I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(01)/D4)
               Hist_N_Xe35_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(02)/D4)
               Hist_N_Sm49_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(03)/D4)
               Hist_Burnup_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(04)/D4)
               Hist_T_Fuel_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(05)/D4)
               Hist_T_Mod_FA (I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(06)/D4)
               Hist_D_Mod_FA (I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(07)/D4)
               Hist_N_I35_FA (I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(08)/D4)
               Hist_N_Nd47_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(09)/D4)
               Hist_N_Nd48_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(10)/D4)
               Hist_N_Nd49_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(11)/D4)
               Hist_N_Pm47_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(12)/D4)
               Hist_N_Ps48_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(13)/D4)
               Hist_N_Pm48_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(14)/D4)
               Hist_N_Pm49_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(15)/D4)
               Hist_N_Sm47_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(16)/D4)
               Hist_N_Sm48_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(17)/D4)
               Hist_N_Gd52_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(18)/D4)
               Hist_N_Gd54_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(19)/D4)
               Hist_N_Gd55_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(20)/D4)
               Hist_N_Gd56_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(21)/D4)
               Hist_N_Gd57_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(22)/D4)
               Hist_N_Gd58_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(23)/D4)
               Hist_N_Gd60_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(24)/D4)
               Hist_N_U34_FA (I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(25)/D4)
               Hist_N_U35_FA (I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(26)/D4)
               Hist_N_U36_FA (I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(27)/D4)
               Hist_N_U37_FA (I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(28)/D4)
               Hist_N_U38_FA (I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(29)/D4)
               Hist_N_Np37_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(30)/D4)
               Hist_N_Np38_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(31)/D4)
               Hist_N_Np39_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(32)/D4)
               Hist_N_Pu38_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(33)/D4)
               Hist_N_Pu39_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(34)/D4)
               Hist_N_Pu40_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(35)/D4)
               Hist_N_Pu41_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(36)/D4)
               Hist_N_Pu42_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(37)/D4)
               Hist_N_Pu43_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(38)/D4)
               Hist_N_Am41_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(39)/D4)
               Hist_N_As42_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(40)/D4)
               Hist_N_Am42_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(41)/D4)
               Hist_N_Am43_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(42)/D4)
               Hist_N_Am44_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(43)/D4)
               Hist_N_Cm42_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(44)/D4)
               Hist_N_Cm43_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(45)/D4)
               Hist_N_Cm44_FA(I_Step, Ixy_1N, Iz) = FoldFactor*(Sum_Real(46)/D4)
            END DO
         END DO
      END DO
     endif

      IF ( Flag_ExDET ) THEN
         DO Iz_Reg = 1, 3
            Hist_DR(I_Step, Iz_Reg) = dr_core(Iz_Reg)
         END DO

         IF ( dr_core(1) + dr_core(3) < 1e-30 ) THEN
            Hist_ASI_DR(I_Step) = D0
         ELSE
            Hist_ASI_DR(I_Step) = (dr_core(1) - dr_core(3)) / (dr_core(1) + dr_core(3))
         END IF
      END IF

      if (flag_indet .and. flag_det_sig) then
         Hist_det_sig_rr(I_Step,:,:)=det_sig_rr
         if (opt_prt_detpow) then
            Hist_det_detpow(I_Step,:,:)=det_detpow
         endif
         if (opt_prt_wprime) then
            Hist_det_wprime(I_Step,:,:)=det_wprime
         endif
      endif

      if (flag_ocean) then
         DO Iz = 1, Nz
            DO Iy_1N = 1, Ny_1N
               DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
                  Hist_N_Xe35_XYZ (I_Step, Ix_1N, Iy_1N, Iz) = N_Xe35_XYZ (Ix_1N, Iy_1N, Iz)
                  Hist_N_Sm49_XYZ (I_Step, Ix_1N, Iy_1N, Iz) = N_Sm49_XYZ (Ix_1N, Iy_1N, Iz)
                  Hist_FFlux_XYZ (I_Step, Ix_1N, Iy_1N, Iz) = FFlux_XYZ (Ix_1N, Iy_1N, Iz)
                  Hist_TFlux_XYZ (I_Step, Ix_1N, Iy_1N, Iz) = TFlux_XYZ (Ix_1N, Iy_1N, Iz)
               END DO
            END DO
         END DO

         DO Iz = 1, Nz
            Iy_in = 1
            DO Iy_1N = 1, Ny_1N
               Ix_in = 1
               DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
                  Ixy_1N = IxIy_1NToIxy_1N(Ix_1N, Iy_1N)

                  DO m = 1, 2
                     IF (I_1Nto4N(Ixy_1N, m) == 0) THEN
                        IF (Ix_1N /= 1) THEN
                           Ix_in = Ix_in + 1
                        END IF
                        CYCLE
                     ELSE
                        Hist_Normal_Power_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = Normal_Power (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_BU_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = BU (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_T_Fuel_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = T_Fuel (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_T_Mod_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = T_Mod (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_N_Xe35_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = N_Xe35 (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_N_Sm49_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = N_Sm49 (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_FFlux_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = Flux (I_1Nto4N(Ixy_1N, m), Iz, 1)
                        Hist_TFlux_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = Flux (I_1Nto4N(Ixy_1N, m), Iz, 2)
                        Ix_in = Ix_in + 1
                     END IF
                  END DO
               END DO

               IF (OPT_Core == 4 .and. Iy_1N == 1) THEN
                  CONTINUE
               ELSE
                  Iy_in = Iy_in + 1
               END IF

               Ix_in = 1

               DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
                  Ixy_1N = IxIy_1NToIxy_1N(Ix_1N, Iy_1N)

                  DO m = 3, 4
                     IF (I_1Nto4N(Ixy_1N, m) == 0) THEN
                        IF (Ix_1N /= 1) THEN
                           Ix_in = Ix_in + 1
                        END IF

                        CYCLE
                     ELSE
                        Hist_Normal_Power_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = Normal_Power (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_BU_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = BU (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_T_Fuel_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = T_FUEL (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_T_Mod_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = T_Mod (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_N_Xe35_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = N_Xe35 (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_N_Sm49_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = N_Sm49 (I_1Nto4N(Ixy_1N, m), Iz)
                        Hist_FFlux_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = Flux (I_1Nto4N(Ixy_1N, m), Iz, 1)
                        Hist_TFlux_XYZ_4N(I_Step, Ix_in, Iy_in, Iz) = Flux (I_1Nto4N(Ixy_1N, m), Iz, 2)
                        Ix_in = Ix_in + 1
                     END IF
                  END DO
               END DO

               Iy_in = Iy_in + 1

            END DO
         END DO
      endif

      DO Iy = 1, Ny
         DO Ix = Ix_Start_y(Iy), Ix_End_y(Iy)
            Hist_Normal_Power_XY_4N(I_Step, Ix, Iy) = Normal_Power_XY_4N(Ix, Iy)
         END DO
      END DO
      DO Iy = 1, Ny
         DO Ix = Ix_Start_y(Iy), Ix_End_y(Iy)
            Hist_BU_XY_4N(I_Step, Ix, Iy) = BU_XY_4N(Ix, Iy)
         END DO
      END DO

      call cal_core_param(i_step)
      call cal_FA_kinf(i_step)

      if (flag_print_kp) then
!         call adjoint_ss
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif
!         call calkp_ss
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif
         do ig=1,n_group_d
            hist_kp_beta   (i_step,ig)=Pbeta(ig)
            hist_kp_lambda (i_step,ig)=Plambda(ig,1)
            hist_kp_zeta   (i_step,ig)=Pzeta(ig,1)
            hist_kp_gentime(i_step)   =iPgenT(1)
         enddo
      endif

      RETURN
      END SUBROUTINE Save_History


!!!#ifdef siarhei_delete 
      SUBROUTINE XYZ_Indexing
      USE Inc_3D
      USE Inc_File
      USE Inc_INP
      USE Inc_Option
      USE Inc_XYZ
      USE Inc_nuclide
      USE Inc_History, only: HIST_BU_XYZ_4N_ANC

      use Inc_TH
      use Inc_Flag, only: flag_tffb,flag_tmfb,flag_THFB ! last - @$^ Siarhei_FR

#ifdef js_r2mpi
      use inc_parallel, only: comm
#endif
      IMPLICIT NONE
      REAL(8) :: FoldFactor, FoldFactor_1N
      REAL(8), DIMENSION(9) :: Sum_Real
      INTEGER :: Ixy_1N, Iz, Ix_1N, Iy_1N
      INTEGER :: Ix, Iy, Ixy
      INTEGER :: FoldCount, m
#ifdef siarhei_fr
      integer ::  ixy_fa
#endif

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [XYZ_Indexing] in Mod_Save'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif



#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      FoldFactor    = D0
      FoldFactor_1N = D0
      FoldCount     = 0
#ifdef siarhei_fr
      ! @$^ Siarhei_FR - used single-core here
      if (flag_THFB) then ! @$^ Siarhei_FR - to avoid errors when TH is off
         do ixy_fa=1,nxy_fa
            ixy=i_fa(ixy_fa)
            if (flag_tffb) then
               do iz=izfuelbot,izfueltop
                  t_fuel(ixy,iz)=tdopl(iz,ixy_fa)*tdopl(iz,ixy_fa)
               end do
            end if
            if (flag_tmfb) then
                  do iz=1,nz
                     t_mod(ixy,iz)=tcool(iz,ixy_fa)+degtok
                     d_mod(ixy,iz)=dcool(iz,ixy_fa)*1d-3
                  end do
            end if
         end do
      end if
#endif


      DO Iy_1N = 1, Ny_1N
         DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
            Ixy_1N = IxIy_1NToIxy_1N(Ix_1N, Iy_1N)
            DO Iz = 1, Nz
               Sum_Real   = D0
               FoldCount  = 0
               DO m = 1, 4
                  IF (I_1Nto4N(Ixy_1N, m) == 0) THEN
                     FoldCount = FoldCount + 1
                     CYCLE
                  ELSE
                     Sum_Real(1) = Sum_Real(1) + Normal_Power (I_1Nto4N(Ixy_1N, m), Iz)

                     Sum_Real(2) = Sum_Real(2) + BU           (I_1Nto4N(Ixy_1N, m), Iz)

                     Sum_Real(3) = Sum_Real(3) + T_Fuel       (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(4) = Sum_Real(4) + T_Mod        (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(5) = Sum_Real(5) + D_Mod        (I_1Nto4N(Ixy_1N, m), Iz)
                     if (flag_ocean) then
                        Sum_Real(6) = Sum_Real(6) + N_Xe35 (I_1Nto4N(Ixy_1N, m), Iz)
                        Sum_Real(7) = Sum_Real(7) + N_Sm49 (I_1Nto4N(Ixy_1N, m), Iz)
                        Sum_Real(8) = Sum_Real(8) + Flux   (I_1Nto4N(Ixy_1N, m), Iz, 1)
                        Sum_Real(9) = Sum_Real(9) + Flux   (I_1Nto4N(Ixy_1N, m), Iz, 2)
                     endif
                  END IF
               END DO

               FoldFactor = D4/(D4 - REAL(FoldCount, 8))
               Normal_Power_XYZ(Ix_1N, Iy_1N, Iz) = FoldFactor*(Sum_Real(1)/D4)
               BU_XYZ(Ix_1N, Iy_1N, Iz) = FoldFactor*(Sum_Real(2)/D4)
               T_Fuel_XYZ(Ix_1N, Iy_1N, Iz) = FoldFactor*(Sum_Real(3)/D4)
               T_Mod_XYZ(Ix_1N, Iy_1N, Iz) = FoldFactor*(Sum_Real(4)/D4)
               D_Mod_XYZ(Ix_1N, Iy_1N, Iz) = FoldFactor*(Sum_Real(5)/D4)
               if (flag_ocean) then
                  N_Xe35_XYZ(Ix_1N, Iy_1N, Iz) = FoldFactor*(Sum_Real(6)/D4)
                  N_Sm49_XYZ(Ix_1N, Iy_1N, Iz) = FoldFactor*(Sum_Real(7)/D4)
                  FFlux_XYZ(Ix_1N, Iy_1N, Iz) = FoldFactor*(Sum_Real(8)/D4)
                  TFlux_XYZ(Ix_1N, Iy_1N, Iz) = FoldFactor*(Sum_Real(9)/D4)
               endif
            END DO
         END DO
      END DO

      DO Iy = 1, Ny
         Do Ix = 1, Nx
            if (IxIytoIxy(Ix,Iy)==0) cycle
            Sum_Real = D0
            DO Iz = IzFuelBot, IzFuelTop
               HIST_BU_XYZ_4N_ANC(Ix,Iy,Iz) = BU (IxIytoIxy(Ix,Iy), Iz)
               Sum_Real(1) = Sum_Real(1) + HIST_BU_XYZ_4N_ANC(Ix,Iy,Iz)* h_z(Iz)
            END DO
            BU_XY_4N(Ix,Iy) = Sum_Real(1) / SUM( h_z(IzFuelBot:IzFuelTop) )
         END DO
      END DO

      DO Iy = 1, Ny
         DO Ix = Ix_Start_y(Iy), Ix_End_y(Iy)
            Ixy = IxIyToIxy(Ix, Iy)
            Sum_Real = D0
            DO Iz = IzFuelBot, IzFuelTop
               Sum_Real(1) = Sum_Real(1) + Normal_Power(Ixy, Iz) * h_z(Iz)
            END DO
            Normal_Power_XY_4N(Ix, Iy) = Sum_Real(1) / SUM( h_z(IzFuelBot:IzFuelTop) )
         END DO
      END DO
      DO Iy_1N = 1, Ny_1N
         DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
            Ixy_1N = IxIy_1NToIxy_1N(Ix_1N, Iy_1N)
            Sum_Real = D0
            DO Iz = IzFuelBot, IzFuelTop
               Sum_Real(1) = Sum_Real(1) + Normal_Power_XYZ (Ix_1N, Iy_1N, Iz) * h_z(Iz)
               Sum_Real(2) = Sum_Real(2) + BU_XYZ           (Ix_1N, Iy_1N, Iz) * h_z(Iz)
               Sum_Real(3) = Sum_Real(3) + T_Fuel_XYZ       (Ix_1N, Iy_1N, Iz) * h_z(Iz)
            END DO
            DO Iz = 1, Nz
               Sum_Real(4) = Sum_Real(4) + T_Mod_XYZ (Ix_1N, Iy_1N, Iz) * h_z(Iz)
               Sum_Real(5) = Sum_Real(5) + D_Mod_XYZ (Ix_1N, Iy_1N, Iz) * h_z(Iz)
            END DO

            Normal_Power_XY(Ix_1N, Iy_1N) = Sum_Real(1) / SUM( h_z(IzFuelBot:IzFuelTop) )
            BU_XY(Ix_1N, Iy_1N) = Sum_Real(2) / SUM( h_z(IzFuelBot:IzFuelTop) )
            T_Fuel_XY(Ix_1N, Iy_1N) = Sum_Real(3) / SUM( h_z(IzFuelBot:IzFuelTop) )
            T_Mod_XY(Ix_1N, Iy_1N) = Sum_Real(4) / SUM( h_z(1:Nz) )
            D_Mod_XY(Ix_1N, Iy_1N) = Sum_Real(5) / SUM( h_z(1:Nz) )
         END DO
      END DO

      DO Iz = IzFuelBot, IzFuelTop
         Sum_Real = D0
         FoldCount = 0
         DO Iy_1N = 1, Ny_1N
            DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
               Ixy_1N = IxIy_1NToIxy_1N(Ix_1N, Iy_1N)
               IF ( I_FARF_1N(Ixy_1N) == 1 ) CYCLE
               DO m = 1, 4
                  IF (I_1Nto4N(Ixy_1N, m) /= 0) THEN

                     Sum_Real(1) = Sum_Real(1) + Normal_Power (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(2) = Sum_Real(2) + BU           (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(3) = Sum_Real(3) + T_Fuel       (I_1Nto4N(Ixy_1N, m), Iz)
                     FoldCount = FoldCount + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         Normal_Power_Z(Iz) = Sum_Real(1)/REAL(FoldCount,8)
         BU_Z(Iz)           = Sum_Real(2)/REAL(FoldCount,8)
         T_Fuel_Z(Iz)       = Sum_Real(3)/REAL(FoldCount,8)
      ENDDO

      DO Iz = 1, Nz
         Sum_Real = D0
         FoldCount = 0
         DO Iy_1N = 1, Ny_1N
            DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
               Ixy_1N = IxIy_1NToIxy_1N(Ix_1N, Iy_1N)
               IF ( I_FARF_1N(Ixy_1N) == 1 ) CYCLE
               DO m = 1, 4
                  IF (I_1Nto4N(Ixy_1N, m) /= 0) THEN
                     Sum_Real(4) = Sum_Real(4) + T_Mod      (I_1Nto4N(Ixy_1N, m), Iz) &
                                               * NodeVolume (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(5) = Sum_Real(5) + D_Mod      (I_1Nto4N(Ixy_1N, m), Iz) &
                                               * NodeVolume (I_1Nto4N(Ixy_1N, m), Iz)
                     Sum_Real(6) = Sum_Real(6) + NodeVolume (I_1Nto4N(Ixy_1N, m), Iz)
                     FoldCount = FoldCount + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         T_Mod_Z(Iz) = Sum_Real(4)/Sum_Real(6)
         D_Mod_Z(Iz) = Sum_Real(5)/Sum_Real(6)
      ENDDO

      RETURN
      END SUBROUTINE XYZ_Indexing
!!!#endif 


      subroutine cal_core_param(i_step)
      use inc_3d, only: keff
      use inc_3d, only: flux
      use inc_maxs, only: nu_maxs_f_3d, maxs_s_3d, maxs_r_3d, maxs_a_3d, maxs_f_3d, d_3d
      use inc_kinetics, only: v_inv
      use inc_history, only: hist_6factor_eps, hist_6factor_p, hist_6factor_eta
      use inc_history, only: hist_6factor_f, hist_6factor_fnl, hist_6factor_tnl, hist_6factor_kinf
      use inc_history, only: hist_core_avg_nu
      use inc_history, only: hist_mfp, hist_mfp_f, hist_mfp_t
      use inc_history, only: hist_n_density, hist_n_density_f, hist_n_density_t
      use inc_history, only: hist_n_speed, hist_n_speed_f, hist_n_speed_t
      use inc_history, only: hist_spectral_idx, hist_n_flux, hist_n_flux_f, hist_n_flux_t
      use inc_history, only: hist_dif_length, hist_dif_length_f, hist_dif_length_t
      use inc_history, only: hist_mig_length
      use inc_history, only: hist_n_lifetime, hist_n_lifetime_f, hist_n_lifetime_t
      use inc_history, only: hist_n_gentime, hist_n_gentime_f, hist_n_gentime_t
      use inc_history, only: hist_n_energy, hist_n_energy_f, hist_n_energy_t
      implicit none
      integer, intent(in) :: i_step
      integer :: ig, iz, ixy, ixy_fa
      real(8) :: rr_v                          ! vol
      real(8) :: rr_fv, rr_fv_1,  rr_fv_2      ! flux * vol
      real(8) :: rr_nffv, rr_nffv_1, rr_nffv_2 ! nu_fission * flux * vol
      real(8) :: rr_sfv_1                      ! s12 * flux * vol
      real(8) :: rr_rfv_1                      ! r * flux * vol
      real(8) :: rr_afv, rr_afv_2              ! a * flux * vol
      real(8) :: rr_ffv                        ! fission * flux * vol
      real(8) :: rr_vfv, rr_vfv_1, rr_vfv_2    ! v_inv * flux * vol
      real(8) :: rr_dfv, rr_dfv_1, rr_dfv_2    ! D * flux * vol
      real(8) :: vtmp, vtmp1, vtmp2
      real(8) :: selfscat1, selfscat2          ! arbitrary self scattering XS

      ! six factor formula : epsilon (fast fission factor)

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [cal_core_param] in Mod_Save'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      rr_nffv=0d0
      rr_nffv_2=0d0
      do iz=1,nz
         do ixy=1,nxy
            do ig=1,n_group
               rr_nffv=rr_nffv+nu_maxs_f_3d(ixy,iz,ig)*flux(ixy,iz,ig)*nodevolume(ixy,iz)
            enddo
            rr_nffv_2=rr_nffv_2+nu_maxs_f_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
         enddo
      enddo
      hist_6factor_eps(i_step)=(rr_nffv)/max(1d-10,rr_nffv_2)

      ! six factor formula : p (resonance escape probability)
      rr_sfv_1=0d0
      rr_rfv_1=0d0
      do iz=1,nz
         do ixy=1,nxy
            rr_sfv_1=rr_sfv_1+maxs_s_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_rfv_1=rr_rfv_1+maxs_r_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
         enddo
      enddo
      hist_6factor_p(i_step)=(rr_sfv_1)/max(1d-10,rr_rfv_1)

      ! six factor formula : eta (reproduction factor)
      rr_nffv_2=0d0
      rr_afv_2=0d0
      do iz=izfuelbot,izfueltop
         do ixy_fa=1,nxy_fa
            ixy=i_fa(ixy_fa)
            rr_nffv_2=rr_nffv_2+nu_maxs_f_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_afv_2=rr_afv_2+maxs_a_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
         enddo
      enddo
      hist_6factor_eta(i_step)=(rr_nffv_2)/max(1d-10,rr_afv_2)

      ! six factor formula : f (thermal utilization factor)
      rr_afv=0d0
      do iz=1,nz
         do ixy=1,nxy
            rr_afv=rr_afv+maxs_a_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
         enddo
      enddo
      hist_6factor_f(i_step)=(rr_afv_2)/max(1d-10,rr_afv)

      ! six factor formula : P_fnl (fast neutron non leakage)
      rr_rfv_1=0d0
      rr_nffv=0d0
      do iz=1,nz
         do ixy=1,nxy
            rr_rfv_1=rr_rfv_1+maxs_r_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            do ig=1,n_group
               rr_nffv=rr_nffv+nu_maxs_f_3d(ixy,iz,ig)*flux(ixy,iz,ig)*nodevolume(ixy,iz)
            enddo
         enddo
      enddo
      hist_6factor_fnl(i_step)=(rr_rfv_1)/max(rr_nffv/keff,1d-10)

      ! six factor formula : P_tnl (thermal neutron non leakage)
      rr_afv_2=0d0
      rr_sfv_1=0d0
      do iz=1,nz
         do ixy=1,nxy
            rr_afv_2=rr_afv_2+maxs_a_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_sfv_1=rr_sfv_1+maxs_s_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
         enddo
      enddo
      hist_6factor_tnl(i_step)=(rr_afv_2)/max(1d-10,rr_sfv_1)

      ! six factor formula : kinf
      hist_6factor_kinf(i_step)=hist_6factor_eps(i_step)*hist_6factor_p(i_step)*hist_6factor_eta(i_step)*hist_6factor_f(i_step)

      ! core average nu
      rr_nffv=0d0
      rr_ffv=0d0
      do iz=1,nz
         do ixy=1,nxy
            do ig=1,n_group
               rr_nffv=rr_nffv+nu_maxs_f_3d(ixy,iz,ig)*flux(ixy,iz,ig)*nodevolume(ixy,iz)
               rr_ffv=rr_ffv+maxs_f_3d(ixy,iz,ig)*flux(ixy,iz,ig)*nodevolume(ixy,iz)
            enddo
         enddo
      enddo
      hist_core_avg_nu(i_step)=(rr_nffv)/max(1d-10,rr_ffv)

      ! mean free path
      selfscat1=0d0 !selfscat1=0.25 ! cm-1
      selfscat2=0d0 !selfscat2=0.80 ! cm-1
      rr_rfv_1=0d0
      rr_afv_2=0d0
      rr_fv_1=0d0
      rr_fv_2=0d0
      rr_fv=0d0
      do iz=1,nz
         do ixy=1,nxy
            rr_rfv_1=rr_rfv_1+(maxs_r_3d(ixy,iz,1)+selfscat1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_afv_2=rr_afv_2+(maxs_a_3d(ixy,iz,2)+selfscat2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_fv_1=rr_fv_1+flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_fv_2=rr_fv_2+flux(ixy,iz,2)*nodevolume(ixy,iz)
            do ig=1,n_group
               rr_fv=rr_fv+flux(ixy,iz,ig)*nodevolume(ixy,iz)
            enddo
         enddo
      enddo
      hist_mfp(i_step)=(rr_fv)/max(1d-10,rr_rfv_1+rr_afv_2)
      hist_mfp_f(i_step)=(rr_fv_1)/max(1d-10,rr_rfv_1)
      hist_mfp_t(i_step)=(rr_fv_2)/max(1d-10,rr_afv_2)

      ! neutron speed & neutron density & spectral index
      rr_fv=0d0
      rr_vfv=0d0
      rr_fv_1=0d0
      rr_vfv_1=0d0
      rr_fv_2=0d0
      rr_vfv_2=0d0
      rr_v=0d0
      do iz=1,nz
         do ixy=1,nxy
            do ig=1,n_group
               rr_fv=rr_fv+flux(ixy,iz,ig)*nodevolume(ixy,iz)
               rr_vfv=rr_vfv+v_inv(ixy,iz,ig)*flux(ixy,iz,ig)*nodevolume(ixy,iz)
            enddo
            rr_fv_1=rr_fv_1+flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_vfv_1=rr_vfv_1+v_inv(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_fv_2=rr_fv_2+flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_vfv_2=rr_vfv_2+v_inv(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_v=rr_v+nodevolume(ixy,iz)
         enddo
      enddo
      hist_n_density(i_step)=rr_vfv/max(rr_v,1d-10)
      hist_n_density_f(i_step)=rr_vfv_1/max(rr_v,1d-10)
      hist_n_density_t(i_step)=rr_vfv_2/max(rr_v,1d-10)
      hist_n_speed(i_step)=rr_fv/max(1d-10,rr_vfv)
      hist_n_speed_f(i_step)=rr_fv_1/max(1d-10,rr_vfv_1)
      hist_n_speed_t(i_step)=rr_fv_2/max(1d-10,rr_vfv_2)
      hist_spectral_idx(i_step)=rr_fv_1/max(1d-10,rr_fv_2)
      hist_n_flux(i_step)=rr_fv/max(1d-10,rr_v)
      hist_n_flux_f(i_step)=rr_fv_1/max(1d-10,rr_v)
      hist_n_flux_t(i_step)=rr_fv_2/max(1d-10,rr_v)

      ! diffusion length & migration length
      rr_dfv=0d0
      rr_dfv_1=0d0
      rr_dfv_2=0d0
      rr_rfv_1=0d0
      rr_afv_2=0d0
      do iz=1,nz
         do ixy=1,nxy
            rr_rfv_1=rr_rfv_1+maxs_r_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_afv_2=rr_afv_2+maxs_a_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_dfv_1=rr_dfv_1+d_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_dfv_2=rr_dfv_2+d_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_fv_1=rr_fv_1+flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_fv_2=rr_fv_2+flux(ixy,iz,2)*nodevolume(ixy,iz)
            do ig=1,n_group
               rr_dfv=rr_dfv+d_3d(ixy,iz,ig)*flux(ixy,iz,ig)*nodevolume(ixy,iz)
               rr_fv=rr_fv+flux(ixy,iz,ig)*nodevolume(ixy,iz)
            enddo
         enddo
      enddo
      hist_dif_length(i_step)=sqrt((rr_dfv)/max(1d-10,rr_rfv_1+rr_afv_2))
      hist_dif_length_f(i_step)=sqrt((rr_dfv_1)/max(1d-10,rr_rfv_1))
      hist_dif_length_t(i_step)=sqrt((rr_dfv_2)/max(1d-10,rr_afv_2))
      hist_mig_length(i_step)=sqrt((rr_dfv_1/max(1d-10,rr_rfv_1))+(rr_dfv_2/max(1d-10,rr_afv_2)))

      ! neutron life-time * neutron generation time
      rr_vfv=0d0
      rr_vfv_1=0d0
      rr_vfv_2=0d0
      rr_sfv_1=0d0
      rr_nffv=0d0
      rr_nffv_1=0d0
      rr_nffv_2=0d0
      do iz=1,nz
         do ixy=1,nxy
            rr_vfv_1=rr_vfv_1+v_inv(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_vfv_2=rr_vfv_2+v_inv(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_sfv_1=rr_sfv_1+maxs_s_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_nffv_1=rr_nffv_1+nu_maxs_f_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_nffv_2=rr_nffv_2+nu_maxs_f_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            do ig=1,n_group
               rr_vfv=rr_vfv+v_inv(ixy,iz,ig)*flux(ixy,iz,ig)*nodevolume(ixy,iz)
               rr_nffv=rr_nffv+nu_maxs_f_3d(ixy,iz,ig)*flux(ixy,iz,ig)*nodevolume(ixy,iz)
            enddo
         enddo
      enddo
      hist_n_lifetime(i_step)=(rr_vfv)/max(1d-10,rr_nffv/keff)
      hist_n_lifetime_f(i_step)=(rr_vfv_1)/max(1d-10,rr_nffv/keff-rr_sfv_1)
      hist_n_lifetime_t(i_step)=(rr_vfv_2)/max(1d-10,rr_sfv_1)
      hist_n_gentime(i_step)=(rr_vfv)/max(1d-10,rr_nffv)
      hist_n_gentime_f(i_step)=(rr_vfv_1)/max(1d-10,rr_nffv_1)
      hist_n_gentime_t(i_step)=(rr_vfv_2)/max(1d-10,rr_nffv_2)

      ! average nuetron energy
      vtmp=hist_n_speed(i_step)/100d0
      vtmp1=hist_n_speed_f(i_step)/100d0
      vtmp2=hist_n_speed_t(i_step)/100d0
      hist_n_energy(i_step)=0.5d0*n_mass_mev*(vtmp/v_light)*(vtmp/v_light)*1d6
      hist_n_energy_f(i_step)=0.5d0*n_mass_mev*(vtmp1/v_light)*(vtmp1/v_light)*1d6
      hist_n_energy_t(i_step)=0.5d0*n_mass_mev*(vtmp2/v_light)*(vtmp2/v_light)*1d6


      return
      end subroutine cal_core_param

      subroutine cal_FA_kinf(i_step)
      use inc_3d,      only: flux
      use inc_maxs,    only: nu_maxs_f_3d, maxs_s_3d, maxs_r_3d, maxs_a_3d, d_3d, kap_maXS_f_3D
      use inc_history, only: HIST_FA_kinf, HIST_FA_kinf_XY, HIST_FA_kinf_XY_4N, HIST_FA_kinf_4N, &
                             HIST_D_1, HIST_D_2, HIST_kF_1, HIST_kF_2, HIST_r_1,                 &
                             HIST_nf_1, HIST_nf_2, HIST_a_1, HIST_a_2, HIST_f_1, HIST_f_2
      implicit none
      integer, intent(in) :: i_step
      integer :: iz, ixy, ixy_fa, m, ixy_1N, Ix_1N, Iy_1N
      integer :: ix, iy
      real(8) :: rr_v                          ! vol
      real(8) :: rr_fv_1,  rr_fv_2             ! flux * vol
      real(8) :: rr_nffv_1, rr_nffv_2          ! nu_fission * flux * vol
      real(8) :: rr_sfv_1                      ! s12 * flux * vol
      real(8) :: rr_rfv_1                      ! r * flux * vol
      real(8) :: rr_afv_2                      ! a * flux * vol
      real(8) :: norm_rr_nffv_1, norm_rr_nffv_2, norm_rr_sfv_1, norm_rr_afv_2, norm_rr_rfv_1
      real(8) :: rr_dfv_1, rr_dfv_2, rr_kffv_1,rr_afv_1,rr_kffv_2


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [cal_FA_kinf] in Mod_Save'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      rr_nffv_1 = 0d0
      rr_nffv_2 = 0d0
      rr_sfv_1  = 0d0
      rr_afv_2  = 0d0
      rr_rfv_1  = 0d0
      rr_fv_1   = 0d0
      rr_fv_2   = 0d0
      norm_rr_nffv_1 = 0d0
      norm_rr_nffv_2 = 0d0
      norm_rr_sfv_1  = 0d0
      norm_rr_afv_2  = 0d0
      norm_rr_rfv_1  = 0d0

      ! ixy_1N
      do ixy_1N=1,nxy_1N
         rr_nffv_1 = 0d0
         rr_nffv_2 = 0d0
         rr_sfv_1  = 0d0
         rr_afv_2  = 0d0
         rr_rfv_1  = 0d0
         rr_fv_1   = 0d0
         rr_fv_2   = 0d0
         do m = 1,4
            ixy = i_1Nto4N(ixy_1N,m)
            if (ixy==0) cycle
            do iz=IzFuelBot,IzFuelTop
               rr_nffv_1=rr_nffv_1+nu_maxs_f_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
               rr_nffv_2=rr_nffv_2+nu_maxs_f_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
               rr_sfv_1 =rr_sfv_1+maxs_s_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
               rr_afv_2 =rr_afv_2+maxs_a_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
               rr_rfv_1 =rr_rfv_1+maxs_r_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
               rr_fv_1  =rr_fv_1+flux(ixy,iz,1)*nodevolume(ixy,iz)
               rr_fv_2  =rr_fv_2+flux(ixy,iz,2)*nodevolume(ixy,iz)
            enddo
         enddo
         norm_rr_nffv_1 = rr_nffv_1 / rr_fv_1
         norm_rr_nffv_2 = rr_nffv_2 / rr_fv_2
         norm_rr_sfv_1  = rr_sfv_1  / rr_fv_1
         norm_rr_afv_2  = rr_afv_2  / rr_fv_2
         norm_rr_rfv_1  = rr_rfv_1  / rr_fv_1
         HIST_FA_kinf(i_step, Ixy_1N) = (norm_rr_nffv_1+norm_rr_nffv_2*norm_rr_sfv_1/norm_rr_afv_2)/norm_rr_rfv_1
      enddo

      do Ix_1N= 1,Nx_1N
         do Iy_1N = 1,Ny_1N
            ixy_1N = IxIy_1NToIxy_1N(Ix_1N, Iy_1N)
            if (ixy_1N==0) cycle
            HIST_FA_kinf_XY(i_step, Ix_1N, Iy_1N) = HIST_FA_kinf(i_step, Ixy_1N)
         enddo
      enddo

      rr_nffv_1 = 0d0
      rr_nffv_2 = 0d0
      rr_sfv_1  = 0d0
      rr_afv_2  = 0d0
      rr_rfv_1  = 0d0
      rr_fv_1   = 0d0
      rr_fv_2   = 0d0
      norm_rr_nffv_1 = 0d0
      norm_rr_nffv_2 = 0d0
      norm_rr_sfv_1  = 0d0
      norm_rr_afv_2  = 0d0
      norm_rr_rfv_1  = 0d0

      ! ixy
      do ixy=1,nxy
         rr_nffv_1 = 0d0
         rr_nffv_2 = 0d0
         rr_sfv_1  = 0d0
         rr_afv_2  = 0d0
         rr_rfv_1  = 0d0
         rr_fv_1   = 0d0
         rr_fv_2   = 0d0
         do iz=IzFuelBot,IzFuelTop
            rr_nffv_1=rr_nffv_1+nu_maxs_f_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_nffv_2=rr_nffv_2+nu_maxs_f_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_sfv_1 =rr_sfv_1+maxs_s_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_afv_2 =rr_afv_2+maxs_a_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_rfv_1 =rr_rfv_1+maxs_r_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_fv_1  =rr_fv_1+flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_fv_2  =rr_fv_2+flux(ixy,iz,2)*nodevolume(ixy,iz)
         enddo
         norm_rr_nffv_1 = rr_nffv_1 / rr_fv_1
         norm_rr_nffv_2 = rr_nffv_2 / rr_fv_2
         norm_rr_sfv_1  = rr_sfv_1  / rr_fv_1
         norm_rr_afv_2  = rr_afv_2  / rr_fv_2
         norm_rr_rfv_1  = rr_rfv_1  / rr_fv_1
         HIST_FA_kinf_4N(i_step, Ixy) = (norm_rr_nffv_1+norm_rr_nffv_2*norm_rr_sfv_1/norm_rr_afv_2)/norm_rr_rfv_1
      enddo

      do Ix= 1,Nx
         do Iy = 1,Ny
            ixy = IxIyToIxy(Ix, Iy)
            if (ixy==0) cycle
            HIST_FA_kinf_XY_4N(i_step, Ix, Iy) = HIST_FA_kinf_4N(i_step, Ixy)
         enddo
      enddo

      ! core
      rr_nffv_1 = 0d0
      rr_nffv_2 = 0d0
      rr_sfv_1  = 0d0
      rr_afv_1  = 0d0
      rr_afv_2  = 0d0
      rr_rfv_1  = 0d0
      rr_fv_1   = 0d0
      rr_fv_2   = 0d0
      norm_rr_nffv_1 = 0d0
      norm_rr_nffv_2 = 0d0
      norm_rr_sfv_1  = 0d0
      norm_rr_afv_2  = 0d0
      norm_rr_rfv_1  = 0d0
      rr_v = 0
      rr_dfv_1 = 0d0
      rr_dfv_2 = 0d0
      rr_rfv_1 = 0d0
      rr_kffv_1= 0d0
      rr_kffv_2= 0d0

      do ixy = 1,Nxy
         do iz = 1,Nz
            rr_afv_2 =rr_afv_2+maxs_a_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_afv_1 =rr_afv_1+maxs_a_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_rfv_1 =rr_rfv_1+maxs_r_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_fv_1  =rr_fv_1+flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_fv_2  =rr_fv_2+flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_dfv_1 =rr_dfv_1+d_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_dfv_2 =rr_dfv_2+d_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
         end do
      end do
      HIST_a_1 (i_step) = rr_afv_1 / rr_fv_1
      HIST_a_2 (i_step) = rr_afv_2 / rr_fv_2
      HIST_D_1 (i_step) = rr_dfv_1 / rr_fv_1
      HIST_D_2 (i_step) = rr_dfv_2 / rr_fv_2
      HIST_r_1 (i_step) = rr_rfv_1 / rr_fv_1

      rr_nffv_1 = 0d0
      rr_nffv_2 = 0d0
      rr_sfv_1  = 0d0
      rr_afv_2  = 0d0
      rr_rfv_1  = 0d0
      rr_fv_1   = 0d0
      rr_fv_2   = 0d0
      norm_rr_nffv_1 = 0d0
      norm_rr_nffv_2 = 0d0
      norm_rr_sfv_1  = 0d0
      norm_rr_afv_2  = 0d0
      norm_rr_rfv_1  = 0d0
      rr_v = 0
      rr_dfv_1 = 0d0
      rr_dfv_2 = 0d0
      rr_rfv_1 = 0d0
      rr_kffv_1= 0d0
      rr_kffv_2= 0d0

      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         do iz=IzFuelBot,IzFuelTop
            rr_nffv_1=rr_nffv_1+nu_maxs_f_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_nffv_2=rr_nffv_2+nu_maxs_f_3d(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_sfv_1 =rr_sfv_1+maxs_s_3d(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_fv_1  =rr_fv_1+flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_fv_2  =rr_fv_2+flux(ixy,iz,2)*nodevolume(ixy,iz)
            rr_v     =rr_v+nodevolume(ixy,iz)
            rr_kffv_1=rr_kffv_1+kap_maXS_f_3D(ixy,iz,1)*flux(ixy,iz,1)*nodevolume(ixy,iz)
            rr_kffv_2=rr_kffv_2+kap_maXS_f_3D(ixy,iz,2)*flux(ixy,iz,2)*nodevolume(ixy,iz)
         enddo
      enddo
      norm_rr_nffv_1 = rr_nffv_1 / rr_fv_1
      norm_rr_nffv_2 = rr_nffv_2 / rr_fv_2
      norm_rr_sfv_1  = rr_sfv_1  / rr_fv_1
      norm_rr_afv_2  = rr_afv_2  / rr_fv_2
      HIST_nf_1(i_step) = norm_rr_nffv_1
      HIST_nf_2(i_step) = norm_rr_nffv_2
      HIST_kf_1(i_step) = rr_kffv_1 / rr_fv_1 * 6.242D+12
      HIST_kf_2(i_step) = rr_kffv_2 / rr_fv_2 * 6.242D+12
      HIST_f_1 (i_step) = rr_fv_1 / rr_v
      HIST_f_2 (i_step) = rr_fv_2 / rr_v

      return
      end subroutine cal_FA_kinf

#ifdef tuan_tr_test 
      subroutine Save_Old_hex
      use inc_solver, only: betap, betap_Old, rvdelt, rvdelt_Old
      use inc_3d
      use inc_depletion
      use inc_kinetics, only: N_Group_d
!      use inc_nuclide
      use inc_th, only: PPower, PPower_Old
      use inc_time
      use inc_transient
      use inc_xs_file, only: I_Tab, I_BU
      use inc_detector
      use inc_crud, only: opt_crud, ind_bu, ind_bu_old

      use Inc_miXS
      use Inc_TPEN, ONLY: fluxf,fluxf_old
      implicit none
      integer :: ixy, iz, ig
      integer :: n, m
      integer :: l0


      PPower_old = PPower
      dT_old = dT
      do ixy=1,nxy
         do iz=1,nz
            i_tab=axialcomp(i_lp_1n(i_4nto1n(ixy)),iz)
            i_tab=new_asym_itab(i_tab,ixy)
#ifdef tuan_tr_test
#else
            do ig=1,n_group
               flux_old(ixy,iz,ig)=flux(ixy,iz,ig)
#ifdef tuan_fr
              fluxf_old(ixy,iz,ig) = fluxf(ixy,iz,ig)
#endif
            enddo
#endif
#ifdef tuan_tr_test
            do ig=1,n_group
               if (flag_transient.and.I_BU==1) then
                   flux_old(ixy,iz,ig)=flux(ixy,iz,ig)
               else
                   flux_old(ixy,iz,ig)=fluxf(ixy,iz,ig)
               endif
               
            enddo
#endif

            if (opt_mode==3) then
               do ig=1,n_group
                  if (.not.flag_transient.and.I_BU==1) then
                     flux_oold(ixy,iz,ig)=flux(ixy,iz,ig)
                  else
                     flux_oold(ixy,iz,ig)=flux_old(ixy,iz,ig)
                  endif
               enddo

               if (.not.flag_transient.and.I_BU==1) then
                  fissrc_oold(ixy,iz)=fissrc(ixy,iz)
               else
                  fissrc_oold(ixy,iz)=fissrc_old(ixy,iz)
               endif
               fissrc_old(ixy,iz)=fissrc(ixy,iz)
            endif

!#ifdef tuan_tr_test
!            do ig=1,n_group
!               if (flag_transient.and.I_BU==1) then
!                   flux_old(ixy,iz,ig)=flux(ixy,iz,ig)
!               else
!                   flux_old(ixy,iz,ig)=fluxf(ixy,iz,ig)
!               endif
!               
!            enddo
!#endif

            BU_old(ixy,iz)=BU(ixy,iz)
            if (opt_crud) then
               Ind_BU_old(ixy,iz)=Ind_BU(ixy,iz)
            endif

            l0=ixy

            if ((flag_transient).or.(.not.flag_transient.and.I_BU==1)) then
               do ig=1,n_group
                  rvdelt_old(ig,l0,iz)=rvdelt(ig,l0,iz)
               enddo
               betap_old(l0,iz)=betap(l0,iz)
            endif

            if (flag_transient) then
               do ig=1,N_Group_d
                  C_d_I_old(ixy,iz,ig)=C_d_I(ixy,iz,ig)
               enddo
            endif


         enddo
      enddo
#ifdef tuan_tr_test
!write(*,*) "WARNING: this part is temporary in transient test, need to remove for depletion"
RETURN
#endif

#ifdef tuan_fr
      N_FR_old (:,:,:) = N_FR(:,:,:)
#endif

#ifdef tuan_fr_tdep
      N_FR_T_old (:,:,:,:) = N_FR_T(:,:,:,:)
#endif

      RETURN
      END SUBROUTINE Save_Old_hex
#endif 

      
      END MODULE Mod_Save


