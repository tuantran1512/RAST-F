
      module Mod_InitCond

      use Inc_Constant
      use Inc_Geometry
      use Inc_RP
      use Mod_Alloc
#ifdef js_mpi
      use inc_parallel, only: comm
#endif

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

      CONTAINS

!!!#ifdef siarhei_delete  ! check later siarhei_rev
      subroutine IniCond_Nodal
      use Mod_Manage
      use Inc_3D
      use Inc_Control
      use Inc_DF
      use Inc_INP
      use Inc_maXS
      use Inc_Option
      use Inc_XS_File
      use Mod_GetSome
      use mod_getnode !, only: get_mtu
!     ! use mod_soldhat, only: init_subaxial_nodal
      use inc_depletion, only: n_bu
      use inc_depletion, only: cycle_bu
      use inc_depletion, only: cycle_day
      use inc_depletion, only: inp_dT
      use inc_depletion, only: inp_dBU
      use inc_depletion, only: inp_ppower
      use inc_depletion, only: inp_core_power
      use inc_depletion, only: inp_fa_power
      use inc_depletion, only: tot_mtu
      use inc_depletion, only: if_buinput
      use inc_flag, only: flag_powhist
      use inc_history, only: hist_ppower
      use inc_th, only: ppower, core_power_100
      use inc_geometry, only: core_n_fa
#ifdef tuan_fr
      use inc_tpen
#endif

#ifdef tuan_tr_test
      use Inc_TH, only: TF_In, TM_In 

      use Inc_File, only: flag_maxs
      use Mod_xsfb !, only: xsfb, xsfbhex, macroxsfbhex_mg, macroxsfbhex
      use Inc_Kinetics, only: beta_d_Tot, N_Group_d
      use Inc_Solver,   only: tbeta, kincomp
      use Inc_Expansion
      use Inc_CR,         only: flag_movecr
      use Inc_maXS
      use Inc_FA
#endif


#ifdef tuan_fr_TherEx
      use Inc_Expansion
      use Mod_GetSome, only: miXS_Func
      USE Mod_XSFB,     ONLY: XSFBHEX_MG
#endif
      implicit none
      integer :: Ixy, Ixy_1N, Iz, Ig, Ixy_FA, Ixy_RF
#ifdef jr_vver
      integer(4) :: ig1
#endif
      integer :: nstep, istep
#ifdef tuan_tr_test
      real(8) :: psil1t, fnorm
      real(8) :: phifc(2)
      integer(4) :: ih,it,ip
      integer(4) :: m,mf
      integer(4) :: l,k
      real(8) :: sum_chi, sum_chid
      integer(4) :: Ig_d
      integer(4) :: icomp, ik
      REAL(8) :: beta0(6)
      DATA beta0/0.0002584d0, 0.00152d0, 0.0013908d0, &
           0.0030704d0, 0.001102d0,0.0002584d0/
      
!      call Geometry_hex
!      call Get_CoreVolume
#endif

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [IniCond_Nodal] in Mod_InitCond'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif


#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif
#ifdef tuan_tr_test
      if (opt_nodal == 4) then
         call inicond_hex
         return
      endif
#endif

      keff = ini_keff

      Flux(:, :, 1)     = ini_Flux * 0.8D0
      Flux(:, :, 2)     = ini_Flux * 0.2D0
      Flux(0, :, :)     = D0
      Flux(Nxy+1, :, :) = D0
#ifdef tuan_fr
      if (.not.allocated(XSset_Hex ))  ALLOCATE(XSset_Hex(1:Nxy, 1:Nz,1:N_Group+6, 1:N_Group))
!      CALL Alloc(XSset_Hex, Nxy, Nz, N_Group + 6, N_Group )
#endif

#ifdef tuan_tr_test
      if (Flag_Card_Transient) then
          beta_d_Tot = D0

          DO Ixy_FA = 1, Nxy_FA
             Ixy = I_FA(Ixy_FA)
             DO Iz = IzFuelBot, IzFuelTop
                icomp = AxialComp( I_LP_1N( I_4Nto1N(Ixy) ), Iz ) 
                ik    = kincomp(icomp)
                DO Ig_d =1 ,N_Group_d 
                   beta_d_Tot(Ixy, Iz) = beta_d_Tot(Ixy, Iz) + tbeta(Ig_d, ik) 
                END DO
             END DO
          END DO
      else
          beta_d_tot(:,:) = sum(beta0)
      endif
#endif

#ifdef jr_vver
      if (opt_nodal == 4) then
         do Ixy = 1, Nxy
            Ixy_1N = I_4Nto1N(Ixy)
            do Iz = 1, Nz
               I_Tab = AxialComp( I_LP_1N( Ixy_1N ), Iz )
               I_SP   = 1
               I_BU   = 1
#ifdef tuan_fr
               if (if_mgxs ) then
                   do Ig = 1, N_Group
!                       maXS_tr_3D_MG    (Ixy, Iz, Ig) = Ref_maXS( I_Tab, Ig )
!                       D_3D_MG          (Ixy, Iz, Ig) = D1 / ( D3*maXS_tr_3D(Ixy, Iz, Ig) )
!                       maXS_a_3D_MG     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (1*N_Group + Ig) )
!                       maXS_f_3D_MG     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (2*N_Group + Ig) )
!                       nu_maXS_f_3D_MG  (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (3*N_Group + Ig) )
!                       kap_maXS_f_3D_MG (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (4*N_Group + Ig) )
!                       maXs_chi_3D_MG   (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (5*N_Group + Ig) )
!                       maXs_chid_3D_MG  (Ixy, Iz, Ig) = maXs_chi_3D (Ixy, Iz, Ig)
!                       maXS_scat_3D_MG  (Ig,:,Ixy,Iz) = Ref_smat( I_Tab, Ig, :)
!                       maXS_r_3D_MG     (Ixy, Iz, Ig) = maXS_a_3D(Ixy, Iz, Ig)
!                       do ig1 = 1, n_group
!                          if (ig1 .ne. ig) maXS_r_3D_MG(Ixy, Iz, Ig) = &
!                             maXS_r_3D_MG(Ixy, Iz, Ig) + maXS_scat_3D_MG(Ig,ig1,Ixy, Iz)
!                       enddo

                       XSset_hex(Ixy, Iz, 1,Ig) = Ref_maXS( I_Tab, Ig )
                       XSset_hex(Ixy, Iz, 2,Ig) = Ref_maXS( I_Tab, (1*N_Group + Ig) )
                       XSset_hex(Ixy, Iz, 3,Ig) = Ref_maXS( I_Tab, (2*N_Group + Ig) )
                       XSset_hex(Ixy, Iz, 4,Ig) = Ref_maXS( I_Tab, (3*N_Group + Ig) )
                       XSset_hex(Ixy, Iz, 5,Ig) = Ref_maXS( I_Tab, (4*N_Group + Ig) )
                       XSset_hex(Ixy, Iz, 6,Ig) = Ref_maXS( I_Tab, (5*N_Group + Ig) )
                       XSset_hex(Ixy, Iz, 6+ig,:) = Ref_smat( I_Tab, Ig, :)
                   enddo
               endif
#endif
#ifdef tuan_tr_test
            sum_chi = 0d0
            sum_chid= 0d0
#endif
               do Ig = 1, N_Group
                  if ( Flag_Card_maXS ) THEN
                     maXS_tr_3D    (Ixy, Iz, Ig) = Ref_maXS( I_Tab, Ig )
                     D_3D          (Ixy, Iz, Ig) = D1 / ( D3*maXS_tr_3D(Ixy, Iz, Ig) )
                     maXS_a_3D     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (1*N_Group + Ig) )
                     maXS_f_3D     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (2*N_Group + Ig) )
                     nu_maXS_f_3D  (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (3*N_Group + Ig) )
                     kap_maXS_f_3D (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (4*N_Group + Ig) )
                     maXs_chi_3D   (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (5*N_Group + Ig) )
                     !if ((Ig.EQ.1).and.(maXs_chi_3D(Ixy,Iz,Ig).EQ.0)) maXs_chi_3D(Ixy,Iz,Ig)=1d0
#ifdef tuan_tr_test
                        sum_chi = sum_chi + maXs_chi_3D   (Ixy, Iz, Ig)
                        if(allocated(ref_chid)) then
                           maXs_chid_3D  (Ixy, Iz, Ig) = Ref_Chid(I_Tab, N_group)
                           sum_chid = sum_chid + maXs_chid_3D (Ixy, Iz, Ig) 
                        else
                           maXs_chid_3D  (Ixy, Iz, Ig) = maXs_chi_3D (Ixy, Iz, Ig)
                        endif
#else 
                        maXs_chid_3D  (Ixy, Iz, Ig) = maXs_chi_3D (Ixy, Iz, Ig)

#endif
                     maXS_scat_3D  (Ig,:,Ixy,Iz) = Ref_smat( I_Tab, Ig, :)
                     maXS_r_3D     (Ixy, Iz, Ig) = maXS_a_3D(Ixy, Iz, Ig)
                     !if (N_Group>2) then
                     do ig1 = 1, n_group
                        if (ig1 .ne. ig) maXS_r_3D(Ixy, Iz, Ig) =  maXS_r_3D(Ixy, Iz, Ig) + maXS_scat_3D(Ig,ig1,Ixy, Iz)
                     enddo
                  ELSE
                     maXS_tr_3D    (Ixy, Iz, Ig) = 0.0d0
                     D_3D          (Ixy, Iz, Ig) = 0.0d0
                     maXS_a_3D     (Ixy, Iz, Ig) = 0.0d0
                     maXS_f_3D     (Ixy, Iz, Ig) = 0.0d0
                     nu_maXS_f_3D  (Ixy, Iz, Ig) = 0.0d0
                     kap_maXS_f_3D (Ixy, Iz, Ig) = 0.0d0
                     maXS_s_3D     (Ixy, Iz, Ig) = 0.0d0
                     maXS_r_3D     (Ixy, Iz, Ig) = 0.0d0
                  endif
                  maXs_s_3D     (Ixy, Iz, 1)  = maXS_scat_3D(1, 2, Ixy, Iz)
                  if ( Flag_ADF .EQV. .TRUE. ) THEN
                     ADF_Lx( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 1, Ig )
                     ADF_Rx( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 2, Ig )
                     ADF_Ly( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 3, Ig )
                     ADF_Ry( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 4, Ig )
                  ELSE
                     ADF_Lx( Ixy, Iz, Ig ) = D1
                     ADF_Rx( Ixy, Iz, Ig ) = D1
                     ADF_Ly( Ixy, Iz, Ig ) = D1
                     ADF_Ry( Ixy, Iz, Ig ) = D1
                  endif
               enddo
#ifdef tuan_tr_test
               if (sum_chi>0d0) then
                  do ig = 1,n_group
                     maXs_chi_3D   (Ixy, Iz, Ig) = maXs_chi_3D   (Ixy, Iz, Ig)/max(sum_chi,1d-30)
                  enddo
               endif
               if (sum_chid>0d0) then
                  do ig = 1,n_group
                     maXs_chid_3D (Ixy, Iz, Ig) = maXs_chid_3D (Ixy, Iz, Ig)/max(sum_chid,1d-30)
                  enddo
               endif
#endif
            enddo
         enddo
      else
         do Ixy = 1, Nxy
            Ixy_1N = I_4Nto1N(Ixy)
            do Iz = 1, Nz
               I_Tab = AxialComp( I_LP_1N( Ixy_1N ), Iz )
               I_SP   = 1
               I_BU   = 1
               do Ig = 1, N_Group
                  if ( Flag_Card_maXS ) THEN
                     maXS_tr_3D    (Ixy, Iz, Ig) = Ref_maXS( I_Tab, Ig )
                     D_3D          (Ixy, Iz, Ig) = D1 / ( D3*maXS_tr_3D(Ixy, Iz, Ig) )
                     maXS_a_3D     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (1*N_Group + Ig) )
                     maXS_f_3D     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (2*N_Group + Ig) )
                     nu_maXS_f_3D  (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (3*N_Group + Ig) )
                     kap_maXS_f_3D (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (4*N_Group + Ig) )
                     maXS_s_3D     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (5*N_Group + Ig) )
                     maXS_r_3D     (Ixy, Iz, Ig) = maXS_a_3D(Ixy, Iz, Ig) + maXS_s_3D(Ixy, Iz, Ig)
                  ELSE
                     maXS_tr_3D    (Ixy, Iz, Ig) = 0.0d0
                     D_3D          (Ixy, Iz, Ig) = 0.0d0
                     maXS_a_3D     (Ixy, Iz, Ig) = 0.0d0
                     maXS_f_3D     (Ixy, Iz, Ig) = 0.0d0
                     nu_maXS_f_3D  (Ixy, Iz, Ig) = 0.0d0
                     kap_maXS_f_3D (Ixy, Iz, Ig) = 0.0d0
                     maXS_s_3D     (Ixy, Iz, Ig) = 0.0d0
                     maXS_r_3D     (Ixy, Iz, Ig) = 0.0d0
                  endif
                  maXS_scat_3D(1, 2, Ixy, Iz) = maXS_s_3D(Ixy, Iz, 1)
                  if ( Flag_ADF .EQV. .TRUE. ) THEN
                     ADF_Lx( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 1, Ig )
                     ADF_Rx( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 2, Ig )
                     ADF_Ly( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 3, Ig )
                     ADF_Ry( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 4, Ig )
                  ELSE
                     ADF_Lx( Ixy, Iz, Ig ) = D1
                     ADF_Rx( Ixy, Iz, Ig ) = D1
                     ADF_Ly( Ixy, Iz, Ig ) = D1
                     ADF_Ry( Ixy, Iz, Ig ) = D1
                  endif
               enddo
            enddo
         enddo
      endif
#else
      do Ixy = 1, Nxy
         Ixy_1N = I_4Nto1N(Ixy)
         do Iz = 1, Nz
            I_Tab = AxialComp( I_LP_1N( Ixy_1N ), Iz )
            I_SP   = 1
            I_BU   = 1
            do Ig = 1, N_Group
               if ( Flag_Card_maXS ) THEN
                  maXS_tr_3D    (Ixy, Iz, Ig) = Ref_maXS( I_Tab, Ig )
                  D_3D          (Ixy, Iz, Ig) = D1 / ( D3*maXS_tr_3D(Ixy, Iz, Ig) )
                  maXS_a_3D     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (1*N_Group + Ig) )
                  maXS_f_3D     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (2*N_Group + Ig) )
                  nu_maXS_f_3D  (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (3*N_Group + Ig) )
                  kap_maXS_f_3D (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (4*N_Group + Ig) )
                  maXS_s_3D     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (5*N_Group + Ig) )
                  maXS_r_3D     (Ixy, Iz, Ig) = maXS_a_3D(Ixy, Iz, Ig) + maXS_s_3D(Ixy, Iz, Ig)
               ELSE
                  maXS_tr_3D    (Ixy, Iz, Ig) = 0.0d0
                  D_3D          (Ixy, Iz, Ig) = 0.0d0
                  maXS_a_3D     (Ixy, Iz, Ig) = 0.0d0
                  maXS_f_3D     (Ixy, Iz, Ig) = 0.0d0
                  nu_maXS_f_3D  (Ixy, Iz, Ig) = 0.0d0
                  kap_maXS_f_3D (Ixy, Iz, Ig) = 0.0d0
                  maXS_s_3D     (Ixy, Iz, Ig) = 0.0d0
                  maXS_r_3D     (Ixy, Iz, Ig) = 0.0d0
               endif
               maXS_scat_3D(1, 2, Ixy, Iz) = maXS_s_3D(Ixy, Iz, 1)
               if ( Flag_ADF .EQV. .TRUE. ) THEN
                  ADF_Lx( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 1, Ig )
                  ADF_Rx( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 2, Ig )
                  ADF_Ly( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 3, Ig )
                  ADF_Ry( Ixy, Iz, Ig ) = ADF_Table( I_Tab, 4, Ig )
               ELSE
                  ADF_Lx( Ixy, Iz, Ig ) = D1
                  ADF_Rx( Ixy, Iz, Ig ) = D1
                  ADF_Ly( Ixy, Iz, Ig ) = D1
                  ADF_Ry( Ixy, Iz, Ig ) = D1
               endif
            enddo
         enddo
      enddo
#endif
#ifdef tuan_tr_test
!      do l=1,nxy
!         do k=1,nz
!            do m = 1,2
!               do mf = mgb(m),mge(m)
!                  fluxf(l,k,mf)=flux(l,k,m)
!               enddo
!            enddo
!         enddo
!      enddo
#endif
      CALL Get_Reg_VF_Hex
       
#ifdef tuan_fr_TherEx
      if (Flag_Card_Thermal_Expansion) then
!          CALL Axial_Expansion_XS
          CALL miXS_func(0, 0)  
          CALL Radial_expansion

!              CALL maXS_func(0, 0)
          
!          IF(miXS_exp) THEN
!              ! update miXS base in the initial Fuel Temp and Coolant Dens.
!              CALL miXS_func(0, 0)
!       	   ! update number density base on the thermal expansion information
!       	   IF(Axial_exp) THEN
!       	       CALL Dens_update(0,0,1)
!       	   ELSE IF (Radial_exp) THEN
!       	   ! this option is not available yet
!       	       CALL Dens_update(0,0,2)
!              END IF
!       	   ! Update maXS base on the new miXS and number density
!              CALL XSFBHEX_MG
!          ELSE IF(maXS_exp) THEN
!              ! update maXS base in the initial Fuel Temp and Coolant Dens.
!              CALL maXS_func(0, 0)
!              ! Update maXS base on the new miXS  
!              CALL XSFBHEX_MG
!          END IF
!          
!          IF(Axial_exp) THEN
!              ! update XS at the mixed area
!              CALL Axial_Expansion_XS
!          ELSE IF (Radial_exp) THEN
!              ! update XS at the mixed area, not available
!!              CALL Radial_Expansion_XS
!          END IF
!          CALL XSFBHEX_MG

!          
      end if
        
#endif

      if (flag_leakage) then
         D_3D_lc = D_3D
         maXS_a_3D_lc = maXS_a_3D
         maXS_s_3D_lc = maXS_s_3D
      endif

      if (Flag_RefDF .and. OPT_Mode == 1) THEN
         do Iz = IzFuelBot, IzFuelTop
            do Ixy_RF = 1, Nxy_RF
               Ixy = I_RF(Ixy_RF)
               do Ig = 1, N_Group
                  if (I_Lx(Ixy) /= 0) THEN
#ifdef tuan_fr_crm                  
                     if ( I_FARF_1N( I_4Nto1N( I_Lx(Ixy) ) ) == 2 ) THEN
#else
                     if ( I_FARF_1N( I_4Nto1N( I_Lx(Ixy) ) ) >= 2 ) THEN
#endif
                        ADF_Lx(Ixy, Iz, Ig) = ADF_Lx(Ixy, Iz, Ig) * ADF_Rx(I_Lx(Ixy), Iz, Ig)
                     endif
                  endif
                  if (I_Rx(Ixy) /= 0) THEN
#ifdef tuan_fr_crm                  
                     if ( I_FARF_1N( I_4Nto1N( I_Rx(Ixy) ) ) == 2 ) THEN
#else
                     if ( I_FARF_1N( I_4Nto1N( I_Rx(Ixy) ) ) >= 2 ) THEN
#endif
                        ADF_Rx(Ixy, Iz, Ig) = ADF_Rx(Ixy, Iz, Ig) * ADF_Lx(I_Rx(Ixy), Iz, Ig)
                     endif
                  endif
                  if (I_Ly(Ixy) /= 0) THEN
#ifdef tuan_fr_crm                  
                     if ( I_FARF_1N( I_4Nto1N( I_Ly(Ixy) ) ) == 2 ) THEN
#else
                     if ( I_FARF_1N( I_4Nto1N( I_Ly(Ixy) ) ) >= 2 ) THEN
#endif
                        ADF_Ly(Ixy, Iz, Ig) = ADF_Ly(Ixy, Iz, Ig) * ADF_Ry(I_Ly(Ixy), Iz, Ig)
                     endif
                  endif
                  if (I_Ry(Ixy) /= 0) THEN
#ifdef tuan_fr_crm                  
                     if ( I_FARF_1N( I_4Nto1N( I_Ry(Ixy) ) ) == 2 ) THEN
#else
                     if ( I_FARF_1N( I_4Nto1N( I_Ry(Ixy) ) ) >= 2 ) THEN
#endif
                        ADF_Ry(Ixy, Iz, Ig) = ADF_Ry(Ixy, Iz, Ig) * ADF_Ly(I_Ry(Ixy), Iz, Ig)
                     endif
                  endif
               enddo
            enddo
         enddo
      endif

      call Get_Source(FisSrc, nu_maXS_f_3D, Flux)

      keff_Old  = keff
      keff_OOld = keff

      FisSrc_Old  = FisSrc
      FisSrc_OOld = FisSrc
      FisSrc_Iout = FisSrc

      do Ixy = 1, Nxy
         do Iz = 1, Nz
            do Ig = 1, N_Group
               Flux_Old ( Ixy, Iz, Ig ) = Flux    ( Ixy, Iz, Ig )
               Flux_OOld( Ixy, Iz, Ig ) = Flux_Old( Ixy, Iz, Ig )
            enddo
         enddo
      enddo

      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         do iz=izfuelbot,izfueltop
            normal_power(ixy,iz)=1d0
         enddo
      enddo

      call alloc_lscoef
      call alloc_lucoef
      call alloc_bicgs
      call alloc_flux
      call alloc_tran
#ifdef jr_vver
      if (opt_nodal == 4) then
         do Ixy = 1,Nxy
            do Iz = 1,Nz
               if (.not.allocated(maXS_r_3D_Old)) allocate(maXS_r_3D_Old(Nxy,Nz,2))
               maXS_r_3D_Old(Ixy,Iz,1) = maXS_r_3D(Ixy,Iz,1)
               maXS_r_3D_Old(Ixy,Iz,2) = maXS_r_3D(Ixy,Iz,2)
            enddo
         enddo
#ifdef tuan_tr_test
         call inicond_hex

#else
         call inicond_hex
#endif         
      endif
#endif

     ! call init_subaxial_nodal
#ifdef tuan_fr

#else
      if (flag_card_maxs) then
         ! get dBU, dT, history ...
         nstep=size(cycle_bu)
         if (.not.allocated( inp_dT           )) allocate( inp_dT           (nstep))
         if (.not.allocated( inp_dBU          )) allocate( inp_dBU          (nstep))
         if (.not.allocated( inp_ppower       )) allocate( inp_ppower       (nstep))
         if (.not.allocated( inp_core_power   )) allocate( inp_core_power   (nstep))
         if (.not.allocated( inp_fa_power     )) allocate( inp_fa_power     (nstep))
         if (.not.allocated( inp_fa_power     )) allocate( inp_fa_power     (nstep))
         inp_dT           =0d0
         inp_dBU          =0d0
         inp_ppower       =0d0
         inp_core_power   =0d0
         inp_fa_power     =0d0

         nstep=n_bu
         do istep=1,nstep
            if (flag_powhist) then
               inp_ppower(istep)=hist_ppower(istep)*1d-2
            else
               inp_ppower(istep)=ppower
            endif
            inp_core_power(istep)=core_power_100*inp_ppower(istep)
            inp_fa_power(istep)=inp_core_power(istep)/real(core_n_fa,8)
         enddo
      endif

#endif

#ifdef tuan_tr_test
!//
!//      if (n_group>2) then
!//         phifc=0d0
!//         do m = 1,2
!//            phifc(m)=1./(mge(m)-mgb(m)+1)
!//         enddo
!//         do m = 1,2
!//            alxr(m)=0
!//            alzl(m)=0
!//            alzr(m)=0
!//            do mf = mgb(m),mge(m)
!//               alxr(m)=alxr(m)+alxrf(mf)*phifc(m)
!//               alzl(m)=alzl(m)+alzlf(mf)*phifc(m)
!//               alzr(m)=alzr(m)+alzrf(mf)*phifc(m)
!//            enddo
!//            reflratf(m)  =(1-2*alxr(m))/(1+2*alxr(m))
!//            reflratzbf(m)=(1-2*alzl(m))/(1+2*alzl(m))
!//            reflratztf(m)=(1-2*alzr(m))/(1+2*alzr(m))
!//         enddo
!//         do ixy = 1,nxy
!//            do iz = 1,nz
!//               do m = 1,2
!//                  do mf = mgb(m),mge(m)
!//                     fluxf(ixy,iz,mf)=phifc(m)
!//                  enddo
!//               enddo
!//            enddo
!//         enddo
!//
!//#ifdef tuan_fr_tdep
!//       DO it = 1,6
!//           tfluxf(it,:,:,:) = fluxf(:,:,:)
!//       END DO
!//#endif
!//      IF (Flag_Card_File) THEN
!//          IF(FLag_miXS) THEN 
!//!              write(*,*) '++++++++++++++++    START miXS Intepolation in  miXS_func    ++++++++++++++++'
!//!              write(*,*) ' '
!//!!              call maXS_func(0, 0)
!//!              DO Iz = IzFuelBot, IzFuelTop
!//!                  DO Ixy_FA = 1, Nxy_FA
!//!                     Ixy = I_FA(Ixy_FA)
!//!                       CALL miXS_func(Ixy,Iz)
!//!                       CALL update_maXS(Ixy,Iz)
!//!                  ENDDO
!//!              ENDDO  
!//!              write(*,*) '++++++++++++++++    END miXS Intepolation in  miXS_func    ++++++++++++++++'
!//
!//          ELSE
!//              call maXS_func(0, 0)
!//          ENDIF
!//      ENDIF
!//         call mgto2g_hex
!//!         call xsfbhex_mg
!//      endif
!//
!//      call Get_Source(FisSrc, nu_maXS_f_3D, Flux)
!//
!//      psil1t=0d0
!//      do iz=1,nz
!//         do ixy=1,nxy
!//            psil1t=psil1t+FisSrc(ixy,iz)
!//         enddo
!//      enddo
!//#ifdef tuan_tr_test
!//      !! NEED TO MODIFY
!//!      Tot_FuelVol = Tot_FuelVol_HEX
!//      fnorm = Tot_FuelVol/psil1t
!//#else
!//      fnorm = Tot_Vol/psil1t
!//#endif
!//      do ig=1,N_Group
!//         do iz=1,nz
!//            do ixy=1,nxy
!//               Flux(ixy,iz,ig)=fnorm*Flux(ixy,iz,ig)
!//            enddo
!//         enddo
!//      enddo
!//      FisSrc=0d0
!//      call Get_Source(FisSrc, nu_maXS_f_3D, Flux)
!//
!//      keff_Old  = keff
!//      keff_OOld = keff
!//      FisSrc_Iout = FisSrc
!//      do Ig = 1, N_Group
!//         do Iz = 1, Nz
!//            do Ixy = 1, Nxy
!//               Flux_Old ( Ixy, Iz, Ig ) = Flux    ( Ixy, Iz, Ig )
!//            enddo
!//         enddo
!//      enddo
!//
!//      do ig=1,2
!//         do iz =1,nz
!//            do ih =1,nassy
!//               hflx(ig,ih,iz)=flux(1,1,ig)
!//            enddo
!//         enddo
!//      enddo
!//      do ig=1,n_group
!//         do iz=1,nz
!//            do ip = 1,ncorn_hex
!//               pflx(ig,ip,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
!//            enddo
!//         enddo
!//      enddo
!//
!//      do iz=1,nz
!//         do ih=1,nassy
!//            do it=1,ntph
!//               do ig=1,n_group
!//                  cnto(ig,it,ih,iz)=0.25*flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
!//               enddo
!//            enddo
!//            do ig=1,n_group
!//               cntzo(ig,1,ih,iz)=0.25*flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
!//               cntzo(ig,2,ih,iz)=0.25*flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
!//            enddo
!//         enddo
!//      enddo
!//      do iz=1,nz
!//         do ih=1,nassy
!//            do ig=1,n_group
!//               fhflx(ig,ih,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)/max(flux(1,1,igc(ig)),1d-30)
!//               fohflx(ig,ih,iz)=fhflx(ig,ih,iz)
!//               hflxf(ig,ih,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
!//               do it=1,6
!//                  aflx(ig,it,ih,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
!//                  fcnto(ig,it,ih,iz)=fhflx(ig,ih,iz)
!//                  focnto(ig,it,ih,iz)=fhflx(ig,ih,iz)
!//               enddo
!//               do it=1,2
!//                  fcntzo(ig,it,ih,iz)=fhflx(ig,ih,iz)
!//                  focntzo(ig,it,ih,iz)=fhflx(ig,ih,iz)
!//               enddo
!//            enddo
!//         enddo
!//      enddo
!//
!//
!//      if (.not.flag_thfb) then
!//         t_mod=tm_in+degtok
!//         t_fuel=tf_in+degtok 
!//      endif
!//      
!//      lupscat = n_group
!//      do ig = 1,n_group
!//         IF(iscatob(ig).LT.ig) lupscat=MIN(lupscat,iscatob(ig)) !! for TPENDriver
!//      enddo
!//
#endif

      return
      end subroutine IniCond_Nodal
!!!#endif 

#ifdef jr_vver
      subroutine inicond_hex
      use Mod_Manage
      use Inc_3D
      use Inc_Control
      use Inc_DF
      use Inc_INP
      use Inc_maXS
      use Inc_Option
      use Mod_GetSome
      use Mod_GetNode
      use Inc_maXS
      use Inc_FA
      use Inc_TPEN
!      use Inc_Solver, only: nkincomp
      use Mod_XSFB, only: mgto2g_hex
      use inc_flag
#ifdef tuan_tr_test
!      use mod_soldhat, only: init_subaxial_nodal
      use inc_depletion, only: n_bu
      use inc_depletion, only: cycle_bu
      use inc_depletion, only: inp_dT
      use inc_depletion, only: inp_dBU
      use inc_depletion, only: inp_ppower
      use inc_depletion, only: inp_core_power
      use inc_depletion, only: inp_fa_power
!      use inc_flag, only: flag_powhist
!      use inc_flag, only: flag_thfb
      use Inc_TH, only: TF_In, TM_In 
      use inc_history, only: hist_ppower
      use inc_th, only: ppower, core_power_100
      use inc_geometry, only: core_n_fa
      use Inc_File, only: flag_maxs
      use Mod_xsfb !, only: xsfb, xsfbhex, macroxsfbhex_mg, macroxsfbhex
      use Inc_Kinetics, only: beta_d_Tot, N_Group_d
      use Inc_Solver,   only: tbeta, kincomp
      use Inc_Expansion
      use Inc_CR,         only: flag_movecr

#endif
      implicit none
      integer :: ixy
      integer :: iz
      integer :: ig
      real(8) :: psil1t, fnorm
      real(8) :: phifc(2)
      integer(4) :: m,mf
      integer(4) :: ih,it,ip
#ifdef tuan_tr_test
      integer(4) :: ig1
      integer :: Ixy_1N, Ixy_FA, Ixy_RF
      integer :: nstep, istep
      integer(4) :: l,k
      real(8) :: sum_chi, sum_chid
      integer(4) :: Ig_d
      integer(4) :: icomp, ik
      REAL(8) :: beta0(6)
      DATA beta0/0.0002584d0, 0.00152d0, 0.0013908d0, &
           0.0030704d0, 0.001102d0,0.0002584d0/
      
      call Geometry_hex
      call Get_CoreVolume
#endif
      !IF (extth .AND. ssinit .AND. .NOT.ssfirst) RETURN
      !IF (extth .AND. rstrt .AND. .NOT.qs1d) RETURN


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [inicond_hex] in Mod_InitCond'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

#ifdef tuan_tr_test
      keff = ini_keff
      
      Flux(:, :, 1)     = ini_Flux * 0.8D0
      Flux(:, :, 2)     = ini_Flux * 0.2D0
      Flux(0, :, :)     = D0
      Flux(Nxy+1, :, :) = D0

      if (.not.allocated(XSset_Hex ))  ALLOCATE(XSset_Hex(1:Nxy, 1:Nz,1:N_Group+6, 1:N_Group))

      if (Flag_Card_Transient) then
          beta_d_Tot = D0

          DO Ixy_FA = 1, Nxy_FA
             Ixy = I_FA(Ixy_FA)
             DO Iz = IzFuelBot, IzFuelTop
                icomp = AxialComp( I_LP_1N( I_4Nto1N(Ixy) ), Iz ) 
                ik    = kincomp(icomp)
                DO Ig_d =1 ,N_Group_d 
                   beta_d_Tot(Ixy, Iz) = beta_d_Tot(Ixy, Iz) + tbeta(Ig_d, ik) 
                END DO
             END DO
          END DO
      else
          beta_d_tot(:,:) = sum(beta0)
      endif

      do Ixy = 1, Nxy
         Ixy_1N = I_4Nto1N(Ixy)
         do Iz = 1, Nz
            I_Tab = AxialComp( I_LP_1N( Ixy_1N ), Iz )
            I_SP   = 1
            I_BU   = 1
            sum_chi = 0d0
            sum_chid= 0d0

            if (if_mgxs ) then
                do Ig = 1, N_Group
                    XSset_hex(Ixy, Iz, 1,Ig) = Ref_maXS( I_Tab, Ig )
                    XSset_hex(Ixy, Iz, 2,Ig) = Ref_maXS( I_Tab, (1*N_Group + Ig) )
                    XSset_hex(Ixy, Iz, 3,Ig) = Ref_maXS( I_Tab, (2*N_Group + Ig) )
                    XSset_hex(Ixy, Iz, 4,Ig) = Ref_maXS( I_Tab, (3*N_Group + Ig) )
                    XSset_hex(Ixy, Iz, 5,Ig) = Ref_maXS( I_Tab, (4*N_Group + Ig) )
                    XSset_hex(Ixy, Iz, 6,Ig) = Ref_maXS( I_Tab, (5*N_Group + Ig) )
                    XSset_hex(Ixy, Iz, 6+ig,:) = Ref_smat( I_Tab, Ig, :)
                enddo
            endif

            do Ig = 1, N_Group
               if ( Flag_Card_maXS ) THEN
                  maXS_tr_3D    (Ixy, Iz, Ig) = Ref_maXS( I_Tab, Ig )
                  D_3D          (Ixy, Iz, Ig) = D1 / ( D3*maXS_tr_3D(Ixy, Iz, Ig) )
                  maXS_a_3D     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (1*N_Group + Ig) )
                  maXS_f_3D     (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (2*N_Group + Ig) ) 
                  nu_maXS_f_3D  (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (3*N_Group + Ig) )  
                  kap_maXS_f_3D (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (4*N_Group + Ig) )
                  maXs_chi_3D   (Ixy, Iz, Ig) = Ref_maXS( I_Tab, (5*N_Group + Ig) )
                  sum_chi = sum_chi + maXs_chi_3D   (Ixy, Iz, Ig)

                  if(allocated(ref_chid)) then
                     maXs_chid_3D  (Ixy, Iz, Ig) = Ref_Chid(I_Tab, N_group)
                     sum_chid = sum_chid + maXs_chid_3D (Ixy, Iz, Ig) 
                  else
                     maXs_chid_3D  (Ixy, Iz, Ig) = maXs_chi_3D (Ixy, Iz, Ig)
                  endif
                  maXS_scat_3D  (Ig,:,Ixy,Iz) = Ref_smat( I_Tab, Ig, :)
                  maXS_r_3D     (Ixy, Iz, Ig) = maXS_a_3D(Ixy, Iz, Ig)
                  !if (N_Group>2) then
                  do ig1 = 1, n_group
                     if (ig1 .ne. ig) maXS_r_3D(Ixy, Iz, Ig) =  maXS_r_3D(Ixy, Iz, Ig) + maXS_scat_3D(Ig,ig1,Ixy, Iz)
                  enddo
               ELSE
                  maXS_tr_3D    (Ixy, Iz, Ig) = 0.0d0 
                  D_3D          (Ixy, Iz, Ig) = 0.0d0 
                  maXS_a_3D     (Ixy, Iz, Ig) = 0.0d0 
                  maXS_f_3D     (Ixy, Iz, Ig) = 0.0d0 
                  nu_maXS_f_3D  (Ixy, Iz, Ig) = 0.0d0 
                  kap_maXS_f_3D (Ixy, Iz, Ig) = 0.0d0 
                  maXS_s_3D     (Ixy, Iz, Ig) = 0.0d0 
                  maXS_r_3D     (Ixy, Iz, Ig) = 0.0d0 
               endif
               maXs_s_3D     (Ixy, Iz, 1)  = maXS_scat_3D(1, 2, Ixy, Iz)
            enddo
            if (sum_chi>0d0) then
               do ig = 1,n_group
                  maXs_chi_3D   (Ixy, Iz, Ig) = maXs_chi_3D   (Ixy, Iz, Ig)/max(sum_chi,1d-30)
               enddo
            endif
            if (sum_chid>0d0) then
               do ig = 1,n_group
                  maXs_chid_3D (Ixy, Iz, Ig) = maXs_chid_3D (Ixy, Iz, Ig)/max(sum_chid,1d-30)
               enddo
            endif
         enddo
      enddo

      do l=1,nxy
         do k=1,nz
            do m = 1,2
               do mf = mgb(m),mge(m)
                  fluxf(l,k,mf)=flux(l,k,m)
               enddo
            enddo
         enddo
      enddo

       
#ifdef tuan_fr_TherEx
      if (Flag_Card_Thermal_Expansion) then
!          CALL Axial_Expansion_XS
          CALL miXS_func(0, 0)  
          CALL Radial_expansion

!              CALL maXS_func(0, 0)
          
!          IF(miXS_exp) THEN
!              ! update miXS base in the initial Fuel Temp and Coolant Dens.
!              CALL miXS_func(0, 0)
!       	   ! update number density base on the thermal expansion information
!       	   IF(Axial_exp) THEN
!       	       CALL Dens_update(0,0,1)
!       	   ELSE IF (Radial_exp) THEN
!       	   ! this option is not available yet
!       	       CALL Dens_update(0,0,2)
!              END IF
!       	   ! Update maXS base on the new miXS and number density
!              CALL XSFBHEX_MG
!          ELSE IF(maXS_exp) THEN
!              ! update maXS base in the initial Fuel Temp and Coolant Dens.
!              CALL maXS_func(0, 0)
!              ! Update maXS base on the new miXS  
!              CALL XSFBHEX_MG
!          END IF
!          
!          IF(Axial_exp) THEN
!              ! update XS at the mixed area
!              CALL Axial_Expansion_XS
!          ELSE IF (Radial_exp) THEN
!              ! update XS at the mixed area, not available
!!              CALL Radial_Expansion_XS
!          END IF
!          CALL XSFBHEX_MG

!          
      end if
        
#endif

      IF (Flag_Card_CR) THEN
!          CALL Get_Reg_VF_Hex
      ENDIF
!//       IF (Flag_Card_File) THEN
!//           IF(FLag_miXS) THEN 
!//               write(*,*) '++++++++++++++++    START miXS Intepolation in  miXS_func    ++++++++++++++++'
!//               write(*,*) ' '
!//               DO Iz = IzFuelBot, IzFuelTop
!//                   DO Ixy_FA = 1, Nxy_FA
!//                      Ixy = I_FA(Ixy_FA)
!// !                       CALL miXS_func(Ixy,Iz)
!// !                       CALL update_maXS(Ixy,Iz)
!//                   ENDDO
!//               ENDDO  
!//               write(*,*) '++++++++++++++++    END miXS Intepolation in  miXS_func    ++++++++++++++++'
!// 
!//           ELSE
!//               call maXS_func(0, 0)
!//           ENDIF
!//       ENDIF
      
!      if (.not.flag_maxs) then
!         call xsfbhex_mg !! = xsfb
!      endif
      if (flag_leakage) then
         D_3D_lc = D_3D 
         maXS_a_3D_lc = maXS_a_3D
         maXS_s_3D_lc = maXS_s_3D
      endif
     

      call Get_Source(FisSrc, nu_maXS_f_3D, Flux)
      
      keff_Old  = keff
      keff_OOld = keff
      
      FisSrc_Old  = FisSrc
      FisSrc_OOld = FisSrc
      FisSrc_Iout = FisSrc
      
      do Ixy = 1, Nxy
         do Iz = 1, Nz
            do Ig = 1, N_Group
               Flux_Old ( Ixy, Iz, Ig ) = Flux    ( Ixy, Iz, Ig )
               Flux_OOld( Ixy, Iz, Ig ) = Flux_Old( Ixy, Iz, Ig )
            enddo
         enddo
      enddo   

      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         do iz=izfuelbot,izfueltop
            normal_power(ixy,iz)=1d0
         enddo
      enddo
      
      call alloc_lscoef
      call alloc_lucoef
#ifdef tuan_tr_test

!#else      
      call alloc_bicgs
      call alloc_flux
      call alloc_tran
#endif
!      do Ixy = 1,Nxy
!         do Iz = 1,Nz
!            if (.not.allocated(maXS_r_3D_Old)) allocate(maXS_r_3D_Old(Nxy,Nz,2))
!            maXS_r_3D_Old(Ixy,Iz,1) = maXS_r_3D(Ixy,Iz,1)
!            maXS_r_3D_Old(Ixy,Iz,2) = maXS_r_3D(Ixy,Iz,2)
!         enddo
!      enddo

!      call init_subaxial_nodal
#ifdef tuan_fr

#else
      if (flag_card_maxs) then
         ! get dBU, dT, history ...
         nstep=size(cycle_bu)
         if (.not.allocated( inp_dT           )) allocate( inp_dT           (nstep))
         if (.not.allocated( inp_dBU          )) allocate( inp_dBU          (nstep))
         if (.not.allocated( inp_ppower       )) allocate( inp_ppower       (nstep))
         if (.not.allocated( inp_core_power   )) allocate( inp_core_power   (nstep))
         if (.not.allocated( inp_fa_power     )) allocate( inp_fa_power     (nstep))
         if (.not.allocated( inp_fa_power     )) allocate( inp_fa_power     (nstep))
         inp_dT           =0d0
         inp_dBU          =0d0
         inp_ppower       =0d0
         inp_core_power   =0d0
         inp_fa_power     =0d0
         
         nstep=n_bu
         do istep=1,nstep
            if (flag_powhist) then
               inp_ppower(istep)=hist_ppower(istep)*1d-2
            else
               inp_ppower(istep)=ppower
            endif
            inp_core_power(istep)=core_power_100*inp_ppower(istep)
            inp_fa_power(istep)=inp_core_power(istep)/real(core_n_fa,8)
         enddo
      endif
#endif

#endif      

      if (n_group>2) then
         phifc=0d0
         do m = 1,2
            phifc(m)=1./(mge(m)-mgb(m)+1)
         enddo
         do m = 1,2
            alxr(m)=0
            alzl(m)=0
            alzr(m)=0
            do mf = mgb(m),mge(m)
               alxr(m)=alxr(m)+alxrf(mf)*phifc(m)
               alzl(m)=alzl(m)+alzlf(mf)*phifc(m)
               alzr(m)=alzr(m)+alzrf(mf)*phifc(m)
            enddo
            reflratf(m)  =(1-2*alxr(m))/(1+2*alxr(m))
            reflratzbf(m)=(1-2*alzl(m))/(1+2*alzl(m))
            reflratztf(m)=(1-2*alzr(m))/(1+2*alzr(m))
         enddo
         do ixy = 1,nxy
            do iz = 1,nz
               do m = 1,2
                  do mf = mgb(m),mge(m)
                     fluxf(ixy,iz,mf)=phifc(m)
                  enddo
               enddo
            enddo
         enddo

#ifdef tuan_fr_tdep
       DO it = 1,6
           tfluxf(it,:,:,:) = fluxf(:,:,:)
       END DO
#endif

      IF (Flag_Card_File) THEN
          IF(FLag_miXS) THEN 
!              write(*,*) '++++++++++++++++    START miXS Intepolation in  miXS_func    ++++++++++++++++'
!              write(*,*) ' '
!!              call maXS_func(0, 0)
!              DO Iz = IzFuelBot, IzFuelTop
!                  DO Ixy_FA = 1, Nxy_FA
!                     Ixy = I_FA(Ixy_FA)
!                       CALL miXS_func(Ixy,Iz)
!                       CALL update_maXS(Ixy,Iz)
!                  ENDDO
!              ENDDO  
!              write(*,*) '++++++++++++++++    END miXS Intepolation in  miXS_func    ++++++++++++++++'

          ELSE
              call maXS_func(0, 0)
          ENDIF
          
      ENDIF
      IF (Flag_Card_CR) THEN
          CALL Get_Reg_VF_Hex
      ENDIF
      
         call mgto2g_hex
!         call xsfbhex_mg
      endif

      call Get_Source(FisSrc, nu_maXS_f_3D, Flux)

      psil1t=0d0
      do iz=1,nz
         do ixy=1,nxy
            psil1t=psil1t+FisSrc(ixy,iz)
         enddo
      enddo
#ifdef tuan_tr_test
      !! NEED TO MODIFY
!      Tot_FuelVol = Tot_FuelVol_HEX
      fnorm = Tot_FuelVol/psil1t
#else
      fnorm = Tot_Vol/psil1t
#endif
      do ig=1,N_Group
         do iz=1,nz
            do ixy=1,nxy
               Flux(ixy,iz,ig)=fnorm*Flux(ixy,iz,ig)
            enddo
         enddo
      enddo
      FisSrc=0d0
      call Get_Source(FisSrc, nu_maXS_f_3D, Flux)

      keff_Old  = keff
      keff_OOld = keff
      FisSrc_Iout = FisSrc
      do Ig = 1, N_Group
         do Iz = 1, Nz
            do Ixy = 1, Nxy
               Flux_Old ( Ixy, Iz, Ig ) = Flux    ( Ixy, Iz, Ig )
            enddo
         enddo
      enddo

      do ig=1,2
         do iz =1,nz
            do ih =1,nassy
               hflx(ig,ih,iz)=flux(1,1,ig)
            enddo
         enddo
      enddo
      do ig=1,n_group
         do iz=1,nz
            do ip = 1,ncorn_hex
               pflx(ig,ip,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
            enddo
         enddo
      enddo

      do iz=1,nz
         do ih=1,nassy
            do it=1,ntph
               do ig=1,n_group
                  cnto(ig,it,ih,iz)=0.25*flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
               enddo
            enddo
            do ig=1,n_group
               cntzo(ig,1,ih,iz)=0.25*flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
               cntzo(ig,2,ih,iz)=0.25*flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
            enddo
         enddo
      enddo
      do iz=1,nz
         do ih=1,nassy
            do ig=1,n_group
               fhflx(ig,ih,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)/max(flux(1,1,igc(ig)),1d-30)
               fohflx(ig,ih,iz)=fhflx(ig,ih,iz)
               hflxf(ig,ih,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
               do it=1,6
                  aflx(ig,it,ih,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
                  fcnto(ig,it,ih,iz)=fhflx(ig,ih,iz)
                  focnto(ig,it,ih,iz)=fhflx(ig,ih,iz)
               enddo
               do it=1,2
                  fcntzo(ig,it,ih,iz)=fhflx(ig,ih,iz)
                  focntzo(ig,it,ih,iz)=fhflx(ig,ih,iz)
               enddo
            enddo
         enddo
      enddo
#ifdef tuan_tr_test

      if (.not.flag_thfb) then
         t_mod=tm_in+degtok
         t_fuel=tf_in+degtok 
      endif
      
      lupscat = n_group
      do ig = 1,n_group
         IF(iscatob(ig).LT.ig) lupscat=MIN(lupscat,iscatob(ig)) !! for TPENDriver
      enddo
#endif
      return
      end subroutine inicond_hex
#endif


!!!#ifdef siarhei_delete  ! check later siarhei_rev
      subroutine IniCond_TH
      use mod_thdriver, only: setth , initth
      use inc_th, only: tf_In, tm_In, flag_th_init_first
      use inc_th, only: flag_th_chanwise, i_chan_1n, chanwise_t_inlet
!      use mod_interface
      use inc_3d, only: t_fuel, t_mod
      use inc_flag
      use inc_inp
#ifdef js_mpc
!      use link_th1d, only: initialize_th1d
#endif
#ifdef js_mpi
      use mod_parallel, only: allreduce
#endif
      implicit none
      integer :: ixy, ixy_fa, iz
      real(8) :: t_inlet

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [IniCond_TH] in Mod_InitCond'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

#ifdef tuan_tr_test
!Write(*,*) 'modify here for transient, need to be check for TH calculation'
          call IniCond_TH_HEX
         return
#endif

#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      if (flag_th_init_first) return

      if (.not.flag_thfb) then
         t_mod=tm_in+degtok
         do ixy_fa=1,nxy_fa
            ixy=i_fa(ixy_fa)
            if (flag_th_chanwise) then
               t_inlet=chanwise_t_inlet(i_chan_1n(i_4nto1n(ixy)))
               t_mod(ixy,iz)=t_inlet+degtok
            endif
            do iz=izfuelbot,izfueltop
               t_fuel(ixy,iz)=tf_in+degtok
            enddo
         enddo
         return
      endif

      if (.not.flag_thfb) return

      call setth
!  #ifdef siarhei_fr
!       if (flag_THFB) then ! @$^ Siarhei_FR
! #endif
         call initth
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
! #ifdef siarhei_fr
!       end if
! #endif

       call inicond_tfuel
!      call interface_th
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
      flag_th_init_first=.true.

      return
      end subroutine IniCond_TH
!!!#endif 
#ifdef tuan_tr_test
      subroutine IniCond_TH_HEX
      use mod_thdriver, only: setth, initth
      use inc_th, only: tf_In, tm_In, flag_th_init_first
      use inc_th, only: flag_th_chanwise, i_chan_1n, chanwise_t_inlet
      use mod_interface
      use inc_3d, only: t_fuel, t_mod
      use inc_flag
      use inc_inp
      use inc_tpen, only: flag_init_thhex
      use Mod_THdriver, only: ssth
      implicit none
      integer :: ixy, ixy_fa, iz
      real(8) :: t_inlet

      if (flag_th_init_first) return

      if (.not.flag_thfb) then
         t_mod=tm_in+degtok
         do ixy_fa=1,nxy_fa
            ixy=i_fa(ixy_fa)
            if (flag_th_chanwise) then
               t_inlet=chanwise_t_inlet(i_chan_1n(i_4nto1n(ixy)))
               t_mod(ixy,iz)=t_inlet+degtok
            endif
            do iz=izfuelbot,izfueltop
               t_fuel(ixy,iz)=tf_in+degtok
            enddo   
         enddo
         return
      endif
       
      if (.not.flag_thfb) return

      call setth
      call initth
      call inicond_tfuel
      call interface_th
      flag_th_init_first=.true.
      flag_init_thhex = .true.
      call ssth
      flag_init_thhex = .false.
      return
      end subroutine IniCond_TH_HEX 
#endif      

      

#ifdef siarhei_delete 
      subroutine IniCond_Depletion
      use Inc_Nuclide
      use Inc_XS_File, only: I_Tab, Ini_NumDen, N_XS_Table
      use Mod_GetNode, only: new_asym_itab
      implicit none
      integer :: Ixy, Iz

      call Alloc( N_U34  , Nxy , Nz )
      call Alloc( N_U35  , Nxy , Nz )
      call Alloc( N_U38  , Nxy , Nz )
      call Alloc( N_U34_Old  , Nxy , Nz )
      call Alloc( N_U35_Old  , Nxy , Nz )
      call Alloc( N_U38_Old  , Nxy , Nz )

      call Alloc( N_Gd52 , Nxy , Nz )
      call Alloc( N_Gd54 , Nxy , Nz )
      call Alloc( N_Gd55 , Nxy , Nz )
      call Alloc( N_Gd56 , Nxy , Nz )
      call Alloc( N_Gd57 , Nxy , Nz )
      call Alloc( N_Gd58 , Nxy , Nz )
      call Alloc( N_Gd60 , Nxy , Nz )
      call Alloc( N_Gd52_Old , Nxy , Nz )
      call Alloc( N_Gd54_Old , Nxy , Nz )
      call Alloc( N_Gd55_Old , Nxy , Nz )
      call Alloc( N_Gd56_Old , Nxy , Nz )
      call Alloc( N_Gd57_Old , Nxy , Nz )
      call Alloc( N_Gd58_Old , Nxy , Nz )
      call Alloc( N_Gd60_Old , Nxy , Nz )

      do Iz = 1, Nz
         do Ixy = 1, Nxy
            I_Tab = AxialComp( I_LP_1N( I_4Nto1N(Ixy) ), Iz )
            I_Tab = new_asym_itab(I_Tab,Ixy)
            N_U34 (Ixy, Iz) = Ini_NumDen( I_Tab, 1  )
            N_U35 (Ixy, Iz) = Ini_NumDen( I_Tab, 2  )
            N_U38 (Ixy, Iz) = Ini_NumDen( I_Tab, 3  )
            N_Gd52(Ixy, Iz) = Ini_NumDen( I_Tab, 4  )
            N_Gd54(Ixy, Iz) = Ini_NumDen( I_Tab, 5  )
            N_Gd55(Ixy, Iz) = Ini_NumDen( I_Tab, 6  )
            N_Gd56(Ixy, Iz) = Ini_NumDen( I_Tab, 7  )
            N_Gd57(Ixy, Iz) = Ini_NumDen( I_Tab, 8  )
            N_Gd58(Ixy, Iz) = Ini_NumDen( I_Tab, 9  )
            N_Gd60(Ixy, Iz) = Ini_NumDen( I_Tab, 10 )

            N_U34_Old (Ixy, Iz) = N_U34 (Ixy, Iz)
            N_U35_Old (Ixy, Iz) = N_U35 (Ixy, Iz)
            N_U38_Old (Ixy, Iz) = N_U38 (Ixy, Iz)
            N_Gd52_Old(Ixy, Iz) = N_Gd52(Ixy, Iz)
            N_Gd54_Old(Ixy, Iz) = N_Gd54(Ixy, Iz)
            N_Gd55_Old(Ixy, Iz) = N_Gd55(Ixy, Iz)
            N_Gd56_Old(Ixy, Iz) = N_Gd56(Ixy, Iz)
            N_Gd57_Old(Ixy, Iz) = N_Gd57(Ixy, Iz)
            N_Gd58_Old(Ixy, Iz) = N_Gd58(Ixy, Iz)
            N_Gd60_Old(Ixy, Iz) = N_Gd60(Ixy, Iz)
         enddo
      enddo
      do I_Tab = 1, N_XS_Table
         if ( Ini_NumDen( I_Tab, 4 ) > D1 ) THEN
            Flag_BP(I_Tab) = .TRUE.
         endif
      enddo

      return
      end subroutine IniCond_Depletion
#endif 


#ifdef siarhei_delete 
      subroutine Jump_in
      use Inc_Nuclide
      use Inc_XS_File, only: I_Tab
      use inc_xs_file, only: bu_ref_nobp,bu_ref_wbp,n_bu_ref,n_bu_ref_nobp,type_tab,N_BU_Ref_wBP
      use inc_option, only: opt_jumpin
      use inc_3d, only: bu,bu_old
      use inc_flag, only : opt_pjumpin
      use Mod_GetNode, only: new_asym_itab
      implicit none
      integer :: Ixy,Iy,Iy_1N,Iz,Isum
      integer(4) :: i_1n,i,i1,i2,ii, nline
      real(8) :: frac1,frac2
      logical(4) :: if_fuel,if_bp,if_exist
      real(8) :: buff_bu

#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      if(.not.opt_jumpin) return

      nline=1
      if (flag_jumpin_2d == 1 .AND. flag_jumpin_4n == 0) then
         call alloc(jump_bu,Nxy_1N,1)
         jump_bu = 0d0
         do Iy_1N = 1, Ny_1N
            if ( Ix_StartFA_y_1N(Iy_1N) == 0 ) cycle
            i1 = IxIy_1NtoIxy_1N( Ix_StartFA_y_1N(Iy_1N), Iy_1N )
            i2 = IxIy_1NtoIxy_1N( Ix_endFA_y_1N  (Iy_1N), Iy_1N )
            jump_bu(i1:i2,1) = inp_jump_bu(1:i2-i1+1,nline)
            nline=nline+1
         enddo
      elseif (flag_jumpin_2d == 0 .AND. flag_jumpin_4n == 0) then
         call alloc(jump_bu,Nxy_1N,Nz)
         jump_bu = 0d0
         do Iz = IzFuelBot, IzFuelTop
            do Iy_1N = 1, Ny_1N
               if ( Ix_StartFA_y_1N(Iy_1N) == 0 ) cycle
               i1 = IxIy_1NtoIxy_1N( Ix_StartFA_y_1N(Iy_1N), Iy_1N )
               i2 = IxIy_1NtoIxy_1N( Ix_endFA_y_1N  (Iy_1N), Iy_1N )
               jump_bu(i1:i2,Iz) = inp_jump_bu(1:i2-i1+1,nline)
               nline=nline+1
            enddo
         enddo
      elseif (flag_jumpin_2d == 1 .AND. flag_jumpin_4n == 1) then
         call alloc(jump_bu,Nxy,1)
         jump_bu = 0d0
         do Iy = 1, Ny
            if ( Ix_StartFA_y_1N( IyToIy_1N(Iy) ) == 0 ) cycle
            i1 = IxIytoIxy( Ix_Start_y(Iy), Iy )
            i2 = IxIytoIxy( Ix_end_y  (Iy), Iy )
            jump_bu(i1:i2,1) = inp_jump_bu(1:i2-i1+1,nline)
            nline=nline+1
         enddo
      elseif (flag_jumpin_2d == 0 .AND. flag_jumpin_4n == 1) then
         call alloc(jump_bu,Nxy,Nz)
         jump_bu = 0d0
         do Iz = IzFuelBot, IzFuelTop
            do Iy = 1, Ny
               if ( Ix_StartFA_y_1N( IyToIy_1N(Iy) ) == 0 ) cycle
               i1 = IxIytoIxy( Ix_Start_y(Iy), Iy )
               i2 = IxIytoIxy( Ix_end_y  (Iy), Iy )
               jump_bu(i1:i2,Iz) = inp_jump_bu(1:i2-i1+1,nline)
               nline=nline+1
            enddo
         enddo
      endif

      write(*,'(a)')'*** Jump in below burnup'
      if     (flag_jumpin_2d == 1 .AND. flag_jumpin_4n == 0) then
         Isum = 0
         do Iy_1N = 1, Ny_1N
            write(*, '(100f8.3)') jump_bu( (Isum+1):(Isum+Nx_y_1N(Iy_1N)), 1)
            Isum = Isum + Nx_y_1N(Iy_1N)
         enddo
      elseif (flag_jumpin_2d == 0 .AND. flag_jumpin_4n == 0) then
         do Iz = 1,Nz
            Isum = 0
            do Iy_1N = 1, Ny_1N
               write(*, '(100f8.3)') jump_bu( (Isum+1):(Isum+Nx_y_1N(Iy_1N)), Iz)
               Isum = Isum + Nx_y_1N(Iy_1N)
            enddo
         enddo
      elseif (flag_jumpin_2d == 1 .AND. flag_jumpin_4n == 1) then
         Isum = 0
         do Iy = 1, Ny
            write(*, '(100f8.3)') jump_bu( (Isum+1):(Isum+Nx_y(Iy)), 1)
            Isum = Isum + Nx_y(Iy)
         enddo
      elseif (flag_jumpin_2d == 0 .AND. flag_jumpin_4n == 1) then
         do Iz = 1, Nz
            Isum = 0
            do Iy = 1, Ny
               write(*, '(100f8.3)') jump_bu( (Isum+1):(Isum+Nx_y(Iy)), Iz)
               Isum = Isum + Nx_y(Iy)
            enddo
         enddo
      endif

      if(allocated(N_U34 )) deallocate(N_U34 )
      if(allocated(N_U35 )) deallocate(N_U35 )
      if(allocated(N_U36 )) deallocate(N_U36 )
      if(allocated(N_U37 )) deallocate(N_U37 )
      if(allocated(N_U38 )) deallocate(N_U38 )
      if(allocated(N_Np37)) deallocate(N_Np37)
      if(allocated(N_Np38)) deallocate(N_Np38)
      if(allocated(N_Np39)) deallocate(N_Np39)
      if(allocated(N_Pu38)) deallocate(N_Pu38)
      if(allocated(N_Pu39)) deallocate(N_Pu39)
      if(allocated(N_Pu40)) deallocate(N_Pu40)
      if(allocated(N_Pu41)) deallocate(N_Pu41)
      if(allocated(N_Pu42)) deallocate(N_Pu42)
      if(allocated(N_Pu43)) deallocate(N_Pu43)
      if(allocated(N_Am41)) deallocate(N_Am41)
      if(allocated(N_As42)) deallocate(N_As42)
      if(allocated(N_Am42)) deallocate(N_Am42)
      if(allocated(N_Am43)) deallocate(N_Am43)
      if(allocated(N_Am44)) deallocate(N_Am44)
      if(allocated(N_Cm42)) deallocate(N_Cm42)
      if(allocated(N_Cm43)) deallocate(N_Cm43)
      if(allocated(N_Cm44)) deallocate(N_Cm44)
      if(allocated(N_I35 )) deallocate(N_I35 )
      if(allocated(N_Xe35)) deallocate(N_Xe35)
      if(allocated(N_Nd47)) deallocate(N_Nd47)
      if(allocated(N_Nd48)) deallocate(N_Nd48)
      if(allocated(N_Nd49)) deallocate(N_Nd49)
      if(allocated(N_Pm47)) deallocate(N_Pm47)
      if(allocated(N_Ps48)) deallocate(N_Ps48)
      if(allocated(N_Pm48)) deallocate(N_Pm48)
      if(allocated(N_Pm49)) deallocate(N_Pm49)
      if(allocated(N_Sm47)) deallocate(N_Sm47)
      if(allocated(N_Sm48)) deallocate(N_Sm48)
      if(allocated(N_Sm49)) deallocate(N_Sm49)
      if(allocated(N_Gd52)) deallocate(N_Gd52)
      if(allocated(N_Gd54)) deallocate(N_Gd54)
      if(allocated(N_Gd55)) deallocate(N_Gd55)
      if(allocated(N_Gd56)) deallocate(N_Gd56)
      if(allocated(N_Gd57)) deallocate(N_Gd57)
      if(allocated(N_Gd58)) deallocate(N_Gd58)
      if(allocated(N_Gd60)) deallocate(N_Gd60)

      if(allocated(N_U34_old))  deallocate(N_U34_old)
      if(allocated(N_U35_old))  deallocate(N_U35_old)
      if(allocated(N_U36_old))  deallocate(N_U36_old)
      if(allocated(N_U37_old))  deallocate(N_U37_old)
      if(allocated(N_U38_old))  deallocate(N_U38_old)
      if(allocated(N_Np37_old)) deallocate(N_Np37_old)
      if(allocated(N_Np38_old)) deallocate(N_Np38_old)
      if(allocated(N_Np39_old)) deallocate(N_Np39_old)
      if(allocated(N_Pu38_old)) deallocate(N_Pu38_old)
      if(allocated(N_Pu39_old)) deallocate(N_Pu39_old)
      if(allocated(N_Pu40_old)) deallocate(N_Pu40_old)
      if(allocated(N_Pu41_old)) deallocate(N_Pu41_old)
      if(allocated(N_Pu42_old)) deallocate(N_Pu42_old)
      if(allocated(N_Pu43_old)) deallocate(N_Pu43_old)
      if(allocated(N_Am41_old)) deallocate(N_Am41_old)
      if(allocated(N_As42_old)) deallocate(N_As42_old)
      if(allocated(N_Am42_old)) deallocate(N_Am42_old)
      if(allocated(N_Am43_old)) deallocate(N_Am43_old)
      if(allocated(N_Am44_old)) deallocate(N_Am44_old)
      if(allocated(N_Cm42_old)) deallocate(N_Cm42_old)
      if(allocated(N_Cm43_old)) deallocate(N_Cm43_old)
      if(allocated(N_Cm44_old)) deallocate(N_Cm44_old)
      if(allocated(N_I35_old))  deallocate(N_I35_old)
      if(allocated(N_Xe35_old)) deallocate(N_Xe35_old)
      if(allocated(N_Nd47_old)) deallocate(N_Nd47_old)
      if(allocated(N_Nd48_old)) deallocate(N_Nd48_old)
      if(allocated(N_Nd49_old)) deallocate(N_Nd49_old)
      if(allocated(N_Pm47_old)) deallocate(N_Pm47_old)
      if(allocated(N_Ps48_old)) deallocate(N_Ps48_old)
      if(allocated(N_Pm48_old)) deallocate(N_Pm48_old)
      if(allocated(N_Pm49_old)) deallocate(N_Pm49_old)
      if(allocated(N_Sm47_old)) deallocate(N_Sm47_old)
      if(allocated(N_Sm48_old)) deallocate(N_Sm48_old)
      if(allocated(N_Sm49_old)) deallocate(N_Sm49_old)
      if(allocated(N_Gd52_old)) deallocate(N_Gd52_old)
      if(allocated(N_Gd54_old)) deallocate(N_Gd54_old)
      if(allocated(N_Gd55_old)) deallocate(N_Gd55_old)
      if(allocated(N_Gd56_old)) deallocate(N_Gd56_old)
      if(allocated(N_Gd57_old)) deallocate(N_Gd57_old)
      if(allocated(N_Gd58_old)) deallocate(N_Gd58_old)
      if(allocated(N_Gd60_old)) deallocate(N_Gd60_old)

      IF (opt_pjumpin ) THEN
         CALL Alloc( BU              , Nxy , Nz )
         CALL Alloc( BU_Old          , Nxy , Nz )
      END IF
      call alloc(N_U34 , Nxy, Nz)
      call alloc(N_U35 , Nxy, Nz)
      call alloc(N_U36 , Nxy, Nz)
      call alloc(N_U37 , Nxy, Nz)
      call alloc(N_U38 , Nxy, Nz)
      call alloc(N_Np37, Nxy, Nz)
      call alloc(N_Np38, Nxy, Nz)
      call alloc(N_Np39, Nxy, Nz)
      call alloc(N_Pu38, Nxy, Nz)
      call alloc(N_Pu39, Nxy, Nz)
      call alloc(N_Pu40, Nxy, Nz)
      call alloc(N_Pu41, Nxy, Nz)
      call alloc(N_Pu42, Nxy, Nz)
      call alloc(N_Pu43, Nxy, Nz)
      call alloc(N_Am41, Nxy, Nz)
      call alloc(N_As42, Nxy, Nz)
      call alloc(N_Am42, Nxy, Nz)
      call alloc(N_Am43, Nxy, Nz)
      call alloc(N_Am44, Nxy, Nz)
      call alloc(N_Cm42, Nxy, Nz)
      call alloc(N_Cm43, Nxy, Nz)
      call alloc(N_Cm44, Nxy, Nz)
      call alloc(N_I35 , Nxy, Nz)
      call alloc(N_Xe35, Nxy, Nz)
      call alloc(N_Nd47, Nxy, Nz)
      call alloc(N_Nd48, Nxy, Nz)
      call alloc(N_Nd49, Nxy, Nz)
      call alloc(N_Pm47, Nxy, Nz)
      call alloc(N_Ps48, Nxy, Nz)
      call alloc(N_Pm48, Nxy, Nz)
      call alloc(N_Pm49, Nxy, Nz)
      call alloc(N_Sm47, Nxy, Nz)
      call alloc(N_Sm48, Nxy, Nz)
      call alloc(N_Sm49, Nxy, Nz)
      call alloc(N_Gd52, Nxy, Nz)
      call alloc(N_Gd54, Nxy, Nz)
      call alloc(N_Gd55, Nxy, Nz)
      call alloc(N_Gd56, Nxy, Nz)
      call alloc(N_Gd57, Nxy, Nz)
      call alloc(N_Gd58, Nxy, Nz)
      call alloc(N_Gd60, Nxy, Nz)

      call alloc(N_U34_old , Nxy, Nz)
      call alloc(N_U35_old , Nxy, Nz)
      call alloc(N_U36_old , Nxy, Nz)
      call alloc(N_U37_old , Nxy, Nz)
      call alloc(N_U38_old , Nxy, Nz)
      call alloc(N_Np37_old, Nxy, Nz)
      call alloc(N_Np38_old, Nxy, Nz)
      call alloc(N_Np39_old, Nxy, Nz)
      call alloc(N_Pu38_old, Nxy, Nz)
      call alloc(N_Pu39_old, Nxy, Nz)
      call alloc(N_Pu40_old, Nxy, Nz)
      call alloc(N_Pu41_old, Nxy, Nz)
      call alloc(N_Pu42_old, Nxy, Nz)
      call alloc(N_Pu43_old, Nxy, Nz)
      call alloc(N_Am41_old, Nxy, Nz)
      call alloc(N_As42_old, Nxy, Nz)
      call alloc(N_Am42_old, Nxy, Nz)
      call alloc(N_Am43_old, Nxy, Nz)
      call alloc(N_Am44_old, Nxy, Nz)
      call alloc(N_Cm42_old, Nxy, Nz)
      call alloc(N_Cm43_old, Nxy, Nz)
      call alloc(N_Cm44_old, Nxy, Nz)
      call alloc(N_I35_old , Nxy, Nz)
      call alloc(N_Xe35_old, Nxy, Nz)
      call alloc(N_Nd47_old, Nxy, Nz)
      call alloc(N_Nd48_old, Nxy, Nz)
      call alloc(N_Nd49_old, Nxy, Nz)
      call alloc(N_Pm47_old, Nxy, Nz)
      call alloc(N_Ps48_old, Nxy, Nz)
      call alloc(N_Pm48_old, Nxy, Nz)
      call alloc(N_Pm49_old, Nxy, Nz)
      call alloc(N_Sm47_old, Nxy, Nz)
      call alloc(N_Sm48_old, Nxy, Nz)
      call alloc(N_Sm49_old, Nxy, Nz)
      call alloc(N_Gd52_old, Nxy, Nz)
      call alloc(N_Gd54_old, Nxy, Nz)
      call alloc(N_Gd55_old, Nxy, Nz)
      call alloc(N_Gd56_old, Nxy, Nz)
      call alloc(N_Gd57_old, Nxy, Nz)
      call alloc(N_Gd58_old, Nxy, Nz)
      call alloc(N_Gd60_old, Nxy, Nz)

      do Iz = 1, Nz
         do Ixy = 1, Nxy
            i_1N= I_4Nto1N(Ixy)
            I_Tab = AxialComp( I_LP_1N( I_1N ), Iz )
            I_Tab = new_asym_itab(I_Tab,Ixy)
            if_fuel=.false.
            if_bp=.false.
            if(Type_Tab( I_Tab ) == 1) then ! Normal fuel
                if_fuel=.true.
                N_BU_Ref  = N_BU_Ref_NoBP
            ELSE if ( (Type_Tab( I_Tab ) == 2) ) THEN
                if_fuel=.true.
                if_bp=.true.
                N_BU_Ref  = N_BU_Ref_wBP(I_Tab)
            endif
            if(if_fuel) then
                if     (flag_jumpin_2d == 1 .AND. flag_jumpin_4n == 0) then
                   buff_bu = jump_bu(i_1n,1)
                elseif (flag_jumpin_2d == 0 .AND. flag_jumpin_4n == 0) then
                   buff_bu = jump_bu(i_1n,iz)
                elseif (flag_jumpin_2d == 1 .AND. flag_jumpin_4n == 1) then
                   buff_bu = jump_bu(Ixy,1)
                elseif (flag_jumpin_2d == 0 .AND. flag_jumpin_4n == 1) then
                   buff_bu = jump_bu(Ixy,Iz)
                endif
                BU(ixy,iz)=buff_bu
                BU_old(ixy,iz)=buff_bu

                frac1=0.0d0
                frac2=0.0d0
                if_exist=.false.
                if(.not.if_bp) then
                   do i=2,n_bu_ref
                       if(BU_Ref_NoBP(I_Tab,i)>buff_bu) then
                           if_exist=.true.
                           frac1=(BU_Ref_NoBP(I_Tab,i)-buff_bu)/(BU_Ref_NoBP(I_Tab,i)-BU_Ref_NoBP(I_Tab,i-1))
                           frac2=1.0d0-frac1
                           i1=i-1
                           i2=i
                           exit
                       endif
                   enddo
                else
                   do i=2,n_bu_ref
                       if(BU_Ref_wBP(i_tab,i)>buff_bu) then
                           if_exist=.true.
                           frac1=(BU_Ref_wBP(i_tab,i)-buff_bu)/(BU_Ref_wBP(i_tab,i)-BU_Ref_wBP(i_tab,i-1))
                           frac2=1.0d0-frac1
                           i1=i-1
                           i2=i
                           exit
                       endif
                   enddo
                endif
                if(.not.if_exist) then
                    write(*,*) 'Jumpin burnup is out of burnup range'
                    write(*,*) iz,ixy,i_1n,i_tab,buff_bu
                    write(*,*) if_fuel,if_bp
                    stop
                endif
                ii=1
                N_U34 (ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_U35 (ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_U36 (ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_U37 (ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_U38 (ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Np37(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Np38(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Np39(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Pu38(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Pu39(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Pu40(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Pu41(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Pu42(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Pu43(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Am41(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_As42(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Am42(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Am43(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Am44(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Cm42(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Cm43(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Cm44(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_I35 (ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Xe35(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Nd47(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Nd48(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Nd49(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Pm47(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Ps48(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Pm48(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Pm49(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Sm47(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Sm48(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Sm49(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                !N_Nd46(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                ii=ii+1
                N_Gd52(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Gd54(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Gd55(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Gd56(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Gd57(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Gd58(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab);ii=ii+1
                N_Gd60(ixy,iz) = frac1*burned_nden(ii,i1,i_tab)+frac2*burned_nden(ii,i2,i_tab)
                if(ii/=42) then
                    write(*,*) 'ii/=42'
                endif
                N_U34_old (ixy,iz) = N_U34 (ixy,iz)
                N_U35_old (ixy,iz) = N_U35 (ixy,iz)
                N_U36_old (ixy,iz) = N_U36 (ixy,iz)
                N_U37_old (ixy,iz) = N_U37 (ixy,iz)
                N_U38_old (ixy,iz) = N_U38 (ixy,iz)
                N_Np37_old(ixy,iz) = N_Np37(ixy,iz)
                N_Np38_old(ixy,iz) = N_Np38(ixy,iz)
                N_Np39_old(ixy,iz) = N_Np39(ixy,iz)
                N_Pu38_old(ixy,iz) = N_Pu38(ixy,iz)
                N_Pu39_old(ixy,iz) = N_Pu39(ixy,iz)
                N_Pu40_old(ixy,iz) = N_Pu40(ixy,iz)
                N_Pu41_old(ixy,iz) = N_Pu41(ixy,iz)
                N_Pu42_old(ixy,iz) = N_Pu42(ixy,iz)
                N_Pu43_old(ixy,iz) = N_Pu43(ixy,iz)
                N_Am41_old(ixy,iz) = N_Am41(ixy,iz)
                N_As42_old(ixy,iz) = N_As42(ixy,iz)
                N_Am42_old(ixy,iz) = N_Am42(ixy,iz)
                N_Am43_old(ixy,iz) = N_Am43(ixy,iz)
                N_Am44_old(ixy,iz) = N_Am44(ixy,iz)
                N_Cm42_old(ixy,iz) = N_Cm42(ixy,iz)
                N_Cm43_old(ixy,iz) = N_Cm43(ixy,iz)
                N_Cm44_old(ixy,iz) = N_Cm44(ixy,iz)
                N_I35_old (ixy,iz) = N_I35 (ixy,iz)
                N_Xe35_old(ixy,iz) = N_Xe35(ixy,iz)
                N_Nd47_old(ixy,iz) = N_Nd47(ixy,iz)
                N_Nd48_old(ixy,iz) = N_Nd48(ixy,iz)
                N_Nd49_old(ixy,iz) = N_Nd49(ixy,iz)
                N_Pm47_old(ixy,iz) = N_Pm47(ixy,iz)
                N_Ps48_old(ixy,iz) = N_Ps48(ixy,iz)
                N_Pm48_old(ixy,iz) = N_Pm48(ixy,iz)
                N_Pm49_old(ixy,iz) = N_Pm49(ixy,iz)
                N_Sm47_old(ixy,iz) = N_Sm47(ixy,iz)
                N_Sm48_old(ixy,iz) = N_Sm48(ixy,iz)
                N_Sm49_old(ixy,iz) = N_Sm49(ixy,iz)
                N_Gd52_old(ixy,iz) = N_Gd52(ixy,iz)
                N_Gd54_old(ixy,iz) = N_Gd54(ixy,iz)
                N_Gd55_old(ixy,iz) = N_Gd55(ixy,iz)
                N_Gd56_old(ixy,iz) = N_Gd56(ixy,iz)
                N_Gd57_old(ixy,iz) = N_Gd57(ixy,iz)
                N_Gd58_old(ixy,iz) = N_Gd58(ixy,iz)
                N_Gd60_old(ixy,iz) = N_Gd60(ixy,iz)
            endif
         enddo
      enddo

      return
      end subroutine Jump_in
#endif 


      subroutine IniCond_TFuel
      use Inc_RP, only: I_Gd_FA
      use Inc_Nuclide
      implicit none
      integer ::  Ixy, Ixy_FA, Iz
      real(8) :: aw_U, aw_O2, aw_Gd2, aw_O3, awt_Gd2O3, awt_UO2, vol


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [IniCond_TFuel] in Mod_InitCond'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call Alloc(I_Gd_FA, Nxy_FA, Nz)

      do Iz=1,Nz
         do Ixy_FA=1,Nxy_FA
            Ixy    = I_FA(Ixy_FA)
            vol    = NodeVolume(Ixy,Iz)
            aw_U   = (N_U34(Ixy,Iz)+N_U35(Ixy,Iz)+N_U38(Ixy,Iz)) * vol
            aw_Gd2 = (N_Gd52(Ixy,Iz)+N_Gd54(Ixy,Iz)+N_Gd55(Ixy,Iz)+N_Gd56(Ixy,Iz) &
                     +N_Gd57(Ixy,Iz)+N_Gd58(Ixy,Iz)+N_Gd60(Ixy,Iz)) * vol

            aw_O2  = 2d0*aw_U
            aw_O3  = aw_Gd2/2d0*3d0

            awt_Gd2O3 = (aw_Gd2+aw_O3)/max(1d-50,aw_U+aw_O2+aw_Gd2+aw_O3)
            awt_UO2   = 1d0-awt_Gd2O3

            I_Gd_FA(Ixy_FA,Iz) = (awt_Gd2O3*M_Gd2O3)/max(1d-50,awt_UO2*M_UO2+awt_Gd2O3*M_Gd2O3)
         enddo
      enddo

      return
      end subroutine IniCond_TFuel


      subroutine GuessIteration
      use Inc_Flag
      use Inc_Option
      use Inc_Control
      use Inc_FluxVar
      use Inc_Solver
      use Inc_3d, only: keff

      implicit none

      ! initialize the iteration parameters

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [GuessIteration] in Mod_InitCond'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      flagr2   = .false.
      flagl2   = .false.
      flaglinf = .false.
      flageig  = .false.
      flagerf  = .false.
      flagth   = .false.

      keff = 1d0

      eshift = 0.04D0

      if ( N_Group > 2 ) THEN
         eshift0 = 1d0
      ELSE
         eshift0 = 1d-1
      endif

      eigvs = keff + eshift0  !ssinit

      reigv = 1d0 / keff   ! reciprocal of eigenvalue

      reigvs  = 1d0 / eigvs
      reigvsd = 0d0
      errl2   = 1d0

      ! transient variable for determine the fix source equation
      rvdelt    = 0d0
      betap     = 1d0
      betap_Old = 1d0

      return
      end subroutine GuessIteration

      end module Mod_InitCond

