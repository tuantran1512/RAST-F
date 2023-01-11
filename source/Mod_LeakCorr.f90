#ifdef siarhei_delete


      MODULE Mod_LeakCorr

      USE Inc_Constant
      USE Inc_Geometry
      USE Inc_RP
      USE Inc_XS_File
      USE Inc_Flag
      use Mod_GetNode, only: new_asym_itab

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE
      real(8),parameter :: tf_max=1500.0d0
      real(8),parameter :: tm_max=620.0d0
      real(8),parameter :: tm_min=280.0d0

      CONTAINS

      SUBROUTINE XSFB_LC

      USE Inc_CR
      USE Inc_Depletion
      USE Inc_INP
      USE Mod_GetSome
      USE Inc_3D
      USE Inc_TH
      USE Inc_maXS
      use inc_xs_file

      IMPLICIT NONE

      INTEGER :: Ixy, Iz


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [XSFB_LC] in Mod_LeakCorr'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      DO Ixy = 1, Nxy
         DO Iz = 1, Nz
            IF ( Flag_Card_maXS .EQV. .TRUE. ) THEN
               CALL maXSsetFB_LC( Ixy, Iz )
            ELSE
               CALL XSsetFB_LC( Ixy, Iz )
            END IF
            IF (Flag_LC_micro ) THEN
               CALL miXSFB_LC( Ixy, Iz )
            ENDIF
            CALL maXSFB_LC( Ixy, Iz )
         END DO
      END DO

      RETURN
      END SUBROUTINE XSFB_LC

      SUBROUTINE maXSsetFB_LC(Ixy, Iz)

      USE Inc_RP, ONLY: I_LP_1N, AxialComp
      USE Inc_XS_File
      USE Mod_XSFB, ONLY: maXSsetFB

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: Ixy, Iz
      INTEGER :: Ixy_1N


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [maXSsetFB_LC] in Mod_LeakCorr'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      XSset  = D0
      HFFset = D1
      Ixy_1N = I_4Nto1N(Ixy)
      I_Tab  = AxialComp( I_LP_1N( Ixy_1N ), Iz )
      I_Tab = new_asym_itab(I_Tab,Ixy)

!js+      IF ( Iz < IzFuelBot ) THEN
!js+         I_Type = 4
!js+      ELSE IF ( Iz > IzFuelTop ) THEN
!js+         I_Type = 5
!js+      ELSE IF ( I_FARF_1N(Ixy_1N) == 1 ) THEN
!js+         I_Type = 3
!js+      ELSE  ! IF ( I_FARF_1N(Ixy_1N) >= 2 ) THEN
!js+         I_Type = 1  ! or 2
!js+      END IF
!js+
!js+      XSset(1)  = Ref_maXS(I_Tab,1)  / leak_ratio(Ixy,Iz,1)
!js+      XSset(3)  = Ref_maXS(I_Tab,3)  * leak_ratio(Ixy,Iz,3)
!js+      XSset(11) = Ref_maXS(I_Tab,11) * leak_ratio(Ixy,Iz,11)

      CALL maXSsetFB(Ixy, Iz)
#ifdef tuan_fr_crm                  
      IF ( I_FARF_1N(Ixy_1N) == 2 ) THEN
#else
      IF ( I_FARF_1N(Ixy_1N) >= 2 ) THEN
#endif
         XSset(1)  = XSset(1)  / leak_ratio(Ixy,Iz,1)
         XSset(3)  = XSset(3)  * leak_ratio(Ixy,Iz,2)
         XSset(11) = XSset(11) * leak_ratio(Ixy,Iz,3)
!         XSset(10) = XSset(10) * leak_ratio(Ixy,Iz,4)
!         XSset(8)  = XSset(8)  * leak_ratio(Ixy,Iz,4)
!         XSset(6)  = XSset(6)  * leak_ratio(Ixy,Iz,4)
      ENDIF

      RETURN
      END SUBROUTINE maXSsetFB_LC

      SUBROUTINE XSsetFB_LC(Ixy, Iz)

      USE Inc_INP, ONLY: Flag_Card_maXS
      USE Inc_RP, ONLY: I_LP_1N, AxialComp
      USE Inc_XS_File
      USE Mod_XSFB, ONLY: XSsetFB

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: Ixy, Iz
      INTEGER :: Ixy_1N


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [XSsetFB_LC] in Mod_LeakCorr'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      Ixy_1N = I_4Nto1N(Ixy)
      I_Tab  = AxialComp( I_LP_1N( Ixy_1N ), Iz )
      I_Type = Type_Tab(I_Tab)
      I_Tab = new_asym_itab(I_Tab,Ixy)

      CALL XSsetFB(Ixy, Iz)
#ifdef tuan_fr_crm                  
      IF ( I_FARF_1N(Ixy_1N) == 2 ) THEN
#else
      IF ( I_FARF_1N(Ixy_1N) >= 2 ) THEN
#endif
         XSset(1)  = XSset(1)  / leak_ratio(Ixy,Iz,1)
         XSset(3)  = XSset(3)  * leak_ratio(Ixy,Iz,2)
         XSset(11) = XSset(11) * leak_ratio(Ixy,Iz,3)
!         XSset(10) = XSset(10) * leak_ratio(Ixy,Iz,4)
!         XSset(8)  = XSset(8)  * leak_ratio(Ixy,Iz,4)
!         XSset(6)  = XSset(6)  * leak_ratio(Ixy,Iz,4)
      ENDIF

      IF ( Flag_Card_maXS .EQV. .FALSE. ) THEN
         XSset( 13:400 ) = XSset( 13:400 ) * barn
      END IF

      RETURN
      END SUBROUTINE XSsetFB_LC

      SUBROUTINE miXSFB_LC(Ixy, Iz)

      USE Inc_DF
      USE Inc_INP, ONLY: Flag_Card_maXS
      USE Inc_Kinetics
      USE Inc_miXS
      USE Inc_Option, ONLY: N_Group
      USE Inc_XS_File, ONLY: XSset

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Ixy, Iz
      INTEGER :: Ig
      INTEGER :: Ixs


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [miXSFB_LC] in Mod_LeakCorr'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF ( Flag_Card_maXS .EQV. .TRUE. ) THEN
         RETURN
      END IF

      Ig  = 0
      Ixs = 13

      DO Ig = 1, N_Group
         miXS_tr_U34 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_U35 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_U36 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_U37 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_U38 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Np37(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Np38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Np39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Pu38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Pu39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Pu40(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Pu41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Pu42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Pu43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Am41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_As42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Am42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Am43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Am44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Cm42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Cm43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Cm44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_tr_U34 (Ixy, Iz, 1) = miXS_tr_U34 (Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_U35 (Ixy, Iz, 1) = miXS_tr_U35 (Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_U36 (Ixy, Iz, 1) = miXS_tr_U36 (Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_U37 (Ixy, Iz, 1) = miXS_tr_U37 (Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_U38 (Ixy, Iz, 1) = miXS_tr_U38 (Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Np37(Ixy, Iz, 1) = miXS_tr_Np37(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Np38(Ixy, Iz, 1) = miXS_tr_Np38(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Np39(Ixy, Iz, 1) = miXS_tr_Np39(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Pu38(Ixy, Iz, 1) = miXS_tr_Pu38(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Pu39(Ixy, Iz, 1) = miXS_tr_Pu39(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Pu40(Ixy, Iz, 1) = miXS_tr_Pu40(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Pu41(Ixy, Iz, 1) = miXS_tr_Pu41(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Pu42(Ixy, Iz, 1) = miXS_tr_Pu42(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Pu43(Ixy, Iz, 1) = miXS_tr_Pu43(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Am41(Ixy, Iz, 1) = miXS_tr_Am41(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_As42(Ixy, Iz, 1) = miXS_tr_As42(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Am42(Ixy, Iz, 1) = miXS_tr_Am42(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Am43(Ixy, Iz, 1) = miXS_tr_Am43(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Am44(Ixy, Iz, 1) = miXS_tr_Am44(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Cm42(Ixy, Iz, 1) = miXS_tr_Cm42(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Cm43(Ixy, Iz, 1) = miXS_tr_Cm43(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Cm44(Ixy, Iz, 1) = miXS_tr_Cm44(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)

      DO Ig = 1, N_Group
         miXS_a_U34 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_U35 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_U36 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_U37 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_U38 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Np37(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Np38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Np39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Pu38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Pu39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Pu40(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Pu41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Pu42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Pu43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Am41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_As42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Am42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Am43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Am44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Cm42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Cm43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Cm44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_a_U34 (Ixy, Iz, 1) = miXS_a_U34 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_U35 (Ixy, Iz, 1) = miXS_a_U35 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_U36 (Ixy, Iz, 1) = miXS_a_U36 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_U37 (Ixy, Iz, 1) = miXS_a_U37 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_U38 (Ixy, Iz, 1) = miXS_a_U38 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Np37(Ixy, Iz, 1) = miXS_a_Np37(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Np38(Ixy, Iz, 1) = miXS_a_Np38(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Np39(Ixy, Iz, 1) = miXS_a_Np39(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Pu38(Ixy, Iz, 1) = miXS_a_Pu38(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Pu39(Ixy, Iz, 1) = miXS_a_Pu39(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Pu40(Ixy, Iz, 1) = miXS_a_Pu40(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Pu41(Ixy, Iz, 1) = miXS_a_Pu41(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Pu42(Ixy, Iz, 1) = miXS_a_Pu42(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Pu43(Ixy, Iz, 1) = miXS_a_Pu43(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Am41(Ixy, Iz, 1) = miXS_a_Am41(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_As42(Ixy, Iz, 1) = miXS_a_As42(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Am42(Ixy, Iz, 1) = miXS_a_Am42(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Am43(Ixy, Iz, 1) = miXS_a_Am43(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Am44(Ixy, Iz, 1) = miXS_a_Am44(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Cm42(Ixy, Iz, 1) = miXS_a_Cm42(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Cm43(Ixy, Iz, 1) = miXS_a_Cm43(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Cm44(Ixy, Iz, 1) = miXS_a_Cm44(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)

      DO Ig = 1, N_Group
         miXS_f_U34 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_U35 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_U36 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_U37 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_U38 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Np37(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Np38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Np39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Pu38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Pu39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Pu40(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Pu41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Pu42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Pu43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Am41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_As42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Am42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Am43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Am44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Cm42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Cm43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_f_Cm44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO

      DO Ig = 1, N_Group
         nu_miXS_f_U34 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_U35 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_U36 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_U37 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_U38 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Np37(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Np38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Np39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Pu38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Pu39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Pu40(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Pu41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Pu42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Pu43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Am41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_As42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Am42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Am43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Am44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Cm42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Cm43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         nu_miXS_f_Cm44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO

      DO Ig = 1, N_Group
         kap_miXS_f_U34 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_U35 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_U36 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_U37 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_U38 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Np37(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Np38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Np39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Pu38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Pu39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Pu40(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Pu41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Pu42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Pu43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Am41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_As42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Am42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Am43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Am44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Cm42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Cm43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         kap_miXS_f_Cm44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO

      DO Ig = 1, N_Group
         miXS_s_U34 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_U35 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_U36 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_U37 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_U38 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Np37(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Np38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Np39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Pu38(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Pu39(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Pu40(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Pu41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Pu42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Pu43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Am41(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_As42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Am42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Am43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Am44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Cm42(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Cm43(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Cm44(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_s_U34 (Ixy, Iz, 1) = miXS_s_U34 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_U35 (Ixy, Iz, 1) = miXS_s_U35 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_U36 (Ixy, Iz, 1) = miXS_s_U36 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_U37 (Ixy, Iz, 1) = miXS_s_U37 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_U38 (Ixy, Iz, 1) = miXS_s_U38 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Np37(Ixy, Iz, 1) = miXS_s_Np37(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Np38(Ixy, Iz, 1) = miXS_s_Np38(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Np39(Ixy, Iz, 1) = miXS_s_Np39(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Pu38(Ixy, Iz, 1) = miXS_s_Pu38(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Pu39(Ixy, Iz, 1) = miXS_s_Pu39(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Pu40(Ixy, Iz, 1) = miXS_s_Pu40(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Pu41(Ixy, Iz, 1) = miXS_s_Pu41(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Pu42(Ixy, Iz, 1) = miXS_s_Pu42(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Pu43(Ixy, Iz, 1) = miXS_s_Pu43(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Am41(Ixy, Iz, 1) = miXS_s_Am41(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_As42(Ixy, Iz, 1) = miXS_s_As42(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Am42(Ixy, Iz, 1) = miXS_s_Am42(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Am43(Ixy, Iz, 1) = miXS_s_Am43(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Am44(Ixy, Iz, 1) = miXS_s_Am44(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Cm42(Ixy, Iz, 1) = miXS_s_Cm42(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Cm43(Ixy, Iz, 1) = miXS_s_Cm43(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Cm44(Ixy, Iz, 1) = miXS_s_Cm44(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)

      DO Ig = 1, N_Group
         miXS_tr_I35 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Xe35(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Nd47(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Nd48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Nd49(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Pm47(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Ps48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Pm48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Pm49(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Sm47(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Sm48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Sm49(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_tr_I35 (Ixy, Iz, 1) = miXS_tr_I35 (Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Xe35(Ixy, Iz, 1) = miXS_tr_Xe35(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Nd47(Ixy, Iz, 1) = miXS_tr_Nd47(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Nd48(Ixy, Iz, 1) = miXS_tr_Nd48(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Nd49(Ixy, Iz, 1) = miXS_tr_Nd49(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Pm47(Ixy, Iz, 1) = miXS_tr_Pm47(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Ps48(Ixy, Iz, 1) = miXS_tr_Ps48(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Pm48(Ixy, Iz, 1) = miXS_tr_Pm48(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Pm49(Ixy, Iz, 1) = miXS_tr_Pm49(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Sm47(Ixy, Iz, 1) = miXS_tr_Sm47(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Sm48(Ixy, Iz, 1) = miXS_tr_Sm48(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Sm49(Ixy, Iz, 1) = miXS_tr_Sm49(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)

      DO Ig = 1, N_Group
         miXS_a_I35 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Xe35(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Nd47(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Nd48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Nd49(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Pm47(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Ps48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Pm48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Pm49(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Sm47(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Sm48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Sm49(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_a_I35 (Ixy, Iz, 1) = miXS_a_I35 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Xe35(Ixy, Iz, 1) = miXS_a_Xe35(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Nd47(Ixy, Iz, 1) = miXS_a_Nd47(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Nd48(Ixy, Iz, 1) = miXS_a_Nd48(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Nd49(Ixy, Iz, 1) = miXS_a_Nd49(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Pm47(Ixy, Iz, 1) = miXS_a_Pm47(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Ps48(Ixy, Iz, 1) = miXS_a_Ps48(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Pm48(Ixy, Iz, 1) = miXS_a_Pm48(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Pm49(Ixy, Iz, 1) = miXS_a_Pm49(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Sm47(Ixy, Iz, 1) = miXS_a_Sm47(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Sm48(Ixy, Iz, 1) = miXS_a_Sm48(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Sm49(Ixy, Iz, 1) = miXS_a_Sm49(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)

      DO Ig = 1, N_Group
         miXS_s_I35 (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Xe35(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Nd47(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Nd48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Nd49(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Pm47(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Ps48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Pm48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Pm49(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Sm47(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Sm48(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Sm49(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_s_I35 (Ixy, Iz, 1) = miXS_s_I35 (Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Xe35(Ixy, Iz, 1) = miXS_s_Xe35(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Nd47(Ixy, Iz, 1) = miXS_s_Nd47(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Nd48(Ixy, Iz, 1) = miXS_s_Nd48(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Nd49(Ixy, Iz, 1) = miXS_s_Nd49(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Pm47(Ixy, Iz, 1) = miXS_s_Pm47(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Ps48(Ixy, Iz, 1) = miXS_s_Ps48(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Pm48(Ixy, Iz, 1) = miXS_s_Pm48(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Pm49(Ixy, Iz, 1) = miXS_s_Pm49(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Sm47(Ixy, Iz, 1) = miXS_s_Sm47(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Sm48(Ixy, Iz, 1) = miXS_s_Sm48(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Sm49(Ixy, Iz, 1) = miXS_s_Sm49(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)

      DO Ig = 1, N_Group
         miXS_tr_Gd52(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Gd54(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Gd55(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Gd56(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Gd57(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Gd58(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_tr_Gd60(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_tr_Gd52(Ixy, Iz, 1) = miXS_tr_Gd52(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Gd54(Ixy, Iz, 1) = miXS_tr_Gd54(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Gd55(Ixy, Iz, 1) = miXS_tr_Gd55(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Gd56(Ixy, Iz, 1) = miXS_tr_Gd56(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Gd57(Ixy, Iz, 1) = miXS_tr_Gd57(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Gd58(Ixy, Iz, 1) = miXS_tr_Gd58(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)
      miXS_tr_Gd60(Ixy, Iz, 1) = miXS_tr_Gd60(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)

      DO Ig = 1, N_Group
         miXS_a_Gd52(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Gd54(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Gd55(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Gd56(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Gd57(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Gd58(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_a_Gd60(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_a_Gd52(Ixy, Iz, 1) = miXS_a_Gd52(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Gd54(Ixy, Iz, 1) = miXS_a_Gd54(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Gd55(Ixy, Iz, 1) = miXS_a_Gd55(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Gd56(Ixy, Iz, 1) = miXS_a_Gd56(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Gd57(Ixy, Iz, 1) = miXS_a_Gd57(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Gd58(Ixy, Iz, 1) = miXS_a_Gd58(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)
      miXS_a_Gd60(Ixy, Iz, 1) = miXS_a_Gd60(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)

      DO Ig = 1, N_Group
         miXS_s_Gd52(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Gd54(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Gd55(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Gd56(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Gd57(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         miXS_s_Gd58(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_s_Gd52(Ixy, Iz, 1) = miXS_s_Gd52(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Gd54(Ixy, Iz, 1) = miXS_s_Gd54(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Gd55(Ixy, Iz, 1) = miXS_s_Gd55(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Gd56(Ixy, Iz, 1) = miXS_s_Gd56(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Gd57(Ixy, Iz, 1) = miXS_s_Gd57(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Gd58(Ixy, Iz, 1) = miXS_s_Gd58(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)
      miXS_s_Gd60(Ixy, Iz, 1) = miXS_s_Gd60(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)

      DO Ig = 1, N_Group
         miXS_tr_B0(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_tr_B0(Ixy, Iz, 1) = miXS_tr_B0(Ixy, Iz, 1) / leak_ratio(Ixy,Iz,1)

      DO Ig = 1, N_Group
         miXS_a_B0(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_a_B0(Ixy, Iz, 1) = miXS_a_B0(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,2)

      DO Ig = 1, N_Group
         miXS_s_B0(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      miXS_s_B0(Ixy, Iz, 1) = miXS_s_B0(Ixy, Iz, 1) * leak_ratio(Ixy,Iz,3)

      miXS_n2n_U35 (Ixy, Iz, 1) = XSset(Ixs);  Ixs = Ixs + 1
      miXS_n2n_U38 (Ixy, Iz, 1) = XSset(Ixs);  Ixs = Ixs + 1
      miXS_n2n_Np37(Ixy, Iz, 1) = XSset(Ixs);  Ixs = Ixs + 1
      miXS_n2n_Pu39(Ixy, Iz, 1) = XSset(Ixs);  Ixs = Ixs + 1

      Y_I35     (Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1
      Y_Xe35    (Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1
      Y_Nd47    (Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1
      Y_Nd48    (Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1
      Y_Nd49    (Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1
      Y_Pm47    (Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1
      Y_Ps48    (Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1
      Y_Pm48    (Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1
      Y_Pm49    (Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1
      Y_Sm49    (Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1
      Y_Xe35_eff(Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1
      Y_Pm49_eff(Ixy, Iz) = XSset(Ixs);  Ixs = Ixs + 1

      DO Ig = 1, N_Group
         ADF_Lx  (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         ADF_Rx  (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         ADF_Ly  (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         ADF_Ry  (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         ADF_Avg (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1

         CDF_LxLy(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         CDF_LxRy(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         CDF_RxLy(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         CDF_RxRy(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO
      DO Ig = 1, N_Group
         v_Inv(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO

      DO Ig = 1, N_Group_d
         beta_d_eff(Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
         lambda_d  (Ixy, Iz, Ig) = XSset(Ixs);  Ixs = Ixs + 1
      END DO

      RETURN
      END SUBROUTINE miXSFB_LC

      SUBROUTINE maXSFB_LC(Ixy, Iz)

      USE Inc_3D , ONLY: Flux
      USE Inc_INP , ONLY: Flag_Card_maXS
      USE Inc_maXS
      USE Inc_miXS
      USE Inc_Nuclide
      USE Inc_Option, ONLY: N_Group
      USE Inc_XS_File, ONLY: XSset

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Ixy, Iz
      INTEGER :: Ig


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [maXSFB_LC] in Mod_LeakCorr'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF ( Flag_Card_maXS .EQV. .TRUE. ) THEN
         maXS_tr_3D(Ixy, Iz, 1) = XSset(1)
         maXS_tr_3D(Ixy, Iz, 2) = XSset(2)

         D_3D(Ixy, Iz, 1) = 1d0 / max(3.0d0*maXS_tr_3D(Ixy, Iz, 1), 1d-50)
         D_3D(Ixy, Iz, 2) = 1d0 / max(3.0d0*maXS_tr_3D(Ixy, Iz, 2), 1d-50)

         maXS_a_3D(Ixy, Iz, 1) = XSset(3)
         maXS_a_3D(Ixy, Iz, 2) = XSset(4)

         maXS_f_3D(Ixy, Iz, 1) = XSset(5)
         maXS_f_3D(Ixy, Iz, 2) = XSset(6)

         nu_maXS_f_3D(Ixy, Iz, 1) = XSset(7)
         nu_maXS_f_3D(Ixy, Iz, 2) = XSset(8)

         kap_maXS_f_3D(Ixy, Iz, 1) = XSset(9)
         kap_maXS_f_3D(Ixy, Iz, 2) = XSset(10)

         maXS_s_3D(Ixy, Iz, 1) = XSset(11)
         maXS_s_3D(Ixy, Iz, 2) = XSset(12)

         maXS_r_3D(Ixy, Iz, 1) = maXS_a_3D(Ixy, Iz, 1) + maXS_s_3D(Ixy, Iz, 1)

         maXS_r_3D(Ixy, Iz, 2) = maXS_a_3D(Ixy, Iz, 2) + maXS_s_3D(Ixy, Iz, 2)

      ELSE
         DO Ig = 1, N_Group
            maXS_tr_3D(Ixy, Iz, Ig) = XSset(Ig)         &
               + N_U34 (Ixy, Iz) * miXS_tr_U34 (Ixy, Iz, Ig)  &
               + N_U35 (Ixy, Iz) * miXS_tr_U35 (Ixy, Iz, Ig)  &
               + N_U36 (Ixy, Iz) * miXS_tr_U36 (Ixy, Iz, Ig)  &
               + N_U37 (Ixy, Iz) * miXS_tr_U37 (Ixy, Iz, Ig)  &
               + N_U38 (Ixy, Iz) * miXS_tr_U38 (Ixy, Iz, Ig)  &
               + N_Np37(Ixy, Iz) * miXS_tr_Np37(Ixy, Iz, Ig)  &
               + N_Np38(Ixy, Iz) * miXS_tr_Np38(Ixy, Iz, Ig)  &
               + N_Np39(Ixy, Iz) * miXS_tr_Np39(Ixy, Iz, Ig)  &
               + N_Pu38(Ixy, Iz) * miXS_tr_Pu38(Ixy, Iz, Ig)  &
               + N_Pu39(Ixy, Iz) * miXS_tr_Pu39(Ixy, Iz, Ig)  &
               + N_Pu40(Ixy, Iz) * miXS_tr_Pu40(Ixy, Iz, Ig)  &
               + N_Pu41(Ixy, Iz) * miXS_tr_Pu41(Ixy, Iz, Ig)  &
               + N_Pu42(Ixy, Iz) * miXS_tr_Pu42(Ixy, Iz, Ig)  &
               + N_Pu43(Ixy, Iz) * miXS_tr_Pu43(Ixy, Iz, Ig)  &
               + N_Am41(Ixy, Iz) * miXS_tr_Am41(Ixy, Iz, Ig)  &
               + N_As42(Ixy, Iz) * miXS_tr_As42(Ixy, Iz, Ig)  &
               + N_Am42(Ixy, Iz) * miXS_tr_Am42(Ixy, Iz, Ig)  &
               + N_Am43(Ixy, Iz) * miXS_tr_Am43(Ixy, Iz, Ig)  &
               + N_Am44(Ixy, Iz) * miXS_tr_Am44(Ixy, Iz, Ig)  &
               + N_Cm42(Ixy, Iz) * miXS_tr_Cm42(Ixy, Iz, Ig)  &
               + N_Cm43(Ixy, Iz) * miXS_tr_Cm43(Ixy, Iz, Ig)  &
               + N_Cm44(Ixy, Iz) * miXS_tr_Cm44(Ixy, Iz, Ig)  &

               + N_I35 (Ixy, Iz) * miXS_tr_I35 (Ixy, Iz, Ig)  &
               + N_Xe35(Ixy, Iz) * miXS_tr_Xe35(Ixy, Iz, Ig)  &
               + N_Nd47(Ixy, Iz) * miXS_tr_Nd47(Ixy, Iz, Ig)  &
               + N_Nd48(Ixy, Iz) * miXS_tr_Nd48(Ixy, Iz, Ig)  &
               + N_Nd49(Ixy, Iz) * miXS_tr_Nd49(Ixy, Iz, Ig)  &
               + N_Pm47(Ixy, Iz) * miXS_tr_Pm47(Ixy, Iz, Ig)  &
               + N_Ps48(Ixy, Iz) * miXS_tr_Ps48(Ixy, Iz, Ig)  &
               + N_Pm48(Ixy, Iz) * miXS_tr_Pm48(Ixy, Iz, Ig)  &
               + N_Pm49(Ixy, Iz) * miXS_tr_Pm49(Ixy, Iz, Ig)  &
               + N_Sm47(Ixy, Iz) * miXS_tr_Sm47(Ixy, Iz, Ig)  &
               + N_Sm48(Ixy, Iz) * miXS_tr_Sm48(Ixy, Iz, Ig)  &
               + N_Sm49(Ixy, Iz) * miXS_tr_Sm49(Ixy, Iz, Ig)  &

               + N_Gd52(Ixy, Iz) * miXS_tr_Gd52(Ixy, Iz, Ig)  &
               + N_Gd54(Ixy, Iz) * miXS_tr_Gd54(Ixy, Iz, Ig)  &
               + N_Gd55(Ixy, Iz) * miXS_tr_Gd55(Ixy, Iz, Ig)  &
               + N_Gd56(Ixy, Iz) * miXS_tr_Gd56(Ixy, Iz, Ig)  &
               + N_Gd57(Ixy, Iz) * miXS_tr_Gd57(Ixy, Iz, Ig)  &
               + N_Gd58(Ixy, Iz) * miXS_tr_Gd58(Ixy, Iz, Ig)  &
               + N_Gd60(Ixy, Iz) * miXS_tr_Gd60(Ixy, Iz, Ig)  &

               + N_B0  (Ixy, Iz) * miXS_tr_B0(Ixy, Iz, Ig)

            D_3D(Ixy, Iz, Ig) = 1d0 / max(3.0d0*maXS_tr_3D(Ixy, Iz, Ig), 1d-50)

            maXS_a_3D(Ixy, Iz, Ig) = XSset(2 + Ig)     &
               + N_U34 (Ixy, Iz) * miXS_a_U34 (Ixy, Iz, Ig)  &
               + N_U35 (Ixy, Iz) * miXS_a_U35 (Ixy, Iz, Ig)  &
               + N_U36 (Ixy, Iz) * miXS_a_U36 (Ixy, Iz, Ig)  &
               + N_U37 (Ixy, Iz) * miXS_a_U37 (Ixy, Iz, Ig)  &
               + N_U38 (Ixy, Iz) * miXS_a_U38 (Ixy, Iz, Ig)  &
               + N_Np37(Ixy, Iz) * miXS_a_Np37(Ixy, Iz, Ig)  &
               + N_Np38(Ixy, Iz) * miXS_a_Np38(Ixy, Iz, Ig)  &
               + N_Np39(Ixy, Iz) * miXS_a_Np39(Ixy, Iz, Ig)  &
               + N_Pu38(Ixy, Iz) * miXS_a_Pu38(Ixy, Iz, Ig)  &
               + N_Pu39(Ixy, Iz) * miXS_a_Pu39(Ixy, Iz, Ig)  &
               + N_Pu40(Ixy, Iz) * miXS_a_Pu40(Ixy, Iz, Ig)  &
               + N_Pu41(Ixy, Iz) * miXS_a_Pu41(Ixy, Iz, Ig)  &
               + N_Pu42(Ixy, Iz) * miXS_a_Pu42(Ixy, Iz, Ig)  &
               + N_Pu43(Ixy, Iz) * miXS_a_Pu43(Ixy, Iz, Ig)  &
               + N_Am41(Ixy, Iz) * miXS_a_Am41(Ixy, Iz, Ig)  &
               + N_As42(Ixy, Iz) * miXS_a_As42(Ixy, Iz, Ig)  &
               + N_Am42(Ixy, Iz) * miXS_a_Am42(Ixy, Iz, Ig)  &
               + N_Am43(Ixy, Iz) * miXS_a_Am43(Ixy, Iz, Ig)  &
               + N_Am44(Ixy, Iz) * miXS_a_Am44(Ixy, Iz, Ig)  &
               + N_Cm42(Ixy, Iz) * miXS_a_Cm42(Ixy, Iz, Ig)  &
               + N_Cm43(Ixy, Iz) * miXS_a_Cm43(Ixy, Iz, Ig)  &
               + N_Cm44(Ixy, Iz) * miXS_a_Cm44(Ixy, Iz, Ig)  &

               + N_I35 (Ixy, Iz) * miXS_a_I35 (Ixy, Iz, Ig)  &
               + N_Xe35(Ixy, Iz) * miXS_a_Xe35(Ixy, Iz, Ig)  &
               + N_Nd47(Ixy, Iz) * miXS_a_Nd47(Ixy, Iz, Ig)  &
               + N_Nd48(Ixy, Iz) * miXS_a_Nd48(Ixy, Iz, Ig)  &
               + N_Nd49(Ixy, Iz) * miXS_a_Nd49(Ixy, Iz, Ig)  &
               + N_Pm47(Ixy, Iz) * miXS_a_Pm47(Ixy, Iz, Ig)  &
               + N_Ps48(Ixy, Iz) * miXS_a_Ps48(Ixy, Iz, Ig)  &
               + N_Pm48(Ixy, Iz) * miXS_a_Pm48(Ixy, Iz, Ig)  &
               + N_Pm49(Ixy, Iz) * miXS_a_Pm49(Ixy, Iz, Ig)  &
               + N_Sm47(Ixy, Iz) * miXS_a_Sm47(Ixy, Iz, Ig)  &
               + N_Sm48(Ixy, Iz) * miXS_a_Sm48(Ixy, Iz, Ig)  &
               + N_Sm49(Ixy, Iz) * miXS_a_Sm49(Ixy, Iz, Ig)  &

               + N_Gd52(Ixy, Iz) * miXS_a_Gd52(Ixy, Iz, Ig)  &
               + N_Gd54(Ixy, Iz) * miXS_a_Gd54(Ixy, Iz, Ig)  &
               + N_Gd55(Ixy, Iz) * miXS_a_Gd55(Ixy, Iz, Ig)  &
               + N_Gd56(Ixy, Iz) * miXS_a_Gd56(Ixy, Iz, Ig)  &
               + N_Gd57(Ixy, Iz) * miXS_a_Gd57(Ixy, Iz, Ig)  &
               + N_Gd58(Ixy, Iz) * miXS_a_Gd58(Ixy, Iz, Ig)  &
               + N_Gd60(Ixy, Iz) * miXS_a_Gd60(Ixy, Iz, Ig)  &

               + N_B0  (Ixy, Iz) * miXS_a_B0(Ixy, Iz, Ig)

            maXS_f_3D(Ixy, Iz, Ig) = XSset(4 + Ig)     &
               + N_U34 (Ixy, Iz) * miXS_f_U34 (Ixy, Iz, Ig)  &
               + N_U35 (Ixy, Iz) * miXS_f_U35 (Ixy, Iz, Ig)  &
               + N_U36 (Ixy, Iz) * miXS_f_U36 (Ixy, Iz, Ig)  &
               + N_U37 (Ixy, Iz) * miXS_f_U37 (Ixy, Iz, Ig)  &
               + N_U38 (Ixy, Iz) * miXS_f_U38 (Ixy, Iz, Ig)  &
               + N_Np37(Ixy, Iz) * miXS_f_Np37(Ixy, Iz, Ig)  &
               + N_Np38(Ixy, Iz) * miXS_f_Np38(Ixy, Iz, Ig)  &
               + N_Np39(Ixy, Iz) * miXS_f_Np39(Ixy, Iz, Ig)  &
               + N_Pu38(Ixy, Iz) * miXS_f_Pu38(Ixy, Iz, Ig)  &
               + N_Pu39(Ixy, Iz) * miXS_f_Pu39(Ixy, Iz, Ig)  &
               + N_Pu40(Ixy, Iz) * miXS_f_Pu40(Ixy, Iz, Ig)  &
               + N_Pu41(Ixy, Iz) * miXS_f_Pu41(Ixy, Iz, Ig)  &
               + N_Pu42(Ixy, Iz) * miXS_f_Pu42(Ixy, Iz, Ig)  &
               + N_Pu43(Ixy, Iz) * miXS_f_Pu43(Ixy, Iz, Ig)  &
               + N_Am41(Ixy, Iz) * miXS_f_Am41(Ixy, Iz, Ig)  &
               + N_As42(Ixy, Iz) * miXS_f_As42(Ixy, Iz, Ig)  &
               + N_Am42(Ixy, Iz) * miXS_f_Am42(Ixy, Iz, Ig)  &
               + N_Am43(Ixy, Iz) * miXS_f_Am43(Ixy, Iz, Ig)  &
               + N_Am44(Ixy, Iz) * miXS_f_Am44(Ixy, Iz, Ig)  &
               + N_Cm42(Ixy, Iz) * miXS_f_Cm42(Ixy, Iz, Ig)  &
               + N_Cm43(Ixy, Iz) * miXS_f_Cm43(Ixy, Iz, Ig)  &
               + N_Cm44(Ixy, Iz) * miXS_f_Cm44(Ixy, Iz, Ig)

            nu_maXS_f_3D(Ixy, Iz, Ig) = XSset(6 + Ig)     &
               + N_U34 (Ixy, Iz) * nu_miXS_f_U34 (Ixy, Iz, Ig)  &
               + N_U35 (Ixy, Iz) * nu_miXS_f_U35 (Ixy, Iz, Ig)  &
               + N_U36 (Ixy, Iz) * nu_miXS_f_U36 (Ixy, Iz, Ig)  &
               + N_U37 (Ixy, Iz) * nu_miXS_f_U37 (Ixy, Iz, Ig)  &
               + N_U38 (Ixy, Iz) * nu_miXS_f_U38 (Ixy, Iz, Ig)  &
               + N_Np37(Ixy, Iz) * nu_miXS_f_Np37(Ixy, Iz, Ig)  &
               + N_Np38(Ixy, Iz) * nu_miXS_f_Np38(Ixy, Iz, Ig)  &
               + N_Np39(Ixy, Iz) * nu_miXS_f_Np39(Ixy, Iz, Ig)  &
               + N_Pu38(Ixy, Iz) * nu_miXS_f_Pu38(Ixy, Iz, Ig)  &
               + N_Pu39(Ixy, Iz) * nu_miXS_f_Pu39(Ixy, Iz, Ig)  &
               + N_Pu40(Ixy, Iz) * nu_miXS_f_Pu40(Ixy, Iz, Ig)  &
               + N_Pu41(Ixy, Iz) * nu_miXS_f_Pu41(Ixy, Iz, Ig)  &
               + N_Pu42(Ixy, Iz) * nu_miXS_f_Pu42(Ixy, Iz, Ig)  &
               + N_Pu43(Ixy, Iz) * nu_miXS_f_Pu43(Ixy, Iz, Ig)  &
               + N_Am41(Ixy, Iz) * nu_miXS_f_Am41(Ixy, Iz, Ig)  &
               + N_As42(Ixy, Iz) * nu_miXS_f_As42(Ixy, Iz, Ig)  &
               + N_Am42(Ixy, Iz) * nu_miXS_f_Am42(Ixy, Iz, Ig)  &
               + N_Am43(Ixy, Iz) * nu_miXS_f_Am43(Ixy, Iz, Ig)  &
               + N_Am44(Ixy, Iz) * nu_miXS_f_Am44(Ixy, Iz, Ig)  &
               + N_Cm42(Ixy, Iz) * nu_miXS_f_Cm42(Ixy, Iz, Ig)  &
               + N_Cm43(Ixy, Iz) * nu_miXS_f_Cm43(Ixy, Iz, Ig)  &
               + N_Cm44(Ixy, Iz) * nu_miXS_f_Cm44(Ixy, Iz, Ig)

            kap_maXS_f_3D(Ixy, Iz, Ig) = XSset(8 + Ig)     &
               + N_U34 (Ixy, Iz) * kap_miXS_f_U34 (Ixy, Iz, Ig)  &
               + N_U35 (Ixy, Iz) * kap_miXS_f_U35 (Ixy, Iz, Ig)  &
               + N_U36 (Ixy, Iz) * kap_miXS_f_U36 (Ixy, Iz, Ig)  &
               + N_U37 (Ixy, Iz) * kap_miXS_f_U37 (Ixy, Iz, Ig)  &
               + N_U38 (Ixy, Iz) * kap_miXS_f_U38 (Ixy, Iz, Ig)  &
               + N_Np37(Ixy, Iz) * kap_miXS_f_Np37(Ixy, Iz, Ig)  &
               + N_Np38(Ixy, Iz) * kap_miXS_f_Np38(Ixy, Iz, Ig)  &
               + N_Np39(Ixy, Iz) * kap_miXS_f_Np39(Ixy, Iz, Ig)  &
               + N_Pu38(Ixy, Iz) * kap_miXS_f_Pu38(Ixy, Iz, Ig)  &
               + N_Pu39(Ixy, Iz) * kap_miXS_f_Pu39(Ixy, Iz, Ig)  &
               + N_Pu40(Ixy, Iz) * kap_miXS_f_Pu40(Ixy, Iz, Ig)  &
               + N_Pu41(Ixy, Iz) * kap_miXS_f_Pu41(Ixy, Iz, Ig)  &
               + N_Pu42(Ixy, Iz) * kap_miXS_f_Pu42(Ixy, Iz, Ig)  &
               + N_Pu43(Ixy, Iz) * kap_miXS_f_Pu43(Ixy, Iz, Ig)  &
               + N_Am41(Ixy, Iz) * kap_miXS_f_Am41(Ixy, Iz, Ig)  &
               + N_As42(Ixy, Iz) * kap_miXS_f_As42(Ixy, Iz, Ig)  &
               + N_Am42(Ixy, Iz) * kap_miXS_f_Am42(Ixy, Iz, Ig)  &
               + N_Am43(Ixy, Iz) * kap_miXS_f_Am43(Ixy, Iz, Ig)  &
               + N_Am44(Ixy, Iz) * kap_miXS_f_Am44(Ixy, Iz, Ig)  &
               + N_Cm42(Ixy, Iz) * kap_miXS_f_Cm42(Ixy, Iz, Ig)  &
               + N_Cm43(Ixy, Iz) * kap_miXS_f_Cm43(Ixy, Iz, Ig)  &
               + N_Cm44(Ixy, Iz) * kap_miXS_f_Cm44(Ixy, Iz, Ig)

            maXS_s_3D(Ixy, Iz, Ig) = XSset(10 + Ig)    &
               + N_U34 (Ixy, Iz) * miXS_s_U34 (Ixy, Iz, Ig)  &
               + N_U35 (Ixy, Iz) * miXS_s_U35 (Ixy, Iz, Ig)  &
               + N_U36 (Ixy, Iz) * miXS_s_U36 (Ixy, Iz, Ig)  &
               + N_U37 (Ixy, Iz) * miXS_s_U37 (Ixy, Iz, Ig)  &
               + N_U38 (Ixy, Iz) * miXS_s_U38 (Ixy, Iz, Ig)  &
               + N_Np37(Ixy, Iz) * miXS_s_Np37(Ixy, Iz, Ig)  &
               + N_Np38(Ixy, Iz) * miXS_s_Np38(Ixy, Iz, Ig)  &
               + N_Np39(Ixy, Iz) * miXS_s_Np39(Ixy, Iz, Ig)  &
               + N_Pu38(Ixy, Iz) * miXS_s_Pu38(Ixy, Iz, Ig)  &
               + N_Pu39(Ixy, Iz) * miXS_s_Pu39(Ixy, Iz, Ig)  &
               + N_Pu40(Ixy, Iz) * miXS_s_Pu40(Ixy, Iz, Ig)  &
               + N_Pu41(Ixy, Iz) * miXS_s_Pu41(Ixy, Iz, Ig)  &
               + N_Pu42(Ixy, Iz) * miXS_s_Pu42(Ixy, Iz, Ig)  &
               + N_Pu43(Ixy, Iz) * miXS_s_Pu43(Ixy, Iz, Ig)  &
               + N_Am41(Ixy, Iz) * miXS_s_Am41(Ixy, Iz, Ig)  &
               + N_As42(Ixy, Iz) * miXS_s_As42(Ixy, Iz, Ig)  &
               + N_Am42(Ixy, Iz) * miXS_s_Am42(Ixy, Iz, Ig)  &
               + N_Am43(Ixy, Iz) * miXS_s_Am43(Ixy, Iz, Ig)  &
               + N_Am44(Ixy, Iz) * miXS_s_Am44(Ixy, Iz, Ig)  &
               + N_Cm42(Ixy, Iz) * miXS_s_Cm42(Ixy, Iz, Ig)  &
               + N_Cm43(Ixy, Iz) * miXS_s_Cm43(Ixy, Iz, Ig)  &
               + N_Cm44(Ixy, Iz) * miXS_s_Cm44(Ixy, Iz, Ig)  &

               + N_I35 (Ixy, Iz) * miXS_s_I35 (Ixy, Iz, Ig)  &
               + N_Xe35(Ixy, Iz) * miXS_s_Xe35(Ixy, Iz, Ig)  &
               + N_Nd47(Ixy, Iz) * miXS_s_Nd47(Ixy, Iz, Ig)  &
               + N_Nd48(Ixy, Iz) * miXS_s_Nd48(Ixy, Iz, Ig)  &
               + N_Nd49(Ixy, Iz) * miXS_s_Nd49(Ixy, Iz, Ig)  &
               + N_Pm47(Ixy, Iz) * miXS_s_Pm47(Ixy, Iz, Ig)  &
               + N_Ps48(Ixy, Iz) * miXS_s_Ps48(Ixy, Iz, Ig)  &
               + N_Pm48(Ixy, Iz) * miXS_s_Pm48(Ixy, Iz, Ig)  &
               + N_Pm49(Ixy, Iz) * miXS_s_Pm49(Ixy, Iz, Ig)  &
               + N_Sm47(Ixy, Iz) * miXS_s_Sm47(Ixy, Iz, Ig)  &
               + N_Sm48(Ixy, Iz) * miXS_s_Sm48(Ixy, Iz, Ig)  &
               + N_Sm49(Ixy, Iz) * miXS_s_Sm49(Ixy, Iz, Ig)  &

               + N_Gd52(Ixy, Iz) * miXS_s_Gd52(Ixy, Iz, Ig)  &
               + N_Gd54(Ixy, Iz) * miXS_s_Gd54(Ixy, Iz, Ig)  &
               + N_Gd55(Ixy, Iz) * miXS_s_Gd55(Ixy, Iz, Ig)  &
               + N_Gd56(Ixy, Iz) * miXS_s_Gd56(Ixy, Iz, Ig)  &
               + N_Gd57(Ixy, Iz) * miXS_s_Gd57(Ixy, Iz, Ig)  &
               + N_Gd58(Ixy, Iz) * miXS_s_Gd58(Ixy, Iz, Ig)  &
               + N_Gd60(Ixy, Iz) * miXS_s_Gd60(Ixy, Iz, Ig)  &

               + N_B0  (Ixy, Iz) * miXS_s_B0(Ixy, Iz, Ig)
         END DO

         IF ( XSset(12) /= D0 ) THEN
            maXS_s_3D(Ixy, Iz, 1) = maXS_s_3D(Ixy, Iz, 1)  &
                                  - ( Flux(Ixy, Iz, 2) / Flux(Ixy, Iz, 1) ) * maXS_s_3D(Ixy, Iz, 2)
         END IF

         maXS_s_3D(Ixy, Iz, 2) = D0

         DO Ig = 1, N_Group
            maXS_r_3D(Ixy, Iz, Ig) = maXS_a_3D(Ixy, Iz, Ig) + maXS_s_3D(Ixy, Iz, Ig)
         END DO

      END IF


      maXS_scat_3D(1, 2, Ixy, Iz) = maXS_s_3D(Ixy, Iz, 1)
      RETURN
      END SUBROUTINE maXSFB_LC

      END MODULE Mod_LeakCorr


#endif
