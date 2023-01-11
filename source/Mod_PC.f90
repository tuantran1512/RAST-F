#ifdef siarhei_delete


      MODULE Mod_PC
      ! [Predictor/Corrector Method]
      ! Previous Burnup Step = n     (Calculated)
      ! Current  Burnup Step = n + 1 (To be Calculated)
      ! maXS_f(n + 1/2) = PC_w*maXS_f_c(n) + (1 - PC_w)*maXS_f_p(n + 1)
      ! Y(n + 1/2)      = PC_w*Y_c(n)      + (1 - PC_w)*Y_p(n + 1)
      ! miXS(n + 1/2)   = PC_w*miXS_c(n)   + (1 - PC_w)*miXS_p(n + 1)
      ! Flux(n + 1/2)   = PC_w*Flux_c(n)   + (1 - PC_w)*Flux_p(n + 1)
      !
      ! maXS_f and Y(Fission Yield) are need in Fission Product Micro Depletion.
      ! miXS is only for Specific Nuclides (Not Residual)
      ! and only for Absorption, Capture, Fission, n2n XS are used.
      !
      ! 1. Get N_p(n + 1) using maXS_f_c(n), Y_c(n), miXS_c(n), Flux_c(n), dT, N_c(n)
      ! 2. maXS(n + 1) Feedback
      ! 3. Get maXS_f_p(n + 1), Y_p(n + 1), miXS_p(n + 1), Flux_p(n + 1)
      !    after Boron Search Iteration
      ! 4. Get maXS_f(n + 1/2), Y(n + 1/2), miXS(n + 1/2), Flux(n + 1/2) using PC_w
      ! 5. Get N_c(n + 1) using maXS_f(n + 1/2), Y(n + 1/2), miXS(n + 1/2), Flux(n + 1/2), dT, N_c(n)
      ! 6. maXS(n + 1) Feedback
      ! 7. Get maXS_f_c(n + 1), Y_c(n + 1), miXS_c(n + 1), Flux_c(n + 1)
      !    after Boron Search Iteration
      ! 8. Go to Next Burnup Step
      !
      ! Step 1 ~ 3 is Predictor Step
      ! Step 4 ~ 5 is Semi Corrector Step
      ! Step 6 ~ 7 is Full Corrector Step

      USE Inc_Constant
      USE Inc_Geometry
      USE Inc_RP
      use Mod_GetNode, only: new_asym_itab

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE PC_Weight_RR
      USE Inc_3D, ONLY: Flux, Flux_Old
      USE Inc_Depletion, ONLY: PC_w
      USE Inc_maXS, ONLY: maXS_f_3D, maXS_f_3D_Old
      USE Inc_miXS
      USE Inc_Option, ONLY: N_Group
      USE Inc_XS_File, ONLY: I_Tab


#ifdef tuan_fr
      USE Inc_TPEN, ONLY: fluxf,fluxf_old
#endif
      IMPLICIT NONE
      INTEGER :: Ixy, Ixy_FA, Iz, Ig

      ! maXS_f Weight
      maXS_f_3D = PC_w*maXS_f_3D_Old + (D1 - PC_w)*maXS_f_3D

      ! Y Weight
      Y_I35      = PC_w*Y_I35_Old      + (D1 - PC_w)*Y_I35
      Y_Xe35     = PC_w*Y_Xe35_Old     + (D1 - PC_w)*Y_Xe35
      Y_Nd47     = PC_w*Y_Nd47_Old     + (D1 - PC_w)*Y_Nd47
      Y_Nd48     = PC_w*Y_Nd48_Old     + (D1 - PC_w)*Y_Nd48
      Y_Nd49     = PC_w*Y_Nd49_Old     + (D1 - PC_w)*Y_Nd49
      Y_Pm47     = PC_w*Y_Pm47_Old     + (D1 - PC_w)*Y_Pm47
      Y_Ps48     = PC_w*Y_Ps48_Old     + (D1 - PC_w)*Y_Ps48
      Y_Pm48     = PC_w*Y_Pm48_Old     + (D1 - PC_w)*Y_Pm48
      Y_Pm49     = PC_w*Y_Pm49_Old     + (D1 - PC_w)*Y_Pm49
      Y_Sm49     = PC_w*Y_Sm49_Old     + (D1 - PC_w)*Y_Sm49
      Y_Xe35_eff = PC_w*Y_Xe35_eff_Old + (D1 - PC_w)*Y_Xe35_eff
      Y_Pm49_eff = PC_w*Y_Pm49_eff_Old + (D1 - PC_w)*Y_Pm49_eff

      ! Flux Weight
      DO Ixy = 1, Nxy
         DO Iz = 1, Nz
            DO Ig = 1, N_Group
               Flux(Ixy, Iz, Ig) = PC_w*Flux_Old(Ixy, Iz, Ig) + (D1 - PC_w)*Flux(Ixy, Iz, Ig)
#ifdef tuan_fr

               Fluxf(Ixy, Iz, Ig) = PC_w*Fluxf_Old(Ixy, Iz, Ig) + (D1 - PC_w)*Fluxf(Ixy, Iz, Ig)
#endif
            END DO
         END DO
      END DO

      ! Heavy Nuclide Absorption miXS
      miXS_a_U34  = PC_w*miXS_a_U34_Old  + (D1 - PC_w)*miXS_a_U34
      miXS_a_U35  = PC_w*miXS_a_U35_Old  + (D1 - PC_w)*miXS_a_U35
      miXS_a_U36  = PC_w*miXS_a_U36_Old  + (D1 - PC_w)*miXS_a_U36
      miXS_a_U37  = PC_w*miXS_a_U37_Old  + (D1 - PC_w)*miXS_a_U37
      miXS_a_U38  = PC_w*miXS_a_U38_Old  + (D1 - PC_w)*miXS_a_U38
      miXS_a_Np37 = PC_w*miXS_a_Np37_Old + (D1 - PC_w)*miXS_a_Np37
      miXS_a_Np38 = PC_w*miXS_a_Np38_Old + (D1 - PC_w)*miXS_a_Np38
      miXS_a_Np39 = PC_w*miXS_a_Np39_Old + (D1 - PC_w)*miXS_a_Np39
      miXS_a_Pu38 = PC_w*miXS_a_Pu38_Old + (D1 - PC_w)*miXS_a_Pu38
      miXS_a_Pu39 = PC_w*miXS_a_Pu39_Old + (D1 - PC_w)*miXS_a_Pu39
      miXS_a_Pu40 = PC_w*miXS_a_Pu40_Old + (D1 - PC_w)*miXS_a_Pu40
      miXS_a_Pu41 = PC_w*miXS_a_Pu41_Old + (D1 - PC_w)*miXS_a_Pu41
      miXS_a_Pu42 = PC_w*miXS_a_Pu42_Old + (D1 - PC_w)*miXS_a_Pu42
      miXS_a_Pu43 = PC_w*miXS_a_Pu43_Old + (D1 - PC_w)*miXS_a_Pu43
      miXS_a_Am41 = PC_w*miXS_a_Am41_Old + (D1 - PC_w)*miXS_a_Am41
      miXS_a_As42 = PC_w*miXS_a_As42_Old + (D1 - PC_w)*miXS_a_As42
      miXS_a_Am42 = PC_w*miXS_a_Am42_Old + (D1 - PC_w)*miXS_a_Am42
      miXS_a_Am43 = PC_w*miXS_a_Am43_Old + (D1 - PC_w)*miXS_a_Am43
      miXS_a_Am44 = PC_w*miXS_a_Am44_Old + (D1 - PC_w)*miXS_a_Am44
      miXS_a_Cm42 = PC_w*miXS_a_Cm42_Old + (D1 - PC_w)*miXS_a_Cm42
      miXS_a_Cm43 = PC_w*miXS_a_Cm43_Old + (D1 - PC_w)*miXS_a_Cm43
      miXS_a_Cm44 = PC_w*miXS_a_Cm44_Old + (D1 - PC_w)*miXS_a_Cm44

      ! Heavy Nuclide Fission miXS
      miXS_f_U34  = PC_w*miXS_f_U34_Old  + (D1 - PC_w)*miXS_f_U34
      miXS_f_U35  = PC_w*miXS_f_U35_Old  + (D1 - PC_w)*miXS_f_U35
      miXS_f_U36  = PC_w*miXS_f_U36_Old  + (D1 - PC_w)*miXS_f_U36
      miXS_f_U37  = PC_w*miXS_f_U37_Old  + (D1 - PC_w)*miXS_f_U37
      miXS_f_U38  = PC_w*miXS_f_U38_Old  + (D1 - PC_w)*miXS_f_U38
      miXS_f_Np37 = PC_w*miXS_f_Np37_Old + (D1 - PC_w)*miXS_f_Np37
      miXS_f_Np38 = PC_w*miXS_f_Np38_Old + (D1 - PC_w)*miXS_f_Np38
      miXS_f_Np39 = PC_w*miXS_f_Np39_Old + (D1 - PC_w)*miXS_f_Np39
      miXS_f_Pu38 = PC_w*miXS_f_Pu38_Old + (D1 - PC_w)*miXS_f_Pu38
      miXS_f_Pu39 = PC_w*miXS_f_Pu39_Old + (D1 - PC_w)*miXS_f_Pu39
      miXS_f_Pu40 = PC_w*miXS_f_Pu40_Old + (D1 - PC_w)*miXS_f_Pu40
      miXS_f_Pu41 = PC_w*miXS_f_Pu41_Old + (D1 - PC_w)*miXS_f_Pu41
      miXS_f_Pu42 = PC_w*miXS_f_Pu42_Old + (D1 - PC_w)*miXS_f_Pu42
      miXS_f_Pu43 = PC_w*miXS_f_Pu43_Old + (D1 - PC_w)*miXS_f_Pu43
      miXS_f_Am41 = PC_w*miXS_f_Am41_Old + (D1 - PC_w)*miXS_f_Am41
      miXS_f_As42 = PC_w*miXS_f_As42_Old + (D1 - PC_w)*miXS_f_As42
      miXS_f_Am42 = PC_w*miXS_f_Am42_Old + (D1 - PC_w)*miXS_f_Am42
      miXS_f_Am43 = PC_w*miXS_f_Am43_Old + (D1 - PC_w)*miXS_f_Am43
      miXS_f_Am44 = PC_w*miXS_f_Am44_Old + (D1 - PC_w)*miXS_f_Am44
      miXS_f_Cm42 = PC_w*miXS_f_Cm42_Old + (D1 - PC_w)*miXS_f_Cm42
      miXS_f_Cm43 = PC_w*miXS_f_Cm43_Old + (D1 - PC_w)*miXS_f_Cm43
      miXS_f_Cm44 = PC_w*miXS_f_Cm44_Old + (D1 - PC_w)*miXS_f_Cm44

      ! Fission Product Absorption miXS
      miXS_a_I35  = PC_w*miXS_a_I35_Old  + (D1 - PC_w)*miXS_a_I35
      miXS_a_Xe35 = PC_w*miXS_a_Xe35_Old + (D1 - PC_w)*miXS_a_Xe35
      miXS_a_Nd47 = PC_w*miXS_a_Nd47_Old + (D1 - PC_w)*miXS_a_Nd47
      miXS_a_Nd48 = PC_w*miXS_a_Nd48_Old + (D1 - PC_w)*miXS_a_Nd48
      miXS_a_Nd49 = PC_w*miXS_a_Nd49_Old + (D1 - PC_w)*miXS_a_Nd49
      miXS_a_Pm47 = PC_w*miXS_a_Pm47_Old + (D1 - PC_w)*miXS_a_Pm47
      miXS_a_Ps48 = PC_w*miXS_a_Ps48_Old + (D1 - PC_w)*miXS_a_Ps48
      miXS_a_Pm48 = PC_w*miXS_a_Pm48_Old + (D1 - PC_w)*miXS_a_Pm48
      miXS_a_Pm49 = PC_w*miXS_a_Pm49_Old + (D1 - PC_w)*miXS_a_Pm49
      miXS_a_Sm47 = PC_w*miXS_a_Sm47_Old + (D1 - PC_w)*miXS_a_Sm47
      miXS_a_Sm48 = PC_w*miXS_a_Sm48_Old + (D1 - PC_w)*miXS_a_Sm48
      miXS_a_Sm49 = PC_w*miXS_a_Sm49_Old + (D1 - PC_w)*miXS_a_Sm49

      ! Burnable Absorber Absorption miXS
      DO Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)

         DO Iz = IzFuelBot, IzFuelTop
            I_Tab = AxialComp( I_LP_1N( I_4Nto1N(Ixy) ), Iz )
            I_Tab = new_asym_itab(I_Tab,Ixy)
            IF ( Flag_BP(I_Tab) .EQV. .FALSE. ) CYCLE

            DO Ig = 1, N_Group
               miXS_a_Gd52(Ixy, Iz, Ig) = PC_w*miXS_a_Gd52_Old(Ixy, Iz, Ig) + (D1 - PC_w)*miXS_a_Gd52(Ixy, Iz, Ig)
               miXS_a_Gd54(Ixy, Iz, Ig) = PC_w*miXS_a_Gd54_Old(Ixy, Iz, Ig) + (D1 - PC_w)*miXS_a_Gd54(Ixy, Iz, Ig)
               miXS_a_Gd55(Ixy, Iz, Ig) = PC_w*miXS_a_Gd55_Old(Ixy, Iz, Ig) + (D1 - PC_w)*miXS_a_Gd55(Ixy, Iz, Ig)
               miXS_a_Gd56(Ixy, Iz, Ig) = PC_w*miXS_a_Gd56_Old(Ixy, Iz, Ig) + (D1 - PC_w)*miXS_a_Gd56(Ixy, Iz, Ig)
               miXS_a_Gd57(Ixy, Iz, Ig) = PC_w*miXS_a_Gd57_Old(Ixy, Iz, Ig) + (D1 - PC_w)*miXS_a_Gd57(Ixy, Iz, Ig)
               miXS_a_Gd58(Ixy, Iz, Ig) = PC_w*miXS_a_Gd58_Old(Ixy, Iz, Ig) + (D1 - PC_w)*miXS_a_Gd58(Ixy, Iz, Ig)
               miXS_a_Gd60(Ixy, Iz, Ig) = PC_w*miXS_a_Gd60_Old(Ixy, Iz, Ig) + (D1 - PC_w)*miXS_a_Gd60(Ixy, Iz, Ig)
            END DO
         END DO
      END DO

      ! Moderator Absorption miXS
      miXS_a_H2O = PC_w*miXS_a_H2O_Old + (D1 - PC_w)*miXS_a_H2O
      miXS_a_B0  = PC_w*miXS_a_B0_Old  + (D1 - PC_w)*miXS_a_B0

      ! (n, 2n) Reaction
      miXS_n2n_U35  = PC_w*miXS_n2n_U35_Old  + (D1 - PC_w)*miXS_n2n_U35
      miXS_n2n_U38  = PC_w*miXS_n2n_U38_Old  + (D1 - PC_w)*miXS_n2n_U38
      miXS_n2n_Np37 = PC_w*miXS_n2n_Np37_Old + (D1 - PC_w)*miXS_n2n_Np37
      miXS_n2n_Pu39 = PC_w*miXS_n2n_Pu39_Old + (D1 - PC_w)*miXS_n2n_Pu39

      miXS_a_B10 = PC_w*miXS_a_B10_Old + (D1 - PC_w)*miXS_a_B10

      RETURN
      END SUBROUTINE PC_Weight_RR


      subroutine PC_Weight_ND(opt)
      use Inc_Depletion, only: PC_w
      use inc_nuclide
      use mod_alloc
#ifdef tuan_fr
      use Inc_miXS,  Only: N_FR, N_FR_Old, N_FR_Predictor
#endif
      implicit none
      logical(1), intent(in) :: opt


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [PC_Weight_ND] in Mod_PC'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.allocated(N_U34_Predictor)) then
         call alloc( N_U34_Predictor , Nxy, Nz )
         call alloc( N_U35_Predictor , Nxy, Nz )
         call alloc( N_U36_Predictor , Nxy, Nz )
         call alloc( N_U37_Predictor , Nxy, Nz )
         call alloc( N_U38_Predictor , Nxy, Nz )
         call alloc( N_Np37_Predictor, Nxy, Nz )
         call alloc( N_Np38_Predictor, Nxy, Nz )
         call alloc( N_Np39_Predictor, Nxy, Nz )
         call alloc( N_Pu38_Predictor, Nxy, Nz )
         call alloc( N_Pu39_Predictor, Nxy, Nz )
         call alloc( N_Pu40_Predictor, Nxy, Nz )
         call alloc( N_Pu41_Predictor, Nxy, Nz )
         call alloc( N_Pu42_Predictor, Nxy, Nz )
         call alloc( N_Pu43_Predictor, Nxy, Nz )
         call alloc( N_Am41_Predictor, Nxy, Nz )
         call alloc( N_As42_Predictor, Nxy, Nz )
         call alloc( N_Am42_Predictor, Nxy, Nz )
         call alloc( N_Am43_Predictor, Nxy, Nz )
         call alloc( N_Am44_Predictor, Nxy, Nz )
         call alloc( N_Cm42_Predictor, Nxy, Nz )
         call alloc( N_Cm43_Predictor, Nxy, Nz )
         call alloc( N_Cm44_Predictor, Nxy, Nz )
         call alloc( N_I35_Predictor , Nxy, Nz )
         call alloc( N_Xe35_Predictor, Nxy, Nz )
         call alloc( N_Nd47_Predictor, Nxy, Nz )
         call alloc( N_Nd48_Predictor, Nxy, Nz )
         call alloc( N_Nd49_Predictor, Nxy, Nz )
         call alloc( N_Pm47_Predictor, Nxy, Nz )
         call alloc( N_Ps48_Predictor, Nxy, Nz )
         call alloc( N_Pm48_Predictor, Nxy, Nz )
         call alloc( N_Pm49_Predictor, Nxy, Nz )
         call alloc( N_Sm47_Predictor, Nxy, Nz )
         call alloc( N_Sm48_Predictor, Nxy, Nz )
         call alloc( N_Sm49_Predictor, Nxy, Nz )
         call alloc( N_Gd52_Predictor, Nxy, Nz )
         call alloc( N_Gd54_Predictor, Nxy, Nz )
         call alloc( N_Gd55_Predictor, Nxy, Nz )
         call alloc( N_Gd56_Predictor, Nxy, Nz )
         call alloc( N_Gd57_Predictor, Nxy, Nz )
         call alloc( N_Gd58_Predictor, Nxy, Nz )
         call alloc( N_Gd60_Predictor, Nxy, Nz )

         N_U34_Predictor  = N_U34
         N_U35_Predictor  = N_U35
         N_U36_Predictor  = N_U36
         N_U37_Predictor  = N_U37
         N_U38_Predictor  = N_U38
         N_Np37_Predictor = N_Np37
         N_Np38_Predictor = N_Np38
         N_Np39_Predictor = N_Np39
         N_Pu38_Predictor = N_Pu38
         N_Pu39_Predictor = N_Pu39
         N_Pu40_Predictor = N_Pu40
         N_Pu41_Predictor = N_Pu41
         N_Pu42_Predictor = N_Pu42
         N_Pu43_Predictor = N_Pu43
         N_Am41_Predictor = N_Am41
         N_As42_Predictor = N_As42
         N_Am42_Predictor = N_Am42
         N_Am43_Predictor = N_Am43
         N_Am44_Predictor = N_Am44
         N_Cm42_Predictor = N_Cm42
         N_Cm43_Predictor = N_Cm43
         N_Cm44_Predictor = N_Cm44
         N_I35_Predictor  = N_I35
         N_Xe35_Predictor = N_Xe35
         N_Nd47_Predictor = N_Nd47
         N_Nd48_Predictor = N_Nd48
         N_Nd49_Predictor = N_Nd49
         N_Pm47_Predictor = N_Pm47
         N_Ps48_Predictor = N_Ps48
         N_Pm48_Predictor = N_Pm48
         N_Pm49_Predictor = N_Pm49
         N_Sm47_Predictor = N_Sm47
         N_Sm48_Predictor = N_Sm48
         N_Sm49_Predictor = N_Sm49
         N_Gd52_Predictor = N_Gd52
         N_Gd54_Predictor = N_Gd54
         N_Gd55_Predictor = N_Gd55
         N_Gd56_Predictor = N_Gd56
         N_Gd57_Predictor = N_Gd57
         N_Gd58_Predictor = N_Gd58
         N_Gd60_Predictor = N_Gd60
      endif
#ifdef tuan_fr
        N_FR_Predictor = N_FR

#endif

      if (.not.opt) then !PREDICTOR
         N_U34  = PC_w*N_U34_Predictor  + (1d0 - PC_w)*N_U34
         N_U35  = PC_w*N_U35_Predictor  + (1d0 - PC_w)*N_U35
         N_U36  = PC_w*N_U36_Predictor  + (1d0 - PC_w)*N_U36
         N_U37  = PC_w*N_U37_Predictor  + (1d0 - PC_w)*N_U37
         N_U38  = PC_w*N_U38_Predictor  + (1d0 - PC_w)*N_U38
         N_Np37 = PC_w*N_Np37_Predictor + (1d0 - PC_w)*N_Np37
         N_Np38 = PC_w*N_Np38_Predictor + (1d0 - PC_w)*N_Np38
         N_Np39 = PC_w*N_Np39_Predictor + (1d0 - PC_w)*N_Np39
         N_Pu38 = PC_w*N_Pu38_Predictor + (1d0 - PC_w)*N_Pu38
         N_Pu39 = PC_w*N_Pu39_Predictor + (1d0 - PC_w)*N_Pu39
         N_Pu40 = PC_w*N_Pu40_Predictor + (1d0 - PC_w)*N_Pu40
         N_Pu41 = PC_w*N_Pu41_Predictor + (1d0 - PC_w)*N_Pu41
         N_Pu42 = PC_w*N_Pu42_Predictor + (1d0 - PC_w)*N_Pu42
         N_Pu43 = PC_w*N_Pu43_Predictor + (1d0 - PC_w)*N_Pu43
         N_Am41 = PC_w*N_Am41_Predictor + (1d0 - PC_w)*N_Am41
         N_As42 = PC_w*N_As42_Predictor + (1d0 - PC_w)*N_As42
         N_Am42 = PC_w*N_Am42_Predictor + (1d0 - PC_w)*N_Am42
         N_Am43 = PC_w*N_Am43_Predictor + (1d0 - PC_w)*N_Am43
         N_Am44 = PC_w*N_Am44_Predictor + (1d0 - PC_w)*N_Am44
         N_Cm42 = PC_w*N_Cm42_Predictor + (1d0 - PC_w)*N_Cm42
         N_Cm43 = PC_w*N_Cm43_Predictor + (1d0 - PC_w)*N_Cm43
         N_Cm44 = PC_w*N_Cm44_Predictor + (1d0 - PC_w)*N_Cm44
         N_I35  = PC_w*N_I35_Predictor  + (1d0 - PC_w)*N_I35
         N_Xe35 = PC_w*N_Xe35_Predictor + (1d0 - PC_w)*N_Xe35
         N_Nd47 = PC_w*N_Nd47_Predictor + (1d0 - PC_w)*N_Nd47
         N_Nd48 = PC_w*N_Nd48_Predictor + (1d0 - PC_w)*N_Nd48
         N_Nd49 = PC_w*N_Nd49_Predictor + (1d0 - PC_w)*N_Nd49
         N_Pm47 = PC_w*N_Pm47_Predictor + (1d0 - PC_w)*N_Pm47
         N_Ps48 = PC_w*N_Ps48_Predictor + (1d0 - PC_w)*N_Ps48
         N_Pm48 = PC_w*N_Pm48_Predictor + (1d0 - PC_w)*N_Pm48
         N_Pm49 = PC_w*N_Pm49_Predictor + (1d0 - PC_w)*N_Pm49
         N_Sm47 = PC_w*N_Sm47_Predictor + (1d0 - PC_w)*N_Sm47
         N_Sm48 = PC_w*N_Sm48_Predictor + (1d0 - PC_w)*N_Sm48
         N_Sm49 = PC_w*N_Sm49_Predictor + (1d0 - PC_w)*N_Sm49
         N_Gd52 = PC_w*N_Gd52_Predictor + (1d0 - PC_w)*N_Gd52
         N_Gd54 = PC_w*N_Gd54_Predictor + (1d0 - PC_w)*N_Gd54
         N_Gd55 = PC_w*N_Gd55_Predictor + (1d0 - PC_w)*N_Gd55
         N_Gd56 = PC_w*N_Gd56_Predictor + (1d0 - PC_w)*N_Gd56
         N_Gd57 = PC_w*N_Gd57_Predictor + (1d0 - PC_w)*N_Gd57
         N_Gd58 = PC_w*N_Gd58_Predictor + (1d0 - PC_w)*N_Gd58
         N_Gd60 = PC_w*N_Gd60_Predictor + (1d0 - PC_w)*N_Gd60
#ifdef tuan_fr

         N_FR = PC_w*N_FR_Predictor + (1d0 - PC_w)*N_FR
#endif

      else !CORRECTOR

         N_U34_Predictor  = N_U34
         N_U35_Predictor  = N_U35
         N_U36_Predictor  = N_U36
         N_U37_Predictor  = N_U37
         N_U38_Predictor  = N_U38
         N_Np37_Predictor = N_Np37
         N_Np38_Predictor = N_Np38
         N_Np39_Predictor = N_Np39
         N_Pu38_Predictor = N_Pu38
         N_Pu39_Predictor = N_Pu39
         N_Pu40_Predictor = N_Pu40
         N_Pu41_Predictor = N_Pu41
         N_Pu42_Predictor = N_Pu42
         N_Pu43_Predictor = N_Pu43
         N_Am41_Predictor = N_Am41
         N_As42_Predictor = N_As42
         N_Am42_Predictor = N_Am42
         N_Am43_Predictor = N_Am43
         N_Am44_Predictor = N_Am44
         N_Cm42_Predictor = N_Cm42
         N_Cm43_Predictor = N_Cm43
         N_Cm44_Predictor = N_Cm44
         N_I35_Predictor  = N_I35
         N_Xe35_Predictor = N_Xe35
         N_Nd47_Predictor = N_Nd47
         N_Nd48_Predictor = N_Nd48
         N_Nd49_Predictor = N_Nd49
         N_Pm47_Predictor = N_Pm47
         N_Ps48_Predictor = N_Ps48
         N_Pm48_Predictor = N_Pm48
         N_Pm49_Predictor = N_Pm49
         N_Sm47_Predictor = N_Sm47
         N_Sm48_Predictor = N_Sm48
         N_Sm49_Predictor = N_Sm49
         N_Gd52_Predictor = N_Gd52
         N_Gd54_Predictor = N_Gd54
         N_Gd55_Predictor = N_Gd55
         N_Gd56_Predictor = N_Gd56
         N_Gd57_Predictor = N_Gd57
         N_Gd58_Predictor = N_Gd58
         N_Gd60_Predictor = N_Gd60
#ifdef tuan_fr
        N_FR_Predictor = N_FR
#endif

      endif

      RETURN
      END SUBROUTINE PC_Weight_ND


      subroutine PC_reset_ND
      use inc_nuclide
      use inc_3d, only: BU, BU_old
#ifdef tuan_fr
      use Inc_miXS, only: N_FR,N_FR_old
#endif

      implicit none


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [PC_reset_ND] in Mod_PC'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      BU = BU_Old

      N_U34  = N_U34_Old
      N_U35  = N_U35_Old
      N_U36  = N_U36_Old
      N_U37  = N_U37_Old
      N_U38  = N_U38_Old
      N_Np37 = N_Np37_Old
      N_Np38 = N_Np38_Old
      N_Np39 = N_Np39_Old
      N_Pu38 = N_Pu38_Old
      N_Pu39 = N_Pu39_Old
      N_Pu40 = N_Pu40_Old
      N_Pu41 = N_Pu41_Old
      N_Pu42 = N_Pu42_Old
      N_Pu43 = N_Pu43_Old
      N_Am41 = N_Am41_Old
      N_As42 = N_As42_Old
      N_Am42 = N_Am42_Old
      N_Am43 = N_Am43_Old
      N_Am44 = N_Am44_Old
      N_Cm42 = N_Cm42_Old
      N_Cm43 = N_Cm43_Old
      N_Cm44 = N_Cm44_Old
      N_I35  = N_I35_Old
      N_Xe35 = N_Xe35_Old
      N_Nd47 = N_Nd47_Old
      N_Nd48 = N_Nd48_Old
      N_Nd49 = N_Nd49_Old
      N_Pm47 = N_Pm47_Old
      N_Ps48 = N_Ps48_Old
      N_Pm48 = N_Pm48_Old
      N_Pm49 = N_Pm49_Old
      N_Sm47 = N_Sm47_Old
      N_Sm48 = N_Sm48_Old
      N_Sm49 = N_Sm49_Old
      N_Gd52 = N_Gd52_Old
      N_Gd54 = N_Gd54_Old
      N_Gd55 = N_Gd55_Old
      N_Gd56 = N_Gd56_Old
      N_Gd57 = N_Gd57_Old
      N_Gd58 = N_Gd58_Old
      N_Gd60 = N_Gd60_Old
#ifdef tuan_fr

      N_FR = N_FR_Old
#endif
      RETURN
      END SUBROUTINE PC_reset_ND

      END MODULE Mod_PC


#endif
