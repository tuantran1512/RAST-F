#ifdef siarhei_delete


      module Mod_Depletion
      use mod_charedit, only: print_msg
      use Inc_Constant
      use Inc_Geometry
      use Inc_RP
      use Mod_GetNode, only: new_asym_itab

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif


      use inc_parallel, only: comm

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none
      real(8) :: nd_cutoff=1d-10
      contains

      subroutine HeavyNuclide
      use Inc_3D, only: Flux
      use Inc_Depletion
      use Inc_miXS
      use Inc_Nuclide
      use Inc_Time, only: dT
!    !  use Mod_Operator, only: Get_InvC
      implicit none
      complex(8) :: Sum_CRAM_Vec(N_Heavy)
      real(8) :: A(N_Heavy, N_Heavy)
      real(8) :: Vec_N(N_Heavy)
      real(8) :: a_U34,  a_U35,  a_U36,  a_U37,  a_U38 , a_Np37, a_Np38, a_Np39, &
               & a_Pu38, a_Pu39, a_Pu40, a_Pu41, a_Pu42, a_Pu43, a_Am41, a_As42, &
               & a_Am42, a_Am43, a_Am44, a_Cm42, a_Cm43, a_Cm44
      real(8) :: f_U34,  f_U35,  f_U36,  f_U37,  f_U38 , f_Np37, f_Np38, f_Np39, &
               & f_Pu38, f_Pu39, f_Pu40, f_Pu41, f_Pu42, f_Pu43, f_Am41, f_As42, &
               & f_Am42, f_Am43, f_Am44, f_Cm42, f_Cm43, f_Cm44
      real(8) :: c_U34,  c_U35,  c_U36,  c_U37,  c_U38 , c_Np37, c_Np38, c_Np39, &
               & c_Pu38, c_Pu39, c_Pu40, c_Pu41, c_Pu42, c_Pu43, c_Am41, c_As42, &
               & c_Am42, c_Am43, c_Am44, c_Cm42, c_Cm43, c_Cm44
      real(8) :: r_U34,  r_U35,  r_U36,  r_U37,  r_U38 , r_Np37, r_Np38, r_Np39, &
               & r_Pu38, r_Pu39, r_Pu40, r_Pu41, r_Pu42, r_Pu43, r_Am41, r_As42, &
               & r_Am42, r_Am43, r_Am44, r_Cm42, r_Cm43, r_Cm44
      real(8) :: d_U34,  d_U35,  d_U36,  d_U37,  d_U38 , d_Np37, d_Np38, d_Np39, &
               & d_Pu38, d_Pu39, d_Pu40, d_Pu41, d_Pu42, d_Pu43, d_Am41, d_As42, &
               & d_Am42, d_Am43, d_Am44, d_Cm42, d_Cm43, d_Cm44
      real(8) :: n2n_U35, n2n_U38, n2n_Np37, n2n_Pu39
      integer :: Ixy, Iz, i

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [HeavyNuclide] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [HeavyNuclide] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.flag_midepl) return

      
#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      if ( abs(dT-0d0) < 1d-10 ) return

      ! decay term
      d_U34  = 8.98324d-14
      d_U35  = 3.12088d-17
      d_U36  = 9.38080d-16
      d_U37  = beta_U37
      d_U38  = 4.91594d-18
      d_Np37 = 1.02643d-14
      d_Np38 = beta_Np38
      d_Np39 = beta_Np39
      d_Pu38 = alpha_Pu38
      d_Pu39 = 9.12756d-13
      d_Pu40 = 3.35990d-12
      d_Pu41 = beta_Pu41
      d_Pu42 = 5.67688d-14
      d_Pu43 = beta_Pu43
      d_Am41 = 5.08172d-11
      d_As42 = 1.20162d-05 !TOdo:js+beta_As42
      d_Am42 = IT_Am42
      d_Am43 = 2.97616d-12
      d_Am44 = beta_Am44
      d_Cm42 = alpha_Cm42
      d_Cm43 = alpha_Cm43
      d_Cm44 = alpha_Cm44

      do Ixy = 1, Nxy
         do Iz = IzFuelBot, IzFuelTop
            if ( I_FARF_1N( I_4Nto1N(Ixy) ) == 1 ) cycle

            a_U34    = dot_product( miXS_a_U34   (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_U35    = dot_product( miXS_a_U35   (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_U36    = dot_product( miXS_a_U36   (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_U37    = dot_product( miXS_a_U37   (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_U38    = dot_product( miXS_a_U38   (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Np37   = dot_product( miXS_a_Np37  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Np38   = dot_product( miXS_a_Np38  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Np39   = dot_product( miXS_a_Np39  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Pu38   = dot_product( miXS_a_Pu38  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Pu39   = dot_product( miXS_a_Pu39  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Pu40   = dot_product( miXS_a_Pu40  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Pu41   = dot_product( miXS_a_Pu41  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Pu42   = dot_product( miXS_a_Pu42  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Pu43   = dot_product( miXS_a_Pu43  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Am41   = dot_product( miXS_a_Am41  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_As42   = dot_product( miXS_a_As42  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Am42   = dot_product( miXS_a_Am42  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Am43   = dot_product( miXS_a_Am43  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Am44   = dot_product( miXS_a_Am44  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Cm42   = dot_product( miXS_a_Cm42  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Cm43   = dot_product( miXS_a_Cm43  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            a_Cm44   = dot_product( miXS_a_Cm44  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_U34    = dot_product( miXS_f_U34   (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_U35    = dot_product( miXS_f_U35   (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_U36    = dot_product( miXS_f_U36   (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_U37    = dot_product( miXS_f_U37   (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_U38    = dot_product( miXS_f_U38   (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Np37   = dot_product( miXS_f_Np37  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Np38   = dot_product( miXS_f_Np38  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Np39   = dot_product( miXS_f_Np39  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Pu38   = dot_product( miXS_f_Pu38  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Pu39   = dot_product( miXS_f_Pu39  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Pu40   = dot_product( miXS_f_Pu40  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Pu41   = dot_product( miXS_f_Pu41  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Pu42   = dot_product( miXS_f_Pu42  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Pu43   = dot_product( miXS_f_Pu43  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Am41   = dot_product( miXS_f_Am41  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_As42   = dot_product( miXS_f_As42  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Am42   = dot_product( miXS_f_Am42  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Am43   = dot_product( miXS_f_Am43  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Am44   = dot_product( miXS_f_Am44  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Cm42   = dot_product( miXS_f_Cm42  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Cm43   = dot_product( miXS_f_Cm43  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            f_Cm44   = dot_product( miXS_f_Cm44  (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            n2n_U35  = dot_product( miXS_n2n_U35 (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            n2n_U38  = dot_product( miXS_n2n_U38 (Ixy,Iz,:), Flux(Ixy,Iz,:) )
            n2n_Np37 = dot_product( miXS_n2n_Np37(Ixy,Iz,:), Flux(Ixy,Iz,:) )
            n2n_Pu39 = dot_product( miXS_n2n_Pu39(Ixy,Iz,:), Flux(Ixy,Iz,:) )

            n2n_U35 = n2n_U38 * 0.71d0
            a_U35  = a_U35  + 2d0 * n2n_U35
            a_U38  = a_U38  + 2d0 * n2n_U38
            a_Np37 = a_Np37 + 2d0 * n2n_Np37
            a_Pu39 = a_Pu39 + 2d0 * n2n_Pu39

            c_U34  = a_U34  - f_U34
            c_U35  = a_U35  - f_U35  - n2n_U35
            c_U36  = a_U36  - f_U36
            c_U37  = a_U37  - f_U37
            c_U38  = a_U38  - f_U38  - n2n_U38
            c_Np37 = a_Np37 - f_Np37 - n2n_Np37
            c_Np38 = a_Np38 - f_Np38
            c_Np39 = a_Np39 - f_Np39
            c_Pu38 = a_Pu38 - f_Pu38
            c_Pu39 = a_Pu39 - f_Pu39 - n2n_Pu39
            c_Pu40 = a_Pu40 - f_Pu40
            c_Pu41 = a_Pu41 - f_Pu41
            c_Pu42 = a_Pu42 - f_Pu42
            c_Pu43 = a_Pu43 - f_Pu43
            c_Am41 = a_Am41 - f_Am41
            c_As42 = a_As42 - f_As42
            c_Am42 = a_Am42 - f_Am42
            c_Am43 = a_Am43 - f_Am43
            c_Am44 = a_Am44 - f_Am44
            c_Cm42 = a_Cm42 - f_Cm42
            c_Cm43 = a_Cm43 - f_Cm43
            c_Cm44 = a_Cm44 - f_Cm44

            r_U34  = a_U34  + d_U34
            r_U35  = a_U35  + d_U35
            r_U36  = a_U36  + d_U36
            r_U37  = a_U37  + d_U37
            r_U38  = a_U38  + d_U38
            r_Np37 = a_Np37 + d_Np37
            r_Np38 = a_Np38 + d_Np38
            r_Np39 = a_Np39 + d_Np39
            r_Pu38 = a_Pu38 + d_Pu38
            r_Pu39 = a_Pu39 + d_Pu39
            r_Pu40 = a_Pu40 + d_Pu40
            r_Pu41 = a_Pu41 + d_Pu41
            r_Pu42 = a_Pu42 + d_Pu42
            r_Pu43 = a_Pu43 + d_Pu43
            r_Am41 = a_Am41 + d_Am41
            r_As42 = a_As42 + d_As42
            r_Am42 = a_Am42 + d_Am42
            r_Am43 = a_Am43 + d_Am43
            r_Am44 = a_Am44 + d_Am44
            r_Cm42 = a_Cm42 + d_Cm42
            r_Cm43 = a_Cm43 + d_Cm43
            r_Cm44 = a_Cm44 + d_Cm44

            A = 0d0
            A(1 , 1 ) = - r_U34
            A(1 , 2 ) = n2n_U35
            A(1 , 9 ) = d_Pu38
            A(2 , 2 ) = - r_U35
            A(2 , 1 ) = c_U34
            A(3 , 3 ) = - r_U36
            A(3 , 2 ) = c_U35
            A(4 , 4 ) = - r_U37
            A(4 , 3 ) = c_U36
            A(4 , 5 ) = n2n_U38
            A(4 , 12) = w_U37 * d_Pu41
            A(5 , 5 ) = - r_U38
            A(5 , 4 ) = c_U37
            A(6 , 6 ) = - r_Np37
            A(6 , 4 ) = d_U37
            A(7 , 7 ) = - r_Np38
            A(7 , 6 ) = c_Np37
            A(7 , 17) = w_Np38 * d_Am42
            A(8 , 8 ) = - r_Np39
            A(8 , 5 ) = c_U38
            A(8 , 7 ) = c_Np38
            A(9 , 9 ) = - r_Pu38
            A(9 , 7 ) = d_Np38
            A(9 , 10) = n2n_Pu39
            A(9 , 20) = d_Cm42
            A(10, 10) = - r_Pu39
            A(10, 8 ) = d_Np39
            A(10, 9 ) = c_Pu38
            A(10, 21) = d_Cm43
            A(11, 11) = - r_Pu40
            A(11, 8 ) = c_Np39
            A(11, 10) = c_Pu39
            A(11, 22) = d_Cm44
            A(12, 12) = - r_Pu41
            A(12, 11) = c_Pu40
            A(13, 13) = - r_Pu42
            A(13, 12) = c_Pu41
            A(13, 16) = w_Pu42 * d_As42
            A(14, 14) = - r_Pu43
            A(14, 13) = c_Pu42
            A(15, 15) = - r_Am41  ! Am-241
            A(15, 12) = w_Am41 * d_Pu41
            A(16, 16) = - r_As42  ! Am-242m
            A(16, 15) = w_As42 * c_Am41
            A(16, 17) = w_As42_IT * d_Am42
            A(17, 17) = - r_Am42  ! Am-242
            A(17, 15) = w_Am42 * c_Am41
            A(18, 18) = - r_Am43
            A(18, 14) = d_Pu43
            A(18, 16) = c_As42
            A(18, 17) = c_Am42
            A(19, 19) = - r_Am44
            A(19, 18) = c_Am43 * 0.06260d0 ! Am243 capture to Am244 (ground) BR=0.06260
            A(20, 20) = - r_Cm42
            A(20, 16) = w_Cm42 * d_As42
            A(21, 21) = - r_Cm43
            A(21, 20) = c_Cm42
            A(22, 22) = - r_Cm44
            A(22, 18) = c_Am43 * (1d0 - 0.06260d0) ! All Am244m beta to Cm244 T_1/2 = 26min
            A(22, 19) = d_Am44
            A(22, 21) = c_Cm43

            Vec_N(1 ) = N_U34_Old  (Ixy,Iz)
            Vec_N(2 ) = N_U35_Old  (Ixy,Iz)
            Vec_N(3 ) = N_U36_Old  (Ixy,Iz)
            Vec_N(4 ) = N_U37_Old  (Ixy,Iz)
            Vec_N(5 ) = N_U38_Old  (Ixy,Iz)
            Vec_N(6 ) = N_Np37_Old (Ixy,Iz)
            Vec_N(7 ) = N_Np38_Old (Ixy,Iz)
            Vec_N(8 ) = N_Np39_Old (Ixy,Iz)
            Vec_N(9 ) = N_Pu38_Old (Ixy,Iz)
            Vec_N(10) = N_Pu39_Old (Ixy,Iz)
            Vec_N(11) = N_Pu40_Old (Ixy,Iz)
            Vec_N(12) = N_Pu41_Old (Ixy,Iz)
            Vec_N(13) = N_Pu42_Old (Ixy,Iz)
            Vec_N(14) = N_Pu43_Old (Ixy,Iz)
            Vec_N(15) = N_Am41_Old (Ixy,Iz)
            Vec_N(16) = N_As42_Old (Ixy,Iz)
            Vec_N(17) = N_Am42_Old (Ixy,Iz)
            Vec_N(18) = N_Am43_Old (Ixy,Iz)
            Vec_N(19) = N_Am44_Old (Ixy,Iz)
            Vec_N(20) = N_Cm42_Old (Ixy,Iz)
            Vec_N(21) = N_Cm43_Old (Ixy,Iz)
            Vec_N(22) = N_Cm44_Old (Ixy,Iz)

            Sum_CRAM_Vec = ( 0d0, 0d0 )
           ! do i = 1, Order_CRAM/2
           !    Sum_CRAM_Vec(:) = Sum_CRAM_Vec(:) &
           !                  & + alpha(i)*MATMUL( Get_invC( dT*A - theta(i)*Mat_I_HN ), Vec_N )
           ! enddo
            Vec_N = alpha_0*Vec_N + 2d0*real(Sum_CRAM_Vec, 8)

            do i=1,22
               if (Vec_N(i)<nd_cutoff) then
                  Vec_N(i)=0d0
               endif
            enddo

            N_U34  (Ixy,Iz) = Vec_N(1 )
            N_U35  (Ixy,Iz) = Vec_N(2 )
            N_U36  (Ixy,Iz) = Vec_N(3 )
            N_U37  (Ixy,Iz) = Vec_N(4 )
            N_U38  (Ixy,Iz) = Vec_N(5 )
            N_Np37 (Ixy,Iz) = Vec_N(6 )
            N_Np38 (Ixy,Iz) = Vec_N(7 )
            N_Np39 (Ixy,Iz) = Vec_N(8 )
            N_Pu38 (Ixy,Iz) = Vec_N(9 )
            N_Pu39 (Ixy,Iz) = Vec_N(10)
            N_Pu40 (Ixy,Iz) = Vec_N(11)
            N_Pu41 (Ixy,Iz) = Vec_N(12)
            N_Pu42 (Ixy,Iz) = Vec_N(13)
            N_Pu43 (Ixy,Iz) = Vec_N(14)
            N_Am41 (Ixy,Iz) = Vec_N(15)
            N_As42 (Ixy,Iz) = Vec_N(16)
            N_Am42 (Ixy,Iz) = Vec_N(17)
            N_Am43 (Ixy,Iz) = Vec_N(18)
            N_Am44 (Ixy,Iz) = Vec_N(19)
            N_Cm42 (Ixy,Iz) = Vec_N(20)
            N_Cm43 (Ixy,Iz) = Vec_N(21)
            N_Cm44 (Ixy,Iz) = Vec_N(22)
         enddo
      enddo

      return
      end subroutine HeavyNuclide


      subroutine FissionProduct(Flag_tr)
      use Inc_3D, only: Flux
      use Inc_Depletion
      use Inc_maXS, only: maXS_f_3D
      use Inc_miXS
      use Inc_Nuclide
      use Inc_Option
      use Inc_Time, only: dT
!     ! use Mod_Operator, only: Get_InvA, Get_InvC
      implicit none
      logical(4), intent(in) :: Flag_tr
      complex(8) :: Sum_CRAM_Mat(N_NdSm, N_NdSm)
      real(8) :: A   (N_NdSm, N_NdSm)
      real(8) :: eAt (N_NdSm, N_NdSm)
      real(8) :: Vec_N(N_NdSm)
      real(8) :: Vec_b(N_NdSm)
      real(8) :: F
      real(8) :: Y_I, Y_X
      real(8) :: a_I, a_X
      real(8) :: r_I, r_X
      real(8) :: p_I, p_X
      real(8) :: e_I, e_X
      real(8) :: b_I, b_X
      real(8) :: I_0, X_0
      real(8) :: Y_Nd7, Y_Nd8, Y_Nd9, Y_Pm7, Y_Ps8, Y_Pm8, Y_Pm9,               Y_Sm9
      real(8) :: a_Nd7, a_Nd8, a_Nd9, a_Pm7, a_Ps8, a_Pm8, a_Pm9, a_Sm7, a_Sm8, a_Sm9
      real(8) :: r_Nd7, r_Nd8, r_Nd9, r_Pm7, r_Ps8, r_Pm8, r_Pm9, r_Sm7, r_Sm8, r_Sm9
      real(8) :: p_Pm9, p_Sm9
      real(8) :: e_Pm9, e_Sm9
      real(8) :: b_Nd7, b_Nd8, b_Nd9, b_Pm7, b_Ps8, b_Pm8, b_Pm9, b_Sm7, b_Sm8, b_Sm9
      real(8) :: Pm9_0, Sm9_0
      integer :: Ixy, Iz, i

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [FissionProduct] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [FissionProduct] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif



#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      if (flag_tr) then
         if ( abs(dT-0d0) < 1d-10 ) return
      endif

      b_I   = beta_I35
      b_X   = beta_Xe35
      b_Nd7 = beta_Nd47
      b_Nd8 = beta_Nd48
      b_Nd9 = beta_Nd49
      b_Pm7 = beta_Pm47
      b_Ps8 = beta_Ps48
      b_Pm8 = beta_Pm48
      b_Pm9 = beta_Pm49
      b_Sm7 = beta_Sm47
      b_Sm8 = beta_Sm48
      b_Sm9 = beta_Sm49

      do Ixy = 1, Nxy
         do Iz = IzFuelBot, IzFuelTop
            if ( I_FARF_1N( I_4Nto1N(Ixy) ) == 1 ) cycle

            F = dot_product( maXS_f_3D(Ixy, Iz, :), Flux(Ixy, Iz, :) )

            if ( OPT_Xe == 0 ) then ! No Xe
               N_I35(Ixy, Iz) = D0
               N_Xe35(Ixy, Iz) = D0
            elseif ( OPT_Xe == 1 ) then ! Equilibrium Xe
               Y_I = Y_I35 (Ixy, Iz)
               Y_X = Y_Xe35(Ixy, Iz)
               a_I = dot_product( miXS_a_I35 (Ixy, Iz, :), Flux(Ixy, Iz, :) )
               a_X = dot_product( miXS_a_Xe35(Ixy, Iz, :), Flux(Ixy, Iz, :) )
               Y_X = Y_Xe35_eff(Ixy, Iz)
               r_I = a_I + b_I
               r_X = a_X + b_X
               p_I = r_I
               p_X = r_X

               N_I35(Ixy, Iz) = Y_I * F / p_I
               N_Xe35(Ixy, Iz) = ( Y_X*F + b_I*N_I35(Ixy, Iz) ) / p_X
            endif

            if ( OPT_SmChain == 1 ) then ! Sm Chain option (1=Default)
               if (OPT_Sm == 0) then ! No Sm
                  N_Nd47(Ixy, Iz) = 0d0
                  N_Nd48(Ixy, Iz) = 0d0
                  N_Nd49(Ixy, Iz) = 0d0
                  N_Pm47(Ixy, Iz) = 0d0
                  N_Ps48(Ixy, Iz) = 0d0
                  N_Pm48(Ixy, Iz) = 0d0
                  N_Pm49(Ixy, Iz) = 0d0
                  N_Sm47(Ixy, Iz) = 0d0
                  N_Sm48(Ixy, Iz) = 0d0
                  N_Sm49(Ixy, Iz) = 0d0
               elseif (OPT_Sm == 1) then ! Eq Sm
                  Y_Nd7 = Y_Nd47(Ixy, Iz)
                  Y_Nd8 = Y_Nd48(Ixy, Iz)
                  Y_Nd9 = Y_Nd49(Ixy, Iz)
                  Y_Pm7 = Y_Pm47(Ixy, Iz)
                  Y_Ps8 = Y_Ps48(Ixy, Iz)
                  Y_Pm8 = Y_Pm48(Ixy, Iz)
                  Y_Pm9 = Y_Pm49(Ixy, Iz)
                  Y_Sm9 = Y_Sm49(Ixy, Iz)

                  a_Nd7 = dot_product( miXS_a_Nd47(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  a_Nd8 = dot_product( miXS_a_Nd48(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  a_Nd9 = dot_product( miXS_a_Nd49(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  a_Pm7 = dot_product( miXS_a_Pm47(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  a_Ps8 = dot_product( miXS_a_Ps48(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  a_Pm8 = dot_product( miXS_a_Pm48(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  a_Pm9 = dot_product( miXS_a_Pm49(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  a_Sm7 = dot_product( miXS_a_Sm47(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  a_Sm8 = dot_product( miXS_a_Sm48(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  a_Sm9 = dot_product( miXS_a_Sm49(Ixy, Iz, :), Flux(Ixy, Iz, :) )

                  r_Nd7 = a_Nd7 + b_Nd7
                  r_Nd8 = a_Nd8 + b_Nd8
                  r_Nd9 = a_Nd9 + b_Nd9
                  r_Pm7 = a_Pm7 + b_Pm7
                  r_Ps8 = a_Ps8 + b_Ps8
                  r_Pm8 = a_Pm8 + b_Pm8
                  r_Pm9 = a_Pm9 + b_Pm9
                  r_Sm7 = a_Sm7 + b_Sm7
                  r_Sm8 = a_Sm8 + b_Sm8
                  r_Sm9 = a_Sm9 + b_Sm9

                  N_Nd47(Ixy, Iz) = Y_Nd7 * F / r_Nd7
                  N_Nd48(Ixy, Iz) = ( Y_Nd8*F + a_Nd7*N_Nd47(Ixy, Iz) ) / r_Nd8
                  N_Nd49(Ixy, Iz) = ( Y_Nd9*F + a_Nd8*N_Nd48(Ixy, Iz) ) / r_Nd9
                  N_Pm47(Ixy, Iz) = ( Y_Pm7*F + b_Nd7*N_Nd47(Ixy, Iz) ) / r_Pm7
                  N_Pm48(Ixy, Iz) = ( Y_Pm8*F + w_Pm48*a_Pm7*N_Pm47(Ixy, Iz) ) / r_Pm8
                  N_Ps48(Ixy, Iz) = ( Y_Ps8*F + w_Ps48*a_Pm7*N_Pm47(Ixy, Iz) + w_Ps48_IT*b_Pm8*N_Pm48(Ixy, Iz) ) / r_Ps8
                  N_Pm49(Ixy, Iz) = ( Y_Pm9*F + b_Nd9*N_Nd49(Ixy, Iz) + a_Ps8*N_Ps48(Ixy, Iz) + a_Pm8*N_Pm48(Ixy, Iz) ) / r_Pm9
                  N_Sm47(Ixy, Iz) = b_Pm7 * N_Pm47(Ixy, Iz) / r_Sm7
                  N_Sm48(Ixy, Iz) = ( b_Ps8*N_Ps48(Ixy, Iz) + w_Sm48*b_Pm8*N_Pm48(Ixy, Iz) + a_Sm7*N_Sm47(Ixy, Iz) ) / r_Sm8
                  N_Sm49(Ixy, Iz) = ( Y_Sm9*F + b_Pm9*N_Pm49(Ixy, Iz) + b_Sm8*N_Sm48(Ixy, Iz) ) / r_Sm9
               endif

            elseif ( OPT_SmChain == 2 ) then ! Sm Chain option (1=Default)

               if ( OPT_Sm == 0 ) then
                  N_Pm49(Ixy, Iz) = D0
                  N_Sm49(Ixy, Iz) = D0
               elseif ( OPT_Sm == 1 ) then
                  Y_Pm9 = Y_Pm49_eff(Ixy, Iz)
                  Y_Sm9 = Y_Sm49    (Ixy, Iz)
                  a_Pm9 = dot_product( miXS_a_Pm49(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  a_Sm9 = dot_product( miXS_a_Sm49(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  r_Pm9 = a_Pm9 + b_Pm9
                  r_Sm9 = a_Sm9 + b_Sm9
                  p_Pm9 = r_Pm9
                  p_Sm9 = r_Sm9
                  N_Pm49(Ixy, Iz) = Y_Pm9 * F / p_Pm9
                  N_Sm49(Ixy, Iz) = ( Y_Sm9*F + b_Pm9*N_Pm49(Ixy, Iz) ) / p_Sm9
               endif
            endif

            if ( ( Flag_tr ) .AND. ( dT /= D0 ) ) then
               F = dot_product( maXS_f_3D(Ixy, Iz, :), Flux(Ixy, Iz, :) )

               if ( OPT_Xe == 2 ) then ! Transient Xe
                  I_0 = N_I35_Old(Ixy, Iz)
                  X_0 = N_Xe35_Old(Ixy, Iz)

                  Y_I = Y_I35 (Ixy, Iz)
                  Y_X = Y_Xe35(Ixy, Iz)
                  Y_X = Y_Xe35_eff(Ixy, Iz)
                  a_I = dot_product( miXS_a_I35 (Ixy, Iz, :), Flux(Ixy, Iz, :) )
                  a_X = dot_product( miXS_a_Xe35(Ixy, Iz, :), Flux(Ixy, Iz, :) )

                  r_I = a_I + b_I
                  r_X = a_X + b_X
                  p_I = r_I
                  p_X = r_X
                  if (p_I*dT>100d0) then
                     e_I = 0d0
                  else
                     e_I = DEXP(- p_I*dT)
                  endif
                  if (p_X*dT>100d0) then
                     e_X = 0d0
                  else
                     e_X = DEXP(- p_X*dT)
                  endif

                  N_I35(Ixy, Iz) = I_0 * e_I + ( Y_I*F / p_I ) * ( D1 - e_I )
                  N_Xe35(Ixy, Iz) = X_0 * e_X + ( D1 - e_X ) / p_X * ( b_I*Y_I*F/p_I + Y_X*F ) &
                                 + ( e_X - e_I ) / ( p_X - p_I ) * ( b_I*Y_I*F/p_I - b_I*I_0 )
               endif

               if ( OPT_SmChain == 1 ) then ! Sm Chain option (1=Default)
                  if ( OPT_Sm == 2 ) then ! Transi ent Sm
                     Y_Nd7 = Y_Nd47(Ixy,Iz)
                     Y_Nd8 = Y_Nd48(Ixy,Iz)
                     Y_Nd9 = Y_Nd49(Ixy,Iz)
                     Y_Pm7 = Y_Pm47(Ixy,Iz)
                     Y_Ps8 = Y_Ps48(Ixy,Iz)
                     Y_Pm8 = Y_Pm48(Ixy,Iz)
                     Y_Pm9 = Y_Pm49(Ixy,Iz)
                     Y_Sm9 = Y_Sm49(Ixy,Iz)

                     a_Nd7 = dot_product( miXS_a_Nd47(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     a_Nd8 = dot_product( miXS_a_Nd48(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     a_Nd9 = dot_product( miXS_a_Nd49(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     a_Pm7 = dot_product( miXS_a_Pm47(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     a_Ps8 = dot_product( miXS_a_Ps48(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     a_Pm8 = dot_product( miXS_a_Pm48(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     a_Pm9 = dot_product( miXS_a_Pm49(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     a_Sm7 = dot_product( miXS_a_Sm47(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     a_Sm8 = dot_product( miXS_a_Sm48(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     a_Sm9 = dot_product( miXS_a_Sm49(Ixy,Iz,:),Flux(Ixy,Iz,:) )

                     ! To Prevent Singula Matrix Cuz there are b_@ = 0
                     if (abs(a_Nd7)<1d-10) a_Nd7=dot_product( miXS_a_Nd47(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     if (abs(a_Nd8)<1d-10) a_Nd8=dot_product( miXS_a_Nd48(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     if (abs(a_Nd9)<1d-10) a_Nd9=dot_product( miXS_a_Nd49(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     if (abs(a_Pm7)<1d-10) a_Pm7=dot_product( miXS_a_Pm47(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     if (abs(a_Ps8)<1d-10) a_Ps8=dot_product( miXS_a_Ps48(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     if (abs(a_Pm8)<1d-10) a_Pm8=dot_product( miXS_a_Pm48(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     if (abs(a_Pm9)<1d-10) a_Pm9=dot_product( miXS_a_Pm49(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     if (abs(a_Sm7)<1d-10) a_Sm7=dot_product( miXS_a_Sm47(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     if (abs(a_Sm8)<1d-10) a_Sm8=dot_product( miXS_a_Sm48(Ixy,Iz,:),Flux(Ixy,Iz,:) )
                     if (abs(a_Sm9)<1d-10) a_Sm9=dot_product( miXS_a_Sm49(Ixy,Iz,:),Flux(Ixy,Iz,:) )

                     r_Nd7 = a_Nd7 + b_Nd7
                     r_Nd8 = a_Nd8 + b_Nd8
                     r_Nd9 = a_Nd9 + b_Nd9
                     r_Pm7 = a_Pm7 + b_Pm7
                     r_Ps8 = a_Ps8 + b_Ps8
                     r_Pm8 = a_Pm8 + b_Pm8
                     r_Pm9 = a_Pm9 + b_Pm9
                     r_Sm7 = a_Sm7 + b_Sm7
                     r_Sm8 = a_Sm8 + b_Sm8
                     r_Sm9 = a_Sm9 + b_Sm9

                     A = D0
                     A(1 , 1 ) = - r_Nd7
                     A(2 , 1 ) = a_Nd7
                     A(2 , 2 ) = - r_Nd8
                     A(3 , 2 ) = a_Nd8
                     A(3 , 3 ) = - r_Nd9
                     A(4 , 1 ) = b_Nd7
                     A(4 , 4 ) = - r_Pm7
                     A(5 , 4 ) = w_Ps48 * a_Pm7
                     A(5 , 5 ) = - r_Ps8
                     A(5 , 6 ) = w_Ps48_IT * b_Pm8
                     A(6 , 4 ) = w_Pm48 * a_Pm7
                     A(6 , 6 ) = - r_Pm8
                     A(7 , 3 ) = b_Nd9
                     A(7 , 5 ) = a_Ps8
                     A(7 , 6 ) = a_Pm8
                     A(7 , 7 ) = - r_Pm9
                     A(8 , 4 ) = b_Pm7
                     A(8 , 8 ) = - r_Sm7
                     A(9 , 5 ) = b_Ps8
                     A(9 , 6 ) = w_Sm48 * b_Pm8
                     A(9 , 8 ) = a_Sm7
                     A(9 , 9 ) = - r_Sm8
                     A(10, 7 ) = b_Pm9
                     A(10, 9 ) = a_Sm8
                     A(10, 10) = - a_Sm9

                     Vec_N(1 ) = N_Nd47_Old(Ixy, Iz)
                     Vec_N(2 ) = N_Nd48_Old(Ixy, Iz)
                     Vec_N(3 ) = N_Nd49_Old(Ixy, Iz)
                     Vec_N(4 ) = N_Pm47_Old(Ixy, Iz)
                     Vec_N(5 ) = N_Ps48_Old(Ixy, Iz)
                     Vec_N(6 ) = N_Pm48_Old(Ixy, Iz)
                     Vec_N(7 ) = N_Pm49_Old(Ixy, Iz)
                     Vec_N(8 ) = N_Sm47_Old(Ixy, Iz)
                     Vec_N(9 ) = N_Sm48_Old(Ixy, Iz)
                     Vec_N(10) = N_Sm49_Old(Ixy, Iz)

                     Vec_b(1 ) = Y_Nd7 * F
                     Vec_b(2 ) = Y_Nd8 * F
                     Vec_b(3 ) = Y_Nd9 * F
                     Vec_b(4 ) = Y_Pm7 * F
                     Vec_b(5 ) = Y_Ps8 * F
                     Vec_b(6 ) = Y_Pm8 * F
                     Vec_b(7 ) = Y_Pm9 * F
                     Vec_b(8 ) = D0
                     Vec_b(9 ) = D0
                     Vec_b(10) = Y_Sm9 * F

                     Sum_CRAM_Mat = ( D0, D0 )
                    ! do i = 1, Order_CRAM/2
                    !    Sum_CRAM_Mat = Sum_CRAM_Mat + alpha(i) * Get_invC( dT*A - theta(i)*Mat_I_FP )
                    ! enddo
                     eAt = D0
                     eAt = alpha_0*real(Mat_I_FP,8) + 2d0*real(Sum_CRAM_Mat, 8)

                   !  Vec_N = MATMUL(eAt,Vec_N)+MATMUL((eAt-real(Mat_I_FP,8)),MATMUL(Get_invA(A),Vec_b))

                     do i=1,10
                        if (Vec_N(i)<nd_cutoff) then
                           Vec_N(i)=0d0
                        endif
                     enddo

                     N_Nd47(Ixy, Iz) = Vec_N(1 )
                     N_Nd48(Ixy, Iz) = Vec_N(2 )
                     N_Nd49(Ixy, Iz) = Vec_N(3 )
                     N_Pm47(Ixy, Iz) = Vec_N(4 )
                     N_Ps48(Ixy, Iz) = Vec_N(5 )
                     N_Pm48(Ixy, Iz) = Vec_N(6 )
                     N_Pm49(Ixy, Iz) = Vec_N(7 )
                     N_Sm47(Ixy, Iz) = Vec_N(8 )
                     N_Sm48(Ixy, Iz) = Vec_N(9 )
                     N_Sm49(Ixy, Iz) = Vec_N(10)
                  endif

               elseif ( OPT_SmChain == 2 ) then ! Sm Chain option (1=Default)

                  if ( OPT_Sm == 2 ) then
                     Pm9_0 = N_Pm49_Old(Ixy, Iz)
                     Sm9_0 = N_Sm49_Old(Ixy, Iz)

                     Y_Pm9 = Y_Pm49_eff(Ixy, Iz)
                     Y_Sm9 = Y_Sm49    (Ixy, Iz)
                     a_Pm9 = dot_product( miXS_a_Pm49(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                     a_Sm9 = dot_product( miXS_a_Sm49(Ixy, Iz, :), Flux(Ixy, Iz, :) )

                     r_Pm9 = a_Pm9 + b_Pm9
                     r_Sm9 = a_Sm9 + b_Sm9
                     p_Pm9 = r_Pm9
                     p_Sm9 = r_Sm9
                     e_Pm9 = DEXP(- p_Pm9*dT)
                     e_Sm9 = DEXP(- p_Sm9*dT)

                     N_Pm49(Ixy, Iz) = Pm9_0 * e_Pm9 + ( Y_Pm9*F / p_Pm9 ) * ( D1 - e_Pm9 )
                     N_Sm49(Ixy, Iz) = Sm9_0 * e_Sm9                                                   &
                      + ( D1 - e_Sm9 ) / p_Sm9 * ( b_Pm9*Y_Pm9*F/p_Pm9 + Y_Sm9*F )                     &
                      + ( e_Sm9 - e_Pm9 ) / ( p_Sm9 - p_Pm9 ) * ( b_Pm9*Y_Pm9*F/p_Pm9 - b_Pm9*Pm9_0 )
                  endif
               endif
            endif
         enddo
      enddo

      return
      end subroutine FissionProduct


      subroutine BurnableAbsorber
      use Inc_3D, only: Flux
      use Inc_Depletion
      use Inc_miXS
      use Inc_Nuclide
      use Inc_Option
      use Inc_Time, only: dT
      use Inc_XS_File, only: I_Tab
!     ! use Mod_Operator, only: Get_InvA, Get_InvC
      implicit none
      real(8) :: Gd152, Gd160, Gdeff
      real(8) :: a_Gd152, a_Gd160, a_Gdeff
      integer :: Ixy, Ixy_FA, Iz

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [BurnableAbsorber] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [BurnableAbsorber] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.flag_midepl) return


#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      if ( abs(dT-0d0)<1d-10 ) return

      do Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)
         do Iz = IzFuelBot, IzFuelTop
            I_Tab = new_asym_itab(AxialComp(I_LP_1N(I_4Nto1N(Ixy)),Iz),Ixy)
            if (.not.Flag_BP(I_Tab)) cycle

            Gd152 = N_Gd52_Old(Ixy,Iz)
            Gd160 = N_Gd60_Old(Ixy,Iz)
            Gdeff = N_Gd58_Old(Ixy,Iz)

            a_Gd152 = dot_product(miXS_a_Gd52(Ixy,Iz,:),Flux(Ixy,Iz,:))
            a_Gd160 = dot_product(miXS_a_Gd60(Ixy,Iz,:),Flux(Ixy,Iz,:))
            a_Gdeff = dot_product(miXS_a_Gd58(Ixy,Iz,:),Flux(Ixy,Iz,:))

            N_Gd52(Ixy,Iz) = Gd152 * exp(-a_Gd152*dT)
            N_Gd60(Ixy,Iz) = Gd160 * exp(-a_Gd160*dT)
            N_Gd58(Ixy,Iz) = Gdeff * exp(-a_Gdeff*dT)
         enddo
      enddo

      if (Opt_Gd == 0) then
         N_Gd52(:,:) = 0d0
         N_Gd54(:,:) = 0d0
         N_Gd55(:,:) = 0d0
         N_Gd56(:,:) = 0d0
         N_Gd57(:,:) = 0d0
         N_Gd58(:,:) = 0d0
         N_Gd60(:,:) = 0d0
      endif

      return
      end subroutine BurnableAbsorber


      subroutine cal_crud_info
      use inc_crud, only: opt_crud
      implicit none

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [cal_crud_info] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [cal_crud_info] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      return
      end subroutine cal_crud_info


    subroutine Save_B0_RRa
    use Inc_Crud, only: B0_RRa, opt_readcrud,n_crudasi
    use inc_3d, only: flux_old,flux
    use inc_depletion, only: flag_pc
    use inc_mixs, only: miXS_a_B0_Old,miXS_a_B0
    implicit none
    integer(4) :: iz,ixy

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Save_B0_RRa] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Save_B0_RRa] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    n_crudasi=0

    if(.not.opt_readcrud) then
        if(.not.Flag_PC) then
            do iz = 1, nz
                do ixy=1,nxy
                    B0_RRa(ixy,iz)=dot_product( miXS_a_B0_Old (Ixy, Iz, :), Flux_Old(Ixy, Iz, :) )
                enddo
            enddo
        else
            do iz = 1, nz
                do ixy=1,nxy
                    B0_RRa(ixy,iz)=dot_product( miXS_a_B0(Ixy, Iz, :), Flux(Ixy, Iz, :) )
                enddo
            enddo
        endif
    endif

    return
    end subroutine Save_B0_RRa

    subroutine Save_crud_info
    use inc_crud, only: opt_crud
    use Inc_Crud, only: n_CrudB0,n_crudb0_old,cruddr,cruddr_old,total_crud,hist_n_CrudB0,hist_n_crudb0_z
#ifdef hj_training
    use Inc_Crud, only: hist_crudt, hist_crudt_z,crud_ifa,crud_nfa
    use Inc_Geometry, only: i_1nto4n
#endif
    use inc_xs_file, only: i_bu
    implicit none
    integer :: ixy,iz
#ifdef hj_training
    integer :: m,i
#endif
    real(8) :: tmp2,vol


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Save_crud_info] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Save_crud_info] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    return
    end subroutine Save_crud_info

    subroutine Crud_Itr_ASIsearch(if_conv)
    use Inc_XYZ, only: ASI
    use Inc_Crud, only: TD_Crud_Asi, TD_Crud_Bnden,n_CRUDASI,ASI_o,ASI_oo,NB_o,NB_oo,cruddr,if_asi_search,eps_crudasi
    use inc_xs_file, only: i_bu
    use inc_parallel, only: comm
    implicit none
    logical(1),intent(inout) :: if_conv
    real(8) :: tmp1
    real(8) :: damp=1.0d0
    real(8) :: slope
    real(8) :: err1,asi_err

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Crud_Itr_ASIsearch] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Crud_Itr_ASIsearch] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    if_conv=.true.

    if(.not.if_asi_search(i_bu)) then
       if(i_bu==1)  then
          TD_Crud_Bnden(i_BU)= 0.0d0
       endif
       if_conv=.true.
       if (comm%if_master) write(*,*) '*********** Skip ASI seach',if_asi_search(i_bu),if_conv
       return
    endif
    NB_oo=NB_o
    NB_o=TD_Crud_Bnden(i_BU)
    ASI_oo=ASI_o
    ASI_o=ASI
    n_CrudASI=n_CrudASI+1
    if(abs(ASI_o-ASI_oo)<1d-30 .or. n_CRUDASI<=2) then
        slope=5D+22
    else
        slope=(NB_o-NB_oo)/(ASI_o-ASI_oo)
    endif
    slope=min(1d+26,slope)
    slope=max(0.0d0,slope)

    tmp1=TD_Crud_Bnden(i_BU)
    TD_Crud_Bnden(i_BU)= &
    (TD_Crud_ASI(i_BU)-ASI)*slope+TD_Crud_Bnden(i_BU)

    !write(*,*) 'Slope =',slope

    TD_Crud_Bnden(i_BU)=damp*TD_Crud_Bnden(i_BU)+(1d0-damp)*tmp1

!    TD_Crud_Bnden(i_BU)=Max(0.0d0,TD_Crud_Bnden(i_BU))
    if(maxval(cruddr)<1d-5) then
        TD_Crud_Bnden(i_BU)=0d0
    endif

    asi_err=abs(ASI_o-ASI_oo)
    err1=abs(ASI-TD_Crud_ASI(i_BU))
    if(((err1<eps_crudasi).or. asi_err<eps_crudasi*0.01d0) .and. n_crudasi>3) then
        if_conv=.true.
    else
        if_conv=.false.
    endif

    if (comm%if_master) write(*,'(a)') '* Update crud infomation - ASI search'
    if (comm%if_master) write(*,'(a)') '----------------------------------------------'
    if (comm%if_master) write(*,'(a,2f13.4)')  '- Target ASI, ASI  =',TD_Crud_ASI(i_BU),ASI
    if (comm%if_master) write(*,'(a,2f13.4,f13.7)') '- Error/Suc error  =',err1,asi_err
    if (comm%if_master) write(*,'(a,2es13.5)')      '- Old/New B0 nden  =',tmp1,TD_Crud_Bnden(i_BU)
    if (comm%if_master) write(*,'(a,es13.5)')  '- Crud B0 #/cm3    =',Td_Crud_Bnden(i_BU)
    if (comm%if_master) write(*,'(a)') '----------------------------------------------'
    !call Get_Cruddr
    call Get_Crud_Nden ! use saved B0_RRa
    call Get_Crud_dtf
    return
    end subroutine Crud_Itr_ASIsearch

    subroutine Crud_Itr_FAASIsearch(if_nodal_conv,if_conv)
    use Inc_Crud, only: itr_FA_BD,itr_fa_asi, n_CRUDASI,&
                        if_asi_search,&
                        FAASI_INP,eps_crudasi,slope_fa_bd
    use inc_xs_file, only: i_bu
    use inc_xyz, only: fa_asi
    use inc_parallel, only: comm
    implicit none
    logical(4),intent(in) :: if_nodal_conv
    logical(1),intent(inout) :: if_conv
    integer(4) :: iy_1n,ix_1n
    real(8) :: tmp1
    real(8) :: damp=1.0d0
    real(8) :: slope
    real(8) :: err1,n_err1,maxerr1,slope_max,slope_min
    real(8) :: asi_err,crud_suc_err

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Crud_Itr_FAASIsearch] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Crud_Itr_FAASIsearch] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    if(.not.if_nodal_conv) then
        return
    endif

    if(.not.if_asi_search(i_bu)) then
       if(i_bu==1)  then
          itr_FA_BD(1,:,:)=0.0d0
       endif
       if_conv=.true.
       return
    endif
    itr_FA_ASI(1,:,:)=FA_ASI(:,:)

    n_CrudASI=n_CrudASI+1
    err1=0.0d0
    n_err1=0.0d0
    maxerr1=0.0d0
    asi_err=0.0
    crud_suc_err=0d0
    slope_fa_bd=0.0d0
    slope_max=0d0
    slope_min=1D+40
    do iy_1n=1,ny_1n
        do ix_1n=1,nx_1n
            if(abs(FA_ASI(ix_1n,iy_1n))>1d-30 .or. abs(FAASI_INP(ix_1n,iy_1n,i_bu))>1d-30) then ! fuel
                tmp1=abs(FA_ASI(ix_1n,iy_1n)-FAASI_INP(ix_1n,iy_1n,i_bu))
                err1=err1+tmp1
                maxerr1=max(maxerr1,tmp1)
                n_err1=n_err1+1
            else
                cycle
            endif
            itr_FA_BD(3,ix_1n,iy_1n)=itr_FA_BD(2,ix_1n,iy_1n)
            itr_FA_BD(2,ix_1n,iy_1n)=itr_FA_BD(1,ix_1n,iy_1n)

            itr_FA_ASI(3,ix_1n,iy_1n)=itr_FA_ASI(2,ix_1n,iy_1n)
            itr_FA_ASI(2,ix_1n,iy_1n)=itr_FA_ASI(1,ix_1n,iy_1n)
            asi_err=max(asi_err,abs(itr_FA_ASI(3,ix_1n,iy_1n)-itr_FA_ASI(2,ix_1n,iy_1n)))

            if(abs(itr_FA_ASI(3,ix_1n,iy_1n)-itr_FA_ASI(2,ix_1n,iy_1n))<1d-30 .or. n_CRUDASI<=2) then
                slope=1D+23
            else
                slope=(itr_FA_BD(2,ix_1n,iy_1n)-itr_FA_BD(3,ix_1n,iy_1n))/(itr_FA_ASI(2,ix_1n,iy_1n)-itr_FA_ASI(3,ix_1n,iy_1n))
            endif
            slope=min(1d+26,slope)
            slope=max(1D+10,slope)
            slope=max(slope,-itr_FA_BD(1,ix_1n,iy_1n)/(FAASI_INP(ix_1n,iy_1n,i_bu)-itr_FA_ASI(1,ix_1n,iy_1n)))
            if(abs(FAASI_INP(ix_1n,iy_1n,i_bu)-itr_FA_ASI(1,ix_1n,iy_1n))< eps_crudasi) then
                slope=0.0d0
            endif

            slope_max=max(slope_max,slope)
            slope_min=min(slope_min,slope)

            slope_fa_bd(ix_1n,iy_1n)=slope
            tmp1=itr_FA_BD(1,ix_1n,iy_1n)
            itr_FA_BD(1,ix_1n,iy_1n)=(FAASI_INP(ix_1n,iy_1n,i_bu)-itr_FA_ASI(1,ix_1n,iy_1n))*slope+itr_FA_BD(1,ix_1n,iy_1n)
            itr_FA_BD(1,ix_1n,iy_1n)=damp*itr_FA_BD(1,ix_1n,iy_1n)+(1d0-damp)*tmp1
            itr_FA_BD(1,ix_1n,iy_1n)=Max(0.0d0,itr_FA_BD(1,ix_1n,iy_1n))
            crud_suc_err= &
              & max(crud_suc_err,abs(itr_FA_BD(1,ix_1n,iy_1n)-itr_FA_BD(2,ix_1n,iy_1n))/max(1.0D-10,itr_FA_BD(1,ix_1n,iy_1n)))
        enddo
    enddo
    err1=err1/max(1d0,n_err1)
    if(((err1<eps_crudasi*0.2d0 .and. maxerr1<eps_crudasi ).or. crud_suc_err<1.0D-5) .and. n_crudasi>3) then
        if_conv=.true.
    else
        if_conv=.false.
    endif

    if (comm%if_master) write(*,'(a)') '* Update crud infomation - FA ASI search'
    if (comm%if_master) write(*,'(a)') '----------------------------------------------'
    if (comm%if_master) write(*,'(a,2f13.4,es13.4,4l5)') '- Avg/Max/Suc error=', &
       & err1,maxerr1,crud_suc_err,err1<eps_crudasi*0.2d0, maxerr1<eps_crudasi,crud_suc_err<1.0D-5,if_conv
    if (comm%if_master) write(*,'(a,3es13.4)') 'Slope=',sum(slope_fa_bd)/dble(n_err1),slope_max,slope_min
    if (comm%if_master) write(*,'(a)') '----------------------------------------------'
    call Get_Crud_Nden ! use saved B0_RRa
    call Get_Crud_dtf
    return
    end subroutine Crud_Itr_FAASIsearch

    subroutine Crud_Itr_CBCCOR
    use Inc_Crud, only: TD_Crud_Bnden,Crud_Bnden
    use inc_xs_file, only: i_bu
    use Inc_TH, only: ppm
    use inc_parallel, only: comm
    implicit none
    real(8) :: tmp1
!    real(8) :: damp=0.7d0
    real(8) :: ppm_Factor


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Crud_Itr_CBCCOR] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Crud_Itr_CBCCOR] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    ppm_factor=max(0.0d0,ppm/800d0)

    tmp1=TD_Crud_Bnden(i_BU)

    TD_Crud_Bnden(i_BU)=Max(0.0d0,Crud_BNDEN*ppm_factor)

!    TD_Crud_Bnden(i_BU)=damp*TD_Crud_Bnden(i_BU)+(1d0-damp)*tmp1


    if (comm%if_master) write(*,'(a)') '* Update crud infomation - CBC correction'
    if (comm%if_master) write(*,'(a)') '----------------------------------------------'
    if (comm%if_master) write(*,'(a,f13.5)')  '- PPM correction   =',ppm_factor
    if (comm%if_master) write(*,'(a,2es13.5)')'- Old/New B0 nden  =',tmp1,TD_Crud_Bnden(i_BU)
    if (comm%if_master) write(*,'(a,es13.5)') '- Crud boron #/cm3 =',Td_Crud_Bnden(i_BU)
    if (comm%if_master) write(*,'(a)') '----------------------------------------------'
    !call Get_Cruddr
    call Get_Crud_Nden ! use saved B0_RRa
    call Get_Crud_dtf
    return
    end subroutine Crud_Itr_CBCCOR

    subroutine Get_Crud_Nden
    use Inc_Crud, only: n_CrudB0,n_crudb0_old,td_Crud_Bnden,&
                        crud_poro,opt_readcrud,opt_crudburn,B0_RRa,    &
                        opt_crudfaasi,itr_fa_bd,burn_factor,opt_crudboa,day_boa,boa_b0,n_day_boa
    use inc_xs_file, only: i_bu
    use inc_geometry, only: ixtoix_1n,iytoiy_1n,ixiytoixy,ny,Ix_Start_y,Ix_end_y
    use Inc_Crud, only: cruddr,cruddr_old,cruddr_thres
    use inc_crud, only: opt_crud_model
    use inc_fa, only: R_Clad, n_pin, h_FA
    use inc_time, only: dt,day
    use inc_parallel, only: comm
    implicit none

    integer(4) :: iz,ixy,ix_1n,iy_1n,ix,iy,i
    real(8) :: vol_Factor
    real(8) :: avg, nn, max1
    real(8) :: RRa,Feed,tmp_dT,tmp1
    integer(4) :: i1,i2
    real(8) :: frac1,frac2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_Crud_Nden] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_Crud_Nden] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    n_CrudB0(:,:)=0.0d0
    nn=0d0
    avg=0d0
    max1=0d0
    if(.not.opt_readcrud) then
        tmp_dT=max(1d-10,dT)
        if(opt_crudfaasi) then
            do iz = 1, nz
                do iy=1,ny
                    iy_1n=iytoiy_1n(iy)
                    do ix=Ix_Start_y(Iy), Ix_end_y(Iy)
                        ix_1n=ixtoix_1n(ix)
                        ixy=IxIyToIxy(ix,iy)
                       ! write(*,*) iz,ix,iy,ixy
                        vol_factor=( (r_clad+cruddr(ixy,iz)*1d-6)**2 &
                           & - (r_clad+cruddr_old(ixy,iz)*1d-6 )**2 )*1d4*pi*n_pin/(h_fa*h_fa)
                        vol_Factor=vol_Factor*(1d0-crud_poro)
                        vol_Factor=max(0d0,vol_factor)
                        !Feed=TD_Crud_Bnden(i_BU)*vol_factor/tmp_dT
                        Feed=itr_FA_BD(1,ix_1n,iy_1n)*vol_factor/tmp_dT
                        RRa=B0_RRa(Ixy, Iz)*burn_factor
                        if(RRa>1d-10 .and. opt_crudburn) then
                            ! Crud boron depletion
                            n_CrudB0(ixy,iz)=Feed/RRa-(Feed/RRa-n_CrudB0_old(ixy,iz))*exp(-RRa*tmp_dT)
                        else
                            ! No boron depletion
                            n_CrudB0(ixy,iz)=Feed*tmp_dT+n_CrudB0_old(ixy,iz)
                        endif
                        n_CrudB0(ixy,iz)=max(0.0d0,itr_FA_BD(1,ix_1n,iy_1n)*vol_factor)
                    enddo
                enddo
            enddo
        elseif(.not.opt_crudfaasi) then
            do iz = 1, nz
                do ixy=1,nxy
                    vol_factor=( (r_clad+cruddr(ixy,iz)*1d-6)**2 - (r_clad+cruddr_old(ixy,iz)*1d-6 )**2 )*1d4*pi*n_pin/(h_fa*h_fa)
                    vol_Factor=vol_Factor*(1d0-crud_poro)
                    vol_Factor=max(0d0,vol_factor)
                    if (opt_crud_model==0) then
                       Feed=TD_Crud_Bnden(i_BU)*vol_factor/tmp_dT
                    elseif (opt_crud_model==1) then
                       if(cruddr(ixy,iz) > cruddr_thres) then
                          Feed=TD_Crud_Bnden(i_BU)*vol_factor/tmp_dT
                       else
                          Feed=0.d0
                       endif
                    endif
                    RRa=B0_RRa(Ixy, Iz)*burn_factor
                    if(RRa>1d-10 .and. opt_crudburn) then
                        ! Crud boron depletion
                        n_CrudB0(ixy,iz)=Feed/RRa-(Feed/RRa-n_CrudB0_old(ixy,iz))*exp(-RRa*tmp_dT)
                    else
                        ! No boron depletion
                        n_CrudB0(ixy,iz)=Feed*tmp_dT+n_CrudB0_old(ixy,iz)
                    endif
                    n_CrudB0(ixy,iz)=max(0.0d0,n_CrudB0(ixy,iz))
                enddo
            enddo
        endif
    elseif(opt_readcrud .and. opt_crudboa) then
        !if_exist=.false.
        !do i=n_day_boa,1,-1
        !    if(day>day_boa(i)) then
        !        if_exist=.true.
        !        exit
        !    endif
        !enddo
        !if(if_exist) then
        !    tmp1=1d0 / 10.7191d0 * NA
        !    !        B0_avg_mass  avo
        !    do iz=1,nz
        !        do ixy=1,nxy
        !            n_crudb0(ixy,iz)=boa_b0(ixy,iz,i)/NodeVolume(Ixy, Iz)*tmp1
        !        enddo
        !    enddo
        !else
        !    n_crudb0=0d0
        !endif

        i1=1
        do i=1,n_day_boa-1
            if(day<=day_boa(i)) then
                exit
            endif
            i1=i
        enddo
        i2=n_day_boa
        do i=n_day_boa,2,-1
            if(day>day_boa(i)) then
                exit
            endif
            i2=i
        enddo
        frac1=(day_boa(i2)-day)/(day_boa(i2)-day_boa(i1))
        frac2=1d0-frac1
        n_crudb0(:,:)=boa_b0(:,:,i1)*frac1+boa_b0(:,:,i2)*frac2
        do iz=1,nz
            do ixy=1,nxy
                n_crudb0(ixy,iz)=max(0.0d0,n_crudb0(ixy,iz))
            enddo
        enddo

        tmp1=1d0 / 10.7191d0 * NA
        !        B0_avg_mass  avo
        do iz=1,nz
            do ixy=1,nxy
                n_crudb0(ixy,iz)=n_crudb0(ixy,iz)/NodeVolume(Ixy, Iz)*tmp1
            enddo
        enddo

    elseif(opt_readcrud) then
        stop 'Indexing was changed. do not use anymore'
        !do ixy_1n=1,nxy_1n
        !    do iz = 1, nz
        !        do m = 1, 4
        !            if (i_1nto4n(ixy_1n, m) == 0) then
        !                cycle
        !            else
        !                n_CrudB0(i_1nto4n(ixy_1n, m),iz)=crud_table(ixy_1n,iz,i_bu)
        !            endif
        !            if(n_CrudB0(i_1nto4n(ixy_1n, m),iz)>0d0) then
        !                nn=nn+1d0
        !                avg=avg+n_CrudB0(i_1nto4n(ixy_1n, m),iz)
        !                max1=max(max1,n_CrudB0(i_1nto4n(ixy_1n, m),iz))
        !            endif
        !        enddo
        !    enddo
        !enddo
    endif

    nn=0d0
    avg=0d0
    max1=0d0
    do iz=1,nz
        do ixy=1,nxy
            if(n_CrudB0(ixy,iz)>0d0) then
                nn=nn+1d0
                avg=avg+n_CrudB0(ixy,iz)
                max1=max(max1,n_CrudB0(ixy,iz))
            endif
        enddo
    enddo

    if (comm%if_master) write(*,'(a,2es13.5)') '- AVG/Max N_B0     =',avg/max(1d0,nn),max1

    return
    end subroutine Get_Crud_Nden


    subroutine Get_SNB_mass
    use Inc_Depletion, only: cycle_BU
    use Inc_Pinvar, only: npin
    use Inc_PinPOW, only: HFF
    use Inc_XS_File, only: I_BU
    use Inc_Crud
    use Inc_Constant, only: DM2
    use Inc_Geometry, only: Nz, I_4Nto1N, h_z  ! axial node height (cm)
    use Inc_RP, only: I_FA  ! ixy_FA to ixy
    use Inc_3D, only: T_Mod, &  ! moderator temperature (K)
                    & D_Mod     ! moderator density (g/cm3)
    use Inc_TH, only: deq,   &  ! equivalent diameter (m)
    !                & nrp4,  &  ! nr + 4, node number at cladding outer surface
    !                & tfuel,tcool,&  ! fuel and coolant temperature (degC)
    !                & htcoef,&  ! heat transfer coefficient (W/m2/K)
                    & rhou,  &  ! coolant mass flux (kg/m2/sec)
                    & zeta,  &  ! heated perimeter (m)
                    & qflux     ! heat flux (W/m2)
    use Inc_Time, only: dT
    use Mod_THfunc, only: fenthal, &  ! water enthalpy (J/kg)
                        & fhcap,   &  ! water specific heat capacity (J/kg/K)
                        & fcond,   &  ! water thermal conductivity (W/m/K)
                        & fdensh      ! water density (kg/m3)
    implicit none
    integer(4) :: Iz, Ixy, Ixy_FA, Ixy_1N
    real(8) :: degC   ! coolant temperature (C)
    real(8) :: hf     ! liquid enthalpy (J/kg)
    real(8), parameter :: hfsat = 1634290.0d0  ! saturation enthalpy of saturated water (J/kg) at 155.11 bar
    real(8), parameter :: hgsat = 2592820.0d0  ! saturation enthalpy of saturated vapor (J/kg) at 155.11 bar
    real(8), parameter :: hfg   = hgsat - hfsat  ! latent heat (J/kg)
    real(8) :: Nu     ! Nusselt Number
    real(8) :: Pe     ! Pe Number
    real(8) :: St     ! St Number
    real(8) :: cpf    ! liquid specific heat (J/kg K)
    real(8) :: kf     ! liquid thermal conductivity (W/m/K)
    real(8) :: ep     ! pumping factor in Lahey model
    real(8) :: chi    ! fraction of wall heat flux used for liquid evaporation
    real(8) :: hcrit  ! critical enthalpy for SNB (J/kg)
    real(8) :: rhof   ! liquid density (kg/m3)
    real(8) :: rhog = 102.04  ! saturated vapor density (kg/m3) at 155.11 bar
    real(8) :: SNB_rate  ! subcooled nucleate boiling rate (kg/sec) for each node
    real(8) :: tot_SNB_rate  ! core subcooled nucleate boiling rate (kg/sec)
    !real(8) :: Aw     ! wall area (m2)
    !real(8) :: V      ! moderator volume (m3)
    !real(8) :: max_form_factor
    integer :: ip, jp
    real(8) :: npin2
    real(8) :: pin_qflux

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_SNB_mass] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_SNB_mass] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif


    if (opt_crud_model/=1) return

    npin2 = npin*npin

    tot_SNB_rate = 0.d0
    do Ixy_FA = 1, Nxy_FA
       Ixy = I_FA(Ixy_FA)
       Ixy_1N = I_4Nto1N(Ixy)
       do Iz = 1, Nz
          do ip = 1, npin
             do jp = 1, npin
                ! set parameters
!                max_form_factor = maxval(HFF(Ixy_1N,Iz,:,:))
                degC = T_Mod(Ixy,Iz) - 273.15  ! (degC)
                rhof = D_Mod(Ixy,Iz)*1.d3  ! (kg/m3)
                hf = fenthal(degC) * HFF(Ixy_1N,Iz,ip,jp) ! (J/kg)
!                hf = fenthal(degC) * max_form_factor ! (J/kg)
                cpf = fhcap(degC)  ! (J/kg/K)
                kf = fcond(degC)  ! (W/m/K)

                ! pumping factor
                ep = rhof*( hfsat - min(hf,hfsat) ) / ( rhog*hfg )

                ! critical enthalpy (Saha-Zuber correlation) (J/kg)
                ! Model used in SPACE code (Ha (2004), modified Saha & Zuber NVG)
                pin_qflux= qflux(Iz,Ixy_FA) * HFF(Ixy_1N,Iz,ip,jp)
                if ( I_BU > 1 ) then
                   Nu = pin_qflux*deq/kf
                   Pe = rhou(Iz,Ixy_FA)*deq*cpf/kf
                   St = Nu/Pe
                   if (Pe > 5.2d4) then
                      hcrit = hfsat - St*Pe**0.124*cpf/0.0287d0
                   else
                      hcrit = hfsat - St*Pe**1.08*cpf/918.525d0
                   endif
!                   if (Pe > 7d4) then
!                      hcrit = hfsat - St*cpf/0.0065d0
!                   else
!                      hcrit = hfsat - Nu*cpf/455d0
!                   endif

                   ! fraction of wall heat flux used for liquid evaporation
                   chi = ( min(hf,hfsat) - hcrit ) / ( (hfsat - hcrit)*(1.d0 + ep) )
                   chi = max(0.d0, chi)

                   ! subcooled nucleate boiling rate for each node (kg/sec)
                   SNB_rate = pin_qflux * zeta * h_z(Iz)*DM2 / hfg * chi / npin2
                   tot_SNB_rate = tot_SNB_rate + SNB_rate
!                   print *, hcrit, hfsat, hf

                   ! SNB mass (kg)
                   SNB_mass(Ixy,Iz) = SNB_mass(Ixy,Iz) + SNB_rate * dT
                   max_SNB_mass = max(max_SNB_mass, SNB_mass(Ixy,Iz))

                endif
             enddo
          enddo
          !write(200416,'( 2(i5,X),4(f15.5,X) )') Ixy, Iz, SNB_rate, SNB_mass(Ixy,Iz), chi, pin_qflux
       enddo
    enddo
    print *,  "Core Boiling (kg/sec)", tot_SNB_rate, cycle_bu(I_BU)
    !write(200416,'(A,f15.5,X,f15.5)') "Core Boiling (kg/sec)", tot_SNB_rate, cycle_bu(I_BU)
    end subroutine Get_SNB_mass


    subroutine Get_Cruddr
    use Inc_Crud, only: SNB_mass, boiling_factor
    use Inc_Crud, only: opt_crud_model
    use Inc_Crud, only: cruddr,cruddr_table,crud_drf,opt_readcrud,Ind_BU,opt_crudboa,day_boa,n_day_boa,boa_dr
#ifdef hj_training
    use Inc_Crud, only: opt_crud_trset, crud_nfa, crud_ifa, bogus_bu
#endif
    use inc_xs_file, only: i_bu
    use inc_geometry, only: nxy_1n,I_1Nto4N
    use inc_time, only: day
    use inc_parallel, only: comm
    implicit none
    integer(4) :: iz,ixy_1n,m,i
    integer(4) :: ixy
    real(8) :: tmp1
    real(8) :: centz(IzFuelBot:IzFuelTop)
    real(8) :: crud_factor(IzFuelBot:IzFuelTop)
    real(8) :: frac1,frac2
    integer(4),parameter :: n_data=17
    real(8) :: crud_h(n_data)
    real(8) :: crud_thick(n_data)
    data crud_h / &
    0.60000, 0.62200, 0.64400, 0.66600, 0.68800, 0.71000, 0.73200, 0.75400, 0.77600, &
      & 0.79800, 0.82000, 0.84200, 0.86400, 0.88600, 0.90800, 0.93000, 0.96500/
    data crud_thick/ &
    0.00000, 0.13562, 0.27125, 0.40687, 0.62387, 0.81374, 1.08499, 1.27486, 1.38336, &
      & 1.43761, 1.49186, 1.49186, 1.46474, 1.35624, 1.08499, 0.67812, 0.00000/


    real(8) :: avg, nn, max1
    real(8) :: avgbu,hzsum,tmp_drf
    integer(4) :: crudz_bot,crudz_top
    integer(4) :: i1,i2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_Cruddr] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_Cruddr] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    cruddr(:,:)=0.0d0
    nn=0d0
    avg=0d0
    max1=0d0
    if(.not.opt_readcrud) then
        centz=0d0
        centz(izfuelbot)=h_z(izfuelbot)*0.5d0
        do Iz = IzFuelBot+1, IzFuelTop
            centz(iz)=centz(iz-1)+h_z(iz-1)*0.5d0+h_z(iz)*0.5d0
        enddo
        centz(:)=centz(:)/sum(h_z(izfuelbot:izfueltop))

        crud_factor=0d0
        do Iz = IzFuelBot, IzFuelTop
            do i=2,n_data
                if( crud_h(i-1) <= centz(iz) .and. centz(iz) < crud_h(i) ) then
                    frac1=centz(iz)-crud_h(i-1)
                    frac2=crud_h(i)-centz(iz)
                    tmp1=1d0/(frac1+frac2)
                    frac1=frac1*tmp1
                    frac2=frac2*tmp1
                    crud_factor(iz)=crud_thick(i-1)*frac2 + crud_thick(i)*frac1
                    exit
                endif
            enddo
        enddo

        centz(:)=centz(:)*sum(h_z(izfuelbot:izfueltop))

        crudz_bot=0
        crudz_top=0
        do Iz = IzFuelBot, IzFuelTop
           if(crud_factor(iz) >0.0d0) then
              crudz_bot=iz
              exit
           endif
        enddo
        do Iz = IzFuelTop,izfuelbot,-1
           if(crud_factor(iz) >0.0d0) then
              crudz_top=iz
              exit
           endif
        enddo

        tmp_drf=crud_drf
#ifdef hj_training
        if(opt_crud_trset) then
           do i = 1,crud_nfa
              ixy_1N = crud_ifa(i)
              avgbu = bogus_bu(i)
              do m = 1, 4
                 if (i_1nto4n(ixy_1n, m) == 0) then
                     cycle
                 else
                    ixy = i_1nto4n(ixy_1n, m)
                    do Iz = crudz_bot, crudz_top
                       if (opt_crud_model==0) then
                          cruddr(ixy,iz)=(crud_factor(iz)*avgbu)*(crud_factor(iz)*avgbu)*tmp_drf
                       elseif (opt_crud_model==1) then
                          cruddr(ixy,iz)= SNB_mass(Ixy,Iz) * boiling_factor
                       endif
                       if(cruddr(ixy,iz)>0d0) then
                          nn=nn+1d0
                          avg=avg+cruddr(ixy,iz)
                          max1=max(max1,cruddr(ixy,iz))
                       endif
                    enddo
                 endif
              enddo
           enddo
        else
#endif
        do ixy=1,nxy
            avgbu=0d0
            hzsum=0d0
            !do Iz = IzFuelBot, IzFuelTop
            do Iz = crudz_bot, crudz_top
                avgbu=avgbu+Ind_Bu(ixy,iz)*h_z(iz)
                hzsum=hzsum+h_z(iz)
            enddo
            avgbu=avgbu/hzsum
            !do Iz = IzFuelBot, IzFuelTop
            do Iz = crudz_bot, crudz_top
                if (opt_crud_model==0) then
                   cruddr(ixy,iz)=(crud_factor(iz)*avgbu)*(crud_factor(iz)*avgbu)*tmp_drf
                elseif (opt_crud_model==1) then
                   cruddr(ixy,iz)= SNB_mass(Ixy,Iz) * boiling_factor
                endif
                if(cruddr(ixy,iz)>0d0) then
                    nn=nn+1d0
                    avg=avg+cruddr(ixy,iz)
                    max1=max(max1,cruddr(ixy,iz))
                endif
            enddo
        enddo
#ifdef hj_training
        endif
#endif

    elseif(opt_readcrud .and. opt_crudboa) then
        i1=1
        do i=1,n_day_boa-1
            if(day<=day_boa(i)) then
                exit
            endif
            i1=i
        enddo
        i2=n_day_boa
        do i=n_day_boa,2,-1
            if(day>day_boa(i)) then
                exit
            endif
            i2=i
        enddo
        frac1=(day_boa(i2)-day)/(day_boa(i2)-day_boa(i1))
        frac2=1d0-frac1
        cruddr(:,:)=boa_dr(:,:,i1)*frac1+boa_dr(:,:,i2)*frac2
        do iz=1,nz
            do ixy=1,nxy
                cruddr(ixy,iz)=max(0.0d0,cruddr(ixy,iz))
            enddo
        enddo

        !if_exist=.false.
        !do i=n_day_boa,1,-1
        !    if(day>day_boa(i)) then
        !        if_exist=.true.
        !        exit
        !    endif
        !enddo
        !if(if_exist) then
        !    cruddr(:,:)=boa_dr(:,:,i)
        !else
        !    cruddr(:,:)=0d0
        !endif
    elseif(opt_readcrud) then
        do ixy_1n=1,nxy_1n
            do iz = 1, nz
                do m = 1, 4
                    if (i_1nto4n(ixy_1n, m) == 0) then
                        cycle
                    else
                        cruddr(i_1nto4n(ixy_1n, m),iz)=cruddr_table(ixy_1n,iz,i_bu)
                    endif

                    if(cruddr(i_1nto4n(ixy_1n, m),iz)>0d0) then
                        nn=nn+1d0
                        avg=avg+cruddr(I_1Nto4N(Ixy_1N, m),iz)
                        max1=max(max1,cruddr(I_1Nto4N(Ixy_1N, m),iz))
                    endif

                enddo
            enddo
        enddo
    endif
    if (comm%if_master) write(*,'(a,2es13.5)') '- AVG/Max CRUD dr  =',avg/max(1d0,nn),max1

    return
    end subroutine Get_Cruddr

    subroutine Get_Crud_DTF
    use inc_Crud, only: cruddtf,cruddr
    use inc_geometry, only: nxy,nz
    use inc_th, only: qwat       ! heat flux [W/m2]
    implicit none
    real(8), parameter :: k_crud = 1.5  ! [W/m/K]
    real(8) :: avg, nn, max1
    integer(4) :: iz,ixy

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_Crud_DTF] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_Crud_DTF] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    cruddtf=0.0d0
!    if(.not.opt_readcrud) then
     nn=0d0
     avg=0d0
     max1=0d0
     do iz=1,nz
         do ixy=1,nxy
             if(cruddr(ixy,iz)>0.0d0) then
                 !cruddtf(ixy,iz)= LHS / ( pi * R_clad * 2d0 * 100.0d0 * 0.015d0 / ( cruddr(ixy,iz) * 1.0D-4 ))
                 cruddtf(ixy,iz) = qwat(ixy,iz) * cruddr(ixy,iz)*1.d-6 / k_crud
                 if(cruddtf(ixy,iz)>0d0) then
                     nn=nn+1d0
                     avg=avg+cruddtf(ixy,iz)
                     max1=max(max1,cruddtf(ixy,iz))
                 endif
             endif
         enddo
     enddo
     !write(*,*) 'Zero cruddtf'
     write(*,'(a,2es13.5)') '- AVG/Max CRUD dTF =',avg/max(1d0,nn),max1
!    elseif(opt_readcrud) then
!        do iz=1,nz
!            do ixy=1,nxy
!                if(cruddr(ixy,iz)>0.0d0) then
!                    cruddtf(ixy,iz)=crtf*1.0d-4*cruddr(ixy,iz)*cruddr(ixy,iz)
!                endif
!            enddo
!        enddo
!    endif
    return
    end subroutine Get_Crud_DTF


      subroutine DEPB10
      use inc_parallel, only: comm
      use Inc_3D
      use Inc_Depletion
      use Inc_INP
      use Inc_maXS
      use Inc_miXS
      use Inc_Nuclide
      use Inc_Option
      use Inc_Time
      use Inc_TH
      use Inc_XS_File
      implicit none
      integer :: Ixy, Iz
      real(8) :: ABSRATE, ASFSUM
      real(8) :: TOTAREA, AREA
      real(8) :: FF, TF
      real(8) :: FASTABS, THABS
      real(8) :: AB10, AB11
      real(8) :: B11DEP
      real(8) :: RRCSB10, RRCSB11
      real(8) :: BORINB11
      real(8) :: WTRB10, WTRB11
      real(8) :: B10, B11
      real(8) :: DELB10B, DELB11B, DELB10M, DELB101, DELB102
      real(8) :: DELB10, DELB11

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [DELB10] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [DEPB10] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.flag_b10dep) return
      if (flag_pc) return ! do not change when corrector

      ! Calculate core total absorption rate if RCSVOL > 0
      if (RCSVOL > 0d0 .AND. Nz > 1) then
         if (dBU > 0d0 .AND. Hour > 0d0) then

            ! initialize
            ABSRATE = 0d0
            TOTAREA = 0d0

            do Ixy = 1, Nxy
               do Iz = 1, Nz

                  ! unpack PBU from BUNOD

                  ! extract fast and thermal flux
                  FF = Flux(Ixy,Iz,1)
                  TF = Flux(Ixy,Iz,2)

                  ! get the desired data from the big sifmic array
                  !FASTABS = miXS_a_Gd55(Ixy,Iz,1)
                  !THABS   = miXS_a_Gd55(Ixy,Iz,2)
                  FASTABS = miXS_a_B10(Ixy,Iz,1)
                  THABS   = miXS_a_B10(Ixy,Iz,2)
                  !write(878,'(a,2i5,4e13.5)') "testtesttest", Ixy,Iz, miXS_a_B10(Ixy,Iz,1) , miXS_a_Gd55(Ixy,Iz,1) , miXS_a_B10(Ixy,Iz,2) , miXS_a_Gd55(Ixy,Iz,2)
!                  if (FASTABS <= 1e-30) then
!                     FASTABS = tmp1
!                  endif
!                  if (THABS <= 1e-30) then
!                     THABS = tmp2
!                  endif
                  ! calculate node area (really volume)
                  if (FASTABS <= 1e-30 .OR. THABS <= 1e-30) then
                     AREA = 0d0
                  else
                     AREA = NodeVolume(Ixy,Iz)
                  endif
!                  AREA = NodeVolume(Ixy,Iz)

                  ! Cumulative absorption rate is area times
                  ! (phi_f * fast B abs + phi_th * thermal B abs * gfactor)
                  ABSRATE = ABSRATE + (FF*FASTABS + &
                                       TF*THABS*GFACT)*AREA
                  TOTAREA = TOTAREA + AREA

               enddo
            enddo
            ! Calculate average absorption rate
            ABSRATE = ABSRATE / TOTAREA
            ASFSUM = ABSRATE * dT
            ! AB10/AB11 = atomic weight of B10/B11
            AB10 = 10.013
            AB11 = 11.009

            ! B10DEP = atom % of B10 in RCS
            ! calculate the atom percent of B11
            ! B10DEP should be from other
            B11DEP = 100.0d0 - B10DEP

            ! calculate the weight % of B10/B11 in coolant
            ! using atom % and the atomic weights
            RRCSB10 = (B10DEP / 100.0d0) * AB10 &
                    / ((B10DEP/100.0d0)*AB10 + (B11DEP/100.0d0)*AB11)
            RRCSB11 = (B11DEP / 100.0d0) * AB11 &
                    / ((B10DEP/100.0d0)*AB10 + (B11DEP/100.0d0)*AB11)

            ! calculate the weight % of B10/B11 in borated water
            ! added to RCS using atom % and atomic weights
            ! WTRB10/WTRB11 = weight % of B10/B11 in the added borated water
            BORINB11 = 100.0d0 - BORINB10
            WTRB10 = (BORINB10 / 100.0d0) * AB10 &
                   / ((BORINB10/100.0d0)*AB10 + (BORINB11/100.0d0)*AB11)
            WTRB11 = (BORINB11 / 100.0d0) * AB11 &
                   / ((BORINB10/100.0d0)*AB10 + (BORINB11/100.0d0)*AB11)

            ! B10/B11 = amount of B10/B11 in coolant at BOS in ppm
            if (PPM > 0d0) then
               B10 = RRCSB10 * PPM
               B11 = RRCSB11 * PPM
            else
               B10 = RRCSB10
               B11 = RRCSB11
            endif

            ! DELB10 = change in the amount of B10 in RCS in ppm
            ! DELB10B = change in the amount of B10 in RCS in ppm
            !           due to borated water addition
            DELB10B = (WTRB10 * BORINFLOW * BORINCONC) &
                    / max(1d-50,(AVGRHO * RCSVOL))

            ! DELB10M = change in the amount of B10 in RCS in ppm
            !           due to makeup water
      ! (i.e. water flow out of RCS replaced by unborated water flow in)
            DELB10M = (B10 * MAKEFLOW) &
                    / max(1d-50,(AVGRHO * RCSVOL))

            ! DELB101 = total change in B10 due to change in flow in/out
            DELB101 = DELB10B - DELB10M

            ! if the DELTIME is long, some of the B10 in the makeup flow
            ! was depleted. As an estimate of this, DELB101 is depleted
            ! over the entire deltime as if all the flow was at BOS.
            ! This value, DELB102 is then averaged with DELB101 to
            ! estimate the actual B10 change due to flow.
            DELB102 = DELB101 * EXPCK(-ASFSUM*CORVOL/RCSVOL)
            DELB10  = (DELB101 + DELB102) / 2.0d0

            ! Finally the depletion of B10 which was in the coolant is
            ! added to the depletion contribution from the flow above
            DELB10 = DELB10 - B10*(1.0d0 - EXPCK(-ASFSUM*CORVOL/RCSVOL))
            ! Calculate the change in amount of B11 in RCS in ppm
            DELB11B = (WTRB11 * BORINFLOW * BORINCONC ) &
                    / max(1d-50,(AVGRHO * RCSVOL))
            DELB10M = (B11 * MAKEFLOW) / max(1d-50,(AVGRHO * RCSVOL))
            DELB11  = DELB11B - DELB10M

            ! Calculate the change in B10 and B11 in ppm
            B10 = B10 + DELB10
            B11 = B11 + DELB11

            ! Calculate B10 in atomic ratio for storage and ppm calc.
            B10DEP = 100.0d0 * AB11 *B10 &
                   / max(1d-50,(AB11 * B10 + AB10 * B11))

         endif

         PPMB10 = PPM * (B10PCT/max(1d-50,B10DEP))

      else

         ! if RCSVOL is still zero, just set adjusted ppm to nominal ppm
         PPMB10 = PPM

      endif
      if (comm%if_master) write(*,*) '--------------------------------------'
      if (comm%if_master) write(*,'(A,F9.2,A6,F9.2,A4)') &
         & 'Before B-10 Depletion : ', PPM*(B10PCT/max(1d-50,BORINB10)), '(PPM)', BORINB10, '(%)'
      if (comm%if_master) write(*,'(A,F9.2,A6,F9.2,A4)') &
         & 'After  B-10 Depletion : ', PPMB10, '(PPM)', B10DEP, '(%)'
      if (comm%if_master) write(*,*) '--------------------------------------'

      return
      end subroutine DEPB10

      FUNCTION EXPCK(X)

      implicit none

      real(8) :: EXPCK
      real(8) :: X

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [EXPCK] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [EXPCK] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (X < -5.0D2) then
         EXPCK = 0d0
      elseif (X > 5.0D2) then
         call print_msg(0,'ERROR OCCUR in EXPCK')
      else
         EXPCK = EXP(X)
      endif

      return
      end FUNCTION EXPCK


      subroutine DEPCR
      use Inc_3D
      use Inc_Depletion
      use Inc_INP
      use Inc_maXS
      use Inc_miXS
      use Inc_Nuclide
      use Inc_Option
      use Inc_Time
      use Inc_TH
      use Inc_XS_File
      use Inc_CR
      implicit none
      integer :: Ixy, Iz, I_CR
      real(8) :: ABSRATE
      real(8) :: FF, TF
      real(8) :: FASTABS, THABS

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [DEPCR] in Mod_Depletion'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [DEPCR] in Mod_Depletion'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      if (.not.flag_crdep) return
      if (flag_pc) return ! do not change when corrector

      CR_mat_frac = 1d0

      if (dBU>1d-10.and.dT>1d-10) then
         ABSRATE = 0d0
         do Ixy = 1, Nxy
            I_CR = I_CR_4N(Ixy)
            if (I_CR /= 0) then
               do Iz = 1, Nz
                  ! extract fast and thermal flux
                  FF = Flux(Ixy,Iz,1)
                  TF = Flux(ixy,Iz,2)

                  ! get the desired date from the big sinfmic array
                  FASTABS = miXS_a_B10(Ixy,Iz,1)
                  THABS   = miXS_a_B10(Ixy,Iz,2)

                  ! Cumulative absorption rate is area times
                  ! (phi_f * fast B abs + phi_th * thermal B abs)
                  ABSRATE = FF*FASTABS + TF*THABS

                  CR_mat_frac(Ixy,Iz) = EXP(-ABSRATE*dT)
               enddo
            endif
         enddo
      endif

      return

      end subroutine DEPCR


      end module Mod_Depletion


#endif
