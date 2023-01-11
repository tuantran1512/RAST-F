      MODULE Mod_GetSome

      USE Inc_Constant
      USE Inc_Geometry
      USE Inc_RP
      USE Mod_Alloc
#ifdef js_mpi
      use inc_parallel, only: comm
#endif

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

      CONTAINS


      SUBROUTINE Get_Avg(Avg, Input, Flag)
      USE Inc_INP
      use mod_charedit, only: print_msg

      IMPLICIT NONE
      REAL(8), DIMENSION(:, :), INTENT(IN) :: Input
      REAL(8), INTENT(OUT) :: Avg
      REAL(8) :: Tot_FlagVol, Tot_Sum
      INTEGER, INTENT(IN) :: Flag
      INTEGER :: Ixy, Iz, Ibot, Itop
      INTEGER :: Check_Even

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_Avg] in Mod_GetSome'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_Avg] in Mod_GetSome'
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

      Ixy  = 0
      Iz   = 0
      Ibot = 0
      Itop = 0

      Check_Even = 0

      Avg         = D0
      Tot_FlagVol = D0
      Tot_Sum     = D0

      IF (MOD((IzFuelBot + IzFuelTop), 2) == 1) THEN
         Check_Even = 1
      END IF

      IF (Flag == 0) THEN
         Ibot = IzFuelBot
         Itop = IzFuelTop
      ELSE IF (Flag == 1) THEN
         Ibot = IzFuelBot
         Itop = (IzFuelBot + IzFuelTop)/2
      ELSE IF (Flag == 2) THEN
         Ibot = (IzFuelBot + IzFuelTop)/2 + Check_Even
         Itop = IzFuelTop
      ELSE
        call print_msg(3,'call Get_Avg', Flag)
      END IF

      DO Ixy = 1, Nxy
         IF ( I_FARF_1N(I_4Nto1N(Ixy)) == 1 ) THEN
            CYCLE
         END IF

         DO Iz = Ibot, Itop

            Tot_FlagVol = Tot_FlagVol +  &
               NodeVolume(Ixy, Iz)

            Tot_Sum = Tot_Sum +  &
               Input(Ixy, Iz) * NodeVolume(Ixy, Iz)
         END DO
      END DO

      Avg = Tot_Sum / Tot_FlagVol
      RETURN
      END SUBROUTINE Get_Avg


      SUBROUTINE Get_Source(Src, maXS, Flux)

      USE Inc_Option , ONLY: N_Group

      IMPLICIT NONE

      REAL(8), DIMENSION(:, :), INTENT(OUT) :: Src
      REAL(8), DIMENSION(:, :, :), INTENT(IN) :: maXS
      REAL(8), DIMENSION(0:Nxy+1, Nz, N_Group), INTENT(IN) :: Flux
      INTEGER :: Ixy, Ixy_FA, Iz, Ig

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_Source] in Mod_GetSome'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_Source] in Mod_GetSome'
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

      Src = D0

      DO Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)

         DO Iz = IzFuelBot, IzFuelTop
            DO Ig = 1, N_Group
               Src(Ixy, Iz) = Src(Ixy, Iz) +  &
                  maXS(Ixy, Iz, Ig) * Flux(Ixy, Iz, Ig) * NodeVolume(Ixy, Iz)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE Get_Source

#ifdef siarhei_delete 
      subroutine get_source_rev(src,maXS,flux)
      use inc_option, only: n_group
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: ixy2iproc, ixyip2ixy
      use mod_parallel, only: allreduce
#endif
      implicit none
      real(8), dimension(:,:), intent(out) :: src
      real(8), dimension(:,:,:), intent(in) :: maXS
      real(8), dimension(0:nxy+1,nz,n_group), intent(in) :: flux
      integer :: ixy, ixy_fa, iz, ig
      integer :: l0

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [get_source_rev] in Mod_GetSome'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      Src=0d0
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         l0=ixy
#ifdef js_mpi
         if (comm%usempi) then
            if (iproc/=ixy2iproc(ixy)) cycle
            l0=ixyip2ixy(ixy,iproc)
         endif
#endif
         do iz=izfuelbot,izfueltop
            do ig=1,n_group
               src(ixy,iz)=src(ixy,iz)+maXS(ig,l0,iz)*flux(ixy,iz,ig)
            enddo
         enddo
      enddo

#ifdef js_mpi
      if (comm%usempi) then
         call allreduce(src,nxy,nz)
      endif
#endif

      return
      end subroutine get_source_rev
#endif 


      subroutine get_asi(ASI,input)
      implicit none
      real(8), intent(in) :: input(:,:)
      real(8), intent(out) :: ASI
      real(8) :: Avg_Bot, Avg_Top
      integer :: ixy, iz, iztb
      real(8) :: h_active, h_botref
      real(8) :: hz_half, h_cumul
      real(8) :: wgt
      real(8) :: rtmp, vtmp, rsum, vsum

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [get_asi] in Mod_GetSome'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_asi] in Mod_GetSome'
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

      if (flag_axial_mirror) then
         call Get_Avg( Avg_Bot, Input, 1)
         call Get_Avg( Avg_Top, Input, 2)
         if (abs(Avg_Bot+Avg_Top)<1d-10) then
            ASI=0d0
         else
            ASI=(Avg_Bot-Avg_Top)/(Avg_Bot+Avg_Top)
         endif
      else
         ! calculate half plane point
         h_active=0d0
         do iz=izfuelbot,izfueltop
            h_active=h_active+gridsize_z(iz)
         enddo
         h_botref=0d0
         do iz=1,izfuelbot-1
            h_botref=h_botref+gridsize_z(iz)
         enddo
         hz_half=h_active/2d0

         ! cal_bot quantity
         rsum=0d0
         vsum=0d0
         h_cumul=0d0
         do iz=izfuelbot,izfueltop
            rtmp=0d0
            vtmp=0d0
            do ixy=1,nxy
               if (i_farf_1n(i_4nto1n(ixy))==1) cycle
               rtmp=rtmp+input(ixy,iz)*nodevolume(ixy,iz)
               vtmp=vtmp+nodevolume(ixy,iz)
            enddo
            h_cumul=h_cumul+gridsize_z(iz)
            if (h_cumul<=hz_half) then
               rsum=rsum+rtmp
               vsum=vsum+vtmp
            else
               wgt=(h_cumul-hz_half)/gridsize_z(iz)
               iztb=iz
               if (abs(wgt-1d0)>1d-10) then
                  rsum=rsum+rtmp*wgt
                  vsum=vsum+vtmp*wgt
               endif
               exit
            endif
         enddo
         if (vsum>1d-10) avg_bot=rsum/vsum

         ! cal_top quantity
         rsum=0d0
         vsum=0d0
         do iz=iztb,izfueltop
            rtmp=0d0
            vtmp=0d0
            do ixy=1,nxy
               if (i_farf_1n(i_4nto1n(ixy))==1) cycle
               rtmp=rtmp+input(ixy,iz)*nodevolume(ixy,iz)
               vtmp=vtmp+nodevolume(ixy,iz)
            enddo
            if (iz>iztb) wgt=1d0
            rsum=rsum+rtmp*wgt
            vsum=vsum+vtmp*wgt
         enddo
         if (vsum>1d-10) avg_top=rsum/vsum

         if (abs(avg_bot+avg_top)<1d-10) then
            ASI=0d0
         else
            ASI=(avg_bot-avg_top)/(avg_bot+avg_top)
         endif
      endif

      return
      end subroutine get_asi


      SUBROUTINE Get_FA_ASI

      use Inc_3D, only: power
      use Inc_xyz, only: FA_ASI

      IMPLICIT NONE

      REAL(8) :: FoldFactor
      REAL(8) :: Sum_Real
      INTEGER :: Ixy_1N, Iz, Ix_1N, Iy_1N
      INTEGER :: FoldCount, m
      real(8),allocatable :: power_1N(:,:,:)
      integer(4) :: iz_up_bot,iz_up_top,iz_down_bot,iz_down_top,check_even
      real(8) :: top_pow,bot_pow

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_FA_ASI] in Mod_GetSome'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_FA_ASI] in Mod_GetSome'
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

      call alloc( power_1n,nx_1n,ny_1n,nz)
      FoldFactor    = D0
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
                     Sum_Real = Sum_Real + Power (I_1Nto4N(Ixy_1N, m), Iz)
                  END IF
               END DO
               FoldFactor = D4/(D4 - REAL(FoldCount, 8))
               power_1n(Ix_1N, Iy_1N, Iz) = FoldFactor*(Sum_Real/D4)
            END DO
         END DO
      END DO

      IF (MOD((IzFuelBot + IzFuelTop), 2) == 1) THEN
         Check_Even=1
      else
         check_Even=0
      END IF

      iz_down_bot= IzFuelBot
      iz_down_top= (IzFuelBot + IzFuelTop)/2
      iz_up_bot  = (IzFuelBot + IzFuelTop)/2 + Check_Even
      iz_up_top  = IzFuelTop

      DO Iy_1N = 1, Ny_1N
         DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
            top_pow=0.0d0
            bot_pow=0.0d0
            do iz=iz_down_bot,iz_down_top
               bot_pow=bot_pow+power_1n(ix_1n,iy_1n,iz)*h_z(iz)
            enddo
            do iz=iz_up_bot,iz_up_top
               top_pow=top_pow+power_1n(ix_1n,iy_1n,iz)*h_z(iz)
            enddo
            FA_ASI(ix_1n,iy_1n)=(bot_pow-top_pow)/max(1d-30,(bot_pow+top_pow))
         END DO
      END DO

      RETURN
      END SUBROUTINE Get_FA_ASI


!!!#ifdef siarhei_delete 
      SUBROUTINE Get_CR_Info

      USE Inc_CR

      IMPLICIT NONE

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_CR_Info] in Mod_GetSome'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      Speed_CR = inch
      OverLap  = 90.0D0*inch

      CALL Alloc( PDIL , N_CR )
      CALL Alloc( PDWL , N_CR )

      RETURN
      END SUBROUTINE Get_CR_Info
!!!#endif 


      SUBROUTINE Get_Reg_VF( Ixy, Iz )

      USE Inc_CR
      USE Inc_XS_File, ONLY: Flag_Reg
      use inc_inp, only: flag_crw_cor, crw_cor
      use inc_inp, only: flag_crw_cor_each, crw_cor_each
      use inc_inp, only: w_cr, w_ucr
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: Ixy, Iz
      INTEGER :: I_CR
      REAL(8) :: AN_Bot, AN_Top
      real(8) :: w_rodded, w_unrodded

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_Reg_VF] in Mod_GetSome'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_Reg_VF] in Mod_GetSome'
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

      I_CR = I_CR_4N(Ixy)
      Reg_VF = D0
      Flag_Reg = .FALSE.

      IF ( Iz /= 1 ) THEN
         AN_Bot = SUM( h_z( 1 : (Iz - 1) ) )
      ELSE
         AN_Bot = D0
      END IF
      AN_Top = AN_Bot + h_z(Iz)

      IF ( I_CR == 0 ) THEN
         Reg_VF(1) = D1
         Flag_Reg(1) = .TRUE.
      ELSE
         IF ( CR_Bot(I_CR) >= AN_Top ) THEN
            Reg_VF(1) = D1
            Flag_Reg(1) = .TRUE.
         ELSE IF ( CR_Top(I_CR) <= AN_Bot ) THEN
            Reg_VF(5) = D1
            Reg_VF(5) = Reg_VF(5) * CR_mat_frac(Ixy,Iz)
            Flag_Reg(5) = .TRUE.
         ELSE IF ( Comp_CR(I_CR) == 1 ) THEN
            IF ( ( CR_Bot(I_CR)                       <= AN_Bot ) .AND.  &
                 ( CR_Bot(I_CR) + Length_CR_Tip(I_CR) >= AN_Top )       ) THEN
               Reg_VF(3) = D1
               Reg_VF(3) = Reg_VF(3) * CR_mat_frac(Ixy,Iz)
               Flag_Reg(3) = .TRUE.
            ELSE IF ( ( CR_Bot(I_CR) + Length_CR_Tip(I_CR) <= AN_Bot ) .AND.  &
                      ( CR_Top(I_CR)                       >= AN_Top )       ) THEN
               Reg_VF(2) = D1
               Reg_VF(2) = Reg_VF(2) * CR_mat_frac(Ixy,Iz)
               Flag_Reg(2) = .TRUE.
            ELSE IF ( ( CR_Bot(I_CR)                       >  AN_Bot ) .AND.  &
                      ( CR_Bot(I_CR) + Length_CR_Tip(I_CR) >= AN_Top )       ) THEN
               Reg_VF(1) = ( CR_Bot(I_CR) - AN_Bot ) / h_z(Iz)
               Reg_VF(3) = D1 - Reg_VF(1)
               Reg_VF(3) = Reg_VF(3) * CR_mat_frac(Ixy,Iz)
               Flag_Reg(1) = .TRUE.
               Flag_Reg(3) = .TRUE.
            ELSE IF ( ( CR_Bot(I_CR)                       <= AN_Bot ) .AND.  &
                      ( CR_Bot(I_CR) + Length_CR_Tip(I_CR) <  AN_Top ) .AND.  &
                      ( CR_Top(I_CR)                       >= AN_Top )       ) THEN
               Reg_VF(3) = ( ( CR_Bot(I_CR) + Length_CR_Tip(I_CR) ) - AN_Bot ) / h_z(Iz)
               Reg_VF(2) = D1 - Reg_VF(3)
               Reg_VF(2) = Reg_VF(2) * CR_mat_frac(Ixy,Iz)
               Reg_VF(3) = Reg_VF(3) * CR_mat_frac(Ixy,Iz)
               Flag_Reg(3) = .TRUE.
               Flag_Reg(2) = .TRUE.
            ELSE IF ( ( CR_Bot(I_CR) + Length_CR_Tip(I_CR) <= AN_Bot ) .AND.  &
                      ( CR_Top(I_CR)                       <  AN_Top )       ) THEN
               Reg_VF(2) = ( CR_Top(I_CR) - AN_Bot ) / h_z(Iz)
               Reg_VF(5) = D1 - Reg_VF(2)
               Reg_VF(5) = Reg_VF(5) * CR_mat_frac(Ixy,Iz)
               Reg_VF(2) = Reg_VF(2) * CR_mat_frac(Ixy,Iz)
               Flag_Reg(2) = .TRUE.
               Flag_Reg(5) = .TRUE.
            ELSE IF ( ( CR_Bot(I_CR)                       >  AN_Bot ) .AND.  &
                      ( CR_Bot(I_CR) + Length_CR_Tip(I_CR) <  AN_Top ) .AND.  &
                      ( CR_Top(I_CR)                       >= AN_Top )       ) THEN
               Reg_VF(1) = ( CR_Bot(I_CR) - AN_Bot ) / h_z(Iz)
               Reg_VF(3) = Length_CR_Tip(I_CR)       / h_z(Iz)
               Reg_VF(2) = D1 - Reg_VF(1) - Reg_VF(3)
               Reg_VF(2) = Reg_VF(2) * CR_mat_frac(Ixy,Iz)
               Reg_VF(3) = Reg_VF(3) * CR_mat_frac(Ixy,Iz)
               Flag_Reg(1) = .TRUE.
               Flag_Reg(3) = .TRUE.
               Flag_Reg(2) = .TRUE.
            ELSE IF ( ( CR_Bot(I_CR)                       <= AN_Bot ) .AND.  &
                      ( CR_Bot(I_CR) + Length_CR_Tip(I_CR) <  AN_Top ) .AND.  &
                      ( CR_Top(I_CR)                       <  AN_Top )       ) THEN
               Reg_VF(3) = ( ( CR_Bot(I_CR) + Length_CR_Tip(I_CR) ) - AN_Bot ) / h_z(Iz)
               Reg_VF(2) = ( Length_CR(I_CR) - Length_CR_Tip(I_CR) )           / h_z(Iz)
               Reg_VF(5) = D1 - Reg_VF(3) - Reg_VF(2)
               Reg_VF(2) = Reg_VF(2) * CR_mat_frac(Ixy,Iz)
               Reg_VF(3) = Reg_VF(3) * CR_mat_frac(Ixy,Iz)
               Reg_VF(5) = Reg_VF(5) * CR_mat_frac(Ixy,Iz)
               Flag_Reg(3) = .TRUE.
               Flag_Reg(2) = .TRUE.
               Flag_Reg(5) = .TRUE.
            ELSE IF ( ( CR_Bot(I_CR)                       >  AN_Bot ) .AND.  &
                      ( CR_Bot(I_CR) + Length_CR_Tip(I_CR) <  AN_Top ) .AND.  &
                      ( CR_Top(I_CR)                       <  AN_Top )       ) THEN
               Reg_VF(1) = ( CR_Bot(I_CR) - AN_Bot )                 / h_z(Iz)
               Reg_VF(3) = Length_CR_Tip(I_CR)                       / h_z(Iz)
               Reg_VF(2) = ( Length_CR(I_CR) - Length_CR_Tip(I_CR) ) / h_z(Iz)
               Reg_VF(5) = D1 - Reg_VF(1) - Reg_VF(3) - Reg_VF(2)

               Reg_VF(2) = Reg_VF(2) * CR_mat_frac(Ixy,Iz)
               Reg_VF(3) = Reg_VF(3) * CR_mat_frac(Ixy,Iz)
               Reg_VF(5) = Reg_VF(5) * CR_mat_frac(Ixy,Iz)

               Flag_Reg(1) = .TRUE.
               Flag_Reg(3) = .TRUE.
               Flag_Reg(2) = .TRUE.
               Flag_Reg(5) = .TRUE.
            ELSE
               PRINT *, "Error, Go to Debug"
               PRINT *, "Check Reg_VF (Full CR)"
               STOP
            END IF
         ELSE !Part strength
            IF ( ( CR_Bot(I_CR) <= AN_Bot ) .AND.  &
                 ( CR_Top(I_CR) >= AN_Top )       ) THEN
               Reg_VF(4) = D1
               Reg_VF(4) = Reg_VF(4) * CR_mat_frac(Ixy,Iz)
               Flag_Reg(4) = .TRUE.
            ELSE IF ( ( CR_Bot(I_CR) >  AN_Bot ) .AND.  &
                      ( CR_Top(I_CR) >= AN_Top )       ) THEN
               Reg_VF(1) = ( CR_Bot(I_CR) - AN_Bot ) / h_z(Iz)
               Reg_VF(4) = D1 - Reg_VF(1)
               Reg_VF(4) = Reg_VF(4) * CR_mat_frac(Ixy,Iz)
               Flag_Reg(1) = .TRUE.
               Flag_Reg(4) = .TRUE.
            ELSE IF ( ( CR_Bot(I_CR) <= AN_Bot ) .AND.  &
                      ( CR_Top(I_CR) <  AN_Top )       ) THEN
               Reg_VF(4) = ( CR_Top(I_CR) - AN_Bot ) / h_z(Iz)
               Reg_VF(5) = D1 - Reg_VF(4)
               Reg_VF(4) = Reg_VF(4) * CR_mat_frac(Ixy,Iz)
               Reg_VF(5) = Reg_VF(5) * CR_mat_frac(Ixy,Iz)
               Flag_Reg(4) = .TRUE.
               Flag_Reg(5) = .TRUE.

            ELSE IF ( ( CR_Bot(I_CR) >  AN_Bot ) .AND.  &
                      ( CR_Top(I_CR) <  AN_Top )       ) THEN
               Reg_VF(1) = ( CR_Bot(I_CR) - AN_Bot ) / h_z(Iz)
               Reg_VF(4) = Length_CR(I_CR)           / h_z(Iz)
               Reg_VF(5) = D1 - Reg_VF(1) - Reg_VF(4)
               Reg_VF(4) = Reg_VF(4) * CR_mat_frac(Ixy,Iz)
               Reg_VF(5) = Reg_VF(5) * CR_mat_frac(Ixy,Iz)
               Flag_Reg(1) = .TRUE.
               Flag_Reg(4) = .TRUE.
               Flag_Reg(5) = .TRUE.
            ELSE
               PRINT *, "Error, Go to Debug"
               PRINT *, "Check Reg_VF (Part CR)"
               STOP
            END IF
         END IF
      END IF

      if (flag_crw_cor) then
         if (flag_crw_cor_each) then
            if (i_cr==0) then
               w_cr=1d0
            else
               w_cr=crw_cor_each(i_cr)
            endif
         else
            w_cr=crw_cor
         endif
         w_rodded=reg_vf(2)+reg_vf(3)+reg_vf(4)
         w_unrodded=reg_vf(1)+reg_vf(5)
         if (w_unrodded>1d-10) then
            w_ucr=(1d0+w_cr*(w_unrodded-1d0))/w_unrodded
         else
            w_ucr=0d0
         endif
         reg_vf(1)=reg_vf(1)*w_ucr
         reg_vf(2)=reg_vf(2)*w_cr
         reg_vf(3)=reg_vf(3)*w_cr
         reg_vf(4)=reg_vf(4)*w_cr
         reg_vf(5)=reg_vf(5)*w_ucr
      endif

      RETURN
      END SUBROUTINE Get_Reg_VF
#ifdef tuan_fr
     SUBROUTINE update_maXS(a,b)
         use Inc_miXS !,     only: miXS_Hex
         USE Inc_RP, ONLY: AxialComp
         use Inc_File,     ONLY: NuNum
         use inc_xs_file, only: XSset_hex
         USE Inc_Option, ONLY: N_Group
! need to remove when kappa fission is correct
        use Inc_Fastdepletion, only: kappaf
        use Inc_Flag, only: flag_miXS
!

         IMPLICIT NONE
         INTEGER :: a,b, Im
         REAL(8) :: Sum_maXS(1:N_Group+6,1:N_Group)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [update_maXS] in Read_File'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         Sum_maXS(:,:) = 0

         if (flag_miXS .AND. (Nunum /= 0)) then
              
             DO Im = 1, NuNum
                   Sum_maXS(:,:) = Sum_maXS(:,:) + miXS_Hex(a,b,Im,:, :)* &
     		      N_FR(a,b,Im)

             END DO
     ! XS_residual need to modify in future
             XSset_Hex(a,b,:,:) = Sum_maXS(:,:) + XS_residual( AxialComp(I_LP_1N( a ), b),  :, :)
! tuan_fr need to remove after update the kappa
!            XSset_Hex(a,b,5,:) = Sum_maXS(5,:)
         end if

     END SUBROUTINE update_maXS

#ifdef tuan_fr_crm 
      SUBROUTINE Get_Reg_VF_HEX
      USE Inc_CR
      USE Inc_RP 
      USE Inc_XS_File 
      use Inc_Geometry, only: I_1Nto4N, I_4Nto1N
      use Inc_Option,        only: n_group
      use Inc_TPEN
      use Inc_maXS, only: maXS_R_3D, maXS_SCAT_3D, maXs_s_3D

      IMPLICIT NONE

      INTEGER :: I_CR, Ixy_1N, N_xs_tmp, N_LP_tmp, i, i2, Ig, Ig1
      INTEGER :: Ix, Isum, Iz, Iz1, Iz2, Iz3, Iz4, Iz5, Iz6
      REAL(8) :: hz1, hz2
      REAL(8) :: W_Vol1, W_Vol2, W_Vol3, W_Vol4, W_Vol5, W_Vol6
      REAL(8) ::  Ref_maXS_tmp(1:4*N_XS_Table,1:6*n_group)
      REAL(8) :: ref_smat_tmp(1:4*N_XS_Table,1:n_group,1:n_group)
      REAL(8) :: F1, F2, F3, F4, F5, F6
      INTEGER(4) :: m,mc,md,mf,ms
      INTEGER(4) :: ixy
      REAL(8) :: flux_r1(N_Group)
      REAL(8) :: flux_u1(N_Group)
      REAL(8) :: flux_r2(N_Group)
      REAL(8) :: flux_u2(N_Group)
      REAL(8) :: WT1(N_Group)
      REAL(8) :: WT2(N_Group)
      REAL(8) :: WT3(N_Group)
      REAL(8) :: WT4(N_Group)
      REAL(8) :: h_fuelbot
      REAL(8) :: h_fueltop


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_Reg_VF_HEX] in Mod_GetSome'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif
      call WATCH('Control SA move','START')

W_Vol1 = 0
W_Vol2 = 0
W_Vol3 = 0
W_Vol4 = 0
W_Vol5 = 0
W_Vol6 = 0


!----------------------------------------
! Determine location in XY directions for CR
      Nxy_surFA(:,:) = 0
      Nxy_CR(:) = 0
      
      DO I_CR = 1, N_CR
      !     
      
          Ix = 0
          Isum = 0
          DO WHILE (Ix < Ix_CR(I_CR) )
              Ix = Ix+1
              IF ( I_FARF_2D(Ix, Iy_CR(I_CR)) == 0 ) THEN
                  CYCLE
              END IF
              Isum = Isum + 1
          
          ENDDO
          Nxy_CR(I_CR) = sum(Nx_y(1:(Iy_CR(I_CR)-1))) + Isum
      
      ! Determine location of axial MATs for CR
          DO Iz = 1, Nz
              
              hz1 = SUM( h_z( 1 : (Iz - 1) ) )
              hz2 = SUM( h_z( 1 : (Iz) ) )
              
              F1 = (BotFol_Bot(I_CR) - hz1)
              F2 = (BotFol_Bot(I_CR) - hz2)
              F3 = (Abs_Bot(I_CR) - hz1)   
              F4 = (Abs_Bot(I_CR) - hz2)   
              F5 = (TopFol_Bot(I_CR) - hz1)   
              F6 = (TopFol_Bot(I_CR) - hz2)   
              
!              IF (hz2 < BotFol_Bot(I_CR)) THEN
              IF (1.0E-6 < F2) THEN
!                  write(*,*) I_comp(Nxy_CR(I_CR), Iz), Nxy_CR(I_CR), Iz, BotFol_Bot(I_CR), hz2
                  DO Ig = 1, N_group
                       XSset_hex(Nxy_CR(I_CR), Iz, 1,Ig) = Ref_maXS(I_comp(Nxy_CR(I_CR), Iz), Ig )
                       XSset_hex(Nxy_CR(I_CR), Iz, 2,Ig) = Ref_maXS(I_comp(Nxy_CR(I_CR), Iz), (1*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 3,Ig) = Ref_maXS(I_comp(Nxy_CR(I_CR), Iz), (2*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 4,Ig) = Ref_maXS(I_comp(Nxy_CR(I_CR), Iz), (3*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 5,Ig) = Ref_maXS(I_comp(Nxy_CR(I_CR), Iz), (4*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 6,Ig) = Ref_maXS(I_comp(Nxy_CR(I_CR), Iz), (5*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 6+ig,:) = Ref_smat( I_comp(Nxy_CR(I_CR), Iz), Ig, :)
                  ENDDO    
                  CYCLE
!              ELSE IF ((hz2 >= BotFol_Bot(I_CR)) .AND. (BotFol_Bot(I_CR)> hz1)) THEN
              ELSE IF ((1.0E-6 >= F2) .AND. (F1 > 1.0E-6) .AND. ((Abs_Bot(I_CR) - BotFol_Bot(I_CR)) > 1.0E-1 )) THEN
                  Iz1 = Iz - 1
                  Iz2 = Iz
                  W_Vol2 = (hz2 - BotFol_Bot(I_CR))/ (hz2 - hz1)
                  IF (W_Vol2 < 0.0001) THEN
                      W_Vol2 = 0D0
                  ELSE IF (W_Vol2 > 0.9999) THEN
                      W_Vol2 = 1D0
                  ENDIF
                  W_Vol1 = 1 - W_Vol2
!              ELSE IF ((hz1> BotFol_Bot(I_CR)) .AND. (hz2 <= Abs_Bot(I_CR))) THEN
              ELSE IF ((1.0E-6 > F1) .AND. (1.0E-6 <= F4)) THEN
                  DO Ig = 1, N_group
                       XSset_hex(Nxy_CR(I_CR), Iz, 1,Ig) = Ref_maXS( BotFol_Mat(I_CR), Ig )
                       XSset_hex(Nxy_CR(I_CR), Iz, 2,Ig) = Ref_maXS( BotFol_Mat(I_CR), (1*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 3,Ig) = Ref_maXS( BotFol_Mat(I_CR), (2*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 4,Ig) = Ref_maXS( BotFol_Mat(I_CR), (3*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 5,Ig) = Ref_maXS( BotFol_Mat(I_CR), (4*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 6,Ig) = Ref_maXS( BotFol_Mat(I_CR), (5*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 6+ig,:) = Ref_smat( BotFol_Mat(I_CR), Ig, :)
                  ENDDO
!              ELSE IF ((hz2 >= (Abs_Bot(I_CR) )) .AND. ((Abs_Bot(I_CR)  )> hz1)) THEN
              ELSE IF ((1.0E-6 > F4) .AND. (F3  > 1.0E-6)) THEN
                  Iz3 = Iz - 1
                  Iz4 = Iz     
                  W_Vol4 = (hz2 - Abs_Bot(I_CR))/ (hz2 - hz1)
                  IF (W_Vol4 < 0.0001) THEN
                      W_Vol4 = 0D0
                  ELSE IF (W_Vol4 > 0.9999) THEN
                      W_Vol4 = 1D0
                  ENDIF
                  W_Vol3 = 1 - W_Vol4
!              ELSE IF ((hz1> Abs_Bot(I_CR)) .AND. (hz2 <= TopFol_Bot(I_CR))) THEN
              ELSE IF ((1.0E-6 > F3) .AND. (1.0E-6 <= F6)) THEN
                  DO Ig = 1, N_group
                       XSset_hex(Nxy_CR(I_CR), Iz, 1,Ig) = Ref_maXS( Abs_Mat(I_CR), Ig )
                       XSset_hex(Nxy_CR(I_CR), Iz, 2,Ig) = Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 3,Ig) = Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 4,Ig) = Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 5,Ig) = Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 6,Ig) = Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 6+ig,:) = Ref_smat( Abs_Mat(I_CR), Ig, :)
                  ENDDO                      
!              ELSE IF ((hz2 > TopFol_Bot(I_CR)+1e-6 ) .AND. (TopFol_Bot(I_CR)+1e-6> hz1)) THEN
              ELSE IF ((1.0E-6 > F6 ) .AND. (F5 > 1.0E-6)) THEN
                  Iz5 = Iz - 1
                  Iz6 = Iz
                  W_Vol6 = (hz2 - TopFol_Bot(I_CR))/ (hz2 - hz1)
                  IF (W_Vol6 < 0.0001) THEN
                      W_Vol6 = 0D0
                  ELSE IF (W_Vol6 > 0.9999) THEN
                      W_Vol6 = 1D0
                  ENDIF
                  W_Vol5 = 1 - W_Vol6
!              ELSE IF (hz1>= TopFol_Bot(I_CR)+1e-6)  THEN
              ELSE IF (1.0E-6 >= F5)  THEN
                  DO Ig = 1, N_group
                       XSset_hex(Nxy_CR(I_CR), Iz, 1,Ig) = Ref_maXS( TopFol_Mat(I_CR), Ig )
                       XSset_hex(Nxy_CR(I_CR), Iz, 2,Ig) = Ref_maXS( TopFol_Mat(I_CR), (1*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 3,Ig) = Ref_maXS( TopFol_Mat(I_CR), (2*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 4,Ig) = Ref_maXS( TopFol_Mat(I_CR), (3*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 5,Ig) = Ref_maXS( TopFol_Mat(I_CR), (4*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 6,Ig) = Ref_maXS( TopFol_Mat(I_CR), (5*N_Group + Ig) )
                       XSset_hex(Nxy_CR(I_CR), Iz, 6+ig,:) = Ref_smat( TopFol_Mat(I_CR), Ig, :)
                  ENDDO
                  
              ENDIF 
              
          ENDDO

      ! for surrouding Fuel SA
              DO i = 1, 6
                  Isum = 0
                  Ix = 0
                  DO WHILE (Ix < Ix_surFA(I_CR,i) )
                      Ix = Ix+1
                      IF ( I_FARF_2D(Ix, Iy_surFA(I_CR,i)) == 0 ) THEN
                          CYCLE
                      END IF
                      Isum = Isum + 1
                  
                  ENDDO
                  Nxy_surFA(I_CR,i) = sum(Nx_y(1:(Iy_surFA(I_CR,i)-1))) + Isum
          ! update XS with SPH_F
                  IF (FLAG_SPH_F) THEN
                      IF (Nxy_surFA(I_CR,i) /= 0) THEN 
                          h_fuelbot = SUM( h_z( 1 : (IzFuelBot-1) ))
                          h_fueltop = SUM( h_z( 1 : (IzFuelTop-1) ))
                          IF ( ( (h_fuelbot - Abs_Bot(I_CR)) <= 1E-6 ) .AND. ((Abs_Bot(I_CR) -h_fueltop) <= 1E-6 )) THEN
                              DO Ig = 1, N_group
                                  XSset_Hex( Nxy_surFA(I_CR,i), Iz4,1, Ig ) = Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Ig)* &
                                     SPH_F(Abs_Mat(I_CR), AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Ig)*W_Vol4 + &
                                     Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Ig)*W_Vol3

                                  XSset_Hex( Nxy_surFA(I_CR,i), Iz4,2, Ig ) = Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),(1*N_Group + Ig))* &
                                     SPH_F(Abs_Mat(I_CR), AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Ig)*W_Vol4 + &
                                     Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),(1*N_Group + Ig))*W_Vol3

                                  XSset_Hex( Nxy_surFA(I_CR,i), Iz4,3, Ig ) = Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),(2*N_Group + Ig))* &
                                     SPH_F(Abs_Mat(I_CR), AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Ig)*W_Vol4 + &
                                     Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),(2*N_Group + Ig))*W_Vol3

                                  XSset_Hex( Nxy_surFA(I_CR,i), Iz4,4, Ig ) = Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),(3*N_Group + Ig))* &
                                     SPH_F(Abs_Mat(I_CR), AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Ig)*W_Vol4 + &
                                     Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),(3*N_Group + Ig))*W_Vol3

                                  XSset_Hex( Nxy_surFA(I_CR,i), Iz4,5, Ig ) = Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),(4*N_Group + Ig))* &
                                     SPH_F(Abs_Mat(I_CR), AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Ig)*W_Vol4 + &
                                     Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),(4*N_Group + Ig))*W_Vol3

                                  XSset_Hex( Nxy_surFA(I_CR,i), Iz4,6, Ig ) = Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),(5*N_Group + Ig))* &
                                     SPH_F(Abs_Mat(I_CR), AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Ig)*W_Vol4 + &
                                     Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),(5*N_Group + Ig))*W_Vol3

                                  XSset_Hex( Nxy_surFA(I_CR,i), Iz4,6 + Ig,: ) = Ref_smat(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Ig,:)* &
                                     SPH_F(Abs_Mat(I_CR), AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),:)*W_Vol4 + &
                                     Ref_smat(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Ig,:)*W_Vol3

!                                  XSset_Hex( Nxy_surFA(I_CR,i), Iz4,Ig, : ) = XSset_Hex( Nxy_surFA(I_CR,i), Iz4,Ig, : )* &
!                                     SPH_F(Abs_Mat(I_CR), AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),:)*W_Vol4 + &
!                                      XSset_Hex( Nxy_surFA(I_CR,i), Iz4,Ig, : )*W_Vol3
!                              XSset_Hex( Nxy_surFA(I_CR,i), Iz4,Ig, : ) = XSset_Hex( Nxy_surFA(I_CR,i), Iz4,Ig, : )*W_Vol4*1 +  XSset_Hex( Nxy_surFA(I_CR,i), Iz4,Ig, : )*W_Vol3
                              ENDDO 
                              i2 = 1
                              DO WHILE (i2 <= (IzFuelTop - Iz4))
                                  DO Ig = 1, N_group
                                      XSset_Hex( Nxy_surFA(I_CR,i), Iz4+i2 ,1, Ig ) =  Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),Ig)* &
                                         SPH_F(Abs_Mat(I_CR),AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),Ig)

                                      XSset_Hex( Nxy_surFA(I_CR,i), Iz4+i2 ,2, Ig ) =  Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),(1*N_Group + Ig))* &
                                         SPH_F(Abs_Mat(I_CR),AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),Ig)

                                      XSset_Hex( Nxy_surFA(I_CR,i), Iz4+i2 ,3, Ig ) =  Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),(2*N_Group + Ig))* &
                                         SPH_F(Abs_Mat(I_CR),AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),Ig)

                                      XSset_Hex( Nxy_surFA(I_CR,i), Iz4+i2 ,4, Ig ) =  Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),(3*N_Group + Ig))* &
                                         SPH_F(Abs_Mat(I_CR),AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),Ig)

                                      XSset_Hex( Nxy_surFA(I_CR,i), Iz4+i2 ,5, Ig ) =  Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),(4*N_Group + Ig))* &
                                         SPH_F(Abs_Mat(I_CR),AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),Ig)

                                      XSset_Hex( Nxy_surFA(I_CR,i), Iz4+i2 ,6, Ig ) =  Ref_maXS(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),(5*N_Group + Ig))* &
                                         SPH_F(Abs_Mat(I_CR),AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),Ig)

                                      XSset_Hex( Nxy_surFA(I_CR,i), Iz4+i2 ,6+Ig,: ) =  Ref_smat(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),Ig,:)* &
                                         SPH_F(Abs_Mat(I_CR),AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),:)

!                                      XSset_Hex( Nxy_surFA(I_CR,i), Iz4+i2 ,ig, : ) = XSset_Hex( Nxy_surFA(I_CR,i), Iz4+i2,Ig, : )* &
!                                          SPH_F(Abs_Mat(I_CR),AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4+i2 ),:)
 !                                 XSset_Hex( Nxy_surFA(I_CR,i), Iz4+i2 ,ig, : ) = XSset_Hex( Nxy_surFA(I_CR,i), Iz4+i2,Ig, : )*1
                                  ENDDO 
                                  i2 = i2 + 1
                              ENDDO
                          END IF
                      ENDIF 
                  ENDIF
              ENDDO
          

! flux weight

    flux_r1(:)=1D0
    flux_u1(:)=1D0
    flux_u2(:)=1D0
    flux_r2(:)=1D0

    WT1 = 1D0
    WT2 = 1D0
    WT3 = 1D0
    WT4 = 1D0

If (flag_decusping_update) then

    flux_r1(:) =(h_z(iz4)*fluxf(Nxy_CR(I_CR),Iz4,:) + h_z( Iz4+1 ) *fluxf(Nxy_CR(I_CR),Iz4+1,:))/SUM( h_z(  Iz4 :   (Iz4+1) ))
    flux_u1(:) =(h_z(iz4)*fluxf(Nxy_CR(I_CR),Iz4,:) + h_z( Iz4-1 ) *fluxf(Nxy_CR(I_CR),Iz4-1,:))/SUM( h_z( (Iz4-1) : Iz4 ))
    WT1 = flux_r1(:)*(SUM( h_z( 1 : (Iz4) ) )- Abs_Bot(I_CR))
    WT2 = flux_u1(:)*( Abs_Bot(I_CR)-SUM( h_z( 1 : (Iz3) ) ))
    if (Iz6 <   Nz) then
        flux_u2(:) =(h_z(iz6)*fluxf(Nxy_CR(I_CR),Iz6,:) + h_z( Iz6+1 ) *fluxf(Nxy_CR(I_CR),Iz6+1,:))/SUM( h_z(  Iz6 :   (Iz6+1) ))
        flux_r2(:) =(h_z(iz6)*fluxf(Nxy_CR(I_CR),Iz6,:) + h_z( Iz6-1 ) *fluxf(Nxy_CR(I_CR),Iz6-1,:))/SUM( h_z( (Iz6-1) : Iz6 ))
        WT3 = flux_u2(:)*(SUM( h_z( 1 : (Iz6) ) )- TopFol_Bot(I_CR))
        WT4 = flux_r2(:)*( TopFol_Bot(I_CR)-SUM( h_z( 1 : (Iz5) ) ))
    endif 

endif

IF (flag_SPH_CR .AND. (Nxy_surFA(I_CR,1)/=0)) THEN

    DO Iz = 1, Nz
        hz1 = SUM( h_z( 1 : (Iz - 1) ) )
        hz2 = SUM( h_z( 1 : (Iz) ) )
        F3 = (Abs_Bot(I_CR) - hz1)   
        F6 = (TopFol_Bot(I_CR) - hz2)   
              
        IF ((1.0E-6 > F3) .AND. (1.0E-6 <= F6)) THEN
            DO Ig = 1, N_group
! note: modify AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz ) in th future 
              XSset_hex(Nxy_CR(I_CR), Iz, 1,Ig) = Ref_maXS( Abs_Mat(I_CR), Ig )* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz ),Abs_Mat(I_CR),Ig)
              XSset_hex(Nxy_CR(I_CR), Iz, 2,Ig) = Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig) )* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz ),Abs_Mat(I_CR),Ig)
              XSset_hex(Nxy_CR(I_CR), Iz, 3,Ig) = Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig) )* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz ),Abs_Mat(I_CR),Ig)
              XSset_hex(Nxy_CR(I_CR), Iz, 4,Ig) = Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig) )* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz ),Abs_Mat(I_CR),Ig)
              XSset_hex(Nxy_CR(I_CR), Iz, 5,Ig) = Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig) )* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz ),Abs_Mat(I_CR),Ig)
              XSset_hex(Nxy_CR(I_CR), Iz, 6,Ig) = Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig) )* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz ),Abs_Mat(I_CR),Ig)
              XSset_hex(Nxy_CR(I_CR), Iz, 6+ig,:) = Ref_smat( Abs_Mat(I_CR), Ig, :)* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz ),Abs_Mat(I_CR),:)
            ENDDO 
        ENDIF
    ENDDO
ENDIF


      ! update XS with volume weight: try 2 way: 1. one big matrix; 2. 7
      ! Replace the XS of Bottom follower, Absorber, and Top follower
      ! in the selected region

      
      ! Interpolate XS using volume weight
If (flag_decusping_update) then
! comment for testing transient, need to remove for normal calculation
     DO Ig = 1, N_group
          IF ((Abs_Bot(I_CR) - BotFol_Bot(I_CR)) > 1.0E-1 ) THEN
              XSset_hex(Nxy_CR(I_CR), Iz2, 1,Ig) = Ref_maXS( BotFol_Mat(I_CR), Ig )* W_Vol2 + Ref_maXS( Rep_MAT, Ig )* W_Vol1
              XSset_hex(Nxy_CR(I_CR), Iz2, 2,Ig) = Ref_maXS( BotFol_Mat(I_CR), (1*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (1*N_Group + Ig)  )* W_Vol1
              XSset_hex(Nxy_CR(I_CR), Iz2, 3,Ig) = Ref_maXS( BotFol_Mat(I_CR), (2*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (2*N_Group + Ig)  )* W_Vol1
              XSset_hex(Nxy_CR(I_CR), Iz2, 4,Ig) = Ref_maXS( BotFol_Mat(I_CR), (3*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (3*N_Group + Ig)  )* W_Vol1
              XSset_hex(Nxy_CR(I_CR), Iz2, 5,Ig) = Ref_maXS( BotFol_Mat(I_CR), (4*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (4*N_Group + Ig)  )* W_Vol1
              XSset_hex(Nxy_CR(I_CR), Iz2, 6,Ig) = Ref_maXS( BotFol_Mat(I_CR), (5*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (5*N_Group + Ig)  )* W_Vol1
              XSset_hex(Nxy_CR(I_CR), Iz2, 6+ig,:) = Ref_smat( BotFol_Mat(I_CR), Ig, :)* W_Vol2 +  Ref_smat( Rep_MAT, Ig, :)* W_Vol1
          ENDIF

          IF(FLAG_SPH_CR) THEN
              
              XSset_hex(Nxy_CR(I_CR), Iz4, 1,Ig) = (Ref_maXS( Abs_Mat(I_CR), Ig )* WT1(Ig)* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Abs_Mat(I_CR),Ig) + &
                  Ref_maXS( BotFol_Mat(I_CR), Ig )* WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 2,Ig) = (Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig) )* WT1(Ig) * &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Abs_Mat(I_CR),Ig)+  &
                  Ref_maXS( BotFol_Mat(I_CR), (1*N_Group + Ig)  )*WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 3,Ig) = (Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig) )* WT1(Ig)* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Abs_Mat(I_CR),Ig) +  &
                  Ref_maXS( BotFol_Mat(I_CR), (2*N_Group + Ig)  )*WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 4,Ig) = (Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig) )* WT1(Ig)* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Abs_Mat(I_CR),Ig) +  &
                  Ref_maXS( BotFol_Mat(I_CR), (3*N_Group + Ig)  )*WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 5,Ig) = (Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig) )* WT1(Ig)* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Abs_Mat(I_CR),Ig) +  &
                  Ref_maXS( BotFol_Mat(I_CR), (4*N_Group + Ig)  )*WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 6,Ig) = (Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig) )* WT1(Ig)* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Abs_Mat(I_CR),Ig) +  &
                  Ref_maXS( BotFol_Mat(I_CR), (5*N_Group + Ig)  )*WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 6+ig,:) = (Ref_smat( Abs_Mat(I_CR), Ig, :)* WT1(Ig)* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz4 ),Abs_Mat(I_CR),:) + &
                  Ref_smat( BotFol_Mat(I_CR), Ig, :)* WT2(Ig))/(WT1(Ig)+WT2(Ig))
                                     
              XSset_hex(Nxy_CR(I_CR), Iz6, 1,Ig) = (Ref_maXS( TopFol_Mat(I_CR), Ig )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), Ig )* WT4(Ig)* &
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz6 ),Abs_Mat(I_CR),Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 2,Ig) = (Ref_maXS( TopFol_Mat(I_CR), (1*N_Group + Ig) )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig)  )* WT4(Ig)*&
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz6 ),Abs_Mat(I_CR),Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 3,Ig) = (Ref_maXS( TopFol_Mat(I_CR), (2*N_Group + Ig) )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig)  )* WT4(Ig)*&
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz6 ),Abs_Mat(I_CR),Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 4,Ig) = (Ref_maXS( TopFol_Mat(I_CR), (3*N_Group + Ig) )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig)  )* WT4(Ig)*&
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz6 ),Abs_Mat(I_CR),Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 5,Ig) = (Ref_maXS( TopFol_Mat(I_CR), (4*N_Group + Ig) )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig)  )* WT4(Ig)*&
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz6 ),Abs_Mat(I_CR),Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 6,Ig) = (Ref_maXS( TopFol_Mat(I_CR), (5*N_Group + Ig) )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig)  )* WT4(Ig)*&
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz6 ),Abs_Mat(I_CR),Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 6+ig,:) = (Ref_smat( TopFol_Mat(I_CR), Ig, :)* WT3(Ig) + &
                  Ref_smat( Abs_Mat(I_CR), Ig, :)* WT4(Ig)*&
                  SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,i)), Iz6 ),Abs_Mat(I_CR),:))/(WT3(Ig)+WT4(Ig))
              
          ELSE
              
              XSset_hex(Nxy_CR(I_CR), Iz4, 1,Ig) = (Ref_maXS( Abs_Mat(I_CR), Ig )* WT1(Ig) + &
                  Ref_maXS( BotFol_Mat(I_CR), Ig )* WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 2,Ig) = (Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig) )* WT1(Ig) +  &
                  Ref_maXS( BotFol_Mat(I_CR), (1*N_Group + Ig)  )*WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 3,Ig) = (Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig) )* WT1(Ig) +  &
                  Ref_maXS( BotFol_Mat(I_CR), (2*N_Group + Ig)  )*WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 4,Ig) = (Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig) )* WT1(Ig) +  &
                  Ref_maXS( BotFol_Mat(I_CR), (3*N_Group + Ig)  )*WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 5,Ig) = (Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig) )* WT1(Ig) +  &
                  Ref_maXS( BotFol_Mat(I_CR), (4*N_Group + Ig)  )*WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 6,Ig) = (Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig) )* WT1(Ig) +  &
                  Ref_maXS( BotFol_Mat(I_CR), (5*N_Group + Ig)  )*WT2(Ig))/(WT1(Ig)+WT2(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz4, 6+ig,:) = (Ref_smat( Abs_Mat(I_CR), Ig, :)* WT1(Ig) + &
                  Ref_smat( BotFol_Mat(I_CR), Ig, :)* WT2(Ig))/(WT1(Ig)+WT2(Ig))
                                     
              XSset_hex(Nxy_CR(I_CR), Iz6, 1,Ig) = (Ref_maXS( TopFol_Mat(I_CR), Ig )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), Ig )* WT4(Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 2,Ig) = (Ref_maXS( TopFol_Mat(I_CR), (1*N_Group + Ig) )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig)  )* WT4(Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 3,Ig) = (Ref_maXS( TopFol_Mat(I_CR), (2*N_Group + Ig) )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig)  )* WT4(Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 4,Ig) = (Ref_maXS( TopFol_Mat(I_CR), (3*N_Group + Ig) )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig)  )* WT4(Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 5,Ig) = (Ref_maXS( TopFol_Mat(I_CR), (4*N_Group + Ig) )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig)  )* WT4(Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 6,Ig) = (Ref_maXS( TopFol_Mat(I_CR), (5*N_Group + Ig) )* WT3(Ig) + &
                  Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig)  )* WT4(Ig))/(WT3(Ig)+WT4(Ig))
              XSset_hex(Nxy_CR(I_CR), Iz6, 6+ig,:) = (Ref_smat( TopFol_Mat(I_CR), Ig, :)* WT3(Ig) + &
                  Ref_smat( Abs_Mat(I_CR), Ig, :)* WT4(Ig))/(WT3(Ig)+WT4(Ig))
          ENDIF
          
          
     ENDDO

ELSE
    
    DO Ig = 1, N_group
         IF ((Abs_Bot(I_CR) - BotFol_Bot(I_CR)) > 1.0E-1 ) THEN
             XSset_hex(Nxy_CR(I_CR), Iz2, 1,Ig) = Ref_maXS( BotFol_Mat(I_CR), Ig )* W_Vol2 + Ref_maXS( Rep_MAT, Ig )* W_Vol1
             XSset_hex(Nxy_CR(I_CR), Iz2, 2,Ig) = Ref_maXS( BotFol_Mat(I_CR), (1*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (1*N_Group + Ig)  )* W_Vol1
             XSset_hex(Nxy_CR(I_CR), Iz2, 3,Ig) = Ref_maXS( BotFol_Mat(I_CR), (2*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (2*N_Group + Ig)  )* W_Vol1
             XSset_hex(Nxy_CR(I_CR), Iz2, 4,Ig) = Ref_maXS( BotFol_Mat(I_CR), (3*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (3*N_Group + Ig)  )* W_Vol1
             XSset_hex(Nxy_CR(I_CR), Iz2, 5,Ig) = Ref_maXS( BotFol_Mat(I_CR), (4*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (4*N_Group + Ig)  )* W_Vol1
             XSset_hex(Nxy_CR(I_CR), Iz2, 6,Ig) = Ref_maXS( BotFol_Mat(I_CR), (5*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (5*N_Group + Ig)  )* W_Vol1
             XSset_hex(Nxy_CR(I_CR), Iz2, 6+ig,:) = Ref_smat( BotFol_Mat(I_CR), Ig, :)* W_Vol2 +  Ref_smat( Rep_MAT, Ig, :)* W_Vol1
         ENDIF
         
         IF(FLAG_SPH_CR .AND. (Nxy_surFA(I_CR,1)/=0)) THEN

             XSset_hex(Nxy_CR(I_CR), Iz4, 1,Ig) = Ref_maXS( Abs_Mat(I_CR), Ig )* W_Vol4* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz4 ),Abs_Mat(I_CR),Ig) + &
                 Ref_maXS( BotFol_Mat(I_CR), Ig )* W_Vol3
             XSset_hex(Nxy_CR(I_CR), Iz4, 2,Ig) = Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig) )* W_Vol4* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz4 ),Abs_Mat(I_CR),Ig)+ &
                 Ref_maXS( BotFol_Mat(I_CR), (1*N_Group + Ig)  )* W_Vol3
             XSset_hex(Nxy_CR(I_CR), Iz4, 3,Ig) = Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig) )* W_Vol4* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz4 ),Abs_Mat(I_CR),Ig)+ &
                 Ref_maXS( BotFol_Mat(I_CR), (2*N_Group + Ig)  )* W_Vol3
             XSset_hex(Nxy_CR(I_CR), Iz4, 4,Ig) = Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig) )* W_Vol4* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz4 ),Abs_Mat(I_CR),Ig)+ &
                 Ref_maXS( BotFol_Mat(I_CR), (3*N_Group + Ig)  )* W_Vol3
             XSset_hex(Nxy_CR(I_CR), Iz4, 5,Ig) = Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig) )* W_Vol4* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz4 ),Abs_Mat(I_CR),Ig)+ &
                 Ref_maXS( BotFol_Mat(I_CR), (4*N_Group + Ig)  )* W_Vol3
             XSset_hex(Nxy_CR(I_CR), Iz4, 6,Ig) = Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig) )* W_Vol4* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz4 ),Abs_Mat(I_CR),Ig)+ &
                 Ref_maXS( BotFol_Mat(I_CR), (5*N_Group + Ig)  )* W_Vol3
             XSset_hex(Nxy_CR(I_CR), Iz4, 6+ig,:) = Ref_smat( Abs_Mat(I_CR), Ig, :)* W_Vol4* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz4 ),Abs_Mat(I_CR),:)+ &
                 Ref_smat( BotFol_Mat(I_CR), Ig, :)* W_Vol3
             
               
             XSset_hex(Nxy_CR(I_CR), Iz6, 1,Ig) = Ref_maXS( TopFol_Mat(I_CR), Ig )* W_Vol6 + &
                 Ref_maXS( Abs_Mat(I_CR), Ig )* W_Vol5* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)),Iz6-1 ),Abs_Mat(I_CR),Ig)
             XSset_hex(Nxy_CR(I_CR), Iz6, 2,Ig) = Ref_maXS( TopFol_Mat(I_CR), (1*N_Group + Ig) )* W_Vol6 + &
                 Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig)  )* W_Vol5* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz6-1  ),Abs_Mat(I_CR),Ig)
             XSset_hex(Nxy_CR(I_CR), Iz6, 3,Ig) = Ref_maXS( TopFol_Mat(I_CR), (2*N_Group + Ig) )* W_Vol6 + &
                 Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig)  )* W_Vol5* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz6-1  ),Abs_Mat(I_CR),Ig)
             XSset_hex(Nxy_CR(I_CR), Iz6, 4,Ig) = Ref_maXS( TopFol_Mat(I_CR), (3*N_Group + Ig) )* W_Vol6 + &
                 Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig)  )* W_Vol5* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz6-1  ),Abs_Mat(I_CR),Ig)
             XSset_hex(Nxy_CR(I_CR), Iz6, 5,Ig) = Ref_maXS( TopFol_Mat(I_CR), (4*N_Group + Ig) )* W_Vol6 + &
                 Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig)  )* W_Vol5* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz6-1  ),Abs_Mat(I_CR),Ig)
             XSset_hex(Nxy_CR(I_CR), Iz6, 6,Ig) = Ref_maXS( TopFol_Mat(I_CR), (5*N_Group + Ig) )* W_Vol6 + &
                 Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig)  )* W_Vol5* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz6-1  ),Abs_Mat(I_CR),Ig)
             XSset_hex(Nxy_CR(I_CR), Iz6, 6+ig,:) = Ref_smat( TopFol_Mat(I_CR), Ig, :)* W_Vol6 + &
                 Ref_smat( Abs_Mat(I_CR), Ig, :)* W_Vol5* &
                 SPH_F(AxialComp( I_LP_1N(Nxy_surFA(I_CR,1)), Iz6-1  ),Abs_Mat(I_CR),:)
        ELSE 
            
            XSset_hex(Nxy_CR(I_CR), Iz4, 1,Ig) = Ref_maXS( Abs_Mat(I_CR), Ig )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), Ig )* W_Vol3
            XSset_hex(Nxy_CR(I_CR), Iz4, 2,Ig) = Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig) )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), (1*N_Group + Ig)  )* W_Vol3
            XSset_hex(Nxy_CR(I_CR), Iz4, 3,Ig) = Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig) )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), (2*N_Group + Ig)  )* W_Vol3
            XSset_hex(Nxy_CR(I_CR), Iz4, 4,Ig) = Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig) )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), (3*N_Group + Ig)  )* W_Vol3
            XSset_hex(Nxy_CR(I_CR), Iz4, 5,Ig) = Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig) )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), (4*N_Group + Ig)  )* W_Vol3
            XSset_hex(Nxy_CR(I_CR), Iz4, 6,Ig) = Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig) )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), (5*N_Group + Ig)  )* W_Vol3
            XSset_hex(Nxy_CR(I_CR), Iz4, 6+ig,:) = Ref_smat( Abs_Mat(I_CR), Ig, :)* W_Vol4 + Ref_smat( BotFol_Mat(I_CR), Ig, :)* W_Vol3
                                   
            XSset_hex(Nxy_CR(I_CR), Iz6, 1,Ig) = Ref_maXS( TopFol_Mat(I_CR), Ig )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), Ig )* W_Vol5
            XSset_hex(Nxy_CR(I_CR), Iz6, 2,Ig) = Ref_maXS( TopFol_Mat(I_CR), (1*N_Group + Ig) )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig)  )* W_Vol5
            XSset_hex(Nxy_CR(I_CR), Iz6, 3,Ig) = Ref_maXS( TopFol_Mat(I_CR), (2*N_Group + Ig) )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig)  )* W_Vol5
            XSset_hex(Nxy_CR(I_CR), Iz6, 4,Ig) = Ref_maXS( TopFol_Mat(I_CR), (3*N_Group + Ig) )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig)  )* W_Vol5
            XSset_hex(Nxy_CR(I_CR), Iz6, 5,Ig) = Ref_maXS( TopFol_Mat(I_CR), (4*N_Group + Ig) )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig)  )* W_Vol5
            XSset_hex(Nxy_CR(I_CR), Iz6, 6,Ig) = Ref_maXS( TopFol_Mat(I_CR), (5*N_Group + Ig) )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig)  )* W_Vol5
            XSset_hex(Nxy_CR(I_CR), Iz6, 6+ig,:) = Ref_smat( TopFol_Mat(I_CR), Ig, :)* W_Vol6 + Ref_smat( Abs_Mat(I_CR), Ig, :)* W_Vol5
        ENDIF
        
    ENDDO
    
ENDIF

!            DO Ig = 1, N_group
!                 IF (BotFol_Bot(I_CR) /= Abs_Bot(I_CR)) THEN
!                     XSset_hex(Nxy_CR(I_CR), Iz2, 1,Ig) = Ref_maXS( BotFol_Mat(I_CR), Ig )* W_Vol2 + Ref_maXS( Rep_MAT, Ig )* W_Vol1
!                     XSset_hex(Nxy_CR(I_CR), Iz2, 2,Ig) = Ref_maXS( BotFol_Mat(I_CR), (1*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (1*N_Group + Ig)  )* W_Vol1
!                     XSset_hex(Nxy_CR(I_CR), Iz2, 3,Ig) = Ref_maXS( BotFol_Mat(I_CR), (2*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (2*N_Group + Ig)  )* W_Vol1
!                     XSset_hex(Nxy_CR(I_CR), Iz2, 4,Ig) = Ref_maXS( BotFol_Mat(I_CR), (3*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (3*N_Group + Ig)  )* W_Vol1
!                     XSset_hex(Nxy_CR(I_CR), Iz2, 5,Ig) = Ref_maXS( BotFol_Mat(I_CR), (4*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (4*N_Group + Ig)  )* W_Vol1
!                     XSset_hex(Nxy_CR(I_CR), Iz2, 6,Ig) = Ref_maXS( BotFol_Mat(I_CR), (5*N_Group + Ig) )* W_Vol2 + Ref_maXS( Rep_MAT, (5*N_Group + Ig)  )* W_Vol1
!                     XSset_hex(Nxy_CR(I_CR), Iz2, 6+ig,:) = Ref_smat( BotFol_Mat(I_CR), Ig, :)* W_Vol2 +  Ref_smat( Rep_MAT, Ig, :)* W_Vol1
!                 ENDIF
!                 XSset_hex(Nxy_CR(I_CR), Iz4, 1,Ig) = Ref_maXS( Abs_Mat(I_CR), Ig )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), Ig )* W_Vol3
!                 XSset_hex(Nxy_CR(I_CR), Iz4, 2,Ig) = Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig) )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), (1*N_Group + Ig)  )* W_Vol3
!                 XSset_hex(Nxy_CR(I_CR), Iz4, 3,Ig) = Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig) )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), (2*N_Group + Ig)  )* W_Vol3
!                 XSset_hex(Nxy_CR(I_CR), Iz4, 4,Ig) = Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig) )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), (3*N_Group + Ig)  )* W_Vol3
!                 XSset_hex(Nxy_CR(I_CR), Iz4, 5,Ig) = Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig) )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), (4*N_Group + Ig)  )* W_Vol3
!                 XSset_hex(Nxy_CR(I_CR), Iz4, 6,Ig) = Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig) )* W_Vol4 + Ref_maXS( BotFol_Mat(I_CR), (5*N_Group + Ig)  )* W_Vol3
!                 XSset_hex(Nxy_CR(I_CR), Iz4, 6+ig,:) = Ref_smat( Abs_Mat(I_CR), Ig, :)* W_Vol4 + Ref_smat( BotFol_Mat(I_CR), Ig, :)* W_Vol3
!                                        
!                 XSset_hex(Nxy_CR(I_CR), Iz6, 1,Ig) = Ref_maXS( TopFol_Mat(I_CR), Ig )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), Ig )* W_Vol5
!                 XSset_hex(Nxy_CR(I_CR), Iz6, 2,Ig) = Ref_maXS( TopFol_Mat(I_CR), (1*N_Group + Ig) )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), (1*N_Group + Ig)  )* W_Vol5
!                 XSset_hex(Nxy_CR(I_CR), Iz6, 3,Ig) = Ref_maXS( TopFol_Mat(I_CR), (2*N_Group + Ig) )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), (2*N_Group + Ig)  )* W_Vol5
!                 XSset_hex(Nxy_CR(I_CR), Iz6, 4,Ig) = Ref_maXS( TopFol_Mat(I_CR), (3*N_Group + Ig) )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), (3*N_Group + Ig)  )* W_Vol5
!                 XSset_hex(Nxy_CR(I_CR), Iz6, 5,Ig) = Ref_maXS( TopFol_Mat(I_CR), (4*N_Group + Ig) )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), (4*N_Group + Ig)  )* W_Vol5
!                 XSset_hex(Nxy_CR(I_CR), Iz6, 6,Ig) = Ref_maXS( TopFol_Mat(I_CR), (5*N_Group + Ig) )* W_Vol6 + Ref_maXS( Abs_Mat(I_CR), (5*N_Group + Ig)  )* W_Vol5
!                 XSset_hex(Nxy_CR(I_CR), Iz6, 6+ig,:) = Ref_smat( TopFol_Mat(I_CR), Ig, :)* W_Vol6 + Ref_smat( Abs_Mat(I_CR), Ig, :)* W_Vol5
!                 
!            ENDDO

      
      ENDDO
      
      call WATCH('Control SA move','END')



!// !----------------------------------------
!// ! update XS
!//       call WATCH('Update XS in CRM','START')
!// 
!//          if(.not. allocated(D_3D_MG         ) ) allocate(D_3D_MG            (nxy,nz,n_group))
!//          if(.not. allocated(maXS_tr_3D_MG   ) ) allocate(maXS_tr_3D_MG      (nxy,nz,n_group))
!//          if(.not. allocated(maXS_a_3D_MG    ) ) allocate(maXS_a_3D_MG       (nxy,nz,n_group))
!//          if(.not. allocated(maXS_f_3D_MG    ) ) allocate(maXS_f_3D_MG       (nxy,nz,n_group))
!//          if(.not. allocated(nu_maXS_f_3D_MG ) ) allocate(nu_maXS_f_3D_MG    (nxy,nz,n_group))
!//          if(.not. allocated(kap_maXS_f_3D_MG) ) allocate(kap_maXS_f_3D_MG   (nxy,nz,n_group))
!//          if(.not. allocated(maXS_s_3D_MG    ) ) allocate(maXS_s_3D_MG       (nxy,nz,n_group))
!//          if(.not. allocated(maXS_r_3D_MG    ) ) allocate(maXS_r_3D_MG       (nxy,nz,n_group))
!//          if(.not. allocated(maXS_chi_3D_MG  ) ) allocate(maXS_chi_3D_MG     (nxy,nz,n_group))
!//          if(.not. allocated(maXS_chid_3D_MG ) ) allocate(maXS_chid_3D_MG    (nxy,nz,n_group))
!// !         if (.not.allocated(xssf))      allocate(xssf(n_group-1,2:n_group,nxy,nz))
!//          if (.not.allocated(maXS_r_3D    )) allocate(maXS_r_3D     ( nxy , nz , n_group ))
!//          if (.not.allocated(maXS_scat_3D )) allocate(maXS_scat_3D  ( n_group, n_group, nxy, nz ))
!//          DO Ixy = 1, Nxy
!//              DO Iz = 1, Nz
!//                   DO Ig = 1, N_Group
!//                       maXS_tr_3D_MG    (Ixy, Iz, Ig) = XSset_Hex( Ixy, Iz,1, Ig )
!//                       D_3D_MG          (Ixy, Iz, Ig) = D1 / ( D3*maXS_tr_3D_MG(Ixy, Iz, Ig) )
!//                       maXS_a_3D_MG     (Ixy, Iz, Ig) = XSset_Hex( Ixy, Iz,2, Ig )
!//                       maXS_f_3D_MG     (Ixy, Iz, Ig) = XSset_Hex( Ixy, Iz,3, Ig )
!//                       nu_maXS_f_3D_MG  (Ixy, Iz, Ig) = XSset_Hex( Ixy, Iz,4, Ig )
!//                       kap_maXS_f_3D_MG (Ixy, Iz, Ig) = XSset_Hex( Ixy, Iz,5, Ig )
!// !                      kap_maXS_f_3D (Ixy, Iz, Ig) = XSset_Hex( Ixy, Iz,5, Ig )
!//                       maXs_chi_3D_MG   (Ixy, Iz, Ig) =XSset_Hex( Ixy, Iz,6, Ig )
!//                       maXs_chid_3D_MG  (Ixy, Iz, Ig) = maXs_chi_3D_MG (Ixy, Iz, Ig)
!//                       maXS_scat_3D  (Ig,:,Ixy,Iz) = XSset_Hex(Ixy,Iz,6+Ig, : )
!//                       maXS_r_3D     (Ixy, Iz, Ig) = maXS_a_3D_MG(Ixy, Iz, Ig)
!//                      do ig1 = 1, n_group
!//                         if (ig1 .ne. ig) maXS_r_3D(Ixy, Iz, Ig) =  maXS_r_3D(Ixy, Iz, Ig) + maXS_scat_3D(Ig,ig1,Ixy, Iz)
!//                      enddo
!//                      maXs_s_3D     (Ixy, Iz, 1)  = maXS_scat_3D(1, 2, Ixy, Iz)
!//                   END DO
!//               END DO
!//          END DO
!// 
!// 
!//       do ixy = 1,Nxy
!//          do md=2,n_group
!//             do m=1,md-1
!//                xssf(m,md,ixy,:)=maXS_scat_3D(m,md,imap(ixy),:)
!//             enddo
!//             do m=md+1,n_group
!//                xssf(m-1,md,ixy,:)=maXS_scat_3D(m,md,imap(ixy),:)
!//             enddo
!//          enddo
!//       enddo
!// !-------------------------------------------------------------------------------
!//       call WATCH('Update XS in CRM','END')


      RETURN
      END SUBROUTINE Get_Reg_VF_HEX
#endif 


#endif



#ifdef tuan_fr

!!!        SUBROUTINE miXS_func(a, b)
!!!         use Inc_File
!!!         USE Inc_XS_File
!!!   !      USE Mod_InitCond
!!!   !      USE Mod_GetNode
!!!         use Inc_miXS
!!!         USE Inc_3D, ONLY: T_Fuel, T_Mod, D_Mod
!!!         USE Inc_RP, ONLY: AxialComp
!!!         use Inc_Option,   only: n_group
!!!   
!!!         IMPLICIT NONE
!!!   
!!!        REAL(8), DIMENSION(:,:), ALLOCATABLE :: W_DC3, W_DC2, W_TF1, W_TF2,W_TF3, W_TF4
!!!        REAL(8), DIMENSION(:,:,:), ALLOCATABLE ::   W_DC1
!!!         REAL(8), DIMENSION(:,:), ALLOCATABLE :: Sum_maXS
!!!   
!!!   
!!!        INTEGER :: Im
!!!        INTEGER :: x, y, x1, y1
!!!        INTEGER :: a,b
!!!        INTEGER :: MAT_ID
!!!        INTEGER :: Ixy, Iz, Ixy_FA
!!!   
!!!   
!!!   
!!!   #ifdef siarhei_plot
!!!            ! @$^ siarhei_fr
!!!            current_sub = &
!!!               '+#+ Entered [miXS_func] in Read_File'
!!!            if(trim(last_saved).ne.trim(current_sub)) then
!!!               write(678,'(A)') trim(current_sub)
!!!               last_saved = current_sub
!!!            end if
!!!            ! =-=-=-=-=-=-=-=
!!!   #endif
!!!   
!!!        IF(.not. allocated(W_DC3    ) ) allocate(W_DC3     (N_Group+6, N_Group ))
!!!        IF(.not. allocated(W_TF1    ) ) allocate(W_TF1     (N_Group+6, N_Group ))
!!!        IF(.not. allocated(W_TF2    ) ) allocate(W_TF2     (N_Group+6, N_Group ))
!!!        IF(.not. allocated(W_DC1    ) ) allocate(W_DC1     (NuNum, N_Group+6, N_Group ))
!!!        IF(.not. allocated(W_TF3    ) ) allocate(W_TF3     (N_Group+6, N_Group ))
!!!        IF(.not. allocated(W_TF4    ) ) allocate(W_TF4     (N_Group+6, N_Group ))
!!!        IF(.not. allocated(W_DC2    ) ) allocate(W_DC2     (N_Group+6, N_Group ))
!!!        IF(.not. allocated(Sum_maXS ) ) allocate(Sum_maXS  (N_Group+6, N_Group ))
!!!        IF (a== 0 .AND. b==0) THEN
!!!           DO MAT_ID = 1, N_XS_Table_tmp
!!!               IF(Fiss(MAT_ID) == 'T') THEN
!!!   ! need to be modified for micro interpolation
!!!                    write(*,*) N_TF(MAT_ID), N_DC(MAT_ID)
!!!                    IF  ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) == 1)) THEN
!!!                        DO Im = 1, NuNum
!!!                            DO Iz = 1, Nz
!!!                                DO Ixy = 1, Nxy
!!!                                    IF (AxialComp( I_LP_1N( Ixy ), Iz ) == MAT_ID) THEN
!!!                                        miXS_Hex(Ixy,Iz,Im,:,:)  = miXS_Hex_ref(Im, MAT_ID, StateID( 1, 1), :, :)
!!!                                    ENDIF
!!!                                END DO 
!!!                            ENDDO 
!!!                            Sum_maXS(:,:) = Sum_maXS(:,:) + miXS_Hex_ref(Im, MAT_ID, StateID( 1, 1), :, :)
!!!                        ENDDO 
!!!                        W_DC2(:,:) = XS_residual( MAT_ID, StateID( 1, 1), :, :) 
!!!                        go to 9412
!!!                    ENDIF
!!!   
!!!   !                 DO x = LBOUND(FuelTemp,dim=1), UBOUND(FuelTemp,dim=1)
!!!   !                     IF  ((UBOUND(FuelTemp,dim=1) > 1).AND.(UBOUND(CoolDens,dim=1) == 1).AND.(FuelTemp_Int <= FuelTemp(x)) ) THEN
!!!                    write(*,*) 'MAT_ID, NuNum:', MAT_ID, NuNum
!!!                    DO x = 1, N_TF(MAT_ID)
!!!                        IF  ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) == 1).AND.(FuelTemp_Int <= FuelTemp(x)) ) THEN
!!!                        	   DO Im = 1, NuNum
!!!                                  W_DC1(Im,:,:)  = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, 1), :, :) +&
!!!                        		       (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
!!!                        			   (miXS_Hex_ref(Im, MAT_ID, StateID( x, 1), :, :) - &
!!!                        			   miXS_Hex_ref(Im, MAT_ID, StateID( x-1, 1), :, :))/ &
!!!                        			   (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!                                       Sum_maXS(:,:) = Sum_maXS(:,:) + W_DC1(Im,:,:)
!!!                              END DO 
!!!                              W_DC2(:,:) = XS_residual( MAT_ID, StateID( x-1, 1), :, :) +&
!!!                        	              (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
!!!                        			   (XS_residual( MAT_ID, StateID( x, 1), :, :) - &
!!!                        			   XS_residual( MAT_ID, StateID( x-1, 1), :, :))/ &
!!!                        			   (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!                              go to 9412
!!!                        END IF 
!!!                        DO y = N_DC(MAT_ID), 1,-1
!!!                            IF ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) > 1).AND.(CoolDens_Int <= CoolDens(y)) ) THEN
!!!                            	   DO Im = 1, NuNum
!!!                                       W_DC1(Im,:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( 1, y+1), :, :) + (CoolDens_Int-CoolDens(y+1))* &
!!!                                              (miXS_Hex_ref(Im, MAT_ID, StateID( 1, y), :, :)- miXS_Hex_ref(Im, MAT_ID, StateID( 1, y+1), :, :))/(CoolDens(y)-CoolDens(y+1))
!!!                                              Sum_maXS(:,:) = Sum_maXS(:,:) + W_DC1(Im,:,:)* &
!!!                                              Dens_Int(MAT_ID,Im)
!!!                                   END DO 
!!!                                   W_DC2(:,:) = XS_residual( MAT_ID, StateID( 1, y+1), :, :) + (CoolDens_Int-CoolDens(y+1))* &
!!!                                         (XS_residual(MAT_ID, StateID( 1, y), :, :)- XS_residual(MAT_ID, StateID( 1, y+1), :, :))/(CoolDens(y)-CoolDens(y+1))
!!!                                   go to 9412
!!!                            
!!!                            ELSE IF ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) > 1) &
!!!                                .AND.(CoolDens_Int <= CoolDens(y)) .AND. (FuelTemp_Int <= FuelTemp(x))) THEN
!!!                                Sum_maXS = 0
!!!                                ! micro xs inter polation
!!!                                W_TF1 = 0
!!!                                W_TF2 = 0
!!!                                Sum_maXS = 0
!!!                                DO Im = 1, NuNum
!!!                                   ! should * Density to be maXS
!!!                                    W_TF1(:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y+1), :, :) +&
!!!                                        (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
!!!                                        (miXS_Hex_ref(Im, MAT_ID, StateID( x, y+1), :, :) - &
!!!                                    miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y+1), :, :))/ &
!!!                                        (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!   
!!!                                    W_TF2(:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y), :, :) +&
!!!                                        (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
!!!                                        (miXS_Hex_ref(Im, MAT_ID, StateID( x, y), :, :) - &
!!!                                    miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y), :, :))/ &
!!!                                        (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!   
!!!                                    W_DC1(Im,:,:) = W_TF1(:,:) + (CoolDens_Int-CoolDens(y+1))* &
!!!                                         (W_TF2(:,:)- W_TF1(:,:))/(CoolDens(y)-CoolDens(y+1))
!!!                                    Sum_maXS(:,:) = Sum_maXS(:,:) + W_DC1(Im,:,:)* &
!!!                                                  Dens_Int(MAT_ID,Im)
!!!   
!!!                                 ENDDO 
!!!                           ! save micro XS for depletion
!!!   
!!!   
!!!                                ! residual macro xs interpolation
!!!                                W_TF3(:,:) = XS_residual( MAT_ID, StateID( x-1, y+1), :, :) +&
!!!                                    (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
!!!                                    (XS_residual( MAT_ID, StateID( x, y+1), :, :) - &
!!!                                XS_residual( MAT_ID, StateID( x-1, y+1), :, :))/ &
!!!                                    (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!   
!!!                                W_TF4(:,:) = XS_residual( MAT_ID, StateID( x-1, y), :, :) +&
!!!                                    (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
!!!                                    (XS_residual( MAT_ID, StateID( x, y), :, :) - &
!!!                                XS_residual( MAT_ID, StateID( x-1, y), :, :))/ &
!!!                                    (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!   
!!!                               W_DC2(:,:) = W_TF3(:,:) + (CoolDens_Int-CoolDens(y+1))* &
!!!                                     (W_TF4(:,:)- W_TF3(:,:))/(CoolDens(y)-CoolDens(y+1))
!!!                                go to 9412
!!!                            END IF
!!!                        END DO 
!!!                    END DO 
!!!   9412             write(*,*) 'XS INTERPOLATION:   DONE', MAT_ID
!!!                    DO Iz = 1, Nz
!!!                        DO Ixy = 1, Nxy
!!!                            IF (AxialComp( I_LP_1N( Ixy ), Iz ) == MAT_ID) THEN
!!!                            
!!!                                XSset_Hex(Ixy,Iz,:,:) = Sum_maXS(:,:) + W_DC2(:,:)
!!!                                
!!!                                IF  ((N_TF(MAT_ID) /= 1).OR.(N_DC(MAT_ID) /= 1)) THEN
!!!                                    DO Im = 1, NuNum
!!!                                        miXS_Hex(Ixy,Iz,Im,:,:) = W_DC1(Im,:,:)
!!!                                    END DO  
!!!                                ENDIF
!!!                                
!!!         ! temporary in for 1 state 221 nuclide
!!!   !                             miXS_Hex(Ixy,Iz,:,:,:) = miXS_Hex_ref(:, AxialComp( I_LP_1N( Ixy ), Iz ), 1, :, :)
!!!   !                             XSset_Hex(Ixy,Iz,:,:)  = XSset_Table_Hex(AxialComp( I_LP_1N( Ixy ), Iz ), 1, :, :)
!!!   
!!!                            ELSE
!!!                                XSset_Hex(Ixy,Iz,:,:) = XSset_Table_Hex(AxialComp( I_LP_1N( Ixy ), Iz ), 1, :, :)
!!!                            END IF
!!!                        END DO 
!!!                    END DO 
!!!               END IF
!!!            END DO 
!!!   
!!!        ELSE
!!!            
!!!   
!!!           MAT_ID = AxialComp( I_LP_1N( a ), b )
!!!           W_TF1 = 0 ! weighting factor of fuel temperature
!!!           W_TF2 = 0 ! weighting factor of fuel temperature
!!!           
!!!           IF  ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) == 1)) THEN
!!!   	         DO Im = 1, NuNum
!!!                    W_DC1(Im,:,:)  = miXS_Hex_ref(Im, MAT_ID, StateID( 1, 1), :, :)
!!!                    Sum_maXS(:,:) = Sum_maXS(:,:) + W_DC1(Im,:,:)
!!!                END DO
!!!                W_DC2(:,:) = XS_residual( MAT_ID, StateID( 1, 1), :, :) 
!!!                go to 9415
!!!           ENDIF
!!!   
!!!   		
!!!           DO x = 1, N_TF(MAT_ID)
!!!   
!!!               IF  ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) == 1).AND.(T_Fuel(a,b) <= FuelTemp(x)) ) THEN
!!!               	   DO Im = 1, NuNum
!!!                         W_DC1(Im,:,:)  = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, 1), :, :) +&
!!!               		       (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
!!!               			   (miXS_Hex_ref(Im, MAT_ID, StateID( x, 1), :, :) - &
!!!               			   miXS_Hex_ref(Im, MAT_ID, StateID( x-1, 1), :, :))/ &
!!!               			   (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!                              Sum_maXS(:,:) = Sum_maXS(:,:) + W_DC1(Im,:,:)
!!!                     END DO
!!!                     W_DC2(:,:) = XS_residual( MAT_ID, StateID( x-1, 1), :, :) +&
!!!               	              (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
!!!               			   (XS_residual( MAT_ID, StateID( x, 1), :, :) - &
!!!               			   XS_residual( MAT_ID, StateID( x-1, 1), :, :))/ &
!!!               			   (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!                     go to 9415
!!!               END IF
!!!   
!!!               DO y = N_DC(MAT_ID), 1,-1
!!!   			
!!!                   IF ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) > 1).AND.(D_Mod(a,b) <= CoolDens(y)) ) THEN
!!!                   	   DO Im = 1, NuNum
!!!                              W_DC1(Im,:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( 1, y+1), :, :) + (D_Mod(a,b)-CoolDens(y+1))* &
!!!                                     (miXS_Hex_ref(Im, MAT_ID, StateID( 1, y), :, :)- miXS_Hex_ref(Im, MAT_ID, StateID( 1, y+1), :, :))/(CoolDens(y)-CoolDens(y+1))
!!!                                     Sum_maXS(:,:) = Sum_maXS(:,:) + W_DC1(Im,:,:)* &
!!!                                     Dens_Int(MAT_ID,Im)
!!!                          END DO
!!!                          W_DC2(:,:) = XS_residual( MAT_ID, StateID( 1, y+1), :, :) + (D_Mod(a,b)-CoolDens(y+1))* &
!!!                                (XS_residual(MAT_ID, StateID( 1, y), :, :)- XS_residual(MAT_ID, StateID( 1, y+1), :, :))/(CoolDens(y)-CoolDens(y+1))
!!!                          go to 9415
!!!                   			
!!!                   ELSE IF ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) > 1) &
!!!                                .AND.(D_Mod(a,b) <= CoolDens(y)) .AND. (T_Fuel(a,b) <= FuelTemp(x))) THEN
!!!                       Sum_maXS = 0
!!!   			        
!!!                       DO Im = 1, NuNum
!!!                       ! should * Density to be maXS
!!!                           W_TF1(:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y+1), :, :) +&
!!!                               (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
!!!                               (miXS_Hex_ref(Im, MAT_ID, StateID( x, y+1), :, :) - &
!!!                           miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y+1), :, :))/ &
!!!                               (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!   			        
!!!                           W_TF2(:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y), :, :) +&
!!!                               (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
!!!                               (miXS_Hex_ref(Im, MAT_ID, StateID( x, y), :, :) - &
!!!                           miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y), :, :))/ &
!!!                               (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!   			        
!!!                           W_DC1(Im,:,:) = W_TF1(:,:) + (D_Mod(a,b)-CoolDens(y+1))* &
!!!                                (W_TF2(:,:)- W_TF1(:,:))/(CoolDens(y)-CoolDens(y+1))
!!!   			        
!!!                           Sum_maXS(:,:) = Sum_maXS(:,:) + W_DC1(Im,:,:)* &
!!!                                  N_FR(a,b,Im)
!!!   			        
!!!                       END DO
!!!   			        
!!!                       ! save micro XS for depletion
!!!                       miXS_Hex(a,b,:,:,:) = W_DC1(:,:,:)
!!!   			        
!!!                       ! residual macro xs interpolation
!!!                       W_TF3(:,:) = XS_residual( MAT_ID, StateID( x-1, y+1), :, :) +&
!!!                           (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
!!!                           (XS_residual( MAT_ID, StateID( x, y+1), :, :) - &
!!!                       XS_residual( MAT_ID, StateID( x-1, y+1), :, :))/ &
!!!                           (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!   			        
!!!                       W_TF4(:,:) = XS_residual( MAT_ID, StateID( x-1, y), :, :) +&
!!!                           (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
!!!                           (XS_residual( MAT_ID, StateID( x, y), :, :) - &
!!!                       XS_residual( MAT_ID, StateID( x-1, y), :, :))/ &
!!!                           (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!!!   			        
!!!                       W_DC2(:,:) = W_TF3(:,:) + (D_Mod(a,b)-CoolDens(y+1))* &
!!!                            (W_TF4(:,:)- W_TF3(:,:))/(CoolDens(y)-CoolDens(y+1))
!!!                       exit
!!!                  END IF
!!!   
!!!               END DO
!!!           END DO
!!!   9415 write(*,*) 'XS INTERPOLATION:   DONE', MAT_ID
!!!   
!!!           XSset_Hex(a,b,:,:) = Sum_maXS(:,:) + W_DC2(:,:)
!!!       END IF
!!!   
!!!        END SUBROUTINE miXS_func


     SUBROUTINE miXS_func(a, b)
      use Inc_File
      USE Inc_XS_File

      use Inc_miXS
      USE Inc_3D, ONLY: T_Fuel, T_Mod, D_Mod
      USE Inc_RP, ONLY: AxialComp
      use Inc_Option,   only: n_group

      IMPLICIT NONE

     REAL(8), DIMENSION(:,:), ALLOCATABLE :: W_DC3, W_DC2, W_TF1, W_TF2,W_TF3, W_TF4
     REAL(8), DIMENSION(:,:,:), ALLOCATABLE ::   W_DC1


     INTEGER :: Im
     INTEGER :: x, y, x1, y1
     INTEGER :: a,b
     INTEGER :: MAT_ID
     INTEGER :: Ixy, Iz, Ixy_FA



#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [miXS_func] in Read_File'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

     IF(.not. allocated(W_DC3    ) ) allocate(W_DC3     (N_Group+6, N_Group ))
     IF(.not. allocated(W_TF1    ) ) allocate(W_TF1     (N_Group+6, N_Group ))
     IF(.not. allocated(W_TF2    ) ) allocate(W_TF2     (N_Group+6, N_Group ))
     IF(.not. allocated(W_DC1    ) ) allocate(W_DC1     (NuNum, N_Group+6, N_Group ))
     IF(.not. allocated(W_TF3    ) ) allocate(W_TF3     (N_Group+6, N_Group ))
     IF(.not. allocated(W_TF4    ) ) allocate(W_TF4     (N_Group+6, N_Group ))
     IF(.not. allocated(W_DC2    ) ) allocate(W_DC2     (N_Group+6, N_Group ))



     IF (a== 0 .AND. b==0) THEN
        DO MAT_ID = 1, N_XS_Table_tmp
            IF(Fiss(MAT_ID) == 'T') THEN
! need to be modified for micro interpolation

                 write(*,*) '----- Composition ID:  ', MAT_ID
                 write(*,*) '----- Number of fuel temperature in XS Table:  ', N_TF(MAT_ID)
                 write(*,*) '----- Number of coolant density in XS Table:   ', N_DC(MAT_ID)

                 IF  ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) == 1)) THEN
                     DO Im = 1, NuNum
                         DO Iz = 1, Nz
                             DO Ixy = 1, Nxy
                                 IF (AxialComp( I_LP_1N( Ixy ), Iz ) == MAT_ID) THEN
                                     miXS_Hex(Ixy,Iz,Im,:,:)  = miXS_Hex_ref(Im, MAT_ID, StateID( 1, 1), :, :)
                                     XS_residual( MAT_ID, :, :)  =  XS_residual_ref( MAT_ID, 1, :, :)
                                 ENDIF
                             END DO 
                         ENDDO 
                     ENDDO 
                     WRITE(*,*) "----- No XS interpolation -----"

                     
                    go to 9412
                 ENDIF

                 DO x = 1, N_TF(MAT_ID)
                     IF  ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) == 1).AND.(FuelTemp_Int <= FuelTemp(x)) ) THEN
                           DO Im = 1, NuNum
                               W_DC1(Im,:,:)  = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, 1), :, :) +&
                                   (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
                                   (miXS_Hex_ref(Im, MAT_ID, StateID( x, 1), :, :) - &
                                   miXS_Hex_ref(Im, MAT_ID, StateID( x-1, 1), :, :))/ &
                                   (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
                           END DO 
                           XS_residual( MAT_ID, :, :)  = XS_residual_ref( MAT_ID, StateID( x-1, 1), :, :) +&
                                      (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
                                   (XS_residual_ref( MAT_ID, StateID( x, 1), :, :) - &
                                   XS_residual_ref( MAT_ID, StateID( x-1, 1), :, :))/ &
                                   (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
                           go to 9412
                     END IF 
                     DO y = N_DC(MAT_ID), 1,-1
                         IF ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) > 1).AND.(CoolDens_Int <= CoolDens(y)) ) THEN
                               DO Im = 1, NuNum
                                    W_DC1(Im,:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( 1, y+1), :, :) + (CoolDens_Int-CoolDens(y+1))* &
                                           (miXS_Hex_ref(Im, MAT_ID, StateID( 1, y), :, :)- miXS_Hex_ref(Im, MAT_ID, StateID( 1, y+1), :, :))/(CoolDens(y)-CoolDens(y+1))
                                END DO 
                                XS_residual( MAT_ID, :, :)  = XS_residual_ref( MAT_ID, StateID( 1, y+1), :, :) + (CoolDens_Int-CoolDens(y+1))* &
                                      (XS_residual_ref(MAT_ID, StateID( 1, y), :, :)- XS_residual_ref(MAT_ID, StateID( 1, y+1), :, :))/(CoolDens(y)-CoolDens(y+1))
                                go to 9412
                         
                         ELSE IF ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) > 1) &
                             .AND.(CoolDens_Int <= CoolDens(y)) .AND. (FuelTemp_Int <= FuelTemp(x))) THEN
                             ! micro xs inter polation
                             W_TF1 = 0
                             W_TF2 = 0
                             DO Im = 1, NuNum
                                ! should * Density to be maXS
                                 W_TF1(:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y+1), :, :) +&
                                     (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
                                     (miXS_Hex_ref(Im, MAT_ID, StateID( x, y+1), :, :) - &
                                 miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y+1), :, :))/ &
                                     (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))

                                 W_TF2(:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y), :, :) +&
                                     (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
                                     (miXS_Hex_ref(Im, MAT_ID, StateID( x, y), :, :) - &
                                 miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y), :, :))/ &
                                     (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))

                                 W_DC1(Im,:,:) = W_TF1(:,:) + (CoolDens_Int-CoolDens(y+1))* &
                                      (W_TF2(:,:)- W_TF1(:,:))/(CoolDens(y)-CoolDens(y+1))

                              ENDDO 
                        ! save micro XS for depletion


                             ! residual macro xs interpolation
                             W_TF3(:,:) = XS_residual_ref( MAT_ID, StateID( x-1, y+1), :, :) +&
                                 (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
                                 (XS_residual_ref( MAT_ID, StateID( x, y+1), :, :) - &
                             XS_residual_ref( MAT_ID, StateID( x-1, y+1), :, :))/ &
                                 (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))

                             W_TF4(:,:) = XS_residual_ref( MAT_ID, StateID( x-1, y), :, :) +&
                                 (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
                                 (XS_residual_ref( MAT_ID, StateID( x, y), :, :) - &
                             XS_residual_ref( MAT_ID, StateID( x-1, y), :, :))/ &
                                 (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))

                            XS_residual( MAT_ID, :, :)  = W_TF3(:,:) + (CoolDens_Int-CoolDens(y+1))* &
                                  (W_TF4(:,:)- W_TF3(:,:))/(CoolDens(y)-CoolDens(y+1))
                             go to 9412
                         END IF
                     END DO 
                 END DO 
9412 continue          !  write(*,*) '----- Finished miXS Interpolation for Composition ID: ', MAT_ID
                 DO Iz = 1, Nz
                     DO Ixy = 1, Nxy
                         IF  ((N_TF(MAT_ID) /= 1).AND.(N_DC(MAT_ID) /= 1)) THEN
                             IF (AxialComp( I_LP_1N( Ixy ), Iz ) == MAT_ID) THEN                     
                                 DO Im = 1, NuNum
                                     miXS_Hex(Ixy,Iz,Im,:,:) = W_DC1(Im,:,:)
                                 END DO  
                             ENDIF
                         ENDIF                                                      
                     END DO 
                 END DO 
            END IF
         END DO 

     ELSE
         

        MAT_ID = AxialComp( I_LP_1N( a ), b )
        W_TF1 = 0 ! weighting factor of fuel temperature
        W_TF2 = 0 ! weighting factor of fuel temperature
        
        IF  ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) == 1)) THEN
             DO Im = 1, NuNum
                 miXS_Hex(a,b,Im,:,:)   = miXS_Hex_ref(Im, MAT_ID, StateID( 1, 1), :, :)
             END DO
             XS_residual( MAT_ID, :, :) = XS_residual_ref( MAT_ID, 1, :, :) 
!             WRITE(*,*) "----- No XS interpolation -----"
             return
        ENDIF

        
        DO x = 1, N_TF(MAT_ID)

            IF  ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) == 1).AND.(T_Fuel(a,b) <= FuelTemp(x)) ) THEN
                   DO Im = 1, NuNum
                      miXS_Hex(a,b,Im,:,:)  = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, 1), :, :) +&
                           (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
                           (miXS_Hex_ref(Im, MAT_ID, StateID( x, 1), :, :) - &
                           miXS_Hex_ref(Im, MAT_ID, StateID( x-1, 1), :, :))/ &
                           (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
                  END DO
                  XS_residual( MAT_ID, :, :) = XS_residual_ref( MAT_ID, StateID( x-1, 1), :, :) +&
                              (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
                           (XS_residual_ref( MAT_ID, StateID( x, 1), :, :) - &
                           XS_residual_ref( MAT_ID, StateID( x-1, 1), :, :))/ &
                           (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
                  go to 9415
            END IF

            DO y = N_DC(MAT_ID), 1,-1
            
                IF ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) > 1).AND.(D_Mod(a,b) <= CoolDens(y)) ) THEN
                       DO Im = 1, NuNum
                           miXS_Hex(a,b,Im,:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( 1, y+1), :, :) + (D_Mod(a,b)-CoolDens(y+1))* &
                                  (miXS_Hex_ref(Im, MAT_ID, StateID( 1, y), :, :)- miXS_Hex_ref(Im, MAT_ID, StateID( 1, y+1), :, :))/(CoolDens(y)-CoolDens(y+1))

                       END DO
                       XS_residual( MAT_ID, :, :) = XS_residual_ref( MAT_ID, StateID( 1, y+1), :, :) + (D_Mod(a,b)-CoolDens(y+1))* &
                             (XS_residual_ref(MAT_ID, StateID( 1, y), :, :)- XS_residual_ref(MAT_ID, StateID( 1, y+1), :, :))/(CoolDens(y)-CoolDens(y+1))
                       go to 9415
                            
                ELSE IF ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) > 1) &
                             .AND.(D_Mod(a,b) <= CoolDens(y)) .AND. (T_Fuel(a,b) <= FuelTemp(x))) THEN
                    
                    DO Im = 1, NuNum
                    ! should * Density to be maXS
                        W_TF1(:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y+1), :, :) +&
                            (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
                            (miXS_Hex_ref(Im, MAT_ID, StateID( x, y+1), :, :) - &
                        miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y+1), :, :))/ &
                            (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
                    
                        W_TF2(:,:) = miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y), :, :) +&
                            (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
                            (miXS_Hex_ref(Im, MAT_ID, StateID( x, y), :, :) - &
                        miXS_Hex_ref(Im, MAT_ID, StateID( x-1, y), :, :))/ &
                            (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
                    
                        miXS_Hex(a,b,Im,:,:) = W_TF1(:,:) + (D_Mod(a,b)-CoolDens(y+1))* &
                             (W_TF2(:,:)- W_TF1(:,:))/(CoolDens(y)-CoolDens(y+1))
                                        
                    END DO
                    
                    ! save micro XS for depletion
!                    miXS_Hex(a,b,:,:,:) = W_DC1(:,:,:)
                    
                    ! residual macro xs interpolation
                    W_TF3(:,:) = XS_residual_ref( MAT_ID, StateID( x-1, y+1), :, :) +&
                        (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
                        (XS_residual_ref( MAT_ID, StateID( x, y+1), :, :) - &
                    XS_residual_ref( MAT_ID, StateID( x-1, y+1), :, :))/ &
                        (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
                    
                    W_TF4(:,:) = XS_residual_ref( MAT_ID, StateID( x-1, y), :, :) +&
                        (SQRT(T_Fuel(a,b))-SQRT(FuelTemp(x-1)))*&
                        (XS_residual_ref( MAT_ID, StateID( x, y), :, :) - &
                    XS_residual_ref( MAT_ID, StateID( x-1, y), :, :))/ &
                        (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
                    
                    XS_residual( MAT_ID, :, :) = W_TF3(:,:) + (D_Mod(a,b)-CoolDens(y+1))* &
                         (W_TF4(:,:)- W_TF3(:,:))/(CoolDens(y)-CoolDens(y+1))
                    exit
               END IF

            END DO
        END DO
9415 continue !write(*,*) '----- Finished miXS Interpolation for Composition ID: ', MAT_ID

    END IF

  END SUBROUTINE miXS_func

     
     
#endif


#ifdef tuan_fr_TherEx 

SUBROUTINE spline(x,y,n,yp1,ypn,y2)

    implicit none

    INTEGER, INTENT(in):: n
    INTEGER:: NMAX
    REAL(8), INTENT(in):: x(1:n),y(1:n)
    REAL, INTENT(in):: yp1,ypn
    REAL(8), INTENT(out):: y2(1:n)
    PARAMETER (NMAX=500)
    !  Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
    !  x1 < x2 < ... < xN , and given values yp1 and ypn for the first derivative of the interpolating function at points 1 and n, respectively, this routine returns an array y2(1:n) of
    !  length n which contains the second derivatives of the interpolating function at the tabulated
    !  points xi. If yp1 and/or ypn are equal to 1 ¡¿ 1030 or larger, the routine is signaled to set
    !  the corresponding boundary condition for a natural spline, with zero second derivative on
    !  that boundary.
    !  Parameter: NMAX is the largest anticipated value of n.
    INTEGER:: i,k
    REAL(8):: p,qn,sig,un,u(NMAX)
	
    if (yp1.gt..99e30) then ! The lower boundary condition is set either to be ¡°natural¡±
        y2(1)=0. 
        u(1)=0.
    else ! or else to have a specified first derivative.
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
	
    do i=2,n-1 
!	This is the decomposition loop of the tridiagonal
!    algorithm. y2 and u are used for temporary
!    storage of the decomposed factors.
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo 
	
    if (ypn.gt..99e30) then ! The upper boundary condition is set either to be ¡°natural¡±
        qn=0. 
        un=0.
    else ! or else to have a specified first derivative.
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    
	do k=n-1,1,-1 ! This is the backsubstitution loop of the tridiagonal algorithm.
	    y2(k)=y2(k)*y2(k+1)+u(k) 
    enddo
    return
END



SUBROUTINE splint(xa,ya,y2a,n,x,y)
    implicit none

    INTEGER, INTENT(in):: n
    REAL(8), INTENT(in):: x,xa(n),y2a(n),ya(n)
    REAL(8), INTENT(out):: y
!       Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
!       xai¡¯s in order), and given the array y2a(1:n), which is the output from spline above,
!       and given a value of x, this routine returns a cubic-spline interpolated value y.
    INTEGER:: k,khi,klo
    REAL(8):: a,b,h
    klo=1 
    khi=n
!	We will find the right place in the table by means of bisection.
!    This is optimal if sequential calls to this routine are at random
!    values of x. If sequential calls are in order, and closely
!    spaced, one would do better to store previous values of
!    klo and khi and test if they remain appropriate on the
!    next call.
    1 if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
            khi=k
        else
            klo=k
        endif
        goto 1
    endif      ! klo and khi now bracket the input value of x.
    h=xa(khi)-xa(klo)
    if (h.eq.0.) pause     ! ¡¯bad xa input in splint¡¯ The xa¡¯s must be distinct.
    a=(xa(khi)-x)/h        ! Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+&
	   ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
    return
END SUBROUTINE

      SUBROUTINE Axial_Expansion_XS

      USE Inc_Expansion
      USE Inc_INP
      USE Inc_maXS
      USE Inc_Geometry
      USE Inc_RP
      USE Inc_Option
      USE Inc_TPEN
      
      
      IMPLICIT NONE
      
      INTEGER :: Ixy_FA, Ixy, i, Ig, ig1
      REAL :: W_Vol1, W_Vol2, XXX
      LOGICAL :: Flag_decusping = .false.
      REAL(8) :: flux_r1(N_Group)
      REAL(8) :: flux_u1(N_Group)
      REAL(8) :: WT1(N_Group)
      REAL(8) :: WT2(N_Group)
      INTEGER :: m, n, k

      
      ALLOCATE( Iz_region  (N_FuelRegion) )
      ALLOCATE( ABC  (N_FuelRegion) )
      ALLOCATE( Fuel_length  (N_FuelRegion) )
      ALLOCATE( Fuel_length_Exp  (N_FuelRegion) )
      m=0
      
      
      DO Ixy_FA = 1, Nxy_FA
          Ixy = I_FA(Ixy_FA)
          DO k=IzFuelBot,IzFuelTop
              m=0
              IF (AxialComp( I_LP_1N(Ixy), k ) /= AxialComp( I_LP_1N(Ixy), k+1 )) THEN
                  m=m+1
                  Iz_region(m) = k
                  DO n = 1,N_FuelType
                      IF (AxialComp( I_LP_1N(Ixy), k )==F_ID(n)) THEN
                           ABC(m) = n
                      ENDIF
                  END DO
              END IF
          END DO
          
      ENDDO
      
      DO Ixy_FA = 1, Nxy_FA
          Ixy = I_FA(Ixy_FA)
          m = 0
          DO k = IzFuelBot,IzFuelTop
              IF (AxialComp( I_LP_1N(Ixy), k ) /= AxialComp( I_LP_1N(Ixy), k+1 )) THEN
                  m=m+1
                  Iz_region(m) = k
                  DO n = 1,N_FuelType
                      IF (AxialComp( I_LP_1N(Ixy), k )==F_ID(n)) THEN
                           ABC(m) = n
                      ENDIF
                  END DO
              END IF
          END DO
              
          DO i = 1, m
              IF (i==1) then
                  Fuel_length(i) = sum(MeshSize_z(IzFuelBot:Iz_region(i)))           
              ELSE
                  Fuel_length(i) = sum(MeshSize_z(Iz_region(i-1):Iz_region(i)))           
              ENDIF 
              Fuel_length_Exp(i) = Fuel_length(i)*((Expect_T_F(ABC(i))-Initial_T_F(ABC(i)))*Exp_Coef(ABC(i))+1)
              
              
              XXX = Fuel_length_Exp(i) - Fuel_length(i)
              W_Vol2 = (MeshSize_z(Iz_region(i)+1)-XXX)/MeshSize_z(Iz_region(i)+1)
              IF (W_Vol2 < 0.0001) THEN
                  W_Vol2 = 0D0
              ELSE IF (W_Vol2 > 0.9999) THEN
                  W_Vol2 = 1D0
              ENDIF
              W_Vol1 = 1 - W_Vol2
              
                    
              IF (Flag_decusping) THEN
              
                  flux_r1(:) =(h_z(Iz_region(i))*fluxf(Ixy,Iz_region(i),:) + h_z( Iz_region(i)+1 ) *fluxf(Ixy,Iz_region(i)+1,:))/SUM( h_z(  Iz_region(i) :   (Iz_region(i)+1) ))
                  flux_u1(:) =(h_z(Iz_region(i)+1)*fluxf(Ixy,Iz_region(i)+1,:) + h_z( Iz_region(i)+2 ) *fluxf(Ixy,Iz_region(i)+2,:))/SUM( h_z( (Iz_region(i)+1) : (Iz_region(i)+2) ))
                  WT2 = flux_r1(:)*(XXX)
                  WT1 = flux_u1(:)*(MeshSize_z(Iz_region(i)+1)-XXX)
                  
                  DO Ig = 1, N_Group
                       maXS_tr_3D_MG    (Ixy, Iz_region(i)+1, Ig) = (maXS_tr_3D_MG    (Ixy, Iz_region(i), Ig)* WT2(Ig) + maXS_tr_3D_MG    (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
                       D_3D_MG          (Ixy, Iz_region(i)+1, Ig) = (D_3D_MG          (Ixy, Iz_region(i), Ig)* WT2(Ig) + D_3D_MG          (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
                       maXS_a_3D_MG     (Ixy, Iz_region(i)+1, Ig) = (maXS_a_3D_MG     (Ixy, Iz_region(i), Ig)* WT2(Ig) + maXS_a_3D_MG     (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
                       maXS_f_3D_MG     (Ixy, Iz_region(i)+1, Ig) = (maXS_f_3D_MG     (Ixy, Iz_region(i), Ig)* WT2(Ig) + maXS_f_3D_MG     (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
                       nu_maXS_f_3D_MG  (Ixy, Iz_region(i)+1, Ig) = (nu_maXS_f_3D_MG  (Ixy, Iz_region(i), Ig)* WT2(Ig) + nu_maXS_f_3D_MG  (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
                       kap_maXS_f_3D_MG (Ixy, Iz_region(i)+1, Ig) = (kap_maXS_f_3D_MG (Ixy, Iz_region(i), Ig)* WT2(Ig) + kap_maXS_f_3D_MG (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
                       maXs_chi_3D_MG   (Ixy, Iz_region(i)+1, Ig) = (maXs_chi_3D_MG   (Ixy, Iz_region(i), Ig)* WT2(Ig) + maXs_chi_3D_MG   (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
                       maXs_chid_3D_MG  (Ixy, Iz_region(i)+1, Ig) = (maXs_chid_3D_MG  (Ixy, Iz_region(i), Ig)* WT2(Ig) + maXs_chid_3D_MG  (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                       maXS_scat_3D (Ig,:,Ixy,Iz_region(i)+1) =  maXS_scat_3D (Ig,:,Ixy,Iz_region(i))* W_Vol2 + maXS_scat_3D (Ig,:,Ixy,Iz_region(i)+1)* W_Vol1
                       maXS_r_3D        (Ixy, Iz_region(i)+1, Ig) =  (maXS_r_3D (Ixy, Iz_region(i), Ig)* WT2(Ig) +  maXS_r_3D (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
                      do ig1 = 1, n_group
                          maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i)+1) =  (maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i))* WT2(Ig) + maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i)+1)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
                  
                         if (ig1 .ne. ig) maXS_r_3D(Ixy, Iz_region(i)+1, Ig) =(maXS_r_3D(Ixy, Iz_region(i), Ig)* WT2(Ig) + maXS_r_3D(Ixy, Iz_region(i)+1, Ig) * WT1(Ig))/(WT1(Ig)+WT2(Ig))
                      enddo
                      maXs_s_3D(Ixy, Iz_region(i)+1, 1)  =  (maXs_s_3D(Ixy, Iz_region(i), 1)* WT2(Ig) +  maXs_s_3D(Ixy, Iz_region(i)+1, 1)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
                   END DO
                      
                  
              ENDIF
              
              
              DO Ig = 1, N_Group
                  maXS_tr_3D_MG    (Ixy, Iz_region(i)+1, Ig) = maXS_tr_3D_MG    (Ixy, Iz_region(i), Ig)* W_Vol1 + maXS_tr_3D_MG    (Ixy, Iz_region(i)+1, Ig)* W_Vol2
                  D_3D_MG          (Ixy, Iz_region(i)+1, Ig) = D_3D_MG          (Ixy, Iz_region(i), Ig)* W_Vol1 + D_3D_MG          (Ixy, Iz_region(i)+1, Ig)* W_Vol2
                  maXS_a_3D_MG     (Ixy, Iz_region(i)+1, Ig) = maXS_a_3D_MG     (Ixy, Iz_region(i), Ig)* W_Vol1 + maXS_a_3D_MG     (Ixy, Iz_region(i)+1, Ig)* W_Vol2
                  maXS_f_3D_MG     (Ixy, Iz_region(i)+1, Ig) = maXS_f_3D_MG     (Ixy, Iz_region(i), Ig)* W_Vol1 + maXS_f_3D_MG     (Ixy, Iz_region(i)+1, Ig)* W_Vol2
                  nu_maXS_f_3D_MG  (Ixy, Iz_region(i)+1, Ig) = nu_maXS_f_3D_MG  (Ixy, Iz_region(i), Ig)* W_Vol1 + nu_maXS_f_3D_MG  (Ixy, Iz_region(i)+1, Ig)* W_Vol2
                  kap_maXS_f_3D_MG (Ixy, Iz_region(i)+1, Ig) = kap_maXS_f_3D_MG (Ixy, Iz_region(i), Ig)* W_Vol1 + kap_maXS_f_3D_MG (Ixy, Iz_region(i)+1, Ig)* W_Vol2
                  maXs_chi_3D_MG   (Ixy, Iz_region(i)+1, Ig) = maXs_chi_3D_MG   (Ixy, Iz_region(i), Ig)* W_Vol1 + maXs_chi_3D_MG   (Ixy, Iz_region(i)+1, Ig)* W_Vol2
                  maXs_chid_3D_MG  (Ixy, Iz_region(i)+1, Ig) = maXs_chid_3D_MG  (Ixy, Iz_region(i), Ig)* W_Vol1 + maXs_chid_3D_MG  (Ixy, Iz_region(i)+1, Ig)* W_Vol2
!                  maXS_scat_3D (Ig,:,Ixy,Iz_region(i)+1) =  maXS_scat_3D (Ig,:,Ixy,Iz_region(i))* W_Vol1 + maXS_scat_3D (Ig,:,Ixy,Iz_region(i)+1)* W_Vol1
                  maXS_r_3D        (Ixy, Iz_region(i)+1, Ig) =  maXS_r_3D (Ixy, Iz_region(i), Ig)* W_Vol1 +  maXS_r_3D (Ixy, Iz_region(i)+1, Ig)* W_Vol2
                  do ig1 = 1, n_group
                      maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i)+1) =  maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i))* W_Vol1 + maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i)+1)* W_Vol2
              
                     if (ig1 .ne. ig) maXS_r_3D(Ixy, Iz_region(i)+1, Ig) =maXS_r_3D(Ixy, Iz_region(i), Ig)* W_Vol1 + maXS_r_3D(Ixy, Iz_region(i)+1, Ig) * W_Vol2
                  enddo
                  maXs_s_3D(Ixy, Iz_region(i)+1, 1)  =  maXs_s_3D(Ixy, Iz_region(i), 1)* W_Vol1 +  maXs_s_3D(Ixy, Iz_region(i)+1, 1)* W_Vol2
              END DO
              
              
              
          ENDDO
      END DO
      
          
      
      
!      DO Ixy_FA = 1, Nxy_FA
!          Ixy = I_FA(Ixy_FA)
!          DO i = 1, N_FuelRegion
!              XXX = Fuel_length_Exp(i) - Fuel_length(i)
!              W_Vol2 = (MeshSize_z(Iz_region(i)+1)-XXX)/MeshSize_z(Iz_region(i)+1)
!              IF (W_Vol2 < 0.0001) THEN
!                  W_Vol2 = 0D0
!              ELSE IF (W_Vol2 > 0.9999) THEN
!                  W_Vol2 = 1D0
!              ENDIF
!              W_Vol1 = 1 - W_Vol2
!              
!                    
!              IF (Flag_decusping) THEN
!              
!                  flux_r1(:) =(h_z(Iz_region(i))*fluxf(Ixy,Iz_region(i),:) + h_z( Iz_region(i)+1 ) *fluxf(Ixy,Iz_region(i)+1,:))/SUM( h_z(  Iz_region(i) :   (Iz_region(i)+1) ))
!                  flux_u1(:) =(h_z(Iz_region(i)+1)*fluxf(Ixy,Iz_region(i)+1,:) + h_z( Iz_region(i)+2 ) *fluxf(Ixy,Iz_region(i)+2,:))/SUM( h_z( (Iz_region(i)+1) : (Iz_region(i)+2) ))
!                  WT2 = flux_r1(:)*(XXX)
!                  WT1 = flux_u1(:)*(MeshSize_z(Iz_region(i)+1)-XXX)
!                  
!                  DO Ig = 1, N_Group
!                       maXS_tr_3D_MG    (Ixy, Iz_region(i)+1, Ig) = (maXS_tr_3D_MG    (Ixy, Iz_region(i), Ig)* WT2(Ig) + maXS_tr_3D_MG    (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                       D_3D_MG          (Ixy, Iz_region(i)+1, Ig) = (D_3D_MG          (Ixy, Iz_region(i), Ig)* WT2(Ig) + D_3D_MG          (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                       maXS_a_3D_MG     (Ixy, Iz_region(i)+1, Ig) = (maXS_a_3D_MG     (Ixy, Iz_region(i), Ig)* WT2(Ig) + maXS_a_3D_MG     (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                       maXS_f_3D_MG     (Ixy, Iz_region(i)+1, Ig) = (maXS_f_3D_MG     (Ixy, Iz_region(i), Ig)* WT2(Ig) + maXS_f_3D_MG     (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                       nu_maXS_f_3D_MG  (Ixy, Iz_region(i)+1, Ig) = (nu_maXS_f_3D_MG  (Ixy, Iz_region(i), Ig)* WT2(Ig) + nu_maXS_f_3D_MG  (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                       kap_maXS_f_3D_MG (Ixy, Iz_region(i)+1, Ig) = (kap_maXS_f_3D_MG (Ixy, Iz_region(i), Ig)* WT2(Ig) + kap_maXS_f_3D_MG (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                       maXs_chi_3D_MG   (Ixy, Iz_region(i)+1, Ig) = (maXs_chi_3D_MG   (Ixy, Iz_region(i), Ig)* WT2(Ig) + maXs_chi_3D_MG   (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                       maXs_chid_3D_MG  (Ixy, Iz_region(i)+1, Ig) = (maXs_chid_3D_MG  (Ixy, Iz_region(i), Ig)* WT2(Ig) + maXs_chid_3D_MG  (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!!                       maXS_scat_3D (Ig,:,Ixy,Iz_region(i)+1) =  maXS_scat_3D (Ig,:,Ixy,Iz_region(i))* W_Vol2 + maXS_scat_3D (Ig,:,Ixy,Iz_region(i)+1)* W_Vol1
!                       maXS_r_3D        (Ixy, Iz_region(i)+1, Ig) =  (maXS_r_3D (Ixy, Iz_region(i), Ig)* WT2(Ig) +  maXS_r_3D (Ixy, Iz_region(i)+1, Ig)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                      do ig1 = 1, n_group
!                          maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i)+1) =  (maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i))* WT2(Ig) + maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i)+1)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                  
!                         if (ig1 .ne. ig) maXS_r_3D(Ixy, Iz_region(i)+1, Ig) =(maXS_r_3D(Ixy, Iz_region(i), Ig)* WT2(Ig) + maXS_r_3D(Ixy, Iz_region(i)+1, Ig) * WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                      enddo
!                      maXs_s_3D(Ixy, Iz_region(i)+1, 1)  =  (maXs_s_3D(Ixy, Iz_region(i), 1)* WT2(Ig) +  maXs_s_3D(Ixy, Iz_region(i)+1, 1)* WT1(Ig))/(WT1(Ig)+WT2(Ig))
!                   END DO
!                      
!                  
!              ENDIF
!              
!              
!              DO Ig = 1, N_Group
!                  maXS_tr_3D_MG    (Ixy, Iz_region(i)+1, Ig) = maXS_tr_3D_MG    (Ixy, Iz_region(i), Ig)* W_Vol2 + maXS_tr_3D_MG    (Ixy, Iz_region(i)+1, Ig)* W_Vol1
!                  D_3D_MG          (Ixy, Iz_region(i)+1, Ig) = D_3D_MG          (Ixy, Iz_region(i), Ig)* W_Vol2 + D_3D_MG          (Ixy, Iz_region(i)+1, Ig)* W_Vol1
!                  maXS_a_3D_MG     (Ixy, Iz_region(i)+1, Ig) = maXS_a_3D_MG     (Ixy, Iz_region(i), Ig)* W_Vol2 + maXS_a_3D_MG     (Ixy, Iz_region(i)+1, Ig)* W_Vol1
!                  maXS_f_3D_MG     (Ixy, Iz_region(i)+1, Ig) = maXS_f_3D_MG     (Ixy, Iz_region(i), Ig)* W_Vol2 + maXS_f_3D_MG     (Ixy, Iz_region(i)+1, Ig)* W_Vol1
!                  nu_maXS_f_3D_MG  (Ixy, Iz_region(i)+1, Ig) = nu_maXS_f_3D_MG  (Ixy, Iz_region(i), Ig)* W_Vol2 + nu_maXS_f_3D_MG  (Ixy, Iz_region(i)+1, Ig)* W_Vol1
!                  kap_maXS_f_3D_MG (Ixy, Iz_region(i)+1, Ig) = kap_maXS_f_3D_MG (Ixy, Iz_region(i), Ig)* W_Vol2 + kap_maXS_f_3D_MG (Ixy, Iz_region(i)+1, Ig)* W_Vol1
!                  maXs_chi_3D_MG   (Ixy, Iz_region(i)+1, Ig) = maXs_chi_3D_MG   (Ixy, Iz_region(i), Ig)* W_Vol2 + maXs_chi_3D_MG   (Ixy, Iz_region(i)+1, Ig)* W_Vol1
!                  maXs_chid_3D_MG  (Ixy, Iz_region(i)+1, Ig) = maXs_chid_3D_MG  (Ixy, Iz_region(i), Ig)* W_Vol2 + maXs_chid_3D_MG  (Ixy, Iz_region(i)+1, Ig)* W_Vol1
!!                  maXS_scat_3D (Ig,:,Ixy,Iz_region(i)+1) =  maXS_scat_3D (Ig,:,Ixy,Iz_region(i))* W_Vol2 + maXS_scat_3D (Ig,:,Ixy,Iz_region(i)+1)* W_Vol1
!                  maXS_r_3D        (Ixy, Iz_region(i)+1, Ig) =  maXS_r_3D (Ixy, Iz_region(i), Ig)* W_Vol2 +  maXS_r_3D (Ixy, Iz_region(i)+1, Ig)* W_Vol1
!                  do ig1 = 1, n_group
!                      maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i)+1) =  maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i))* W_Vol2 + maXS_scat_3D (Ig,ig1,Ixy,Iz_region(i)+1)* W_Vol1
!              
!                     if (ig1 .ne. ig) maXS_r_3D(Ixy, Iz_region(i)+1, Ig) =maXS_r_3D(Ixy, Iz_region(i), Ig)* W_Vol2 + maXS_r_3D(Ixy, Iz_region(i)+1, Ig) * W_Vol1
!                  enddo
!                  maXs_s_3D(Ixy, Iz_region(i)+1, 1)  =  maXs_s_3D(Ixy, Iz_region(i), 1)* W_Vol2 +  maXs_s_3D(Ixy, Iz_region(i)+1, 1)* W_Vol1
!              END DO
!              
!              
!          END DO
!      END DO
      
      END SUBROUTINE Axial_Expansion_XS


      SUBROUTINE Radial_expansion

      USE Inc_Expansion
      USE Inc_INP
      USE Inc_maXS
      USE Inc_Geometry
      USE Inc_RP
      USE Inc_Option
      USE Inc_TPEN
      use Inc_miXS, ONLY: N_FR, ZAID
      use Inc_File,     ONLY: NuNum
!       USE Read_File, ONLY:  update_maXS
      
      IMPLICIT NONE
      
      INTEGER :: Ixy_FA, Ixy
      REAL(8) :: R1, R2,F, W_Fuel, W_Cool
      INTEGER :: m, n, k, i,j,Iz

      Radial_exp = .true.
      
      R1 = GridSize_x(1)
      if (.not.allocated(N_FR_Exp )) Allocate(N_FR_Exp( Nxy, Nz, NuNum ))
      DO Ixy_FA = 1, Nxy_FA
          Ixy = I_FA(Ixy_FA)
          
          R2 = GridSize_x(1) * ((Expect_T_F(1)-Initial_T_F(1))*Exp_Coef(1)+1)
          
          W_Fuel = R1/R2
          W_Cool = (R1/R2)*((R2/R1)**2-F)/(1-F)
          DO i = 1, NuNum
              DO j = 1,  N_CooMAT
                  IF(ZAID(i) == CooMAT_ID(j)) THEN
                      N_FR(Ixy,1,i) = N_FR(Ixy,1,i)  * W_Cool       
                      GO TO 2101
                  END IF 
              ENDDO
              N_FR(Ixy,1,i) = N_FR(Ixy,1,i)* W_Fuel
2101                    write(*,*) 'Expeansion ----'
          ENDDO
! update_maXS wih new density

!         DO Ixy_FA = 1, Nxy_FA
!            Ixy = I_FA(Ixy_FA)

            DO Iz = IzFuelBot, IzFuelTop
                call  update_maXS(Ixy,Iz)
            END DO
!         END DO
      END DO
      
      
      
      END SUBROUTINE Radial_expansion

      
#endif


!#ifdef tuan_fr_TherEx 


SUBROUTINE maXS_func(a, b)
      use Inc_File
      USE Inc_XS_File
!      USE Mod_InitCond
!      USE Mod_GetNode
      use Inc_miXS
      USE Inc_3D, ONLY: T_Fuel, T_Mod, D_Mod
      USE Inc_RP, ONLY: AxialComp
      use Inc_Option,   only: n_group
      use splines, only: spline3
      use Inc_TH, ONLY: tfuel,dcool,tdopl
      
      
    implicit none
	
	INTEGER,INTENT(in)  :: a,b
    INTEGER             :: x,y
    INTEGER             :: Ixy, Iz, a_xy
    INTEGER             :: MAT_ID
    INTEGER             :: Ig1,Ig2

    REAL, DIMENSION(:,:), ALLOCATABLE :: W_TF1, W_TF2
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE ::   W_DC1
    
    LOGICAL :: Flag_linear  = .false.
    LOGICAL :: Flag_cubic = .false.
    

     REAL(8):: dummy1(1:N_Group) 
     REAL(8):: dummy2(1:N_Group) 
     REAL(8):: dummy3(1:2) 
    
    IF(.not. allocated(W_TF1)) allocate( W_TF1(N_Group+6, N_Group ))
    IF(.not. allocated(W_TF2)) allocate( W_TF2(N_Group+6, N_Group ))
    IF(.not. allocated(W_DC1)) allocate( W_DC1(N_XS_Table_tmp, N_Group+6, N_Group ))

!	W_DC1(:,:,;) = 0
!     write(*,*) FuelTemp(1:N_TF(1)),XSset_Table_Hex(  1, StateID( 1:N_TF(1), 1),1,1)
!     pause
    IF (a== 0 .AND. b==0) THEN
        DO MAT_ID = 1, N_XS_Table_tmp
            ! consider only fuel assembly
            IF(Fiss(MAT_ID) == 'T') THEN
                IF  ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) == 1)) THEN
                     W_DC1(MAT_ID,:,:) = XSset_Table_Hex( MAT_ID, StateID( 1, 1), :, :) 
                     go to 9418
                ENDIF

        	    DO x = 1, N_TF(MAT_ID)
                    
                    IF ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) == 1).AND.(FuelTemp_Int <= FuelTemp(x)) ) THEN

!                        IF (Flag_cubic) Then
!                             DO Ig1 = 1, N_Group
!     						    DO Ig2 = 1, N_Group + 6
!                                      call spline(FuelTemp(1:N_TF(MAT_ID)),XSset_Table_Hex(  MAT_ID, StateID( 1:N_TF(MAT_ID), 1), Ig2, Ig1),N_TF(MAT_ID),1e30,1e30,dummy1(1:N_TF(MAT_ID)))
!                                      call splint(FuelTemp(1:N_TF(MAT_ID)),XSset_Table_Hex(  MAT_ID, StateID( 1:N_TF(MAT_ID), 1), Ig2, Ig1),dummy1(1:N_TF(MAT_ID)),N_TF(MAT_ID),FuelTemp_Int,dummy2(1) )
!                                      W_DC1(MAT_ID,Ig2,Ig1)  = dummy2(1)
!     					        END DO
!                             END DO
!                        ELSE
                            W_DC1(MAT_ID,:,:)  = XSset_Table_Hex(  MAT_ID, StateID( x-1, 1), :, :) +&
                      		  (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
                      		  (XSset_Table_Hex(  MAT_ID, StateID( x, 1), :, :) - &
                      		  XSset_Table_Hex(  MAT_ID, StateID( x-1, 1), :, :))/ &
                      		  (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
!                        ENDIF
                                                
                        go to 9418
                    END IF 
                    
                    DO y =  N_DC(MAT_ID), 1,-1
                        IF ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) > 1).AND.(CoolDens_Int <= CoolDens(y)) ) THEN
                            
!                            IF (Flag_cubic) Then
!                                 DO Ig1 = 1, N_Group
!     					     	    DO Ig2 = 1, N_Group + 6
!                                          call spline(CoolDens(1:N_DC(MAT_ID)),XSset_Table_Hex(  MAT_ID, StateID( 1, 1:N_DC(MAT_ID)), Ig2, Ig1),N_DC(MAT_ID),1e30,1e30,dummy1(1:N_DC(MAT_ID)))
!                                          call splint(CoolDens(1:N_DC(MAT_ID)),XSset_Table_Hex(  MAT_ID, StateID( 1, 1:N_DC(MAT_ID)), Ig2, Ig1),dummy1(1:N_DC(MAT_ID)),N_DC(MAT_ID),CoolDens_Int,dummy2(1) )
!                                          W_DC1(MAT_ID,Ig2,Ig1)  = dummy2(1)
!     					            END DO
!                                 END DO
!                            ELSE
                                W_DC1(MAT_ID,:,:) = XSset_Table_Hex(  MAT_ID, StateID( 1, y+1), :, :) + (CoolDens_Int-CoolDens(y+1))* &
                                       (XSset_Table_Hex(  MAT_ID, StateID( 1, y), :, :)- XSset_Table_Hex(  MAT_ID, StateID( 1, y+1), :, :))/(CoolDens(y)-CoolDens(y+1))
!                            ENDIF
                             go to 9418
                        
        		        ELSE IF ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) > 1) &
                             .AND.(CoolDens_Int <= CoolDens(y)) .AND. (FuelTemp_Int <= FuelTemp(x))) THEN
					        W_TF1 = 0
                            W_TF2 = 0
!                            IF (Flag_cubic) Then
!                                 DO Ig1 = 1, N_Group
!     					     	    DO Ig2 = 1, N_Group + 6
!                                      call spline(FuelTemp(1:N_TF(MAT_ID)),XSset_Table_Hex(  MAT_ID, StateID( 1:N_TF(MAT_ID), y+1), Ig2, Ig1),N_TF(MAT_ID),1e30,1e30,dummy1(1:N_TF(MAT_ID)))
!                                      call splint(FuelTemp(1:N_TF(MAT_ID)),XSset_Table_Hex(  MAT_ID, StateID( 1:N_TF(MAT_ID), y+1), Ig2, Ig1),dummy1(1:N_TF(MAT_ID)),N_TF(MAT_ID),FuelTemp_Int,dummy2(1) )
!                                      W_TF1(Ig2,Ig1)  = dummy2(1)
!                                      call spline(FuelTemp(1:N_TF(MAT_ID)),XSset_Table_Hex(  MAT_ID, StateID( 1:N_TF(MAT_ID), y ), Ig2, Ig1),N_TF(MAT_ID),1e30,1e30,dummy1(1:N_TF(MAT_ID)))
!                                      call splint(FuelTemp(1:N_TF(MAT_ID)),XSset_Table_Hex(  MAT_ID, StateID( 1:N_TF(MAT_ID), y ), Ig2, Ig1),dummy1(1:N_TF(MAT_ID)),N_TF(MAT_ID),FuelTemp_Int,dummy2(2) )
!                                      W_TF2(Ig2,Ig1)  = dummy2(2)
!     					            END DO
!                                 END DO
!                                 DO Ig1 = 1, N_Group
!     					     	    DO Ig2 = 1, N_Group + 6
!                                        dummy3(1) = W_TF2(Ig2,Ig1)
!                                        dummy3(2) = W_TF1(Ig2,Ig1)
!                                        call spline(CoolDens(y:y+1),dummy3(1:2),2,1e30,1e30,dummy1(1:2))
!                                        call splint(CoolDens(y:y+1),dummy3(1:2),dummy1(1:2),2,CoolDens_Int,dummy2(3) )
!                                        W_DC1(MAT_ID,Ig2,Ig1)  = dummy2(3)
!     					            END DO
!                                 END DO
                                 
!                            ELSE
                            
                                W_TF1(:,:) = XSset_Table_Hex( MAT_ID, StateID( x-1, y+1), :, :) +&
                                    (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
                                    (XSset_Table_Hex( MAT_ID, StateID( x, y+1), :, :) - &
                                XSset_Table_Hex( MAT_ID, StateID( x-1, y+1), :, :))/ &
                                    (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
                                
                                W_TF2(:,:) = XSset_Table_Hex( MAT_ID, StateID( x-1, y), :, :) +&
                                    (SQRT(FuelTemp_Int)-SQRT(FuelTemp(x-1)))*&
                                    (XSset_Table_Hex( MAT_ID, StateID( x, y), :, :) - &
                                XSset_Table_Hex( MAT_ID, StateID( x-1, y), :, :))/ &
                                    (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
                                
                                W_DC1(MAT_ID,:,:) = W_TF1(:,:) + (CoolDens_Int-CoolDens(y+1))* &
                                     (W_TF2(:,:)- W_TF1(:,:))/(CoolDens(y)-CoolDens(y+1))
!                            ENDIF
                            

                            go to 9418
        		 	   END IF
        		    
        		    ENDDO
        	    ENDDO
            9418 continue !write(*,*) 'maXS INTERPOLATION:   DONE', MAT_ID
            END IF
        ENDDO
        	
        DO Ixy = 1, Nxy
            DO Iz = 1, Nz
        	    IF (Fiss(AxialComp( I_LP_1N( Ixy ), Iz )) == 'T') THEN
                    XSset_Hex(Ixy,Iz,:,:) = W_DC1(AxialComp( I_LP_1N( Ixy ), Iz ),:,:) 
        		ELSE
                    XSset_Hex(Ixy,Iz,:,:) = XSset_Table_Hex(AxialComp( I_LP_1N( Ixy ), Iz ), 1, :, :)
        			
        		ENDIF
        	ENDDO
        ENDDO
	   
    ELSE
        a_xy =I_FA (a)
        MAT_ID = AxialComp( I_LP_1N( a_xy ), b )
        W_TF1 = 0 ! weighting factor of fuel temperature
        W_TF2 = 0 ! weighting factor of fuel temperature
        IF(Fiss(MAT_ID) == 'T') THEN
        ! need to be modified for micro interpolation
           IF  ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) == 1)) THEN
                W_DC1(MAT_ID,:,:) = XSset_Table_Hex( MAT_ID, StateID( 1, 1), :, :) 
                go to 9413
           ENDIF
            
           DO x = 1, N_TF(MAT_ID)

               IF  ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) == 1).AND.((tdopl(b,a)**2) <= FuelTemp(x)) ) THEN
                     W_DC1(MAT_ID,:,:)  = XSset_Table_Hex(  MAT_ID, StateID( x-1, 1), :, :) +&
               		  (SQRT((tdopl(b,a)**2))-SQRT(FuelTemp(x-1)))*&
               		  (XSset_Table_Hex(  MAT_ID, StateID( x, 1), :, :) - &
               		  XSset_Table_Hex(  MAT_ID, StateID( x-1, 1), :, :))/ &
               		  (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
                     go to 9413
               END IF 
                    
               
               
               DO y =  N_DC(MAT_ID), 1,-1
                   IF ((N_TF(MAT_ID) == 1).AND.(N_DC(MAT_ID) > 1).AND.(dcool(b,a)*1d-3 <= CoolDens(y)) ) THEN
                           W_DC1(MAT_ID,:,:) = XSset_Table_Hex( MAT_ID, StateID( 1, y+1), :, :) + (dcool(b,a)*1d-3 -CoolDens(y+1))* &
                                  (XSset_Table_Hex( MAT_ID, StateID( 1, y), :, :)- XSset_Table_Hex( MAT_ID, StateID( 1, y+1), :, :))/(CoolDens(y)-CoolDens(y+1))

                          go to 9413
                   
                   ELSE IF ((N_TF(MAT_ID) > 1).AND.(N_DC(MAT_ID) > 1) &
                       .AND.(dcool(b,a)*1d-3  <= CoolDens(y)) .AND. ((tdopl(b,a)**2) <= FuelTemp(x))) THEN
                       W_TF1(:,:) = XSset_Table_Hex( MAT_ID, StateID( x-1, y+1), :, :) +&
                           (SQRT((tdopl(b,a)**2))-SQRT(FuelTemp(x-1)))*&
                           (XSset_Table_Hex( MAT_ID, StateID( x, y+1), :, :) - &
                       XSset_Table_Hex( MAT_ID, StateID( x-1, y+1), :, :))/ &
                           (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
         
                       W_TF2(:,:) = XSset_Table_Hex( MAT_ID, StateID( x-1, y), :, :) +&
                           (SQRT((tdopl(b,a)**2))-SQRT(FuelTemp(x-1)))*&
                           (XSset_Table_Hex( MAT_ID, StateID( x, y), :, :) - &
                       XSset_Table_Hex( MAT_ID, StateID( x-1, y), :, :))/ &
                           (SQRT(FuelTemp(x))-SQRT(FuelTemp(x-1)))
         
                       W_DC1(MAT_ID,:,:) = W_TF1(:,:) + (dcool(b,a)*1d-3 -CoolDens(y+1))* &
                            (W_TF2(:,:)- W_TF1(:,:))/(CoolDens(y)-CoolDens(y+1))
                       go to 9413
            	   END IF
               
               ENDDO
           ENDDO
		   
        9413 continue !write(*,*) 'maXS INTERPOLATION:   DONE', MAT_ID
		XSset_Hex(a_xy,b,:,:) = W_DC1(MAT_ID,:,:)
        END IF

	ENDIF

	
END SUBROUTINE 

#ifdef tuan_fr_TherEx 

SUBROUTINE Dens_update(x,y,a)
    ! number density of nuclides N_FR
    USE Inc_miXS, ONLY: N_FR, N_FR_old
!	USE Thermal_INP
	
    implicit none
	
	INTEGER,INTENT(in)  :: x,y,a
    INTEGER             :: i,j
    REAL(8)             :: exp_value
           
    IF(a==1) THEN	
	! exp_value = delta T*thermal expansion coefficient
	    N_FR(x,y,:) =  N_FR(x,y,:)/exp_value
	ELSE IF (a==2) THEN
	    N_FR(x,y,:) =  N_FR(x,y,:)/exp_value/exp_value
		
    ELSE IF (a==3) THEN
	    N_FR(x,y,:) =  N_FR(x,y,:)/exp_value/exp_value/exp_value
	   
    ENDIF
	
	
END SUBROUTINE 


#endif 

      subroutine Get_ADF
      use Inc_DF
      use Inc_INP, only: Flag_Card_maXS
      use Inc_Option, only: N_Group
      use Inc_XS_File, only: Flag_RefDF, flag_1nadf
      implicit none
      integer :: Ixy, Ixy_FA, Ixy_RF, Ixy_1N, Iz, Ig, I_Tab, m

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_ADF] in Mod_GetSome'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_ADF] in Mod_GetSome'
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

      if (Flag_Card_maXS) return

      Buff_ADF_Lx=ADF_Lx
      Buff_ADF_Rx=ADF_Rx
      Buff_ADF_Ly=ADF_Ly
      Buff_ADF_Ry=ADF_Ry
      ADF_Lz=1d0
      ADF_Rz=1d0

      do Iz=1,(IzFuelBot-1)
         do Ixy=1,Nxy
            Ixy_1N=I_4Nto1N(Ixy)
            do Ig=1,N_Group
               if (I_Rz(Iz)/= 0) then
                  if (I_FARF_1N_3D(Ixy_1N,I_Rz(Iz))>=2) then
                     ADF_Rz(Ixy,Iz,Ig)=Buff_ADF_Lx(Ixy,Iz,Ig)
                  endif
               endif
               ADF_Lx(Ixy,Iz,Ig)=1d0
               ADF_Rx(Ixy,Iz,Ig)=1d0
               ADF_Ly(Ixy,Iz,Ig)=1d0
               ADF_Ry(Ixy,Iz,Ig)=1d0
            enddo
         enddo
      enddo

      do Iz=(IzFuelTop+1),Nz
         do Ixy=1,Nxy
            Ixy_1N=I_4Nto1N(Ixy)
            do Ig=1,N_Group
               if (I_Lz(Iz)/=0) then
                  if (I_FARF_1N_3D(Ixy_1N,I_Lz(Iz))>=2) then
                     ADF_Lz(Ixy,Iz,Ig)=Buff_ADF_Lx(Ixy,Iz,Ig)
                  endif
               endif
               ADF_Lx(Ixy,Iz,Ig)=1d0
               ADF_Rx(Ixy,Iz,Ig)=1d0
               ADF_Ly(Ixy,Iz,Ig)=1d0
               ADF_Ry(Ixy,Iz,Ig)=1d0
            enddo
         enddo
      enddo

      do Iz=IzFuelBot,IzFuelTop
         do Ixy_RF=1,Nxy_RF
            Ixy=I_RF(Ixy_RF)
            do Ig=1,N_Group
               if (I_Lx(Ixy)==0) then
                  ADF_Lx(Ixy,Iz,Ig)=1d0
               else
#ifdef tuan_fr_crm                  
                  if (I_FARF_1N(I_4Nto1N(I_Lx(Ixy)))==2) then
#else                      
                  if (I_FARF_1N(I_4Nto1N(I_Lx(Ixy)))>=2) then
#endif
                     ADF_Lx(Ixy,Iz,Ig)=Buff_ADF_Lx(Ixy,Iz,Ig)
                  else
                     ADF_Lx(Ixy,Iz,Ig)=1d0
                  endif
               endif
               if (I_Rx(Ixy)==0) then
                  ADF_Rx(Ixy,Iz,Ig)=1d0
               else
#ifdef tuan_fr_crm                  
                  if (I_FARF_1N(I_4Nto1N(I_Rx(Ixy)))==2) then
#else
                  if (I_FARF_1N(I_4Nto1N(I_Rx(Ixy)))>=2) then
#endif
                     ADF_Rx(Ixy,Iz,Ig)=Buff_ADF_Lx(Ixy,Iz,Ig)
                  else
                     ADF_Rx(Ixy,Iz,Ig)=1d0
                  endif
               endif
               if (I_Ly(Ixy)==0) then
                  ADF_Ly(Ixy,Iz,Ig)=1d0
               else
#ifdef tuan_fr_crm                  
                  if (I_FARF_1N(I_4Nto1N(I_Ly(Ixy)))==2) then
#else
                  if (I_FARF_1N(I_4Nto1N(I_Ly(Ixy)))>=2) then
#endif
                     ADF_Ly(Ixy,Iz,Ig)=Buff_ADF_Lx(Ixy,Iz,Ig)
                  else
                     ADF_Ly(Ixy,Iz,Ig)=1d0
                  endif
               endif
               if (I_Ry(Ixy)==0) then
                  ADF_Ry(Ixy,Iz,Ig)=1d0
               else
#ifdef tuan_fr_crm                  
                  if (I_FARF_1N(I_4Nto1N(I_Ry(Ixy)))==2) then
#else
                  if (I_FARF_1N(I_4Nto1N(I_Ry(Ixy)))>=2) then
#endif
                     ADF_Ry(Ixy,Iz,Ig)=Buff_ADF_Lx(Ixy,Iz,Ig)
                  else
                     ADF_Ry(Ixy,Iz,Ig)=1d0
                  endif
               endif
            enddo
         enddo
      enddo

      do Iz=IzFuelBot,IzFuelTop
         do Ixy_FA=1,Nxy_FA
            Ixy=I_FA(Ixy_FA)
            Ixy_1N=I_4Nto1N(Ixy)
            I_Tab=AxialComp(I_LP_1N(Ixy_1N),Iz)
            do Ig=1,N_Group
               do m=1,4
                  if (Ixy==I_1Nto4N(Ixy_1N,m)) then
                     select case (m)
                     case (1)
                        ADF_Rx(Ixy,Iz,Ig)=Buff_ADF_Lx(Ixy,Iz,Ig)
                        ADF_Ry(Ixy,Iz,Ig)=Buff_ADF_Ly(Ixy,Iz,Ig)
                        ADF_Lx(Ixy,Iz,Ig)=Buff_ADF_Rx(Ixy,Iz,Ig)
                        ADF_Ly(Ixy,Iz,Ig)=Buff_ADF_Ry(Ixy,Iz,Ig)
                        if (.not.flag_4n1fa) then
                           if (flag_1nadf) then
                              ADF_Rx(Ixy,Iz,Ig)=Buff_ADF_Rx(Ixy,Iz,Ig)
                              ADF_Ry(Ixy,Iz,Ig)=Buff_ADF_Ry(Ixy,Iz,Ig)
                              ADF_Lx(Ixy,Iz,Ig)=Buff_ADF_Lx(Ixy,Iz,Ig)
                              ADF_Ly(Ixy,Iz,Ig)=Buff_ADF_Ly(Ixy,Iz,Ig)
                           else
                              ADF_Rx(Ixy,Iz,Ig)=Buff_ADF_Rx(Ixy,Iz,Ig)
                              ADF_Ry(Ixy,Iz,Ig)=Buff_ADF_Rx(Ixy,Iz,Ig)
                              ADF_Lx(Ixy,Iz,Ig)=Buff_ADF_Rx(Ixy,Iz,Ig)
                              ADF_Ly(Ixy,Iz,Ig)=Buff_ADF_Rx(Ixy,Iz,Ig)
                           endif
                        endif
                     case (2)
                        ADF_Lx(Ixy,Iz,Ig)=Buff_ADF_Ly(Ixy,Iz,Ig)
                        ADF_Ry(Ixy,Iz,Ig)=Buff_ADF_Lx(Ixy,Iz,Ig)
                        ADF_Ly(Ixy,Iz,Ig)=Buff_ADF_Rx(Ixy,Iz,Ig)
                        ADF_Rx(Ixy,Iz,Ig)=Buff_ADF_Ry(Ixy,Iz,Ig)
                     case (3)
                        ADF_Rx(Ixy,Iz,Ig)=Buff_ADF_Ly(Ixy,Iz,Ig)
                        ADF_Ly(Ixy,Iz,Ig)=Buff_ADF_Lx(Ixy,Iz,Ig)
                        ADF_Lx(Ixy,Iz,Ig)=Buff_ADF_Ry(Ixy,Iz,Ig)
                        ADF_Ry(Ixy,Iz,Ig)=Buff_ADF_Rx(Ixy,Iz,Ig)
                     case (4)
                        ADF_Lx(Ixy,Iz,Ig)=Buff_ADF_Lx(Ixy,Iz,Ig)
                        ADF_Ly(Ixy,Iz,Ig)=Buff_ADF_Ly(Ixy,Iz,Ig)
                        ADF_Rx(Ixy,Iz,Ig)=Buff_ADF_Rx(Ixy,Iz,Ig)
                        ADF_Ry(Ixy,Iz,Ig)=Buff_ADF_Ry(Ixy,Iz,Ig)
                     end select
                  endif
               enddo
            enddo
         enddo
      enddo

      if (.not.Flag_RefDF) return

      do Iz=IzFuelBot,IzFuelTop
         do Ixy_RF=1,Nxy_RF
            Ixy=I_RF(Ixy_RF)
            do Ig=1,N_Group
               if (I_Lx(Ixy) /= 0) then
#ifdef tuan_fr_crm                  
                  if (I_FARF_1N(I_4Nto1N(I_Lx(Ixy)))==2) then
#else
                  if (I_FARF_1N(I_4Nto1N(I_Lx(Ixy)))>=2) then
#endif
                     ADF_Lx(Ixy,Iz,Ig)=ADF_Lx(Ixy,Iz,Ig)*ADF_Rx(I_Lx(Ixy),Iz,Ig)
                  endif
               endif
               if (I_Rx(Ixy) /= 0) then
#ifdef tuan_fr_crm                  
                  if (I_FARF_1N(I_4Nto1N(I_Rx(Ixy)))==2) then
#else
                  if (I_FARF_1N(I_4Nto1N(I_Rx(Ixy)))>=2) then
#endif
                     ADF_Rx(Ixy,Iz,Ig)=ADF_Rx(Ixy,Iz,Ig)*ADF_Lx(I_Rx(Ixy),Iz,Ig)
                  endif
               endif
               if (I_Ly(Ixy) /= 0) then
#ifdef tuan_fr_crm                  
                  if (I_FARF_1N(I_4Nto1N(I_Ly(Ixy)))==2) then
#else
                  if (I_FARF_1N(I_4Nto1N(I_Ly(Ixy)))>=2) then
#endif
                     ADF_Ly(Ixy,Iz,Ig)=ADF_Ly(Ixy,Iz,Ig)*ADF_Ry(I_Ly(Ixy),Iz,Ig)
                  endif
               endif
               if (I_Ry(Ixy) /= 0) then
#ifdef tuan_fr_crm                  
                  if (I_FARF_1N(I_4Nto1N(I_Ry(Ixy)))==2) then
#else
                  if (I_FARF_1N(I_4Nto1N(I_Ry(Ixy)))>=2) then
#endif
                     ADF_Ry(Ixy,Iz,Ig)=ADF_Ry(Ixy,Iz,Ig)*ADF_Ly(I_Ry(Ixy),Iz,Ig)
                  endif
               endif
            enddo
         enddo
      enddo

      return
      end subroutine Get_ADF


      SUBROUTINE Get_POW
      USE Inc_3D, ONLY: Flux, Power
      USE Inc_maXS, ONLY: kap_maXS_f_3D
      USE Inc_Option, ONLY: N_Group
#ifdef tuan_fr
      USE Inc_TPEN, ONLY: fluxf,  kap_maXS_f_3D_MG, hflxf,imapsol
!      USE Mod_TPENDrive, ONLY: p2r
#endif ! commented @$^ Siarhei_FR
    !  use Inc_DecayHeat, only: use_decay_heat
!    !  use Mod_DecayHeat, only: update_decay_heat
!    !  use Mod_DecayHeat, only: power_plus_decayheat
      IMPLICIT NONE
      INTEGER :: Ixy, Ixy_FA, Iz, Ig

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_POW] in Mod_GetSome'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_POW] in Mod_GetSome'
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

!!!#ifdef tuan_fr
!!!        real(8) :: value1(1:nxy,1:nz,1:N_Group)
!!!
!!!
!!!#endif ! commented @$^ Siarhei_FR

      Power = D0

!     ! if (use_decay_heat) call update_decay_heat
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 

      DO Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)

         DO Iz = IzFuelBot, IzFuelTop
            DO Ig = 1, 2 ! N_Group
               Power(Ixy, Iz) = Power(Ixy, Iz) +  &
                  kap_maXS_f_3D(Ixy, Iz, Ig) * Flux(Ixy, Iz, Ig)
!                  kap_maXS_f_3D_MG(Ixy, Iz, Ig) &
!                   * hflxf(Ig,imapsol(Ixy),Iz)
!                  *Fluxf(Ixy, Iz, Ig) 
!                  *Fluxf(Ixy, Iz, Ig) 

            END DO

         END DO
      END DO

!     ! if (use_decay_heat) call power_plus_decayheat
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 

      RETURN
      END SUBROUTINE Get_POW


      SUBROUTINE Get_LinPOW

      USE Inc_3D, ONLY: Power
      USE Inc_TH, ONLY: LinPOW

      IMPLICIT NONE

      INTEGER :: Ixy, Ixy_FA, Iz

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_LinPOW] in Mod_GetSome'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_LinPOW] in Mod_GetSome'
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

      DO Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)

         DO Iz = IzFuelBot, IzFuelTop
            LinPOW(Ixy, Iz) = Power(Ixy, Iz) * NodeVolume(Ixy, Iz) / h_z(Iz) * DP2
            ! DP2    :: W/cm   => W/m
         END DO
      END DO

      RETURN
      END SUBROUTINE Get_LinPOW


#ifdef siarhei_delete 
      SUBROUTINE Get_FoldAvg(Avg, Input, I_N, I_R, FlagRot, Flag)
      USE Inc_INP
      IMPLICIT NONE
      REAL(8), DIMENSION(:, :, :), INTENT(IN) :: Input
      REAL(8), DIMENSION(:), INTENT(OUT) :: Avg
      INTEGER, INTENT(IN) :: Flag, FlagRot, I_N, I_R
      INTEGER :: Iz

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_FoldAvg] in Mod_GetSome'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      Avg = D0

      DO Iz = 1, Nz
         SELECT CASE(Flag)
         CASE(10)
            Avg(Iz) = Input(I_N, I_R, Iz)
         CASE(7)
            SELECT CASE(FlagRot)
            CASE(1)
               Avg(Iz) = 0.5_8 * (Input(I_N, I_R, Iz) + Input(I_N, I_R - 2, Iz))
            CASE(2)
               Avg(Iz) = 0.5_8 * (Input(I_N, I_R, Iz) + Input(I_N, I_R + 1, Iz))
            CASE(3)
               Avg(Iz) = 0.5_8 * (Input(I_N, I_R, Iz) + Input(I_N, I_R - 1, Iz))
            CASE(4)
               Avg(Iz) = 0.5_8 * (Input(I_N, I_R, Iz) + Input(I_N, I_R + 2, Iz))
            END SELECT
         CASE(6)
            SELECT CASE(FlagRot)
            CASE(1)
               Avg(Iz) = 0.5_8 * (Input(I_N, I_R, Iz) + Input(I_N, I_R - 1, Iz))
            CASE(2)
               Avg(Iz) = 0.5_8 * (Input(I_N, I_R, Iz) + Input(I_N, I_R - 2, Iz))
            CASE(3)
               Avg(Iz) = 0.5_8 * (Input(I_N, I_R, Iz) + Input(I_N, I_R + 2, Iz))
            CASE(4)
               Avg(Iz) = 0.5_8 * (Input(I_N, I_R, Iz) + Input(I_N, I_R + 1, Iz))
            END SELECT
         CASE(4)
            Avg(Iz) = Sum(Input(I_N, 1:4, Iz)) / D4
         END SELECT
      END DO

      RETURN
      END SUBROUTINE Get_FoldAvg
#endif 


      subroutine get_normal_power
      use inc_3d, only: power
      use inc_3d, only: avg_power
      use inc_3d, only: normal_power
      implicit none
      integer :: ixy_fa, ixy, iz

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [get_normal_power] in Mod_GetSome'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_normal_power] in Mod_GetSome'
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

      call get_avg(avg_power,power,0)
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         do iz=izfuelbot,izfueltop
            normal_power(ixy,iz)=power(ixy,iz)/max(1d-10,avg_power)
         enddo
      enddo

      return
      end subroutine get_normal_power


      subroutine get_burnup(opt)
      use inc_time, only: dt
      use inc_rp, only: i_farf_1n
      use inc_3d, only: mtu
      use inc_3d, only: power
      use inc_3d, only: BU
      use inc_3d, only: BU_old
      use inc_3d, only: BU_Predictor
      use inc_crud, only: opt_crud
      use inc_crud, only: ind_bu
      use inc_crud, only: ind_bu_old
      use inc_crud, only: ind_bu_predictor
      use inc_depletion, only: PC_w
      use inc_depletion, only: opt_pc
      implicit none
      logical(1), intent(in) :: opt
      integer :: ixy, ixy_1n, iz

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [get_burnup] in Mod_GetSome'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_burnup] in Mod_GetSome'
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

      if (dT<1d-10) return

      select case (opt_pc)
      case (0) ! original Full PC
         do ixy=1,nxy
            ixy_1n=i_4nto1n(ixy)
            if (i_farf_1n(ixy_1n)==1) cycle
            do iz=izfuelbot,izfueltop
               if (.not.opt) then
                  bu_predictor(ixy,iz)=bu(ixy,iz) &
                            & +power(ixy,iz)*nodevolume(ixy,iz) &
                            & *1d-9*(dT/86400d0)/mtu(ixy,iz)/2d0
                  bu(ixy,iz)=bu(ixy,iz) &
                            & +power(ixy,iz)*nodevolume(ixy,iz) &
                            & *1d-9*(dT/86400d0)/mtu(ixy,iz)

               else
                  bu(ixy,iz)=bu_predictor(ixy,iz) &
                            & +power(ixy,iz)*nodevolume(ixy,iz) &
                            & *1d-9*(dT/86400d0)/mtu(ixy,iz)/2d0
 
                endif
            enddo
         enddo

      case (1) ! No PC
         do ixy=1,nxy
            ixy_1n=i_4nto1n(ixy)
            if (i_farf_1n(ixy_1n)==1) cycle
            do iz=izfuelbot,izfueltop
               bu(ixy,iz)=bu_old(ixy,iz) &
                         & +power(ixy,iz)*nodevolume(ixy,iz) &
                         & *1d-9*(dT/86400d0)/mtu(ixy,iz)

            enddo
         enddo

      case (2) ! Semi PC
         do ixy=1,nxy
            ixy_1n=i_4nto1n(ixy)
            if (i_farf_1n(ixy_1n)==1) cycle
            do iz=izfuelbot,izfueltop
               bu(ixy,iz)=bu_old(ixy,iz) &
                         & +power(ixy,iz)*nodevolume(ixy,iz) &
                         & *1d-9*(dT/86400d0)/mtu(ixy,iz)

            enddo
         enddo

         if (.not.opt) then !PREDICTOR
            BU = PC_w*BU_Predictor + (1d0-PC_w)*BU

         else !CORRECTOR - save BU_predictor
            BU_Predictor = BU

         endif

      case (3) ! Full PC new
         do ixy=1,nxy
            ixy_1n=i_4nto1n(ixy)
            if (i_farf_1n(ixy_1n)==1) cycle
            do iz=izfuelbot,izfueltop
               bu(ixy,iz)=bu_old(ixy,iz) &
                         & +power(ixy,iz)*nodevolume(ixy,iz) &
                         & *1d-9*(dT/86400d0)/mtu(ixy,iz)

            enddo
         enddo

         if (.not.opt) then !PREDICTOR
            BU = PC_w*BU_Predictor + (1d0-PC_w)*BU

         else !CORRECTOR - save BU_predictor
            BU_Predictor = BU

         endif
      end select
      return
      end subroutine get_burnup


      subroutine LevelBalancing
      use Inc_3D, only: Flux, FisSrc, FisSrc_Iout, Power
      use Inc_TH, only: LevFactor
#ifdef tuan_fr
      use Inc_TPEN, only: fluxf,imapsol,hflxf, aflx
#endif
#ifdef tuan_fr_tdep
      use Inc_TPEN, only: tfluxf
#endif
      use inc_option, only: n_group
      implicit none
      integer :: Ixy, Iz, Ig
#ifdef tuan_fr_tdep
      integer :: It
#endif

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [LevelBalancing] in Mod_GetSome'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [LevelBalancing] in Mod_GetSome'
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

      do Ixy = 1, Nxy
         if ( I_FARF_1N( I_4Nto1N(Ixy) ) == 1 ) cycle
         do Iz = IzFuelBot, IzFuelTop
            Power(Ixy, Iz) = Power(Ixy, Iz) * LevFactor
            FisSrc     (Ixy, Iz) = FisSrc     (Ixy, Iz) * LevFactor
            FisSrc_Iout(Ixy, Iz) = FisSrc_Iout(Ixy, Iz) * LevFactor
         enddo
      enddo

     do Ig = 1, N_Group
          do Iz = 1, Nz
              do Ixy = 1, Nxy
               Flux(Ixy, Iz, Ig) = Flux(Ixy, Iz, Ig) * LevFactor

#ifdef tuan_fr
               Fluxf(Ixy, Iz, Ig) = hflxf(Ig,imapsol(Ixy),Iz) * LevFactor
#endif

#ifdef tuan_fr_tdep
              if(.not. allocated(tfluxf) ) allocate(tfluxf( 1:6,1: Nxy, 1:Nz, 1:N_Group))
            do It = 1,6

                  tfluxf(It,Ixy,Iz,Ig) = aflx(Ig,It,imapsol(Ixy),iz) * LevFactor
            enddo
#endif
        end do
         enddo
     enddo

      return
      end subroutine LevelBalancing

      END MODULE Mod_GetSome
