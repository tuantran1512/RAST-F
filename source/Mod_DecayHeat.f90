#ifdef siarhei_delete

      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! This module contains functions for Decay Heat calculation =-=
      ! option for Transient scenarios. If enabled, it calculates =-=
      ! current decay heat fraction of power and adds it to the =-=-=
      ! actual value of power, by that replacing some of fission  =-=
      ! power.  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


      module Mod_DecayHeat
      use Inc_DecayHeat

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

      contains

      ! set initial values for Decay_Heat variables
      subroutine initialize_decay_heat(Nxy,Nz,N_Group)
      use Inc_DecayHeat
      implicit none
      integer                             :: Nxy, Nz, N_Group

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [initialize_decay_heat] in Mod_DecayHeat'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [initialize_decay_heat] in Mod_DecayHeat'
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

      if (.not.allocated(Flux_dh_old)) allocate(Flux_dh_old(Nxy,Nz,N_Group))
      if (.not.allocated(Decay_Power)) allocate(Decay_Power(Nxy,Nz))
      if (.not.allocated(dep_total_DH)) allocate(dep_total_DH(Nxy,Nz))
      if (.not.allocated(DH_precurs)) allocate(DH_precurs(Nxy,Nz,N_Group,6))
      tn          =  NONE
      Flux_dh_old =  NONE
      tn_Old      =  NONE
      Decay_Power =  NONE
      DH_precurs =   NONE
      C_tn =         NONE
      C_tn_Old =     NONE
      dC_tn =        NONE
      dC_tn_Old =    NONE
      dP_tn =        NONE
      dP_tn_Old =    NONE
      curr_time_dh = NONE

      if (beta_user == 0.0) beta_user = beta_dunn

      end subroutine initialize_decay_heat

      subroutine update_decay_heat()
      use Inc_DecayHeat
      use Inc_Depletion, only: Tot_MTU
      use Inc_Geometry, only: NodeVolume,Nxy,Nz,Opt_Core
      use Inc_maXS, only: kap_maXS_f_3D
      use Inc_3D, only: Flux
      use Inc_Option, only: N_Group
      use Inc_Transient,  only: Flag_Transient
      use Inc_Time, only: sec,dt_Tr

      implicit none
      integer        :: i,k,l,m

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [update_decay_heat] in Mod_DecayHeat'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [update_decay_heat] in Mod_DecayHeat'
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

      if (init_decay_heat) then
         call initialize_decay_heat(Nxy,Nz,N_Group)
         write(*,*) '*** Decay Heat Initialized ***'
         write(*,*) '******************************'
         init_decay_heat = .false.
      end if
      curr_time_dh_old = curr_time_dh
      if (Tot_MTU==0) then
         Tot_MTU = 63.568
      end if
      curr_time_dh = init_time_dh
      ! curr_time_dh = Cycle_Bu(1)*Tot_MTU*24*3600/Core_Power_100*1.D9
      ! if (Cycle_Bu(1).eq.0) curr_time_dh = 600
      Decay_Power = NONE
      total_DH = NONE
      if (.not.Flag_Transient) DH_precurs = NONE
      do k = 1, Nxy
         do l = 1, Nz
            do m=1,N_Group
               if (kap_maXS_f_3D(k,l,m).ne.0) then
                  do i=1,6
                     C_tn_Old(i) = DH_precurs(k,l,m,i)
                  end do
                  if (Flag_Transient) then
                     tn = sec
                     if ((sec - dt_Tr) < 0) then
                        tn_Old = 0
                     else
                        tn_Old = tn - dt_Tr
                     end if
                  else
                     tn = curr_time_dh
                     tn_Old = curr_time_dh_old
                  end if
                  P_dec_heat = NONE
                  call get_decay_heat(Flux(k,l,m),Flux_dh_old(k,l,m),kap_maXS_f_3D(k,l,m))
                  Flux_dh_old(k,l,m) = Flux(k,l,m)
                  do i=1,6
                     DH_precurs(k,l,m,i) = C_tn(i)
                     C_tn(i) = NONE
                  end do
                  Decay_Power(k,l) = Decay_Power(k,l) + P_dec_heat
               end if
            end do
            total_DH = total_DH + Decay_Power(k,l) * NodeVolume(k,l)
         end do
      end do
      if (.not.Flag_Transient) curr_time_dh = NONE
      if (OPT_Core == 4) total_DH  = total_DH * 4
      if (.not.Flag_Transient) dep_total_DH = Decay_Power

      end subroutine update_decay_heat

      ! subroutine for updating the current value of fission power
      ! only needed for plotting as an auxiliary plotting data
      subroutine update_fission_power
      use Inc_DecayHeat
      use Inc_3D, only: Power
      use Inc_Geometry, only: NodeVolume,Nxy,Nz

      implicit none
      integer                             :: k, l

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [update_fission_power] in Mod_DecayHeat'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [update_fission_power] in Mod_DecayHeat'
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

      total_FP = NONE
      do k = 1, Nxy
         do l = 1, Nz
            total_FP = total_FP + (Power(k,l)-Decay_Power(k,l))*NodeVolume(k,l)
         end do
      end do

      end subroutine update_fission_power

      ! main subroutine for updating precursors and DH power (1 node)
      subroutine get_decay_heat(Flux,Flux_Old,kappainp)
      use Inc_DecayHeat

      implicit none
      real(8)                             :: Flux, Flux_Old
      integer                             :: i
      real(8),intent(in)                  :: kappainp
      real(8)                             :: kappa

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [get_decay_heat] in Mod_DecayHeat'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_decay_heat] in Mod_DecayHeat'
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

      ! Update concentrations of precursors, find DH power fractions
      dC_tn_Old      = dC_tn
      dP_tn_Old      = dP_tn
      P_dec_heat     = NONE
      kappa = kappainp / (1-beta_user)
      do i = 1,6
         C_tn(i) = C_tn_Old(i)*DEXP(-1*lam6_user(i)*(tn-tn_Old)) + &
            (kappa*0.5*(Flux+Flux_Old)*beta6_user(i) / &
            lam6_user(i))*(1-DEXP(-1*lam6_user(i)*(tn-tn_Old)))
         if (kappa.ne.0) then
            dC_tn(i) = C_tn(i) / (kappa*Flux*beta6_user(i)/ &
               lam6_user(i))
            dP_tn(i) = C_tn(i)*lam6_user(i)
         else
            dC_tn(i) = 0.0
            dP_tn(i) = 0.0
         end if
         P_dec_heat = P_dec_heat + dP_tn(i)
      end do

      end subroutine get_decay_heat

      ! add the value of decay heat to Power(Nxy,Nz)
      subroutine power_plus_decayheat()
      use Inc_DecayHeat, only: Decay_Power,use_decay_heat
      use Inc_Geometry, only: Nxy,Nz
      use Inc_3D, only: Power
      implicit none
      integer        :: k,l

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [power_plus_decayheat] in Mod_DecayHeat'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [power_plus_decayheat] in Mod_DecayHeat'
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

      if (use_decay_heat) then
         do k = 1, Nxy
            do l = 1, Nz
               Power(k,l) = Power(k,l) + Decay_Power(k,l)
            end do
         end do
      end if

      end subroutine power_plus_decayheat

      ! remove the value of decay heat from Power(Nxy,Nz)
      subroutine power_minus_decayheat()
      use Inc_DecayHeat, only: Decay_Power,use_decay_heat
      use Inc_Geometry, only: Nxy,Nz
      use Inc_3D, only: Power
      implicit none
      integer        :: k,l

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [power_minus_decayheat] in Mod_DecayHeat'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [power_minus_decayheat] in Mod_DecayHeat'
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

      if (use_decay_heat) then
         do k = 1, Nxy
            do l = 1, Nz
              Power(k,l) = Power(k,l) - Decay_Power(k,l)
            end do
         end do
      end if

      end subroutine power_minus_decayheat

      ! concentrations and Decay Heat Power for the depletion step that is
      ! considered as the beginning of given transient calculation.
      subroutine get_initial_decay_heat(time_0,time_1,Nxy, Nz, N_Group)
      use Inc_DecayHeat
      use Inc_3D,                      only: Flux, Flux_Old
      use Inc_maXS,                    only: kap_maXS_f_3D
      implicit none
      integer                             :: time_0 ! zero time of cycle
      integer                             :: time_1 ! time of transient calc begin
      integer                             :: d_time = 1 ! delta t = 1 sec
      integer                             :: t
      integer                             :: k,l,m
      integer                             :: Nxy, Nz, N_Group
      integer                             :: i

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [get_initial_decay_heat] in Mod_DecayHeat'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_initial_decay_heat] in Mod_DecayHeat'
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

      do t = (time_0 + d_time),time_1,d_time
         tn = t
         tn_Old = t - d_time

         Decay_Power = NONE
         total_DH = NONE
         do k = 1, Nxy
            do l = 1, Nz
               do m=1,N_Group
                  if (kap_maXS_f_3D(k,l,m).ne.0) then
                     do i=1,6
                        C_tn_Old(i) = DH_precurs(k,l,m,i)
                     end do
                     P_dec_heat = NONE
                     call get_decay_heat(Flux(k,l,m),Flux_Old(k,l,m),kap_maXS_f_3D(k,l,m))
                     do i=1,6
                        DH_precurs(k,l,m,i) = C_tn(i)
                        C_tn(i) = NONE
                     end do
                     Decay_Power(k,l) = Decay_Power(k,l) + P_dec_heat
                  end if
               end do
            end do
         end do
      end do

      end subroutine get_initial_decay_heat

      ! deallocate arrays and return kappa to the original state
      ! (without Decay Heat)
      subroutine quit_decay_heat
      use Inc_DecayHeat
      implicit none

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [quit_decay_heat] in Mod_DecayHeat'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [quit_decay_heat] in Mod_DecayHeat'
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

      use_decay_heat = .false.

      if (allocated(Flux_dh_old)) deallocate(Flux_dh_old)
      if (allocated(Decay_Power)) deallocate(Decay_Power)
      if (allocated(dep_total_DH)) deallocate(dep_total_DH)
      if (allocated(DH_precurs)) deallocate(DH_precurs)

      write(*,*) '*** Decay Heat Data Deallocated ***'
      write(*,*) '***********************************'

      end subroutine quit_decay_heat

      end module Mod_DecayHeat


#endif
