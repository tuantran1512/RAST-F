      PROGRAM RAST_K
      use Mod_Watch
      use Inc_Option, only: OPT_Mode
      use Mod_Manage, only: Initializing , Set_Memory ! check Set_Memory later
      use Mod_InitCond, only: IniCond_Nodal, IniCond_TH
    !  use Mod_InitCond, only: jump_in
      use mod_charedit, only: print_msg
      use Mod_Finalize, only: Finalize_All
      use Read_I, only: Read_In
      use Run_Steady, only: Search_Eigval !, Depletion
#ifdef tuan_fr
      use Run_Steady, only: Depletion_Hex
#endif
#ifdef siarhei_tr_hex
      !use Run_Steady , only: Search_CR_Criticality
#endif
#ifdef tuan_tr_test
      use Run_Transient, only: Transient_hex ! @$^ needs modification Siarhei_FR
#endif 
    !  use Run_Steady, only: LoadFollow_QS, Depletion_QS
!!      use Run_Transient, only: Transient ! @$^ needs modification Siarhei_FR
      use Write_O, only: Write_MAP, Write_SUM, Write_OUT  !, Write_ANC ! tmp_removal
!!      use mod_snf, only: write_snf
      use inc_parallel, only: comm
      use mod_parallel, only: init_comm
#ifdef js_mpi
      use mod_parallel, only: final_parallel
#endif

#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction
     ! use Mod_PinPow_Hex , only: deallocate_pin_pow_hex
     ! use Inc_PinPow_Hex
      ! -=-=-=-=-=-=-=-=-
#endif

#ifdef jr_vver
      use inc_option, only: opt_nodal
      use mod_getnode, only: geometry_hex
#endif

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

#ifdef siarhei_plot
      ! @$^ siarhei_fr

      last_saved = ' '
      current_sub = '----'
      open(678,file='called_routines_in_order.dat',status='replace',&
         access='sequential',form='formatted',action='write')

      current_sub = '^#^ Entered [RAST-K] in Main'
      if(trim(last_saved).ne.trim(current_sub)) then
         write(678,'(A)') trim(current_sub)
         last_saved = current_sub
      end if
      ! =-=-=-=-=-=-=-=
#endif

      call init_comm

      call WATCH('Total simulation','START')
      call WATCH('Initialize', 'START')
      call Initializing  ! Initializing for Static Data
      call Read_In       ! Read Input Files
      call Set_Memory    ! Dynamic Allocation for Dynamic Data ! check later (needs detailed revision)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
#ifdef jr_vver
      if (opt_nodal == 4) then
#ifdef tuan_tr_test

#else
         CALL Geometry_hex
#endif         
      endif
#endif
      call IniCond_Nodal ! Set Initial Condition for Nodal Variables ! check later siarhei_rev
          !              dummy_filler = 1 ! @$^ siarhei_plot 

      call IniCond_TH    ! Set Initial Condition for T/H Variables ! check later siarhei_rev
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
!      call jump_in
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif            
      call WATCH('Initialize', 'END')

      select case (opt_mode)
      case (1)                ! EG
         call Search_EigVal   ! Steady State Eigenvalue Search
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
      case (2)                ! ST
!         call Depletion       ! Steady State Depletion Calculation Driver
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif   
      case (3)                ! TR
!         call Depletion       ! Initial Condition Calculation (Steady State)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
#ifdef tuan_tr_test       
            call Depletion_Hex       ! Steady State Depletion Calculation Driver without Critical Search
!            call Search_CR_Criticality
!            call Search_EigVal
            call Transient_hex       ! Steady State Depletion Calculation Driver without Critical Search
#endif
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
      case (4)                ! DP
#ifdef tuan_fr
         call Depletion_Hex       ! Steady State Depletion Calculation Driver without Critical Search
           !            dummy_filler = 1 ! @$^ siarhei_plot 

write(*,*) '*** END DEPLETION ... '
#else
!         call Depletion       ! Steady State Depletion Calculation Driver without Critical Search
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif   
#endif
      case (5)                ! LS
!         call Depletion       ! Steady State Depletion Calculation Driver
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
 
#endif   
!         call LoadFollow_QS   ! Quasi Steady State Load Follow Calculation Driver with Ciritical Search using CR Position
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif    
      case (6)                ! QS
!         call Depletion       ! Steady State Depletion Calculation Driver
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
!         call Depletion_QS    ! Quasi Steady State Depletion Calculation Driver with Ciritical Search using Boron
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
      case (7)                ! PS
!         call Depletion       ! Steady State Depletion Calculation Driver
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
!         !call PowerSearch_QS ! Quasi Steady State Depletion Calculation Driver with Ciritical Search using Power Level
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
      case (8)                ! LT
         !call LoadFollow_Tr  ! Transiet State Load Follow Calculation Driver
      end select

      if (comm%if_master) then
         call WATCH('Write','START')
         call Write_MAP
!        ! call Write_ANC ! commented @$^ Siarhei_FR
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
!         call Write_SUM
         call Write_OUT
!        ! call Write_SNF ! commented @$^ Siarhei_FR
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         call WATCH('Write','END')
      endif


      call print_msg(0,'*** END RAST-K 2.2 ...')
      call WATCH('Total simulation','END')

      if (comm%if_master) then
         call print_msg(0,'  - Elapsed Time')
         do iclock=1,n_clock
            if (time(iclock)>0.0) then
               call print_message_float(watch_name(iclock),time(iclock),'sec')
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
            endif
         enddo
      endif
      call finalize_all
#ifdef js_mpi
!      call final_parallel
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
#endif

#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction
     ! call deallocate_pin_pow_hex
      ! -=-=-=-=-=-=-=-=-
#endif
      close(678) ! @$^ siarhei_fr debugging file
     ! pause ! @$^ Siarhei checking print output
      END PROGRAM RAST_K
