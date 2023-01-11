#ifdef siarhei_delete


      module Mod_Function
      use Inc_Constant
      use Inc_File
      use Inc_Flag
      use Inc_Geometry
      use Inc_RP
      use Inc_Option

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none
      contains

      subroutine search_10ppm
      use inc_xs_file, only: i_bu
      use inc_depletion, only: n_bu
      use inc_th, only: ppm
      use inc_history, only: hist_ppm
      use inc_depletion, only: cycle_bu
      use inc_depletion, only: if_get_interpol
      use inc_depletion, only: if_conv_10ppm
      use inc_depletion, only: targetppm
      use mod_charedit, only: print_msg
      use mod_branch, only: backup_data
      use mod_branch, only: reset_data
      use inc_depletion, only: cycle_day
      use inc_depletion, only: inp_dT
      use inc_depletion, only: inp_dBU
      use inc_depletion, only: inp_ppower
      use inc_depletion, only: inp_core_power
      use inc_depletion, only: inp_fa_power
      use inc_depletion, only: tot_mtu
      implicit none
      real(8) :: x1,x2,y1,y2,dx,dy

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [search_10ppm] in Mod_Function'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [search_10ppm] in Mod_Function'
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

      if (abs(ppm-targetppm)<1d-2) then
         if_conv_10ppm=.true.
      else
         if (.not. if_get_interpol) then
            if ((PPM-targetppm)<0d0) then
               if_get_interpol = .true.
            else
               if (I_BU >= N_BU) then
                  call backup_data
                  if (I_BU == 1) then
                     Cycle_BU(I_BU+1) = Cycle_BU(I_BU) + 0.5d0
                     call print_msg(0,'* Find BU for ',targetppm,'ppm CBC')
                     call print_msg(0,'----------------------------------------------')
                     call print_msg(0,' > Next Cycle BU :', Cycle_BU(I_BU+1),'GWd/t')
                     call print_msg(0,'----------------------------------------------')
                  else
                     x1=Cycle_BU(I_BU-1)
                     x2=Cycle_BU(I_BU)
                     y1=Hist_PPM(I_BU-1)
                     y2=Hist_PPM(I_BU)
                     dx=x2-x1
                     dy=y2-y1
                     Cycle_BU(I_BU+1) = (targetppm+dy/dx*x1-y1)/(dy/dx)
                     call print_msg(0,'* Find BU for ',targetppm,'ppm CBC')
                     call print_msg(0,'----------------------------------------------')
                     call print_msg(0,' Old:',x1,'GWd/t,',y1,'ppm')
                     call print_msg(0,' New:',x2,'GWd/t,',y2,'ppm')
                     call print_msg(0,' > Next Cycle BU :', Cycle_BU(I_BU+1),'GWd/t')
                     call print_msg(0,'----------------------------------------------')
                  endif
                  N_BU=I_BU+1
               endif
            endif
         endif
         if (if_get_interpol) then
            x1=Cycle_BU(I_BU-1)
            x2=Cycle_BU(I_BU)
            y1=Hist_PPM(I_BU-1)
            y2=Hist_PPM(I_BU)
            dx=x2-x1
            dy=y2-y1
            Cycle_BU(I_BU) = (targetppm+dy/dx*x1-y1)/(dy/dx)
            call print_msg(0,'* Find BU for ',targetppm,'ppm CBC')
            call print_msg(0,'----------------------------------------------')
            call print_msg(0,' Old:',x1,'GWd/t,',y1,'ppm')
            call print_msg(0,' New:',x2,'GWd/t,',y2,'ppm')
            call print_msg(0,' > Next Cycle BU :', Cycle_BU(I_BU),'GWd/t')
            call print_msg(0,'----------------------------------------------')
            call reset_data
            I_BU=I_BU-1
         endif
      endif

      inp_ppower     (I_BU+1) = inp_ppower(I_BU)
      inp_core_power (I_BU+1) = inp_core_power(I_BU)
      inp_fa_power   (I_BU+1) = inp_fa_power(I_BU)
      inp_dBU        (I_BU+1) = cycle_BU(I_BU+1)-cycle_BU(I_BU)
      inp_dT         (I_BU+1) = inp_dBU(I_BU+1)*tot_mtu*86400d0*1d+9/inp_core_power(I_BU+1)
      cycle_day      (I_BU+1) = cycle_day(I_BU)+inp_dT(I_BU+1)/86400d0

      return
      end subroutine search_10ppm


      end module Mod_Function


#endif
