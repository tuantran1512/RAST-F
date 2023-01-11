!#ifdef siarhei_delete

!#elif tuan_tr_test    

      MODULE Mod_Interface

      USE Inc_Constant
      USE Inc_Geometry
      USE Inc_RP
      use mod_charedit, only: print_msg
	  use ieee_arithmetic ! @$^ siarhei_fr for using nvfortran (ieee_is_nan is not working)
#ifdef js_mpi
      use inc_parallel, only: comm
#endif

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

      CONTAINS

      subroutine backup_TH
      use inc_3d, only: t_mod_old_xsfb, t_mod
      use inc_3d, only: d_mod_old_xsfb, d_mod
      use inc_3d, only: t_fuel_old_xsfb, t_fuel
#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif
      T_mod_Old_XSFB =T_mod
      d_mod_Old_XSFB =d_mod
      T_fuel_Old_XSFB=T_fuel

      return
      end subroutine backup_TH


      subroutine check_chth(opt_damp)
      use inc_3d, only: t_mod_old_xsfb, t_mod
      use inc_3d, only: d_mod_old_xsfb, d_mod
      use inc_3d, only: t_fuel_old_xsfb, t_fuel
      implicit none
      logical(1),intent(in) :: opt_damp
      real(8) :: maxerr
      integer(4) ::  ixy_fa,ixy,iz
      real(8),parameter :: eps_del_t=0.1d0 ! 0.5 K
      real(8) :: damp=0.4d0
#ifdef js_r2mpi

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [check_chth] in Mod_Interface'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.comm%usempi.and..not.comm%if_master) return
#endif
      if (.not.opt_core==1) return
      if (opt_damp) then
         damp=0.4d0
      else
         damp=1.0d0
      endif
      maxerr=0d0
      DO Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)
         DO Iz = 1, Nz
            maxerr=max(maxerr,abs(T_mod(ixy,iz)-T_mod_Old_XSFB(ixy,iz)))
            maxerr=max(maxerr,abs(T_fuel(ixy,iz)-T_fuel_Old_XSFB(ixy,iz)))
         enddo
      enddo
      if(maxerr>eps_del_t) then
          !write(*,'(a,f10.2)') '* Updated TH. Max dT:', maxerr

          T_mod = (1d0 - damp) * T_mod_Old_XSFB  + damp * T_mod
          d_mod = (1d0 - damp) * d_mod_Old_XSFB  + damp * d_mod
          T_fuel= (1d0 - damp) * T_fuel_Old_XSFB + damp * T_fuel
      else
          !write(*,'(a,f10.2)') '* Skip TH.    Max dT:', maxerr
          T_mod = T_mod_Old_XSFB
          d_mod = d_mod_Old_XSFB
          T_fuel= T_fuel_Old_XSFB
      endif
      return
      end subroutine check_chth

      subroutine damp_th(opt_damp)
      USE Inc_3D, only: t_fuel, t_fuel_old_xsfb
      USE Inc_3D, only: t_mod, t_mod_old_xsfb
      USE Inc_3D, only: d_mod, d_mod_old_xsfb
      implicit none
      logical(4),intent(in) :: opt_damp
      real(8) :: maxerr
      integer(4) ::  ixy_fa,ixy,iz
      real(8),parameter :: eps_del_t=0.1d0 ! 0.5 K
      real(8) :: damp=0.4d0
#ifdef js_r2mpi

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [damp_th] in Mod_Interface'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.comm%usempi.and..not.comm%if_master) return
#endif
      if (opt_damp) then
         damp=0.4d0
      else
         damp=1.0d0
      endif
      maxerr=0d0
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         do iz=1,nz
            maxerr=max(maxerr,abs(t_mod(ixy,iz)-t_mod_old_xsfb(ixy,iz)))
            maxerr=max(maxerr,abs(t_fuel(ixy,iz)-t_fuel_old_xsfb(ixy,iz)))
         enddo
      enddo
      if (maxerr>eps_del_t) then
          t_mod = (1d0-damp)*t_mod_old_xsfb  + damp*t_mod
          d_mod = (1d0-damp)*d_mod_old_xsfb  + damp*d_mod
          t_fuel= (1d0-damp)*t_fuel_old_xsfb + damp*t_fuel
      else
          t_mod = t_mod_old_xsfb
          d_mod = d_mod_old_xsfb
          t_fuel= t_fuel_old_xsfb
      endif

      return
      end subroutine damp_th

      subroutine add_crud_dtf
      USE Inc_3D
      use inc_Crud, only: cruddtf
      implicit none
      integer(4) ::  ixy_fa,ixy,iz


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [add_crud_dtf] in Mod_Interface'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      DO Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)
         DO Iz = IzFuelBot, IzFuelTop
            T_fuel(ixy,iz)=T_fuel(ixy,iz)+cruddtf(ixy,iz)
         enddo
      enddo

      end subroutine add_crud_dtf


      subroutine interface_th
      use inc_xs_file, only: if_branch
      use inc_3d, only: t_fuel, t_mod, d_mod
      use inc_th, only: tcool, dcool, tdopl
      use inc_th, only: tm_in, tm_out, tfmax, tfuel
      use inc_th, only: ppower
      use mod_getsome, only: get_avg
      use inc_depletion, only: avgrho, flag_b10dep
      use inc_flag, only: flag_tmfb, flag_tffb

#ifdef js_mpi
      use inc_parallel, only: iproc
      use inc_parallel, only: ixy_faip2ixy_fa, iproc2nxy_fa, ixy_fa2iproc
      use mod_parallel, only: barrier, allreduce
#endif
      implicit none
      integer :: ixy, ixy_fa, ixy_rf, iz
      integer :: ith_bk
      real(8) :: tmpvol
#ifdef js_mpi
      real(8), allocatable :: tfmax_array(:)
#endif
#ifdef js_r2mpi

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [interface_th] in Mod_Interface'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      t_fuel=0d0
      t_mod=0d0
      d_mod=0d0

      tfmax=0d0
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         ith_bk=ixy_fa
#ifdef js_mpi
         if (comm%usempi) then
            if (iproc/=ixy_fa2iproc(ixy_fa)) cycle
            ith_bk=ixy_faip2ixy_fa(ixy_fa,iproc)
         endif
#endif
         do iz=izfuelbot,izfueltop
            if (tfmax<tfuel(1,iz,ith_bk)) then
               tfmax=tfuel(1,iz,ith_bk)
            endif
         enddo
      enddo
#ifdef js_mpi
      if (comm%usempi) then
         allocate(tfmax_array(1:comm%n_procs))
         tfmax_array=0d0
         tfmax_array(iproc)=tfmax
         call allreduce(tfmax_array,comm%n_procs)
         tfmax=maxval(tfmax_array)
         deallocate(tfmax_array)
      endif
#endif

      tm_out=0d0
      tmpvol=0d0
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         ith_bk=ixy_fa
#ifdef js_mpi
         if (comm%usempi) then
            if (iproc/=ixy_fa2iproc(ixy_fa)) cycle
            ith_bk=ixy_faip2ixy_fa(ixy_fa,iproc)
         endif
#endif
         tm_out=tm_out+tcool(nz,ith_bk)*nodevolume(ixy,nz)
         tmpvol=tmpvol+nodevolume(ixy,nz)
      enddo
#ifdef js_mpi
      if (comm%usempi) then
         call allreduce(tm_out)
         call allreduce(tmpvol)
      endif
#endif
      tm_out=tm_out/max(1d-10,tmpvol)

      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         ith_bk=ixy_fa
#ifdef js_mpi
         if (comm%usempi) then
            if (iproc/=ixy_fa2iproc(ixy_fa)) cycle
            ith_bk=ixy_faip2ixy_fa(ixy_fa,iproc)
         endif
#endif
         if (flag_tffb) then
            do iz=izfuelbot,izfueltop
               t_fuel(ixy,iz)=tdopl(iz,ith_bk)*tdopl(iz,ith_bk)
            enddo
         endif
         if (flag_tmfb) then
            do iz=1,nz
               t_mod(ixy,iz)=tcool(iz,ith_bk)+degtok
               d_mod(ixy,iz)=dcool(iz,ith_bk)*1d-3
            enddo
         endif
      enddo
      do ixy_rf=1,nxy_rf
         ixy=i_rf(ixy_rf)
         ith_bk=nxy_fa+1
#ifdef js_mpi
         if (comm%usempi) then
            if (.not.comm%if_master) cycle
            ith_bk=iproc2nxy_fa(iproc)+1
         endif
#endif
         if (flag_tffb) then
            do iz=izfuelbot,izfueltop
               t_fuel(ixy,iz)=0d0
            enddo
         endif
         if (flag_tmfb) then
            do iz=1,nz
               t_mod(ixy,iz)=tcool(iz,ith_bk)+degtok
               d_mod(ixy,iz)=dcool(iz,ith_bk)*1d-3
            enddo
         endif
      enddo

#ifdef js_mpi
      if (comm%usempi) then
         call barrier
         call allreduce(t_fuel,nxy,nz)
         call allreduce(t_mod,nxy,nz)
         call allreduce(d_mod,nxy,nz)
      endif
#endif

      if (flag_b10dep) then
         call get_avg(avgrho,d_mod,0)
      endif

#ifdef js_mpi
      if (comm%usempi) call barrier
#endif

      if (tfmax>2865d0) then
         call print_msg(2,'Fuel melting (Centerline Temp. > 2865C)')
      endif
      if ((if_branch/=4).and.(tm_in<280d0.or.tm_out>340d0)) then
         call print_msg(2,'Coolant Temp. is not hot state (280C < Tm < 340C)')
      endif

      ! check solution is Nan
      if (ppower>1d-3) then
         do iz=1,nz
            do ixy=1,nxy
               if (ieee_is_nan(t_fuel(ixy,iz))) then
                  call print_msg(3,'Solution diverged ... t_fuel')
                  call print_msg(3,'Error occur at ixy,iz',ixy,iz)
                  stop
               endif
               if (ieee_is_nan(t_mod(ixy,iz))) then
                  call print_msg(3,'Solution diverged ... t_mod')
                  call print_msg(3,'Error occur at ixy,iz',ixy,iz)
                  stop
               endif
            enddo
         enddo
      endif

      return
      end subroutine interface_th


      subroutine restore_force_th
      use inc_inp, only: flag_force_ftctr
      use inc_inp, only: flag_force_mtctr
      use inc_3d, only: t_fuel, t_fuel_bk, t_fuel_ss
      use inc_3d, only: t_mod, t_mod_bk, t_mod_ss
      use inc_geometry, only: nxy, nz
      use mod_alloc
      implicit none

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [restore_force_th] in Mod_Interface'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (flag_force_ftctr) then
         call alloc(t_fuel_bk,nxy,nz)
         t_fuel_bk=t_fuel
         call alloc(t_fuel_ss,nxy,nz)
         t_fuel_ss=t_fuel
      endif
      if (flag_force_mtctr) then
         call alloc(t_mod_bk,nxy,nz)
         t_mod_bk=t_mod
         call alloc(t_mod_ss,nxy,nz)
         t_mod_ss=t_mod
      endif
      return
      end subroutine restore_force_th

      subroutine adjust_force_th
      use inc_inp, only: flag_force_ftctr, force_ftctr, force_dftctr
      use inc_inp, only: flag_force_mtctr, force_mtctr, force_dmtctr
      use inc_3d, only: t_fuel, t_fuel_bk, t_fuel_ss
      use inc_3d, only: t_mod, t_mod_bk, t_mod_ss
      use inc_geometry, only: nxy, nz
      implicit none
      integer :: l,k

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [adjust_force_th] in Mod_Interface'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (flag_force_ftctr) then
         t_fuel_bk=t_fuel
         do k=1,nz
            do l=1,nxy
               t_fuel(l,k)=(t_fuel_bk(l,k)-t_fuel_ss(l,k)) &
                         &*force_ftctr+force_dftctr+t_fuel_ss(l,k)
            enddo
         enddo
      endif
      if (flag_force_mtctr) then
         t_mod_bk=t_mod
         do k=1,nz
            do l=1,nxy
               t_mod(l,k)=(t_mod_bk(l,k)-t_mod_ss(l,k)) &
                         &*force_mtctr+force_dmtctr+t_mod_ss(l,k)
            enddo
         enddo
      endif
      return
      end subroutine adjust_force_th

      subroutine return_force_th
      use inc_inp, only: flag_force_ftctr
      use inc_inp, only: flag_force_mtctr
      use inc_3d, only: t_fuel, t_fuel_bk
      use inc_3d, only: t_mod, t_mod_bk
      implicit none

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [return_force_th] in Mod_Interface'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (flag_force_ftctr) then
         t_fuel=t_fuel_bk
      endif
      if (flag_force_mtctr) then
         t_mod=t_mod_bk
      endif
      return
      end subroutine return_force_th

#ifdef js_mpc
      subroutine mp_backup_th
      use inc_th, only: tfuel, tfuel_bk
      use inc_th, only: tdopl, tdopl_bk
      use inc_th, only: tcool, tcool_bk
      use inc_th, only: dcool, dcool_bk
      use inc_th, only: hcool, hcool_bk
      use inc_th, only: nr
      use inc_th, only: nzth
      use inc_th, only: nchan
      use inc_3d, only: t_fuel, t_fuel_old_xsfb
      use inc_3d, only: t_mod, t_mod_old_xsfb
      use inc_3d, only: d_mod, d_mod_old_xsfb
      use mod_alloc
      implicit none
#ifdef js_r2mpi

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mp_backup_th] in Mod_Interface'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.comm%usempi.and..not.comm%if_master) return
#endif
      if (.not.allocated(tfuel_bk)) allocate(tfuel_bk(nr+5,nzth,nchan))
      if (.not.allocated(tdopl_bk)) allocate(tdopl_bk(0:nzth,0:nchan+1))
      if (.not.allocated(tcool_bk)) allocate(tcool_bk(0:nzth,0:nchan+1))
      if (.not.allocated(dcool_bk)) allocate(dcool_bk(0:nzth,0:nchan+1))
      if (.not.allocated(hcool_bk)) allocate(hcool_bk(nzth,nchan+1))
      tfuel_bk=tfuel
      tdopl_bk=tdopl
      tcool_bk=tcool
      dcool_bk=dcool
      hcool_bk=hcool
      t_mod_old_xsfb=t_mod
      d_mod_old_xsfb=d_mod
      t_fuel_old_xsfb=t_fuel
      return
      end subroutine mp_backup_th
#endif
#ifdef js_mpc
      subroutine mp_return_th
      use inc_th, only: tfuel, tfuel_bk
      use inc_th, only: tdopl, tdopl_bk
      use inc_th, only: tcool, tcool_bk
      use inc_th, only: dcool, dcool_bk
      use inc_th, only: hcool, hcool_bk
      use inc_3d, only: t_fuel, t_fuel_old_xsfb
      use inc_3d, only: t_mod, t_mod_old_xsfb
      use inc_3d, only: d_mod, d_mod_old_xsfb
      implicit none
#ifdef js_r2mpi

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mp_return_th] in Mod_Interface'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.comm%usempi.and..not.comm%if_master) return
#endif
      tfuel=tfuel_bk
      tdopl=tdopl_bk
      tcool=tcool_bk
      dcool=dcool_bk
      hcool=hcool_bk
      t_mod=t_mod_old_xsfb
      d_mod=d_mod_old_xsfb
      t_fuel=t_fuel_old_xsfb
      return
      end subroutine mp_return_th
#endif

      END MODULE Mod_Interface

!#endif
