#ifdef siarhei_delete


      module Mod_SSnodal
      use Inc_Constant
      use Inc_Flag
      use Inc_Geometry
      use Inc_RP
      use Inc_Option
      use inc_parallel, only: comm
#ifdef js_mpi
      use mod_parallel, only: bcast, barrier
#endif

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

      contains

      subroutine NodalSolver_St
      use Inc_maXS
      use Inc_FluxVar
      use Mod_Wielandt
      use Mod_SolLS
      use Inc_Lscoef
      use Mod_SolLU
!      use Mod_Soldhat
      use Mod_THdriver
      use Mod_SolFDM
      use Mod_InitCond, only: GuessIteration
!      use Mod_Interface
      use Inc_3D, only: Flux, Power, Avg_Power, keff
      use Inc_Control
      use Inc_Option
      use mod_charedit, only: print_msg
      use Inc_TH, only: Avg_CorePower, LevFactor, tfmax
      use Inc_XS_File, only: I_BU, flag_leakage
!      use Mod_LeakCorr, only: XSFB_LC
      use Mod_XSFB, only: XSFB
      use mod_getsome, only: levelbalancing
      use mod_getsome, only: get_asi
      use mod_getsome, only: get_avg
      use mod_getsome, only: get_fa_asi
      use mod_getsome, only: get_linpow
      use mod_getsome, only: get_pow
      use mod_getsome, only: get_source
      use mod_getsome, only: get_normal_power
!      use Mod_Depletion, only: FissionProduct
      use Inc_Crud, only: opt_crud, opt_crudasi, opt_crudcbc, opt_crudfaasi
      use Inc_XYZ, only: ASI
!      use Mod_Depletion, only: Crud_Itr_ASIsearch,Crud_itr_CBCCOR,Crud_Itr_FAASIsearch
      use inc_inp, only: flag_force_asi, force_asi, sigcor, flag_force_ao
      use inc_branch, only: flag_stack
      use inc_extsrc, only: flag_extsrc, extsrc
!      use mod_extsrc, only: set_fsrc_extsrc, write_extsrc
      use inc_inp, only: flag_ftnth_detail, switch_saveth
      use inc_detector, only: flag_det_pow_sol, DET_BU
      use inc_depletion, only: cycle_bu
      use inc_flag, only: flag_rod
      use inc_cr, only: cr_bot, rodmove
#ifdef js_mpc
!      use link_th1d, only: run_th1d_solver
      use inc_multiphysics, only: flag_mpc
      use inc_multiphysics, only: opt_mpc
!      use mod_pinpow, only: pinpow_s
#ifdef js_ctf
!      use link_ctf, only: run_ctf_solver
#endif
#ifdef js_frap
!      use link_frap, only: run_fcn_solver
#endif
#endif
      implicit none
      integer :: Iout
      integer :: k, l, m, ioutbeg, ioutend, nskipnodal
      integer :: icy_CMFD, icy_TH, iin
      real(8) :: reigvdel, fs
      save ioutbeg,nskipnodal
      data ioutbeg,nskipnodal/1,1/
      real(8) :: keff_0,keff_1,ppm_0,ppm_1,err_old
      integer(4) :: i_xs_upd
      real(8) :: ppm_old=0d0, ppm_old2=0d0, ppm_tmp=0d0
      integer(4) :: nupd_ppm
      real(8) :: k_old=0d0, k_old2=0d0
      real(8) :: dkdppm=0d0
      logical(1) :: ifconv_ppm = .false.
      logical(1) :: ifoscil_ppm = .false.
      real(8) :: eps_ppm
      real(8) :: avg_normal_power_old=0d0
      real(8) :: avg_normal_power_old2=0d0
      real(8) :: dpdppm
      real(8) :: dppm
      real(8) :: dppm_max=5d0
      logical(1) :: ifupdate_dhat=.false.
      logical(1) :: flag_crudasi_conv
      real(8) :: dampfac
      logical(1) :: flag_conv_asi_first=.true.
      real(8) :: asi_ori, asi_cal, asi_inp
      real(8) :: mulfac, mulfac_old
      integer :: check_even
      logical(1) :: exit_signal
      integer :: th_signal=0
      real(8) :: tmp_tm_avg=0.0d0
      integer(4) :: i_tavg_itr=0
      logical(1) :: flag_tavg_conv
      integer :: nupd_cr
      real(8) :: hcr_max
      real(8) :: hcr_min
      real(8) :: dcr
      real(8) :: cr_old
      real(8) :: cr_old2
      real(8) :: dkdcr=0d0
      logical(1) :: ifconv_cr
#ifdef js_mpc
      logical(1) :: flag_mpcall
      integer :: nmp_precond=0
      integer :: nmp_call
#endif
      ioutbeg=1
      Ioutend=Iout_Max
#ifdef js_mpc
      nmp_call=0
#endif

#ifdef js_r2mpi
     call barrier
     if (comm%if_master) then
#endif
      call XSFB

      nupd_ppm = 1
      dnw = 0d0
      dnn = 0d0
      dnb = 0d0
      call GuessIteration
      call FdmCoupl
      call setls
      call ilufac2d
      ifnodal  = .true.               ! true: need nodal calculation update
      nodal    = .true.
      iflsupd  = .false.              ! true: need update linear system
      icy_CMFD = - Period_CMFD
      icy_TH   = - Period_TH
      flagl2   = .false.
      flaglinf = .false.
      flageig  = .false.
      flagerf  = .false.
#ifdef js_mpc
      flag_mpcall=.true.
      if (flag_mpc) flag_mpcall=.false.
#endif
      eps_ppm=1d-2
      !if (opt_core==1) eps_ppm=1d-1
      ifconv_ppm=.false.
      Flag_Conv_PowLev=.true.
      dampfac=1d+2
      !if (opt_core==1) dampfac=1d+3

      if (flag_rod) then
         ifconv_cr=.false.
         nupd_cr=1
         hcr_min=sum(h_z(1:izfuelbot-1))
         hcr_max=sum(h_z(1:nz))
         dcr=10d0
         cr_old=0d0
      else
         ifconv_cr=.true.
      endif


         flag_crudasi_conv=.true.

      if (opt_findtavg.and.flag_thfb.and..not.flag_doing_branch) then
         i_tavg_itr=0
         itr_tavg(:)=tm_in
         itr_fa_flowrate(:)=fa_massflow
         flag_tavg_conv=.false.
         if(opt_findtavg_init.and.i_bu>1) then
           flag_tavg_conv=.true.
         endif
      else
         flag_tavg_conv=.true.
      endif

      keff_0=keff
      keff_1=keff_0
      ppm_0=ppm
      ppm_1=ppm_0
      i_xs_upd=0
      err_old=1.0

      ! set initial dppm
      dppm = dppm_max

#ifdef js_r2mpi
     endif
#endif

      do Iout=Ioutbeg,Ioutend
         exit_signal=.false.
#ifdef js_r2mpi
        call barrier
        if (comm%if_master) then
#endif
         icy_CMFD=icy_CMFD+1
         if (icy_CMFD>=Period_CMFD) icy_CMFD=mod(icy_CMFD,Period_CMFD)
         icy_TH=icy_TH+1
         if (icy_TH>=Period_TH) icy_TH=mod(icy_TH,Period_TH)
         iflsupd=.false.  ! ls update judgment
         ifnlupd=.false.  ! non-linear system update judgment
         reigvdel=reigv-reigvs

         call Get_Source(FisSrc,nu_maXS_f_3D,Flux)
         if (OPT_Mode<7) then
            do k=1,Nz
               do l=1,Nxy
                  fs=FisSrc(l,k)*reigvdel
                  do m=1,N_Group
                     src(m,l,k)=maxs_chi_3d(l,k,m)*fs
                  enddo
                  if (flag_extsrc) src(1,l,k)=src(1,l,k)+extsrc(l,k)
               enddo
            enddo
         else  ! FSP
            do k=1,Nz
               do l=1,Nxy
                  src(1,l,k)=0d0
                  src(2,l,k)=0d0
               enddo
            enddo
         endif
         call WATCH('SS CMFD','START')

         if (flag_det_pow_sol.and.(DET_BU-Cycle_BU(I_bu)<1d-10)) then
            call SolveLinearSystem_det(iin,.FALSE.)
         elseif (flag_det_pow_sol)then
            flag_det_pow_sol=.false.
         else
            call SolveLinearSystem(iin,.FALSE.)
         endif
         call WATCH('SS CMFD','END')

         call wiel(icy_CMFD)

         ! output iteration table
         if (.not.flag_stack) then
            if (comm%if_master) write(*,601) iout,iin,keff,flageig,rerrl2,flagl2,errlinf,flaglinf,domr
         endif

         ! exit outer if convergence achieved
         !damping
         if (ifconv_ppm.and.flageig.and..not.flag_extsrc) then
            if (rerrl2<=EPS_Global*dampfac) flagl2=.true.
            if (errlinf<=EPS_Local*dampfac) flaglinf=.true.
         endif
         flagneut=flageig.and.flagl2.and.flaglinf
         if (flag_extsrc) flagneut=flagneut.and.ifconv_ppm.and.ifupdate_dhat
#ifdef js_mpc
         flagneut=flagneut.and.flag_mpcall
#endif
         skip=.false.

         if (.not.nodal.or.icy_CMFD==1) then
            if (flagneut.and.Flag_Conv_PowLev.and.flag_crudasi_conv.and.flag_tavg_conv) then
               if (opt_mode==2) then
                  if (abs(keff-k_target)<EPS_Critical) then
                     skip=.true.
                     exit_signal=.true.
                  endif
               elseif (opt_mode/=1.and.flag_rod) then
                  if (ifconv_cr) then
                     skip=.true.
                     exit_signal=.true.
                  endif
               else
                  if (flag_force_asi) then
                     call get_asi (asi_cal, power)
                     if (flag_conv_asi_first) then
                        asi_ori=asi_cal
                        flag_conv_asi_first=.false.
                        check_even=0
                        if (mod(nz,2)==1) check_even=1
                        mulfac_old=1d0
                        mulfac=mulfac_old
                     endif
                     asi_inp=force_asi
                     if (abs(asi_cal-asi_inp)>1d-3.and.iout<=1000) then
                        ! mac_sig_res_abs correction
                        mulfac_old=mulfac
                        mulfac=1d0-(asi_cal-asi_inp)/10d0
                        mulfac=mulfac*mulfac_old
                        do k=1,nz/2
                           sigcor(k)=2d0-mulfac
                        enddo
                        do k=nz/2+1+check_even,nz
                           sigcor(k)=mulfac
                        enddo
                        call xsfb
                        call setls
                        call ilufac2d
                        if (flag_force_ao) then
                           call print_msg(0,'* Correct for target AO  : ',-asi_cal,'  => ',-asi_inp)
                        else
                           call print_msg(0,'* Correct for target ASI : ', asi_cal,'  => ', asi_inp)
                        endif
                     else
                        if (iout>1000) then
                           call print_msg(0,'* Reach limit AO/ASI search iteration')
                        endif
                        if (flag_force_ao) then
                           call print_msg(0,'* Initial AO  : ',-asi_ori,'  Final AO  : ',-asi_cal,'  Target AO  : ',-asi_inp)
                        else
                           call print_msg(0,'* Initial ASI : ', asi_ori,'  Final ASI : ', asi_cal,'  Target ASI : ', asi_inp)
                        endif
                        skip=.true.
                        exit_signal=.true.
                     endif
                  else
                     skip=.true.
                     exit_signal=.true.
                  endif
               endif
            endif
         endif
#ifdef js_r2mpi
        endif
        call barrier
        call bcast(exit_signal,comm%i_master)
#endif
         if (exit_signal) exit

#ifdef js_r2mpi
        if (comm%if_master) then
#endif
         ! invoke nodal update
         ifxsupd=.false.
         ifnlupd=(icy_CMFD==0).or.(iout==3).or.(flagneut)
         if (nodal.and.ifnlupd.and.ifnodal) then
            if (.not.flag_extsrc.or.rerrl2<1d-3) then
               call WATCH('SS Nodal','START')
               call updatedhat
               call WATCH('SS Nodal','END')
               iflsupd=.true.
               ifupdate_dhat=.true.
            endif
         endif

         ! invoke T/H update
         ifnlupd=(icy_TH==0).or.(iout==3)
         th_signal=0
         if (ifnlupd) then
            call Get_POW
            call Get_Avg(Avg_Power,Power,0)
            if (flag_extsrc) then
               flag_conv_powlev=flaglinf
            else
               LevFactor = Avg_CorePower/Avg_Power
               call LevelBalancing
               !js+if ((abs(keff_Old-k_target)<0.001d0).and.(abs(LevFactor-D1)>EPS_PowLev)) then
               if (abs(LevFactor-1d0)>EPS_PowLev) then
                  Flag_Conv_PowLev=.false.
               else
                  Flag_Conv_PowLev=.true.
               endif
            endif
            call Get_POW
            call Get_LinPOW
            call get_normal_power
            if (flag_thfb) then
#ifdef js_mpc
               if (flag_mpc) then
                  ! precondition
                  if (nmp_precond<30) then
                     if (flageig.and.flagl2.and.flaglinf) then
                        nmp_precond=30
                     else
                        nmp_precond=nmp_precond+1
                     endif
                     th_signal=1
                     goto 1921
                  endif
                  if (iout<20.or.iout>300) goto 1921
                  nmp_call=nmp_call+1
                  th_signal=opt_mpc+1
                  flag_mpcall=.true.
1921              continue
               else
                  th_signal=1
                  if (.not.flag_th_original) then
                     th_signal=th_signal+10
                     nmp_call=nmp_call+1
                  endif
               endif
#else
               th_signal=1
#endif
               ifxsupd=.true.  ! Invoke XSFB
            endif
         endif

         if (flag_ftnth_detail) switch_saveth=.true.
#ifdef js_r2mpi
        endif
        call barrier
        call bcast(th_signal,comm%i_master)
#endif
         select case (th_signal)
         case (1)
            call WATCH('TH Feedback','START')
            call ssth
            call backup_TH
            call Interface_TH
            call check_chth(flag_tavg_conv) !Damping
            call WATCH('TH Feedback','END')
#ifdef js_mpc
         case (11)
            if (flag_th1d_pbp) then
               call WATCH('Pin Power Recon','START')
               call pinpow_s
               call WATCH('Pin Power Recon','END')
            endif
            call WATCH('TH Feedback','START')
            call run_th1d_solver
            call print_msg(0,' Run Next TH1D  ... Done', nmp_call)
            call backup_TH
            call Interface_TH
            call check_chth(flag_tavg_conv) !Damping
            call WATCH('TH Feedback','END')
#ifdef js_ctf
         case (2)
            call WATCH('Pin Power Recon','START')
            call pinpow_s
            call WATCH('Pin Power Recon','END')
            call print_msg(0,' Run Next CTF  ... ', nmp_call)
            call WATCH('CTF','START')
            call run_ctf_solver()
            call WATCH('CTF','END')
            call print_msg(0,' Run Next CTF  ... Done', nmp_call)
            call WATCH('TH Feedback','START')
            call backup_TH
            call Interface_TH
            call damp_th(.true.)
            call WATCH('TH Feedback','END')
#endif
#ifdef js_frap
         case (3)
            call WATCH('Pin Power Recon','START')
            call pinpow_s
            call WATCH('Pin Power Recon','END')
            call print_msg(0,' Run Next FRAPCON  ... ', nmp_call)
            call WATCH('FRAPI','START')
            call run_fcn_solver(.false.,2,1)
            call WATCH('FRAPI','END')
            call print_msg(0,' Run Next FRAPCON  ... Done', nmp_call)
            call WATCH('TH Feedback','START')
            call backup_TH
            call Interface_TH
            call damp_th(.true.)
            call WATCH('TH Feedback','END')
#ifdef js_ctf
         case (4)
            call WATCH('Pin Power Recon','START')
            call pinpow_s
            call WATCH('Pin Power Recon','END')
            call print_msg(0,' Run Next CTF  ... ', nmp_call)
            call WATCH('CTF','START')
            call run_ctf_solver()
            call WATCH('CTF','END')
            call print_msg(0,' Run Next CTF  ... Done', nmp_call)
            call print_msg(0,' Run Next FRAPCON  ... ', nmp_call)
            call WATCH('FRAPI','START')
            call run_fcn_solver(.false.,1,2)
            call WATCH('FRAPI','END')
            call print_msg(0,' Run Next FRAPCON  ... Done', nmp_call)
            call WATCH('TH Feedback','START')
            call backup_TH
            call Interface_TH
            call damp_th(.true.)
            call WATCH('TH Feedback','END')
#endif
#endif
#endif
         end select

#ifdef js_r2mpi
        if (comm%if_master) then
         if (flag_ftnth_detail) switch_saveth=.false.
#endif
         if (ifnlupd) then
            if (OPT_Mode/=4.and.opt_mode/=1.and..not.flag_rod) then
               if (rerrl2<1d-2) then
                  if (flag_extsrc.and.rerrl2<1d-3) then
                     if (nupd_ppm==1) then
                        ppm_old=ppm
                        avg_normal_power_old=avg_power/avg_corepower
                        ppm=ppm+dppm_max
                     else
                        ppm_old2=ppm_old
                        ppm_old=ppm
                        avg_normal_power_old2=avg_normal_power_old
                        avg_normal_power_old=avg_power/avg_corepower
                        dpdppm=(avg_normal_power_old-avg_normal_power_old2)/(ppm_old-ppm_old2)
                        dppm=(1d0-avg_normal_power_old)/dpdppm
                        if (abs(dppm)>dppm_max) then
                           if (dppm>0d0) then
                              ppm=ppm_old+dppm_max
                           else
                              ppm=ppm_old-dppm_max
                           endif
                        else
                           ppm=ppm_old+dppm
                        endif
                     endif
                    ! check osillation & reduce dppm_max (KHNP CRI)
                    if ((PPM-ppm_old==dppm_max).and.(ppm_old-ppm_old2==-dppm_max).and.(nupd_ppm>3)) then
                       ifoscil_ppm=.true.
                    else
                       ifoscil_ppm=.false.
                    endif
                    if (ifoscil_ppm) then
                       dppm_max = dppm_max * 0.5d0
                       if (comm%if_master) write(*,'(a,f10.2,a)') '* Decrease Max dPPM: ', dppm_max,' ppm'
                    endif
                     ifconv_ppm=.false.
                     if (abs(ppm-ppm_old)<1d-2) ifconv_ppm=.true.
                     if (comm%if_master) write(*,'(a,f10.2,a)') &
                        & '* Update boron concentration: ', ppm, ' ppm'
                     if (comm%if_master) write(*,'(a,2es15.5)') &
                        & '* dPPM/ dPdPPM: ', dppm, dpdppm
                     if (comm%if_master) write(*,'(a,2es15.5)') &
                        & '* Power / Target [W]: ',avg_power*tot_fuelvol, avg_corepower*tot_fuelvol
                  else
                     if (.not.flag_extsrc) then
                        ppm_old2=ppm_old
                        ppm_old=ppm
                        k_old2=k_old
                        k_old=keff
                        ppm_tmp=PPM+(keff-k_target)*1.0d+4
                        if ( (abs(ppm_tmp-ppm_old)>abs(ppm_old-ppm_old2))     &
                           & .and. ((ppm_tmp-ppm_old)*(ppm_old-ppm_old2)<0d0) &
                           & .and. (nupd_ppm>3) ) then
                           ifoscil_ppm=.true.
                        else
                           ifoscil_ppm=.false.
                        endif
                        if (ifoscil_ppm) then
                           call print_msg(1,'Oscillation occurred during CBC search')
                           dkdppm=(k_old-k_old2)/(ppm_old-ppm_old2)
                           PPM=(k_target-k_old+dkdppm*ppm_old)/dkdppm
                        else
                           PPM=PPM+(keff-k_target)*1.0d+4
                        endif
                        ifconv_ppm=.false.
                        if (abs(ppm-ppm_old)<eps_ppm) ifconv_ppm=.true.
                        call print_msg(0,'* Update boron concentration: ', PPM,' ppm')
                     endif
                  endif
                  ifxsupd=.true.  ! Invoke XSFB
                  nupd_ppm=nupd_ppm+1
               endif
               if (opt_findtavg.and.rerrl2<1d-3) then
                  call Get_Avg(tmp_tm_avg,T_Mod,0)
                  call change_massflow(tmp_tm_avg,i_tavg_itr,flag_tavg_conv)
               endif

            elseif (rerrl2<1d-4.and.opt_mode/=1.and.flag_rod) then
               if (.not.ifconv_cr) then
                  cr_old2=cr_old
                  cr_old=cr_bot(rodmove)
                  k_old2=k_old
                  k_old=keff
                  if (nupd_cr==1) then
                     dcr=10d0
                  elseif (nupd_cr==2) then
                     dcr=5d0
                     if (k_old>k_old2) dcr=-10d0
                  else
                     if (abs(cr_old-cr_old2)>1d-3) then
                        dkdcr=(k_old-k_old2)/(cr_old-cr_old2)
                        dcr=(k_target-k_old)/dkdcr
                     else
                        dcr=0.001
                     endif
                  endif
                  dcr=min(10d0,max(-10d0,dcr))
                  cr_bot(rodmove)=cr_bot(rodmove)+dcr
                  if (cr_bot(rodmove)<hcr_min) then
                     cr_bot(rodmove)=hcr_min
                  elseif (cr_bot(rodmove)>hcr_max) then
                     cr_bot(rodmove)=hcr_max
                  endif
                  call print_msg(0,'* Update control rod',rodmove,' position =',cr_bot(rodmove))
                  if (abs(keff-k_target)<1d-4) ifconv_cr=.true.
                  ifxsupd=.true.
                  nupd_cr=nupd_cr+1
               else
                  if (abs(keff-k_target)>1d-4) ifconv_cr=.false.
               endif
            endif
         endif

         if ((OPT_Xe/=2).or.(OPT_Sm/=2)) then
            !js+call FissionProduct(.false.)
            if (ifnlupd) then
               call FissionProduct(.false.)
               ifxsupd=.true.
            endif
         endif

         ! update xsec and others associated with xsec change
         if (ifxsupd) then
            if (Flag_Leakage) then
               call print_msg(0,'* Leakage Correction *')
               call WATCH('SS Nodal','START')
               call updatedhat
               call WATCH('SS Nodal','END')
               call XSFB_LC
               call FdmCoupl
               call setls
            else
               call XSFB
            endif
            call Get_Source(FisSrc,nu_maXS_f_3D,Flux)
            call FdmCoupl
            iflsupd=.true.
         endif

         if (iflsupd) then
            ! update linear system
            call setls
            ! update preconditioner
            call ilufac2d
         else
            ! correct diagonal components of the coefficient matrix
            reigvdel=reigvs-reigvsd
            do k=1,nz
               do l=1,nxy
                  am(1,l,k)=am(1,l,k)-af(1,l,k)*reigvdel
                  af2(l,k)=af(2,l,k)*reigvs
               enddo
            enddo
            reigvsd=reigvs
         endif
#ifdef js_r2mpi
        endif
#endif
      enddo

      call get_asi(asi,power)
      if (flag_extsrc) call write_extsrc

#ifdef js_r2mpi
     call barrier
#endif

601   format(i7,i4,f12.5,l2,2(1p,e13.4,l2),0p,f8.4,f10.2)
      return
      end subroutine NodalSolver_St


      end module Mod_SSnodal


#endif
