!#ifdef siarhei_delete
#ifdef tuan_tr_test

      MODULE Mod_TRnodal
      use Inc_Solver
      use Inc_Control
      use Mod_SolFDM
      use Mod_Manage
      use Mod_SolLS
      use Mod_SolLU
      use Mod_THdriver
      use Inc_FluxVar
      use Inc_Lscoef, only : af
      use Inc_maXS, only: maxs_chid_3d
      use Inc_Adjoint
      use Inc_Constant
      use Inc_Geometry
      use Inc_Flag
      use Inc_RP
      use Inc_3D
      use Inc_FA
      use Inc_INP
      use Inc_Kinetics
      use Inc_Option
      use Inc_TH, only: PPower
      use Inc_Time
      use Inc_Transient
      use Mod_GetSome, only: Get_Source

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

      CONTAINS

      subroutine settran
      implicit none


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [settran] in Mod_TRnodal'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      cetak    = 0.5d0 ! central differencing as the default
      cetaf    = 0.5d0 ! theta for heat conduction in fuel
      cetac    = 1.0d0 ! theta for heat convection in coolan
      nordprec = 2d0   ! second order precursor integration

      return
      end subroutine settran


      subroutine inittran
      use Inc_3D
      use Inc_TH, only: ppower_0
      use Inc_CR
      use Inc_Kinetics
      use inc_mixs, only: miXS_a_Xe35, miXS_a_Sm49
      use inc_nuclide, only: N_Xe35, N_Sm49
      implicit none
      real(8) :: tmp1, vol
      real(8) :: psi0, sumd
      integer :: k, l, m, ik, icomp, Ixy, Ixy_FA, Iz, Ig_d
      integer :: nxyfb, nzfb


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [inittran] in Mod_TRnodal'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      rhoadj = D0
      call Alloc_FBRHO

      ! assign transient iteration control parameters
      Period_TH    = 5          ! Number of T/H updates per nodal update 5
      Period_CMFD  = 1          ! Number of iterations to be performed before the first nonlinear update 1
      Period_NonL  = 5          ! Nonlinear update cycle 5

      ! assign nodewise kinetics parameters
      if (Flag_Card_maXS) then
         beta_d_Tot=0d0
         do k=1,Nz
            do l=1,Nxy
               icomp=I_Comp(l,k)
               ik=kincomp(icomp)
               do m=1,N_Group
                  v_Inv(l,k,m)=1d0/tvelo(m,ik)
               enddo
               do Ig_d=1,N_Group_d
                  beta_d_Tot(l,k)=beta_d_Tot(l,k)+tbeta(Ig_d,ik)
                  beta_d_eff(l,k,Ig_d)=tbeta(Ig_d,ik)
                  lambda_d(l,k,Ig_d)=lambda_d_I(Ig_d,ik)
               enddo
            enddo
         enddo
      else ! .not. flag_card_maxs
         call alloc_tran
         do Ixy_FA=1,Nxy_FA
            Ixy=I_FA(Ixy_FA)
            do Iz=IzFuelBot,IzFuelTop
               do Ig_d=1,N_Group_d
                  beta_d_Tot(Ixy,Iz)=beta_d_Tot(Ixy,Iz)+beta_d_eff(Ixy,Iz,Ig_d)
               enddo
            enddo
         enddo
      endif

      if (flag_card_maxs) then
         if (flag_force_beff) then
            beta_d_tot=0d0
            do m=1,n_group_d
               do k=1,nz
                  do l=1,nxy
                     beta_d_eff(l,k,m)=beta_d_eff(l,k,m)*force_beff(m)
                     beta_d_tot(l,k)=beta_d_tot(l,k)+beta_d_eff(l,k,m)
                  enddo
               enddo
            enddo
         endif
         if (flag_force_lamb) then
            do m=1,n_group_d
               do k=1,nz
                  do l=1,nxy
                     lambda_d(l,k,m)=lambda_d(l,k,m)*force_lamb(m)
                  enddo
               enddo
            enddo
         endif
         if (flag_force_velo) then
            do m=1,n_group
               do k=1,nz
                  do l=1,nxy
                     v_inv(l,k,m)=v_inv(l,k,m)*force_velo(m)
                  enddo
               enddo
            enddo
         endif
      endif

      ! initialize frequently used constants
      reigv=1d0/keff
      cetak0=cetak
      cetak=1d0
      cetakb=1d0-cetak
      cetakr=cetakb/cetak
      cetafb=1d0-cetaf
      cetacb=1d0-cetac
      cetadelt=dT
      if (Flag_THFB) then
         cetacr = cetacb / cetac
         kgap   = cetaf  * h_Gap * delr
         kgapb  = cetafb * h_Gap * delr
         tmp1   = h_Gap  * h_Clad * (4d0-h_Clad/R_Gap)*R_Fuel / R_Gap
         kgap2  = cetaf  * h_Gap * h_Clad * R_Fuel / R_Gap
         kgap4  = cetaf  * tmp1
         kgap4b = cetafb * tmp1
      endif

      ! calculate steady state precursor concentration
      af           = reigv * af
      nu_maXs_f_3D = reigv * nu_maXs_f_3D

      call Get_Source(FisSrc,nu_maXS_f_3D,Flux)
      FisSrc_Old  = FisSrc
      FisSrc_OOld = FisSrc

      do Ixy_FA=1,Nxy_FA
         Ixy=I_FA(Ixy_FA)
         do Iz=IzFuelBot,IzFuelTop
            do Ig_d=1,N_Group_d
               if (abs(lambda_d(Ixy,Iz,Ig_d))>1d-10) then
                  C_d_I(Ixy,Iz,Ig_d)=beta_d_eff(Ixy,Iz,Ig_d) &
                      &  /lambda_d(Ixy,Iz,Ig_d)*FisSrc_Old(Ixy,Iz)
               endif
            enddo
         enddo
      enddo

      psi0=ppower_0
      sumv0=0d0
      sumd=0d0
      do k=1,Nz
         do l=1,Nxy
            vol=NodeVolume(l,k)
            sumd=sumd+FisSrc(l,k)*flux_adj(l,k,1)
            do m=1,N_Group
               sumv0=sumv0+Flux(l,k,m)*v_inv(l,k,m)*flux_adj(l,k,m)*vol
            enddo
         enddo
      enddo
      sumd=1d0/sumd
      sumv0=sumv0/psi0

      ! save initial feedback variables (time=0.0)
      nxyfb = nchan
      nzfb  = nzth
      init_tFuel  = T_Fuel
      init_tMod   = T_Mod
      init_dMod   = D_Mod
      init_ppm    = PPM
      if (n_cr>0) init_CR_Bot = Abs_Bot
      if (n_cr>0) init_CR_Top = TopFol_Bot
      init_miXS_a_Xe35 = miXS_a_Xe35
      init_miXS_a_Sm49 = miXS_a_Sm49
      init_N_Xe35 = N_Xe35
      init_N_Sm49 = N_Sm49
      do k=1,Nz
         do l=1,Nxy+Nx+Ny
            do m=1,N_Group
               dnw0(m,l,k)=dnw(m,l,k)
               dnn0(m,l,k)=dnn(m,l,k)
            enddo
         enddo
      enddo
      do k=1,Nz+1
         do l=1,Nxy
            do m=1,N_Group
               dnb0(m,l,k)=dnb(m,l,k)
            enddo
         enddo
      enddo

      keff=1d0
      reigv=1d0
      reigvs=1d0

      return
      end subroutine inittran


      subroutine precurs
      use Inc_RP
      use Inc_Kinetics, only: N_Group_d
      use Inc_maXS, only: nu_maXS_f_3D
      use Inc_Option, only: N_Group
      use Inc_Transient, only: C_d_I
      implicit none
      integer :: k, l, lf, ik, m, i
      real(8) :: vol, psin


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [precurs] in Mod_TRnodal'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif


      do k=IzFuelBot,IzFuelTop
         do lf=1,Nxy_FA
            l=nodef(lf)
            vol=NodeVolume(l,k)
            psin=0d0
            do m=1,N_Group
               psin=psin+nu_maXs_f_3D(l,k,m)*Flux(l,k,m)
            enddo
            psin=psin*vol
            do i=1,N_Group_d
               C_d_I(l,k,i) = C_d_I(l,k,i)    * cappa(l,k,i)  &
                          & + FisSrc_OOld(l,k)* omegam(i,l,k) &
                          & + FisSrc_Old(l,k) * omega0(i,l,k) &
                          & + psin            * omegap(i,l,k)
            enddo
         enddo
      enddo

      return
      end subroutine precurs


      subroutine settfsp(first)
      use Inc_BiCG
      use Inc_Solver
      use Inc_Lscoef, only: af
      use Inc_FixSrc_Tr
      use Inc_maXS
      use Inc_Time
      use inc_extsrc, only: flag_extsrc, extsrc
      implicit none
      integer :: i, k, l, m, ik
      real(8) :: ramprec
      real(8) :: vol
      logical(4) :: first


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [settfsp] in Mod_TRnodal'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif
#ifdef tuan_tr_test
         call settfsp_hex(first)
         return
#endif
!      if (abs(cetak-1d0)>1d-10) then
!         call AxB(Flux,flux_add)
!      endif
!
!      do k=1,Nz
!         do l=1,Nxy
!            if (I_Comp(l,k)==0) cycle
!            ik=kincomp(I_Comp(l,k))
!            vol=NodeVolume(l,k)
!            rvdtdvol(1)=vol*rvdelt_Old(1,l,k)
!            rvdtdvol(2)=vol*rvdelt_Old(2,l,k)
!            do m=1,N_Group
!               if ( (Flux_Old(l,k,m)>0d0).and.(Flux_OOld(l,k,m)>0d0).and.&
!                    (Flux_Old(l,k,m)/Flux_OOld(l,k,m)<20d0) ) then !prevent excessive extrapolation
!                  u_omega(m,l,k)=log(abs(Flux_Old(l,k,m)/Flux_OOld(l,k,m)))/dT_Old
!               else
!                  u_omega(m,l,k)=0d0
!               endif
!               flux_extra(m,l,k)=dexp(u_omega(m,l,k)*dT)
!            enddo
!
!            if (.not.first) then
!               do m=1,N_Group
!                  rvdelt(m,l,k)=v_Inv(l,k,m)/cetadelt+u_omega(m,l,k)*v_Inv(l,k,m)
!                  rvdtvol(m)=vol*v_Inv(l,k,m)/cetadelt
!               enddo
!               if (ltolfa(l)==0.or.k>IzFuelTop.or.k<IzFuelBot) then
!                  cetakrl=0d0
!               else
!                  cetakrl=cetakr
!               endif
!            else
!               do m=1,N_Group
!                  if (ltolfa(l)==0.or.k>IzFuelTop.or.k<IzFuelBot) then
!                     cetakrl=0d0
!                     rcdeltl=1d0/dT
!                  else
!                     cetakrl=cetakr
!                     rcdeltl=1d0/(dT*cetak)
!                  endif
!                  rvdelt(m,l,k)=rcdeltl*v_Inv(l,k,m)
!                  rvdtvol(m)=vol*rvdelt(m,l,k)
!               enddo
!            endif
!
!            omegalm=0d0
!            omegal0=0d0
!            omegalp=0d0
!            if (abs(af(2,l,k))>1d-10) then
!               if (nordprec==1) then
!                  do i=1,N_Group_d
!                     omegam(i,l,k)=0d0
!                     omega0(i,l,k)=beta_d_eff(l,k,i)*(capbrldt(l,k,i)-cappa(l,k,i))
!                     omegap(i,l,k)=beta_d_eff(l,k,i)*(1-capbrldt(l,k,i))
!                     omegalm=0d0
!                     omegal0=omegal0+omega0(i,l,k)
!                     omegalp=omegalp+omegap(i,l,k)
!                  enddo
!               else
!                  do i=1,N_Group_d
!                     omegam(i,l,k)=beta_d_eff(l,k,i)/lambda_d(l,k,i) &
!                        & *rldtgp1(l,k,i)*(2*capbrldt(l,k,i)-gamma*cappap1(l,k,i))
!                     omega0(i,l,k)=beta_d_eff(l,k,i)/lambda_d(l,k,i) &
!                        & *(rldt(l,k,i)*(cappap1(l,k,i)+capbrldt2(l,k,i)*rgamma)-cappa(l,k,i))
!                     omegap(i,l,k)=beta_d_eff(l,k,i)/lambda_d(l,k,i) &
!                        & *(1d0-rldtgp1(l,k,i)*(2d0+capbrldt2(l,k,i)*rgamma))
!                     omegalm=omegalm+omegam(i,l,k)*lambda_d(l,k,i)
!                     omegal0=omegal0+omega0(i,l,k)*lambda_d(l,k,i)
!                     omegalp=omegalp+omegap(i,l,k)*lambda_d(l,k,i)
!                  enddo
!               endif
!            endif
!
!            ! establish effective source
!            ! calculate g*phi term
!            omegalpd=betap_Old(l,k)+beta_d_Tot(l,k)-1d0
!            sd=omegalpd*FisSrc_Old(l,k)
!            flux_add(l,k,1)=flux_add(l,k,1)-rvdtdvol(1) &
!               & *Flux_Old(l,k,1)+maxs_chid_3d(l,k,1)*sd
!            flux_add(l,k,2)=flux_add(l,k,2)-rvdtdvol(2) &
!               & *Flux_Old(l,k,2)+maxs_chid_3d(l,k,2)*sd
!            sd=0d0
!            sdn=0d0
!            do i=1,N_Group_d
!               ramprec=lambda_d(l,k,i)*C_d_I(l,k,i)
!               sd=sd+ramprec
!               sdn=sdn+cappa(l,k,i)*ramprec
!            enddo
!
!            src_tr_fac1(1,l,k)=Flux_Old(l,k,1)*rvdtvol(1) &
!               & +cetakrl*(sd-flux_add(l,k,1))
!            src_tr_fac1(2,l,k)=Flux_Old(l,k,2)*rvdtvol(2) &
!               & -cetakrl*flux_add(l,k,2)
!            src_tr_fac2(1,l,k)=sdn+omegalm*FisSrc_OOld(l,k) &
!               & +omegal0*FisSrc_Old(l,k)
!            src_tr_fac2(2,l,k)=0d0
!            src(1,l,k)=flux_extra(1,l,k)*(src_tr_fac1(1,l,k)-cetakrl*u_omega(1,l,k) &
!               & *v_Inv(l,k,1)*vol*Flux_Old(l,k,1))+src_tr_fac2(1,l,k)
!            src(2,l,k)=flux_extra(2,l,k)*(src_tr_fac1(2,l,k)-cetakrl*u_omega(2,l,k) &
!               & *v_Inv(l,k,2)*vol*Flux_Old(l,k,2))
!            if (flag_extsrc) then
!               src(1,l,k)=src(1,l,k)+extsrc(l,k)*(1d0+cetakrl*flux_extra(1,l,k))
!            endif
!            betap(l,k)=1d0-beta_d_Tot(l,k)+omegalp
!         enddo
!      enddo
!
!      ifnodal=.true.

      return
      end subroutine settfsp


#ifdef tuan_tr_test
      subroutine tr_outer_hex
      use Inc_BiCG
!      use Mod_Soldhat
      use Inc_maXS, only: maXS_a_3D, maXS_r_3D, maXS_s_3D
      use Mod_GetSome
      use Mod_Interface
      use Mod_XSFB
      !use inc_extsrc, only: flag_extsrc
      use Mod_SolTPEN, only: DtilHex, setlshex, &
                             SolveLinearSystem_hex
      use Mod_SStpen,  only: IluFacHex
      use Mod_TPENDrive
      use Inc_TPEN
      implicit none
      integer(4) :: th_signal
      logical(1) :: exit_signal
      logical(1) :: ifth
      logical(1) :: notifth
      logical(1) :: notifnodal
      logical(1) :: chkth
      logical(1) :: chknodal
      logical(1) :: flagnodal
      integer(4) :: nth, nnodal
      integer(4) :: nthmax=10
      integer(4) :: icy, icynodal
      integer(4) :: iout, iin
      integer(4) :: k, l, m

      flagnodal=.false. 

      if (.not.Flag_THFB) flagth=.true.
      ifth=.not.flagth
      notifth=flagth
      nth=0
      if (.not.nodal) ifnodal=.false.
      notifnodal=.not.ifnodal
      nnodal=0

      ! begin loop
      icy      = -Period_CMFD
      icynodal = -Period_CMFD
      
      Iout_Max = 500 ! for test compard with PARCS 
      do Iout=1,Iout_Max
        exit_signal=.false.
         iflsupd=.false.
         ifxsupd=.false.
         ifnlupd=.false.

         call WATCH('TR CMFD','START') 
         call SolveLinearSystem_hex(iin,.true.)
         call WATCH('TR CMFD','END') 

         do m=1,N_Group
            do k=1,Nz
               do l=1,Nxy
                  if (Flux(l,k,m)<=0) then
                     Flux(l,k,m) = 1d-30
                  endif
               enddo
            enddo
         enddo
         
!         call get_fissrc 
         call Get_Source(FisSrc,nu_maXS_f_3D,Flux)

         icy=icy+1
         icynodal=icynodal+1
         if ((icy>0).and.(r2<=EPS_ERF)) flagerf=.true.
         ! check convergence
         chkth=flagth
         if (notifth.or.(nth/=0.and.icy==1)) chkth=.true.
         chknodal=.false.
         if (.not.nodal) chknodal=.true.
         if (nnodal/=0) then
            chknodal=.true.
         endif
         if(flag_thfb) then
             if (chkth.and.chknodal) then
                if (flagr2.and.flagth.and.flagl2.and.flaglinf) exit_signal=.true.
                if (flagr2.and.flagth) exit_signal=.true.
                if (flagr2.and.flagl2.and.flaglinf) exit_signal=.true.
             endif
         else

             if (chknodal) then
                if (flagr2.and.flagth.and.flagl2.and.flaglinf) exit_signal=.true.
                if (flagr2.and.flagth) exit_signal=.true.
!                if (flagr2) exit_signal=.true.
                if (flagr2.and.flagl2.and.flaglinf) exit_signal=.true.
!write(*,*) 'come here:  12   3 2', exit_signal, flagr2, flagl2, flaglinf

             endif
             
         endif
         
         if (exit_signal) exit
!write(*,*) 'come here:  12   4'

         ! invoke non-linear updates
         ifnlupd=(mod(icy,Period_NonL)==0).or.(flagerf).or.(flagr2)
         th_signal=0
         if (ifnlupd) then
            if (flag_thfb) then
               th_signal=1
               if (nth>nthmax) th_signal=0
            endif
            if ((th_signal>0) .and. ifth) then
               call Get_POW
               call Get_LinPOW
               call trtf_hex
               ifxsupd=.true.
               nth=nth+1
               icy=0
            endif
         endif

         if (ifnlupd) then
            ! invoke nodal update
            flagnodal=.false.
            if (nodal.and.(ifnodal.or.((mod(nth,Period_TH)==0).and.(iout>Period_CMFD)))) then
               call WATCH('TR Nodal','START')
               call TpenDriver
               call WATCH('TR Nodal','END')
               if (N_Group>2) ifxsupd=.true.
               nnodal=nnodal+1
               icynodal=0
               flagnodal=.true.
               iflsupd=.true.

            endif

            ! update cross sections
            if (ifxsupd) then
               
               call adjust_force_th
               if (Flag_Card_maXS) then
!                  CALL XSFBHEX_MG
                  CALL maXSFBHEX
               else
                  call XSFBHEX
               endif

               call return_force_th
                call Get_Source(FisSrc,nu_maXS_f_3D,Flux)

!               call get_fissrc 
!                      call Get_Source(FisSrc,nu_maXS_f_3D,Flux)

               call FdmCoupl
               iflsupd=.true.
            endif

            ! setup new linear system
            call setlshex(1)

            call ilufachex

            icy=0
         endif
        if (n_group > 2) call mphihmg
      enddo

      return
      end subroutine tr_outer_hex


      SUBROUTINE mphihmg
      use inc_tpen
      use inc_3d, only: flux
      use inc_maxs
      use inc_geometry, only: gridsize_z, gridsize_x, nodevolume

      !! for multigroup calculation
      IMPLICIT NONE
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:)   :: cntz0i
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:)   :: cntz1i
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:)   :: hflxfn
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:,:) :: cnti
      REAL(8)    :: rrt3h 
      REAL(8)    :: dfdm 
      REAL(8)    :: dhatm 
      REAL(8)    :: bt 
      REAL(8)    :: oneflxm 
      REAL(8)    :: w 
      REAL(8)    :: oneflxn
      REAL(8)    :: phis
      REAL(8)    :: cjn 
      REAL(8)    :: cntoavg 
      REAL(8)    :: ttt 
      REAL(8)    :: flx1 
      REAL(8)    :: flx2 
      REAL(8)    :: cntoavg1
      REAL(8)    :: cntoavg2
      REAL(8)    :: hzr 
      REAL(8)    :: tt1 
      REAL(8)    :: tt2 
      REAL(8)    :: arear 
      REAL(8)    :: vol 
      REAL(8)    :: sumlk 
      REAL(8)    :: sumss 
      INTEGER(4) :: ig 
      INTEGER(4) :: iz 
      INTEGER(4) :: ih 
      INTEGER(4) :: it 
      INTEGER(4) :: nn 
      INTEGER(4) :: isfc 
      INTEGER(4) :: ig2 
      INTEGER(4) :: ind 
      INTEGER(4) :: iz1 
      INTEGER(4) :: iz2
      INTEGER(4) :: id1 
      INTEGER(4) :: id2 
      INTEGER(4) :: neigdn 
      INTEGER(4) :: neigup
      LOGICAL(1) :: ifbcref(2)
      REAL(8)    :: rt3
      REAL(8)    :: hside 
      REAL(8)    :: sqrt3 
      REAL(8)    :: rsqrt3 
      REAL(8)    :: hexarea

      hexarea = gridsize_x(1)*gridsize_x(1)/(3d0**(1d0/2d0))/2d0/2d0*6d0

      IF(.NOT.ALLOCATED(cnti)) THEN
         ALLOCATE(cnti(n_group,ntph),cntz0i(n_group),cntz1i(n_group),hflxfn(n_group))
         cntz0i=0
         cntz1i=0
         hflxfn=0
         cnti=0
      ENDIF

      DO ig=1,2
         ifbcref(ig)=.FALSE.
         IF(alxr(ig).EQ.0) ifbcref(ig)=.TRUE.
      ENDDO

      rt3=1.73205080756888
      sqrt3=SQRT(3.d0)
      rsqrt3=1/sqrt3
      hside=gridsize_x(1)*rsqrt3
      rrt3h=1./rt3/hside

      ! surface flux cal. from phis=w*phin+(1-w)phim+beta*(phim+phin)/2
      DO iz=1,nz
         DO ih=1,nassy
            DO it=1,6
               nn=neignd(it,ih)
               isfc=neigsfc(it,ih)
               DO ig=1,2
                  dfdm=dfd(ig,isfc,iz)
                  dhatm=wtdhat(it,ih)*dhat(ig,isfc,iz)
                  bt=betaphis(ig,isfc,iz)
                  oneflxm=flux(imap(ih),iz,ig)
                  IF(nn.NE.0) THEN
                     w=dfdm/d_3d(imap(ih),iz,ig)*0.5
                     oneflxn=flux(imap(nn),iz,ig)
                     phis=w*oneflxn+(1.-w)*oneflxm+bt*(oneflxm+oneflxn)*0.5
                     cjn=-((dfdm+dhatm)*oneflxn-(dfdm-dhatm)*oneflxm) &
                          /rt3/hside
                  ELSE
                     cjn=(dfdm-dhatm)*oneflxm/rt3/hside
                     phis=bt*oneflxm
                  ENDIF
                  cntoavg=(phis+2.*cjn)*0.25
                  DO ig2=mgb(ig),mge(ig)
                     ttt=fcnto(ig2,it,ih,iz)
                     cnto(ig2,it,ih,iz)=ttt*cntoavg
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
! update z-direction out-current
      IF(nz.GT.1) THEN
         DO iz=1,nz+1
            ind=neigsndz(3,iz)
            IF(ind.EQ.12) THEN
               iz1=neigsndz(1,iz)
               iz2=neigsndz(2,iz)
               id1=neigsndz(4,iz)
               id2=neigsndz(5,iz)
               hzr=gridsize_z(iz2)/gridsize_z(iz1)
               DO ih=1,nassy
                  DO ig=1,2
                     dfdm=dfdz(ig,ih,iz)
                     dhatm=dhatz(ig,ih,iz)
                     bt=betaphisz(ig,ih,iz)
                     flx1=flux(imap(ih),iz1,ig)
                     flx2=flux(imap(ih),iz2,ig)
                     w=dfdm/d_3d(imap(ih),iz1,ig)/(1.+hzr)
                     phis=w*flx2+(1.-w)*flx1+bt*(flx1+flx2)*0.5
                     cjn=-((dfdm+dhatm)*flx2-(dfdm-dhatm)*flx1)*2. &
                          /(gridsize_z(iz1)+gridsize_z(iz2))
                     cntoavg1=(phis+2.*cjn)*0.25
                     cntoavg2=(phis-2.*cjn)*0.25
                     DO ig2=mgb(ig),mge(ig)
                        tt1=fcntzo(ig2,id1,ih,iz1)
                        tt2=fcntzo(ig2,id2,ih,iz2)
                        cntzo(ig2,id1,ih,iz1)=tt1*cntoavg1
                        cntzo(ig2,id2,ih,iz2)=tt2*cntoavg2
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               iz1=neigsndz(ind,iz)
               id1=neigsndz(ind+3,iz)
               IF(iz.EQ.1) THEN
                  DO ig=1,2 ! n_group
                     alphaz(ig)=alzl(ig)
                  ENDDO
               ELSE
                  DO ig=1,2 ! n_group
                     alphaz(ig)=alzr(ig)
                  ENDDO
               ENDIF
               DO ih=1,nassy
                  DO ig=1,2 ! n_group
                     dfdm=dfdz(ig,ih,iz)
                     dhatm=dhatz(ig,ih,iz)
                     bt=betaphisz(ig,ih,iz)
                     IF(ind.EQ.2) dhatm=-dhatm
                     flx1=flux(imap(ih),iz1,ig)
                     cjn=(dfdm-dhatm)*flx1/gridsize_z(iz1)
                     phis=bt*flx1
                     cntoavg1=(phis+2.*cjn)*0.25
                     DO ig2=mgb(ig),mge(ig)
                        tt1=fcntzo(ig2,id1,ih,iz1)
                        cntzo(ig2,id1,ih,iz1)=tt1*cntoavg1
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
      DO iz=1,nz
         arear=hside*gridsize_z(iz)
         DO ih=1,nassy
            DO it=1,6
               nn=neignd(it,ih)
               IF(nn.EQ.0) THEN
                  DO ig=1,n_group
                     cnti(ig,it)=cnto(ig,it,ih,iz)*reflratf(ig)
                  ENDDO
               ELSE
                  DO ig=1,n_group
                     cnti(ig,it)=cnto(ig,neigjin(it,ih),nn,iz)
                  ENDDO
               ENDIF
            ENDDO
! preparation for axial solver(cnti,srczn,pbflx)
            neigdn=neigz(1,iz)
            neigup=neigz(2,iz)
            IF(neigdn.NE.0) THEN
               DO ig=1,n_group
                  cntz0i(ig)=cntzo(ig,2,ih,neigdn)
               ENDDO
            ELSE
               DO ig=1,n_group
                  cntz0i(ig)=reflratzbf(ig)*cntzo(ig,1,ih,iz)
               ENDDO
            ENDIF
            IF(neigup.NE.0) THEN
               DO ig=1,n_group
                  cntz1i(ig)=cntzo(ig,1,ih,neigup)
               ENDDO
            ELSE
               DO ig=1,n_group
                  cntz1i(ig)=reflratztf(ig)*cntzo(ig,2,ih,iz)
               ENDDO
            ENDIF
            vol=NodeVolume(imap(ih),iz)
            DO ig=1,2 !n_group
               DO ig2=mgb(ig),mge(ig)
                  hflxfn(ig2)=fhflx(ig2,ih,iz)*flux(imap(ih),iz,ig)
               ENDDO
            ENDDO
!
            DO ig=1,n_group
               sumlk=0
               DO it=1,6
                  sumlk=sumlk+cnto(ig,it,ih,iz)-cnti(ig,it)
               ENDDO
               sumlk=sumlk*arear*wtass(ih)
               sumlk=sumlk+(cntzo(ig,1,ih,iz)-cntz0i(ig) &
                    +cntzo(ig,2,ih,iz)-cntz1i(ig))*hexarea*wtass(ih)
               sumss=0
               DO ig2=iscatib(ig),ig-1
                  sumss=sumss+xssf(ig2,ig,ih,iz)*hflxfn(ig2)
               ENDDO
               DO ig2=ig+1,iscatie(ig)
                  sumss=sumss+xssf(ig2-1,ig,ih,iz)*hflxfn(ig2)
               ENDDO
               aphif(ig,ih,iz)=sumlk &
                    +(maXs_r_3D_MG(imap(ih),iz,ig)*hflxfn(ig)-sumss)*vol
            ENDDO
         ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE mphihmg

      subroutine settfsp_hex(first)
      use Inc_BiCG
      use Inc_Solver
      use Inc_Lscoef, only: af
      use Inc_FixSrc_TR
      use Inc_maXS
      use Inc_Time
      use inc_extsrc, only: flag_extsrc, extsrc
!      use inc_tpen, only: imapsol, imap
!      use inc_tpen, only: velo_tr
      use inc_tpen
!      use Inc_transient, only: flag_exp
      implicit none
      integer :: i, k, l, m, ik
      real(8) :: ramprec
      real(8) :: vol
      logical(4) :: first

      integer(4) :: ixy
!      real(8) :: cetakrl, rcdeltl
      real(8) :: phifl

      !! === INDEX ============ !!
      !! af : hex index         !!
      !! ====================== !!


      if (abs(cetak-1d0)>1d-10) then
         call AxBHEX(Flux,flux_add)
      endif
!      do k=1,Nz
      do k=IzFuelBot, IzFuelTop
         do ixy = 1,nxy   !! hex index
            l = imap(ixy) !! square index
            if (I_Comp(l,k)==0) cycle
            if ( I_FARF_1N( I_4Nto1N(l) ) /= 2 ) cycle
            
            ik=kincomp(I_Comp(l,k))                                    
            vol=NodeVolume(l,k)
            rvdtdvol(1)=vol*rvdelt_Old(1,l,k)
            rvdtdvol(2)=vol*rvdelt_Old(2,l,k)
            do m=1,2
               if ( (Flux_Old(l,k,m)>0d0).and.(Flux_OOld(l,k,m)>0d0).and.&
                    (Flux_Old(l,k,m)/Flux_OOld(l,k,m)<20d0) ) then !prevent excessive extrapolation
                  u_omega(m,l,k)=log(abs(Flux_Old(l,k,m)/Flux_OOld(l,k,m)))/dT_Old
               else
                  u_omega(m,l,k)=0d0
               endif
               flux_extra(m,l,k)=dexp(u_omega(m,l,k)*dT)
            enddo

!#ifdef jr_opt
            if ((flag_exp) .and. .not.first) then
!#else
!            if (.not.first) then
!#endif
!#ifdef jr_opt
               do m=1,2
                  if (n_group > 2) then
                     rvdelt(m,l,k)=1d0/velo_tr(l,k,m)/cetadelt+u_omega(m,l,k)*1d0/velo_tr(l,k,m)
                     rvdtvol(m)=vol*1d0/velo_tr(l,k,m)/cetadelt
                  else
                     rvdelt(m,l,k)=v_Inv(l,k,m)/cetadelt+u_omega(m,l,k)*v_Inv(l,k,m)
                     rvdtvol(m)=vol*v_Inv(l,k,m)/cetadelt
                  endif
               enddo
               if (ltolfa(l)==0.or.k>IzFuelTop.or.k<IzFuelBot) then
                  cetakrl=0d0
               else
                  cetakrl=cetakr
               endif
            else
!#endif
               do m=1,2
                  if (ltolfa(l)==0.or.k>IzFuelTop.or.k<IzFuelBot) then
                     cetakrl=0d0
                     rcdeltl=1d0/dT
                  else
                     cetakrl=cetakr
                     rcdeltl=1d0/(dT*cetak)
                  endif
                  if (n_group > 2) then
                     rvdelt(m,l,k)=rcdeltl*1d0/velo_tr(l,k,m)
                     rvdtvol(m)=vol*rvdelt(m,l,k)
                  else
                     rvdelt(m,l,k)=rcdeltl*v_Inv(l,k,m)
                     rvdtvol(m)=vol*rvdelt(m,l,k)
                  endif
               enddo

!#ifdef jr_opt
            endif
!#endif
            
            omegalm=0d0
            omegal0=0d0
            omegalp=0d0
            if (abs(af(2,imapsol(l),k))>1d-10) then
               if (nordprec==1) then
                  do i=1,N_Group_d
                     omegam(i,l,k)=0d0
                     omega0(i,l,k)=beta_d_eff(l,k,i)*(capbrldt(l,k,i)-cappa(l,k,i))
                     omegap(i,l,k)=beta_d_eff(l,k,i)*(1-capbrldt(l,k,i))          
                     omegalm=0d0
                     omegal0=omegal0+omega0(i,l,k)
                     omegalp=omegalp+omegap(i,l,k)
                  enddo
               else
                  do i=1,N_Group_d
                     omegam(i,l,k)=beta_d_eff(l,k,i)/lambda_d(l,k,i) &
                        & *rldtgp1(l,k,i)*(2*capbrldt(l,k,i)-gamma*cappap1(l,k,i))
                     omega0(i,l,k)=beta_d_eff(l,k,i)/lambda_d(l,k,i) &
                        & *(rldt(l,k,i)*(cappap1(l,k,i)+capbrldt2(l,k,i)*rgamma)-cappa(l,k,i))
                     omegap(i,l,k)=beta_d_eff(l,k,i)/lambda_d(l,k,i) &
                        & *(1d0-rldtgp1(l,k,i)*(2d0+capbrldt2(l,k,i)*rgamma))
                     omegalm=omegalm+omegam(i,l,k)*lambda_d(l,k,i)
                     omegal0=omegal0+omega0(i,l,k)*lambda_d(l,k,i)
                     omegalp=omegalp+omegap(i,l,k)*lambda_d(l,k,i)
                  enddo
               endif
            endif
      
            ! establish effective source
            ! calculate g*phi term
            omegalpd=betap_Old(l,k)+beta_d_Tot(l,k)-1d0
            sd=omegalpd*FisSrc_Old(l,k)
            flux_add(l,k,1)=flux_add(l,k,1)-rvdtdvol(1) &
               & *Flux_Old(l,k,1)+maxs_chid_3d(l,k,1)*sd
            flux_add(l,k,2)=flux_add(l,k,2)-rvdtdvol(2) &
               & *Flux_Old(l,k,2)+maxs_chid_3d(l,k,2)*sd
            sd=0d0
            sdn=0d0
            do i=1,N_Group_d
               ramprec=lambda_d(l,k,i)*C_d_I(l,k,i)
               sd=sd+ramprec
               sdn=sdn+cappa(l,k,i)*ramprec
            enddo
            
!#ifdef jr_opt
            if (.not.flag_exp) then
                sdnt = sdn+omegalm*FisSrc_OOld(l,k)+omegal0*FisSrc_Old(l,k) 
                do m = 1,2
                   src(m,l,k) = rvdtvol(m)*flux_old(l,k,m) &
                      +cetakrl*(maXS_chid_3D(l,k,m)*sd-flux_add(l,k,m)) &
                      +maXS_chid_3D(l,k,m)*sdnt         
                enddo
                if (n_group > 2) then
                   spnt=(1-beta_d_Tot(l,k))*FisSrc_Old(l,k)
                   do m = 1,n_group
                      phifl=flux_old(l,k,igc(m))*fluxf(l,k,m)*nodevolume(l,k)
                      rvdeltf(m,l,k) = v_inv(l,k,m)*rcdeltl 
                      srcf(m,l,k)    = rvdeltf(m,l,k)*phifl &
                          -aphif(m,ixy,k)*cetakrl &
                          +maXs_chi_3D_MG(l,k,m)*spnt*cetakrl &
                          +maXs_chid_3D_MG(l,k,m)*(sdnt+sd*cetakrl) 
                   enddo
                endif
            else
!#endif
                src_tr_fac1(1,l,k)=Flux_Old(l,k,1)*rvdtvol(1) &
                   & +cetakrl*(sd-flux_add(l,k,1))
                src_tr_fac1(2,l,k)=Flux_Old(l,k,2)*rvdtvol(2) &
                   & -cetakrl*flux_add(l,k,2)
                src_tr_fac2(1,l,k)=sdn+omegalm*FisSrc_OOld(l,k) &
                   & +omegal0*FisSrc_Old(l,k)
                src_tr_fac2(2,l,k)=0d0
                if (n_group > 2) then
                   src(1,l,k)=flux_extra(1,l,k)*(src_tr_fac1(1,l,k)-cetakrl*u_omega(1,l,k) &
                      & *1d0/velo_tr(l,k,1)*vol*Flux_Old(l,k,1))+src_tr_fac2(1,l,k)
                   src(2,l,k)=flux_extra(2,l,k)*(src_tr_fac1(2,l,k)-cetakrl*u_omega(2,l,k) &
                      & *1d0/velo_tr(l,k,2)*vol*Flux_Old(l,k,2))
                else
                   src(1,l,k)=flux_extra(1,l,k)*(src_tr_fac1(1,l,k)-cetakrl*u_omega(1,l,k) &
                      & *v_Inv(l,k,1)*vol*Flux_Old(l,k,1))+src_tr_fac2(1,l,k)
                   src(2,l,k)=flux_extra(2,l,k)*(src_tr_fac1(2,l,k)-cetakrl*u_omega(2,l,k) &
                      & *v_Inv(l,k,2)*vol*Flux_Old(l,k,2))
                endif
!#ifdef jr_opt
            endif
!#endif
            if (flag_extsrc) then
               src(1,l,k)=src(1,l,k)+extsrc(l,k)*(1d0+cetakrl*flux_extra(1,l,k))
            endif
            betap(l,k)=1d0-beta_d_Tot(l,k)+omegalp
         enddo
      enddo
      ifnodal=.true.

      return
      end subroutine settfsp_hex

#endif
#ifdef tuan_tr_test

      subroutine inittran_hex
      use Inc_3D
      use Inc_TH, only: ppower_0
      use Inc_CR
      use Inc_Kinetics
      use inc_tpen
      use inc_xs_file, only: Ref_maXS

      implicit none
      real(8) :: tmp1, vol
      real(8) :: psi0, sumd
      integer :: k, l, m, ik, icomp, Ixy, Ixy_FA, Iz, Ig_d
      integer :: nxyfb, nzfb
      integer(4) :: ig
      real(8) :: rvdeltl
      integer(4) :: mf

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [inittran] in Mod_TRnodal'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      rhoadj = D0
      call Alloc_FBRHO
      if (n_group > 2) then
         if (.not.allocated(velo_tr)) allocate(velo_tr(nxy,nz,2))
         velo_tr = 0d0
      endif

      ! assign transient iteration control parameters
      Period_TH    = 5          ! Number of T/H updates per nodal update 5
      Period_CMFD  = 1          ! Number of iterations to be performed before the first nonlinear update 1
      Period_NonL  = 5          ! Nonlinear update cycle 5
    !! for comparison with PARCS
      EPS_ERF = 1.D-4 

      ! assign nodewise kinetics parameters
!      if (Flag_Card_maXS) then
         beta_d_Tot=0d0
         DO k = IzFuelBot, IzFuelTop
            DO Ixy_FA = 1, Nxy_FA
            l = I_FA(Ixy_FA)
               icomp=I_Comp(l,k)
               ik=kincomp(icomp)
               do m=1,N_Group
                  v_Inv(l,k,m)=1d0/tvelo(m,ik)
               enddo

               if (n_group > 2) then
                  do m = 1,2
                     rvdeltl = 0 
                     do mf = mgb(m),mge(m)
                        rvdeltl=rvdeltl+v_Inv(l,k,mf)*fluxf(l,k,mf)
                     enddo
                     velo_tr(l,k,m) = 1/rvdeltl
                  enddo
               endif

               if (nu_maXs_f_3D(l,k,2) .NE. 0) then
                  do Ig_d=1,N_Group_d
                     beta_d_Tot(l,k)=beta_d_Tot(l,k)+tbeta(Ig_d,ik)
                     beta_d_eff(l,k,Ig_d)=tbeta(Ig_d,ik)
                     lambda_d(l,k,Ig_d)=lambda_d_I(Ig_d,ik)
                  enddo
               endif
            enddo
         enddo
!      else ! .not. flag_card_maxs
!         call alloc_tran
!         do Ixy_FA=1,Nxy_FA
!            Ixy=I_FA(Ixy_FA)
!            do Iz=IzFuelBot,IzFuelTop
!               do Ig_d=1,N_Group_d
!                  beta_d_Tot(Ixy,Iz)=beta_d_Tot(Ixy,Iz)+beta_d_eff(Ixy,Iz,Ig_d)
!               enddo
!            enddo
!         enddo
!      endif

      if (flag_card_maxs) then
         if (flag_force_beff) then
            beta_d_tot=0d0
            do m=1,n_group_d
               do k=1,nz
                  do l=1,nxy
                     beta_d_eff(l,k,m)=beta_d_eff(l,k,m)*force_beff(m)
                     beta_d_tot(l,k)=beta_d_tot(l,k)+beta_d_eff(l,k,m)
                  enddo
               enddo
            enddo
         endif
         if (flag_force_lamb) then
            do m=1,n_group_d
               do k=1,nz
                  do l=1,nxy
                     lambda_d(l,k,m)=lambda_d(l,k,m)*force_lamb(m)
                  enddo
               enddo
            enddo
         endif
         if (flag_force_velo) then
            do m=1,n_group
               do k=1,nz
                  do l=1,nxy
                     v_inv(l,k,m)=v_inv(l,k,m)*force_velo(m)
                  enddo
               enddo
            enddo
         endif
      endif

      ! initialize frequently used constants
      reigv=1d0/keff
      cetak0=cetak
      cetak=1d0
      cetakb=1d0-cetak
      cetakr=cetakb/cetak
      cetafb=1d0-cetaf
      cetacb=1d0-cetac
      cetadelt=dT
      if (Flag_THFB) then
         cetacr = cetacb / cetac
         kgap   = cetaf  * h_Gap * delr
         kgapb  = cetafb * h_Gap * delr
         tmp1   = h_Gap  * h_Clad * (4d0-h_Clad/R_Gap)*R_Fuel / R_Gap
         kgap2  = cetaf  * h_Gap * h_Clad * R_Fuel / R_Gap
         kgap4  = cetaf  * tmp1
         kgap4b = cetafb * tmp1
      endif

      ! calculate steady state precursor concentration
      do k = 1,nz
         do l = 1,nxy
            flux_oold(l,k,1) = flux(l,k,1)
            flux_oold(l,k,2) = flux(l,k,2)
            flux_old(l,k,1) = flux(l,k,1)
            flux_old(l,k,2) = flux(l,k,2)
!            write(1002,'(2I5,2E13.6)') l,k, flux(l,k,1),flux(l,k,2)

         enddo
      enddo
      af           = reigv * af
      nu_maXs_f_3D = reigv * nu_maXs_f_3D
      !! make nu to adjust core critical
!      if (Flag_Card_maXS) then
         do ig = 1, n_group
!            Ref_maXS(:, (3*N_Group + ig))=Ref_maXS(:, (3*N_Group + ig))*reigv
            XSset_Hex( :, :,4, Ig )=XSset_Hex( :, :,4, Ig )*reigv
         enddo
         nufactor = reigv
!      endif
      

      call Get_Source(FisSrc,nu_maXS_f_3D,Flux)
!      call get_fissrc 

      FisSrc_Old  = FisSrc
      FisSrc_OOld = FisSrc

      do Ixy_FA=1,Nxy_FA
         Ixy=I_FA(Ixy_FA)
         do Iz=IzFuelBot,IzFuelTop
            do Ig_d=1,N_Group_d
               if (abs(lambda_d(Ixy,Iz,Ig_d))>1d-10) then
                  C_d_I(Ixy,Iz,Ig_d)=beta_d_eff(Ixy,Iz,Ig_d) &
                      &  /lambda_d(Ixy,Iz,Ig_d)*FisSrc_Old(Ixy,Iz)
               endif
            enddo
         enddo
      enddo

      psi0=ppower_0
      sumv0=0d0
      sumd=0d0
!      do k=1,Nz
!         do l=1,Nxy
     DO k = IzFuelBot, IzFuelTop
         DO Ixy_FA = 1, Nxy_FA
            l = I_FA(Ixy_FA)
            vol=NodeVolume(l,k)
            sumd=sumd+FisSrc(l,k)*flux_adj(l,k,1)
            do m=1,2
               sumv0=sumv0+Flux(l,k,m)*1/velo_tr(l,k,m)*flux_adj(l,k,m)*vol
            enddo
         enddo
      enddo
      sumd=1d0/sumd
      sumv0=sumv0/psi0

      ! save initial feedback variables (time=0.0)
      nxyfb = nchan
      nzfb  = nzth
      t_Fuel0  = T_Fuel
      t_Mod0   = T_Mod
      d_Mod0   = D_Mod
#ifdef tuan_tr_test
if(flag_thfb) then
      !! === For calculate reactivity === !!
      if (.not.allocated(tdopl0)) allocate(tdopl0(0:nz,0:nxy_fa))
      if (.not.allocated(tcool0)) allocate(tcool0(0:nz,0:nxy_fa))
      if (.not.allocated(dcool0)) allocate(dcool0(0:nz,0:nxy_fa))
      do Iz = 1,Nz
         do Ixy_Fa = 1,nxy_fa
            tdopl0(iz,ixy_fa) = tdopl(iz,ixy_fa)
            tcool0(iz,ixy_fa) = tcool(iz,ixy_fa)
            dcool0(iz,ixy_fa) = dcool(iz,ixy_fa)
         enddo
      enddo
endif
#endif      
      
      if (n_cr>0)  CR_Bot0 = Abs_Bot
      if (n_cr>0)  CR_Top0 = TopFol_Bot
      if (n_cr>0)  BotFol_Bot0 = BotFol_Bot
!      do k=1,Nz
!         do l=1,Nxy+Nx+Ny
!            do m=1,N_Group
!               dnw0(m,l,k)=dnw(m,l,k)
!               dnn0(m,l,k)=dnn(m,l,k)
!            enddo
!         enddo
!      enddo
!      do k=1,Nz+1
!         do l=1,Nxy
!            do m=1,N_Group
!               dnb0(m,l,k)=dnb(m,l,k)
!            enddo
!         enddo
!      enddo
      !! feedback rho, hexagonal geometry
      !! this parameter has relationship with reactivity calculation.
      dhat0 =dhat
      dhatz0=dhatz
      keff=1d0
      reigv=1d0
      reigvs=1d0
      iptype=1 !! for TPEN solver (has relationship with TPENBC subroutine)

      return
      end subroutine inittran_hex


      subroutine settran_hex
      implicit none


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [settran] in Mod_TRnodal'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      cetak    = 0.5d0 ! central differencing as the default
      cetaf    = 0.5d0 ! theta for heat conduction in fuel
      cetac    = 1.0d0 ! theta for heat convection in coolan
      nordprec = 2d0   ! second order precursor integration

      return
      end subroutine settran_hex
#endif

      END MODULE Mod_TRnodal

#endif
