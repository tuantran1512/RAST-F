!#ifdef siarhei_delete    
#ifdef tuan_tr_test    

      MODULE Mod_CalReact

      use Inc_Constant
      use Inc_Geometry
      use Inc_maXS
      use Inc_3D
      use Inc_Solver
      use Inc_Kinetics
      use Inc_Control
      use Inc_Lscoef
      use Inc_Adjoint
      use Inc_Transient
      use Inc_miXS
      use Inc_Nuclide
      use Inc_File, only: W_PKD
      use Mod_SolLS, only: AxBHEX
      use Inc_FluxVar
      use Inc_RP
      use Inc_Option
      use Inc_TH, only: ppower
      use Inc_Time, only: sec
      use Inc_Utile, only: div_line2, div_line
#ifdef tuan_tr_test
      use Inc_Flag, only:  flag_thfb
#endif

      
#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

      CONTAINS

      subroutine calrho(first)
      use Mod_GetSome
      implicit none
      logical(4) :: first
      integer :: k, l, m, i
      real(8) :: sumn, sumd, betatavg, betatadj, psil2, sumv, vol
      real(8) :: rvdtvol, rhoadjd, pcfact, temp
      integer :: next=1
      integer :: KX,MX,NNX,KY,MY,NNY,KZ,MZ,NNZ
      real(8) :: BRANGE(2),ERANGE(2),HRANGE(2)
      integer :: IX, IY, IZ, N, NN, LL
      real(8) :: T(50), tmp_beta(6)=0d0, tmp_lamb(6)=0d0
      real(8) :: tmp_BU, tmp_ENR, tmp_H2U, tmp_VOL
      integer :: IFTY

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [calrho] in Mod_CalReact'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [calrho] in Mod_CalReact'
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

      call AxBHEX(flux,flux_add)

      sumn=0d0
      sumd=0d0
      sumv=0d0
      betatavg=0d0
      betatadj=0d0
      psil2=0d0
      do i=1,N_Group_d
         Pzeta(i,next)=0d0
         Plambda(i,next)=0d0
         Pbeta(i)=0d0
      enddo
      do k=1,nz
         do l=1,nxy
            vol=NodeVolume(l,k)
            do m=1,N_group
               rvdtvol=vol*rvdelt(m,l,k)
               flux_add(l,k,m)=flux_add(l,k,m)-rvdtvol*flux(l,k,m) &
                  & +maxs_chi_3d(l,k,m)*(betap(l,k)-1d0)*FisSrc(l,k)
               sumn=sumn+flux_add(l,k,m)*flux_adj(l,k,m)
               sumd=sumd+maxs_chi_3d(l,k,m)*FisSrc(l,k)*flux_adj(l,k,m)
               sumv=sumv+flux(l,k,m)*v_inv(l,k,m)*flux_adj(l,k,m)*vol
               betatadj=betatadj+beta_d_tot(l,k)*FisSrc(l,k) &
                  & *flux_adj(l,k,m)*maxs_chid_3d(l,k,m)
               psil2=psil2+FisSrc(l,k)*flux_adj(l,k,m)*maxs_chi_3d(l,k,m)
            enddo
            betatavg=betatavg+beta_d_tot(l,k)*FisSrc(l,k)
            do i=1,N_Group_d
               temp=flux_adj(l,k,1)*C_d_I(l,k,i)
               Pzeta(i,next)=Pzeta(i,next)+temp
               Plambda(i,next)=Plambda(i,next)+temp*lambda_d(l,k,i)
               Pbeta(i)=Pbeta(i)+flux_adj(l,k,1)*FisSrc(l,k)*beta_d_eff(l,k,i)
            enddo
         enddo
      enddo
      betatavg=betatadj/psil2

      rhoadjd = rhoadj
      sumd=1d0/sumd
      rhoadj=-sumn*sumd/betatavg
      iPgenT(next)=sumv*sumd

      pcfact=ppower/sumv*sumv0
      do i=1,N_Group_d
         Pbeta(i)=Pbeta(i)*sumd
         Plambda(i,next)=Plambda(i,next)/Pzeta(i,next)
      enddo

      if (flag_print_kp_anc) then
         KX=5; MX=0; NNX=4
         KY=5; MY=0; NNY=4
         KZ=2; MZ=0; NNZ=1
         BRANGE(1) = 1.0D-3; BRANGE(2) = 6.0D+3
         ERANGE(1) = 0.5D+0; ERANGE(2) = 6.0D+0
         HRANGE(1) = 4.0D+0; HRANGE(2) = 5.0D+0

         tmp_BU = 0d0; tmp_VOL = 0d0
         do l=1,Nxy
            do k=IzFuelBot,IzFuelTop
               tmp_BU=tmp_BU+BU(l,k)*h_x(l)*h_y(l)*h_z(k)
               tmp_VOL=tmp_VOL+h_x(l)*h_y(l)*h_z(k)
            enddo
         enddo
         tmp_BU=tmp_BU/tmp_VOL

         tmp_ENR=4.5D0

         IFTY=int((tmp_BU+1d-1),4)
         IFTY=MAX(IFTY,1)
         tmp_H2U=1D0*COEF_TYPEF(IFTY)/tmp_ENR

         tmp_BU =MAX(BRANGE(1),MIN(tmp_BU ,BRANGE(2)))
         tmp_ENR=MAX(ERANGE(1),MIN(tmp_ENR,ERANGE(2)))
         tmp_H2U=MAX(HRANGE(1),MIN(tmp_H2U,HRANGE(2)))

         T=0d0
         N=1
         do IZ=MZ,NNZ
            do IY=MY,NNY
               do IX=MX,NNX
                  T(N)=(tmp_BU**IX)*(tmp_ENR**IY)*(tmp_H2U**IZ)
                  N=N+1
               enddo
            enddo
         enddo
         NN=KX*KY*KZ
         do LL=1,NN
            tmp_beta(1)=tmp_beta(1)+T(LL)*COEF_FB1(LL)
            tmp_beta(2)=tmp_beta(2)+T(LL)*COEF_FB2(LL)
            tmp_beta(3)=tmp_beta(3)+T(LL)*COEF_FB3(LL)
            tmp_beta(4)=tmp_beta(4)+T(LL)*COEF_FB4(LL)
            tmp_beta(5)=tmp_beta(5)+T(LL)*COEF_FB5(LL)
            tmp_beta(6)=tmp_beta(6)+T(LL)*COEF_FB6(LL)
            tmp_lamb(1)=tmp_lamb(1)+T(LL)*COEF_FL1(LL)
            tmp_lamb(2)=tmp_lamb(2)+T(LL)*COEF_FL2(LL)
            tmp_lamb(3)=tmp_lamb(3)+T(LL)*COEF_FL3(LL)
            tmp_lamb(4)=tmp_lamb(4)+T(LL)*COEF_FL4(LL)
            tmp_lamb(5)=tmp_lamb(5)+T(LL)*COEF_FL5(LL)
            tmp_lamb(6)=tmp_lamb(6)+T(LL)*COEF_FL6(LL)
         enddo
      endif

      if (first) then
         if (flag_print_kp_anc) then
            write(W_PKD,*) 'ANC Kinetics Parameter ******************'
            write(W_PKD,'(A10,F6.3,A6,F6.3,A6,F6.3)') ' Input BU=',tmp_BU, ', ENR=',tmp_ENR, ', H2U=',tmp_H2U
            write(W_PKD,*) 'Group      beta      lambda'
            do LL=1,6
               write(W_PKD,'(I4, 2F12.6)') LL, tmp_beta(LL), tmp_lamb(LL)
            enddo
            write(W_PKD,*) 'ANC Kinetics Parameter ******************'
            write(W_PKD,*) ''
         endif

         write(W_PKD,*) Div_Line2(1:70)
         write(W_PKD,*) "Point Kinetics Data"
         write(W_PKD,*) ""
         write(W_PKD,*) "[ Reactivity Coefficients ]"
         write(W_PKD,*) "[ Kinetics Parameters ]"
         write(W_PKD,*) "[ Time Dependent Reactivity ]"
         write(W_PKD,*) ""

         call fdbkrho(betatavg,sumd,rhoadj,.true.)

         write(W_PKD,*) Div_Line2(1:70)
         write(W_PKD,*) "[ Kinetics Parameters ]"
         write(W_PKD,*) ""
         write(W_PKD,*) "  Group    Beta      Lambda[1/sec]"
         do i=1,N_Group_d
            write(W_PKD,'(I6,f12.6,f14.6)') i, Pbeta(i), Plambda(i,next)
         end do
         write(W_PKD,'(A,f10.6)') "   Total", sum(Pbeta)
         write(W_PKD,*) ""

         write(W_PKD,*) Div_Line2(1:70)
         write(W_PKD,*) "[ Time Dependent Reactivity ]"
         write(W_PKD,*) ""
         !write(W_PKD,*) "    Time[sec]       rho[$]     Gen. Time[sec]     plevel"
         write(W_PKD,'(4a13,a20,a80)') "Time [s]","rho [$]","plevel","Gen. T [s]", &
                      & "beta(1:n_group_d)","lambda(1:n_group_d)"
      endif
      !WRITE(W_PKD,'(1p,20e15.6)') Sec, rhoadj, iPgenT(next), ppower
      WRITE(W_PKD,'(1p,16e13.5)') Sec, rhoadj, ppower/pcfact, iPgenT(next), &
                                & Pbeta(:), Plambda(:,next)

      call fdbkrho(betatavg,sumd,rhoadj,.false.)

      return
      end subroutine calrho


      subroutine fdbkrho(betatavg,sumd,rhoadjn,first)
      use Inc_TH
      use Inc_File
      use Inc_CR
      use Inc_INP, only: Flag_Card_maXS
      use Mod_Alloc
      implicit none
      logical(4) :: first
      real(8) :: toutavgt, tfmaxt
      real(8) :: rhon, rhoadjn, rhod, betatavg, sumd
      integer :: l, k, m
      real(8) :: tfrho
      real(8) :: dmrho
      real(8) :: ppmrho
      real(8) :: crrho
      real(8) :: sumrho
      real(8) :: xesmrho
      real(8) :: leakrho
      real(8) :: nullrho
      real(8) :: tmp_delta
      logical(1), save :: flag_alloc_first=.true.

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [fdbkrho] in Mod_CalReact'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fdbkrho] in Mod_CalReact'
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

      if (flag_alloc_first) then
         flag_alloc_first=.false.
         WRITE(W_RHO,'(1x,13a15)') 'time', 'sumrho', 'tfrho', 'dmrho', &
                                 & 'ppmrho','crrho','xesmrho', 'leakrho', 'nullrho'
      endif

      if (.not.allocated(bk_CR_Bot)) call alloc (bk_CR_Bot ,n_cr)
      if (.not.allocated(bk_CR_Top)) call alloc (bk_CR_Top ,n_cr)
      if (.not.allocated(bk_tFuel )) call alloc (bk_tFuel  ,nxy,nz)
      if (.not.allocated(bk_tMod  )) call alloc (bk_tMod   ,nxy,nz)
      if (.not.allocated(bk_dMod  )) call alloc (bk_dMod   ,nxy,nz)
      if (.not.allocated(bk_amd   )) call alloc (bk_amd    ,N_group,nxy,nz)
      if (.not.allocated(bk_afd   )) call alloc (bk_afd    ,N_group,nxy,nz)
      if (.not.allocated(bk_af2d  )) call alloc (bk_af2d   ,nxy,nz)
      if (.not.allocated(bk_scatd )) call alloc (bk_scatd  ,nxy,nz)
      if (.not.allocated(bk_ccwd  )) call alloc (bk_ccwd   ,N_group,nxy,nz)
      if (.not.allocated(bk_cced  )) call alloc (bk_cced   ,N_group,nxy,nz)
      if (.not.allocated(bk_ccnd  )) call alloc (bk_ccnd   ,N_group,nxy,nz)
      if (.not.allocated(bk_ccsd  )) call alloc (bk_ccsd   ,N_group,nxy,nz)
      if (.not.allocated(bk_ccbd  )) call alloc (bk_ccbd   ,N_group,nxy,nz)
      if (.not.allocated(bk_cctd  )) call alloc (bk_cctd   ,N_group,nxy,nz)
      if (.not.allocated(bk_dnwd  )) call alloc (bk_dnwd   ,N_group,nxy+nx+ny,nz)
      if (.not.allocated(bk_dnnd  )) call alloc (bk_dnnd   ,N_group,nxy+nx+ny,nz)
      if (.not.allocated(bk_dnbd  )) call alloc (bk_dnbd   ,N_group,nxy,nz+1)
      if (.not.allocated(bk_xsxead)) call alloc (bk_xsxead ,N_group,nxy,nz)
      if (.not.allocated(bk_xssmad)) call alloc (bk_xssmad ,N_group,nxy,nz)
      if (.not.allocated(bk_rnxed )) call alloc (bk_rnxed  ,nxy,nz)
      if (.not.allocated(bk_rnsmd )) call alloc (bk_rnsmd  ,nxy,nz)

      toutavgt = toutavg
      tfmaxt   = tfmax

      bk_tFuel = T_Fuel
      bk_tMod  = T_Mod
      if (Flag_Card_maXS) then
         bk_dMod=D_Mod
      endif
      bk_ppm    = PPM
      if (n_cr>0) bk_CR_Bot = CR_Bot
      if (n_cr>0) bk_CR_Top = CR_Top

      ! Xe/Sm number density <- need to be added
      do k=1,nz
         do l=1,nxy
            do m=1, N_group
               bk_xsxead(m,l,k)=miXS_a_Xe35(l,k,m)
               bk_xssmad(m,l,k)=miXS_a_Sm49(l,k,m)
            enddo
            bk_rnxed(l,k)=N_Xe35(l,k)
            bk_rnsmd(l,k)=N_Sm49(l,k)
         enddo
      enddo

      if (first) then
         rhocoef = 0

         ! Compute boron reactivity coeficient
         tmp_delta = 10d0
         PPM = PPM + tmp_delta
         rhocoef(1) = (fbrho(betatavg,sumd)-rhoadjn)/tmp_delta
         PPM = bk_ppm

         ! Compute coolant density coeficient
         if (Flag_Card_maXS) then
             tmp_delta = 1.d-3
             d_mod = d_mod + tmp_delta
             rhocoef(2) = fbrho(betatavg,sumd)-rhoadjn
             d_mod = bk_dMod
         end if

         ! Compute coolant temperature coeficient
         t_mod = t_mod + tmp_delta
         rhocoef(3) = (fbrho(betatavg,sumd)-rhoadjn)/tmp_delta
         t_mod = bk_tMod

         ! Compute doppler temperature coeficient
         tmp_delta = 10.d0
         t_fuel = t_fuel + tmp_delta
         rhocoef(4) = (fbrho(betatavg,sumd)-rhoadjn)/tmp_delta
         T_Fuel = bk_tFuel

         ! Compute xe/sm component of reactivity
         do k=1,nz
            do l=1,nxy
               N_Xe35(l,k)=N_Xe35(l,k)+1d+10
            enddo
         enddo
         rhocoef(5)=(fbrho(betatavg,sumd)-rhoadjn)*1d-10
         do k=1,nz
            do l=1,nxy
               N_Xe35(l,k)=bk_rnxed(l,k)
            enddo
         enddo

         do k=1,nz
            do l=1,nxy
               N_Sm49(l,k)=N_Sm49(l,k)+1d+12
            enddo
         enddo
         rhocoef(6)=(fbrho(betatavg,sumd)-rhoadjn)*1d-12
         do k=1,nz
            do l=1,nxy
               N_Sm49(l,k)=bk_rnsmd(l,k)
            enddo
         enddo

         write(W_PKD,*) Div_Line2(1:70)
         write(W_PKD,*) "[ Reactivity Coefficients ]"
         write(W_PKD,'(A,Es12.5,A)') "  - Boron reactivity worth : ", rhocoef(1), " [$/ppm]"
         if (Flag_Card_maXS) then
            write(W_PKD,'(A,Es12.5,A)') "  - Coolant density worth  : ", rhocoef(2), " [$/kgm^-3]"
         endif
         write(W_PKD,'(A,Es12.5,A)') "  - Coolant temp. worth    : ", rhocoef(3), " [$/K]"
         write(W_PKD,'(A,Es12.5,A)') "  - doppler temp. worth    : ", rhocoef(4), " [$/K]"
         write(W_PKD,'(A,Es12.5,A)') "  - Xenon worth            : ", rhocoef(5), " [$/cm^-3]"
         write(W_PKD,'(A,Es12.5,A)') "  - Samarium worth         : ", rhocoef(6), " [$/cm^-3]"
         write(W_PKD,*) ""
         return
      endif

      do k=1,nz
         do l=1,nxy
            do m=1,N_group
               bk_amd (m,l,k) = am(m,l,k)
               bk_afd (m,l,k) = af(m,l,k)
               bk_ccwd(m,l,k) = ccw(m,l,k)
               bk_cced(m,l,k) = cce(m,l,k)
               bk_ccnd(m,l,k) = ccn(m,l,k)
               bk_ccsd(m,l,k) = ccs(m,l,k)
               bk_ccbd(m,l,k) = ccb(m,l,k)
               bk_cctd(m,l,k) = cct(m,l,k)
            enddo
            bk_af2d (l,k) = af2(l,k)
            bk_scatd(l,k) = scat(l,k)
         enddo
      enddo
      do k=1,nz
         do l=1,nxy+nx+ny
            do m=1,N_group
               bk_dnwd(m,l,k)=dnw(m,l,k)
               bk_dnnd(m,l,k)=dnn(m,l,k)
            enddo
         enddo
      enddo
      do k=1,nz+1
         do l=1,nxy
            do m=1,N_group
               bk_dnbd(m,l,k)=dnb(m,l,k)
            enddo
         enddo
      enddo

      rhon = rhoadjn

      ! compute "nodal leakage" component of reactivity.
      do k=1,nz
         do l=1,nxy+nx+ny
            do m=1,N_group
               dnw(m,l,k)=dnw0(m,l,k)
               dnn(m,l,k)=dnn0(m,l,k)
            enddo
         enddo
      enddo
      do k=1,nz+1
         do l=1,nxy
            do m=1,N_group
               dnb(m,l,k)=dnb0(m,l,k)
            enddo
         enddo
      enddo
      rhod = rhon
      rhon = fbrho(betatavg,sumd)
      leakrho = rhod-rhon

      ! Compute control component of reactivity
      if (n_cr>0) CR_Bot = init_CR_Bot
      if (n_cr>0) CR_Top = init_CR_Top
      rhod   = rhon
      rhon   = fbrho(betatavg,sumd)
      crrho  = rhod - rhon

      ! Compute boron component of reactivity
      ppm    = init_ppm
      rhod   = rhon
      rhon   = fbrho(betatavg,sumd)
      ppmrho = rhod - rhon

      !! Compute moderator density/temperature component of reactivity
      t_mod = init_tMod
      if (Flag_Card_maXS) then
         d_mod = init_dMod
      endif
      rhod  = rhon
      rhon  = fbrho(betatavg,sumd)
      dmrho = rhod - rhon

      ! compute doppler temperature component of reactivity
      t_fuel = init_tFuel
      rhod  = rhon
      rhon  = fbrho(betatavg,sumd)
      tfrho = rhod - rhon

      ! compute xe/sm component of reactivity
      do k=1,nz
         do l=1,nxy
            do m=1,N_group
               miXS_a_Xe35(l,k,m)=init_miXS_a_Xe35(l,k,m)
               miXS_a_Sm49(l,k,m)=init_miXS_a_Sm49(l,k,m)
            enddo
            N_Xe35(l,k)=init_N_Xe35(l,k)
            N_Sm49(l,k)=init_N_Sm49(l,k)
         enddo
      enddo
      rhod    = rhon
      rhon    = fbrho(betatavg,sumd)
      xesmrho = rhod-rhon
      nullrho = rhon

      sumrho = tfrho + dmrho + ppmrho + crrho + xesmrho + leakrho + nullrho

      ! print reactivity components
      WRITE(W_RHO,601) sec, sumrho, tfrho, dmrho, ppmrho, &
                  & crrho, xesmrho, leakrho, nullrho
601   FORMAT(1x,f15.8,9(2x,1p,e13.6))

      ! reset matrix and feedback variables to current time step values
      T_Fuel = bk_tFuel
      T_Mod  = bk_tMod
      if (Flag_Card_maXS) then
         D_Mod = bk_dMod
      endif
      PPM = bk_ppm
      if (n_cr>0) CR_Bot = bk_CR_Bot
      if (n_cr>0) CR_Top = bk_CR_Top

      do k=1,nz
         do l=1,nxy
            do m=1,N_group
               miXS_a_Xe35(l,k,m)=bk_xsxead(m,l,k)
               miXS_a_Sm49(l,k,m)=bk_xssmad(m,l,k)
            enddo
            N_Xe35(l,k)=bk_rnxed(l,k)
            N_Sm49(l,k)=bk_rnsmd(l,k)
         enddo
      enddo

      do k=1,nz
         do l=1,nxy
            do m=1,N_group
               am (m,l,k) = bk_amd (m,l,k)
               af (m,l,k) = bk_afd (m,l,k)
               ccw(m,l,k) = bk_ccwd(m,l,k)
               cce(m,l,k) = bk_cced(m,l,k)
               ccn(m,l,k) = bk_ccnd(m,l,k)
               ccs(m,l,k) = bk_ccsd(m,l,k)
               ccb(m,l,k) = bk_ccbd(m,l,k)
               cct(m,l,k) = bk_cctd(m,l,k)
            enddo
            af2 (l,k) = bk_af2d(l,k)
            scat(l,k) = bk_scatd(l,k)
         enddo
      enddo
      do k=1,nz
         do l=1,nxy+nx+ny
            do m=1,N_group
               dnw(m,l,k)=bk_dnwd(m,l,k)
               dnn(m,l,k)=bk_dnnd(m,l,k)
            enddo
         enddo
      enddo
      do k=1,nz+1
         do l=1,nxy
            do m=1,N_group
               dnb(m,l,k)=bk_dnbd(m,l,k)
            enddo
         enddo
      enddo

      toutavg=toutavgt
      tfmax=tfmaxt

      return
      end subroutine fdbkrho


      function fbrho(betatavg, sumd) result(fbrho_return)
      use Mod_SolLS, only: setls
      use Mod_XSFB, only: XSFBHEX
      use Mod_SolFDM, only: FdmCoupl
      use Inc_maXS, only: nu_maXS_f_3D
!      use mod_interface, only: adjust_force_th, return_force_th
      implicit none
      integer :: k, l
      real(8) :: sumn, sumd, betatavg, rvdtvol(N_group), vol, psixs
      real(8) :: fbrho_return

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [fbrho] in Mod_CalReact'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fbrho] in Mod_CalReact'
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
#ifdef tuan_tr_test

#else
      call adjust_force_th
#endif      
      call XSFBHEX
#ifdef tuan_tr_test

#else
      call return_force_th
#endif      

      call FdmCoupl
      call setls
#ifdef tuan_tr_test

#else    
      call AxB(flux,flux_add)
#endif
      sumn=0
      do k=1,nz
         do l=1,nxy
            vol=nodevolume(l,k)
            rvdtvol(1)=vol*rvdelt(1,l,k)
            rvdtvol(2)=vol*rvdelt(2,l,k)
            psixs=(nu_maXS_f_3D(l,k,1)*flux(l,k,1)+nu_maXS_f_3D(l,k,2)*flux(l,k,2))*vol
            flux_add(l,k,1)=flux_add(l,k,1)-rvdtvol(1)*flux(l,k,1)+(betap(l,k)-1d0)*psixs
            flux_add(l,k,2)=flux_add(l,k,2)-rvdtvol(2)*flux(l,k,2)
            sumn=sumn+flux_add(l,k,1)*flux_adj(l,k,1)+flux_add(l,k,2)*flux_adj(l,k,2)
         enddo
      enddo
      fbrho_return=-sumn*sumd/betatavg

      return
      end function fbrho

#ifdef tuan_tr_test
      subroutine calrho_hex(first)
      use Mod_GetSome
      use Mod_SolLS, only: AxBHEX 
      use Inc_TPEN
      implicit none
      logical(4) :: first
      integer :: k, l, m, i
      real(8) :: sumn, sumd, betatavg, betatadj, psil2, sumv, vol
      real(8) :: rvdtvol, rhoadjd, pcfact, temp
      integer :: next=1
      integer :: KX,MX,NNX,KY,MY,NNY,KZ,MZ,NNZ
      real(8) :: BRANGE(2),ERANGE(2),HRANGE(2)
      integer :: IX, IY, IZ, N, NN, LL
      real(8) :: T(50), tmp_beta(6)=0d0, tmp_lamb(6)=0d0
      real(8) :: tmp_BU, tmp_ENR, tmp_H2U, tmp_VOL
      integer :: IFTY

      call AxBHex(flux,flux_add)
      sumn=0d0
      sumd=0d0
      sumv=0d0
      betatavg=0d0
      betatadj=0d0
      psil2=0d0
      do i=1,N_Group_d
         Pzeta(i,next)=0d0
         Plambda(i,next)=0d0
         Pbeta(i)=0d0
      enddo

      do k=1,nz
         do l=1,nxy
            !! === INDEX ============= !!
            !! flux_add : square index !! 
            !! flux     : square index !!
            !! rvdelt   : square index !!
            !! FisSrc   : square index !!
            !! betap    : square index !!
            !! ======================= !!
            vol=NodeVolume(l,k)
            do m=1,2
               rvdtvol=vol*rvdelt(m,l,k)
               flux_add(l,k,m)=flux_add(l,k,m)-rvdtvol*flux(l,k,m) &
                  & +maxs_chi_3d(l,k,m)*(betap(l,k)-1d0)*FisSrc(l,k)
               sumn=sumn+flux_add(l,k,m)*flux_adj(l,k,m)
! sumd = equation (9.68)               
               sumd=sumd+maxs_chi_3d(l,k,m)*FisSrc(l,k)*flux_adj(l,k,m)
! sumv = numerator of equation (9.70)
               sumv=sumv+flux(l,k,m)*v_inv(l,k,m)*flux_adj(l,k,m)*vol
! betatadj = numerator of equation (9.73)
               betatadj=betatadj+beta_d_tot(l,k)*FisSrc(l,k) &
                  & *flux_adj(l,k,m)*maxs_chid_3d(l,k,m)
!psil2 = equation (9.68) 
               psil2=psil2+FisSrc(l,k)*flux_adj(l,k,m)*maxs_chi_3d(l,k,m)
            enddo
            betatavg=betatavg+beta_d_tot(l,k)*FisSrc(l,k)
            do i=1,N_Group_d
               temp=flux_adj(l,k,1)*C_d_I(l,k,i)
! Pzeta = equation (9.71)               
               Pzeta(i,next)=Pzeta(i,next)+temp
! Plambda = numerator of equation (9.71)
               Plambda(i,next)=Plambda(i,next)+temp*lambda_d(l,k,i)
! Pbeta = numerator of equation (9.73)
               Pbeta(i)=Pbeta(i)+flux_adj(l,k,1)*FisSrc(l,k)*beta_d_eff(l,k,i)
            enddo
         enddo
      enddo
!
      betatavg=betatadj/psil2
      
      rhoadjd = rhoadj
      sumd=1d0/sumd
!      
      rhoadj=-sumn*sumd/betatavg
      iPgenT(next)=sumv*sumd
   
      pcfact=ppower/sumv*sumv0
 !      Pbeta
      do i=1,N_Group_d
         Pbeta(i)=Pbeta(i)*sumd
         Plambda(i,next)=Plambda(i,next)/Pzeta(i,next)
      enddo
   
      if (flag_print_kp_anc) then
         KX=5; MX=0; NNX=4
         KY=5; MY=0; NNY=4
         KZ=2; MZ=0; NNZ=1
         BRANGE(1) = 1.0D-3; BRANGE(2) = 6.0D+3
         ERANGE(1) = 0.5D+0; ERANGE(2) = 6.0D+0
         HRANGE(1) = 4.0D+0; HRANGE(2) = 5.0D+0
   
         tmp_BU = 0d0; tmp_VOL = 0d0
         do l=1,Nxy
            do k=IzFuelBot,IzFuelTop
               tmp_BU=tmp_BU+BU(l,k)*h_x(l)*h_y(l)*h_z(k)
               tmp_VOL=tmp_VOL+h_x(l)*h_y(l)*h_z(k)
            enddo
         enddo
         tmp_BU=tmp_BU/tmp_VOL
   
         tmp_ENR=4.5D0
   
         IFTY=int((tmp_BU+1d-1),4)
         IFTY=MAX(IFTY,1)
         tmp_H2U=1D0*COEF_TYPEF(IFTY)/tmp_ENR
         
         tmp_BU =MAX(BRANGE(1),MIN(tmp_BU ,BRANGE(2)))
         tmp_ENR=MAX(ERANGE(1),MIN(tmp_ENR,ERANGE(2)))
         tmp_H2U=MAX(HRANGE(1),MIN(tmp_H2U,HRANGE(2)))
   
         T=0d0
         N=1
         do IZ=MZ,NNZ
            do IY=MY,NNY
               do IX=MX,NNX
                  T(N)=(tmp_BU**IX)*(tmp_ENR**IY)*(tmp_H2U**IZ)
                  N=N+1
               enddo
            enddo
         enddo
         NN=KX*KY*KZ
         do LL=1,NN
            tmp_beta(1)=tmp_beta(1)+T(LL)*COEF_FB1(LL)
            tmp_beta(2)=tmp_beta(2)+T(LL)*COEF_FB2(LL)
            tmp_beta(3)=tmp_beta(3)+T(LL)*COEF_FB3(LL)
            tmp_beta(4)=tmp_beta(4)+T(LL)*COEF_FB4(LL)
            tmp_beta(5)=tmp_beta(5)+T(LL)*COEF_FB5(LL)
            tmp_beta(6)=tmp_beta(6)+T(LL)*COEF_FB6(LL)
            tmp_lamb(1)=tmp_lamb(1)+T(LL)*COEF_FL1(LL)
            tmp_lamb(2)=tmp_lamb(2)+T(LL)*COEF_FL2(LL)
            tmp_lamb(3)=tmp_lamb(3)+T(LL)*COEF_FL3(LL)
            tmp_lamb(4)=tmp_lamb(4)+T(LL)*COEF_FL4(LL)
            tmp_lamb(5)=tmp_lamb(5)+T(LL)*COEF_FL5(LL)
            tmp_lamb(6)=tmp_lamb(6)+T(LL)*COEF_FL6(LL)
         enddo
      endif
  
      if (first) then
         if (flag_print_kp_anc) then
            write(W_PKD,*) 'ANC Kinetics Parameter ******************'
            write(W_PKD,'(A10,F6.3,A6,F6.3,A6,F6.3)') ' Input BU=',tmp_BU, ', ENR=',tmp_ENR, ', H2U=',tmp_H2U
            write(W_PKD,*) 'Group      beta      lambda'
            do LL=1,6
               write(W_PKD,'(I4, 2F12.6)') LL, tmp_beta(LL), tmp_lamb(LL)
            enddo
            write(W_PKD,*) 'ANC Kinetics Parameter ******************'
            write(W_PKD,*) ''
         endif
   
         write(W_PKD,*) Div_Line2(1:70)
         write(W_PKD,*) "Point Kinetics Data"
         write(W_PKD,*) ""
         write(W_PKD,*) "[ Reactivity Coefficients ]"
         write(W_PKD,*) "[ Kinetics Parameters ]"
         write(W_PKD,*) "[ Time Dependent Reactivity ]"
         write(W_PKD,*) ""
         
         call fdbkrho_hex(betatavg,sumd,rhoadj,.true.) 
   
         write(W_PKD,*) Div_Line2(1:70)
         write(W_PKD,*) "[ Kinetics Parameters ]"
         write(W_PKD,*) ""
         write(W_PKD,*) "  Group    Beta      Lambda[1/sec]"
         do i=1,N_Group_d
            write(W_PKD,'(I6,f12.6,f14.6)') i, Pbeta(i), Plambda(i,next)
         end do
         write(W_PKD,'(A,f10.6)') "   Total", sum(Pbeta)
         write(W_PKD,*) ""
         
         write(W_PKD,*) Div_Line2(1:70)
         write(W_PKD,*) "[ Time Dependent Reactivity ]"
         write(W_PKD,*) ""
         !write(W_PKD,*) "    Time[sec]       rho[$]     Gen. Time[sec]     plevel"
         write(W_PKD,'(4a13,a20,a80)') "Time [s]","rho [$]","plevel","Gen. T [s]", &
                      & "beta(1:n_group_d)","lambda(1:n_group_d)"
      endif
      !WRITE(W_PKD,'(1p,20e15.6)') Sec, rhoadj, iPgenT(next), ppower
      WRITE(W_PKD,'(1p,16e13.5)') Sec, rhoadj, ppower/pcfact, iPgenT(next), &
                                & Pbeta(:), Plambda(:,next)
   
      call fdbkrho_hex(betatavg,sumd,rhoadj,.false.)
   
      return
      end subroutine calrho_hex

      subroutine fdbkrho_hex(betatavg,sumd,rhoadjn,first)  
      use Inc_TH
      use Inc_File
      use Inc_CR
      use Inc_INP, only: Flag_Card_maXS
      use Mod_Alloc
      use inc_tpen

      implicit none
      logical(4) :: first
      real(8) :: toutavgt, tfmaxt
      real(8) :: rhon, rhoadjn, rhod, betatavg, sumd
      integer :: l, k, m
      real(8) :: tfrho
      real(8) :: dmrho
      real(8) :: ppmrho
      real(8) :: crrho
      real(8) :: sumrho
      real(8) :: xesmrho
      real(8) :: leakrho
      real(8) :: nullrho
      real(8) :: tmp_delta
      logical(1), save :: flag_alloc_first=.true.
      integer(4) :: ixy, iz, ixy_fa
      !integer(4) :: ih, it, isfc

   
      if (flag_alloc_first) then
         flag_alloc_first=.false.
         IF(flag_thfb) then
             WRITE(W_RHO,'(1x,13a15)') 'time', 'sumrho', 'tfrho', 'dmrho',&
                                 & 'crrho', 'leakrho', 'nullrho'
         ELSE
             WRITE(W_RHO,'(1x,11a15)') 'time', 'sumrho', &
                                 & 'crrho', 'leakrho', 'nullrho'
             
         ENDIF
         
      endif
      
      if (.not.allocated(bk_CR_Bot))     call alloc (bk_CR_Bot ,n_cr)
      if (.not.allocated(bk_CR_Top))     call alloc (bk_CR_Top ,n_cr)
      if (.not.allocated(bk_BotFol_Bot)) call alloc (bk_BotFol_Bot ,n_cr)
      if (.not.allocated(bk_tFuel )) call alloc (bk_tFuel  ,nxy,nz)
      if (.not.allocated(bk_tMod  )) call alloc (bk_tMod   ,nxy,nz)
      if (.not.allocated(bk_dMod  )) call alloc (bk_dMod   ,nxy,nz)
      if (.not.allocated(bk_amd   )) call alloc (bk_amd    ,N_group,nxy,nz)
      if (.not.allocated(bk_afd   )) call alloc (bk_afd    ,n_group,nxy,nz)
      if (.not.allocated(bk_af2d  )) call alloc (bk_af2d   ,nxy,nz)
      if (.not.allocated(bk_scatd )) call alloc (bk_scatd  ,nxy,nz)
      if (.not.allocated(bk_ccwd  )) call alloc (bk_ccwd   ,N_group,nxy,nz)
      if (.not.allocated(bk_cced  )) call alloc (bk_cced   ,N_group,nxy,nz)
      if (.not.allocated(bk_ccnd  )) call alloc (bk_ccnd   ,N_group,nxy,nz)
      if (.not.allocated(bk_ccsd  )) call alloc (bk_ccsd   ,N_group,nxy,nz)
      if (.not.allocated(bk_ccbd  )) call alloc (bk_ccbd   ,N_group,nxy,nz)
      if (.not.allocated(bk_cctd  )) call alloc (bk_cctd   ,N_group,nxy,nz)
      if (.not.allocated(bk_dnwd  )) call alloc (bk_dnwd   ,N_group,nxy+nx+ny,nz)
      if (.not.allocated(bk_dnnd  )) call alloc (bk_dnnd   ,N_group,nxy+nx+ny,nz)
      if (.not.allocated(bk_dnbd  )) call alloc (bk_dnbd   ,N_group,nxy,nz+1)
!      
      if (.not.allocated(bk_xsxead)) call alloc (bk_xsxead ,N_group,nxy,nz)
      if (.not.allocated(bk_xssmad)) call alloc (bk_xssmad ,N_group,nxy,nz)
!      
      if (.not.allocated(bk_rnxed )) call alloc (bk_rnxed  ,nxy,nz)
      if (.not.allocated(bk_rnsmd )) call alloc (bk_rnsmd  ,nxy,nz)
      
      toutavgt = toutavg
      tfmaxt   = tfmax

      !! === square === !!
      !!@@ bk_tFuel = T_Fuel
      !!@@ bk_tMod  = T_Mod
      !! === hex === !!
      if (flag_thfb) then
          bk_tFuel = 0d0
          bk_tMod  = 0d0
          do Iz = IzFuelBot, IzFuelTop
             do Ixy_FA = 1,nxy_FA
                Ixy = I_FA(Ixy_FA)
                if (I_FARF_1N(Ixy) == 2) then
                   bk_tFuel(Ixy,Iz) = tdopl(Iz,Ixy_FA)
                   bk_tMod (Ixy,Iz) = tcool(Iz,Ixy_FA)
                   bk_dMod (Ixy,Iz) = dcool(Iz,Ixy_FA)
                   !! Initilization of temperture (Fuel, Moderator)
                   T_Fuel  (Ixy,Iz) = tdopl(Iz,Ixy_FA)*tdopl(Iz,Ixy_FA) !-273.15
                   T_Mod   (Ixy,Iz) = tcool(Iz,Ixy_FA)+273.15
                endif
             enddo
          enddo
      endif
      
          
!      if (Flag_Card_maXS) then
!         bk_dMod=D_Mod
!#ifdef jr_tr
!         bk_dmod=dcool
!#endif
!      endif
      bk_ppm    = PPM
      if (n_cr>0) bk_CR_Bot = Abs_Bot
      if (n_cr>0) bk_CR_Top = TopFol_Bot
      if (n_cr>0) bk_BotFol_Bot = BotFol_Bot
      
      ! Xe/Sm number density <- need to be added
!      do k=1,nz
!         do l=1,nxy
!            do m=1,2
!               bk_xsxead(m,l,k)=miXS_a_Xe35(l,k,m)
!               bk_xssmad(m,l,k)=miXS_a_Sm49(l,k,m)
!            enddo
!            bk_rnxed(l,k)=N_Xe35(l,k)
!            bk_rnsmd(l,k)=N_Sm49(l,k)
!         enddo
!      enddo
      
      if (first) then
         rhocoef = 0
   
         ! Compute boron reactivity coeficient
!         tmp_delta = 10d0
!         PPM = PPM + tmp_delta
!         rhocoef(1) = (fbrho_hex(betatavg,sumd)-rhoadjn)/tmp_delta
!         PPM = bk_ppm
         
         ! Compute coolant density coeficient
         if (flag_thfb) then

             if (Flag_Card_maXS) then
                tmp_delta = 1.d-3
                !d_mod = d_mod + tmp_delta
                dcool = dcool + tmp_delta
                rhocoef(2) = fbrho_hex(betatavg,sumd)-rhoadjn
                !end if
                do Iz = IzFuelBot, IzFuelTop
                   do Ixy_FA = 1,nxy_FA
                      Ixy = I_FA(Ixy_FA)
                      if (I_FARF_1N(Ixy) == 2) then
                         dcool(Iz,Ixy_FA) = dcool0(Iz,Ixy_FA) 
                         !!dcool(Iz,Ixy_FA) = bk_dMod (Ixy,Iz)
                      endif
                   enddo
                enddo
             endif
             
             ! Compute coolant temperature coeficient
             !!t_mod = t_mod + tmp_delta
             tmp_delta = 10.d0
             tcool = tcool + tmp_delta
             rhocoef(3) = (fbrho_hex(betatavg,sumd)-rhoadjn)/tmp_delta
             !! t_mod = bk_tMod
             !! === hex === !!
             do Iz = IzFuelBot, IzFuelTop
                do Ixy_FA = 1,nxy_FA
                   Ixy = I_FA(Ixy_FA)
                   if (I_FARF_1N(Ixy) == 2) then
                      tcool(Iz,Ixy_FA) = tcool0 (Iz,Ixy_FA)
                      !t_mod(Ixy,Iz)    = tcool0 (Iz,Ixy_FA)
                      !!tcool(Iz,Ixy_FA) = bk_tMod (Ixy,Iz)
                   endif
                enddo
             enddo
             
             ! Compute doppler temperature coeficient
             tmp_delta = 10.d0
             !! t_fuel = t_fuel + tmp_delta
             do k = IzFuelBot, IzFuelTop
                do l = 1,nxy_fa
                   tdopl(k,l) = dsqrt(tdopl(k,l)*tdopl(k,l)+10d0) !tdopl + tmp_delta !!
                enddo
             enddo
             rhocoef(4) = (fbrho_hex(betatavg,sumd)-rhoadjn)/tmp_delta
             !! T_Fuel = bk_tFuel
             !! === hex === !!
             do Iz = IzFuelBot, IzFuelTop
                do Ixy_FA = 1,nxy_FA
                   Ixy = I_FA(Ixy_FA)
                   if (I_FARF_1N(Ixy) == 2) then
                      tdopl(Iz,Ixy_FA) = tdopl0(Iz,Ixy_FA) 
                      !t_fuel(Ixy,Iz)   = tdopl0(Iz,Ixy_FA)
                      !!tdopl(Iz,Ixy_FA) = bk_tFuel(Ixy,Iz)
                   endif
                enddo
             enddo
             
             
             ! Compute xe/sm component of reactivity
!             do k=1,nz
!                do l=1,nxy
!                   N_Xe35(l,k)=N_Xe35(l,k)+1d+10
!                enddo
!             enddo
!             rhocoef(5)=(fbrho_hex(betatavg,sumd)-rhoadjn)*1d-10
!             do k=1,nz
!                do l=1,nxy
!                   N_Xe35(l,k)=bk_rnxed(l,k)
!                enddo
!             enddo
!             
!             do k=1,nz
!                do l=1,nxy
!                   N_Sm49(l,k)=N_Sm49(l,k)+1d+12
!                enddo
!             enddo
!             rhocoef(6)=(fbrho_hex(betatavg,sumd)-rhoadjn)*1d-12
!             do k=1,nz
!                do l=1,nxy
!                   N_Sm49(l,k)=bk_rnsmd(l,k)
!                enddo
!             enddo
             
             write(W_PKD,*) Div_Line2(1:70)
             write(W_PKD,*) "[ Reactivity Coefficients ]"
!             write(W_PKD,'(A,Es12.5,A)') "  - Boron reactivity worth : ", rhocoef(1), " [$/ppm]"
             if (Flag_Card_maXS) then
                write(W_PKD,'(A,Es12.5,A)') "  - Coolant density worth  : ", rhocoef(2), " [$/kgm^-3]"
             endif
             write(W_PKD,'(A,Es12.5,A)') "  - Coolant temp. worth    : ", rhocoef(3), " [$/K]"
             write(W_PKD,'(A,Es12.5,A)') "  - doppler temp. worth    : ", rhocoef(4), " [$/K]"
!             write(W_PKD,'(A,Es12.5,A)') "  - Xenon worth            : ", rhocoef(5), " [$/cm^-3]"
!             write(W_PKD,'(A,Es12.5,A)') "  - Samarium worth         : ", rhocoef(6), " [$/cm^-3]"
             write(W_PKD,*) ""
             return
         
         else
             return
         endif
      endif
   
      !! hexagonal geometry
      dcmatd=dcmat
      cmatd =cmat
      dhatd =dhat
      dhatzd=dhatz
      bk_afd=af
      !! hexagonal geometry
   
      rhon = rhoadjn

      ! compute "nodal leakage" component of reactivity.
      !! hexagonal geometry
      rhod = fbrho_hex(betatavg,sumd)
      dhat = dhat0
      dhatz = dhatz0

      !! hexagonal geometry
      rhod = rhon
      rhon = fbrho_hex(betatavg,sumd)
      leakrho = rhod-rhon
     
      ! Compute control component of reactivity
!Write(*,*) "Tuan: need to modify for control rod reactivity: fdbkrho_hex"
      if (n_cr>0) Abs_Bot = CR_Bot0
      if (n_cr>0) TopFol_Bot = CR_Top0
      if (n_cr>0) BotFol_Bot = BotFol_Bot0

      rhod   = rhon
      rhon   = fbrho_hex(betatavg,sumd)
      crrho  = rhod - rhon
      nullrho = rhon

      ! Compute boron component of reactivity
!      ppm    = ppm0
!      rhod   = rhon
!      rhon   = fbrho_hex(betatavg,sumd)
!      ppmrho = rhod - rhon
  
      !! Compute moderator density/temperature component of reactivity
      !! === For lattice === !!
      !t_mod = init_tMod
      !if (Flag_Card_maXS) then
      !   d_mod = init_dMod
      !endif
      !! =================== !!
      if (flag_thfb) then
          do Iz = IzFuelBot, IzFuelTop
             do Ixy_FA = 1,nxy_FA
                Ixy = I_FA(Ixy_FA)
                if (I_FARF_1N(Ixy) == 2) then
                   tcool(Iz,Ixy_FA) = tcool0 (Iz,Ixy_FA)
                   dcool(Iz,Ixy_FA) = dcool0 (Iz,Ixy_FA)
                   !t_mod(Ixy,Iz)    = tcool0 (Iz,Ixy_FA)
                endif
             enddo
          enddo
          rhod  = rhon
          rhon  = fbrho_hex(betatavg,sumd)
          dmrho = rhod - rhon
          
          ! compute doppler temperature component of reactivity
          !! === For lattice === !!
          !t_fuel = init_tFuel
          !! =================== !!
          do iz = IzFuelBot, IzFuelTop
             do ixy_fa = 1,nxy_fa
                ixy = i_fa(ixy_fa)
                if (i_farf_1n(ixy) == 2) then
                   tdopl(iz,ixy_fa) = tdopl0 (iz,ixy_fa)
                endif
             enddo
          enddo
          rhod  = rhon
          rhon  = fbrho_hex(betatavg,sumd)
          tfrho = rhod - rhon
          nullrho = rhon

          ! compute xe/sm component of reactivity
!          do k=1,nz
!             do l=1,nxy
!                do m=1,2
!                   miXS_a_Xe35(l,k,m)=xsxea0(l,k,m)
!                   miXS_a_Sm49(l,k,m)=xssma0(l,k,m)
!                enddo
!                N_Xe35(l,k)=rnxe0(l,k)
!                N_Sm49(l,k)=rnsm0(l,k)
!             enddo
!          enddo
!          rhod    = rhon
!          rhon    = fbrho_hex(betatavg,sumd)
!          xesmrho = rhod-rhon
      endif
      
      if (flag_thfb) then

          sumrho = tfrho + dmrho + crrho + leakrho + nullrho
          
          ! print reactivity components
          WRITE(W_RHO,601) sec, sumrho, tfrho, dmrho,  &
                      & crrho,  leakrho, nullrho
    !#ifdef jr_tr
    !      WRITE(*,601) sec, sumrho, tfrho, dmrho, ppmrho, &
    !                  & crrho, xesmrho, leakrho, nullrho
    !#endif
601 FORMAT(1x,f15.8,7(2x,1p,e13.6))
      else
          sumrho = crrho +  leakrho + nullrho
          
          ! print reactivity components
          WRITE(W_RHO,602) sec, sumrho,  &
                      & crrho, leakrho, nullrho
    !#ifdef jr_tr
    !      WRITE(*,601) sec, sumrho, tfrho, dmrho, ppmrho, &
    !                  & crrho, xesmrho, leakrho, nullrho
    !#endif
602 FORMAT(1x,f15.8,5(2x,1p,e13.6))
      endif
      
      ! reset matrix and feedback variables to current time step values
      !! === For lattice ====== !!
      !!@@ T_Fuel = bk_tFuel
      !!@@ T_Mod  = bk_tMod
      !!@@ if (Flag_Card_maXS) then
      !!@@    D_Mod = bk_dMod
      !!@@ endif
      !! ====================== !!
      if(flag_thfb) then
          
          do iz = IzFuelBot, IzFuelTop
             do ixy_fa = 1,nxy_fa
                ixy = i_fa(ixy_fa)
                if (i_farf_1n(ixy) == 2) then
                   tdopl(iz,ixy_fa) = bk_tFuel(ixy,iz) 
                   tcool(iz,ixy_fa) = bk_tMod (ixy,iz)
                   dcool(iz,ixy_fa) = bk_dMod (ixy,iz)
                endif
             enddo
          enddo
!          PPM = bk_ppm

!          do k=1,nz
!             do l=1,nxy
!                do m=1,2
!                   miXS_a_Xe35(l,k,m)=bk_xsxead(m,l,k)
!                   miXS_a_Sm49(l,k,m)=bk_xssmad(m,l,k)
!                enddo
!                N_Xe35(l,k)=bk_rnxed(l,k)
!                N_Sm49(l,k)=bk_rnsmd(l,k)
!             enddo
!          enddo
          toutavg=toutavgt
          tfmax=tfmaxt
      endif
      if (n_cr>0) Abs_Bot = bk_CR_Bot
      if (n_cr>0) TopFol_Bot = bk_CR_Top
      if (n_cr>0) BotFol_Bot = bk_BotFol_Bot
          
   
      !! hexagonal geometry   
      dcmat = dcmatd
      cmat  = cmatd
      dhat  = dhatd
      dhatz = dhatzd
      af    = bk_afd
      !! hexagonal geometry   


      return
      end subroutine fdbkrho_hex


      function fbrho_hex(betatavg, sumd) result(fbrho_return)
      use Mod_SolLS,   only: AxBHEX 
      use Mod_SolTPEN, only: setlshex
      use Mod_XSFB !,    only: XSFBHEX, macroXSFBHEX, macroXSFBHEX_MG
      use Mod_SolFDM,  only: FdmCoupl
      use Inc_maXS,    only: nu_maXS_f_3D
#ifdef tuan_tr_test

#else
      use mod_interface, only: adjust_force_th, return_force_th
#endif
      use Inc_tpen
      use Inc_INP,     only: Flag_Card_maXS
      implicit none
      integer :: k, l
      real(8) :: sumn, sumd, betatavg, rvdtvol(N_group), vol, psixs
      real(8) :: fbrho_return
#ifdef tuan_tr_test

#else
      call adjust_force_th
#endif
      if (Flag_Card_maXS) then
          call maXSFBHEX
      else
         call XSFBHEX
      endif
#ifdef tuan_tr_test

#else
      call return_force_th
#endif
      call FdmCoupl
      call setlshex(1)
      call AxBHEX(flux,flux_add)
     
      sumn=0
      do k=1,nz
         do l=1,nxy
            !! === INDEX ============= !!
            !! rvdelt   : square index !!
            !! flux     : square index !!
            !! flux_add : square index !!
            !! ======================= !!
            vol=nodevolume(l,k)
            rvdtvol(1)=vol*rvdelt(1,l,k)
            rvdtvol(2)=vol*rvdelt(2,l,k)
            psixs=(nu_maXS_f_3D(l,k,1)*flux(l,k,1)+nu_maXS_f_3D(l,k,2)*flux(l,k,2))*vol
            flux_add(l,k,1)=flux_add(l,k,1)-rvdtvol(1)*flux(l,k,1)+(betap(l,k)-1d0)*psixs
            flux_add(l,k,2)=flux_add(l,k,2)-rvdtvol(2)*flux(l,k,2)
            sumn=sumn+flux_add(l,k,1)*flux_adj(l,k,1)+flux_add(l,k,2)*flux_adj(l,k,2)

         enddo
      enddo
      fbrho_return=-sumn*sumd/betatavg
     
      return
      end function fbrho_hex

#endif      
      
      END MODULE Mod_CalReact


#endif
