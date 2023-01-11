#ifdef jr_vver
      MODULE Mod_SStpen
      USE Inc_Constant
      USE Inc_Flag
      USE Inc_Geometry
      USE Inc_RP
      USE Inc_Option
      use inc_parallel, only: comm
      USE Inc_Lscoef
      USE Inc_Control
      USE Inc_TPEN
      USE Inc_maXS
      USE Mod_SolFDM
      USE Mod_Save
      USE Mod_SolTPEN
      USE Mod_TPENDrive
      USE Mod_SolLU
!      USE Mod_Soldhat
      USE Mod_SolLS
#ifdef js_mpi
      use mod_parallel, only: bcast, barrier
#endif
#ifdef tuan_tr_test
      USE Inc_INP, only: Flag_Card_maXS, Flag_Card_CR
#endif

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE TpenSolver_St
      !USE Inc_maXS
      USE Inc_FluxVar
      USE Mod_Wielandt
      USE Mod_SolLS
      !USE Inc_Lscoef
      USE Mod_SolLU
!      USE Mod_Soldhat
      USE Mod_THdriver
      USE Mod_SolFDM
      USE Mod_InitCond, only: GuessIteration
!      USE Mod_Interface
      USE Inc_3D, ONLY: Flux, Power
      !USE Inc_Control
      USE Inc_Option
      use mod_charedit, only: print_msg
      USE Inc_TH, ONLY: Avg_CorePower, LevFactor, tfmax
      USE Inc_XS_File, only: flag_leakage
!      USE Mod_LeakCorr, ONLY: XSFB_LC
      USE Mod_XSFB !, ONLY: XSFBHEX_MG,mgto2g_hex,XSFBHEX
     ! USE Mod_GetSome, only: levelbalancing
      USE Mod_GetSome
!    !  USE Mod_Depletion, ONLY: FissionProduct
      use Inc_Crud, only: opt_crud, opt_crudasi
      !use Inc_XYZ, only: ASI
!     ! USE Mod_Depletion, ONLY: Crud_Itr_ASIsearch,Crud_itr_CBCCOR,Crud_Itr_FAASIsearch
!      use mod_extsrc, only: set_fsrc_extsrc, fatal_error


#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction
      use Mod_PinPow_Hex, only: calculate_pin_power_hex
      use Inc_PinPow_Hex, only: pin_pow_hex_needed
      ! -=-=-=-=-=-=-=-=-
#endif

      IMPLICIT NONE
      !
      integer :: ixy_fa
      INTEGER :: Iout, Ixy, Iz
      INTEGER :: ioutbeg, ioutend, nskipnodal
      INTEGER :: icy_CMFD, icy_TH
!      REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: Flux_temp
      !REAL(8) :: dxs, dxsn
      SAVE ioutbeg,nskipnodal
      DATA ioutbeg,nskipnodal/1,1/
      logical(1) :: flag_crudasi_conv
!      CALL XSFB

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [TpenSolver_St] in Mod_SStpen'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif
#ifdef tuan_tr_test
! Tuan: Need to modify XSFB
      CALL XSFBHEX
#endif
      
      write(*,'(A)') '* Leakage Correction1 *'
      !if ( flag_thfb ) then
      !  CALL hexth
      !else
#ifdef tuan_fr
!         DO Iz = IzFuelBot, IzFuelTop
!             DO Ixy_FA = 1, Nxy_FA
!                 Ixy = I_FA(Ixy_FA)
!                 call  update_maXS(Ixy,Iz)
!             END DO
!         END DO
!         CALL XSFBHEX_MG
#endif
       
      
      dnw = D0
      dnn = D0
      dnb = D0
      CALL GuessIteration
      CALL FdmCoupl
      CALL setls
      CALL ilufac2d

      ifnodal  = .TRUE.
      nodal    = .TRUE.
      iflsupd  = .FALSE.
      ioutbeg  = 1
      Ioutend  = Iout_Max
      icy_CMFD = - Period_CMFD
      icy_TH   = - Period_TH
      flagl2   = .FALSE.
      flaglinf = .FALSE.
      flageig  = .FALSE.
      flagerf  = .FALSE.
      Flag_Conv_PowLev = .TRUE.

         flag_crudasi_conv=.True.
      nintot = 0


      DO Iout = 1,1 !Ioutbeg, Ioutend
         !call xs_tpen !!!!!! NEED TO CHECK
         call eignv_tpen
!            dummy_filler = 1 ! @$^ siarhei_plot 
         !call scarpget('phi', Flux(0:Nxy-1,:,:))
         !call scarpget('eigv', keff)
#ifdef tuan_fr
         CALL Get_POW
#ifdef tuan_tr_test
         if(opt_mode == 3) CALL Get_LinPOW
#endif

         CALL Get_Avg(Avg_Power, Power, 0)
         LevFactor = Avg_CorePower/Avg_Power

         if(opt_mode /= 3) call LevelBalancing
#ifdef tuan_tr_test
         if(opt_mode == 3) then
             CALL Get_POW
             CALL Get_LinPOW
             CALL Get_Avg(Avg_Power, Power, 0)
         endif
         
#endif
         CALL Get_normal_power
#endif

!#ifdef jr_vver
!         CALL Get_POW
!         CALL Get_LinPOW
!         CALL Get_Avg(Avg_Power, Power, 0)
!
!         DO Ixy_FA = 1, Nxy_FA
!            Ixy = I_FA(Ixy_FA)
!
!            DO Iz = IzFuelBot, IzFuelTop
!               Normal_Power(Ixy, Iz) = Power(Ixy, Iz) / Avg_Power
!            END DO
!         END DO
!#endif
         CALL XYZ_Indexing
!            dummy_filler = 1 ! @$^ siarhei_plot 


#ifdef EFFICIENT
#else
            ifxsupd = .TRUE.  ! Invoke XSFB
#endif
         IF ( (OPT_Xe /= 2) .OR. (OPT_Sm /= 2) ) THEN
#ifdef BRCAL
           ! CALL FissionProduct(.FALSE.)

#ifdef EFFICIENT
            if(ifnlupd) ifxsupd = .TRUE.  ! Invoke XSFB
#endif
#else
          !  IF ( .NOT. ( (OPT_Mode == 4))) then
          !     CALL FissionProduct(.FALSE.)
          !  END IF
#endif
         END IF
         IF ( ifxsupd ) THEN
#ifdef check_test
            ifcbcupd=.true.
#endif

!#ifdef Leak_Corr
!            IF (Flag_Leakage) THEN
!               write(*,'(A)') '* Leakage Correction *'
!!               CALL updatedhat
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif
!!               CALL XSFB_LC
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif
!               CALL FdmCoupl
!               CALL setls
!            ELSE
!               CALL XSFB
!            ENDIF
#ifdef tuan_tr_test
            if (Flag_Card_maXS) then
!                 CALL XSFBHEX_MG
                 CALL maXSFBHEX
            else
!               call XSFBHEX
            endif
#endif

            CALL Get_Source(FisSrc, nu_maXS_f_3D, Flux)
            CALL FdmCoupl

            iflsupd = .TRUE.
         END IF


      END DO

#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction
      ! run pin power reconstruction if needed
      ! at this point, corner/surface/etc fluxes are saved from Mod_TPENDrive
      if (pin_pow_hex_needed) then
         call calculate_pin_power_hex
      end if
     ! write(*,*) 'SStpen Siarhei pause'
     ! pause ! Siarhei pause (temp)

      ! -=-=-=-=-=-=-=-=-
#endif

      RETURN
      END SUBROUTINE TpenSolver_St

!!!#ifdef siarhei_delete ! ! check later siarhei_rev
      SUBROUTINE eignv_tpen
      USE Mod_SolTPEN

!!!#ifdef siarhei_ppr
!!!      ! Siarhei pin power reconstruction
!!!      use Inc_PinPow_Hex
!!!      use Mod_PinPow_Hex
!!!      ! -=-=-=-=-=-=-=-=-
!!!#endif

      IMPLICIT NONE
!      call p_init ! NEED TO CHECK
      call DtilHex !is included in subroutine, FdmCoupl
      call setlshex(1)
      call ilufachex
      call outerss

      RETURN
      END SUBROUTINE eignv_tpen
!!!#endif 

      SUBROUTINE outerss
      USE Inc_3D!,       ONLY: FisSrc_Iout, FisSrc, Flux, keff
      USE Inc_Fluxvar,  ONLY: src
      USE Mod_Wielandt, ONLY: wiel_hex
      USE Mod_XSFB!,     ONLY: XSFB,XSFBHEX,XSFBHEX_MG
      USE Inc_TH!,   ONLY: ppm
#ifdef tuan_frth
!      USE Inc_Control, ONLY:  Period_TH, Period_CMFD
!      use Mod_Interface
      USE Mod_THdriver
      USE Mod_GetSome
#endif
#ifdef tuan_tr_test
      use Mod_GetSome, only: Get_Reg_VF
#endif
#ifdef tuan_fr_crm
     USE Inc_CR
#endif

      IMPLICIT NONE
! outer iteration for the steady state calculation
!      CHARACTER(len=79) :: amesg
!      CHARACTER(len=65) ::  amesg65
!      LOGICAL :: iflsupd,ifnlupd,ifxsupd, skip
!      LOGICAL :: ifnlupd,ifxsupd, skip
      LOGICAL,save :: firstout=.true.
      INTEGER(4) :: k, l, m, iout, ioutbeg, nskipnodal
      INTEGER(4) :: ioutend, icy, iin, lf
      REAL(8) ::  dxs, dxsn, reigvdel, fs
      SAVE ioutbeg,nskipnodal
      DATA ioutbeg,nskipnodal/1,1/
      INTEGER(4) :: mf
!
      INTEGER(4) :: ind
#ifdef tuan_frth
      INTEGER :: icy_CMFD, icy_TH

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [outerss] in Mod_SStpen'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

! update          Period_TH 
      Period_TH   = 5
      icy_CMFD = - Period_CMFD
      icy_TH   = - Period_TH
#endif
!      INTEGER(4) :: iptype
!
! reduce the convergence criteria if performing SS initialization
!
!     if (extth .AND. .NOT.tran) ssinit = .TRUE. !! need to modify
!     when add the extth, please modify this ...
!
      IF (extth .AND. ssinit) THEN                            
!         dxs=0
!         DO k=1,nz
!            DO l=1,nxy
!               ind = imap(l)
!               DO m=1,2
!               !DO m=1,n_group
!                  IF (maxs_r_3d(ind,k,m).NE.0.0) THEN
!                     dxsn=ABS(maxs_r_3d_old(ind,k,m)/maxs_r_3d(ind,k,m)-1)
!                     dxs=MAX(dxs,dxsn)
!                  ELSE
!                     dxs=1.0
!                  ENDIF
!               ENDDO
!            ENDDO
!         ENDDO
         ! epstf = epstft
         IF(nskipnodal.GE.nodalcy) THEN
         !IF(dxs.GT.epstft .OR. nskipnodal.GE.nodalcy) THEN
            ifnodal=.TRUE.
            ioutend=ioutbeg+nodalcy-1 !temp
            nskipnodal=0
         ELSE
            ifnodal=.FALSE.
            ioutend=ioutbeg
            nskipnodal=nskipnodal+1
         ENDIF
         IF(firstout)THEN
            icy=-ninitoutt
            firstout=.FALSE.
         ELSE
            icy=0
         ENDIF
      ELSE
         ifnodal=.TRUE.
         iflsupd = .FALSE.
         ioutbeg=1
         ioutend=iout_max
         icy=-ninitoutt
         flagl2=.false.
         flaglinf=.false.
         flageig=.false.
         flagerf=.false.
         flagth=.false.
      ENDIF
!               WRITE(amesg,599)
!
      iptype = 0 ! If you adds this option, please change this value
      DO iout=ioutbeg,ioutend
         icy=icy+1
         iflsupd=.false.
         ifnlupd=.false.
         reigvdel=reigv-reigvs
#ifdef tuan_frth
         if (flag_THFB) then
             icy_CMFD = icy_CMFD + 1
             icy_TH   = icy_TH + 1
             if (icy_TH>=Period_TH) icy_TH=mod(icy_TH,Period_TH)
             if (icy_CMFD>=Period_CMFD) icy_CMFD=mod(icy_CMFD,Period_CMFD)
         endif
#endif
         !IF(iptype.EQ.0) THEN
         DO k=1,nz
            DO l=1,nxy
               ind=imap(l)
               fs=FisSrc(ind,k)*reigvdel
               src(1,ind,k)=maxs_chi_3d(ind,k,1)*fs  !mg
               src(2,ind,k)=maxs_chi_3d(ind,k,2)*fs
            ENDDO
         ENDDO
         CALL SolveLinearSystem_hex(iin,.FALSE.)
         DO l=1,nxy
            ind=imap(l)
            DO k=1,nz
               DO m=1,2
               !DO m=1,n_group
                  IF(Flux(ind,k,m) .LE. 0)  THEN
                     Flux(ind,k,m) = 1.0e-30
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
!          WRITE(*,*) 'MAX MIN FLUX: ', maxval(Flux), minval(Flux)
         nintot=nintot+iin
         CALL wiel_hex(icy)
!         WRITE(*,*) '----------------------------------------------------------------------------------------------------'
!         WRITE(*,599)
         WRITE(*,601) iout,iin,keff,flageig,erreig,flagl2           &
              ,             rerrl2,flaglinf,domr, icy !,ppm
!         WRITE(*,*) ' '
!         CALL message(false,popt(2),popt(14),amesg)
!
!         WRITE(*,601) iout,iin,effk,flageig,rerrl2,flagl2           &
!              ,             errlinf,flaglinf,domr,ppm
! exit outer if convergence achieved
         flagneut=flageig .AND. flagl2 .AND. flaglinf
         skip = FALSE
! removed by Tuan         
!         IF (extth.AND.ssinit) THEN
!            IF (ifnodal) THEN
!               IF (icy.EQ.1.AND.flagneut) THEN
!                  skip = TRUE
!                  EXIT
!               END IF
!            ELSE
!               IF (flagneut) THEN
!                  skip = TRUE
!                  EXIT
!               END IF
!            ENDIF
!         ELSE
!            IF(.NOT.nodal .OR. icy.EQ.1) THEN
!               IF(flagneut) THEN
!                  skip = TRUE
!                  EXIT
!               END IF
!
!            ENDIF
!         ENDIF
!
            IF( icy.EQ.1) THEN
               IF(flagneut) THEN
!                   write(*,*) 'come here and see:'
!                   pause
                  skip = TRUE
                  EXIT
               END IF

            ENDIF
! .. invoke nodal update
         ifxsupd=.false.
         ifnlupd=(icy.EQ.nodalcy) .OR. (iout.EQ.ninitoutt) .OR. flagneut
!wRITE(1995,*) 'ITERATION      ', iout
         IF(nodal .AND. ifnlupd .AND. ifnodal) THEN
!
! .. rod cusping correction
            CALL decusp
! .. nodal kernel
            do l=1,nxy
               do k=1,nz
                  do m = 1,2
                     do mf = mgb(m),mge(m)
                        fluxf(l,k,mf)=flux(l,k,m)
                     enddo
                  enddo
               enddo
            enddo
            CALL TpenDriver
            if_cycle2  = .true.
! set zero(1.0e-30) for negative flux values
            DO l=1,nxy
               ind=imap(l)
               DO k=1,nz
                  DO m=1,2
                  !DO m=1,n_group
                     IF(flux(ind,k,m) .LE. 0)  THEN
                        flux(ind,k,m) = 1d-30
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
!  end of change
            !IF(multigroup) ifxsupd=.true.
            ! multigroup = .true.
            IF(.true.) ifxsupd=.true.
            iflsupd=.true.
            icy=0

         ENDIF

#ifdef tuan_frth
         IF( ifnlupd .and. flag_thfb) THEN

!        ifnlupd=(icy_TH==0).or.(iout==3)
!         th_signal=0
!         if (ifnlupd) then
            call Get_POW
!            call Get_LinPOW
            call Get_Avg(Avg_Power,Power,0)
            LevFactor = Avg_CorePower/Avg_Power
            if(opt_mode /= 3) call LevelBalancing


            !js+if ((abs(keff_Old-k_target)<0.001d0).and.(abs(LevFactor-D1)>EPS_PowLev)) then
            if (abs(LevFactor-1d0)>EPS_PowLev) then
               Flag_Conv_PowLev=.false.
            else
               Flag_Conv_PowLev=.true.
            endif

!            if(opt_mode == 3) then
                call Get_POW
                call Get_LinPOW
                call get_normal_power
!            endif
            
            if (flag_thfb) then

!            th_signal=1
                ifxsupd=.true.  ! Invoke XSFB
#ifdef tuan_tr_test
            CALL ssth_hex
#else
            CALL ssth
#endif

            endif

!         endif
         ENDIF
!
!!         if (flag_ftnth_detail) switch_saveth=.true.
!
!!         select case (th_signal)
!!         case (1)
!            call WATCH('TH Feedback','START')


!            call backup_TH
!            dummy_filler = 1 ! @$^ siarhei_plot 
!            call Interface_TH
!            dummy_filler = 1 ! @$^ siarhei_plot 
!!            call check_chth(flag_tavg_conv) !Damping
!            dummy_filler = 1 ! @$^ siarhei_plot 
!            call WATCH('TH Feedback','END')
!!         end select
!
!
!!         if (ifnlupd) then
!!            if (OPT_Mode/=4.and.opt_mode/=1.and..not.flag_rod) then
!!                  ifxsupd=.true.  ! Invoke XSFB
!!               if (opt_findtavg.and.rerrl2<1d-3) then
!!                  call Get_Avg(tmp_tm_avg,T_Mod,0)
!!                  call change_massflow(tmp_tm_avg,i_tavg_itr,flag_tavg_conv)
!            dummy_filler = 1 ! @$^ siarhei_plot 
!!               endif
!!            endif
!!         endif
#endif


!         lupscat = n_group
         IF(lupscat.EQ.1) ifxsupd=.true.
! .. update xsec and others associated with xsec change
         
         IF(ifxsupd) THEN
            !CALL XSFB    !! 56 pcm
#ifdef tuan_fr_crm
             if (flag_decusping .and. flagl2) then
                 flag_decusping_update = .TRUE.
                 call Get_Reg_VF_HEX
             
             endif 
#endif        

#ifdef tuan_tr_test
            CALL maXSFBHEX
#else
            CALL XSFBHEX_MG
#endif
            CALL DtilHex
            CALL DtilHex
            DO k=kfs,kfe
               DO lf=1,nfuel
                  l=nodef_hex(lf)
                  ind=imap(l)
                  FisSrc(ind,k)=0.
                  DO m=1,2
                  !DO m=1,n_group
                     FisSrc(ind,k)=FisSrc(ind,k)+nu_maXs_f_3D(ind,k,m)*flux(ind,k,m)
                  ENDDO
                  FisSrc(ind,k)=FisSrc(ind,k)*nodevolume(ind,k)
               ENDDO
            ENDDO
            iflsupd=.true.
         ENDIF

         !ELSEIF(scatto1)THEN
         !   DO k=1,nz
         !      DO l=1,nxy
         !         xss(l,k)=xsso(1,l,k)-xsso(2,l,k)*phi(2,l,k)/phi(1,l,k)
         !         xst(1,l,k)=xsa(1,l,k)+xss(l,k)
         !      ENDDO
         !   ENDDO
         !   iflsupd=true
         !ENDIF
!
! .. update linear system
         IF(iflsupd) THEN
            CALL setlshex(1)
!
! .. update preconditioner
!            CALL ilufac2d
            
            CALL IluFacHex
         ELSE
!
! correct diagonal components of the coefficient matrix
            reigvdel=reigvs-reigvsd
            IF(.not.if_hexgometry) THEN
               DO k=1,nz
                  DO l=1,nxy
                     am(1,l,k)=am(1,l,k)-af(1,l,k)*reigvdel
                     af2(l,k)=af(2,l,k)*reigvs
                  ENDDO
               ENDDO
            ELSE
               CALL setlshex(2)
            ENDIF
            reigvsd=reigvs
         ENDIF
      ENDDO  !do iout
      IF(.NOT.skip) iout=iout-1
      ioutbeg=iout+1
      nouttot=iout
!599   FORMAT("  Itr Nin    k-eff        Error k-eff      F.S error"      &
!           ,      "   Dom. R     PPM")
!600   FORMAT(79('-'))
!601      FORMAT(i7,i4,f12.5,l2,2(1p,e13.4,l2),0p,f8.4)
601      FORMAT(i7,i4,f12.5,l2,2(1p,e13.4,l2),0p,f8.4,i4)
      !ENDIF
      RETURN
      END SUBROUTINE outerss

      SUBROUTINE decusp

      IMPLICIT NONE

!      LOGICAL,SAVE :: not3anm=.TRUE.
!      REAL(8)      :: jnet,janm(n_group),phiavgr(n_group,3),phiavgrd(n_group,3)
!      REAL(8)      :: cofbdz(n_group)
!      REAL(8)      :: rl(n_group),rlp(n_group), rm1(n_group),rmp1(n_group)
!      REAL(8)      :: rm2(n_group),rmp2(n_group), df(n_group),dfp(n_group)
!      REAL(8)      :: fwf, tb, white, black, rhz, hzk, hzkb, fw, x
!      REAL(8)      :: rf, hzkt, hxy, difl, diflb, hzr, phiintd, phiint
!      REAL(8)      :: dc1, dc2, dnbd
      REAL(8), ALLOCATABLE :: cofbd(:)
!      REAL(8)      :: hxi, hxiw, diflw, albf, hyj
!      INTEGER(4)   :: nprnodes, iprod, la, kp, kpm1, kpp1, lp, l
!      INTEGER(4)   :: k, kf, m, ngray, ireg, ind, iz1, iz2, ig
!      INTEGER(4)   :: nfb, nft, i, j, kl, klp, ki, kip, iflag
!      INTEGER(4)   :: iz, it, isfc, is1, is2, ih, lw, le, ln, ls, ic
      REAL(8)      :: rt3, sqrt3, rsqrt3, hside

      INTEGER(4)   :: nprod, idecusp


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [decusp] in Mod_SStpen'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      nprod = 1  !! add control rod bank
      idecusp = 0

      sqrt3=SQRT(3.d0)
      rt3=1.73205080756888
      rsqrt3=1/sqrt3
      hside=gridsize_x(1)*rsqrt3

      !fwf(x)=x*(0.3364*x+0.1718)+0.4917
!! alxrf(1) is zero, alxrf(2) is zero
      IF(.NOT.ALLOCATED(cofbd))THEN
         ALLOCATE(cofbd(2))
         !IF(hex)THEN
         cofbd(1)=rt3*alxrf(1)*hside
         cofbd(2)=rt3*alxrf(2)*hside
         !ENDIF
      ENDIF
!
!      IF(fmfd.AND.idecusp.GT.0) THEN
!         WRITE(ioscn,*) 'not clearly checked and verified with fmfd'
!         WRITE(ioscn,*) 'idecusp=',idecusp
!      ENDIF
!
! ==== Need to modify for control rod operation
!!!!!      IF(nprod*idecusp /= 0 .AND. .NOT.scram) THEN
!!!!!         tb=dclock()
!!!!!         ntdecusp=ntdecusp+1
!!!!!         nfine=10
!!!!!         ntfine=3*nfine
!!!!!         ffine=nfine
!!!!!         nfreg(1)=nfine
!!!!!         nfreg(4)=nfine
!!!!!         nprnodes=0
!!!!!         singlenode = .TRUE.
!!!!!         DO iprod=1,nprod
!!!!!            la=ipprod(1,iprod)
!!!!!            kp=ipprod(2,iprod)
!!!!!            kpm1=kp-1
!!!!!            kpp1=kp+1
!!!!!            DO lp=laptr(la-1)+1,laptr(la)
!!!!!               l=latol(lp)
!!!!!               nprnodes=nprnodes+1
!!!!!               gray=zero
!!!!!               DO ic = 1,max_cr_cmp
!!!!!                  gray=gray+crbdens(la,kp,ic)
!!!!!               ENDDO
!!!!!               IF(qs1d .OR. idecusp.EQ.3) THEN
!!!!!                  fw=fwf(gray)
!!!!!                  fweight(1,l)=fw
!!!!!                  fweight(2,l)=fw
!!!!!                  CYCLE
!!!!!               ENDIF
!!!!!
!!!!!               fwgray(1)=gray*fweight(1,l)
!!!!!                  fweight(2,l)=fw
!!!!!                  CYCLE
!!!!!               ENDIF
!!!!!
!!!!!               fwgray(1)=gray*fweight(1,l)
!!!!!               fwgray(2)=gray*fweight(2,l)
!!!!!!
!!!!!! determine fine mesh sizes
!!!!!               ngray=NINT(gray*nfine)
!!!!!               ngray=MAX(1,ngray) !allow at least one mesh for unrodded region
!!!!!               ngray=MIN(nfine-1,ngray) !allow at least one mesh for rodded rg
!!!!!               IF (lwropt.EQ.0 .OR. lwropt.EQ.3) THEN                                          !bwrcr
!!!!!                  nfreg(2)=nfine-ngray
!!!!!                  nfreg(3)=ngray
!!!!!                  white=nfreg(2)
!!!!!                  black=nfreg(3)
!!!!!               ELSE                                                           !bwrcr
!!!!!                  nfreg(2)=ngray                                               !bwrcr
!!!!!                  nfreg(3)=nfine-ngray                                         !bwrcr
!!!!!                  white=nfreg(3)                                               !bwrcr
!!!!!                  black=nfreg(2)                                               !bwrcr
!!!!!               ENDIF                                                          !bwrcr
!!!!!               rhz=1/hz(kp)
!!!!!               hzf(1)=hz(kpm1)/ffine*rhz
!!!!!               IF (lwropt.EQ.0 .OR. lwropt.EQ.3) THEN                         !bwrcr
!!!!!                  hzf(2)=(1-gray)/white
!!!!!                  hzf(3)=gray/black
!!!!!               ELSE                                                           !bwrcr
!!!!!                  hzf(2)=gray/black                                            !bwrcr
!!!!!                  hzf(3)=(1-gray)/white                                        !bwrcr
!!!!!               ENDIF                                                          !bwrcr
!!!!!               hzf(4)=hz(kpp1)/ffine*rhz
!!!!!!
!!!!!!
!!!!!               hzf(1)=hzf(1)*hz(kp-1)
!!!!!               hzf(2)=hzf(2)*hz(kp)
!!!!!               hzf(3)=hzf(3)*hz(kp)
!!!!!               hzf(4)=hzf(4)*hz(kp+1)
!!!!!!
!!!!!! setup 1d linear system for the three node problem
!!!!!               CALL set3node(l,kp)
!!!!!!
!!!!!! solve the three node problem
!!!!!               CALL sol3node
!!!!!!
!!!!!! obtain average flux in the rodded node
!!!!!               k=nfine
!!!!!               DO ireg=2,3
!!!!!                  phiavgr(1,ireg)=0
!!!!!                  phiavgr(2,ireg)=0
!!!!!                  DO kf=1,nfreg(ireg)
!!!!!                     k=k+1
!!!!!                     phiavgr(1,ireg)=phiavgr(1,ireg)+phi1d(1,k)
!!!!!                     phiavgr(2,ireg)=phiavgr(2,ireg)+phi1d(2,k)
!!!!!                  ENDDO
!!!!!               ENDDO
!!!!!               phiavgr(1,1)=(phiavgr(1,2)*hzf(2)+phiavgr(1,3)*hzf(3))         &
!!!!!                    /hz(kp)
!!!!!               phiavgr(2,1)=(phiavgr(2,2)*hzf(2)+phiavgr(2,3)*hzf(3))         &
!!!!!                    /hz(kp)
!!!!!               phiavgrd(1,1)=phiavgr(1,1)
!!!!!               phiavgrd(2,1)=phiavgr(2,1)
!!!!!               IF (lwropt.EQ.0 .OR. lwropt.EQ.3) THEN                         !bwrcr
!!!!!                  rf=hzf(3)/hz(kp)
!!!!!                  fweight(1,l)=phiavgr(1,3)*rf/phiavgr(1,1)/gray
!!!!!                  fweight(2,l)=phiavgr(2,3)*rf/phiavgr(2,1)/gray
!!!!!               ELSE                                                           !bwrcr
!!!!!                  rf=hzf(2)/hz(kp)                                             !bwrcr
!!!!!                  fweight(1,l)=phiavgr(1,2)*rf/phiavgr(1,1)/gray               !bwrcr
!!!!!                  fweight(2,l)=phiavgr(2,2)*rf/phiavgr(2,1)/gray               !bwrcr
!!!!!               ENDIF                                                          !bwrcr
!!!!!               ksta=kp
!!!!!               kfin=kp
!!!!!               lsta=l
!!!!!               lfin=l
!!!!!               CALL xsecfb
!!!!!               psi(l,kp)=volnode(l,kp)*(xsnf(1,l,kp)*phi(1,l,kp)              &
!!!!!                    + xsnf(2,l,kp)*phi(2,l,kp))
!!!!!               IF(hex)THEN
!!!!!                  iz=kp
!!!!!                  DO it=1,6
!!!!!                     isfc=neigsfc(it,l)
!!!!!                     ind=neigsnd(3,isfc)
!!!!!                     IF(ind.EQ.12) THEN
!!!!!                        is1=neigsnd(1,isfc)
!!!!!                        is2=neigsnd(2,isfc)
!!!!!                        DO ig=1,ng
!!!!!                           dc1=xsd(ig,is1,iz)
!!!!!                           dc2=xsd(ig,is2,iz)
!!!!!                           dfd(ig,isfc,iz)=2.*dc1*dc2/(dc1+dc2)
!!!!!                        ENDDO
!!!!!                     ELSE
!!!!!                        is1=neigsnd(ind,isfc)
!!!!!                        DO ig=1,ng
!!!!!                           dc1=xsd(ig,is1,iz)
!!!!!                           dfd(ig,isfc,iz)=dc1*2.*cofbd(ig)/(cofbd(ig)+2.*dc1)
!!!!!                        ENDDO
!!!!!                     ENDIF
!!!!!                  ENDDO
!!!!!
!!!!!                  DO iz=kp,kpp1
!!!!!                     iz1=neigsndz(1,iz)
!!!!!                     iz2=neigsndz(2,iz)
!!!!!                     hzr=hz(iz2)/hz(iz1)
!!!!!                     ih=l
!!!!!                     DO ig=1,ng
!!!!!                        dc1=xsd(ig,ih,iz1)
!!!!!                        dc2=xsd(ig,ih,iz2)
!!!!!                        dfdz(ig,ih,iz)=dc2*(1.+hzr)/(hzr+dc2/dc1)
!!!!!                     ENDDO
!!!!!                  ENDDO
!!!!!               ENDIF
!!!!!!
!!!!!               IF(idecusp.EQ.1) CYCLE
!!!!!               DO k=1,ntfine
!!!!!                  phi1dd(1,k)=phi1d(1,k)
!!!!!                  phi1dd(2,k)=phi1d(2,k)
!!!!!               ENDDO
!!!!!!
!!!!!! solve the homegeneous problem
!!!!!               gray=0
!!!!!               fwgray(1)=0
!!!!!               fwgray(2)=0
!!!!!!
!!!!!! setup 1d linear system for the three node problem
!!!!!               CALL set3homo(l,kp)
!!!!!!
!!!!!! solve the three node problem
!!!!!               CALL sol3node
!!!!!!
!!!!!! obtain average flux in the rodded node
!!!!!               k=nfine
!!!!!               DO ireg=2,3
!!!!!                  phiavgr(1,ireg)=0
!!!!!                  phiavgr(2,ireg)=0
!!!!!                  DO kf=1,nfreg(ireg)
!!!!!                     k=k+1
!!!!!                     phiavgr(1,ireg)=phiavgr(1,ireg)+phi1d(1,k)
!!!!!                     phiavgr(2,ireg)=phiavgr(2,ireg)+phi1d(2,k)
!!!!!                  ENDDO
!!!!!               ENDDO
!!!!!               phiavgr(1,1)=(phiavgr(1,2)*hzf(2)+phiavgr(1,3)*hzf(3))         &
!!!!!                    /hz(kp)
!!!!!               phiavgr(2,1)=(phiavgr(2,2)*hzf(2)+phiavgr(2,3)*hzf(3))         &
!!!!!                    /hz(kp)
!!!!! !
!!!!! ! calculate adf's
!!!!!               nfb=nfine
!!!!!               nft=nfb+1
!!!!!               DO m=1,ng
!!!!!                  phiintd=0.5*(phi1dd(m,nft)+phi1dd(m,nfb))/phiavgrd(m,1)
!!!!!                  phiint=0.5*(phi1d(m,nft)+phi1d(m,nfb))/phiavgr(m,1)
!!!!!                  adfb(m,l,kp)=phiintd/phiint
!!!!!               ENDDO
!!!!!               nfb=2*nfine
!!!!!               nft=nfb+1
!!!!!               DO m=1,ng
!!!!!                  phiintd=0.5*(phi1dd(m,nft)+phi1dd(m,nfb))/phiavgrd(m,1)
!!!!!                  phiint=0.5*(phi1d(m,nft)+phi1d(m,nfb))/phiavgr(m,1)
!!!!!                  adft(m,l,kp)=phiintd/phiint
!!!!!               ENDDO
!!!!!               adfb(1,l,kpm1)=1
!!!!!               adfb(2,l,kpm1)=1
!!!!!               adft(1,l,kpm1)=1
!!!!!               adft(2,l,kpm1)=1
!!!!!               adfb(1,l,kpp1)=1
!!!!!               adfb(2,l,kpp1)=1
!!!!!               adft(1,l,kpp1)=1
!!!!!               adft(2,l,kpp1)=1
!!!!!            ENDDO
!!!!!         ENDDO
!!!!! !
!!!!!         ksta=1
!!!!!         kfin=nz
!!!!!         lsta=1
!!!!!         lfin=nxy
!!!!!         singlenode=.FALSE.
!!!!! !        tdecusp=tdecusp+dclock()-tb
!!!!!         IF (popt(2)) THEN                                                 !v21m14
!!!!!            WRITE(ioscn,*)'Rod Cusping Correction...',ntdecusp,nprnodes
!!!!!            WRITE(ioutp,*)'Rod Cusping Correction...',ntdecusp,nprnodes     !v21m14
!!!!!         ENDIF                                                             !v21m14
!!!!!      ENDIF

      RETURN
      END SUBROUTINE decusp

      SUBROUTINE IluFacHex
      IMPLICIT NONE
!
! perform incomplete LU factorization for the 3D coefficient matrices
!
      REAL(8) :: rinv(4),rmlt(4), rdet
      INTEGER(4) :: k, ir, l, ifr, ig, ir2, ipos
      INTEGER(4) :: ia1(4),ia2(4),ib1(4),ib2(4)
      DATA ia1/1,1,3,3/
      DATA ia2/2,2,4,4/
      DATA ib1/1,2,1,2/
      DATA ib2/3,4,3,4/
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: cftm
!
      !ALLOCATE(cftm(4, 7, nassy))

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [IluFacHex] in Mod_SStpen'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      !ALLOCATE(cftm(4, 7, nassy))
      ALLOCATE(cftm(N_group*N_group, 7, nassy))
      cftm=0d0
!
      DO k=1,nz
         DO l=1,nxy
            DO ir=1,7
               DO ig=1,4
                  cftm(ig,ir,l)=0
               ENDDO
            ENDDO
            DO ir=1,ipntr(0,l)
               cftm(1,ir,l)=-cmat(1,ir,l,k)
               cftm(4,ir,l)=-cmat(2,ir,l,k)
            ENDDO
            DO ig=1,4
               cftm(ig,7,l)=dcmat(ig,l,k)
            ENDDO
         ENDDO
!
         DO l=2,nxy
            DO ir=1,ilubnd(l)
               ifr=ipntr(ir,l)
               rdet=1/(cftm(1,7,ifr)*cftm(4,7,ifr)                 &
                    -cftm(2,7,ifr)*cftm(3,7,ifr))
               rinv(1)=cftm(4,7,ifr)*rdet
               rinv(2)=-cftm(2,7,ifr)*rdet
               rinv(3)=-cftm(3,7,ifr)*rdet
               rinv(4)=cftm(1,7,ifr)*rdet
               DO ig=1,4
                  rmlt(ig)=cftm(ia1(ig),ir,l)*rinv(ib1(ig))        &
                       +cftm(ia2(ig),ir,l)*rinv(ib2(ig))
               ENDDO
               DO ig=1,4
                  cftm(ig,ir,l)=rmlt(ig)
               ENDDO
               DO ir2=ir+1,ipntr(0,l)
                  ipos=iastopnt(ipntr(ir2,l),ifr)
                  IF(ipos.NE.0) THEN
                     DO ig=1,4
                        cftm(ig,ir2,l)=cftm(ig,ir2,l)              &
                             -rmlt(ia1(ig))*cftm(ib1(ig),ipos,ifr) &
                             -rmlt(ia2(ig))*cftm(ib2(ig),ipos,ifr)
                     ENDDO
                  ENDIF
               ENDDO
               ipos=iastopnt(l,ifr)
               DO ig=1,4
                  cftm(ig,7,l)=cftm(ig,7,l)                         &
                       -rmlt(ia1(ig))*cftm(ib1(ig),ipos,ifr)        &
                       -rmlt(ia2(ig))*cftm(ib2(ig),ipos,ifr)
               ENDDO
            ENDDO
         ENDDO
!
         DO l=1,nxy
            DO ir=1,6
               xlufac(1,ir,l,k)=cftm(1,ir,l)
               xlufac(2,ir,l,k)=cftm(4,ir,l)
            ENDDO
            rdet=1/(cftm(1,7,l)*cftm(4,7,l)                         &
                 -cftm(2,7,l)*cftm(3,7,l))
            delinv(1,l,k)=cftm(4,7,l)*rdet
            delinv(2,l,k)=-cftm(2,7,l)*rdet
            delinv(3,l,k)=-cftm(3,7,l)*rdet
            delinv(4,l,k)=cftm(1,7,l)*rdet
         ENDDO
      ENDDO
!
      DEALLOCATE(cftm)
!
      RETURN
      END SUBROUTINE IluFacHex

      END MODULE
#endif
