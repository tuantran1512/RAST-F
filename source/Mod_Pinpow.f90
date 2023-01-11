#ifdef siarhei_delete


      MODULE Mod_Pinpow

      USE Inc_FluxVar
      USE Inc_Control
      USE Inc_PinVar
      USE Inc_Solver
      USE Inc_Constant
      USE Inc_Geometry
      USE Inc_RP
      USE Inc_3D, ONLY: Flux
      USE Inc_DF
      USE Inc_maXS
      USE Inc_Option, ONLY: N_Group
      USE Inc_PinPOW
      use inc_branch, only: flag_coef


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE
      INTEGER,ALLOCATABLE,SAVE :: ipjptol(:,:)

      ! contains the routines for pin power reconstruction

      CONTAINS

      SUBROUTINE PinPOW

!      USE Mod_Soldhat, ONLY: updcurnxy

      ! calculate pin power distribution for the requested fuel assemblies
      IMPLICIT NONE
      INTEGER :: lfap, lfa

      ! generate corner point index and determime neigboring nodes and corners

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [PinPOW] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF ( Flag_First_pinpow .EQV. .FALSE. ) THEN
         Flag_First_pinpow = .TRUE.
         ALLOCATE(ipjptol(npin,npin))
         ipjptol = 0
         CALL initcorn
      END IF

      ! update the currents
      CALL updcurnxy
      CALL afenparam
      CALL setlscorn

      ! calculate corner point fluxes first
      CALL cornflx

      ! calculate detector response
      DO lfap = 1, Nxy_FA_1N
         lfa = pploc(lfap)
         CALL homoflx(lfa,npin,ipjptol)
         CALL ffmult(lfa,ipjptol)
      END DO

      if (.not.flag_coef) CALL calppeak

      RETURN
      END SUBROUTINE PinPOW


      subroutine pinpow_f ! after one step calculation, save pbu
      use mod_xsfb, only: hffsetfb, hfffb
      implicit none
      integer :: ixy_fa, ixy, ixy_1n, iz


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [pinpow_f] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.flag_pinpow) return

      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         ixy_1n=i_4nto1n(ixy)
         do iz=izfuelbot,izfueltop
            call hffsetfb(ixy,iz)
            call hfffb(ixy,ixy_1n,iz)
         enddo
      enddo
      flag_savepbu=.true.
      call pinpow
      flag_savepbu=.false.

      return
      end subroutine pinpow_f


      subroutine pinpow_s ! don't save pbu, usually for coupled calculation
!      use mod_soldhat, only: updcurnxy
      use mod_getsome, only: get_pow, get_linpow, get_avg
      use mod_save, only: xyz_indexing
      use mod_xsfb, only: hffsetfb, hfffb, xsfb
      use inc_3d, only: normal_power, power, avg_power
      use inc_th, only: core_power_100
      use inc_fa, only: h_fa
#ifdef js_r2mpi
      use inc_parallel, only: comm
#endif
      implicit none
      integer :: lfap, lfa, l
      integer :: ixy_fa, ixy, iz, ixy_1n
      integer :: ip, jp, ia, ja, k
      real(8) :: powsum, pownorm, tmpvol

#ifdef js_r2mpi

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [pinpow_s] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      if (.not.flag_first_pinpow) call xsfb
      call get_pow
      call get_linpow
      call get_avg(avg_power,power,0)
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         do iz=izfuelbot,izfueltop
            normal_power(ixy,iz)=power(ixy,iz)/max(1d-10,avg_power)
         enddo
      enddo
      call xyz_indexing
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         ixy_1n=i_4nto1n(ixy)
         do iz=izfuelbot,izfueltop
            call hffsetfb(ixy,iz)
            call hfffb(ixy,ixy_1n,iz)
         enddo
      enddo

      if (.not.flag_first_pinpow) then
         allocate(ipjptol(npin,npin))
         ipjptol=0
         call initcorn
      endif
      call updcurnxy
      call afenparam
      call setlscorn
      call cornflx
      do lfap=1,nxy_fa_1n
         lfa=pploc(lfap)
         call homoflx(lfa,npin,ipjptol)
         call ffmult(lfa,ipjptol)
      enddo

      if (.not.flag_first_pinpow) then
         flag_first_pinpow=.true.

         pownorm=core_power_100
         if (opt_core==4.and.nxy/=1) pownorm=pownorm/4d0
         powsum=0d0
         do lfap=1,nxy_fa_1n
            lfa=pploc(lfap)
            ixy_1n=i_fa_1n(lfa)
            l  = lfatol(lfaptr(lfa-1)+1)
            ia = Ix_4Nto1N(ltox(l))
            ja = Iy_4Nto1N(ltoy(l))
            do k=izfuelbot,izfueltop
               do ip=1,npin
                  do jp=1,npin
                     pinpow_3d(ia,ja,k,ip,jp)=hff(ixy_1n,k,ip,jp)
                     powsum=powsum+pinpow_3d(ia,ja,k,ip,jp)
                  enddo
               enddo
            enddo
         enddo
         pinpow_3d=pinpow_3d/max(1d-10,powsum)*pownorm
         do lfap=1,nxy_fa_1n
            lfa=pploc(lfap)
            l  = lfatol(lfaptr(lfa-1)+1)
            ia = Ix_4Nto1N(ltox(l))
            ja = Iy_4Nto1N(ltoy(l))
            do k=izfuelbot,izfueltop
               tmpvol=meshsize_z(k)*(h_fa/npin)**2d0
               do ip=1,npin
                  do jp=1,npin
                     pinpow_3d(ia,ja,k,ip,jp)=pinpow_3d(ia,ja,k,ip,jp)/max(1d-10,tmpvol)
                  enddo
               enddo
            enddo
         enddo
      endif

      call calppeak

      return
      end subroutine pinpow_s

      SUBROUTINE initcorn

      IMPLICIT NONE

      ! determine the corner point geometry paramters and initialize corner flux

      INTEGER :: ldiag(4),lneigh(4)
      INTEGER :: ic,j ,lc, ibeg, iend, i
      INTEGER :: ie, js, iw, jn, id, ic1, ic2, ic3, idp1
      INTEGER :: l, nneigh, k, in, lcp
      REAL(8) :: wf
      LOGICAL :: notifprev

      ! determine the corner point index

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [initcorn] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      j      = 1
      lc     = 0
      rsqrt2 = D1 / DSQRT(D2)

      ! -- at the north boundary
      IF ( .NOT. ( BC_Ly == 1 .AND. MeshSize_y(j) /= MeshSize_y(j+1) ) .OR. ( lc == 0 ) ) THEN
         ibeg = Ix_Start_y(j)
         iend = Ix_End_y(j) + 1
         IF ( ( BC_Lx == 1 ) .AND. ( MeshSize_x(ibeg) /= MeshSize_x(ibeg+1) ) ) THEN
            ibeg = ibeg + 1
         END IF
         nxcs(j) = ibeg
         nxce(j) = iend
         DO i = ibeg, iend
            lc = lc + 1
            nodec(i, j) = lc
            lctox(lc)   = i
            lctoy(lc)   = j
         END DO
      END IF

      ! -- interior rows
      DO j = 2, ny
         ibeg = Ix_Start_y(j)
         iend = Ix_End_y(j) + 1
         IF ( ( BC_Lx == 1 ) .AND. ( MeshSize_x(ibeg) /= MeshSize_x(ibeg+1) ) ) THEN
            ibeg = ibeg + 1
         END IF
         IF ( Ix_Start_y(j) > Ix_Start_y(j-1) ) THEN
            ibeg = Ix_Start_y(j-1)
         END IF
         IF ( Ix_End_y(j) < Ix_End_y(j-1) ) THEN
            iend = Ix_End_y(j-1) + 1
         END IF
         nxcs(j) = ibeg
         nxce(j) = iend
         DO i = ibeg, iend
            lc = lc + 1
            nodec(i, j) = lc
            lctox(lc)   = i
            lctoy(lc)   = j
         END DO
      END DO

      ! -- at the south boundary
      j    = ny + 1
      ibeg = Ix_Start_y(ny)
      iend = Ix_End_y(ny) + 1
      IF ( ( BC_Lx == 1 ) .AND. ( MeshSize_x(ibeg) /=  MeshSize_x(ibeg+1) ) ) THEN
         ibeg = ibeg + 1
      END IF
      nxcs(j) = ibeg
      nxce(j) = iend
      DO i = ibeg, iend
         lc = lc + 1

         nodec(i, j) = lc
         lctox(lc)   = i
         lctoy(lc)   = j
      END DO
      ncorn = lc

      ! process reflective boundary condition
      j = 1
      IF ( BC_Ly == 1 ) THEN
         js = j + 1
         IF ( MeshSize_y(j) /= MeshSize_y(js) ) THEN ! corner at ther center of the node
            DO i = nxcs(js), nxce(js)
               nodec(i, j) = nodec(i, js)
            END DO
         END IF
      END IF

      i = 1
      IF ( BC_Lx == 1 ) THEN
         ie = i + 1
         IF ( MeshSize_x(i) /= MeshSize_x(ie) ) THEN ! corner at ther center of the node
            DO j = 1, ny + 1
               nodec(i, j) = nodec(ie, j)
            END DO
         END IF
         j = 1
         IF ( BC_Ly == 1 ) THEN
            IF ( ( MeshSize_y(j) /= MeshSize_y(js) ) .AND. ( MeshSize_x(i) /= MeshSize_x(ie) ) ) THEN
               nodec(i, j) = nodec(ie, js)
            END IF
         END IF
      END IF

      ! determine the neighboring nodes and corner coupling
      DO lc = 1, ncorn
         i = lctox(lc)
         j = lctoy(lc)
         iw = i - 1
         ie = i + 1
         jn = j - 1
         js = j + 1
         ldiag(1) = nodec(iw, jn)
         ldiag(2) = nodec(ie, jn)
         ldiag(3) = nodec(ie, js)
         ldiag(4) = nodec(iw, js)
         lneigh(1) = nodec(iw, j )
         lneigh(2) = nodec(i , jn)
         lneigh(3) = nodec(ie, j )
         lneigh(4) = nodec(i , js)
         lcn(1, lc) = nodel(iw, jn)
         lcn(2, lc) = nodel(i , jn)
         lcn(3, lc) = nodel(i , j )
         lcn(4, lc) = nodel(iw, j )
         ic = 0
         notifprev = .TRUE.
         DO id = 1, 4
            IF ( ( lcn(id,lc) > 0 ) .AND. ( lcn(id,lc) <= nxy ) ) THEN
               IF ( ic /=0 ) THEN
                  ic = ic - 1
               END IF
               ic1 = ic + 1
               ic2 = ic + 2
               ic3 = ic + 3
               ic  = ic + 3
               notifprev = .TRUE.
            ELSE
               lcn(id, lc) = 0
               IF ( ic /= 0 .AND. notifprev ) THEN
                  ic = ic + 1
                  notifprev = .FALSE.
               END IF
               CYCLE
            END IF
            lcc(ic1, lc) = lneigh(id)
            lcc(ic2, lc) = ldiag(id)
            idp1 = id + 1
            IF ( id == 4 ) THEN
               idp1 = 1
               IF ( ( MOD(ic3, 2) /= 0 ) .AND. ( ic3 <= 8 ) ) THEN
                  lcc(ic3, lc) = lneigh(idp1)
               ELSE
                  ic3 = ic3 - 1
               END IF
            ELSE
               lcc(ic3, lc) = lneigh(idp1)
            END IF
         END DO
         nneighc(lc) = ic3
      END DO
      DO l = 1, nxy
         i = ltox(l)
         j = ltoy(l)
         lcnw(l) = nodec(i  , j  )
         lcsw(l) = nodec(i  , j+1)
         lcne(l) = nodec(i+1, j  )
         lcse(l) = nodec(i+1, j+1)
      END DO

      ! initialize corner fluxes
      ! -- at the boundaries
      j = 1
      DO lc = 1, ncorn
         i = lctox(lc)
         j = lctoy(lc)
         ldiag(1) = nodel(i-1,j-1)
         ldiag(2) = nodel(i,j-1)
         ldiag(3) = nodel(i,j)
         ldiag(4) = nodel(i-1,j)
         nneigh = 0
         DO id = 1, 4
            IF ( ldiag(id).GT.0 .AND. ldiag(id).LE.nxy) THEN
               nneigh=nneigh+1
               lneigh(nneigh)=ldiag(id)
            END IF
         END DO
         wf=nneigh
         wf=1/wf
         DO k=1,nz
            phicorn(1,lc,k)=0
            phicorn(2,lc,k)=0
            DO in=1,nneigh
               l=lneigh(in)
               phicorn(1,lc,k)=phicorn(1,lc,k)+flux(l,k,1)
               phicorn(2,lc,k)=phicorn(2,lc,k)+flux(l,k,2)
            END DO
            phicorn(1,lc,k)=phicorn(1,lc,k)*wf
            phicorn(2,lc,k)=phicorn(2,lc,k)*wf
         END DO
      END DO

      ! zero out at the external boundary
      IF( BC_Ly /= 1) THEN
         DO i=1,Nx
            j=Iy_Start_x(i)
            lc=nodec(i,j)
            lcp=nodec(i+1,j)
            IF(phicorn(1,lc,1).NE.0.) THEN
               DO k=1,nz
                  phicorn(1,lc,k)=0
                  phicorn(2,lc,k)=0
               END DO
            END IF
            IF(phicorn(1,lcp,1).NE.0.) THEN
               DO k=1,nz
                  phicorn(1,lcp,k)=0
                  phicorn(2,lcp,k)=0
               END DO
            END IF
         END DO
      END IF
      IF( BC_Ry /= 1 ) THEN
         DO i=1,Nx
            j=Iy_End_x(i)+1
            lc=nodec(i,j)
            lcp=nodec(i+1,j)
            IF(phicorn(1,lc,1).NE.0.) THEN
               DO k=1,nz
                  phicorn(1,lc,k)=0
                  phicorn(2,lc,k)=0
               END DO
            END IF
            IF(phicorn(1,lcp,1).NE.0.) THEN
               DO k=1,nz
                  phicorn(1,lcp,k)=0
                  phicorn(2,lcp,k)=0
               END DO
            END IF
         END DO
      END IF
      IF( BC_Lx /= 1) THEN
         DO j=1,ny
            i=Ix_Start_y(j)
            lc=nodec(i,j)
            lcp=nodec(i,j+1)
            IF(phicorn(1,lc,1).NE.0.) THEN
               DO k=1,nz
                  phicorn(1,lc,k)=0
                  phicorn(2,lc,k)=0
               END DO
            END IF
            IF(phicorn(1,lcp,1).NE.0.) THEN
               DO k=1,nz
                  phicorn(1,lcp,k)=0
                  phicorn(2,lcp,k)=0
               END DO
            END IF
         END DO
      END IF
      IF( BC_Rx /= 1 ) THEN
         DO j=1,ny
            i=Ix_End_y(j)+1
            lc=nodec(i,j)
            lcp=nodec(i,j+1)
            IF(phicorn(1,lc,1).NE.0.) THEN
               DO k=1,nz
                  phicorn(1,lc,k)=0
                  phicorn(2,lc,k)=0
               END DO
            END IF
            IF(phicorn(1,lcp,1).NE.0.) THEN
               DO k=1,nz
                  phicorn(1,lcp,k)=0
                  phicorn(2,lcp,k)=0
               END DO
            END IF
         END DO
      END IF

      RETURN
      END SUBROUTINE initcorn


      SUBROUTINE afenparam

      USE Inc_3D, ONLY: FisSrc
      USE Inc_Transient, ONLY: Flag_Transient

      IMPLICIT NONE

      ! calculate anm parameters which depends only on node properties

      REAL(8) :: kah,muh,kaht,muht,katemp
      REAL(8) :: xstp(N_Group),seff(N_Group)
      REAL(8) :: fdum, eps
      LOGICAL(1) :: ifcrik, ifcrim
      DATA ifcrik,ifcrim,fdum/.FALSE.,.FALSE.,1.0D0/
      DATA eps /1.D-4/
      INTEGER :: k, kp, l
      REAL(8) :: rkeff, hh, hht, rhz, rvol, akr, rxsd1, sr1, sr2
      REAL(8) :: sf1, b, c, bsqrt, rxss

      ! calculate afen parameters

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [afenparam] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      rkeff=reigv
      IF ( Flag_4N1FA ) THEN
         hh = Half*MeshSize_x(1)
      ELSE
         if (Nxy==1) then
            hh = MeshSize_x(1)
         else
            hh = Half*MeshSize_x(1)
         endif
      END IF

      hht=hh*(1.0D0/sqrt(2.0D0))
      DO k=1,nz
         kp=k+1
         rhz=1/MeshSize_z(k)
         DO l=1,nxy
            ! calculate effective source
            rvol=1.0D0/nodevolume(l,k)
            IF( Flag_Transient .EQV. .TRUE. ) THEN
               seff(1)=src(1,l,k)*rvol-rvdelt(1,l,k)*flux(l,k,1) &
                    -(1.0D0-betap(l,k))*FisSrc(l,k)*rvol         &
                    -(j_net_z_3d(1,l,kp)-j_net_z_3d(1,l,k))*rhz  &
                    -vr(l,k,1)*rvol
               seff(2)=src(2,l,k)*rvol-rvdelt(2,l,k)*flux(l,k,2) &
                    -(j_net_z_3d(2,l,kp)-j_net_z_3d(2,l,k))*rhz  &
                    -vr(l,k,2)*rvol
            ELSE
               seff(1)=-(j_net_z_3d(1,l,kp)-j_net_z_3d(1,l,k))*rhz
               seff(2)=-(j_net_z_3d(2,l,kp)-j_net_z_3d(2,l,k))*rhz
            END IF
            ! correct total cross section for effective source
            xstp(1)=maxs_r_3d(l,k,1)-seff(1)/flux(l,k,1)
            xstp(2)=maxs_r_3d(l,k,2)-seff(2)/flux(l,k,2)

            akr=(nu_maxs_f_3d(l,k,1)+nu_maxs_f_3d(l,k,2)*maxs_s_3d(l,k,1)/max(1d-10,xstp(2)))/max(1d-10,xstp(1))*rkeff
            akratio(l,k)=akr
            rxsd1=1/max(1d-10,d_3d(l,k,1))
            sr1=xstp(1)*rxsd1
            sr2=xstp(2)/max(1d-10,d_3d(l,k,2))
            sf1=nu_maxs_f_3d(l,k,1)*rxsd1
            b=0.5*(sr1+sr2-rkeff*sf1)
            bsqrt=SQRT(max(1d-10,b*b-(1-akr)*sr1*sr2))  !the number is radical should be non-negative if keff*xsnf(1)*xsnf(2)/xsd(1)/xsd(2)>=0
            c=b-bsqrt
            IF(c<0)THEN
               kflag(l,k)=.TRUE.
               akappa(l,k)=SQRT(-c)
            ELSE
               kflag(l,k)=.false.
               akappa(l,k)=SQRT(c)
            END IF

            if (akappa(l,k)<1d-30) akappa(l,k)=1d-20

            amu(l,k)=DSQRT(dabs(bsqrt+b))  ! a flag for  bsqrt+b<0  should be added in a better fix
            rxss=1/max(1d-10,maxs_s_3d(l,k,1))
            ar(l,k)=(d_3d(l,k,2)*(bsqrt-b)+xstp(2))*rxss
            as(l,k)=(-d_3d(l,k,2)*(bsqrt+b)+xstp(2))*rxss
            ardet(l,k)=1/max(1d-10,ar(l,k)-as(l,k))

            IF((akappa(l,k) < epsanm)) THEN
               ifcrik=.true.
               bctka(l,k)=epsanm
            ELSE
               ifcrik=.false.
               bctka(l,k)=akappa(l,k)
            END IF
            IF(ABS(amu(l,k)) < epsanm) THEN
               ifcrim=.true.
               bctmu(l,k)=epsanm
            ELSE
               ifcrim=.false.
               bctmu(l,k)=amu(l,k)
            END IF
            ! -- calculate sn and cn functions
            kah=akappa(l,k)*hh
            muh=amu(l,k)*hh
            katemp=akappa(l,k)
            CALL fsncnp(kflag(l,k),ifcrik,hh,katemp,kah,bctka(l,k) &
                 , asnka(l,k),asnkatl(l,k),acnka(l,k),acnkatl(l,k))
            IF(ifcrim) THEN
               asnmu(l,k)=hh+(muh**2)*hh/6+(muh**4)*hh/120
               asnmutl(l,k)=hh+((bctmu(l,k)*hh)**2)*hh/6 + &
                    ((bctmu(l,k)*hh)**4)*hh/120
            ELSE
               asnmu(l,k)=SINH(muh)/amu(l,k)
               asnmutl(l,k)=SINH(bctmu(l,k)*hh)/bctmu(l,k)
            END IF
            acnmu(l,k)=SQRT(1+(amu(l,k)*asnmu(l,k))**2)
            acnmutl(l,k)=SQRT(1+(bctmu(l,k)*asnmutl(l,k))**2)

            ! calculate cross terms
            kaht=akappa(l,k)*hht
            muht=amu(l,k)*hht
            CALL fsncnp(kflag(l,k),ifcrik,hh,katemp,kaht,bctka(l,k), &
                 asnkat(l,k),asnkattl(l,k),acnkat(l,k),acnkattl(l,k))
            katemp=hht*bctmu(l,k)
            IF(ifcrim) THEN
               asnmut(l,k)=hh+(muht**2)*hh/6+(muht**4)*hh/120
               asnmuttl(l,k)=hh+(katemp**2)*hh/6+(katemp**4)*hh/120
            ELSE
               asnmut(l,k)=SINH(muht)*hh/muht
               asnmuttl(l,k)=SINH(katemp)*hh/katemp
            END IF
            acnmut(l,k)=SQRT(1+(muht*asnmut(l,k)/hh)**2)
            acnmuttl(l,k)=SQRT(1+(asnmuttl(l,k)*katemp/hh)**2)
         END DO
      END DO
      RETURN
      END SUBROUTINE afenparam


      SUBROUTINE setlscorn

      IMPLICIT NONE

      ! setup the linear system for corner flux calculation
      INTEGER :: i, j, k, l, m, ic, ic1, ic2, ic3
      INTEGER :: in, idc0, idc1, idc3
      INTEGER :: la, lc, lp, lfa, ipr, lpbeg, lpend, nneigh, nodesfa
      REAL(8) :: alb, anum, alba, cdfmt, curka
      REAL(8) :: curmu, dsign
      REAL(8) :: r, s, d_1, d_2, hh, rd1, rd2, rdet, rdetrs
      DIMENSION curka(2),curmu(2)
      DIMENSION cdfmt(N_Group),idc0(4),idc1(4),idc3(4)
      DATA idc0/3,4,1,2/
      DATA idc1/4,1,2,3/
      DATA idc3/2,3,4,1/
      LOGICAL notifprev


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [setlscorn] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF ( Flag_4N1FA ) THEN
         hh = Half*MeshSize_x(1)
      ELSE
         hh = MeshSize_x(1)
      END IF

      DO k = IzFuelBot, IzFuelTop
         DO lc=1,ncorn
            DO ic=0,nneighc(lc)
               cc11(ic,lc,k)=0
               cc12(ic,lc,k)=0
               cc21(ic,lc,k)=0
               cc22(ic,lc,k)=0
            END DO
         END DO
         DO l=1,nxy
            dsign=1
            IF(kflag(l,k)) dsign=-1
            CALL cornercoef(hh,dsign,akappa(l,k),bctka(l,k)            &
                 , asnka(l,k),asnkatl(l,k),acnka(l,k),acnkatl(l,k)     &
                 , asnkat(l,k),asnkattl(l,k),acnkat(l,k),acnkattl(l,k) &
                 , curka,cornfka(:,l,k))
            dsign=1
            CALL cornercoef(hh,dsign,amu(l,k),bctmu(l,k)                &
                 , asnmu(l,k),asnmutl(l,k),acnmu(l,k),acnmutl(l,k)      &
                 , asnmut(l,k),asnmuttl(l,k),acnmut(l,k),acnmuttl(l,k)  &
                 , curmu,cornfmu(:,l,k))

            r=ar(l,k)
            s=as(l,k)
            rdet=ardet(l,k)
            rd1=rdet/max(1d-10,d_3d(l,k,1))
            rd2=rdet/max(1d-10,d_3d(l,k,2))
            d_1=d_3d(l,k,1)
            d_2=d_3d(l,k,2)
            DO in=1,2
               anum=r*curka(in)-s*curmu(in)
               cur11(in,l,k)=rdet*anum
               anum=curka(in)-curmu(in)
               cur12(in,l,k)=-r*s*d_1*rd2*anum
               cur21(in,l,k)=d_2*rd1*anum
               anum=s*curka(in)-r*curmu(in)
               cur22(in,l,k)=-rdet*anum
            END DO
         END DO
      END DO

      ! assign corner pointer discontinutity factors
      CDF_LxLy = D1

      DO k=izfuelbot,izfueltop
         DO l=1,nxy
            DO m=1,N_Group
               cdf(1,m,l,k) = D1 / CDF_RxRy(l, k, m)  ! Copy edge CDF first
               cdf(2,m,l,k) = D1 / CDF_RxRy(l, k, m)  ! Copy edge CDF first
               cdf(3,m,l,k) = D1 / CDF_RxRy(l, k, m)  ! Copy edge CDF first
               cdf(4,m,l,k) = D1 / CDF_RxRy(l, k, m)  ! Copy edge CDF first
            END DO
         END DO
      END DO

      ! assign cdf at the middle of each assembly
      IF(npfa.EQ.4) THEN
         DO k=izfuelbot,izfueltop
            ipr=k            ! not to use izpr for only coarse mesh
            DO lfa = 1, Nxy_FA_1N
               l=lfatol(lfaptr(lfa))
               la=ltola(l)
               lpbeg=lfaptr(lfa-1)+1
               lpend=lfaptr(lfa)
               nodesfa=lpend-lpbeg+1

               ! Soo, make adfmt 1 / FA symmetric line CDF

               DO m = 1,N_Group
                  cdfmt(m) = D1 / CDF_RxLy(l, k, m)  ! RxLy or LxRy. Same.
               END DO

               ! Soo, note iquad(l)
               ! ---------
               ! | 1 | 2 |
               ! ---------
               ! | 4 | 3 |
               ! ---------

               ! Soo, note cdf()
               !  1     2
               !   -----
               !   |   |
               !   -----
               !  4     3

               DO lp = lpbeg, lpend
                  l = lfatol(lp)

                  IF(iquad(l).EQ.1) THEN
                     DO m=1,N_Group
                        ! cdf(3,m,l,k) = D1
                        ! cdf(2,m,l,k) = D1 / CDF_RxRy(l, k, m)
                        ! cdf(4,m,l,k) = D1 / CDF_LxLy(l, k, m)
                        cdf(3,m,l,k) = D1 ! Center edge
                        cdf(2,m,l,k) = cdfmt(m)
                        cdf(4,m,l,k) = cdfmt(m)
                     END DO
                  ELSE IF(iquad(l).EQ.2) THEN
                     DO m=1,N_Group
                        ! cdf(4,m,l,k) = D1
                        ! cdf(1,m,l,k) = D1 / CDF_RxLy(l, k, m)
                        ! cdf(3,m,l,k) = D1 / CDF_LxRy(l, k, m)
                        cdf(4,m,l,k) = D1 ! Center edge
                        cdf(1,m,l,k) = cdfmt(m)
                        cdf(3,m,l,k) = cdfmt(m)
                     END DO
                  ELSE IF (iquad(l).EQ.3) THEN
                     DO m=1,N_Group
                        ! cdf(1,m,l,k) = D1
                        ! cdf(2,m,l,k) = D1 / CDF_RxRy(l, k, m)
                        ! cdf(4,m,l,k) = D1 / CDF_LxLy(l, k, m)
                        cdf(1,m,l,k) = D1 ! Center edge
                        cdf(2,m,l,k) = cdfmt(m)
                        cdf(4,m,l,k) = cdfmt(m)
                     END DO
                  ELSE
                     DO m=1,N_Group
                        ! cdf(2,m,l,k) = D1
                        ! cdf(1,m,l,k) = D1 / CDF_RxLy(l, k, m)
                        ! cdf(3,m,l,k) = D1 / CDF_LxRy(l, k, m)
                        cdf(2,m,l,k) = 1d0 ! Center edge
                        cdf(1,m,l,k) = cdfmt(m)
                        cdf(3,m,l,k) = cdfmt(m)
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END IF
      alb=0
      IF( BC_Rx == 3 .OR. BC_Ry == 3 ) alb=cinf
      IF( BC_Rx == 2 .OR. BC_Ry == 2 ) alb=0.5

      DO lc=1,ncorn
         nneigh=nneighc(lc)
         ic=0
         notifprev=.true.
         DO in=1,4
            l=lcn(in,lc)
            IF(l.NE.0) THEN
               IF(ic.NE.0) ic=ic-1
               ic1=ic+1
               ic2=ic+2
               ic3=ic+3
               IF(in.EQ.4) THEN
                  IF(MOD(ic3,2).EQ.0 .OR. ic3.GT.8) ic3=1
               END IF
               ic=ic+3
               notifprev=.true.
            ELSE
               IF(ic.NE.0 .AND. notifprev) THEN
                  ic=ic+1
                  notifprev=.false.
               END IF
               CYCLE
            END IF
            DO k=izfuelbot,izfueltop
               r=ar(l,k)
               s=as(l,k)
               d_1=d_3d(l,k,1)
               d_2=d_3d(l,k,2)
               rdet=ardet(l,k)
               rdetrs=rdet*r*s

               anum=r*cornfka(0,l,k)-s*cornfmu(0,l,k)
               cc11(0,lc,k)=cc11(0,lc,k)+rdet*anum*d_1                   &
                    *cdf(idc0(in),1,l,k)
               anum=cornfka(0,l,k)-cornfmu(0,l,k)
               cc12(0,lc,k)=cc12(0,lc,k)-rdetrs*anum*d_1                 &
                    *cdf(idc0(in),2,l,k)
               cc21(0,lc,k)=cc21(0,lc,k)+rdet*anum*d_2                   &
                    *cdf(idc0(in),1,l,k)
               anum=s*cornfka(0,l,k)-r*cornfmu(0,l,k)
               cc22(0,lc,k)=cc22(0,lc,k)-rdet*anum*d_2                   &
                    *cdf(idc0(in),2,l,k)

               anum=r*cornfka(2,l,k)-s*cornfmu(2,l,k)
               cc11(ic2,lc,k)=cc11(ic2,lc,k)+rdet*anum*d_1               &
                    *cdf(in,1,l,k)
               anum=cornfka(2,l,k)-cornfmu(2,l,k)
               cc12(ic2,lc,k)=cc12(ic2,lc,k)-rdetrs*anum*d_1             &
                    *cdf(in,2,l,k)
               cc21(ic2,lc,k)=cc21(ic2,lc,k)+rdet*anum*d_2               &
                    *cdf(in,1,l,k)
               anum=s*cornfka(2,l,k)-r*cornfmu(2,l,k)
               cc22(ic2,lc,k)=cc22(ic2,lc,k)-rdet*anum*d_2               &
                    *cdf(in,2,l,k)

               anum=r*cornfka(1,l,k)-s*cornfmu(1,l,k)
               cc11(ic1,lc,k)=cc11(ic1,lc,k)+rdet*anum*d_1               &
                    *cdf(idc1(in),1,l,k)
               cc11(ic3,lc,k)=cc11(ic3,lc,k)+rdet*anum*d_1               &
                    *cdf(idc3(in),1,l,k)
               anum=cornfka(1,l,k)-cornfmu(1,l,k)
               cc12(ic1,lc,k)=cc12(ic1,lc,k)-rdetrs*anum*d_1             &
                    *cdf(idc1(in),2,l,k)
               cc12(ic3,lc,k)=cc12(ic3,lc,k)-rdetrs*anum*d_1             &
                    *cdf(idc3(in),2,l,k)
               cc21(ic1,lc,k)=cc21(ic1,lc,k)+rdet*anum*d_2               &
                    *cdf(idc1(in),1,l,k)
               cc21(ic3,lc,k)=cc21(ic3,lc,k)+rdet*anum*d_2               &
                    *cdf(idc3(in),1,l,k)
               anum=s*cornfka(1,l,k)-r*cornfmu(1,l,k)
               cc22(ic1,lc,k)=cc22(ic1,lc,k)-rdet*anum*d_2               &
                    *cdf(idc1(in),2,l,k)
               cc22(ic3,lc,k)=cc22(ic3,lc,k)-rdet*anum*d_2               &
                    *cdf(idc3(in),2,l,k)
            END DO
         END DO

         ! incorporate boundary condition
         IF(nneigh.NE.8) THEN
            i=lctox(lc)
            j=lctoy(lc)
            IF((i.EQ.1 .AND. BC_Lx == 1 ) .OR.                         &
                 (j.EQ.1 .AND. BC_Ly == 1 )) THEN
               alba=0
               IF(i.EQ.Nx+1 .OR. j.EQ.ny+1) alba=alb
            ELSE
               alba=alb*2
            END IF
            DO k=izfuelbot,izfueltop
               cc11(0,lc,k)=cc11(0,lc,k)+alba
               cc22(0,lc,k)=cc22(0,lc,k)+alba
            END DO
         END IF

         ! store the inverse of the 2x2 diagonal block element
         DO k=izfuelbot,izfueltop
            rdet=1/max(1d-10,cc11(0,lc,k)*cc22(0,lc,k)-cc12(0,lc,k)*cc21(0,lc,k))
            anum=cc11(0,lc,k)
            cc11(0,lc,k)=rdet*cc22(0,lc,k)
            cc12(0,lc,k)=-rdet*cc12(0,lc,k)
            cc21(0,lc,k)=-rdet*cc21(0,lc,k)
            cc22(0,lc,k)=rdet*anum
         END DO
      END DO

      RETURN
      END SUBROUTINE setlscorn


      SUBROUTINE cornercoef(hh,dsign,ka,bct,S,stl,C,ctl,St,sttl,Ct,cttl,curf,cornf)

      IMPLICIT NONE

      ! calculate coefficients of the corner point balance equation

      REAL(8) :: S,C,hh,ka,St,Ct,kat,kah,dsign, bct
      REAL(8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t12,t13,Xt,C2t,S2t
      REAL(8) :: z1,z2,z3,z4,z5,rsqrt2, c2ttl, cttl, s2ttl, sttl
      REAL(8) :: xttls, xttlc, beta, gama, ctl, stl
      REAL(8) :: curf(2),cornf(0:2)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [cornercoef] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      rsqrt2 = D1 / DSQRT(D2)

      kah=ka*hh
      kat=rsqrt2*ka
      C2ttl=Cttl*Ct
      s2ttl=sttl*st
      c2t=ct*ct
      S2t=St*St
      Xt=St*Ct
      xttls=ct*sttl
      xttlc=st*cttl
      beta=hh*xttlc/xttls
      gama=(ka/bct)**2
      t1=1/(S-C*beta)
      t2=1/(hh*c2t*s-ctl*s2t*gama)/D2
      t3=1/Xttls
      t4=C2ttl+dsign*bct*bct*s2ttl/D2
      t5=ct/2/st
      t6=0.5*C*t1
      t7=0.5*S*kah*dsign*ka*Xt*t2
      t8=0.5*stl*ka*ka*dsign*S2t*t2
      t9=0.5*C*hh*t1*t3*t4
      t10=t7-t8
      t12=t10+t5
      t13=t6-t9
      cornf(0)=t12+t13
      cornf(1)=t10-t5
      cornf(2)=t12-t13
      z1=0.5*C*beta*t1
      z2=hh*stl*C2t*t2
      z3=Ctl*kah*Xt*ka*t2/bct/bct
      z4=hh*s*t1*t3*t4/2
      z5=z4-z1
      curf(1)=-z2+z3-z5
      curf(2)=z2-z3-z5

      RETURN
      END SUBROUTINE cornercoef

      SUBROUTINE cornflx

      IMPLICIT NONE

      ! determine the corner point fluxes based on node average flux distributions
      INTEGER :: j, l, le, ls, k, i

      ! incorporate reflective bc on current

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [cornflx] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF( BC_Lx == 1 .AND. MeshSize_x(1).NE.MeshSize_x(2)) THEN
         DO j=1,ny
            l=nodel(Ix_Start_y(j),j)
            le=nodel(Ix_Start_y(j)+1,j)
            DO k=1,nz
               j_net_x_3d(1,l,k)=-j_net_x_3d(1,le,k)
               j_net_x_3d(2,l,k)=-j_net_x_3d(2,le,k)
            END DO
         END DO
      END IF
      IF( BC_Ly == 1 .AND. MeshSize_y(1).NE.MeshSize_y(2)) THEN
         j=1
         DO i=1,Nx
            l=nodel(i,j)
            ls=nodel(i,j+1)
            DO k=1,nz
               j_net_y_3d(1,l,k)=-j_net_y_3d(1,ls,k)
               j_net_y_3d(2,l,k)=-j_net_y_3d(2,ls,k)
            END DO
         END DO
      END IF
      ! loop over planes
      DO k=izfuelbot,izfueltop
         CALL cpbrhs(k)
         CALL solcorn( k, phicorn(:,1:ncorn,k), cc11(:,1:ncorn,k), cc12(:,1:ncorn,k), cc21(:,1:ncorn,k), cc22(:,1:ncorn,k) )
      END DO

      RETURN
      END SUBROUTINE cornflx


      SUBROUTINE cpbrhs(k)

      IMPLICIT NONE

      ! determine the rhs of the CPB linear system containing sources associated
      REAL(8) :: curw(N_Group),cure(N_Group),curn(N_Group),curs(N_Group)
      REAL(8) :: curnear(N_Group),curfar(N_Group)
      LOGICAL :: notifn,notifw,ifsymx,ifsymy
      INTEGER :: k, lc, jp1, le, ls, lc1, lc2, lc3
      INTEGER :: lc4, m, i, j, l


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [cpbrhs] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      DO lc=1,ncorn
         cpbsrc(1,lc)=0
         cpbsrc(2,lc)=0
      END DO
      ifsymx=( BC_Lx == 1 ) .AND. (MeshSize_x(1).NE.MeshSize_x(2))
      ifsymy=( BC_Ly == 1 ) .AND. (MeshSize_y(1).NE.MeshSize_y(2))
      DO j=1,ny
         jp1=j+1
         notifn=.true.
         IF(ifsymy .AND. j.EQ.1) notifn=.false.
         DO i=Ix_Start_y(j),Ix_End_y(j)
            l=nodel(i,j)
            le=nodel(i+1,j)
            ls=nodel(i,jp1)
            lc1=lcnw(l)
            lc2=lcne(l)
            lc3=lcse(l)
            lc4=lcsw(l)
            notifw=.true.
            IF(ifsymx .AND. i.EQ.Ix_Start_y(j)) notifw=.false.
            DO m=1,N_Group
               curw(m)=j_net_x_3d(m,l,k)
               cure(m)=j_net_x_3d(m,le,k)
               curn(m)=j_net_y_3d(m,l,k)
               curs(m)=j_net_y_3d(m,ls,k)
            END DO
            ! contribution to the south-east corner
            DO m=1,N_Group
               curnear(m)=-cure(m)-curs(m)
               curfar(m)=-curw(m)-curn(m)
            END DO
            CALL fillcpb(cpbsrc(:,lc3),l,k,curnear,curfar)
            ! contribution to the north-east corner
            IF(notifn) THEN
               DO m=1,N_Group
                  curnear(m)=-cure(m)+curn(m)
                  curfar(m)=-curw(m)+curs(m)
               END DO
               CALL fillcpb(cpbsrc(:,lc2),l,k,curnear,curfar)
            END IF
            ! contribution to the north-west corner
            IF(notifn .AND. notifw) THEN
               DO m=1,N_Group
                  curnear(m)=curw(m)+curn(m)
                  curfar(m)=cure(m)+curs(m)
               END DO
               CALL fillcpb(cpbsrc(:,lc1),l,k,curnear,curfar)
            END IF
            ! contribution to the south-west corner
            IF(notifw) THEN
               DO m=1,N_Group
                  curnear(m)=curw(m)-curs(m)
                  curfar(m)=cure(m)-curn(m)
               END DO
               CALL fillcpb(cpbsrc(:,lc4),l,k,curnear,curfar)
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE cpbrhs


      SUBROUTINE fillcpb(avec,l,k,curnear,curfar)

      IMPLICIT NONE

      REAL(8) :: avec(N_Group),curnear(N_Group),curfar(N_Group)
      INTEGER :: l, k

      avec(1)=avec(1)                                                 &
           +cur11(1,l,k)*curnear(1)+cur11(2,l,k)*curfar(1)            &
           +cur12(1,l,k)*curnear(2)+cur12(2,l,k)*curfar(2)
      avec(2)=avec(2)                                                 &
           +cur21(1,l,k)*curnear(1)+cur21(2,l,k)*curfar(1)            &
           +cur22(1,l,k)*curnear(2)+cur22(2,l,k)*curfar(2)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fillcpb] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      RETURN
      END SUBROUTINE fillcpb


      SUBROUTINE solcorn(k,phic,a11,a12,a21,a22)


      IMPLICIT NONE

      ! perform iteration to determine the corner fluxes at a plane
      LOGICAL :: notconv
      INTEGER :: i,j,k,l,lc,ln,id,in,itr,nneigh,ldiag(4),lneigh(4)
      REAL(8) :: eps,err1,err2,rhs1,rhs2,sumd,phicd1,phicd2
      REAL(8),ALLOCATABLE,SAVE :: phicmss(:,:)
      REAL(8) :: wf,sumn
      REAL(8) :: phic(N_Group,1:ncorn),a11(0:8,1:ncorn),a12(0:8,1:ncorn),a21(0:8,1:ncorn),a22(0:8,1:ncorn)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [solcorn] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF(.NOT.ALLOCATED(phicmss)) THEN
         ALLOCATE(phicmss(N_Group,ncorn))
         phicmss=0
      END IF
      notconv=.true.
      itr=0
      ! initial guess
      DO lc=1,ncorn
         i=lctox(lc)
         j=lctoy(lc)
         ldiag(1)=nodel(i-1,j-1)
         ldiag(2)=nodel(i,j-1)
         ldiag(3)=nodel(i,j)
         ldiag(4)=nodel(i-1,j)
         nneigh=0
         DO id=1,4
            IF(ldiag(id).GT.0 .AND. ldiag(id).LE.nxy) THEN
               nneigh=nneigh+1
               lneigh(nneigh)=ldiag(id)
            END IF
         END DO
         wf=nneigh
         wf=1/wf
         phicmss(1,lc)=0
         phicmss(2,lc)=0
         DO in=1,nneigh
            l=lneigh(in)
            phicmss(1,lc)=phicmss(1,lc)+flux(l,k,1)
            phicmss(2,lc)=phicmss(2,lc)+flux(l,k,2)
         END DO
         phicmss(1,lc)=phicmss(1,lc)*wf
         phicmss(2,lc)=phicmss(2,lc)*wf
      END DO
      DO lc=1,ncorn
         phic(1,lc)=phicmss(1,lc)
         phic(2,lc)=phicmss(2,lc)
      END DO

      DO WHILE(notconv)
         sumd=0
         sumn=0
         DO lc=1,ncorn
            rhs1=cpbsrc(1,lc)
            rhs2=cpbsrc(2,lc)
            DO in=1,nneighc(lc)
               ln=lcc(in,lc)
               rhs1=rhs1-(a11(in,lc)*phic(1,ln)+a12(in,lc)*phic(2,ln))
               rhs2=rhs2-(a21(in,lc)*phic(1,ln)+a22(in,lc)*phic(2,ln))
            END DO
            phicd1=phic(1,lc)
            phicd2=phic(2,lc)
            phic(1,lc)=rhs1*a11(0,lc)+rhs2*a12(0,lc)
            phic(2,lc)=rhs1*a21(0,lc)+rhs2*a22(0,lc)
            err1=phic(1,lc)-phicd1
            err2=phic(2,lc)-phicd2
            sumn=sumn+err1*err1+err2*err2
            sumd=sumd+phic(1,lc)*phicd1+phic(2,lc)*phicd2
         END DO
         errl2=SQRT(ABS(sumn/max(1d-10,sumd)))
         itr=itr+1
         IF(errl2.LT.EPS_Global*0.01) notconv=.false.
         IF(itr.GT.Iout_Max) RETURN
      END DO

      eps=1/cinf
      DO lc=1,ncorn
         IF(ABS(phic(1,lc)).LT.eps .OR. phic(1,lc).LT.0) phic(1,lc)=0
         IF(ABS(phic(2,lc)).LT.eps .OR. phic(2,lc).LT.0) phic(2,lc)=0
      END DO

      RETURN
      END SUBROUTINE solcorn


      SUBROUTINE homoflx(idfa,npin,ipjptol)

      ! calculate the homogeneouS intra-nodal flux distribution
      IMPLICIT NONE

      REAL(8) :: psic(ng2,4) !1-nw,2-ne,3-se,4-sw
      REAL(8) :: psidi(ng2,4) !1-xl,2-xr,3-yl,4-yr
      INTEGER :: lc(4)
      INTEGER :: lfa, idfa, ibegl, iendl, neut, iendr, istepl
      INTEGER :: ibegr, istepr, inode, nnodesfa, lp, l
      INTEGER :: i, j, le, ls, k, m, jbeg, jend, jstep
      INTEGER :: ibeg, iend, istep, id, lcid, ic, ir, jr, npin, ipin
      INTEGER :: ipjptol(npin,npin) !dmm
      REAL(8) :: hh, fnpin, hpin, rdet, r, s, rd1, rd2
      REAL(8) :: psic1, psic2, psic3, psic4, psixs, psixd, psiys, psiyd
      REAL(8) :: s2toka, s2toka8, hhk, sumf1, sumf2, anum
      REAL(8) :: fdum
      REAL(8) :: fka(ng2),rka(ng2),dsign(ng2)
      REAL(8) :: sn(ng2+2),cn(ng2+2),snt(ng2+2),cnt(ng2+2)
      REAL(8) :: s2t(ng2),c2t(ng2),xt(ng2),hhka(ng2)
      REAL(8) :: snx(ng2+2),cnx(ng2+2),stx(ng2+2),ctx(ng2+2)
      REAL(8) :: snx0(ng2+2),cnx0(ng2+2),stx0(ng2+2),ctx0(ng2+2)
      REAL(8) :: sny(ng2+2),cny(ng2+2),sty(ng2+2),cty(ng2+2)
      REAL(8) :: cnh(ng2+2),snh(ng2+2),snhx(ng2+2),snhy(ng2+2)
      REAL(8) :: cnht(ng2+2),snht(ng2+2),snhxt(ng2+2),snhyt(ng2+2)
      REAL(8) :: betax(ng2+2),betaxt(ng2+2),betay(ng2+2),betayt(ng2+2)
      REAL(8) :: snxa(ng2+2),cnxa(ng2+2),stxa(ng2+2),ctxa(ng2+2)
      REAL(8) :: snya(ng2+2),cnya(ng2+2),stya(ng2+2),ctya(ng2+2)
      REAL(8) :: w(8,ng2+4),gama(ng2),beta(ng2),bct2(ng2)
      REAL(8) :: a(8,ng2),fxyggg(ng2)
      REAL(8) :: phibar(ng2),psibar(ng2) !node-average
      LOGICAL(1) :: ifoddpin, kflagl(ng2), ifcri(ng2)

      ! determine modal fluxes and currents

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [homoflx] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF ( Flag_4N1FA ) THEN
         hh = Half*MeshSize_x(1)
      ELSE
         hh = MeshSize_x(1)
      END IF

      lfa=idfa
      fnpin=npin
      ibegl=1
      fdum=1.0   ! taken from the definition
      IF(npfa.EQ.1) THEN
         iendl=npin
      ELSE
         neut=2
         IF(MOD(npin,2).EQ.1) THEN
            iendl=(npin+1)/neut
            iendr=iendl
         ELSE
            iendl=npin/neut
            iendr=iendl+1
         END IF
      END IF
      istepl=1
      ibegr=npin
      istepr=-1
      IF(MOD(npin,2).EQ.0) THEN
         ifoddpin=.false.
      ELSE
         ifoddpin=.true.
      END IF
      inode=0
      nnodesfa=lfaptr(lfa)-lfaptr(lfa-1)
      DO lp=lfaptr(lfa-1)+1,lfaptr(lfa)
         inode=inode+1
         l=lfatol(lp)
         i=ltox(l)
         j=ltoy(l)
         lc(1)=lcnw(l)
         lc(2)=lcne(l)
         lc(3)=lcse(l)
         lc(4)=lcsw(l)
         le=nodel(i+1,j)
         ls=nodel(i,j+1)
         DO k=izfuelbot,izfueltop
            hpin=pinpitch(l,k)
            rdet=ardet(l,k)
            r=ar(l,k)
            s=as(l,k)
            sn(1)=asnka(l,k)
            sn(2)=asnmu(l,k)
            sn(3)=asnkatl(l,k)
            sn(4)=asnmutl(l,k)
            cn(1)=acnka(l,k)
            cn(2)=acnmu(l,k)
            cn(3)=acnkatl(l,k)
            cn(4)=acnmutl(l,k)
            snt(1)=asnkat(l,k)
            snt(2)=asnmut(l,k)
            snt(3)=asnkattl(l,k)
            snt(4)=asnmuttl(l,k)
            cnt(1)=acnkat(l,k)
            cnt(2)=acnmut(l,k)
            cnt(3)=acnkattl(l,k)
            cnt(4)=acnmuttl(l,k)
            hhka(1)=hh*akappa(l,k)
            hhka(2)=hh*amu(l,k)
            fka(1)=akappa(l,k)
            fka(2)=amu(l,k)
            kflagl(2)=.false.
            IF(kflag(l,k)) THEN
               dsign(1)=-1
               kflagl(1)=.true.
            ELSE
               dsign(1)=1
               kflagl(1)=.false.
            END IF
            dsign(2)=1
            IF((fka(1).LT.epsanm)) THEN
               ifcri(1)=.true.
            ELSE
               ifcri(1)=.false.
            END IF
            IF(fka(2).LT.epsanm) THEN
               ifcri(2)=.true.
            ELSE
               ifcri(2)=.false.
            END IF
            bct2(1)=bctka(l,k)
            bct2(2)=bctmu(l,k)
            DO m=1,N_Group
               rka(m)=1/max(1d-10,fka(m))
               s2t(m)=snt(m)*snt(m)
               c2t(m)=cnt(m)*cnt(m)
               xt(m)=snt(m)*cnt(m)
               beta(m)=hh*snt(m)*cnt(m+2)/max(1d-10,cnt(m))/max(1d-10,snt(m+2))
               gama(m)=(fka(m)/max(1d-10,bct2(m)))**2
            END DO
            DO id=1,4
               lcid=lc(id)
               psic(1,id)=rdet*(phicorn(1,lcid,k)*cdf(id,1,l,k)         &
                    -s*phicorn(2,lcid,k)*cdf(id,2,l,k))
               psic(2,id)=rdet*(r*phicorn(2,lcid,k)*cdf(id,2,l,k)       &
                    -phicorn(1,lcid,k)*cdf(id,1,l,k))
            END DO
            rd1=-1/max(1d-10,d_3d(l,k,1))
            rd2=-1/max(1d-10,d_3d(l,k,2))
            psidi(1,1)=rdet*(j_net_x_3d(1,l,k)*rd1-s*j_net_x_3d(2,l,k)*rd2)
            psidi(2,1)=rdet*(r*j_net_x_3d(2,l,k)*rd2-j_net_x_3d(1,l,k)*rd1)
            psidi(1,2)=rdet*(j_net_x_3d(1,le,k)*rd1-s*j_net_x_3d(2,le,k)*rd2)
            psidi(2,2)=rdet*(r*j_net_x_3d(2,le,k)*rd2-j_net_x_3d(1,le,k)*rd1)
            psidi(1,3)=rdet*(j_net_y_3d(1,l,k)*rd1-s*j_net_y_3d(2,l,k)*rd2)
            psidi(2,3)=rdet*(r*j_net_y_3d(2,l,k)*rd2-j_net_y_3d(1,l,k)*rd1)
            psidi(1,4)=rdet*(j_net_y_3d(1,ls,k)*rd1-s*j_net_y_3d(2,ls,k)*rd2)
            psidi(2,4)=rdet*(r*j_net_y_3d(2,ls,k)*rd2-j_net_y_3d(1,ls,k)*rd1)
            DO m=1,N_Group
               psic1=0.25*(psic(m,3)-psic(m,4)-psic(m,1)+psic(m,2))
               psic2=0.25*(psic(m,3)-psic(m,4)+psic(m,1)-psic(m,2))
               psic3=0.25*(psic(m,3)+psic(m,4)-psic(m,1)-psic(m,2))
               psic4=0.25*(psic(m,3)+psic(m,4)+psic(m,1)+psic(m,2))
               psixs=0.5*(psidi(m,2)+psidi(m,1))
               psixd=0.5*(psidi(m,2)-psidi(m,1))*dsign(m)
               psiys=0.5*(psidi(m,4)+psidi(m,3))
               psiyd=0.5*(psidi(m,4)-psidi(m,3))*dsign(m)
               a(6,m)=psic2/max(1d-10,s2t(m))
               if (abs(sn(m)-beta(m)*cn(m))<1d-10) then
                  a(1,m)=(psic1-beta(m)*psixs)/max(1d-10,(sn(m)-beta(m)*cn(m)))
               else
                  a(1,m)=(psic1-beta(m)*psixs)/(sn(m)-beta(m)*cn(m))
               endif
               if (abs(sn(m)-beta(m)*cn(m))<1d-10) then
                  a(3,m)=(psic3-beta(m)*psiys)/max(1d-10,(sn(m)-beta(m)*cn(m)))
               else
                  a(3,m)=(psic3-beta(m)*psiys)/(sn(m)-beta(m)*cn(m))
               endif
               a(5,m)=(psic1-sn(m)*a(1,m))/max(1d-10,cnt(m+2))/max(1d-10,snt(m))
               a(7,m)=(psic3-sn(m)*a(3,m))/max(1d-10,cnt(m+2))/max(1d-10,snt(m))
               s2toka=s2t(m)*gama(m)
               a(8,m)=(2*hh*sn(m)*psic4-2*hh*cn(m+2)/max(1d-10,(bct2(m)**2))* &
                    (psixd+psiyd))/max(1d-10,(hh*2*sn(m)*c2t(m)-2*cn(m+2)*s2toka))
               s2toka8=s2t(m)*a(8,m)*gama(m)
               a(2,m)=(2*hh*psixd/max(1d-10,(bct2(m)**2))-s2toka8)/2/max(1d-10,sn(m))/hh
               a(4,m)=(2*hh*psiyd/max(1d-10,(bct2(m)**2))-s2toka8)/2/max(1d-10,sn(m))/hh
               psibar(m)=(a(2,m)+a(4,m))*sn(m)/max(1d-10,hhka(m))                  &
                    +2*a(8,m)*s2t(m)/max(1d-10,(hhka(m)*hhka(m)))
            END DO
            phibar(1)=r*psibar(1)+s*psibar(2)
            phibar(2)=psibar(1)+psibar(2)

            IF (abs(xoffset(l,k)-0d0)<1d-10) THEN
               DO m=1,N_Group+2
                  snxa(m)=sn(m)
                  cnxa(m)=cn(m)
                  stxa(m)=snt(m)
                  ctxa(m)=cnt(m)
               END DO
            ELSE
               ! -- calculate sn0 and cn0 functions
               DO m=1,N_Group
                  hhk=fka(m)*(hh-xoffset(l,k))
                  CALL fsncnp(kflagl(m),ifcri(m),(hh-xoffset(l,k))    &
                       , fka(m),hhk,bct2(m),snxa(m),snxa(m+2)      &
                       , cnxa(m),cnxa(m+2))
                  hhk=hhk*rsqrt2
                  CALL fsncnp(kflagl(m),ifcri(m),(hh-xoffset(l,k))  &
                       , fka(m),hhk,bct2(m),stxa(m),stxa(m+2)    &
                       , ctxa(m),ctxa(m+2))
               END DO
            END IF
            IF (abs(yoffset(l,k)-0d0)<1d-10) THEN
               DO m=1,N_Group+2
                  snya(m)=sn(m)
                  cnya(m)=cn(m)
                  stya(m)=snt(m)
                  ctya(m)=cnt(m)
               END DO
            ELSE
               ! -- calculate sn0 and cn0 functions
               DO m=1,N_Group
                  hhk=fka(m)*(hh-yoffset(l,k))
                  CALL fsncnp(kflagl(m),ifcri(m),(hh-yoffset(l,k))   &
                       , fka(m),hhk,bct2(m),snya(m),snya(m+2)   &
                       , cnya(m),cnya(m+2))
                  hhk=hhk*rsqrt2
                  CALL fsncnp(kflagl(m),ifcri(m),(hh-yoffset(l,k))  &
                       , fka(m),hhk,bct2(m),stya(m),stya(m+2)   &
                       , ctya(m),ctya(m+2))
               END DO
            END IF

            IF(npfa.EQ.1) THEN
               ibeg=ibegl
               iend=iendl
               istep=istepl
               jbeg=ibegl
               jend=iendl
               jstep=istepl
               DO m=1,N_Group+2
                  snx0(m)=-snxa(m)
                  cnx0(m)=cnxa(m)
                  stx0(m)=-stxa(m)
                  ctx0(m)=ctxa(m)
                  sny(m)=-snya(m)
                  cny(m)=cnya(m)
                  sty(m)=-stya(m)
                  cty(m)=ctya(m)
               END DO
            ELSEIF(npfa.EQ.4) THEN
               IF(iquad(l).EQ.1 .OR. iquad(l).EQ. 4) THEN
                  ibeg=ibegl
                  iend=iendl
                  istep=istepl
                  DO m=1,N_Group+2
                     snx0(m)=-snxa(m)
                     cnx0(m)=cnxa(m)
                     stx0(m)=-stxa(m)
                     ctx0(m)=ctxa(m)
                  END DO
               ELSE
                  ibeg=ibegr
                  iend=iendr
                  istep=istepr
                  DO m=1,N_Group+2
                     snx0(m)=snxa(m)
                     cnx0(m)=cnxa(m)
                     stx0(m)=stxa(m)
                     ctx0(m)=ctxa(m)
                  END DO
               END IF
               IF(iquad(l).EQ.1 .OR. iquad(l).EQ.2) THEN
                  jbeg=ibegl
                  jend=iendl
                  jstep=istepl
                  DO m=1,N_Group+2
                     sny(m)=-snya(m)
                     cny(m)=cnya(m)
                     sty(m)=-stya(m)
                     cty(m)=ctya(m)
                  END DO
               ELSE
                  jbeg=ibegr
                  jend=iendr
                  jstep=istepr
                  DO m=1,N_Group+2
                     sny(m)=snya(m)
                     cny(m)=cnya(m)
                     sty(m)=stya(m)
                     cty(m)=ctya(m)
                  END DO
               END IF
            END IF

            DO m=1,N_Group
               hhk=fka(m)*hpin
               rka(m)=fka(m)*fka(m)*dsign(m)/2
               CALL fsncnp(kflagl(m),ifcri(m),hpin, fka(m),hhk,bct2(m),snh(m),snh(m+2), cnh(m),cnh(m+2))
               betax(m)=(cnh(m)-1)*dsign(m)
               betax(m+2)=(cnh(m+2)-1)*dsign(m)
               betay(m)=betax(m)
               betay(m+2)=betax(m+2)
               snhx(m)=snh(m)
               snhx(m+2)=snh(m+2)
               snhy(m)=snh(m)
               snhy(m+2)=snh(m+2)
               hhk=hhk*rsqrt2
               CALL fsncnp(kflagl(m),ifcri(m),hpin, fka(m),hhk,bct2(m),snht(m),snht(m+2), cnht(m),cnht(m+2))
               betaxt(m)=(cnht(m)-1)*dsign(m)
               betaxt(m+2)=(cnht(m+2)-1)*dsign(m)
               betayt(m)=betaxt(m)
               betayt(m+2)=betaxt(m+2)
               snhxt(m)=snht(m)
               snhxt(m+2)=snht(m+2)
               snhyt(m)=snht(m)
               snhyt(m+2)=snht(m+2)
               IF(ibeg.EQ.npin) THEN
                  snhx(m)=-snhx(m)
                  snhx(m+2)=-snhx(m+2)
                  snhxt(m)=-snhxt(m)
                  snhxt(m+2)=-snhxt(m+2)
                  betax(m)=-betax(m)
                  betax(m+2)=-betax(m+2)
                  betaxt(m)=-betaxt(m)
                  betaxt(m+2)=-betaxt(m+2)
               END IF
               IF(jbeg.EQ.npin) THEN
                  snhy(m)=-snhy(m)
                  snhy(m+2)=-snhy(m+2)
                  snhyt(m)=-snhyt(m)
                  snhyt(m+2)=-snhyt(m+2)
                  betay(m)=-betay(m)
                  betay(m+2)=-betay(m+2)
                  betayt(m)=-betayt(m)
                  betayt(m+2)=-betayt(m+2)
               END IF
               w(1,m)=betax(m)*a(1,m)/max(1d-10,(fka(m)**2))
               w(1,m+2)=snh(m+2)*a(2,m)
               w(2,m+2)=betax(m+2)*a(2,m)*dsign(m)
               w(2,m)=snh(m)*a(1,m)
               w(3,m)=betay(m)*a(3,m)/max(1d-10,(fka(m)**2))
               w(3,m+2)=snh(m+2)*a(4,m)
               w(4,m+2)=betay(m+2)*a(4,m)*dsign(m)
               w(4,m)=snh(m)*a(3,m)
               w(5,m)=snht(m)*snht(m+2)*a(5,m)
               w(5,m+2)=snht(m)*betayt(m)*a(6,m)*2/max(1d-10,(fka(m)**2))+snht(m)*a(8,m)*dsign(m)*betaxt(m)
               if (abs(rka(m))<1d-10) then
                  w(5,m+4)=betaxt(m+2)*betayt(m)*a(7,m)/max(1d-10,rka(m))
               else
                  w(5,m+4)=betaxt(m+2)*betayt(m)*a(7,m)/rka(m)
               endif
               w(6,m)=snht(m)*betayt(m+2)*a(5,m)*dsign(m)
               w(6,m+2)=snht(m)*snht(m)*a(6,m)+betaxt(m)*betayt(m)*a(8,m)
               w(6,m+4)=betaxt(m+2)*snht(m)*a(7,m)*dsign(m)
               if (abs(rka(m))<1d-10) then
                  w(7,m)=betayt(m+2)*betaxt(m)*a(5,m)/max(1d-10,rka(m))
               else
                  w(7,m)=betayt(m+2)*betaxt(m)*a(5,m)/rka(m)
               endif
               if (abs(rka(m))<1d-10) then
                  w(7,m+2)=snht(m)*betaxt(m)*a(6,m)/max(1d-10,rka(m))*dsign(m)+dsign(m)*snht(m)*betayt(m)*a(8,m)
               else
                  w(7,m+2)=snht(m)*betaxt(m)*a(6,m)/rka(m)*dsign(m)+dsign(m)*snht(m)*betayt(m)*a(8,m)
               endif
               w(7,m+4)=snht(m+2)*snht(m)*a(7,m)
               if (abs(rka(m))<1d-10) then
                  w(8,m)=betaxt(m)*snht(m)*a(5,m)/max(1d-10,rka(m))*dsign(m)
               else
                  w(8,m)=betaxt(m)*snht(m)*a(5,m)/rka(m)*dsign(m)
               endif
               w(8,m+2)=betaxt(m)*betayt(m)*a(6,m)/max(1d-10,(rka(m)**2))+a(8,m)*snht(m)*snht(m)
               if (abs(rka(m))<1d-10) then
                  w(8,m+4)=snht(m+2)*a(7,m)*betayt(m)/max(1d-10,rka(m))*dsign(m)
               else
                  w(8,m+4)=snht(m+2)*a(7,m)*betayt(m)/rka(m)*dsign(m)
               endif
               DO ic=1,4
                  w(ic,m)=w(ic,m)/hpin
                  w(ic,m+2)=w(ic,m+2)/hpin
                  w(ic+4,m)=w(ic+4,m)/hpin/hpin
                  w(ic+4,m+2)=w(ic+4,m+2)/hpin/hpin
                  w(ic+4,m+4)=w(ic+4,m+4)/hpin/hpin
               END DO
            END DO
            sumf1=0
            sumf2=0
            DO j=jbeg,jend,jstep
               DO m=1,N_Group+2
                  snx(m)=snx0(m)
                  cnx(m)=cnx0(m)
                  stx(m)=stx0(m)
                  ctx(m)=ctx0(m)
               END DO
               DO i=ibeg,iend,istep
                  DO m=1,N_Group
                     fxyggg(m)=cnx(m)*w(1,m)+cnx(m+2)*w(1,m+2)+    &
                          snx(m)*w(2,m)+snx(m+2)*w(2,m+2)+    &
                          cny(m)*w(3,m)+cny(m+2)*w(3,m+2)+    &
                          sny(m)*w(4,m)+sny(m+2)*w(4,m+2)+    &
                          stx(m)*cty(m+2)*w(5,m)+            &
                          stx(m)*cty(m)*w(5,m+2)+            &
                          stx(m+2)*cty(m)*w(5,m+4)+         &
                          stx(m)*sty(m+2)*w(6,m)+           &
                          stx(m)*sty(m)*w(6,m+2)+           &
                          stx(m+2)*sty(m)*w(6,m+4)+         &
                          ctx(m)*sty(m+2)*w(7,m)+           &
                          ctx(m)*sty(m)*w(7,m+2)+           &
                          ctx(m+2)*sty(m)*w(7,m+4)+         &
                          ctx(m)*cty(m+2)*w(8,m)+           &
                          ctx(m)*cty(m)*w(8,m+2)+           &
                          ctx(m+2)*cty(m)*w(8,m+4)
                  END DO
                  phihom(i,j,lfa,k,1)=r*fxyggg(1)+s*fxyggg(2)
                  phihom(i,j,lfa,k,2)=fxyggg(1)+fxyggg(2)
                  ipjptol(i,j)=l
                  sumf1=sumf1+phihom(i,j,lfa,k,1)
                  sumf2=sumf2+phihom(i,j,lfa,k,2)
                  ! update sn and cn functions in the x-direction
                  DO m=1,N_Group
                     anum=snx(m)*dsign(m)*fka(m)*fka(m)
                     snx(m)=snx(m)*cnh(m)+cnx(m)*snhx(m)
                     cnx(m)=cnx(m)*cnh(m)+anum*snhx(m)
                     anum=snx(m+2)*dsign(m)*bct2(m)*bct2(m)
                     snx(m+2)=snx(m+2)*cnh(m+2)+cnx(m+2)*snhx(m+2)
                     cnx(m+2)=cnx(m+2)*cnh(m+2)+anum*snhx(m+2)
                     anum=stx(m)*dsign(m)*fka(m)*fka(m)/2
                     stx(m)=stx(m)*cnht(m)+ctx(m)*snhxt(m)
                     ctx(m)=ctx(m)*cnht(m)+anum*snhxt(m)
                     anum=stx(m+2)*dsign(m)*bct2(m)*bct2(m)/2
                     stx(m+2)=stx(m+2)*cnht(m+2)+ctx(m+2)*snhxt(m+2)
                     ctx(m+2)=ctx(m+2)*cnht(m+2)+anum*snhxt(m+2)
                  END DO
               END DO
               ! update sn and cn functions in the y-direction
               DO m=1,N_Group
                  anum=sny(m)*dsign(m)*fka(m)*fka(m)
                  sny(m)=sny(m)*cnh(m)+cny(m)*snhy(m)
                  cny(m)=cny(m)*cnh(m)+anum*snhy(m)
                  anum=sny(m+2)*dsign(m)*bct2(m)*bct2(m)
                  sny(m+2)=sny(m+2)*cnh(m+2)+cny(m+2)*snhy(m+2)
                  cny(m+2)=cny(m+2)*cnh(m+2)+anum*snhy(m+2)
                  anum=sty(m)*dsign(m)*fka(m)*fka(m)/2
                  sty(m)=sty(m)*cnht(m)+cty(m)*snhyt(m)
                  cty(m)=cty(m)*cnht(m)+anum*snhyt(m)
                  anum=sty(m+2)*dsign(m)*bct2(m)*bct2(m)/2
                  sty(m+2)=sty(m+2)*cnht(m+2)+cty(m+2)*snhyt(m+2)
                  cty(m+2)=cty(m+2)*cnht(m+2)+anum*snhyt(m+2)
               END DO
            END DO
         END DO  ! k = IzFuelBot, IzFuelTop
      END DO
      ! expand homoflux
      IF(npfa.EQ.4 .AND. nnodesfa.NE.4) THEN
         ipin=mod(npin,2)
         if(ipin==1)then
             IF(nnodesfa.EQ.1) THEN !at the center assembly
                DO k=izfuelbot,izfueltop
                   DO j=ibegr,iendr,istepr
                      DO i=ibegl,iendl-1
                         ir=npin-i+1
                         phihom(i,j,lfa,k,1)=phihom(ir,j,lfa,k,1)
                         phihom(i,j,lfa,k,2)=phihom(ir,j,lfa,k,2)
                         ipjptol(i,j)=ipjptol(ir,j)
                      END DO
                   END DO
                   DO j=ibegl,iendl-1
                      jr=npin-j+1
                      DO i=1,npin
                         phihom(i,j,lfa,k,1)=phihom(i,jr,lfa,k,1)
                         phihom(i,j,lfa,k,2)=phihom(i,jr,lfa,k,2)
                         ipjptol(i,j)=ipjptol(i,jr)
                      END DO
                   END DO
                END DO
             ELSEIF(I_Lx_P2R(l).EQ.l) THEN !on y axis
                DO k=izfuelbot,izfueltop
                   DO j=1,npin
                      DO i=ibegl,iendl-1
                         ir=npin-i+1
                         phihom(i,j,lfa,k,1)=phihom(ir,j,lfa,k,1)
                         phihom(i,j,lfa,k,2)=phihom(ir,j,lfa,k,2)
                         ipjptol(i,j)=ipjptol(ir,j)
                      END DO
                   END DO
                END DO
             ELSEIF(I_Ly_P2R(l).EQ.l) THEN !on x axis
                DO k=izfuelbot,izfueltop
                   DO j=ibegl,iendl-1
                      jr=npin-j+1
                      DO i=1,npin
                         phihom(i,j,lfa,k,1)=phihom(i,jr,lfa,k,1)
                         phihom(i,j,lfa,k,2)=phihom(i,jr,lfa,k,2)
                         ipjptol(i,j)=ipjptol(i,jr)
                      END DO
                   END DO
                END DO
             END IF
         else  ! for even number of pins
             IF(nnodesfa.EQ.1) THEN !at the center assembly
                DO k=izfuelbot,izfueltop
                   DO j=ibegr,iendr,istepr
                      DO i=ibegl,iendl
                         ir=npin-i+1
                         phihom(i,j,lfa,k,1)=phihom(ir,j,lfa,k,1)
                         phihom(i,j,lfa,k,2)=phihom(ir,j,lfa,k,2)
                         ipjptol(i,j)=ipjptol(ir,j)
                      END DO
                   END DO
                   DO j=ibegl,iendl
                      jr=npin-j+1
                      DO i=1,npin
                         phihom(i,j,lfa,k,1)=phihom(i,jr,lfa,k,1)
                         phihom(i,j,lfa,k,2)=phihom(i,jr,lfa,k,2)
                         ipjptol(i,j)=ipjptol(i,jr)
                      END DO
                   END DO
                END DO
             ELSEIF(I_Lx_P2R(l).EQ.l) THEN !on y axis
                DO k=izfuelbot,izfueltop
                   DO j=1,npin
                      DO i=ibegl,iendl
                         ir=npin-i+1
                         phihom(i,j,lfa,k,1)=phihom(ir,j,lfa,k,1)
                         phihom(i,j,lfa,k,2)=phihom(ir,j,lfa,k,2)
                         ipjptol(i,j)=ipjptol(ir,j)
                      END DO
                   END DO
                END DO
             ELSEIF(I_Ly_P2R(l).EQ.l) THEN !on x axis
                DO k=izfuelbot,izfueltop
                   DO j=ibegl,iendl
                      jr=npin-j+1
                      DO i=1,npin
                         phihom(i,j,lfa,k,1)=phihom(i,jr,lfa,k,1)
                         phihom(i,j,lfa,k,2)=phihom(i,jr,lfa,k,2)
                         ipjptol(i,j)=ipjptol(i,jr)
                      END DO
                   END DO
                END DO
             END IF
         END IF
      END IF

      RETURN
      END SUBROUTINE homoflx


      SUBROUTINE ffmult(lfa,ipjptol)

      USE Inc_FA
      USE Inc_3D, ONLY: Avg_Power
      USE Inc_TH, ONLY: PPower
      USE Inc_Time
      USE Inc_XYZ
      USE Inc_PinPOW
#ifdef JR_SRCTRM
      use inc_SNF, only: flag_snf_pin, phihom_M, SNF_k
#endif

      IMPLICIT NONE

      ! multiplies homogeneous flux by form functions to obtain pinwise powers
      ! control rod insertion is taken care of by using weighted form functions
      INTEGER :: ipjptol(npin,npin)
      INTEGER :: l, ia, ja, nnodes, jp, ip, k, Ig, la, lfa
      INTEGER :: lptr, ipmax, jpmax, icomp
      INTEGER :: Ixy_1N, I_LP
      REAL(8) :: powtot(2), avgpinp(2), avgasyp(2), Buff_Sum(2)
      REAL(8) :: fnodes, hzk
      REAL(8) :: fnorm, fnorm1, fnorm2
      REAL(8) :: rnorm, pmax
#ifdef JR_SRCTRM
      REAL(8) :: temp1
      REAL(8) :: temp_sum_1, temp_sum_2, temp_mid_1, temp_mid_2
      REAL(8) :: temp_pin_sum_1, temp_pin_sum_2, temp_pin_mid_1, temp_pin_mid_2
      INTEGER :: pin_num, m
#endif

      ! determine fuel assembly coordinate and number of nodes in the fa

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [ffmult] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      l  = lfatol(lfaptr(lfa-1)+1)
      ia = Ix_4Nto1N(ltox(l))
      ja = Iy_4Nto1N(ltoy(l))

      nnodes = lfaptr(lfa)-lfaptr(lfa-1)
      fnodes = real(nnodes,8)

      Ixy_1N = I_FA_1N( lfa    )
      I_LP   = I_LP_1N( Ixy_1N )

      if (dT<1e-30) then
         if (sum(pinbu_3d(:,:,:,:,:))<1d-10) then
            PinBU_2D(ia,ja,:,:) = BU_XY(ia,ja)
            do k=izfuelbot,izfueltop
               PinBU_3D(ia,ja,k,:,:) = BU_XYZ(ia,ja,k)
            enddo
         endif
      endif
      ! reset the axially integrated radial pin power
      powvalr = D0
      DO k = izfuelbot, izfueltop
         powtot = D0
         hzk = MeshSize_z(k) / hactive
         l   = lfatol( lfaptr( lfa - 1 ) + 1 )
         la  = ltola(l)
         icomp = I_Comp(l, k)
         DO jp = 1, npin
            DO ip = 1, npin
               l = ipjptol(ip, jp)
               DO Ig = 1, N_Group
                  if (phihom(ip,jp,lfa,k,ig)<1d-10) phihom(ip,jp,lfa,k,ig)=1d-10
                  powval(ip,jp,k,Ig) = phihom(ip,jp,lfa,k,Ig) * kap_maxs_f_3d(l,k,Ig) * HFF(Ixy_1N,k,ip,jp)
                  powtot(Ig) = powtot(Ig) + powval(ip,jp,k,Ig)
               END DO
            END DO
         END DO
         DO Ig = 1, N_Group
            avgpinp(Ig) = powtot(ig) / N_Pin
         END DO

         Buff_Sum = D0
         DO lptr = lfaptr(lfa-1) + 1, lfaptr(lfa)
            l = lfatol(lptr)
            DO Ig = 1, N_Group
               Buff_Sum(Ig) = Buff_Sum(Ig) + flux(l,k,Ig)*kap_maxs_f_3d(l,k,Ig)
            END DO
         END DO

         DO Ig = 1, N_Group
            avgasyp(Ig) = PPower * Buff_Sum(Ig) / max(1d-10,fnodes)  !pp@1
         END DO
         fnorm = SUM(avgasyp)
         rnorm = D1 / max(1d-10,fnorm)
         fnorm1 = avgasyp(1) / max(1d-10,avgpinp(1))
         fnorm2 = avgasyp(2) / max(1d-10,avgpinp(2))
         Buff_Sum = D0
         pmax     = D0
         DO jp = 1, npin
            DO ip = 1, npin
               powvalr(ip, jp, k) = SUM( powval(ip,jp,k,:) )

               Buff_Sum(1) = Buff_Sum(1) + (powvalr(ip, jp, k))

               IF ( powvalr(ip, jp, k) > pmax ) THEN
                  pmax  = powvalr(ip, jp, k)
                  ipmax = ip
                  jpmax = jp
               END IF
            END DO
         END DO

         ! Store 3-D Peak Values
         pppeak(lfa, k) = pmax
         Peak_X(lfa, k) = ip
         Peak_Y(lfa, k) = jp

         ! IntraFA Normalization in a k Plane
         fnorm = Buff_Sum(1) / (real(npin,8)*real(npin,8)) !js+ N_Pin !wrong way but best now..
         rnorm = D1 / max(1d-10,fnorm)

         DO jp = 1, npin
            DO ip = 1, npin
               PinPOW_3D(ia,ja,k,ip,jp) = Normal_Power_XYZ(ia,ja,k)*Avg_Power*powvalr(ip,jp,k)*rnorm
               ! unit [W/cm3]           = [-]                      * [W/cm3] * [-]            * [-]
               if (pinpow_3d(ia,ja,k,ip,jp)<0d0) pinpow_3d(ia,ja,k,ip,jp)=1d-10
               PinQ_3D(ia,ja,k,ip,jp) = PinPOW_3D(ia,ja,k,ip,jp)*(h_FA/real(npin,8))**D2*MeshSize_z(k)
               ! unit [W]             = [W/cm3]                 * [cm2]         * [cm]
               if (flag_savepbu) then
                  PinBU_3D(ia,ja,k,ip,jp)=PinBU_3D(ia,ja,k,ip,jp) &
                                        &+dT*PinQ_3D(ia,ja,k,ip,jp)/max(1d-10,PinMTU_3D(I_LP,k))*DM9/86400.D0
               endif
            END DO
         END DO
      END DO  ! k = IzFuelBot, IzFuelTop

#ifdef JR_SRCTRM
      IF (Flag_SNF_pin) THEN
         temp1 = 1d+24
         pin_num = 0
         temp_pin_mid_1 = 0d0
         temp_pin_mid_2 = 0d0
         temp_pin_sum_1 = 0d0
         temp_pin_sum_2 = 0d0
         DO jp = 1, npin
            DO ip = 1, npin
               if (HFF(Ixy_1N,SNF_k,ip,jp)==0d0) CYCLE
               pin_num = pin_num+1
               temp_pin_sum_1 = temp_pin_sum_1+phihom(ip,jp,Ixy_1N,SNF_k,1)
               temp_pin_sum_2 = temp_pin_sum_2+phihom(ip,jp,Ixy_1N,SNF_k,2)
            ENDDO
         ENDDO
         temp_pin_mid_1 = temp_pin_sum_1 / real(pin_num)
         temp_pin_mid_2 = temp_pin_sum_2 / real(pin_num)
         temp_sum_1 = 0d0
         temp_sum_2 = 0d0
         temp_mid_1 = 0d0
         temp_mid_2 = 0d0
         DO m = 1, 4
            temp_sum_1 = (temp_sum_1 + Flux(I_1Nto4N(Ixy_1N,m),SNF_k,1))
            temp_sum_2 = (temp_sum_2 + Flux(I_1Nto4N(Ixy_1N,m),SNF_k,2))
         ENDDO
         temp_mid_1 =  temp_sum_1 / 4d0
         temp_mid_2 =  temp_sum_2 / 4d0
         pin_num = 0
         DO jp = 1, npin
            DO ip = 1, npin
               pin_num = pin_num+1
            ENDDO
         ENDDO
         if (temp1 == 0d0) temp1 = 1d0
         DO jp = 1, int(npin)
            DO ip = 1, int(npin)
               write(*,*) ip, jp, phihom_M(ip,jp,Ixy_1N,SNF_k,1), phihom_M(ip,jp,Ixy_1N,SNF_k,2)
            ENDDO
         ENDDO
      END IF
#endif


      Buff_Sum = D0
      pmax     = D0
      DO jp = 1, npin
         DO ip = 1, npin
            DO k = IzFuelBot, IzFuelTop
               powvalr(ip, jp, 0) = powvalr(ip, jp, 0) + powvalr(ip, jp, k)
            END DO
            Buff_Sum(1) = Buff_Sum(1) + powvalr(ip, jp, 0)
            IF ( powvalr(ip, jp, 0) > pmax ) THEN
               pmax  = powvalr(ip, jp, 0)
               ipmax = ip
               jpmax = jp
            END IF
         END DO
      END DO

      ! Store 2-D Peak Values
      pppeak(lfa, 0) = pmax
      Peak_X(lfa, 0) = ip
      Peak_Y(lfa, 0) = jp
      fnorm = Buff_Sum(1) / (npin*npin) !js+ N_Pin !wrong way but best now..
      rnorm = D1 / max(1d-10,fnorm)
      DO jp = 1, npin
         DO ip = 1, npin
            PinPOW_2D(ia,ja,ip,jp) = Normal_Power_XY(ia,ja)*Avg_Power*powvalr(ip,jp,0)*rnorm
            PinQ_2D(ia,ja,ip,jp) = PinPOW_2D(ia,ja,ip,jp)*(h_FA/npin)**D2*hactive
            if (flag_savepbu) then
               PinBU_2D(ia,ja,ip,jp)=PinBU_2D(ia,ja,ip,jp) &
                                   &+dT*PinQ_2D(ia,ja,ip,jp)/max(1d-10,PinMTU_2D(I_LP))*DM9/86400.D0
            endif
         END DO
      END DO

      RETURN
      END SUBROUTINE ffmult


      SUBROUTINE fsncnp(kf,nf,h,ka,x,bct,fsn,fsntl,fcn,fcntl)

      IMPLICIT NONE

      LOGICAL(1) :: kf, nf
      REAL(8) :: fsn,fsntl,fcn,fcntl,ka,x,bct,h,xsinh


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fsncnp] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF(nf) THEN                                                  ! crit.node
         IF(kf) THEN
            fsn=h-x*x*h/6+(x**4)*h/120
            fsntl=h-((x/max(ka,1d-30)*bct)**2)*h/6+((x/max(ka,1d-30)*bct)**4)*h/120      ! crit.node
            fcn=COS(x)                                              ! crit.node
            fcntl=COS(x/max(ka,1d-30)*bct)
         ELSE
            fsn=h+x*x*h/6+(x**4)*h/120
            fsntl=h+((x/max(ka,1d-30)*bct)**2)*h/6+((x/max(ka,1d-30)*bct)**4)*h/120      ! crit.node
            fcn=SQRT(1+(fsn*x/h)**2)
            fcntl=SQRT(1+(fsntl*(x/h/max(ka,1d-30)*bct))**2)
         END IF
      ELSE                                                         ! crit.node
         IF(kf) THEN
            fsn=SIN(x)*h/max(x,1d-30)
            fcn=COS(x)
         ELSE
            xsinh=SINH(x)
            fsn = xsinh*h/max(x,1d-30)
            fcn=SQRT(1+xsinh**2)
         END IF
         fcntl=fcn
         fsntl=fsn
      END IF                                                        ! crit.node
      RETURN
      END SUBROUTINE fsncnp


      SUBROUTINE calppeak
      USE Inc_FA
      USE Inc_Time
      USE Inc_XYZ
      USE mod_Alloc
      IMPLICIT NONE
      INTEGER :: l, ia, ja, jp, ip, k
      REAL(8), allocatable :: tmp_Fq(:,:,:,:,:), tmp_Fr(:,:,:,:)
      INTEGER :: Ifap, Ifa
      real(8) :: tmp_pp
      real(8) :: tmp_pv
      real(8), allocatable :: tmp_Fz(:)
      real(8) :: tmp_pz
      real(8) :: tmp_az
      real(8), allocatable :: tmp_FdelH(:,:,:,:)
      integer :: tmp_pi
      real(8), allocatable :: Fxy_val(:,:,:)
      integer :: nz_active, kk
      real(8),allocatable :: accumul_H(:), percent_H(:)
      real(8) :: tot_h, tol_izfxy
      logical(1), save :: flag_first_izfxy=.true.
      integer, save :: izfxybot, izfxytop
      INTEGER :: Ix_1N, Iy_1N

      !js+ modifying for calculation of pin peaking factor Fxy, Fr, Fz, Fq, FdH
      ! Fq = (Max of Power Density of Pin) / (Average Power Density of pin)
      ! Average Power Density

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [calppeak] in Mod_Pinpow'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      allocate(tmp_Fq(Nx_1N, Ny_1N, Nz, npin, npin))
      tmp_pp=0d0
      tmp_pv=0d0
      tmp_pi=0
      DO Ifap = 1, Nxy_FA_1N
         Ifa = pploc(Ifap)
         l  = lfatol(lfaptr(Ifa-1)+1)
         ia = Ix_4Nto1N(ltox(l))
         ja = Iy_4Nto1N(ltoy(l))
         DO k = izfuelbot, izfueltop
            DO jp = 1, npin
               DO ip = 1, npin
                   if(abs(PinPOW_3D(ia, ja, k, ip, jp)-0d0)>1d-10) tmp_pi=tmp_pi+1
                   tmp_pp=tmp_pp + PinPOW_3D(ia, ja, k, ip, jp)
               END DO
            END DO
         END DO
      END DO
      tmp_pp = tmp_pp/max(1d-10,real(tmp_pi,8))
      tmp_Fq = PinPOW_3D/max(1d-10,tmp_pp)
      Fq_val = maxval(tmp_Fq)
      Fq_loc = maxloc(tmp_Fq)
      deallocate(tmp_Fq)

      ! Fr = (Max of Axially integrated Power Density of Pin) / (Average Power Density of pin)
      allocate(tmp_Fr( Nx_1N, Ny_1N, npin, npin))
      tmp_pi=0
      DO Ifap = 1, Nxy_FA_1N
         Ifa = pploc(Ifap)
         l  = lfatol(lfaptr(Ifa-1)+1)
         ia = Ix_4Nto1N(ltox(l))
         ja = Iy_4Nto1N(ltoy(l))
         DO jp = 1, npin
            DO ip = 1, npin
               if(abs(PinPOW_2D(ia, ja, ip, jp)-0d0)>1d-10) tmp_pi=tmp_pi+1
               tmp_pp = tmp_pp + PinPow_2D(ia, ja, ip, jp)
            END DO
         END DO
      END DO
      tmp_pp=tmp_pp/max(1d-10,real(tmp_pi,8))
      tmp_Fr=PinPow_2d/max(1d-10,tmp_pp)
      Fr_val = maxval(tmp_Fr)
      Fr_loc = maxloc(tmp_Fr)
      deallocate(tmp_Fr)

      ! Fz = (Max of Radially integrated Power Density of pin) / (Average Power Density of pin)
      allocate(tmp_Fz(Nz))
      tmp_Fz = 0d0
      tmp_pz = 0d0
      DO k = izfuelbot, izfueltop
         tmp_pi=0
         tmp_pp=0d0
         tmp_az=0d0
         DO Ifap = 1, Nxy_FA_1N
            Ifa = pploc(Ifap)
            l  = lfatol(lfaptr(Ifa-1)+1)
            ia = Ix_4Nto1N(ltox(l))
            ja = Iy_4Nto1N(ltoy(l))
            DO jp = 1, npin
               DO ip = 1, npin
                  if(abs(PinPOW_3D(ia, ja, k, ip, jp)-0d0)>1d-10) tmp_pi=tmp_pi+1
                  tmp_pp = tmp_pp + PinPow_3D(ia, ja, k, ip, jp) * (h_FA/npin)**2d0 * MeshSize_z(k)
                  tmp_az = tmp_az + (h_FA/npin)**2d0 * MeshSize_z(k)
               END DO
            END DO
         END DO
         tmp_Fz(k) = tmp_pp / max(1d-10,tmp_az)
         tmp_pz = tmp_pz + tmp_Fz(k)
      END DO
      tmp_pz = tmp_pz / (izfueltop - izfuelbot + 1)
      tmp_Fz = tmp_Fz / max(1d-10,tmp_pz)
      Fz_ave = sum(tmp_Fz)/(izfueltop - izfuelbot + 1)
      Fz_val = maxval(tmp_Fz)
      Fz_loc = maxloc(tmp_Fz)
      deallocate(tmp_Fz)

      ! FdH = (Integral of LinPow of Peak Rod) / (Average rod power)
      allocate(tmp_FdelH( Nx_1N, Ny_1N, npin, npin))
      tmp_pi=0
      tmp_pp=0d0
      DO Ifap = 1, Nxy_FA_1N
         Ifa = pploc(Ifap)
         l  = lfatol(lfaptr(Ifa-1)+1)
         ia = Ix_4Nto1N(ltox(l))
         ja = Iy_4Nto1N(ltoy(l))
         DO jp = 1, npin
            DO ip = 1, npin
               if ((abs(PinQ_2D(ia,ja,ip,jp)-0d0)>1d-10)) tmp_pi=tmp_pi+1
               tmp_pp = tmp_pp + PinQ_2D(ia, ja, ip, jp)
            END DO
         END DO
      END DO
      tmp_pp = tmp_pp/max(1d-10,real(tmp_pi,8))
      tmp_FdelH = PinQ_2D/max(1d-10,tmp_pp)
      FdH_val = maxval(tmp_FdelH)
      FdH_val_XY = 0d0
      DO Ifap = 1, Nxy_FA_1N
         Ifa = pploc(Ifap)
         l  = lfatol(lfaptr(Ifa-1)+1)
         ia = Ix_4Nto1N(ltox(l))
         ja = Iy_4Nto1N(ltoy(l))
         FdH_val_XY(ia,ja) = maxval(tmp_FdelH(ia,ja,:,:))
      END DO
      FdH_loc = maxloc(tmp_FdelH)
      deallocate(tmp_FdelH)

      ! Fxy
      !define active axial plane for fxy calculation
      if (flag_first_izfxy) then
         flag_first_izfxy=.false.
         nz_active=izfueltop-izfuelbot+1
         allocate(accumul_H(nz_active))
         tot_h=0d0
         kk=1
         do k=izfuelbot,izfueltop
            if (kk==1) then
               accumul_H(kk)=GridSize_z(k)
            else
               accumul_H(kk)=accumul_H(kk-1)+GridSize_z(k)
            endif
            kk=kk+1
         enddo
         tot_h=accumul_H(nz_active)
         allocate(percent_H(nz_active))
         do k=1,nz_active
            percent_H(k)=accumul_H(k)/tot_h
         enddo
         tol_izfxy=0.1d0 ! ignore TOP/BOTTOM 10% axial plane
         do k=nz_active,1,-1
            izfxytop=k
            if (percent_H(k)<1d0-tol_izfxy) exit
         enddo
         izfxybot=nz_active-izfxytop+1
         deallocate(accumul_H)
         deallocate(percent_H)
         if (izfxybot>izfxytop) then
            izfxybot=izfuelbot
            izfxytop=izfueltop
         endif
      endif

      allocate(Fxy_val( Nx_1N, Ny_1N, Nz))
      Fxy_val=0d0
      DO k = izfxybot, izfxytop
         tmp_pp = 0d0
         tmp_pi = 0
         DO Ifap = 1, Nxy_FA_1N
            Ifa = pploc(Ifap)
            l  = lfatol(lfaptr(Ifa-1)+1)
            ia = Ix_4Nto1N(ltox(l))
            ja = Iy_4Nto1N(ltoy(l))
            DO jp = 1, npin
               DO ip = 1, npin
                  if ((abs(PinPow_2D(ia,ja,ip,jp)-0d0)>1d-10)) tmp_pi=tmp_pi+1
                  tmp_pp = tmp_pp + PinPow_3D(ia, ja, k, ip, jp)
               END DO
            END DO
         END DO
         tmp_pp = tmp_pp/max(1d-10,real(tmp_pi,8))
         DO Ifap = 1, Nxy_FA_1N
            Ifa = pploc(Ifap)
            l  = lfatol(lfaptr(Ifa-1)+1)
            ia = Ix_4Nto1N(ltox(l))
            ja = Iy_4Nto1N(ltoy(l))
            Fxy_val(ia,ja,k)=maxval(PinPow_3D(ia,ja,k,:,:))/max(1d-10,tmp_pp)
         END DO
      END DO
      max_Fxy_val = maxval(Fxy_val)
      max_Fxy_loc = maxloc(Fxy_val)
      DO Iy_1N=1,Ny_1N
         DO Ix_1N=1,Nx_1N
            Fxy_val_XY(Ix_1N, Iy_1N) = maxval(Fxy_val(Ix_1N,Iy_1N,:))
         END DO
      END DO
      deallocate(Fxy_val)

      ! overwrite the peaking factor from wk for GRP
      Fq  = Fq_val
      Fxy = max_Fxy_val
      Fr  = Fr_val
      Fdh = FdH_val

      ! Find Max PinBU_2D, PinPOW_2D, PinPOW_3D with location
      max_PinBU_2D  = maxval(PinBU_2D)
      max_PinPOW_2D = maxval(PinPOW_2D)
      max_PinPOW_3D = maxval(PinPOW_3D)
      max_PinQ_2D = maxval(PinQ_2D)
      max_PinQ_3D = maxval(PinQ_3D)
      max_PinBU_2D_loc  = maxloc(PinBU_2D)
      max_PinPOW_2D_loc = maxloc(PinPOW_2D)
      max_PinPOW_3D_loc = maxloc(PinPOW_3D)
      max_PinQ_2D_loc = maxloc(PinQ_2D)
      max_PinQ_3D_loc = maxloc(PinQ_3D)

      !js+ modifying for calculation of pin peaking factor Fxy, Fr, Fz, Fq, FdH
      END SUBROUTINE calppeak


      END MODULE Mod_Pinpow


#endif
