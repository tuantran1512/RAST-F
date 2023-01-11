#ifdef jr_vver

      MODULE Mod_TPENDrive

      USE Inc_TPEN
      USE Inc_Option,   ONLY: n_group
      USE Inc_Geometry!, ONLY: GridSize_x
      USE Inc_Control!,  ONLY: rerrl2
      USE Inc_3D
      USE Inc_maXS
      USE Inc_Kinetics, ONLY: beta_d_Tot

#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction
     ! use Mod_PinPow_Hex
      use Inc_PinPow_Hex
      ! -=-=-=-=-=-=-=-=-
#endif


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE TpenDriver

      IMPLICIT NONE!
! drive TPEN calculation to update radial and axial nodal
! coupling coefficicents!
      LOGICAL fsweep
      CHARACTER(len=65) :: amesg
      REAL(8) :: rtpen, errl2tpen1, rtpend
      REAL(8) :: rtpendd
      INTEGER(4) :: isweep0, isweep, nswpmin, nsweep
      DATA isweep0/0/
!
      
      REAL(8) :: epsr2=0.001
      fsweep = .TRUE.

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [TpenDriver] in Mod_TPENDrive'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      ntnodal=ntnodal+1
      !tlap=dclock()
!
! set tpen boundary condition
!     !! CALL TpenBc
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
! 40 pcm difference between TpenBc and TpenBcMg2
!      IF (Flag_transient) THEN
!          CALL TpenBcMg2
!      ELSE
          CALL TpenBc
!      ENDIF
      
!      CALL TpenBc
!
      nsweep=2*n_group
      errl2tpen1=rerrl2
    !!  if (n_group > 2) then
         rtpen=residtpen_mg()
    !!  else
    !!     rtpen=residtpen()
    !!  endif
      rtpend=rtpen
      DO isweep=1,10
!         timptst=dclock()
! solve for point flux
       !!  if (n_group>2) then
            CALL SolpFlxMg

       !!  else
       !!     CALL SolpFlx
       !!  endif
!         timpted=dclock()
         !tpoint=tpoint+timpted-timptst
         !! twogroup=.not.multigroup
         !if(twogroup) then
         !call SwpTpen(err,1,nassy,1)
         !else
       !!  if (n_group>2) then
            call SwpTpenMg(1,nassy,1,1,nz,1)

       !!  else
       !!     call SwpTpen(1,nassy,1,1,nz,1)
       !!  endif
         !endif
!
! check convergence of residual and exit sweep if needed

         rtpendd=rtpend
         rtpend=rtpen
       !!  if (n_group > 2) then
            rtpen=residtpen_mg()
       !!  else
       !!     rtpen=residtpen()
       !!  endif
         IF( n_group .LE. 12 ) THEN
            nswpmin=2
         ELSE
            nswpmin=4
         ENDIF
         IF(isweep.GT.nswpmin) THEN
            IF(rtpen.LT.rtpend .AND. rtpend.LT.rtpendd) THEN
               IF(rtpen/rtpend .LT. 1.1*rtpend/rtpendd) THEN
                  fsweep = .FALSE.
                  EXIT
               END IF
               IF(iptype.EQ.1 .AND. rtpen.LT.epsr2) THEN
                  fsweep = .FALSE.
                  EXIT
               END IF
            ENDIF
            IF(rtpen.GT.rtpend .AND. rtpend.GT.rtpendd) THEN
               fsweep = .FALSE.
               EXIT
            END IF
         ENDIF
!         write(*,*) 'isweep',isweep
!         timpened=dclock()
         !ttpen=ttpen+timpened-timpted
      ENDDO
      IF(fsweep) isweep=isweep-1
      ntsweep=ntsweep+isweep
!
! update CMFD coupling coefficients
    !!  if (n_group>2) then
         CALL UpdDhatMg
    !!  else
    !!     CALL UpdDhat
    !!  endif
!      tnodal=tnodal+dclock()-tlap
      WRITE(*,'(a,2i5)') "TPEN Nodal update...",ntnodal,ntsweep
!      WRITE(amesg,'(a,2i5)') "TPEN Nodal update...",ntnodal,ntsweep
      !CALL message(.TRUE. ,popt(2),.TRUE.,amesg)
!stop

      RETURN
      END SUBROUTINE TpenDriver

      SUBROUTINE UpdDhatMg
!
! Update CMFD coefficient(dhat and beta) from Nodal Solution
!
      IMPLICIT NONE
      INTEGER(4) :: k,l,m,ig,ih,it,md,mf,ind,id1,id2,is1,is2,ig2,iz,iz1,iz2,mvz,isfc
      REAL(8) :: w,hzr,cnts,dfdm,cnto1,cnto2,flxs,flx1,flx2,fdflxs,sumcnt,sumflx
      REAL(8) :: xsd1,xssup,xsslocal
      DIMENSION mvz(2)
      DATA mvz/2,1/
      REAL(8) :: rt3, sqrt3, rsqrt3, hside
      REAL(8) :: hexarea
!      integer(4) :: ndivhs


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [UpdDhatMg] in Mod_TPENDrive'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      sqrt3=SQRT(3.d0)
      rt3=1.73205080756888
      rsqrt3=1d0/sqrt3
      hside=gridsize_x(1)*rsqrt3
      hexarea=1.5d0*gridsize_x(1)*hside
!
!
!
! variables for nem solver in z-direction
! atleak : transverse leakage for axial direction
! atleaku,atleakl : transverse leakage for upper and lower node
! xsdu,xsdl : xsd(diffusion coefficients) for upper and lower node
!
!
! update D(diffusion coefficient) for CMFD
      DO iz=1,nz
         DO ih=1,nassy
            DO ig=1,2
               sumflx=0
               DO ig2=mgb(ig),mge(ig)
                  sumflx=sumflx+hflxf(ig2,ih,iz)
               ENDDO
               hflx(ig,ih,iz)=sumflx
               flux(imap(ih),iz,ig)=hflx(ig,ih,iz)  !temp
               d_3d(imap(ih),iz,ig)=0
               DO ig2=mgb(ig),mge(ig)
                  fohflx(ig2,ih,iz)=fhflx(ig2,ih,iz)
                  fhflx(ig2,ih,iz)=hflxf(ig2,ih,iz)/sumflx
                  fluxf(imap(ih),iz,ig2)=fhflx(ig2,ih,iz)
                  d_3d(imap(ih),iz,ig)=d_3d(imap(ih),iz,ig) &
                       +1/(3*maXS_tr_3D_MG(imap(ih),iz,ig2))*fluxf(imap(ih),iz,ig2)
               ENDDO
               maXS_tr_3D(imap(ih),iz,ig)=1/(d_3d(imap(ih),iz,ig)*3)
            ENDDO
            DO it=1,6
               DO ig=1,2
                  sumcnt=0
                  DO ig2=mgb(ig),mge(ig)
                     sumcnt=sumcnt+cnto(ig2,it,ih,iz)
                  ENDDO
                  DO ig2=mgb(ig),mge(ig)
                     focnto(ig2,it,ih,iz)=fcnto(ig2,it,ih,iz)
                     fcnto(ig2,it,ih,iz)=cnto(ig2,it,ih,iz)/sumcnt
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      IF(nz.GT.1) THEN
         DO iz=1,nz
            DO ih=1,nassy
               DO it=1,2
                  DO ig=1,2
                     sumcnt=0
                     DO ig2=mgb(ig),mge(ig)
                        sumcnt=sumcnt+cntzo(ig2,it,ih,iz)
                     ENDDO
                     DO ig2=mgb(ig),mge(ig)
                        focntzo(ig2,it,ih,iz)=fcntzo(ig2,it,ih,iz)
                        fcntzo(ig2,it,ih,iz)=cntzo(ig2,it,ih,iz)/sumcnt
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      DO isfc=1,nxsfc
         ind=neigsnd(3,isfc)
         IF(ind.EQ.12) THEN
            is1=neigsnd(1,isfc)
            is2=neigsnd(2,isfc)
            id1=neigsnd(4,isfc)
            id2=neigsnd(5,isfc)
            DO iz=1,nz
               DO ig=1,2
                  cnto1=0
                  cnto2=0
                  DO ig2=mgb(ig),mge(ig)
                     cnto1=cnto1+cnto(ig2,id1,is1,iz)
                     cnto2=cnto2+cnto(ig2,id2,is2,iz)
                  ENDDO
                  flxs=2.*(cnto1+cnto2)
                  cnts=cnto1-cnto2
                  flx1=hflx(ig,is1,iz)
                  flx2=hflx(ig,is2,iz)
                  dfdm=dfd(ig,isfc,iz)
                  xsd1=d_3d(imap(is1),iz,ig)
                  w=dfdm/xsd1/2.
                  fdflxs=(1.-w)*flx1+w*flx2
                  dhat(ig,isfc,iz)= &
                       -(rt3*hside*cnts+dfdm*(flx2-flx1))/(flx2+flx1)
                  betaphis(ig,isfc,iz)=2.*(flxs-fdflxs)/(flx1+flx2)
               ENDDO
            ENDDO
         ELSE
            is1=neigsnd(ind,isfc)
            id1=neigsnd(ind+3,isfc)
            DO iz=1,nz
               DO ig=1,2
                  cnto1=0
                  cnto2=0
                  DO ig2=mgb(ig),mge(ig)
                     cnto1=cnto1+cnto(ig2,id1,is1,iz)
                     cnto2=cnto2+reflratf(ig2)*cnto(ig2,id1,is1,iz) ! y
                  ENDDO
                  flxs=2.*(cnto1+cnto2)
                  cnts=cnto1-cnto2
                  flx1=hflx(ig,is1,iz)
                  dfdm=dfd(ig,isfc,iz)
                  dhat(ig,isfc,iz)=(dfdm*flx1-rt3*hside*cnts)/flx1
                  IF(ind.EQ.2) dhat(ig,isfc,iz)=-dhat(ig,isfc,iz)
                  betaphis(ig,isfc,iz)=flxs/flx1
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      IF(nz /= 1) THEN
         !! update D in z-direction for CMFD
         DO iz=1,nz+1
            ind=neigsndz(3,iz)
            IF(ind.EQ.12) THEN
               iz1=neigsndz(1,iz)
               iz2=neigsndz(2,iz)
               id1=neigsndz(4,iz)
               id2=neigsndz(5,iz)
               hzr=hz(iz2)/hz(iz1)
               DO ih=1,nassy
                  DO ig=1,2
                     cnto1=0
                     cnto2=0
                     DO ig2=mgb(ig),mge(ig)
                        cnto1=cnto1+cntzo(ig2,id1,ih,iz1)
                        cnto2=cnto2+cntzo(ig2,id2,ih,iz2)
                     ENDDO
                     flxs=2.*(cnto1+cnto2)
                     cnts=cnto1-cnto2
                     flx1=hflx(ig,ih,iz1)
                     flx2=hflx(ig,ih,iz2)
                     dfdm=dfdz(ig,ih,iz)
                     xsd1=d_3d(imap(ih),iz1,ig)
                     w=dfdm/xsd1/(1.+hzr)
                     fdflxs=(1.-w)*flx1+w*flx2
                     dhatz(ig,ih,iz)=-((hz(iz1)+hz(iz2))*cnts*0.5 &
                          +dfdm*(flx2-flx1))/(flx2+flx1)
                     betaphisz(ig,ih,iz)=2.*(flxs-fdflxs)/(flx1+flx2)
!                     WRITE(1995,*) ig,ih,iz, ind,  dhatz(ig,ih,iz),betaphisz(ig,ih,iz)
                  ENDDO
               ENDDO
            ELSE
               iz1=neigsndz(ind,iz)
               id1=neigsndz(ind+3,iz)
               IF(iz.EQ.1) THEN
                  DO ig=1,n_group
                     reflratzf(ig)=reflratzbf(ig)
                  ENDDO
               ELSE
                  DO ig=1,n_group
                     reflratzf(ig)=reflratztf(ig)
                  ENDDO
               ENDIF
               DO ih=1,nassy
                  DO ig=1,2
                     cnto1=0
                     cnto2=0
                     DO ig2=mgb(ig),mge(ig)
                        cnto1=cnto1+cntzo(ig2,id1,ih,iz1)
                        cnto2=cnto2+reflratzf(ig2)*cntzo(ig2,id1,ih,iz1)
                     ENDDO
                     flxs=2.*(cnto1+cnto2)
                     cnts=cnto1-cnto2
                     flx1=hflx(ig,ih,iz1)
                     dfdm=dfdz(ig,ih,iz)
                     dhatz(ig,ih,iz)=(dfdm*flx1-hz(iz1)*cnts)/flx1
                     IF(ind.EQ.2) dhatz(ig,ih,iz)=-dhatz(ig,ih,iz)
                     betaphisz(ig,ih,iz)=flxs/flx1
!                     WRITE(1995,*) ig,ih,iz, ind,  dhatz(ig,ih,iz),betaphisz(ig,ih,iz)

                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      END IF
      !! update old abs xsec
      DO k=1,nz
         DO l=1,nxy
            xsslocal=0
            DO md=mgb(2),mge(2)
               DO m=mgb(1),mge(1)
                  xsslocal=xsslocal+xssf(m,md,l,k)*fluxf(imap(l),k,m)
               ENDDO
            ENDDO
            xssup=0
            DO md=mgb(1)+1,mge(1)
               DO m=mgb(2),mge(2)
                  xssup=xssup+xssf(m-1,md,l,k)*fluxf(imap(l),k,m)
               ENDDO
            ENDDO
            xssup=xssup*flux(imap(l),k,2)/flux(imap(l),k,1)
            DO m=1,2
               xstd(m,l,k)=0
               DO mf=mgb(m),mge(m)
                  xstd(m,l,k)=xstd(m,l,k)+maXS_a_3D_MG(imap(l),k,mf)*fluxf(imap(l),k,mf)
               ENDDO
            ENDDO
            xstd(1,l,k)=xstd(1,l,k)+xsslocal-xssup
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE UpdDhatMg

      SUBROUTINE SwpTpenMg(irfrom,irto,irstp,iafrom,iato,iastp)

#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction
      use Inc_PinPow_Hex
      ! -=-=-=-=-=-=-=-=-
#endif

      IMPLICIT NONE!
! control NEM Solver for z-direction!
      INTEGER(4) :: m,ig,ih,it,iz,nn,md,irfrom,irto,irstp,icg
      INTEGER(4) :: iafrom,iato,iastp,jupstt,neigup,neigdn
      REAL(8) :: sum,rhz,atleakm,areavol,cntohom,cntihet
      REAL(8) :: cntihom,cntohet,errl2tpen,flxl2tpen
! variables for nem solver in z-direction
! atleak : transverse leakage for axial direction
! atleaku,atleakl : transverse leakage for upper and lower node
! xsdu,xsdl : xsd(diffusion coefficients) for upper and lower node
! cntz0i,cntz1i : incoming currents from lower and upper surface
!
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: atleakl,atleaku
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: xsdu,xsdl,cntz0i,cntz1i,chieff
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:,:) :: xssmg,cnti,pbflx, srczn
      REAL(8) :: mv(6), hzm, hzu, hzl, xsdm
      LOGICAL ifbcref(100),ifbcrefzb(100),ifbcrefzt(100)
      SAVE areavol
      DATA mv/4,5,6,1,2,3/
      REAL(8) :: rt3, sqrt3, rsqrt3, hside
      REAL(8) :: hexarea
!      integer(4) :: ndivhs

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SwpTpenMg] in Mod_TPENDrive'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      sqrt3=SQRT(3.d0)
      rt3=1.73205080756888
      rsqrt3=1d0/sqrt3
      hside=gridsize_x(1)*rsqrt3
      hexarea=1.5d0*gridsize_x(1)*hside
!
!

      IF(.NOT.ALLOCATED(atleakl)) THEN
         areavol=2.*rt3/9.D0/hside
         ALLOCATE(atleakl(n_group),atleaku(n_group),xsdu(n_group),xsdl(n_group) &
              ,          cntz0i(n_group),cntz1i(n_group),chieff(n_group))
         atleakl=0
         atleaku=0
         xsdu=0
         xsdl=0
         cntz0i=0
         cntz1i=0
         chieff=0
         ALLOCATE(xssmg(n_group,n_group),cnti(n_group,ntph),pbflx(n_group,ntph) &
              ,          srczn(n_group,ntph))
         xssmg=0
         cnti=0
         pbflx=0
         srczn=0
      ENDIF
 
      DO ig=1,n_group
         ifbcref(ig)=.FALSE.
         ifbcrefzb(ig)=.FALSE.
         ifbcrefzt(ig)=.FALSE.
         IF(alxrf(ig).EQ.0) ifbcref(ig)=.TRUE.
         IF(alzlf(ig).EQ.0) ifbcrefzb(ig)=.TRUE.
         IF(alzrf(ig).EQ.0) ifbcrefzt(ig)=.TRUE.
      ENDDO
      !! average transverse leakage cal. for z-direction
      DO ih=1,nassy
         DO iz=1,nz
            DO ig=1,n_group
               atleak(ig,ih,iz)=-sefve(ig,ih,iz)
               srcz(ig,ih,iz)=sefve(ig,ih,iz)
            ENDDO
         ENDDO
      ENDDO

      IF(nassy /= 1) THEN
         DO ih=1,nassy
            DO iz=1,nz
               DO ig=1,n_group
                  sum=0
                  DO it=1,6
                     nn=neignd(it,ih)
                     IF(nn.EQ.0) THEN
                        sum=sum+(1.-reflratf(ig))*cnto(ig,it,ih,iz)
                     ELSE
                        sum=sum+cnto(ig,it,ih,iz)-cnto(ig,neigjin(it,ih),nn,iz)
                     ENDIF
                  ENDDO
                  atleak(ig,ih,iz)=atleak(ig,ih,iz)+sum*areavol
               ENDDO
            ENDDO
         ENDDO
      !! Axial Source for TPEN cal.
      END IF

      IF(nz /= 1) THEN
         DO iz=1,nz
            rhz=1./hz(iz)
            IF(iz.EQ.1) THEN
               DO ih=1,nassy
                  DO ig=1,n_group
                     srcz(ig,ih,iz)=srcz(ig,ih,iz) &
                          +((reflratzbf(ig)-1.)*cntzo(ig,1,ih,iz) &
                          +cntzo(ig,1,ih,iz+1)-cntzo(ig,2,ih,iz))*rhz
                  ENDDO
               ENDDO
            ELSE
               IF(iz.EQ.nz) THEN
                  DO ih=1,nassy
                     DO ig=1,n_group
                        srcz(ig,ih,iz)=srcz(ig,ih,iz) &
                             +(cntzo(ig,2,ih,iz-1)-cntzo(ig,1,ih,iz) &
                             +(reflratztf(ig)-1.)*cntzo(ig,2,ih,iz))*rhz

                     ENDDO
                  ENDDO
               ELSE
                  DO ih=1,nassy
                     DO ig=1,n_group
                        srcz(ig,ih,iz)=srcz(ig,ih,iz) &
                             +(cntzo(ig,2,ih,iz-1)+cntzo(ig,1,ih,iz+1) &
                             -cntzo(ig,1,ih,iz)-cntzo(ig,2,ih,iz))*rhz

                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      !! solve from irfrom to irto in radial direction
      !!   and from iafrom to iato in axial direction
      END IF

      errl2tpen=0
      flxl2tpen=0
      DO iz=iafrom,iato,iastp
         DO ih=irfrom,irto,irstp
            !! preparation for radial solver(cnti,srczn,pbflx)
            IF(nassy.GT.1) THEN
               DO it=1,6
                  nn=neignd(it,ih)
                  IF(nn.EQ.0) THEN
                     DO ig=1,n_group
                        cnti(ig,it)=cnto(ig,it,ih,iz)*reflratf(ig)
                        IF(ifbcref(ig)) THEN
                           srczn(ig,it)=srcz(ig,ih,iz)
                        ELSE
                           srczn(ig,it)=0.
                        ENDIF
                     ENDDO
                  ELSE
                     DO ig=1,n_group
                        cnti(ig,it)=cnto(ig,neigjin(it,ih),nn,iz)
                        srczn(ig,it)=srcz(ig,nn,iz)
                     ENDDO
                  ENDIF
               ENDDO
               DO it=1,6
                  nn=neigpt(it,ih)
                  DO ig=1,n_group
                     cntihet=cnti(ig,it)
                     cntohet=cnto(ig,it,ih,iz)
                     if (xsadf(ig,it,ih,iz) .EQ. 0) then
                         cnti(ig,it)=0.5*((cntohet+cntihet)/1d0 &
                              -(cntohet-cntihet))
                     else
                         cnti(ig,it)=0.5*((cntohet+cntihet)/xsadf(ig,it,ih,iz) &
                              -(cntohet-cntihet))
                     endif
                     pbflx(ig,it)=pflx(ig,nn,iz)/adfpt(ig,it,ih,iz)
                  ENDDO
               ENDDO
            ENDIF
            !! preparation for axial solver(cnti,srczn,pbflx)
            hzm=hz(iz)
            hzu=hzm
            hzl=hzm
            IF(nz.GT.1) THEN
               neigdn=neigz(1,iz)
               neigup=neigz(2,iz)
               DO ig=1,n_group
                  atleakm=atleak(ig,ih,iz)
                  atleakl(ig)=atleakm
                  atleaku(ig)=atleakm
                  ! xsdf = xsd = d_3d
                  xsdm=d_3d_MG(imap(ih),iz,ig)
                  xsdl(ig)=xsdm
                  xsdu(ig)=xsdm
               ENDDO
               IF(neigdn.NE.0) THEN
                  DO ig=1,n_group
                     cntz0i(ig)=cntzo(ig,2,ih,neigdn)
                     atleakl(ig)=atleak(ig,ih,neigdn)
                     xsdl(ig)=d_3d_MG(imap(ih),neigdn,ig)
                     hzl=hz(neigdn)
                  ENDDO
               ELSE
                  DO ig=1,n_group
                     cntz0i(ig)=reflratzbf(ig)*cntzo(ig,1,ih,iz)
                     IF(.NOT.ifbcrefzb(ig)) atleakl(ig)=0
                  ENDDO
               ENDIF
               IF(neigup.NE.0) THEN
                  DO ig=1,n_group
                     cntz1i(ig)=cntzo(ig,1,ih,neigup)
                     atleaku(ig)=atleak(ig,ih,neigup)
                     xsdu(ig)=d_3d_MG(imap(ih),neigup,ig)
                     hzu=hz(neigup)
                  ENDDO
               ELSE
                  DO ig=1,n_group
                     cntz1i(ig)=reflratztf(ig)*cntzo(ig,2,ih,iz)
                     IF(.NOT.ifbcrefzt(ig)) atleaku(ig)=0
                  ENDDO
               ENDIF
            ENDIF

            DO ig=1,n_group
               icg=igc(ig)
               cntihet=cntz0i(ig)
               cntohet=cntzo(ig,1,ih,iz)
               cntz0i(ig)=0.5*((cntohet+cntihet)/adfb(icg,ih,iz) &
                    -(cntohet-cntihet))
               cntihet=cntz1i(ig)
               cntohet=cntzo(ig,2,ih,iz)
               cntz1i(ig)=0.5*((cntohet+cntihet)/adft(icg,ih,iz) &
                    -(cntohet-cntihet))
            ENDDO
            DO md=2,n_group
               DO m=1,md-1
                  xssmg(m,md)=xssf(m,md,ih,iz)
               ENDDO
               DO m=md+1,n_group  !temp
                  xssmg(m,md)=xssf(m-1,md,ih,iz)
               ENDDO
            ENDDO
            DO m=1,n_group
               chieff(m)=maXs_chi_3D_MG(imap(ih),iz,m) &
                    +beta_d_tot(imap(ih),iz)*(maXs_chid_3D_MG(imap(ih),iz,m)-maXs_chi_3D_MG(imap(ih),iz,m))
            ENDDO
            jupstt=lupscat

            !! xstf = xst = maXs_r_3D
            !! xsnff = xsnf = nu_maXs_f_3D
            !! xsdf = xsd = d_3d
!CALL CPU_TIME ( time_begin_2_61 )
            CALL OneTpenMg(nassy,nz,jupstt &
                 ,    iscatib,iscatie &
                 ,    maXs_r_3D_MG(imap(ih),iz,:),nu_maXs_f_3D_MG(imap(ih),iz,:),xssmg,d_3d_MG(imap(ih),iz,:) &
                 ,    chieff,hside,hzm,hzu,hzl &
                 ,    pbflx,cnti,cntz0i,cntz1i,sefve(:,ih,iz),srczn &
                 ,    atleaku,atleakl,reigv &
                 ,    aflx(:,:,ih,iz),xmom(:,:,ih,iz),ymom(:,:,ih,iz) &
                 ,    zmom1(:,ih,iz),zmom2(:,ih,iz) &
                 ,    cnto(:,:,ih,iz),cntzo(:,1,ih,iz),cntzo(:,2,ih,iz) &
                 ,    hflxf(:,ih,iz))
 !CALL CPU_TIME ( time_end_2_61 )
 !WRITE (*,*) 'Time of operation was 2_61 ', time_end_2_61 - time_begin_2_61, ' seconds'

            DO it=1,6
               DO ig=1,n_group
                  cntihom=cnti(ig,it)
                  cntohom=cnto(ig,it,ih,iz)
                  if (n_group>2) xsadf(ig,it,ih,iz)=1d0
                  cnto(ig,it,ih,iz)= &
                       0.5*((cntohom+cntihom)*xsadf(ig,it,ih,iz) &
                       +(cntohom-cntihom))
               ENDDO
            ENDDO
            DO ig=1,n_group
               icg=igc(ig)
               cntihom=cntz0i(ig)
               cntohom=cntzo(ig,1,ih,iz)
               cntzo(ig,1,ih,iz)=0.5*((cntohom+cntihom)*adfb(ig,ih,iz) &
                    +(cntohom-cntihom))
               cntihom=cntz1i(ig)
               cntohom=cntzo(ig,2,ih,iz)
               cntzo(ig,2,ih,iz)=0.5*((cntohom+cntihom)*adft(ig,ih,iz) &
                    +(cntohom-cntihom))
            ENDDO

         END DO
      END DO


      RETURN
      END SUBROUTINE SwpTpenMg

!//      SUBROUTINE OneTpenMg(nassy,nz,jupstt,iscatib,iscatie &
!//     ,    xst,xsnf,xss,xsd,chi,sz,hz,hzu,hzl,pflbt,cntit &
!//     ,    cntz0i,cntz1i,seff,srczn, dtlu,dtll,rxkeff, aflxt &
!//     ,    xmomt,ymomt,zmomt1,zmomt2,cntot,cntz0o,cntz1o,hflx)
!//
!//#ifdef siarhei_ppr
!//      ! Siarhei pin power reconstruction
!//      use Inc_PinPow_Hex
!//      ! -=-=-=-=-=-=-=-=-
!//#endif
!//
!//      IMPLICIT NONE
!//      REAL(8) :: rxkeff, sz, hz, hzu, hzl
!//      REAL(8) :: xst(n_group),xsnf(n_group),xsd(n_group),chi(n_group), &
!//                 xss(n_group,n_group)
!//      REAL(8) :: pflbt(n_group,6),cntit(n_group,6),cntz0i(n_group)
!//      REAL(8) :: cntz1i(n_group),seff(n_group),srczn(n_group,6), &
!//                 dtlu(n_group),dtll(n_group)
!//      REAL(8) :: cntot(n_group,6),cntz0o(n_group),cntz1o(n_group)
!//      REAL(8) :: hflx(n_group),aflxt(n_group,6),xmomt(n_group,6),ymomt(n_group,6)
!//      REAL(8) :: zmomt1(n_group),zmomt2(n_group)
!//      REAL(8) :: asrc(6),xsrc(6),ysrc(6)
!//      REAL(8) ::  sflb(6),sfle(6), tdt1(6),tdt2(6),tdt3(6)
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: tem1,tem2,hflxold
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: raxy,ca1,csfi2,rptj
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: cbd1,cbd2,cbd3,rbd,rsfsm
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: aax,czj1,czj2,ctl1i,cmom1
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: cmom2,ccpt1,ccpt2
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: pflcbd,sj1bd,sj0bd,csj
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: smm1bd,smm2bd,csmm1,csmm2
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: csmm2hfl,csmm2sjs,csmm1sj
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: xax,scgr,scgrb, cntzisum,cntisum,pflbsum
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:,:) :: asrcbd,xsrcbd,ysrcbd
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:,:) :: aflxbd,xmombd,ymombd,cntim
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:,:) :: pflbm,ccsf,sflebd,ssrcbbd
!//      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:,:) :: ccgr,rcgr
!//      REAL(8) :: cnto(6),aflx(6),xmom(6),ymom(6)
!//      REAL(8) :: onethree, onenine, onesix, one12, oneten, rsz
!//      REAL(8) :: rsz2, rfac1, rfac2, rfac3, areavol, rhz, rhz2, rhzl
!//      REAL(8) :: rhzu, faccsf, faccpt, faccbd, facmx, bt, gam
!//      REAL(8) :: rsf, csf1, rpt, cpt1, cpt2, cpt3, ttt, fac1, rsfi
!//      REAL(8) :: csfi1, rpti, cpti1,fac2, fac3, fac4, rsi1, rsi2
!//      REAL(8) :: rsi3, rsi4, t1, t2, t3, t4, ccpt3, ccpt4, q1, q2, q4
!//      REAL(8) :: tt, q3, q5, q7, ccsf11, ccsf12, ccsf13, ccsf14, ccsf15
!//      REAL(8) :: ccsf16, ccsf17, datj1, datj2, datj3, datc1, datc2, csf4
!//      REAL(8) :: csfi3, cpt4, gamz, rdet, rzj1, rzj2, ctl2i, rbtmb, rbtmu
!//      REAL(8) :: tl0bd, tl1bd, btz, dat, sumfs, sumcnti, sumcntzi, tfsrc
!//      REAL(8) :: tfsrcx, tfsrcy, ttt1, ttt2, ttt3, ttt4, tfsrcz1, tfsrcz2
!//      REAL(8) :: dat0, ad1, ad2, ad3, pflc, tfsrcz, xd1, xd2, xd3, yd1
!//      REAL(8) :: yd2, yd3, hflxr, dat1, smm1, smm2, sbal, sj1, sj0, sjsum
!//      REAL(8) :: dtlm, phit, phib, flxzsum, cntzosum, dat2, sumcnto
!//      REAL(8) :: sumscat, sumsrc, sumtot, sumcntzo, sumrhs, sumlhs
!//      REAL(8) :: errimbal, fcgr
!//      INTEGER(4) :: nassy, nz, jupstt, ii, j, i1, j1, iter, jstt, j2, i3
!//      INTEGER(4) :: iscatib(n_group),iscatie(n_group), nitrsrc, it, j3
!//      LOGICAL(1) ::  if2d,if1d,if3d, skip
!//      SAVE onenine,onesix,one12,oneten,onethree,rsz,rsz2,areavol &
!//           ,    rfac1,rfac2,rfac3,if2d,if1d,if3d
!//!
!//      INTEGER(4) :: mp1(6),mp2(6),mp3(6),mp4(6),mp5(6)
!//      DATA mp1/2,3,4,5,6,1/
!//      DATA mp2/3,4,5,6,1,2/
!//      DATA mp3/4,5,6,1,2,3/
!//      DATA mp4/5,6,1,2,3,4/
!//      DATA mp5/6,1,2,3,4,5/
!//      REAL(8) :: rt3, sqrt3, rsqrt3, hside
!//      REAL(8) :: hexarea
!//!      integer(4) :: ndivhs
!//
!//
!//#ifdef siarhei_plot
!//         ! @$^ siarhei_fr
!//         current_sub = &
!//            '+#+ Entered [OneTpenMg] in Mod_TPENDrive'
!//         if(trim(last_saved).ne.trim(current_sub)) then
!//            write(678,'(A)') trim(current_sub)
!//            last_saved = current_sub
!//         end if
!//         ! =-=-=-=-=-=-=-=
!//#endif
!//
!//      sqrt3=SQRT(3.d0)
!//      rt3=1.73205080756888
!//      rsqrt3=1d0/sqrt3
!//      hside=gridsize_x(1)*rsqrt3
!//      hexarea=1.5d0*gridsize_x(1)*hside
!//
!//      skip = .FALSE.
!//!
!//      IF(.NOT.ALLOCATED(tem1)) THEN
!//         onethree=1/3.0D0
!//         onenine=onethree*onethree
!//         onesix=0.5*onethree
!//         one12=0.5*onesix
!//         oneten=0.1D0
!//         rsz=1./sz
!//         rsz2=rsz*rsz
!//         rfac1=1/540.0D0
!//         rfac2=1/3240.0D0
!//         rfac3=1/360.0D0
!//         areavol=2.*rt3*onenine*rsz
!//      
!//         if1d=.FALSE.
!//         if2d=.FALSE.
!//         if3d=.FALSE.
!//         IF(nassy.EQ.1) if1d=.TRUE.
!//         IF(nz.EQ.1) if2d=.TRUE.
!//         IF(.NOT.if1d.AND..NOT.if2d) if3d=.TRUE.
!//         ALLOCATE(tem1(n_group),tem2(n_group),hflxold(n_group) &
!//              ,          asrcbd(n_group,6),xsrcbd(n_group,6),ysrcbd(n_group,6) &
!//              ,          raxy(n_group),ca1(n_group),csfi2(n_group),rptj(n_group) &
!//              ,          cbd1(n_group),cbd2(n_group),cbd3(n_group),rbd(n_group),rsfsm(n_group) &
!//              ,          aax(n_group),czj1(n_group),czj2(n_group) &
!//              ,          ctl1i(n_group),cmom1(n_group),cmom2(n_group) &
!//              ,          ccpt1(n_group),ccpt2(n_group),ccsf(9,n_group) &
!//              ,          pflcbd(n_group),sflebd(6,n_group),ssrcbbd(6,n_group) &
!//              ,          aflxbd(6,n_group),xmombd(6,n_group),ymombd(6,n_group) &
!//              ,          sj1bd(n_group),sj0bd(n_group),csj(n_group) &
!//              ,          smm1bd(n_group),smm2bd(n_group),csmm1(n_group),csmm2(n_group) &
!//              ,          csmm2hfl(n_group),csmm2sjs(n_group),csmm1sj(n_group),xax(n_group) &
!//              ,          cntim(6,n_group),pflbm(6,n_group) &
!//              ,          scgr(n_group-jupstt+1),ccgr(n_group-jupstt+1,n_group-jupstt+1) &
!//              ,          rcgr(n_group-jupstt+1,n_group-jupstt+1),scgrb(n_group-jupstt+1) &
!//              ,          cntzisum(n_group),cntisum(n_group),pflbsum(n_group))
!//         tem1=0
!//         tem2=0
!//         hflxold=0
!//         raxy=0
!//         ca1=0
!//         csfi2=0
!//         rptj=0
!//         cbd1=0
!//         cbd2=0
!//         cbd3=0
!//         rbd=0
!//         rsfsm=0
!//         aax=0
!//         czj1=0
!//         czj2=0
!//         ctl1i=0
!//         cmom1=0
!//         cmom2=0
!//         ccpt1=0
!//         ccpt2=0
!//         pflcbd=0
!//         sj1bd=0
!//         sj0bd=0
!//         csj=0
!//         smm1bd=0
!//         smm2bd=0
!//         csmm1=0
!//         csmm2=0
!//         csmm2hfl=0
!//         csmm2sjs=0
!//         csmm1sj=0
!//         xax=0
!//         scgr=0
!//         scgrb=0
!//         cntzisum=0
!//         cntisum=0
!//         pflbsum=0
!//         asrcbd=0
!//         xsrcbd=0
!//         ysrcbd=0
!//         aflxbd=0
!//         xmombd=0
!//         ymombd=0
!//         cntim=0
!//         pflbm=0
!//         ccsf=0
!//         sflebd=0
!//         ssrcbbd=0
!//         ccgr=0
!//         rcgr=0
!//!         first=.FALSE.
!//      ENDIF
!//      rhz=1./hz
!//      rhz2=rhz*rhz
!//      rhzl=1./hzl
!//      rhzu=1./hzu
!//      faccsf=-190.*onenine*rhz
!//      faccpt=-5.*onethree*rhz
!//      faccbd=-80.*onenine*rhz
!//      facmx=rhz/54.
!//      
!//      DO ii=1,6
!//         DO j=1,n_group
!//            cntim(ii,j)=cntit(j,ii)
!//            pflbm(ii,j)=pflbt(j,ii)
!//            asrcbd(j,ii)=0
!//            xsrcbd(j,ii)=0
!//            ysrcbd(j,ii)=0
!//         ENDDO
!//      ENDDO
!//      DO j=1,n_group
!//         hflxold(j)=hflx(j)
!//         tem1(j)=seff(j)
!//         tem2(j)=seff(j)
!//         cntzisum(j)=cntz0i(j)+cntz1i(j)
!//         cntisum(j)=cntim(1,j)+cntim(2,j)+cntim(3,j)+cntim(4,j) &
!//              +cntim(5,j)+cntim(6,j)
!//         pflbsum(j)=pflbm(1,j)+pflbm(2,j)+pflbm(3,j)+pflbm(4,j) &
!//              +pflbm(5,j)+pflbm(6,j)
!//      ENDDO
!//     
!//      DO j=1,n_group
!//         IF(.NOT.if1d) THEN
!//            bt=xsd(j)*rsz2
!//            gam=rt3*sz/xsd(j)
!//            raxy(j)=1./(xst(j)+80.*bt)
!//            ca1(j)=32.*bt*raxy(j)
!//            rsf=1./(40.*ca1(j)-24.)
!//            csf1=5.*ca1(j)*rsf
!//            rpt=1./(5.*ca1(j)-6.)
!//            cpt1=-5.*onesix*ca1(j)*rpt
!//            cpt3=-3.*cpt1
!//            cpt2=-cpt3+rpt
!//            rbd(j)=4./(80.*ca1(j)-48.-gam)
!//            cbd1(j)=5.*ca1(j)*rbd(j)
!//            cbd2(j)=rbd(j)-cbd1(j)
!//            DO ii=1,6
!//               ssrcbbd(ii,j)=rbd(j)*(-0.5*(pflbm(ii,j)+pflbm(mp1(ii),j)) &
!//                    -gam*cntim(ii,j))
!//            ENDDO
!//            ttt=onesix*ca1(j)
!//            DO ii=1,6
!//               tdt3(ii)=ttt*pflbm(ii,j)
!//            ENDDO
!//            DO ii=1,6
!//               i1=mp1(ii)
!//               aflxbd(ii,j)=-tdt3(ii)-tdt3(i1)
!//               xmombd(ii,j)=onesix*(tdt3(ii)+tdt3(i1))
!//               ymombd(ii,j)=0.5*(tdt3(ii)-tdt3(i1))
!//            ENDDO
!//      !   delete the boundary surfaces
!//            fac1=csf1*cbd1(j)
!//            rsfi=1./(1.-2.*fac1)
!//            csfi1=rsfi*(csf1-fac1)
!//            csfi2(j)=rsfi*(rsf-2.*csf1*cbd2(j))
!//            rpti=1./(1.-6.*cpt2*cbd2(j))
!//            cpti1=rpti*(cpt3-2.*cpt2*cbd1(j))
!//      !   delele the inner surfaces
!//            fac1=onesix/(-1.+csfi1)
!//            fac2=onesix/(1.+csfi1)
!//            fac3=onesix/(-1.+2.*csfi1)
!//            fac4=onesix/(1.+2.*csfi1)
!//            rsi1=-2.*(fac1-fac2)-fac3+fac4
!//            rsi2=fac1+fac2+fac3+fac4
!//            rsi3=fac1-fac2-fac3+fac4
!//            rsi4=-2.*(fac1+fac2)+fac3+fac4
!//            rsfsm(j)=rsi1+2.*(rsi2+rsi3)+rsi4
!//            rptj(j)=1./(1.-6.*cpti1*rsfsm(j)*csfi2(j))
!//      !
!//            t1=cpti1*rsfsm(j)*rsfi
!//            t2=rbd(j)*(rpti*cpt2-2.*t1*csf1)
!//            t3=t1*rsf
!//            ccpt1(j)=10.*(t2+2.*t3)
!//            ccpt2(j)=15.*(rpti*rpt+4.*(t2-t3))
!//            ccpt3=-rpti*cpt1-t1*0.25-3.*t3+t2
!//            ccpt4=gam*t2
!//            pflcbd(j)=ccpt3*pflbsum(j)+ccpt4*cntisum(j)
!//      !
!//            t1=rsfi*rsi1
!//            t2=rsfi*rsi2
!//            t3=rsfi*rsi3
!//            t4=rsfi*rsi4
!//            tt=csf1*rbd(j)
!//            q4=-rsf+tt
!//            q1=10.*q4
!//            ccsf(1,j)=(t1+t2)*q1
!//            ccsf(2,j)=(t2+t3)*q1
!//            ccsf(3,j)=(t3+t4)*q1
!//            q2=60.*tt
!//            q3=-30.*rsf
!//            q7=q2-q3
!//            ccsf(4,j)=(t1+t2)*q7
!//            ccsf(5,j)=(t2+t3)*q7
!//            ccsf(6,j)=(t3+t4)*q7
!//            ccsf(7,j)=(t1-t2)*q3
!//            ccsf(8,j)=(t2-t3)*q3
!//            ccsf(9,j)=(t3-t4)*q3
!//            q5=2.*csf1-rsf+q4
!//            ccsf11=q4*t1+q5*t2
!//            ccsf12=0.5*q5*(t1+t3)+q4*t2
!//            ccsf13=0.5*q5*(t2+t4)+q4*t3
!//            ccsf14=q5*t3+q4*t4
!//            tt=tt*gam
!//            ccsf15=tt*(t1+t2)
!//            ccsf16=tt*(t2+t3)
!//            ccsf17=tt*(t3+t4)
!//
!//            DO ii=1,3
!//               datj1=cntim(ii,j)+cntim(mp5(ii),j)
!//               datj2=ccsf16*(cntim(mp1(ii),j)+cntim(mp4(ii),j))
!//               datj3=cntim(mp2(ii),j)+cntim(mp3(ii),j)
!//               datc1=pflbm(mp1(ii),j)+pflbm(mp5(ii),j)
!//               datc2=pflbm(mp2(ii),j)+pflbm(mp4(ii),j)
!//               sflebd(ii,j)=ccsf11*pflbm(ii,j)+ccsf12*datc1 &
!//                    +ccsf13*datc2+ccsf14*pflbm(mp3(ii),j) &
!//                    +ccsf15*datj1+datj2+ccsf17*datj3
!//               sflebd(mp3(ii),j)=ccsf11*pflbm(mp3(ii),j)+ccsf12*datc2 &
!//                    +ccsf13*datc1+ccsf14*pflbm(ii,j) &
!//                    +ccsf15*datj3+datj2+ccsf17*datj1
!//            ENDDO
!//         ENDIF
!//!
!//         IF(if3d) THEN
!//            csf4=faccsf*raxy(j)*rsf
!//            cbd3(j)=faccbd*raxy(j)*rbd(j)
!//            csfi3=rsfi*(csf4-2.*csf1*cbd3(j))
!//            cpt4=rpti*(faccpt*raxy(j)*rpt-6.*cpt2*cbd3(j))
!//            cpt4=rptj(j)*(cpt4-6.*cpti1*rsfsm(j)*csfi3)
!//            csf4=rsfsm(j)*(csfi3-csfi2(j)*cpt4)
!//            cbd3(j)=cbd3(j)-cbd2(j)*cpt4-2.*cbd1(j)*csf4
!//            aax(j)=rhz*raxy(j)+ca1(j)*(2.*csf4+cbd3(j)-onesix*cpt4)
!//            aax(j)=0.5*aax(j)
!//            xax(j)=-onesix*ca1(j)*(cbd3(j)-csf4-onethree*cpt4) &
!//                 +facmx*raxy(j)
!//         ENDIF
!//!
!//         IF(if1d) THEN
!//            aax(j)=0.5*rhz/xst(j)
!//         ENDIF
!//         IF(.NOT.if2d) THEN
!//            gamz=hz/xsd(j)
!//            fac1=10.*aax(j)+8.+0.25*gamz
!//            fac2=10.*aax(j)+2.
!//            rdet=1./(fac1*fac1-fac2*fac2)
!//            rzj1=rdet*fac1
!//            rzj2=rdet*fac2
!//            sj1bd(j)=gamz*(rzj1*cntz1i(j)-rzj2*cntz0i(j))
!//            sj0bd(j)=gamz*(rzj1*cntz0i(j)-rzj2*cntz1i(j))
!//            csj(j)=10.*(rzj1-rzj2)
!//            czj1(j)=rzj1+rzj2
!//            czj2(j)=rzj1-rzj2
!//         ENDIF
!//         IF(if1d) THEN
!//            ctl1i(j)=0
!//            ctl2i=0
!//            smm1bd(j)=0
!//            smm2bd(j)=0
!//         ELSE
!//            fac2=105.*areavol*cbd3(j)*czj2(j)
!//            rbtmb=1./(rhzl+rhz)
!//            rbtmu=1./(rhzu+rhz)
!//            tl0bd=dtll(j)*rhzl*rbtmb
!//            tl1bd=dtlu(j)*rhzu*rbtmu
!//            smm1bd(j)=-onesix*(tl1bd-tl0bd)
!//            smm2bd(j)=oneten*(tl1bd+tl0bd)
!//            csmm1(j)=-onesix*rhz*(rbtmu-rbtmb)
!//            csmm2(j)=oneten*(-2.+rhz*(rbtmu+rbtmb))
!//            ctl1i(j)=-csmm1(j)*fac2
!//            ctl2i=-csmm2(j)*fac2
!//         ENDIF
!//         IF(.NOT.if2d) THEN
!//            btz=xsd(j)*rhz2
!//            dat=14.*btz*(1.+2.*aax(j))
!//            cmom2(j)=xst(j)+140.*btz-70.*dat*czj2(j)+ctl2i
!//            cmom1(j)=xst(j)+60.*btz*(1.-5.*czj1(j))
!//            csmm2hfl(j)=28.*btz
!//            csmm2sjs(j)=-14.*btz*(1.+2.*aax(j))
!//            csmm1sj(j)=10.*btz
!//         ENDIF
!//      ENDDO
!//!
!//      sumfs=0.
!//      sumcnti=0
!//      sumcntzi=0
!//      DO j=1,n_group
!//         sumfs=sumfs+xsnf(j)*hflx(j)
!//      ENDDO
!//      sumfs=rxkeff*sumfs
!//      DO j=1,n_group-jupstt+1
!//         scgrb(j)=0
!//      ENDDO
!//      IF(.NOT.if1d) THEN
!//         DO j=1,n_group
!//            tem1(j)=tem1(j)+rhz*cntzisum(j)
!//            sumcnti=sumcnti+cntisum(j)
!//         ENDDO
!//         DO j1=jupstt,n_group
!//            j=j1-jupstt+1
!//            scgrb(j)=scgrb(j)+areavol*cntisum(j1)
!//         ENDDO
!//         DO ii=1,6
!//            tfsrc=0.
!//            tfsrcx=0.
!//            tfsrcy=0.
!//            DO j=1,n_group
!//               tfsrc=tfsrc+xsnf(j)*aflxt(j,ii)
!//               tfsrcx=tfsrcx+xsnf(j)*xmomt(j,ii)
!//               tfsrcy=tfsrcy+xsnf(j)*ymomt(j,ii)
!//            ENDDO
!//            tfsrc=rxkeff*tfsrc
!//            tfsrcx=rxkeff*tfsrcx
!//            tfsrcy=rxkeff*tfsrcy
!//            DO j=1,n_group
!//               ttt1=srczn(j,mp1(ii))+srczn(j,mp5(ii))
!//               ttt2=srczn(j,mp2(ii))+srczn(j,mp4(ii))
!//               ttt3=srczn(j,mp1(ii))-srczn(j,mp5(ii))
!//               ttt4=srczn(j,mp2(ii))-srczn(j,mp4(ii))
!//               asrcbd(j,ii)=chi(j)*tfsrc+tem1(j)+(83.*srczn(j,ii)+17.*ttt1 &
!//                    -37.*ttt2-43.*srczn(j,mp3(ii)))*rfac1
!//               xsrcbd(j,ii)=chi(j)*tfsrcx+(-60.*tem1(j)+59.*srczn(j,ii) &
!//                    +14.*ttt1-10.*ttt2-7.*srczn(j,mp3(ii)))*rfac2
!//               ysrcbd(j,ii)=chi(j)*tfsrcy-(9.*ttt3+ttt4)*rfac3
!//            ENDDO
!//         ENDDO
!//      ENDIF
!//!
!//      IF(.NOT.if2d) THEN
!//         DO j=1,n_group
!//            tem2(j)=tem2(j)+areavol*cntisum(j)
!//            sumcntzi=sumcntzi+cntzisum(j)
!//         ENDDO
!//         DO j1=jupstt,n_group
!//            j=j1-jupstt+1
!//            scgrb(j)=scgrb(j)+rhz*cntzisum(j1)
!//         ENDDO
!//         tfsrcz1=0.
!//         tfsrcz2=0.
!//         DO j=1,n_group
!//            tfsrcz1=tfsrcz1+xsnf(j)*zmomt1(j)
!//            tfsrcz2=tfsrcz2+xsnf(j)*zmomt2(j)
!//         ENDDO
!//         tfsrcz1=rxkeff*tfsrcz1
!//         tfsrcz2=rxkeff*tfsrcz2
!//         tfsrcz=sumfs
!//         DO j=1,n_group
!//            smm1bd(j)=smm1bd(j)+chi(j)*tfsrcz1
!//            smm2bd(j)=smm2bd(j)+chi(j)*tfsrcz2
!//         ENDDO
!//      ENDIF
!//!
!//      nitrsrc=3
!//      DO iter=1,nitrsrc
!//         IF(iter.EQ.1) THEN
!//            jstt=1
!//         ELSE
!//            jstt=jupstt
!//         ENDIF
!//         DO j=jstt,n_group
!//            DO ii=1,6
!//               asrc(ii)=asrcbd(j,ii)
!//               xsrc(ii)=xsrcbd(j,ii)
!//               ysrc(ii)=ysrcbd(j,ii)
!//               DO j2=iscatib(j),iscatie(j)
!//                  asrc(ii)=asrc(ii)+xss(j2,j)*aflxt(j2,ii)
!//                  xsrc(ii)=xsrc(ii)+xss(j2,j)*xmomt(j2,ii)
!//                  ysrc(ii)=ysrc(ii)+xss(j2,j)*ymomt(j2,ii)
!//               ENDDO
!//               asrc(ii)=raxy(j)*asrc(ii)
!//               xsrc(ii)=raxy(j)*xsrc(ii)
!//               ysrc(ii)=raxy(j)*ysrc(ii)
!//            ENDDO
!//!
!//            pflc=rptj(j)*(pflcbd(j) &
!//                 +ccpt1(j)*(asrc(1)+asrc(2)+asrc(3)+asrc(4)+asrc(5)+asrc(6)) &
!//                 +ccpt2(j)*(xsrc(1)+xsrc(2)+xsrc(3)+xsrc(4)+xsrc(5)+xsrc(6)))
!//            !! calculate surface flux
!//            dat0=rsfsm(j)*csfi2(j)*pflc
!//            DO ii=1,3
!//               i3=mp3(ii)
!//               ad1=asrc(ii)+asrc(mp5(ii))
!//               ad2=ccsf(2,j)*(asrc(mp1(ii))+asrc(mp4(ii)))
!//               ad3=asrc(mp2(ii))+asrc(i3)
!//               xd1=xsrc(ii)+xsrc(mp5(ii))
!//               xd2=ccsf(5,j)*(xsrc(mp1(ii))+xsrc(mp4(ii)))
!//               xd3=xsrc(mp2(ii))+xsrc(i3)
!//               yd1=ysrc(ii)-ysrc(mp5(ii))
!//               yd2=ccsf(8,j)*(ysrc(mp1(ii))-ysrc(mp4(ii)))
!//               yd3=ysrc(mp2(ii))-ysrc(i3)
!//               sfle(ii)=ccsf(1,j)*ad1+ad2+ccsf(3,j)*ad3 &
!//                    +ccsf(4,j)*xd1+xd2+ccsf(6,j)*xd3 &
!//                    +ccsf(7,j)*yd1+yd2+ccsf(9,j)*yd3-dat0+sflebd(ii,j)
!//               sfle(i3)=ccsf(1,j)*ad3+ad2+ccsf(3,j)*ad1 &
!//                    +ccsf(4,j)*xd3+xd2+ccsf(6,j)*xd1 &
!//                    -ccsf(7,j)*yd3-yd2-ccsf(9,j)*yd1-dat0+sflebd(i3,j)
!//            ENDDO
!//#ifdef siarhei_ppr
!//            ! Siarhei pin power reconstruction
!//            if (pin_pow_hex_needed.and.coord_ppr.gt.0) then
!//               ! triangle inner boundary flux | alloc(6,n_group,nassy,nz)
!//               do ii = 1,6
!//                  inner_b_flx_h(ii,j,coord_ppr,coord_iz) = sfle(ii)
!//               end do
!//            end if
!//            ! -=-=-=-=-=-=-=-=-
!//#endif
!//
!//!  calculate boundary surface flux
!//            dat0=cbd2(j)*pflc
!//            DO ii=1,6
!//               sflb(ii)=ssrcbbd(ii,j)-10.*rbd(j)*(asrc(ii)+6.*xsrc(ii)) &
!//                    -cbd1(j)*(sfle(ii)+sfle(mp1(ii)))-dat0
!//               cnto(ii)=0.5*sflb(ii)-cntim(ii,j)
!//            ENDDO
!//#ifdef siarhei_ppr
!//            ! Siarhei pin power reconstruction
!//            if (pin_pow_hex_needed.and.coord_ppr.gt.0) then
!//               ! outer surface boundary flux | alloc(6,n_group,nassy,nz)
!//               do ii = 1,6
!//                  outer_b_flx_h(ii,j,coord_ppr,coord_iz) = sflb(ii)
!//               end do
!//            end if
!//            ! -=-=-=-=-=-=-=-=-
!//#endif
!//
!//            !!  calculate node average flux for each small triangle
!//            hflxr=0
!//            DO ii=1,6
!//               tdt1(ii)=ca1(j)*sfle(ii)
!//               tdt2(ii)=ca1(j)*sflb(ii)
!//            ENDDO
!//            dat1=onesix*ca1(j)*pflc
!//            DO ii=1,6
!//               i1=mp1(ii)
!//               aflx(ii)=aflxbd(ii,j) &
!//                    +asrc(ii)+tdt1(ii)+tdt1(i1)+tdt2(ii)-dat1
!//               xmom(ii)=xmombd(ii,j)+xsrc(ii) &
!//                    +one12*(2.*tdt2(ii)-tdt1(ii)-tdt1(i1)-4.*dat1)
!//               ymom(ii)=ymombd(ii,j)+ysrc(ii)-0.25*(tdt1(i1)-tdt1(ii))
!//               hflxr=hflxr+aflx(ii)
!//            ENDDO
!//            hflxr=onesix*hflxr
!//            IF(if2d) THEN
!//               hflx(j)=hflxr
!//               DO ii=1,6
!//                  aflxt(j,ii)=aflx(ii)
!//                  xmomt(j,ii)=xmom(ii)
!//                  ymomt(j,ii)=ymom(ii)
!//                  cntot(j,ii)=cnto(ii)
!//               ENDDO
!//               CYCLE
!//            ENDIF
!//            !!
!//            !!cc  calculate coefficients for axial nem                  ccc
!//            !!cc     phi_h=aax*(J_b+J_t)+hsrcax                         ccc
!//            !!cc     sum of phi_bd=csf4*(J_b+J_t)+bsrcax                ccc
!//            !!
!//            dat=2.*aax(j)*cntzisum(j)
!//            hflxr=hflxr+dat
!//            DO ii=1,6
!//               aflx(ii)=aflx(ii)+dat
!//            ENDDO
!//
!//            smm1=smm1bd(j)
!//            smm2=smm2bd(j)
!//            DO j2=iscatib(j),iscatie(j)
!//               smm1=smm1+xss(j2,j)*zmomt1(j2)
!//               smm2=smm2+xss(j2,j)*zmomt2(j2)
!//            ENDDO
!//            IF(if1d) THEN
!//               sbal=chi(j)*tfsrcz+seff(j)
!//               DO j2=iscatib(j),iscatie(j)
!//                  sbal=sbal+xss(j2,j)*hflx(j2)
!//               ENDDO
!//               hflxr=sbal/xst(j)+4.*aax(j)*cntzisum(j)
!//            ENDIF
!//            !! calculation for boundary surface flux
!//            ttt=csj(j)*hflxr
!//            sj1=sj1bd(j)+ttt
!//            sj0=sj0bd(j)+ttt
!//            sjsum=sj0+sj1
!//            !! Transverse Leakage and coeff's at interface
!//            IF(.NOT.if1d) THEN
!//               dtlm=areavol*(cnto(1)+cnto(2)+cnto(3) &
!//                    +cnto(4)+cnto(5)+cnto(6))-tem2(j) &
!//                    -1.5*areavol*cbd3(j)*(sjsum-2.*cntzisum(j))
!//               smm1=smm1+csmm1(j)*dtlm
!//               smm2=smm2+csmm2(j)*dtlm
!//            ENDIF
!//            !!  z-moments
!//            zmomt2(j)=(smm2+csmm2sjs(j)*sjsum+csmm2hfl(j)*hflxr)/cmom2(j)
!//            zmomt1(j)=(smm1+csmm1sj(j)*(sj1-sj0)-ctl1i(j)*zmomt2(j))/cmom1(j)
!//            !!  calculate outgoing axial currents
!//            fac1=15.*czj1(j)*zmomt1(j)
!//            fac2=35.*czj2(j)*zmomt2(j)
!//            phit=fac1-fac2+sj1
!//            phib=-fac1-fac2+sj0
!//            cntz1o(j)=0.5*phit-cntz1i(j)
!//            cntz0o(j)=0.5*phib-cntz0i(j)
!//            flxzsum=phit+phib
!//            cntzosum=cntz1o(j)+cntz0o(j)
!//            !!  calculate hex-avg. flux
!//            hflx(j)=hflxr-aax(j)*flxzsum
!//            IF(if1d) CYCLE
!//            !!  calculate outgoing radial currents
!//            fac1=0.5*cbd3(j)*cntzosum
!//            dat=aax(j)*flxzsum
!//            dat2=xax(j)*cntzosum
!//            DO ii=1,6
!//               cntot(j,ii)=cnto(ii)-fac1
!//               aflxt(j,ii)=aflx(ii)-dat
!//               xmomt(j,ii)=xmom(ii)+dat2
!//               ymomt(j,ii)=ymom(ii)
!//            ENDDO
!//!            write(*,*) 'hflxr 4:', hflx(j)
!//!            write(*,*) 'aflx 4: ', (aflxt(j,1)+aflxt(j,2)+aflxt(j,3)&
!//!                           +aflxt(j,4)+aflxt(j,5)+aflxt(j,6))/6
!//
!//#ifdef siarhei_ppr
!//            ! Siarhei pin power reconstruction
!//            if (pin_pow_hex_needed.and.coord_ppr.gt.0) then
!//               ! fluxes for calculating expansion coefficients | alloc(6,n_group,nassy,nz)
!//               do ii = 1,6
!//                  Fl_momx(ii,j,coord_ppr,coord_iz) = xmomt(j,ii)
!//                  Fl_momy(ii,j,coord_ppr,coord_iz) = ymomt(j,ii)
!//                  Fl_avg(ii,j,coord_ppr,coord_iz) = aflxt(j,ii)
!//               end do
!//            end if
!//            ! -=-=-=-=-=-=-=-=-
!//#endif
!//
!//         ENDDO
!//         !!   check residual
!//         sumcnto=0
!//         IF(.NOT.if1d) THEN
!//            DO it=1,6
!//               DO j=1,n_group
!//                  sumcnto=sumcnto+cntot(j,it)
!//               ENDDO
!//            ENDDO
!//         ENDIF
!//         sumcntzo=0
!//         sumscat=0
!//         sumsrc=0
!//         sumtot=0
!//         DO j=1,n_group
!//            sumcntzo=sumcntzo+cntz0o(j)+cntz1o(j)
!//            sumsrc=sumsrc+seff(j)
!//            sumtot=sumtot+xst(j)*hflx(j)
!//            DO j2=iscatib(j),iscatie(j)
!//               sumscat=sumscat+xss(j2,j)*hflx(j2)
!//            ENDDO
!//         ENDDO
!//         sumrhs=sumsrc+sumfs+areavol*sumcnti+rhz*sumcntzi
!//         sumlhs=sumtot-sumscat+areavol*sumcnto+rhz*sumcntzo
!//         errimbal=ABS((sumrhs-sumlhs)/sumtot)
!//         IF(errimbal.LT.0.05*epsl2) THEN
!//            skip = .TRUE.
!//            EXIT
!//         END IF
!//         !!  coarse group acceleration
!//         DO j1=jupstt,n_group
!//            j=j1-jupstt+1
!//            scgr(j)=scgrb(j)+seff(j1)+chi(j1)*sumfs
!//            DO j2=iscatib(j1),jupstt-1
!//               scgr(j)=scgr(j)+xss(j2,j1)*hflx(j2)
!//            ENDDO
!//            DO j2=jupstt,n_group
!//               ccgr(j2-jupstt+1,j)=-xss(j2,j1)*hflx(j2)
!//            ENDDO
!//            ccgr(j,j)=xst(j1)*hflx(j1)+ccgr(j,j)
!//         ENDDO
!//         IF(.NOT.if1d) THEN
!//            DO j1=jupstt,n_group
!//               j=j1-jupstt+1
!//               ccgr(j,j)=ccgr(j,j) &
!//                    +areavol*(cntot(j1,1)+cntot(j1,2)+cntot(j1,3) &
!//                    +cntot(j1,4)+cntot(j1,5)+cntot(j1,6))
!//            ENDDO
!//         ENDIF
!//         IF(.NOT.if2d) THEN
!//            DO j1=jupstt,n_group
!//               j=j1-jupstt+1
!//               ccgr(j,j)=ccgr(j,j)+rhz*(cntz0o(j1)+cntz1o(j1))
!//            ENDDO
!//         ENDIF
!//         CALL Invag(n_group-jupstt+1,ccgr,rcgr)
!//         DO j=1,n_group-jupstt+1
!//            fcgr=0
!//            DO j2=1,n_group-jupstt+1
!//               fcgr=fcgr+rcgr(j2,j)*scgr(j2)
!//            ENDDO
!//            j3=j+jupstt-1
!//            hflx(j3)=fcgr*hflx(j3)
!//            DO ii=1,6
!//               cntot(j3,ii)=fcgr*cntot(j3,ii)
!//               aflxt(j3,ii)=fcgr*aflxt(j3,ii)
!//               xmomt(j3,ii)=fcgr*xmomt(j3,ii)
!//               ymomt(j3,ii)=fcgr*ymomt(j3,ii)
!//            ENDDO
!//         ENDDO
!//      ENDDO
!//      IF(.NOT.skip) iter=iter-1
!//
!//      RETURN
!//      END SUBROUTINE OneTpenMg


      SUBROUTINE OneTpenMg(nassy,nz,jupstt,iscatib,iscatie &
     ,    xst,xsnf,xss,xsd,chi,sz,hz,hzu,hzl,pflbt,cntit &
     ,    cntz0i,cntz1i,seff,srczn, dtlu,dtll,rxkeff, aflxt &
     ,    xmomt,ymomt,zmomt1,zmomt2,cntot,cntz0o,cntz1o,hflx)

      IMPLICIT NONE

      REAL(8) :: rxkeff, sz, hz, hzu, hzl
      REAL(8) :: xst(n_group),xsnf(n_group),xsd(n_group),chi(n_group), &
                 xss(n_group,n_group)
      REAL(8) :: pflbt(n_group,6),cntit(n_group,6),cntz0i(n_group)
      REAL(8) :: cntz1i(n_group),seff(n_group),srczn(n_group,6), &
                 dtlu(n_group),dtll(n_group)
      REAL(8) :: cntot(n_group,6),cntz0o(n_group),cntz1o(n_group)
      REAL(8) :: hflx(n_group),aflxt(n_group,6),xmomt(n_group,6),ymomt(n_group,6)
      REAL(8) :: zmomt1(n_group),zmomt2(n_group)
      REAL(8) :: asrc(6),xsrc(6),ysrc(6)
      REAL(8) ::  sflb(6),sfle(6), tdt1(6),tdt2(6),tdt3(6)
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: tem1,tem2,hflxold
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: raxy,ca1,csfi2,rptj
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: cbd1,cbd2,cbd3,rbd,rsfsm
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: aax,czj1,czj2,ctl1i,cmom1
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: cmom2,ccpt1,ccpt2
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: pflcbd,sj1bd,sj0bd,csj
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: smm1bd,smm2bd,csmm1,csmm2
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: csmm2hfl,csmm2sjs,csmm1sj
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: xax,scgr,scgrb, cntzisum,cntisum,pflbsum
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:,:) :: asrcbd,xsrcbd,ysrcbd
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:,:) :: aflxbd,xmombd,ymombd,cntim
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:,:) :: pflbm,ccsf,sflebd,ssrcbbd
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:,:) :: ccgr,rcgr
      REAL(8) :: cnto(6),aflx(6),xmom(6),ymom(6)
      REAL(8) :: onethree, onenine, onesix, one12, oneten, rsz
      REAL(8) :: rsz2, rfac1, rfac2, rfac3, areavol, rhz, rhz2, rhzl
      REAL(8) :: rhzu, faccsf, faccpt, faccbd, facmx, bt, gam
      REAL(8) :: rsf, csf1, rpt, cpt1, cpt2, cpt3, ttt, fac1, rsfi
      REAL(8) :: csfi1, rpti, cpti1,fac2, fac3, fac4, rsi1, rsi2
      REAL(8) :: rsi3, rsi4, t1, t2, t3, t4, ccpt3, ccpt4, q1, q2, q4
      REAL(8) :: tt, q3, q5, q7, ccsf11, ccsf12, ccsf13, ccsf14, ccsf15
      REAL(8) :: ccsf16, ccsf17, datj1, datj2, datj3, datc1, datc2, csf4
      REAL(8) :: csfi3, cpt4, gamz, rdet, rzj1, rzj2, ctl2i, rbtmb, rbtmu
      REAL(8) :: tl0bd, tl1bd, btz, dat, sumfs, sumcnti, sumcntzi, tfsrc
      REAL(8) :: tfsrcx, tfsrcy, ttt1, ttt2, ttt3, ttt4, tfsrcz1, tfsrcz2
      REAL(8) :: dat0, ad1, ad2, ad3, pflc, tfsrcz, xd1, xd2, xd3, yd1
      REAL(8) :: yd2, yd3, hflxr, dat1, smm1, smm2, sbal, sj1, sj0, sjsum
      REAL(8) :: dtlm, phit, phib, flxzsum, cntzosum, dat2, sumcnto
      REAL(8) :: sumscat, sumsrc, sumtot, sumcntzo, sumrhs, sumlhs
      REAL(8) :: errimbal, fcgr
      INTEGER(4) :: nassy, nz, jupstt, ii, j, i1, j1, iter, jstt, j2, i3
      INTEGER(4) :: iscatib(n_group),iscatie(n_group), nitrsrc, it, j3
      LOGICAL(1) ::  if2d,if1d,if3d, skip
      SAVE onenine,onesix,one12,oneten,onethree,rsz,rsz2,areavol &
           ,    rfac1,rfac2,rfac3,if2d,if1d,if3d
!
      INTEGER(4) :: mp1(6),mp2(6),mp3(6),mp4(6),mp5(6)
      DATA mp1/2,3,4,5,6,1/
      DATA mp2/3,4,5,6,1,2/
      DATA mp3/4,5,6,1,2,3/
      DATA mp4/5,6,1,2,3,4/
      DATA mp5/6,1,2,3,4,5/
      REAL(8) :: rt3, sqrt3, rsqrt3, hside
      REAL(8) :: hexarea
#ifdef jr_vver
      REAL(8) :: epsl2=1d-5 !! Inc_Control.f90
#endif

      sqrt3=SQRT(3.d0)
      rt3=1.73205080756888
      rsqrt3=1d0/sqrt3
      hside=gridsize_x(1)*rsqrt3
      hexarea=1.5d0*gridsize_x(1)*hside

      skip = .FALSE.
!
      IF(.NOT.ALLOCATED(tem1)) THEN
         onethree=1/3.0D0
         onenine=onethree*onethree
         onesix=0.5*onethree
         one12=0.5*onesix
         oneten=0.1D0
         rsz=1./sz
         rsz2=rsz*rsz
         rfac1=1/540.0D0
         rfac2=1/3240.0D0
         rfac3=1/360.0D0
         areavol=2.*rt3*onenine*rsz
      
         if1d=.FALSE.
         if2d=.FALSE.
         if3d=.FALSE.
         IF(nassy.EQ.1) if1d=.TRUE.
         IF(nz.EQ.1) if2d=.TRUE.
         IF(.NOT.if1d.AND..NOT.if2d) if3d=.TRUE.
         ALLOCATE(tem1(n_group),tem2(n_group),hflxold(n_group) &
              ,          asrcbd(n_group,6),xsrcbd(n_group,6),ysrcbd(n_group,6) &
              ,          raxy(n_group),ca1(n_group),csfi2(n_group),rptj(n_group) &
              ,          cbd1(n_group),cbd2(n_group),cbd3(n_group),rbd(n_group),rsfsm(n_group) &
              ,          aax(n_group),czj1(n_group),czj2(n_group) &
              ,          ctl1i(n_group),cmom1(n_group),cmom2(n_group) &
              ,          ccpt1(n_group),ccpt2(n_group),ccsf(9,n_group) &
              ,          pflcbd(n_group),sflebd(6,n_group),ssrcbbd(6,n_group) &
              ,          aflxbd(6,n_group),xmombd(6,n_group),ymombd(6,n_group) &
              ,          sj1bd(n_group),sj0bd(n_group),csj(n_group) &
              ,          smm1bd(n_group),smm2bd(n_group),csmm1(n_group),csmm2(n_group) &
              ,          csmm2hfl(n_group),csmm2sjs(n_group),csmm1sj(n_group),xax(n_group) &
              ,          cntim(6,n_group),pflbm(6,n_group) &
              ,          scgr(n_group-jupstt+1),ccgr(n_group-jupstt+1,n_group-jupstt+1) &
              ,          rcgr(n_group-jupstt+1,n_group-jupstt+1),scgrb(n_group-jupstt+1) &
              ,          cntzisum(n_group),cntisum(n_group),pflbsum(n_group))
         tem1=0
         tem2=0
         hflxold=0
         raxy=0
         ca1=0
         csfi2=0
         rptj=0
         cbd1=0
         cbd2=0
         cbd3=0
         rbd=0
         rsfsm=0
         aax=0
         czj1=0
         czj2=0
         ctl1i=0
         cmom1=0
         cmom2=0
         ccpt1=0
         ccpt2=0
         pflcbd=0
         sj1bd=0
         sj0bd=0
         csj=0
         smm1bd=0
         smm2bd=0
         csmm1=0
         csmm2=0
         csmm2hfl=0
         csmm2sjs=0
         csmm1sj=0
         xax=0
         scgr=0
         scgrb=0
         cntzisum=0
         cntisum=0
         pflbsum=0
         asrcbd=0
         xsrcbd=0
         ysrcbd=0
         aflxbd=0
         xmombd=0
         ymombd=0
         cntim=0
         pflbm=0
         ccsf=0
         sflebd=0
         ssrcbbd=0
         ccgr=0
         rcgr=0
!         first=.FALSE.
      ENDIF
      rhz=1./hz
      rhz2=rhz*rhz
      rhzl=1./hzl
      rhzu=1./hzu
      faccsf=-190.*onenine*rhz
      faccpt=-5.*onethree*rhz
      faccbd=-80.*onenine*rhz
      facmx=rhz/54.
      
      DO ii=1,6
         DO j=1,n_group
            cntim(ii,j)=cntit(j,ii)
            pflbm(ii,j)=pflbt(j,ii)
            asrcbd(j,ii)=0
            xsrcbd(j,ii)=0
            ysrcbd(j,ii)=0
         ENDDO
      ENDDO
      DO j=1,n_group
         hflxold(j)=hflx(j)
         tem1(j)=seff(j)
         tem2(j)=seff(j)
         cntzisum(j)=cntz0i(j)+cntz1i(j)
         cntisum(j)=cntim(1,j)+cntim(2,j)+cntim(3,j)+cntim(4,j) &
              +cntim(5,j)+cntim(6,j)
         pflbsum(j)=pflbm(1,j)+pflbm(2,j)+pflbm(3,j)+pflbm(4,j) &
              +pflbm(5,j)+pflbm(6,j)
      ENDDO
     
      DO j=1,n_group
         IF(.NOT.if1d) THEN
            bt=xsd(j)*rsz2
            gam=rt3*sz/xsd(j)
            raxy(j)=1./(xst(j)+80.*bt)
            ca1(j)=32.*bt*raxy(j)
            rsf=1./(40.*ca1(j)-24.)
            csf1=5.*ca1(j)*rsf
            rpt=1./(5.*ca1(j)-6.)
            cpt1=-5.*onesix*ca1(j)*rpt
            cpt3=-3.*cpt1
            cpt2=-cpt3+rpt
            rbd(j)=4./(80.*ca1(j)-48.-gam)
            cbd1(j)=5.*ca1(j)*rbd(j)
            cbd2(j)=rbd(j)-cbd1(j)
            DO ii=1,6
               ssrcbbd(ii,j)=rbd(j)*(-0.5*(pflbm(ii,j)+pflbm(mp1(ii),j)) &
                    -gam*cntim(ii,j))
            ENDDO
            ttt=onesix*ca1(j)
            DO ii=1,6
               tdt3(ii)=ttt*pflbm(ii,j)
            ENDDO
            DO ii=1,6
               i1=mp1(ii)
               aflxbd(ii,j)=-tdt3(ii)-tdt3(i1)
               xmombd(ii,j)=onesix*(tdt3(ii)+tdt3(i1))
               ymombd(ii,j)=0.5*(tdt3(ii)-tdt3(i1))
            ENDDO
      !   delete the boundary surfaces
            fac1=csf1*cbd1(j)
            rsfi=1./(1.-2.*fac1)
            csfi1=rsfi*(csf1-fac1)
            csfi2(j)=rsfi*(rsf-2.*csf1*cbd2(j))
            rpti=1./(1.-6.*cpt2*cbd2(j))
            cpti1=rpti*(cpt3-2.*cpt2*cbd1(j))
      !   delele the inner surfaces
            fac1=onesix/(-1.+csfi1)
            fac2=onesix/(1.+csfi1)
            fac3=onesix/(-1.+2.*csfi1)
            fac4=onesix/(1.+2.*csfi1)
            rsi1=-2.*(fac1-fac2)-fac3+fac4
            rsi2=fac1+fac2+fac3+fac4
            rsi3=fac1-fac2-fac3+fac4
            rsi4=-2.*(fac1+fac2)+fac3+fac4
            rsfsm(j)=rsi1+2.*(rsi2+rsi3)+rsi4
            rptj(j)=1./(1.-6.*cpti1*rsfsm(j)*csfi2(j))
      !
            t1=cpti1*rsfsm(j)*rsfi
            t2=rbd(j)*(rpti*cpt2-2.*t1*csf1)
            t3=t1*rsf
            ccpt1(j)=10.*(t2+2.*t3)
            ccpt2(j)=15.*(rpti*rpt+4.*(t2-t3))
            ccpt3=-rpti*cpt1-t1*0.25-3.*t3+t2
            ccpt4=gam*t2
            pflcbd(j)=ccpt3*pflbsum(j)+ccpt4*cntisum(j)
      !
            t1=rsfi*rsi1
            t2=rsfi*rsi2
            t3=rsfi*rsi3
            t4=rsfi*rsi4
            tt=csf1*rbd(j)
            q4=-rsf+tt
            q1=10.*q4
            ccsf(1,j)=(t1+t2)*q1
            ccsf(2,j)=(t2+t3)*q1
            ccsf(3,j)=(t3+t4)*q1
            q2=60.*tt
            q3=-30.*rsf
            q7=q2-q3
            ccsf(4,j)=(t1+t2)*q7
            ccsf(5,j)=(t2+t3)*q7
            ccsf(6,j)=(t3+t4)*q7
            ccsf(7,j)=(t1-t2)*q3
            ccsf(8,j)=(t2-t3)*q3
            ccsf(9,j)=(t3-t4)*q3
            q5=2.*csf1-rsf+q4
            ccsf11=q4*t1+q5*t2
            ccsf12=0.5*q5*(t1+t3)+q4*t2
            ccsf13=0.5*q5*(t2+t4)+q4*t3
            ccsf14=q5*t3+q4*t4
            tt=tt*gam
            ccsf15=tt*(t1+t2)
            ccsf16=tt*(t2+t3)
            ccsf17=tt*(t3+t4)

            DO ii=1,3
               datj1=cntim(ii,j)+cntim(mp5(ii),j)
               datj2=ccsf16*(cntim(mp1(ii),j)+cntim(mp4(ii),j))
               datj3=cntim(mp2(ii),j)+cntim(mp3(ii),j)
               datc1=pflbm(mp1(ii),j)+pflbm(mp5(ii),j)
               datc2=pflbm(mp2(ii),j)+pflbm(mp4(ii),j)
               sflebd(ii,j)=ccsf11*pflbm(ii,j)+ccsf12*datc1 &
                    +ccsf13*datc2+ccsf14*pflbm(mp3(ii),j) &
                    +ccsf15*datj1+datj2+ccsf17*datj3
               sflebd(mp3(ii),j)=ccsf11*pflbm(mp3(ii),j)+ccsf12*datc2 &
                    +ccsf13*datc1+ccsf14*pflbm(ii,j) &
                    +ccsf15*datj3+datj2+ccsf17*datj1
            ENDDO
         ENDIF
!
         IF(if3d) THEN
            csf4=faccsf*raxy(j)*rsf
            cbd3(j)=faccbd*raxy(j)*rbd(j)
            csfi3=rsfi*(csf4-2.*csf1*cbd3(j))
            cpt4=rpti*(faccpt*raxy(j)*rpt-6.*cpt2*cbd3(j))
            cpt4=rptj(j)*(cpt4-6.*cpti1*rsfsm(j)*csfi3)
            csf4=rsfsm(j)*(csfi3-csfi2(j)*cpt4)
            cbd3(j)=cbd3(j)-cbd2(j)*cpt4-2.*cbd1(j)*csf4
            aax(j)=rhz*raxy(j)+ca1(j)*(2.*csf4+cbd3(j)-onesix*cpt4)
            aax(j)=0.5*aax(j)
            xax(j)=-onesix*ca1(j)*(cbd3(j)-csf4-onethree*cpt4) &
                 +facmx*raxy(j)
         ENDIF
!
         IF(if1d) THEN
            aax(j)=0.5*rhz/xst(j)
         ENDIF
         IF(.NOT.if2d) THEN
            gamz=hz/xsd(j)
            fac1=10.*aax(j)+8.+0.25*gamz
            fac2=10.*aax(j)+2.
            rdet=1./(fac1*fac1-fac2*fac2)
            rzj1=rdet*fac1
            rzj2=rdet*fac2
            sj1bd(j)=gamz*(rzj1*cntz1i(j)-rzj2*cntz0i(j))
            sj0bd(j)=gamz*(rzj1*cntz0i(j)-rzj2*cntz1i(j))
            csj(j)=10.*(rzj1-rzj2)
            czj1(j)=rzj1+rzj2
            czj2(j)=rzj1-rzj2
         ENDIF
         IF(if1d) THEN
            ctl1i(j)=0
            ctl2i=0
            smm1bd(j)=0
            smm2bd(j)=0
         ELSE
            fac2=105.*areavol*cbd3(j)*czj2(j)
            rbtmb=1./(rhzl+rhz)
            rbtmu=1./(rhzu+rhz)
            tl0bd=dtll(j)*rhzl*rbtmb
            tl1bd=dtlu(j)*rhzu*rbtmu
            smm1bd(j)=-onesix*(tl1bd-tl0bd)
            smm2bd(j)=oneten*(tl1bd+tl0bd)
            csmm1(j)=-onesix*rhz*(rbtmu-rbtmb)
            csmm2(j)=oneten*(-2.+rhz*(rbtmu+rbtmb))
            ctl1i(j)=-csmm1(j)*fac2
            ctl2i=-csmm2(j)*fac2
         ENDIF
         IF(.NOT.if2d) THEN
            btz=xsd(j)*rhz2
            dat=14.*btz*(1.+2.*aax(j))
            cmom2(j)=xst(j)+140.*btz-70.*dat*czj2(j)+ctl2i
            cmom1(j)=xst(j)+60.*btz*(1.-5.*czj1(j))
            csmm2hfl(j)=28.*btz
            csmm2sjs(j)=-14.*btz*(1.+2.*aax(j))
            csmm1sj(j)=10.*btz
         ENDIF
      ENDDO
!
      sumfs=0.
      sumcnti=0
      sumcntzi=0
      DO j=1,n_group
         sumfs=sumfs+xsnf(j)*hflx(j)
      ENDDO
      sumfs=rxkeff*sumfs
      DO j=1,n_group-jupstt+1
         scgrb(j)=0
      ENDDO
      IF(.NOT.if1d) THEN
         DO j=1,n_group
            tem1(j)=tem1(j)+rhz*cntzisum(j)
            sumcnti=sumcnti+cntisum(j)
         ENDDO
         DO j1=jupstt,n_group
            j=j1-jupstt+1
            scgrb(j)=scgrb(j)+areavol*cntisum(j1)
         ENDDO
         DO ii=1,6
            tfsrc=0.
            tfsrcx=0.
            tfsrcy=0.
            DO j=1,n_group
               tfsrc=tfsrc+xsnf(j)*aflxt(j,ii)
               tfsrcx=tfsrcx+xsnf(j)*xmomt(j,ii)
               tfsrcy=tfsrcy+xsnf(j)*ymomt(j,ii)
            ENDDO
            tfsrc=rxkeff*tfsrc
            tfsrcx=rxkeff*tfsrcx
            tfsrcy=rxkeff*tfsrcy
            DO j=1,n_group
               ttt1=srczn(j,mp1(ii))+srczn(j,mp5(ii))
               ttt2=srczn(j,mp2(ii))+srczn(j,mp4(ii))
               ttt3=srczn(j,mp1(ii))-srczn(j,mp5(ii))
               ttt4=srczn(j,mp2(ii))-srczn(j,mp4(ii))
               asrcbd(j,ii)=chi(j)*tfsrc+tem1(j)+(83.*srczn(j,ii)+17.*ttt1 &
                    -37.*ttt2-43.*srczn(j,mp3(ii)))*rfac1
               xsrcbd(j,ii)=chi(j)*tfsrcx+(-60.*tem1(j)+59.*srczn(j,ii) &
                    +14.*ttt1-10.*ttt2-7.*srczn(j,mp3(ii)))*rfac2
               ysrcbd(j,ii)=chi(j)*tfsrcy-(9.*ttt3+ttt4)*rfac3
            ENDDO
         ENDDO
      ENDIF
!
      IF(.NOT.if2d) THEN
         DO j=1,n_group
            tem2(j)=tem2(j)+areavol*cntisum(j)
            sumcntzi=sumcntzi+cntzisum(j)
         ENDDO
         DO j1=jupstt,n_group
            j=j1-jupstt+1
            scgrb(j)=scgrb(j)+rhz*cntzisum(j1)
         ENDDO
         tfsrcz1=0.
         tfsrcz2=0.
         DO j=1,n_group
            tfsrcz1=tfsrcz1+xsnf(j)*zmomt1(j)
            tfsrcz2=tfsrcz2+xsnf(j)*zmomt2(j)
         ENDDO
         tfsrcz1=rxkeff*tfsrcz1
         tfsrcz2=rxkeff*tfsrcz2
         tfsrcz=sumfs
         DO j=1,n_group
            smm1bd(j)=smm1bd(j)+chi(j)*tfsrcz1
            smm2bd(j)=smm2bd(j)+chi(j)*tfsrcz2
         ENDDO
      ENDIF
!
      nitrsrc=3
      DO iter=1,nitrsrc
         IF(iter.EQ.1) THEN
            jstt=1
         ELSE
            jstt=jupstt
         ENDIF
         DO j=jstt,n_group
            DO ii=1,6
               asrc(ii)=asrcbd(j,ii)
               xsrc(ii)=xsrcbd(j,ii)
               ysrc(ii)=ysrcbd(j,ii)
               DO j2=iscatib(j),iscatie(j)
                  asrc(ii)=asrc(ii)+xss(j2,j)*aflxt(j2,ii)
                  xsrc(ii)=xsrc(ii)+xss(j2,j)*xmomt(j2,ii)
                  ysrc(ii)=ysrc(ii)+xss(j2,j)*ymomt(j2,ii)
               ENDDO
               asrc(ii)=raxy(j)*asrc(ii)
               xsrc(ii)=raxy(j)*xsrc(ii)
               ysrc(ii)=raxy(j)*ysrc(ii)
            ENDDO
!
            pflc=rptj(j)*(pflcbd(j) &
                 +ccpt1(j)*(asrc(1)+asrc(2)+asrc(3)+asrc(4)+asrc(5)+asrc(6)) &
                 +ccpt2(j)*(xsrc(1)+xsrc(2)+xsrc(3)+xsrc(4)+xsrc(5)+xsrc(6)))
            !! calculate surface flux
            dat0=rsfsm(j)*csfi2(j)*pflc
            DO ii=1,3
               i3=mp3(ii)
               ad1=asrc(ii)+asrc(mp5(ii))
               ad2=ccsf(2,j)*(asrc(mp1(ii))+asrc(mp4(ii)))
               ad3=asrc(mp2(ii))+asrc(i3)
               xd1=xsrc(ii)+xsrc(mp5(ii))
               xd2=ccsf(5,j)*(xsrc(mp1(ii))+xsrc(mp4(ii)))
               xd3=xsrc(mp2(ii))+xsrc(i3)
               yd1=ysrc(ii)-ysrc(mp5(ii))
               yd2=ccsf(8,j)*(ysrc(mp1(ii))-ysrc(mp4(ii)))
               yd3=ysrc(mp2(ii))-ysrc(i3)
               sfle(ii)=ccsf(1,j)*ad1+ad2+ccsf(3,j)*ad3 &
                    +ccsf(4,j)*xd1+xd2+ccsf(6,j)*xd3 &
                    +ccsf(7,j)*yd1+yd2+ccsf(9,j)*yd3-dat0+sflebd(ii,j)
               sfle(i3)=ccsf(1,j)*ad3+ad2+ccsf(3,j)*ad1 &
                    +ccsf(4,j)*xd3+xd2+ccsf(6,j)*xd1 &
                    -ccsf(7,j)*yd3-yd2-ccsf(9,j)*yd1-dat0+sflebd(i3,j)
            ENDDO
            !!  calculate boundary surface flux
            dat0=cbd2(j)*pflc
            DO ii=1,6
               sflb(ii)=ssrcbbd(ii,j)-10.*rbd(j)*(asrc(ii)+6.*xsrc(ii)) &
                    -cbd1(j)*(sfle(ii)+sfle(mp1(ii)))-dat0
               cnto(ii)=0.5*sflb(ii)-cntim(ii,j)
            ENDDO
            !!  calculate node average flux for each small triangle
            hflxr=0
            DO ii=1,6
               tdt1(ii)=ca1(j)*sfle(ii)
               tdt2(ii)=ca1(j)*sflb(ii)
            ENDDO
            dat1=onesix*ca1(j)*pflc
            DO ii=1,6
               i1=mp1(ii)
               aflx(ii)=aflxbd(ii,j) &
                    +asrc(ii)+tdt1(ii)+tdt1(i1)+tdt2(ii)-dat1
               xmom(ii)=xmombd(ii,j)+xsrc(ii) &
                    +one12*(2.*tdt2(ii)-tdt1(ii)-tdt1(i1)-4.*dat1)
               ymom(ii)=ymombd(ii,j)+ysrc(ii)-0.25*(tdt1(i1)-tdt1(ii))
               hflxr=hflxr+aflx(ii)
            ENDDO
            hflxr=onesix*hflxr
            IF(if2d) THEN
               hflx(j)=hflxr
               DO ii=1,6
                  aflxt(j,ii)=aflx(ii)
                  xmomt(j,ii)=xmom(ii)
                  ymomt(j,ii)=ymom(ii)
                  cntot(j,ii)=cnto(ii)
               ENDDO
               CYCLE
            ENDIF
            !!
            !!cc  calculate coefficients for axial nem                  ccc
            !!cc     phi_h=aax*(J_b+J_t)+hsrcax                         ccc
            !!cc     sum of phi_bd=csf4*(J_b+J_t)+bsrcax                ccc
            !!
            dat=2.*aax(j)*cntzisum(j)
            hflxr=hflxr+dat
            DO ii=1,6
               aflx(ii)=aflx(ii)+dat
            ENDDO

            smm1=smm1bd(j)
            smm2=smm2bd(j)
            DO j2=iscatib(j),iscatie(j)
               smm1=smm1+xss(j2,j)*zmomt1(j2)
               smm2=smm2+xss(j2,j)*zmomt2(j2)
            ENDDO
            IF(if1d) THEN
               sbal=chi(j)*tfsrcz+seff(j)
               DO j2=iscatib(j),iscatie(j)
                  sbal=sbal+xss(j2,j)*hflx(j2)
               ENDDO
               hflxr=sbal/xst(j)+4.*aax(j)*cntzisum(j)
            ENDIF
            !! calculation for boundary surface flux
            ttt=csj(j)*hflxr
            sj1=sj1bd(j)+ttt
            sj0=sj0bd(j)+ttt
            sjsum=sj0+sj1
            !! Transverse Leakage and coeff's at interface
            IF(.NOT.if1d) THEN
               dtlm=areavol*(cnto(1)+cnto(2)+cnto(3) &
                    +cnto(4)+cnto(5)+cnto(6))-tem2(j) &
                    -1.5*areavol*cbd3(j)*(sjsum-2.*cntzisum(j))
               smm1=smm1+csmm1(j)*dtlm
               smm2=smm2+csmm2(j)*dtlm
            ENDIF
            !!  z-moments
            zmomt2(j)=(smm2+csmm2sjs(j)*sjsum+csmm2hfl(j)*hflxr)/cmom2(j)
            zmomt1(j)=(smm1+csmm1sj(j)*(sj1-sj0)-ctl1i(j)*zmomt2(j))/cmom1(j)
            !!  calculate outgoing axial currents
            fac1=15.*czj1(j)*zmomt1(j)
            fac2=35.*czj2(j)*zmomt2(j)
            phit=fac1-fac2+sj1
            phib=-fac1-fac2+sj0
            cntz1o(j)=0.5*phit-cntz1i(j)
            cntz0o(j)=0.5*phib-cntz0i(j)
            flxzsum=phit+phib
            cntzosum=cntz1o(j)+cntz0o(j)
            !!  calculate hex-avg. flux
            hflx(j)=hflxr-aax(j)*flxzsum
            IF(if1d) CYCLE
            !!  calculate outgoing radial currents
            fac1=0.5*cbd3(j)*cntzosum
            dat=aax(j)*flxzsum
            dat2=xax(j)*cntzosum
            DO ii=1,6
               cntot(j,ii)=cnto(ii)-fac1
               aflxt(j,ii)=aflx(ii)-dat
               xmomt(j,ii)=xmom(ii)+dat2
               ymomt(j,ii)=ymom(ii)
            ENDDO
         ENDDO
         !!   check residual
         sumcnto=0
         IF(.NOT.if1d) THEN
            DO it=1,6
               DO j=1,n_group
                  sumcnto=sumcnto+cntot(j,it)
               ENDDO
            ENDDO
         ENDIF
         sumcntzo=0
         sumscat=0
         sumsrc=0
         sumtot=0
         DO j=1,n_group
            sumcntzo=sumcntzo+cntz0o(j)+cntz1o(j)
            sumsrc=sumsrc+seff(j)
            sumtot=sumtot+xst(j)*hflx(j)
            DO j2=iscatib(j),iscatie(j)
               sumscat=sumscat+xss(j2,j)*hflx(j2)
            ENDDO
         ENDDO
         sumrhs=sumsrc+sumfs+areavol*sumcnti+rhz*sumcntzi
         sumlhs=sumtot-sumscat+areavol*sumcnto+rhz*sumcntzo
         errimbal=ABS((sumrhs-sumlhs)/sumtot)
         IF(errimbal.LT.0.05*epsl2) THEN
            skip = .TRUE.
            EXIT
         END IF
         !!  coarse group acceleration
         DO j1=jupstt,n_group
            j=j1-jupstt+1
            scgr(j)=scgrb(j)+seff(j1)+chi(j1)*sumfs
            DO j2=iscatib(j1),jupstt-1
               scgr(j)=scgr(j)+xss(j2,j1)*hflx(j2)
            ENDDO
            DO j2=jupstt,n_group
               ccgr(j2-jupstt+1,j)=-xss(j2,j1)*hflx(j2)
            ENDDO
            ccgr(j,j)=xst(j1)*hflx(j1)+ccgr(j,j)
         ENDDO
         IF(.NOT.if1d) THEN
            DO j1=jupstt,n_group
               j=j1-jupstt+1
               ccgr(j,j)=ccgr(j,j) &
                    +areavol*(cntot(j1,1)+cntot(j1,2)+cntot(j1,3) &
                    +cntot(j1,4)+cntot(j1,5)+cntot(j1,6))
            ENDDO
         ENDIF
         IF(.NOT.if2d) THEN
            DO j1=jupstt,n_group
               j=j1-jupstt+1
               ccgr(j,j)=ccgr(j,j)+rhz*(cntz0o(j1)+cntz1o(j1))
            ENDDO
         ENDIF
         CALL Invag(n_group-jupstt+1,ccgr,rcgr)
         DO j=1,n_group-jupstt+1
            fcgr=0
            DO j2=1,n_group-jupstt+1
               fcgr=fcgr+rcgr(j2,j)*scgr(j2)
            ENDDO
            j3=j+jupstt-1
            hflx(j3)=fcgr*hflx(j3)
            DO ii=1,6
               cntot(j3,ii)=fcgr*cntot(j3,ii)
               aflxt(j3,ii)=fcgr*aflxt(j3,ii)
               xmomt(j3,ii)=fcgr*xmomt(j3,ii)
               ymomt(j3,ii)=fcgr*ymomt(j3,ii)
            ENDDO
         ENDDO
      ENDDO
      IF(.NOT.skip) iter=iter-1

      RETURN
      END SUBROUTINE OneTpenMg

      
      SUBROUTINE Invag(n,a,r)
      IMPLICIT NONE!
!  inverse all group nxn matrix!
      INTEGER(4) :: n, j, j2, j3
      REAL(8) :: a(n,n), r(n,n), t(n,n), rdet, tem
!

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Invag] in Mod_TPENDrive'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      DO j=1,n
         DO j2=1,n
            t(j2,j)=a(j2,j)
            r(j2,j)=0.
         ENDDO
         r(j,j)=1.
      ENDDO
!
      DO j=1,n

         IF (t(j,j) .EQ. 0.0) THEN
            rdet = 1.0e30
         ELSE
            rdet = 1./t(j,j)
         ENDIF

         DO j2=j+1,n
            t(j2,j)=rdet*t(j2,j)
         ENDDO
         DO j2=1,j
            r(j2,j)=rdet*r(j2,j)
         ENDDO
         DO j2=1,j-1
            tem=-t(j,j2)
            DO j3=j+1,n
               t(j3,j2)=t(j3,j2)+tem*t(j3,j)
            ENDDO
            DO j3=1,j
               r(j3,j2)=r(j3,j2)+tem*r(j3,j)
            ENDDO
         ENDDO
         DO j2=j+1,n
            tem=-t(j,j2)
            DO j3=j+1,n
               t(j3,j2)=t(j3,j2)+tem*t(j3,j)
            ENDDO
            DO j3=1,j
               r(j3,j2)=r(j3,j2)+tem*r(j3,j)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE Invag


      SUBROUTINE SolpFlxMg

#ifdef siarhei_ppr
         ! Siarhei pin power reconstruction
         use Inc_PinPow_Hex
         ! -=-=-=-=-=-=-=-=-
#endif

      IMPLICIT NONE!
! point flux solver using CPB relation!
      INTEGER(4) :: ig,ih,it,ip,iz,mp1,mp2,mp3,mp4,mp5,nn,ipt,iter,nmem
      REAL(8) :: dat1,dat2,dat3,wts1,wts2,wts3,wts4,pbdvl,phipt,bsflx,dathflx
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: psrc,pcoef,psrctem, xsdwtp(:,:)
      REAL(8) :: xsdwt
      DIMENSION bsflx(6),phipt(6)
      DIMENSION mp1(6),mp2(6),mp3(6),mp4(6),mp5(6)
      DATA mp1/2,3,4,5,6,1/
      DATA mp2/3,4,5,6,1,2/
      DATA mp3/4,5,6,1,2,3/
      DATA mp4/5,6,1,2,3,4/
      DATA mp5/6,1,2,3,4,5/
!

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SolpFlxMg] in Mod_TPENDrive'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF(.NOT.ALLOCATED(psrc)) THEN
         nmem=n_group*6
         ALLOCATE(psrc(ncorn_hex),pcoef(ncorn_hex),psrctem(ncorn_hex))
         psrc=0
         pcoef=0
         psrctem=0
         ALLOCATE(xsdwtp(3,nassy))
         xsdwtp=0
      ENDIF
!
      DO iz=1,nz
         DO ig=1,n_group
            DO ip=1,nxpnt
               pcoef(ip)=0
               psrc(ip)=0.
               pflxt(ip)=pflx(ig,ip,iz)
            ENDDO
            DO ih=1,nassy
               ! xsdf = xsd = d_3d
               xsdwt=wtass(ih)*d_3d_Mg(imap(ih),iz,ig)
               xsdwtp(1,ih)=chlval(1)*xsdwt
               xsdwtp(2,ih)=chlval(2)*xsdwt
               xsdwtp(3,ih)=chlval(3)*xsdwt
               wts1=xsdwt*chlval(4)
               wts2=xsdwt*chlval(5)
               wts3=xsdwt*chlval(6)
               wts4=xsdwt*chlval(7)
               DO it=1,6
                  nn=neignd(it,ih)
                  IF(nn.EQ.0) THEN
                     bsflx(it)=2.*(1.+reflratf(ig))*cnto(ig,it,ih,iz)
                  ELSE
                     bsflx(it)=2.*(cnto(ig,it,ih,iz) &
                          +cnto(ig,neigjin(it,ih),nn,iz))
                  ENDIF
                  if (xsadf(ig,it,ih,iz) .NE. 0) then
                      bsflx(it)=bsflx(it)/xsadf(ig,it,ih,iz)
                  else
                      bsflx(it)=bsflx(it)/1d0 !xsadf(ig,it,ih,iz)
                  endif
               ENDDO
               dathflx=wts4*hflxf(ig,ih,iz)
               DO it=1,3
                  nn=neigpt(it,ih)
                  dat1=bsflx(it)+bsflx(mp5(it))
                  dat2=bsflx(mp2(it))+bsflx(mp3(it))
                  dat3=wts3*(bsflx(mp1(it))+bsflx(mp4(it)))
                  pbdvl=pbdv(ig,nn)*wtass(ih)
                  psrc(nn)=psrc(nn)+wts1*dat1+wts2*dat2+dat3+dathflx
                  pcoef(nn)=pcoef(nn)+xsdwt/adfpt(ig,it,ih,iz)+pbdvl
                  nn=neigpt(mp3(it),ih)
                  pbdvl=pbdv(ig,nn)*wtass(ih)
                  psrc(nn)=psrc(nn)+wts1*dat2+wts2*dat1+dat3+dathflx
                  pcoef(nn)=pcoef(nn)+xsdwt/adfpt(ig,mp3(it),ih,iz)+pbdvl
               ENDDO
            ENDDO

            DO iter=1,3
               DO ip=1,nxpnt
                  psrctem(ip)=psrc(ip)
               ENDDO
               DO ih=1,nassy
                  DO it=1,6
                     phipt(it)=pflxt(neigpt(it,ih))/adfpt(ig,it,ih,iz)
                  ENDDO
                  DO it=1,3
                     dat1=phipt(mp1(it))+phipt(mp5(it))
                     dat2=phipt(mp2(it))+phipt(mp4(it))
                     ipt=neigpt(it,ih)
                     psrctem(ipt)=psrctem(ipt) &
                          +xsdwtp(3,ih)*phipt(mp3(it)) &
                          +xsdwtp(1,ih)*dat1+xsdwtp(2,ih)*dat2
                     ipt=neigpt(mp3(it),ih)
                     psrctem(ipt)=psrctem(ipt) &
                          +xsdwtp(3,ih)*phipt(it) &
                          +xsdwtp(1,ih)*dat2+xsdwtp(2,ih)*dat1
                  ENDDO
               ENDDO
               DO ip=1,nxpnt
                  pflxt(ip)=0.7d0*psrctem(ip)/pcoef(ip)+0.3d0*pflxt(ip)
               ENDDO
            ENDDO
            DO ip=1,nxpnt
               pflx(ig,ip,iz)=pflxt(ip)
#ifdef siarhei_ppr
               ! Siarhei pin power reconstruction
               if (pin_pow_hex_needed) &
                  corner_flux_h(ip,ig,iz) = pflxt(ip)
                  ! fluxes for calculating expansion coefficients
                  ! | alloc(6,n_group,nassy,nz)
               ! -=-=-=-=-=-=-=-=-
#endif

            ENDDO
         ENDDO
      ENDDO
!stop

      RETURN
      END SUBROUTINE SolpFlxMg

      FUNCTION residtpen_mg()

      IMPLICIT NONE
      REAL(8) :: residtpen_mg
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: cntz0i, cntz1i, hflxfn, cnti(:,:)
!      REAL(8),POINTER,SAVE,DIMENSION(:,:,:) :: hflxtmp            !(n_group,nassy,nz)
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:,:,:) :: hflxtmp
      REAL(8) :: totresid, totfsrc, arear, vol, sumlkr, sumlkz, sumss
      REAL(8) :: sumfs, sumes, sumtxs, sumlk, resid
      INTEGER(4) :: iz, ih, it, ig, nn, neigdn, neigup, ig2
      REAL(8) :: rt3, sqrt3, rsqrt3, hside
      REAL(8) :: hexarea
!      integer(4) :: ndivhs


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [residtpen_mg] in Mod_TPENDrive'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      sqrt3=SQRT(3.d0)
      rt3=1.73205080756888
      rsqrt3=1d0/sqrt3
      hside=gridsize_x(1)*rsqrt3
      hexarea=1.5d0*gridsize_x(1)*hside
!
      IF(.NOT.ALLOCATED(cnti)) THEN
          ALLOCATE(cnti(n_group,ntph),cntz0i(n_group),cntz1i(n_group),hflxfn(n_group))
          allocate(hflxtmp(n_group,nassy,nz))
          cnti=0
          cntz0i=0
          cntz1i=0
          hflxfn=0
          ! multigroup = .true. (when, .not.FMFD)
          !IF(multigroup) THEN
          hflxtmp=0d0
      ENDIF
      do ig =1,n_group
         do ih = 1,nassy
            do iz = 1, nz
               hflxtmp(ig,ih,iz)=hflxf(ig,ih,iz)
            enddo
         enddo
      enddo
          !ELSE
          !   hflxtmp=>hflx
          !ENDIF

      totresid=0
      totfsrc=0

      DO iz=1,nz
         arear=hside*hz(iz)
         DO ih=1,nassy
            IF(nassy.GT.1) THEN
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
            ENDIF
            IF(nz.GT.1) THEN
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
            ENDIF
            vol=nodevolume(imap(ih),iz)
            sumlkr=0
            sumlkz=0
            sumss=0
            sumfs=0
            sumes=0
            sumtxs=0
            DO ig=1,n_group
               DO it=1,6
                  sumlkr=sumlkr+cnto(ig,it,ih,iz)-cnti(ig,it)
               ENDDO
               sumlkz=sumlkz+cntzo(ig,1,ih,iz)-cntz0i(ig) &
                    +cntzo(ig,2,ih,iz)-cntz1i(ig)
               DO ig2=iscatib(ig),ig-1
                  sumss=sumss+xssf(ig2,ig,ih,iz)*hflxtmp(ig2,ih,iz)
               ENDDO
               DO ig2=ig+1,iscatie(ig)
                  sumss=sumss+xssf(ig2-1,ig,ih,iz)*hflxtmp(ig2,ih,iz)
               ENDDO
               ! xsnff = xsnf = nu_maXs_f_3D
               sumfs=sumfs+nu_maXs_f_3D_mg(imap(ih),iz,ig)*hflxtmp(ig,ih,iz)
               sumes=sumes+sefve(ig,ih,iz)
               ! xstf = xst = maxs_r_3d
               sumtxs=sumtxs+maXs_r_3D_mg(imap(ih),iz,ig)*hflxtmp(ig,ih,iz)
            ENDDO
            sumlk=(sumlkr*arear+sumlkz*hexarea)*wtass(ih)
            resid=(sumfs*reigv+sumss+sumes-sumtxs)*vol-sumlk
            totresid=totresid+resid*resid
            sumfs=sumfs*vol
            totfsrc=totfsrc+sumfs*sumfs
         ENDDO
      ENDDO
      residtpen_mg=SQRT(totresid/totfsrc)
!
      RETURN
      END FUNCTION residtpen_mg

!#ifdef siarhei_delete 
      SUBROUTINE TpenBc

      IMPLICIT NONE
! Update Nodal Solution from CMFD results
      INTEGER(4) :: k,l,m
      REAL(8) :: sumf,fnorm
      REAL(8) :: vol,sump
      DIMENSION sumf(2),sump(2),fnorm(2)
      REAL(8) :: rt3, sqrt3, rsqrt3, hside
!      integer(4) :: ndivhs

      sqrt3=SQRT(3.d0)
      rt3=1.73205080756888
      rsqrt3=1d0/sqrt3
      hside=gridsize_x(1)*rsqrt3
!
      DO m=1,2
         sumf(m)=0
         sump(m)=0
      ENDDO
      DO k=1,nz
         DO l=1,nxy
            vol=nodevolume(imap(l),k)
            DO m=1,2
               sumf(m)=sumf(m)+hflx(m,l,k)*vol
               sump(m)=sump(m)+flux(imap(l),k,m)*vol
            ENDDO
         ENDDO
      ENDDO
      !stop
      DO m=1,2
         fnorm(m)=sump(m)/sumf(m)
      ENDDO
 
! multigroup = .true. (= .not.FMFD)
!      IF(multigroup) THEN
      CALL TpenBcMg(fnorm)
!      ELSE

      RETURN
      END SUBROUTINE TpenBc
!#endif 

      SUBROUTINE TpenBcMg(fnorm)
      USE Inc_Transient, Only: Flag_Transient
      USE Inc_Solver, only: betap, rvdelt
#ifdef tuan_tr_test
      use inc_tpen, only: rvdeltf
#endif
      
      IMPLICIT NONE
!
! Update Nodal Solution from CMFD results
!
      INTEGER(4) :: k,l,m,nn,iz2
      INTEGER(4) :: ig,ih,ip,it,iz,id1,id2,ig2,iz1,igf,ind,isfc
      REAL(8)    :: bt,cjn,hzr,flx1,flx2,rvol,dfdm,dhatm,betamw,cntoavg,cntoavg1,cntoavg2
      REAL(8)    :: w,wf,ttt,tt1,tt2,ttt2,phis,rrt3h,onewf,oneflxm,oneflxn,fnorm
      DIMENSION fnorm(2)
      LOGICAL ifbcref(2)
      REAL(8)    :: rt3, sqrt3, rsqrt3, hside
!      integer(4) :: ndivhs

      sqrt3=SQRT(3.d0)
      rt3=1.73205080756888
      rsqrt3=1d0/sqrt3
      hside=gridsize_x(1)*rsqrt3
!
!
      DO iz=1,nz
         DO ih=1,nassy
            DO ig=1,2
               hflx(ig,ih,iz)=flux(imap(ih),iz,ig)
            ENDDO
         ENDDO
         DO ip=1,nxpnt
            DO ig=1,2
               DO igf=mgb(ig),mge(ig)
                  pflx(igf,ip,iz)=pflx(igf,ip,iz)*fnorm(ig)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      DO ig=1,2
         ifbcref(ig)=.FALSE.
         IF(alxr(ig).EQ.0) ifbcref(ig)=.TRUE.
      ENDDO
!
      IF(iptype.EQ.0) THEN
         wf=1.15D0
      ELSE
         wf=1
      ENDIF

!
      onewf=1.-wf
      rrt3h=1./rt3/hside
      DO iz=1,nz
         DO ih=1,nassy
            DO ig=1,2
               DO ig2=mgb(ig),mge(ig)
                  ttt=wf*fhflx(ig2,ih,iz)+(1.-wf)*fohflx(ig2,ih,iz)
                  ttt2=ttt*hflx(ig,ih,iz)/hflxf(ig2,ih,iz)
                  hflxf(ig2,ih,iz)=ttt*hflx(ig,ih,iz)
                  DO it=1,6
                     aflx(ig2,it,ih,iz)=aflx(ig2,it,ih,iz)*ttt2
                     xmom(ig2,it,ih,iz)=xmom(ig2,it,ih,iz)*ttt2
                     ymom(ig2,it,ih,iz)=ymom(ig2,it,ih,iz)*ttt2
                  ENDDO
                  zmom1(ig2,ih,iz)=zmom1(ig2,ih,iz)*ttt2
                  zmom2(ig2,ih,iz)=zmom2(ig2,ih,iz)*ttt2
               ENDDO
            ENDDO
         ENDDO
      ENDDO
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
                  oneflxm=hflx(ig,ih,iz)
                  IF(nn.NE.0) THEN
                     w=dfdm/d_3d(imap(ih),iz,ig)/2.
                     oneflxn=hflx(ig,nn,iz)
                     phis=w*oneflxn+(1.-w)*oneflxm+0.5*bt*(oneflxm+oneflxn)
                     cjn=-((dfdm+dhatm)*oneflxn-(dfdm-dhatm)*oneflxm) &
                          /rt3/hside
                     cntoavg=0.25*(phis+2.*cjn)
                  ELSE
                     cjn=(dfdm-dhatm)*oneflxm/rt3/hside
                     IF(bt.EQ.0) THEN
                        phis=oneflxm/(1+0.5*gridsize_x(1)*alxr(ig)/d_3d(imap(ih),iz,ig)) ! y
                     ELSE
                        phis=bt*oneflxm
                     ENDIF
                     cntoavg=0.25*(phis+2.*cjn)
                  ENDIF
                  DO ig2=mgb(ig),mge(ig)
                     ttt=wf*fcnto(ig2,it,ih,iz)+(1.-wf)*focnto(ig2,it,ih,iz)
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
               hzr=hz(iz2)/hz(iz1)
               DO ih=1,nassy
                  DO ig=1,2
                     dfdm=dfdz(ig,ih,iz)
                     dhatm=dhatz(ig,ih,iz)
                     bt=betaphisz(ig,ih,iz)
                     flx1=hflx(ig,ih,iz1)
                     flx2=hflx(ig,ih,iz2)
                     w=dfdm/d_3d(imap(ih),iz1,ig)/(1.+hzr)
                     phis=w*flx2+(1.-w)*flx1+0.5*bt*(flx1+flx2)
                     cjn=-((dfdm+dhatm)*flx2-(dfdm-dhatm)*flx1)*2. &
                          /(hz(iz1)+hz(iz2))
                     cntoavg1=0.25*(phis+2.*cjn)
                     cntoavg2=0.25*(phis-2.*cjn)
                     DO ig2=mgb(ig),mge(ig)
                        tt1=wf*fcntzo(ig2,id1,ih,iz1) &
                             +(1.-wf)*focntzo(ig2,id1,ih,iz1)
                        tt2=wf*fcntzo(ig2,id2,ih,iz2) &
                             +(1.-wf)*focntzo(ig2,id2,ih,iz2)
                        cntzo(ig2,id1,ih,iz1)=tt1*cntoavg1
                        cntzo(ig2,id2,ih,iz2)=tt2*cntoavg2
!                     write(1995,*) 'tt1,tt2, cntzo: ',tt1, tt2, cntzo(ig2,id1,ih,iz1), cntzo(ig2,id2,ih,iz2)

                     ENDDO
!                     write(1995,*) ig, ih, iz,iz1, iz2,dfdm, dhatm,bt, flx1,flx2, cjn, phis, cntoavg1

                  ENDDO
               ENDDO
            ELSE
               iz1=neigsndz(ind,iz)
               id1=neigsndz(ind+3,iz)
               IF(iz.EQ.1) THEN
                  DO ig=1,2
                     alphaz(ig)=alzl(ig)
                  ENDDO
               ELSE
                  DO ig=1,2
                     alphaz(ig)=alzr(ig)
                  ENDDO
               ENDIF
               DO ih=1,nassy
                  DO ig=1,2
                     dfdm=dfdz(ig,ih,iz)
                     dhatm=dhatz(ig,ih,iz)
                     bt=betaphisz(ig,ih,iz)
                     IF(ind.EQ.2) dhatm=-dhatm
                     flx1=hflx(ig,ih,iz1)
                     cjn=(dfdm-dhatm)*flx1/hz(iz1)
                     IF(bt.EQ.0) THEN
                        phis=flx1/(1+0.25*hz(iz1)/d_3d(imap(ih),iz1,ig)*alphaz(ig))
                     ELSE
                        phis=bt*flx1
                     ENDIF
                     cntoavg1=0.25*(phis+2.*cjn)
                     DO ig2=mgb(ig),mge(ig)
                        tt1=wf*fcntzo(ig2,id1,ih,iz1) &
                             +(1.-wf)*focntzo(ig2,id1,ih,iz1)
                        cntzo(ig2,id1,ih,iz1)=tt1*cntoavg1
!                        write(1995,*) 'tt1, cntzo: ',tt1, cntzo(ig2,id1,ih,iz1)
                     ENDDO
!                     write(1995,*) ig, ih, iz,dfdm, dhatm,bt, flx1, cjn, phis, cntoavg1
                     
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!stop
!
      IF(.NOT.(iptype.EQ.1 .AND. Flag_Transient)) RETURN
! calculate effective trasient fixed source
   if (n_group > 2) then
      DO k=1,nz
         DO l=1,nxy
            rvol=1d0/nodevolume(imap(l),k)
            betamw=1-betap(imap(l),k)
            DO m=1,n_group
#ifdef tuan_tr_test
               !! === INDEX ===================== !!
               !! sefve  : hex index              !!
               !! src    : square index           !!
               !! rvdelt : square index           !!
               !! FisSrc : square index           !!
               !! =============================== !!
                  !! for test with parcs 
                  sefve(m,l,k)=srcf(m,imap(l),k)*rvol
                  sefve(m,l,k)=srcf(m,imap(l),k)*rvol &
                       -rvdeltf(m,imap(l),k)*hflxf(m,l,k) &
                       -betamw*maxs_chid_3d_mg(imap(l),k,m)*FisSrc(imap(l),k)*rvol
#else                
               sefve(m,l,k)=srcf(m,l,k)*rvol
               sefve(m,l,k)=srcf(m,l,k)*rvol &
                    -rvdelt(m,imap(l),k)*hflxf(m,l,k) &
                    -betamw*maxs_chid_3d_mg(imap(l),k,m)*FisSrc(imap(l),k)*rvol
#endif
            ENDDO
         ENDDO
      ENDDO
   else
      DO k=1,nz
         DO l=1,nxy
            rvol=1d0/nodevolume(imap(l),k)
            betamw=1-betap(imap(l),k)
            DO m=1,n_group
               sefve(m,l,k)=srcf(m,l,k)*rvol
               sefve(m,l,k)=srcf(m,l,k)*rvol &
                    -rvdelt(m,imap(l),k)*hflxf(m,l,k) &
                    -betamw*maxs_chid_3d_mg(imap(l),k,m)*FisSrc(imap(l),k)*rvol
            ENDDO
         ENDDO
      ENDDO
   endif
!stop
      RETURN
      END SUBROUTINE TpenBcMg


      SUBROUTINE TpenBcMg2!(fnorm)
      USE Inc_Transient, Only: Flag_Transient
      USE Inc_Solver, only: betap, rvdelt

      IMPLICIT NONE!
! Update Nodal Solution from CMFD results!
      INTEGER(4) :: k,l,m,nn,iz2
      INTEGER(4) :: ig,ih,ip,it,iz,id1,id2,ig2,iz1,igf,ind,isfc
      REAL(8)    :: bt,cjn,hzr,flx1,flx2,rvol,dfdm,dhatm,betamw,cntoavg,cntoavg1,cntoavg2
      REAL(8)    :: w,wf,ttt,tt1,tt2,ttt2,phis,rrt3h,onewf,oneflxm,oneflxn,fnorm
      DIMENSION fnorm(2)
      LOGICAL ifbcref(n_group) !!(2)
      REAL(8)    :: rt3, sqrt3, rsqrt3, hside
!      integer(4) :: ndivhs


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [TpenBcMg!] in Mod_TPENDrive'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      fnorm(:) = 1d0 !! @$^

      sqrt3=SQRT(3.d0)
      rt3=1.73205080756888
      rsqrt3=1d0/sqrt3
      hside=gridsize_x(1)*rsqrt3
!
!
      DO iz=1,nz
         DO ih=1,nassy
            DO ig = 1,n_group!!1,2
               hflx(ig,ih,iz)=flux(imap(ih),iz,ig)
            ENDDO
         ENDDO
       !!  DO ip=1,nxpnt
       !!     DO ig=1,2
       !!        DO igf=mgb(ig),mge(ig)
       !!           pflx(igf,ip,iz)=pflx(igf,ip,iz)*fnorm(ig)
       !!        ENDDO
       !!     ENDDO
       !!  ENDDO
      ENDDO
!
      DO ig = 1,n_group !!ig=1,2
         ifbcref(ig)=.FALSE.
         IF(alxr(ig).EQ.0) ifbcref(ig)=.TRUE.
      ENDDO
!
      IF(iptype.EQ.0) THEN
         wf=1.15D0
      ELSE
         wf=1
      ENDIF

!
      onewf=1.-wf
      rrt3h=1./rt3/hside
      DO iz=1,nz
         DO ih=1,nassy
            DO ig=1,2 !!checked
               DO ig2=mgb(ig),mge(ig)
                  ttt=wf*fhflx(ig2,ih,iz)+(1.-wf)*fohflx(ig2,ih,iz)
                  ttt2=ttt*hflx(ig,ih,iz)/hflxf(ig2,ih,iz)
                  hflxf(ig2,ih,iz)=ttt*hflx(ig,ih,iz)
                  DO it=1,6
                     aflx(ig2,it,ih,iz)=aflx(ig2,it,ih,iz)*ttt2
                     xmom(ig2,it,ih,iz)=xmom(ig2,it,ih,iz)*ttt2
                     ymom(ig2,it,ih,iz)=ymom(ig2,it,ih,iz)*ttt2
                  ENDDO
                  zmom1(ig2,ih,iz)=zmom1(ig2,ih,iz)*ttt2
                  zmom2(ig2,ih,iz)=zmom2(ig2,ih,iz)*ttt2
               ENDDO
            ENDDO
         ENDDO
      ENDDO
! surface flux cal. from phis=w*phin+(1-w)phim+beta*(phim+phin)/2
      DO iz=1,nz
         DO ih=1,nassy
            DO it=1,6
               nn=neignd(it,ih)
               isfc=neigsfc(it,ih)
               DO ig=1,2 !!checked
                  dfdm=dfd(ig,isfc,iz)
                  dhatm=wtdhat(it,ih)*dhat(ig,isfc,iz)
                  bt=betaphis(ig,isfc,iz)
                  oneflxm=hflx(ig,ih,iz)
                  IF(nn.NE.0) THEN
                     w=dfdm/d_3d(imap(ih),iz,ig)/2.
                     oneflxn=hflx(ig,nn,iz)
                     phis=w*oneflxn+(1.-w)*oneflxm+0.5*bt*(oneflxm+oneflxn)
                     cjn=-((dfdm+dhatm)*oneflxn-(dfdm-dhatm)*oneflxm) &
                          /rt3/hside
                     cntoavg=0.25*(phis+2.*cjn)
                  ELSE
                     cjn=(dfdm-dhatm)*oneflxm/rt3/hside
                     IF(bt.EQ.0) THEN
                        phis=oneflxm/(1+0.5*gridsize_x(1)*alxr(ig)/d_3d(imap(ih),iz,ig)) ! y
                     ELSE
                        phis=bt*oneflxm
                     ENDIF
                     cntoavg=0.25*(phis+2.*cjn)
                  ENDIF
                  DO ig2=mgb(ig),mge(ig)
                     ttt=wf*fcnto(ig2,it,ih,iz)+(1.-wf)*focnto(ig2,it,ih,iz)
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
               hzr=hz(iz2)/hz(iz1)
               DO ih=1,nassy
                  DO ig=1,2 !!checked
                     dfdm=dfdz(ig,ih,iz)
                     dhatm=dhatz(ig,ih,iz)
                     bt=betaphisz(ig,ih,iz)
                     flx1=hflx(ig,ih,iz1)
                     flx2=hflx(ig,ih,iz2)
                     w=dfdm/d_3d(imap(ih),iz1,ig)/(1.+hzr)
                     phis=w*flx2+(1.-w)*flx1+0.5*bt*(flx1+flx2)
                     cjn=-((dfdm+dhatm)*flx2-(dfdm-dhatm)*flx1)*2. &
                          /(hz(iz1)+hz(iz2))
                     cntoavg1=0.25*(phis+2.*cjn)
                     cntoavg2=0.25*(phis-2.*cjn)
                     DO ig2=mgb(ig),mge(ig)
                        tt1=wf*fcntzo(ig2,id1,ih,iz1) &
                             +(1.-wf)*focntzo(ig2,id1,ih,iz1)
                        tt2=wf*fcntzo(ig2,id2,ih,iz2) &
                             +(1.-wf)*focntzo(ig2,id2,ih,iz2)
                        cntzo(ig2,id1,ih,iz1)=tt1*cntoavg1
                        cntzo(ig2,id2,ih,iz2)=tt2*cntoavg2
!                     write(1995,*) 'tt1,tt2, cntzo: ',tt1, tt2, cntzo(ig2,id1,ih,iz1), cntzo(ig2,id2,ih,iz2)

                     ENDDO
!                     write(1995,*) ig, ih, iz,iz1, iz2,dfdm, dhatm,bt, flx1,flx2, cjn, phis, cntoavg1
                  ENDDO
               ENDDO
            ELSE
               iz1=neigsndz(ind,iz)
               id1=neigsndz(ind+3,iz)
               IF(iz.EQ.1) THEN
                  DO ig = 1,n_group !!ig=1,2
                     alphaz(ig)=alzl(ig)
                  ENDDO
               ELSE
                  DO ig = 1,n_group !!ig=1,2
                     alphaz(ig)=alzr(ig)
                  ENDDO
               ENDIF
               DO ih=1,nassy
                  DO ig=1,2 !!checked
                     dfdm=dfdz(ig,ih,iz)
                     dhatm=dhatz(ig,ih,iz)
                     bt=betaphisz(ig,ih,iz)
                     IF(ind.EQ.2) dhatm=-dhatm
                     flx1=hflx(ig,ih,iz1)
                     cjn=(dfdm-dhatm)*flx1/hz(iz1)
                     IF(bt.EQ.0) THEN
                        phis=flx1/(1+0.25*hz(iz1)/d_3d(imap(ih),iz1,ig)*alphaz(ig))
                     ELSE
                        phis=bt*flx1
                     ENDIF
                     cntoavg1=0.25*(phis+2.*cjn)
                     DO ig2=mgb(ig),mge(ig)
                        tt1=wf*fcntzo(ig2,id1,ih,iz1) &
                             +(1.-wf)*focntzo(ig2,id1,ih,iz1)
                        cntzo(ig2,id1,ih,iz1)=tt1*cntoavg1
!                        write(1995,*) 'tt1, cntzo: ',tt1, cntzo(ig2,id1,ih,iz1)

                     ENDDO
!                     write(1995,*) ig, ih, iz,dfdm, dhatm,bt, flx1, cjn, phis, cntoavg1

                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!stop
!
      IF(.NOT.(iptype.EQ.1 .AND. Flag_Transient)) RETURN
! calculate effective trasient fixed source
   if (n_group > 2) then
      DO k=1,nz
         DO l=1,nxy
            rvol=1d0/nodevolume(imap(l),k)
            betamw=1-betap(imap(l),k)
            DO m=1,n_group
#ifdef tuan_tr_test
               !! === INDEX ===================== !!
               !! sefve  : hex index              !!
               !! src    : square index           !!
               !! rvdelt : square index           !!
               !! FisSrc : square index           !!
               !! =============================== !!
                  !! for test with parcs 
                  sefve(m,l,k)=srcf(m,imap(l),k)*rvol
                  sefve(m,l,k)=srcf(m,imap(l),k)*rvol &
                       -rvdeltf(m,imap(l),k)*hflxf(m,l,k) &
                       -betamw*maxs_chid_3d_mg(imap(l),k,m)*FisSrc(imap(l),k)*rvol
#else                
               sefve(m,l,k)=srcf(m,l,k)*rvol
               sefve(m,l,k)=srcf(m,l,k)*rvol &
                    -rvdelt(m,imap(l),k)*hflxf(m,l,k) &
                    -betamw*maxs_chid_3d_mg(imap(l),k,m)*FisSrc(imap(l),k)*rvol
#endif
            ENDDO
         ENDDO
      ENDDO
   else
      DO k=1,nz
         DO l=1,nxy
            rvol=1d0/nodevolume(imap(l),k)
            betamw=1-betap(imap(l),k)
            DO m=1,n_group
               sefve(m,l,k)=srcf(m,l,k)*rvol
               sefve(m,l,k)=srcf(m,l,k)*rvol &
                    -rvdelt(m,imap(l),k)*hflxf(m,l,k) &
                    -betamw*maxs_chid_3d_mg(imap(l),k,m)*FisSrc(imap(l),k)*rvol
            ENDDO
         ENDDO
      ENDDO
   endif
!stop
      RETURN
      END SUBROUTINE TpenBcMg2


      
      
#ifdef tuan_fr

#ifdef siarhei_delete 
!    subroutine scarpinit
!        USE Inc_TPEN
!        implicit none
!        integer :: i, j, n, ind
!        allocate(imapf(nxy_hex), imapb(nxy_hex))
!
!        imapf = 0
!        imapb  = 0
!
!            n = 1
!            do j = -1,ny_hex+2
!                do i = -3,nx_hex+4
!                    ind = nodel_hex(i,j)
!                    if (ind > 0) then
!                        imapb(n) = ind
!                        n = n + 1
!                    endif
!                enddo
!            enddo
!
!            do i = 1, nxy_hex
!                n = imapb(i)
!                imapf(n) = i
!            enddo
!
!    end subroutine scarpinit
#endif 

#ifdef siarhei_delete 
    subroutine p2r(a,b)
        ! Reodering of XY nodes from SCARP to RAST-K
        USE Inc_TPEN
        use Inc_Option, only: N_Group
        implicit none
        integer :: i, j, k, l
        real(8), INTENT(out) :: a(:,:,:)
        real(8), INTENT(in) :: b(:,:,:)


        do k = 1, N_Group
            do i = 1, nassy
                l = imapb(i)
                a(i,:,k) = b(k,l,:)
            enddo
        enddo
    end subroutine p2r
#endif 

#endif


      END MODULE Mod_TPENDrive
#endif
