#ifdef jr_vver

      MODULE Mod_SolTPEN

      USE Inc_Constant
      USE Inc_Geometry
      USE Inc_3D
      USE Inc_Control
      USE Inc_maXS, ONLY: nu_maXS_f_3D
      USE Inc_Option, ONLY: N_Group
      USE Inc_TPEN
      USE Inc_Kinetics, ONLY: beta_d_tot
      USE Mod_SolLS

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif
      
      CONTAINS

!!!#ifdef siarhei_delete 
      SUBROUTINE SolveLinearSystem_hex(iin,checkfiss)

      IMPLICIT NONE

      INTEGER:: iin
      LOGICAL :: checkfiss
      
#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SolveLinearSystem_hex] in Mod_SolTPEN'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif
#ifdef tuan_tr_test
      call BicgStab_hex(iin,checkfiss)
#else
      CALL BicgStab(iin,checkfiss)
#endif
      END SUBROUTINE SolveLinearSystem_hex
!!!#endif 

      SUBROUTINE setlshex(iroute)

      USE Mod_Alloc
      USE Inc_Lscoef
      USE Inc_FluxVar
      USE Inc_RP
      USE Inc_maXS, ONLY: maXS_r_3D, nu_maXS_f_3D,      &
                          maxs_chi_3d, maxs_chid_3d,  &
                          maxs_s_3d
      USE Inc_Transient
      USE Inc_Solver
      USE Inc_Kinetics, ONLY: beta_d_tot

      
      IMPLICIT NONE

      INTEGER(4) :: i, k, l, m
      INTEGER(4) :: ih, ig, ir, it, iz, isfc, ir2, ifr, iroute
      REAL(8) :: chi1, chi2, chid, chip, tmpm(6)
      REAL(8) :: vol
      REAL(8) :: r9hs2
      REAL(8) :: rt3, sqrt3, rsqrt3, hside
      REAL(8) :: hexarea
      INTEGER(4) :: ind
!      integer(4) :: ndivhs


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [setlshex] in Mod_SolTPEN'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      sqrt3=SQRT(3.d0)
      rt3=1.73205080756888
      rsqrt3=1/sqrt3
      hside=gridsize_x(1)*rsqrt3
      hexarea=1.5d0*gridsize_x(1)*hside
!
!
      r9hs2=2./9.d0/(hside*hside)      
      IF(iroute.EQ.2) then
         Do iz = 1,nz
            Do ih =1,nxy_1n
               ind=imap(ih)
               vol=nodevolume(ind,iz) 
               chi=betap(ind,iz)*reigvs
               ! betat = beta_d_tot
               chip=(1-beta_d_tot(ind,iz))*reigvs
               chid=(betap(ind,iz)+beta_d_tot(ind,iz)-1)*reigvs
               chi1=chip*maxs_chi_3d(ind,iz,1)+chid*maxs_chid_3d(ind,iz,1)
               chi2=chip*maxs_chi_3d(ind,iz,2)+chid*maxs_chid_3d(ind,iz,2)
               dcmat(1,ih,iz)=(maXs_r_3d(ind,iz,1)-chi1*nu_maXs_f_3D(ind,iz,1) &  
                    +dsum(1,ih,iz))*vol
               dcmat(2,ih,iz)=-chi1*nu_maXs_f_3D(ind,iz,2)*vol
               dcmat(3,ih,iz)=-(maXS_s_3D(ind, Iz, 1)+chi2*nu_maXs_f_3D(ind,iz,1))*vol
               dcmat(4,ih,iz)= (maXs_r_3d(ind,iz,2)-chi2*nu_maXs_f_3D(ind,iz,2) &
                    +dsum(2,ih,iz))*vol
            ENDDO
         ENDDO
         RETURN
      ENDIF

      DO iz=1,nz
         o1: DO ih=1,nxy_1n
            DO it=1,6
               isfc=neigsfc(it,ih)
               DO ig=1,2
                  cmat(ig,it,ih,iz)=r9hs2 &
                       *(dfd(ig,isfc,iz)+wtdhat(it,ih)*dhat(ig,isfc,iz))
               ENDDO
            ENDDO
            IF(nz.EQ.1) CYCLE o1
            isfc=neigsfcz(1,iz)
            DO ig=1,2
               cmat(ig,7,ih,iz)=rhzbar2(1,iz) &
                    *(dfdz(ig,ih,isfc)-dhatz(ig,ih,iz))
            ENDDO
            isfc=neigsfcz(2,iz)
            DO ig=1,2
               cmat(ig,8,ih,iz)=rhzbar2(2,iz) &
                    *(dfdz(ig,ih,isfc)+dhatz(ig,ih,iz+1))
            ENDDO
         ENDDO o1
      ENDDO
!
! Dsum Calculation
      DO iz=1,nz
         DO ih=1,nxy_1n
            DO ig=1,2
               dsum(ig,ih,iz)=0
               DO it=1,6
                  isfc=neigsfc(it,ih)
                  dsum(ig,ih,iz)=dsum(ig,ih,iz) &
                       +dfd(ig,isfc,iz)-wtdhat(it,ih)*dhat(ig,isfc,iz)
               ENDDO
               dsum(ig,ih,iz)=dsum(ig,ih,iz)*r9hs2
               IF(nz.GT.1) THEN
                  isfc=neigsfcz(1,iz)
                  dsum(ig,ih,iz)=dsum(ig,ih,iz) &
                       +rhzbar2(1,iz)*(dfdz(ig,ih,isfc)+dhatz(ig,ih,iz))
                  isfc=neigsfcz(2,iz)
                  dsum(ig,ih,iz)=dsum(ig,ih,iz) &
                       +rhzbar2(2,iz)*(dfdz(ig,ih,isfc)-dhatz(ig,ih,iz+1))
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
! condense cmat on radial direction

      DO iz=1,nz
         DO ih=1,nxy_1n
            DO ir2=1,ineigcond(0,0,ih)
               ifr=ineigcond(ir2,0,ih)
               dsum(1,ih,iz)=dsum(1,ih,iz)-cmat(1,ifr,ih,iz)
               dsum(2,ih,iz)=dsum(2,ih,iz)-cmat(2,ifr,ih,iz)
            ENDDO
            DO ig=1,2
               DO ir=1,ipntr(0,ih)
                  tmpm(ir)=0
                  DO ir2=1,ineigcond(0,ir,ih)
                     ifr=ineigcond(ir2,ir,ih)
                     tmpm(ir)=tmpm(ir)+cmat(ig,ifr,ih,iz)
                  ENDDO
               ENDDO
               DO ir=1,ipntr(0,ih)
                  cmat(ig,ir,ih,iz)=tmpm(ir)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO iz=1,nz
         DO ih=1,nxy_1n
            ind=imap(ih)
            chi=betap(ind,iz)*reigvs
            chip=(1-beta_d_tot(ind,iz))*reigvs
            chid=(betap(ind,iz)+beta_d_tot(ind,iz)-1)*reigvs
            chi1=chip*maxs_chi_3d(ind,iz,1)+chid*maxs_chid_3d(ind,iz,1)
            chi2=chip*maxs_chi_3d(ind,iz,2)+chid*maxs_chid_3d(ind,iz,2)
            dcmat(1,ih,iz)=maxs_r_3d(ind,iz,1) &
                 -chi1*nu_maxs_f_3d(ind,iz,1)+dsum(1,ih,iz)
            dcmat(2,ih,iz)=-chi1*nu_maxs_f_3d(ind,iz,2)
            dcmat(3,ih,iz)=-chi2*nu_maxs_f_3d(ind,iz,1)-maXS_s_3D(ind, Iz, 1)
            dcmat(4,ih,iz)=maxs_r_3d(ind,iz,2)-chi2*nu_maxs_f_3d(ind,iz,2) &
                 +dsum(2,ih,iz)
         ENDDO
      ENDDO
!
      DO k=1,nz
         DO l=1,nxy
            ind=imap(l)
            vol=nodevolume(ind,k)
            af(1,l,k)=nu_maxs_f_3d(ind,k,1)*vol
            af(2,l,k)=nu_maxs_f_3d(ind,k,2)*vol
            dcmat(1,l,k)=dcmat(1,l,k)+rvdelt(1,ind,k)
            dcmat(4,l,k)=dcmat(4,l,k)+rvdelt(2,ind,k)
            DO i=1,4
               dcmat(i,l,k)=dcmat(i,l,k)*vol
            ENDDO
            DO it=1,ipntr(0,l)
               DO m=1,2
                  cmat(m,it,l,k)=cmat(m,it,l,k)*vol
               ENDDO
            ENDDO
            DO it=7,8
               DO m=1,2
                  cmat(m,it,l,k)=cmat(m,it,l,k)*vol
               ENDDO
            ENDDO            
         ENDDO
      ENDDO
      DO l=1,nxy
         DO m=1,2
            cmat(m,7,l,1)=0
            cmat(m,8,l,nz)=0
         ENDDO
      ENDDO

      l=1
      k=1
!
!801   FORMAT(2i4,3f12.4,3x,1p,2e14.6,3x,0p,8f12.4)
      RETURN
      END SUBROUTINE setlshex

      SUBROUTINE premat
      IMPLICIT NONE
      INTEGER(4) :: ii, ih, ir, ir2, imin, inext, inum, irfr, irtmp(6), irto
!
!   Generation of information for Condensation of Matrix cmat
!   ipntr(0,ia) : The number of elements in cmat(ig,nn,ia,iz)
!

! radial direction condensation

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [premat] in Mod_SolTPEN'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      DO ih=1,nassy
         ipntr(0,ih)=0
         DO ir=0,6
            ineigcond(0,ir,ih)=0
         ENDDO
         o1: DO irfr=1,6
            inum=neignd(irfr,ih)
            IF(inum.NE.0) THEN
               IF(inum.EQ.ih) THEN
                  inext=ineigcond(0,0,ih)+1
                  ineigcond(inext,0,ih)=irfr
                  ineigcond(0,0,ih)=inext
               ELSE
                  DO irto=1,ipntr(0,ih)
                     IF(inum.EQ.ipntr(irto,ih)) THEN
                        inext=ineigcond(0,irto,ih)+1
                        ineigcond(inext,irto,ih)=irfr
                        ineigcond(0,irto,ih)=inext
                        CYCLE o1
                     ENDIF
                  ENDDO
                  inext=ipntr(0,ih)+1
                  ipntr(inext,ih)=inum
                  ineigcond(0,inext,ih)=1
                  ineigcond(1,inext,ih)=irfr
                  ipntr(0,ih)=inext
               ENDIF
            ENDIF
         ENDDO o1
         DO ir=1,ipntr(0,ih)
            imin=ir
            DO ir2=ir+1,ipntr(0,ih)
               IF(ipntr(ir2,ih).LT.ipntr(imin,ih)) imin=ir2
            ENDDO
            IF(imin.NE.ir) THEN
               DO ir2=1,ineigcond(0,ir,ih)
                  irtmp(ir2)=ineigcond(ir2,ir,ih)
               ENDDO
               DO ir2=1,ineigcond(0,imin,ih)
                  ineigcond(ir2,ir,ih)=ineigcond(ir2,imin,ih)
               ENDDO
               DO ir2=1,ineigcond(0,ir,ih)
                  ineigcond(ir2,imin,ih)=irtmp(ir2)
               ENDDO
               ii=ineigcond(0,ir,ih)
               ineigcond(0,ir,ih)=ineigcond(0,imin,ih)
               ineigcond(0,imin,ih)=ii
               ii=ipntr(ir,ih)
               ipntr(ir,ih)=ipntr(imin,ih)
               ipntr(imin,ih)=ii
            ENDIF
         ENDDO
         ilubnd(ih)=0
         DO ir=1,ipntr(0,ih)
            IF(ipntr(ir,ih).LT.ih) ilubnd(ih)=ir
            iastopnt(ipntr(ir,ih),ih)=ir
         ENDDO
         iastopnt(ih,ih)=7
      ENDDO

      RETURN
      END SUBROUTINE premat

      SUBROUTINE DtilHex
      use inc_maXS, only: d_3d
      IMPLICIT NONE
!
! Calculate D_tilde for hex that is defined by ordinary FDM
      REAL(8) :: cofbd(n_group),cofbdz(n_group), dc1, dc2, hzr
      INTEGER(4) :: isfc, ind, is1, is2, iz, ig, iz1, iz2, ih
      real(8) :: rt3, sqrt3, rsqrt3, hside
!      integer(4) :: ndivhs


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [DtilHex] in Mod_SolTPEN'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      sqrt3=SQRT(3.d0)
      rt3=1.73205080756888
      rsqrt3=1/sqrt3
      hside=gridsize_x(1)*rsqrt3

      cofbd(1)=rt3*alxr(1)*hside
      cofbd(2)=rt3*alxr(2)*hside


      DO isfc=1,nxsfc
         ind=neigsnd(3,isfc)
         IF(ind.EQ.12) THEN
            is1=neigsnd(1,isfc)
            is2=neigsnd(2,isfc)
            DO iz=1,nz
               DO ig=1,2
                  dc1=D_3D(imap(is1),iz,ig)
                  dc2=D_3D(imap(is2),iz,ig)
                  dfd(ig,isfc,iz)=2.*dc1*dc2/(dc1+dc2)
               ENDDO
            ENDDO
         ELSE
            is1=neigsnd(ind,isfc)
            DO iz=1,nz
               DO ig=1,2
                  dc1=D_3D(imap(is1),iz,ig)
                  dfd(ig,isfc,iz)=dc1*2.*cofbd(ig)/(cofbd(ig)+2.*dc1)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      IF(nz.EQ.1) RETURN   !if 2D
      DO iz=1,nz+1
         ind=neigsndz(3,iz)
         IF(ind.EQ.12) THEN
            iz1=neigsndz(1,iz)
            iz2=neigsndz(2,iz)
            hzr=hz(iz2)/hz(iz1)
            DO ih=1,nassy
               DO ig=1,2
                  dc1=D_3D(imap(ih),iz1,ig)
                  dc2=D_3D(imap(ih),iz2,ig)
                  dfdz(ig,ih,iz)=dc2*(1.+hzr)/(hzr+dc2/dc1)
               ENDDO
            ENDDO
         ELSE
            iz1=neigsndz(ind,iz)
            IF(iz.EQ.1) THEN
               cofbdz(1)=hz(iz1)*alzl(1)
               cofbdz(2)=hz(iz1)*alzl(2)
            ELSE
               cofbdz(1)=hz(iz1)*alzr(1)
               cofbdz(2)=hz(iz1)*alzr(2)
            ENDIF
            DO ih=1,nassy
               DO ig=1,2
                  dc1=D_3D(imap(ih),iz1,ig)
                  dfdz(ig,ih,iz)=dc1*2.*cofbdz(ig)/(cofbdz(ig)+2.*dc1)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
!
!
      RETURN
      END SUBROUTINE DtilHex

      END MODULE Mod_SolTPEN

#endif
