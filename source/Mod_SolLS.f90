
      MODULE Mod_SolLS

      USE Inc_Constant
      USE Inc_Geometry
      USE Inc_3D
      USE Inc_Control
      USE Inc_maXS, ONLY: nu_maXS_f_3D
      USE Inc_Option, ONLY: N_Group
#ifdef jr_vver
      USE Inc_TPEN, ONLY: if_hexgometry, imap
#endif

      CONTAINS

#ifdef siarhei_delete 
      SUBROUTINE SolveLinearSystem(iin,checkfiss)


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE
      INTEGER:: iin
      LOGICAL :: checkfiss

      CALL BicgStab(iin,checkfiss)

      END SUBROUTINE SolveLinearSystem
#endif 


      SUBROUTINE BicgStab(iin,checkfiss)

      USE Inc_Lscoef
      USE Inc_FluxVar
      USE Inc_BiCG
      use Inc_ExtSrc, only: flag_ExtSrc
      IMPLICIT NONE
      ! solves a system of linear equations by preconditioned conjugate gradient method
      INTEGER, INTENT(out) :: iin
      INTEGER :: k, l, Ixy, Iz
      REAL(8) :: errl2t,r2t,err,psipsidt,errlinft
      REAL(8) :: r20epserf, crhod
      REAL(8) :: r0v, pts, ptt, vol
      LOGICAL, INTENT(in) :: checkfiss

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [BicgStab] in Mod_SolLS'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call AxBHex( Flux, flux_add )
      

      
      r20 = 0
      b2  = 0
      DO Iz = 1, Nz
         DO Ixy = 1, Nxy
            vr0(Ixy, Iz, 1) = src(1, Ixy, Iz) - flux_add(Ixy, Iz, 1)
            vr0(Ixy, Iz, 2) = src(2, Ixy, Iz) - flux_add(Ixy, Iz, 2)
            vr (Ixy, Iz, 1) = vr0(Ixy, Iz, 1)
            vr (Ixy, Iz, 2) = vr0(Ixy, Iz, 2)
            r20 = r20 + vr0(Ixy, Iz, 1)**2 + vr0(Ixy, Iz, 2)**2
            b2  = b2  + src(1, Ixy, Iz)**2 + src(2, Ixy, Iz)**2
            vp(Ixy, Iz, 1) = D0
            vp(Ixy, Iz, 2) = D0
            vv(Ixy, Iz, 1) = D0
            vv(Ixy, Iz, 2) = D0
            if(isnan(vr0(Ixy, Iz, 1))) then
                write(*,*) 'vr0 NaN 1', Ixy, Iz, src(1, Ixy, Iz),flux_add(Ixy, Iz, 1)
                stop
            endif 
            
            if(isnan(vr0(Ixy, Iz, 2))) then
                write(*,*) 'vr0 NaN 2', Ixy, Iz
                stop
            endif
            
            if(isnan(vr(Ixy, Iz, 1))) then
                write(*,*) 'vr NaN1 ', Ixy, Iz
                stop
            endif 
            
            if(isnan(vr(Ixy, Iz, 2))) then
                write(*,*) 'vr NaN2', Ixy, Iz
                stop
            endif 
            
         END DO
      END DO
      r20epserf = r20 * EPS_Residual * EPS_Residual
      r20 = SQRT(r20)
      b2  = SQRT(b2)
      calpha = 1
      crho   = 1
      comega = 1
      DO Iin = 1, Iin_Max
         crhod = crho
         crho  = scalp(vr0, vr)
         cbeta = crho*calpha / (crhod*comega)

         DO k = 1, Nz
            DO l = 1, Nxy
               vp(l, k, 1) = vr(l, k, 1) + cbeta *( vp(l, k, 1) - comega*vv(l, k, 1) )
               vp(l, k, 2) = vr(l, k, 2) + cbeta *( vp(l, k, 2) - comega*vv(l, k, 2) )
            END DO
         END DO

         CALL MInvHex(vp, vy)
         call AxBHex( vy, vv )

         r0v    = scalp(vr0, vv)
         calpha = crho / r0v

         DO k=1,Nz
            DO l=1,Nxy
               vs(l, k, 1) = vr(l, k, 1) - calpha*vv(l, k, 1)
               vs(l, k, 2) = vr(l, k, 2) - calpha*vv(l, k, 2)
            END DO
         END DO
         CALL MInvHex(vs, vz)
         call AxBHex( vz, vt )


         pts = scalp(vs, vt)
         ptt = scalp(vt, vt)
         comega = pts / ptt
         IF ( checkfiss ) THEN
            flagr2= .FALSE.
            flagl2= .FALSE.
            flaglinf= .FALSE.
            flagerf= .FALSE.
            errl2t=0
            r2t=0
            psipsidt=0
            errlinft=0
            DO k = 1, Nz
               DO l = 1,Nxy
                  vol = NodeVolume(l,k)
                  Flux(l, k, 1) = Flux(l, k, 1) + calpha*vy(l, k, 1) + comega*vz(l, k, 1)
                  Flux(l, k, 2) = Flux(l, k, 2) + calpha*vy(l, k, 2) + comega*vz(l, k, 2)
                  vr(l, k, 1) = vs(l, k, 1) - comega*vt(l, k, 1)
                  vr(l, k, 2) = vs(l, k, 2) - comega*vt(l, k, 2)
                  r2t = r2t + vr(l, k, 1)**2 + vr(l, k, 2)**2
                  IF ( nu_maXs_f_3D(l, k, 2) /= D0 ) THEN
                     FisSrc_Iout(l, k) = FisSrc(l,k)
                     FisSrc(l, k) = ( nu_maXs_f_3D(l, k, 1)*Flux(l, k, 1) + nu_maXs_f_3D(l, k, 2)*Flux(l, k, 2) ) * vol
                     err = FisSrc(l, k) - FisSrc_Iout(l, k)
                     errlinft =MAX( errlinft, ABS( err/FisSrc(l,k) ) )
                     errl2t   = errl2t + err**2
                     psipsidt = psipsidt + FisSrc(l, k)*FisSrc_Iout(l, k)
                  END IF
               END DO
            END DO
            

            
            errl2=SQRT(errl2t/ABS(psipsidt))
            errlinf=errlinft
            r2=SQRT(r2t)
            r2ob2=r2/b2
            r2=r2/r20
! write(*,*) "aaaa, convergence :  ", r2ob2,EPS_Residual, errl2,eps_global,errlinf, eps_local

            IF( r2ob2   <= EPS_Residual ) flagr2   = .TRUE.
            IF( errl2   <= EPS_Global   ) flagl2   = .TRUE.
            IF( errlinf <= EPS_Local    ) flaglinf = .TRUE.
            IF( r2 <= EPS_ERF ) THEN
               flagerf = .TRUE.
               RETURN
            END IF                  
         ELSE
            r2t = 0
            DO k = 1, Nz
               DO l = 1, Nxy
                  vol = NodeVolume(l, k)
                  Flux(l, k, 1) = Flux(l, k, 1) + calpha*vy(l, k, 1) + comega*vz(l, k, 1)
                  Flux(l, k, 2) = Flux(l, k, 2) + calpha*vy(l, k, 2) + comega*vz(l, k, 2)
                  vr(l, k, 1) = vs(l, k, 1) - comega*vt(l, k, 1)
                  vr(l, k, 2) = vs(l, k, 2) - comega*vt(l, k, 2)
                  r2t = r2t + vr(l, k, 1)**2 + vr(l, k, 2)**2
               END DO
            END DO
            IF(r2t.LT.r20epserf)RETURN
         END IF
      END DO
      Iin = Iin_Max

      ! Negative flux correction
      ! skip negative flux correction for transient analysis with external source (KHNP CRI)
      if (.not.flag_extsrc) then
         do k=1,Nz
            do l=1,Nxy
               if (Flux(l,k,1)<0d0) Flux(l,k,1)=1d-30
               if (Flux(l,k,2)<0d0) Flux(l,k,2)=1d-30
            enddo
         enddo
      endif
      
      RETURN
      END SUBROUTINE BicgStab

#ifdef siarhei_delete 
      SUBROUTINE AxB(b,ab)

      USE Inc_Lscoef

      IMPLICIT NONE
      ! Matrix-vector product
      INTEGER :: k, kp1, km1, ln, ls, Ixy
      REAL(8), INTENT(in) :: b(0:Nxy+1, Nz, N_Group)
      REAL(8), INTENT(out) :: ab(Nxy, Nz, N_Group)    !dmm

      DO k=1,Nz
         kp1=k+1
         km1=k-1
         IF(kp1.GT.Nz) kp1=Nz
         IF(km1.LT.1) km1=1
         DO Ixy = 1, Nxy
            ln = I_Ly_P2R(Ixy)
            ls = I_Ry_P2R(Ixy)
            ab(Ixy, k, 1) = am ( 1, Ixy, k ) * b( Ixy    , k  , 1 )  &
                          + ccw( 1, Ixy, k ) * b( Ixy - 1, k  , 1 )  &
                          + cce( 1, Ixy, k ) * b( Ixy + 1, k  , 1 )  &
                          + ccn( 1, Ixy, k ) * b( ln     , k  , 1 )  &
                          + ccs( 1, Ixy, k ) * b( ls     , k  , 1 )  &
                          + ccb( 1, Ixy, k ) * b( Ixy    , km1, 1 )  &
                          + cct( 1, Ixy, k ) * b( Ixy    , kp1, 1 )  &
                          - af2   ( Ixy, k ) * b( Ixy    , k  , 2 )
            ab(Ixy, k, 2) = am ( 2, Ixy, k ) * b( Ixy    , k  , 2 )  &
                          + ccw( 2, Ixy, k ) * b( Ixy - 1, k  , 2 )  &
                          + cce( 2, Ixy, k ) * b( Ixy + 1, k  , 2 )  &
                          + ccn( 2, Ixy, k ) * b( ln     , k  , 2 )  &
                          + ccs( 2, Ixy, k ) * b( ls     , k  , 2 )  &
                          + ccb( 2, Ixy, k ) * b( Ixy    , km1, 2 )  &
                          + cct( 2, Ixy, k ) * b( Ixy    , kp1, 2 )  &
                          - scat  ( Ixy, k ) * b( Ixy    , k  , 1 )
         END DO
      END DO

      END SUBROUTINE AxB
#endif 

#ifdef siarhei_delete 
      SUBROUTINE Minv(b,x)

      USE Inc_Lscoef
      USE Inc_BiCG

      IMPLICIT NONE
      ! solve Mx=b for x given b
      REAL(8) :: x(0:Nxy+1, Nz, N_Group), b(Nxy, Nz, N_Group)
      INTEGER :: l, k, kp1

      ! forward solve
      DO l=-Nx+1,Nxy+Nx !dmm
         s(1,l)=zero
         s(2,l)=zero
      END DO
      DO k=1,Nz-1
         DO l=1,Nxy
            b0(1,l)=b(l, k, 1)-ccb(1,l,k)*s(1,l)
            b0(2,l)=b(l, k, 2)-ccb(2,l,k)*s(2,l)
         END DO
!         CALL sol2d(k,b0,s)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         DO l=1,Nxy
            x(l, k, 1) = s(1, l)
            x(l, k, 2) = s(2, l)
         END DO
      END DO
      ! on the top plane
      k=Nz
      DO l=1,Nxy
         b0(1, l) = b(l, k, 1) - ccb(1,l,k)*s(1,l)
         b0(2, l) = b(l, k, 2) - ccb(2,l,k)*s(2,l)
      END DO
!      CALL sol2d(k,b0,s)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
      DO l=1,Nxy
         x(l, k, 1) = s(1, l)
         x(l, k, 2) = s(2, l)
      END DO
      ! backward
      DO k = Nz-1,1,-1
         kp1=k+1
         DO l=1,Nxy
            b0(1, l) = x(l, kp1, 1)*cct(1,l,k)
            b0(2, l) = x(l, kp1, 2)*cct(2,l,k)
         END DO
!         CALL sol2d(k,b0,s)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         DO l = 1, Nxy
            x(l, k, 1) = x(l, k, 1) - s(1, l)
            x(l, k, 2) = x(l, k, 2) - s(2, l)
         END DO
      END DO

      RETURN

      END SUBROUTINE Minv
#endif 

#ifdef siarhei_delete 
      SUBROUTINE sol2d(k,b,x)

      USE Inc_BiCG
      USE Inc_Lscoef

      IMPLICIT NONE
      ! solve a 2d problem using precalculated LU factors
      INTEGER :: i,j,k,l,jp1,ls,lout
      REAL(8) :: b(N_Group,Nxy),x(N_Group,-Nx+1:Nxy+Nx) !dmm

      !  forward solve
      j=1
      DO i=1,Nx !dmm
         s1dl(1,i)=0
         s1dl(2,i)=0
      END DO
      lout=0
      DO j=1,Ny
         l=lout
         DO i=Ix_Start_y(j),Ix_End_y(j)
            l=l+1
            b01d(1,i)=b(1,l)-ccn(1,l,k)*s1dl(1,i)
            b01d(2,i)=b(2,l)-ccn(2,l,k)*s1dl(2,i)
         END DO
         ! solve 1d problem
!         CALL sol1d(j,k,b01d,s1dl)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         l=lout
         DO i=Ix_Start_y(j),Ix_End_y(j)
            l=l+1
            x(1,l)=s1dl(1,i)
            x(2,l)=s1dl(2,i)
         END DO
         lout=lout+nrnx(j)
      END DO
      !  backward solve
      lout=Nxy-nrnx(Ny)
      jp1=Ny
      DO j=Ny-1,1,-1
         l=lout
         DO i=Ix_End_y(j),Ix_Start_y(j),-1
            ls=nodel(i,jp1)
            b01d(1,i)=x(1,ls)*ccs(1,l,k)
            b01d(2,i)=x(2,ls)*ccs(2,l,k)
            l=l-1
         END DO
!         CALL sol1d(j,k,b01d,s1dl)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         l=lout
         DO i=Ix_End_y(j),Ix_Start_y(j),-1
            x(1,l)=x(1,l)-s1dl(1,i)
            x(2,l)=x(2,l)-s1dl(2,i)
            l=l-1
         END DO
         lout=lout-nrnx(j)
         jp1=j
      END DO

      RETURN
      END SUBROUTINE sol2d
#endif 

#ifdef siarhei_delete 
      SUBROUTINE sol1d(irow,k,b,x)

      USE Inc_Lscoef
      USE Inc_BiCG

      IMPLICIT NONE
      !  solve 1D problem using predetermined LU factors
      INTEGER :: i, k, l, im1, ip1, irow, ibeg, iend
      REAL(8) :: b1i, b2i
      REAL(8) :: x(N_Group,Nx),b(N_Group,Nx) !dmm

      !  forward substitution
      ibeg=Ix_Start_y(irow)
      iend=Ix_End_y(irow)
      l=nodel(ibeg,irow)
      i=ibeg
      y(1,i)=delinv(1,l,k)*b(1,i)+delinv(2,l,k)*b(2,i)
      y(2,i)=delinv(3,l,k)*b(1,i)+delinv(4,l,k)*b(2,i)
      im1=i
      DO i=ibeg+1,iend
         l=l+1
         b1i=b(1,i)-al(1,l,k)*y(1,im1)-al(2,l,k)*y(2,im1)
         b2i=b(2,i)-al(3,l,k)*y(1,im1)-al(4,l,k)*y(2,im1)
         y(1,i)=delinv(1,l,k)*b1i+delinv(2,l,k)*b2i
         y(2,i)=delinv(3,l,k)*b1i+delinv(4,l,k)*b2i
         im1=i
      END DO
      !  backward substitution
      x(1,iend)=y(1,iend)
      x(2,iend)=y(2,iend)
      ip1=iend
      DO i=iend-1,ibeg,-1
         l=l-1
         x(1,i)=y(1,i)-deliau(1,l,k)*x(1,ip1)-deliau(2,l,k)*x(2,ip1)
         x(2,i)=y(2,i)-deliau(3,l,k)*x(1,ip1)-deliau(4,l,k)*x(2,ip1)
         ip1=i
      END DO

      RETURN
      END SUBROUTINE sol1d
#endif 

      FUNCTION scalp(x, y)

      USE Inc_Option, ONLY: N_Group

      IMPLICIT NONE
      ! scalar product of two vectors
      INTEGER :: k, l
      REAL(8) :: scalp, x(Nxy, Nz, N_Group), y(Nxy, Nz, N_Group), psum


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [scalp] in Mod_SolLS'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      psum = 0d0
      DO k = 1, Nz
         DO l = 1, Nxy
            psum = psum + x(l, k, 1)*y(l, k, 1) + x(l, k, 2)*y(l, k, 2)
         END DO
      END DO
      scalp = psum

      RETURN
      END FUNCTION scalp

#ifdef siarhei_delete 
      SUBROUTINE SolveLinearSystem_mg(ig,iin)

      IMPLICIT NONE
      INTEGER :: iin
      INTEGER :: ig

!      CALL BicgStab_mg(ig,iin)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 

      END SUBROUTINE SolveLinearSystem_mg
#endif 

#ifdef siarhei_delete 
      SUBROUTINE BicgStab_mg(ig,iin)

      USE Inc_Lscoef
      USE Inc_FluxVar
      USE Inc_BiCG

      IMPLICIT NONE
      ! solves a system of linear equations by preconditioned conjugate gradient method
      INTEGER, INTENT(out) :: iin
      INTEGER :: k, l, Ixy, Iz
      INTEGER :: ig
      REAL(8) :: r2t
      REAL(8) :: r20epserf,crhod
      REAL(8) :: r0v, pts, ptt, vol

!      CALL AxB_mg(ig, Flux, flux_add )
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 

      r20 = 0
      b2  = 0
      DO Iz = 1, Nz
         DO Ixy = 1, Nxy
            vr0(Ixy, Iz, ig) = src(ig, Ixy, Iz) - flux_add(Ixy, Iz, ig)
            if(abs(vr0(Ixy,Iz,ig))<1e-30) vr0(Ixy,Iz,ig)=1e-30
            vr(Ixy, Iz, ig) = vr0(Ixy, Iz, ig)
            r20 = r20 + vr0(Ixy, Iz, ig)*vr0(Ixy, Iz, ig)
            !b2  = b2  + src(ig, Ixy, Iz)*src(ig, Ixy, Iz)
            vp(Ixy, Iz, ig) = 0D0
            vv(Ixy, Iz, ig) = 0D0
         END DO
      END DO
      r20epserf = r20 * EPS_Residual * EPS_Residual
      r20 = SQRT(r20)
      b2  = SQRT(b2)
      calpha = 1
      crho   = 1
      comega = 1
      DO iin = 1, Iin_Max
         crhod = crho
!         crho  = scalp_mg(ig,vr0, vr)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         cbeta = crho*calpha / (crhod*comega)
         DO k = 1, Nz
            DO l = 1, Nxy
               vp(l, k, ig) = vr(l, k, ig) + cbeta *( vp(l, k, ig) - comega*vv(l, k, ig) )
            END DO
         END DO
!         CALL MInv_mg(ig,vp, vy)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
!         CALL AxB_mg(ig,vy, vv)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
!         r0v    = scalp_mg(ig,vr0, vv)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         calpha = crho / r0v
         DO k=1,Nz
            DO l=1,Nxy
               vs(l, k, ig) = vr(l, k, ig) - calpha*vv(l, k, ig)
               if(abs(vs(l,k,ig)).LE.1e-30) vs(l,k,ig)=1e-30       ! for vv=0 case
            ENDDO
         ENDDO
!         CALL MInv_mg(ig, vs, vz)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
!         CALL AxB_mg(ig, vz, vt)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
!         pts = scalp_mg(ig,vs,vt)
            dummy_filler = 1 ! @$^ siarhei_plot 
!         ptt = scalp_mg(ig,vt,vt)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         comega = pts / ptt
         r2t = 0
         DO k = 1, Nz
            DO l = 1, Nxy
               vol = NodeVolume(l, k)
               Flux(l,k,ig) = Flux(l,k,ig) + calpha*vy(l,k,ig) + comega*vz(l,k,ig)
               vr(l,k,ig) = vs(l,k,ig) - comega*vt(l,k,ig)
               r2t = r2t + vr(l,k,ig)*vr(l,k,ig)
            ENDDO
         ENDDO
         !IF(r2t.LT.r20epserf .and. iin>=3) then
         IF(r2t.LT.r20epserf)THEN
            RETURN
         ENDIF
      ENDDO
      !Iin = Iin_Max

      RETURN
      END SUBROUTINE BicgStab_mg
#endif 

#ifdef siarhei_delete 
      SUBROUTINE AxB_mg(ig,b,ab)

      USE Inc_Lscoef

      IMPLICIT NONE
      ! Matrix-vector product
      INTEGER :: k, kp1, km1, ln, ls, Ixy
      INTEGER :: ig
      REAL(8), INTENT(in) :: b(0:Nxy+1, Nz, N_Group)
      REAL(8), INTENT(out) :: ab(Nxy, Nz, N_Group)    !dmm

      DO k=1,Nz
         kp1=k+1
         km1=k-1
         IF(kp1.GT.Nz) kp1=Nz
         IF(km1.LT.1) km1=1
         DO Ixy = 1, Nxy
            ln = I_Ly_P2R(Ixy)
            ls = I_Ry_P2R(Ixy)
            ab(Ixy,k,ig) = am (ig,Ixy,k)*b(Ixy,  k,  ig)+ccw(ig,Ixy,k)*b(Ixy-1,k,  ig) &
                         + cce(ig,Ixy,k)*b(Ixy+1,k,  ig)+ccn(ig,Ixy,k)*b(ln   ,k,  ig) &
                         + ccs(ig,Ixy,k)*b(ls   ,k,  ig)+ccb(ig,Ixy,k)*b(Ixy  ,km1,ig) &
                         + cct(ig,Ixy,k)*b(Ixy  ,kp1,ig)
         END DO
      ENDDO

      END SUBROUTINE AxB_mg
#endif 

#ifdef siarhei_delete 
      SUBROUTINE Minv_mg(ig,b,x)

      USE Inc_Lscoef
      USE Inc_BiCG

      IMPLICIT NONE
      ! solve Mx=b for x given b
      REAL(8) :: x(0:Nxy+1, Nz, N_Group), b(Nxy, Nz, N_Group)
      INTEGER :: l, k, kp1
      INTEGER :: ig

      ! forward solve
      DO l=-Nx+1,Nxy+Nx
         s(ig,l)=0d0
      ENDDO
      DO k=1,Nz-1
         DO l=1,Nxy
            b0(ig,l)=b(l,k,ig)-ccb(ig,l,k)*s(ig,l)
         ENDDO
!         CALL sol2d_mg(ig,k,b0,s)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         DO l=1,Nxy
            x(l,k,ig) = s(ig,l)
         ENDDO
      ENDDO
      ! on the top plane
      k=Nz
      DO l=1,Nxy
         b0(ig,l) = b(l,k,ig) - ccb(ig,l,k)*s(ig,l)
      ENDDO
!      CALL sol2d_mg(ig,k,b0,s)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
      DO l=1,Nxy
         x(l,k,ig) = s(ig,l)
      ENDDO
      ! backward
      DO k = Nz-1,1,-1
         kp1=k+1
         DO l=1,Nxy
            b0(ig,l) = x(l,kp1,ig)*cct(ig,l,k)
         ENDDO
!         CALL sol2d_mg(ig,k,b0,s)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         DO l = 1, Nxy
            x(l,k,ig) = x(l,k,ig) - s(ig,l)
         ENDDO
      ENDDO

      RETURN

      END SUBROUTINE Minv_mg
#endif 

#ifdef siarhei_delete 
      SUBROUTINE sol2d_mg(ig,k,b,x)

      USE Inc_BiCG
      USE Inc_Lscoef

      IMPLICIT NONE
      ! solve a 2d problem using precalculated LU factors
      INTEGER :: i,j,k,l,jp1,ls,lout
      INTEGER :: ig
      REAL(8) :: b(N_Group,Nxy),x(N_Group,-Nx+1:Nxy+Nx) !dmm

      !  forward solve
      j=1
      DO i=1,Nx !dmm
         s1dl(ig,i)=0
      ENDDO
      lout=0
      DO j=1,Ny
         l=lout
         DO i=Ix_Start_y(j),Ix_End_y(j)
            l=l+1
            b01d(ig,i)=b(ig,l)-ccn(ig,l,k)*s1dl(ig,i)
         ENDDO
         ! solve 1d problem
!         CALL sol1d_mg(ig,j,k,b01d,s1dl)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         l=lout
         DO i=Ix_Start_y(j),Ix_End_y(j)
            l=l+1
            x(ig,l)=s1dl(ig,i)
         ENDDO
         lout=lout+nrnx(j)
      ENDDO
      !  backward solve
      lout=Nxy-nrnx(Ny)
      jp1=Ny
      DO j=Ny-1,1,-1
         l=lout
         DO i=Ix_End_y(j),Ix_Start_y(j),-1
            ls=nodel(i,jp1)
            b01d(ig,i)=x(ig,ls)*ccs(ig,l,k)
            l=l-1
         ENDDO
!         CALL sol1d_mg(ig,j,k,b01d,s1dl)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         l=lout
         DO i=Ix_End_y(j),Ix_Start_y(j),-1
            x(ig,l)=x(ig,l)-s1dl(ig,i)
            l=l-1
         ENDDO
         lout=lout-nrnx(j)
         jp1=j
      ENDDO

      RETURN
      END SUBROUTINE sol2d_mg
#endif 

#ifdef siarhei_delete 
      SUBROUTINE sol1d_mg(ig,irow,k,b,x)

      USE Inc_Lscoef
      USE Inc_BiCG

      IMPLICIT NONE
      !  solve 1D problem using predetermined LU factors
      INTEGER :: i, k, l, im1, ip1, irow, ibeg, iend
      INTEGER :: ig
      REAL(8) :: b1i
      REAL(8) :: x(N_Group,Nx),b(N_Group,Nx) !dmm

      !  forward substitution
      ibeg=Ix_Start_y(irow)
      iend=Ix_End_y(irow)
      l=nodel(ibeg,irow)
      i=ibeg
      y(ig,i)=delinv(ig,l,k)*b(ig,i)
      im1=i
      DO i=ibeg+1,iend
         l=l+1
         b1i=b(ig,i)-al(ig,l,k)*y(ig,im1)
         y(ig,i)=delinv(ig,l,k)*b1i
         im1=i
      ENDDO
      !  backward substitution
      x(ig,iend)=y(ig,iend)
      ip1=iend
      DO i=iend-1,ibeg,-1
         l=l-1
         x(ig,i)=y(ig,i)-deliau(ig,l,k)*x(ig,ip1)
         ip1=i
      ENDDO

      RETURN
      END SUBROUTINE sol1d_mg
#endif 

#ifdef siarhei_delete 
      FUNCTION scalp_mg(ig,x,y)

      USE Inc_Option , ONLY: N_Group

      IMPLICIT NONE
      ! scalar product of two vectors
      INTEGER :: k, l
      INTEGER :: ig
      REAL(8) :: scalp_mg, x(Nxy, Nz, N_Group), y(Nxy, Nz, N_Group), psum

      psum = 0d0
      DO k = 1, Nz
         DO l = 1, Nxy
            psum = psum + x(l,k,ig)*y(l,k,ig)
         ENDDO
      ENDDO
!      scalp_mg = psum
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 

      RETURN
      END FUNCTION scalp_mg
#endif 


      SUBROUTINE setls

      USE Mod_Alloc
      USE Inc_Lscoef
      USE Inc_FluxVar
      USE Inc_RP
      USE Inc_maXS, ONLY: maXS_r_3D, nu_maXS_f_3D, maXS_s_3D
      USE Inc_Transient
      USE Inc_Solver

      IMPLICIT NONE
      ! linear system setup for cartesian geometry
      INTEGER :: i, j, k, l, m, kt, kb, le, ln, ls, js, lw
      REAL(8) :: hxi, hxy, hyz, hxz, hyj, hzk, vol

      ! begin the generation of the coefficient matrix

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [setls] in Mod_SolLS'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      DO k=1,Nz
         kt=k+1
         hzk = MeshSize_z(k)
         DO j=1,Ny
            hyj = MeshSize_y(j)
            hyz = hzk*hyj
            js=j+1
            DO i=Ix_Start_y(j),Ix_End_y(j)
               hxi = MeshSize_x(i)
               hxy = hxi*hyj
               hxz = hxi*hzk
               l   = nodel(i,j)
               le  = nodel(i+1,j)
               ls  = nodel(i,js)
               vol = NodeVolume(l,k)
               DO m=1,N_Group
                  ccw(m,l,k)=(-dfw(m,l,k)+dnw(m,l,k))*hyz
                  ccn(m,l,k)=(-dfn(m,l,k)+dnn(m,l,k))*hxz
                  ccb(m,l,k)=(-dfb(m,l,k)+dnb(m,l,k))*hxy
                  cce(m,l,k)=(-dfw(m,le,k)-dnw(m,le,k))*hyz
                  ccs(m,l,k)=(-dfn(m,ls,k)-dnn(m,ls,k))*hxz
                  cct(m,l,k)=(-dfb(m,l,kt)-dnb(m,l,kt))*hxy
                  am(m,l,k)=(maXs_r_3D(l, k, m)+rvdelt(m,l,k))*vol  &
                       +(dfw(m,l,k)+dnw(m,l,k))*hyz                 &
                       +(dfn(m,l,k)+dnn(m,l,k))*hxz                 &
                       +(dfb(m,l,k)+dnb(m,l,k))*hxy                 &
                       +(dfw(m,le,k)-dnw(m,le,k))*hyz               &
                       +(dfn(m,ls,k)-dnn(m,ls,k))*hxz               &
                       +(dfb(m,l,kt)-dnb(m,l,kt))*hxy
                  af(m,l,k)=betap(l,k)*nu_maXs_f_3D(l, k, m)*vol
               END DO
               scat(l, k) = maXS_s_3D(l, k, 1)*vol
            END DO
         END DO
      END DO
      ! zero out unnecessary boundary coupling coefficients
      kb=1
      kt=Nz
      DO l=1,Nxy
         DO m=1,N_Group
             ccb(m,l,kb)=zero
             cct(m,l,kt)=zero
         END DO
      END DO
      DO k=1,Nz
         DO i=1,Nx
            ln=nodel(i,Iy_Start_x(i))
            ls=nodel(i,Iy_End_x(i))
            DO m=1,N_Group
                ccn(m,ln,k)=zero
                ccs(m,ls,k)=zero
            END DO
         END DO
         DO j=1,Ny
            lw=nodel(Ix_Start_y(j),j)
            le=nodel(Ix_End_y(j),j)
            DO m=1,N_Group
                ccw(m,lw,k)=zero
                cce(m,le,k)=zero
            END DO
         END DO
      END DO
      ! linear system for the fsp-like problem: note regivs=0 for chebyshev
      DO k=1,Nz
         DO l=1,Nxy
            am(1,l,k)=am(1,l,k)-af(1,l,k)*reigvs
            af2(l,k)=af(2,l,k)*reigvs
         END DO
      END DO

      END SUBROUTINE setls

#ifdef siarhei_delete 
      SUBROUTINE SolveLinearSystem_det(iin,checkfiss)
      IMPLICIT NONE
      INTEGER:: iin
      LOGICAL :: checkfiss


!      CALL BicgStab_det(iin,checkfiss)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
      END SUBROUTINE SolveLinearSystem_det
#endif 


#ifdef siarhei_delete 
      SUBROUTINE BicgStab_det(iin,checkfiss)
      USE Inc_Lscoef
      USE Inc_FluxVar
      USE Inc_BiCG
      USE Inc_3D
      USE Inc_Detector, only: DET_C, flag_det_mk_matrix, DET_R, flag_det_mk_ainv
      IMPLICIT NONE
      ! solves a system of linear equations by preconditioned conjugate gradient method
      INTEGER, INTENT(out) :: iin
      INTEGER :: k, l, Ixy, Iz
      REAL(8) :: errl2t,r2t,err,psipsidt,errlinft
      REAL(8) :: r20epserf, crhod
      REAL(8) :: r0v, pts, ptt, vol
      LOGICAL, INTENT(in) :: checkfiss

      flag_det_mk_ainv = .true.
      flag_det_mk_matrix = .true.
!      CALL AxB_det( Flux, flux_add )
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
      r20 = 0
      b2  = 0
      DO Iz = 1, Nz
         DO Ixy = 1, Nxy
            vr0(Ixy, Iz, 1) = (DET_C(Ixy, Iz, 1) - flux_add(Ixy, Iz, 1))
            vr0(Ixy, Iz, 2) = (DET_C(Ixy, Iz, 2) - flux_add(Ixy, Iz, 2))
            vr (Ixy, Iz, 1) = vr0(Ixy, Iz, 1)
            vr (Ixy, Iz, 2) = vr0(Ixy, Iz, 2)
            r20 = r20 + vr0(Ixy, Iz, 1)**2 + vr0(Ixy, Iz, 2)**2
            b2  = b2  + DET_C(Ixy, Iz, 1)**2 + DET_C(Ixy, Iz, 2)**2
            vp(Ixy, Iz, 1) = D0
            vp(Ixy, Iz, 2) = D0
            vv(Ixy, Iz, 1) = D0
            vv(Ixy, Iz, 2) = D0
         END DO
      END DO
      r20epserf = r20 * EPS_Residual * EPS_Residual
      r20 = SQRT(r20)
      b2  = SQRT(b2)
      calpha = 1
      crho   = 1
      comega = 1
      DO Iin = 1, Iin_Max
         crhod = crho
         crho  = scalp(vr0, vr)
         cbeta = crho*calpha / (crhod*comega)
         DO k = 1, Nz
            DO l = 1, Nxy
               vp(l, k, 1) = vr(l, k, 1) + cbeta *( vp(l, k, 1) - comega*vv(l, k, 1) )
               vp(l, k, 2) = vr(l, k, 2) + cbeta *( vp(l, k, 2) - comega*vv(l, k, 2) )
            END DO
         END DO
!         CALL RINV_CHOLESKY(DET_R,Nxy*Nz*N_group,vp, vy)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         flag_det_mk_matrix = .false.
!         CALL AxB_det(vy, vv)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         r0v    = scalp(vr0, vv)
         calpha = crho / r0v
         DO k=1,Nz
            DO l=1,Nxy
               vs(l, k, 1) = vr(l, k, 1) - calpha*vv(l, k, 1)
               vs(l, k, 2) = vr(l, k, 2) - calpha*vv(l, k, 2)
            END DO
         END DO
!         CALL RINV_CHOLESKY(DET_R,Nxy*Nz*N_group,vs,vz)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
!         CALL AxB_det(vz, vt)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         pts = scalp(vs, vt)
         ptt = scalp(vt, vt)
         comega = pts / ptt
         IF ( checkfiss ) THEN
            flagr2= .FALSE.
            flagl2= .FALSE.
            flaglinf= .FALSE.
            flagerf= .FALSE.
            errl2t=0
            r2t=0
            psipsidt=0
            errlinft=0
            DO k = 1, Nz
               DO l = 1,Nxy
                  vol = NodeVolume(l,k)
                  Flux(l, k, 1) = Flux(l, k, 1) + calpha*vy(l, k, 1) + comega*vz(l, k, 1)
                  Flux(l, k, 2) = Flux(l, k, 2) + calpha*vy(l, k, 2) + comega*vz(l, k, 2)
                  vr(l, k, 1) = vs(l, k, 1) - comega*vt(l, k, 1)
                  vr(l, k, 2) = vs(l, k, 2) - comega*vt(l, k, 2)
                  r2t = r2t + vr(l, k, 1)**2 + vr(l, k, 2)**2
                  IF ( nu_maXs_f_3D(l, k, 2) /= D0 ) THEN
                     FisSrc_Iout(l, k) = FisSrc(l,k)
                     FisSrc(l, k) = ( nu_maXs_f_3D(l, k, 1)*Flux(l, k, 1) + nu_maXs_f_3D(l, k, 2)*Flux(l, k, 2) ) * vol
                     err = FisSrc(l, k) - FisSrc_Iout(l, k)
                     errlinft =MAX( errlinft, ABS( err/FisSrc(l,k) ) )
                     errl2t   = errl2t + err**2
                     psipsidt = psipsidt + FisSrc(l, k)*FisSrc_Iout(l, k)
                  END IF
               END DO
            END DO
            errl2=SQRT(errl2t/ABS(psipsidt))
            errlinf=errlinft
            r2=SQRT(r2t)
            r2ob2=r2/b2
            r2=r2/r20
            IF( r2ob2   <= EPS_Residual ) flagr2   = .TRUE.
            IF( errl2   <= EPS_Global   ) flagl2   = .TRUE.
            IF( errlinf <= EPS_Local    ) flaglinf = .TRUE.
            IF( r2 <= EPS_ERF ) THEN
               flagerf = .TRUE.
               RETURN
            END IF
         ELSE
            r2t = 0
            DO k = 1, Nz
               DO l = 1, Nxy
                  vol = NodeVolume(l, k)
                  Flux(l, k, 1) = Flux(l, k, 1) + calpha*vy(l, k, 1) + comega*vz(l, k, 1)
                  Flux(l, k, 2) = Flux(l, k, 2) + calpha*vy(l, k, 2) + comega*vz(l, k, 2)
                  vr(l, k, 1) = vs(l, k, 1) - comega*vt(l, k, 1)
                  vr(l, k, 2) = vs(l, k, 2) - comega*vt(l, k, 2)
                  r2t = r2t + vr(l, k, 1)**2 + vr(l, k, 2)**2
               END DO
            END DO
            IF(abs(r2t).LT.abs(r20epserf))RETURN
         END IF
      END DO
      Iin = Iin_Max
      do k=1,Nz
         do l=1,Nxy
            if (Flux(l,k,1)<1d-30) Flux(l,k,1)=1d-30
            if (Flux(l,k,2)<1d-30) Flux(l,k,2)=1d-30
         enddo
      enddo
      RETURN
      END SUBROUTINE BicgStab_det
#endif 

#ifdef siarhei_delete 
      SUBROUTINE AxB_det(b,ab)

      USE Inc_Lscoef
      USE Inc_Detector,   only: DET_A, DET_F, DET_R, DET_D, DET_C,   &
                               DET_NPOW, DET_NPOW_Ixy, DET_NPOW_Iz, DET_NPOW_POW,  &
                               flag_det_mk_matrix
      USE Inc_Fluxvar,    only: src
      USE Inc_MaXS,       only: kap_maXS_f_3D
      IMPLICIT NONE
      ! Matrix-vector product
      INTEGER :: k, kp1, km1, ln, ls, Ixy
      REAL(8), INTENT(in) :: b(0:Nxy+1, Nz, N_Group)
      REAL(8), INTENT(out) :: ab(Nxy, Nz, N_Group)   !dmm for 3D VER
      INTEGER(4) :: row_p, col_p, tmp_ixy, tmp_iz, tmp_ig
      real(8)    :: mid
      integer(4) :: row_p_1, ixy_b, iz_b, ig_b, i, j, Iz, Ig
      integer(4) :: r, rr
      integer(4), allocatable, dimension(:,:) :: map_A
      integer(4), allocatable, dimension(:,:) :: map_AT
      integer(4), allocatable, dimension(:,:) :: n_nonzero_A
      integer(4), allocatable, dimension(:,:) :: n_nonzero_AT

      !! INFOMARION ---------------------------------------- !!
      !     MATRIX DET_A, DET_D, DET_R                        !
      !     [DET_R]= ([DET_A][DET_D])*([DET_A])               !
      !                               ([DET_D])               !
      !     Ig > Ixy > Iz  row1. (1,1,1)   row5. (1,1,2)      !
      !      (Ixy,Iz,Ig)   row2. (1,2,1)   row6. (1,2,2)      !
      !                    row3. (2,1,1)   row7. (2,1,2)      !
      !                    row4. (2,2,1)   row8. (2,2,2)      !
      !! --------------------------------------------------- !!
      IF (flag_det_mk_matrix) THEN
       DET_A = 0d0
       DET_F = 0d0
       DO k=1,Nz
          kp1=k+1
          km1=k-1
          IF(kp1.GT.Nz) kp1=Nz
          IF(km1.LT.1) km1=1
          DO Ixy = 1, Nxy
             ln = I_Ly_P2R(Ixy)
             ls = I_Ry_P2R(Ixy)
             ab(Ixy, k, 1) = am ( 1, Ixy, k ) * b( Ixy    , k  , 1 )  &
                           + ccw( 1, Ixy, k ) * b( Ixy - 1, k  , 1 )  &
                           + cce( 1, Ixy, k ) * b( Ixy + 1, k  , 1 )  &
                           + ccn( 1, Ixy, k ) * b( ln     , k  , 1 )  &
                           + ccs( 1, Ixy, k ) * b( ls     , k  , 1 )  &
                           + ccb( 1, Ixy, k ) * b( Ixy    , km1, 1 )  &
                           + cct( 1, Ixy, k ) * b( Ixy    , kp1, 1 )  &
                           - af2   ( Ixy, k ) * b( Ixy    , k  , 2 )
             ab(Ixy, k, 2) = am ( 2, Ixy, k ) * b( Ixy    , k  , 2 )  &
                           + ccw( 2, Ixy, k ) * b( Ixy - 1, k  , 2 )  &
                           + cce( 2, Ixy, k ) * b( Ixy + 1, k  , 2 )  &
                           + ccn( 2, Ixy, k ) * b( ln     , k  , 2 )  &
                           + ccs( 2, Ixy, k ) * b( ls     , k  , 2 )  &
                           + ccb( 2, Ixy, k ) * b( Ixy    , km1, 2 )  &
                           + cct( 2, Ixy, k ) * b( Ixy    , kp1, 2 )  &
                           - scat  ( Ixy, k ) * b( Ixy    , k  , 1 )
             !! Make [M-reigev*F] MATRIX
             ! ===== group 1
             row_p = Nz*(Ixy-1)+Nxy*Nz*(1-1)+k
             ! am
             tmp_ixy  = Ixy
             tmp_iz   = k
             tmp_ig   = 1
             col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
             DET_A (row_p, col_p) = DET_A (row_p, col_p) +  am(1,Ixy,k) !+af(1,Ixy,k)*reigvs
             ! ccw
             tmp_ixy  = Ixy-1
             if (tmp_ixy > 0) then
                tmp_iz   = k
                tmp_ig   = 1
                col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
                DET_A (row_p, col_p) = DET_A (row_p, col_p) +  ccw(1,Ixy,k)
             endif
             ! cce
             tmp_ixy  = Ixy+1
             if (tmp_ixy < Nxy+1) then
                tmp_iz   = k
                tmp_ig   = 1
                col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
                DET_A (row_p, col_p) = DET_A (row_p, col_p) +  cce(1,Ixy,k)
             endif
             ! ccn
             tmp_ixy  = ln
             if (tmp_ixy > 0) then
                tmp_iz   = k
                tmp_ig   = 1
                col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
                DET_A (row_p, col_p) = DET_A (row_p, col_p) +  ccn(1,Ixy,k)
             endif
             ! ccs
             tmp_ixy  = ls
             if (tmp_ixy < Nxy+1) then
                tmp_iz   = k
                tmp_ig   = 1
                col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
                DET_A (row_p, col_p) = DET_A (row_p, col_p) +  ccs(1,Ixy,k)
             endif
             ! ccb
             tmp_ixy  = Ixy
             tmp_iz   = km1
             tmp_ig   = 1
             col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
             DET_A (row_p, col_p) = DET_A (row_p, col_p) +  ccb(1,Ixy,k)
             ! cct
             tmp_ixy  = Ixy
             tmp_iz   = kp1
             tmp_ig   = 1
             col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
             DET_A (row_p, col_p) = DET_A (row_p, col_p) +  cct(1,Ixy,k)
             ! -af2
             tmp_ixy  = Ixy
             tmp_iz   = k
             tmp_ig   = 2
             col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
             DET_A (row_p, col_p) = DET_A (row_p, col_p) -  af2(Ixy,k)

             ! ===== group 2
             row_p = Nz*(Ixy-1)+Nxy*Nz*(2-1)+k
             ! am
             tmp_ixy  = Ixy
             tmp_iz   = k
             tmp_ig   = 2
             col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
             DET_A (row_p, col_p) = DET_A (row_p, col_p) +  am(2,Ixy,k)
             ! ccw
             tmp_ixy  = Ixy-1
             if (tmp_ixy > 0) then
                tmp_iz   = k
                tmp_ig   = 2
                col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
                DET_A (row_p, col_p) = DET_A (row_p, col_p) +  ccw(2,Ixy,k)
             endif
             ! cce
             tmp_ixy  = Ixy+1
             if (tmp_ixy < Nxy+1) then
                tmp_iz   = k
                tmp_ig   = 2
                col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
                DET_A (row_p, col_p) = DET_A (row_p, col_p) +  cce(2,Ixy,k)
             endif
             ! ccn
             tmp_ixy  = ln
             if (tmp_ixy > 0) then
                tmp_iz   = k
                tmp_ig   = 2
                col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
                DET_A (row_p, col_p) = DET_A (row_p, col_p) +  ccn(2,Ixy,k)
             endif
             ! ccs
             tmp_ixy  = ls
             if (tmp_ixy < Nxy+1) then
                tmp_iz   = k
                tmp_ig   = 2
                col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
                DET_A (row_p, col_p) = DET_A (row_p, col_p) +  ccs(2,Ixy,k)
             endif
             ! ccb
             tmp_ixy  = Ixy
             tmp_iz   = km1
             tmp_ig   = 2
             col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
             DET_A (row_p, col_p) = DET_A (row_p, col_p) +  ccb(2,Ixy,k)
             ! cct
             tmp_ixy  = Ixy
             tmp_iz   = kp1
             tmp_ig   = 2
             col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
             DET_A (row_p, col_p) = DET_A (row_p, col_p) +  cct(2,Ixy,k)
             ! -scat
             tmp_ixy  = Ixy
             tmp_iz   = k
             tmp_ig   = 1
             col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
             DET_A (row_p, col_p) = DET_A (row_p, col_p) -  scat(Ixy,k)
          END DO
       END DO
       !! MAKE DETECTOR MATRIX
       IF (allocated(DET_D)) deallocate(DET_D)
       allocate(DET_D(DET_NPOW,Nxy*Nz*N_Group))
       DET_D = 0d0
       i = 1
       Do i = 1, DET_NPOW
          tmp_ixy  = DET_NPOW_Ixy(i)
          tmp_iz   = DET_NPOW_Iz(i)
          tmp_ig   = 1
          col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
          DET_D(i, col_p) = kap_maXS_f_3D(tmp_ixy, tmp_Iz, tmp_ig)
          tmp_ig   = 2
          col_p = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
          DET_D(i, col_p) = kap_maXS_f_3D(tmp_Ixy, tmp_Iz, tmp_ig)
       ENDDO
       !! MAKE R MATRIX
       DET_R = 0d0
!!!! ================ DET_R cal with index ======================================= !!!
       allocate(map_A(Nxy*Nz*N_group,Nxy*Nz*N_group))
       allocate(n_nonzero_A(Nxy*Nz*N_group,1))
       allocate(map_AT(Nxy*Nz*N_group,Nxy*Nz*N_group))
       allocate(n_nonzero_AT(Nxy*Nz*N_group,1))
       map_A = 0
       n_nonzero_A = 0
       map_AT = 0
       n_nonzero_AT = 0
       do j = 1,Nxy*Nz*N_group
          do i = 1,Nxy*Nz*N_group
             if (DET_A(i,j) .EQ. 0) cycle
             n_nonzero_A(j,1) = n_nonzero_A(j,1)+1
             map_A(n_nonzero_A(j,1),j)=i
             n_nonzero_AT(i,1) = n_nonzero_AT(i,1)+1
             map_AT(i,n_nonzero_AT(i,1))=j
          enddo
       enddo
       ! diagonal term
       do j = 1,Nxy*Nz*N_group
          do i = 1,(n_nonzero_A(j,1))
             k = map_A(i,j)
             do r = 1,(n_nonzero_AT(k,1))
                rr = map_AT(k,r)
                DET_R(rr,j)=DET_R(rr,j)+DET_A(k,rr)*DET_A(k,j)
             enddo
          enddo
       enddo
       deallocate(map_A)
       deallocate(map_AT)
       deallocate(n_nonzero_A)
       deallocate(n_nonzero_AT)
      ENDIF
      ab = 0d0
      do tmp_ig = 1,2
         Do tmp_ixy = 1, Nxy
            DO tmp_iz = 1, Nz
               col_p   = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
               do ig = 1, 2
                  Do Ixy = 1, Nxy
                     Do Iz = 1, Nz
                        row_p_1   = Nz*(ixy-1)+Nxy*Nz*(ig-1)+iz  ! group 1
                        if (DET_R(row_p_1,col_P) .EQ. 0) cycle
                        ab(Ixy,Iz,ig) = ab(Ixy,Iz,ig) + DET_R(row_p_1,col_P)*b(tmp_Ixy,tmp_Iz,tmp_ig)
                     ENDDO
                  ENDDO
               enddo
            ENDDO
         ENDDO
      enddo
!! === ab calculate end ===================================================== !!
      !! SRC UPdate
      DET_C = 0d0
      Do Ig = 1, N_group
         Do Ixy = 1, Nxy
            Do Iz = 1,Nz
               tmp_ixy = Ixy
               tmp_iz  = Iz
               tmp_ig  = ig
               row_p   = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz  ! group 1
               mid     = 0d0
               Do Ig_b  = 1, N_group
                  DO Ixy_b = 1, Nxy
                     DO Iz_b = 1, Nz
                        tmp_ixy = Ixy_b
                        tmp_iz  = Iz_b
                        tmp_ig  = Ig_b
                        col_p   = Nz*(tmp_ixy-1)+Nxy*Nz*(tmp_ig-1)+tmp_iz
                        mid     = mid + DET_A(col_p,row_p)*src(Ig_b,Ixy_b,Iz_b)
                     ENDDO
                  ENDDO
               ENDDO
               DO i = 1, DET_NPOW
                  mid = mid + DET_D(i,row_p)*DET_NPOW_POW(i) ! group 1
               ENDDO
               DET_C(Ixy,Iz,ig) = mid
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE AxB_det
#endif 

#ifdef siarhei_delete 
      SUBROUTINE RINV_CHOLESKY_R(a,n,ss,ff)
      !! Cholesky decomposition
      USE Inc_Detector, only: DET_RINV, &
                              flag_det_mk_ainv
      IMPlICIT NONE
      INTEGER(4) :: n, i, j, k
      INTEGER(4) :: ixy, iz, ig
      REAL(8) :: mid
      real(8) :: max_val
      double precision a(n,n) !, c(n,n)
      double precision L(n,n), U(n) !, b(n), d(n), x(n)
      double precision e(n,n), y(n), r(n,n), ae(n), z(n), v(n)
      double precision v1m(n),v2m(n)
      REAL(8) :: ff(0:Nxy+1, Nz, N_Group), ss(Nxy, Nz, N_Group)

      e=0d0
      u=0d0
      r=0d0
      ae=0d0
      y=0d0
      v1m=0d0
      v2m=0d0
      z=0d0
      v=0d0
      max_val=0d0
      l=0d0

      !! L  --------------------------
      if (flag_det_mk_ainv) then
         do k = 1,n
            do i = 1,n
               if (a(k,i) == 0d0) cycle
               e(k,i) = a(k,i)
            enddo
         enddo
         max_val = maxval(a)
         do k = 1,n
            if (e(k,k) < 0d0) cycle
            e(k,k) = e(k,k)**(0.5d0)
            L(k,k) = e(k,k)
            do j = (k+1),n
               if (e(k,j) == 0d0) cycle
               if (e(j,k) == 0d0) cycle
               e(k,j)=e(k,j)/(e(k,k)**0.5d0)
               if (a(k,j) < (0.02d0*max_val)) cycle
               L(k,j)=e(k,j)
            enddo
            do i = (k+1),n
               do j = (k+1),n
                  if (a(k,i) < (0.02d0*max_val)) cycle
                  if (a(k,j) < (0.02d0*max_val)) cycle
                  if (a(i,k) < (0.02d0*max_val)) cycle
                  if (a(j,k) < (0.02d0*max_val)) cycle
                  if (a(i,j) < (0.02d0*max_val)) cycle
                  if (a(j,i) < (0.02d0*max_val)) cycle
                  if (e(k,i) == 0d0) cycle
                  if (e(i,k) == 0d0) cycle
                  if (e(k,j) == 0d0) cycle
                  if (e(j,k) == 0d0) cycle
                  e(i,j)=e(i,j)-e(k,i)*e(k,j)
               enddo
            enddo
         enddo
         do i = 1,n
            r(i,:) = L(i,:)
            DET_RINV(i,:) = L(i,:)
         enddo
         flag_det_mk_ainv = .false.
      else
         do i = 1,n
            r(i,:) = DET_RINV(i,:)
         enddo
      endif
      !! Initialization
      do ixy=1,Nxy
         do iz=1,Nz
            do ig=1,N_group
               mid = Nz*(ixy-1)+Nxy*Nz*(ig-1)+iz
               y(int(mid))=ss(ixy,iz,ig)
            enddo
         enddo
      enddo
      !y = ss
      !z = ff
      !! y = M*z = (R'*R)*z
      !! v = R*z
      !! y = R'*v
      !! v = inv(R')*y

      !! solve v by forward substitution
      v(1)=y(1)/R(1,1)
      do i = 2,n
         mid = y(i)
         do j = 1, i-1
            if (R(j,i) == 0) cycle
            mid = mid - R(j,i)*v(j)
         enddo
         v(i) = mid/R(i,i)
      enddo

      !! solve z by backward substitution
      z(n) = v(n)/R(n,n)
      do i = n-1,1,-1
         !mid = 0d0
         mid = v(i)
         do j = i+1,n
            if (R(i,j) == 0) cycle
            mid = mid - R(i,j)*z(j)
         enddo
         z(i) = mid / R(i,i)
      enddo

      !! update flux
      do ixy = 1,Nxy
         do iz = 1,Nz
            do ig = 1,N_group
               mid = Nz*(ixy-1)+Nxy*Nz*(ig-1)+iz
               ff(ixy,iz,ig) = v(int(mid))
            enddo
         enddo
      enddo

      END SUBROUTINE RINV_CHOLESKY_R
#endif 

#ifdef siarhei_delete 
      SUBROUTINE RINV_CHOLESKY(a,n,ss,ff)
      !! Cholesky decomposition
      USE Inc_Detector, only: DET_RINV, &
                              flag_det_mk_ainv
      IMPlICIT NONE
      INTEGER(4) :: n, i, j, k
      INTEGER(4) :: ixy, iz, ig
      REAL(8) :: mid, mid_2
      real(8) :: max_val
      double precision a(n,n) !, c(n,n)
      double precision L(n,n), U(n) !, b(n), d(n), x(n)
      double precision e(n,n), y(n), r(n,n), ae(n), z(n), v(n)
      double precision v1m(n),v2m(n)
      REAL(8) :: ff(0:Nxy+1, Nz, N_Group), ss(Nxy, Nz, N_Group)

      e=0d0
      u=0d0
      r=0d0
      ae=0d0
      y=0d0
      v1m=0d0
      v2m=0d0
      z=0d0
      v=0d0
      max_val=0d0
      l=0d0

      !! L  --------------------------
      if (flag_det_mk_ainv) then
         max_val = maxval(a)
         L(1,1) = a(1,1)**0.5d0
         do k = 2,n
            if (a(k,1) > (0.02d0*max_val)) L(k,1) = a(k,1)/L(1,1)
            do i = 2,(k-1)
               mid=0d0
               if (a(k,i) < (0.02d0*max_val)) cycle
               if (a(i,k) < (0.02d0*max_val)) cycle
               do j = 1,(i-1)
                  if (L(i,j) == 0) cycle
                  if (L(k,j) == 0) cycle
                  mid=mid+L(i,j)*L(k,j)
               enddo
               L(k,i)=1d0/L(1,1)*(a(k,i)-mid)
            enddo
            mid_2=0d0
            do j = 1,(k-1)
               if (L(k,j) == 0) cycle
               mid_2=mid_2+L(k,j)*L(k,j)
            enddo
            if (a(k,k) < mid_2) cycle
            L(k,k)=(a(k,k)-mid_2)**0.5d0
         enddo
         do i = 1,n
              r(:,i) = L(i,:)
              DET_RINV(:,i) = L(i,:)
         enddo
         flag_det_mk_ainv = .false.
      else
         do i = 1,n
            r(:,i) = DET_RINV(:,i)
         enddo
      endif

      !! Initialization
      do ixy=1,Nxy
         do iz=1,Nz
            do ig=1,N_group
               mid = Nz*(ixy-1)+Nxy*Nz*(ig-1)+iz
               y(int(mid))=ss(ixy,iz,ig)
            enddo
         enddo
      enddo
      !y = ss
      !z = ff
      !! y = M*z = (R'*R)*z
      !! v = R*z
      !! y = R'*v
      !! v = inv(R')*y

      !! solve v by forward substitution
      v(1)=y(1)/R(1,1)
      do i = 2,n
         mid = y(i)
         do j = 1, i-1
            if (R(j,i) == 0) cycle
            mid = mid - R(j,i)*v(j)
         enddo
         v(i) = mid/R(i,i)
      enddo

      !! solve z by backward substitution
      z(n) = v(n)/R(n,n)
      do i = n-1,1,-1
         !mid = 0d0
         mid = v(i)
         do j = i+1,n
            if (R(i,j) == 0) cycle
            mid = mid - R(i,j)*z(j)
         enddo
         z(i) = mid / R(i,i)
      enddo

      !! update flux
      do ixy = 1,Nxy
         do iz = 1,Nz
            do ig = 1,N_group
               mid = Nz*(ixy-1)+Nxy*Nz*(ig-1)+iz
               ff(ixy,iz,ig) = v(int(mid))
            enddo
         enddo
      enddo

      END SUBROUTINE RINV_CHOLESKY
#endif 

#ifdef siarhei_delete 
      SUBROUTINE AINV_QR(a,n,ss,ff)
      USE Inc_Detector, only: DET_R, DET_RINV, &
                              flag_det_mk_ainv
      IMPlICIT NONE
      INTEGER(4) :: n, i, j, k
      INTEGER(4) :: ixy, iz, ig
      REAL(8) :: mid
      real(8) :: max_val
      real(8) :: norm
      double precision a(n,n) !, c(n,n)
      double precision U(n) !, b(n), d(n), x(n)
      double precision e(n,n), y(n), r(n,n), ae(n), q(n,n), z(n), v(n)
      double precision v1m(n),v2m(n)
      REAL(8) :: ff(0:Nxy+1, Nz, N_Group), ss(Nxy, Nz, N_Group)

      e=0d0
      u=0d0
      r=0d0
      ae=0d0
      y=0d0
      v1m=0d0
      v2m=0d0
      z=0d0
      v=0d0
      max_val=0d0

      !call cpu_time(timer(2))
      !! QR --------------------------
      max_val = maxval(DET_R)
      if (flag_det_mk_ainv) then
         do i = 1,n
            mid = 0d0
            do k = 1,n
               if (a(i,k) == 0) cycle
               !if (i /= k) cycle
               mid = mid + a(i,k)*a(i,k)
            enddo
            norm = mid**0.5d0
            !! row
            r(i,i) = norm
            q(i,:)   = a(i,:)/r(i,i)
            do j = i+1,n
               if (DET_R(i,j) == 0) cycle
               r(i,j)= sum(q(i,:)*a(:,j))
               if (r(i,j) < (0.02d0*max_val)) then
                  r(i,j) = 0d0
                  cycle
               endif
               do k = 1,n
                  a(k,j)=a(k,j)-q(i,k)*r(i,j)
               enddo
            enddo
         enddo
         do i = 1,n
            DET_RINV(:,i) = r(:,i)
         enddo
      else
         do i = 1,n
            r(:,i) = DET_RINV(:,i)
         enddo
      endif
      !! Initialization
      do ixy=1,Nxy
         do iz=1,Nz
            do ig=1,N_group
               mid = Nz*(ixy-1)+Nxy*Nz*(ig-1)+iz
               y(int(mid))=ss(ixy,iz,ig)
            enddo
         enddo
      enddo
      !y = ss
      !z = ff
      !! y = M*z = (R'*R)*z
      !! v = R*z
      !! y = R'*v
      !! v = inv(R')*y

      !! solve v by forward substitution
      v(1)=y(1)/R(1,1)
      do i = 2,n
         mid = y(i)
         do j = 1, i-1
            if (R(j,i) == 0) cycle
            mid = mid - R(j,i)*v(j)
         enddo
         v(i) = mid/R(i,i)
      enddo

      !! solve z by backward substitution
      z(n) = v(n)/R(n,n)
      do i = n-1,1,-1
         mid = v(i)
         do j = i+1,n
            if (R(i,j) == 0) cycle
            mid = mid - R(i,j)*z(j)
         enddo
         z(i) = mid / R(i,i)
      enddo

      !! update flux
      do ixy = 1,Nxy
         do iz = 1,Nz
            do ig = 1,N_group
               mid = Nz*(ixy-1)+Nxy*Nz*(ig-1)+iz
               ff(ixy,iz,ig) = z(int(mid))
            enddo
         enddo
      enddo

      END SUBROUTINE AINV_QR
#endif 

#ifdef siarhei_delete 
      SUBROUTINE AINV_SGS(n,ss,ff)
      USE Inc_Detector, only: DET_R
      IMPlICIT NONE
      INTEGER(4) :: n, i, j, k
      INTEGER(4) :: ixy, iz, ig
      REAL(8) :: mid, sum
      real(8) :: max_val
      double precision U(n) !L(n,n), U(n) !, b(n), d(n), x(n)
      double precision e(n,n), y(n), r(n,n), ae(n), z(n), v(n)
      double precision v1m(n),v2m(n)
      REAL(8) :: ff(0:Nxy+1, Nz, N_Group), ss(Nxy, Nz, N_Group)
      integer(4),dimension(:,:),allocatable :: n_nonzero_R
      real(8),dimension(:,:),allocatable    :: rcd_R
      integer(4),dimension(:,:),allocatable :: map_R

      allocate(n_nonzero_R(n,1))
      allocate(map_R(n,n))
      allocate(rcd_R(n,n))
      n_nonzero_R = 0 ; map_R = 0 ; rcd_R = 0d0

      e=0d0
      u=0d0
      r=0d0
      ae=0d0
      y=0d0
      v1m=0d0
      v2m=0d0
      z=0d0
      v=0d0
      max_val=0d0

      do ixy = 1,Nxy
         do iz = 1,Nz
            do ig = 1,N_group
               mid = Nz*(ixy-1)+Nxy*Nz*(ig-1)+iz
               z(int(mid))=ff(ixy,iz,ig)
               y(int(mid))=ss(ixy,iz,ig)
            enddo
         enddo
      enddo
      !! DET_R * z = y
      !! QR --------------------------
      do j = 1,n
         do i = 1,n
            if (DET_R(i,j) .EQ. 0) cycle
            n_nonzero_R(j,1)=n_nonzero_R(j,1)+1
            map_R(n_nonzero_R(j,1),j)=i
         enddo
      enddo
      do j = 1,n
         sum = 0d0
         do i = 1,n_nonzero_R(j,1)
            k = map_R(i,j)
            sum = sum + rcd_R(k,j)*z(k)
         enddo
         z(j) = (y(j)-sum)/R(j,j)
      enddo
      do j = 1,n
         do i = 1,n_nonzero_R(j,1)
            k = map_R(i,j)
            rcd_R(i,j) = DET_R(k,j)
         enddo
      enddo
      !! update flux
      do ixy = 1,Nxy
         do iz = 1,Nz
            do ig = 1,N_group
               mid = Nz*(ixy-1)+Nxy*Nz*(ig-1)+iz
               ff(ixy,iz,ig) = z(int(mid))
            enddo
         enddo
      enddo

      END SUBROUTINE AINV_SGS
#endif 

#ifdef jr_vver
       SUBROUTINE AxBHex(b,ab)
       USE Inc_TPEN

       IMPLICIT NONE
! Matrix-vector product for tpen
       INTEGER(4) :: k, kp1, km1, l, is, ln
       REAL(8), INTENT(in) :: b(0:nxy+1,nz,n_group)
       REAL(8), INTENT(out) :: ab(nxy,nz,n_group)
!
!! INDEXING
!  b = FLUX(Nxy,Nz,N_group), ab = FLUX_add(Nxy,Nz,N_group)
!  Sample::  b(imap(l),k,1) or b(imap(l),k,2)
!
!

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [AxBHex] in Mod_SolLS'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif
       DO k=1,nz
          kp1=k+1
          km1=k-1
          IF(kp1.GT.nz) kp1=nz
          IF(km1.LT.1) km1=1
          DO l=1,nxy
             ab(imap(l),k,1)=dcmat(1,l,k)*b(imap(l),k,1)+dcmat(2,l,k)*b(imap(l),k,2) &
                  -cmat(1,7,l,k)*b(imap(l),km1,1)-cmat(1,8,l,k)*b(imap(l),kp1,1)
             ab(imap(l),k,2)=dcmat(3,l,k)*b(imap(l),k,1)+dcmat(4,l,k)*b(imap(l),k,2) &
                  -cmat(2,7,l,k)*b(imap(l),km1,2)-cmat(2,8,l,k)*b(imap(l),kp1,2)
             DO is=1,ipntr(0,l)
                ln=ipntr(is,l)
                ab(imap(l),k,1)=ab(imap(l),k,1)-cmat(1,is,l,k)*b(imap(ln),k,1)
                ab(imap(l),k,2)=ab(imap(l),k,2)-cmat(2,is,l,k)*b(imap(ln),k,2)
             ENDDO
          ENDDO
       ENDDO
 !write(*,*) 'dcmat: ', dcmat(1,1,5), dcmat(2,1,5), dcmat(3,1,5), dcmat(4,1,5)
 !write(*,*) 'cmat: ', cmat(1,7,1,5), cmat(1,8,1,5), cmat(2,7,1,5), cmat(2,8,1,5)
 
       RETURN
       END SUBROUTINE AxBHex


       SUBROUTINE MinvHex(b,x_mid)
       USE Inc_TPEN
       USE Inc_Lscoef
!
       IMPLICIT NONE
! solve Mx=b for x given b
       REAL(8)    :: x_mid(0:nxy+1,nz,n_group),b(nxy,nz,n_group)
       INTEGER(4) :: l, k, kp1
       !REAL(8),ALLOCATABLE,DIMENSION(:,:) :: dumrs_mid
       REAL(8),ALLOCATABLE,DIMENSION(:,:,:) :: x 

       !if (.not.allocated(dumrs_mid)) allocate(dumrs_mid(nxy,n_group))
       if (.not.allocated(x)) allocate(x(n_group,0:nxy+1,nz))
       x=0d0
       !dumrs_mid=0d0

!

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [MinvHex] in Mod_SolLS'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

!
!! Indexing
       do k=1,nz
          do l=1,nxy
             x(1,l,k)=x_mid(imap(l),k,1)
             x(2,l,k)=x_mid(imap(l),k,2)
          enddo
       enddo
!
! forward solve
!
!   for bottom plane
       DO l=1,nxy
          dumrv(1,l)=b(imap(l),1,1)
          dumrv(2,l)=b(imap(l),1,2)
       ENDDO
       ! for debug :: r2 origin; delinv, x 
       !! delinv(1,l,k)
       !! delinv(2,l,k)
       !! delinv(3,l,k)
       !! delinv(4,l,k)
       CALL sol2dhex(delinv(:,:,1),xlufac(:,:,:,1),dumrv,x(:,1:,1))
!
       DO k=2,nz
          DO l=1,nxy
             dumrv(1,l)=b(imap(l),k,1)+cmat(1,7,l,k)*x(1,l,k-1)
             dumrv(2,l)=b(imap(l),k,2)+cmat(2,7,l,k)*x(2,l,k-1)
          ENDDO
!
          CALL sol2dhex(delinv(:,:,k),xlufac(:,:,:,k),dumrv,x(:,1:,k))
!
       ENDDO
!
! backward solve
       DO k=nz-1,1,-1
          kp1=k+1
          DO l=1,nxy
             dumrv(1,l)=-cmat(1,8,l,k)*x(1,l,kp1)
             dumrv(2,l)=-cmat(2,8,l,k)*x(2,l,kp1)
          ENDDO
!
          CALL sol2dhex(delinv(:,:,k),xlufac(:,:,:,k),dumrv,dumrs)
!
          DO l=1,nxy
             x(1,l,k)=x(1,l,k)-dumrs(1,l)
             x(2,l,k)=x(2,l,k)-dumrs(2,l)
          ENDDO
       ENDDO
       !! for indexing
       do k=1,nz
          do l=1,nxy
             x_mid(l,k,1) = x(1,imapsol(l),k)
             x_mid(l,k,2) = x(2,imapsol(l),k)
          enddo
       enddo
       !! --- 
       !stop
      
       deallocate(x)


       RETURN
       END SUBROUTINE MinvHex

       SUBROUTINE sol2dhex(df,xlu,b,x)
       !! Indexing
       !! 1) df(4,nxy)
       !! 2) xlu(n_group,6,nxy)
       !! 3) b(n_group,nxy)
       !! 4) x(n_group,nxy)
       USE Inc_TPEN
       USE Inc_Lscoef
       USE Inc_BiCG

       IMPLICIT NONE
! solve a 2d problem using precalculated LU factors
!  df : diagonal term D of preconditioner
!  xlu : off-diagonal term
       REAL(8) :: df(4,nxy),xlu(n_group,6,nassy)
       REAL(8) :: x(n_group,nxy),b(n_group,nxy)
       REAL(8) :: sb(n_group)
       INTEGER(4) :: l, j, ifr
!  forward solve

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [sol2dhex] in Mod_SolLS'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

!  forward solve
       DO l=1,nxy
          sb(1)=0
          sb(2)=0
          DO j=1,ilubnd(l)
             ifr=ipntr(j,l)
             sb(1)=sb(1)+xlu(1,j,l)*b(1,ifr)
             sb(2)=sb(2)+xlu(2,j,l)*b(2,ifr)
          ENDDO
          b(1,l)=b(1,l)-sb(1)
          b(2,l)=b(2,l)-sb(2)
       ENDDO
!
!  backward solve
       DO l=nxy,1,-1
          sb(1)=b(1,l)
          sb(2)=b(2,l)
          DO j=ilubnd(l)+1,ipntr(0,l)
             ifr=ipntr(j,l)
             sb(1)=sb(1)-xlu(1,j,l)*x(1,ifr)
             sb(2)=sb(2)-xlu(2,j,l)*x(2,ifr)
          ENDDO
          x(1,l)=df(1,l)*sb(1)+df(2,l)*sb(2)
          x(2,l)=df(3,l)*sb(1)+df(4,l)*sb(2)
       ENDDO
!
       RETURN
       END SUBROUTINE sol2dhex

#endif

#ifdef tuan_tr_test
      SUBROUTINE BicgStab_HEX(iin,checkfiss)

      USE Inc_Lscoef
      USE Inc_FluxVar
      USE Inc_BiCG
      use Inc_ExtSrc, only: flag_ExtSrc
!#ifdef jr_tr_mg
!      use inc_tpen, only: imap
!#endif
      IMPLICIT NONE

      ! solves a system of linear equations by preconditioned conjugate gradient method
      INTEGER, INTENT(out) :: iin
      INTEGER :: k, l, Ixy, Iz
      REAL(8) :: errl2t,r2t,err,psipsidt,errlinft
      REAL(8) :: r20epserf, crhod
      REAL(8) :: r0v, pts, ptt, vol
      LOGICAL, INTENT(in) :: checkfiss
!#ifdef jr_tr_mg
!      integer(4) :: ig
!#endif

      call AxBHex( Flux, flux_add )
      r20 = 0
      b2  = 0
      DO Iz = 1, Nz
         DO Ixy = 1, Nxy
            vr0(Ixy, Iz, 1) = src(1, Ixy, Iz) - flux_add(Ixy, Iz, 1)
            vr0(Ixy, Iz, 2) = src(2, Ixy, Iz) - flux_add(Ixy, Iz, 2)
            vr (Ixy, Iz, 1) = vr0(Ixy, Iz, 1)
            vr (Ixy, Iz, 2) = vr0(Ixy, Iz, 2)
            r20 = r20 + vr0(Ixy, Iz, 1)**2 + vr0(Ixy, Iz, 2)**2
            b2  = b2  + src(1, Ixy, Iz)**2 + src(2, Ixy, Iz)**2
            vp(Ixy, Iz, 1) = D0
            vp(Ixy, Iz, 2) = D0
            vv(Ixy, Iz, 1) = D0
            vv(Ixy, Iz, 2) = D0
            if(isnan(vr0(Ixy, Iz, 1))) then
                write(*,*) 'vr0 NaN 1', Ixy, Iz, src(1, Ixy, Iz),flux_add(Ixy, Iz, 1)
            write(*,*) 'nufissiton',   nu_maXs_f_3D( Ixy, Iz, 1),nu_maXs_f_3D( Ixy, Iz, 2)

                stop
            endif 
            
            if(isnan(vr0(Ixy, Iz, 2))) then
                write(*,*) 'vr0 NaN 2', Ixy, Iz
            write(*,*) 'nufissiton',   nu_maXs_f_3D( Ixy, Iz, 1),nu_maXs_f_3D( Ixy, Iz, 2)

                stop
            endif
            
            if(isnan(vr(Ixy, Iz, 1))) then
                write(*,*) 'vr NaN1 ', Ixy, Iz
            write(*,*) 'nufissiton',   nu_maXs_f_3D( Ixy, Iz, 1),nu_maXs_f_3D( Ixy, Iz, 2)

                stop
            endif 
            
            if(isnan(vr(Ixy, Iz, 2))) then
                write(*,*) 'vr NaN2', Ixy, Iz
            write(*,*) 'nufissiton',   nu_maXs_f_3D( Ixy, Iz, 1),nu_maXs_f_3D( Ixy, Iz, 2)

                stop
            endif 
            
         END DO
      END DO
      r20epserf = r20 * EPS_Residual * EPS_Residual
      r20 = SQRT(r20)
      b2  = SQRT(b2)
      calpha = 1
      crho   = 1
      comega = 1
      DO Iin = 1, Iin_Max
         crhod = crho
         crho  = scalp(vr0, vr)
         if (crho == 0) then
             write(*,*) ' crho = 0',Iin
!             stop
         endif 
         
         cbeta = crho*calpha / (crhod*comega)
         DO k = 1, Nz
            DO l = 1, Nxy
               vp(l, k, 1) = vr(l, k, 1) + cbeta *( vp(l, k, 1) - comega*vv(l, k, 1) )
               vp(l, k, 2) = vr(l, k, 2) + cbeta *( vp(l, k, 2) - comega*vv(l, k, 2) )
               if(isnan( vp(l, k, 1))) then
                   write(*,*) 'vr(l, k, 1) NaN 1',k,l,Iin,vr(l, k, 1),cbeta, crho,calpha , crhod,comega
                   write(*,*) 'vs, comega,vt',  vs(l, k, 1),  comega,vt(l, k, 1)

!                   stop
               endif 
            END DO
         END DO
         CALL MInvHex(vp, vy)
         call AxBHex( vy, vv )
         r0v    = scalp(vr0, vv)
         calpha = crho / r0v
         DO k=1,Nz
            DO l=1,Nxy
               vs(l, k, 1) = vr(l, k, 1) - calpha*vv(l, k, 1)
               vs(l, k, 2) = vr(l, k, 2) - calpha*vv(l, k, 2)
            END DO
         END DO
         CALL MInvHex(vs, vz)
         call AxBHex( vz, vt )
         pts = scalp(vs, vt)
         ptt = scalp(vt, vt)
         comega = pts / ptt
!write(*,*) 'BICGSTAB information:'
!write(*,*) '     Iteration: ', Iin
!write(*,*) '     Alpha: ', crho, r0v, calpha
!write(*,*) '     Omega: ', pts , ptt,comega
!write(*,*) '     max and min vy 1: ', maxval(vy(:,:,1)), minval(vy(:,:,1)), sum(vy(:,:,1))
!write(*,*) '     max and min vy 2: ', maxval(vy(:,:,2)), minval(vy(:,:,2)), sum(vy(:,:,2))
!write(*,*) '     max and min vz 1: ', maxval(vz(:,:,1)), minval(vz(:,:,1)) , sum(vz(:,:,1))
!write(*,*) '     max and min vz 2: ', maxval(vz(:,:,2)), minval(vz(:,:,2)) , sum(vz(:,:,2))
         IF ( checkfiss ) THEN
            flagr2= .FALSE.
            flagl2= .FALSE.
            flaglinf= .FALSE.
            flagerf= .FALSE.
            errl2t=0
            r2t=0
            psipsidt=0
            errlinft=0
            DO k = 1, Nz
               DO l = 1,Nxy
                  vol = NodeVolume(l,k)
                  Flux(l, k, 1) = Flux(l, k, 1) + calpha*vy(l, k, 1) + comega*vz(l, k, 1)
                  Flux(l, k, 2) = Flux(l, k, 2) + calpha*vy(l, k, 2) + comega*vz(l, k, 2)
                  vr(l, k, 1) = vs(l, k, 1) - comega*vt(l, k, 1)
                  vr(l, k, 2) = vs(l, k, 2) - comega*vt(l, k, 2)
                  r2t = r2t + vr(l, k, 1)**2 + vr(l, k, 2)**2
                  IF ( nu_maXs_f_3D(l, k, 2) /= D0 ) THEN
                     FisSrc_Iout(l, k) = FisSrc(l,k)
                     FisSrc(l, k) = ( nu_maXs_f_3D(l, k, 1)*Flux(l, k, 1) + nu_maXs_f_3D(l, k, 2)*Flux(l, k, 2) ) * vol
                     err = FisSrc(l, k) - FisSrc_Iout(l, k)
                     errlinft =MAX( errlinft, ABS( err/FisSrc(l,k) ) )
                     errl2t   = errl2t + err**2
                     psipsidt = psipsidt + FisSrc(l, k)*FisSrc_Iout(l, k)
                  END IF
               END DO
            END DO
            errl2=SQRT(errl2t/ABS(psipsidt))
            errlinf=errlinft
            r2=SQRT(r2t)
            r2ob2=r2/b2
            r2=r2/r20
!write(*,*) "aaaa, convergence :  ", r2ob2,EPS_Residual, errl2,eps_global,errlinf, eps_local

            IF( r2ob2   <= EPS_Residual ) flagr2   = .TRUE.
            IF( errl2   <= EPS_Global   ) flagl2   = .TRUE.
            IF( errlinf <= EPS_Local    ) flaglinf = .TRUE.
            IF( r2 <= EPS_ERF ) THEN
               flagerf = .TRUE.
               RETURN
            END IF                  
         ELSE
            r2t = 0
            DO k = 1, Nz
               DO l = 1, Nxy
                  vol = NodeVolume(l, k)
                  Flux(l, k, 1) = Flux(l, k, 1) + calpha*vy(l, k, 1) + comega*vz(l, k, 1)
                  Flux(l, k, 2) = Flux(l, k, 2) + calpha*vy(l, k, 2) + comega*vz(l, k, 2)
                  vr(l, k, 1) = vs(l, k, 1) - comega*vt(l, k, 1)
                  vr(l, k, 2) = vs(l, k, 2) - comega*vt(l, k, 2)
                  r2t = r2t + vr(l, k, 1)**2 + vr(l, k, 2)**2
                  if (isnan(Flux(l,k,1))) then
                      write(*,*) 'l, k', l,k
                      write(*,*)  calpha, comega
                      write(*,*)  'vy(l, k, 1)', vy(l, k, 1)
                      write(*,*)  'vz(l, k, 1)', vz(l, k, 1)
!                      stop
                  endif
                  
!                  if(vr(l, k, 1) == 0) then
!                      write(*,*) 'vr = 0',l, k, vs(l, k, 1) , comega,vt(l, k, 1)
!                  endif
                  
               END DO
            END DO
            IF(r2t.LT.r20epserf)then
               RETURN
            endif
         END IF
      END DO
      Iin = Iin_Max

      ! Negative flux correction
      ! skip negative flux correction for transient analysis with external source (KHNP CRI)
      if (.not.flag_extsrc) then
         do k=1,Nz
            do l=1,Nxy
               if (Flux(l,k,1)<0d0) Flux(l,k,1)=1d-30
               if (Flux(l,k,2)<0d0) Flux(l,k,2)=1d-30
            enddo
         enddo
      endif

      RETURN
      END SUBROUTINE BicgStab_HEX
#endif

      END MODULE Mod_SolLS
