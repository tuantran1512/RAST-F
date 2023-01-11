
      MODULE Mod_SolLU

      USE Inc_Constant
      USE Inc_Lscoef
      USE Inc_Geometry


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

      CONTAINS

      SUBROUTINE ilufac2d

      IMPLICIT NONE

      ! perform incomplete LU factorization for the 2D coefficient matrices
      INTEGER :: k, j, l, i, jm1, ln, lnm1, lnp1


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [ilufac2d] in Mod_SolLU'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      DO k=1,Nz
         j=1
         l=0
         DO i=Ix_Start_y(j),Ix_End_y(j)
            l=l+1
            del(1,i)=am(1,l,k)
            del(2,i)=-af2(l,k)
            del(3,i)=-scat(l,k)
            del(4,i)=am(2,l,k)
            al(1,i,k)=ccw(1,l,k)
            al(2,i,k)=zero
            al(3,i,k)=zero
            al(4,i,k)=ccw(2,l,k)
            au(1,i)=cce(1,l,k)
            au(2,i)=zero
            au(3,i)=zero
            au(4,i)=cce(2,l,k)
         END DO
         DO j=2,Ny
            jm1=j-1
            ! obtain incomplete lu factor for the 1d matrix of the row
            CALL Factor1d(jm1,k)
            ! obtain the inverse of the 1d matrix
            CALL Abi1d(jm1,k)
            DO i=Ix_Start_y(j),Ix_End_y(j)
               l=nodel(i,j)
               ln=nodel(i,jm1)
               IF(ln.GT.0 .AND. ln.LE.Nxy) THEN
                  del(1,i)=  am(1,l,k)-ccn(1,l,k)*ainvd(1,i)*ccs(1,ln,k)
                  del(2,i)=  -af2(l,k)-ccn(1,l,k)*ainvd(2,i)*ccs(2,ln,k)
                  del(3,i)= -scat(l,k)-ccn(2,l,k)*ainvd(3,i)*ccs(1,ln,k)
                  del(4,i)=  am(2,l,k)-ccn(2,l,k)*ainvd(4,i)*ccs(2,ln,k)
               ELSE
                  del(1,i)=  am(1,l,k)
                  del(2,i)=  -af2(l,k)
                  del(3,i)= -scat(l,k)
                  del(4,i)=  am(2,l,k)
               END IF
               lnm1=nodel(i-1,jm1)
               IF(i.NE.Ix_Start_y(j) .AND. lnm1.GT.0 .AND. lnm1.LE.Nxy) THEN
                  al(1,l,k)=ccw(1,l,k)                                  &
                       -ccn(1,l,k)*ainvl(1,i)*ccs(1,lnm1,k)
                  al(2,l,k)=-ccn(1,l,k)*ainvl(2,i)*ccs(2,lnm1,k)
                  al(3,l,k)=-ccn(2,l,k)*ainvl(3,i)*ccs(1,lnm1,k)
                  al(4,l,k)=ccw(2,l,k)                                  &
                       -ccn(2,l,k)*ainvl(4,i)*ccs(2,lnm1,k)
               ELSE
                  al(1,l,k)=ccw(1,l,k)
                  al(2,l,k)=zero
                  al(3,l,k)=zero
                  al(4,l,k)=ccw(2,l,k)
               END IF
               lnp1=nodel(i+1,jm1)
               IF(i.NE.Ix_End_y(j) .AND. lnp1.GT.0 .AND. lnp1.LE.Nxy) THEN
                  au(1,i)=cce(1,l,k)                                    &
                       -ccn(1,l,k)*ainvu(1,i)*ccs(1,lnp1,k)
                  au(2,i)=-ccn(1,l,k)*ainvu(2,i)*ccs(2,lnp1,k)
                  au(3,i)=-ccn(2,l,k)*ainvu(3,i)*ccs(1,lnp1,k)
                  au(4,i)=cce(2,l,k)                                    &
                       -ccn(2,l,k)*ainvu(4,i)*ccs(2,lnp1,k)
               ELSE
                  au(1,i)=cce(1,l,k)
                  au(2,i)=zero
                  au(3,i)=zero
                  au(4,i)=cce(2,l,k)
               END IF
            END DO
         END DO
         ! obtain incomplete lu factor for the 1d matrix of the last row
         CALL Factor1d(Ny,k)
      END DO

      RETURN
      END SUBROUTINE ilufac2d

      SUBROUTINE abi1d(irow,k)

      IMPLICIT NONE

      ! approximate block inverse from the LU factors
      INTEGER :: irow, k
      INTEGER :: m, l, i , lp1, mp1
      REAL(8) :: al1, al2, al3, al4
      REAL(8) :: au1, au2, au3, au4


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [abi1d] in Mod_SolLU'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      m=Ix_End_y(irow)
      l=nodel(m,irow)
      DO i=1,4
         ainvd(i,m)=delinv(i,l,k)
      END DO
      DO i=Ix_End_y(irow)-1,Ix_Start_y(irow),-1
         lp1=l
         l=l-1
         mp1=m
         m=m-1
         ! lower part of the inverse
         al1=ainvd(1,mp1)*al(1,lp1,k)+ainvd(2,mp1)*al(3,lp1,k)
         al2=ainvd(1,mp1)*al(2,lp1,k)+ainvd(2,mp1)*al(4,lp1,k)
         al3=ainvd(3,mp1)*al(1,lp1,k)+ainvd(4,mp1)*al(3,lp1,k)
         al4=ainvd(3,mp1)*al(2,lp1,k)+ainvd(4,mp1)*al(4,lp1,k)
         ainvl(1,mp1)=-al1*delinv(1,l,k)-al2*delinv(3,l,k)
         ainvl(2,mp1)=-al1*delinv(2,l,k)-al2*delinv(4,l,k)
         ainvl(3,mp1)=-al3*delinv(1,l,k)-al4*delinv(3,l,k)
         ainvl(4,mp1)=-al3*delinv(2,l,k)-al4*delinv(4,l,k)
         ! upper part of the inverse
         au1=delinv(1,l,k)*au(1,m)+delinv(2,l,k)*au(3,m)
         au2=delinv(1,l,k)*au(2,m)+delinv(2,l,k)*au(4,m)
         au3=delinv(3,l,k)*au(1,m)+delinv(4,l,k)*au(3,m)
         au4=delinv(3,l,k)*au(2,m)+delinv(4,l,k)*au(4,m)
         ainvu(1,m)=-au1*ainvd(1,mp1)-au2*ainvd(3,mp1)
         ainvu(2,m)=-au1*ainvd(2,mp1)-au2*ainvd(4,mp1)
         ainvu(3,m)=-au3*ainvd(1,mp1)-au4*ainvd(3,mp1)
         ainvu(4,m)=-au3*ainvd(2,mp1)-au4*ainvd(4,mp1)
         ! diagonal part
         ainvd(1,m)=delinv(1,l,k)-au1*ainvl(1,mp1)-au2*ainvl(3,mp1)
         ainvd(2,m)=delinv(2,l,k)-au1*ainvl(2,mp1)-au2*ainvl(4,mp1)
         ainvd(3,m)=delinv(3,l,k)-au3*ainvl(1,mp1)-au4*ainvl(3,mp1)
         ainvd(4,m)=delinv(4,l,k)-au3*ainvl(2,mp1)-au4*ainvl(4,mp1)
      END DO

      END SUBROUTINE abi1d

      SUBROUTINE factor1d(irow,k)

      IMPLICIT NONE

      ! incomplete factorization of block-tridiagonal matices
      INTEGER :: irow, k, i, l, im1, lm1
      REAL(8) :: ald1, ald2, ald3, ald4, f


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [factor1d] in Mod_SolLU'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      i=Ix_Start_y(irow)
      l=nodel(i,irow)
      ! inverse of a 2x2 matrix
      f=1/(del(1,i)*del(4,i)-del(2,i)*del(3,i))
      delinv(1,l,k)=del(4,i)*f
      delinv(2,l,k)=-del(2,i)*f
      delinv(3,l,k)=-del(3,i)*f
      delinv(4,l,k)=del(1,i)*f
      ! calc. inv(del)*u for later use in backsub
      deliau(1,l,k)=delinv(1,l,k)*au(1,i)+delinv(2,l,k)*au(3,i)
      deliau(2,l,k)=delinv(1,l,k)*au(2,i)+delinv(2,l,k)*au(4,i)
      deliau(3,l,k)=delinv(3,l,k)*au(1,i)+delinv(4,l,k)*au(3,i)
      deliau(4,l,k)=delinv(3,l,k)*au(2,i)+delinv(4,l,k)*au(4,i)
      im1=i
      DO i=Ix_Start_y(irow)+1,Ix_End_y(irow)
         lm1=l
         l=l+1
         ald1=al(1,l,k)*delinv(1,lm1,k)+al(2,l,k)*delinv(3,lm1,k)
         ald2=al(1,l,k)*delinv(2,lm1,k)+al(2,l,k)*delinv(4,lm1,k)
         ald3=al(3,l,k)*delinv(1,lm1,k)+al(4,l,k)*delinv(3,lm1,k)
         ald4=al(3,l,k)*delinv(2,lm1,k)+al(4,l,k)*delinv(4,lm1,k)
         del(1,i)=del(1,i)-ald1*au(1,im1)-ald2*au(3,im1)
         del(2,i)=del(2,i)-ald1*au(2,im1)-ald2*au(4,im1)
         del(3,i)=del(3,i)-ald3*au(1,im1)-ald4*au(3,im1)
         del(4,i)=del(4,i)-ald3*au(2,im1)-ald4*au(4,im1)
         f=1/(del(1,i)*del(4,i)-del(2,i)*del(3,i))
         delinv(1,l,k)=del(4,i)*f
         delinv(2,l,k)=-del(2,i)*f
         delinv(3,l,k)=-del(3,i)*f
         delinv(4,l,k)=del(1,i)*f
         ! calc. inv(del)*u for later use in backsub
         deliau(1,l,k)=delinv(1,l,k)*au(1,i)+delinv(2,l,k)*au(3,i)
         deliau(2,l,k)=delinv(1,l,k)*au(2,i)+delinv(2,l,k)*au(4,i)
         deliau(3,l,k)=delinv(3,l,k)*au(1,i)+delinv(4,l,k)*au(3,i)
         deliau(4,l,k)=delinv(3,l,k)*au(2,i)+delinv(4,l,k)*au(4,i)
         im1=i
      END DO

      END SUBROUTINE factor1d

      END MODULE Mod_SolLU
