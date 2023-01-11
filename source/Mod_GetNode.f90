
      MODULE Mod_GetNode

      USE Inc_INP
      USE Inc_Constant
      USE Inc_Geometry
      USE Inc_RP
      USE Mod_Alloc


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE Get_RadialIndex
      !
      ! Below Description is for "1/4 Core" and "4N/1FA" Case
      ! Row    = Y-Axis, Left  = (-) Direction
      ! Column = X-Axis, Right = (+) Direction
      !
      !     X    [1N/1FA]
      !     -->
      ! Y | 1   2   3   4
      !   V
      !     5   6   7   8
      !
      !     9   10  11
      !
      !     12  13
      !
      !      X   [4N/1FA]
      !      -->
      !      --------------------------
      ! Y | |1 | 2   3 | 4   5 | 6   7 |
      !   V |--|-------|-------|-------|
      !     |8 | 9   10| 11  12| 13  14|
      !     |15| 16  17| 18  19| 20  21|
      !     |--|-------|-------|-------
      !     |22| 23  24| 25  26|
      !     |27| 28  29| 30  31|
      !     |--|-------|-------
      !     |32| 33  34|           (Ixy = 16) = (Ix = 2, Iy = 3)
      !     |35| 36  37|           I_4Nto1N(16) = 6
      !      ----------
      IMPLICIT NONE

      INTEGER :: Ixy
      INTEGER :: Ix
      INTEGER :: Iy
      INTEGER :: Ix_1N
      INTEGER :: Iy_1N
      INTEGER :: Ixy_1N
      INTEGER :: i
      INTEGER :: j
      INTEGER :: Isum
      INTEGER :: Sum_Int
      INTEGER :: Sum_Int2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_RadialIndex] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_RadialIndex] in Mod_GetNode'
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

      CALL Alloc( Nx_y, Ny )

      CALL Alloc( Ixy_Start_y , Ny )
      CALL Alloc( Ixy_End_y   , Ny )

      CALL Alloc( Ix_Start_y_1N   , Ny_1N )
      CALL Alloc( Ix_StartFA_y_1N , Ny_1N )
      CALL Alloc( Ix_End_y_1N     , Ny_1N )
      CALL Alloc( Ix_EndFA_y_1N   , Ny_1N )
      CALL Alloc( Iy_Start_x_1N   , Nx_1N )
      CALL Alloc( Iy_StartFA_x_1N , Nx_1N )
      CALL Alloc( Iy_End_x_1N     , Nx_1N )
      CALL Alloc( Iy_EndFA_x_1N   , Nx_1N )

      ! To Prevent Sweeping from 0 to 0
      !Ix_Start_y_1N   = - 1
      !Ix_StartFA_y_1N = - 1
      !Ix_End_y_1N     = - 1
      Ix_EndFA_y_1N   = - 1
      !Iy_Start_x_1N   = - 1
      !Iy_StartFA_x_1N = - 1
      !Iy_End_x_1N     = - 1
      Iy_EndFA_x_1N   = - 1

      CALL Alloc( Ix_Start_y , Ny + 1 )
      CALL Alloc( Ix_End_y   , Ny + 1 )
      CALL Alloc( Iy_Start_x , Nx + 1 )
      CALL Alloc( Iy_End_x   , Nx + 1 )


      DO Iy = 1, Ny
         Isum = 0

         DO Ix = 1, Nx
            IF ( I_FARF_2D(Ix, Iy) == 0 ) THEN
               CYCLE
            END IF

            Isum = Isum + 1
         END DO

         Nx_y(Iy) = Isum
      END DO

      Nxy  = SUM(Nx_y)
      Nxyz = Nxy*Nz

      CALL Alloc( IxIyToIxy       , Nx    , Ny    )
      CALL Alloc( IxIyToIxy_1N    , Nx    , Ny    )
      CALL Alloc( IxIy_1NToIxy    , Nx_1N , Ny_1N )
      CALL Alloc( IxIy_1NToIxy_1N , Nx_1N , Ny_1N )
      CALL Alloc( IxToIx_1N       , Nx            )
      CALL Alloc( IyToIy_1N       ,         Ny    )

      Ixy = 0

      DO Iy = 1, Ny
         DO Iy_1N = 1, Ny_1N
            IF ( Iy <= SUM( Mesh_y(1:Iy_1N) ) ) THEN
               EXIT
            END IF
         END DO

         IyToIy_1N(Iy) = Iy_1N

         IF ( Iy == 1 ) THEN
            DO Ix = 1, Nx
               DO Ix_1N = 1, Nx_1N
                  IF ( Ix <= SUM( Mesh_x(1:Ix_1N) ) ) THEN
                     EXIT
                  END IF
               END DO

               IxToIx_1N(Ix) = Ix_1N
            END DO
         END IF

         DO Ix = 1, Nx
            IF ( I_FARF_2D(Ix, Iy) == 0 ) THEN
               CYCLE
            END IF

            Ixy   = Ixy + 1
            Ix_1N = IxToIx_1N(Ix)

            Isum = Ix_1N

            DO i = 1, Ix_1N
               IF ( I_FARF_1N_2D(i, Iy_1N) == 0 ) THEN
                  Isum = Isum - 1
               END IF
            END DO

            IxIyToIxy   (Ix, Iy) = Ixy
            IxIyToIxy_1N(Ix, Iy) = SUM( Nx_y_1N( 1 : Iy_1N - 1 ) ) + Isum

            IxIy_1NToIxy   (Ix_1N, Iy_1N) = IxIyToIxy   (Ix, Iy)
            IxIy_1NToIxy_1N(Ix_1N, Iy_1N) = IxIyToIxy_1N(Ix, Iy)
         END DO
      END DO

      call alloc( ixy_1ntoix_1n, nxy_1n)
      call alloc( ixy_1ntoiy_1n, nxy_1n)
      do iy_1n=1,ny_1n
         do ix_1n=1,nx_1n
            if (i_farf_1n_2d(ix_1n,iy_1n)<=0) cycle
            ixy_1n=ixiy_1ntoixy_1n(ix_1n,iy_1n)
            ixy_1ntoix_1n(ixy_1n)=ix_1n
            ixy_1ntoiy_1n(ixy_1n)=iy_1n
         enddo
      enddo

      call alloc(ixytoix,nxy)
      ixytoix=0
      call alloc(ixytoiy,nxy)
      ixytoiy=0
      do iy=1,ny
         do ix=1,nx
            if (i_farf_2d(ix,iy)==0) cycle
               ixy=ixiytoixy(ix,iy)
               ixytoix(ixy)=ix
               ixytoiy(ixy)=iy
         enddo
      enddo

      Isum = 0

      DO Iy = 1, Ny
         Ixy_Start_y(Iy) = Isum + 1
         Isum            = Isum + Nx_y(Iy)
         Ixy_End_y(Iy)   = Isum
      END DO

      DO Iy_1N = 1, Ny_1N
         DO Ix_1N = 1, Nx_1N
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) /= 0 ) THEN
               Ix_Start_y_1N(Iy_1N) = Ix_1N

               EXIT
            END IF
         END DO

         DO Ix_1N = Ix_Start_y_1N(Iy_1N), Nx_1N
#ifdef tuan_fr_crm
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) == 2 ) THEN
               Ix_StartFA_y_1N(Iy_1N) = Ix_1N
#else
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) > 1 ) THEN
               Ix_StartFA_y_1N(Iy_1N) = Ix_1N
#endif
               
               EXIT
            END IF
         END DO

         DO Ix_1N = Nx_1N, 1, - 1
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) /= 0 ) THEN
               Ix_End_y_1N(Iy_1N) = Ix_1N

               EXIT
            END IF
         END DO

         DO Ix_1N = Ix_End_y_1N(Iy_1N), 1, - 1
#ifdef tuan_fr_crm
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) == 2) THEN
               Ix_EndFA_y_1N(Iy_1N) = Ix_1N

#else
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) > 1) THEN
               Ix_EndFA_y_1N(Iy_1N) = Ix_1N

#endif

               EXIT
            END IF
         END DO
      END DO

      DO Ix_1N = 1, Nx_1N
         DO Iy_1N = 1, Ny_1N
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) /= 0 ) THEN
               Iy_Start_x_1N(Ix_1N) = Iy_1N

               EXIT
            END IF
         END DO

         DO Iy_1N = Iy_Start_x_1N(Ix_1N), Ny_1N
#ifdef tuan_fr_crm
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) == 2 ) THEN
               Iy_StartFA_x_1N(Ix_1N) = Iy_1N
#else
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) > 1 ) THEN
               Iy_StartFA_x_1N(Ix_1N) = Iy_1N
#endif


               EXIT
            END IF
         END DO

         DO Iy_1N = Ny_1N, 1, - 1
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) /= 0 ) THEN
               Iy_End_x_1N(Ix_1N) = Iy_1N

               EXIT
            END IF
         END DO

         DO Iy_1N = Iy_End_x_1N(Ix_1N), 1, - 1
#ifdef tuan_fr_crm
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) == 2) THEN
               Iy_EndFA_x_1N(Ix_1N) = Iy_1N
#else
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) > 1) THEN
               Iy_EndFA_x_1N(Ix_1N) = Iy_1N
#endif
               EXIT
            END IF
         END DO
      END DO

      Iy = 0

      DO Iy_1N = 1, Ny_1N
         Sum_Int  = 0
         Sum_Int2 = 0

         DO Ix_1N = 1, ( Ix_Start_y_1N(Iy_1N) - 1 )
            Sum_Int = Sum_Int + Mesh_x(Ix_1N)
         END DO

         DO Ix_1N = Nx_1N, ( Ix_End_y_1N(Iy_1N) + 1 ), - 1
            Sum_Int2 = Sum_Int2 + Mesh_x(Ix_1N)
         END DO

         DO i = 1, Mesh_y(Iy_1N)
            Iy = Iy + 1

            Ix_Start_y(Iy) = Sum_Int + 1
            Ix_End_y  (Iy) = Nx - Sum_Int2
         END DO
      END DO

      Ix_Start_y(Ny + 1) = Ix_Start_y(Ny)
      Ix_End_y  (Ny + 1) = Ix_End_y  (Ny) !js+ + 1

      Ix = 0

      DO Ix_1N = 1, Nx_1N
         Sum_Int  = 0
         Sum_Int2 = 0

         DO Iy_1N = 1, ( Iy_Start_x_1N(Ix_1N) - 1 )
            Sum_Int = Sum_Int + Mesh_y(Iy_1N)
         END DO

         DO Iy_1N = Ny_1N, ( Iy_End_x_1N(Ix_1N) + 1 ), - 1
            Sum_Int2 = Sum_Int2 + Mesh_y(Iy_1N)
         END DO

         DO j = 1, Mesh_x(Ix_1N)
            Ix = Ix + 1

            Iy_Start_x(Ix) = Sum_Int + 1
            Iy_End_x  (Ix) = Ny - Sum_Int2
         END DO
      END DO

      Iy_Start_x(Nx + 1) = Iy_Start_x(Nx)
      Iy_End_x  (Nx + 1) = Iy_End_x  (Nx) + 1

      RETURN
      END SUBROUTINE Get_RadialIndex


      SUBROUTINE Get_AdjacentNodeIndex

      IMPLICIT NONE

      INTEGER :: Ixy, Ix, Iy, Iz

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_AdjacentNodeIndex] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_AdjacentNodeIndex] in Mod_GetNode'
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

      CALL Alloc( I_Rx, Nxy )
      CALL Alloc( I_Lx, Nxy )
      CALL Alloc( I_Ry, Nxy )
      CALL Alloc( I_Ly, Nxy )
      CALL Alloc( I_Rz, Nz  )
      CALL Alloc( I_Lz, Nz  )

      IF ( Nxy == 1 ) THEN
         GOTO 999
      END IF

      Ix  = 1
      Iy  = 1
      Ixy = IxIyToIxy(Ix, Iy)

      IF ( Ixy /= 0 ) THEN
         IF ( Nx >= 2 ) THEN
            I_Rx(Ixy) = IxIyToIxy(Ix + 1, Iy)
         END IF

         IF ( Ny >= 2 ) THEN
            I_Ry(Ixy) = IxIyToIxy(Ix, Iy + 1)
         END IF
      END IF

      Ix  = Nx
      Iy  = 1
      Ixy = IxIyToIxy(Ix, Iy)

      IF ( Ixy /= 0 ) THEN
         IF ( Nx >= 2 ) THEN
            I_Lx(Ixy) = IxIyToIxy(Ix - 1, Iy)
         END IF

         IF ( Ny >= 2 ) THEN
            I_Ry(Ixy) = IxIyToIxy(Ix, Iy + 1)
         END IF
      END IF

      Ix  = 1
      Iy  = Ny
      Ixy = IxIyToIxy(Ix, Iy)

      IF ( Ixy /= 0 ) THEN
         IF ( Nx >= 2 ) THEN
            I_Rx(Ixy) = IxIyToIxy(Ix + 1, Iy)
         END IF

         IF ( Ny >= 2 ) THEN
            I_Ly(Ixy) = IxIyToIxy(Ix, Iy - 1)
         END IF
      END IF

      Ix  = Nx
      Iy  = Ny
      Ixy = IxIyToIxy(Ix, Iy)

      IF ( Ixy /= 0 ) THEN
         IF ( Nx >= 2 ) THEN
            I_Lx(Ixy) = IxIyToIxy(Ix - 1, Iy)
         END IF

         IF ( Ny >= 2 ) THEN
            I_Ly(Ixy) = IxIyToIxy(Ix, Iy - 1)
         END IF
      END IF

      Ix = 1

      IF ( Ny >= 3 ) THEN
         DO Iy = 2, Ny - 1
            Ixy = IxIyToIxy(Ix, Iy)

            IF ( Ixy /= 0 ) THEN
               I_Rx(Ixy) = IxIyToIxy(Ix + 1, Iy    )
               I_Ly(Ixy) = IxIyToIxy(Ix    , Iy - 1)
               I_Ry(Ixy) = IxIyToIxy(Ix    , Iy + 1)
            END IF
         END DO
      END IF

      Ix = Nx

      IF ( Ny >= 3 ) THEN
         DO Iy = 2, Ny - 1
            Ixy = IxIyToIxy(Ix, Iy)

            IF ( Ixy /= 0 ) THEN
               I_Lx(Ixy) = IxIyToIxy(Ix - 1, Iy    )
               I_Ly(Ixy) = IxIyToIxy(Ix    , Iy - 1)
               I_Ry(Ixy) = IxIyToIxy(Ix    , Iy + 1)
            END IF
         END DO
      END IF

      Iy = 1

      IF ( Nx >= 3 ) THEN
         DO Ix = 2, Nx - 1
            Ixy = IxIyToIxy(Ix, Iy)

            IF ( Ixy /= 0 ) THEN
               I_Lx(Ixy) = IxIyToIxy(Ix - 1, Iy    )
               I_Rx(Ixy) = IxIyToIxy(Ix + 1, Iy    )
               I_Ry(Ixy) = IxIyToIxy(Ix    , Iy + 1)
            END IF
         END DO
      END IF

      Iy = Ny

      IF ( Nx >= 3 ) THEN
         DO Ix = 2, Nx - 1
            Ixy = IxIyToIxy(Ix, Iy)

            IF ( Ixy /= 0 ) THEN
               I_Lx(Ixy) = IxIyToIxy(Ix - 1, Iy    )
               I_Rx(Ixy) = IxIyToIxy(Ix + 1, Iy    )
               I_Ly(Ixy) = IxIyToIxy(Ix    , Iy - 1)
            END IF
         END DO
      END IF

      IF ( ( Nx >= 3 ) .AND. ( Ny >= 3 ) ) THEN
         DO Iy = 2, Ny - 1
            DO Ix = 2, Nx - 1
               Ixy = IxIyToIxy(Ix, Iy)

               IF ( Ixy /= 0 ) THEN
                  I_Lx(Ixy) = IxIyToIxy(Ix - 1, Iy    )
                  I_Rx(Ixy) = IxIyToIxy(Ix + 1, Iy    )
                  I_Ly(Ixy) = IxIyToIxy(Ix    , Iy - 1)
                  I_Ry(Ixy) = IxIyToIxy(Ix    , Iy + 1)
               END IF
            END DO
         END DO
      END IF

      DO Iz = 1, Nz
         I_Rz(Iz) = Iz + 1
         I_Lz(Iz) = Iz - 1
      END DO

      I_Lz(1)  = 0
      I_Rz(Nz) = 0

999   RETURN
      END SUBROUTINE Get_AdjacentNodeIndex


      SUBROUTINE Get_NodeID

      IMPLICIT NONE

      INTEGER :: Ixy, Ix, Iy

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_NodeID] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_NodeID] in Mod_GetNode'
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

      CALL Alloc( N_ID , Nxy )

      IF ( Nxy == 1 ) THEN
        N_ID(1) = - 5

        GOTO 999
      END IF

      DO Iy = 1, Ny
         DO Ix = Ix_Start_y(Iy), Ix_End_y(Iy)
            Ixy = IxIyToIxy(Ix, Iy)

            IF ((I_Lx(Ixy)==0).AND.(I_Rx(Ixy)/=0).AND.(I_Ly(Ixy)==0).AND.(I_Ry(Ixy)/=0)) THEN
               N_ID(Ixy) = 1
            ELSE IF ((I_Lx(Ixy)/=0).AND.(I_Rx(Ixy)/=0).AND.(I_Ly(Ixy)==0).AND.(I_Ry(Ixy)/=0)) THEN
               N_ID(Ixy) = 2
            ELSE IF ((I_Lx(Ixy)/=0).AND.(I_Rx(Ixy)==0).AND.(I_Ly(Ixy)==0).AND.(I_Ry(Ixy)/=0)) THEN
               N_ID(Ixy) = 3
            ELSE IF ((I_Lx(Ixy)==0).AND.(I_Rx(Ixy)/=0).AND.(I_Ly(Ixy)/=0).AND.(I_Ry(Ixy)/=0)) THEN
               N_ID(Ixy) = 4
            ELSE IF ((I_Lx(Ixy)/=0).AND.(I_Rx(Ixy)/=0).AND.(I_Ly(Ixy)/=0).AND.(I_Ry(Ixy)/=0)) THEN
               N_ID(Ixy) = 5
            ELSE IF ((I_Lx(Ixy)/=0).AND.(I_Rx(Ixy)==0).AND.(I_Ly(Ixy)/=0).AND.(I_Ry(Ixy)/=0)) THEN
               N_ID(Ixy) = 6
            ELSE IF ((I_Lx(Ixy)==0).AND.(I_Rx(Ixy)/=0).AND.(I_Ly(Ixy)/=0).AND.(I_Ry(Ixy)==0)) THEN
               N_ID(Ixy) = 7
            ELSE IF ((I_Lx(Ixy)/=0).AND.(I_Rx(Ixy)/=0).AND.(I_Ly(Ixy)/=0).AND.(I_Ry(Ixy)==0)) THEN
               N_ID(Ixy) = 8
            ELSE IF ((I_Lx(Ixy)/=0).AND.(I_Rx(Ixy)==0).AND.(I_Ly(Ixy)/=0).AND.(I_Ry(Ixy)==0)) THEN
               N_ID(Ixy) = 9
            ELSE IF ((I_Lx(Ixy)/=0).AND.(I_Rx(Ixy)==0).AND.(I_Ly(Ixy)/=0).AND.(I_Ry(Ixy)==0)) THEN
               N_ID(Ixy) = - 1
            ELSE IF ((I_Lx(Ixy)==0).AND.(I_Rx(Ixy)==0).AND.(I_Ly(Ixy)/=0).AND.(I_Ry(Ixy)==0)) THEN
               N_ID(Ixy) = - 2
            ELSE IF ((I_Lx(Ixy)==0).AND.(I_Rx(Ixy)/=0).AND.(I_Ly(Ixy)/=0).AND.(I_Ry(Ixy)==0)) THEN
               N_ID(Ixy) = - 3
            ELSE IF ((I_Lx(Ixy)/=0).AND.(I_Rx(Ixy)==0).AND.(I_Ly(Ixy)==0).AND.(I_Ry(Ixy)==0)) THEN
               N_ID(Ixy) = - 4
            ELSE IF ((I_Lx(Ixy)==0).AND.(I_Rx(Ixy)==0).AND.(I_Ly(Ixy)==0).AND.(I_Ry(Ixy)==0)) THEN
               N_ID(Ixy) = - 5
            ELSE IF ((I_Lx(Ixy)==0).AND.(I_Rx(Ixy)/=0).AND.(I_Ly(Ixy)==0).AND.(I_Ry(Ixy)==0)) THEN
               N_ID(Ixy) = - 6
            ELSE IF ((I_Lx(Ixy)/=0).AND.(I_Rx(Ixy)==0).AND.(I_Ly(Ixy)==0).AND.(I_Ry(Ixy)/=0)) THEN
               N_ID(Ixy) = - 7
            ELSE IF ((I_Lx(Ixy)==0).AND.(I_Rx(Ixy)==0).AND.(I_Ly(Ixy)==0).AND.(I_Ry(Ixy)/=0)) THEN
               N_ID(Ixy) = - 8
            ELSE IF ((I_Lx(Ixy)==0).AND.(I_Rx(Ixy)/=0).AND.(I_Ly(Ixy)==0).AND.(I_Ry(Ixy)/=0)) THEN
               N_ID(Ixy) = - 9
            ELSE IF ((I_Lx(Ixy)==0).AND.(I_Rx(Ixy)==0).AND.(I_Ly(Ixy)/=0).AND.(I_Ry(Ixy)/=0)) THEN
               N_ID(Ixy) = 456
            ELSE IF ((I_Lx(Ixy)/=0).AND.(I_Rx(Ixy)/=0).AND.(I_Ly(Ixy)==0).AND.(I_Ry(Ixy)==0)) THEN
               N_ID(Ixy) = 258
            END IF
         END DO
      END DO

999   RETURN
      END SUBROUTINE Get_NodeID


      SUBROUTINE Get_I_4Nto1N

      IMPLICIT NONE

      INTEGER :: Ixy, Ix, Iy, Ix_1N, Iy_1N

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_I_4Nto1N] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_I_4Nto1N] in Mod_GetNode'
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

      CALL Alloc( I_4Nto1N , Nxy )

      DO Iy = 1, Ny
         Iy_1N = IyToIy_1N(Iy)
         DO Ix = Ix_Start_y(Iy), Ix_End_y(Iy)
            Ix_1N = IxToIx_1N(Ix)
            Ixy   = IxIyToIxy(Ix, Iy)
            I_4Nto1N(Ixy) = IxIyToIxy_1N(Ix, Iy)
         END DO
      END DO

      RETURN
      END SUBROUTINE Get_I_4Nto1N


      SUBROUTINE Get_I_1Nto4N
      !      -------
      ! Y | | 1 | 2 |   1: Lx + Ly   2: Rx + Ly
      !   V |---|---|
      !     | 3 | 4 |   3: Lx + Ry   4: Rx + Ry
      !      -------

      IMPLICIT NONE

      INTEGER :: Ixy, Ix, Iy, Ixy_1N, Ix_1N, Iy_1N
      INTEGER :: i, j
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_4Nin1FA

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_I_1Nto4N] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_I_1Nto4N] in Mod_GetNode'
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

      CALL Alloc( I_1Nto4N , Nxy_1N , 4 )
      CALL Alloc( Ix_4Nto1N , Nx )
      CALL Alloc( Iy_4Nto1N , Ny )
      CALL Alloc( I_4Nin1FA , Nxy_1N )

      DO Iy = 1, Ny
         Iy_1N = IyToIy_1N(Iy)
         DO Ix = Ix_Start_y(Iy), Ix_End_y(Iy)
            Ix_1N  = IxToIx_1N(Ix)
            Ixy    = IxIyToIxy(Ix, Iy)
            Ixy_1N = IxIyToIxy_1N(Ix, Iy)
            IF ( OPT_Core == 4 ) THEN
               IF ( Nxy == 1) THEN
                  I_4Nin1FA(1) = 4
               ELSE IF ( (Ix == 1)        .AND. (Iy ==1)         .AND.  &
                         (Mesh_x(1) == 1) .AND. (Mesh_x(2) == 2) .AND.  &
                         (GridSize_x(1) / GridSize_x(2) == 0.5D0 )     ) THEN
                  I_4Nin1FA(Ixy_1N) = 4  ! Ixy_1N == 1
                  Flag_HalfCenter = .TRUE.
               ELSE IF ( (Flag_HalfCenter .EQV. .TRUE.) .AND. (Iy ==1) ) THEN
                  IF ( I_4Nin1FA(Ixy_1N) == 0 ) THEN
                     I_4Nin1FA(Ixy_1N) = 3
                  ELSE
                     I_4Nin1FA(Ixy_1N) = 4
                  END IF
               ELSE IF ( (Flag_HalfCenter .EQV. .TRUE.) .AND. (Ix ==1) ) THEN
                  IF ( I_4Nin1FA(Ixy_1N) == 0 ) THEN
                     I_4Nin1FA(Ixy_1N) = 2
                  ELSE
                     I_4Nin1FA(Ixy_1N) = 4
                  END IF
               ELSE
                  I_4Nin1FA(Ixy_1N) = I_4Nin1FA(Ixy_1N) + 1
               END IF
            ELSE
               I_4Nin1FA(Ixy_1N) = I_4Nin1FA(Ixy_1N) + 1
            END IF
            I_1Nto4N(Ixy_1N, I_4Nin1FA(Ixy_1N)) = Ixy
         END DO
      END DO

      Ix = 0
      Iy = 0

      DO Ix_1N = 1, Nx_1N
         DO i = 1, Mesh_x(Ix_1N)
            Ix = Ix + 1
            Ix_4Nto1N(Ix) = Ix_1N
         END DO
      END DO
      DO Iy_1N = 1, Ny_1N
         DO j = 1, Mesh_y(Iy_1N)
            Iy = Iy + 1
            Iy_4Nto1N(Iy) = Iy_1N
         END DO
      END DO

      RETURN
      END SUBROUTINE Get_I_1Nto4N


      SUBROUTINE Get_FAandRF_Index
      USE Inc_CR, ONLY: I_CR_1N, I_CR_4N
      USE Inc_PinPOW
      USE Inc_XS_File

#ifdef tuan_fr
      use Inc_Option,        only: n_group
#endif
#ifdef js_mpi
      use inc_parallel, only: comm, iproc, node2iproc
      use mod_parallel, only: allreduce
#endif
      IMPLICIT NONE
      INTEGER :: Ixy, Ixy_FA, Ixy_RF, Ixy_1N
      integer :: Ixy_FA_1N, Ix_1N, Iy_1N
      integer :: Iz, I_LP
      INTEGER :: m, Icount, i, ii
      integer :: ibot, itop
      real(8) :: htop, hbot
      LOGICAL(1) :: Flag_LoopExit
      LOGICAL(1), SAVE :: Flag_First

#ifdef js_mpi
      integer, allocatable :: tmp_farf_1n_3d(:,:)
#endif

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_FAandRF_Index] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_FAandRF_Index] in Mod_GetNode'
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

      IF ( Flag_First .EQV. .FALSE. ) THEN
         Flag_First = .TRUE.

         DO Iy_1N = 1, Ny_1N
            DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
#ifdef tuan_fr_crm                
               IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) /= 2 ) THEN
#else
               IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) == 1 ) THEN
#endif
                  Nxy_RF_1N = Nxy_RF_1N + 1
               END IF
            END DO
         END DO

         Nxy_FA_1N = Nxy_1N - Nxy_RF_1N

         CALL Alloc( I_CR_1N , Nxy_1N )
         CALL Alloc( I_CR_4N , Nxy    )
         CALL Alloc( I_FA_1N , Nxy_FA_1N )
         CALL Alloc( I_FARF_1N , Nxy_1N )

         IF (Nxy_RF_1N >= 1) THEN
            CALL Alloc( I_RF_1N , Nxy_RF_1N )
         ELSE
            CALL Alloc( I_RF_1N , Nxy_FA_1N )
         END IF

         IF ( OPT_Core == 4 ) THEN
            IF ( Flag_HalfCenter ) THEN
               Core_N_FARF  = Nxy_1N   * 4 - 3 &
                  & - 2*( Ix_End_y_1N  (1) - Ix_Start_y_1N  (1) ) &
                  & - 2*( Iy_End_x_1N  (1) - Iy_Start_x_1N  (1) )
               Core_N_FA   = Nxy_FA_1N * 4 - 3 &
                  & - 2*( Ix_EndFA_y_1N(1) - Ix_StartFA_y_1N(1) ) &
                  & - 2*( Iy_EndFA_x_1N(1) - Iy_StartFA_x_1N(1) )
            ELSE
               Core_N_FARF = Nxy_1N    * 1
               Core_N_FA   = Nxy_FA_1N * 1
            END IF
         ELSE
            Core_N_FARF = Nxy_1N
            Core_N_FA   = Nxy_FA_1N
         END IF

         i  = 0
         ii = 0
         DO Iy_1N = 1, Ny_1N
            DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
               Ixy_1N = IxIy_1NToIxy_1N(Ix_1N, Iy_1N)
               IF ( Ixy_1N == 0 ) CYCLE
#ifdef tuan_fr_crm
               IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) /= 2 ) THEN
#else
               IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) == 1 ) THEN
#endif
                  i = i + 1
                  I_RF_1N(i)        = Ixy_1N
                  I_FARF_1N(Ixy_1N) = 1
#ifdef tuan_fr_crm
               ELSE IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) == 2 ) THEN
#else
               ELSE IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) >= 2 ) THEN
#endif
                  ii = ii + 1
                  I_FA_1N(ii)       = Ixy_1N
                  I_FARF_1N(Ixy_1N) = 2
               END IF
            END DO
         END DO

         DO Ixy_FA_1N = 1, Nxy_FA_1N
            Icount = 0
            DO m = 1, 4
               IF (I_1Nto4N(I_FA_1N(Ixy_FA_1N), m) == 0) THEN
                  Icount = Icount + 1
               END IF
            END DO
            Nxy_FA = Nxy_FA + 4 - Icount
         END DO

         Nxy_RF = Nxy - Nxy_FA

         CALL Alloc( I_FA , Nxy_FA )
         CALL Alloc( I_RF , Nxy_RF )

         Ixy_FA = 0
         Ixy_RF = 0
         DO Ixy = 1, Nxy
            Ixy_1N = I_4Nto1N(Ixy)
#ifdef tuan_fr_crm
            IF ( I_FARF_1N(Ixy_1N) == 2 ) THEN
#else
            IF ( I_FARF_1N(Ixy_1N) >= 2 ) THEN
#endif
               Ixy_FA = Ixy_FA + 1
               I_FA(Ixy_FA) = Ixy
            ELSE
               Ixy_RF = Ixy_RF + 1
               I_RF(Ixy_RF) = Ixy
            END IF
         END DO

!         call alloc(ixytoifa,nxy)
 !        Allocate (ixytoifa(1:nxy))
 !        ixytoifa=0
 !        do ixy=1,nxy
 !           do ixy_fa=1,nxy_fa
 !              if (ixy==i_fa(ixy_fa)) then
 !                 ixytoifa(ixy)=ixy_fa
 !              endif
 !           enddo
 !        enddo

      ELSE
         CALL Alloc( I_FARF_1N_3D , Nxy_1N , Nz )
         I_SP = 1
         I_BU = 1
         DO Ixy = 1, Nxy
            Ixy_1N = I_4Nto1N(Ixy)
            I_LP   = I_LP_1N( Ixy_1N )
            IF ( I_FARF_1N( Ixy_1N ) == 1 ) THEN
               I_FARF_1N_3D( Ixy_1N, : ) = 1
            ELSE
               DO Iz = 1, Nz
                  I_Tab = AxialComp( I_LP, Iz )
                  I_Tab = new_asym_itab(I_Tab,Ixy)
                  IF ( Flag_Card_maXS ) THEN
#ifdef tuan_fr
                         IF ( sum(Ref_maXS( I_Tab, (4*n_group+1):(5*n_group))) /= D0 )  THEN
                            I_FARF_1N_3D( Ixy_1N, Iz ) = 2
                         ELSE
                            I_FARF_1N_3D( Ixy_1N, Iz ) = 1
                         END IF
#else
                     IF ( ( Ref_maXS( I_Tab, 7) /= D0 ) .OR.  &
                          ( Ref_maXS( I_Tab, 8) /= D0 )      ) THEN
                        I_FARF_1N_3D( Ixy_1N, Iz ) = 2
                     ELSE
                        I_FARF_1N_3D( Ixy_1N, Iz ) = 1
                     END IF
#endif
                  ELSE
                     IF ( ( XSset_Table( 7, I_SP, I_BU, I_Tab) /= 0d0 ) .OR.  &
                          ( XSset_Table( 8, I_SP, I_BU, I_Tab) /= 0d0 )      ) THEN
                        I_FARF_1N_3D( Ixy_1N, Iz ) = 2
                     ELSE
                        I_FARF_1N_3D( Ixy_1N, Iz ) = 1
                     END IF

                  END IF
               END DO
            END IF
         END DO
#ifdef js_mpi
         if (comm%usempi) then
            call alloc(tmp_farf_1n_3d,nxy_1n,nz)
            do iz=1,nz
               do ixy=1,nxy
                  ixy_1n=i_4nto1n(ixy)
                  if (iproc/=node2iproc(ixy,iz)) cycle
                  tmp_farf_1n_3d(ixy_1n,iz)=i_farf_1n_3d(ixy_1n,iz)
               enddo
            enddo
            call allreduce(tmp_farf_1n_3d,nxy_1n,nz)
            i_farf_1n_3d=tmp_farf_1n_3d
            deallocate (tmp_farf_1n_3d)
         endif
#endif
         Flag_LoopExit = .FALSE.
         DO Iz = 1, Nz
            IF ( Flag_LoopExit .EQV. .TRUE. ) THEN
               Flag_LoopExit = .FALSE.
               EXIT
            END IF

            DO Ixy = 1, Nxy
               Ixy_1N = I_4Nto1N(Ixy)
               IF ( I_FARF_1N_3D( Ixy_1N, Iz ) >= 2 ) THEN
                  IzFuelBot = Iz
                  Flag_LoopExit = .TRUE.
                  EXIT
               END IF
            END DO
         END DO

         ! cal bot axial size
         ibot=izfuelbot
         itop=(izfuelbot+izfueltop)/2
         hbot=0d0
         do iz=ibot,itop
            hbot=hbot+gridsize_z(iz)
         enddo
         ! cal top axial size
         ibot=(izfuelbot+izfueltop)/2+mod(izfuelbot+izfueltop,2)
         itop=izfueltop
         htop=0d0
         do iz=ibot,itop
            htop=htop+gridsize_z(iz)
         enddo
         ! check core axial mirror
         if (abs(hbot-htop)<1d-10) then
            flag_axial_mirror=.true.
         else
            flag_axial_mirror=.false.
         endif

         Flag_LoopExit = .FALSE.
         DO Iz = Nz, 1, - 1
            IF ( Flag_LoopExit .EQV. .TRUE. ) THEN
               Flag_LoopExit = .FALSE.
               EXIT
            END IF

            DO Ixy = 1, Nxy
               Ixy_1N = I_4Nto1N(Ixy)
               IF ( I_FARF_1N_3D( Ixy_1N, Iz ) >= 2 ) THEN
                  IzFuelTop = Iz
                  Flag_LoopExit = .TRUE.
                  EXIT
               END IF
            END DO
         END DO

         if (flag_auto_geom) then
            GridSize_z(izfuelbot:izfueltop)=coef_the_z*GridSize_z(izfuelbot:izfueltop)
            MeshSize_z(izfuelbot:izfueltop) = GridSize_z(izfuelbot:izfueltop)
            call Get_NodeSize
!            call Get_NodeVolume
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         endif

         IF ( ALLOCATED( PinMTU_3D ) ) THEN
            DO Ixy = 1, Nxy
               Ixy_1N = I_4Nto1N(Ixy)
               I_LP   = I_LP_1N( Ixy_1N )
               DO Iz = IzFuelBot, IzFuelTop
                  PinMTU_3D(I_LP, Iz) = PinMTU_2D(I_LP) / (IzFuelTop - IzFuelBot + 1)
               END DO
            END DO
         END IF
         CALL SetGeometry
      END IF

      RETURN
      END SUBROUTINE Get_FAandRF_Index


      subroutine Get_MTU
      use Inc_Depletion, only: Tot_MTU
      use Inc_3D, only: MTU
      use inc_nuclide, only: hmmass
      use inc_pinpow, only: PinMTU_2D
      use inc_rp, only: N_LP, index_LP
      use inc_fa, only: n_pin
      use inc_pinpow,only: PinMTU_3D
      use inc_xs_file, only: asym_tab, asym2itab
      use mod_charedit, only: print_msg
      use inc_depletion, only: n_bu
      use inc_depletion, only: cycle_bu
      use inc_depletion, only: cycle_day
      use inc_depletion, only: inp_dT
      use inc_depletion, only: inp_dBU
      use inc_depletion, only: inp_ppower
      use inc_depletion, only: inp_core_power
      use inc_depletion, only: inp_fa_power
      use inc_depletion, only: tot_mtu
      use inc_depletion, only: if_buinput
      use inc_flag, only: flag_powhist
      use inc_history, only: hist_ppower
      use inc_th, only: ppower, core_power_100
      use inc_geometry, only: core_n_fa
      implicit none
      real(8) :: Nin1FA
      integer :: Ixy, Iz, itab
      integer :: Nz_Fuel
      integer :: i,i_lp, ia
      integer :: Ixy_1N
      integer :: i_asym, i_aset, isub
      real(8) :: hmmsub
      integer :: nstep, istep

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_MTU] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_MTU] in Mod_GetNode'
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

      call Alloc( MTU , Nxy , Nz )
      if ( Flag_4N1FA ) then
         Nin1FA=4d0
      else
         Nin1FA=1d0
      endif
      Nz_Fuel=IzFuelTop-IzFuelBot+1
      if (sum(hmmass)>1E-10) then
         call print_msg(0,'Use HM mass from XS file...')
         LP_MTU=0.0d0
         do ia=1,N_LP
            i=index_LP(ia)
            do iz=1,nz
               if (.not.allocated(asym_tab)) then
                  lp_mtu(i)=lp_mtu(i)+hmmass(axialcomp(i,iz))*gridsize_z(iz)
               else
                  itab=axialcomp(i,iz)
                  i_asym=asym_tab(itab,1)
                  if (i_asym > 0) then
                     i_aset=asym_tab(itab,2)
                     hmmsub=0d0
                     do isub=1,4
                        hmmsub=hmmsub+hmmass( asym2itab(isub,i_aset) )
                     enddo
                     hmmsub=hmmsub/4d0
                     lp_mtu(i)=lp_mtu(i)+hmmsub*gridsize_z(iz)
                  else
                     lp_mtu(i)=lp_mtu(i)+hmmass(itab)*gridsize_z(iz)
                  endif
               endif
            enddo
         enddo
         call print_msg(0,' No.','   FA','          MTU')
         do ia=1,n_LP
            i=index_LP(ia)
            if(lp_mtu(i)<1d-10) then
                lp_mtu(i)=0d0
            endif
            call print_msg(0,i,trim(Name_LP(i)),lp_mtu(i))
         enddo
      endif

      do ia = 1, N_LP
         I_LP=index_LP(ia)
         PinMTU_2D(I_LP) = LP_MTU(I_LP) / N_Pin
      enddo
      if (allocated(PinMTU_3D)) then
         do Ixy = 1, Nxy
            Ixy_1N = I_4Nto1N(Ixy)
            I_LP   = I_LP_1N( Ixy_1N )
            do Iz = IzFuelBot, IzFuelTop
               PinMTU_3D(I_LP, Iz) = PinMTU_2D(I_LP) / (IzFuelTop - IzFuelBot + 1)
            enddo
         enddo
      endif

      do ixy=1,nxy
         if (i_farf_1n(i_4nto1n(ixy))==1) cycle
         do Iz = IzFuelBot, IzFuelTop
            if(sum(hmmass)>1E-10) then
               mtu(ixy,iz)=hmmass(axialcomp(i_lp_1n(i_4nto1n(ixy)),iz))*gridsize_z(iz)
               if (flag_4n1fa) mtu(ixy,iz)=mtu(ixy,iz)/4d0
            else
               MTU(Ixy, Iz) = LP_MTU( I_LP_1N( I_4Nto1N(Ixy) ) ) / ( Nin1FA*Nz_Fuel )
            endif
         enddo
      enddo

      Tot_MTU=SUM(MTU)
      if (OPT_Core==4) then
         Tot_MTU=SUM(MTU)*4d0
      endif

      call print_msg(0,' ')
      call print_msg(0,'Total Loading: ',Tot_MTU, ' [MT]')
      call print_msg(0,' ')


      ! get dBU, dT, history ...
      nstep=size(cycle_bu)
      if (.not.allocated( inp_dT           )) allocate( inp_dT           (nstep))
      if (.not.allocated( inp_dBU          )) allocate( inp_dBU          (nstep))
      if (.not.allocated( inp_ppower       )) allocate( inp_ppower       (nstep))
      if (.not.allocated( inp_core_power   )) allocate( inp_core_power   (nstep))
      if (.not.allocated( inp_fa_power     )) allocate( inp_fa_power     (nstep))
      if (.not.allocated( inp_fa_power     )) allocate( inp_fa_power     (nstep))
      inp_dT           =0d0
      inp_dBU          =0d0
      inp_ppower       =0d0
      inp_core_power   =0d0
      inp_fa_power     =0d0

      nstep=n_bu
      do istep=1,nstep
         if (flag_powhist) then
            inp_ppower(istep)=hist_ppower(istep)*1d-2
         else
            inp_ppower(istep)=ppower
         endif
         inp_core_power(istep)=core_power_100*inp_ppower(istep)
         inp_fa_power(istep)=inp_core_power(istep)/real(core_n_fa,8)
         if (istep>1) then
            if (if_buinput) then
               inp_dBU(istep)=cycle_bu(istep)-cycle_bu(istep-1)
               inp_dT(istep)=inp_dBU(istep)*tot_mtu*86400d0*1d+9/inp_core_power(istep)
               cycle_day(istep)=cycle_day(istep-1)+inp_dT(istep)/86400d0
            else
               inp_dT(istep)=(cycle_day(istep)-cycle_day(istep-1))*86400d0
               inp_dBU(istep)=inp_dT(istep)*inp_core_power(istep)/(tot_mtu*86400d0*1d+9)
               if (inp_dT(istep)<1d-10) inp_dT(istep)=0d0
               if (inp_dBU(istep)<1d-10) inp_dBU(istep)=0d0
               cycle_bu(istep)=cycle_bu(istep-1)+inp_dBU(istep)
            endif
         endif
      enddo

      return
      end subroutine Get_MTU


      SUBROUTINE Get_NodeSize

      IMPLICIT NONE
      INTEGER :: Ixy, Ix, Iy, Iz

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_NodeSize] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_NodeSize] in Mod_GetNode'
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

      if(.not.allocated(h_x)) CALL Alloc( h_x, Nxy )
      if(.not.allocated(h_y)) CALL Alloc( h_y, Nxy )
      if(.not.allocated(h_z)) CALL Alloc( h_z, Nz  )

      DO Iy = 1, Ny
         DO Ix = Ix_Start_y(Iy), Ix_End_y(Iy)
            Ixy = IxIyToIxy(Ix, Iy)
            h_x(Ixy) = MeshSize_x(Ix)
            h_y(Ixy) = MeshSize_y(Iy)
         END DO
      END DO
      DO Iz = 1, Nz
         h_z(Iz) = MeshSize_z(Iz)
      END DO


      RETURN
      END SUBROUTINE Get_NodeSize


#ifdef siarhei_delete 
      SUBROUTINE Get_NodeVolume

      IMPLICIT NONE
      INTEGER :: Ixy, Iz

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_NodeVolume] in Mod_GetNode'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(.not.allocated(NodeVolume)) CALL Alloc( NodeVolume, Nxy, Nz )

      DO Ixy = 1, Nxy
         DO Iz = 1, Nz
            NodeVolume(Ixy, Iz) = h_x(Ixy) * h_y(Ixy) * h_z(Iz)
         END DO
      END DO



      RETURN
      END SUBROUTINE Get_NodeVolume
#endif 


      SUBROUTINE Get_CoreVolume
      USE Inc_TH, ONLY: Avg_CorePower, Core_Power !!!!!Avg_CorePower0
#ifdef tuan_fr
!    use MTU from input
     use Inc_miXS, ONLY: Dens_Int, ZAID
     use Inc_File,     ONLY: NuNum
     USE Inc_RP, ONLY: AxialComp, I_LP_1N
     USE Inc_Depletion , ONLY: Tot_MTU
     use Inc_FA, ONLY: N_Pin, R_Fuel

#endif

      IMPLICIT NONE
      INTEGER :: Ixy, Iz

#ifdef tuan_fr
      REAL(8):: Sum_Ndens
      INTEGER :: i, Ixy_FA
#endif

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_CoreVolume] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_CoreVolume] in Mod_GetNode'
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

      DO Iz = 1, Nz
         DO Ixy = 1, Nxy
            Tot_Vol = Tot_Vol + NodeVolume(Ixy, Iz)
            IF ( I_FARF_1N_3D( I_4Nto1N(Ixy), Iz ) == 1 ) CYCLE
            Tot_FuelVol = Tot_FuelVol + NodeVolume(Ixy, Iz)
         END DO
      END DO

      Core_Height = SUM( MeshSize_z( IzFuelBot : IzFuelTop ) )
      IF ( OPT_Core == 4 ) THEN
         Tot_Vol     = Tot_Vol     * D4
         Tot_FuelVol = Tot_FuelVol * D4
      END IF
#ifdef tuan_fr
      Tot_FuelVol =0d0

      DO Iz = IzFuelBot,IzFuelTop
         DO Ixy_FA = 1, Nxy_FA
            Ixy=I_FA(Ixy_FA)
            Tot_FuelVol = Tot_FuelVol+N_Pin*(3.141592653589793238*(R_Fuel*100)**2)*MeshSize_z(iz)
         ENDDO
      ENDDO
      
!        need to remove when fraction of fuel in hete/homo are given in input or xs file
!        this volume is for SMLFR only
!         Tot_FuelVol = Tot_FuelVol * 0.447949860           ! tuan_fr
!         Tot_FuelVol = Tot_FuelVol * 0.1          ! tuan_fr
!         Tot_FuelVol = Tot_FuelVol            ! tuan_fr
#endif
!         Tot_FuelVol = 1792456         ! tuan_fr
!      Avg_CorePower = Core_Power / Tot_FuelVol
      !!!!!Avg_CorePower0 = Core_Power / Tot_FuelVol



      RETURN
      END SUBROUTINE Get_CoreVolume


      SUBROUTINE Get_I_Comp

      IMPLICIT NONE
      INTEGER :: Ixy, Iz, itab

      ALLOCATE (I_Comp (1:Nxy,1:Nz))

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_I_Comp] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_I_Comp] in Mod_GetNode'
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

!      CALL Alloc( I_Comp, Nxy, Nz )
      DO Ixy = 1, Nxy
          DO Iz = 1, Nz
             I_Comp(Ixy, Iz) = AxialComp( I_LP_1N( I_4Nto1N(Ixy) ), Iz )
             itab = I_Comp(Ixy, Iz)
             itab = new_asym_itab(itab,Ixy)
             I_Comp(Ixy, Iz) = itab
          END DO
      END DO

      RETURN
      END SUBROUTINE Get_I_Comp


!!!#ifdef siarhei_delete 
      SUBROUTINE Get_I_CR_4N(Buff_Char3)
      USE Inc_CR, ONLY: I_CR_1N, I_CR_4N, N_CR

      IMPLICIT NONE
      CHARACTER(3), DIMENSION(N_CR), INTENT(IN) ::  Buff_Char3
      INTEGER :: Ixy, Ixy_1N, m, I_CR

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_I_CR_4N] in Mod_GetNode'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      DO Ixy = 1, Nxy
         Ixy_1N       = I_4Nto1N(Ixy)
         I_CR         = I_CR_1N(Ixy_1N)
         I_CR_4N(Ixy) = I_CR
      END DO

      DO Ixy = 1, Nxy
         Ixy_1N = I_4Nto1N(Ixy)
         I_CR   = I_CR_1N(Ixy_1N)
         IF ( I_CR /= 0 ) THEN
            IF ( Buff_Char3(I_CR)(2:3) == '12' ) THEN
               DO m = 1, 4
                  IF ( Ixy == I_1Nto4N( Ixy_1N, m ) ) THEN
                     EXIT
                  END IF
               END DO
               SELECT CASE (m)
               CASE (1)
                  IF ( I_Lx(Ixy) /= 0 ) THEN
                     IF ( I_CR_1N( I_4Nto1N( I_Lx(Ixy) ) ) /= 0 ) THEN
                        PRINT *, "Error, Check 12 Finger Type CR"
                        STOP
                     ELSE
                        I_CR_4N( I_Lx(Ixy) ) = I_CR
                     END IF
                  END IF
                  IF ( I_Ly(Ixy) /= 0 ) THEN
                     IF ( I_CR_1N( I_4Nto1N( I_Ly(Ixy) ) ) /= 0 ) THEN
                        PRINT *, "Error, Check 12 Finger Type CR"
                        STOP
                     ELSE
                        I_CR_4N( I_Ly(Ixy) ) = I_CR
                     END IF
                  END IF

               CASE (2)
                  IF ( I_Rx(Ixy) /= 0 ) THEN
                     IF ( I_CR_1N( I_4Nto1N( I_Rx(Ixy) ) ) /= 0 ) THEN
                        PRINT *, "Error, Check 12 Finger Type CR"
                        STOP
                     ELSE
                        I_CR_4N( I_Rx(Ixy) ) = I_CR
                     END IF
                  END IF
                  IF ( I_Ly(Ixy) /= 0 ) THEN
                     IF ( I_CR_1N( I_4Nto1N( I_Ly(Ixy) ) ) /= 0 ) THEN
                        PRINT *, "Error, Check 12 Finger Type CR"
                        STOP
                     ELSE
                        I_CR_4N( I_Ly(Ixy) ) = I_CR
                     END IF
                  END IF

               CASE (3)
                  IF ( I_Lx(Ixy) /= 0 ) THEN
                     IF ( I_CR_1N( I_4Nto1N( I_Lx(Ixy) ) ) /= 0 ) THEN
                        PRINT *, "Error, Check 12 Finger Type CR"
                        STOP
                     ELSE
                        I_CR_4N( I_Lx(Ixy) ) = I_CR
                     END IF
                  END IF
                  IF ( I_Ry(Ixy) /= 0 ) THEN
                     IF ( I_CR_1N( I_4Nto1N( I_Ry(Ixy) ) ) /= 0 ) THEN
                        PRINT *, "Error, Check 12 Finger Type CR"
                        STOP
                     ELSE
                        I_CR_4N( I_Ry(Ixy) ) = I_CR
                     END IF
                  END IF

               CASE (4)
                  IF ( I_Rx(Ixy) /= 0 ) THEN
                     IF ( I_CR_1N( I_4Nto1N( I_Rx(Ixy) ) ) /= 0 ) THEN
                        PRINT *, "Error, Check 12 Finger Type CR"
                        STOP
                     ELSE
                        I_CR_4N( I_Rx(Ixy) ) = I_CR
                     END IF
                  END IF
                  IF ( I_Ry(Ixy) /= 0 ) THEN
                     IF ( I_CR_1N( I_4Nto1N( I_Ry(Ixy) ) ) /= 0 ) THEN
                        PRINT *, "Error, Check 12 Finger Type CR"
                        STOP
                     ELSE
                        I_CR_4N( I_Ry(Ixy) ) = I_CR
                     END IF
                  END IF
               END SELECT
            END IF
         END IF
      END DO

      RETURN
      END SUBROUTINE Get_I_CR_4N
!!!#endif 

#ifdef siarhei_delete 
      subroutine get_asym2itab
      use Inc_XS_File, only: asym_tab, asym2itab, N_XS_Table
      implicit none
      integer :: n_iset, i
      integer :: i_asym, i_aset

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [get_asym2itab] in Mod_GetNode'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      n_iset=maxval(asym_tab(:,2))
      if (.not. allocated(asym2itab)) then
         call alloc(asym2itab,4,n_iset)
         asym2itab=0
      endif

      do i=1,N_XS_Table
         i_asym=asym_tab(i,1)
         i_aset=asym_tab(i,2)
         if (i_asym > 0) then
            asym2itab(i_asym,i_aset)=i
         endif
      enddo

      return
      end subroutine get_asym2itab
#endif 


      function new_asym_itab(inp1,inp2)
      use Inc_XS_File, only: asym_tab, asym2itab, asym_rot
      use Inc_Geometry, only: I_1Nto4N, I_4Nto1N
      implicit none
      integer :: new_asym_itab
      integer :: inp1, inp2, inp3, inp4, inp5
      integer :: i1, i2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [new_asym_itab] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [new_asym_itab] in Mod_GetNode'
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

      inp3 = I_4Nto1N(inp2)  ! Ixy_1N = I_4Nto1N(Ixy)
      inp4 = I_LP_1N(inp3)   ! I_LP_1N => LP
      inp5 = asym_rot(inp4)  ! asym_rot => 0/1/2/3

      if (.not. allocated(asym_tab)) then
         new_asym_itab = inp1
         return
      endif
      if (.not. allocated(asym2itab)) then
         new_asym_itab = inp1
         return
      endif

      if (asym_tab(inp1,1) > 0) then
         i2=asym_tab(inp1,2)
         if (inp5 == 0) then
            if     ( I_1Nto4N(inp3, 1) == inp2 ) then; i1=1;
            elseif ( I_1Nto4N(inp3, 2) == inp2 ) then; i1=2;
            elseif ( I_1Nto4N(inp3, 3) == inp2 ) then; i1=3;
            elseif ( I_1Nto4N(inp3, 4) == inp2 ) then; i1=4;
            endif
         elseif (inp5==1) then
            if     ( I_1Nto4N(inp3, 1) == inp2 ) then; i1=2;
            elseif ( I_1Nto4N(inp3, 2) == inp2 ) then; i1=4;
            elseif ( I_1Nto4N(inp3, 3) == inp2 ) then; i1=1;
            elseif ( I_1Nto4N(inp3, 4) == inp2 ) then; i1=3;
            endif
         elseif (inp5==2) then
            if     ( I_1Nto4N(inp3, 1) == inp2 ) then; i1=4;
            elseif ( I_1Nto4N(inp3, 2) == inp2 ) then; i1=3;
            elseif ( I_1Nto4N(inp3, 3) == inp2 ) then; i1=2;
            elseif ( I_1Nto4N(inp3, 4) == inp2 ) then; i1=1;
            endif
         elseif (inp5==3) then
            if     ( I_1Nto4N(inp3, 1) == inp2 ) then; i1=3;
            elseif ( I_1Nto4N(inp3, 2) == inp2 ) then; i1=1;
            elseif ( I_1Nto4N(inp3, 3) == inp2 ) then; i1=4;
            elseif ( I_1Nto4N(inp3, 4) == inp2 ) then; i1=2;
            endif
         endif
         new_asym_itab = asym2itab(i1,i2)
      else
         new_asym_itab = inp1
      endif

      return
      end function new_asym_itab


      SUBROUTINE SetGeometry
      USE Inc_Geometry
      USE Inc_RP
      USE Inc_Option
      USE Mod_Alloc
      USE Inc_PinVar
      USE Inc_DF, ONLY: cdf
     ! use mod_manage, only: alloc_pin
      use Inc_PinPOW, only: flag_pinpow
#ifdef js_mpi
      use inc_parallel, only: comm
#endif
      IMPLICIT NONE
      INTEGER:: i,j,k,lptr
      INTEGER:: l,la,lfa,jm1,jp1,km1
      INTEGER:: i0,j0,iw,js,jn,ie
      INTEGER:: iam
      INTEGER:: nfnodes,nanodes,nperassy
      INTEGER:: ndomxy
      INTEGER:: nx2,nx2m1
      REAL(8):: hr
      INTEGER:: nnodesfa,inode,lp
      INTEGER,ALLOCATABLE :: laptr(:),laptrt(:),lfaptrt(:)
      INTEGER :: Ix_1N, Iy_1N

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [SetGeometry] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [SetGeometry] in Mod_GetNode'
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

      ! radial node numbering for the core region
      if (.not.allocated(lfaptr)) then
         call alloc0(lfaptr,0,Nxy_1N)
         call alloc0(nodelfa,1,Nx_1N, 1,Ny_1N)
      endif
      call alloc0(laptr,0,Nxy_1N)
      call alloc0(laptrt,0,Nxy_1N)
      call alloc0(lfaptrt,0,Nxy_1N)

      nodelfa = 0
      la=0
      lfa=0
      nfnodes=0
      nanodes=0
      laptr(0)=0
      lfaptr(0)=0
      npfa=1

      DO Iy_1N = 1, Ny_1N
         DO Ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
            la = la + 1
            nperassy=Mesh_x(Ix_1N)*Mesh_y(Iy_1N)
            npfa=MAX(npfa,nperassy)
            nanodes=nanodes+Mesh_x(Ix_1N)*Mesh_y(Iy_1N)
            laptr(la)=nanodes
            IF(Ix_1N.GE.Ix_StartFA_y_1N(Iy_1N) .AND. Ix_1N.LE. Ix_EndFA_y_1N(Iy_1N)) THEN
               lfa=lfa+1
               nodelfa(Ix_1N,Iy_1N)=lfa
               nfnodes=nfnodes+nperassy
               lfaptr(lfa)=nfnodes
            END IF
         END DO
      END DO

      DO la = 1, Nxy_1N
         laptrt(la) = laptr(la-1)
      END DO
      DO lfa = 1, Nxy_FA_1N
         lfaptrt(lfa) = lfaptr(lfa-1)
      END DO

      if (.not.allocated(I_Lx_P2R)) then
         call alloc(I_Lx_P2R,Nxy)
         call alloc(I_Rx_P2R,Nxy)
         call alloc(I_Ly_P2R,Nxy)
         call alloc(I_Ry_P2R,Nxy)
         call alloc0(znode,0,Nz+1)
         call alloc0(zcent,0,Nz+1)
         call alloc0(zb,0,Nz+1)
      endif

      nx2=Nx+Ny
      nx2m1=nx2-1

      if (.not.allocated(ltola)) then
         call alloc0(ltola,0,Nxy)
         call alloc (latol,Nxy)
         call alloc0(ltox,-nx2m1,Nxy+nx2+1)
         call alloc0(ltoy,-nx2m1,Nxy+nx2+1)
         call alloc0(nodel,-3,Nx+4,-1,Ny+2)
         call alloc (nodef,Nxy)
         call alloc0(ltolf,0,Nxy)
         call alloc0(ltolfa,0,Nxy)
         call alloc (lfatol,Nxy)
      endif

      l=0
      Nxy_FA_P2R=0
      nodel=0
      DO j=1,Ny
         nrnx(j)=Ix_End_y(j)-Ix_Start_y(j)+1
         Iy_1N=Iy_4Nto1N(j)
         DO i=Ix_Start_y(j),Ix_End_y(j)
            l=l+1
            Ix_1N=Ix_4Nto1N(i)
            nodel(i,j)=l
            la = IxIy_1NToIxy_1N(Ix_1N,Iy_1N)
            ltola(l)=la
            laptrt(la)=laptrt(la)+1
            latol(laptrt(la))=l
            ltox(l)=i
            ltoy(l)=j
#ifdef tuan_fr_crm
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) == 2) THEN
#else
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) > 1) THEN
#endif
               Nxy_FA_P2R=Nxy_FA_P2R+1
               nodef(Nxy_FA_P2R)=l
               ltolf(l)=Nxy_FA_P2R
               lfa=nodelfa(Ix_1N,Iy_1N)
               ltolfa(l)=lfa
               lfaptrt(lfa)=lfaptrt(lfa)+1
               lfatol(lfaptrt(lfa))=l
            ELSE
               ltolfa(l)=0
            END IF
         END DO
      END DO

      ! determine neighboring nodes
      DO j=1,Ny
          jm1=j-1
          jp1=j+1
          DO i=Ix_Start_y(j),Ix_End_y(j)
             l=nodel(i,j)
             I_Lx_P2R(l)=nodel(i-1,j)
             I_Rx_P2R(l)=nodel(i+1,j)
             I_Ly_P2R(l)=nodel(i,jm1)
             I_Ry_P2R(l)=nodel(i,jp1)
             IF(I_Lx_P2R(l).EQ.0) I_Lx_P2R(l)=l
             IF(I_Rx_P2R(l).EQ.0) I_Rx_P2R(l)=l
             IF(I_Ly_P2R(l).EQ.0) I_Ly_P2R(l)=l
             IF(I_Ry_P2R(l).EQ.0) I_Ry_P2R(l)=l
          END DO
      END DO

      ! calculate node volumes and axial coordinate of top of each plane
      hr=MeshSize_x(1)
      znode(0)=0
      km1=0
      DO k=1,Nz
         znode(k)=znode(km1)+ MeshSize_z(k)
         km1=k
      END DO

      ! for axial plots
      zb(0)=0
      DO k=1,Nz
         zb(k)=znode(k)
         zcent(k)=0.5*(zb(k-1)+zb(k))
      END DO
      hactive = D0
      DO k=IzFuelBot,IzFuelTop
         hactive=hactive+ MeshSize_z(k)
         DO lfa = 1, Nxy_FA_1N
            DO lptr=lfaptr(lfa-1)+1,lfaptr(lfa)
               l=lfatol(lptr)
            END DO
         END DO
      END DO

      if (izfuelbot<=1) then
         hbot_active = 0d0
      else
         hbot_active = sum(MeshSize_z(1:IzFuelBot-1))
      endif
      htop_active = sum(MeshSize_z(1:IzFuelTop))
      ! fill a nonzero for the meshes outside the problem geometry
      MeshSize_x(0)=1.
      MeshSize_x(Nx+1)=1.
      MeshSize_y(0)=1.
      MeshSize_y(Ny+1)=1.
      MeshSize_z(0)=1
      MeshSize_z(Nz+1)=1

      ! establish node number correspondence
      DO j=1,Ny
         DO i=1,Ix_Start_y(j)-1
            nodel(i,j)=0
         END DO
         DO i=Ix_End_y(j)+1,Nx
            nodel(i,j)=0
         END DO
      END DO
      ! west boundary ordering
      l=-Nx+1
      iw=0
      DO j=Ny,1,-1
         i0=Ix_Start_y(j)-1
         l=l-1
         nodel(iw,j)=l
         nodel(i0,j)=l
         ltox(l)=i0
         ltoy(l)=j
      END DO
      ! north boundary ordering
      l=1
      jn=0
      DO i=Nx,1,-1
         j0=Iy_Start_x(i)-1
         l=l-1
         nodel(i,jn)=l
         nodel(i,j0)=l
         ltox(l)=i
         ltoy(l)=j0
      END DO
      ! east boundary ordering
      l=Nxy+Nx
      ie=Nx+1
      DO j=1,Ny
         i0=Ix_End_y(j)+1
         l=l+1
         nodel(ie,j)=l
         nodel(i0,j)=l
         ltox(l)=i0
         ltoy(l)=j
      END DO
      ! south boundary ordering
      l=Nxy
      js=Ny+1
      DO i=1,Nx
         j0=Iy_End_x(i)+1
         l=l+1
         nodel(i,js)=l
         nodel(i,j0)=l
         ltox(l)=i
         ltoy(l)=j0
      END DO
      ! south-esat corner added for corner flux calculation
      l=Nxy+Nx+Ny+1
      nodel(Nx+1,Ny+1)=l
      ltox(l)=Nx+1
      ltoy(l)=Ny+1
      ndomx=1
      ndomy=1
      ndomz=1
      ndomxy=ndomx*ndomy
      ndom=ndomz*ndomxy
      iam=0
      idom=iam+1
      idomxy=MOD(iam,ndomxy)
      idomx=MOD(idomxy,ndomx)+1
      idomy=idomxy/ndomx+1
      idomz=iam/ndomxy+1
      idomz=1 !temp for DEC alpha
      idomxy=idomxy+1

      ! set ncorner and iquad for pin power reconstuction
      ncorn=1
      DO i=Ix_Start_y(1),Ix_End_y(1)
          ncorn=ncorn+1
      END DO
      DO j=1,ny
          ncorn=ncorn+1
          if ( (Ix_End_y(j)-Ix_Start_y(j))<(Ix_End_y(j+1)-Ix_Start_y(j+1)) ) then
             DO i=Ix_Start_y(j+1),Ix_End_y(j+1)
                 ncorn=ncorn+1
             END DO
          else
             DO i=Ix_Start_y(j),Ix_End_y(j)
                 ncorn=ncorn+1
             END DO
          endif
      END DO

      if (.not.allocated(iquad)) then
#ifdef jr_vver
         allocate(iquad(0:nxy))
#else
         call alloc(iquad,nxy)
#endif
      endif

      DO lfa = 1, Nxy_FA_1N
         nnodesfa=lfaptr(lfa)-lfaptr(lfa-1)
         inode=0
         DO lp=lfaptr(lfa-1)+1,lfaptr(lfa)
            inode=inode+1
            l=lfatol(lp)
            IF(nnodesfa.EQ.1) THEN
               IF(ltox(l).EQ.1) THEN
                  IF(ltoy(l).EQ.1) THEN
                     iquad(l)=3
                  ELSE
                     iquad(l)=2
                  END IF
               ELSE
                  IF(ltoy(l).EQ.1) THEN
                     iquad(l)=4
                  ELSE
                     iquad(l)=1
                  END IF
               END IF
            ELSEIF(nnodesfa.EQ.2) THEN
               IF(ltox(l).EQ.1) THEN
                  IF(inode.EQ.1) THEN
                     iquad(l)=2
                  ELSE
                     iquad(l)=3
                  END IF
               ELSEIF(ltox(l).EQ.Nx) THEN
                  IF(inode.EQ.1) THEN
                     iquad(l)=1
                  ELSE
                     iquad(l) = 3
                  END IF
               ELSEIF(ltoy(l).EQ.1) THEN
                  IF(inode.EQ.1) THEN
                     iquad(l)=4
                  ELSE
                     iquad(l)=3
                  END IF
               ELSE
                  IF(inode.EQ.1) THEN
                     iquad(l)=1
                  ELSE
                     iquad(l)=2
                  END IF
               END IF
            ELSE
               IF(inode.LT.3) THEN
                  iquad(l)=inode
               ELSE
                  iquad(l)=7-inode
               END IF
            END IF
         END DO
      END DO

      if (flag_pinpow) then
         call alloc(cdf, 4, n_group, nxy, nz )
         cdf=1d0
!         call alloc_pin
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         i = 0
         do Iy_1N=1,Ny_1N
            do Ix_1N=1,Ny_1N
#ifdef tuan_fr_crm
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) == 2) THEN
#else
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) > 1) THEN
#endif
                  i=i+1
                  pploc(i)=nodelfa(Ix_1N,Iy_1N)
               endif
            enddo
         enddo

         if (npfa==1) then
            if (nxy==1) then
               pinpitch = h_x(1)*D2 / npin
               pin_area = (h_x(1)*D2/npin)*(h_x(1)*D2/npin)
            else
               pinpitch = h_x(1) / npin
               pin_area = (h_x(1)/npin)*(h_x(1)/npin)
            endif
         else
            pinpitch = h_x(1)*D2 / npin
            pin_area = (h_x(1)*D2/npin)*(h_x(1)*D2/npin)
         endif
      endif



      RETURN
      END SUBROUTINE SetGeometry

#ifdef js_mpi
#ifdef siarhei_delete 
      subroutine set_mpi_nodel
      use inc_geometry, only: nxy
      use inc_geometry, only: ltox, ltoy
      use inc_parallel, only: iproc
      use inc_parallel, only: ixy2iproc, ixyip2ixy
      use inc_parallel, only: nodelmpi, nnodelmpi
      implicit none
      integer :: l, i, j, is, ie, js, je
      integer :: l0, l3, flagn

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [set_mpi_nodel] in Mod_GetNode'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      flagn=0
      do l=1,nxy
         i=ltox(l)
         j=ltoy(l)
         if (iproc/=ixy2iproc(l)) cycle
         if (flagn==0) then
            is=i
            ie=i
            js=j
            je=j
         endif
         flagn=flagn+1
         l0=ixyip2ixy(l,iproc)
         if (i<is) is=i
         if (i>ie) ie=i
         if (j<js) js=j
         if (j>je) je=j
      enddo

      allocate(nodelmpi(-1:nx+1,-1:ny+1))
      nodelmpi=0

      l3=0
      do j=js-1,je+1
         do i=is-1,ie+1
            l3=l3+1
            nodelmpi(i,j)=l3 ! global i,j -> local l
         enddo
      enddo
      nnodelmpi=l3

      return
      end subroutine set_mpi_nodel
#endif 
#endif

#ifdef jr_vver
      SUBROUTINE Geometry_hex
      USE Inc_TPEN
!      USE Inc_Control,  only: ninitoutt
!      USE Mod_XSFB,     only: mgto2g_hex
      !USE Inc_3D,       only: flux
      IMPLICIT NONE

!      INTEGER(4),ALLOCATABLE,DIMENSION(:,:), SAVE :: iradconf
!      INTEGER(4),ALLOCATABLE,DIMENSION(:)   :: lfaptrt
!      INTEGER(4),ALLOCATABLE,DIMENSION(:)   :: laptrt


!      INTEGER(4),ALLOCATABLE,DIMENSION(:)   :: laptr
      
!      INTEGER(4),ALLOCATABLE,DIMENSION(:), SAVE   :: ncrbasy
      INTEGER(4) :: i !, ndataf
      INTEGER(4) :: iat, k
      INTEGER(4) :: km1, j, ipr, la, lfa
      INTEGER(4) :: l, icomp
      REAL(8)    :: xsecnf, val
!      REAL(8), ALLOCATABLE, SAVE :: hasyx(:),hasyy(:) !temp
!      CHARACTER(len=4) :: astr
      LOGICAL(1)       :: skip30
      !LOGICAL(1), SAVE :: gfirst = TRUE
!      INTEGER(4) :: mf

      INTEGER(4) :: nxskip
!      LOGICAL(1) :: Inp_Ref = .false.
!      REAL(8)    :: hr, hrsq
      INTEGER(4) :: iam = 1

#ifdef tuan_tr_test
      integer(4) :: lptr
#else
      INTEGER(4),ALLOCATABLE,DIMENSION(:)   :: laptr
#endif

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Geometry_hex] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Geometry_hex] in Mod_GetNode'
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


      nxp1=0
      nyp1=0




! === NodeVolume_hex === !
      call Get_NodeVolume_hex
!
#ifdef tuan_tr_test
      npr=Nz
#else
      
#endif
!
! === saninput === !
      call hex_map_make
!
! initialze assembly weight array
      val=1
      wtass=val
      npr=Nz

      call scanlay
      call readlay
      call orderhex
!
! generate planar regions
      npr=1
      izpr(0)=0
      izpr(1)=1
      DO k=2,nz
         km1=k-1
         DO iat=1,nat
            IF(iasypr(km1,iat).NE.iasypr(k,iat)) THEN
               npr=npr+1
               EXIT
            ENDIF
         ENDDO
         izpr(k)=npr
      ENDDO
!
      DO k=1,nz
         IF(izpr(k-1).NE.izpr(k)) THEN
            ipr=izpr(k)
            DO la=1,nassy
               iat=iasytype(la)
               iprcomp(la,ipr)=iasypr(k,iat)
            ENDDO
         ENDIF
      ENDDO

      DO k=1,nz
         DO la=1,nassy
            iat=iasytype(la)
            xsid(la,k)=iasypr(k,iat)
         ENDDO
      ENDDO
!
! determine fuel assemblies
      if(.not.allocated(laptr))      allocate(laptr(0:Nxy_1N))
      if(.not.allocated(lfaptr_hex)) allocate(lfaptr_hex(0:Nxy_FA_1N))
      if(.not.allocated(latol_hex))  allocate(latol_hex(nxy))
      if(.not.allocated(ltola_hex))  allocate(ltola_hex(0:nxy))
      if(.not.allocated(ltolfa_hex)) allocate(ltolfa_hex(0:nxy))
      if(.not.allocated(lfatol_hex)) allocate(lfatol_hex(nxy))
      if(.not.allocated(nodef_hex))  allocate(nodef_hex(nxy))
      if(.not.allocated(nxfas))      allocate(nxfas(Ny_1N))
      if(.not.allocated(nxfae))      allocate(nxfae(Ny_1N))
      if(.not.allocated(nxas))       allocate(nxas(Ny_1N))
      if(.not.allocated(nxae))       allocate(nxae(Ny_1N))

! determine fuel assemblies
      lfa=0
      l=0
      laptr(0)=0
      lfaptr_hex(0)=0
      wtasssum=0d0
      DO la=1,nassy
         l=l+1
         latol_hex(la)=l
         ltola_hex(l)=la
         laptr(la)=la
         IF(iffuel(iasytype(la))) THEN
            lfa=lfa+1
            ltolfa_hex(l)=lfa
            lfaptr_hex(lfa)=lfa
            lfatol_hex(lfa)=l
            nodef_hex(lfa)=l
            wtasssum=wtasssum+wtass(la)
         ENDIF
      ENDDO
! end of change
      nfuelfa=lfa
      nfuel=nfuelfa
      nchan_hex=nfuelfa
      IF( (.NOT.fdbk).OR.extth ) THEN
         nchan_hex=nxy
      ELSE
         nchan_hex=nfuelfa
      ENDIF
!
      nxskip = 2
      DO j=1,ny
         DO i=nxs(j),nxe(j),nxskip
            lfa=ltolfa_hex(nodel_hex(i,j))
            IF((ltolfa_hex(nodel_hex(i-nxskip,j)).EQ.0) .AND. &
                 lfa.NE.0) nxfas(j)=i
            IF((ltolfa_hex(nodel_hex(i+nxskip,j)).EQ.0) .AND. &
                 lfa.NE.0) nxfae(j)=i
            nxas(j)=nxs(j)
            nxae(j)=nxe(j)
         ENDDO
      ENDDO
!
!
      xsecnf = 0d0
!
! determine starting and ending fuel plane numbers
      skip30 = FALSE
      la=ltola_hex(lfatol_hex(1)) !first fuel assy
      kfs=1
      kfe=nz
      DO k=1,nz
         icomp=iprcomp(la,izpr(k))
         IF(icomp.EQ.0) THEN
            skip30 = TRUE
            EXIT
         END IF
         IF(xsecnf.NE.0) THEN
            kfs=k
            EXIT
         ENDIF
      ENDDO
      IF(.NOT.skip30) THEN
         DO k=nz,1,-1
            icomp=iprcomp(la,izpr(k))
            IF(xsecnf.NE.0) THEN
               kfe=k
               EXIT
            ENDIF
         ENDDO
      END IF
!  establish bank number to assembly number correspondence
!      IF(nrcrb.NE.0.AND.hexbyring) CALL assignbank
!
! initial axial flux shape
!
! -- flat flux initial guess using 4:1 fast to thermal flux ratio
!
! for axial plots
      !IF(flxlevel.EQ.0.) flxlevel=1
      !DO k=1,nz
      !   axf(1,k)=0.8*flxlevel
      !   axf(2,k)=0.2*flxlevel
      !ENDDO
!
! calculate node volumes and axial coordinate of top of each plane
!     check subroutine Get_NodeVolume_hex
!
!
!
! fill a nonzero for the meshes outside the problem geometry
      nxp1=nx_hex+1
      nyp1=ny_hex+1
      hx(0)=one
      hx(nxp1)=one
      hy(0)=one
      hy(nyp1)=one
      hz(0)=one
      hz(nzp1)=one
!
!
! determine domain structure, dummies for serial runs
      ndomx_hex=1
      ndomy_hex=1
      ndomz_hex=1
      ndomxy_hex=ndomx_hex*ndomy_hex
      ndom_hex=ndomz_hex*ndomxy_hex
      !IF(ndom.NE.ncpus) THEN
      !   WRITE(*,*)' Domain Configuration Inconsistent'
      !   STOP ' Domain Configuration Inconsistent'
      !ENDIF
      idom_hex=iam+1
      idomxy_hex=MOD(iam,ndomxy_hex)
      idomx_hex=MOD(idomxy_hex,ndomx_hex)+1
      idomy_hex=idomxy_hex/ndomx_hex+1
      idomz_hex=iam/ndomxy_hex+1
      idomz_hex=1 !temp for DEC alpha
      idomxy_hex=idomxy_hex+1
!
      nxyb=nfuelfa
      ltolb=ltolfa
!
      ! nmzbi is zero.
      !IF(nmzbi>0)THEN
      !   IF(nmzbi<nz)THEN
      !      bndlpow=.TRUE.
      !      ktokb=0
      !      kb=0
      !      kze=0
      !      DO k=1,nmzbi
      !         kzs=kze+1
      !         kze=kze+nmbz(k)
      !         IF(kzs<=kfe.AND.kze>=kfs)THEN
      !            kb=kb+1
      !            DO kz=kzs,kze
      !               ktokb(kz)=kb
      !            ENDDO
      !         ENDIF
      !      ENDDO
      !      nzb=kb
      !      kbe=kb
      !      kbs=1
      !   ENDIF
      !   DEALLOCATE(nmbz)
      !ELSE
!
      DO k=1,nz
         ktokb(k)=k
      ENDDO
      nzb=nz
      kbs=kfs
      kbe=kfe
!
      !ENDIF
      call rearrange_index_hex

#ifdef tuan_tr_test
      !! volume
      Tot_FuelVol_HEX=0
      do k = IzFuelBot, IzFuelTop
         do lfa = 1,Nxy_FA_1N
            do lptr = lfaptr_hex(lfa-1)+1,lfaptr_hex(lfa)
               l = lfatol_hex(lptr)
               Tot_FuelVol_HEX = Tot_FuelVol_HEX + nodevolume(imap(l),k)
            enddo
         enddo
      enddo 
#endif


      ninitoutt = 3

!      nxp1=0
!      nyp1=0
!
!
!
!
!! === NodeVolume_hex === !
!      call Get_NodeVolume_hex
!!
!#ifdef tuan_tr_test
!      npr=Nz
!#endif
!! === saninput === !
!      call hex_map_make
!!
!! initialze assembly weight array
!      val=1
!      wtass=val
!      npr=Nz
!
!      call scanlay
!      call readlay
!      call orderhex
!!
!! generate planar regions
!      npr=1
!      izpr(0)=0
!      izpr(1)=1
!      DO k=2,nz
!         km1=k-1
!         DO iat=1,nat
!            IF(iasypr(km1,iat).NE.iasypr(k,iat)) THEN
!               npr=npr+1
!               EXIT
!            ENDIF
!         ENDDO
!         izpr(k)=npr
!      ENDDO
!!
!      DO k=1,nz
!         IF(izpr(k-1).NE.izpr(k)) THEN
!            ipr=izpr(k)
!            DO la=1,nassy
!               iat=iasytype(la)
!               iprcomp(la,ipr)=iasypr(k,iat)
!            ENDDO
!         ENDIF
!      ENDDO
!
!      DO k=1,nz
!         DO la=1,nassy
!            iat=iasytype(la)
!            xsid(la,k)=iasypr(k,iat)
!         ENDDO
!      ENDDO
!!
!! determine fuel assemblies
!      ALLOCATE(laptr(0:Nxy_1N))
!      ALLOCATE(lfaptr_hex(0:Nxy_FA_1N))
!      ALLOCATE(latol_hex(nxy))
!      ALLOCATE(ltola_hex(0:nxy))
!      ALLOCATE(ltolfa_hex(0:nxy))
!      ALLOCATE(lfatol_hex(nxy))
!      ALLOCATE(nodef_hex(nxy))
!      ALLOCATE(nxfas(Ny_1N))
!      ALLOCATE(nxfae(Ny_1N))
!      ALLOCATE(nxas(Ny_1N))
!      ALLOCATE(nxae(Ny_1N))
!
!! determine fuel assemblies
!      lfa=0
!      l=0
!      laptr(0)=0
!      lfaptr_hex(0)=0
!      wtasssum=0d0
!      DO la=1,nassy
!         l=l+1
!         latol_hex(la)=l
!         ltola_hex(l)=la
!         laptr(la)=la
!         IF(iffuel(iasytype(la))) THEN
!            lfa=lfa+1
!            ltolfa_hex(l)=lfa
!!            lfaptr_hex(lfa)=lfa
!            lfatol_hex(lfa)=l
!            nodef_hex(lfa)=l
!            wtasssum=wtasssum+wtass(la)
!         ENDIF
!      ENDDO
!! end of change
!      nfuelfa=lfa
!      nfuel=nfuelfa
!      nchan_hex=nfuelfa
!      IF( (.NOT.fdbk).OR.extth ) THEN
!         nchan_hex=nxy
!      ELSE
!         nchan_hex=nfuelfa
!      ENDIF
!!
!      nxskip = 2
!      DO j=1,ny
!         DO i=nxs(j),nxe(j),nxskip
!            lfa=ltolfa_hex(nodel_hex(i,j))
!            IF((ltolfa_hex(nodel_hex(i-nxskip,j)).EQ.0) .AND. &
!                 lfa.NE.0) nxfas(j)=i
!            IF((ltolfa_hex(nodel_hex(i+nxskip,j)).EQ.0) .AND. &
!                 lfa.NE.0) nxfae(j)=i
!            nxas(j)=nxs(j)
!            nxae(j)=nxe(j)
!         ENDDO
!      ENDDO
!!
!!
!      xsecnf = 0d0
!!
!! determine starting and ending fuel plane numbers
!      skip30 = FALSE
!      la=ltola_hex(lfatol_hex(1)) !first fuel assy
!      kfs=1
!      kfe=nz
!      DO k=1,nz
!         icomp=iprcomp(la,izpr(k))
!         IF(icomp.EQ.0) THEN
!            skip30 = TRUE
!            EXIT
!         END IF
!         IF(xsecnf.NE.0) THEN
!            kfs=k
!            EXIT
!         ENDIF
!      ENDDO
!      IF(.NOT.skip30) THEN
!         DO k=nz,1,-1
!            icomp=iprcomp(la,izpr(k))
!            IF(xsecnf.NE.0) THEN
!               kfe=k
!               EXIT
!            ENDIF
!         ENDDO
!      END IF
!!  establish bank number to assembly number correspondence
!!      IF(nrcrb.NE.0.AND.hexbyring) CALL assignbank
!!
!! initial axial flux shape
!!
!! -- flat flux initial guess using 4:1 fast to thermal flux ratio
!!
!! for axial plots
!      !IF(flxlevel.EQ.0.) flxlevel=1
!      !DO k=1,nz
!      !   axf(1,k)=0.8*flxlevel
!      !   axf(2,k)=0.2*flxlevel
!      !ENDDO
!!
!! calculate node volumes and axial coordinate of top of each plane
!!     check subroutine Get_NodeVolume_hex
!!
!!
!!
!! fill a nonzero for the meshes outside the problem geometry
!      nxp1=nx_hex+1
!      nyp1=ny_hex+1
!      hx(0)=one
!      hx(nxp1)=one
!      hy(0)=one
!      hy(nyp1)=one
!      hz(0)=one
!      hz(nzp1)=one
!!
!!
!! determine domain structure, dummies for serial runs
!      ndomx_hex=1
!      ndomy_hex=1
!      ndomz_hex=1
!      ndomxy_hex=ndomx_hex*ndomy_hex
!      ndom_hex=ndomz_hex*ndomxy_hex
!      !IF(ndom.NE.ncpus) THEN
!      !   WRITE(*,*)' Domain Configuration Inconsistent'
!      !   STOP ' Domain Configuration Inconsistent'
!      !ENDIF
!      idom_hex=iam+1
!      idomxy_hex=MOD(iam,ndomxy_hex)
!      idomx_hex=MOD(idomxy_hex,ndomx_hex)+1
!      idomy_hex=idomxy_hex/ndomx_hex+1
!      idomz_hex=iam/ndomxy_hex+1
!      idomz_hex=1 !temp for DEC alpha
!      idomxy_hex=idomxy_hex+1
!!
!      nxyb=nfuelfa
!      ltolb=ltolfa
!!
!      ! nmzbi is zero.
!      !IF(nmzbi>0)THEN
!      !   IF(nmzbi<nz)THEN
!      !      bndlpow=.TRUE.
!      !      ktokb=0
!      !      kb=0
!      !      kze=0
!      !      DO k=1,nmzbi
!      !         kzs=kze+1
!      !         kze=kze+nmbz(k)
!      !         IF(kzs<=kfe.AND.kze>=kfs)THEN
!      !            kb=kb+1
!      !            DO kz=kzs,kze
!      !               ktokb(kz)=kb
!      !            ENDDO
!      !         ENDIF
!      !      ENDDO
!      !      nzb=kb
!      !      kbe=kb
!      !      kbs=1
!      !   ENDIF
!      !   DEALLOCATE(nmbz)
!      !ELSE
!!
!      DO k=1,nz
!         ktokb(k)=k
!      ENDDO
!      nzb=nz
!      kbs=kfs
!      kbe=kfe
!!
!      !ENDIF
!      call rearrange_index_hex
!
!      ninitoutt = 3

      RETURN
      END SUBROUTINE Geometry_hex

      SUBROUTINE rearrange_index_hex
      USE Inc_TPEN
#ifdef tuan_tr_test
      USE Inc_Flag, only: Flag_THFB 
#endif
      
      IMPLICIT NONE
      INTEGER(4) :: i, j, n, ind
!      INTEGER(4), ALLOCATABLE, DIMENSION(:) :: imap_mid
!      INTEGER(4) :: ixy, iz
#ifdef tuan_tr_test
      INTEGER(4) :: ixy, ixy_fa
#endif
      
#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [rearrange_index_hex] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [rearrange_index_hex] in Mod_GetNode'
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

      if(.not.allocated(imap)) allocate(imap(nxy))
      if(.not.allocated(imapsol)) allocate(imapsol(nxy))
      imap=0
      imapsol=0

      n = 1
! modified by Tuan
!      do j = -1,ny_hex+2
!          do i = -3,nx_hex+4
!              ind = nodel_hex(i,j)
!              if (ind > 0) then
!                  imapsol(n) = ind
!                  n = n + 1
!              endif
!          enddo
!      enddo
      do j = 1,ny_hex+4
          do i = 1,nx_hex+8
              ind = nodel_hex(i,j)
              if (ind > 0) then
                  imapsol(n) = ind
                  n = n + 1
              endif
          enddo
      enddo

      do i = 1, nxy
          n = imapsol(i)
          imap(n) = i
      enddo
#ifdef tuan_tr_test
      if (flag_thfb) then
         if (.not.allocated(ixytoifa)) allocate(ixytoifa(nxy))
         ixytoifa = 0
         do ixy_fa = 1,nxy_fa
            ixy = i_fa(ixy_fa) 
            ixytoifa(ixy) = ixy_fa
         enddo
      endif
#endif
      
      
      RETURN
      END SUBROUTINE rearrange_index_hex

      SUBROUTINE Get_NodeVolume_hex
      USE Inc_TPEN, only: ndivhs, ntph
      IMPLICIT NONE
      INTEGER :: ixy, iz
      real(8) :: hexarea

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Get_NodeVolume_hex] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_NodeVolume_hex] in Mod_GetNode'
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

      if(.not.allocated(NodeVolume)) CALL Alloc( NodeVolume, Nxy, Nz )
      hexarea = GridSize_x(1)*GridSize_x(1)/(3d0**(1d0/2d0))/2d0/2d0*6d0 ! hexa_geometry
      do ixy = 1,nxy
         do iz = 1, nz
            !NodeVolume(ixy,iz) = hexarea*h_z(iz)
            NodeVolume(ixy,iz) = hexarea*MeshSize_z(iz)
         enddo
      enddo

      ndivhs=1
      ntph=6*ndivhs*ndivhs
      !nxskip=2

      RETURN
      END SUBROUTINE Get_NodeVolume_hex

      SUBROUTINE orderhex
      use Inc_Option, only: N_Group
      USE Inc_TPEN!, only: isymang, isolang, isymtype, isymmetry, nring,     &
                  !        icordys, icordye, icordxs, icordxe, hex_map,      &
                  !        ltox_hex, ltoy_hex, nxs, nxe, nys, nye, nrnx_hex, &
                  !        nassy, pbdv
      USE Mod_SolTPEN, only: premat

#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction
      use Mod_PinPow_Hex!, only: allocate_pinpow_hex
      use Inc_PinPow_Hex, only: pin_pow_hex_needed,saved_nodel_hex
      use Inc_File,       only: Name_INP, Len_INP
      ! -=-=-=-=-=-=-=-=-
#endif

      IMPLICIT NONE
      real(8) :: value(100)
      real(8) :: fac, fac2
      integer(4) :: iy, ix, i0, j0, inum, imvalx
      integer(4) :: icy, imval, ittt, i, j, l, item, m
      integer(4) :: nbdy, icel, isfc
      integer(4) :: id, iz, it, ih, iptmy, i2, iptyou, i3, idif
      integer(4) :: mp(2)
      data mp/2,1/
      logical(4) :: skip
      real(8) :: rt3, sqrt3, rsqrt3, hside
!      integer(4) :: ndivhs
      integer(4) :: nnz

#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction
      character(200)             :: filename
      integer                    :: ci,cj
      character(20)              :: print_name_tmp
      ! Tuan added
      INTEGER(4) :: n, ind



#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [orderhex] in Mod_GetNode'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      filename = Name_INP(1:Len_INP)//'_hff_input.dat'
      ! -=-=-=-=-=-=-=-=-
#endif

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [orderhex] in Mod_GetNode'
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

      nring = (Nx_1N/2)+mod(Nx_1N,2)
      !nxfc  = 8*(nring+1)
      !nyfc  = 4*(nring+1)
      !npr   = nz
      !if(isolang.EQ.360) then
      !do iy=icordys,-1
      !   do ix=icordxs,icordxe
      !      hex_map(ix,iy)=0
      !   enddo
      !enddo
      !do ix=icordxs,-1
      !   hex_map(ix,0)=0
      !enddo
      !endif
      i0=2*nring+1
      j0=nring+1

      if(i0.NE.1 .AND. icordxe/4.NE.nring) i0=i0-1

      inum=0
      imvalx=MOD(icordxe,4)
      do icy=1,3
         do iy=icordys,icordye,2
            imval=mod(iy,4)
            if(imvalx.eq.0) then
               if(imval.eq.0) then
                  ittt=icordxs
               else
                  ittt=icordxs-6
               endif
            else
               if(imval.eq.0) then
                  ittt=icordxs-6
               else
                  ittt=icordxs
               endif
            endif
            ittt=ittt+(icy-1)*4
! y-location on natural number assembly coordinate
            j=iy/2+j0
            do ix=ittt,icordxe,12
               if(hex_map(ix,iy).ne.0) then
                  inum=inum+1
! x-location on natural number assembly coordinate
                  i=ix/2+i0
                  l=inum
                  nodel_hex(i,j)=l
                  ltox_hex(l)=i
                  ltoy_hex(l)=j
                  if(hex_map(ix-4,iy).eq.0) nxs(j)=i
                  if(hex_map(ix+4,iy).eq.0) nxe(j)=i
                  iasytype(inum)=hex_map(ix,iy)
                  hex_map(ix,iy)=inum
               endif
            enddo
         enddo
      enddo
      !stop

#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction | save nodel_hex
      ! save FA arrangements in the core
      ci = 1
      cj = 1
! modified by Tuan - error in here:      Subscript #2 of the array NODEL_HEX has value -1 which is less than the lower bound of 1  
      do iy = -1,(ny_hex+2)
         ci = 1
         do ix = -3,(nx_hex+4)
            saved_nodel_hex(ci,cj) = nodel_hex(ix,iy)
            ci = ci + 1
         end do
         cj = cj + 1
      end do
      ! -=-=-=-=-=-=-=-=-
      if(.not.allocated(imapsol)) allocate(imapsol(nxy))
      n=1
      ind = 0
      do j = -1,ny_hex+2
          do i = -3,nx_hex+4
              ind = nodel_hex(i,j)
              if (ind > 0) then
                  imapsol(n) = ind
                  n = n + 1
              endif
          enddo
      enddo

#endif

      do j=1,ny
         nrnx_hex(j)=nxe(j)-nxs(j)+1
      enddo
      nassy=inum
! === Surface ordering === !
      inum=0
      do iy=icordys-1,icordye+1
         do ix=icordxs-2,icordxe+2
            item=iasurf(ix,iy)
            if(item.NE.0) then
               inum=inum+1
               iasurf(ix,iy)=inum
            endif
         enddo
      enddo
      nxsfc=inum
! === Point ordering === !
! set value at boundary point and point numbering
      do m=1,n_group
         if (if_hexbc) then
            value(m)=4964.d0*rt3*hside*hexbc/74165.d0
         else
            value(m)=4964.d0*rt3*hside*albedo(BC_LX)/74165.d0
         endif
      enddo
      if(hside.LT.5.0.OR.ndivhs.GE.2) then
         fac=1.
         fac2=1.
      else
         fac=1.-rt3/2.
         fac2=1.
      endif
      inum=0
      DO iy=icordys-1,icordye+1
         DO ix=icordxs-2,icordxe+2
            nbdy=iapoint(ix,iy)
            IF(nbdy.NE.0) THEN
               inum=inum+1
               IF(nbdy.EQ.1) THEN
                  DO m=1,n_group
                     pbdv(m,inum)=value(m)*fac2
                  ENDDO
               ENDIF
               IF(nbdy.EQ.2) THEN
                  DO m=1,n_group
                     pbdv(m,inum)=value(m)*fac*0.5
                  ENDDO
               ENDIF
               iapoint(ix,iy)=inum
            ENDIF
         ENDDO
      ENDDO
      nxpnt=inum
!     generate neighbor information
      DO iy=-nyfc,nyfc
         DO ix=-nxfc,nxfc
            layh(ix,iy)=0
            layp(ix,iy)=0
            lays(ix,iy)=0
         ENDDO
      ENDDO
      DO iy=icordys-1,icordye+1
         DO ix=icordxs-2,icordxe+2
            layh(ix,iy)=hex_map(ix,iy)
            layp(ix,iy)=iapoint(ix,iy)
            lays(ix,iy)=iasurf(ix,iy)
         ENDDO
      ENDDO
!     neighbor node linkage
!     1) initialization
      DO iy=icordys,icordye
         DO ix=icordxs,icordxe
            icel=hex_map(ix,iy)
            IF(icel.NE.0) THEN
               neignd(1,icel)=layh(ix-4,iy)
               neignd(2,icel)=layh(ix-2,iy-2)
               neignd(3,icel)=layh(ix+2,iy-2)
               neignd(4,icel)=layh(ix+4,iy)
               neignd(5,icel)=layh(ix+2,iy+2)
               neignd(6,icel)=layh(ix-2,iy+2)
               neigpt(1,icel)=layp(ix-2,iy+1)
               neigpt(2,icel)=layp(ix-2,iy-1)
               neigpt(3,icel)=layp(ix,iy-1)
               neigpt(4,icel)=layp(ix+2,iy-1)
               neigpt(5,icel)=layp(ix+2,iy+1)
               neigpt(6,icel)=layp(ix,iy+1)
               neigsfc(1,icel)=lays(ix-2,iy)
               neigsfc(2,icel)=lays(ix-1,iy-1)
               neigsfc(3,icel)=lays(ix+1,iy-1)
               neigsfc(4,icel)=lays(ix+2,iy)
               neigsfc(5,icel)=lays(ix+1,iy+1)
               neigsfc(6,icel)=lays(ix-1,iy+1)
               neigjin(1,icel)=4
               neigjin(2,icel)=5
               neigjin(3,icel)=6
               neigjin(4,icel)=1
               neigjin(5,icel)=2
               neigjin(6,icel)=3
               wtdhat(1,icel)=-1.
               wtdhat(2,icel)=-1.
               wtdhat(3,icel)=-1.
               wtdhat(4,icel)=1.
               wtdhat(5,icel)=1.
               wtdhat(6,icel)=1.
            ENDIF
         ENDDO
      ENDDO
!
      DO iy=icordys,icordye,2
         DO ix=icordxs-2,icordxe+2,2
            isfc=iasurf(ix,iy)
            IF(isfc.NE.0) THEN
               neigsnd(1,isfc)=layh(ix-2,iy)
               neigsnd(2,isfc)=layh(ix+2,iy)
               neigsnd(4,isfc)=4
               neigsnd(5,isfc)=1
            ENDIF
         ENDDO
      ENDDO
      DO iy=icordys-1,icordye+1,2
         DO ix=icordxs-1,icordxe+1,2
            isfc=iasurf(ix,iy)
            IF(isfc.NE.0) THEN
               neigsnd(1,isfc)=layh(ix-1,iy-1)
               neigsnd(2,isfc)=layh(ix+1,iy+1)
               neigsnd(4,isfc)=5
               neigsnd(5,isfc)=2
               IF(neigsnd(1,isfc).EQ.0) THEN
                  neigsnd(1,isfc)=layh(ix+1,iy-1)
                  neigsnd(4,isfc)=6
               ENDIF
               IF(neigsnd(2,isfc).EQ.0) THEN
                  neigsnd(2,isfc)=layh(ix-1,iy+1)
                  neigsnd(5,isfc)=3
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      DO isfc=1,nxsfc
         neigsnd(3,isfc)=12
         DO id=1,2
            IF(neigsnd(id,isfc).EQ.0) neigsnd(3,isfc)=mp(id)
         ENDDO
      ENDDO
!
      DO iz=1,nz
         neigz(1,iz)=iz-1
         neigz(2,iz)=iz+1
      ENDDO
      DO iz=1,nz
         neigsfcz(1,iz)=iz
         neigsfcz(2,iz)=iz+1
      ENDDO
      neigz(1,1)=0
      neigz(2,nz)=0
      DO iz=1,nz+1
         neigsndz(1,iz)=iz-1
         neigsndz(2,iz)=iz
         neigsndz(3,iz)=12
         neigsndz(4,iz)=2
         neigsndz(5,iz)=1
      ENDDO
      neigsndz(1,1)=0
      neigsndz(3,1)=2
      neigsndz(2,nz+1)=0
      neigsndz(3,nz+1)=1
!
!! Geometry
      rhzbar2=0d0
      if (nz /= 1)then
         do iz=1,nz
            do it=1,2
               nnz=neigz(it,iz)
               if (nnz.EQ.0) then
                  rhzbar2(it,iz)=1./hz(iz)/hz(iz)
               else
                  rhzbar2(it,iz)=2./hz(iz)/(hz(iz)+hz(nnz))
               endif
            enddo
         enddo
      endif
!
! === set parameters for point flux solver(pntslv.F)
      DO ih=1,nassy
         DO it=1,6
            iptmy=neigpt(it,ih)
            DO i2=1,6
               IF(i2.NE.it) THEN
                  iptyou=neigpt(i2,ih)
                  skip = .FALSE.
                  DO i3=1,neignpt(iptmy)
                     IF(neigppt(i3,iptmy).EQ.iptyou) THEN
                        skip = .TRUE.
                        EXIT
                     END IF
                  ENDDO
                  IF(.NOT.skip) THEN
                     neignpt(iptmy)=neignpt(iptmy)+1
                     neigppt(neignpt(iptmy),iptmy)=iptyou
                  END IF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DO ih=1,nassy
         DO it=1,6
            iptmy=neigpt(it,ih)
            DO i2=1,6
               IF(i2.NE.it) THEN
                  idif=i2-it
                  IF(idif.LT.0) idif=idif+6
                  iptyou=neigpt(i2,ih)
                  skip = .FALSE.
                  DO i3=1,neignpt(iptmy)
                     IF(neigppt(i3,iptmy).EQ.iptyou) THEN
                        imatid(idif,it,ih)=i3
                        skip = .TRUE.
                        EXIT
                     ENDIF
                  ENDDO
                  IF(.NOT.skip) WRITE(*,*)"Error when searching position of point"
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
! generation of core edit data
!
      icoreys=icordye
      icoreye=icordys
      icoref=icordxe
      DO iy=icordys,icordye,2
         icorexs(iy)=icordxe
         icorexe(iy)=icordxs
         DO ix=icordxs,icordxe,2
            IF(hex_map(ix,iy).NE.0) THEN
               IF(iy.LT.icoreys) icoreys=iy
               IF(iy.GT.icoreye) icoreye=iy
               IF(ix.LT.icorexs(iy)) icorexs(iy)=ix
               IF(ix.GT.icorexe(iy)) icorexe(iy)=ix
               IF(ix.LT.icoref) icoref=ix
            ENDIF
         ENDDO
      ENDDO
!
! generation of core edit data
!
      icoreyps=icordye+2
      icoreype=icordys-2
      icorepf=icordxe+2
      DO iy=icordys-1,icordye+1,2
         icorexps(iy)=icordxe+2
         icorexpe(iy)=icordxs-2
         DO ix=icordxs-2,icordxe+2,2
            IF(iapoint(ix,iy).NE.0) THEN
               IF(iy.LT.icoreyps) icoreyps=iy
               IF(iy.GT.icoreype) icoreype=iy
               IF(ix.LT.icorexps(iy)) icorexps(iy)=ix
               IF(ix.GT.icorexpe(iy)) icorexpe(iy)=ix
               IF(ix.LT.icorepf) icorepf=ix
            ENDIF
         ENDDO
      ENDDO
!
! generation of core edit data
!
      icoreyss=icordye+1
      icoreyse=icordys-1
      icoresf=icordxe+1
      DO iy=icordys-1,icordye+1
         icorexss(iy)=icordxe+1
         icorexse(iy)=icordxs-1
         DO ix=icordxs-2,icordxe+2
            IF(iasurf(ix,iy).NE.0) THEN
               IF(iy.LT.icoreyss) icoreyss=iy
               IF(iy.GT.icoreyse) icoreyse=iy
               IF(ix.LT.icorexss(iy)) icorexss(iy)=ix
               IF(ix.GT.icorexse(iy)) icorexse(iy)=ix
               IF(ix.LT.icoresf) icoresf=ix
            ENDIF
         ENDDO
      ENDDO
!
! generation of core edit data for fuel assembly
!
      icoreyfs=icordye
      icoreyfe=icordys
      icoreff=icordxe
      DO iy=icordys,icordye,2
         icorexfs(iy)=icordxe
         icorexfe(iy)=icordxs
         DO ix=icordxs,icordxe,2
            IF(hex_map(ix,iy).NE.0) THEN
               IF(iffuel(iasytype(hex_map(ix,iy)))) THEN
                  IF(iy.LT.icoreyfs) icoreyfs=iy
                  IF(iy.GT.icoreyfe) icoreyfe=iy
                  IF(ix.LT.icorexfs(iy)) icorexfs(iy)=ix
                  IF(ix.GT.icorexfe(iy)) icorexfe(iy)=ix
                  IF(ix.LT.icoreff) icoreff=ix
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
! Generation of information for Condensation of Matrix cmat
!
      CALL premat
!

#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction
      if (pin_pow_hex_needed) then
         call read_hff_hex(filename)
      end if
      ! -=-=-=-=-=-=-=-=-
#endif

      RETURN
      END SUBROUTINE orderhex

      REAL(8) FUNCTION albedo(ibc)

      INTEGER(4) :: ibc

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [albedo] in Mod_GetNode'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      select case(ibc)
      case (1)
         albedo=0.0 ! reflector boundary condition
      case (2)
         albedo=0.5 ! zero incoming current
      case (3)
         albedo=1.0 ! zero flux
      end select

      RETURN
      END FUNCTION

      SUBROUTINE readlay
      USE Inc_TPEN!, only: icordys, icordye, nxrow, nring, hex_map_index, &
                  !        icordxs, icordxe, hex_map
      IMPLICIT NONE

      !INTEGER(4), INTENT(in) :: ifile, isyminp
      INTEGER(4) :: iy1, nhalfs, ir
      INTEGER(4) :: nxcol, ic, nhalfcol, ix1, imodv
      INTEGER(4) :: iy2, iy, ix

      !INTEGER(4) :: nn
      INTEGER(4) :: k

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [readlay] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [readlay] in Mod_GetNode'
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

      iy1 = 0
      icordys = 0
      icordye = 2*nxrow-2
      !if (isyminp.EQ.360) then
      nhalfs=nxrow/2
      iy1=-2*nhalfs
      icordys=iy1
      icordye=-iy1
      icordxs=0
      icordxe=0
      ir=0
      k=1

      do while (ir.LT.nxrow)
         ir=ir+1
         nxcol=sum(hex_map_index(ir,:))
         nhalfcol=nxcol/2
         ix1=-4*nhalfcol
         imodv=mod(nxcol,2)
         if(imodv.EQ.0) ix1=ix1+2
         do ic=1,nxcol
            if(ix1.LT.icordxs) icordxs=ix1
            if(ix1.GT.icordxe) icordxe=ix1
            iy2=(ix1+iy1)/2
            if(iy2.GT.icordye) icordye=iy2
            hex_map(ix1,iy1)=hex_map_mid(ir,ic)
            ix1=ix1+4
         enddo
         iy1=iy1+2
         continue
      enddo
      icordxs=-icordxe
      icordys=-icordye

! surface and point generation for 360 degree
      do iy=icordys,icordye
         do ix=icordxs,icordxe
            if (hex_map(ix,iy).NE.0) then
               iapoint(ix-2,iy-1)=iapoint(ix-2,iy-1)+1
               iapoint(ix,iy-1)=iapoint(ix,iy-1)+1
               iapoint(ix+2,iy-1)=iapoint(ix+2,iy-1)+1
               iapoint(ix-2,iy+1)=iapoint(ix-2,iy+1)+1
               iapoint(ix,iy+1)=iapoint(ix,iy+1)+1
               iapoint(ix+2,iy+1)=iapoint(ix+2,iy+1)+1
               iasurf(ix-1,iy-1)=1
               iasurf(ix+1,iy-1)=1
               iasurf(ix-2,iy)=1
               iasurf(ix+2,iy)=1
               iasurf(ix-1,iy+1)=1
               iasurf(ix+1,iy+1)=1
            endif
         enddo
      enddo
#ifdef tuan_tr_test
      if (if_hex_cr) then
         call hex_crbank
      endif
#endif
      RETURN
      END SUBROUTINE readlay

      SUBROUTINE scanlay
      USE Inc_TPEN
      USE Inc_Option,   ONLY: n_group
      !USE Inc_maXS, ONLY: nu_maXS_f_3D
!      USE Inc_Control,  ONLY: nmultht,nupdcyt,ntnodal
      USE Inc_3D,       ONLY: flux
      USE Inc_Solver,   ONLY: betap
      USE Inc_maXS
      USE Inc_Geometry, ONLY: Gridsize_z
      USE Inc_XS_File,  ONLY: N_XS_Table
#ifdef tuan_tr_test
      use Inc_XS_File,  only: ref_smat
#endif
      
      IMPLICIT NONE
! scan lay(loading pattern) card and generates maximum value for radial grid
      !LOGICAL, SAVE :: sfirst = TRUE
      !CHARACTER(len=65) :: amesg
      !CHARACTER(len=4) ::  astr
      INTEGER(4) :: i,ic, ir,ix
      INTEGER(4) :: ibegh, iendh, jbegh, jendh
      INTEGER(4) :: inums, imodv,inump, inumh, ix1, iy, iy1, iy2
      INTEGER(4) :: nhalfs,nxcol,ionerow,nhalfcol
      DIMENSION ionerow(1024)
      INTEGER(4) :: la, Ix_1N, Iy_1N, j
      !INTEGER(4) :: nring
      INTEGER(4) :: iz, I_Tab
      !LOGICAL(1) :: tmp_flag

      INTEGER(4) :: ig,  ih, ip
      REAL(8)    :: rdat
      INTEGER(4) :: it
!      INTEGER(4) :: n
      INTEGER(4) :: m, md
      INTEGER(4) :: mf, ms

!      INTEGER(4) :: n1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [scanlay] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [scanlay] in Mod_GetNode'
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

      nring = (Nx_1N/2)+mod(Nx_1N,2)
      nxrow=2*nring+1
      nat=0 !initialize number of assembly types to zero
! read radial core configuration
      iy1=0
      icordys=0
      icordye=2*nxrow-2
      nhalfs=nxrow/2
      iy1=-2*nhalfs
      icordys=iy1
      icordye=-iy1
      icordxs=0
      icordxe=0
      ir=0
      nassyinp=0

!!
      nxrow=Nx_1N
!      n1=0
      out: DO WHILE(ir.LT.nxrow)
         ir=ir+1
         nxcol=sum(hex_map_index(ir,:))
         !READ(oneline,*) (ionerow(ic),ic=1,nxcol)
         nassyinp=nassyinp+nxcol
         DO ic=1,nxcol
            ionerow(ic)=hex_map_mid(ir,ic)
            nat=MAX(nat,hex_map_mid(ir,ic))
         ENDDO
         !IF(isyminp.EQ.360) THEN
         nhalfcol=nxcol/2
         ix1=-4*nhalfcol
         imodv=MOD(nxcol,2)
         IF(imodv.EQ.0) ix1=ix1+2
         !ENDIF
         DO ic=1,nxcol
            IF(ix1.LT.icordxs) icordxs=ix1
            IF(ix1.GT.icordxe) icordxe=ix1
            iy2=(ix1+iy1)/2
            IF(iy2.GT.icordye) icordye=iy2
            layh(ix1,iy1)=ionerow(ic)
            !if (layh(ix1,iy1).NE.0) n1=n1+1
            ix1=ix1+4
         ENDDO
         iy1=iy1+2
      ENDDO out
      nxrow=ir
      icordxs=-icordxe
      icordys=-icordye
!
! generate load patterb for 360 degree
!
!
! surface and point generation for 360 degree
      DO iy=icordys,icordye
         DO ix=icordxs,icordxe
            IF(layh(ix,iy).NE.0) THEN
               layp(ix-2,iy-1)=1
               layp(ix,iy-1)=1
               layp(ix+2,iy-1)=1
               layp(ix-2,iy+1)=1
               layp(ix,iy+1)=1
               layp(ix+2,iy+1)=1
               lays(ix-1,iy-1)=1
               lays(ix+1,iy-1)=1
               lays(ix-2,iy)=1
               lays(ix+2,iy)=1
               lays(ix-1,iy+1)=1
               lays(ix+1,iy+1)=1
            ENDIF
         ENDDO
      ENDDO
! hex node generation and numbering
      !! isolang .LT. 360
      !DO iy=icordys-1,-1
      !   DO ix=icordxs-2,icordxe+2
      !      layh(ix,iy)=0
      !      layp(ix,iy)=0
      !      lays(ix,iy)=0
      !   ENDDO
      !ENDDO
      !DO ix=icordxs-2,-1
      !   layh(ix,0)=0
      !ENDDO
!
      inumh=0
      inums=0
      inump=0
!
      nx_hex=4*nring+1
      ny_hex=nx_hex
!
      ibegh=0
      iendh=0
      jbegh=0
      jendh=0
!
      DO iy=icordys-1,icordye+1
         DO ix=icordxs-2,icordxe+2
            IF(layh(ix,iy).NE.0) THEN
               inumh=inumh+1
               ibegh=MIN(ibegh,ix)
               iendh=MAX(iendh,ix)
               jbegh=MIN(jbegh,iy)
               jendh=MAX(jendh,iy)
            ENDIF
            IF(layp(ix,iy).NE.0) THEN
               inump=inump+1
            ENDIF
            IF(lays(ix,iy).NE.0) THEN
               inums=inums+1
            ENDIF
         ENDDO
      ENDDO
      nassy=inumh
      ncorn_hex=inump
      nsurf=inums
!
      if (.not.allocated(pbdv))      allocate(pbdv(n_group,ncorn_hex))
      if (.not.allocated(neignd))    allocate(neignd(6,nassy))
      if (.not.allocated(neigpt))    allocate(neigpt(6,nassy))
      if (.not.allocated(neigsfc))   allocate(neigsfc(6,nassy))
      if (.not.allocated(neigjin))   allocate(neigjin(6,nassy))
      if (.not.allocated(wtdhat))    allocate(wtdhat(ntph,nassy))
      if (.not.allocated(neigsnd))   allocate(neigsnd(5,nsurf))
      if (.not.allocated(neigz))     allocate(neigz(2,nz))
      if (.not.allocated(neigsfcz))  allocate(neigsfcz(2,nz))
      if (.not.allocated(neigsndz))  allocate(neigsndz(5,nz+1))
      if (.not.allocated(imatid))    allocate(imatid(5,6,nassy))
      !if (.not.allocated(inptr))     allocate(inptr(0:6,nassy))
      if (.not.allocated(neignpt))   allocate(neignpt(ncorn_hex))
      if (.not.allocated(neigppt))   allocate(neigppt(12,ncorn_hex))
      if (.not.allocated(iffuel))    allocate(iffuel(nat))
      if (.not.allocated(ineigcond)) allocate(ineigcond(0:6,0:6,nassy))
      if (.not.allocated(ipntr))     allocate(ipntr(0:6,nassy))
      if (.not.allocated(ilubnd))    allocate(ilubnd(nassy))
      if (.not.allocated(iastopnt))  allocate(iastopnt(nassy,nassy))
      if (.not.allocated(iasypr))    allocate(iasypr(Nz,nat))
      if (.not.allocated(xsid))      allocate(xsid(1:Nxy,Nz))
      if (.not.allocated(hx))        allocate(hx(0:nx_hex+1))
      if (.not.allocated(hy))        allocate(hy(0:ny_hex+1))
      if (.not.allocated(hz))        allocate(hz(0:nz+1))
      if (.not.allocated(ktokb))     allocate(ktokb(nz))
      if (.not.allocated(cmat))      allocate(cmat(N_group,8,nassy,nz))
      if (.not.allocated(dcmat))     allocate(dcmat(n_group*n_group,nassy,nz))
      if (.not.allocated(cmatd))     allocate(cmatd(n_group,8,nassy,nz))
      if (.not.allocated(dcmatd))    allocate(dcmatd(n_group*n_group,nassy,nz))
      if (.not.allocated(xlufac))    allocate(xlufac(n_group,6,nassy,nz))
      if (.not.allocated(alxrf))     allocate(alxrf(n_group))
      if (.not.allocated(alxr))      allocate(alxr(n_group))
      if (.not.allocated(hflx))      allocate(hflx(n_group,nassy,nz))
      if (.not.allocated(pflx))      allocate(pflx(n_group,ncorn_hex,nz))
      if (.not.allocated(mge))       allocate(mge(n_group))
      if (.not.allocated(mgb))       allocate(mgb(n_group))
      if (.not.allocated(dfd))       allocate(dfd(n_group,nsurf,nz))
      if (.not.allocated(dfdz))      allocate(dfdz(n_group,nassy,nz+1))
      if (.not.allocated(dhat))      allocate(dhat(n_group,nsurf,nz))
      if (.not.allocated(dhatz))     allocate(dhatz(n_group,nassy,nz+1))
      if (.not.allocated(betaphis))  allocate(betaphis(n_group,nsurf,nz))
      if (.not.allocated(betaphisz)) allocate(betaphisz(n_group,nassy,nz+1))
      if (.not.allocated(betaphisd)) allocate(betaphisd(n_group,nsurf,nz))
      if (.not.allocated(betaphiszd))allocate(betaphiszd(n_group,nassy,nz+1))
      if (.not.allocated(fhflx))     allocate(fhflx(n_group,nassy,nz))
      if (.not.allocated(fohflx))    allocate(fohflx(n_group,nassy,nz))
      if (.not.allocated(hflxf))     allocate(hflxf(n_group,nassy,nz))
      if (.not.allocated(aflx))      allocate(aflx(n_group,ntph,nassy,nz))
      if (.not.allocated(xmom))      allocate(xmom(n_group,ntph,nassy,nz))
      if (.not.allocated(ymom))      allocate(ymom(n_group,ntph,nassy,nz))
      if (.not.allocated(zmom1))     allocate(zmom1(n_group,nassy,nz))
      if (.not.allocated(zmom2))     allocate(zmom2(n_group,nassy,nz))
      if (.not.allocated(fcnto))     allocate(fcnto(n_group,6,nassy,nz))
      if (.not.allocated(focnto))    allocate(focnto(n_group,6,nassy,nz))
      if (.not.allocated(cnto))      allocate(cnto(n_group,ntph,nassy,nz))
      if (.not.allocated(fcntzo))    allocate(fcntzo(n_group,2,nassy,nz))
      if (.not.allocated(focntzo))   allocate(focntzo(n_group,2,nassy,nz))
      if (.not.allocated(cntzo))     allocate(cntzo(n_group,2,nassy,nz))
      if (.not.allocated(alzl))      allocate(alzl(n_group))
      if (.not.allocated(alzlf))     allocate(alzlf(n_group))
      if (.not.allocated(alphaz))    allocate(alphaz(n_group))
      if (.not.allocated(alzr))      allocate(alzr(n_group))
      if (.not.allocated(alzrf))     allocate(alzrf(n_group))
      if (.not.allocated(reflratf))  allocate(reflratf(n_group))
      if (.not.allocated(reflratzbf))allocate(reflratzbf(n_group))
      if (.not.allocated(reflratztf))allocate(reflratztf(n_group))
      if (.not.allocated(reflratzf)) allocate(reflratzf(n_group))
      if (.not.allocated(pflxt))     allocate(pflxt(ncorn_hex))
      if (.not.allocated(xsadf))     allocate(xsadf(n_group,6,nxy,nz))
      if (.not.allocated(atleak))    allocate(atleak(n_group,nassy,nz))
      if (.not.allocated(srcz))      allocate(srcz(n_group,nassy,nz))
      if (.not.allocated(igc))       allocate(igc(n_group))
      if (.not.allocated(adft))      allocate(adft(n_group,nxy,nz))
      if (.not.allocated(adfb))      allocate(adfb(n_group,nxy,nz))
      if (.not.allocated(xssf))      allocate(xssf(n_group-1,2:n_group,nxy,nz))
      if (.not.allocated(dsum))      allocate(dsum(n_group,nassy,nz))
      if (.not.allocated(rhzbar2))   allocate(rhzbar2(2,nz))
      if (.not.allocated(xstd))      allocate(xstd(n_group,nxy,nz))
      if (.not.allocated(adfpt))     allocate(adfpt(n_group,6,nxy,nz))
      if (.not.allocated(iscatib))   allocate(iscatib(n_group))
      if (.not.allocated(iscatie))   allocate(iscatie(n_group))
      if (.not.allocated(iscatob))   allocate(iscatob(n_group))
      if (.not.allocated(iscatoe))   allocate(iscatoe(n_group))
      if (.not.allocated(betap))     allocate(betap(nxy,nz))
      if (.not.allocated(dumrv))     allocate(dumrv(n_group,nassy))
      if (.not.allocated(dumrs))     allocate(dumrs(n_group,nassy))
      if (.not.allocated(fluxf))     allocate(fluxf(nxy,nz,n_group))
#ifdef tuan_fr
      if (.not.allocated(fluxf_old)) allocate(fluxf_old(nxy,nz,n_group))
#endif

#ifdef tuan_tr_test
      if (.not.allocated(dhatd))     allocate(dhatd(n_group,nsurf,nz))
      if (.not.allocated(dhatz0))    allocate(dhatzd(n_group,nassy,nz+1))
      if (.not.allocated(dhat0))     allocate(dhat0(n_group,nsurf,nz))
      if (.not.allocated(dhatz0))    allocate(dhatz0(n_group,nassy,nz+1))
      if (.not.allocated(rvdeltf))   allocate(rvdeltf(n_group,nxy,nz))
      rvdeltf = 0d0
      if (.not.allocated(aphif))     allocate(aphif(n_group,nxy,nz))
      aphif=0d0
#endif

#ifdef tuan_fr_tdep
       if (.not.allocated(tfluxf))     allocate(tfluxf(6,nxy,nz,n_group))
            tfluxf = 0

#endif
      dumrv=0d0
      dumrs=0d0
!! Geometry
      rhzbar2=0d0
      fluxf=0d0
!      if (nz /= 1)then
!         do iz=1,nz
!            do it=1,2
!               nnz=neigz(it,iz)
!            enddo
!         enddo
!      endif
!! XS
      xssf=0d0
      !do ixy = 1,Nxy
      !   do md=2,n_group
      !      do m=1,md-1
      !         xssf(m,md,ixy,:)=maXS_scat_3D(m,md,imap(ixy),:)
      !      enddo
      !      if ((md+1)<n_group) then
      !      do m=md+1,n_group
      !         xssf(m-1,md,i,:)=maXS_scat_3D(m,md,imap(ixy),:)
      !      enddo
      !      endif
      !   enddo
      !enddo
!! define parameter
!      nodalcy =nmultht*nupdcyt ! assume nodal calcualtion
!! scattering (in-scattering)
!! --- base two-group:: need to modify
      iscatie=0
      iscatib=0
      iscatib(1)=1
      iscatie(1)=0
      do ig = 2,n_group
         iscatib(ig)=1
         iscatie(ig)=n_group-1
      enddo
      do ig = 1,n_group-1
         iscatob(ig)=ig+1
         iscatoe(ig)=n_group
      enddo
      iscatob(n_group)=n_group
      iscatoe(n_group)=n_group-1
!! if you need multigroup calculation, please modify
      if (n_group>2) then
         iscatib(2:n_group)=n_group
         iscatie(1:n_group-1)=2
         iscatib(1)=1
         iscatie(1)=0
         !DO Ixy = 1,Nxy
         !   DO Iz=1,Nz
         !      DO md=2,n_group
         !         DO m=1,md-1
         !            IF( maXs_scat_3D(m,md,Ixy,Iz).NE.0) THEN
         !               iscatib(md)=MIN(m,iscatib(md))
         !               EXIT
         !            ENDIF
         !         ENDDO
         !         iscatie(md)=md-1
         !         DO m=n_group,md+1,-1
         !            ms=m-1
         !            IF(maXs_scat_3D(ms,md,Ixy,Iz).NE.0) THEN
         !               iscatie(md)=MAX(m,iscatie(md))
         !               EXIT
         !            ENDIF
         !         ENDDO
         !      ENDDO
         !   ENDDO
         !ENDDO
         DO I_Tab = 1,N_XS_Table
            DO md=2,n_group
               DO m=1,md-1
                  IF( .TRUE. ) THEN
                  !IF( Ref_smat(I_Tab,m,md).NE.0) THEN
                     iscatib(md)=MIN(m,iscatib(md))
                     EXIT
                  ENDIF
               ENDDO
               iscatie(md)=md-1
               DO m=n_group,md+1,-1
                  ms=m-1
                  IF( .TRUE. ) THEN
                  !IF( Ref_smat(I_Tab,ms,md).NE.0) THEN
                     iscatie(md)=MAX(m,iscatie(md))
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
#ifdef tuan_tr_test
!          write(*,*) "it comes here: scanlay"
           do m = 1, n_group
              iscatob(m) = m+1
              do md = 2,m-1
                 ms = m - 1
                 if (Ref_smat(I_Tab,ms,md).NE.0) then 
                    iscatob(m) = min(md,iscatob(m))
                    exit
                 endif
              enddo
              do md = n_group, m+1, -1
                 if (Ref_smat(I_Tab,ms,md).NE.0) then
                    iscatoe(m) = max(md,iscatoe(m))
                    exit
                 endif
              enddo
           enddo
#endif
         ENDDO
      endif
!! initialization
      betap=1d0
      adfpt=1d0
      xstd=0d0
      dsum=0d0
      hflx=0d0
      pflx=0d0
      fhflx=0d0
      fohflx=0d0
      hflxf=0d0
      aflx=0d0
      xmom=0d0
      ymom=0d0
      zmom1=0d0
      zmom2=0d0
      fcnto=0d0
      focnto=0d0
      cnto=0d0
      fcntzo=0d0
      focntzo=0d0
      cntzo=0d0
      alzl=0d0
      alzlf=0d0
      alphaz=0d0
      alzr=0d0
      alzrf=0d0
      alxr=0d0
      alxrf=0d0
      reflratzf=0d0
      reflratf=0d0
      reflratzbf=0d0
      reflratztf=0d0
      pflxt=0d0
      xsadf=1d0     !! ??? check the value ???
      xsadf(:,:,:,:) =1d0
      do ig = 1,n_group
         do it = 1,6
            do ih = 1,nxy
               do iz = 1,nz
                  xsadf(ig,it,ih,iz)=1d0
               enddo
            enddo
         enddo
      enddo
      atleak=0d0
      srcz=0d0
      adft=1d0
      adfb=1d0
!! initialization
      rdat=1./74165.d0
      chlval(1)=-22230.*rdat
      chlval(2)=2590.*rdat
      chlval(3)=-2187.*rdat
      chlval(4)=58968.*rdat
      chlval(5)=-600.*rdat
      chlval(6)=4704.*rdat
      chlval(7)=-10512.*rdat
!! igc ! Need to modify: current mode for two-group problem
      igc(1) = 1
      igc(2) = 2
!! albedo
      do ig=1,n_group
         alzlf(ig)=albedo(bc_lz)
         alzl(ig) =albedo(bc_lz)
         alzrf(ig)=albedo(bc_rz)
         alzr(ig) =albedo(bc_rz)
         if(if_hexbc) then
            alxrf(ig)=hexbc
            alxr(ig) =hexbc
         else
            alxrf(ig)=albedo(bc_lx)
            alxr(ig) =albedo(bc_lx)
         endif
         reflratf(ig)  =(1-2*alxrf(ig))/(1+2*alxrf(ig))
         reflratzbf(ig)=(1-2*alzlf(ig))/(1+2*alzlf(ig))
         reflratztf(ig)=(1-2*alzrf(ig))/(1+2*alzrf(ig))
      enddo
!! albedo
      do ig=1,n_group
         do iz =1,nz
            do ih =1,nassy
               hflx(ig,ih,iz)=flux(1,1,ig)
            enddo
         enddo
      enddo
!stop
      mge(1) = 1 ; mge(2) = 2
      mgb(1) = 1 ; mgb(2) = 2
      ! mge(1), mge(2)
      ! mgb(1), mgb(2)
      !! for multi-group
      mgb(1) = 1 ; mgb(2) = n_group/2+1
      mge(1) = n_group/2 ; mge(2) = n_group
! === multi-group
! vlaue: igc, mge, mgb
! 2-group :: ng2=2
      do m = 1,2
         do mf = mgb(m), mge(m)
            igc(mf)=m
         enddo
      enddo

#ifdef tuan_tr_test

#else
! modify value: pflx
! pflx(ig,ip,iz)=phi0f(igc(ig))
! === multi-group
      do ig=1,n_group
         do iz=1,nz
            do ip = 1,ncorn_hex
               pflx(ig,ip,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
            enddo
         enddo
      enddo

      do iz=1,nz
         do ih=1,nassy
            do it=1,ntph
               do ig=1,n_group
                  cnto(ig,it,ih,iz)=0.25*flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
               enddo
            enddo
            do ig=1,n_group
               cntzo(ig,1,ih,iz)=0.25*flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
               cntzo(ig,2,ih,iz)=0.25*flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
            enddo
         enddo
      enddo
      do iz=1,nz
         do ih=1,nassy
            do ig=1,n_group
               fhflx(ig,ih,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)/max(flux(1,1,ig),1d-30)
               fohflx(ig,ih,iz)=fhflx(ig,ih,iz)
               hflxf(ig,ih,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
               do it=1,6
                  aflx(ig,it,ih,iz)=flux(1,1,igc(ig))/(mge(igc(ig))-mgb(igc(ig))+1)
                  fcnto(ig,it,ih,iz)=fhflx(ig,ih,iz)
                  focnto(ig,it,ih,iz)=fhflx(ig,ih,iz)
               enddo
               do it=1,2
                  fcntzo(ig,it,ih,iz)=fhflx(ig,ih,iz)
                  focntzo(ig,it,ih,iz)=fhflx(ig,ih,iz)
               enddo
            enddo
         enddo
      enddo
!stop
#endif

      betaphis=0d0
      betaphisz=0d0
      betaphisd=0d0
      betaphiszd=0d0
      dhat=0d0
      dhatz=0d0
#ifdef tuan_tr_test
      dhatd=0d0
      dhatzd=0d0
      dhat0=0d0
      dhatz0=0d0
#endif
      dfd=0d0
      dfdz=0d0
      ntnodal=0
!      alxr=0d0
!      alxrf=0d0
      xlufac = 0d0
      cmat = 0d0
      dcmat = 0d0
      cmatd = 0d0
      dcmatd = 0d0
      ktokb = 0
      hz = 0d0
      hy = 0d0
      hx = 0d0
      pbdv = 0d0
      neignd = 0
      neigpt  = 0
      neigsfc = 0
      neigjin = 0
      wtdhat = 0d0
      neigsnd = 0
      neigz = 0
      neigsfcz = 0
      neigsndz = 0
      imatid = 0
      !inptr = 0
      neignpt = 0
      neigppt = 0
      iffuel = .false.
      ineigcond = 0
      ipntr = 0
      ilubnd = 0
      iastopnt = 0
      iasypr = 0
      xsid = 0
!! Geometry
      do iz=1,nz
         hz(iz)=Gridsize_z(iz)
      enddo
!
! generation of iffuel
      DO j=1,Ny
         Iy_1N=Iy_4Nto1N(j)
         DO i=Ix_Start_y(j),Ix_End_y(j)
            Ix_1N=Ix_4Nto1N(i)
            la = IxIy_1NToIxy_1N(Ix_1N,Iy_1N)
#ifdef tuan_fr_crm
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) == 2) THEN
#else
            IF ( I_FARF_1N_2D(Ix_1N, Iy_1N) > 1) THEN
#endif
               iffuel(I_LP_1N(la))=.true.
            END IF
         END DO
      END DO
!
      nasyx=(iendh-ibegh)/2+1
      nasyy=(jendh-jbegh)/2+1
      nx_hex=nasyx
      ny_hex=nasyy
      nxy_hex=nassy
      nchan_hex=nassy
      npr=nz
!
      RETURN
      END SUBROUTINE scanlay

      SUBROUTINE hex_map_make
      USE Inc_TPEN
      USE Inc_Option,    ONLY: N_group
!      USE Inc_Transient, ONLY: Flag_Transient
#ifdef tuan_tr_test
      use Inc_Option, only: OPT_Mode
#endif
#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction
      use Inc_PinPow_Hex, only: saved_nodel_hex
      ! -=-=-=-=-=-=-=-=-
#endif

      IMPLICIT NONE

      INTEGER(4) :: nn, nr
      INTEGER(4) :: i, j
!      INTEGER(4) :: nxfc, nyfc, nring
      INTEGER(4) :: nx2m1, nx2
      INTEGER(4) :: a,b,k,n
      INTEGER(4) :: ix, iy
      INTEGER(4) :: iy1
!      INTEGER(4) :: nmax
      !LOGICAL(1) :: flag

      INTEGER(4),ALLOCATABLE,DIMENSION(:,:) :: hex_map_mid_b


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [hex_map_make] in Mod_GetNode'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [hex_map_make] in Mod_GetNode'
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

!    make map hex index
!    example: 
!   1 1 0
!   1 1 1
!   0 1 1
!
      nn = Nx_1N
      nr = (Nx_1N+1)/2
      IF (.not. allocated(hex_map_index)) ALLOCATE(hex_map_index(nn,nn))
      do i = 1,nn
         do j = 1,nn
            if ((i<nr) .and. (j<nr+i)) then
               hex_map_index(i,j) = 1
            else if ((i>nr) .and. (j>i-nr)) then
               hex_map_index(i,j) = 1
            else if (i==nr) then
               hex_map_index(i,j) = 1
            else
               hex_map_index(i,j) = 0
            endif
         enddo
      enddo
!      stop

      nring = (Nx_1N/2)+mod(Nx_1N,2)
      nxfc  = 8*(nring+1)
      nyfc  = 4*(nring+1)
      nx2=Nx+Ny
      nx2m1=nx2-1
      nzp1 = Nz+1
      nx_hex=4*nring+1
      ny_hex=nx_hex
#ifdef tuan_fr
      npr =Nz
#endif

      if (.not.allocated(chi))       allocate(chi(n_group))
      if (.not.allocated(hex_map))   allocate(hex_map(-nxfc:nxfc,-nyfc:nyfc))
      if (.not.allocated(hex_map_mid_b))allocate(hex_map_mid_b(nn,nn))
      if (.not.allocated(hex_map_mid))allocate(hex_map_mid(nn,nn))
      if (.not.allocated(iapoint))   allocate(iapoint(-nxfc:nxfc,-nyfc:nyfc))
      if (.not.allocated(iasurf))    allocate(iasurf(-nxfc:nxfc,-nyfc:nyfc))
      if (.not.allocated(ltox_hex))  allocate(ltox_hex(-nx2m1:nxy+nx2+1))
      if (.not.allocated(ltoy_hex))  allocate(ltoy_hex(-nx2m1:nxy+nx2+1))
! modified by Tuan
!      if (.not.allocated(nodel_hex)) allocate(nodel_hex(1:nx_hex+8,1:ny_hex+4))
      if (.not.allocated(nodel_hex)) allocate(nodel_hex(-3:nx_hex+4,-1:ny_hex+2))
!      CALL dmalloc0(nodel_hex,-3,nx_hex+4,-1,ny_hex+2)  !radial node index
!      CALL dmalloc0(nodel_hex,1,nx+8,1,ny+4)  !radial node index
      if (.not.allocated(nxs))       allocate(nxs(ny_hex+1))
      if (.not.allocated(nxe))       allocate(nxe(ny_hex+1))
      if (.not.allocated(nys))       allocate(nys(nx_hex+1))
      if (.not.allocated(nye))       allocate(nye(nx_hex+1))
      if (.not.allocated(iasytype))  allocate(iasytype(0:nxy_1N))
      if (.not.allocated(nrnx_hex))  allocate(nrnx_hex(ny+1))
      if (.not.allocated(layh))      allocate(layh(-nxfc:nxfc,-nyfc:nyfc))
      if (.not.allocated(layp))      allocate(layp(-nxfc:nxfc,-nyfc:nyfc))
      if (.not.allocated(lays))      allocate(lays(-nxfc:nxfc,-nyfc:nyfc))
      if (.not.allocated(icorexs))   allocate(icorexs(-nyfc:nyfc))
      if (.not.allocated(icorexe))   allocate(icorexe(-nyfc:nyfc))
      if (.not.allocated(icorexps))  allocate(icorexps(-nyfc:nyfc))
      if (.not.allocated(icorexpe))  allocate(icorexpe(-nyfc:nyfc))
      if (.not.allocated(icorexss))  allocate(icorexss(-nyfc:nyfc))
      if (.not.allocated(icorexse))  allocate(icorexse(-nyfc:nyfc))
      if (.not.allocated(icorexfs))  allocate(icorexfs(-nyfc:nyfc))
      if (.not.allocated(icorexfe))  allocate(icorexfe(-nyfc:nyfc))
      if (.not.allocated(wtass))     allocate(wtass(nxy_1N))
      if (.not.allocated(wtasssum))  allocate(wtasssum(nxy_1N))
      if (.not.allocated(izpr))      allocate(izpr(0:nzp1))
      if (.not.allocated(iprcomp))   allocate(iprcomp(0:Nxy_1N,npr))
      if (.not.allocated(ltolb))     allocate(ltolb(0:Nxy_1N))

      if (.not.allocated(sefve))     allocate(sefve(n_group,-nx2m1:nxy,0:nz))
      if (.not.allocated(xschi))     allocate(xschi(Nxy,Nz,N_group))
      if (.not.allocated(xschid))    allocate(xschid(Nxy,Nz,N_group))
      ! Need to modify for transient calculation
      !if (Flag_Transient) then
      !   if (.not.allocated(srcf))   allocate(srcf(n_group,nxy,nz))
      !   srcf=0d0
      !   srcf=src
      !endif
#ifdef tuan_tr_test
      if (OPT_Mode == 3) then
         if (.not.allocated(srcf))   allocate(srcf(n_group,nxy,nz))
         srcf=0d0
         !srcf=src
      endif
#endif
#ifdef siarhei_ppr
      ! Siarhei pin power reconstruction | allocate positive (x,y) nodel_hex
      if (.not.allocated(saved_nodel_hex)) allocate(saved_nodel_hex(nx_hex+4+4,ny_hex+2+2))
      saved_nodel_hex = 0d0
      ! -=-=-=-=-=-=-=-=-
#endif

      hex_map_mid_b=0
      chi=0d0
      sefve=0d0
      ltolb = 0
      hex_map=0
      hex_map_mid=0
      iapoint=0
      iasurf=0
      ltox_hex=0
      ltoy_hex=0
      nodel_hex=0
      nxs=0
      nxe=0
      nys=0
      nye=0
      iasytype = 0
      nrnx_hex = 0
      layh = 0
      layp = 0
      lays = 0
      icorexs = 0
      icorexe = 0
      icorexps = 0
      icorexpe = 0
      icorexss = 0
      icorexse = 0
      icorexfs = 0
      icorexfe = 0
      wtass = 0d0
      wtasssum = 0d0
      izpr = 0
      iprcomp = 0
!! === xschi
      xschi(:,:,1)=1d0
      xschi(:,:,2)=0d0
      xschid(:,:,1)=1d0
      xschid(:,:,2)=0d0
!! === hex_map
!! I_LP_1N = hex_map_mid_b = hex_map_mid
      k=1
      do ix=1,nn
         iy1=1
         do iy=1,nn
            a=i_farf_1n_2d(iy, ix)
            b=hex_map_index(ix, iy)
            if (b == 1) then
                if (a==0) then
                    hex_map_mid_b(ix,iy1)= 0
                    hex_map_mid(ix,iy1)= 0
                    iy1=iy1+1
                else
                    n = i_lp_1n(k)
                    hex_map_mid_b(ix,iy1)= n
                    hex_map_mid(ix,iy1)= n
                    iy1=iy1+1
                endif
            endif
            if (a > 0) k = k + 1
         enddo
      enddo

      RETURN
      END SUBROUTINE hex_map_make
#endif


#ifdef tuan_tr_test
      SUBROUTINE hex_crbank
      !! FOR FAST REACTOR TRANSIENT CALCULATION 
      use Inc_TPEN
      use Inc_CR
      use Inc_Geometry

      IMPLICIT NONE
      INTEGER(4) :: i, ival, iy1, nhalfs, ir
      INTEGER(4) :: nxcol, nhalfcol, ix1, imodv
      INTEGER(4) :: iy, ix
      INTEGER(4) :: ii
      INTEGER(4) :: Ixy_1N_S, Ixy_1N_E
      INTEGER(4) :: Iy_1N
      INTEGER(4) :: istrod, icrb, ipend, icr, ncrbasy

     
      !! === ALLOCATE === !!
      if (.not.allocated(layh))    allocate(layh(-nxfc:nxfc,-nyfc:nyfc))
      if (.not.allocated(laybank)) allocate(laybank(-nxfc:nxfc,-nyfc:nyfc))
      !! === ALLOCATE === !!

      ival=0
      layh=ival
      iy1=0
      nxrow=Nx_1N
      nhalfs = nxrow/2
      iy1=-2*nhalfs
      ir=0
      ii = 1

      laybank = 0
      !! Full core model
      !! === control rod bank position === !!
      do Iy_1N = 1,Ny_1N
         IF (Ix_Start_y_1N(Iy_1N) == 0) cycle 
         Ixy_1N_S    = IxIy_1NToIxy_1N( Ix_Start_y_1N(Iy_1N), Iy_1N )
         Ixy_1N_E    = IxIy_1NToIxy_1N( Ix_End_y_1N  (Iy_1N), Iy_1N )
         nxcol = Ixy_1N_E - Ixy_1N_S + 1
         ir = ir+1
         nhalfcol=nxcol/2
         ix1=-4*nhalfcol
         imodv=mod(nxcol,2)
         if (imodv==0) ix1=ix1+2
         do i = Ixy_1N_S,Ixy_1N_E
            laybank(ix1,iy1) = I_CR_1N(i)
            ix1 = ix1+4
         enddo
         iy1 = iy1+2
      enddo

      do iy=icordys,-1
         do ix = icordxs, icordxe
            laybank(ix,iy)=0
         enddo
      enddo
      do ix=icordxs,-1
         laybank(ix,0)=0
      enddo
      !! === control rod bank position === !!

      !! === assign crbank === !!
      istrod = 0
      if (.not.allocated(lstroda))  allocate(lstroda(0:N_CR))
      if (.not.allocated(lstrodb))  allocate(lstrodb(0:N_CR))
      if (.not.allocated(lcrbptr))  allocate(lcrbptr(0:N_CR))
      if (.not.allocated(lcrbtola)) allocate(lcrbtola(nxy))
!      if (.not.allocated(laptr))    allocate(laptr(nxy))
      do iy=icordys,icordye,2
         do ix=icordxs,icordxe,2
            icrb=laybank(ix,iy)
            if(icrb<0)then
               istrod=istrod+1
               icrb=-icrb
               laybank(ix,iy)=icrb
               lstroda(istrod)=hex_map(ix,iy)
               lstrodb(-icrb)=.true.
            endif
         enddo
      enddo
      nstrod=istrod
      ipend=0
      lcrbptr(0)=0
      do icr=1,n_cr
         ncrbasy=0
         do iy=icordys,icordye,2
            do ix=icordxs,icordxe,2
               icrb=laybank(ix,iy)
               if(icr.eq.icrb) then
                  ncrbasy=ncrbasy+1
                  lcrbtola(ipend+ncrbasy)=hex_map(ix,iy)
               endif
            enddo
         enddo
         ipend=ipend+ncrbasy
         lcrbptr(icr)=ipend
      enddo
      !! === assing crbank === !!

      RETURN
      END SUBROUTINE hex_crbank

#endif

      END MODULE Mod_GetNode
