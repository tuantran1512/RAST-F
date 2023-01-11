#ifdef siarhei_delete


      MODULE Mod_Operator


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

      CONTAINS

      RECURSIVE FUNCTION Get_Factorial(N) RESULT(fac)

      IMPLICIT NONE

      REAL(8), INTENT(IN) :: N
      REAL(8) :: fac

      !IF ( abs(N-0d0)<1d-10 ) THEN

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_Factorial] in Mod_Operator'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF ( N==0d0 ) THEN
         fac = 1.0D0
      ELSE
         fac = N*Get_Factorial(N - 1.0D0)
      END IF

      RETURN
      END FUNCTION Get_Factorial


      FUNCTION Get_DetA(A) RESULT(D)

      IMPLICIT NONE

      REAL(8), DIMENSION(:, :), INTENT(IN) :: A
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: X
      REAL(8) :: D
      REAL(8) :: Quotient
      INTEGER :: N
      INTEGER :: i, j, m


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_DetA] in Mod_Operator'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      i = 0
      j = 0
      m = 0
      D = 1.0D0
      N = SIZE(A, 1)

      ALLOCATE(X(N, N))
      X = 0.0D0
      X = A
      DO i = 1, N
         D = D*X(i, i)  ! It Use Transformation Properties of Determinant
                        ! while Elementary Row Operation.
                        ! This function make matrix "A" to Upper Triangular Matrix
                        ! which has Reading 1 in Diagonal Term.
                        ! Therefore, its Diagonal Product is "1" or "0" and it's Determinat of "A"
         IF (i == N) EXIT  ! Gauss Elimination
         IF ( X(i,i)==0d0 ) THEN
            D = 0.0D0
            WRITE(*, '(A)') "*** Input Matrix is Singular Matrix!!"
            DEALLOCATE(X)
            STOP  ! If any Diagonal Term is "0", then det(A) = 0
         END IF

         Quotient = X(i, i)
         DO j = i, N
            X(i, j) = X(i, j)/Quotient  ! (det => det*A(i, i))
         END DO

         DO m = (i + 1), N
            Quotient = X(m, i)
            DO j = i, N
               X(m, j) = X(m, j) - Quotient*X(i, j)  ! (det => det)
            END DO
         END DO
      END DO

      DEALLOCATE(X)

      RETURN
      END FUNCTION Get_DetA

      FUNCTION Get_InvA(A) RESULT(invA)

      IMPLICIT NONE

      REAL(8), DIMENSION(:, :), INTENT(IN) :: A
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: invA, X
      REAL(8) :: detA, Quotient
      INTEGER :: N
      INTEGER :: i, j, m


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_InvA] in Mod_Operator'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      i    = 0
      j    = 0
      m    = 0

      N = SIZE(A, 1)

      ALLOCATE(invA(N, N))
      ALLOCATE(X(N, 2*N))  ! X = [A | eye]
      invA = 0.0D0
      X    = 0.0D0

      detA = Get_DetA(A)

      IF ( detA==0d0 ) THEN
      !IF ( abs(detA-0d0)<1d-10 ) THEN
         WRITE(*, '(A)') "*** Input Matrix is Singular Matrix!!"
         DEALLOCATE(X)
         STOP
      END IF

      DO i = 1, N
         DO j = 1, N
            X(i, j) = A(i, j)
         END DO
      END DO
      DO i = 1, N
         X(i, (i + N)) = 1.0D0
      END DO

      DO i = 1, N
         Quotient = X(i, i)
         DO j = i, 2*N
            X(i, j) = X(i, j)/Quotient
         END DO
         IF (i == N) EXIT
         DO m = (i + 1), N
            Quotient = X(m, i)
            DO j = i, 2*N
               X(m, j) = X(m, j) - Quotient*X(i, j)
            END DO
         END DO
      END DO

      ! This Procedure Makes Matrix X_Left => Reduced Echelon Form,
      ! Then X_Right(Identity Matrix) => inv(A)
      DO i = N, 2, - 1
         DO m = (i - 1), 1, - 1
            Quotient = X(m, i)
            DO j = 1, 2*N
               X(m, j) = X(m, j) - Quotient*X(i, j)
            END DO
         END DO
      END DO

      DO i = 1, N
         DO j = 1, N
            invA(i, j) = X(i, (j + N))
         END DO
      END DO

      DEALLOCATE(X)

      RETURN
      END FUNCTION Get_InvA


      FUNCTION Get_InvD(D) RESULT(invD)

      IMPLICIT NONE

      REAL(8), DIMENSION(:, :), INTENT(IN) :: D
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: invD
      REAL(8) :: detD
      INTEGER :: N
      INTEGER :: i


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_InvD] in Mod_Operator'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      i    = 0
      N = SIZE(D, 1)
      ALLOCATE(invD(N, N))
      invD = 0.0D0
      detD = Get_DetA(D)

      !IF ( abs(detD-0d0)<1d-10 ) THEN
      IF ( detD==0d0 ) THEN
         WRITE(*, '(A)') "*** Input Matrix is Singular Matrix!!"
         STOP
      END IF
      DO i = 1, N
         invD(i, i) = 1.0D0/D(i, i)
      END DO

      RETURN
      END FUNCTION Get_InvD


      FUNCTION Get_InvC(C) RESULT(invC)

      IMPLICIT NONE

      COMPLEX(8), DIMENSION(:, :), INTENT(IN) :: C
      COMPLEX(8), DIMENSION(:, :), ALLOCATABLE :: invC, X
      COMPLEX(8) :: Quotient
      INTEGER :: N
      INTEGER :: i, j, m


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Get_InvC] in Mod_Operator'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      i    = 0
      j    = 0
      m    = 0
      N = SIZE(REAL(C, 8), 1)
      ALLOCATE(invC(N, N))
      ALLOCATE(X(N, 2*N))  ! X = [C | eye]
      invC = (0.0D0, 0.0D0)
      X    = (0.0D0, 0.0D0)

      DO i = 1, N
         DO j = 1, N
            X(i, j) = C(i, j)
         END DO
      END DO
      DO i = 1, N
         X(i, (i + N)) = (1.0D0, 0.0D0)
      END DO
      DO i = 1, N
         Quotient = X(i, i)
         DO j = i, 2*N
            X(i, j) = X(i, j)/Quotient
         END DO
         IF (i == N) EXIT
         DO m = (i + 1), N
            Quotient = X(m, i)
            DO j = i, 2*N
               X(m, j) = X(m, j) - Quotient*X(i, j)
            END DO
         END DO
      END DO

      ! This Procedure Makes Matrix X_Left => Reduced Echelon Form,
      ! Then X_Right(Identity Matrix) => inv(C)
      DO i = N, 2, - 1
         DO m = (i - 1), 1, - 1
            Quotient = X(m, i)
            DO j = 1, 2*N
               X(m, j) = X(m, j) - Quotient*X(i, j)
            END DO
         END DO
      END DO
      DO i = 1, N
         DO j = 1, N
            invC(i, j) = X(i, (j + N))
         END DO
      END DO

      DEALLOCATE(X)

      RETURN
      END FUNCTION Get_InvC

      END MODULE Mod_Operator

#endif
