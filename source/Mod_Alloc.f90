
      MODULE Mod_Alloc

      USE Inc_Constant
      USE Inc_Utile

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif

      IMPLICIT NONE

      INTEGER :: N_SGL_DIGITS, N_DBL_DIGITS, N_INT_ORDER
      INTEGER :: NBFT, NBDT, NBIT
      PARAMETER (N_SGL_DIGITS=6)
      PARAMETER (N_DBL_DIGITS=15)
      PARAMETER (N_INT_ORDER=8)
      PARAMETER (NBFT=SELECTED_REAL_KIND(N_SGL_DIGITS))
      PARAMETER (NBDT=SELECTED_REAL_KIND(N_DBL_DIGITS))
      PARAMETER (NBIT=SELECTED_INT_KIND(N_INT_ORDER))

      INTERFACE Alloc
         MODULE PROCEDURE Alloc_D1
         MODULE PROCEDURE Alloc_I1
         MODULE PROCEDURE Alloc_C1
         MODULE PROCEDURE Alloc_L1
         MODULE PROCEDURE Alloc_D2
         MODULE PROCEDURE Alloc_I2
         MODULE PROCEDURE Alloc_C2
         MODULE PROCEDURE Alloc_L2
         MODULE PROCEDURE Alloc_D3
         MODULE PROCEDURE Alloc_I3
         MODULE PROCEDURE Alloc_C3
         MODULE PROCEDURE Alloc_L3
         MODULE PROCEDURE Alloc_D4
         MODULE PROCEDURE Alloc_I4
         MODULE PROCEDURE Alloc_C4
         MODULE PROCEDURE Alloc_L4
         MODULE PROCEDURE Alloc_D5
         MODULE PROCEDURE Alloc_I5
         MODULE PROCEDURE Alloc_C5
         MODULE PROCEDURE Alloc_L5
         MODULE PROCEDURE Alloc_D6
         MODULE PROCEDURE Alloc_I6
         MODULE PROCEDURE Alloc_C6
         MODULE PROCEDURE Alloc_L6
         MODULE PROCEDURE Alloc_D7
         MODULE PROCEDURE Alloc_I7
         MODULE PROCEDURE Alloc_C7
         MODULE PROCEDURE Alloc_L7
      END INTERFACE

      INTERFACE Alloc0
         MODULE PROCEDURE Alloc0_D1
         MODULE PROCEDURE Alloc0_I1
         MODULE PROCEDURE Alloc0_C1
         MODULE PROCEDURE Alloc0_L1
         MODULE PROCEDURE Alloc0_D2
         MODULE PROCEDURE Alloc0_I2
         MODULE PROCEDURE Alloc0_C2
         MODULE PROCEDURE Alloc0_L2
         MODULE PROCEDURE Alloc0_D3
         MODULE PROCEDURE Alloc0_I3
         MODULE PROCEDURE Alloc0_C3
         MODULE PROCEDURE Alloc0_L3
         MODULE PROCEDURE Alloc0_D4
         MODULE PROCEDURE Alloc0_I4
         MODULE PROCEDURE Alloc0_C4
         MODULE PROCEDURE Alloc0_L4
         MODULE PROCEDURE Alloc0_D5
         MODULE PROCEDURE Alloc0_I5
         MODULE PROCEDURE Alloc0_C5
         MODULE PROCEDURE Alloc0_L5
      END INTERFACE

      INTERFACE dmalloc
         MODULE PROCEDURE mallocf1
         MODULE PROCEDURE malloci1
         MODULE PROCEDURE mallocl1
         MODULE PROCEDURE mallocf2
         MODULE PROCEDURE malloci2
         MODULE PROCEDURE mallocl2
         MODULE PROCEDURE mallocf3
         MODULE PROCEDURE malloci3
         MODULE PROCEDURE mallocf4
         MODULE PROCEDURE malloci4
         MODULE PROCEDURE mallocf5
         MODULE PROCEDURE malloci5
         MODULE PROCEDURE mallocf7
         MODULE PROCEDURE malloci7
      END INTERFACE

      INTERFACE dmalloc0
         MODULE PROCEDURE mallocf01
         MODULE PROCEDURE malloci01
         MODULE PROCEDURE mallocl01
         MODULE PROCEDURE mallocf02
         MODULE PROCEDURE malloci02
         MODULE PROCEDURE mallocl02
         MODULE PROCEDURE mallocf03
         MODULE PROCEDURE malloci03
         MODULE PROCEDURE mallocf04
         MODULE PROCEDURE malloci04
      END INTERFACE

      CONTAINS

      SUBROUTINE Alloc_D1(a, n1)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: n1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_D1] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1))
      a = D0

      nbyte=nbyte+SIZE(a)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_I1(a, n1)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: n1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_I1] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_C1(a, n1)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: n1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_C1] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1))
      a = ''

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_L1(a, n1)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: n1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_L1] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_D2(a, n1, n2)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: a(:, :)
      INTEGER, INTENT(IN) :: n1, n2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_D2] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2))
      a = D0

      nbyte=nbyte+SIZE(a)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_I2(a, n1, n2)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:, :)
      INTEGER, INTENT(IN) :: n1, n2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_I2] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_C2(a, n1, n2)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:, :)
      INTEGER, INTENT(IN) :: n1, n2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_C2] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2))
      a = ''

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_L2(a, n1, n2)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:, :)
      INTEGER, INTENT(IN) :: n1, n2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_L2] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_D3(a, n1, n2, n3)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_D3] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3))
      a = D0

      nbyte=nbyte+SIZE(a)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_I3(a, n1, n2, n3)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_I3] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_C3(a, n1, n2, n3)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_C3] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3))
      a = ''

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_L3(a, n1, n2, n3)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_L3] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_D4(a, n1, n2, n3, n4)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_D4] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4))
      a = D0

      nbyte=nbyte+SIZE(a)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_I4(a, n1, n2, n3, n4)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_I4] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_C4(a, n1, n2, n3, n4)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_C4] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4))
      a = ''

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_L4(a, n1, n2, n3, n4)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_L4] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_D5(a, n1, n2, n3, n4, n5)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_D5] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5))
      a = D0

      nbyte=nbyte+SIZE(a)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_I5(a, n1, n2, n3, n4, n5)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_I5] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_C5(a, n1, n2, n3, n4, n5)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_C5] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5))
      a = ''

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_L5(a, n1, n2, n3, n4, n5)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_L5] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_D6(a, n1, n2, n3, n4, n5, n6)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5, n6

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_D6] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5, n6))
      a = D0

      nbyte=nbyte+SIZE(a)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_I6(a, n1, n2, n3, n4, n5, n6)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5, n6

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_I6] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5, n6))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_C6(a, n1, n2, n3, n4, n5, n6)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5, n6

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_C6] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5, n6))
      a = ''

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_L6(a, n1, n2, n3, n4, n5, n6)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5, n6

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_L6] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5, n6))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_D7(a, n1, n2, n3, n4, n5, n6, n7)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5, n6, n7

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_D7] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5, n6, n7))
      a = D0

      nbyte=nbyte+SIZE(a)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_I7(a, n1, n2, n3, n4, n5, n6, n7)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5, n6, n7

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_I7] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5, n6, n7))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_C7(a, n1, n2, n3, n4, n5, n6, n7)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5, n6, n7

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_C7] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5, n6, n7))
      a = ''

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc_L7(a, n1, n2, n3, n4, n5, n6, n7)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :, :, :)
      INTEGER, INTENT(IN) :: n1, n2, n3, n4, n5, n6, n7

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc_L7] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(n1, n2, n3, n4, n5, n6, n7))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_D1(a, nb1, ne1)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: nb1, ne1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_D1] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1))
      a = D0

      nbyte=nbyte+SIZE(a)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_I1(a, nb1, ne1)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: nb1, ne1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_I1] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_C1(a, nb1, ne1)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: nb1, ne1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_C1] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      ALLOCATE(a(nb1:ne1))
      a = ''

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_L1(a, nb1, ne1)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: nb1, ne1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_L1] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_D2(a, nb1, ne1, nb2, ne2)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: a(:, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_D2] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2))
      a = D0

      nbyte=nbyte+SIZE(a)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_I2(a, nb1, ne1, nb2, ne2)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_I2] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_C2(a, nb1, ne1, nb2, ne2)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_C2] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2))
      a = ''

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_L2(a, nb1, ne1, nb2, ne2)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_L2] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_D3(aa, nb1, ne1, nb2, ne2, nb3, ne3)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: aa(:, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_D3] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      ALLOCATE(aa(nb1:ne1, nb2:ne2, nb3:ne3))
      aa = D0

      nbyte=nbyte+SIZE(aa)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_I3(a, nb1, ne1, nb2, ne2, nb3, ne3)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_I3] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2, nb3:ne3))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_C3(a, nb1, ne1, nb2, ne2, nb3, ne3)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_C3] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2, nb3:ne3))
      a = ''

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_L3(a, nb1, ne1, nb2, ne2, nb3, ne3)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_L3] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2, nb3:ne3))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_D4(a, nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_D4] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2, nb3:ne3, nb4:ne4))
      a = D0

      nbyte=nbyte+SIZE(a)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_I4(a, nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_I4] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2, nb3:ne3, nb4:ne4))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_C4(a, nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_C4] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2, nb3:ne3, nb4:ne4))
      a = ''

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_L4(a, nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_L4] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2, nb3:ne3, nb4:ne4))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_D5(a, nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4, nb5, ne5)
      REAL(8), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4, nb5, ne5

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_D5] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2, nb3:ne3, nb4:ne4, nb5:ne5))
      a = D0

      nbyte=nbyte+SIZE(a)*8d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_I5(a, nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4, nb5, ne5)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4, nb5, ne5

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_I5] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2, nb3:ne3, nb4:ne4, nb5:ne5))
      a = 0

      nbyte=nbyte+SIZE(a)*4d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_C5(a, nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4, nb5, ne5)
      CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4, nb5, ne5

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_C5] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2, nb3:ne3, nb4:ne4, nb5:ne5))
      a = ''

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE


      SUBROUTINE Alloc0_L5(a, nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4, nb5, ne5)
      LOGICAL(1), ALLOCATABLE, INTENT(INOUT) :: a(:, :, :, :, :)
      INTEGER, INTENT(IN) :: nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4, nb5, ne5

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Alloc0_L5] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(Allocated(a)) deallocate(a)
      ALLOCATE(a(nb1:ne1, nb2:ne2, nb3:ne3, nb4:ne4, nb5:ne5))
      a = .FALSE.

      nbyte=nbyte+SIZE(a)*1d0
      RETURN
      END SUBROUTINE

      SUBROUTINE mallocf1(a,n1)
            IMPLICIT NONE
            INTEGER :: n1
            REAL(8),POINTER :: a(:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocf1] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1))
            a=zero
            nbyte=nbyte+SIZE(a)*8
      END SUBROUTINE mallocf1

      SUBROUTINE malloci1(a,n1)
            IMPLICIT NONE
            INTEGER :: n1
            INTEGER,POINTER :: a(:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [malloci1] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1))
            a=0
            nbyte=nbyte+SIZE(a)*4
      END SUBROUTINE malloci1
!
      SUBROUTINE mallocl1(a,n1)
            IMPLICIT NONE
            INTEGER :: n1
            LOGICAL,POINTER :: a(:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocl1] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1))
            a=.FALSE.
            nbyte=nbyte+SIZE(a)
      END SUBROUTINE mallocl1
!
      SUBROUTINE mallocf2(a,n1,n2)
            IMPLICIT NONE
            INTEGER :: n1, n2
            REAL(8),POINTER :: a(:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocf2] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1,n2))
            a=zero
            nbyte=nbyte+SIZE(a)*8
      END SUBROUTINE mallocf2
!
!
      SUBROUTINE malloci2(a,n1,n2)
            IMPLICIT NONE
            INTEGER :: n1, n2
            INTEGER,POINTER :: a(:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [malloci2] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1,n2))
            a=0
            nbyte=nbyte+SIZE(a)*4
      END SUBROUTINE malloci2
!
      SUBROUTINE mallocl2(a,n1,n2)
            IMPLICIT NONE
            INTEGER :: n1, n2
            LOGICAL,POINTER :: a(:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocl2] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1,n2))
            a=.FALSE.
            nbyte=nbyte+SIZE(a)
      END SUBROUTINE mallocl2
!
      SUBROUTINE mallocf3(a,n1,n2,n3)
            IMPLICIT NONE
            INTEGER :: n1, n2, n3
            REAL(8),POINTER :: a(:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocf3] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1,n2,n3))
            a=zero
            nbyte=nbyte+SIZE(a)*8
      END SUBROUTINE mallocf3
!
!
      SUBROUTINE malloci3(a,n1,n2,n3)
            IMPLICIT NONE
            INTEGER :: n1, n2, n3
            INTEGER,POINTER :: a(:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [malloci3] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1,n2,n3))
            a=0
            nbyte=nbyte+SIZE(a)*4
      END SUBROUTINE malloci3
!
      SUBROUTINE mallocf4(a,n1,n2,n3,n4)
            IMPLICIT NONE
            INTEGER :: n1, n2, n3, n4
            REAL(8),POINTER :: a(:,:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocf4] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1,n2,n3,n4))
            a=zero
            nbyte=nbyte+SIZE(a)*8
      END SUBROUTINE mallocf4
!
!
      SUBROUTINE malloci4(a,n1,n2,n3,n4)
            IMPLICIT NONE
            INTEGER :: n1, n2, n3, n4
            INTEGER,POINTER :: a(:,:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [malloci4] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1,n2,n3,n4))
            a=0
            nbyte=nbyte+SIZE(a)*4
      END SUBROUTINE malloci4
!
      SUBROUTINE mallocf5(a,n1,n2,n3,n4,n5)
            IMPLICIT NONE
            INTEGER :: n1, n2, n3, n4, n5
            REAL(8),POINTER :: a(:,:,:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocf5] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1,n2,n3,n4,n5))
            a=zero
            nbyte=nbyte+SIZE(a)*8
      END SUBROUTINE mallocf5
!
!
      SUBROUTINE malloci5(a,n1,n2,n3,n4,n5)
            IMPLICIT NONE
            INTEGER :: n1, n2, n3, n4, n5
            INTEGER,POINTER :: a(:,:,:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [malloci5] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1,n2,n3,n4,n5))
            a=0
            nbyte=nbyte+SIZE(a)*4
      END SUBROUTINE malloci5
!
      SUBROUTINE mallocf7(a,n1,n2,n3,n4,n5,n6,n7)
            IMPLICIT NONE
            INTEGER :: n1, n2, n3, n4, n5, n6, n7
            REAL(8),POINTER :: a(:,:,:,:,:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocf7] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1,n2,n3,n4,n5,n6,n7))
            a=zero
            nbyte=nbyte+SIZE(a)*8
      END SUBROUTINE mallocf7
!
!
      SUBROUTINE malloci7(a,n1,n2,n3,n4,n5,n6,n7)
            IMPLICIT NONE
            INTEGER :: n1, n2, n3, n4, n5, n6, n7
            INTEGER,POINTER :: a(:,:,:,:,:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [malloci7] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(n1,n2,n3,n4,n5,n6,n7))
            a=0
            nbyte=nbyte+SIZE(a)*4
      END SUBROUTINE malloci7
!
      SUBROUTINE mallocf01(a,nb1,ne1)
            IMPLICIT NONE
            INTEGER :: nb1, ne1
            REAL(8),POINTER :: a(:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocf01] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(nb1:ne1))
            a=zero
            nbyte=nbyte+SIZE(a)*8
      END SUBROUTINE mallocf01
!
!
      SUBROUTINE malloci01(a,nb1,ne1)
            IMPLICIT NONE
            INTEGER :: nb1, ne1
            INTEGER,POINTER :: a(:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [malloci01] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(nb1:ne1))
            a=0
            nbyte=nbyte+SIZE(a)*4
      END SUBROUTINE malloci01
!
      SUBROUTINE mallocl01(a,nb1,ne1)
            IMPLICIT NONE
            INTEGER :: nb1, ne1
            LOGICAL,POINTER :: a(:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocl01] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(nb1:ne1))
            a=.FALSE.
            nbyte=nbyte+SIZE(a)
      END SUBROUTINE mallocl01
!
      SUBROUTINE mallocf02(a,nb1,ne1,nb2,ne2)
            IMPLICIT NONE
            INTEGER :: nb1, ne1, nb2, ne2
            REAL(8),POINTER :: a(:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocf02] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(nb1:ne1,nb2:ne2))
            a=zero
            nbyte=nbyte+SIZE(a)*8
      END SUBROUTINE mallocf02
!
!
      SUBROUTINE malloci02(a,nb1,ne1,nb2,ne2)
            IMPLICIT NONE
            INTEGER :: nb1, ne1, nb2, ne2
            INTEGER,POINTER :: a(:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [malloci02] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(nb1:ne1,nb2:ne2))
            a=0
            nbyte=nbyte+SIZE(a)*4
      END SUBROUTINE malloci02
!
      SUBROUTINE mallocl02(a,nb1,ne1,nb2,ne2)
            IMPLICIT NONE
            INTEGER :: nb1, ne1, nb2, ne2
            LOGICAL,POINTER :: a(:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocl02] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(nb1:ne1,nb2:ne2))
            a=.FALSE.
            nbyte=nbyte+SIZE(a)
      END SUBROUTINE mallocl02
!
      SUBROUTINE mallocf03(a,nb1,ne1,nb2,ne2,nb3,ne3)
            IMPLICIT NONE
            INTEGER :: nb1, ne1, nb2, ne2, nb3, ne3
            REAL(8),POINTER :: a(:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocf03] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(nb1:ne1,nb2:ne2,nb3:ne3))
            a=zero
            nbyte=nbyte+SIZE(a)*8
      END SUBROUTINE mallocf03
!
!
      SUBROUTINE malloci03(a,nb1,ne1,nb2,ne2,nb3,ne3)
            IMPLICIT NONE
            INTEGER :: nb1, ne1, nb2, ne2, nb3, ne3
            INTEGER,POINTER :: a(:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [malloci03] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(nb1:ne1,nb2:ne2,nb3:ne3))
            a=0
            nbyte=nbyte+SIZE(a)*4
      END SUBROUTINE malloci03
!
      SUBROUTINE mallocf04(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4)
            IMPLICIT NONE
            INTEGER :: nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4
            REAL(8),POINTER :: a(:,:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mallocf04] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4))
            a=zero
            nbyte=nbyte+SIZE(a)*8
      END SUBROUTINE mallocf04
!
!
      SUBROUTINE malloci04(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4)
            IMPLICIT NONE
            INTEGER :: nb1, ne1, nb2, ne2, nb3, ne3, nb4, ne4
            INTEGER,POINTER :: a(:,:,:,:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [malloci04] in Mod_Alloc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

            ALLOCATE(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4))
            a=0
            nbyte=nbyte+SIZE(a)*4
      END SUBROUTINE malloci04

      END MODULE Mod_Alloc

