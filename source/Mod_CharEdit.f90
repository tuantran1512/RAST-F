
      MODULE Mod_CharEdit

      USE Inc_Constant
      USE Inc_Geometry
      USE Inc_RP
      use inc_parallel, only: comm

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif


      IMPLICIT NONE

      interface print_msg
         module procedure print_msg_a
         module procedure print_msg_aa
         module procedure print_msg_aaa
         module procedure print_msg_afa
         module procedure print_msg_afaf
         module procedure print_msg_afafa
         module procedure print_msg_afafaf
         module procedure print_msg_afafafa
         module procedure print_msg_aiai
         module procedure print_msg_aiaia
         module procedure print_msg_ae
         module procedure print_msg_aee
         module procedure print_msg_ai
         module procedure print_msg_aii
         module procedure print_msg_aiaf
         module procedure print_msg_iae
         module procedure print_msg_iff
      end interface print_msg

      CONTAINS

      SUBROUTINE CharToUpper(I_Char)

      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: I_Char
      INTEGER, PARAMETER :: INDXA = 97, IDNXZ = 122
      INTEGER :: Len_Char, ASCII, i

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [CharToUpper] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      Len_Char = LEN_TRIM(I_Char)
      DO i = 1, Len_Char
         IF ( I_Char(i:i) == ' ' ) THEN
            CYCLE
         END IF
         ASCII = ICHAR( I_Char(i:i) )
         IF ( ( ASCII >= INDXA ) .AND. ( ASCII <= IDNXZ ) ) THEN
            I_Char(i:i) = CHAR( ASCII - 32 )
         END IF
      END DO

      RETURN
      END SUBROUTINE CharToUpper


      SUBROUTINE RemoveSubCard(I_Char, I_SubCard)

      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: I_Char
      CHARACTER(*), INTENT(IN) :: I_SubCard
      INTEGER :: Len_Char, Len_SubCard, i
      CHARACTER(1024) :: Buff_I_Char

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [RemoveSubCard] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      Buff_I_Char = I_Char

      CALL CharToUpper(Buff_I_Char)

      Len_Char    = LEN_TRIM(I_Char)
      Len_SubCard = LEN_TRIM(I_SubCard)

      DO i = 1, Len_Char
         IF ( Buff_I_Char( i : i + Len_SubCard - 1 ) == I_SubCard ) THEN
            I_Char( 1 : LEN(I_Char) - i - Len_SubCard + 1 ) = I_Char( i + Len_SubCard : LEN(I_Char) )
            EXIT
         END IF
      END DO

      Len_Char = LEN_TRIM(I_Char)
      DO i = 1, Len_Char
         IF ( I_Char(i:i) /= ' ' ) THEN
            I_Char( 1 : Len_Char - i + 1 ) = I_Char( i : Len_Char )
            I_Char( Len_Char - i + 2 : Len_Char ) = ''
            EXIT
         END IF
      END DO


      Len_Char = LEN_TRIM(I_Char)
      DO i = 1, Len_Char
         IF ( I_Char(i:i) == '!' ) THEN
            I_Char(i:Len_Char) = ''
            EXIT
         END IF
      END DO

      RETURN
      END SUBROUTINE RemoveSubCard


      SUBROUTINE Count_Field( OneLine, N_Field )

      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: OneLine
      INTEGER, INTENT(INOUT) :: N_Field
      INTEGER :: Len_Char, i, j, k
      LOGICAL(1) :: Flag_LoopExit

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Count_Field] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      Len_Char = LEN_TRIM(OneLine)

      IF ( Len_Char == 0 ) THEN
         RETURN
      END IF

      i = 0
      DO WHILE (.TRUE.)
         i = i + 1
         IF ( OneLine(i:i) == ',' ) THEN
            OneLine( i + 2 : Len_Char + 1 ) = OneLine( i + 1 : Len_Char )
            OneLine( i + 1 : i + 1 ) = ' '
            Len_Char = Len_Char + 1
            i = i + 1
         END IF
         IF ( i + 1 >= Len_Char ) THEN
            EXIT
         END IF
      END DO

      i = 0
      RemBlank: DO WHILE (.TRUE.)
         Flag_LoopExit = .TRUE.
         DO WHILE (.TRUE.)
            i = i + 1
            IF ( OneLine( i : i + 1 ) == '  ' ) THEN
               OneLine( i : Len_Char - 1 ) = OneLine( i + 1 : Len_Char )
               OneLine( Len_Char : Len_Char ) = ' '
               Len_Char = Len_Char - 1
               i = i + 1
               Flag_LoopExit = .FALSE.
            END IF
            IF ( i + 1 >= Len_Char ) THEN
               EXIT
            END IF
         END DO
         Len_Char = LEN_TRIM(OneLine)
         i = 0
         IF ( Flag_LoopExit ) THEN
            EXIT RemBlank
         END IF
      END DO RemBlank

      i = 0
      DO WHILE (.TRUE.)
         i = i + 1
         IF ( OneLine( i : i + 1 ) == ' ,' ) THEN
            IF ( i /= 1 ) THEN
               IF ( OneLine( i - 1 : i - 1 ) == ',' ) THEN
                  i = i + 1
                  CYCLE
               END IF
            END IF
            OneLine( i : Len_Char - 1 ) = OneLine( i + 1 : Len_Char )
            OneLine( Len_Char : Len_Char ) = ' '
            Len_Char = Len_Char - 1
            i = i + 1
            Flag_LoopExit = .FALSE.
         END IF
         IF ( i + 1 >= Len_Char ) THEN
            EXIT
         END IF
      END DO

      Len_Char = LEN_TRIM(OneLine)
      N_Field = 0
      IF ( ( OneLine(1:1) /= ',' ) .AND. ( OneLine(1:1) /= ' ' ) ) THEN
         N_Field = 1
      END IF
      DO i = 1, Len_Char
         IF ( OneLine(i:i) == ' ' ) THEN
            N_Field = N_Field + 1
         END IF
      END DO
      IF ( OneLine(1:1) == ',' ) THEN
         N_Field = N_Field + 1
         OneLine( 3 : Len_Char + 2 ) = OneLine( 1 : Len_Char )
         OneLine( 1 : 2 ) = '-1'
         Len_Char = Len_Char + 2
      END IF

      IF ( Len_Char > 1 ) THEN
         IF ( OneLine( Len_Char : Len_Char ) == ',' ) THEN
            N_Field = N_Field + 1
            OneLine( Len_Char + 1 : Len_Char + 3 ) = ' -1'
            Len_Char = Len_Char + 3
         END IF
      END IF

      i = 0
      DO WHILE (.TRUE.)
         i = i + 1
         IF ( OneLine( i : i + 2 ) == ', ,' ) THEN
            OneLine( i + 4 : Len_Char + 2 ) = OneLine( i + 2 : Len_Char )
            OneLine( i + 2 : i + 3 ) = '-1'
            Len_Char = Len_Char + 2
         END IF
         IF ( i + 2 >= Len_Char ) THEN
            EXIT
         END IF
      END DO

      DO i = 1, Len_Char
         IF ( OneLine( i : i ) == '*' ) THEN
            j = i
            DO WHILE (.TRUE.)
               j = j - 1
               IF ( j == 0 ) THEN
                  EXIT
               END IF
               IF ( OneLine( j : j ) == ' ' ) THEN
                  EXIT
               END IF
            END DO
            READ( OneLine( j + 1 : i - 1 ), * ) k
            N_Field = N_Field + k - 1
         END IF
      END DO

      RETURN
      END SUBROUTINE Count_Field


      SUBROUTINE RemoveGivenChar(I_Char, I_GivenChar)

      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: I_Char
      CHARACTER(*), INTENT(IN) :: I_GivenChar
      INTEGER :: Len_Char, Len_GivenChar, i

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [RemoveGivenChar] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      Len_Char      = LEN_TRIM(I_Char)
      Len_GivenChar = LEN_TRIM(I_GivenChar)

      DO i = 1, ( Len_Char - Len_GivenChar + 1 )
         IF ( I_Char( i : i + Len_GivenChar - 1 ) == I_GivenChar ) THEN
            I_Char( i : i + Len_GivenChar - 1 ) = ''
         END IF
      END DO

      RETURN
      END SUBROUTINE RemoveGivenChar

      ! -----------------------------------------------------------------------
      subroutine print_msg_a(i_tag,str1)
      implicit none
      character(*), intent(in) :: str1
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''
      if (comm%if_master) then

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_a] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1
      endif
100   format(a,a)
      return
      end subroutine print_msg_a

      subroutine print_msg_aa(i_tag,str1,str2)
      implicit none
      character(*), intent(in) :: str1
      character(*), intent(in) :: str2
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_aa] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,str2
      endif
100   format(a,a,a)
      return
      end subroutine print_msg_aa

      subroutine print_msg_aaa(i_tag,str1,str2,str3)
      implicit none
      character(*), intent(in) :: str1
      character(*), intent(in) :: str2
      character(*), intent(in) :: str3
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_aaa] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,str2,str3
      endif
100   format(a,a,a,a)
      return
      end subroutine print_msg_aaa

      subroutine print_msg_afa(i_tag,str1,real1,str2)
      implicit none
      character(*), intent(in) :: str1
      real(8), intent(in) :: real1
      character(*), intent(in) :: str2
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_afa] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         if (real1<1000d0) then
            write(*,100) trim(tag)//' ',str1,real1,str2
         else
            write(*,101) trim(tag)//' ',str1,real1,str2
         endif
      endif
100   format(a,a,f8.3,a)
101   format(a,a,f10.2,a)
      return
      end subroutine print_msg_afa

      subroutine print_msg_afaf(i_tag,str1,real1,str2,real2)
      implicit none
      character(*), intent(in) :: str1
      real(8), intent(in) :: real1
      character(*), intent(in) :: str2
      real(8), intent(in) :: real2
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_afaf] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,real1,str2,real2
      endif
100   format(a,a,f8.3,a,f8.3)
      return
      end subroutine print_msg_afaf

      subroutine print_msg_afafa(i_tag,str1,real1,str2,real2,str3)
      implicit none
      character(*), intent(in) :: str1
      real(8), intent(in) :: real1
      character(*), intent(in) :: str2
      real(8), intent(in) :: real2
      character(*), intent(in) :: str3
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_afafa] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,real1,str2,real2,str3
      endif
100   format(a,a,f8.3,a,f8.3,a)
      return
      end subroutine print_msg_afafa

      subroutine print_msg_afafaf(i_tag,str1,real1,str2,real2,str3,real3)
      implicit none
      character(*), intent(in) :: str1
      real(8), intent(in) :: real1
      character(*), intent(in) :: str2
      real(8), intent(in) :: real2
      character(*), intent(in) :: str3
      real(8), intent(in) :: real3
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_afafaf] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,real1,str2,real2,str3,real3
      endif
100   format(a,a,f8.3,a,f8.3,a,f8.3)
      return
      end subroutine print_msg_afafaf

      subroutine print_msg_afafafa(i_tag,str1,real1,str2,real2,str3,real3,str4)
      implicit none
      character(*), intent(in) :: str1
      real(8), intent(in) :: real1
      character(*), intent(in) :: str2
      real(8), intent(in) :: real2
      character(*), intent(in) :: str3
      real(8), intent(in) :: real3
      character(*), intent(in) :: str4
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_afafafa] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,real1,str2,real2,str3,real3,str4
      endif
100   format(a,a,f10.6,a,f10.3,a,f10.4,a)
      return
      end subroutine print_msg_afafafa

      subroutine print_msg_aiai(i_tag,str1,int1,str2,int2)
      implicit none
      character(*), intent(in) :: str1
      integer, intent(in) :: int1
      character(*), intent(in) :: str2
      integer, intent(in) :: int2
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_aiai] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         if (int2<10000) then
            write(*,100) trim(tag)//' ',str1,int1,str2,int2
         else
            write(*,200) trim(tag)//' ',str1,int1,str2,int2
         endif
      endif
100   format(a,a,i5,a,i5)
200   format(a,a,i5,a,i9)
      return
      end subroutine print_msg_aiai

      subroutine print_msg_aiaia(i_tag,str1,int1,str2,int2,str3)
      implicit none
      character(*), intent(in) :: str1
      integer, intent(in) :: int1
      character(*), intent(in) :: str2
      integer, intent(in) :: int2
      character(*), intent(in) :: str3
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_aiaia] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,int1,str2,int2,str3
      endif
100   format(a,a,i5,a,i5,a)
      return
      end subroutine print_msg_aiaia

      subroutine print_msg_ae(i_tag,str1,real1)
      implicit none
      character(*), intent(in) :: str1
      real(8), intent(in) :: real1
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_ae] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,real1
      endif
100   format(a,a,es13.5)
      return
      end subroutine print_msg_ae

      subroutine print_msg_aee(i_tag,str1,real1,real2)
      implicit none
      character(*), intent(in) :: str1
      real(8), intent(in) :: real1
      real(8), intent(in) :: real2
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_aee] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,real1,real2
      endif
100   format(a,a,es13.5,es13.5)
      return
      end subroutine print_msg_aee

      subroutine print_msg_ai(i_tag,str1,int1)
      implicit none
      character(*), intent(in) :: str1
      integer, intent(in) :: int1
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_ai] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,int1
      endif
100   format(a,a," ",i0)
      return
      end subroutine print_msg_ai

      subroutine print_msg_aii(i_tag,str1,int1,int2)
      implicit none
      character(*), intent(in) :: str1
      integer, intent(in) :: int1
      integer, intent(in) :: int2
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_aii] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,int1,int2
      endif
100   format(a,a," ",i0," ",i0)
      return
      end subroutine print_msg_aii

      subroutine print_msg_aiaf(i_tag,str1,int1,str2,real1)
      implicit none
      character(*), intent(in) :: str1
      integer, intent(in) :: int1
      character(*), intent(in) :: str2
      real(8), intent(in) :: real1
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_aiaf] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',str1,int1,str2,real1
      endif
100   format(a,a," ",i0," ",a," ",f10.2)
      return
      end subroutine print_msg_aiaf

      subroutine print_msg_iae(i_tag,int1,str1,real1)
      implicit none
      integer, intent(in) :: int1
      character(*), intent(in) :: str1
      real(8), intent(in) :: real1
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_iae] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',int1,str1,real1
      endif
100   format(a,i4,a5,es13.5)
      return
      end subroutine print_msg_iae

      subroutine print_msg_iff(i_tag,int1,real1,real2)
      implicit none
      integer, intent(in) :: int1
      real(8), intent(in) :: real1
      real(8), intent(in) :: real2
      integer, intent(in) :: i_tag
      character(len=20) :: tag=''


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_msg_iff] in Mod_CharEdit'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (comm%if_master) then
         select case (i_tag)
         case (0) ! Note
            tag=''
         case (1) ! Note
            tag='*** Note:'
         case (2) ! Warning
            tag='*** Warning:'
         case (3) ! Error
            tag='*** Error:'
         case default
            write(*,'(a,i0)') 'No tag in PRINT_MSG'
            stop
         end select
         write(*,100) trim(tag)//' ',int1,real1,real2
      endif
100   format(a,i10,f10.6,f10.3)
      return
      end subroutine print_msg_iff
      ! -----------------------------------------------------------------------

      END MODULE Mod_CharEdit