#ifdef siarhei_delete

      module mod_utility


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

      interface
      integer function c_chdir(path) bind(C,name="chdir")
      use iso_c_binding
      character(kind=c_char) :: path(*)
      end function c_chdir
      end interface

      interface c2c
      module procedure real_to_char
      module procedure int4_to_char
      module procedure int8_to_char
      end interface c2c

      interface get_digit
      module procedure int4_get_digit
      end interface get_digit

      contains

      integer function newunit(unit)
         integer, intent(out), optional :: unit
         integer, parameter :: LUN_MIN=10, LUN_MAX=1000
         logical :: opened
         integer :: lun
         newunit=-1
         do lun=LUN_MIN,LUN_MAX
           inquire(unit=lun,opened=opened)
           if (.not. opened) then
             newunit=lun
             exit
           end if
         end do
         if (present(unit)) unit=newunit
      end function newunit

      function replace_text (s,text,rep)  result(outs)

         character(*)        :: s,text,rep
         character(len(s)+100) :: outs     ! provide outs with extra 100 char len
         integer             :: i, nt, nr

         outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
         do
            i = index(outs,text(:nt)) ; IF (i == 0) exit
            outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
         enddo

      endfunction replace_text

      subroutine chdir0(path,err)
      use iso_c_binding
      implicit none
      character(*) :: path
      integer, optional, intent(out) :: err
      integer :: loc_err
      loc_err=c_chdir(trim(path)//c_null_char)

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [chdir0] in mod_utility'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (present(err)) err=loc_err
      end subroutine chdir0

      function real_to_char(num, digit_int) result(char_out)
      implicit none
      character(len=30) :: char
      character(len=30) :: form
      character(len=30) :: digit
      character(len=:), allocatable :: char_out
      integer :: n_out
      real(8) :: num
      integer, optional :: digit_int


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [real_to_char] in mod_utility'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      write(digit,*) abs(digit_int)
      form=""
      if (digit_int>=0) then
         form="(f30."//trim(adjustl(digit))//")"
      else
         form="(Es30."//trim(adjustl(digit))//")"
      endif

      write(char,form) num
      char=trim(adjustl(char))

      n_out=len(trim(char))
      allocate(character(n_out)::char_out)
      char_out=char(1:n_out)

      end function real_to_char

      function int4_to_char(num, digit) result(char_out)
      implicit none
      character(len=30) :: char
      character(len=30) :: form
      character(len=:), allocatable :: char_out
      integer :: n_out
      integer :: num
      integer :: i
      integer, optional :: digit
      integer :: digit0
      integer :: digit1


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [int4_to_char] in mod_utility'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      form="(i30)"
      write(char,form) num
      char=trim(adjustl(char))

      if (present(digit)) then
         digit0=abs(digit)
         digit1=get_digit(num)
         do i=1,digit0-digit1
            if (digit>0) then
               char="0"//trim(char)
            else
               char=" "//trim(char)
            endif
         enddo
      endif

      n_out=len(trim(char))
      allocate(character(n_out)::char_out)
      char_out=char(1:n_out)

      end function int4_to_char

      function int8_to_char(num, digit) result(char_out)
      implicit none
      character(len=30) :: char
      character(len=30) :: form
      character(len=:), allocatable :: char_out
      integer :: n_out
      integer(8) :: num
      integer :: i
      integer, optional :: digit
      integer :: digit0
      integer :: digit1


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [int8_to_char] in mod_utility'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      form="(i30)"
      write(char,form) num
      char=trim(adjustl(char))

      if (present(digit)) then
         digit0=abs(digit)
         digit1=get_digit(int(num,4))
         do i=1,digit0-digit1
            if (digit>0) then
               char="0"//trim(char)
            else
               char=" "//trim(char)
            endif
         enddo
      endif

      n_out=len(trim(char))
      allocate(character(n_out)::char_out)
      char_out=char(1:n_out)

      end function int8_to_char


      function int4_get_digit(num) result(digit0)
      implicit none
      integer :: digit0
      integer :: num
      integer :: t


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [int4_get_digit] in mod_utility'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      t=num
      digit0=1
      do
         if (t<10) then
            exit
          endif
          t=t/10
          digit0=digit0+1
      enddo

      end function int4_get_digit

      subroutine mkdir(folder, clear)
      implicit none
      character(len=*), intent(in) :: folder
      logical, intent(in) :: clear

      !js+! window version
      !js+if (clear) then
      !js+   call system_call ('rmdir /s/q  '//trim(folder)//" 1>NUL 2>NUL")
      !js+   call system_call ('del  '//trim(folder)//" 1>NUL 2>NUL")
      !js+endif
      !js+call system_call ('mkdir '//trim(folder)//" 1>NUL 2>NUL")

      ! linux version

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [mkdir] in mod_utility'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (clear) then
         call system_call('rm -rf '//trim(folder))!//" 1>/dev/null 2>/dev/null")
      endif
      call system_call('mkdir -p '//trim(folder))!//" 1>/dev/null 2>/dev/null")


      return
      end subroutine mkdir

      subroutine system_call(comm)
      implicit none
      character(len=*) :: comm

      !js+! window version
      !js+call execute_command_line(trim(comm))

      ! linux gnu compiler

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [system_call] in mod_utility'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call execute_command_line(trim(comm))

      !js+! linux intel compiler
      !js+call systemqq(trim(comm))

      return
      end subroutine system_call

      subroutine get_pwd(pwd)
      implicit none
      character(len=*) :: pwd
      character(len=10) :: key
      character(len=300) :: name


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_pwd] in mod_utility'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      name=""
      key='PWD'

      ! linux gnu
      call get_environment_variable(key,name)

      !! linux intel
      !call getenvqq(key,name)

      if (name=="") pwd="N/A"
      pwd=name

      return
      end subroutine get_pwd

      end module mod_utility



#endif
