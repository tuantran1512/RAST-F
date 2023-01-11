      module Mod_Watch
      integer(4),parameter :: n_max_time=40
      integer(4)           :: tictoc_time=0
      real(8)              :: time(n_max_time)=0d0
      integer(4)           :: time_start(n_max_time)=0
      integer(4)           :: time_rate(n_max_time)=0
      integer(4)           :: time_end(n_max_time)  =0
      character(len=18)    :: watch_name(n_max_time)=''
      integer(4),parameter :: max_length_time=18
      integer(4)           :: nnn_time=0
      integer(4)           :: n_msg(1:3)=0
      integer :: n_clock, iclock
      end module Mod_Watch

      subroutine WATCH(what,START_or_END)
      use Mod_Watch
      use mod_charedit, only: chartoupper, print_msg
      character(*) :: what
      character(*) :: START_or_END
      integer(4) :: i,ii

      if(len(what)>30) then
         call print_msg(3,'Too long watch name ',what)
         stop
      endif
      do i=1,tictoc_time
         if (trim(adjustl(what))==trim(adjustl(watch_name(i)))) then
            ii=i
            exit
         endif
      enddo
      if(tictoc_time==0 .or. i==tictoc_time+1) then
         tictoc_time=tictoc_time+1
         if(tictoc_time>n_max_time) then
            call print_msg(3,'tictoc_time exceeds n_max_time')
            stop
         endif
         watch_name(tictoc_time)=trim(adjustl(what))
         ii=tictoc_time
      endif

      call chartoupper(START_OR_END)
      if (trim(adjustl(START_OR_END))=='START') then
         call SYSTEM_CLOCK(time_start(ii),time_rate(ii))
      elseif(trim(adjustl(START_OR_END))=='END') then
         call SYSTEM_CLOCK(time_end(ii))
         time(ii)=time(ii)+dble(time_end(ii)-time_start(ii))/max(1d-30,dble(time_rate(ii)))
      else
         call print_msg(3,'#WATCH: START or END?')
         stop
      endif

      if (n_clock<ii) n_clock=ii

      return
      end subroutine WATCH

!#ifdef siarhei_delete 
      subroutine print_message_float(string,print_data,print_unit)
      use inc_file, only: w_anc, flag_anc

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none
      character(*),intent(in) :: string
      character(*),intent(in) :: print_unit
      real(8),intent(in)      :: print_data
      character(len=41)       :: string2
      integer :: target_fname
      logical(1) :: fname_on=.false.
      logical(1), save :: flag_first_print_time=.true.

      write(string2,'(a41)') string // ':'
      string2=adjustl(trim(string2))
      write(*,'(4x,"* ",a41,f19.3,a,a)') string2, print_data,' ',print_unit

      if (flag_anc) then
         target_fname=w_anc
         fname_on=.true.
      endif

      if (fname_on) then
         if (flag_first_print_time) then
            write(target_fname,'(a)') ' *** SIMULATION TIME *** '
            flag_first_print_time=.false.
         endif
         write(target_fname,'(4x,"* ",a41,f19.3,a,a)') string2, print_data,' ',print_unit
      endif

      end subroutine print_message_float
!#endif 

