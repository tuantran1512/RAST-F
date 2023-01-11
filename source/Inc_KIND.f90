
      module inc_kind

      integer, parameter :: l1k=selected_int_kind(2)       ! logical(1)
      integer, parameter :: i1k=selected_int_kind(2)       ! integer(1)
      integer, parameter :: i4k=selected_int_kind(9)       ! integer(4)
      integer, parameter :: i8k=selected_int_kind(18)      ! integer(8)
      integer, parameter :: r4k=selected_real_kind(6,37)   ! real(4)
      integer, parameter :: r8k=selected_real_kind(15,307) ! real(8)

      end module inc_kind
