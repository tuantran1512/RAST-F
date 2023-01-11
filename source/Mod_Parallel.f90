      module mod_parallel
      use inc_parallel


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

#ifdef js_mpi
      interface bcast
         module procedure bcast_real8
         module procedure bcast_real8_1d
         module procedure bcast_real8_2d
         module procedure bcast_real8_3d
         module procedure bcast_real8_4d
         module procedure bcast_real4_1d
         module procedure bcast_int4
         module procedure bcast_int4_1d
         module procedure bcast_int4_2d
         module procedure bcast_int4_3d
         module procedure bcast_int4_4d
         module procedure bcast_int1
         module procedure bcast_int1_1d
         module procedure bcast_int1_2d
         module procedure bcast_int1_3d
         module procedure bcast_int1_4d
         module procedure bcast_logical4
         module procedure bcast_logical4_1d
         module procedure bcast_logical1
         module procedure bcast_logical1_1d
      end interface bcast

      interface allreduce
         module procedure allreduce_int4_0d
         module procedure allreduce_int4_1d
         module procedure allreduce_int4_2d
         module procedure allreduce_int4_3d
         module procedure allreduce_int4_4d
         module procedure allreduce_int4_5d
         module procedure allreduce_int4_6d
         !module procedure allreduce_int1_4d
         !module procedure allreduce_int2_4d
         module procedure allreduce_real8_0d
         module procedure allreduce_real8_1d
         module procedure allreduce_real8_2d
         module procedure allreduce_real8_3d
         module procedure allreduce_real8_4d
         module procedure allreduce_real8_5d
         module procedure allreduce_real8_6d
         module procedure allreduce_real8_7d
         !module procedure allreduce_real4_0d
         !module procedure allreduce_real4_1d
         !module procedure allreduce_real4_2d
         !module procedure allreduce_real4_3d
         !module procedure allreduce_real4_4d
         !module procedure allreduce_real4_5d
         !module procedure allreduce_real4_6d
      end interface allreduce

      !interface scatterv
      !   module procedure scatterv_real8_1d
      !   module procedure scatterv_real8_2d
      !   !module procedure scatterv_real4_1d
      !   !module procedure scatterv_real4_2d
      !end interface scatterv

      interface reduce
         module procedure reduce_real8_1d
         module procedure reduce_real8_2d
         module procedure reduce_real8_3d
         module procedure reduce_real8_4d
         !module procedure reduce_real4_1d
      end interface reduce

      !interface gather
      !   module procedure gather_real8
      !   module procedure gather_real8_1d
      !   module procedure gather_real8_2d
      !end interface gather
#endif

      contains

      subroutine init_comm
      implicit none


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [init_comm] in mod_parallel'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      comm%n_procs = 1
      comm%i_proc = 0
      comm%i_master = 0
      comm%if_master = .true.
      comm%usempi = .false.
#ifdef js_mpi
!      call init_parallel
            dummy_filler = 1 ! @$^ siarhei_plot 
      iproc=comm%i_proc+1
#endif
      return
      end subroutine init_comm

#ifdef js_mpi
#ifdef siarhei_delete 
      subroutine init_parallel
      implicit none

      call mpi_init(error_mpi)
      comm%i_comm = mpi_comm_world
      comm%tag = 1

      call mpi_comm_size(comm%i_comm, comm%n_procs, error_mpi)
      call mpi_comm_rank(comm%i_comm, comm%i_proc, error_mpi)

      if (comm%i_proc == 0) then
         comm%if_master = .true.
      else
         comm%if_master = .false.
      endif

      if (comm%n_procs>1) then
         comm%usempi=.true.
      endif

      return
      end subroutine init_parallel
#endif 


#ifdef siarhei_delete 
      subroutine final_parallel
      implicit none

!      call barrier
            dummy_filler = 1 ! @$^ siarhei_plot 
      call mpi_finalize(error_mpi)

      return
      end subroutine final_parallel
#endif 
#endif

      subroutine end_parallel
#ifdef js_mpi
      include 'mpif.h'
#endif
#ifdef js_mpi
      error_mpi=0
      status(1:mpi_status_size)=0
#endif
      n_ompths=0
      end subroutine end_parallel

#ifdef js_mpi
      ! ------------------------------------------------------------------------
#ifdef siarhei_delete 
      subroutine barrier
      implicit none
      call mpi_barrier(mpi_comm_world,error_mpi)
      return
      end subroutine barrier
#endif 
      ! ------------------------------------------------------------------------

      ! ------------------------------------------------------------------------
#ifdef siarhei_delete 
      subroutine bcast_real8(tmp1,i_proc)
      implicit none
      real(8), intent(inout) :: tmp1
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, 1, mpi_real8, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_real8
#endif 

#ifdef siarhei_delete 
      subroutine bcast_real8_1d(tmp1,n1,i_proc)
      implicit none
      real(8), intent(inout) :: tmp1(:)
      integer, intent(in) :: n1
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1, mpi_real8, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_real8_1d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_real8_2d(tmp1,n1,n2,i_proc)
      implicit none
      real(8), intent(inout) :: tmp1(:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1*n2, mpi_real8, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_real8_2d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_real8_3d(tmp1,n1,n2,n3,i_proc)
      implicit none
      real(8), intent(inout) :: tmp1(:,:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1*n2*n3, mpi_real8, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_real8_3d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_real8_4d(tmp1,n1,n2,n3,n4,i_proc)
      implicit none
      real(8), intent(inout) :: tmp1(:,:,:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: n4
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1*n2*n3*n4, mpi_real8, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_real8_4d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_real4_1d(tmp1,n1,i_proc)
      implicit none
      real(4), intent(inout) :: tmp1(:)
      integer, intent(in) :: n1
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1, mpi_real4, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_real4_1d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_int4(tmp1,i_proc)
      implicit none
      integer(4), intent(inout) :: tmp1
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, 1, mpi_integer4, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_int4
#endif 

#ifdef siarhei_delete 
      subroutine bcast_int4_1d(tmp1,n1,i_proc)
      implicit none
      integer(4), intent(inout) :: tmp1(:)
      integer, intent(in) :: n1
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1, mpi_integer4, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_int4_1d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_int4_2d(tmp1,n1,n2,i_proc)
      implicit none
      integer(4), intent(inout) :: tmp1(:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1*n2, mpi_integer4, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_int4_2d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_int4_3d(tmp1,n1,n2,n3,i_proc)
      implicit none
      integer(4), intent(inout) :: tmp1(:,:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1*n2*n3, mpi_integer4, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_int4_3d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_int4_4d(tmp1,n1,n2,n3,n4,i_proc)
      implicit none
      integer(4), intent(inout) :: tmp1(:,:,:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: n4
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1*n2*n3*n4, mpi_integer4, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_int4_4d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_int1(tmp1,i_proc)
      implicit none
      integer(1), intent(inout) :: tmp1(:)
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, mpi_integer1, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_int1
#endif 

#ifdef siarhei_delete 
      subroutine bcast_int1_1d(tmp1,n1,i_proc)
      implicit none
      integer(1), intent(inout) :: tmp1(:)
      integer, intent(in) :: n1
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1, mpi_integer1, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_int1_1d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_int1_2d(tmp1,n1,n2,i_proc)
      implicit none
      integer(1), intent(inout) :: tmp1(:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1*n2, mpi_integer1, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_int1_2d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_int1_3d(tmp1,n1,n2,n3,i_proc)
      implicit none
      integer(1), intent(inout) :: tmp1(:,:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1*n2*n3, mpi_integer1, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_int1_3d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_int1_4d(tmp1,n1,n2,n3,n4,i_proc)
      implicit none
      integer(1), intent(inout) :: tmp1(:,:,:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: n4
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1*n2*n3*n4, mpi_integer1, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_int1_4d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_logical4(tmp1,i_proc)
      implicit none
      logical(4), intent(inout) :: tmp1
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, 1, mpi_logical, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_logical4
#endif 

#ifdef siarhei_delete 
      subroutine bcast_logical4_1d(tmp1,n1,i_proc)
      implicit none
      logical(4), intent(inout) :: tmp1(:)
      integer, intent(in) :: n1
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1, mpi_logical, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_logical4_1d
#endif 

#ifdef siarhei_delete 
      subroutine bcast_logical1(tmp1,i_proc)
      implicit none
      logical(1), intent(inout) :: tmp1
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, 1, mpi_logical, i_proc, mpi_comm_world, error_mpi)
      !js+call mpi_bcast(tmp1, 1, mpi_logical1, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_logical1
#endif 

#ifdef siarhei_delete 
      subroutine bcast_logical1_1d(tmp1,n1,i_proc)
      implicit none
      logical(1), intent(inout) :: tmp1(:)
      integer, intent(in) :: n1
      integer, intent(in) :: i_proc
      call mpi_bcast(tmp1, n1, mpi_logical, i_proc, mpi_comm_world, error_mpi)
      return
      end subroutine bcast_logical1_1d
#endif 
      ! ------------------------------------------------------------------------

      ! ------------------------------------------------------------------------
#ifdef siarhei_delete 
      subroutine allreduce_int4_0d(dat)
      implicit none
      integer(4) :: dat
      integer(4) :: buff
      buff=dat
      call mpi_allreduce(buff, dat, 1, mpi_integer4, mpi_sum, mpi_comm_world, error_mpi)
      return
      end subroutine allreduce_int4_0d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_int4_1d(dat, n1)
      implicit none
      integer(4) :: dat(:)
      integer(4), allocatable :: buff(:)
      integer :: n1
      allocate(buff(1:n1))
      buff=dat
      call mpi_allreduce(buff, dat, n1, mpi_integer4, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_int4_1d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_int4_2d(dat, n1, n2)
      implicit none
      integer(4) :: dat(:,:)
      integer(4), allocatable :: buff(:,:)
      integer :: n1
      integer :: n2
      allocate(buff(1:n1,1:n2))
      buff=dat
      call mpi_allreduce(buff, dat, n1*n2, mpi_integer4, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_int4_2d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_int4_3d(dat, n1, n2, n3)
      implicit none
      integer(4) :: dat(:,:,:)
      integer(4), allocatable :: buff(:,:,:)
      integer :: n1
      integer :: n2
      integer :: n3
      allocate(buff(1:n1,1:n2,1:n3))
      buff=dat
      call mpi_allreduce(buff, dat, n1*n2*n3, mpi_integer4, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_int4_3d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_int4_4d(dat, n1, n2, n3, n4)
      implicit none
      integer(4) :: dat(:,:,:,:)
      integer(4), allocatable :: buff(:,:,:,:)
      integer :: n1
      integer :: n2
      integer :: n3
      integer :: n4
      allocate(buff(1:n1,1:n2,1:n3,1:n4))
      buff=dat
      call mpi_allreduce(buff, dat, n1*n2*n3*n4, mpi_integer4, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_int4_4d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_int4_5d(dat, n1, n2, n3, n4, n5)
      implicit none
      integer(4) :: dat(:,:,:,:,:)
      integer(4), allocatable :: buff(:,:,:,:,:)
      integer :: n1
      integer :: n2
      integer :: n3
      integer :: n4
      integer :: n5
      allocate(buff(1:n1,1:n2,1:n3,1:n4,1:n5))
      buff=dat
      call mpi_allreduce(buff, dat, n1*n2*n3*n4*n5, mpi_integer4, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_int4_5d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_int4_6d(dat, n1, n2, n3, n4, n5, n6)
      implicit none
      integer(4) :: dat(:,:,:,:,:,:)
      integer(4), allocatable :: buff(:,:,:,:,:,:)
      integer :: n1
      integer :: n2
      integer :: n3
      integer :: n4
      integer :: n5
      integer :: n6
      allocate(buff(1:n1,1:n2,1:n3,1:n4,1:n5,1:n6))
      buff=dat
      call mpi_allreduce(buff, dat, n1*n2*n3*n4*n5*n6, mpi_integer4, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_int4_6d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_real8_0d(dat)
      implicit none
      real(8) :: dat
      real(8) :: buff
      buff=dat
      call mpi_allreduce(buff, dat, 1, mpi_real8, mpi_sum, mpi_comm_world, error_mpi)
      return
      end subroutine allreduce_real8_0d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_real8_1d(dat, n1)
      implicit none
      real(8) :: dat(:)
      real(8), allocatable :: buff(:)
      integer :: n1
      allocate(buff(1:n1))
      buff=dat
      call mpi_allreduce(buff, dat, n1, mpi_real8, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_real8_1d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_real8_2d(dat, n1, n2)
      implicit none
      real(8) :: dat(:,:)
      real(8), allocatable :: buff(:,:)
      integer :: n1
      integer :: n2
      allocate(buff(1:n1,1:n2))
      buff=dat
      call mpi_allreduce(buff, dat, n1*n2, mpi_real8, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_real8_2d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_real8_3d(dat, n1, n2, n3)
      implicit none
      real(8) :: dat(:,:,:)
      real(8), allocatable :: buff(:,:,:)
      integer :: n1
      integer :: n2
      integer :: n3
      allocate(buff(1:n1,1:n2,1:n3))
      buff=dat
      call mpi_allreduce(buff, dat, n1*n2*n3, mpi_real8, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_real8_3d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_real8_4d(dat, n1, n2, n3, n4)
      implicit none
      real(8) :: dat(:,:,:,:)
      real(8), allocatable :: buff(:,:,:,:)
      integer :: n1
      integer :: n2
      integer :: n3
      integer :: n4
      allocate(buff(1:n1,1:n2,1:n3,1:n4))
      buff=dat
      call mpi_allreduce(buff, dat, n1*n2*n3*n4, mpi_real8, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_real8_4d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_real8_5d(dat, n1, n2, n3, n4, n5)
      implicit none
      real(8) :: dat(:,:,:,:,:)
      real(8), allocatable :: buff(:,:,:,:,:)
      integer :: n1
      integer :: n2
      integer :: n3
      integer :: n4
      integer :: n5
      allocate(buff(1:n1,1:n2,1:n3,1:n4,1:n5))
      buff=dat
      call mpi_allreduce(buff, dat, n1*n2*n3*n4*n5, mpi_real8, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_real8_5d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_real8_6d(dat, n1, n2, n3, n4, n5, n6)
      implicit none
      real(8) :: dat(:,:,:,:,:,:)
      real(8), allocatable :: buff(:,:,:,:,:,:)
      integer :: n1
      integer :: n2
      integer :: n3
      integer :: n4
      integer :: n5
      integer :: n6
      allocate(buff(1:n1,1:n2,1:n3,1:n4,1:n5,1:n6))
      buff=dat
      call mpi_allreduce(buff, dat, n1*n2*n3*n4*n5*n6, mpi_real8, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_real8_6d
#endif 

#ifdef siarhei_delete 
      subroutine allreduce_real8_7d(dat, n1, n2, n3, n4, n5, n6, n7)
      implicit none
      real(8) :: dat(:,:,:,:,:,:,:)
      real(8), allocatable :: buff(:,:,:,:,:,:,:)
      integer :: n1
      integer :: n2
      integer :: n3
      integer :: n4
      integer :: n5
      integer :: n6
      integer :: n7
      allocate(buff(1:n1,1:n2,1:n3,1:n4,1:n5,1:n6,1:n7))
      buff=dat
      call mpi_allreduce(buff, dat, n1*n2*n3*n4*n5*n6*n7, mpi_real8, mpi_sum, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine allreduce_real8_7d
#endif 
      ! ------------------------------------------------------------------------

      ! ------------------------------------------------------------------------
      !subroutine gather_real8OC_1d(i_proc,sendbuf,recvbuf)
      !
      !end subroutine scatterv_real8_1d
      ! ------------------------------------------------------------------------

      ! ------------------------------------------------------------------------
      !subroutine scatterv_real8_1d(i_proc,sendbuf,recvbuf)
      !
      !end subroutine scatterv_real8_1d
      ! ------------------------------------------------------------------------

      ! ------------------------------------------------------------------------
#ifdef siarhei_delete 
      subroutine reduce_real8_1d(dat,n1,i_proc)
      implicit none
      real(8), intent(inout) :: dat(:)
      real(8), allocatable :: buff(:)
      integer, intent(in) :: n1
      integer, intent(in) :: i_proc
      allocate(buff(1:n1))
      buff=dat
      call mpi_reduce(buff, dat, n1, mpi_real8, mpi_sum, i_proc, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine reduce_real8_1d
#endif 

#ifdef siarhei_delete 
      subroutine reduce_real8_2d(dat,n1,n2,i_proc)
      implicit none
      real(8), intent(inout) :: dat(:,:)
      real(8), allocatable :: buff(:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: i_proc
      allocate(buff(1:n1,1:n2))
      buff=dat
      call mpi_reduce(buff, dat, n1*n2, mpi_real8, mpi_sum, i_proc, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine reduce_real8_2d
#endif 

#ifdef siarhei_delete 
      subroutine reduce_real8_3d(dat,n1,n2,n3,i_proc)
      implicit none
      real(8), intent(inout) :: dat(:,:,:)
      real(8), allocatable :: buff(:,:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: i_proc
      allocate(buff(1:n1,1:n2,1:n3))
      buff=dat
      call mpi_reduce(buff, dat, n1*n2*n3, mpi_real8, mpi_sum, i_proc, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine reduce_real8_3d
#endif 

#ifdef siarhei_delete 
      subroutine reduce_real8_4d(dat,n1,n2,n3,n4,i_proc)
      implicit none
      real(8), intent(inout) :: dat(:,:,:,:)
      real(8), allocatable :: buff(:,:,:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: n4
      integer, intent(in) :: i_proc
      allocate(buff(1:n1,1:n2,1:n3,1:n4))
      buff=dat
      call mpi_reduce(buff, dat, n1*n2*n3*n4, mpi_real8, mpi_sum, i_proc, mpi_comm_world, error_mpi)
      deallocate(buff)
      return
      end subroutine reduce_real8_4d
#endif 
      ! ------------------------------------------------------------------------
#endif

#ifdef siarhei_delete 
      subroutine set_parallel

      implicit none

      if (comm%usempi) comm%usempi=.false.
      return


      return
      end subroutine set_parallel
#endif 


      end module mod_parallel
