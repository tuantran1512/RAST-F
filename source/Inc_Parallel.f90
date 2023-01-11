      module inc_parallel
#ifdef js_mpi
      include 'mpif.h'
#endif
      type :: comm_type
         integer :: n_procs=1
         integer :: n_proc=0
         integer :: i_proc=0
         integer :: i_master=0
         integer :: i_comm=0
         integer :: tag=0
         logical(1) :: if_master=.true.
         logical(1) :: usempi=.false.
      end type comm_type

      type (comm_type) :: comm
#ifdef js_mpi
      integer :: error_mpi=0
      integer :: status(1:mpi_status_size)=0
      integer :: iproc
#endif
      integer :: n_ompths=0

      type :: arrange_type
         integer, allocatable :: n4pxs(:)
         integer, allocatable :: n4pxs_fa(:)
         integer, allocatable :: n4pxs_rf(:)

         integer, allocatable :: iproc2ixs(:)
         integer, allocatable :: iproc2n(:)
         integer, allocatable :: mat_ipixs(:,:)
         integer, allocatable :: matn_ipixs(:,:)

         ! for TH
         integer, allocatable :: nxy_fa4proc(:)

         ! for neutronics
         integer, allocatable :: nxy4proc(:)
      end type arrange_type
      type (arrange_type) :: arr

      integer, allocatable :: node2iproc(:,:)
      integer, allocatable :: mat_ipixs(:,:)
      integer, allocatable :: matn_ipixs(:,:)

      integer, allocatable :: ixyiz2ixyz(:,:)
      integer, allocatable :: ixyz2ixy(:)
      integer, allocatable :: ixyz2iz(:)
      integer, allocatable :: iproc2nxyz(:)
      integer, allocatable :: ixyzip2ixyz(:,:)

      ! for TH
      integer, allocatable :: ixy_fa2iproc(:)
      integer, allocatable :: iproc2nxy_fa(:)
      integer, allocatable :: ixy_faip2ixy_fa(:,:)

      ! for neutronics
      integer, allocatable :: ixy2iproc(:)
      integer, allocatable :: iproc2nxy(:)
      integer, allocatable :: ixyip2ixy(:,:)
      integer, allocatable :: ixy4p_l2g(:)

      ! nodel
      integer, allocatable :: nodelmpi(:,:)
      integer :: nnodelmpi

      real(8), allocatable :: buff_ndnxynz(:,:,:)
      real(8), allocatable :: buff_nxynz(:,:)
      real(8), allocatable :: buff_nxynzng(:,:,:)
      real(8), allocatable :: buff_nxynzmg(:,:,:)
      real(8), allocatable :: buff_ngngnxynz(:,:,:,:)

      end module inc_parallel
