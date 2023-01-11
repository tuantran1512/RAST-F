
      module Inc_DF
      implicit none

      real(8), allocatable :: ADF_Lx      (:,:,:)
      real(8), allocatable :: ADF_Rx      (:,:,:)
      real(8), allocatable :: ADF_Ly      (:,:,:)
      real(8), allocatable :: ADF_Ry      (:,:,:)
      real(8), allocatable :: ADF_Lz      (:,:,:)
      real(8), allocatable :: ADF_Rz      (:,:,:)
      real(8), allocatable :: ADF_Avg     (:,:,:)
      real(8), allocatable :: Buff_ADF_Lx (:,:,:)
      real(8), allocatable :: Buff_ADF_Rx (:,:,:)
      real(8), allocatable :: Buff_ADF_Ly (:,:,:)
      real(8), allocatable :: Buff_ADF_Ry (:,:,:)

      real(8), allocatable :: cdf (:,:,:,:)

      real(8), allocatable :: CDF_LxLy (:,:,:)
      real(8), allocatable :: CDF_LxRy (:,:,:)
      real(8), allocatable :: CDF_RxLy (:,:,:)
      real(8), allocatable :: CDF_RxRy (:,:,:)

      ! Axial DF
      logical(1) :: flag_axialdf=.false. !default
      real(8), allocatable :: corr_ADF_Lz(:,:,:)
      real(8), allocatable :: corr_ADF_Rz(:,:,:)
      real(8), allocatable :: phisurf_inpz(:,:,:)
      real(8), allocatable :: phisurf_subz(:,:,:)
      real(8), allocatable :: L_1subz_3D(:,:,:)
      real(8), allocatable :: L_2subz_3D(:,:,:)
      integer, allocatable :: nz_subz(:)
      integer :: nsubz
      integer, allocatable :: inp2sub_beg(:)
      integer, allocatable :: inp2sub_end(:)
      integer, allocatable :: sub2inp(:)
      real(8), allocatable :: hsubz(:)

      end module Inc_DF
