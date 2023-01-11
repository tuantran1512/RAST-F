      module Inc_Lscoef
      implicit none

      ! alloc_lscoef
      real(8), allocatable :: am    (:,:,:)
      real(8), allocatable :: amcc  (:,:,:)
      real(8), allocatable :: af    (:,:,:)
      real(8), allocatable :: af2   (:,:)
      real(8), allocatable :: scat  (:,:)
      real(8), allocatable :: scatv (:,:,:,:)
      real(8), allocatable :: ccw   (:,:,:)
      real(8), allocatable :: cce   (:,:,:)
      real(8), allocatable :: ccn   (:,:,:)
      real(8), allocatable :: ccs   (:,:,:)
      real(8), allocatable :: ccb   (:,:,:)
      real(8), allocatable :: cct   (:,:,:)
      real(8), allocatable :: dfw   (:,:,:)
      real(8), allocatable :: dfn   (:,:,:)
      real(8), allocatable :: dfb   (:,:,:)
      real(8), allocatable :: dnw   (:,:,:)
      real(8), allocatable :: dnn   (:,:,:)
      real(8), allocatable :: dnb   (:,:,:)

      ! alloc_lu
      real(8), allocatable :: del    (:,:)
      real(8), allocatable :: au     (:,:)
      real(8), allocatable :: ainvl  (:,:)
      real(8), allocatable :: ainvu  (:,:)
      real(8), allocatable :: ainvd  (:,:)
      real(8), allocatable :: delinv (:,:,:)
      real(8), allocatable :: al     (:,:,:)
      real(8), allocatable :: deliau (:,:,:)

      end module Inc_Lscoef
