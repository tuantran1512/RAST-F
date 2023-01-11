      module Inc_PinVar
      implicit none

      integer :: npin
      integer :: ncorn             !number of corner points
      real(8) :: rsqrt2

      integer, ALLOCATABLE :: pploc   (:)
      integer, ALLOCATABLE :: nodec   (:,:)
      integer, ALLOCATABLE :: nxcs    (:)
      integer, ALLOCATABLE :: nxce    (:)
      integer, ALLOCATABLE :: lctox   (:)
      integer, ALLOCATABLE :: lctoy   (:)
      integer, ALLOCATABLE :: nneighc (:)
      integer, ALLOCATABLE :: lcc     (:,:)
      integer, ALLOCATABLE :: lcn     (:,:)
      integer, ALLOCATABLE :: lcnw    (:)
      integer, ALLOCATABLE :: lcsw    (:)
      integer, ALLOCATABLE :: lcne    (:)
      integer, ALLOCATABLE :: lcse    (:)
      integer, ALLOCATABLE :: iquad   (:)
      real(8), ALLOCATABLE :: asnmut   (:,:)
      real(8), ALLOCATABLE :: asnmutl  (:,:)
      real(8), ALLOCATABLE :: asnmuttl (:,:)
      real(8), ALLOCATABLE :: acnmut   (:,:)
      real(8), ALLOCATABLE :: acnmutl  (:,:)
      real(8), ALLOCATABLE :: acnmuttl (:,:)
      real(8), ALLOCATABLE :: asnkat   (:,:)
      real(8), ALLOCATABLE :: asnkatl  (:,:)
      real(8), ALLOCATABLE :: asnkattl (:,:)
      real(8), ALLOCATABLE :: acnkat   (:,:)
      real(8), ALLOCATABLE :: acnkatl  (:,:)
      real(8), ALLOCATABLE :: acnkattl (:,:)
      real(8), ALLOCATABLE :: bctka    (:,:)
      real(8), ALLOCATABLE :: bctmu    (:,:)
      real(8), ALLOCATABLE :: cc11     (:,:,:)
      real(8), ALLOCATABLE :: cc12     (:,:,:)
      real(8), ALLOCATABLE :: cc21     (:,:,:)
      real(8), ALLOCATABLE :: cc22     (:,:,:)
      real(8), ALLOCATABLE :: cur11    (:,:,:)
      real(8), ALLOCATABLE :: cur12    (:,:,:)
      real(8), ALLOCATABLE :: cur21    (:,:,:)
      real(8), ALLOCATABLE :: cur22    (:,:,:)
      real(8), ALLOCATABLE :: cornfka  (:,:,:)
      real(8), ALLOCATABLE :: cornfmu  (:,:,:)
      real(8), ALLOCATABLE :: cpbsrc   (:,:)
      real(8), ALLOCATABLE :: phicorn  (:,:,:)
      real(8), ALLOCATABLE :: phihom   (:,:,:,:,:)
      real(8), ALLOCATABLE :: pppeak   (:,:)
      real(8), ALLOCATABLE :: Peak_X   (:,:)
      real(8), ALLOCATABLE :: Peak_Y   (:,:)
      real(8), ALLOCATABLE :: powvalr  (:,:,:)
      real(8), ALLOCATABLE :: powval   (:,:,:,:)
      real(8), ALLOCATABLE :: RevNorm  (:,:)
      real(8), ALLOCATABLE :: pinpitch (:,:)
      real(8), ALLOCATABLE :: xoffset  (:,:)
      real(8), ALLOCATABLE :: yoffset  (:,:)

      real(8) :: pin_area

      logical(1) :: flag_pin_inform=.false.
      integer, allocatable :: pinmaps(:,:,:)    !(size_x, size_y, as type) : only one size is considered
      real(8), allocatable :: pin_inform(:,:,:) !(var type, pin type, as type)

      END MODULE Inc_PinVar
