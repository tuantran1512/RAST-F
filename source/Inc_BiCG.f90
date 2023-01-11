
      MODULE Inc_BiCG
      implicit none

      real(8) :: calpha
      real(8) :: cbeta
      real(8) :: crho
      real(8) :: comega

      real(8), allocatable :: vr0  (:,:,:)
      real(8), allocatable :: vp   (:,:,:)
      real(8), allocatable :: vv   (:,:,:)
      real(8), allocatable :: vs   (:,:,:)
      real(8), allocatable :: vt   (:,:,:)
      real(8), allocatable :: vy   (:,:,:)
      real(8), allocatable :: vz   (:,:,:)
      real(8), allocatable :: s    (:,:)
      real(8), allocatable :: b0   (:,:)
      real(8), allocatable :: s1dl (:,:)
      real(8), allocatable :: b01d (:,:)
      real(8), allocatable :: y    (:,:)

      real(8) :: b2
      real(8) :: r2
      real(8) :: r20
      real(8) :: r2ob2

      end module Inc_BiCG
