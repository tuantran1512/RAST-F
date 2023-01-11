
      module Inc_Kinetics

      implicit none

      integer :: N_Group_d
      real(8), allocatable :: beta_d_Tot (:,:)
      real(8), allocatable :: beta_d     (:,:,:)
      real(8), allocatable :: beta_d_eff (:,:,:)
      real(8), allocatable :: lambda_d   (:,:,:)
      real(8), allocatable :: v_Inv      (:,:,:)

      end module Inc_Kinetics
