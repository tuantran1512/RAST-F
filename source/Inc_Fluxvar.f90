
      module Inc_FluxVar
      implicit none

      real(8), allocatable :: flux_add   (:,:,:)
      real(8), allocatable :: vr         (:,:,:)

      real(8), allocatable :: src        (:,:,:)
      real(8), allocatable :: src_tr     (:,:,:)

      real(8), allocatable :: j_net_x_3D (:,:,:)
      real(8), allocatable :: j_net_y_3D (:,:,:)
      real(8), allocatable :: j_net_z_3D (:,:,:)

      real(8), allocatable :: L_0x_3D    (:,:,:)
      real(8), allocatable :: L_1x_3D    (:,:,:)
      real(8), allocatable :: L_2x_3D    (:,:,:)
      real(8), allocatable :: L_0y_3D    (:,:,:)
      real(8), allocatable :: L_1y_3D    (:,:,:)
      real(8), allocatable :: L_2y_3D    (:,:,:)
      real(8), allocatable :: L_0z_3D    (:,:,:)
      real(8), allocatable :: L_1z_3D    (:,:,:)
      real(8), allocatable :: L_2z_3D    (:,:,:)

      end module Inc_FluxVar
