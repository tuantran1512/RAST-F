
      module Inc_XYZ
      implicit none

      real(8) :: ASI

      real(8), allocatable :: normal_power_z   (:)
      real(8), allocatable :: normal_power_xy  (:,:)
      real(8), allocatable :: normal_power_xyz (:,:,:)
      real(8), allocatable :: bu_z             (:)
      real(8), allocatable :: bu_xy            (:,:)
      real(8), allocatable :: bu_xyz           (:,:,:)
      real(8), allocatable :: t_fuel_z         (:)
      real(8), allocatable :: t_fuel_xy        (:,:)
      real(8), allocatable :: t_fuel_xyz       (:,:,:)
      real(8), allocatable :: t_mod_z          (:)
      real(8), allocatable :: t_mod_xy         (:,:)
      real(8), allocatable :: t_mod_xyz        (:,:,:)
      real(8), allocatable :: d_mod_z          (:)
      real(8), allocatable :: d_mod_xy         (:,:)
      real(8), allocatable :: d_mod_xyz        (:,:,:)


      real(8), allocatable ::  FA_ASI (:,:)

      real(8), allocatable :: N_Xe35_XYZ (:,:,:)
      real(8), allocatable :: N_Sm49_XYZ (:,:,:)
      real(8), allocatable :: FFlux_XYZ  (:,:,:)
      real(8), allocatable :: TFlux_XYZ  (:,:,:)

      real(8), allocatable :: fFlux_adj_XYZ (:,:,:)
      real(8), allocatable :: tFlux_adj_XYZ (:,:,:)
      real(8), allocatable :: fFlux_adj_XY  (:,:)
      real(8), allocatable :: tFlux_adj_XY  (:,:)

      real(8), allocatable :: Normal_Power_XY_4N (:,:)
      real(8), allocatable :: BU_XY_4N           (:,:)

      real(8), allocatable :: Fxy_val_XY (:,:)
      real(8), allocatable :: FdH_val_XY (:,:)

      end module Inc_XYZ

