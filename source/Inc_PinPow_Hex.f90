   module Inc_PinPow_Hex

      implicit none
      ! variables for reading HFF data
      real(8),allocatable                    :: HFF_hex(:,:,:) ! 2D HFF plux FA index
      real(8),allocatable                    :: zero_HFF_hex(:,:) ! all zero elements (empty)
      real(8),allocatable                    :: inp_HFF_hex(:,:,:) ! input format of HFF
      integer                                :: radNum ! number of pins in radial direction (incl central tube)
      integer                                :: N_FA_parts ! FA sliced in how many parts (6 - triangle, 2 - half, 1 - full hex)
      real(8)                                :: radPtch ! fuel pin pitch

      ! variables for geometry/flux assignment
      real(8),allocatable,dimension(:,:,:,:,:) :: xy_Flux_hex !(rows,columns,nassy,n_group,Nz)
      integer                                :: coord_ppr ! imap(ih) - radial location
      integer                                :: coord_iz ! iz - axial location
      real(8)                                :: triangle_height_hex
      real(8)                                :: center_x_triangle_hex
      integer,allocatable,dimension(:,:)     :: saved_nodel_hex ! save (x,y) node config
      integer,allocatable,dimension(:,:)     :: saved_nodel_hex_adj ! removed zeros between FAs
      real(8),allocatable,dimension(:,:)     :: x_coord_fa_hex
      real(8),allocatable,dimension(:,:)     :: y_coord_fa_hex
      real(8),allocatable,dimension(:,:)     :: p_coord_fa_hex
      real(8),allocatable,dimension(:,:)     :: u_coord_fa_hex
      real(8),allocatable,dimension(:,:)     :: rad_coord_fa_hex
      real(8),allocatable,dimension(:,:)     :: sin_coord_fa_hex
      real(8),allocatable,dimension(:,:)     :: ang_coord_fa_hex
      real(8),allocatable,dimension(:,:,:)   :: x_y_u_p_coord_fa_hex
      integer,allocatable,dimension(:,:)     :: location_in_hex

    !  ! remove later
    !  integer                                :: nassy = 1 !<- remove later (temp)
    !  integer                                :: n_group = 1 !<- remove later (temp)

     ! ! variables for setting up Ф(x,y) smooth flux function | alloc(6,n_group,nassy)
     ! real(8),allocatable,dimension(:,:,:)   :: C_g0, A_gx, A_gy
     ! real(8),allocatable,dimension(:,:,:)   :: B_gx, B_gu, B_gp
     ! real(8),allocatable,dimension(:,:,:)   :: C_gx, C_gu, C_gp

      ! fluxes for calculating expansion coefficients | alloc(6,n_group,nassy,Nz)
      real(8),allocatable,dimension(:,:,:,:) :: Fl_momx,Fl_momy,Fl_avg
      real(8),allocatable,dimension(:,:,:)   :: corner_flux_h
      real(8),allocatable,dimension(:,:,:,:) :: inner_b_flx_h,outer_b_flx_h
      real(8),allocatable,dimension(:,:,:)   :: central_flux_h

      ! fluxes for calculating expansion coefficients | alloc(6,n_group,nassy)
      ! averaged over Z-axis
      real(8),allocatable,dimension(:,:,:) :: Fl_momx_av,Fl_momy_av,Fl_avg_av
      real(8),allocatable,dimension(:,:)   :: corner_flux_h_av
      real(8),allocatable,dimension(:,:,:) :: inner_b_flx_h_av,outer_b_flx_h_av
      real(8),allocatable,dimension(:,:)   :: central_flux_h_av

      real(8),allocatable,dimension(:,:,:)   :: Fl_c_x,Fl_c_p,Fl_c_u ! maybe not needed
      real(8),allocatable,dimension(:,:,:)   :: Fl_s_x,Fl_s_p,Fl_s_u ! maybe not needed

      logical                                :: pin_pow_hex_needed = .false.

      ! calculation of pin power
                                                ! (xpin,ypin,nassy,Nz)
      ! Fq -axial mesh rod linear power/average core axial mesh rod linear power [kW/cm]/[kW/cm]
      ! FdH - integral rod linear power/core average intergral rod linear power [kW/rod]/[kW/rod]
      ! Fxy - same as Fq but only for one axial layer (layer-averaged, not core-averaged denominator)
      ! Fz - average average axial power at elevation Z (Core_Power(Z)/Volume(Z))
      real(8),allocatable,dimension(:,:,:,:)   :: hom_pow_hex ! total MG power in each pin loc
      real(8),allocatable,dimension(:,:,:,:)   :: xy_pow_hex ! ^ prev value multiplied by HFF
      real(8),allocatable,dimension(:,:,:,:)   :: lin_pow_hex ! ^ prev value in W/cm
      real(8),allocatable,dimension(:,:,:,:)   :: xy_Fq_hex ! Fq factor for each pin mesh
      real(8),allocatable,dimension(:,:,:,:)   :: xy_Fxy_hex ! Fxy factor for each pin mesh
      real(8),allocatable,dimension(:)         :: xy_Fz_hex ! Fz factor for each Z-layer
      real(8),allocatable,dimension(:,:,:)     :: integral_pin_pow_hex ! rod power for FdH
      real(8),allocatable,dimension(:,:,:)     :: xy_FdH_hex ! FdH factor for each pin

   end module Inc_PinPow_Hex
