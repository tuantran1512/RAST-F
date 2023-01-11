
      module Inc_maXS
      implicit none

      real(8), allocatable :: D_3D          (:,:,:)
      real(8), allocatable :: maXS_tr_3D    (:,:,:)
      real(8), allocatable :: maXS_a_3D     (:,:,:)
      real(8), allocatable :: maXS_f_3D     (:,:,:)
      real(8), allocatable :: nu_maXS_f_3D  (:,:,:)
      real(8), allocatable :: kap_maXS_f_3D (:,:,:)
      real(8), allocatable :: maXS_s_3D     (:,:,:)
      real(8), allocatable :: maXS_r_3D     (:,:,:)
      real(8), allocatable :: maXS_f_3D_old (:,:,:)
      real(8), allocatable :: maXS_chi_3D   (:,:,:)
      real(8), allocatable :: maXS_chid_3D  (:,:,:)
      real(8), allocatable :: maXS_scat_3D  (:,:,:,:)
#ifdef jr_vver
      real(8), allocatable :: maXS_r_3D_old (:,:,:)
#endif


#ifdef siarhei_tr_hex
      real(8), allocatable :: D_3D_CR          (:,:,:)
      real(8), allocatable :: maXS_tr_3D_CR    (:,:,:)
      real(8), allocatable :: maXS_a_3D_CR     (:,:,:)
      real(8), allocatable :: maXS_f_3D_CR     (:,:,:)
      real(8), allocatable :: nu_maXS_f_3D_CR  (:,:,:)
      real(8), allocatable :: kap_maXS_f_3D_CR (:,:,:)
      real(8), allocatable :: maXS_s_3D_CR     (:,:,:)
      real(8), allocatable :: maXS_r_3D_CR     (:,:,:)
      real(8), allocatable :: maXS_f_3D_old_CR (:,:,:)
      real(8), allocatable :: maXS_chi_3D_CR   (:,:,:)
      real(8), allocatable :: maXS_chid_3D_CR  (:,:,:)
      real(8), allocatable :: maXS_scat_3D_CR  (:,:,:,:)
#endif

      real(8), allocatable :: D_3D_lc      (:,:,:)
      real(8), allocatable :: maXS_a_3D_lc (:,:,:)
      real(8), allocatable :: maXS_s_3D_lc (:,:,:)
      real(8), allocatable :: buck2_lc     (:,:,:)
      real(8), allocatable :: dbuck2_lc    (:,:,:)

      real(8) :: nu
      real(8) :: kappa

      end module Inc_maXS
