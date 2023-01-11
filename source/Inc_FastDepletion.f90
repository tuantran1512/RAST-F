! constants types and variables for FastDepletion.f90 module
! and for fast reactor depletion in general

   module Inc_FastDepletion


      implicit none

      type :: reaction_product_data
         character(10)           :: what_reaction
         character(10)           :: product_name
         integer                 :: coord ! product coordinate
         real(8)                 :: probability
      end type reaction_product_data

      type :: nuclide_wiki
         character(10)           :: name
         integer                 :: ID
         real(8)                 :: decay_const
         integer                 :: fission_branch
         real(8),allocatable     :: fission_yield(:)
         integer                 :: len_neutron
         integer                 :: len_decay
         type(reaction_product_data),allocatable :: neutron(:)
         type(reaction_product_data),allocatable :: decay(:)
      end type nuclide_wiki

      type(nuclide_wiki),allocatable :: all_names(:)

      logical                    :: block1 = .false.
      logical                    :: block2 = .false.
      logical                    :: block3 = .false.
      logical                    :: headline = .false.
      integer                    :: number_hm
      integer                    :: number_fp
      integer                    :: numNuk,numHM,numFP,maxDPerN,numDpath
      integer                    :: maxNPerN,numNpath,numFPYset
      character(10),allocatable  :: hm_names(:)
      character(10),allocatable  :: fp_names(:)
      integer,allocatable        :: fission_flag(:)

      complex(8),allocatable     :: cram_alpha(:)
      complex(8),allocatable     :: cram_theta(:)
      complex(8),allocatable     :: sum_CRAM_vec(:)
      real(8)                    :: cram_alpha_0
      integer                    :: order_cram


      real(8),allocatable        :: A_matrix(:,:)
      real(8),allocatable        :: vec_N_old(:)
      complex(8),allocatable     :: mat_I_all(:,:)

      real(8), dimension(:,:,:,:), allocatable ::xs_f
      real(8), dimension(:,:,:,:), allocatable ::xs_n2n
      real(8), dimension(:,:,:,:), allocatable ::xs_a

      real(8), allocatable         :: kappaf(:)
      real(8),allocatable        :: all_concentrations(:,:,:,:) 
#ifdef tuan_fr_tdep
      real(8),allocatable        :: A_matrix_t(:,:,:)
      real(8),allocatable        :: NDens_T(:,:,:,:,:)

#endif
   end module Inc_FastDepletion
