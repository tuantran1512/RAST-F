#ifdef JR_SRCTRM
      MODULE Inc_SRC

      INTEGER(4)                                  :: n_neutron_grp=16
      INTEGER(4)                                  :: n_gamma_grp=39
      INTEGER(4)                                  :: n_nuclide=0
      INTEGER(4)                                  :: actp=903
      INTEGER(4)                                  :: acti=174
      INTEGER(4)                                  :: fp=1149
      INTEGER(4)                                  :: tota=2226
      INTEGER(4),ALLOCATABLE                      :: iso6(:)
      INTEGER(4),ALLOCATABLE                      :: ORG_ID(:)
      REAL(8),ALLOCATABLE                         :: gram_den(:)
      REAL(8),ALLOCATABLE                         :: radio(:)                             ! nuclide-wise radioactivity in curies
      REAL(8),ALLOCATABLE                         :: t_power(:)                           ! nuclide-wise thermal power in watts
      REAL(8),ALLOCATABLE                         :: g_power(:)                           ! nuclide-wise gamma power in watts
      REAL(8),ALLOCATABLE                         :: a_nsrc(:)                            ! nuclide-wise alpha-n neutron source in neutrons/sec
      REAL(8),ALLOCATABLE                         :: s_nsrc(:)                            ! nuclide-wise spontaneous fission neutron source in neutrons/sec
      REAL(8),ALLOCATABLE                         :: d_nsrc(:)                            ! nuclide-wise delayed neutron source in neutrons/sec
      REAL(8),ALLOCATABLE                         :: a_nspec(:)                           ! material-wise alpha-n neutron spectra in neutrons/sec (16group)
      REAL(8),ALLOCATABLE                         :: s_nspec(:)                           ! material-wise spontaneous fission neutron spectra in neutrons/sec (16group)
      REAL(8),ALLOCATABLE                         :: t_nspec(:)                           ! material-wise total neutron spectra in neutrons/sec (16group)
      REAL(8),ALLOCATABLE                         :: d_nspec(:)                           ! material-wise delayed neutron spectra in neutrons/sec (16group)
      REAL(8),ALLOCATABLE                         :: gespec(:)                            ! material-wise gamma energy release rates in MeV/sec (38group)
      REAL(8),ALLOCATABLE                         :: gspec(:)                             ! material-wise gamma spectra in photons/sec (38group)
      REAL(8),ALLOCATABLE                         :: pseb(:)                              ! photon spectra energy boundaries  (MeV)
      REAL(8),ALLOCATABLE                         :: nseb(:)                              ! neutron spectra energy boundaries
      REAL(8),ALLOCATABLE                         :: t_nspec_wdn(:)
      REAL(8)                                     :: ttnspec_wdn=0d0
      REAL(8)                                     :: tdh=0d0
      REAL(8)                                     :: tac=0d0
      REAL(8)                                     :: tgp=0d0
      REAL(8)                                     :: tgrd=0d0
      REAL(8)                                     :: tansrc=0d0
      REAL(8)                                     :: tsnsrc=0d0
      REAL(8)                                     :: tdnsrc=0d0
      REAL(8)                                     :: tanspec=0d0
      REAL(8)                                     :: tsnspec=0d0
      REAL(8)                                     :: tdnspec=0d0
      REAL(8)                                     :: ttnspec=0d0
      REAL(8)                                     :: tgespec=0d0
      REAL(8)                                     :: tgspec=0d0

      END MODULE Inc_SRC
#endif
