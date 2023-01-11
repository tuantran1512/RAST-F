#ifdef jr_vver
      MODULE Inc_TPEN
      IMPLICIT NONE
! ============================================================================= !
! NOTICE
! ============================================================================= !
! 1. Nxy = Nxy_1N (in hexagonal geometry)
!
      LOGICAL(1) :: if_cycle2 = .false.

#ifdef tuan_tr_test
      LOGICAL(1) :: flag_jr=.false.
      REAL(8)    :: nufactor = 1d0
!      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:) :: Fluxf_old
      real(8)                                   :: powlin !! for test compared with parcs
      real(8)                                   :: plevel0
      !!real(8)                                   :: hac !! for test compared with parcs
      LOGICAL(1)                                :: flag_iter_hex = .false. !! jr: for TH Feedback
      LOGICAL(1)                                :: flag_print_nu = .false.
      INTEGER(4)                                :: cr_opt = 1 ! 1: parcs solver, 2: rast-k solver
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: laybank 
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: lstroda 
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: laptr 
      LOGICAL(1), ALLOCATABLE, DIMENSION(:)     :: lstrodb 
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: lcrbptr 
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: lcrbtola 
      INTEGER(4)                                :: nstrod
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:) :: aphif 
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:) :: dXS_dCR_ADF_hex
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:) :: dXS_dCR_CDF_hex
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:)   :: ADF_hex
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:)   :: CDF_hex
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: reflratz
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: dhatd
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: dhatzd
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: dhat0
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: dhatz0
      REAL(8), ALLOCATABLE :: velo_tr (:,:,:) 
      REAL(8), ALLOCATABLE :: rvdeltf(:,:,:)
      LOGICAL(1)                                :: if_hex_cr = .false. ! maXs feedback, CR
      REAL(8)                                   :: Tot_FuelVol_HEX
      LOGICAL(1)                                :: flag_init_thhex = .true.
      LOGICAL(1)                                :: flag_exp=.true.

#endif    
!
! ============================================================================= !
! GetNode.f90
! ============================================================================= !
      INTEGER(4)                                :: nintot
      LOGICAL(1)                                :: if_hexgometry = .false.
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: chi
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: imap
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: imapsol ! for solver
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:) :: Fluxf
#ifdef tuan_fr
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:) :: Fluxf_old
#endif

#ifdef tuan_fr_tdep
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: tfluxf
#endif
      LOGICAL(1)                                :: if_hexbc = .false.
      REAL(8)                                   :: hexbc
      LOGICAL(1)                                :: if_mgxs=.false.
      !! === imap ===
      ! 1) nodevolume
      ! 2) maXs_
      ! 3) Flux
      ! 4) betap
      ! 5) beta_d_tot
      ! 6) rvdelt
      ! 7) src
      !! ============
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:) :: xschi, xschid
      INTEGER(4)                                :: iptype
      LOGICAL(1)                                :: flag_hex_geom = .false.
      INTEGER(4)                                :: lupscat
      INTEGER(4)                                :: ndivhs
      INTEGER(4)                                :: ntph
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: wtass
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: wtasssum
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: izpr
      INTEGER(4)                                :: nzp1
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: iasypr
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: iprcomp
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: xsid(:,:)
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: lfaptr_hex
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: latol_hex
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: ltola_hex
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: ltolfa_hex
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: ltolb
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: lfatol_hex
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: nodef_hex
      INTEGER(4)                                :: nfuelfa
      INTEGER(4)                                :: nfuel
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: nxfas, nxfae, nxas, nxae
      INTEGER(4)                                :: kfs, kfe
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: hx
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: hy
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: hz
      INTEGER(4)                                :: ndomx_hex, ndomy_hex, ndomz_hex
      INTEGER(4)                                :: ndomxy_hex
      INTEGER(4)                                :: ndom_hex
      INTEGER(4)                                :: idom_hex
      INTEGER(4)                                :: idomxy_hex
      INTEGER(4)                                :: idomx_hex
      INTEGER(4)                                :: idomy_hex
      INTEGER(4)                                :: idomz_hex
      INTEGER(4)                                :: nxyb
      INTEGER(4)                                :: nmzbi
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: ktokb
      INTEGER(4)                                :: nzb
      INTEGER(4)                                :: kbs
      INTEGER(4)                                :: kbe
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: alxr, alxrf
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: alzl, alzlf
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: reflratf
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: reflratzbf
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: reflratztf
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: reflratzf
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: alphaz
      REAL(8)   , ALLOCATABLE, DIMENSION(:)     :: alzr, alzrf
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:) :: sefve
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:) :: dsum
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:)   :: rhzbar2
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: igc
      INTEGER(4)                                :: ntsweep
#ifdef tuan_fr
    integer, allocatable ::  imapb(:), imapf(:)

#endif

      ! === subroutine orderhex === !
      INTEGER(4)                                :: isymang   = 360 ! Current: fixed version (whole core), Need to add option
      INTEGER(4)                                :: isolang   = 360 ! Current: fixed version (whole core), Need to add option
      INTEGER(4)                                :: isymtype  = 1
      INTEGER(4)                                :: isymmetry = 1
      INTEGER(4)                                :: nring
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: ltox, ltoy
! modified by Tuan
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: nodel_hex
!      INTEGER(4), POINTER, DIMENSION(:):: nodel_hex(:,:)      !(-3:nx+4,-1:ny+2) !radial node index
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: nxs, nxe, nys, nye
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: iasytype
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: nrnx_hex
      INTEGER(4)                                :: nassy
      INTEGER(4)                                :: nxsfc
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:)   :: pbdv ! boundary setting value from net J=alpha flux for all points
      INTEGER(4)                                :: nxpnt
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: neignd
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: neigpt
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: neigsfc
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: neigjin
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:)   :: wtdhat
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: neigsnd
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: neigz
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: neigsfcz
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: neigsndz
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:,:) :: imatid
      INTEGER(4)                                :: icoreys
      INTEGER(4)                                :: icoreye
      INTEGER(4)                                :: icoref
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: icorexs
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: icorexe
      INTEGER(4)                                :: icoreyps
      INTEGER(4)                                :: icoreype
      INTEGER(4)                                :: icorepf
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: icorexps
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: icorexpe
      INTEGER(4)                                :: icoreyss
      INTEGER(4)                                :: icoreyse
      INTEGER(4)                                :: icoresf
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: icorexss
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: icorexse
      INTEGER(4)                                :: icoreyfs
      INTEGER(4)                                :: icoreyfe
      INTEGER(4)                                :: icoreff
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: icorexfs
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: icorexfe
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: neignpt
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: neigppt
      LOGICAL, ALLOCATABLE, DIMENSION(:)        :: iffuel
      LOGICAL, ALLOCATABLE, DIMENSION(:)        :: iffuelhex
      INTEGER(4)                                :: nxp1, nyp1
      ! === subroutine orderhex === !

      ! === subroutine readlay === !
      INTEGER(4)                                :: icordxs
      INTEGER(4)                                :: icordxe
      INTEGER(4)                                :: icordys
      INTEGER(4)                                :: icordye
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: hex_map
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: hex_map_mid
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: iapoint
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: iasurf
      INTEGER(4)                                :: nxrow    ! = Nx_1N
      ! === subroutine readlay === !

      ! === subroutine scanlay === !
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: layh
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: layp
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: lays
      INTEGER(4)                                :: nat      ! number of assembly types
      INTEGER(4)                                :: nassyinp ! number of assemblies specified in the input
      INTEGER(4)                                :: nx_hex, ny_hex, nxy_hex
      INTEGER(4)                                :: nsurf
      INTEGER(4)                                :: nasyx
      INTEGER(4)                                :: nasyy
!      INTEGER(4)                            :: nassy
      INTEGER(4)                                :: nchan_hex
      INTEGER(4)                                :: npr
      ! === subroutine scanlay === !

      ! === subroutine hex_map_make === !
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)   :: hex_map_index ! = pattern in scarpkit.f90
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: ltox_hex
      INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: ltoy_hex
      INTEGER(4)                                :: ncorn_hex
      INTEGER(4)                                :: nxfc, nyfc
      ! === subroutine hex_map_make === !

! ============================================================================= !
! Mod_SolTPEN.f90
! ============================================================================= !
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)     :: ipntr
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:,:)   :: ineigcond
      INTEGER(4), ALLOCATABLE, DIMENSION(:)       :: ilubnd
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:)     :: iastopnt
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: cmat
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: dcmat
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: cmatd
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: dcmatd
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: xlufac
      LOGICAL(1)                                  :: extth=.false. ! default option, please check this option
      LOGICAL(1)                                  :: fdbk =.false. ! if you use THFB, please turn on this option
      LOGICAL(1)                                  :: ssinit=.false.
      INTEGER(4)                                  :: nodalcy=3

! ============================================================================= !
! Mod_TPENDrive.f90
! ============================================================================= !
!      REAL(8), DIMENSION(:, :), ALLOCATABLE ::  beta_d_Tot
      INTEGER(4)                                  :: ntnodal
      INTEGER(4)                                  :: ninitoutt
      REAL(8)                                     :: tlap
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: hflx
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: pflx
      INTEGER(4), ALLOCATABLE, DIMENSION(:)       :: mgb, mge ! if you needs multigroup calculation, please check those values
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: dfd
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: dfdz
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: dhat
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: dhatz
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: betaphis
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: betaphisz
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: betaphisd
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: betaphiszd
!      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: betap
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: fhflx
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: fohflx
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: hflxf
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: aflx
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: xmom
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: ymom
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: zmom1
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: zmom2
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: fcnto
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: focnto
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: cnto
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: fcntzo
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: focntzo
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: cntzo
#ifdef tuan_tr_test
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: srcf !! transient calculation, srcf => src
#else
      REAL(8)   , POINTER,     DIMENSION(:,:,:)   :: srcf !! transient calculation, srcf => src
#endif    
      REAL(8)   , ALLOCATABLE, DIMENSION(:)       :: pflxt
      REAL(8)                                     :: chlval(7)
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: xsadf
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: atleak
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: srcz
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: adft !! have relationship with control rod
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: adfb !! have relationship with control rod
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: xssf !! need to check
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:)   :: xstd !! need to check
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:,:,:) :: adfpt !! have relationship with control rod
      INTEGER(4), ALLOCATABLE, DIMENSION(:)       :: iscatib
      INTEGER(4), ALLOCATABLE, DIMENSION(:)       :: iscatie
      INTEGER(4), ALLOCATABLE, DIMENSION(:)       :: iscatob
      INTEGER(4), ALLOCATABLE, DIMENSION(:)       :: iscatoe
      INTEGER(4)                                  :: nouttot
!
! ============================================================================= !
! Mod_SolLS.f90
! ============================================================================= !
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:)     :: dumrv
      REAL(8)   , ALLOCATABLE, DIMENSION(:,:)     :: dumrs

! ============================================================================= !
! Mod_XSFB.f90
! ============================================================================= !
      REAL(8) :: ppmref = 1.0
      REAL(8) :: tmref  = 290
      REAL(8) :: dmref  = 0
      REAL(8) :: tfref  = 480
      ! =====================
      REAL(8), ALLOCATABLE ::D_3D_MG          (:,:,:)
      REAL(8), ALLOCATABLE ::maXS_tr_3D_MG    (:,:,:)
      REAL(8), ALLOCATABLE ::maXS_a_3D_MG     (:,:,:)
      REAL(8), ALLOCATABLE ::maXS_f_3D_MG     (:,:,:)
      REAL(8), ALLOCATABLE ::nu_maXS_f_3D_MG  (:,:,:)
      REAL(8), ALLOCATABLE ::kap_maXS_f_3D_MG (:,:,:)
      REAL(8), ALLOCATABLE ::maXS_s_3D_MG     (:,:,:)
      REAL(8), ALLOCATABLE ::maXS_r_3D_MG     (:,:,:)
      REAL(8), ALLOCATABLE ::maXS_chi_3D_MG   (:,:,:)
      REAL(8), ALLOCATABLE ::maXS_chid_3D_MG  (:,:,:)
      REAL(8), ALLOCATABLE ::maXS_scat_3D_MG  (:,:,:,:)

      END MODULE Inc_TPEN
#endif
