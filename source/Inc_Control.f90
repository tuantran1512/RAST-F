
      module Inc_Control
      implicit none

      real(8) :: EPS_Critical
      real(8) :: EPS_keff_St
      real(8) :: EPS_Linf_St
      real(8) :: EPS_Linf_Tr
      real(8) :: EPS_PowLev
      real(8) :: EPS_Flux
      real(8) :: EPS_TF
      real(8) :: EPS_TM
      real(8) :: EPS_Xe
      real(8) :: EPS_Sm

      real(8) :: EPS_Global
      real(8) :: EPS_Local
      real(8) :: EPS_Residual
      real(8) :: EPS_ERF
      real(8) :: EPS_DoppTF

      logical(1) :: Flag_Conv_PowLev
      logical(1) :: Flag_Conv_Critical

      integer :: Iout_Max
      integer :: Iin_Max
      integer :: I_CriSearch_Max
      integer :: I_XSFB_Max
      integer :: I_CRCH_Max
#ifdef K1LF
      integer :: I_BRCH_Max
#endif
      integer :: Period_TH
      integer :: Period_NonL
      integer :: Period_CMFD

      real(8) :: ini_keff
      real(8) :: ini_Flux
      real(8) :: ini_jOut
      real(8) :: ini_TL
      real(8) :: k_target = 1d0

      real(8) :: eigvd
      real(8) :: eigvs
      real(8) :: eshift
      real(8) :: eshift0
      real(8) :: reigv
      real(8) :: reigvd
      real(8) :: reigvs
      real(8) :: reigvsd
      real(8) :: domr

      real(8) :: errlinf
      real(8) :: errl2
      real(8) :: errl2d
      real(8) :: erreig
      real(8) :: rerrl2

      logical :: flagr2
      logical :: flagl2
      logical :: flaglinf
      logical :: flagerf
      logical :: flageig
      logical :: ifnodal
      logical :: nodal
      logical :: flagth

      logical :: iflsupd
      logical :: ifnlupd
      logical :: ifxsupd
      logical :: flagneut
      logical :: skip


#ifdef jr_vver
      real(8) :: epsl2=1d-5
#endif

#ifdef siarhei_plot
      character(100)    :: last_saved,current_sub
#endif
      integer(2)        :: dummy_filler = 0

      end module Inc_Control
