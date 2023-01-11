! ============================================================================
! module Inc_ExtSrc
! ============================================================================
!> @brief
!> variables for transient analysis with external source
!> @date
!> 2020.03.10
      module Inc_ExtSrc

      ! cooling time of burned fuel [sec]
      real(8) :: cool_time

      ! external source multiplier
      real(8) :: extsrc_multiplier

      ! flag for external source
      logical :: flag_ExtSrc = .false.

      ! decay constant [#/sec] (ENDF/B-VII.1 decay data)
      real(8), parameter :: &
      & dkcon_U34  = 8.95297E-14, &
      & dkcon_U35  = 3.12298E-17, &
      & dkcon_U36  = 9.38495E-16, &
      & dkcon_U38  = 4.91933E-18, &
      & dkcon_Np37 = 1.02517E-14, &
      & dkcon_Pu38 = 2.50622E-10, &
      & dkcon_Pu39 = 9.11636E-13, &
      & dkcon_Pu40 = 3.35003E-12, &
      & dkcon_Pu41 = 1.53811E-09, &
      & dkcon_Pu42 = 5.88475E-14, &
      & dkcon_Am41 = 5.08080E-11, &
      & dkcon_As42 = 1.55883E-10, &
      & dkcon_Am43 = 2.98230E-12, &
      & dkcon_Cm42 = 4.92361E-08, &
      & dkcon_Cm43 = 7.55311E-10, &
      & dkcon_Cm44 = 1.21367E-09


      ! sf decay branching ratio (ENDF/B-VII.1 decay data)
      real(8), parameter :: &
      & sf_br_U34  = 0.164E-10, &
      & sf_br_U35  = 0.700E-10, &
      & sf_br_U36  = 0.940E-09, &
      & sf_br_U38  = 0.545E-06, &
      & sf_br_Np37 = 0.210E-11, &
      & sf_br_Pu38 = 0.185E-08, &
      & sf_br_Pu39 = 0.300E-11, &
      & sf_br_Pu40 = 0.575E-07, &
      & sf_br_Pu41 = 0.240E-15, &
      & sf_br_Pu42 = 0.554E-05, &
      & sf_br_Am41 = 0.430E-11, &
      & sf_br_As42 = 0.470E-10, &
      & sf_br_Am43 = 0.370E-10, &
      & sf_br_Cm42 = 0.637E-07, &
      & sf_br_Cm43 = 0.530E-10, &
      & sf_br_Cm44 = 0.137E-05


      ! neutron yield per sf (ENDF/B-VII.1 decay data)
      real(8), parameter :: &
      & nu_sf_U34  = 1.810, &
      & nu_sf_U35  = 1.860, &
      & nu_sf_U36  = 1.910, &
      & nu_sf_U38  = 2.010, &
      & nu_sf_Np37 = 2.050, &
      & nu_sf_Pu38 = 2.220, &
      & nu_sf_Pu39 = 2.160, &
      & nu_sf_Pu40 = 2.160, &
      & nu_sf_Pu41 = 2.250, &
      & nu_sf_Pu42 = 2.150, &
      & nu_sf_Am41 = 2.270, &
      & nu_sf_As42 = 2.340, &
      & nu_sf_Am43 = 2.420, &
      & nu_sf_Cm42 = 2.520, &
      & nu_sf_Cm43 = 2.680, &
      & nu_sf_Cm44 = 2.690


      ! alpha decay branching ratio (ENDF/B-VII.1 decay data)
      real(8), parameter :: &
      & alpa_br_U34  = 1.00000E+00, &
      & alpa_br_U35  = 1.00000E+00, &
      & alpa_br_U36  = 1.00000E+00, &
      & alpa_br_U38  = 1.00000E+00, &
      & alpa_br_Np37 = 1.00000E+00, &
      & alpa_br_Pu38 = 1.00000E+00, &
      & alpa_br_Pu39 = 9.99537E-01, &
      & alpa_br_Pu40 = 1.00000E+00, &
      & alpa_br_Pu41 = 2.45000E-05, &
      & alpa_br_Pu42 = 9.99995E-01, &
      & alpa_br_Am41 = 1.00000E+00, &
      & alpa_br_As42 = 4.59000E-03, &
      & alpa_br_Am43 = 1.00000E+00, &
      & alpa_br_Cm42 = 1.00000E+00, &
      & alpa_br_Cm43 = 9.97100E-01, &
      & alpa_br_Cm44 = 9.99999E-01



      ! natural abundance (fraction) of O16, O17, O18
      real(8), parameter :: &
      & NA_O16  = 0.997620, &
      & NA_O17  = 0.000380, &
      & NA_O18  = 0.002000


      ! tau = (a,n) integral production cross section based on ORIGEN (a,n) library
      ! Ref: D.P. Griesheimer, et al., In-Line (a,n) Source Sampling Methodology for
      !      Monte Carlo Radiation Transport Simulations, M&C 2017
      real(8), parameter :: &
      & tau_O17_U34  = 3.81059E+08, &
      & tau_O17_U35  = 2.90138E+08, &
      & tau_O17_U36  = 3.12903E+08, &
      & tau_O17_U38  = 2.58142E+08, &
      & tau_O17_Np37 = 3.87748E+08, &
      & tau_O17_Pu38 = 6.47889E+08, &
      & tau_O17_Pu39 = 4.99916E+08, &
      & tau_O17_Pu40 = 5.02410E+08, &
      & tau_O17_Pu41 = 0.00000E+00, &
      & tau_O17_Pu42 = 4.12971E+08, &
      & tau_O17_Am41 = 6.49502E+08, &
      & tau_O17_As42 = 0.00000E+00, &
      & tau_O17_Am43 = 5.53842E+08, &
      & tau_O17_Cm42 = 9.47577E+08, &
      & tau_O17_Cm43 = 8.07766E+08, &
      & tau_O17_Cm44 = 7.91616E+08


      real(8), parameter :: &
      & tau_O18_U34  = 8.71868E+08, &
      & tau_O18_U35  = 5.42141E+08, &
      & tau_O18_U36  = 6.60296E+08, &
      & tau_O18_U38  = 4.37421E+08, &
      & tau_O18_Np37 = 8.84609E+08, &
      & tau_O18_Pu38 = 1.41473E+09, &
      & tau_O18_Pu39 = 1.11660E+09, &
      & tau_O18_Pu40 = 1.11991E+09, &
      & tau_O18_Pu41 = 0.00000E+00, &
      & tau_O18_Pu42 = 9.47786E+08, &
      & tau_O18_Am41 = 1.41812E+09, &
      & tau_O18_As42 = 0.00000E+00, &
      & tau_O18_Am43 = 1.21849E+09, &
      & tau_O18_Cm42 = 2.03362E+09, &
      & tau_O18_Cm43 = 1.74769E+09, &
      & tau_O18_Cm44 = 1.71435E+09


      real(8), parameter :: &
      & sqrtZ_O16  = 2.82843, &
      & sqrtZ_O17  = 2.82843, &
      & sqrtZ_O18  = 2.82843, &
      & sqrtZ_I35  = 7.28011, &
      & sqrtZ_Xe35 = 7.34847, &
      & sqrtZ_Nd47 = 7.74597, &
      & sqrtZ_Nd48 = 7.74597, &
      & sqrtZ_Nd49 = 7.74597, &
      & sqrtZ_Pm47 = 7.81025, &
      & sqrtZ_Ps48 = 7.81025, &
      & sqrtZ_Pm48 = 7.81025, &
      & sqrtZ_Pm49 = 7.81025, &
      & sqrtZ_Sm47 = 7.87401, &
      & sqrtZ_Sm48 = 7.87401, &
      & sqrtZ_Sm49 = 7.87401, &
      & sqrtZ_Gd52 = 8.00000, &
      & sqrtZ_Gd54 = 8.00000, &
      & sqrtZ_Gd55 = 8.00000, &
      & sqrtZ_Gd56 = 8.00000, &
      & sqrtZ_Gd57 = 8.00000, &
      & sqrtZ_Gd58 = 8.00000, &
      & sqrtZ_Gd60 = 8.00000, &
      & sqrtZ_U34  = 9.59166, &
      & sqrtZ_U35  = 9.59166, &
      & sqrtZ_U36  = 9.59166, &
      & sqrtZ_U37  = 9.59166, &
      & sqrtZ_U38  = 9.59166, &
      & sqrtZ_Np37 = 9.64365, &
      & sqrtZ_Np38 = 9.64365, &
      & sqrtZ_Np39 = 9.64365, &
      & sqrtZ_Pu38 = 9.69536, &
      & sqrtZ_Pu39 = 9.69536, &
      & sqrtZ_Pu40 = 9.69536, &
      & sqrtZ_Pu41 = 9.69536, &
      & sqrtZ_Pu42 = 9.69536, &
      & sqrtZ_Pu43 = 9.69536, &
      & sqrtZ_Am41 = 9.74679, &
      & sqrtZ_As42 = 9.74679, &
      & sqrtZ_Am42 = 9.74679, &
      & sqrtZ_Am43 = 9.74679, &
      & sqrtZ_Am44 = 9.74679, &
      & sqrtZ_Cm42 = 9.79796, &
      & sqrtZ_Cm43 = 9.79796, &
      & sqrtZ_Cm44 = 9.79796


      ! volume-integrated external source (lambda*N) [#/sec]
      real(8) :: tot_ExtSrc
      real(8) :: tot_ExtSrc_sf
      real(8) :: tot_ExtSrc_an
      real(8), dimension(:,:), allocatable :: ExtSrc   ! ( nxy, nz )
      real(8), dimension(:), allocatable :: ExtSrc_Z   ! ( nz )
      real(8), dimension(:), allocatable :: ExtSrc_XY  ! ( nxy )
      real(8), dimension(:), allocatable :: ExtSrc_XY_1N  ! ( nxy_1N )

      end module Inc_ExtSrc
