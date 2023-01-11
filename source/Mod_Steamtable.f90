#ifdef siarhei_delete

module mod_steamtable


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

implicit none

contains

!(IAPWS-IF97) from CTF  :HJ KIM
!Thermodynamic Properties of Water and Steam

!----------------------------------------------------------------------------
function get_water_density(flag,inp1,inp2) result(outp)
!  flag == 1
!       input 1 :pressure in MPa
!             2 :temperature in K
!       output  :water density g/cc

!  flag == 2
!       input 1 :pressure in MPa
!             2 :specific enthalpy in kJ/kg
!       output  :water temperature in K

    implicit none

    integer(4), intent(in)    ::  flag
    real(8),    intent(in)    ::  inp1
    real(8),    intent(in)    ::  inp2
    real(8)                   ::  outp
    real(8)                   ::  temp


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_water_density] in mod_steamtable'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    if(flag .eq. 1) then
        outp = water_PT2D(inp1, inp2)
    else if(flag .eq. 2) then
        temp = water_PH2T(inp1, inp2)
        outp = water_PT2D(inp1, temp)
    else
        ! stop simulation
        write(*,*)"Error in get_water_density: flag is invalid"
    endif

end function
!!----------------------------------------------------------------------------


!!----------------------------------------------------------------------------
function water_PH2T(P,H) result(T)
!input  P         :pressure in MPa
!       H         :specific enthalpy in kJ/kg
!output T         :water temperature in K

    implicit none

    real(8), intent(in)    ::  P
    real(8), intent(in)    ::  H
    real(8)                ::  T
    real(8), parameter     ::  a1 = -0.23872489924521D+03
    real(8), parameter     ::  a2 =  0.40421188637945D+03
    real(8), parameter     ::  a3 =  0.11349746881718D+03
    real(8), parameter     ::  a4 = -0.58457616048039D+01
    real(8), parameter     ::  a5 = -0.15285482413140D-03
    real(8), parameter     ::  a6 = -0.10866707695377D-05
    real(8), parameter     ::  a7 = -0.13391744872602D+02
    real(8), parameter     ::  a8 =  0.43211039183559D+02
    real(8), parameter     ::  a9 = -0.54010067170506D+02
    real(8), parameter     ::  a10=  0.30535892203916D+02
    real(8), parameter     ::  a11= -0.65964749423638D+01
    real(8), parameter     ::  a12=  0.93965400878363D-02
    real(8), parameter     ::  a13=  0.11573647505340D-06
    real(8), parameter     ::  a14= -0.25858641282073D-04
    real(8), parameter     ::  a15= -0.40644363084799D-08
    real(8), parameter     ::  a16=  0.66456186191635D-07
    real(8), parameter     ::  a17=  0.80670734103027D-10
    real(8), parameter     ::  a18= -0.93477771213947D-12
    real(8), parameter     ::  a19=  0.58265442020601D-14
    real(8), parameter     ::  a20= -0.15020185953503D-16
    real(8)                ::  hn
    real(8)                ::  h2
    real(8)                ::  h4
    real(8)                ::  h6
    real(8)                ::  h10


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [water_PH2T] in mod_steamtable'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

     hn = 1.0d0+h*4.d-04
     h2 = hn*hn
     h4 = h2*h2
     h6 = h4*h2
     h10= h6*h4

     T =                                                      &
        a1             + hn*(a2          + hn*(a3          +  &
        h4*(a4         + h10*h6*(a5      + h10*(a6         +  &
        p*(a13         + p*(a15          + p*(a17          +  &
        p*(a18         + p*(a19          + p*a20)))))))))) +  &
        p*(a7          + hn*(a8          + hn*(a9          +  &
        hn*(a10        + hn*(a11         + h6*(a12         +  &
        p*(a14         + p*a16)))))))

end function
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
function water_PT2D(P,T) result(D)
!input  P         :pressure in MPa
!       T         :temperature in K
!output rho_water :water density g/cc

    implicit none

    real(8), intent(in)  :: P
    real(8), intent(in)  :: T
    real(8)              :: D
    real(8), parameter   :: a8 =-0.79068310787773D-08
    real(8), parameter   :: a9 = 0.16949507886565D-07
    real(8), parameter   :: a10= 0.53021235478366D-06
    real(8), parameter   :: a11= 0.90824711621636D-06
    real(8), parameter   :: a12= 0.60983184277677D-06
    real(8), parameter   :: a13= 0.14752738052287D-08
    real(8), parameter   :: a14= 0.26348204437581D-07
    real(8), parameter   :: a15= 0.16753299313106D-07
    real(8), parameter   :: a16=-0.26614606756583D-08
    real(8), parameter   :: a17= 0.24649255061299D-09
    real(8), parameter   :: a18= 0.40593624756496D-19
    real(8), parameter   :: a19= 0.26535353478691D-08
    real(8), parameter   :: a20= 0.23680051381070D-09
    real(8), parameter   :: a21= 0.71369114266351D-13
    real(8), parameter   :: a22= 0.25045010666356D-09
    real(8), parameter   :: a23= 0.72784546444320D-10
    real(8), parameter   :: a24= 0.16017135514411D-16
    real(8), parameter   :: a25= 0.56562757086698D-10
    real(8), parameter   :: a26= 0.28443854062250D-12
    real(8), parameter   :: a27= 0.38920900760265D-13
    real(8), parameter   :: a28= 0.40317346616717D-21
    real(8), parameter   :: a29=-0.92975593753127D-23
    real(8), parameter   :: a30=-0.21323943802988D-25
    real(8), parameter   :: a31= 0.10007510864939D-25
    real(8), parameter   :: a32=-0.15777067572480D-26
    real(8), parameter   :: a33= 0.83571296309235D-28
    real(8)              :: tsich
    real(8)              :: targ
    real(8)              :: targ2
    real(8)              :: targ4
    real(8)              :: targ8
    real(8)              :: targinv
    real(8)              :: targinv2
    real(8)              :: targinv3
    real(8)              :: targinv4
    real(8)              :: targinv7
    real(8)              :: targinv8
    real(8)              :: parg
    real(8)              :: parg2
    real(8)              :: parg3
    real(8)              :: parg6
    real(8)              :: pt
    real(8)              :: vpt1n


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [water_PT2D] in mod_steamtable'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

     tsich=1386.d0/1.222d0
     if ((T .GT. 0.D0) .AND. (DABS(T-TSICH) .GT. 1.D-12)) then
         targ = 1386.d0/t-1.222d0
         targ2 = targ*targ
         targ4 = targ2*targ2
         targ8 = targ4*targ4
         targinv = 1.d0/targ
         targinv2 = targinv*targinv
         targinv3 = targinv2*targinv
         targinv4 = targinv3*targinv
         targinv7 = targinv4*targinv3
         targinv8 = targinv7*targinv

         parg = 7.1d0-p*0.60496067755596D-01
         parg2 = parg*parg
         parg3 = parg2*parg
         parg6 = parg3*parg3

         pt = parg*targinv

         vpt1n=((((((((((((((pt*a33 + a32)*pt + a31)*pt + a30)*           &
              parg6*targinv7 + a29)*parg2*targinv2 + a28)*                &
              parg6*parg6*parg*targinv8*targinv8*targinv2 + a26)*         &
              targinv*targinv4 + a27)*parg3*targinv + targinv3*a25)*      &
              parg + a22)*targinv3 + a23)*targinv2 + targ8*targ2*a24)*    &
              parg + targ4*targ2*a21 + a20 + targinv4*a19)*parg +         &
              (targ8*targ8*a18 + targ2*a17 + a16)*targ + a15 +            &
              targinv3*a14)*parg + (targ2*a13 + a12)*targ + a11 +         &
              (targinv2*a8 + a9)*targinv7 + targinv*a10)*t

         D = 1/vpt1n
         D = D * 1D-3
     else
        !Stop simulation
        write(*,*) "Out of Range"
     endif
endfunction
!!----------------------------------------------------------------------------

!!----------------------------------------------------------------------------
function water_PT2H(P,T) result(H)
!input  P         :pressure in MPa
!       T         :temperature in K
!output H         :water enthalpy kJ/kg

    implicit none

    real(8), intent(in)   :: P
    real(8), intent(in)   :: T
    real(8)               :: H
    real(8), parameter    :: a01=-0.18720692775139D+03
    real(8), parameter    :: a02= 0.54083364671138D+03
    real(8), parameter    :: a04= 0.21656306556572D+04
    real(8), parameter    :: a05=-0.12255145525730D+04
    real(8), parameter    :: a06= 0.30266937911228D+03
    real(8), parameter    :: a07=-0.42516429081128D+02
    real(8), parameter    :: a08= 0.25975485679232D+01
    real(8), parameter    :: a09=-0.16303507737913D+01
    real(8), parameter    :: a10= 0.27182613947704D+01
    real(8), parameter    :: a11= 0.12147472571259D+02
    real(8), parameter    :: a13=-0.13971601220485D+02
    real(8), parameter    :: a14=-0.10139813560979D+00
    real(8), parameter    :: a15= 0.90547896843533D+00
    real(8), parameter    :: a17= 0.30487803863262D-01
    real(8), parameter    :: a18=-0.84709309503347D-02
    real(8), parameter    :: a19=-0.79051996435264D-11
    real(8), parameter    :: a20= 0.81058711826910D-01
    real(8), parameter    :: a22=-0.32702156038568D-05
    real(8), parameter    :: a23= 0.71724465059049D-02
    real(8), parameter    :: a24= 0.83376808703816D-03
    real(8), parameter    :: a25=-0.91740466143441D-09
    real(8), parameter    :: a26= 0.20734169140086D-02
    real(8), parameter    :: a27= 0.89603964175208D-05
    real(8), parameter    :: a28= 0.66877530790508D-06
    real(8), parameter    :: a29= 0.12755771455453D-13
    real(8), parameter    :: a30=-0.28710377452428D-15
    real(8), parameter    :: a31=-0.64016099916299D-18
    real(8), parameter    :: a32= 0.29806124175367D-18
    real(8), parameter    :: a33=-0.46640228230284D-19
    real(8), parameter    :: a34= 0.24531669269267D-20
    real(8)               :: tsich
    real(8)               :: targ
    real(8)               :: targ2
    real(8)               :: targ4
    real(8)               :: targ8
    real(8)               :: targinv
    real(8)               :: targinv2
    real(8)               :: targinv4
    real(8)               :: targinv8
    real(8)               :: targinv9
    real(8)               :: parg
    real(8)               :: parg2
    real(8)               :: parg4
    real(8)               :: parg8
    real(8)               :: pt
    real(8)               :: hpt1n


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [water_PT2H] in mod_steamtable'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    tsich=1386.d0/1.222d0
    if ((T .GT. 0.D0) .AND. (DABS(T-TSICH) .GT. 1.D-12)) then
        targ = 1386.d0/t-1.222d0
        targ2 = targ*targ
        targ4 = targ2*targ2
        targ8 = targ4*targ4
        targinv = 1.d0/targ
        targinv2 = targinv*targinv
        targinv4 = targinv2*targinv2
        targinv8 = targinv4*targinv4
        targinv9 = targinv8*targinv

        parg = 7.1d0-p*0.60496067755596D-01
        parg2 = parg*parg
        parg4 = parg2*parg2
        parg8 = parg4*parg4

        pt = parg*targinv

        hpt1n = ((((((((((((pt*a34 + a33)*pt + a32)*pt + a31)*                &
                parg4*parg2*targinv9 + targinv2*a30)*parg2 +                  &
                a29)*parg8*parg4*parg*targinv8*targinv8*targinv8 +            &
                targinv*a28)*targinv2 + targinv8*a27)*parg2*parg*             &
                targinv4 + targinv9*a26)*parg + targinv2*(targinv*            &
                a24 + targinv4*a23) + targ8*targ*a25)*parg + targ4*targ*      &
                a22 + targinv4*targinv*a20)*parg + targ8*targ8*a19 +          &
                targ2*a18 + targinv4*a15 + a17)*parg + targinv8*              &
                (targinv2*a09 + a10) + targ2*a14 + targinv2*a11 + a13)*parg + &
                (targ2*a08 + targ*a07 + a06)*targ2 +                          &
                targinv2*(targinv*a01 + a02) + targ*a05 + a04
        H = hpt1n
    else
       !Stop simulation
       H = 0d0
       write(*,*) "Out of Range"
    endif

endfunction
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
function adjust_rhol_for_boron(rhol,boron_conc) result(rho_adj)

! input   rhol       :water density in g/cc
!         boron_conc :boron concentration in ppm
!
! output  rho_adj    :adjusted water density with boron in g/cc
    implicit none

    real(8), intent(in) :: rhol       ! g/cc
    real(8), intent(in) :: boron_conc ! ppm
    real(8)             :: rho_adj    ! g/cc
    real(8), parameter :: rho_boron=1.435D0  ! g/cc ! Density of boric acid solute H3BO3
    real(8), parameter :: cbm_tol=0.0d0      ! boron tolerance for correction [ppm]

    real(8) :: v, boric_acid_ppm


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [adjust_rhol_for_boron] in mod_steamtable'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    if (boron_conc>cbm_tol) then ! Adjust density to account for presence of boron
       ! Convert boron ppm to boric acid ppm
       boric_acid_ppm = boron_conc/0.1746
       v = (1.0-boric_acid_ppm/1.0E6)/rhol+boric_acid_ppm/1.0E6/rho_boron
       rho_adj=1.0/v
    else
       rho_adj=rhol
    end if

end function adjust_rhol_for_boron
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
function water_PT2Cp(P,T) result(Cp)
!input   P     : pressure in MPa.
!        T     : temperature in K.
!output  Cp    : isochoric heat capacity in kJ/(kg * K)
!                isochoric (process in which the pressure stays constant)

implicit none

    real(8), intent(in)   ::  P
    real(8), intent(in)   ::  T
    real(8)               ::  Cp
    real(8), parameter    ::  e01=  0.14632971213167D0*6.D0
    real(8), parameter    ::  e02= -0.84548187169114D0*2.D0
    real(8), parameter    ::  e05= -0.95791963387872D0*2.D0
    real(8), parameter    ::  e06=  0.15772038513228D0*6.D0
    real(8), parameter    ::  e07= -0.16616417199501D-01*12.D0
    real(8), parameter    ::  e08=  0.81214629983568D-03*20.D0
    real(8), parameter    ::  e09=  0.28319080123804D-03*90.D0
    real(8), parameter    ::  e10= -0.60706301565874D-03*56.D0
    real(8), parameter    ::  e11= -0.18990068218419D-01*2.D0
    real(8), parameter    ::  e14= -0.52838357969930D-04*6.D0
    real(8), parameter    ::  e15= -0.47184321073267D-03*12.D0
    real(8), parameter    ::  e18= -0.44141845330846D-05*6.D0
    real(8), parameter    ::  e19= -0.72694996297594D-15*272.D0
    real(8), parameter    ::  e20= -0.31679644845054D-04*20.D0
    real(8), parameter    ::  e22= -0.85205128120103D-09*30.D0
    real(8), parameter    ::  e23= -0.22425281908000D-05*30.D0
    real(8), parameter    ::  e24= -0.65171222895601D-06*6.D0
    real(8), parameter    ::  e25= -0.14341729937924D-12*90.D0
    real(8), parameter    ::  e26= -0.40516996860117D-06*72.D0
    real(8), parameter    ::  e27= -0.12734301741641D-08*132.D0
    real(8), parameter    ::  e28= -0.17424871230634D-09*42.D0
    real(8), parameter    ::  e29= -0.68762131295531D-18*870.D0
    real(8), parameter    ::  e30=  0.14478307828521D-19*992.D0
    real(8), parameter    ::  e31=  0.26335781662795D-22*1482.D0
    real(8), parameter    ::  e32= -0.11947622640071D-22*1560.D0
    real(8), parameter    ::  e33=  0.18228094581404D-23*1640.D0
    real(8), parameter    ::  e34= -0.93537087292458D-25*1722.D0
    real(8), parameter    ::  c15= -0.47184321073267D-03*2.D0
    real(8), parameter    ::  c16= -0.30001780793026D-03*2.D0
    real(8), parameter    ::  c17=  0.47661393906987D-04*2.D0
    real(8), parameter    ::  c18= -0.44141845330846D-05*2.D0
    real(8), parameter    ::  c19= -0.72694996297594D-15*2.D0
    real(8), parameter    ::  c20= -0.31679644845054D-04*6.D0
    real(8), parameter    ::  c21= -0.28270797985312D-05*6.D0
    real(8), parameter    ::  c22= -0.85205128120103D-09*6.D0
    real(8), parameter    ::  c23= -0.22425281908000D-05*12.D0
    real(8), parameter    ::  c24= -0.65171222895601D-06*12.D0
    real(8), parameter    ::  c25= -0.14341729937924D-12*12.D0
    real(8), parameter    ::  c26= -0.40516996860117D-06*20.D0
    real(8), parameter    ::  c27= -0.12734301741641D-08*56.D0
    real(8), parameter    ::  c28= -0.17424871230634D-09*56.D0
    real(8), parameter    ::  c29= -0.68762131295531D-18*420.D0
    real(8), parameter    ::  c30=  0.14478307828521D-19*506.D0
    real(8), parameter    ::  c31=  0.26335781662795D-22*812.D0
    real(8), parameter    ::  c32= -0.11947622640071D-22*870.D0
    real(8), parameter    ::  c33=  0.18228094581404D-23*930.D0
    real(8), parameter    ::  c34= -0.93537087292458D-25*992.D0
    real(8), parameter    ::  b09=  0.28319080123804D-03
    real(8), parameter    ::  b10= -0.60706301565874D-03
    real(8), parameter    ::  b11= -0.18990068218419D-01
    real(8), parameter    ::  b12= -0.32529748770505D-01
    real(8), parameter    ::  b13= -0.21841717175414D-01
    real(8), parameter    ::  b14= -0.52838357969930D-04
    real(8), parameter    ::  b15= -0.47184321073267D-03*2.D0
    real(8), parameter    ::  b16= -0.30001780793026D-03*2.D0
    real(8), parameter    ::  b17=  0.47661393906987D-04*2.D0
    real(8), parameter    ::  b18= -0.44141845330846D-05*2.D0
    real(8), parameter    ::  b19= -0.72694996297594D-15*2.D0
    real(8), parameter    ::  b20= -0.31679644845054D-04*3.D0
    real(8), parameter    ::  b21= -0.28270797985312D-05*3.D0
    real(8), parameter    ::  b22= -0.85205128120103D-09*3.D0
    real(8), parameter    ::  b23= -0.22425281908000D-05*4.D0
    real(8), parameter    ::  b24= -0.65171222895601D-06*4.D0
    real(8), parameter    ::  b25= -0.14341729937924D-12*4.D0
    real(8), parameter    ::  b26= -0.40516996860117D-06*5.D0
    real(8), parameter    ::  b27= -0.12734301741641D-08*8.D0
    real(8), parameter    ::  b28= -0.17424871230634D-09*8.D0
    real(8), parameter    ::  b29= -0.68762131295531D-18*21.D0
    real(8), parameter    ::  b30=  0.14478307828521D-19*23.D0
    real(8), parameter    ::  b31=  0.26335781662795D-22*29.D0
    real(8), parameter    ::  b32= -0.11947622640071D-22*30.D0
    real(8), parameter    ::  b33=  0.18228094581404D-23*31.D0
    real(8), parameter    ::  b34= -0.93537087292458D-25*32.D0
    real(8), parameter    ::  f09= -0.28319080123804D-03*9.D0
    real(8), parameter    ::  f10=  0.60706301565874D-03*7.D0
    real(8), parameter    ::  f11=  0.18990068218419D-01*1.D0
    real(8), parameter    ::  f13= -0.21841717175414D-01
    real(8), parameter    ::  f14= -0.52838357969930D-04*3.D0
    real(8), parameter    ::  f15=  0.47184321073267D-03*6.D0
    real(8), parameter    ::  f17=  0.47661393906987D-04*2.D0
    real(8), parameter    ::  f18= -0.44141845330846D-05*6.D0
    real(8), parameter    ::  f19= -0.72694996297594D-15*34.D0
    real(8), parameter    ::  f20=  0.31679644845054D-04*12.D0
    real(8), parameter    ::  f22= -0.85205128120103D-09*18.D0
    real(8), parameter    ::  f23=  0.22425281908000D-05*20.D0
    real(8), parameter    ::  f24=  0.65171222895601D-06*8.D0
    real(8), parameter    ::  f25= -0.14341729937924D-12*40.D0
    real(8), parameter    ::  f26=  0.40516996860117D-06*40.D0
    real(8), parameter    ::  f27=  0.12734301741641D-08*88.D0
    real(8), parameter    ::  f28=  0.17424871230634D-09*48.D0
    real(8), parameter    ::  f29=  0.68762131295531D-18*609.D0
    real(8), parameter    ::  f30= -0.14478307828521D-19*713.D0
    real(8), parameter    ::  f31= -0.26335781662795D-22*1102.D0
    real(8), parameter    ::  f32=  0.11947622640071D-22*1170.D0
    real(8), parameter    ::  f33= -0.18228094581404D-23*1240.D0
    real(8), parameter    ::  f34=  0.93537087292458D-25*1312.D0
    real(8), parameter    ::  tn= 1386.D0
    real(8), parameter    ::  pn= 16.53D0
    real(8), parameter    ::  pnq= 1.D0/pn
    real(8), parameter    ::  r= 0.461526D0
    real(8)               ::  tsich
    real(8)               ::  tau
    real(8)               ::  taun
    real(8)               ::  taun2
    real(8)               ::  taun3
    real(8)               ::  taun6
    real(8)               ::  taun12
    real(8)               ::  tauninv
    real(8)               ::  tauninv2
    real(8)               ::  tauninv3
    real(8)               ::  tauninv7
    real(8)               ::  pin
    real(8)               ::  pin2
    real(8)               ::  pin3
    real(8)               ::  pin6
    real(8)               ::  gatt1
    real(8)               ::  gapp1
    real(8)               ::  gap1
    real(8)               ::  gapt1


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [water_PT2Cp] in mod_steamtable'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    tsich=1386.d0/1.222d0
    if ((T .GT. 0.D0) .AND. (DABS(T-TSICH) .GT. 1.D-12)) then
        tau= tn/t
        taun= tau - 1.222d0
        taun2= taun*taun
        taun3= taun2*taun
        taun6= taun3*taun3
        taun12= taun6*taun6
        tauninv= 1.d0/taun
        tauninv2= tauninv*tauninv
        tauninv3= tauninv2*tauninv
        tauninv7= tauninv3*tauninv3*tauninv

        pin= 7.1d0 - p*pnq
        pin2= pin*pin
        pin3= pin2*pin
        pin6= pin3*pin3

        gatt1= tauninv2*(tauninv2*(e01 + taun*(e02 + taun3*(e05 +  &
               taun*(e06 + taun*(e07 + taun*e08)))) +              &
               pin*  (tauninv7*(e09 + taun2*(e10 + taun6*(e11 +    &
                     taun2*taun2*e14))))) +                        &
               pin2* tauninv3*(e15 + taun6*(e18 +                  &
                     taun12*taun2*e19) +                           &
               pin*  tauninv*(e20 + taun6*taun3*taun*e22 +         &
               pin*  tauninv*(e23 + taun3*(e24 + taun12*e25) +     &
               pin*  tauninv3*(e26 +                               &
               pin3*(tauninv3*(e27 + taun3*taun2*e28) +            &
               pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(e29 +     &
               pin2* tauninv2*(e30 +                               &
               pin6* tauninv7*(e31 +                               &
               pin*  tauninv*(e32 +                                &
               pin*  tauninv*(e33 +                                &
               pin*  tauninv*e34))) ))) ))) ))

        gapp1= tauninv3*(c15 + taun3*(c16 + taun*(c17 +            &
               taun2*(c18 + taun12*taun2*c19))) +                  &
               pin*  tauninv*(c20 + taun3*taun*(c21 + taun6*c22) + &
               pin*  tauninv*(c23 + taun3*(c24 + taun12*c25) +     &
               pin*  tauninv3*(c26 +                               &
               pin3* (tauninv3*(c27 + taun3*taun2*c28) +           &
               pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(c29 +     &
               pin2* tauninv2*(c30 +                               &
               pin6* tauninv7*(c31 +                               &
               pin*  tauninv*(c32 +                                &
               pin*  tauninv*(c33 +                                &
               pin*  tauninv*c34))) ))) ))) )

        gap1= -(tauninv7*(tauninv2*b09 + b10 + taun6*(b11 +        &
              taun*(b12 + taun*(b13 + taun2*b14)))) +              &
              pin*  tauninv3*(b15 + taun3*(b16 + taun*(b17 +       &
                    taun2*(b18 + taun12*taun2*b19))) +             &
              pin*  tauninv*(b20 + taun3*taun*(b21 + taun6*b22) +  &
              pin*  tauninv*(b23 + taun3*(b24 + taun12*b25) +      &
              pin*  tauninv3*(b26 +                                &
              pin3* (tauninv3*(b27 + taun3*taun2*b28) +            &
              pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(b29 +      &
              pin2* tauninv2*(b30 +                                &
              pin6* tauninv7*(b31 +                                &
              pin*  tauninv*(b32 +                                 &
              pin*  tauninv*(b33 +                                 &
              pin*  tauninv*b34))) ))) ))) ))

        gapt1= -(tauninv*(tauninv7*(tauninv2*f09 + f10 + taun6*(f11 + &
               taun2*(f13 + taun2*f14))) +                            &
               pin*  tauninv3*(f15 + taun3*taun*(f17 +                &
                     taun2*(f18 + taun12*taun2*f19)) +                &
               pin*  tauninv*(f20 + taun6*taun3*taun*f22 +            &
               pin*  tauninv*(f23 + taun3*(f24 + taun12*f25) +        &
               pin*  tauninv3*(f26 +                                  &
               pin3*(tauninv3*(f27 + taun3*taun2*f28) +               &
               pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(f29 +        &
               pin2* tauninv2*(f30 +                                  &
               pin6* tauninv7*(f31 +                                  &
               pin*  tauninv*(f32 +                                   &
               pin*  tauninv*(f33 +                                   &
               pin*  tauninv*f34))) ))) ))) )))

        Cp =   r*(-tau*tau*gatt1+(gap1-tau*gapt1)*(gap1-tau*gapt1)/gapp1)
    else
       !stop simulation
       write(*,*) "Out of Range in function PT2Cp"
       Cp =0.D0
    endif

endfunction
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
FUNCTION water_PT2Cv(P,T) result(Cv)
!input   P     : pressure in MPa.
!        T     : temperature in K.
!output  Cv    : isobaric heat capacity in kJ/(kg * K)
!                isobaric (process in which the volume stays constant)
    implicit none

    real(8), intent(in)    ::  P
    real(8), intent(in)    ::  T
    real(8)                ::  Cv
    real(8), parameter     ::  e01=  0.14632971213167D0*6.D0
    real(8), parameter     ::  e02= -0.84548187169114D0*2.D0
    real(8), parameter     ::  e05= -0.95791963387872D0*2.D0
    real(8), parameter     ::  e06=  0.15772038513228D0*6.D0
    real(8), parameter     ::  e07= -0.16616417199501D-01*12.D0
    real(8), parameter     ::  e08=  0.81214629983568D-03*20.D0
    real(8), parameter     ::  e09=  0.28319080123804D-03*90.D0
    real(8), parameter     ::  e10= -0.60706301565874D-03*56.D0
    real(8), parameter     ::  e11= -0.18990068218419D-01*2.D0
    real(8), parameter     ::  e14= -0.52838357969930D-04*6.D0
    real(8), parameter     ::  e15= -0.47184321073267D-03*12.D0
    real(8), parameter     ::  e18= -0.44141845330846D-05*6.D0
    real(8), parameter     ::  e19= -0.72694996297594D-15*272.D0
    real(8), parameter     ::  e20= -0.31679644845054D-04*20.D0
    real(8), parameter     ::  e22= -0.85205128120103D-09*30.D0
    real(8), parameter     ::  e23= -0.22425281908000D-05*30.D0
    real(8), parameter     ::  e24= -0.65171222895601D-06*6.D0
    real(8), parameter     ::  e25= -0.14341729937924D-12*90.D0
    real(8), parameter     ::  e26= -0.40516996860117D-06*72.D0
    real(8), parameter     ::  e27= -0.12734301741641D-08*132.D0
    real(8), parameter     ::  e28= -0.17424871230634D-09*42.D0
    real(8), parameter     ::  e29= -0.68762131295531D-18*870.D0
    real(8), parameter     ::  e30=  0.14478307828521D-19*992.D0
    real(8), parameter     ::  e31=  0.26335781662795D-22*1482.D0
    real(8), parameter     ::  e32= -0.11947622640071D-22*1560.D0
    real(8), parameter     ::  e33=  0.18228094581404D-23*1640.D0
    real(8), parameter     ::  e34= -0.93537087292458D-25*1722.D0
    real(8), parameter     ::  tn= 1386.0D0
    real(8), parameter     ::  pn= 16.53D0
    real(8), parameter     ::  pnq= 1.D0/pn
    real(8), parameter     ::  r= 0.461526D0
    real(8)                ::  tsich
    real(8)                ::  tau
    real(8)                ::  taun
    real(8)                ::  taun2
    real(8)                ::  taun3
    real(8)                ::  taun6
    real(8)                ::  taun12
    real(8)                ::  tauninv
    real(8)                ::  tauninv2
    real(8)                ::  tauninv3
    real(8)                ::  tauninv7
    real(8)                ::  pin
    real(8)                ::  pin2
    real(8)                ::  pin3
    real(8)                ::  pin6


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [water_PT2Cv] in mod_steamtable'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    tsich=1386.d0/1.222d0
    if ((T .GT. 0.D0) .AND. (DABS(T-TSICH) .GT. 1.D-12)) then
        tau= tn/t
        taun= tau - 1.222D0
        taun2= taun*taun
        taun3= taun2*taun
        taun6= taun3*taun3
        taun12= taun6*taun6
        tauninv= 1.d0/taun
        tauninv2= tauninv*tauninv
        tauninv3= tauninv2*tauninv
        tauninv7= tauninv3*tauninv3*tauninv

        pin= 7.1d0 - p*pnq
        pin2= pin*pin
        pin3= pin2*pin
        pin6= pin3*pin3

        Cv =   (tauninv2*(tauninv2*(e01 + taun*(e02 + taun3*(e05 +  &
               taun*(e06 + taun*(e07 + taun*e08)))) +               &
               pin*  (tauninv7*(e09 + taun2*(e10 + taun6*(e11 +     &
                     taun2*taun2*e14))))) +                         &
               pin2* tauninv3*(e15 + taun6*(e18 +                   &
                     taun12*taun2*e19) +                            &
               pin*  tauninv*(e20 + taun6*taun3*taun*e22 +          &
               pin*  tauninv*(e23 + taun3*(e24 + taun12*e25) +      &
               pin*  tauninv3*(e26 +                                &
               pin3*(tauninv3*(e27 + taun3*taun2*e28) +             &
               pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(e29 +      &
               pin2* tauninv2*(e30 +                                &
               pin6* tauninv7*(e31 +                                &
               pin*  tauninv*(e32 +                                 &
               pin*  tauninv*(e33 +                                 &
               pin*  tauninv*e34))))))))))))*(-tau*tau*r)
    else
       !stop simulation
       write(*,*) "Out of range in PT2Cv"
       Cv =0.D0
    endif

endfunction
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------

function water_PT2DVDP(P,T) result(DVDP)
!P   is pressure in MPa
!T   is temperature in K
!DVDP is in m^3/(kg*MPa)

    implicit none

    real(8), intent(in)   ::  P
    real(8), intent(in)   ::  T
    real(8)               ::  DVDP
    real(8), parameter    ::  c15= -0.47184321073267D-03*2.D0
    real(8), parameter    ::  c16= -0.30001780793026D-03*2.D0
    real(8), parameter    ::  c17=  0.47661393906987D-04*2.D0
    real(8), parameter    ::  c18= -0.44141845330846D-05*2.D0
    real(8), parameter    ::  c19= -0.72694996297594D-15*2.D0
    real(8), parameter    ::  c20= -0.31679644845054D-04*6.D0
    real(8), parameter    ::  c21= -0.28270797985312D-05*6.D0
    real(8), parameter    ::  c22= -0.85205128120103D-09*6.D0
    real(8), parameter    ::  c23= -0.22425281908000D-05*12.D0
    real(8), parameter    ::  c24= -0.65171222895601D-06*12.D0
    real(8), parameter    ::  c25= -0.14341729937924D-12*12.D0
    real(8), parameter    ::  c26= -0.40516996860117D-06*20.D0
    real(8), parameter    ::  c27= -0.12734301741641D-08*56.D0
    real(8), parameter    ::  c28= -0.17424871230634D-09*56.D0
    real(8), parameter    ::  c29= -0.68762131295531D-18*420.D0
    real(8), parameter    ::  c30=  0.14478307828521D-19*506.D0
    real(8), parameter    ::  c31=  0.26335781662795D-22*812.D0
    real(8), parameter    ::  c32= -0.11947622640071D-22*870.D0
    real(8), parameter    ::  c33=  0.18228094581404D-23*930.D0
    real(8), parameter    ::  c34= -0.93537087292458D-25*992.D0
    real(8), parameter    ::  a06=-0.33032641670203D-04*2.D0
    real(8), parameter    ::  a07=-0.18948987516315D-03*2.D0
    real(8), parameter    ::  a08=-0.39392777243355D-02*2.D0
    real(8), parameter    ::  a09=-0.43797295650573D-01*2.D0
    real(8), parameter    ::  a10=-0.26674547914087D-04*2.D0
    real(8), parameter    ::  a11= 0.20481737692309D-07*6.D0
    real(8), parameter    ::  a12= 0.43870667284435D-06*6.D0
    real(8), parameter    ::  a13=-0.32277677238570D-04*6.D0
    real(8), parameter    ::  a14=-0.15033924542148D-02*6.D0
    real(8), parameter    ::  a15=-0.40668253562649D-01*6.D0
    real(8), parameter    ::  a16=-0.78847309559367D-09*12.D0
    real(8), parameter    ::  a17= 0.12790717852285D-07*12.D0
    real(8), parameter    ::  a18= 0.48225372718507D-06*12.D0
    real(8), parameter    ::  a19= 0.22922076337661D-05*20.D0
    real(8), parameter    ::  a20=-0.16714766451061D-10*30.D0
    real(8), parameter    ::  a21=-0.21171472321355D-02*30.D0
    real(8), parameter    ::  a22=-0.23895741934104D02*30.D0
    real(8), parameter    ::  a23=-0.59059564324270D-17*42.D0
    real(8), parameter    ::  a24=-0.12621808899101D-05*42.D0
    real(8), parameter    ::  a25=-0.38946842435739D-01*42.D0
    real(8), parameter    ::  a26= 0.11256211360459D-10*56.D0
    real(8), parameter    ::  a27=-0.82311340897998D01*56.D0
    real(8), parameter    ::  a28= 0.19809712802088D-07*72.D0
    real(8), parameter    ::  a29= 0.10406965210174D-18*90.D0
    real(8), parameter    ::  a30=-0.10234747095929D-12*90.D0
    real(8), parameter    ::  a31=-0.10018179379511D-08*90.D0
    real(8), parameter    ::  a32=-0.80882908646985D-10*240.D0
    real(8), parameter    ::  a33= 0.10693031879409D0*240.D0
    real(8), parameter    ::  a34=-0.33662250574171D0*306.D0
    real(8), parameter    ::  a35= 0.89185845355421D-24*380.D0
    real(8), parameter    ::  a36= 0.30629316876232D-12*380.D0
    real(8), parameter    ::  a37=-0.42002467698208D-05*380.D0
    real(8), parameter    ::  a38=-0.59056029685639D-25*420.D0
    real(8), parameter    ::  a39= 0.37826947613457D-05*462.D0
    real(8), parameter    ::  a40=-0.12768608934681D-14*506.D0
    real(8), parameter    ::  a41= 0.73087610595061D-28*552.D0
    real(8), parameter    ::  a42= 0.55414715350778D-16*552.D0
    real(8), parameter    ::  a43=-0.94369707241210D-06*552.D0
    real(8), parameter    ::  r01=  0.10658070028513D01
    real(8), parameter    ::  r09= -0.12654315477714D01
    real(8), parameter    ::  r10= -0.11524407806681D01
    real(8), parameter    ::  r11=  0.88521043984318D0
    real(8), parameter    ::  r12= -0.64207765181607D0
    real(8), parameter    ::  r13=  0.38493460186671D0*2.D0
    real(8), parameter    ::  r14= -0.85214708824206D0*2.D0
    real(8), parameter    ::  r15=  0.48972281541877D01*2.D0
    real(8), parameter    ::  r16= -0.30502617256965D01*2.D0
    real(8), parameter    ::  r17=  0.39420536879154D-01*2.D0
    real(8), parameter    ::  r18=  0.12558408424308D0*2.D0
    real(8), parameter    ::  r19= -0.27999329698710D0*3.D0
    real(8), parameter    ::  r20=  0.13899799569460D01*3.D0
    real(8), parameter    ::  r21= -0.20189915023570D01*3.D0
    real(8), parameter    ::  r22= -0.82147637173963D-02*3.D0
    real(8), parameter    ::  r23= -0.47596035734923D0*3.D0
    real(8), parameter    ::  r24=  0.43984074473500D-01*4.D0
    real(8), parameter    ::  r25= -0.44476435428739D0*4.D0
    real(8), parameter    ::  r26=  0.90572070719733D0*4.D0
    real(8), parameter    ::  r27=  0.70522450087967D0*4.D0
    real(8), parameter    ::  r28=  0.10770512626332D0*5.D0
    real(8), parameter    ::  r29= -0.32913623258954D0*5.D0
    real(8), parameter    ::  r30= -0.50871062041158D0*5.D0
    real(8), parameter    ::  r31= -0.22175400873096D-01*6.D0
    real(8), parameter    ::  r32=  0.94260751665092D-01*6.D0
    real(8), parameter    ::  r33=  0.16436278447961D0*6.D0
    real(8), parameter    ::  r34= -0.13503372241348D-01*7.D0
    real(8), parameter    ::  r35= -0.14834345352472D-01*8.D0
    real(8), parameter    ::  r36=  0.57922953628084D-03*9.D0
    real(8), parameter    ::  r37=  0.32308904703711D-02*9.D0
    real(8), parameter    ::  r38=  0.80964802996215D-04*10.D0
    real(8), parameter    ::  r39= -0.16557679795037D-03*10.D0
    real(8), parameter    ::  r40= -0.44923899061815D-04*11.D0
    real(8), parameter    ::  s01=  0.10658070028513D01
    real(8), parameter    ::  s13=  0.38493460186671D0*2.D0
    real(8), parameter    ::  s14= -0.85214708824206D0*2.D0
    real(8), parameter    ::  s15=  0.48972281541877D01*2.D0
    real(8), parameter    ::  s16= -0.30502617256965D01*2.D0
    real(8), parameter    ::  s17=  0.39420536879154D-01*2.D0
    real(8), parameter    ::  s18=  0.12558408424308D0*2.D0
    real(8), parameter    ::  s19= -0.27999329698710D0*6.D0
    real(8), parameter    ::  s20=  0.13899799569460D01*6.D0
    real(8), parameter    ::  s21= -0.20189915023570D01*6.D0
    real(8), parameter    ::  s22= -0.82147637173963D-02*6.D0
    real(8), parameter    ::  s23= -0.47596035734923D0*6.D0
    real(8), parameter    ::  s24=  0.43984074473500D-01*12.D0
    real(8), parameter    ::  s25= -0.44476435428739D0*12.D0
    real(8), parameter    ::  s26=  0.90572070719733D0*12.D0
    real(8), parameter    ::  s27=  0.70522450087967D0*12.D0
    real(8), parameter    ::  s28=  0.10770512626332D0*20.D0
    real(8), parameter    ::  s29= -0.32913623258954D0*20.D0
    real(8), parameter    ::  s30= -0.50871062041158D0*20.D0
    real(8), parameter    ::  s31= -0.22175400873096D-01*30.D0
    real(8), parameter    ::  s32=  0.94260751665092D-01*30.D0
    real(8), parameter    ::  s33=  0.16436278447961D0*30.D0
    real(8), parameter    ::  s34= -0.13503372241348D-01*42.D0
    real(8), parameter    ::  s35= -0.14834345352472D-01*56.D0
    real(8), parameter    ::  s36=  0.57922953628084D-03*72.D0
    real(8), parameter    ::  s37=  0.32308904703711D-02*72.D0
    real(8), parameter    ::  s38=  0.80964802996215D-04*90.D0
    real(8), parameter    ::  s39= -0.16557679795037D-03*90.D0
    real(8), parameter    ::  s40= -0.44923899061815D-04*110.D0
    real(8), parameter    ::  H10= -0.36668082266957D-05
    real(8), parameter    ::  H11=  0.17887679267013D-06, R1=461.526D-6
    real(8), parameter    ::  R=0.461526D3
    real(8)               ::  TN
    real(8)               ::  PN
    real(8)               ::  tau
    real(8)               ::  TAU2
    real(8)               ::  taun
    real(8)               ::  taun2
    real(8)               ::  taun3
    real(8)               ::  taun6
    real(8)               ::  taun12
    real(8)               ::  tauninv
    real(8)               ::  tauninv2
    real(8)               ::  tauninv3
    real(8)               ::  tauninv7
    real(8)               ::  PI
    real(8)               ::  pin
    real(8)               ::  pin2
    real(8)               ::  pin3
    real(8)               ::  pin6
    real(8)               ::  gapp
    real(8)               ::  GA0
    real(8)               ::  DVDPI


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [water_PT2DVDP] in mod_steamtable'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    TN= 1386.D0
    PN= 16.53D0

    tau= tn/t
    TAU2=TAU*TAU
    taun= tau - 1.222d0
    taun2= taun*taun
    taun3= taun2*taun
    taun6= taun3*taun3
    taun12= taun6*taun6
    tauninv= 1.d0/taun
    tauninv2= tauninv*tauninv
    tauninv3= tauninv2*tauninv
    tauninv7= tauninv3*tauninv3*tauninv

    PI=P/PN
    pin= 7.1d0 - p/pn
    pin2= pin*pin
    pin3= pin2*pin
    pin6= pin3*pin3

    gapp= tauninv3*(c15 + taun3*(c16 + taun*(c17 +            &
          taun2*(c18 + taun12*taun2*c19))) +                  &
          pin*  tauninv*(c20 + taun3*taun*(c21 + taun6*c22) + &
          pin*  tauninv*(c23 + taun3*(c24 + taun12*c25) +     &
          pin*  tauninv3*(c26 +                               &
          pin3* (tauninv3*(c27 + taun3*taun2*c28) +           &
          pin6*pin6*pin*tauninv7*tauninv7*tauninv7*(c29 +     &
          pin2* tauninv2*(c30 +                               &
          pin6* tauninv7*(c31 +                               &
          pin*  tauninv*(c32 +                                &
          pin*  tauninv*(c33 +                                &
          pin*  tauninv*c34))) ))) ))) )

    GA0=0.D0

    DVDPI=R*T/P/P*PN*(PI*PI*GAPP-GA0)
    DVDP=DVDPI/PN*1.D-6
    !dvdp in m^3/(kg*MPa)

endfunction
!----------------------------------------------------------------------------


!IAPWS Industrial Formulation 2008
!----------------------------------------------------------------------------
function water_PDT2Eta(Di,T) result(Eta)
! T   is temperature in K
! Di  is density in g/cc
! Eta is viscoity in Pa-s
    implicit none

    real(8), intent(in) :: Di
    real(8), intent(in) :: T
    real(8)             :: Eta
    real(8),parameter,dimension(4) ::  &! Parameters for eq. 11 in table 1.
    A=[1.67752D0,2.20462D0,0.6366564D0,-0.241605D0]
    real(8), parameter, dimension(21) ::  &! Parameters for eq. 12 in table 2.
     H=[ 5.20094D-1,  8.50895D-2,  -1.08374D0, -2.89555D-1,          &
     2.22531D-1,  9.99115D-1,   1.88797D0,   1.26613D0,  1.20573D-1, &
    -2.81378D-1, -9.06851D-1, -7.72479D-1, -4.89837D-1, -2.57040D-1, &
     1.61913D-1,  2.57399D-1, -3.25372D-2,  6.98452D-2,              &
     8.72102D-3, -4.35673D-3, -5.93264D-4]
    integer(4), parameter, dimension(21) :: & ! Exponential orders for eq. 12 in table 2.
        I=[0,1,2,3,0,1,2,3,5,0,1,2,3,4,0,1,0,3,4,3,5]
    integer(4), parameter, dimension(21) :: &
        J=[0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,4,4,5,6,6]
    real(8), parameter ::  T_ref=647.10D0       ! Tempeature [K]
    real(8), parameter ::  D_ref=322.0D0        ! Density [kg/m^3]
    real(8)    :: Tr,Th,D,Dr,Dh
    real(8)    :: u0,u1  ! u = u0*u1*u2 where u2=1 for ind. applications
    integer(4) :: ii


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [water_PDT2Eta] in mod_steamtable'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    D = Di * 1D+3 ! g/cc -> kg/m3
    ! Normalize Values
    Tr=T/T_ref
    Dr=D/D_ref

    ! Initialize Variables
    u0=0.0D0
    u1=0.0D0
    Eta=0.0D0

    ! Calculate u0 using equation 11
    do ii = 1, size(A,1)
      u0=u0+A(ii)/(Tr**(ii-1))
    end do
    u0=100.0D0*sqrt(Tr)/u0

    ! Calculate u1 using equation 12
    do ii = 1, size(H,1)
      Th=(1.0D0/Tr-1.0D0)**I(ii)
      Dh=(Dr-1.0D0)**J(ii)
      u1=u1+Th*Dh*H(ii)
    end do
    u1=exp(Dr*u1)

    ! Calculate Dynamic Viscosity using equation 10
    Eta=u0*u1*1D-6 ! Convert from uPa-s to Pa-s
    return

end function
!----------------------------------------------------------------------------


!IAPWS Industrial Formulation 2011
!----------------------------------------------------------------------------
function water_PDT2C(P,Di,T) result(C)
!input   P     : pressure in MPa
!        Di    : density in g/cc
!        T     : temperature in K
!output  C     : conductivity in W/m-K

    implicit none

    real(8), intent(in) :: P,Di,T
    real(8)             :: C
    real(8),parameter,dimension(5) ::    & ! Table 1. Coefficients for  eq. 11
    L0=[2.443221D-3,1.323095D-2,6.770357D-3,-3.454586D-3,4.096266D-4]

    real(8),parameter,dimension(5,6) ::  & ! Table 2. Coefficients for eq. 17
        L1=reshape([   1.60397357D0,    2.33771842D0,   2.19650529D0, &
                      -1.21051378D0,    -2.7203370D0,                 &
                     -0.646013523D0,   -2.78843778D0,  -4.54580785D0, &
                       1.60812989D0,    4.57586331D0,                 &
                      0.111443906D0,    1.53616167D0,   3.55777244D0, &
                     -0.621178141D0,   -3.18369245D0,                 &
                      0.102997357D0,  -0.463045512D0,  -1.40944978D0, &
                     0.0716373224D0,     1.1168348D0,                 &
                    -0.0504123634D0,  0.0832827019D0,  0.275418278D0, &
                              0.0D0,   -0.19268305D0,                 &
                    0.00609859258D0, -0.0071901245D0,-0.0205938816D0, &
                              0.0D0,   0.012913842D0],                &
                     shape(L1))

    real(8),parameter,dimension(5,6) ::   &  ! Table 4. Coefficients for eq. 25
        A=reshape(                                                      &
        [ 6.53786807199516D0,  6.52717759281799D0, 5.35500529896124D0,  & ! i=1
          1.55225959906681D0,  1.11999926419994D0,                      &
         -5.61149954923348D0, -6.30816983387575D0,-3.96415689925446D0,  & ! i=2
         0.464621290821181D0, 0.595748562571649D0,                      &
          3.39624167361325D0,  8.08379285492595D0, 8.91990208918795D0,  & ! i=3
          8.93237374861479D0,  9.88952565078920D0,                      &
         -2.27492629730878D0, -9.82240510197603D0,-12.0338729505790D0,  & ! i=4
         -11.0321960061126D0, -10.3255051147040D0,                      &
          10.2631854662709D0,  12.1358413791395D0, 9.19494865194302D0,  & ! i=5
          6.16780999933360D0,  4.66861294457414D0,                      &
          1.97815050331519D0, -5.54349664571295D0,-2.16866274479712D0,  & ! i=6
        -0.965458722086812D0,-0.503243546373828D0]                      &
        , shape(A) )

    real(8), parameter, dimension(5) ::  & ! Ranges for array A from Table 4 in Dimensionless Density
        B=[0.0D0,0.310559006D0,0.776397516D0,1.242236025D0,1.863354037D0]
    ! Parameters for k2

    real(8), parameter ::  & ! Parameters from table 3
        lambda=177.8514D0, qd=1.0/(0.40D-0),  T_crit=1.5D0,  &
        Gamma_o=0.06D0, u_g=0.630D0/1.239D0, xi_o=0.13D0,    &
        pi = 4.*atan(1.) ! This is just pi(numerical const.)

    real(8), parameter ::   & ! Reference Values at Critical Point
       T_ref=647.10D0,      & ! Tempeature [K]
       D_ref=322.0D0,       & ! Density [kg/m^3]
       P_ref=22.064D0,      & ! Pressure [MPa]
       R=0.46151805D0,      & ! specific gas constant [kJ/kg-K]
       u_ref=1.0D-6           ! Viscosity [Pa-s]

    ! Variables used in calculations
    real(8) :: Tr,Th,D,Dr,Dh ! Normalized and Manipulated Temperature and Densities
    real(8) :: zeta1, zeta2, Chi,xi,y,z,mu,cp,cv,cp_cv ! Calculated suport variables
    real(8) :: k0, k1, k2 ! k=(k0*k1+k2) Base equation of calculation

    ! Index Counters
    integer(4) :: i,j,jj

    ! Call regsopt if needed
    ! Warning! regsopt might change the intent of ireg to in,out
    ! Make sure IREG is reset after this function is called
    ! D -> kg/m3

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [water_PDT2C] in mod_steamtable'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

    D = Di * 1D+3
    ! Normalize Values
    Tr=T/T_ref
    Dr=D/D_ref

    ! Initialize Variables
    zeta2=0.0
    k0=0.0
    k1=0.0
    k2=0.0
    C=0.0

    !--- Calculate k0 using eq. 16
    do i = 1, size(L0)
      k0=k0+L0(i)/(Tr**(i-1))
    end do
    k0=sqrt(Tr)/k0

    !--- Calculate k1 using eq. 17
    do i=1, size(L1,1)
      Th=(1.0D0/Tr-1.0D0)**(i-1)
      do j=1, size(L1,2)
        Dh=(Dr-1.0D0)**(j-1)
        k1=k1+Th*Dh*L1(i,j)
      end do
    end do
    k1=exp(Dr*k1)

    !--- Calculate k2 using industrial formulation

    ! Select the correct jj value for array A, from table 6
    jj=size(B)
    do j=1, ( size(B,1) -1 )
      if( (Dr>B(j)) .and. (Dr<=B(j+1)) ) then
        jj=j
      end if
    end do

    ! Compute zeta2 equation 25
    do i=1, size(A,2)
      zeta2=zeta2+(A(jj,i)*Dr**(i-1))
    end do
    zeta2=1.0D0/zeta2

    ! Compute dimensionless isothermal compressibility using eq. 24
    ! Evaluated as dv/dP * D^2 *P_ref/D_ref
    ! where v is specific volume, D is density, and P is pressure
    zeta1=-D*D*P_ref/D_ref*(water_PT2DVDP(P,T))  ! dvdp in m^3/(kg*MPa)
    if(zeta1<0.0D0) zeta1=0.0D0
    if(zeta1>1D13) zeta1=1D13

    ! Calculate Chi from equation 23
    Chi=Dr*(zeta1-zeta2*T_crit/Tr)
    if(Chi<0.0D0) Chi=0.0D0

    xi=xi_o*(Chi/Gamma_o)**(u_g)  ! Calculate xi from equation
    cp=water_PT2Cp(P,T)           ! Isobaric Specific Heat  [kJ/kg-K]
    cv=water_PT2Cv(P,T)           ! Isochoric Specific Heat [kJ/kg-K]
    y=qd*xi                       ! Dimensionless value for chord length [nm/nm]
    cp_cv=cp/cv                   ! Ratio of specific heats
    mu=water_PDT2Eta(Di,T)      ! Dynamic Viscosity   [Pa-s]
    cp=cp/R                       ! Reference to R
    mu=mu/u_ref                   ! Reference to ref. visc.

    ! Compute z according to equation
    if(y<1.2D-7) then
      z=0.0D0
    else
      z=2/(pi*y)*( ((1-1.0D0/cp_cv)*atan(y)+y/cp_cv)-  &
        (1-exp(-1.0D0/(1.0D0/y+y*y/(3.0D0*Dr*Dr)) ) ) )
    end if
    if(z<0.0D0) then
      z=0.0D0
    end if

    ! Compute k2
    k2=0

    ! Combine correlations to evaluate thermal conductivity
    C=(k0*k1+k2)*1.0D-3 ! Conversion factor from mW to W

    return

end function
!----------------------------------------------------------------------------


endmodule mod_steamtable

#endif
