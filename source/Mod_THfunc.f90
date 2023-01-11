
      MODULE Mod_THfunc

      USE Inc_Constant
      USE Inc_TH


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE
#ifdef tuan_frth
        real(8)      :: x_Pu                                                ! weight fraction Pu/(Pu+U)
        real(8)      :: x_O2M                                               ! oxygen excess (U.Pu)O2+x
        real(8)      :: TD                                                  ! Theory density
        real(8)      :: mol_mass                                            ! Mol mass
        real(8)      :: t_melting                                           ! melting temperature
        real(8)      :: x_Zr                                                !
        real(8)      :: weight_Zr                                           !
        real(8)      :: weight_Pu                                                !
#endif
      CONTAINS

#ifdef tuan_frth

!!!#ifdef siarhei_delete 
    subroutine linear_interpolation (x, y, x_point, y_point)

        real(8), intent(in)  :: x(:)
        real(8), intent(in)  :: y(:)
        real(8), intent(in)  :: x_point
        real(8), intent(out) :: y_point

        integer  :: i

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [linear_interpolation] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif


        do i = LBOUND(x, dim=1), UBOUND(x, dim=1)
            if (x_point <= x(i) + 1.0D-09)  then
                y_point = y(i)
                exit

            else if (x_point < x(i) - 1.0D-09)  then
                y_point = y(i-1) + (x_point-x(i-1)) * ((y(i)-y(i-1))/(x(i)-x(i-1)))
                exit
            end if
        end do

    end subroutine linear_interpolation
!!!#endif 

!!!#ifdef siarhei_delete 
   subroutine FuelProperty
!        real(8)      :: x_Pu                                                ! weight fraction Pu/(Pu+U)
!        real(8)      :: x_O2M                                               ! oxygen excess (U.Pu)O2+x
!        real(8)      :: TD                                                  ! Theory density
!        real(8)      :: mol_mass                                            ! Mol mass
!        real(8)      :: t_melting                                           ! melting temperature
!        real(8)      :: x_Zr                                                !
!        real(8)      :: weight_Zr                                           !
!        real(8)      :: weight_Pu                                                !

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [FuelProperty] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif
   
      select case(fuel_type)
        case(1)                                                                 ! U
            t_melting = 1405.0D0
            mol_mass  = 238.0D0

        case(2)                                                                 ! Pu
            t_melting = 913.0D0
            mol_mass  = 244.0D0

        case(3)                                                                 ! Pu-Zr
            weight_Zr = 0.618D0
!            if (PRESENT(weight_Zr))  then
!                weight_Zr = weight_Zr
!            end if

            weight_Pu = 1.0D0 - weight_Zr

            x_Zr = (weight_Zr*244.0D0) / (91.22D0 + weight_Zr*(244.0D0-91.22D0))
            x_Pu = 1.0D0 - x_Zr

            t_melting = x_Zr*2128.0D0 + x_Pu*913.0D0
            mol_mass  = x_Zr*91.22D0 + x_Pu*244.0D0
        case(4)                                                                 ! IAEA-UO2
            t_melting = 3120.0D0
            mol_mass  = 270.3D0
            TD        = 0.95D0

            x_Pu  = 0.0D0
            x_O2M = 2.0D0

        case(5)                                                                 ! IAEA-PuO2
            t_melting = 2663.0D0
            mol_mass  = 276.045D0
            TD        = 0.95D0

            x_Pu  = 1.0D0
            x_O2M = 2.0D0

        case(6)                                                                 ! IAEA-MOX (U0.8Pu0.2)O2+x
            t_melting = 3023.0D0
            mol_mass  = 271.2D0
            TD        = 0.95D0

            x_Pu  = 0.2D0*276.045D0 / (0.2D0*276.045D0 + 0.8D0*270.3D0)
            x_O2M = 2.0D0

!            if (PRESENT(x_Pu))  then
!                x_Pu  = x_Pu*276.045D0 / (x_Pu*276.045D0 + (1.0D0-x_Pu)*270.3D0)
!            end if
!            if (PRESENT(x_O2M))  then
!                x_O2M  = x_O2M
!            end if

        case(7)                                                                 ! PARCS-UO2, same as 1
            t_melting = 3120.0D0
            mol_mass  = 270.3D0
            TD        = 0.95D0

            x_Pu  = 0.0D0
            x_O2M = 2.0D0

        case(8)                                                                 ! PARCS-MOX, same as 3
            t_melting = 3023.0D0
            mol_mass  = 271.2D0
            TD        = 0.95D0

            x_Pu  = 0.2D0*276.045D0 / (0.2D0*276.045D0 + 0.8D0*270.3D0)
            x_O2M = 2.0D0

!        case(6)                                                                 ! COBRA-MATPRO-UO2, same as 1
!            t_melting = 3120.0D0
!            mol_mass  = 270.3D0
!            TD        = 0.95D0
!
!            x_Pu  = 0.0D0
!            x_O2M = 2.0D0
!
!        case(7)                                                                 ! COBRA-NEA-UO2, same as 1
!            t_melting = 3120.0D0
!            mol_mass  = 270.3D0
!            TD        = 0.95D0
!
!            x_Pu  = 0.0D0
!            x_O2M = 2.0D0
!
        case(9)                                                                 ! RELAP5-UO2, same as 1
            t_melting = 3120.0D0
            mol_mass  = 270.3D0
            TD        = 0.95D0

            x_Pu  = 0.0D0
            x_O2M = 2.0D0

        case(10)                                                                 ! RELAP5-MOX (U1-c.Puc)O2+x, same as 3
            t_melting = 3023.0D0
            mol_mass  = 271.2D0
            TD        = 0.95D0

            x_Pu  = 0.2D0
            x_O2M = 2.0D0

!            if (PRESENT(x_Pu))  then
!                x_Pu  = x_Pu*276.045D0 / (x_Pu*276.045D0 + (1.0D0-x_Pu)*270.3D0)
!            end if
!            if (PRESENT(x_O2M))  then
!                x_O2M  = x_O2M
!            end if

        case(11)                                                                ! OECD/NEA beam trip, same as 1
            t_melting = 3120.0D0
            mol_mass  = 270.3D0
            TD        = 0.95D0

            x_Pu  = 0.0D0
            x_O2M = 2.0D0

!        case default

        case(12)                                                                ! OECD/NEA beam trip, same as 1
            t_melting = 3123.0D0
            mol_mass  = 252.3D0


      end select

   end subroutine FuelProperty
!!!#endif 

!!!#ifdef siarhei_delete 
   subroutine CladProperty
   
   
#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [CladProperty] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

        select case(clad_type)
        case(1)                                                                 ! IAEA, 316
            t_melting = 1703.0D0
            mol_mass  = 55.9354D0

        case(2)                                                                 ! RELAP5, 304, as 1
            t_melting = 1703.0D0
            mol_mass  = 55.9354D0

        case(3)                                                                 ! HT9, as 1
            t_melting = 1703.0D0
            mol_mass  = 55.9354D0

        case(4)                                                                 ! OECD/NEA beam trip , as 1
            t_melting = 1703.0D0
            mol_mass  = 55.9354D0
        case(5)                                                                 ! 15-15Ti
            t_melting = 1673.0D0
            mol_mass  = 55.9354D0
        end select

   end subroutine CladProperty
!!!#endif 


    function f_lm_temp_search(h, t0) result (t)
        implicit none
        integer :: count = 0
        real(8) :: h, t0
        real(8) :: ent, cap, delta, t
        t = t0

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [f_lm_temp_search] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

        delta = 1.D+9
        do while (abs(delta) > 1.D-5)
            ent = fenthal(t) ! siarhei_rev
!            dummy_filler = 1 ! @$^ siarhei_plot 
            cap = fhcap(t)
            delta = (h - ent) / cap
            t = t + delta
            count = count + 1
            if (count > 100) exit
        enddo
        return
    end function f_lm_temp_search


#endif


      FUNCTION fkf(t, ifa, i_z, fkf_th1d)
      use Inc_Solver
      use Inc_Option,    ONLY: OPT_Tburn
      use Inc_3D,        ONLY: BU
      use Inc_RP,        ONLY: I_Gd_FA
      use Inc_RP,        ONLY: I_FA
      use Inc_MATFB
      use inc_inp, only: fuelk_con

      IMPLICIT NONE

      ! thermal conductivity of uo2 in w/m-C, t in K
      INTEGER(4) :: ifa, i_z, i
      REAL(8) :: fkf, t, wogd, lht, rht, fkf_th1d, frap_fuelk_bu_ratio, fkf0
      real(8) :: a, b, c, d

#ifdef tuan_frth
        real(8)  :: key(2)
        real(8)  :: value(2)
        real(8)  :: wdeltap                                                 ! wdeltap, the weight fraction of delta-pu phase
        real(8)  :: vdeltap                                                 ! valphap, the volume fraction of delta-pu phase
        real(8)  :: valphaz                                                 ! vdeltaz, the volume fraction of alpha-zr phase
        real(8)  :: density_Zr(2,8), density_Pu(2,11)
        real(8)  :: rhozr,rhopu,rhodeltap,rhoalphaz
        real(8)  :: cd,cdalphazr,cddeltapu,wz,wp
        real(8)  :: degree
        real(8)      :: TD                                                  ! theory density
        real(8)      :: x_Pu                                                ! weight fraction Pu/(Pu+U)
        real(8)      :: x_O2M                                               ! oxygen excess (U.Pu)O2+x
        real(8)  :: e_th, Cv, t_1, t_2
        real(8)  :: Cv_U, Cv_Pu

        integer  :: j
#endif

      ! FRAPCON fuel_k_cond BU ratio

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fkf] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      SELECT CASE(OPT_Tburn)
      CASE(2,3,4)

      ! At BU point
      i = I_FA(ifa)
      FKB_h = 1d0/(1d0 + 396d0 * exp(-FKB_Q/t))
      wogd = I_Gd_FA(ifa, i_z)
      lht = FKB_A + FKB_sA * wogd + FKB_B * t + 0.00187d0 * BU(i, i_z) &
      + (1d0 - 0.9d0 * exp(-0.04d0 * (BU(i, i_z)))) * 0.038d0 * BU(i, i_z)**(0.28d0) * FKB_h
      rht = FKB_E/(t * t) * exp(-FKB_F/t)
      fkf = 1/lht + rht
      fkf_th1d = akfuel(0)+t*(akfuel(1)+t*(akfuel(2)+t*akfuel(3))) &
          +akfuel(4)/(t-akfuel(5))
      fkf = fkf_th1d + fuelk_con * (fkf - fkf_th1d)

      ! At BU=0
      i = I_FA(ifa)
      FKB_h = 1d0/(1d0 + 396d0 * exp(-FKB_Q/t))
      wogd = I_Gd_FA(ifa, i_z)
      lht = FKB_A + FKB_sA * wogd + FKB_B * t + 0.00187d0 * 0d0 &
      + (1d0 - 0.9d0 * exp(-0.04d0 * (0d0))) * 0.038d0 * 0d0**(0.28d0) * FKB_h
      rht = FKB_E/(t * t) * exp(-FKB_F/t)
      fkf0 = 1/lht + rht
      fkf_th1d = akfuel(0)+t*(akfuel(1)+t*(akfuel(2)+t*akfuel(3))) &
          +akfuel(4)/(t-akfuel(5))
      fkf0 = fkf_th1d + fuelk_con * (fkf0 - fkf_th1d)
      frap_fuelk_bu_ratio = fkf/fkf0
      END SELECT

      IF (OPT_Tburn == 1) THEN
         i = I_FA(ifa)
         FKB_h = 1d0/(1d0 + 396d0 * exp(-FKB_Q/t))
         wogd = I_Gd_FA(ifa, i_z)
         lht = FKB_A + FKB_sA * wogd + FKB_B * t + 0.00187d0 * BU(i, i_z) &
            + (1d0 - 0.9d0 * exp(-0.04d0 * (BU(i, i_z)))) * 0.038d0 * BU(i, i_z)**(0.28d0) * FKB_h
         rht = FKB_E/(t * t) * exp(-FKB_F/t)
         fkf = 1/lht + rht
         fkf_th1d = akfuel(0)+t*(akfuel(1)+t*(akfuel(2)+t*akfuel(3))) &
             +akfuel(4)/(t-akfuel(5))
         fkf = fkf_th1d + fuelk_con * (fkf - fkf_th1d)
      ELSE IF(OPT_Tburn == 2) THEN
         a = 0.015d0
         b = 1.54835d-4
         c = 3.3d+8
         d = 7040d0
         fkf = 1/(a+b*t) + c/(t*t)*exp(-d/t)
         fkf = frap_fuelk_bu_ratio * fkf
      ELSE IF(OPT_Tburn == 3) THEN
         a = 0.01544d0
         b = 1.5d-4
         c = 3.3d+8
         d = 7420d0
         fkf = 1/(a+b*t) + c/(t*t)*exp(-d/t)
         fkf = frap_fuelk_bu_ratio * fkf
      ELSE IF(OPT_Tburn == 4) THEN
         a = 0.02d0
         b = 1.51946d-4
         c = 3.3d+8
         d = 7658.92119d0
         fkf = 1/(a+b*t) + c/(t*t)*exp(-d/t)
         fkf = frap_fuelk_bu_ratio * fkf
      ELSE
         fkf=akfuel(0)+t*(akfuel(1)+t*(akfuel(2)+t*akfuel(3))) &
             +akfuel(4)/(t-akfuel(5))
      ENDIF

#ifdef tuan_frth
      call FuelProperty

      select case(fuel_type)
        case(1)                                                                 ! U
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if           
            fkf = 22.0D0 + 0.023D0*(t-273.15D0)
            
        case(2)                                                                 ! Pu
         
        if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
        end if
            ! step-wise approximation
            if (t <= 395.0D0)  then
                key   = [273.15D0, 395.0D0]
                value = [5.2D0, 6.6D0]
                call linear_interpolation (key, value, t, fkf)
            else if (t <= 479.0D0)  then
                key   = [395.0D0, 479.0D0]
                value = [7.87D0, 8.67D0]
                call linear_interpolation (key, value, t, fkf)
            else if (t <= 592.0D0)  then
                key   = [479.0D0, 592.0D0]
                value = [8.97D0, 10.5D0]
                call linear_interpolation (key, value, t, fkf)
            else if (t <= 724.0D0)  then
                key   = [592.0D0, 724.0D0]
                value = [10.97D0, 12.1D0]
                call linear_interpolation (key, value, t, fkf)
            else if (t <= 749.0D0)  then
                fkf = 7.72D0
            else
                fkf = 12.35D0
            end if
        
        case(3)                                                                 ! Pu-Zr
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if

            density_Zr(1,1) = 500.0D0;   density_Zr(2,1) = 6.546D0
            density_Zr(1,2) = 600.0D0;   density_Zr(2,2) = 6.532D0
            density_Zr(1,3) = 700.0D0;   density_Zr(2,3) = 6.518D0
            density_Zr(1,4) = 800.0D0;   density_Zr(2,4) = 6.503D0
            density_Zr(1,5) = 900.0D0;   density_Zr(2,5) = 6.488D0
            density_Zr(1,6) = 1000.0D0;  density_Zr(2,6) = 6.471D0
            density_Zr(1,7) = 1100.0D0;  density_Zr(2,7) = 6.456D0
            density_Zr(1,8) = 1135.0D0;  density_Zr(2,8) = 6.450D0
            
            density_Pu(1,1) = 500.0D0;  density_Pu(2,1) = 17.075D0
            density_Pu(1,2) = 550.0D0;  density_Pu(2,2) = 16.991D0
            density_Pu(1,3) = 593.1D0;  density_Pu(2,3) = 15.843D0              ! extend to 540 K;
            density_Pu(1,4) = 600.0D0;  density_Pu(2,4) = 15.846D0
            density_Pu(1,5) = 700.0D0;  density_Pu(2,5) = 15.884D0
            density_Pu(1,6) = 736.0D0;  density_Pu(2,6) = 15.898D0
            density_Pu(1,7) = 736.0D0;  density_Pu(2,7) = 15.921D0
            density_Pu(1,8) = 755.7D0;  density_Pu(2,8) = 15.935D0
            density_Pu(1,9) = 755.7D0;  density_Pu(2,9) = 16.444D0
            density_Pu(1,10) = 800.0D0; density_Pu(2,10) = 16.369D0
            density_Pu(1,11) = 900.0D0; density_Pu(2,11) = 16.201D0
            
            wz=0.65
            wp=0.35
            temperature: if (540.0D0<=t .AND. t<=870.0D0)  then
                
                ! delta-Pu phase
                if (0.0D0<=wz .AND. wz<=0.45D0)  then
                    a = (1.0D0-SQRT(1.0D0-wp))*3.225D0 + SQRT(1.0D0-wp)*(29.469D0-118.811D0*wp+88.893D0*wp**2)
                    b = (1.0D0-SQRT(1.0D0-wp))*0.0296D0 + SQRT(1.0D0-wp)*(0.0117D0-0.00716D0*wp)
                    c = SQRT(1.0D0-wp)*1.922D-5
                   fkf = a+b*t+c*t**2
                
                ! alpha-Zr phase
                else if (0.8D0<=wz .AND. wz<=1.0D0)  then  
                    a = (1.0D0-SQRT(1.0D0-wz))*8.853D0 + SQRT(1.0D0-wz)*(30.57D0-82.301D0*wz+71.456D0*wz**2)
                    b = (1.0D0-SQRT(1.0D0-wz))*7.082D-3 + SQRT(1.0D0-wz)*(0.01895D0-0.02453D0*wz)
                    c = (1.0D0-SQRT(1.0D0-wz))*2.533D-6 + SQRT(1.0D0-wz)*8.111D-6
                    d = (1.0D0-SQRT(1.0D0-wz))*2.992D3
                    fkf = a+b*t+c*t**2+d/t
                    
                ! delta-pu + alpha-zr phase
                else if (0.45D0<wz .AND. wz<0.8D0)  then                
                    ! bruggeman model -- aaa fuels handbook from ANL
                    ! delta-pu phase contains 45wt% zr and 55wt% pu, alpha-zr phase contains 80wt% zr and 20wt% pu
                    ! according to lever rule
                    wdeltap = (wz-0.8D0)/(0.45D0-0.8D0)
                    rhozr = 6.488D0
                    rhopu = 16.201D0
                    
                    call linear_interpolation (density_Zr(1,:), density_Zr(2,:), t, rhozr)
                    call linear_interpolation (density_Pu(1,:), density_Pu(2,:), t, rhopu)
            
                    rhodeltap = 1.0D0 / (0.45D0/rhozr + 0.55D0/rhopu)
                    rhoalphaz = 1.0D0 / (0.8D0/rhozr + 0.2D0/rhopu)
                    vdeltap = rhoalphaz*wdeltap / (rhodeltap - (rhodeltap-rhoalphaz)*wdeltap)
                    valphaz = 1.0D0 - vdeltap
            
                    ! delta-pu phase, 45wt% zr and 55wt% pu
                    a = (1.0D0-SQRT(1.0D0-0.55D0))*3.225D0 + SQRT(1.0D0-0.55D0)*(29.469D0-118.811D0*0.55D0+88.893D0*0.55D0**2)
                    b = (1.0D0-SQRT(1.0D0-0.55D0))*0.0296D0 + SQRT(1.0D0-0.55D0)*(0.0117D0-0.00716D0*0.55D0)
                    c = SQRT(1.0D0-0.55D0)*1.922D-5
                    cddeltapu = a+b*t+c*t**2
            
                    ! alpha-zr phase, 80wt% zr and 20wt% pu
                    a = (1.0D0-SQRT(1.0D0-0.8D0))*8.853D0 + SQRT(1.0D0-0.8D0)*(30.57D0-82.301D0*0.8D0+71.456D0*0.8D0**2)
                    b = (1.0D0-SQRT(1.0D0-0.8D0))*7.082D-3 + SQRT(1.0D0-0.8D0)*(0.01895D0-0.02453D0*0.8D0)
                    c = (1.0D0-SQRT(1.0D0-0.8D0))*2.533D-6 + SQRT(1.0D0-0.8D0)*8.111D-6
                    d = (1.0D0-SQRT(1.0D0-0.8D0))*2.992D3
                    cdalphazr = a+b*t+c*t**2+d/t
            
                    a = (3.0D0*vdeltap-1.0D0)*cddeltapu+(3.0D0*valphaz-1.0D0)*cdalphazr
            
                    fkf = (a+SQRT(a**2+8.0D0*cdalphazr*cddeltapu))/4.0D0
                end if
                
            else temperature
                a = (1.0D0-SQRT(1.0D0-wz))*8.853D0 + SQRT(1.0D0-wz)*(wz*(-98.806D0+147.895D0*wz-26.883D0*wz**2)+(1.0D0-wz)*9.507D0)
                b = (1.0D0-SQRT(1.0D0-wz))*7.082D-3 + SQRT(1.0D0-wz)*(wz*(0.0512D0-0.0601D0*wz)+(1.0D0-wz)*0.0184D0)
                c = (1.0D0-SQRT(1.0D0-wz))*2.533D-6 + SQRT(1.0D0-wz)*(wz*8.699D-6)
                d = (1.0D0-SQRT(1.0D0-wz))*2.992D3
                
                fkf =  a+b*t+c*t**2+d/t
                
            end if temperature
            
	case(4)                                                                 ! IAEA-UO2
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if            
            t = t / 1000.0D0
            
            fkf = 100.0D0/(7.5408D0+17.692D0*t+3.6142D0*t**2) + (6400.0D0/t**2.5D0)*EXP(-16.35D0/t)
            
        case(5)                                                                 ! IAEA-PuO2
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
            TD        = 0.95D0
            fkf = 8.441D0 - 7.445D-3*t + 2.236D-6*t**2
            fkf = fkf * (1.0D0 - 2.5D0*((1.0D0-TD)/TD))
        
        case(6) 
	   if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
!            x_O2M = x_O2M                                            ! IAEA-MOX (U0.8Pu0.2)O2+x
            x_O2M = 2.0D0
            degree = 2.0D0 - x_O2M
            fkf = 1.0D0/(1.528D0*SQRT(degree+0.00931D0)-0.1055D0+2.885D-4*t) + 76.38D-12*t**3
        
        case(7)                                                                 ! PARCS-UO2  
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if            
            fkf = 1.05D0+t*(0.0D0+t*(0.0D0+t*0.0D0))+2150.0D0/(t-73.15D0)
            
        case(8)        
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if 
            fkf = 1.05D0+t*(0.0D0+t*(0.0D0+t*0.0D0))+2150.0D0/(t-73.15D0)
            fkf = fkf * 0.90D0
            

        case(9)                                                                 ! RELAP5-UO2
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if            
            A    = 0.339D0
            B    = 0.06867D0
            e_th = 1.0D-5*t - 3.0D-3 + 4.0D-2*EXP(-6.9D-20/(1.38D-23*t))
            Cv   = (296.7D0*535.285D0**2*EXP(535.285D0/t)) / (t**2*(EXP(535.285D0/t)-1)**2)
            D    = 0.95D0
            
            if (t < 1364.0D0)  then
                t_1 = 6.5D0 - 0.00469D0*t
            else if (t > 1834.0D0)  then
                t_1 = -1.0D0
            else 
                t_1 = 0.10284D0 + (-0.1D0 - 0.10284D0) * (t-1364.0D0) / (1834.0D0-1364.0D0)
            end if
            
            if (t < 1800.0D0)  then
                t_2 = t
            else if (t > 2300.0D0)  then
                t_2 = 2050.0D0
            else 
                t_2 = 1800.0D0 + (2050.0D0-1800.0D0) * (t-1800.0D0) / (2300.0D0-1800.0D0)
            end if
            
            fkf = (D/(1.0D0+t_1*(1.0D0-D))) * (Cv/(A+B*t_2)/(1.0D0+3*e_th)) + 5.2997D-3*t*EXP(-13358.0D0/t)*(1.0D0+0.169D0*((13358.0D0/t)+2.0D0)**2)
            
        case(10)                                                                 ! RELAP5-MOX (U1-c.Puc)O2+x
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
            x_Pu  = 0.2D0
            x_O2M = 2.0D0
            TD        = 0.95D0
            x_Pu  = x_Pu*276.045D0 / (x_Pu*276.045D0 + (1.0D0-x_Pu)*270.3D0)
            x_O2M  = x_O2M
 
            A    = 0.339D0 + 12.6D0*ABS(2.0D0-x_O2M)
            B    = 0.06867D0 * (1.0D0 + 0.6238D0*x_Pu)
            e_th = 1.0D-5*t - 3.0D-3 + 4.0D-2*EXP(-6.9D-20/(1.38D-23*t))
            D    = TD
        
            Cv_U   = (296.7D0*535.285D0**2*EXP(535.285D0/t)) / (t**2*(EXP(535.285D0/t)-1)**2)
            Cv_Pu  = (347.4D0*571.000D0**2*EXP(571.000D0/t)) / (t**2*(EXP(571.000D0/t)-1)**2)
            Cv = Cv_Pu*x_Pu + Cv_U*(1-x_Pu)
            
            if (t < 1364.0D0)  then
                t_1 = 6.5D0 - 0.00469D0*t
            else if (t > 1834.0D0)  then
                t_1 = -1.0D0
            else 
                t_1 = 0.10284D0 + (-0.1D0 - 0.10284D0) * (t-1364.0D0) / (1834.0D0-1364.0D0)
            end if
            
            if (t < 1800.0D0)  then
                t_2 = t
            else if (t > 2300.0D0)  then
                t_2 = 2050.0D0
            else 
                t_2 = 1800.0D0 + (2050.0D0-1800.0D0) * (t-1800.0D0) / (2300.0D0-1800.0D0)
            end if
            
            fkf = (D/(1.0D0+t_1*(1.0D0-D))) * (Cv/(A+B*t_2)/(1.0D0+3*e_th)) + 5.2997D-3*t*EXP(-13358.0D0/t)*(1.0D0+0.169D0*((13358.0D0/t)+2.0D0)**2)
            
        case(11)                                                                ! OECD/NEA beam trip
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
            fkf = 1.0/(0.042+2.71D-4*t) + 6.9D-11*t**3
!        case(12)                                                                 ! COBRA-MATPRO-UO2
!
!            t = t - CKELVIN
!            t = 32.0D0 + (9.0D0/5.0D0)*t
!            
!            fkf = COBRA_fuel_conductivity (t, real(1.0, 8))
!            fkf = fkf * 6.23D3
!            
!        case(13)                                                                 ! COBRA-NEA-UO2
!           
!            t = t - CKELVIN
!            t = 32.0D0 + (9.0D0/5.0D0)*t
!            
!            fkf = COBRA_fuel_conductivity (t, real(-1.0, 8))
!            fkf = fkf * 6.23D3
!       
        case(12)                                                                 ! UN
            if (t > t_melting)  then
                t  = t_melting
                write(*,*) '              *** WARNING ***  '
                write(*,*) '         Excess melting temperature'
                write(*,*) 'Melting temperature: ', t_melting
                write(*,*) 'Fuel temperature: ', t
                write(*,*) 'Set Fuel temperature = Melting temperature '
            
            end if           
            fkf = 1.41D0 * t**0.39
            
        end select


#endif
      RETURN
      END FUNCTION fkf



      FUNCTION fkc(j, t, ifa, i_z, qc)
      use Inc_Option,    ONLY: OPT_OXcal
      use Inc_3D,        ONLY: BU
      use Inc_Geometry
      use Inc_FA,        ONLY: R_Clad
      use Inc_RP,        ONLY: I_FA
      use Inc_MATFB

      IMPLICIT NONE
      ! thermal conductivity of Zr in w/m-C, t in K
      INTEGER(4) :: i, j, ifa, i_z
      REAL(8) :: fkc,t,afkc,rox_bu,f_h,fox_h,dt,RO_Clad,qflx,qc


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fkc] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (OPT_OXcal == 1) then
         if (j == nrp4) then
            i = I_FA(ifa)
            afkc = OX_c1 - t * (OX_c2 - t * (OX_c3 - t * OX_c4))
            f_h = (sum(MeshSize_z(IzFuelBot:i_z)) - MeshSize_z(i_z)/2d0) / Core_Height
            fox_h = -49.826d0 * f_h**6d0 + 121.55d0 * f_h**5d0 - 113.69d0 * f_h**4d0 &
               + 51.347d0 * f_h**3d0 -10.243d0 * f_h**2d0 + 1.0476d0 * f_h + 0.1413d0
            rox_bu = (0.0005d0 * BU(i, i_z)**3d0 - 0.0089d0 * BU(i, i_z)**2d0 &
               + 0.225d0 * BU(i, i_z) + 6.7479d0) * 1d-6
            RO_Clad = R_Clad + rox_bu
            qflx = qc * afp / zeta
            !js+dt = -qflx*fox_h * R_Clad * log(RO_Clad/R_Clad)/afkc
            dt = -qflx*fox_h * RO_Clad * log(RO_Clad/R_Clad)/afkc
            fkc=akclad(0)+(t+dt)*(akclad(1)+(t+dt)*(akclad(2)+(t+dt)*akclad(3)))
         else
            fkc=akclad(0)+t*(akclad(1)+t*(akclad(2)+t*akclad(3)))
         end if
      elseif (OPT_OXcal == 2) then ! FeCrAl: C36M
         fkc=8.187d0+t*(1.368d-2-t*9.184d-7)
      else
         fkc=akclad(0)+t*(akclad(1)+t*(akclad(2)+t*akclad(3)))
      end if

#ifdef tuan_frth
      call CladProperty
      select case(clad_type)
      case(1)  
          if (t  > t_melting)  then
              t  = t_melting
              !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get conductivity, uppper exceed')
              !  call a_warning%print (FILES%TH_WARNING)
          end if
      
          fkc = 9.248D0 + 0.01571D0*t
          
      case(2)                                                                 ! RELAP5, 304
          if (t  > t_melting)  then
              t  = t_melting
              !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get conductivity, uppper exceed')
              !  call a_warning%print (FILES%TH_WARNING)
          end if
          
          if (300.0D0<=t .AND. t<= 1671.0D0)  then
              fkc = 7.58D0 + 0.0189D0*t
          else if (t <= 1727.0D0)  then
              fkc = 610.9393D0 - 0.342176D0*t
          else 
              fkc = 20.0D0
          end if

	  case(3) 
        ! HT-9 	  
           if (500.0D0<=t .AND. t<=1027.0D0) then
                   fkc = 17.622D0+2.428D-2*t-1.696D-5*t**2
               else
                   fkc = 12.027D0+1.218D-2*t
           end if

      case(4)                                                                 ! OECD/NEA beam trip 
          if (t > t_melting)  then
              t  = t_melting
              !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get conductivity, uppper exceed')
              !  call a_warning%print (FILES%TH_WARNING)
          end if
          
          fkc = 15.4767D0 + t*3.448D-3
      
      case(5)                                                                 ! 15-15Ti
          if (t > t_melting)  then
              t  = t_melting
              !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get density, upper exceed')
              !  call a_warning%print (FILES%TH_WARNING)
          end if    
!            
          fkc = 13.95D0+1.163D-2*t
!      case(6)                                                                 ! IAEA metal-Zr
!          if (TC_In + 273.15 > 2000.0D0)  then
!              t  = 2000.0D0
!              !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
!              !  call a_warning%print (FILES%TH_WARNING)
!          else 
!              t = TC_In + 273.15
!          end if
!          
!          fkc = 8.8527D0+7.0820D-3*t+2.5329D-6*t**2+2.9918D3*t**(-1)
!          
!      case(7)                                                                 ! IAEA- Zr + 1% Nb (E-110)
!          if (TC_In + 273.15 > 1600.0D0)  then
!              t  = 1600.0D0
!              !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
!              !  call a_warning%print (FILES%TH_WARNING)
!          else 
!              t = TC_In + 273.15
!          end if
!          
!          if (300.0D0<=t .AND. t<=1100.0D0)  then
!              fkc = 23.5D0 - 0.0192D0*t + 1.68D-5*t**2
!          else
!              fkc = 1.5D0 + 0.02D0*t
!          end if
!      
!      case(8)                                                                 ! IAEA- Zr + 2.5% Nb (E-125)
!          if (TC_In + 273.15 > 1100.0D0)  then
!              t  = 1100.0D0
!              !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
!              !  call a_warning%print (FILES%TH_WARNING)
!          else 
!              t = TC_In + 273.15
!          end if
!          
!          fkc = 14.0D0 + 0.0115D0*t
!      
!      case(9)                                                                 ! PARCS Zr-alloy
!          if (TC_In + 273.15 > t_melting)  then
!              t  = t_melting
!              !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
!              !  call a_warning%print (FILES%TH_WARNING)
!          else 
!              t = TC_In + 273.15
!          end if
!          
!          fkc = 7.51D0+t*(2.09D-2+t*(-1.45D-5+t*7.67D-9))
!      
!      case(10)                                                                 ! MATPRO Zr-alloy
!          if (TC_In + 273.15 > t_melting)  then
!              t  = t_melting
!              !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
!              !  call a_warning%print (FILES%TH_WARNING)
!          else 
!              t = TC_In + 273.15
!          end if
!          t = t - CKELVIN
!          t = 32.0D0 + (9.0D0/5.0D0)*t
!          
!          this%capacity = COBRA_clad_capacity (t, real(1.0, KREAL))
!          this%capacity = this%capacity * 4186.3D0
!          
!      case(11)                                                                 ! NEA Zr-alloy
!          if (TC_In + 273.15 > t_melting)  then
!              t  = t_melting
!              !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
!              !  call a_warning%print (FILES%TH_WARNING)
!          else 
!              t = TC_In + 273.15
!          end if
!          t = t - CKELVIN
!          t = 32.0D0 + (9.0D0/5.0D0)*t
!          
!          this%capacity = COBRA_clad_capacity (t, real(-1.0, KREAL))
!          this%capacity = this%capacity * 4186.3D0
      
      
      end select
 


#endif
      RETURN
      END FUNCTION fkc


      FUNCTION ftfavg(x,nr,dr2)

      IMPLICIT NONE

      ! fuel volume avg. temperature
      REAL(8) :: ftfavg, x(*), dr2
      REAL(8) :: area , xavg, a
      INTEGER :: i, nrm1, l, nr


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [ftfavg] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      xavg=zero
      i=1
      a=(x(i+1)-x(i))/dr2
      xavg=(7*x(i)+x(i+1))*0.25d0
      DO l=2,nr
         area=2d0/3d0*l*x(l+1)+44d0/3d0*i*x(l)+2d0/3d0*(i-1d0)*x(l-1)
         xavg=xavg+area
         i=l
      END DO
      nrm1=nr-1
      area=(16d0/3d0*nrm1+4d0+5d0/24d0)*x(nr+1)+(10d0/3d0*nrm1+2.25d0)*x(nr) &
           -(2d0/3d0*nrm1+11d0/24d0)*x(nr-1)
      ftfavg=(xavg+area)*dr2/(8.*(nr*nr*dr2))

      RETURN
      END FUNCTION ftfavg


!!!#ifdef siarhei_delete  ! siarhei_rev
      FUNCTION fenthalf(tcel)

      IMPLICIT NONE

      REAL(8) :: fenthalf, t, tcel

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fenthalf] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      t=tcel+273.15d0
      fenthalf=(162.3d0+t*(0.1519d0+t*(-7.970d-5+t*1.601d-8)))*t ! siarhei_rev
!            dummy_filler = 1 ! @$^ siarhei_plot 

      RETURN
      END FUNCTION fenthalf
!!!#endif 


!!!#ifdef siarhei_delete ! siarhei_rev
      FUNCTION frhocpf(t)
      IMPLICIT NONE
      ! volumetric heat capacity of uo2 in J/m^3-C, t in K
      REAL(8) :: frhocpf,t
#ifdef tuan_frth
      REAL(8) :: capacity_Zr,capacity_Pu
      REAL(8) :: Cp_U, Cp_Pu, Cp_PuO2, Cp_UO2
      REAL(8) :: rho, k_1, k_2, k_3, theta, ED, R
      REAL(8) :: key(2), value(2)

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [frhocpf] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call FuelProperty ! siarhei_rev
!            dummy_filler = 1 ! @$^ siarhei_plot 

#endif

      frhocpf=arcpfuel(0)+t*(arcpfuel(1)+t*(arcpfuel(2)+t*arcpfuel(3))) ! siarhei_rev
!            dummy_filler = 1 ! @$^ siarhei_plot 

#ifdef tuan_frth
      select case(fuel_type)
	  case(1)
              if (t > t_melting)  then
                  t  = t_melting
              end if
              
              if (293.0D0<=t .AND. t<=942.0D0)  then
                  frhocpf = 104.82D0 + 5.3686D-3*t + 10.1823D-5*t**2
              else if (t<=1049.0D0)  then
                  frhocpf = 176.4D0
              else
                  frhocpf = 156.8D0
              end if
 
	  case(2)
              if (t > t_melting)  then
                  t  = t_melting
              end if
              
              ! step-wise approximation
              if (t <= 395.0D0)  then
                  key   = [273.15D0, 395.0D0]
                  value = [32.0D0, 34.3D0]
                  call linear_interpolation (key, value, t, frhocpf)
              else if (t <= 479.0D0)  then
                  key   = [395.0D0, 479.0D0]
                  value = [34.3D0, 36.0D0]
                  call linear_interpolation (key, value, t, frhocpf)
              else if (t <= 592.0D0)  then
                  key   = [479.0D0, 592.0D0]
                  value = [34.8D0, 39.8D0]
                  call linear_interpolation (key, value, t, frhocpf)
              else if (t <= 724.0D0)  then
                  frhocpf = 37.7D0
              else if (t <= 749.0D0)  then
                  frhocpf = 37.4D0
              else
                  frhocpf = 35.0D0
              end if
              
              frhocpf = frhocpf * 1000.0D0 / mol_mass
 
          case(3)
	    if (t > t_melting)  then
                t  = t_melting
            end if
      ! Pu-Zr
            capacity_Zr = 0.0D0
            if (298.15D0<=t .AND. t<=1135.0D0)  then
                capacity_Zr = 22.839D0+9.091D-3*t-2.132D4*t**(-2)
            else 
                capacity_Zr = 12.885D0+9.976D-3*t+5.518D6*t**(-2)
            end if
            capacity_Zr = capacity_Zr * 1000.0D0 / 91.22D0
            
            capacity_Pu = 0.0D0
            if (298.15D0<=t .AND. t<=397.6D0)  then
                capacity_Pu = 18.126D0 + 4.482D-2*t
            else if (t <= 487.9D0)  then
                capacity_Pu = 27.416D0 + 1.306D-2*t
            else if (t <= 593.1D0)  then
                capacity_Pu = 22.023D0 + 2.296D-2*t
            else if (t <= 736.0D0)  then
                capacity_Pu = 28.478D0 + 1.081D-2*t
            else if (t <= 755.7D0)  then
                capacity_Pu = 35.560D0
            else if (t <= 913.0D0)  then
                capacity_Pu = 33.720D0
            else
                capacity_Pu = 33.720D0
            end if
            capacity_Pu = capacity_Pu * 1000.0D0 / 244.0D0
            
            frhocpf = capacity_Zr*0.81 + capacity_Pu*0.19

          case(4)                                                                 ! IAEA-UO2
            if (t  > t_melting)  then
                t  = t_melting

            end if
            t = t / 1000.0D0
            
            frhocpf = 52.1743D0+87.951D0*t-82.2411D0*t**2+31.542D0*t**3-2.6334D0*t**4-0.71391D0*t**(-2)
            frhocpf = frhocpf * 1000.0D0/ mol_mass
            
        case(5)                                                                 ! IAEA-PuO2
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
            
            frhocpf = -4.243D-6*t**2 + 2.366D-3*t + 293.1D0
        
        case(6)                                                                 ! IAEA-MOX (U0.8Pu0.2)O2+x
            if (t  > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
            t = t / 1000.0D0
            
            Cp_U = 52.1743D0+87.951D0*t-82.2411D0*t**2+31.542D0*t**3-2.6334D0*t**4-0.71391D0*t**(-2)
            Cp_U = Cp_U * 1000.0D0/ mol_mass
            
            Cp_Pu = -4.243D-6*t**2 + 2.366D-3*t + 293.1D0
        
            frhocpf = 0.2D0*Cp_Pu + 0.8D0*Cp_U
        
        case(7)                                                                 ! PARCS-UO2
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if

            rho = 10282.0D0
            rho = 10412.0D0
            
            frhocpf = 162.3D0*rho+t*(0.3038D0*rho+t*(-2.391D-4*rho+t*6.404D-8*rho))
            frhocpf = frhocpf / rho
            
        case(8)                                                                 ! PARCS-MOX
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if

            rho = 11080.0D0
            
            frhocpf = 162.3D0*rho+t*(0.3038D0*rho+t*(-2.391D-4*rho+t*6.404D-8*rho))
            frhocpf = frhocpf / rho
            
!        case(6)                                                                 ! COBRA-MATPRO-UO2
!            if (TF_In + 273.15 > t_melting)  then
!                t  = t_melting
!                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
!                !  call a_warning%print (FILES%TH_WARNING)
!            else 
!                t = TF_In + 273.15
!            end if
!            t = t - CKELVIN
!            t = 32.0D0 + (9.0D0/5.0D0)*t
!            
!            frhocpf = COBRA_fuel_capacity (t, real(1.0, KREAL))
!            frhocpf = frhocpf * 4186.3D0
!            
!        case(7)                                                                 ! COBRA-NEA-UO2
!            if (TF_In + 273.15 > t_melting)  then
!                t  = t_melting
!                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
!                !  call a_warning%print (FILES%TH_WARNING)
!            else 
!                t = TF_In + 273.15
!            end if
!            t = t - CKELVIN
!            t = 32.0D0 + (9.0D0/5.0D0)*t
!            
!            frhocpf = COBRA_fuel_capacity (t, real(-1.0, KREAL))
!            frhocpf = frhocpf * 4186.3D0
!            
        case(9)                                                                 ! RELAP5-UO2
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
            
            frhocpf = (296.7D0*535.285D0**2*EXP(535.285D0/t)) / (t**2*(EXP(535.285D0/t)-1)**2) + 2.43D-2*t    &
                &   + ((2*8.745D7*1.577D5) / (2*8.3143D0*t**2)) * EXP(-1.577D5/(8.3143D0*t))
            
        case(10)                                                                 ! RELAP5-MOX (U1-c.Puc)O2+x
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
            
            k_1   = 296.7D0
            k_2   = 2.43D-2
            k_3   = 8.745D7
            theta = 535.285D0
            ED    = 1.577D5
            R     = 8.3143D0
            Cp_U = (k_1*theta**2*EXP(theta/t)) / (t**2*(EXP(theta/t)-1)**2) + k_2*t    &
                &   + ((x_O2M*k_3*ED) / (2*R*t**2)) * EXP(-ED/(R*t))
                
            k_1   = 347.4D0
            k_2   = 3.95D-4
            k_3   = 3.860D7
            theta = 571.000D0
            ED    = 1.967D5
            R     = 8.3143D0
            Cp_Pu = (k_1*theta**2*EXP(theta/t)) / (t**2*(EXP(theta/t)-1)**2) + k_2*t    &
                &   + ((x_O2M*k_3*ED) / (2*R*t**2)) * EXP(-ED/(R*t))
                
            frhocpf = (x_Pu*276.045D0*Cp_Pu + (1-x_Pu)*270.3D0*Cp_U) / (x_Pu*276.045D0 + (1-x_Pu)*270.3D0)
            
        case(11)                                                                ! OECD/NEA beam trip
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if

            Cp_UO2 = 81.825 + 0.78695*t - 1.1552E-3*t**2 + 9.9037E-7*t**3 - 5.1982E-10*t**4 + 1.5241E-13*t**5 - 1.7906E-17*t**6
            Cp_PuO2 = -4.9236E6/(t**2) + 240.89 + 0.32556*t - 3.5398E-4*t**2 + 1.512E-7*t**3 - 1.9707E-11*t**4
            frhocpf = (214.65*Cp_UO2 + 55.56*Cp_PuO2) / 270.21

      case(12)
!write(*,*) 'tuan1: ', frhocpf
              if (t > t_melting)  then                     ! UN
                  t  = t_melting
              end if
!                 from IAEA              
                  frhocpf = (0.2029D0*(365.7/t)**2*exp(365.7/t)/(exp(365.7/t)-1)**&
                      2+3.766D-5*t+1.048D9*exp(-18081/t)/t**2)*1000.0D0
!write(*,*) 'frhocpf1: ', frhocpf
                  frhocpf = 252*0.001*(51.14*(365.7/t)*exp(365.7/t)/(exp(365.7/t)-1)**&
                      2+9.491D-3*t+2.642D11*exp(-18081/t)/t**2)
            
!write(*,*) 'frhocpf2: ', frhocpf
        end select
 

#endif

      RETURN
      END FUNCTION frhocpf
!!!#endif 

!!!#ifdef siarhei_delete 
      FUNCTION frhocpc(t)

      IMPLICIT NONE

      ! volumetric heat capacity of Zr in J/m^3-C, t in K
      REAL(8) :: frhocpc,t

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [frhocpc] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      frhocpc=arcpclad(0)+t*(arcpclad(1)+t*(arcpclad(2)+t*arcpclad(3)))
!            dummy_filler = 1 ! @$^ siarhei_plot 
#ifdef tuan_frth
        call CladProperty
 
        select case(clad_type)
        case(1)                                                                 ! IAEA, 316
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
            
            frhocpc = 462.0D0 + 0.134D0*t
            
        case(2)                                                                 ! RELAP5, 304
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
            
            if (300.0D0<=t .AND. t<= 1671.0D0)  then
                frhocpc = 326.0D0 - 0.242D0*t + 3.71D0*(t**0.719D0)
            else 
                frhocpc = 691.98D0
            end if
            
        case(3)                                                                 ! HT9, same as 1
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
        
            frhocpc = 460.0D0 + 0.134D0*t
            
        case(4)                                                                 ! OECD/NEA beam trip, same as 1
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
        
            frhocpc = 620.0
        case(5)                                                                 ! 15-15Ti
            if (t > t_melting)  then
                t  = t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            end if
            
            frhocpc = 431.0D0+0.177D0*t + 8.72D-5*t**2
            
        end select


! end added by Tuan      
      ! HT-9
!       frhocpc= 460.0D0 + 0.134D0*t

#endif

      RETURN
      END FUNCTION frhocpc
!!!#endif 


      FUNCTION fcond(t)

      IMPLICIT NONE

      ! cubic poly for thermal conductivity as a func of temperature
      REAL(8) :: fcond, t
#ifdef tuan_frth
      REAL(8) :: tk


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fcond] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      tk = t + 273.15d0

#endif

      IF ( fluid_type == 0 ) THEN
         fcond = 8.9182016D-1 - 2.1996892D-3*t + 9.9347652D-6*t*t - 2.0862471D-8*t*t*t
      ELSE IF( fluid_type == 1 ) THEN
         fcond = 1.1075040D+0 - 4.8694050D-3*t + 1.7068352D-5*t*t - 2.6096304D-8*t*t*t
#ifdef tuan_frth
      ELSE IF( fluid_type == 2 ) THEN
 ! need to check the unit (Cecius or Kelvin)
!            fcond = 3.6100d0 + 1.517d-2 * tk - 1.741d-6 * tk**2
            fcond = 3.284+0.01617*tk - 0.000002305*tk**2   
      ELSE IF( fluid_type == 3 ) THEN
 ! modified by Tuan
            fcond = 9.2d0 + 0.011d0 * tk
 ! end modified 
      ELSE IF( fluid_type == 4 ) THEN
        ! need to check the unit (Cecius or Kelvin)
!    		fcond = 92.948d0 - 5.809d-2*(tk-273.15d0) + 1.1727d-5*(tk-273.15d0)**2
         fcond =  110.0-0.0648*tk+0.0000116*tk**2 ! SARAX perporty table

#endif

      ELSE
         WRITE(*,*) "Fluid Type: 0 for h2o, 1 for d2o"
         STOP "Fluid Type: 0 for h2o, 1 for d2o"
      END IF

      RETURN
      END FUNCTION fcond


!!!#ifdef siarhei_delete 
      FUNCTION fdens(t)
      use Inc_XS_File, only: if_th
!     ! use MOD_STEAMTABLE, only: get_water_density
      IMPLICIT NONE
      ! cubic poly for density as a func of temperature
      REAL(8) ::  fdens, t
      REAL(8) :: t_f                    !F
      REAL(8) :: pressure_psia          ! psia
      REAL(8) :: rhof                   ! lbm/ft^3
#ifdef tuan_frth
      REAL(8) :: tk

      tk = t + 273.15d0

#endif
#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fdens] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF (fluid_type == 0) THEN
         select case(if_th)
         case(1)
            ! cubic polynomial for density as a function of temperature, max err=0.0692%
            ! 15.5 Mpa,  280 < T < 340, rho in Kg/M^3, T in C
            fdens=5.9901166d+03+t*(-5.1618182d+01+t*(1.7541848d-01+t*(-2.0613054d-04)))
!            dummy_filler = 1 ! @$^ siarhei_plot 
         case(2)
            t_f = t * 1.8d0 + 32.0d0
            pressure_psia = Core_Pressure*14.5037738d0
            IF (t_f <= 450.0d0) THEN
                rhof = 65.43d0 - t_f * 0.03486d0
            ELSE
                rhof = -810.063d0 + (((-1.06824d-8 * t_f + 2.29413d-5) * t_f &
                    - 0.0184872d0) * t_f + 6.56174d0) * t_f
            END IF
            fdens = rhof*16.018453258153393d0 ! lbm/ft3 to kg/m^3
!            dummy_filler = 1 ! @$^ siarhei_plot 
         case(3)
            if(t>280.0d0) then
                fdens=5.9901166d+03+t*(-5.1618182d+01+t*(1.7541848d-01+t*(-2.0613054d-04)))
!            dummy_filler = 1 ! @$^ siarhei_plot 
            else
                fdens=-8.198d-21*t**8 + 9.589d-18*t**7 -5.009d-15*t**6 +1.517d-12*t**5 &
                    -3.077d-10*t**4 +4.535d-8*t**3-6.697d-6*t**2 -8.916d-6*t + 1.008d0
!            dummy_filler = 1 ! @$^ siarhei_plot 
                fdens=fdens*1.0d+3
!            dummy_filler = 1 ! @$^ siarhei_plot 
            endif
         case(4)
            fdens = 1d0
!            dummy_filler = 1 ! @$^ siarhei_plot 
           ! fdens=get_water_density(1,Core_Pressure*0.10d0,t+273.15d0)*1d+3 ! P in MPa, T in K
!            dummy_filler = 1 ! @$^ siarhei_plot 
         end select
      ELSEIF(fluid_type == 1) THEN
         fdens=-1.2368792d-04*t*t*t+9.5084035d-02*t*t-2.6208319d-01*t+3.4368724d+03
!            dummy_filler = 1 ! @$^ siarhei_plot 
#ifdef tuan_frth
      ELSEIF(fluid_type == 2) THEN
!            fdens = 11096.0d0 - 1.3236d0 * tk
            fdens =  11065.0-1.293*tk

      ELSEIF(fluid_type == 3) THEN
            fdens = 11367.0d0 - 1.1944d0 * tk
            !fdens =  11441.0d0 - 1.2795*tk

      ELSEIF(fluid_type == 4) THEN
            fdens = 9.50076d2 - 2.2976d-1*(tk-273.15d0) - 1.46049d-5*(tk-273.15d0)**2 + 5.63788d-9*(tk-273.15d0)**3
            !fdens = 1014.0d0-0.235d0*tk

#endif

      ELSE
         WRITE(*,*)'Fluid Type: 0 for h2o, 1 for d2o '
         STOP 'Fluid Type: 0 for h2o, 1 for d2o '
      END IF
      RETURN
      END FUNCTION fdens
!!!#endif 


      FUNCTION fdensh(h)

      IMPLICIT NONE

      ! cubic poly for density as a func of enthalpy
      REAL(8) :: fdensh, h, hk


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fdensh] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      hk=0.001d0*h
      IF(fluid_type == 0) THEN
         fdensh=1.4532039d+03                                              &
              +hk*(-1.1243975d+00+hk*(7.4502004d-04+hk*(-2.3216531d-07)))
      ELSEIF(fluid_type == 1) THEN
         fdensh=-2.5859687d-07*hk*hk*hk+3.8133349d-04*hk*hk                &
              -5.4411527d-01*hk+1.1525456d+03
#ifdef tuan_frth
      ELSEIF(fluid_type == 2) THEN
         fdensh = fdens( ftemp(h) )
!            dummy_filler = 1 ! @$^ siarhei_plot 

      ELSEIF(fluid_type == 3) THEN
         fdensh = fdens( ftemp(h) )
!            dummy_filler = 1 ! @$^ siarhei_plot 

      ELSEIF(fluid_type == 4) THEN
         fdensh = fdens( ftemp(h) )
!            dummy_filler = 1 ! @$^ siarhei_plot 

#endif
      ELSE
         WRITE(*,*)'Fluid Type: 0 for h2o, 1 for d2o '
         STOP 'Fluid Type: 0 for h2o, 1 for d2o '
      END IF

      RETURN
      END FUNCTION fdensh


!!!#ifdef siarhei_delete 
      FUNCTION fenthal(t)

      IMPLICIT NONE

      ! cubic poly for enthalpy as a func of temperature
      REAL(8) :: y
      REAL(8) :: fenthal, t
#ifdef tuan_frth
      REAL(8) :: tk

      tk = t + 273.15d0

#endif
#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fenthal] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF(fluid_type == 0) THEN
         y=-5.9301427d+03                                                  &
              +t*(6.5488800d+01+t*(-2.1237562d-01+t*(2.4941725d-04)))
      ELSEIF(fluid_type == 1) THEN
         y=1.2046906d-04*t*t*t-9.3015736d-02*t*t                           &
              +2.8651125d+01*t-2.5956445d+03
#ifdef tuan_frth
      ELSEIF(fluid_type == 2) THEN
!         y = (159.0d0 * tk - 1.3600d-2 * tk**2 + 7.120d-9/3.0d0 * tk**3)/1.d3
         y = 164.8 * (tk-398.0)-0.0197*(tk**2-398.0**2)+0.000004167*(tk**3-398.0**3) &
             + 456000*((tk**(-1)-398.0**(-1)))


      ELSEIF(fluid_type == 3) THEN
         y = (175.1d0 * tk - 4.961d-2/2.0d0 * tk**2 + 1.985d-5/3.0d0 * tk**3 &
          - 2.099d-9/4.0d0 * tk**4 + 1.524d6 / tk)/1.d3
       !  y= 176.2*(t-600.0)-0.024615*(t**2-600.0**2)+0.000005147*(t**3-600.0**3) + 1524000*((t**(-1)-600.0**(-1))) 

      ELSEIF(fluid_type == 4) THEN
         y = (-1.3792162d+5 +1.4370158d+3*(tk-273.15d0) -0.290302d0*&
         (tk-273.15d0)**2 +1.5427196d-4*(tk-273.15d0)**3)/1.d3
       !  y = 1658*(t-371.0) - 0.42395*(t**2-371.0**2) + 0.0001484667*(t**3-371.0**3) + 3001000*((t**(-1)-371.0**(-1)))

#endif
      ELSE
         WRITE(*,*)'Fluid Type: 0 for h2o, 1 for d2o '
         STOP 'Fluid Type: 0 for h2o, 1 for d2o '
      END IF
#ifdef tuan_fr

      fenthal=y
#else
      fenthal=y*1000d0
#endif
!            dummy_filler = 1 ! @$^ siarhei_plot 

      RETURN
      END FUNCTION fenthal
!!!#endif 


      FUNCTION fhcap(t)

      IMPLICIT NONE

      ! quartic poly for heat capacity as a func of temperature
      REAL(8) :: y
      REAL(8) :: fhcap, t
#ifdef tuan_frth
      REAL(8) :: tk


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fhcap] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      tk = t + 273.15d0

#endif

      IF ( fluid_type == 0 ) THEN
         y = 3.0455749D+3 - 4.0684599D+1*t + 2.0411250D-1*t*t - 4.5526705D-4*t*t*t + 3.8115453D-7*t*t*t*t
      ELSE IF ( fluid_type == 1 ) THEN
         y = 6.2739180D+2 - 8.9935937D+0*t + 4.8717545D-2*t*t - 1.1747742D-4*t*t*t + 1.0658597D-7*t*t*t*t
#ifdef tuan_frth
      ELSE IF ( fluid_type == 2 ) THEN
!         y = (159.0d0 - 2.720d-2 * tk + 7.120d-6 * tk**2)/1000.D0  
           y = (164.8-0.0394*tk+0.0000125*tk**2-456000*tk**(-2))   


      ELSE IF ( fluid_type == 3 ) THEN
          y = 175.1d0 - 4.961d-2 * tk + 1.985d-5 * tk**2 - 2.099d-9 * tk**3 - 1.524d6 * tk**(-2)
!         y = (175.1d0 - 4.961d-2 * tk + 1.985d-5 * tk**2 - 2.099d-9 * tk**3 &
!		 - 1.524d6 * tk**(-2))/1000.D0
      ELSE IF ( fluid_type == 4 ) THEN
         y = (1.4370158d+3 -0.290302d0*2.0d0*tk +1.5427196d-4*3.0d0*tk**2)/1000.D0
!y = 1.4370158d+3 -0.290302d0*2.0d0*(t-273.15d0) +1.5427196d-4*3.0d0*(t-273.15d0)**2
! y = -3001000*t**(-2)+1658-0.8479*t+0.0004454*t**2
#endif

      ELSE
         WRITE(*,*) "Fluid Type: 0 for h2o, 1 for d2o"
         STOP "Fluid Type: 0 for h2o, 1 for d2o"
      END IF
#ifdef tuan_fr
      fhcap = y 
#else
      fhcap = y * 1000.D0
#endif

      RETURN
      END FUNCTION fhcap



      FUNCTION fhtcoef( t, deq, rhou )

      IMPLICIT NONE

      ! wall-to-coolant heat transfer coeffcient
      REAL(8) :: fhtcoef, deq, rhou, t
      REAL(8) :: k, mu, cp, pr, re


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fhtcoef] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      k  = fcond(t)
      mu = fvisco(t)
      cp = fhcap(t)
      pr = cp * mu / k
      re = deq * rhou / mu

      fhtcoef = 0.023D0*k/deq*(pr**0.4D0)*(re**0.8D0)
! write(*,*) "This heat transfer coeficient is for LBE"
#ifdef tuan_frth
! Lyon Martinelli correlation from TRACE
!       fhtcoef = (4.8D0+0.025D08*(pr*re)**0.8d0)*k/deq
! need to modify, pitch-to-diameter ratio, v
!      fhtcoef = 0.047*(1-dexp(-3.8D0*(1.130-1)))*((pr*re)**0.77+250)*k/deq
! Ushakov  correlation
      fhtcoef = 7.55*1.13D0-20D0/(1.13D0**13)+(0.041D0*(pr*re)**(0.56+0.19*1.13D0)/1.13D0**2) *k/deq
! Graber correlation
!      fhtcoef = (0.25D0+6.2*1.13D0+(0.032*1.13D0-0.007D0)*(pr*re)**(0.8D0-0.024d0*1.130D0)) *k/deq
! Subbotin correlation
!      fhtcoef = (0.58D0 *(2*sqrt(3.0)*1.130D0**2/3.141592-1)**0.55)*(pr*re)**0.45 *k/deq
!    new
!      fhtcoef = (2.75+0.02*(pr*re)**0.8) *k/deq
#endif
      RETURN
      END FUNCTION fhtcoef

      FUNCTION ftemp(h)

      IMPLICIT NONE

      ! cubic poly for temperature as a func of enthalpy
      REAL(8) :: ftemp, h, hk


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [ftemp] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      hk=0.001*h
      IF(fluid_type == 0) THEN
         ftemp=1.4851739d+02                                               &
              +hk*(-1.2764991d-01+hk*(3.0781294d-04+hk*(-9.5429959d-08)))
      ELSEIF(fluid_type == 1) THEN
         ftemp=-9.7160887d-08*hk*hk*hk+1.4489550d-04*hk*hk                 &
              +1.4472728d-01*hk+1.2268248d+02
#ifdef tuan_frth
      ELSEIF(fluid_type == 2) THEN
         ftemp = 450.d0
         ftemp = f_lm_temp_search(h, ftemp)
      ELSEIF(fluid_type == 3) THEN
         ftemp = 500.d0
         ftemp = f_lm_temp_search(h, ftemp)
      ELSEIF(fluid_type == 4) THEN
         ftemp = 500.d0
         ftemp = f_lm_temp_search(h, ftemp)

#endif
      ELSE
         WRITE(*,*)'Fluid Type: 0 for h2o, 1 for d2o '
         STOP 'Fluid Type: 0 for h2o, 1 for d2o '
      END IF

      RETURN
      END FUNCTION ftemp


      FUNCTION fvisco(t)
      IMPLICIT NONE

      ! cubic poly for viscosity as a func of temperature
      REAL(8) :: fvisco, t
#ifdef tuan_frth
      REAL(8) :: tk


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fvisco] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      tk = t + 273.15d0

#endif

      IF( fluid_type == 0 ) THEN
         fvisco = 9.0836878D-4 - 7.4542195D-6*t + 2.3658072D-8*t*t - 2.6398601D-11*t*t*t
      ELSE IF( fluid_type == 1 ) THEN
         fvisco = 7.8299421D-4 - 6.1211876D-6*t + 1.9507509D-8*t*t - 2.2490154D-11*t*t*t
#ifdef tuan_frth
      ELSE IF( fluid_type == 2 ) THEN
         fvisco = 4.940d-4    * exp(754.1d0/tk)

      ELSE IF( fluid_type == 3 ) THEN
         fvisco = 4.550d-4    * exp(1069.0d0/tk)

      ELSE IF( fluid_type == 4 ) THEN
         fvisco = 3.241903d-3 * exp(508.07d0/tk-0.4925d0*log(tk))


#endif
      ELSE
         WRITE(*,*) "Fluid Type: 0 for h2o, 1 for d2o"
         STOP "Fluid Type: 0 for h2o, 1 for d2o"
      END IF

      RETURN
      END FUNCTION fvisco


!!!#ifdef siarhei_delete 
      FUNCTION findh(t,h)

      IMPLICIT NONE

      ! find enthalpy for given temperature, input h is the initial guess
      ! use fenthal instead of this unless it is really necessary
      INTEGER :: i
      REAL(8) :: t,h,x1,x2,xn,y1,y2,findh,slope

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [findh] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      x1=1.01d0*h
      x2=h
      y1=ftemp(x1)
      y2=ftemp(x2)
      DO i=1,20
         IF(ABS(x2-x1).LT.0.01d0) EXIT
         slope=(y2-y1)/(x2-x1)
         xn=x2+(t-y2)/slope
         x1=x2
         y1=y2
         x2=xn
         y2=ftemp(x2)
      END DO
      findh=xn
!            dummy_filler = 1 ! @$^ siarhei_plot 

      RETURN
      END FUNCTION findh
!!!#endif 


!!!#ifdef siarhei_delete 
      FUNCTION fdens_slb(p,t,opt)
!     ! use MOD_STEAMTABLE, only: get_water_density
      IMPLICIT NONE
      REAL(8) :: fdens_slb, t, p
      integer :: opt

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fdens_slb] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      select case(opt)
      case(0) !th1d
          fdens_slb=5.9901166d+03+t*(-5.1618182d+01+t*(1.7541848d-01+t*(-2.0613054d-04)))
!            dummy_filler = 1 ! @$^ siarhei_plot 
      case(1) !ctf
         fdens_slb = 1d0
!            dummy_filler = 1 ! @$^ siarhei_plot 
         ! fdens_slb=get_water_density(1, p*0.10d0, t+273.15d0)*1.0d+3 ! P in MPa, T in K
!            dummy_filler = 1 ! @$^ siarhei_plot 
      end select

      RETURN
      END FUNCTION fdens_slb
!!!#endif 


!!!#ifdef siarhei_delete 
      real(8) function fenthrise(ftemp, ftref, fbu, OMratio, facmolt)
      !function to compute enthalpy rise from reference state(initial state)
      ! input
      !    ftref  : fuel temperature in [K] at reference state
      !    ttemp  : fuel temperature in [K] at current state
      !    fbu    : fuel burnup in [MWd/MTU]
      !    OMRatio: ratio of fuel oxygen atom to U and Pu atom
      !    facmolt: fraction of molten fuel
      ! output
      !    fenthrise: enthalpy rise relative to reference state in [cal/g]
      real(8), intent(in) :: ftemp
      real(8), intent(in) :: ftref
      real(8), intent(in) :: fbu
      real(8), intent(in),optional :: OMRatio
      real(8), intent(in),optional :: facmolt
      real(8) :: omrate = 2.0d0
      real(8) :: frcmlt = 0.0d0
      real(8) :: fenth_cur
      real(8) :: fenth_ref

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fenthrise] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(present(OMRatio)) omrate = OMratio
      if(present(facmolt)) frcmlt = facmolt

      fenth_cur = fenthalpy(ftemp=ftemp, fbu=fbu, OMratio=omrate, facmolt=frcmlt)
!            dummy_filler = 1 ! @$^ siarhei_plot 
      fenth_ref = fenthalpy(ftemp=ftref, fbu=fbu, OMratio=omrate, facmolt=frcmlt)
!            dummy_filler = 1 ! @$^ siarhei_plot 
      fenthrise = fenth_cur - fenth_ref
!            dummy_filler = 1 ! @$^ siarhei_plot 

      end function fenthrise
!!!#endif 


!!!#ifdef siarhei_delete 
      real(8) function fenthalpy(ftemp, fbu, OMratio, facmolt)
      ! function to compute fuel enthalpy relative to 0 K
      ! function from fraptran, data from matpro-11
      ! input
      !    ttemp  : fuel temperature in [K] at current state
      !    fbu    : fuel burnup in [MWd/MTU]
      !    OMRatio: ratio of fuel oxygen atom to U and Pu atom
      !    facmolt: fraction of molten fuel
      ! output
      !    fenthrise: enthalpy rise relative to reference state in [cal/g]
      real(8), intent(in) :: ftemp
      real(8), intent(in) :: fbu
      real(8), intent(in),optional :: OMRatio
      real(8), intent(in),optional :: facmolt
      real(8), parameter :: hefus = 27.4d4
      real(8), parameter :: fcpmol = 503.0d0
      real(8), parameter :: cg2jk = 4186.8d0 !unit conversion cal/g to k/j
      ! The following statements contain constants from matpro-11
      real(8), parameter :: c1u = 296.7d0
      real(8), parameter :: c2u = 2.43d-2
      real(8), parameter :: c3u = 8.745d7
      real(8), parameter :: thu = 535.285d0
      real(8), parameter :: edu = 1.577d5
      real(8), parameter :: c1pu = 347.4d0
      real(8), parameter :: c2pu = 3.95d-4
      real(8), parameter :: c3pu = 3.860d7
      real(8), parameter :: thpu = 571.0d0
      real(8), parameter :: edpu = 1.967d5
      real(8) :: tx
      real(8) :: omrate = 2.0d0
      real(8) :: frcmlt = 0.0d0
      real(8) :: ftmelt ! melting temperature func. of burnup


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fenthalpy] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if(present(OMRatio)) omrate = OMratio
      if(present(facmolt)) frcmlt = facmolt
      ftmelt = 3113.0d0 - 5d0*fbu/10000d0

      tx = MIN(ftemp,ftmelt)

      fenthalpy = cpdt(c1u, thu, c2u, omrate, edu, tx, c3u)
!            dummy_filler = 1 ! @$^ siarhei_plot 
      if(ftemp .gt. ftmelt - 2d0) then
         fenthalpy = fenthalpy + hefus * frcmlt
!            dummy_filler = 1 ! @$^ siarhei_plot 
         if(ftemp .gt. ftmelt + 2d0) then
            fenthalpy = fenthalpy + (ftemp - ftmelt) * frcmlt
!            dummy_filler = 1 ! @$^ siarhei_plot 
         endif
      endif

      fenthalpy = fenthalpy/cg2jk
!            dummy_filler = 1 ! @$^ siarhei_plot 

      end function fenthalpy
!!!#endif 


!!!#ifdef siarhei_delete 
      real(8) FUNCTION cpdt (c1, th, c2, otm, ed, t, c3)
      ! integral of the fuel specific heat with respect to temperature
      real(8), intent(in) :: c1, th, c2, otm, ed, t, c3

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [cpdt] in Mod_THfunc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      cpdt = c1 * th * (1.0d0 / (exp(th / t) - 1.0d0)) + c2 * t ** 2 / 2.0d0 + &
        &    c3 * otm * exp(-ed / (t * 8.314d0)) / 2.0d0
!            dummy_filler = 1 ! @$^ siarhei_plot 

      END FUNCTION cpdt
!!!#endif 


      END MODULE Mod_THfunc
