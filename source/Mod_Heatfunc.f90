#ifdef siarhei_delete

module mod_heatfunc

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

implicit none

contains

! From CTF old version : JS (180512)
! contains pure procedures for calculation of heat transfer coefficients

function calc_qchf_Biasi (de,ireg,alfl,alfv,g,h55,hspv,p,xa) result (qchf)
   !---------------------------------------------------------------------------
   ! COBRA-TF standard CHF correlation (Biasi)
   ! -------------------------------------------
   ! --> for forced convection DNB including ramps to
   !     - annular film dryout (for high void fractions)
   !     - pool boiling DNB (for low mass flow rates)
   ! --> using SI units: g,cm,s
   !
   ! Arguments ----------------------------------------------------------------
   !   de   - Hydraulid diameter     [ft]
   !   alfl - Liquid volume fraction [-]
   !   alfv - Vapor volume fraction  [-]
   !   g    - Total mass flux        [lbm/{hr ft^2}]
   !   h55  - (See "fluidprop")
   !   hspv - Single phase vapor HTC [btu/{hr ft^2 F}]
   !   p    - Pressure               [psi]
   !   xa   - Quality                [-]
   !   ireg - Flow regime
   real(8), intent(in) :: de,alfl,alfv,g,h55,hspv,p,xa
   integer, intent(in) :: ireg
   real(8)             :: qchf
   ! Local variables ----------------------------------------------------------
   real(8) :: dm,en,fp,gm,gm6,gmz,hp,pm,qchf1,qchf2,qchfz,ramp
   !-------------------------------------------------------  --------------------

   ! Convert from US units to SI units
   en = 0.4d0
   gm = 0.0001356d0 * g
   ! In case of no post CHF conditions (flow regime <> 6)
   if (ireg /= 6) then
      ! In case of annular film dryout
      if (alfl < 0.1d0) then
         qchf = h55
         if (alfl < 0.01d0) then
            qchf = h55 * max(0.2d0,200d0*(alfl-0.005d0))
         end if
      else
         gmz = max(30d0,gm)
         dm  = 30.48d0 * de
         pm  = max(1.38d0,0.06897d0*p)
         if (dm < 1.0d0)  then
            en = 0.6d0
         endif
         ! Low quality equation
         fp = 0.7249d0 + 0.099d0 * pm * exp(-0.032d0*pm)
         gm6 = 1d0 / gmz**0.1667d0
         qchf1 = 2.764e+07 * gm6 * (1.468d0*fp*gm6-xa) / dm**en
         ! High quality equation
         hp = -1.159d0 + 0.149d0 * pm * exp(-0.019d0*pm) &
              + 8.99d0 * pm / (10d0 + pm*pm)
         qchf2 = 15.048e+07 * hp * (1d0-xa) / (dm**en * gmz**0.6d0)
         ! Convert from SI units to US units
         qchf1 = qchf1 * 3.171e-01
         qchf2 = qchf2 * 3.171e-01
         ! Use the maximum of the low quality equation, of the high
         ! quality equation, or of the fix value 90000
         qchf = max(qchf1,qchf2,90000d0)
         ! Low flowrate - ramp to Zuber CHF (Pool boiling DNB)
         qchfz = max(0.2d0,(1d0-alfv)) * h55
         if (gm < 30d0) then
            qchf = max(qchf,qchfz)
         endif
         ! High void fraction - ramp to annular film dryout
         if (alfl < 0.15d0) then
             qchf = 20d0*(qchf*(alfl-0.1) + h55*(0.15d0-alfl))
         endif
      end if
   else
      ! In case of post CHF conditions (flow regime 6)
      ! Low flowrate - use Zuber CHF
      qchf = h55
   end if
   ! Ramp to single phase vapor
   ramp = min(1d0,100d0*(1d0-alfv))
   qchf = ramp * qchf + (1d0-ramp) * 50d0 * hspv

end function calc_qchf_Biasi

function calc_qchf_W3 (de,hin,hf,p,g,xaeq) result (qchf)
   !---------------------------------------------------------------------------
   ! W3 correlation for forced convection DNB
   !
   ! Arguments ----------------------------------------------------------------
   real(8), intent(in) :: de    ! Hydraulic diameter  [ft]
   real(8), intent(in) :: hin   ! Inlet enthalpy      [btu/lbm]
   real(8), intent(in) :: hf    ! Saturation enthalpy [btu/lmb]
   real(8), intent(in) :: p     ! Pressure            [psi]
   real(8), intent(in) :: g     ! Total mass flux     [lbm/{hr ft^2}]
   real(8), intent(in) :: xaeq  ! Equilibrium quality [-]
   real(8)             :: qchf
   ! Local variables ----------------------------------------------------------
   real(8) :: c4,c5,c6,c6a,c7,c8,c9,qdnb
   !---------------------------------------------------------------------------

   c4  = 2.022d0 - 0.4302d0 * (p/1000d0)
   c5  = 0.1722d0 - 0.0984d0 * (p/1000d0)
   c6  = EXP(18.177d0*xaeq-(4.129d0*xaeq*(p/1000d0)))
   c6a = 1.157d0 - 0.869d0 * xaeq
   c7  = (0.1484d0 - 1.596d0 * xaeq + 0.1729d0 &
         * xaeq * abs(xaeq))*(g/1.E+06) + 1.037d0
   c8  = 0.2664d0 + (0.8357d0 * exp(-(3.151d0*12d0*de)))
   c9  = 0.8258d0 + 0.000794d0 * (hf-hin)
   qdnb = 1.E+06 * ((c4+c5*c6) * c6a * c7 * c8 * c9)
   ! For now the Tong F Factor is assumed to be constant (1.1)
   qchf = qdnb / 1.1d0 ! [BTU/hr-ft^2]

end function calc_qchf_W3
!js+
!js+ real function calc_Nusselt_spv (re,pr) result(nu)
!js+    use DittusBoelter
!js+    !---------------------------------------------------------------------------
!js+    ! Calculates Nusselt number for single-phase vapor
!js+    !
!js+    ! Arguments ----------------------------------------------------------------
!js+    real, intent(in) :: re, pr ! Reynolds' number and Prandtl number
!js+    !---------------------------------------------------------------------------
!js+    if (re <= 25200.0) then
!js+       ! Maximum of:
!js+       ! - Nu = 10 for laminar flow
!js+       ! - Wong and Hochreiter FLECHT-SEASET correlation
!js+       nu = max(10.0, 0.0797 * re**0.6774 * pr**(1./3.))
!js+    else
!js+       ! Dittus-Boelter correlation
!js+       ! VUQ PIRT
!js+       ! nu = 0.023 * re**0.8 * pr**0.4
!js+       nu = k_db_1 * re**k_db_2 * pr**k_db_3
!js+    end if
!js+ end function calc_Nusselt_spv
!js+
!js+ real function calc_spvmd (area,de,cgam,dsd,sdn,alfe,ddrop,ffric,rhol,rhov,uv,&
!js+       vcont,vrel)
!js+    !---------------------------------------------------------------------------
!js+    ! Calculates multiplier for single-phase vapor heat transfer coefficient to
!js+    ! take into account enhanced heat transfer due to dispersed droplets.
!js+    !
!js+    ! Modules ------------------------------------------------------------------
!js+
!js+    ! Arguments ----------------------------------------------------------------
!js+    !  area  - Continuity cell flow area          [ft^2]
!js+    !  de    - Hydraulic diameter                 [ft]
!js+    !  cgam  - Entrained velocity                 [ft/s]
!js+    !  dsd   - Diameter of small droplets         [ft]
!js+    !  sdn   - Small droplet number density       [1/ft^3]
!js+    !  alfe  - Entrained liquid volume fraction   [-]
!js+    !  ddrop - Droplet diameter                   [ft]
!js+    !  ffric - Friction factor                    [-]
!js+    !  rhol  - Liquid density                     [lbm/ft^3]
!js+    !  rhov  - Vapor density                      [lbm/ft^3]
!js+    !  uv    - Vapor viscosity                    [lbm/{hr ft}]
!js+    !  vcont - Vapor velocity                     [ft/s]
!js+    !  vrel  - Droplet velocity relative to vapor [ft/s]
!js+    real, intent(in) :: area,de,cgam,dsd,sdn,alfe,ddrop,ffric,rhol,rhov,uv,&
!js+          vcont,vrel
!js+    ! Local variables ----------------------------------------------------------
!js+    real :: aesd,cdrop,cdvr2,red,spvms,vrmax,vsd
!js+    !---------------------------------------------------------------------------
!js+
!js+    ! Compute pressure effect on drop repulsion
!js+    ! (coefp is never used in the code !)
!js+    !coefp = 96.56 * (sigma*ddrop*rhof)**0.25 / sqrt(uf)
!js+    red = max(2.0,t_hr_s*rhov*vrel*ddrop/uv)
!js+    cdrop = (24./red) * (1.0 + 0.15*red**0.687)
!js+    ! Add effect from small drop field
!js+    spvms = 0.0
!js+    if (dsd >= 1.0e-05) then
!js+       cdvr2 = 1.333 * (rhol-rhov) * 32.2 * dsd / rhov
!js+       vrmax = sqrt(cdvr2/0.45)
!js+       vsd = max(cgam,vcont-vrmax,1.0)
!js+       aesd = sdn * PI * dsd**3 / (6. * vsd * area)
!js+       spvms = 0.75 * aesd * de * cdvr2 / (dsd * ffric * vcont**2)
!js+    end if
!js+    calc_spvmd = sqrt(1.0 + 0.75 * alfe * de * cdrop / &
!js+               (ddrop*ffric) *(vrel/vcont)**2 + spvms)
!js+    calc_spvmd = min(10.,calc_spvmd)
!js+
!js+ end function calc_spvmd
!js+
!js+ subroutine calc_Tchf (alfl,alfv,h11,hf,hl,hspl,hspv,p,Tf,Tl,qchf,Tchf)
!js+    !---------------------------------------------------------------------------
!js+    ! Calculation of the CHF temperature 'tchf'
!js+    ! (iterative procedure to find the wall temperature at which the
!js+    !  heat flux according to the Chen nucleate boiling correlation is
!js+    !  equal to the critical heat flux)
!js+    !
!js+    ! Arguments ----------------------------------------------------------------
!js+    real, intent(in)  :: alfl ! Liquid void at the top of this cell
!js+    real, intent(in)  :: alfv ! Vapor void at the top of this cell
!js+    real, intent(in)  :: h11  ! Fluid properties table needed by the chen
!js+                              ! correlation
!js+    real, intent(in)  :: hf   ! Liquid saturation enthalpy
!js+    real, intent(in)  :: hl   ! Liquid enthalpy
!js+    real, intent(in)  :: hspl ! Liquid single phase HTC [BTU/hr-ft^2-F]
!js+    real, intent(in)  :: hspv ! Vapor single phase HTC [BTU/hr-ft^2-F]
!js+    real, intent(in)  :: p    ! Pressure [psia]
!js+    real, intent(in)  :: tf   ! Saturation temperature [F]
!js+    real, intent(in)  :: tl   ! Liquid temperature [F]
!js+    real, intent(in out) :: qchf ! The critical heat flux [BTU/hr-ft^2]
!js+    real, intent(out) :: Tchf ! The temperature at CHF [F]
!js+    ! Local variables ----------------------------------------------------------
!js+    integer :: niter
!js+    real    :: dhtc,dTf,dTl,eps,hnb,htcl1,htcv1,qcmin,qcmax,qnb,qnew,&
!js+               Tcmin,Tcmax
!js+    real, parameter :: convergence_criteria=0.01
!js+    !---------------------------------------------------------------------------
!js+
!js+    ! Use liquid and vapor modifiers
!js+    htcl1 = hspl
!js+    htcv1 = hspv
!js+    if (alfv > 0.999) htcl1 = htcl1 * max(0.0,(0.9999-alfv)/0.0009)
!js+    if (alfv < 0.05)  htcv1 = htcv1 * max(0.0,(alfv-0.01)/0.04)
!js+
!js+
!js+    ! Calculate the boundaries on the CHF temperature
!js+    Tcmin = Tf + 20.
!js+    Tcmax = min(Tf+200.0,705.3)
!js+
!js+    ! Now calculate the corresponding boundaries on the critical heat flux
!js+    ! Corresponding to Tmin
!js+    dTl   = Tcmin - Tl
!js+    dTf   = Tcmin - Tf
!js+    if (ihtc == 0) then
!js+       call htc_chen(dtf,h11,alfl,hnb,dhtc,qnb)
!js+    else
!js+       call htc_thom(p,dTf,dTl,hl,hf,alfl,qnb,hnb,dhtc)
!js+    end if
!js+    qnew  = qnb + htcl1 * dTl + htcv1 * dTf
!js+    qcmin = qnew
!js+
!js+    ! Corresponding to Tmax
!js+    dTl   = Tcmax - Tl
!js+    dTf   = Tcmax - Tf
!js+    if (ihtc == 0) then
!js+       call htc_chen(dtf,h11,alfl,hnb,dhtc,qnb)
!js+    else
!js+       call htc_thom(p,dTf,dTl,hl,hf,alfl,qnb,hnb,dhtc)
!js+    end if
!js+    qnew  = qnb + htcl1 * dTl + htcv1 * dTf
!js+    qcmax = qnew
!js+
!js+    ! Before we even search for Tchf, check to see if the calculated
!js+    ! q"chf is outside of the boundaries we have just calculated
!js+    if (qchf >= qcmax) then
!js+       qchf = qcmax
!js+       Tchf = Tcmax
!js+       return
!js+    else if (qchf <= qcmin) then
!js+       qchf = qcmin
!js+       Tchf = Tcmin
!js+       return
!js+    end if
!js+
!js+    ! Iteratively find Tchf, initializing it to the midpoint of Tmin and Tmax
!js+    Tchf = 0.5*(Tcmin+Tcmax)
!js+    niter = 0
!js+    ! Iterate to find CHF temperature
!js+    do niter = 1, 100
!js+       dTl   = Tchf - Tl
!js+       dTf   = Tchf - Tf
!js+       if (ihtc == 0) then
!js+          call htc_chen(dtf,h11,alfl,hnb,dhtc,qnb)
!js+       else
!js+          call htc_thom(p,dTf,dTl,hl,hf,alfl,qnb,hnb,dhtc)
!js+       end if
!js+       qnew  = qnb + htcl1 * dTl + htcv1 * dTf
!js+
!js+       if (qchf >= qnew)  Tcmin = Tchf
!js+       if (qchf <= qnew)  Tcmax = Tchf
!js+
!js+       ! Finished if error sufficiently small
!js+       eps = abs(qchf-qnew) / qchf
!js+       if (eps < convergence_criteria) return
!js+
!js+       tchf = 0.5 * (Tcmin + Tcmax)
!js+    end do
!js+
!js+    ! If convergence fails -----------------------------------------------------
!js+    if (ittyout == 0) then
!js+       write(itty,*) 'Too many iteration to calculate critical',&
!js+             'heat flux temperature Tchf !! Residual EPS = ',eps
!js+    end if
!js+    msg = ""
!js+    msg(1) = 'Too many iteration to calculate critical heat '
!js+    msg(2) = 'flux temperature Tchf !! Residual EPS = '
!js+    write(msg(3),*) eps
!js+    call print_error
!js+
!js+ end subroutine calc_Tchf
!js+
!js+ subroutine calc_view_factors (de,alfe,alfl,ddrop,p,Tv,fwd,fwg)
!js+    ! Arguments ----------------------------------------------------------------
!js+    real, intent(in)  :: de,    &! Hydraulic diameter [ft]
!js+                         alfe,  &! Entrained droplet volume fraction
!js+                         alfl,  &! Continuous liquid volume fraction
!js+                         ddrop, &! Droplet diameter   [ft]
!js+                         p,     &! Pressure           [psi]
!js+                         Tv      ! Vapor temperature  [F]
!js+    real, intent(out) :: fwd,   &! Wall to droplet view factor
!js+                         fwg     ! Wall to vapor view factor
!js+    ! Local variables ----------------------------------------------------------
!js+    real :: af1,av1,ef,eg,egef,ew,Tvr,r1,r2,r3,r4
!js+    !---------------------------------------------------------------------------
!js+
!js+    af1  = 1.11 * (alfe/ddrop + alfl/de)
!js+    ef   = 1. - exp(-0.85*af1*de)
!js+    ef   = max(0.001,min(0.75,ef))
!js+    Tvr  = Tv + 460.
!js+    av1  = (p/14.7) * (5.6*(1000./Tvr)**2 - 0.3*(1000./Tvr)**4)
!js+    eg   = 1.0 - exp(-0.85*av1*de)
!js+    eg   = max(0.001,min(0.75,eg))
!js+    ew   = 0.6
!js+    egef = 1.0 / (1.0 - eg*ef)
!js+    r1   = egef * (1.0-eg) / eg
!js+    r2   = egef * (1.0-ef) / ef
!js+    r3   = egef + (1.0-ew) / ew
!js+    r4   = 1.0 + r3/r1 + r3/r2
!js+    fwg  = 1.0 / (r1*r4)
!js+    fwd  = 1.0 / (r2*r4)
!js+    fwg  = 1.713e-09 * fwg
!js+    fwd  = min(7.71e-10,1.713e-09*fwd)
!js+
!js+ end subroutine calc_view_factors
!js+
!js+ subroutine film_props (de,hg,p,reu,T,hsp,k,rho,u)
!js+    !---------------------------------------------------------------------------
!js+    ! Calculates vapor properties at film temperature (twall + tvap) / 2
!js+    ! (Created from part 3 of the "prop" subroutine)
!js+    !
!js+    ! Arguments ----------------------------------------------------------------
!js+    real, intent(in)    :: de,  &! Hydraulic diameter   [ft]
!js+                           hg,  &! Sat. vapor enthalpy  [btu/lbm]
!js+                           p,   &! Pressure             [psi]
!js+                           reu   ! Reynolds number of vapor divided by dynamic
!js+                                 ! viscosity            [{hr ft}/lbm]
!js+    real, intent(inout) :: T     ! Film temperature     [F]
!js+    real, intent(out)   :: hsp, &! Single-phase htc     [btu/{hr ft^2 F}]
!js+                           k,   &! Thermal conductivity [btu/{hr ft F}]
!js+                           rho, &! Density              [lbm/ft^3]
!js+                           u     ! Dynamic viscosity    [lbm/{hr ft}]
!js+    ! Local variables ----------------------------------------------------------
!js+    real :: cp, hv, nu, pr, re
!js+    !---------------------------------------------------------------------------
!js+
!js+    !jwm288 - The calculation is made in two stages (t->h & h->t) as artifact
!js+    !         of how CTF's former hgas,tgas, and transp subroutines functioned.
!js+    !         This should be simplified, but will slightly change results.
!js+    call vapor_props(p=p,T_in=t,h=hv,rho=rho)
!js+    hv = max(hv,hg)
!js+    call vapor_props(p=p,h_in=hv,rho_in=rho,cp=cp,k=k,T=T,u=u)
!js+
!js+    ! Prandtl number
!js+    pr = cp * u / k
!js+    ! Reynolds number
!js+    re = reu / u
!js+    ! Calculate Nusselt number
!js+    nu = calc_Nusselt_spv(re,pr)
!js+    ! Calculate heat transfer coefficient
!js+    hsp = nu * k / de
!js+
!js+ end subroutine film_props
!js+
!js+ subroutine hot_wall(nrad,alfa,alfl,coefd,de,fwd,fwg,gle,h33,h55,hfg,&
!js+    hg,pl,qchf,scbmod,Tchf,Tf,Tl,Tg,Tmin,Twall,hspv,htcb,htcl,htclr,&
!js+    htcv,htcvr,imode)
!js+    !---------------------------------------------------------------------------
!js+    ! Calculate HTCs when the CHF temperature has been exceeded
!js+    !
!js+    ! Modules ------------------------------------------------------------------
!js+    use fluidprops, only: vapor_props
!js+    ! Arguments ----------------------------------------------------------------
!js+    integer, intent(in)  :: nrad
!js+    real,    intent(in)  :: alfa,alfl,coefd,de,fwd,fwg,gle,h33,h55,hfg,hg,pl,&
!js+                            qchf,scbmod,Tchf,Tf,Tl,Tg,Tmin,Twall,hspv
!js+    real,    intent(out) :: htcb,htcl,htclr,htcv,htcvr
!js+    integer, intent(out) :: imode
!js+    ! Local variables ----------------------------------------------------------
!js+    real :: dtf,dtg,dtl,effn,fwet,hfb,hfbl,qdffb,qiafb,qiliq,qivap,qsde,qtb,&
!js+            qtbd,ramp,rampl,rampq,rfb,tcfb,Tf4,Tfb,Tg4,tw4,ufb
!js+    !---------------------------------------------------------------------------
!js+
!js+    ! Initializations
!js+    htcl = 0.0
!js+    htcb = 0.0
!js+    dTl = twall - Tl ! [F]
!js+    dTf = twall - Tf ! [F]
!js+    dTg = twall - Tg ! [F]
!js+    if (dTl == 0.0)  dTl = 1.0e-6
!js+    if (dTf == 0.0)  dTf = 1.0e-6
!js+    if (dTg == 0.0)  dTg = 1.0e-6
!js+    htclr = 0.0
!js+    htcvr = 0.0
!js+
!js+    ! Vapor heat transfer
!js+    htcv = hspv
!js+    imode = 2
!js+
!js+    ! Film boiling heat transfer
!js+    ! Drop deposition heat transfer
!js+    qsde = min(hfg*coefd*gle,qchf)
!js+    effn = exp(1.-((Twall+460.)/(Tf+460.))**2)
!js+    htcb = qsde * effn / dTf
!js+    ! Wall - drop radiation (Sun, Dix, Tien)
!js+    if (nrad == 0) then
!js+       Tw4 = (Twall+460.)**4
!js+       Tf4 = (Tf+460.)**4
!js+       htclr = fwd * (tw4 - tf4)
!js+       htcl = htclr / dTl
!js+       ! Wall - steam radiation
!js+       Tg4 = (Tg + 460.)**4
!js+       htcvr = fwg * (Tw4 - Tg4)
!js+       htcv = htcv + htcvr / dTg
!js+    endif
!js+    imode = 9
!js+
!js+    ! Check alfa (vapor volume fraction = void fraction)
!js+    if (alfa <= 0.95) then
!js+       ! Dispersed flow film boiling (dffb) heat flux
!js+       qdffb = max(0.0,htcv*dTg)
!js+       ! Wall - core radiation
!js+       ! Inverted annular film boiling (iafb) heat flux
!js+       ! Evaluate vapor properties for Bromley correlation.
!js+       Tfb = 0.5 * (Twall + Tf)
!js+       call vapor_props(p=pl,T_in=Tfb,h=hfb,k=tcfb,rho=rfb,u=ufb)
!js+       hfbl = max(hfb,hg)
!js+       hfbl = h33 * Tcfb**0.75 * (rfb/(ufb*dTf*de))**0.25 * de**0.172
!js+       qiafb = hfbl * dTf
!js+       qtb = alfl * h55 * effn
!js+
!js+       ! Use maximum
!js+       if ((qdffb+htcb*dTf) <= (qiafb+qTb)) then
!js+          qiafb = max(0.0,qiafb-qdffb)
!js+          qiliq = qiafb
!js+          qivap = 0.0
!js+          if (dtg >= 0.0) then
!js+             rampq = 1.11 * alfa * sqrt(dTg/dTf)
!js+             qiliq = (1.0 - rampq) * qiafb
!js+             qivap = rampq * qiafb
!js+          endif
!js+
!js+          ! Check alfa (vapor volume fraction = void fraction)
!js+          if (alfa >= 0.4) then
!js+             ! Ramp between inverted annular film boiling (iafb) and
!js+             ! dispersed flow film boiling (dffb)
!js+             ! (The heat flux is linearely interpolated between the
!js+             !  values for iafb and dffb with void fraction.)
!js+             ramp = (0.95 - alfa) / 0.5
!js+             htcb = max(htcb,ramp*qTb/dTf)
!js+             htcb = htcb + ramp * qiliq / dTf
!js+             htcv = htcv + ramp * qivap / dTg
!js+             imode = 8
!js+          else
!js+             ! Inverted annular film boiling (modified Bromley correl.)
!js+             htcb = (qtb + qiliq) / dTf
!js+             htcv = htcv + qivap / dTg
!js+             imode = 7
!js+          endif
!js+       endif
!js+    endif
!js+
!js+
!js+    ! Check if wall temperature exceeds minimum film boiling
!js+    ! temperature
!js+    if (Twall > (Tmin+1.0)) then
!js+       return
!js+    end if
!js+
!js+    ! Transition boiling
!js+    ! (based on Bjornards fraction of wall wettable model)
!js+    fwet = ((Twall-Tmin) / (Tchf-Tmin))**2
!js+    ! This ramp was added to simulate dryout.
!js+    rampl = MIN(1.0,MAX(0.0,400.0*(alfl-0.0025)))
!js+    qtbd = fwet * qchf * rampl
!js+    htcb = scbmod * qtbd / dTf
!js+    htcl = htcl + (1.0-scbmod) * qtbd / max(1.,dTl)
!js+    imode = 6
!js+
!js+ end subroutine hot_wall

end module mod_heatfunc

#endif
