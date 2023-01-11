
      MODULE Mod_THdriver
      USE Inc_Solver
      USE Inc_Constant
      USE Inc_Flag
      USE Inc_Geometry
      USE Inc_RP
      USE Inc_Control
      USE Inc_TH
      USE Inc_Time
      USE Inc_Option
      USE Inc_MATFB
!     ! use mod_heatfunc
      use inc_xs_file

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE
      INTEGER, DIMENSION(:), ALLOCATABLE :: nthx, nthy
      REAL(8):: dzcetadt,dr2odt,tw2odt,cetacr,kgap,kgapb,kgap2,kgap4,kgap4b

      CONTAINS

      SUBROUTINE setth
      IMPLICIT NONE
      INTEGER:: i
      REAL(8) :: rhof, rhoc


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [setth] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      nzth=Nz
#ifdef siarhei_fr
      write(*,*) 'Here to initialize TH Fast Reactor (init cond)' ! @$^ Siarhei_FR
      if (fluid_type.gt.1) if_th = 3
#else
      fluid_type = 0
#endif
      
      DO i=1, Nx_1N ! T/H Mesh Same As Neutronic Mesh
         nthx(i) = Mesh_x(i)
      END DO
      DO i=1,Ny_1N ! T/H Mesh Same As Neutronic Mesh
         nthy(i) = Mesh_y(i)
      END DO
      nzth=Nz
      ! fuel conductivity correlation in w/m-C, t in K
      akfuel(0)=1.05d0   !constant term
      akfuel(1)=0d0      !1-st order
      akfuel(2)=0d0      !2-nd order
      akfuel(3)=0d0      !3-rd order
      akfuel(4)=2150d0   !-1-th order
      akfuel(5)=73.15d0  !reference in the -1 th term
      akclad(0)=7.51d0   !constant term
      akclad(1)=2.09e-2  !1-st order
      akclad(2)=-1.45e-5 !2-nd order
      akclad(3)=7.67e-9  !3-rd order
      ! volumetric heat capacity of uo2 in J/m^3-C, t in K
      rhof=10282.D0      !density of fuel kg/m^3
      arcpfuel(0)=162.3d0*rhof
      arcpfuel(1)=0.3038d0*rhof
      arcpfuel(2)=-2.391e-4*rhof
      arcpfuel(3)=6.404e-8*rhof
      rhoc=6600.D0              !density of Zr
      arcpclad(0)=252.54d0*rhoc
      arcpclad(1)=0.11474d0*rhoc
      arcpclad(2)=0d0
      arcpclad(3)=0d0
      
      END SUBROUTINE setth

!!!#ifdef siarhei_delete 
      subroutine initth
      use Mod_THfunc, only: fdens, fenthal
      use mod_manage, only: alloc_th
      use Inc_Geometry
      use Inc_FA
      use Mod_GetSome, only: Get_POW, Get_LinPOW
!     ! use Mod_H2O, only: H2O_T2D, H2O_T2H
      use Inc_XS_File, only: if_th
      implicit none
      integer, allocatable, save :: ineut(:), jneut(:)
      integer :: ia, ja, k, lc, npthy, ji, jb, iabeg, iaend
      integer :: npthx, ii, ib, j, jn, i, in, l
      integer :: nchanp1, kth, ir
      integer :: xx, yy
      real(8) :: chanvf, hin
      real(8) :: rhoin, rhouin, rhohuin, uin
      real(8) :: ptot
      real(8) :: t_inlet, g_inlet

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [initth] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

#ifdef tuan_tr_test
      if (flag_THFB) then
         call initth_hex
         return
      endif
#endif

      Ptot = 0
      allocate(ltochan(0:Nxy))
      allocate(lchantol(1:Nxy))
      allocate(lchanptr(0:Nxy))
      ltochan=0
      lchantol=0
      lchanptr=0

      ! determine neutronic node structure
      if (.not.allocated(ineut)) then
         allocate(ineut(0:Nx_1N))
         allocate(jneut(0:Ny_1N))
         ineut(0)=0
         jneut(0)=0
         do ia=1,nx_1n
            ineut(ia)=ineut(ia-1)+mesh_x(ia)
         enddo
         do ja=1,ny_1n
            jneut(ja)=jneut(ja-1)+mesh_y(ja)
         enddo
         if (nzth==0) nzth=nz

         ! channel assignment
         nchan=0
         lc=0
         do ja=1,Ny_1N
            npthy=mesh_y(ja)/nthy(ja)
            do ji=1,nthy(ja)
               jb=jneut(ja-1)+(ji-1)*npthy
               iabeg=ix_startfa_y_1n(ja)
               iaend=ix_endfa_y_1n(ja)
               do ia=iabeg,iaend
                  npthx=Mesh_x(ia)/nthx(ia)
                  do ii=1,nthx(ia)
                     nchan=nchan+1
                     ib=ineut(ia-1)+(ii-1)*npthx
                     j=jb
                     do jn=1,npthy
                        j=j+1
                        i=ib
                        do in=1,npthx
                           i=i+1
                           l=nodel(i,j)
                           ltochan(l)=nchan
                           lc=lc+1
                           lchantol(lc)=l
                        enddo
                        lchanptr(nchan)=lc
                     enddo
                  enddo
               enddo
            enddo
         enddo
         lchanptr(0)=0

         call alloc_th
         ! channel to node correspondence for the radial reflector
         nchanp1=nchan+1

         ! initialize constants
         yy=1
         do ja=1,ny_1n
            if (nthy(ja)>yy) yy=nthy(ja)
         enddo
         xx=1
         do ja=1,nx_1n
            if (nthx(ja)>xx) xx=nthy(ja)
         enddo
         chanvf=1d0/real(xx*yy,8)         !channel volume fraction per assembly

         if (opt_core==4.and.nxy==1) chanvf=0.25d0
         chanvf_save=chanvf

         acf         = ( (h_FA*DM2)**2-PI_P2R*(N_Pin*R_Clad**2+N_GT*R_GT_OClad**2) )*chanvf  ! coolant flow area
         afp         = N_Pin * PI_P2R * R_Fuel**2 * chanvf                                        ! total fuel pellet area
         xi          = 2*pi_P2R*(N_Pin*R_Clad+N_GT*R_GT_OClad)*chanvf                             ! wetted perimeter
         zeta        = N_Pin*2*pi_P2R*R_Clad*chanvf                                               ! heated perimeter
         FA_MassFlow = FA_MassFlow * chanvf
         zetap       = zeta/acf                                                                   ! heated perimeter density
         deq         = 4*acf/xi                                                                   ! equivalent diameter
         fracdf      = 1 - Frac_Gamma                                                             ! fuel heat deposit fraction

         ! radial nodalization (in fuel pin)
         nrp1=nr+1
         nrp2=nr+2
         nrp3=nr+3
         nrp4=nr+4
         nrp5=nr+5
         delr=r_fuel/real(nr,8)
         delr2=delr*delr
         delrw=0.5d0*h_clad
         delrw2=delrw*delrw
         tworm=h_clad/(r_gap+0.5d0*h_clad)
         do i=1,nrp1
            r(i)=delr*(i-1)
         enddo
         r(nrp2)=r_gap
         r(nrp3)=r_gap+delrw
         r(nrp4)=r_clad
         kgap=h_gap*delr
         kgap2=h_gap*h_clad*r_fuel/r_gap
         kgap4=h_gap*h_clad*(4d0-h_clad/r_gap)*r_fuel/r_gap

         ! assign inlet condition to all nodes
         tdopin=sqrt(TF_In+CKELVIN)

         if (if_th==4) then
            din= 1d0 !H2O_T2D(TM_In + DegToK)*1d3
         else
            din=fdens(TM_In)
!            dummy_filler = 1 ! @$^ siarhei_plot 
         endif

         do l=1,nchan
            do k=1,nzth
               do ir=1,nrp1
                  tfuel(ir,k,l)=TF_In
                  tdopl(k,l)=tdopin
               enddo
               do ir=nrp2,nrp4
                  tfuel(ir,k,l)=TF_In
               enddo
               tcool(k,l)=TM_In
               dcool(k,l)=din
               if (flag_th_chanwise) then
                  t_inlet=chanwise_t_inlet(i_chan_1n(i_4nto1n(i_fa(l))))
                  tcool(k,l)=t_inlet
                  if (if_th==4) then
                     dcool(k,l)= 1d0 !H2O_T2D(t_inlet+degtok)*1d3
                  else
                     dcool(k,l)=fdens(t_inlet)
!            dummy_filler = 1 ! @$^ siarhei_plot 
                  endif
               endif
            enddo
         enddo

         ! assign reflector t/h condition
         do kth=1,nzth
            dcool(kth,nchanp1)=din
            tcool(kth,nchanp1)=TM_In
            tdopl(kth,nchanp1)=tdopin
         enddo

         ! initialize volume and junction variables
         if (if_th==4) then
            hin=1d0 !H2O_T2H(TM_In+degtok)
            rhoin=1d0 !H2O_T2D(TM_In+degtok)*1d3
         else
            hin=fenthal(TM_In)
!            dummy_filler = 1 ! @$^ siarhei_plot 
            rhoin=fdens(TM_In)
!            dummy_filler = 1 ! @$^ siarhei_plot 
         endif
         rhouin=FA_MassFlow/acf
         rhohuin=rhouin*hin
         uin=rhouin/rhoin

         do l=1,nchan
            k=0
            if (flag_th_chanwise) then
               t_inlet=chanwise_t_inlet(i_chan_1n(i_4nto1n(i_fa(l))))
               g_inlet=chanwise_g_inlet(i_chan_1n(i_4nto1n(i_fa(l))))*chanvf_save
               if (if_th==4) then
                  hin= 1d0 !H2O_T2H(t_inlet+degtok)
                  rhoin=1d0 !H2O_T2D(t_inlet+degtok)
               else
                  hin=fenthal(t_inlet)
!            dummy_filler = 1 ! @$^ siarhei_plot 
                  rhoin=fdens(t_inlet)
!            dummy_filler = 1 ! @$^ siarhei_plot 
               endif
               rhouin=g_inlet/acf
               rhohuin=rhouin*hin
               uin=rhouin/rhoin
            endif
            rhou(k,l)=rhouin
            rhohu(k,l)=rhohuin
!write(*,*) 'tuan in initth',k,l, rhohu(k,l), rhohuin
            u(k,l)=uin
            ud(k,l)=uin
            do k=1,nzth
               hcool(k,l)=hin
               rhou(k,l)=rhouin
               rhohu(k,l)=rhohuin
               u(k,l)=uin
               ud(k,l)=uin
            enddo
         enddo

      endif

      end subroutine initth
!!!#endif 


!!!#ifdef siarhei_delete 
      SUBROUTINE trtf
      USE Inc_FA, ONLY: h_Clad
#ifdef js_r2mpi
      use inc_parallel, only: comm
#endif
#ifdef tuan_tr_test
      use inc_tpen, only: powlin, plevel0, imap
      use inc_3d, only: power!, Avg_power
#endif
      IMPLICIT NONE
      INTEGER :: Ixy, Ixy_FA, Iz,l
      REAL(8) :: qfn,qprime
      REAL(8) :: qf

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [trtf] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

!      tdoplmax=0
!      dr2odt=delr2/dT
!      tw2odt = h_Clad**2 / dT
!
!      do l = 1,nxy
!         ixy = imap(l)
!         Ixy_FA = ixytoifa(ixy) 
!         if (I_FARF_1N(ixy) .NE. 2) cycle
!         do Iz = IzFuelBot, IzFuelTop

      do Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)
         do Iz = IzFuelBot, IzFuelTop
#ifdef tuan_tr_test
            if (opt_mode ==3) then
                qprime = power(ixy, iz)*powlin*plevel0 
            else
                qprime = LinPOW(Ixy, Iz)
            endif
            

#else
            qprime = LinPOW(Ixy, Iz)
#endif            
            qfn    = fracdf * qprime / afp
            qf     = cetafb * qvol(Iz, Ixy_FA) + cetaf * qfn
            CALL tfcaltr( Iz, Ixy_FA, tcool(Iz, Ixy_FA), htcoef(Iz, Ixy_FA), qf, .FALSE. )
!            dummy_filler = 1 ! @$^ siarhei_plot 
         enddo
      enddo

      flagth=.false.
      if (tdoplmax<1d-4) flagth=.true.

      RETURN
      END SUBROUTINE trtf
!!!#endif 


#ifdef tuan_tr_test
      SUBROUTINE trtf_hex
      USE Inc_FA, ONLY: h_Clad
      use inc_tpen, only: powlin, plevel0, imap
      use inc_3d, only: power!, Avg_power

      IMPLICIT NONE
      INTEGER :: Ixy, Ixy_FA, Iz,l
      REAL(8) :: qfn,qprime
      REAL(8) :: qf

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [trtf] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif
      tdoplmax=0
      dr2odt=delr2/dT
      tw2odt = h_Clad**2 / dT

      do l = 1,nxy
         ixy = imap(l)
         Ixy_FA = ixytoifa(ixy) 
         if (I_FARF_1N(ixy) .NE. 2) cycle
         do Iz = IzFuelBot, IzFuelTop
            qprime = power(ixy, iz)*powlin*plevel0 
            qfn    = fracdf * qprime / afp
            qf     = cetafb * qvol(Iz, Ixy_FA) + cetaf * qfn
            CALL tfcaltr_hex( Iz, Ixy_FA, tcool(Iz, Ixy_FA), htcoef(Iz, Ixy_FA), qf, .FALSE. )
         enddo
      enddo
      flagth=.false.
      if (tdoplmax<1d-4) flagth=.true.
      
      RETURN
      END SUBROUTINE trtf_hex
#endif 


!!!#ifdef siarhei_delete 
      SUBROUTINE trth
      USE Inc_FluxVar
      USE Mod_THfunc !, ONLY: fhtcoef, ftemp
      USE Inc_FA , ONLY: h_Clad
!     ! USE MOD_H2O,  ONLY: H2O_H2T
      USE Inc_XS_File, ONLY: if_th
#ifdef js_r2mpi
      use inc_parallel, only: comm
#endif
      IMPLICIT NONE
      INTEGER :: Ixy, Ixy_FA, Iz
      REAL(8) :: rhouin,rhohuin,rhouinn,rhohuinn,qc, tb, qf
      REAL(8) :: hout, tout,qfn,qeffn,qprime,dz

#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [trth] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      tfmax   = 0
      dr2odt  = delr2/dT
      tw2odt  = h_Clad**2 / dT
      toutavg = 0
      
      do Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)
         rhouin   = rhou(0, Ixy_FA)
         rhohuin  = rhohu(0, Ixy_FA)
         rhouinn  = rhouin
         rhohuinn = rhohuin
         do Iz = 1, Nz
            dz=meshsize_z(Iz)*1d-2
            qprime = LinPOW(Ixy, Iz)
            qfn=fracdf*qprime/afp
            qf=cetafb*qvol(Iz, Ixy_FA)+cetaf*qfn
            qvol(Iz, Ixy_FA)=qfn
            htcoef(Iz, Ixy_FA)=fhtcoef(tcool(Iz, Ixy_FA),deq,rhou(Iz, Ixy_FA))
            tb=tcool(Iz, Ixy_FA)
            if (iz<izfuelbot.or.iz>izfueltop) then
               tfuel(:,iz,ixy_fa)=tcool(iz,ixy_fa)
            endif
            CALL tfcaltr(Iz,Ixy_FA,tb,htcoef(Iz, Ixy_FA),qf,.TRUE.)
            tfmax=MAX(tfmax,tfuel(1,Iz ,Ixy_FA))
            qflux(iz,ixy_fa)=htcoef(Iz, Ixy_FA)*(tfuel(nrp4,Iz ,Ixy_FA)-tb)
            dzcetadt=dz/(cetac*dT)
            qc    = Frac_Gamma * qprime / acf            
            qeffn = qflux(iz,ixy_fa) * zetap + qc
            qc=qeffn+cetacr*qeff(Iz, Ixy_FA)            
            qeff(Iz, Ixy_FA)=qeffn
            CALL tccaltr(Iz,Ixy_FA,rhouin,rhohuin,rhouinn,rhohuinn,qc)
         enddo
         hout = rhohu(IzFuelTop, Ixy_FA) / rhou(IzFuelTop, Ixy_FA)
         if (if_th==4) then
            tout = 1d0 ! H2O_H2T(hout)-degtok
         else
            tout = ftemp(hout)
         endif  
         if (idomz==ndomz) toutavg=toutavg+tout
      enddo
      
      toutavg = toutavg / float(nchan)

      RETURN
      END SUBROUTINE trth
!!!#endif 


!!!#ifdef siarhei_delete 
      SUBROUTINE tfcaltr(k,l,tb,hwall,qf,tfupd)
      USE Mod_THfunc !, ONLY: fkf, fkc, ftfavg,frhocpf, frhocpc
      USE Inc_FA, ONLY: R_Clad, h_Clad
      use inc_inp, only: gaph_con
      USE Inc_3D, ONLY: BU
      USE Inc_FA, ONLY: R_FUEL, R_GAP
      USE Inc_Option, ONLY: OPT_BUh
      LOGICAL :: tfupd
      REAL(8) :: kconv,kconv1,kconv4,kconv4b
      REAL(8) :: kfi,kmr,kml,kmrb,kmlb
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: kf,kfb,kfm,kfmb
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: x,rhocp,ad,al,au,b
      REAL(8) :: tb,hwall,qf
      REAL(8) :: tdopld,tdoplold,aldi,alphab,alpha,ri,qfd2,fkf_th1d
      INTEGER :: i,l,k,m,im1,ip1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [tfcaltr] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.allocated(kf)) then
         allocate(kf(nrp5))
         allocate(kfb(nrp5))
         allocate(kfm(nrp5))
         allocate(kfmb(nrp5))
         allocate(x(nrp5))
         allocate(rhocp(nrp5))
         allocate(ad(nrp5))
         allocate(al(nrp5))
         allocate(au(nrp5))
         allocate(b(nrp5))
      endif

      ! update coefficients
      tdopld=tdopl(k,l)
      do i=1,nrp4
         x(i)=tfuel(i,k,l)+CKELVIN
      enddo
      do i=1,nrp1
         kfi=fkf(x(i), l, k, fkf_th1d)
         kf(i)=cetaf*kfi
         kfb(i)=cetafb*kfi
         rhocp(i)=frhocpf(x(i))*dr2odt
      enddo
      ! coefficient at the middle points
      do i=1,nr
         kfm(i)=0.5d0*(kf(i)+kf(i+1))
         kfmb(i)=0.5d0*(kfb(i)+kfb(i+1))
      enddo
      do i=nrp2,nrp4
         kfi=fkc(i, x(i), l, k, qf)
         kf(i)=cetaf*kfi
         kfb(i)=cetafb*kfi
         rhocp(i)=frhocpc(x(i))*tw2odt
      enddo
      m=nrp3
      kmr=0.5d0*(kf(m)+kf(m+1))
      kml=0.5d0*(kf(m)+kf(m-1))
      kmrb=0.5d0*(kfb(m)+kfb(m+1))
      kmlb=0.5d0*(kfb(m)+kfb(m-1))
      ! setup linear system
      do i=1,nrp4
         x(i)=tfuel(i,k,l)
      enddo
      qfd2=qf*delr2
      m=1
      ad(m)=rhocp(m)+4*kf(m)
      au(m)=-4*kf(m)
      b(m)=qfd2+x(m)*rhocp(m)-4*kfb(m)*x(m)+4*kfb(m)*x(m+1)
      i=m
      do m=2,nr
         ri=1/float(i)
         ad(m)=rhocp(m)+(kfm(i)+kfm(m)+0.5d0*(kfm(m)-kfm(i))*ri)
         al(m)=-kfm(i)*(1d0-0.5d0*ri)
         au(m)=-kfm(m)*(1d0+0.5d0*ri)
         b(m)=qfd2+x(m)*rhocp(m)                                        &
              -(kfmb(i)+kfmb(m)+0.5d0*(kfmb(m)-kfmb(i))*ri)*x(m)           &
              +kfmb(i)*(1d0-0.5d0*ri)*x(m-1)                                 &
              +kfmb(m)*(1d0+0.5d0*ri)*x(m+1)
         i=m
      enddo

      if (OPT_BUh == 1) then
         GK_BU = BU(I_FA(l), k)
         if (GK_BU <= GK_fit_BU_cut) then
            h_gap = GK_a1 * GK_BU**6 + GK_b1 * GK_BU**5 + GK_c1 * GK_BU**4 &
               + GK_d1 * GK_BU**3 + GK_e1 * GK_BU**2 + GK_f1 * GK_BU + GK_g1
         else
            h_gap = GK_a2 * GK_BU**6 + GK_b2 * GK_BU**5 + GK_c2 * GK_BU**4 &
               + GK_d2 * GK_BU**3 + GK_e2 * GK_BU**2 + GK_f2 * GK_BU + GK_g2
         endif
         h_gap = h_gap*gaph_con
         kgap  = h_Gap*delr*cetaf
         kgapb = h_Gap*delr*cetafb
         kgap2 = h_Gap*h_Clad*R_Fuel/R_Gap*cetaf
         kgap4 = h_Gap*h_Clad*(4d0-h_Clad/R_Gap)*R_Fuel/R_Gap*cetaf
         kgap4b= h_Gap*h_Clad*(4d0-h_Clad/R_Gap)*R_Fuel/R_Gap*cetafb
      endif

      m=nrp1
      alpha=kgap*(1d0-kf(m-1)/kf(m))
      alphab=alpha/cetaf*cetafb
      ad(m)=rhocp(m)+2*(kf(m)+kgap*(1d0+0.5d0/nr))+alpha
      al(m)=-2*kf(m)
      au(m)=-2*kgap*(1+0.5/nr)-alpha
      b(m)=qfd2+x(m)*rhocp(m)-(2*(kfb(m)+kgapb*(1+0.5/nr))+alphab)*x(m) &
           +2*kfb(m)*x(m-1)+(2*kgapb*(1+0.5/nr)+alphab)*x(m+1)
      m=nrp2
      alpha=2*kgap2*(kf(m+1)/kf(m)-1)
      alphab=alpha/cetaf*cetafb
      ad(m)=rhocp(m)+8*kf(m)+kgap4-alpha
      al(m)=-kgap4+alpha
      au(m)=-8*kf(m)
      b(m)=x(m)*rhocp(m)-(8*kfb(m)+kgap4b)*x(m)+kgap4b*x(m-1)           &
           +8*kfb(m)*x(m+1)+x(m)*alphab-x(m-1)*alphab
      m=nrp3
      ad(m)=rhocp(m)+4*(kmr+kml)+tworm*(kmr-kml)
      al(m)=-kml*(4-tworm)
      au(m)=-kmr*(4+tworm)
      b(m)=x(m)*rhocp(m)-(4*(kmrb+kmlb)+tworm*(kmrb-kmlb))*x(m)         &
           +kmlb*(4-tworm)*x(m-1)+kmrb*(4+tworm)*x(m+1)
      m=nrp4
      kconv1 = cetaf * hwall * h_Clad
      alpha=2*kconv1*(1-kf(m-1)/kf(m))
      alphab=alpha/cetaf*cetafb
      kconv  = hwall * h_Clad * ( 4 + h_Clad/R_Clad )
      kconv4 = cetaf * kconv
      kconv4b=cetafb*kconv
      ad(m)=rhocp(m)+8*kf(m)+kconv4+alpha
      al(m)=-8*kf(m)
      b(m)=x(m)*rhocp(m)-(8*kfb(m)+kconv4b+alphab)*x(m)+8*kfb(m)*x(m-1) &
           +(kconv+alpha/cetaf)*tb
      ! solve the tridiagonal system by gauss elimination
      im1=1
      do i=2,nrp4
         aldi=al(i)/ad(im1)
         ad(i)=ad(i)-aldi*au(im1)
         b(i)=b(i)-aldi*b(im1)
         im1=i
      enddo
      i=nrp4
      ip1=nrp4
      x(i)=b(i)/ad(i)
      do i=nrp3,1,-1
         x(i)=(b(i)-au(i)*x(ip1))/ad(i)
         ip1=i
      enddo
      if (tfupd) then
         do i=1,nrp4
            tfuel(i,k,l)=x(i)
         enddo
      endif
      tdoplold=tdopl(k,l)
      tfuel(nrp5,k,l)=ftfavg(x,nr,delr2)
      
      if ( OPT_TFcal == 1 ) then
         ! BE1; Volume Avg
         tdopl(k,l) = tfuel(nrp5,k,l)   
      elseif ( OPT_TFcal == 2 ) then
         ! BE2; BE1*0.92 + Surf*0.08
         tdopl(k,l) = wfcl*tfuel(nrp5,k,l) + wfsurf*x(nrp1)
      elseif ( OPT_TFcal == 3 ) then
         ! NEA; Center*0.3 + Surf*0.7
         tdopl(k,l) = wfcl*x(1) + wfsurf*x(nrp1)  
      elseif ( (OPT_TFcal >= 4) .and. (OPT_TFcal <= 8) ) then
         ! BE2; BE1*0.92 + Surf*0.08
         tdopl(k,l) = wfcl*tfuel(nrp5,k,l) + wfsurf*x(nrp1)
      else
          STOP "CHECK OPT_TFcal !!"
      endif
      tdopl(k,l) = DSQRT( tdopl(k,l) + CKELVIN )
      tdoplmax=max(tdoplmax,abs(1d0-tdoplold/tdopl(k,l)))

      RETURN
      END SUBROUTINE tfcaltr
      
#ifdef tuan_tr_test 
      SUBROUTINE tfcaltr_hex(k,l,tb,hwall,qf,tfupd)
      USE Mod_THfunc !, ONLY: fkf, fkc, ftfavg,frhocpf, frhocpc
      USE Inc_FA, ONLY: R_Clad, h_Clad
      use inc_inp, only: gaph_con
      USE Inc_3D, ONLY: BU
      USE Inc_FA, ONLY: R_FUEL, R_GAP
      USE Inc_Option, ONLY: OPT_BUh
      LOGICAL :: tfupd
      REAL(8) :: kconv,kconv1,kconv4,kconv4b
      REAL(8) :: kfi,kmr,kml,kmrb,kmlb
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: kf,kfb,kfm,kfmb
      REAL(8),ALLOCATABLE,SAVE,DIMENSION(:) :: x,rhocp,ad,al,au,b
      REAL(8) :: tb,hwall,qf
      REAL(8) :: tdopld,tdoplold,aldi,alphab,alpha,ri,qfd2,fkf_th1d
      INTEGER :: i,l,k,m,im1,ip1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [tfcaltr] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.allocated(kf)) then
         allocate(kf(nrp5))
         allocate(kfb(nrp5))
         allocate(kfm(nrp5))
         allocate(kfmb(nrp5))
         allocate(x(nrp5))
         allocate(rhocp(nrp5))
         allocate(ad(nrp5))
         allocate(al(nrp5))
         allocate(au(nrp5))
         allocate(b(nrp5))
      endif

      ! update coefficients
      tdopld=tdopl(k,l)
      do i=1,nrp4
         x(i)=tfuel(i,k,l)+CKELVIN
      enddo
      do i=1,nrp1
         kfi=fkf(x(i), l, k, fkf_th1d)
         kf(i)=cetaf*kfi
         kfb(i)=cetafb*kfi
         rhocp(i)=frhocpf(x(i))*dr2odt
      enddo
      ! coefficient at the middle points
      do i=1,nr
         kfm(i)=0.5d0*(kf(i)+kf(i+1))
         kfmb(i)=0.5d0*(kfb(i)+kfb(i+1))
      enddo
      do i=nrp2,nrp4
         kfi=fkc(i, x(i), l, k, qf)
         kf(i)=cetaf*kfi
         kfb(i)=cetafb*kfi
         rhocp(i)=frhocpc(x(i))*tw2odt
      enddo
      m=nrp3
      kmr=0.5d0*(kf(m)+kf(m+1))
      kml=0.5d0*(kf(m)+kf(m-1))
      kmrb=0.5d0*(kfb(m)+kfb(m+1))
      kmlb=0.5d0*(kfb(m)+kfb(m-1))
      ! setup linear system
      do i=1,nrp4
         x(i)=tfuel(i,k,l)
      enddo
      qfd2=qf*delr2
      m=1
      ad(m)=rhocp(m)+4*kf(m)
      au(m)=-4*kf(m)
      b(m)=qfd2+x(m)*rhocp(m)-4*kfb(m)*x(m)+4*kfb(m)*x(m+1)
      i=m
      do m=2,nr
         ri=1/float(i)
         ad(m)=rhocp(m)+(kfm(i)+kfm(m)+0.5d0*(kfm(m)-kfm(i))*ri)
         al(m)=-kfm(i)*(1d0-0.5d0*ri)
         au(m)=-kfm(m)*(1d0+0.5d0*ri)
         b(m)=qfd2+x(m)*rhocp(m)                                        &
              -(kfmb(i)+kfmb(m)+0.5d0*(kfmb(m)-kfmb(i))*ri)*x(m)           &
              +kfmb(i)*(1d0-0.5d0*ri)*x(m-1)                                 &
              +kfmb(m)*(1d0+0.5d0*ri)*x(m+1)
         i=m
      enddo

      if (OPT_BUh == 1) then
         GK_BU = BU(I_FA(l), k)
         if (GK_BU <= GK_fit_BU_cut) then
            h_gap = GK_a1 * GK_BU**6 + GK_b1 * GK_BU**5 + GK_c1 * GK_BU**4 &
               + GK_d1 * GK_BU**3 + GK_e1 * GK_BU**2 + GK_f1 * GK_BU + GK_g1
         else
            h_gap = GK_a2 * GK_BU**6 + GK_b2 * GK_BU**5 + GK_c2 * GK_BU**4 &
               + GK_d2 * GK_BU**3 + GK_e2 * GK_BU**2 + GK_f2 * GK_BU + GK_g2
         endif
         h_gap = h_gap*gaph_con
         kgap  = h_Gap*delr*cetaf
         kgapb = h_Gap*delr*cetafb
         kgap2 = h_Gap*h_Clad*R_Fuel/R_Gap*cetaf
         kgap4 = h_Gap*h_Clad*(4d0-h_Clad/R_Gap)*R_Fuel/R_Gap*cetaf
         kgap4b= h_Gap*h_Clad*(4d0-h_Clad/R_Gap)*R_Fuel/R_Gap*cetafb
      endif
      
      m=nrp1
      alpha=kgap*(1d0-kf(m-1)/kf(m))
      alphab=alpha/cetaf*cetafb
      ad(m)=rhocp(m)+2*(kf(m)+kgap*(1d0+0.5d0/nr))+alpha
      al(m)=-2*kf(m)
      au(m)=-2*kgap*(1+0.5/nr)-alpha
      b(m)=qfd2+x(m)*rhocp(m)-(2*(kfb(m)+kgapb*(1+0.5/nr))+alphab)*x(m) &
           +2*kfb(m)*x(m-1)+(2*kgapb*(1+0.5/nr)+alphab)*x(m+1)
      m=nrp2
      alpha=2*kgap2*(kf(m+1)/kf(m)-1)
      alphab=alpha/cetaf*cetafb
      ad(m)=rhocp(m)+8*kf(m)+kgap4-alpha
      al(m)=-kgap4+alpha
      au(m)=-8*kf(m)
      b(m)=x(m)*rhocp(m)-(8*kfb(m)+kgap4b)*x(m)+kgap4b*x(m-1)           &
           +8*kfb(m)*x(m+1)+x(m)*alphab-x(m-1)*alphab
      m=nrp3
      ad(m)=rhocp(m)+4*(kmr+kml)+tworm*(kmr-kml)
      al(m)=-kml*(4-tworm)
      au(m)=-kmr*(4+tworm)
      b(m)=x(m)*rhocp(m)-(4*(kmrb+kmlb)+tworm*(kmrb-kmlb))*x(m)         &
           +kmlb*(4-tworm)*x(m-1)+kmrb*(4+tworm)*x(m+1)
      m=nrp4
      kconv1 = cetaf * hwall * h_Clad
      alpha=2*kconv1*(1-kf(m-1)/kf(m))
      alphab=alpha/cetaf*cetafb
      kconv  = hwall * h_Clad * ( 4 + h_Clad/R_Clad )
      kconv4 = cetaf * kconv
      kconv4b=cetafb*kconv
      ad(m)=rhocp(m)+8*kf(m)+kconv4+alpha
      al(m)=-8*kf(m)
      b(m)=x(m)*rhocp(m)-(8*kfb(m)+kconv4b+alphab)*x(m)+8*kfb(m)*x(m-1) &
           +(kconv+alpha/cetaf)*tb
      ! solve the tridiagonal system by gauss elimination
      im1=1
      do i=2,nrp4
         aldi=al(i)/ad(im1)
         ad(i)=ad(i)-aldi*au(im1)
         b(i)=b(i)-aldi*b(im1)
         im1=i
      enddo
      i=nrp4
      ip1=nrp4
      x(i)=b(i)/ad(i)
      do i=nrp3,1,-1
         x(i)=(b(i)-au(i)*x(ip1))/ad(i)
         ip1=i
      enddo
      if (tfupd) then
         do i=1,nrp4
            tfuel(i,k,l)=x(i)
         enddo
      endif
      tdoplold=tdopl(k,l)
      tfuel(nrp5,k,l)=ftfavg(x,nr,delr2)
      
      if ( OPT_TFcal == 1 ) then
         ! BE1; Volume Avg
         tdopl(k,l) = tfuel(nrp5,k,l)   
      elseif ( OPT_TFcal == 2 ) then
         ! BE2; BE1*0.92 + Surf*0.08
         tdopl(k,l) = wfcl*tfuel(nrp5,k,l) + wfsurf*x(nrp1)
      elseif ( OPT_TFcal == 3 ) then
         ! NEA; Center*0.3 + Surf*0.7
         tdopl(k,l) = wfcl*x(1) + wfsurf*x(nrp1)  
      elseif ( (OPT_TFcal >= 4) .and. (OPT_TFcal <= 8) ) then
         ! BE2; BE1*0.92 + Surf*0.08
         tdopl(k,l) = wfcl*tfuel(nrp5,k,l) + wfsurf*x(nrp1)
      else
          STOP "CHECK OPT_TFcal !!"
      endif
      tdopl(k,l) = DSQRT( tdopl(k,l) + CKELVIN )
      tdoplmax=max(tdoplmax,abs(1d0-tdoplold/tdopl(k,l)))

      RETURN
      END SUBROUTINE tfcaltr_hex
#endif 


!!!#ifdef siarhei_delete 
      SUBROUTINE tccaltr(k,l,rhouin,rhohuin,rhouinn,rhohuinn,qc)
      USE Mod_THfunc !, ONLY: fdensh, ftemp
!     ! USE MOD_H2O,  ONLY: H2O_H2D, H2O_H2T
      USE Inc_XS_File, ONLY: if_th
      IMPLICIT NONE
      ! transient coolant temperature calculation
      INTEGER :: k,l
      REAL(8) :: rhouin,rhohuin,rhouinn,rhohuinn,qc,hr,hl,rhor,rhol
      REAL(8) :: rhoh,toutn,rhohn,rhooutn,rhououtn,rhohuoutn,rhouout,rhohuout
      REAL(8) :: rhon,hn,tn,houtn,sqterm,tmpterm
      REAL(8) :: a,b,c,dh,dz,eps,rho,alpha,beta,biga,bigb
      REAL(8) :: delh,delh1,delh2,delh3,delta,delrhou,drhodh,drhohdh,delrhohu
      REAL(8) :: h,hinn,hout,uavg,uout,uoutn,uoutd
      REAL(8) :: gamma1
      DATA eps/0.001/

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [tccaltr] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      h=hcool(k,l)
      rho=dcool(k,l)
      rhouout=rhou(k,l)
      rhohuout=rhohu(k,l)
      uoutd=ud(k,l)
      uout=u(k,l)
      dz=meshsize_z(k)*1d-2
      ! calculate derivatives
      hr=(1d0+eps)*h
      hl=(1d0-eps)*h
      dh=hr-hl
      if (if_th==4) then
         rhor=1d0 !H2O_H2D(hr)*1d3 
         rhol=1d0 !H2O_H2D(hl)*1d3
      else
         rhor=fdensh(hr)
         rhol=fdensh(hl)
      endif  
      drhodh=(rhor-rhol)/dh
      drhohdh=(rhor*hr-rhol*hl)/dh
      ! calculate coefficients
      delrhou=rhouout-rhouin
      delrhohu=rhohuout-rhohuin
      hinn=rhohuinn/rhouinn
      alpha=-dzcetadt*drhodh
      beta=rhouinn-cetacr*delrhou
      tmpterm=2d0*h-hinn
      gamma1=alpha*tmpterm+2d0*beta
      delta=beta*tmpterm
      a=2d0*alpha
      b=gamma1+dzcetadt*drhohdh
      c=cetacr*delrhohu+delta-rhohuinn-dz*qc
      sqterm=SQRT(b*b-4d0*a*c)
      delh1=0.5d0*(-b+sqterm)/a
      delh2=0.5d0*(-b-sqterm)/a
      hout=rhohuout/rhouout
      biga=dzcetadt*drhohdh+alpha*hout
      bigb=dz*qc-rhouinn*hout+rhohuinn+cetacr*(delrhou*hout-delrhohu)
      delh3=bigb/biga
      if (abs(delh1-delh3)<abs(delh2-delh3)) then
         delh=delh1
      else
         delh=delh2
      endif
      hn=h+delh
      if (if_th==4) then
         rhon=1d0 !H2O_H2D(hn)*1d3
         tn=1d0 !H2O_H2T(hn)-degtok
      else
         rhon=fdensh(hn)
         tn=ftemp(hn)
      endif  
      rhoh=rho*h
      rhohn=rhon*hn
      rhououtn=rhouinn-dzcetadt*(rhon-rho)-cetacr*delrhou
      rhohuoutn=rhohuinn-dzcetadt*(rhohn-rhoh)-cetacr*delrhohu+dz*qc
      houtn=rhohuoutn/rhououtn
      if (if_th==4) then
         toutn=1d0 !H2O_H2T(houtn)-degtok
         rhooutn=1d0 !H2O_H2D(houtn)*1d3
      else
         toutn=ftemp(houtn)
         rhooutn=fdensh(houtn)
      endif  
      uoutn=rhououtn/rhooutn
      uavg=cetac*(uoutn+cetacr*uout)
      if (((uoutd-uout)*(uoutn-uout)<-1d-8*uavg*uavg).and.(cetac<1d0)) then
         uoutn=uavg
         rhououtn=rhooutn*uoutn
         rhohuoutn=rhououtn*houtn
      endif
      ! update variables
      dcool(k,l)=rhon
      tcool(k,l)=tn
      hcool(k,l)=hn
      ud(k,l)=u(k,l)
      u(k,l)=uoutn
      rhouin=rhou(k,l)
      rhohuin=rhohu(k,l)
      rhou(k,l)=rhououtn
      rhohu(k,l)=rhohuoutn
      rhouinn=rhououtn
      rhohuinn=rhohuoutn

      RETURN
      END SUBROUTINE tccaltr
!!!#endif 

#ifdef tuan_tr_test 

      SUBROUTINE tccaltr_hex(k,l,rhouin,rhohuin,rhouinn,rhohuinn,qc)
      USE Mod_THfunc !, ONLY: fdensh, ftemp
!     ! USE MOD_H2O,  ONLY: H2O_H2D, H2O_H2T
      USE Inc_XS_File, ONLY: if_th
      IMPLICIT NONE
      ! transient coolant temperature calculation
      INTEGER :: k,l
      REAL(8) :: rhouin,rhohuin,rhouinn,rhohuinn,qc,hr,hl,rhor,rhol
      REAL(8) :: rhoh,toutn,rhohn,rhooutn,rhououtn,rhohuoutn,rhouout,rhohuout
      REAL(8) :: rhon,hn,tn,houtn,sqterm,tmpterm
      REAL(8) :: a,b,c,dh,dz,eps,rho,alpha,beta,biga,bigb
      REAL(8) :: delh,delh1,delh2,delh3,delta,delrhou,drhodh,drhohdh,delrhohu
      REAL(8) :: h,hinn,hout,uavg,uout,uoutn,uoutd
      REAL(8) :: gamma1
      DATA eps/0.001/

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [tccaltr] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      h=hcool(k,l)
      rho=dcool(k,l)
      rhouout=rhou(k,l)
      rhohuout=rhohu(k,l)
      uoutd=ud(k,l)
      uout=u(k,l)
      dz=meshsize_z(k)*1d-2
      ! calculate derivatives
      hr=(1d0+eps)*h
      hl=(1d0-eps)*h
      dh=hr-hl
      if (if_th==4) then
         rhor=1d0 !H2O_H2D(hr)*1d3 
         rhol=1d0 !H2O_H2D(hl)*1d3
      else
         rhor=fdensh(hr)
         rhol=fdensh(hl)
      endif  
      drhodh=(rhor-rhol)/dh
      drhohdh=(rhor*hr-rhol*hl)/dh
      ! calculate coefficients
      delrhou=rhouout-rhouin
      delrhohu=rhohuout-rhohuin
      hinn=rhohuinn/rhouinn
      alpha=-dzcetadt*drhodh
      beta=rhouinn-cetacr*delrhou
      tmpterm=2d0*h-hinn
      gamma1=alpha*tmpterm+2d0*beta
      delta=beta*tmpterm
      a=2d0*alpha
      b=gamma1+dzcetadt*drhohdh
      c=cetacr*delrhohu+delta-rhohuinn-dz*qc
      sqterm=SQRT(b*b-4d0*a*c)
      delh1=0.5d0*(-b+sqterm)/a
      delh2=0.5d0*(-b-sqterm)/a
      hout=rhohuout/rhouout
      biga=dzcetadt*drhohdh+alpha*hout
      bigb=dz*qc-rhouinn*hout+rhohuinn+cetacr*(delrhou*hout-delrhohu)
      delh3=bigb/biga
      if (abs(delh1-delh3)<abs(delh2-delh3)) then
         delh=delh1
      else
         delh=delh2
      endif
      hn=h+delh
      if (if_th==4) then
         rhon=1d0 !H2O_H2D(hn)*1d3
         tn=1d0 !H2O_H2T(hn)-degtok
      else
         rhon=fdensh(hn)
         tn=ftemp(hn)
      endif  
      rhoh=rho*h
      rhohn=rhon*hn
      rhououtn=rhouinn-dzcetadt*(rhon-rho)-cetacr*delrhou
      rhohuoutn=rhohuinn-dzcetadt*(rhohn-rhoh)-cetacr*delrhohu+dz*qc
      houtn=rhohuoutn/rhououtn
      if (if_th==4) then
         toutn=1d0 !H2O_H2T(houtn)-degtok
         rhooutn=1d0 !H2O_H2D(houtn)*1d3
      else
         toutn=ftemp(houtn)
         rhooutn=fdensh(houtn)
      endif  
      uoutn=rhououtn/rhooutn
      uavg=cetac*(uoutn+cetacr*uout)
      if (((uoutd-uout)*(uoutn-uout)<-1d-8*uavg*uavg).and.(cetac<1d0)) then
         uoutn=uavg
         rhououtn=rhooutn*uoutn
         rhohuoutn=rhououtn*houtn
      endif
      ! update variables
      dcool(k,l)=rhon
      tcool(k,l)=tn
      hcool(k,l)=hn
      ud(k,l)=u(k,l)
      u(k,l)=uoutn
      rhouin=rhou(k,l)
      rhohuin=rhohu(k,l)
      rhou(k,l)=rhououtn
      rhohu(k,l)=rhohuoutn
      rhouinn=rhououtn
      rhohuinn=rhohuoutn

      RETURN
      END SUBROUTINE tccaltr_hex
#endif 


!!#ifdef siarhei_delete 
      subroutine ssth
      use Mod_THfunc
!     ! USE MOD_H2O,  ONLY: H2O_H2D, H2O_H2T
      USE Inc_XS_File, ONLY: if_th
#ifdef js_r2mpi
      use inc_parallel, only: comm
#endif
      implicit none
      integer :: ixy,ixy_fa,iz
      real(8) :: qc,qprime,qeffnew,rhohuin,rhohuout
      real(8) :: tfmaxt,toutavgt,fnchan,qf,rhouin
      real(8) :: hout,tout
      real(8) :: rtmp, rvol

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [ssth] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

#ifdef tuan_tr_test
!      if (if_hexgeometry) then
         call ssth_hex
         return
!      endif
#endif

#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif
      tfmax    = 0
      toutavg  = 0
      tdoplmax = 0
      ! coolant temperature calculation
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         rhouin=rhou(0,ixy_fa)
         rhohuin=rhohu(0,ixy_fa)
         do iz=1,nz
            qprime=LinPOW(ixy,iz) ! W/m
            qf=fracdf*qprime/afp
            qvol(iz,ixy_fa)=qf
            qc=Frac_Gamma*qprime/acf
            qflux(iz,ixy_fa)=qf*afp/zeta
            qeffnew=qflux(iz,ixy_fa)*zetap+qc
            qeff(iz,ixy_fa)=qeffnew
            rhohuout=rhohuin+qeffnew*(meshsize_z(iz)*1d-2)
            hcool(iz,ixy_fa)=0.5d0*(rhohuout+rhohuin)/rhouin
            if (if_th==4) then
               tcool(iz,ixy_fa)=1 !H2O_H2T(hcool(iz,ixy_fa))-degtok
               dcool(iz,ixy_fa)=1 !H2O_H2D(hcool(iz,ixy_fa))*1d3
            else
               tcool(iz,ixy_fa)=ftemp(hcool(iz,ixy_fa))
               dcool(iz,ixy_fa)=fdensh(hcool(iz,ixy_fa))
            endif
            rhohuin=rhohuout
            rhohu(iz,ixy_fa)=rhohuout
         enddo
         hout=rhohu(nzth,ixy_fa)/rhou(nzth,ixy_fa)
         if (if_th==4) then
            tout=1 !H2O_H2T(hout)-degtok
         else
            tout=ftemp(hout)
         endif
         if (idomz==ndomz) toutavg=toutavg+tout
      enddo

      ! coolant temperature for reflector
      do iz=1,nz
         rtmp=0d0
         rvol=0d0
         do ixy_fa=1,nxy_fa
            ixy=i_fa(ixy_fa)
            !rtmp=rtmp+tcool(iz,ixy_fa)*nodevolume(ixy,iz)
            rtmp=rtmp+hcool(iz,ixy_fa)*nodevolume(ixy,iz)
            rvol=rvol+nodevolume(ixy,iz)
         enddo
         !tcool(iz,nxy_fa+1)=rtmp/max(1d-10,rvol)
         hcool(iz,nxy_fa+1)=rtmp/max(1d-10,rvol)
         tcool(iz,nxy_fa+1)=ftemp(hcool(iz,nxy_fa+1))
         dcool(iz,nxy_fa+1)=fdensh(hcool(iz,nxy_fa+1))
      enddo

      ! fuel temperature calculation
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         do iz=izfuelbot,izfueltop
            if (opt_tfutab) then
               call get_tfutab(iz,ixy_fa,tcool(iz,ixy_fa))
!            dummy_filler = 1 ! @$^ siarhei_plot 
               tfmax=max(tfmax,tfuel(1,iz,ixy_fa))
            else
               qprime=LinPOW(ixy,iz)
               qf=fracdf*qprime/afp
               ! heat transfer coefficient
               htcoef(iz,ixy_fa)=fhtcoef(tcool(iz,ixy_fa),deq,rhou(iz,ixy_fa))
               call tfcalss(iz,ixy_fa,tcool(iz,ixy_fa),htcoef(iz,ixy_fa),qf)
               tfmax=max(tfmax,tfuel(1,iz,ixy_fa))
            endif
         enddo
      enddo

      flagth=.false.
      if (tdoplmax<EPS_DoppTF) flagth=.true.
      tfmaxt=tfmax
      toutavgt=toutavg
      fnchan=nchan
      tfmax=tfmaxt
      fnchan=fnchan/float(ndomz)
      toutavg=toutavgt/fnchan

      return
      end subroutine ssth
!!!#endif 


      subroutine tfcalss(k,l,tb,hwall,qf)
      use mod_charedit, only: print_msg
      use mod_thfunc !, only: fkf, fkc, ftfavg
      use inc_inp, only: gaph_con
      use inc_fa, only: R_Clad, h_Clad
      use mod_alloc
      use inc_3d, only: BU
      use inc_fa, only: R_FUEL, R_GAP
      use inc_option, only: OPT_BUh
      use inc_inp, only: I_ATF_TCD_FA
      use inc_option, only: OPT_TCD
      implicit none
      integer :: i, l, k, m, itr, im1, ip1
      real(8) :: kconv, kconv1
      real(8) :: kmr, kml
      real(8), allocatable,save :: kf(:)
      real(8), allocatable,save :: kfb(:)
      real(8), allocatable,save :: kfm(:)
      real(8), allocatable,save :: kfmb(:)
      real(8), allocatable,save :: x(:)
      real(8), allocatable,save :: xd(:)
      real(8), allocatable,save :: ad(:)
      real(8), allocatable,save :: al(:)
      real(8), allocatable,save :: au(:)
      real(8), allocatable,save :: b(:)
      real(8) :: tb, hwall, qf
      real(8) :: tdoplold, errtf, qfd2, ri, alpha, aldi, fkf_th1d
      real(8) :: EPS_DoppTF
      data EPS_DoppTF/0.01D0/


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [tfcalss] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (OPT_TCD(1)==9) then
         OPT_Tburn=I_ATF_TCD_FA(l)
      endif

      if (.not.allocated(kf)) then
         call alloc( kf   , nrp5 )
         call alloc( kfb  , nrp5 )
         call alloc( kfm  , nrp5 )
         call alloc( kfmb , nrp5 )
         call alloc( x    , nrp5 )
         call alloc( xd   , nrp5 )
         call alloc( ad   , nrp5 )
         call alloc( al   , nrp5 )
         call alloc( au   , nrp5 )
         call alloc( b    , nrp5 )
      endif

      errtf=CINF
      itr=0

      do while (errtf>EPS_DoppTF)
         itr=itr+1
         ! update coefficients
         do i=1,nrp4
            x(i)=tfuel(i,k,l)+CKELVIN
         enddo
         do i=1,nrp1
            kf(i)=fkf(x(i), l, k, fkf_th1d)
         enddo
                  
         ! coefficient at the middle points
         do i=1,nr
            kfm(i)=0.5d0*(kf(i)+kf(i+1))
         enddo
         do i=nrp2,nrp4
            kf(i)=fkc(i, x(i), l, k ,qf)
         enddo
      
         m=nrp3
         kmr=0.5d0*(kf(m)+kf(m+1))
         kml=0.5d0*(kf(m)+kf(m-1))

         ! setup linear system
         do i=1,nrp4
            x(i)=tfuel(i,k,l)
            xd(i)=x(i)
         enddo
         qfd2=qf*delr2
         m=1
         ad(m)=4d0*kf(m)
         au(m)=-4d0*kf(m)
         b(m)=qfd2
         i=m
         do m=2,nr
            ri=1/float(i)
            ad(m)=kfm(i)+kfm(m)+0.5d0*(kfm(m)-kfm(i))*ri
            al(m)=-kfm(i)*(1d0-0.5d0*ri)
            au(m)=-kfm(m)*(1d0+0.5d0*ri)
            b(m)=qfd2
            i=m
         enddo
               
         if (OPT_BUh == 1) then
            GK_BU = BU(I_FA(l), k)
            if (GK_BU <= GK_fit_BU_cut) then
               h_gap = GK_a1 * GK_BU**6 + GK_b1 * GK_BU**5 + GK_c1 * GK_BU**4 &
                  + GK_d1 * GK_BU**3 + GK_e1 * GK_BU**2 + GK_f1 * GK_BU + GK_g1
            else
               h_gap = GK_a2 * GK_BU**6 + GK_b2 * GK_BU**5 + GK_c2 * GK_BU**4 &
                  + GK_d2 * GK_BU**3 + GK_e2 * GK_BU**2 + GK_f2 * GK_BU + GK_g2
            endif
            h_gap = h_gap*gaph_con
            kgap  = h_Gap*delr
            kgap2 = h_Gap*h_Clad*R_Fuel/R_Gap
            kgap4 = h_Gap*h_Clad*(4-h_Clad/R_Gap)*R_Fuel/R_Gap
         endif

         m=nrp1
         alpha=kgap*(1-kf(m-1)/kf(m))
         ad(m)=2*(kf(m)+kgap*(1+0.5/nr))+alpha
         al(m)=-2*kf(m)
         au(m)=-2*kgap*(1+0.5/nr)-alpha
         b(m)=qfd2
         m=nrp2
         alpha=2*kgap2*(kf(m+1)/kf(m)-1)
         ad(m)=8*kf(m)+kgap4-alpha
         al(m)=-kgap4+alpha
         au(m)=-8*kf(m)
         b(m)=0
         m=nrp3
         ad(m)=4*(kmr+kml)+tworm*(kmr-kml)
         al(m)=-kml*(4-tworm)
         au(m)=-kmr*(4+tworm)
         b(m)=0
         m=nrp4
         kconv1 = hwall * h_Clad
         alpha=2*kconv1*(1-kf(m-1)/kf(m))
         kconv = hwall * h_Clad * ( 4+ h_Clad/R_Clad )
         ad(m)=8*kf(m)+kconv+alpha
         al(m)=-8*kf(m)
         b(m)=(kconv+alpha)*tb
         ! solve the tridiagonal system by gauss elimination
         im1=1
         do i=2,nrp4
            aldi=al(i)/ad(im1)
            ad(i)=ad(i)-aldi*au(im1)
            b(i)=b(i)-aldi*b(im1)
            im1=i
         enddo
         i=nrp4
         ip1=nrp4
         x(i)=b(i)/ad(i)
         do i=nrp3,1,-1
            x(i)=(b(i)-au(i)*x(ip1))/ad(i)
            ip1=i
         enddo
         errtf=0
         do i=1,nrp4
            tfuel(i,k,l)=x(i)
            errtf=MAX(errtf,ABS(x(i)-xd(i)))
         enddo
         if (itr>1000) then
            call print_msg(2,'Fuel temp. calculation is not converged')
            exit
         endif
      enddo
            
      tfuel(nrp5,k,l)=ftfavg(x,nr,delr2)
      tdoplold=tdopl(k,l)

      if (OPT_TFcal==1) then
         ! BE1; Volume Avg
         tdopl(k,l)=tfuel(nrp5,k,l)   
      elseif (OPT_TFcal==2) then
         ! BE2; BE1*0.92 + Surf*0.08
         tdopl(k,l)=wfcl*tfuel(nrp5,k,l)+wfsurf*x(nrp1)
      elseif (OPT_TFcal==3) then
         ! NEA; Center*0.3 + Surf*0.7
         tdopl(k,l)=wfcl*x(1)+wfsurf*x(nrp1)  
      elseif ((OPT_TFcal>=4).and.(OPT_TFcal<=8)) then
         ! BE2; BE1*0.92 + Surf*0.08
         tdopl(k,l)=wfcl*tfuel(nrp5,k,l)+wfsurf*x(nrp1)
      else
          stop "CHECK OPT_TFcal !!"
      endif
      
      tdopl(k,l)=sqrt(tdopl(k,l)+CKELVIN)
      tdoplmax=max(tdoplmax,abs(1d0-tdoplold/tdopl(k,l)))

      return
      end subroutine tfcalss

!!!#ifdef siarhei_delete 
      subroutine get_tfutab(iz,ixy_fa,tmo)
      use inc_3d, only: bu, power, avg_power
      implicit none
      integer :: iz, ixy_fa, ixy
      integer :: ibu, jbu, ipow, jpow
      real(8) :: rpf, x, x1, x2, y1, y2
      real(8) :: table_tfu_corr(5)
      real(8) :: sol, tabtfu, tfu, tmo

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_tfutab] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      ixy=i_fa(ixy_fa)

      ! get relative power fraction
      rpf=power(ixy,iz)/avg_power*PPower

      ! find BU points
      jbu=1
      do ibu=1,10
         if (BU(ixy,iz)>table_bu(ibu)) jbu=ibu
      enddo
      table_tfu_corr=0d0
      do ipow=1,5
         x=BU(ixy,iz)
         x1=table_bu(jbu)
         x2=table_bu(jbu+1)
         y1=table_tfu(ipow,jbu)
         y2=table_tfu(ipow,jbu+1)
         table_tfu_corr(ipow)=(y1*(x-x2)-y2*(x-x1))/(x1-x2)
      enddo

      ! Lagrange interpolation
      tabtfu=0d0
      do ipow=1,5
         sol=table_tfu_corr(ipow)
         do jpow=1,5
            if (ipow==jpow) cycle
            sol=sol*(rpf-table_pow(jpow))/(table_pow(ipow)-table_pow(jpow))
         enddo
         tabtfu=tabtfu+sol
      enddo

      ! get TFU
      tmo=tmo ! C
      tfu=tmo+degtok+segtfu(1)+(segtfu(2)+tabtfu)*rpf+segtfu(3)*rpf*rpf ! K

      ! assign tfuel & tdopl
      tfuel(1:14,iz,ixy_fa)=tfu
      tdopl(iz,ixy_fa)=sqrt(tfu)
      return
      end subroutine get_tfutab
!!!#endif 


!!!#ifdef siarhei_delete 
      subroutine change_massflow(tmp_tm_Avg,i_tavg_itr,flag_tavg_conv)
      use inc_flag, only: flag_tmmhist
      use inc_history, only: hist_tmm_inp
      use inc_xs_file, only: i_bu
      use mod_charedit, only: print_msg
      implicit none
      real(8),intent(in) :: tmp_tm_avg
      integer(4),intent(inout) :: i_tavg_itr
      logical(1),intent(inout) :: flag_tavg_conv
      real(8) :: temp
      real(8) :: slope
      real(8) :: dtt

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [change_massflow] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (opt_findtavg_init.and.i_bu>1) then
         core_massflow=FA_massflow/chanvf_save*Core_N_FA
         call print_msg(0,'* Skip finding mass flow rate:',core_massflow,' kg/sec')
         return
      endif
      if (flag_tmmhist.and.i_bu==1) then
         if (abs(target_tavg-hist_tmm_inp(i_bu))>0.001d0) then
            call print_msg(0,'Error in hist tmm')
         endif
      endif

      if (flag_tmmhist.and.i_bu>1) then
         target_tavg=hist_tmm_inp(i_bu)
      endif
      temp=tmp_tm_Avg-degtok
      i_tavg_itr=i_tavg_itr+1

      if (abs(temp-target_tavg)<0.1d0) then
         flag_tavg_conv=.true.
      else
         flag_tavg_conv=.false.
      endif

      ! backup for extrapolation
      itr_tavg(1)=itr_tavg(2)
      itr_tavg(2)=itr_tavg(3)
      itr_fa_flowrate(1)=itr_fa_flowrate(2)
      itr_fa_flowrate(2)=itr_fa_flowrate(3)

      itr_tavg(3)=temp               ! current one
      itr_fa_flowrate(3)=fa_massflow ! current one
    
      if (abs((itr_tavg(3)-itr_tavg(2)))<1.0d-4) then
         dtt=1.0D-4
      else
         dtt=(itr_tavg(3)-itr_tavg(2))
      endif
      slope=(itr_fa_flowrate(3)-itr_fa_flowrate(2))/dtt !(itr_tavg(3)-itr_tavg(2))
      slope=min(-1d-3,slope)
      slope=max(-100d0,slope)
      fa_massflow=slope*(target_tavg-itr_tavg(3))+itr_fa_flowrate(3)
      core_massflow=FA_massflow/chanvf_save*Core_N_FA

      call print_msg(0,'* Update core mass flow rate: ',core_massflow,' kg/sec')
      if (core_massflow<0d0) then
         call print_msg(3,'Negative core mass flow rate')
         stop
      elseif (core_massflow/inp_massflow<0.5d0) then
         call print_msg(3,'Diverge core mass flow rate')
         stop
      elseif (core_massflow/inp_massflow>1.7d0) then
         call print_msg(3,'Diverge core mass flow rate')
         stop
      endif

      call restore_massflow

      return
      end subroutine change_massflow
!!!#endif 


!!!#ifdef siarhei_delete 
      subroutine restore_massflow
      use Mod_THfunc !, only: fdens, fenthal
!     ! use MOD_H2O, only: H2O_T2D, H2O_T2H
      use Inc_XS_File, only: if_th
      implicit none
      integer(4) :: ixy,iz
      REAL(8) :: hin, rhoin, rhouin, rhohuin,uin

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [restore_massflow] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (if_th==4) then
         hin=1d0 !H2O_T2H(TM_In + DegToK)
         rhoin=1d0 !H2O_T2D(TM_In + DegToK)*1d3
      else
         hin    = fenthal(TM_In)
         rhoin  = fdens(TM_In)         
      endif  

      rhouin  = FA_MassFlow / acf
      rhohuin = rhouin * hin
      uin     = rhouin / rhoin

      do ixy = 1, nchan
         Iz = 0
         rhou (iz,ixy) = rhouin
         rhohu(iz,ixy) = rhohuin
         u    (iz,ixy) = uin
         ud   (iz,ixy) = uin
         do iz = 1, nzth
            rhou (iz,ixy) = rhouin
            rhohu(iz,ixy) = rhohuin
            u    (iz,ixy) = uin
            ud   (iz,ixy) = uin
         enddo
      enddo
      return
      end subroutine restore_massflow
!!!#endif 


!!!#ifdef siarhei_delete 
      subroutine restore_tm_in
      use Mod_THfunc !, only: fdens, fenthal
      use mod_manage, only: alloc_th
!     ! USE MOD_H2O,  ONLY: H2O_T2D
      USE Inc_XS_File, ONLY: if_th
      implicit none
      integer :: k, l, nchanp1, kth

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [restore_tm_in] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif
      ! assign inlet condition to all nodes
      if (if_th==4) then
         din=1d0 !H2O_T2D(TM_In + DegToK)*1d3
      else
         din=fdens(TM_In)
      endif  
      do l=1,nchan
         do k=1,nzth
            tcool(k,l)=TM_In
            dcool(k,l)=din
         enddo
      enddo
      nchanp1=nchan+1
      ! assign reflector t/h condition
      do kth=1,nzth
         dcool(kth,nchanp1)=din
         tcool(kth,nchanp1)=TM_In
      enddo

      return
      end subroutine restore_tm_in
!!!#endif 

!!!#ifdef siarhei_delete 
      subroutine restore_tm_in_slb(dp,dtm_in,ichan,flag_rf)
      use Mod_THfunc !, only: fdens_slb, fenthal, fdens
      use mod_manage, only: alloc_th
      use inc_transient
!     ! use MOD_H2O, only: H2O_T2D, H2O_T2H
      implicit none
      integer :: k, l
      integer :: nchanp1
      real(8) :: dtm_in, p, dp
      integer :: ichan, flag_rf

      nchanp1=nchan+1

      do k=1,nzth
         do l=1,Nxy_FA
            if (ichan == i_chan_fa_4n(l)) then
               p=dp+Core_Pressure
               tcool(k,l)=tcool(k,l)+dtm_in
!               din=fdens_slb(p,tcool(k,l),opt_ctf)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
               dcool(k,l)=din
            endif
         enddo
         if (flag_rf==1) then
            p=dp+Core_Pressure
            tcool(k,nchanp1)=tcool(k,nchanp1)+dtm_in
!            din=fdens_slb(p,tcool(k,nchanp1),opt_ctf)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
            dcool(k,nchanp1)=din
         endif
      enddo

      return
      end subroutine restore_tm_in_slb
!!!#endif 


      subroutine cal_dnbr
      use inc_geometry, only: izfuelbot, izfueltop, nz, nxy
      use inc_geometry, only: i_4nto1n
      use inc_rp, only: nxy_fa, i_fa
      use mod_thfunc !, only: fenthal
!     ! use mod_heatfunc, only: calc_qchf_w3
      use mod_alloc
!     ! use MOD_H2O, only: H2O_T2D, H2O_T2H
      use Inc_XS_File, only: if_th
      implicit none
      integer :: ixy_fa, ixy, iz, ixyp, izp
      real(8) :: de, hin, hw, hv, h, pres, go, qual
      real(8) :: t_inlet, g_inlet


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [cal_dnbr] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.allocated(qchf)) call alloc ( qchf , nxy , nz )
      if (.not.allocated(qwat)) call alloc ( qwat , nxy , nz )
      if (.not.allocated(dnbr)) call alloc ( dnbr , nxy , nz )

      do iz=izfuelbot,izfueltop
         do ixy_fa=1,nxy_fa
            ixy=i_fa(ixy_fa)
            de=deq*3.28084d0                       ! hydraulic diameter  ! m -> ft
            if (flag_th_chanwise) then
               t_inlet=chanwise_t_inlet(i_chan_1n(i_4nto1n(ixy)))
               g_inlet=chanwise_g_inlet(i_chan_1n(i_4nto1n(ixy)))*chanvf_save
            else
               t_inlet=tm_in
               g_inlet=fa_massflow
            endif
            if (if_th==4) then
               hin=1 !H2O_T2H(t_inlet+degtok)*0.0004299226d0 ! inlet enthalpy      ! J/kg -> Btu/lbm
            else
               hin=fenthal(t_inlet)*0.0004299226d0       ! inlet enthalpy      ! J/kg -> Btu/lbm
!            dummy_filler = 1 ! @$^ siarhei_plot 
            endif
            ! https://www.tlv.com/global/TI/calculator/steam-table-pressure.html#
            ! hw : specific enthalpy of saturated water (only depends on pressure, ref 155.11 bar)
            ! hv : specific enthalpy of saturated vapor (only depends on pressure, ref 155.11 bar)
            hw   = 1634290.0d0*0.0004299226d0          ! saturation enthalpy ! J/kg -> Btu/lbm
            hv   = 2592820.0d0*0.0004299226d0          ! saturation enthalpy ! J/kg -> Btu/lbm
            if (if_th==4) then
               h = 1 !H2O_T2H(tcool(iz,ixy_fa) + DegToK)*0.0004299226d0
            else
               h = fenthal(tcool(iz,ixy_fa))*0.0004299226d0
!            dummy_filler = 1 ! @$^ siarhei_plot 
            endif
            pres = core_pressure*14.5038d0             ! pressure            ! bar -> psi
            go   = g_inlet/acf*737.3376d0              ! mass flux           ! kg/sec-m2 -> lbm/hr-ft2
            ! quality = (h - hw_sat) / (hv_sat - hw_sat)
            qual = (h-hw)/(hv-hw)                      ! equilibrium quality ! [-]

           ! qchf(ixy,iz) = calc_qchf_W3(de,hin,hw,pres,go,qual)
!            dummy_filler = 1 ! @$^ siarhei_plot 
            qchf(ixy,iz) = qchf(ixy,iz)/0.31721046d0                        ! Btu/hr-ft2 -> W/m2
            qwat(ixy,iz) = fracdf*LinPOW(ixy,iz)/zeta                       ! heat flux W/m2
            if (qwat(ixy,iz)>1d-3) then
               dnbr(ixy,iz) = qchf(ixy,iz)/qwat(ixy,iz)
            else
               dnbr(ixy,iz) = 1000
            endif
         enddo
      enddo

      mdnbr=999.999
      ixyp=1
      izp=1
      do iz=izfuelbot,izfueltop
         do ixy_fa=1,nxy_fa
            ixy=i_fa(ixy_fa)
            if (dnbr(ixy,iz)<mdnbr.and.dnbr(ixy,iz)>1d-10) then
               mdnbr=dnbr(ixy,iz)
               ixyp=ixy
               izp=iz
            endif
         enddo
      enddo
      mqchf=qchf(ixyp,izp)/10000d0
      mqwat=qwat(ixyp,izp)/10000d0

      return
      end subroutine cal_dnbr

#ifdef tuan_tr_test
      subroutine initth_hex
      use Mod_THfunc, only: fdens, fenthal
      use mod_manage, only: alloc_th
      use Inc_Geometry
      use Inc_FA
      use Mod_GetSome, only: Get_POW, Get_LinPOW
!     ! use Mod_H2O, only: H2O_T2D, H2O_T2H
      use Inc_XS_File, only: if_th
      use inc_tpen, only: powlin, plevel0

      implicit none
      integer, allocatable, save :: ineut(:), jneut(:)
      integer :: ia, ja, k, lc, npthy, ji, jb, iabeg, iaend
      integer :: npthx, ii, ib, j, jn, i, in, l
      integer :: nchanp1, kth, ir
      integer :: xx, yy
      real(8) :: chanvf, hin
      real(8) :: rhoin, rhouin, rhohuin, uin
      real(8) :: ptot
      real(8) :: t_inlet, g_inlet
      real(8) :: powfa
      real(8) :: hac


      if_th = 1 ! using rast-k origin steam table

      Ptot = 0
      allocate(ltochan(0:Nxy))
      allocate(lchantol(1:Nxy))
      allocate(lchanptr(0:Nxy))
      ltochan=0
      lchantol=0
      lchanptr=0

      ! determine neutronic node structure
      if (.not.allocated(ineut)) then
         allocate(ineut(0:Nx_1N))
         allocate(jneut(0:Ny_1N))
         ineut(0)=0
         jneut(0)=0
         do ia=1,nx_1n
            ineut(ia)=ineut(ia-1)+mesh_x(ia)
         enddo
         do ja=1,ny_1n
            jneut(ja)=jneut(ja-1)+mesh_y(ja)
         enddo
         if (nzth==0) nzth=nz

         ! channel assignment
         nchan=0
         lc=0
         do ja=1,Ny_1N
            npthy=mesh_y(ja)/nthy(ja)
            do ji=1,nthy(ja)
               jb=jneut(ja-1)+(ji-1)*npthy
               iabeg=ix_startfa_y_1n(ja)
               iaend=ix_endfa_y_1n(ja)
               do ia=iabeg,iaend
                  npthx=Mesh_x(ia)/nthx(ia)
                  do ii=1,nthx(ia)
                     nchan=nchan+1
                     ib=ineut(ia-1)+(ii-1)*npthx
                     j=jb
                     do jn=1,npthy
                        j=j+1
                        i=ib
                        do in=1,npthx
                           i=i+1
                           l=nodel(i,j)
                           ltochan(l)=nchan
                           lc=lc+1
                           lchantol(lc)=l
                        enddo
                        lchanptr(nchan)=lc
                     enddo
                  enddo
               enddo
            enddo
         enddo
         lchanptr(0)=0
         
         call alloc_th
         ! channel to node correspondence for the radial reflector
         nchanp1=nchan+1

         ! initialize constants
         yy=1
         do ja=1,ny_1n
            if (nthy(ja)>yy) yy=nthy(ja)
         enddo
         xx=1
         do ja=1,nx_1n
            if (nthx(ja)>xx) xx=nthy(ja)
         enddo
         chanvf=1d0/real(xx*yy,8)         !channel volume fraction per assembly

         if (opt_core==4.and.nxy==1) chanvf=0.25d0
         chanvf_save=chanvf

         acf         = ( 0.5*1.73205080757*(h_FA*DM2)**2 -PI*(N_Pin * R_Clad**2) )*chanvf  ! coolant flow area
         afp         = N_Pin * PI * R_Fuel**2 * chanvf                                        ! total fuel pellet area
         xi          = 2*PI*(N_Pin*R_Clad)*chanvf                             ! wetted perimeter
         zeta        = N_Pin*2*PI*R_Clad*chanvf                                               ! heated perimeter
         FA_MassFlow = FA_MassFlow * chanvf
         zetap       = zeta/acf                                                                   ! heated perimeter density
         deq         = 4*acf/xi                                                                   ! equivalent diameter
         fracdf      = 1 - Frac_Gamma                                                             ! fuel heat deposit fraction
         if (Opt_mode==3) then
             hac         = 0d0
             do k = IzFuelBot, IzFuelTop
                hac = hac + gridsize_z(k)/100d0
             enddo
             powfa       = Core_Power_100/Nxy_FA
             write(*,*) 'WARNING: need to check here: initth_hex'         
             powlin      = powfa / hac*chanvf !! for test !! need to remove 
             plevel0 = PPower
         endif

         ! radial nodalization (in fuel pin)
         nrp1=nr+1
         nrp2=nr+2
         nrp3=nr+3
         nrp4=nr+4
         nrp5=nr+5
         delr=r_fuel/real(nr,8)
         delr2=delr*delr
         delrw=0.5d0*h_clad
         delrw2=delrw*delrw
         tworm=h_clad/(r_gap+0.5d0*h_clad)
         do i=1,nrp1
            r(i)=delr*(i-1)
         enddo
         r(nrp2)=r_gap
         r(nrp3)=r_gap+delrw
         r(nrp4)=r_clad
         kgap=h_gap*delr
         kgap2=h_gap*h_clad*r_fuel/r_gap
         kgap4=h_gap*h_clad*(4d0-h_clad/r_gap)*r_fuel/r_gap

         ! assign inlet condition to all nodes
         tdopin=sqrt(TF_In+CKELVIN)
         if (if_th==4) then
            din= 1d0 !H2O_T2D(TM_In + DegToK)*1d3
         else
            din=fdens(TM_In)
         endif  
         do l=1,nchan
            do k=1,nzth
               do ir=1,nrp1
                  tfuel(ir,k,l)=TF_In
                  tdopl(k,l)=tdopin
               enddo
               do ir=nrp2,nrp4
                  tfuel(ir,k,l)=TF_In
               enddo
               tcool(k,l)=TM_In
               dcool(k,l)=din
               if (flag_th_chanwise) then
                  t_inlet=chanwise_t_inlet(i_chan_1n(i_4nto1n(i_fa(l))))
                  tcool(k,l)=t_inlet
                  if (if_th==4) then
                     dcool(k,l)= 1d0 ! H2O_T2D(t_inlet+degtok)*1d3
                  else
                     dcool(k,l)=fdens(t_inlet)
                  endif
               endif
            enddo
         enddo
         
         ! assign reflector t/h condition
         do kth=1,nzth
            dcool(kth,nchanp1)=din
            tcool(kth,nchanp1)=TM_In
            tdopl(kth,nchanp1)=tdopin
         enddo

         ! initialize volume and junction variables
         if (if_th==4) then
            hin=1d0 !H2O_T2H(TM_In+degtok)
            rhoin=1d0 !H2O_T2D(TM_In+degtok)*1d3
         else
            hin=fenthal(TM_In)
            rhoin=fdens(TM_In)
         endif  
         rhouin=FA_MassFlow/acf
         rhohuin=rhouin*hin
         uin=rhouin/rhoin

         do l=1,nchan
            k=0
            if (flag_th_chanwise) then
               t_inlet=chanwise_t_inlet(i_chan_1n(i_4nto1n(i_fa(l))))
               g_inlet=chanwise_g_inlet(i_chan_1n(i_4nto1n(i_fa(l))))*chanvf_save
               if (if_th==4) then
                  hin= 1d0 !H2O_T2H(t_inlet+degtok)
                  rhoin=1d0! H2O_T2D(t_inlet+degtok)
               else
                  hin=fenthal(t_inlet)
                  rhoin=fdens(t_inlet)
               endif
               rhouin=g_inlet/acf
               rhohuin=rhouin*hin
               uin=rhouin/rhoin
            endif
            rhou(k,l)=rhouin
            rhohu(k,l)=rhohuin
            u(k,l)=uin
            ud(k,l)=uin
            do k=1,nzth
               hcool(k,l)=hin
               rhou(k,l)=rhouin
               rhohu(k,l)=rhohuin
               u(k,l)=uin
               ud(k,l)=uin
            enddo
         enddo

      endif

      end subroutine initth_hex

      subroutine ssth_hex
      use Mod_THfunc
      USE Inc_XS_File, ONLY: if_th
      use inc_tpen!, only: powlin, plevel0,flag_init_thhex
      use inc_3D, only: power, avg_power
      implicit none
      integer :: ixy,ixy_fa,iz
      real(8) :: qc,qprime,qeffnew,rhohuin,rhohuout
      real(8) :: tfmaxt,toutavgt,fnchan,qf,rhouin
      real(8) :: hout,tout
      real(8) :: rtmp, rvol
      integer(4) :: l,k
      real(8)    :: hac
      real(8)    :: relpow_init
      REAL(8) :: powfa 

      hac = 0d0
      powfa = Core_Power_100/max(nxy_fa,1)



#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [ssth] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      tfmax    = 0
      toutavg  = 0
      tdoplmax = 0


!      Tot_Vol = Tot_FuelVol !! NEED TO MODIFY

      if (opt_mode ==3 .and. flag_init_thhex) then
         LinPOW = 0d0
         do k = IzFuelBot, IzFuelTop
            hac = hac + Gridsize_z(k)
         enddo
         relpow_init = Core_Power_100/max(nxy_fa,1)/(hac*1E-2)*PPower
         do ixy_fa=1,nxy_fa
            ixy = i_fa(ixy_fa)
            do iz = IzFuelBot, IzFuelTop
               LinPOW(ixy,iz) = relpow_init 
            enddo
         enddo
      endif
!
!
      ! coolant temperature calculation
!      do l = 1,nxy
!         ixy=imap(l)
!         ixy_fa = ixytoifa(ixy)
!         if (I_FARF_1N(ixy) .NE. 2) cycle
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         rhouin=rhou(0,ixy_fa)
         rhohuin=rhohu(0,ixy_fa)
         do iz=1,nz
            if(opt_mode==3 ) then
                if (flag_init_thhex) then
                  qprime=LinPOW(ixy,iz) ! W/m
                else
                  qprime = power(ixy, iz)/Avg_power*powlin*ppower 
                
                endif
            else
                qprime=LinPOW(ixy,iz) ! W/m
            endif

            qf=fracdf*qprime/afp
            qvol(iz,ixy_fa)=qf
            qc=Frac_Gamma*qprime/acf
            qflux(iz,ixy_fa)=qf*afp/zeta
            qeffnew=qflux(iz,ixy_fa)*zetap+qc
            qeff(iz,ixy_fa)=qeffnew
            rhohuout=rhohuin+qeffnew*(meshsize_z(iz)*1d-2)
            hcool(iz,ixy_fa)=0.5d0*(rhohuout+rhohuin)/rhouin
            if (if_th==4) then
               tcool(iz,ixy_fa)=1d0 !H2O_H2T(hcool(iz,ixy_fa))-degtok
               dcool(iz,ixy_fa)=1d0 !H2O_H2D(hcool(iz,ixy_fa))*1d3
            else
               tcool(iz,ixy_fa)=ftemp(hcool(iz,ixy_fa))
               dcool(iz,ixy_fa)=fdensh(hcool(iz,ixy_fa))
            endif
            rhohuin=rhohuout
            rhohu(iz,ixy_fa)=rhohuout
         enddo
         hout=rhohu(nzth,ixy_fa)/rhou(nzth,ixy_fa)
         if (if_th==4) then
            tout=1d0 !H2O_H2T(hout)-degtok
         else
            tout=ftemp(hout)
         endif
         if (idomz==ndomz) toutavg=toutavg+tout
      enddo

      ! coolant temperature for reflector
      do iz=1,nz
         rtmp=0d0
         rvol=0d0
         do ixy_fa=1,nxy_fa
            ixy=i_fa(ixy_fa)
            !rtmp=rtmp+tcool(iz,ixy_fa)*nodevolume(ixy,iz)
            rtmp=rtmp+hcool(iz,ixy_fa)*nodevolume(ixy,iz)
            rvol=rvol+nodevolume(ixy,iz)
         enddo
         !tcool(iz,nxy_fa+1)=rtmp/max(1d-10,rvol)
         hcool(iz,nxy_fa+1)=rtmp/max(1d-10,rvol)
         tcool(iz,nxy_fa+1)=ftemp(hcool(iz,nxy_fa+1))
         dcool(iz,nxy_fa+1)=fdensh(hcool(iz,nxy_fa+1))
      enddo

      ! fuel temperature calculation
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         do iz=izfuelbot,izfueltop
            if (opt_tfutab) then
               call get_tfutab(iz,ixy_fa,tcool(iz,ixy_fa))
               tfmax=max(tfmax,tfuel(1,iz,ixy_fa))
            else
               qprime=LinPOW(ixy,iz)
               qf=fracdf*qprime/afp    
               ! heat transfer coefficient
               htcoef(iz,ixy_fa)=fhtcoef(tcool(iz,ixy_fa),deq,rhou(iz,ixy_fa))
               call tfcalss(iz,ixy_fa,tcool(iz,ixy_fa),htcoef(iz,ixy_fa),qf)
               tfmax=max(tfmax,tfuel(1,iz,ixy_fa))
            endif
         enddo
      enddo

      flagth=.false.
      if (tdoplmax<EPS_DoppTF) flagth=.true.
      tfmaxt=tfmax
      toutavgt=toutavg
      fnchan=nchan
      tfmax=tfmaxt
      fnchan=fnchan/float(ndomz)
      toutavg=toutavgt/fnchan
      write(*,601) tfmax,toutavg,tdoplmax,flagth

601   FORMAT(" Max. Tf=",f7.1,"C,  Avg. Outlet Temp.=",f6.2,"C"         &
           , ", Max. Doppler Change=",1p,e9.2,l2)

      return
      end subroutine ssth_hex


      SUBROUTINE trth_hex
      USE Inc_FluxVar
      USE Mod_THfunc !, ONLY: fhtcoef, ftemp
      USE Inc_FA , ONLY: h_Clad
!     ! USE MOD_H2O,  ONLY: H2O_H2T
      USE Inc_XS_File, ONLY: if_th
      use inc_tpen, only: imap, powlin, plevel0
      use Inc_3D, only: power

      IMPLICIT NONE
      INTEGER :: Ixy, Ixy_FA, Iz
      REAL(8) :: rhouin,rhohuin,rhouinn,rhohuinn,qc, tb, qf
      REAL(8) :: hout, tout,qfn,qeffn,qprime,dz
      integer(4) :: l!,k

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [trth] in Mod_THdriver'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      tfmax   = 0
      dr2odt  = delr2/dT
      tw2odt  = h_Clad**2 / dT
      toutavg = 0

!      do l = 1,nxy
!         ixy = imap(l)
!         if (I_FARF_1N(ixy) .NE. 2) cycle
!         Ixy_FA = ixytoifa(ixy) 
      do ixy_fa=1,nxy_fa
         ixy=i_fa(ixy_fa)
         
         rhouin   = rhou(0, Ixy_FA)
         rhohuin  = rhohu(0, Ixy_FA)
         rhouinn  = rhouin
         rhohuinn = rhohuin
         do Iz = IzFuelBot,IzFuelTop  
            dz=meshsize_z(Iz)*1d-2
            qprime = power(ixy, iz)*powlin*plevel0 
            qfn=fracdf*qprime/afp
            qf=cetafb*qvol(Iz, Ixy_FA)+cetaf*qfn
            qvol(Iz, Ixy_FA)=qfn
            htcoef(Iz, Ixy_FA)=fhtcoef(tcool(Iz, Ixy_FA),deq,rhou(Iz, Ixy_FA))
            tb=tcool(Iz, Ixy_FA)
            if (iz<izfuelbot.or.iz>izfueltop) then
               tfuel(:,iz,ixy_fa)=tcool(iz,ixy_fa)
            endif
            CALL tfcaltr_hex(Iz,Ixy_FA,tb,htcoef(Iz, Ixy_FA),qf,.TRUE.)
            tfmax=MAX(tfmax,tfuel(1,Iz ,Ixy_FA))
            qflux(iz,ixy_fa)=htcoef(Iz, Ixy_FA)*(tfuel(nrp4,Iz ,Ixy_FA)-tb)
            dzcetadt=dz/(cetac*dT)
            qc    = Frac_Gamma * qprime / acf            
            qeffn = qflux(iz,ixy_fa) * zetap + qc
            qc=qeffn+cetacr*qeff(Iz, Ixy_FA)            
            qeff(Iz, Ixy_FA)=qeffn
            CALL tccaltr_hex(Iz,Ixy_FA,rhouin,rhohuin,rhouinn,rhohuinn,qc)
         enddo
         hout = rhohu(IzFuelTop, Ixy_FA) / rhou(IzFuelTop, Ixy_FA)
         if (if_th==4) then
            tout = 1 !H2O_H2T(hout)-degtok
         else
            tout = ftemp(hout)
         endif
         if (idomz==ndomz) toutavg=toutavg+tout
      enddo

      toutavg = toutavg / float(nchan)

      RETURN
      END SUBROUTINE trth_hex
#endif


      END MODULE Mod_THdriver

