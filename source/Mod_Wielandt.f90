
      MODULE Mod_Wielandt
      use ieee_arithmetic ! @$^ siarhei_fr for using nvfortran (ieee_is_nan is not working)


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

      CONTAINS

#ifdef siarhei_delete 
      subroutine wiel(icy)
      use Inc_FluxVar
      use Inc_Lscoef
      use Mod_SolLS
      use Inc_Geometry
      use Inc_3D
      use Inc_Control
      use Inc_Option, ONLY: N_Group
      use mod_charedit, only: print_msg
      use inc_extsrc, only: flag_extsrc
!      use mod_extsrc, only: set_fsrc_extsrc
      implicit none
      integer :: k,l,m,lf,icy
      real(8) :: err,eigshft,errlinfs
      real(8) :: gamman,gammad,gamma
      real(8) :: sumf,summ, psi_max

      if (flag_extsrc) then
!         call set_fsrc_extsrc
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         return
      endif

      flagl2=.false.
      flaglinf=.false.
      flageig=.false.

      ! compute new fission source and corresponding integral quantities
      errl2d=errl2
      errlinf=0
      gamman=0
      gammad=0
      errl2=0
      psi_max=0.0

      do k=IzFuelBot,IzFuelTop
         do lf=1,Nxy_FA_P2R
            l=nodef(lf)
            psi_max=max(psi_max,abs(FisSrc_Iout(l,k)))
         enddo
      enddo

      do k=IzFuelBot,IzFuelTop
         do lf=1,Nxy_FA_P2R
            l=nodef(lf)
            FisSrc_Iout(l,k)=FisSrc(l,k)
            FisSrc(l,k)=0d0
            do m=1,N_Group
               FisSrc(l,k)=FisSrc(l,k)+af(m,l,k)*Flux(l,k,m)
            enddo
            err=FisSrc(l,k)-FisSrc_Iout(l,k)
            errlinfs=errlinf
            if (FisSrc(l,k)>psi_max*0.05d0) then
               errlinf=max(errlinfs,abs(err/FisSrc(l,k)))
            endif
            gammad=gammad+FisSrc_Iout(l,k)*FisSrc(l,k)
            gamman=gamman+FisSrc(l,k)*FisSrc(l,k)
            errl2=errl2+err*err
         enddo
      enddo

      ! compute new eigenvalue
      eigvd = keff
      if (icy<=1) then
!         call AxB(Flux,flux_add)
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         sumf=0d0
         summ=0d0
         do k=1,Nz
            do l=1,Nxy
               sumf=sumf+FisSrc(l,k)
               summ=summ+FisSrc(l,k)*reigvs  !N_Group,temp
               do m=1,N_Group
                  summ=summ+flux_add(l,k,m)
               enddo
            enddo
         enddo
         keff=sumf/summ
      else
         gamma=gammad/gamman
         keff=1d0/(reigv*gamma+(1d0-gamma)*reigvs)
      endif
      reigv=1d0/keff
      erreig=abs(keff-eigvd)

      if (erreig<EPS_keff_St) flageig=.true.

      ! compute fission source errors and estimate the dominance ratio :
      gammad=abs(gammad)
      rerrl2=sqrt(errl2/gammad)
      domr=sqrt(errl2/errl2d)
      if (domr>10) domr=10
      if (rerrl2<EPS_Global) flagl2=.true.
      if (errlinf<EPS_Local) flaglinf=.true.

      ! shift eigenvalue
      if (icy<0) then
         eigshft=eshift0
      else
         eigshft=eshift
      endif
      eigvs=keff+eigshft
      reigvsd=reigvs
      reigvs=1d0/eigvs

      ! check solution is NaN
      do k=1,nz
         do l=1,nxy
            if (ieee_is_nan(fissrc(l,k))) then
               call print_msg(3,'Solution diverged ... psi')
               stop
            endif
         enddo
      enddo
      if (ieee_is_nan(keff)) then
         call print_msg(3,'Solution diverged ... keff')
         stop
      endif

      return
      end subroutine wiel
#endif 

#ifdef jr_vver
      subroutine wiel_hex(icy)
      use Inc_FluxVar
      use Inc_Lscoef
      use Mod_SolLS
      use Inc_Geometry
      use Inc_3D
      use Inc_Control
      use Inc_Option, ONLY: N_Group
      use mod_charedit, only: print_msg
      use inc_extsrc, only: flag_extsrc
!      use mod_extsrc, only: set_fsrc_extsrc
      use Inc_TPEN, only: if_hexgometry,nodef_hex
      implicit none
      integer :: k,l,m,lf,icy
      real(8) :: err,eigshft,errlinfs
      real(8) :: gamman,gammad,gamma
      real(8) :: sumf,summ, psi_max


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [wiel_hex] in Mod_Wielandt'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (flag_extsrc) then
!         call set_fsrc_extsrc
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
         return
      endif

      flagl2=.false.
      flaglinf=.false.
      flageig=.false.

      ! compute new fission source and corresponding integral quantities
      errl2d=errl2
      errlinf=0
      gamman=0
      gammad=0
      errl2=0
      psi_max=0.0

#if tuan_fr_crm
      do k=IzFuelBot,IzFuelTop
#else
      do k=1,Nz
#endif
         do lf=1,Nxy_FA_P2R
            l=imap(nodef_hex(lf))
            psi_max=max(psi_max,abs(FisSrc_Iout(l,k)))
         enddo
      enddo
#if tuan_fr_crm
      do k=IzFuelBot,IzFuelTop
#else
      do k=1,Nz
#endif
         do lf=1,Nxy_FA_P2R
            l=imap(nodef_hex(lf))
            FisSrc_Iout(l,k)=FisSrc(l,k)
            FisSrc(l,k)=0d0
            do m=1,2 !N_Group
               FisSrc(l,k)=FisSrc(l,k)+af(m,nodef_hex(lf),k)*Flux(l,k,m)
            enddo
            if (ieee_is_nan(fissrc(l,k))) then
                write(*,*) 'l, k', l,k, nodef_hex(lf),k
               call print_msg(3,'Solution diverged2 ... psi')
            write(*,*)  k, l, lf, 'af: ',af(:,nodef_hex(lf),k)
            write(*,*)  k, l, lf, 'flux: ',Flux(l,k,:)
!               stop
            endif
!            write(*,*)  k, l, lf, 'af: ',af(:,nodef_hex(lf),k)
!            write(*,*)  k, l, lf, 'flux: ',Flux(l,k,:)
            err=FisSrc(l,k)-FisSrc_Iout(l,k)
            errlinfs=errlinf
            if (FisSrc(l,k)>psi_max*0.05d0) then
               errlinf=max(errlinfs,abs(err/FisSrc(l,k)))
            endif
            gammad=gammad+FisSrc_Iout(l,k)*FisSrc(l,k)
            gamman=gamman+FisSrc(l,k)*FisSrc(l,k)
            errl2=errl2+err*err
         enddo
      enddo
!write(*,*) gamman, gammad

      ! compute new eigenvalue
      eigvd = keff
      if (icy<=1) then
         call AxBhex(Flux,flux_add)
         sumf=0d0
         summ=0d0
         do k=1,Nz
            do l=1,Nxy
               sumf=sumf+FisSrc(l,k)
               summ=summ+FisSrc(l,k)*reigvs  !N_Group,temp
               do m=1,N_Group
                  summ=summ+flux_add(l,k,m)
               enddo
            enddo
         enddo
         keff=sumf/summ
      else
         gamma=gammad/gamman
         keff=1d0/(reigv*gamma+(1d0-gamma)*reigvs)
      endif

      reigv=1d0/keff
      erreig=abs(keff-eigvd)

      if (erreig<EPS_keff_St) flageig=.true.

      ! compute fission source errors and estimate the dominance ratio :
      gammad=abs(gammad)

      rerrl2=sqrt(errl2/gammad)

      domr=sqrt(errl2/errl2d)

      if (domr>10) domr=10

      if (rerrl2<EPS_Global) flagl2=.true.

      if (errlinf<EPS_Local) flaglinf=.true.

      ! shift eigenvalue
      if (icy<0) then
         eigshft=eshift0
      else
         eigshft=eshift
      endif

      eigvs=keff+eigshft

      reigvsd=reigvs

      reigvs=1d0/eigvs

      ! check solution is NaN
      do k=1,nz
         do l=1,nxy
            if (ieee_is_nan(fissrc(l,k))) then
                write(*,*) 'l, k', l,k
               call print_msg(3,'Solution diverged2 ... psi')
               stop
            endif
         enddo
      enddo
      if (ieee_is_nan(keff)) then
         call print_msg(3,'Solution diverged ... keff')
         stop
      endif

      return
      end subroutine wiel_hex

#endif

      END MODULE Mod_Wielandt
