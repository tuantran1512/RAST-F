#ifdef siarhei_delete

! ============================================================================
! module Mod_ExtSrc
! ============================================================================
!> @brief
!> set external source density for transient analysis
!> @date
!> 2020.03.10
      module Mod_ExtSrc
!      use Mod_Depletion, only: HeavyNuclide, BurnableAbsorber, FissionProduct
      use Inc_Option, only: N_Group
      use Inc_Time, only: dT ! dT = depletion time step [sec]
      use Inc_3D, only: Flux, Flux_Old, FisSrc, BU, MTU
      use Inc_ExtSrc
      use Inc_Nuclide
      use Inc_RP, only: I_FARF_1N_3D
      use Inc_Geometry, only: izFuelTop, izFuelBot, nz, nxy, NodeVolume, ny, nxy_1N, &
                            & Ix_Start_y, Ix_End_y, IxIytoIxy, Ix_4Nto1N, Iy_4Nto1N, &
                            & ixytoix, ixytoiy, IxIy_1NToIxy_1N, IxIyToIxy_1N, &
                            & Ix_Start_y_1N, Ix_End_y_1N, Ny_1N, Nx_1N
      use inc_parallel, only: comm

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

      contains

      !===============================================================================
      ! Set node-wise volume-integrated external source
      !===============================================================================
      subroutine Set_ExtSrc
      implicit none
      integer :: ix, iy, ixy, iz, ix_1N, iy_1N, ixy_1N
      real(8) :: tot_volume
      real(8) :: tot_N_U34,  tot_N_U35,  tot_N_U36,  tot_N_U37,  tot_N_U38,  &
      &          tot_N_Np37, tot_N_Np38, tot_N_Np39, tot_N_Pu38, tot_N_Pu39, tot_N_Pu40, &
      &          tot_N_Pu41, tot_N_Pu42, tot_N_Pu43, tot_N_Am41, tot_N_As42, tot_N_Am42, &
      &          tot_N_Am43, tot_N_Am44, tot_N_Cm42, tot_N_Cm43, tot_N_Cm44, &
      &          tot_N_I35,  tot_N_Xe35, tot_N_Nd47, tot_N_Nd48, tot_N_Nd49, &
      &          tot_N_Pm47, tot_N_Ps48, tot_N_Pm48, tot_N_Pm49, tot_N_Sm47, &
      &          tot_N_Sm48, tot_N_Sm49
      real(8), allocatable :: buf_Flux(:,:,:)
      real(8), allocatable :: buf_Flux_Old(:,:,:)
      real(8) :: buf_dT
      real(8) :: tmp_sf     !> temp. spontaneous fission source
      real(8) :: tmp_an     !> temp. (a,n) source
      real(8) :: tmp_omega  !> slowing down parameter
      real(8) :: tmp_N_Gd   !> temporary number density Gd
      real(8) :: tmp_N_U    !> temporary number density U
      real(8) :: tmp_N_O    !> temporary number density O
      real(8) :: tmp_N_O16  !> temporary number density O16
      real(8) :: tmp_N_O17  !> temporary number density O17
      real(8) :: tmp_N_O18  !> temporary number density O18
      real(8) :: tmp_N_fiss  !> temporary number density of fission
      real(8), parameter :: convert_GWD_to_nfiss = 2.69633130E+24  !> unit conversion factor for GWD to number of fissions

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Set_ExtSrc] in Mod_ExtSrc'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Set_ExtSrc] in Mod_ExtSrc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

#ifdef js_r2mpi
      if (.not.comm%usempi.and..not.comm%if_master) return
#endif

      if (.not.flag_extsrc) return

      ! allocate ExtSrc
      if (allocated(ExtSrc)) then
         !call fatal_error("ExtSrc is already allocated.")
      else
         allocate( ExtSrc       (nxy,nz) )
         allocate( ExtSrc_Z     (nz)     )
         allocate( ExtSrc_XY    (nxy)    )
         allocate( ExtSrc_XY_1N (nxy_1n) )
      endif

      ! initialize ExtSrc
      tot_ExtSrc = 0.d0
      tot_ExtSrc_sf = 0.d0
      tot_ExtSrc_an = 0.d0
      ExtSrc = 0.d0
      ExtSrc_Z = 0.d0
      ExtSrc_XY = 0.d0
      ExtSrc_XY_1N = 0.d0

      ! for debug
      tot_N_U34  = 0.d0
      tot_N_U35  = 0.d0
      tot_N_U36  = 0.d0
      tot_N_U37  = 0.d0
      tot_N_U38  = 0.d0
      tot_N_Np37 = 0.d0
      tot_N_Np38 = 0.d0
      tot_N_Np39 = 0.d0
      tot_N_Pu38 = 0.d0
      tot_N_Pu39 = 0.d0
      tot_N_Pu40 = 0.d0
      tot_N_Pu41 = 0.d0
      tot_N_Pu42 = 0.d0
      tot_N_Pu43 = 0.d0
      tot_N_Am41 = 0.d0
      tot_N_As42 = 0.d0  ! Am-242m
      tot_N_Am42 = 0.d0  ! Am-242
      tot_N_Am43 = 0.d0
      tot_N_Am44 = 0.d0
      tot_N_Cm42 = 0.d0
      tot_N_Cm43 = 0.d0
      tot_N_Cm44 = 0.d0
      tot_N_I35  = 0.d0
      tot_N_Xe35 = 0.d0
      tot_N_Nd47 = 0.d0
      tot_N_Nd48 = 0.d0
      tot_N_Nd49 = 0.d0
      tot_N_Pm47 = 0.d0
      tot_N_Ps48 = 0.d0
      tot_N_Pm48 = 0.d0
      tot_N_Pm49 = 0.d0
      tot_N_Sm47 = 0.d0
      tot_N_Sm48 = 0.d0
      tot_N_Sm49 = 0.d0

      tot_volume = 0.d0

      ! backup flux and time step
      if(.not.allocated(buf_Flux_Old)) allocate(buf_Flux_Old(nxy, nz, n_group))
      if(.not.allocated(buf_Flux))     allocate(buf_Flux(nxy, nz, n_group))
      buf_Flux_Old = Flux_Old
      buf_Flux     = Flux
      buf_dT       = dT

      ! set flux as zero for cooling period
      Flux_Old = 0.d0
      Flux = 0.d0
      dT = cool_time

      ! depletion calculation with zero flux for cooling period
      call HeavyNuclide
      call BurnableAbsorber
      call FissionProduct(.TRUE.)

      ! roll-back flux and time step data
      Flux_Old = buf_Flux_Old
      Flux     = buf_Flux
      dT       = buf_dT
      deallocate(buf_Flux_Old)
      deallocate(buf_Flux)

      ! set ExtSrc
      do iz=izFuelBot,izFuelTop
         do ixy=1,nxy
            ix=ixytoix(ixy)
            iy=ixytoiy(ixy)
            ix_1N=Ix_4Nto1N(ix)
            iy_1N=Iy_4Nto1N(iy)
            ixy_1N=IxIytoIxy_1N(ix,iy)

            if ( I_FARF_1N_3D(ixy_1N, iz) == 0 .or. I_FARF_1N_3D(ixy_1N, iz) == 1 ) then
               ! Non Fuel Region
               tmp_sf = 0.d0
               tmp_an = 0.d0
            else
               ! Fuel Region
               ! spontaneous fission
               tmp_sf=( nu_sf_U34  * sf_br_U34 * dkcon_U34  * N_U34 (ixy,iz) + &
                     &  nu_sf_U35  * sf_br_U35 * dkcon_U35  * N_U35 (ixy,iz) + &
                     &  nu_sf_U36  * sf_br_U36 * dkcon_U36  * N_U36 (ixy,iz) + &
                     &  nu_sf_U38  * sf_br_U38 * dkcon_U38  * N_U38 (ixy,iz) + &
                     &  nu_sf_Np37 * sf_br_Np37* dkcon_Np37 * N_Np37(ixy,iz) + &
                     &  nu_sf_Pu38 * sf_br_Pu38* dkcon_Pu38 * N_Pu38(ixy,iz) + &
                     &  nu_sf_Pu39 * sf_br_Pu39* dkcon_Pu39 * N_Pu39(ixy,iz) + &
                     &  nu_sf_Pu40 * sf_br_Pu40* dkcon_Pu40 * N_Pu40(ixy,iz) + &
                     &  nu_sf_Pu41 * sf_br_Pu41* dkcon_Pu41 * N_Pu41(ixy,iz) + &
                     &  nu_sf_Pu42 * sf_br_Pu42* dkcon_Pu42 * N_Pu42(ixy,iz) + &
                     &  nu_sf_Am41 * sf_br_Am41* dkcon_Am41 * N_Am41(ixy,iz) + &
                     &  nu_sf_As42 * sf_br_As42* dkcon_As42 * N_As42(ixy,iz) + &
                     &  nu_sf_Am43 * sf_br_Am43* dkcon_Am43 * N_Am43(ixy,iz) + &
                     &  nu_sf_Cm42 * sf_br_Cm42* dkcon_Cm42 * N_Cm42(ixy,iz) + &
                     &  nu_sf_Cm43 * sf_br_Cm43* dkcon_Cm43 * N_Cm43(ixy,iz) + &
                     &  nu_sf_Cm44 * sf_br_Cm44* dkcon_Cm44 * N_Cm44(ixy,iz) ) &
                     &  * NodeVolume(ixy,iz) * extsrc_multiplier

               ! O16, O17, O18
               tmp_N_U = N_U34(ixy,iz)  + N_U35(ixy,iz)  + N_U36(ixy,iz)  + N_U38 (ixy,iz) + &
                       & N_Np37(ixy,iz) + N_Pu38(ixy,iz) + N_Pu39(ixy,iz) + N_Pu40(ixy,iz) + &
                       & N_Pu41(ixy,iz) + N_Pu42(ixy,iz) + N_Am41(ixy,iz) + N_As42(ixy,iz) + &
                       & N_Am43(ixy,iz) + N_Cm42(ixy,iz) + N_Cm43(ixy,iz) + N_Cm44(ixy,iz)

               tmp_N_Gd = N_Gd52(ixy,iz) +  N_Gd54(ixy,iz) + N_Gd55(ixy,iz) + N_Gd56(ixy,iz) + &
                        & N_Gd57(ixy,iz) +  N_Gd58(ixy,iz) + N_Gd60(ixy,iz)

               tmp_N_fiss = convert_GWD_to_nfiss * BU(ixy,iz) * MTU(ixy,iz) / NodeVolume(ixy,iz)

               tmp_N_O = 2*tmp_N_U + 1.5*tmp_N_Gd + 2*tmp_N_fiss
               tmp_N_O16 = tmp_N_O * NA_O16
               tmp_N_O17 = tmp_N_O * NA_O17
               tmp_N_O18 = tmp_N_O * NA_O18

               ! calculate slowing down parameter
               tmp_omega = (   tmp_N_O16      * sqrtZ_O16  + &
                             & tmp_N_O17      * sqrtZ_O17  + &
                             & tmp_N_O18      * sqrtZ_O18  + &
                             & N_I35 (ixy,iz) * sqrtZ_I35  + &
                             & N_Xe35(ixy,iz) * sqrtZ_Xe35 + &
                             & N_Nd47(ixy,iz) * sqrtZ_Nd47 + &
                             & N_Nd48(ixy,iz) * sqrtZ_Nd48 + &
                             & N_Nd49(ixy,iz) * sqrtZ_Nd49 + &
                             & N_Pm47(ixy,iz) * sqrtZ_Pm47 + &
                             & N_Ps48(ixy,iz) * sqrtZ_Ps48 + &
                             & N_Pm48(ixy,iz) * sqrtZ_Pm48 + &
                             & N_Pm49(ixy,iz) * sqrtZ_Pm49 + &
                             & N_Sm47(ixy,iz) * sqrtZ_Sm47 + &
                             & N_Sm48(ixy,iz) * sqrtZ_Sm48 + &
                             & N_Sm49(ixy,iz) * sqrtZ_Sm49 + &
                             & N_Gd52(ixy,iz) * sqrtZ_Gd52 + &
                             & N_Gd54(ixy,iz) * sqrtZ_Gd54 + &
                             & N_Gd55(ixy,iz) * sqrtZ_Gd55 + &
                             & N_Gd56(ixy,iz) * sqrtZ_Gd56 + &
                             & N_Gd57(ixy,iz) * sqrtZ_Gd57 + &
                             & N_Gd58(ixy,iz) * sqrtZ_Gd58 + &
                             & N_Gd60(ixy,iz) * sqrtZ_Gd60 + &
                             & N_U34 (ixy,iz) * sqrtZ_U34  + &
                             & N_U35 (ixy,iz) * sqrtZ_U35  + &
                             & N_U36 (ixy,iz) * sqrtZ_U36  + &
                             & N_U37 (ixy,iz) * sqrtZ_U37  + &
                             & N_U38 (ixy,iz) * sqrtZ_U38  + &
                             & N_Np37(ixy,iz) * sqrtZ_Np37 + &
                             & N_Np38(ixy,iz) * sqrtZ_Np38 + &
                             & N_Np39(ixy,iz) * sqrtZ_Np39 + &
                             & N_Pu38(ixy,iz) * sqrtZ_Pu38 + &
                             & N_Pu39(ixy,iz) * sqrtZ_Pu39 + &
                             & N_Pu40(ixy,iz) * sqrtZ_Pu40 + &
                             & N_Pu41(ixy,iz) * sqrtZ_Pu41 + &
                             & N_Pu42(ixy,iz) * sqrtZ_Pu42 + &
                             & N_Pu43(ixy,iz) * sqrtZ_Pu43 + &
                             & N_Am41(ixy,iz) * sqrtZ_Am41 + &
                             & N_Am42(ixy,iz) * sqrtZ_Am42 + &
                             & N_As42(ixy,iz) * sqrtZ_As42 + &
                             & N_Am43(ixy,iz) * sqrtZ_Am43 + &
                             & N_Am44(ixy,iz) * sqrtZ_Am44 + &
                             & N_Cm42(ixy,iz) * sqrtZ_Cm42 + &
                             & N_Cm43(ixy,iz) * sqrtZ_Cm43 + &
                             & N_Cm44(ixy,iz) * sqrtZ_Cm44   )


               tmp_omega = 1.d0/(1.866E13 * tmp_omega)

               ! (a,n) reaction
               tmp_an = (&
                         & dkcon_U34 *alpa_br_U34 *N_U34(ixy,iz)  * ( tmp_N_O17*tau_O17_U34  + tmp_N_O18*tau_O18_U34  ) + &
                         & dkcon_U35 *alpa_br_U35 *N_U35(ixy,iz)  * ( tmp_N_O17*tau_O17_U35  + tmp_N_O18*tau_O18_U35  ) + &
                         & dkcon_U36 *alpa_br_U36 *N_U36(ixy,iz)  * ( tmp_N_O17*tau_O17_U36  + tmp_N_O18*tau_O18_U36  ) + &
                         & dkcon_U38 *alpa_br_U38 *N_U38(ixy,iz)  * ( tmp_N_O17*tau_O17_U38  + tmp_N_O18*tau_O18_U38  ) + &
                         & dkcon_Np37*alpa_br_Np37*N_Np37(ixy,iz) * ( tmp_N_O17*tau_O17_Np37 + tmp_N_O18*tau_O18_Np37 ) + &
                         & dkcon_Pu38*alpa_br_Pu38*N_Pu38(ixy,iz) * ( tmp_N_O17*tau_O17_Pu38 + tmp_N_O18*tau_O18_Pu38 ) + &
                         & dkcon_Pu39*alpa_br_Pu39*N_Pu39(ixy,iz) * ( tmp_N_O17*tau_O17_Pu39 + tmp_N_O18*tau_O18_Pu39 ) + &
                         & dkcon_Pu40*alpa_br_Pu40*N_Pu40(ixy,iz) * ( tmp_N_O17*tau_O17_Pu40 + tmp_N_O18*tau_O18_Pu40 ) + &
                         & dkcon_Pu41*alpa_br_Pu41*N_Pu41(ixy,iz) * ( tmp_N_O17*tau_O17_Pu41 + tmp_N_O18*tau_O18_Pu41 ) + &
                         & dkcon_Pu42*alpa_br_Pu42*N_Pu42(ixy,iz) * ( tmp_N_O17*tau_O17_Pu42 + tmp_N_O18*tau_O18_Pu42 ) + &
                         & dkcon_Am41*alpa_br_Am41*N_Am41(ixy,iz) * ( tmp_N_O17*tau_O17_Am41 + tmp_N_O18*tau_O18_Am41 ) + &
                         & dkcon_As42*alpa_br_As42*N_As42(ixy,iz) * ( tmp_N_O17*tau_O17_As42 + tmp_N_O18*tau_O18_As42 ) + &
                         & dkcon_Am43*alpa_br_Am43*N_Am43(ixy,iz) * ( tmp_N_O17*tau_O17_Am43 + tmp_N_O18*tau_O18_Am43 ) + &
                         & dkcon_Cm42*alpa_br_Cm42*N_Cm42(ixy,iz) * ( tmp_N_O17*tau_O17_Cm42 + tmp_N_O18*tau_O18_Cm42 ) + &
                         & dkcon_Cm43*alpa_br_Cm43*N_Cm43(ixy,iz) * ( tmp_N_O17*tau_O17_Cm43 + tmp_N_O18*tau_O18_Cm43 ) + &
                         & dkcon_Cm44*alpa_br_Cm44*N_Cm44(ixy,iz) * ( tmp_N_O17*tau_O17_Cm44 + tmp_N_O18*tau_O18_Cm44 )   &
                         & ) * tmp_omega * NodeVolume(ixy,iz) * extsrc_multiplier
            endif

            ! calculate external source
            ExtSrc(ixy,iz) = tmp_sf + tmp_an
            ExtSrc_Z(iz) = ExtSrc_Z(iz) + ExtSrc(ixy,iz)
            ExtSrc_XY(ixy) = ExtSrc_XY(ixy) + ExtSrc(ixy,iz)
            ExtSrc_XY_1N(ixy_1N) = ExtSrc_XY_1N(ixy_1N) + ExtSrc(ixy,iz)
            tot_ExtSrc = tot_ExtSrc + ExtSrc(ixy,iz)
            tot_ExtSrc_sf = tot_ExtSrc_sf + tmp_sf
            tot_ExtSrc_an = tot_ExtSrc_an + tmp_an
         enddo
      enddo

      return
      end subroutine Set_ExtSrc


      !===============================================================================
      ! Set fission source for external source problem
      !===============================================================================
      subroutine Set_Fsrc_ExtSrc
      use Inc_FluxVar
      use Inc_Lscoef
      use Mod_SolLS
      use Inc_Geometry
      use Inc_3D
      use Inc_Control
      implicit none
      integer :: k,l,m,lf
      real(8) :: err,errlinfs,psi_max,gamman,gammad

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Set_Fsrc_ExtSrc] in Mod_ExtSrc'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Set_Fsrc_ExtSrc] in Mod_ExtSrc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      flagl2  = .FALSE.
      flaglinf= .FALSE.

      ! setting for fixed source problem
      keff = 1.d0
      reigv = 1.d0    ! 1/keff
      reigvs = 1.d0   ! 1/(keff + delk)
      reigvsd = 1.d0  ! reigvs at previous outer iteration

      ! do ont update eigenvalue (keff is fixed as 1.0)
      flageig = .TRUE. ! eigenvalue convergence flag
      erreig  = 0.d0   ! eigenvalue error

      ! compute new fission source and corresponding integral quantities
      errl2d=errl2
      errlinf=0
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

      ! compute fission source errors and estimate the dominance ratio
      gammad=abs(gammad)
      rerrl2=sqrt(errl2/gammad)
      domr=SQRT(errl2/errl2d)
      if (domr.GT.10) domr=10
      if (rerrl2<=EPS_Global) flagl2=.true.
      if (errlinf<=EPS_Local) flaglinf=.true.

      return
      end subroutine Set_Fsrc_ExtSrc


      !===============================================================================
      ! Write steady-state
      !===============================================================================
      subroutine Write_ExtSrc
      use Inc_File, only: Len_INP, Name_INP
      implicit none
      character(len=:), allocatable :: fname ! output file for external source distribution
      integer :: len_fname                   ! length of output file name
      integer :: unit_fname = 20191204
      real(8), allocatable :: FisSrc_2D(:)
      integer :: ix, iy, ixy, iz, ix_1N, iy_1N, ixy_1N


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Write_ExtSrc] in Mod_ExtSrc'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Write_ExtSrc] in Mod_ExtSrc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      ! write external source distribution
      len_fname = Len_INP + 4
      allocate( character(len=len_fname) :: fname )
      fname(1        :Len_INP  ) = Name_INP( 1: Len_INP )
      fname(Len_INP+1:Len_INP+4) = ".ext"

      open(unit_fname, file=trim(fname), action="write", status="replace")
      write(unit_fname,*) "out1 :: Total External Source [#/sec]", tot_ExtSrc
      write(unit_fname,*) "        Spontaneous Fission   [#/sec]", tot_ExtSrc_sf
      write(unit_fname,*) "        Alpha N Reaction      [#/sec]", tot_ExtSrc_an

      write(unit_fname,*) "out2 :: Axial External Source [#/sec]"
      do iz = izFuelBot, izFuelTop
         write(unit_fname,*) iz, ExtSrc_Z(iz)
      enddo

      write(unit_fname,*) "out3 :: Assembly-Wise External Source [#/sec]"
      do iy_1N = 1, Ny_1N
         do ix_1N = Ix_Start_y_1N(Iy_1N), Ix_End_y_1N(Iy_1N)
            ixy_1N = IxIy_1NToIxy_1N(Ix_1N,Iy_1N)
            write(unit_fname,*) ix_1N, iy_1N, ixy_1N, ExtSrc_XY_1N(ixy_1N)
         enddo
      enddo

      write(unit_fname,*) "out4 :: 2-D External Source [#/sec]"
      do ixy = 1, nxy
         ix = ixytoix(ixy)
         iy = ixytoiy(ixy)
         write(unit_fname,*) ix, iy, ExtSrc_XY(ixy)
      enddo

      write(unit_fname,*) "out5 :: 3-D External Source [#/sec]"
      do iz = izFuelBot, izFuelTop
         do ixy = 1, nxy
            write(unit_fname, '(2i5, 10ES13.5)') iz, ixy, NodeVolume(ixy,iz), ExtSrc(ixy,iz), N_U35(ixy,Iz), N_U38(ixy,Iz), N_Pu38(ixy,Iz), N_Pu39(ixy,Iz), N_Pu40(ixy,Iz), N_Pu42(ixy,Iz), N_Cm42(ixy,Iz), N_Cm44(ixy,Iz)
         enddo
      enddo


      write(unit_fname,*) "out6 :: 2-D Fission Source [#/sec]"
      allocate( FisSrc_2D(nxy) )
      FisSrc_2D = 0.d0
      do ixy = 1, nxy
         ix = ixytoix(ixy)
         iy = ixytoiy(ixy)
         do iz = 1, nz
           FisSrc_2D(ixy) = FisSrc_2D(ixy) + FisSrc(ixy,iz)
         enddo
         write(unit_fname,*) ix, iy, FisSrc_2D(ixy)
      enddo
      close(unit_fname)

      return
      end subroutine Write_ExtSrc


      !===============================================================================
      ! Prints error message and stops.
      !===============================================================================
      subroutine fatal_error(message)
      implicit none
      character(len=*), intent(in) :: message

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [fatal_error] in Mod_ExtSrc'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fatal_error] in Mod_ExtSrc'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      write(*,*) message
      stop
      end subroutine


      end module Mod_ExtSrc


#endif
