#ifdef siarhei_delete


      module Mod_LiteVer
      use Inc_XS_File
      use Inc_Geometry
      use Inc_3D
      use Inc_RP
      use Inc_Nuclide
      use Inc_Constant
      use mod_alloc

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none
      real(8), parameter :: tf_max=1500.0d0
      real(8), parameter :: tm_max=620.0d0
      real(8), parameter :: tm_min=280.0d0

      contains

      subroutine MacroXs_Table_alloc()
      implicit none


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [MacroXs_Table_alloc] in Mod_LiteVer'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.(allocated(macro_xs_table))) then
          allocate(macro_xs_table(Nxy, Nz))
          allocate(macroxsset(57))
      endif

      return
      end subroutine MacroXs_Table_alloc


      subroutine MacroXs_Table_dealloc()
      implicit none
      integer :: ix, iy, iz, ixy


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [MacroXs_Table_dealloc] in Mod_LiteVer'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (allocated(macro_xs_table)) then
         do ix=1,nx
            do iy=1,ny
               ixy=IxIyToIxy(ix,iy)
               if (ixy==0) then
                  cycle
               endif
               do iz=1,nz
                  deallocate(macro_xs_table(ixy,iz)%xs)
                  deallocate(macro_xs_table(ixy,iz)%hff)
               enddo
            enddo
         enddo
         deallocate(macro_xs_table)
         deallocate(macroXSset)
      endif

      return
      end subroutine MacroXs_Table_dealloc


      subroutine Micro_Xstable_2_Macro_Xstable()
      use mod_getnode, only: new_asym_itab
      implicit none
      integer :: ix, iy, iz, ixy, ixy_1n
      integer :: i_table, i_type, i_matrix, N_matrix
      integer :: pref1, pref2, pref
      integer :: pbrch1, pbrch2, pbrch
      real(8) :: burnup, wref1, wref2, wbrch1, wbrch2
      real(8) :: MacroXs_ref1(1:57) , MacroXs_ref2(1:57)
      real(8) :: MacroXs_brch1(1:57) , MacroXs_brch2(1:57)
      real(8) :: MacroXs_brref1(1:57) , MacroXs_brref2(1:57)
      real(8) :: n_density(1:42)
      real(8), allocatable :: bu_ref(:)
      real(8), allocatable :: bu_brch(:)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Micro_Xstable_2_Macro_Xstable] in Mod_LiteVer'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      do ix=1,nx
         do iy=1,ny
            ixy=IxIyToIxy(ix,iy)
            if (ixy==0) cycle
            do iz=1,nz
               burnup=bu(ixy,iz)
               i_table=I_Comp(ixy,iz)
               ixy_1n=i_4nto1n(ixy)
               i_table=new_asym_itab(axialcomp(i_lp_1n(ixy_1n),iz),ixy)
               i_type=Type_Tab(i_table)

               if (i_type==1) then
                  allocate(bu_ref(size(BU_Ref_NoBP(i_table,:))))
                  allocate(bu_brch(size(BU_Brch_NoBP(i_table,:))))
                  bu_ref(:)=BU_Ref_NoBP(i_table,:)
                  bu_brch(:)=BU_Brch_NoBP(i_table,:)
                  N_matrix=N_SP_FA
               elseif (i_type==2) then
                  allocate(bu_ref(size(BU_Ref_wBP(i_table,:))))
                  allocate(bu_brch(size(BU_Brch_wBP(i_table,:))))
                  bu_ref(:)=BU_Ref_wBP(i_table,:)
                  bu_brch(:)=BU_Brch_wBP(i_table,:)
                  N_matrix=N_SP_FA
               else
                  allocate(bu_ref(1))
                  allocate(bu_brch(1))
                  bu_ref(:)=0d0
                  bu_brch(:)=0d0
                  N_matrix = N_SP_RF
               endif

               allocate(macro_xs_table(ixy,iz)%xs(57,N_matrix))
               allocate(macro_xs_table(ixy,iz)%hff(XS_file_NbyN,XS_File_NbyN,N_matrix))
               macro_xs_table(ixy,iz)%xs(:,:)=0d0
               macro_xs_table(ixy,iz)%hff(:,:,:)=0d0

               do i_matrix=1,N_matrix ! number of state
                  if (i_type<=2) then
                     if (burnup>1d-10) then
                        pref1=1
                        pref2=size(bu_ref)
                        do while(.true.)
                           pref=(pref1+pref2)/2
                           if (bu_ref(pref)<burnup.and.burnup<=bu_ref(pref+1)) then
                              exit
                           elseif (burnup>bu_ref(pref+1)) then
                              pref1=pref
                           else
                              pref2=pref
                           endif
                        enddo
                        wref1=(bu_ref(pref+1)-burnup)/(bu_ref(pref+1)-bu_ref(pref))
                        wref2=(burnup-bu_ref(pref))/(bu_ref(pref+1)-bu_ref(pref))

                        pbrch1=1
                        pbrch2=size(bu_brch)
                        do while(.true.)
                           pbrch=(pbrch1+pbrch2)/2
                           if (bu_brch(pbrch)<burnup.and.burnup<=bu_brch(pbrch+1)) then
                              exit
                           elseif (burnup>bu_brch(pbrch+1)) then
                              pbrch1=pbrch
                           else
                              pbrch2=pbrch
                           endif
                        enddo
                        wbrch1=(bu_brch(pbrch+1)-burnup)/(bu_brch(pbrch+1)-bu_brch(pbrch))
                        wbrch2=(burnup-bu_brch(pbrch))/(bu_brch(pbrch+1)-bu_brch(pbrch))
                     else
                        pref=1
                        wref1=1d0
                        wref2=0d0
                        pbrch=1
                        wbrch1=1d0
                        wbrch2=0d0
                     endif

                     n_density( 1)=N_U34 (ixy,iz)
                     n_density( 2)=N_U35 (ixy,iz)
                     n_density( 3)=N_U36 (ixy,iz)
                     n_density( 4)=N_U37 (ixy,iz)
                     n_density( 5)=N_U38 (ixy,iz)
                     n_density( 6)=N_Np37(ixy,iz)
                     n_density( 7)=N_Np38(ixy,iz)
                     n_density( 8)=N_Np39(ixy,iz)
                     n_density( 9)=N_Pu38(ixy,iz)
                     n_density(10)=N_Pu39(ixy,iz)
                     n_density(11)=N_Pu40(ixy,iz)
                     n_density(12)=N_Pu41(ixy,iz)
                     n_density(13)=N_Pu42(ixy,iz)
                     n_density(14)=N_Pu43(ixy,iz)
                     n_density(15)=N_Am41(ixy,iz)
                     n_density(16)=N_As42(ixy,iz)
                     n_density(17)=N_Am42(ixy,iz)
                     n_density(18)=N_Am43(ixy,iz)
                     n_density(19)=N_Am44(ixy,iz)
                     n_density(20)=N_Cm42(ixy,iz)
                     n_density(21)=N_Cm43(ixy,iz)
                     n_density(22)=N_Cm44(ixy,iz)
                     n_density(23)=N_I35 (ixy,iz)
                     n_density(24)=N_Xe35(ixy,iz)
                     n_density(25)=N_Nd47(ixy,iz)
                     n_density(26)=N_Nd48(ixy,iz)
                     n_density(27)=N_Nd49(ixy,iz)
                     n_density(28)=N_Pm47(ixy,iz)
                     n_density(29)=N_Ps48(ixy,iz)
                     n_density(30)=N_Pm48(ixy,iz)
                     n_density(31)=N_Pm49(ixy,iz)
                     n_density(32)=N_Sm47(ixy,iz)
                     n_density(33)=N_Sm48(ixy,iz)
                     n_density(34)=N_Sm49(ixy,iz)
                     n_density(35)=N_Gd52(ixy,iz)
                     n_density(36)=N_Gd54(ixy,iz)
                     n_density(37)=N_Gd55(ixy,iz)
                     n_density(38)=N_Gd56(ixy,iz)
                     n_density(39)=N_Gd57(ixy,iz)
                     n_density(40)=N_Gd58(ixy,iz)
                     n_density(41)=N_Gd60(ixy,iz)
                     n_density(42)=N_B0  (ixy,iz)

                     call cal_MacXs_by_MicXs( Xsset_Table_Base(:,pref  ,i_table), MacroXs_ref1, n_density)
                     call cal_MacXs_by_MicXs( Xsset_Table_Base(:,pref+1,i_table), MacroXs_ref2, n_density)
                     call cal_MacXs_by_MicXs( Xsset_Table(:,i_matrix,pbrch  ,i_table), MacroXs_brch1, n_density)
                     call cal_MacXs_by_MicXs( Xsset_Table(:,i_matrix,pbrch+1,i_table), MacroXs_brch2, n_density)
                     call cal_MacXs_by_MicXs( Xsset_Table(:,1,pbrch  ,i_table), MacroXs_brref1, n_density)
                     call cal_MacXs_by_MicXs( Xsset_Table(:,1,pbrch+1,i_table), MacroXs_brref2, n_density)

                     macro_xs_table(ixy,iz)%xs(:,i_matrix)    = wref1  * MacroXs_ref1(:)   &
                                                            & + wref2  * MacroXs_ref2(:)   &
                                                            & - wbrch1 * MacroXs_brref1(:) &
                                                            & - wbrch2 * MacroXs_brref2(:) &
                                                            & + wbrch1 * MacroXs_brch1(:)  &
                                                            & + wbrch2 * MacroXs_brch2(:)
                     macro_xs_table(ixy,iz)%hff(:,:,i_matrix) = wref1  * HFF_Table_Base(i_table,pref  ,:,:)      &
                                                            & + wref2  * HFF_Table_Base(i_table,pref+1,:,:)      &
                                                            & - wbrch1 * HFF_Table(i_table,1,pref  ,:,:)         &
                                                            & - wbrch2 * HFF_Table(i_table,1,pref+1,:,:)         &
                                                            & + wbrch1 * HFF_Table(i_table,i_matrix,pbrch  ,:,:) &
                                                            & + wbrch2 * HFF_Table(i_table,i_matrix,pbrch+1,:,:)

                  else ! reflector
                     macro_xs_table(ixy,iz)%xs(1:12,i_matrix)  = Xsset_table(1:12,i_matrix,1,i_table)
                     macro_xs_table(ixy,iz)%xs(13,i_matrix)    = 0.0D+0
                     macro_xs_table(ixy,iz)%xs(14:57,i_matrix) = Xsset_table(401:444,i_matrix,1,i_table)
                     macro_xs_table(ixy,iz)%HFF(:,:,i_matrix)  = 0.0D+0
                  endif

               enddo
               deallocate(bu_ref)
               deallocate(bu_brch)
            enddo
         enddo
      enddo

      !js+deallocate(XSset_Table)
      !js+deallocate(XSset_Table_Base)
      !js+deallocate(HFF_Table)
      !js+deallocate(HFF_Table_Base)

      return
      end subroutine Micro_Xstable_2_Macro_Xstable


      subroutine cal_MacXs_by_MicXs(microXS,macroXS,n_density)
      implicit none
      real(8),intent(in)  :: microXS(1:444), n_density(42)
      real(8),intent(out) :: macroXS(1:57)
      integer             :: ig, ixs


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [cal_MacXs_by_MicXs] in Mod_LiteVer'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      macroXS=0d0
      ixs=13

      ! HM nuclide
      ! transport cross section
      do ig=1,2
         macroXS(ig)    = macroXS(ig)    + n_density(1)  * microXS(ixs);  ixs= ixs+1   !  U-234
         macroXS(ig)    = macroXS(ig)    + n_density(2)  * microXS(ixs);  ixs= ixs+1   !  U-235
         macroXS(ig)    = macroXS(ig)    + n_density(3)  * microXS(ixs);  ixs= ixs+1   !  U-236
         macroXS(ig)    = macroXS(ig)    + n_density(4)  * microXS(ixs);  ixs= ixs+1   !  U-237
         macroXS(ig)    = macroXS(ig)    + n_density(5)  * microXS(ixs);  ixs= ixs+1   !  U-238
         macroXS(ig)    = macroXS(ig)    + n_density(6)  * microXS(ixs);  ixs= ixs+1   ! Np-237
         macroXS(ig)    = macroXS(ig)    + n_density(7)  * microXS(ixs);  ixs= ixs+1   ! Np-238
         macroXS(ig)    = macroXS(ig)    + n_density(8)  * microXS(ixs);  ixs= ixs+1   ! Np-239
         macroXS(ig)    = macroXS(ig)    + n_density(9)  * microXS(ixs);  ixs= ixs+1   ! Pu-238
         macroXS(ig)    = macroXS(ig)    + n_density(10) * microXS(ixs);  ixs= ixs+1   ! Pu-239
         macroXS(ig)    = macroXS(ig)    + n_density(11) * microXS(ixs);  ixs= ixs+1   ! Pu-240
         macroXS(ig)    = macroXS(ig)    + n_density(12) * microXS(ixs);  ixs= ixs+1   ! Pu-241
         macroXS(ig)    = macroXS(ig)    + n_density(13) * microXS(ixs);  ixs= ixs+1   ! Pu-242
         macroXS(ig)    = macroXS(ig)    + n_density(14) * microXS(ixs);  ixs= ixs+1   ! Pu-243
         macroXS(ig)    = macroXS(ig)    + n_density(15) * microXS(ixs);  ixs= ixs+1   ! Am-241
         macroXS(ig)    = macroXS(ig)    + n_density(16) * microXS(ixs);  ixs= ixs+1   ! Am-242
         macroXS(ig)    = macroXS(ig)    + n_density(17) * microXS(ixs);  ixs= ixs+1   ! Am-242m
         macroXS(ig)    = macroXS(ig)    + n_density(18) * microXS(ixs);  ixs= ixs+1   ! Am-243
         macroXS(ig)    = macroXS(ig)    + n_density(19) * microXS(ixs);  ixs= ixs+1   ! Am-244
         macroXS(ig)    = macroXS(ig)    + n_density(20) * microXS(ixs);  ixs= ixs+1   ! Cm-242
         macroXS(ig)    = macroXS(ig)    + n_density(21) * microXS(ixs);  ixs= ixs+1   ! Cm-243
         macroXS(ig)    = macroXS(ig)    + n_density(22) * microXS(ixs);  ixs= ixs+1   ! Cm-244
      enddo
      ! absorption cross section
      do ig=1,2
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(1)  * microXS(ixs);  ixs= ixs+1   !  U-234
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(2)  * microXS(ixs);  ixs= ixs+1   !  U-235
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(3)  * microXS(ixs);  ixs= ixs+1   !  U-236
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(4)  * microXS(ixs);  ixs= ixs+1   !  U-237
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(5)  * microXS(ixs);  ixs= ixs+1   !  U-238
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(6)  * microXS(ixs);  ixs= ixs+1   ! Np-237
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(7)  * microXS(ixs);  ixs= ixs+1   ! Np-238
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(8)  * microXS(ixs);  ixs= ixs+1   ! Np-239
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(9)  * microXS(ixs);  ixs= ixs+1   ! Pu-238
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(10) * microXS(ixs);  ixs= ixs+1   ! Pu-239
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(11) * microXS(ixs);  ixs= ixs+1   ! Pu-240
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(12) * microXS(ixs);  ixs= ixs+1   ! Pu-241
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(13) * microXS(ixs);  ixs= ixs+1   ! Pu-242
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(14) * microXS(ixs);  ixs= ixs+1   ! Pu-243
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(15) * microXS(ixs);  ixs= ixs+1   ! Am-241
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(16) * microXS(ixs);  ixs= ixs+1   ! Am-242
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(17) * microXS(ixs);  ixs= ixs+1   ! Am-242m
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(18) * microXS(ixs);  ixs= ixs+1   ! Am-243
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(19) * microXS(ixs);  ixs= ixs+1   ! Am-244
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(20) * microXS(ixs);  ixs= ixs+1   ! Cm-242
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(21) * microXS(ixs);  ixs= ixs+1   ! Cm-243
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(22) * microXS(ixs);  ixs= ixs+1   ! Cm-244
      enddo
      ! fission cross section
      do ig=1,2
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(1)  * microXS(ixs);  ixs= ixs+1   !  U-234
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(2)  * microXS(ixs);  ixs= ixs+1   !  U-235
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(3)  * microXS(ixs);  ixs= ixs+1   !  U-236
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(4)  * microXS(ixs);  ixs= ixs+1   !  U-237
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(5)  * microXS(ixs);  ixs= ixs+1   !  U-238
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(6)  * microXS(ixs);  ixs= ixs+1   ! Np-237
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(7)  * microXS(ixs);  ixs= ixs+1   ! Np-238
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(8)  * microXS(ixs);  ixs= ixs+1   ! Np-239
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(9)  * microXS(ixs);  ixs= ixs+1   ! Pu-238
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(10) * microXS(ixs);  ixs= ixs+1   ! Pu-239
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(11) * microXS(ixs);  ixs= ixs+1   ! Pu-240
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(12) * microXS(ixs);  ixs= ixs+1   ! Pu-241
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(13) * microXS(ixs);  ixs= ixs+1   ! Pu-242
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(14) * microXS(ixs);  ixs= ixs+1   ! Pu-243
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(15) * microXS(ixs);  ixs= ixs+1   ! Am-241
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(16) * microXS(ixs);  ixs= ixs+1   ! Am-242
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(17) * microXS(ixs);  ixs= ixs+1   ! Am-242m
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(18) * microXS(ixs);  ixs= ixs+1   ! Am-243
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(19) * microXS(ixs);  ixs= ixs+1   ! Am-244
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(20) * microXS(ixs);  ixs= ixs+1   ! Cm-242
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(21) * microXS(ixs);  ixs= ixs+1   ! Cm-243
         macroXS(ig+4)  = macroXS(ig+4)  + n_density(22) * microXS(ixs);  ixs= ixs+1   ! Cm-244
      enddo
      ! nu-fission cross section
      do ig=1,2
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(1)  * microXS(ixs);  ixs= ixs+1   !  U-234
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(2)  * microXS(ixs);  ixs= ixs+1   !  U-235
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(3)  * microXS(ixs);  ixs= ixs+1   !  U-236
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(4)  * microXS(ixs);  ixs= ixs+1   !  U-237
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(5)  * microXS(ixs);  ixs= ixs+1   !  U-238
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(6)  * microXS(ixs);  ixs= ixs+1   ! Np-237
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(7)  * microXS(ixs);  ixs= ixs+1   ! Np-238
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(8)  * microXS(ixs);  ixs= ixs+1   ! Np-239
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(9)  * microXS(ixs);  ixs= ixs+1   ! Pu-238
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(10) * microXS(ixs);  ixs= ixs+1   ! Pu-239
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(11) * microXS(ixs);  ixs= ixs+1   ! Pu-240
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(12) * microXS(ixs);  ixs= ixs+1   ! Pu-241
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(13) * microXS(ixs);  ixs= ixs+1   ! Pu-242
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(14) * microXS(ixs);  ixs= ixs+1   ! Pu-243
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(15) * microXS(ixs);  ixs= ixs+1   ! Am-241
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(16) * microXS(ixs);  ixs= ixs+1   ! Am-242
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(17) * microXS(ixs);  ixs= ixs+1   ! Am-242m
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(18) * microXS(ixs);  ixs= ixs+1   ! Am-243
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(19) * microXS(ixs);  ixs= ixs+1   ! Am-244
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(20) * microXS(ixs);  ixs= ixs+1   ! Cm-242
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(21) * microXS(ixs);  ixs= ixs+1   ! Cm-243
         macroXS(ig+6)  = macroXS(ig+6)  + n_density(22) * microXS(ixs);  ixs= ixs+1   ! Cm-244
      enddo
      ! kappa-fission cross section
      do ig=1,2
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(1)  * microXS(ixs);  ixs= ixs+1   !  U-234
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(2)  * microXS(ixs);  ixs= ixs+1   !  U-235
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(3)  * microXS(ixs);  ixs= ixs+1   !  U-236
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(4)  * microXS(ixs);  ixs= ixs+1   !  U-237
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(5)  * microXS(ixs);  ixs= ixs+1   !  U-238
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(6)  * microXS(ixs);  ixs= ixs+1   ! Np-237
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(7)  * microXS(ixs);  ixs= ixs+1   ! Np-238
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(8)  * microXS(ixs);  ixs= ixs+1   ! Np-239
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(9)  * microXS(ixs);  ixs= ixs+1   ! Pu-238
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(10) * microXS(ixs);  ixs= ixs+1   ! Pu-239
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(11) * microXS(ixs);  ixs= ixs+1   ! Pu-240
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(12) * microXS(ixs);  ixs= ixs+1   ! Pu-241
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(13) * microXS(ixs);  ixs= ixs+1   ! Pu-242
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(14) * microXS(ixs);  ixs= ixs+1   ! Pu-243
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(15) * microXS(ixs);  ixs= ixs+1   ! Am-241
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(16) * microXS(ixs);  ixs= ixs+1   ! Am-242
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(17) * microXS(ixs);  ixs= ixs+1   ! Am-242m
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(18) * microXS(ixs);  ixs= ixs+1   ! Am-243
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(19) * microXS(ixs);  ixs= ixs+1   ! Am-244
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(20) * microXS(ixs);  ixs= ixs+1   ! Cm-242
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(21) * microXS(ixs);  ixs= ixs+1   ! Cm-243
         macroXS(ig+8)  = macroXS(ig+8)  + n_density(22) * microXS(ixs);  ixs= ixs+1   ! Cm-244
      enddo
      ! scattering cross section
      do ig=1,2
         macroXS(ig+10) = macroXS(ig+10) + n_density(1)  * microXS(ixs);  ixs= ixs+1   !  U-234
         macroXS(ig+10) = macroXS(ig+10) + n_density(2)  * microXS(ixs);  ixs= ixs+1   !  U-235
         macroXS(ig+10) = macroXS(ig+10) + n_density(3)  * microXS(ixs);  ixs= ixs+1   !  U-236
         macroXS(ig+10) = macroXS(ig+10) + n_density(4)  * microXS(ixs);  ixs= ixs+1   !  U-237
         macroXS(ig+10) = macroXS(ig+10) + n_density(5)  * microXS(ixs);  ixs= ixs+1   !  U-238
         macroXS(ig+10) = macroXS(ig+10) + n_density(6)  * microXS(ixs);  ixs= ixs+1   ! Np-237
         macroXS(ig+10) = macroXS(ig+10) + n_density(7)  * microXS(ixs);  ixs= ixs+1   ! Np-238
         macroXS(ig+10) = macroXS(ig+10) + n_density(8)  * microXS(ixs);  ixs= ixs+1   ! Np-239
         macroXS(ig+10) = macroXS(ig+10) + n_density(9)  * microXS(ixs);  ixs= ixs+1   ! Pu-238
         macroXS(ig+10) = macroXS(ig+10) + n_density(10) * microXS(ixs);  ixs= ixs+1   ! Pu-239
         macroXS(ig+10) = macroXS(ig+10) + n_density(11) * microXS(ixs);  ixs= ixs+1   ! Pu-240
         macroXS(ig+10) = macroXS(ig+10) + n_density(12) * microXS(ixs);  ixs= ixs+1   ! Pu-241
         macroXS(ig+10) = macroXS(ig+10) + n_density(13) * microXS(ixs);  ixs= ixs+1   ! Pu-242
         macroXS(ig+10) = macroXS(ig+10) + n_density(14) * microXS(ixs);  ixs= ixs+1   ! Pu-243
         macroXS(ig+10) = macroXS(ig+10) + n_density(15) * microXS(ixs);  ixs= ixs+1   ! Am-241
         macroXS(ig+10) = macroXS(ig+10) + n_density(16) * microXS(ixs);  ixs= ixs+1   ! Am-242
         macroXS(ig+10) = macroXS(ig+10) + n_density(17) * microXS(ixs);  ixs= ixs+1   ! Am-242m
         macroXS(ig+10) = macroXS(ig+10) + n_density(18) * microXS(ixs);  ixs= ixs+1   ! Am-243
         macroXS(ig+10) = macroXS(ig+10) + n_density(19) * microXS(ixs);  ixs= ixs+1   ! Am-244
         macroXS(ig+10) = macroXS(ig+10) + n_density(20) * microXS(ixs);  ixs= ixs+1   ! Cm-242
         macroXS(ig+10) = macroXS(ig+10) + n_density(21) * microXS(ixs);  ixs= ixs+1   ! Cm-243
         macroXS(ig+10) = macroXS(ig+10) + n_density(22) * microXS(ixs);  ixs= ixs+1   ! Cm-244
      enddo

      ! Fission Production nuclide
      ! transport cross section
      do ig=1,2
         macroXS(ig)    = macroXS(ig)    + n_density(23) * microXS(ixs);  ixs= ixs+1   !  I-135
         macroXS(ig)    = macroXS(ig)    + n_density(24) * microXS(ixs);  ixs= ixs+1   ! Xe-135
         macroXS(ig)    = macroXS(ig)    + n_density(25) * microXS(ixs);  ixs= ixs+1   ! Nd-147
         macroXS(ig)    = macroXS(ig)    + n_density(26) * microXS(ixs);  ixs= ixs+1   ! Nd-148
         macroXS(ig)    = macroXS(ig)    + n_density(27) * microXS(ixs);  ixs= ixs+1   ! Nd-149
         macroXS(ig)    = macroXS(ig)    + n_density(28) * microXS(ixs);  ixs= ixs+1   ! Pm-147
         macroXS(ig)    = macroXS(ig)    + n_density(29) * microXS(ixs);  ixs= ixs+1   ! Pm-148m
         macroXS(ig)    = macroXS(ig)    + n_density(30) * microXS(ixs);  ixs= ixs+1   ! Pm-148
         macroXS(ig)    = macroXS(ig)    + n_density(31) * microXS(ixs);  ixs= ixs+1   ! Pm-149
         macroXS(ig)    = macroXS(ig)    + n_density(32) * microXS(ixs);  ixs= ixs+1   ! Sm-147
         macroXS(ig)    = macroXS(ig)    + n_density(33) * microXS(ixs);  ixs= ixs+1   ! Sm-148
         macroXS(ig)    = macroXS(ig)    + n_density(34) * microXS(ixs);  ixs= ixs+1   ! Sm-149
      enddo
      ! absorption cross section
      do ig=1,2
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(23) * microXS(ixs);  ixs= ixs+1   !  I-135
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(24) * microXS(ixs);  ixs= ixs+1   ! Xe-135
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(25) * microXS(ixs);  ixs= ixs+1   ! Nd-147
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(26) * microXS(ixs);  ixs= ixs+1   ! Nd-148
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(27) * microXS(ixs);  ixs= ixs+1   ! Nd-149
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(28) * microXS(ixs);  ixs= ixs+1   ! Pm-147
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(29) * microXS(ixs);  ixs= ixs+1   ! Pm-148m
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(30) * microXS(ixs);  ixs= ixs+1   ! Pm-148
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(31) * microXS(ixs);  ixs= ixs+1   ! Pm-149
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(32) * microXS(ixs);  ixs= ixs+1   ! Sm-147
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(33) * microXS(ixs);  ixs= ixs+1   ! Sm-148
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(34) * microXS(ixs);  ixs= ixs+1   ! Sm-149
      enddo
      ! scattering cross section
      do ig=1,2
         macroXS(ig+10) = macroXS(ig+10) + n_density(23) * microXS(ixs);  ixs= ixs+1   !  I-135
         macroXS(ig+10) = macroXS(ig+10) + n_density(24) * microXS(ixs);  ixs= ixs+1   ! Xe-135
         macroXS(ig+10) = macroXS(ig+10) + n_density(25) * microXS(ixs);  ixs= ixs+1   ! Nd-147
         macroXS(ig+10) = macroXS(ig+10) + n_density(26) * microXS(ixs);  ixs= ixs+1   ! Nd-148
         macroXS(ig+10) = macroXS(ig+10) + n_density(27) * microXS(ixs);  ixs= ixs+1   ! Nd-149
         macroXS(ig+10) = macroXS(ig+10) + n_density(28) * microXS(ixs);  ixs= ixs+1   ! Pm-147
         macroXS(ig+10) = macroXS(ig+10) + n_density(29) * microXS(ixs);  ixs= ixs+1   ! Pm-148m
         macroXS(ig+10) = macroXS(ig+10) + n_density(30) * microXS(ixs);  ixs= ixs+1   ! Pm-148
         macroXS(ig+10) = macroXS(ig+10) + n_density(31) * microXS(ixs);  ixs= ixs+1   ! Pm-149
         macroXS(ig+10) = macroXS(ig+10) + n_density(32) * microXS(ixs);  ixs= ixs+1   ! Sm-147
         macroXS(ig+10) = macroXS(ig+10) + n_density(33) * microXS(ixs);  ixs= ixs+1   ! Sm-148
         macroXS(ig+10) = macroXS(ig+10) + n_density(34) * microXS(ixs);  ixs= ixs+1   ! Sm-149
      enddo

      ! BP nuclide
      ! transport cross section
      do ig=1,2
         macroXS(ig)    = macroXS(ig)    + n_density(35) * microXS(ixs);  ixs= ixs+1   ! Gd-152
         macroXS(ig)    = macroXS(ig)    + n_density(36) * microXS(ixs);  ixs= ixs+1   ! Gd-154
         macroXS(ig)    = macroXS(ig)    + n_density(37) * microXS(ixs);  ixs= ixs+1   ! Gd-155
         macroXS(ig)    = macroXS(ig)    + n_density(38) * microXS(ixs);  ixs= ixs+1   ! Gd-156
         macroXS(ig)    = macroXS(ig)    + n_density(39) * microXS(ixs);  ixs= ixs+1   ! Gd-157
         macroXS(ig)    = macroXS(ig)    + n_density(40) * microXS(ixs);  ixs= ixs+1   ! Gd-158
         macroXS(ig)    = macroXS(ig)    + n_density(41) * microXS(ixs);  ixs= ixs+1   ! Gd-160
      enddo
      ! absorption cross section
      do ig=1,2
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(35) * microXS(ixs);  ixs= ixs+1   ! Gd-152
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(36) * microXS(ixs);  ixs= ixs+1   ! Gd-154
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(37) * microXS(ixs);  ixs= ixs+1   ! Gd-155
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(38) * microXS(ixs);  ixs= ixs+1   ! Gd-156
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(39) * microXS(ixs);  ixs= ixs+1   ! Gd-157
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(40) * microXS(ixs);  ixs= ixs+1   ! Gd-158
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(41) * microXS(ixs);  ixs= ixs+1   ! Gd-160
      enddo
      ! scattering cross section
      do ig=1,2
         macroXS(ig+10) = macroXS(ig+10) + n_density(35) * microXS(ixs);  ixs= ixs+1   ! Gd-152
         macroXS(ig+10) = macroXS(ig+10) + n_density(36) * microXS(ixs);  ixs= ixs+1   ! Gd-154
         macroXS(ig+10) = macroXS(ig+10) + n_density(37) * microXS(ixs);  ixs= ixs+1   ! Gd-155
         macroXS(ig+10) = macroXS(ig+10) + n_density(38) * microXS(ixs);  ixs= ixs+1   ! Gd-156
         macroXS(ig+10) = macroXS(ig+10) + n_density(39) * microXS(ixs);  ixs= ixs+1   ! Gd-157
         macroXS(ig+10) = macroXS(ig+10) + n_density(40) * microXS(ixs);  ixs= ixs+1   ! Gd-158
         macroXS(ig+10) = macroXS(ig+10) + n_density(41) * microXS(ixs);  ixs= ixs+1   ! Gd-160
      enddo

      ! B10 nuclide
      ! transport cross section
      do ig=1,2
         macroXS(ig)    = macroXS(ig)    + n_density(42) * microXS(ixs);  ixs= ixs+1   !  B-10
      enddo
      ! absorption cross section
      do ig=1,2
         macroXS(ig+2)  = macroXS(ig+2)  + n_density(42) * microXS(ixs);  ixs= ixs+1   !  B-10
      enddo
      ! scattering cross section
      do ig=1,2
         macroXS(ig+10) = macroXS(ig+10) + n_density(42) * microXS(ixs);  ixs= ixs+1   !  B-10
      enddo

      macroXS(13) = macroXS(13) + n_density(2)  * microXS(ixs);  ixs= ixs+1   !  U-235, (n,2n), g1
      macroXS(13) = macroXS(13) + n_density(5)  * microXS(ixs);  ixs= ixs+1   !  U-238, (n,2n), g1
      macroXS(13) = macroXS(13) + n_density(6)  * microXS(ixs);  ixs= ixs+1   ! Np-237, (n,2n), g1
      macroXS(13) = macroXS(13) + n_density(10) * microXS(ixs);  ixs= ixs+1   ! Pu-239, (n,2n), g1

      macroXS(1:13)=macroXS(1:13)*barn
      ixs=1
      macroXS(1)  = macroXS(1)  + microXS(ixs); ixs = ixs + 1   ! residual transport     g1
      macroXS(2)  = macroXS(2)  + microXS(ixs); ixs = ixs + 1   ! residual transport     g2
      macroXS(3)  = macroXS(3)  + microXS(ixs); ixs = ixs + 1   ! residual absorption    g1
      macroXS(4)  = macroXS(4)  + microXS(ixs); ixs = ixs + 1   ! residual absorption    g2
      macroXS(5)  = macroXS(5)  + microXS(ixs); ixs = ixs + 1   ! residual fission       g1
      macroXS(6)  = macroXS(6)  + microXS(ixs); ixs = ixs + 1   ! residual fission       g2
      macroXS(7)  = macroXS(7)  + microXS(ixs); ixs = ixs + 1   ! residual nu-fission    g1
      macroXS(8)  = macroXS(8)  + microXS(ixs); ixs = ixs + 1   ! residual nu-fission    g2
      macroXS(9)  = macroXS(9)  + microXS(ixs); ixs = ixs + 1   ! residual kappa-fission g1
      macroXS(10) = macroXS(10) + microXS(ixs); ixs = ixs + 1   ! residual kappa-fission g2
      macroXS(11) = macroXS(11) + microXS(ixs); ixs = ixs + 1   ! residual scattering    g1
      macroXS(12) = macroXS(12) + microXS(ixs); ixs = ixs + 1   ! residual scattering    g2

      macroXS(14:57)=microXS(401:444) ! FY, ADF, CDF, KP

      return
      end subroutine cal_MacXs_by_MicXs


      subroutine maXS2solver(ixy,iz)
      use Inc_maXS
      use Inc_miXS
      use Inc_Option, only: N_Group
      use Inc_Detector
      use INC_DF
      use INC_Kinetics
      implicit none
      integer, intent(in) :: ixy, iz
      integer :: ig, ixs


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [maXS2solver] in Mod_LiteVer'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      maXS_tr_3D    (ixy,iz,1) = MacroXSset(1)
      maXS_tr_3D    (ixy,iz,2) = MacroXSset(2)
      D_3D          (ixy,iz,1) = 1d0 / max(3.0d0*maXS_tr_3D(ixy,iz,1),1d-30)
      D_3D          (ixy,iz,2) = 1d0 / max(3.0d0*maXS_tr_3D(ixy,iz,2),1d-30)
      maXS_a_3D     (ixy,iz,1) = MacroXSset(3)
      maXS_a_3D     (ixy,iz,2) = MacroXSset(4)
      maXS_f_3D     (ixy,iz,1) = MacroXSset(5)
      maXS_f_3D     (ixy,iz,2) = MacroXSset(6)
      nu_maXS_f_3D  (ixy,iz,1) = MacroXSset(7)
      nu_maXS_f_3D  (ixy,iz,2) = MacroXSset(8)
      kap_maXS_f_3D (ixy,iz,1) = MacroXSset(9)
      kap_maXS_f_3D (ixy,iz,2) = MacroXSset(10)
      maXS_s_3D     (ixy,iz,1) = MacroXSset(11)
      maXS_s_3D     (ixy,iz,2) = MacroXSset(12)
      maXS_r_3D     (ixy,iz,1) = maXS_a_3D(ixy,iz,1) + maXS_s_3D(ixy,iz,1)
      maXS_r_3D     (ixy,iz,2) = maXS_a_3D(ixy,iz,2) + maXS_s_3D(ixy,iz,2)
      maXS_s_3D     (ixy,iz,2) = 0d0
      do ig=1,N_Group
         maXS_r_3D(ixy,iz,ig) = maXS_a_3D(ixy,iz,ig) + maXS_s_3D(ixy,iz,ig)
      enddo
      maXS_scat_3D(1,2,ixy,iz) = maXS_s_3D(ixy,iz,1)

      ixs=14
      Y_I35     (ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1
      Y_Xe35    (ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1
      Y_Nd47    (ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1
      Y_Nd48    (ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1
      Y_Nd49    (ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1
      Y_Pm47    (ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1
      Y_Ps48    (ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1
      Y_Pm48    (ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1
      Y_Pm49    (ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1
      Y_Sm49    (ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1
      Y_Xe35_eff(ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1
      Y_Pm49_eff(ixy,iz) = MacroXSset(ixs);  ixs = ixs + 1

      ixs=26
      do ig=1,N_Group
         ADF_Lx  (ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
         ADF_Rx  (ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
         ADF_Ly  (ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
         ADF_Ry  (ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
         ADF_Avg (ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
         CDF_LxLy(ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
         CDF_LxRy(ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
         CDF_RxLy(ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
         CDF_RxRy(ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
      enddo

      ixs=44
      do ig=1,N_Group
         v_Inv(ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
      enddo
      do ig=1,N_Group_d
         beta_d_eff(ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
         lambda_d  (ixy,iz,ig) = MacroXSset(ixs);  ixs=ixs+1
      enddo

      return
      end subroutine maXS2solver


      subroutine macroxsfb
      use Inc_INP
      use Mod_GetSome, only: Get_Reg_VF, Get_ADF
      use Inc_TH
      use Inc_maXS
      use inc_flag, only: flag_tr_macroxs_1d
      implicit none
      integer :: ixy, iz


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [macroxsfb] in Mod_LiteVer'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call WATCH('maXS Feedback','START')
      do ixy=1,nxy
         do iz=1,nz
            call Get_Reg_VF(ixy,iz)
            if (flag_tr_macroxs_1d) then
               call MacroXSsetFB_1D(ixy,iz)
            else
               call MacroXSsetFB(ixy,iz)
            endif
            call maXS2solver(ixy,iz)
         enddo
      enddo
      call Get_ADF
      call WATCH('maXS Feedback','END')

      return
      end subroutine macroxsfb


      subroutine macroxssetfb(Ixy, Iz)
      use Inc_3D, only: T_Fuel, T_Mod
      use Inc_TH, only: PPM, TM_In
      use inc_th, only: flag_th_chanwise, i_chan_1n, chanwise_t_inlet
      use Inc_CR, only: Reg_VF
      use Inc_RP, only: I_LP_1N, AxialComp
      use Inc_Detector
      use Inc_INP, only: flag_force_tm, force_tm, force_dtm
      use Inc_INP, only: flag_force_tf, force_tf, force_dtf
      use Inc_INP, only: flag_force_tdop, force_tdop, force_dtdop
      use Inc_INP, only: flag_force_ftc, force_ftc, flag_force_mtc, force_mtc
      use Mod_GetNode, only: new_asym_itab
      implicit none
      integer, intent(in) :: Ixy, Iz
      integer :: Ixy_1N
      integer(4) :: i
      real(8) :: tm_in_bak, chg_Tm, chg_tf, tmp_tdop
      real(8) :: tf_bnd, tm_bnd


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [macroxssetfb] in Mod_LiteVer'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (flag_th_chanwise) then
         tm_in_bak=chanwise_t_inlet(i_chan_1n(i_4nto1n(ixy)))+degtok
      else
         tm_in_bak=Tm_in+DegToK
      endif

      if (flag_force_tm) then
         chg_Tm=tm_in_bak+(T_Mod(Ixy,Iz)-tm_in_bak)*force_tm+force_dtm
      else
         chg_Tm=T_Mod(Ixy,Iz)
      endif
      if (flag_force_tf) then
         chg_Tf=t_fuel(ixy,iz)*force_tf+force_dtf
      elseif (flag_force_tdop) then
         tmp_tdop=sqrt(t_fuel(ixy,iz))
         tmp_tdop=tmp_tdop*force_tdop+force_dtdop
         chg_Tf=tmp_tdop*tmp_tdop
      else
         chg_Tf=t_fuel(ixy,iz)
      endif
      if (flag_force_ftc) chg_tf=chg_tf+force_ftc
      if (flag_force_mtc) chg_tm=chg_tm+force_mtc

      tf_bnd=min(tf_max,chg_tf)
      tm_bnd=max(tm_min,min(tm_max,chg_Tm))

      W_SP(:,Ixy,Iz) = 0d0
      W_SP_1 = 0d0

      Ixy_1N = I_4Nto1N(Ixy)
      I_Tab = new_asym_itab(AxialComp(I_LP_1N(Ixy_1N),Iz),Ixy)
      I_Type = Type_Tab(I_Tab)

      if_state=0
      if (if_branch==4) WC_TF=0d0
      if (if_branch==4) WC_TM=0d0
      W_PPM=0d0
      W_TM=0d0
      W_TF=0d0

      ! calculate the weight factors for PPM
      W_PPM(1) = (PPM       -Var_PPM(2))*(PPM       -Var_PPM(3))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(1)-Var_PPM(2))*(Var_PPM(1)-Var_PPM(3))*(Var_PPM(1)-Var_PPM(4)) )
      W_PPM(2) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(3))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(2)-Var_PPM(1))*(Var_PPM(2)-Var_PPM(3))*(Var_PPM(2)-Var_PPM(4)) )
      W_PPM(3) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(2))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(3)-Var_PPM(1))*(Var_PPM(3)-Var_PPM(2))*(Var_PPM(3)-Var_PPM(4)) )
      W_PPM(4) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(2))*(PPM       -Var_PPM(3)) / &
               ( (Var_PPM(4)-Var_PPM(1))*(Var_PPM(4)-Var_PPM(2))*(Var_PPM(4)-Var_PPM(3)) )
      if ((i_type<3) .and. (if_branch==4)) then
         if ((tm_bnd>=Var_ref_TM(1)).and.(tf_bnd>=Var_ref_TF(1))) then
            if_state=1
            W_TF(1) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(2)))*(sqrt(tf_bnd)       -sqrt(Var_ref_TF(3))) / &
                    ( (sqrt(Var_ref_TF(1))-sqrt(Var_ref_TF(2)))*(sqrt(Var_ref_TF(1))-sqrt(Var_ref_TF(3))) )
            W_TF(2) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(1)))*(sqrt(tf_bnd)       -sqrt(Var_ref_TF(3))) / &
                    ( (sqrt(Var_ref_TF(2))-sqrt(Var_ref_TF(1)))*(sqrt(Var_ref_TF(2))-sqrt(Var_ref_TF(3))) )
            W_TF(3) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(1)))*(sqrt(tf_bnd)       -sqrt(Var_ref_TF(2))) / &
                    ( (sqrt(Var_ref_TF(3))-sqrt(Var_ref_TF(1)))*(sqrt(Var_ref_TF(3))-sqrt(Var_ref_TF(2))) )
            W_TF(4) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(3))) / &
                    ( (sqrt(Var_ref_TF(2))-sqrt(Var_ref_TF(3))) )
            W_TF(5) = 1d0 - W_TF(4)
            ! calculate the weight factors for TM
            W_TM(1) = (tm_bnd       -Var_ref_TM(2))*(tm_bnd       -Var_ref_TM(3)) / &
                    ( (Var_ref_TM(1)-Var_ref_TM(2))*(Var_ref_TM(1)-Var_ref_TM(3)) )
            W_TM(2) = (tm_bnd       -Var_ref_TM(1))*(tm_bnd       -Var_ref_TM(3)) / &
                    ( (Var_ref_TM(2)-Var_ref_TM(1))*(Var_ref_TM(2)-Var_ref_TM(3)) )
            W_TM(3) = (tm_bnd       -Var_ref_TM(1))*(tm_bnd       -Var_ref_TM(2)) / &
                    ( (Var_ref_TM(3)-Var_ref_TM(1))*(Var_ref_TM(3)-Var_ref_TM(2)) )
            W_TM(4) = (tm_bnd       -Var_ref_TM(3)) / &
                    ( (Var_ref_TM(1)-Var_ref_TM(3)) )
            W_TM(5) = 1d0 - W_TM(4)
         elseif (((tm_bnd<=Var_TM(9)).and.(tf_bnd<=Var_TF(2))).and.((tm_bnd>=Var_TM(1)-5d0).and.(tf_bnd>=Var_TF(1)-5d0))) then
            if(tm_bnd<=400d0) then
               if_state=2
               WC_TF(1)= (sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                       ( (sqrt(Var_TF(1))-sqrt(Var_TF(2))) )
               WC_TF(2)= (sqrt(tf_bnd)   -sqrt(Var_TF(1))) / &
                       ( (sqrt(Var_TF(2))-sqrt(Var_TF(1))) )
               WC_TM(1) = (tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(3))*(tm_bnd   -Var_TM(4)) / &
                        ( (Var_TM(1)-Var_TM(2))*(Var_TM(1)-Var_TM(3))*(Var_TM(1)-Var_TM(4)) )
               WC_TM(2) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(3))*(tm_bnd   -Var_TM(4)) / &
                        ( (Var_TM(2)-Var_TM(1))*(Var_TM(2)-Var_TM(3))*(Var_TM(2)-Var_TM(4)) )
               WC_TM(3) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(4)) / &
                        ( (Var_TM(3)-Var_TM(1))*(Var_TM(3)-Var_TM(2))*(Var_TM(3)-Var_TM(4)) )
               WC_TM(4) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(3)) / &
                        ( (Var_TM(4)-Var_TM(1))*(Var_TM(4)-Var_TM(2))*(Var_TM(4)-Var_TM(3)) )
           elseif(tm_bnd>400d0 .and. tm_bnd<=500d0) then
               if_state=3
               WC_TF(1)= (sqrt(tf_bnd)   -sqrt(Var_TF(2)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(1))-sqrt(Var_TF(2)))*(sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
               WC_TF(2)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(2))-sqrt(Var_TF(1)))*(sqrt(Var_TF(2))-sqrt(Var_TF(3))) )
               WC_TF(3)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                       ( (sqrt(Var_TF(3))-sqrt(Var_TF(1)))*(sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
               WC_TM(1) = (tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(6))*(tm_bnd   -Var_TM(7)) / &
                        ( (Var_TM(4)-Var_TM(5))*(Var_TM(4)-Var_TM(6))*(Var_TM(4)-Var_TM(7)) )
               WC_TM(2) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(6))*(tm_bnd   -Var_TM(7)) / &
                        ( (Var_TM(5)-Var_TM(4))*(Var_TM(5)-Var_TM(6))*(Var_TM(5)-Var_TM(7)) )
               WC_TM(3) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(7)) / &
                        ( (Var_TM(6)-Var_TM(4))*(Var_TM(6)-Var_TM(5))*(Var_TM(6)-Var_TM(7)) )
               WC_TM(4) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(6)) / &
                        ( (Var_TM(7)-Var_TM(4))*(Var_TM(7)-Var_TM(5))*(Var_TM(7)-Var_TM(6)) )
           elseif(tm_bnd>500d0) then
              if_state=4
              WC_TF(1)= (sqrt(tf_bnd)   -sqrt(Var_TF(2)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                      ( (sqrt(Var_TF(1))-sqrt(Var_TF(2)))*(sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
              WC_TF(2)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                      ( (sqrt(Var_TF(2))-sqrt(Var_TF(1)))*(sqrt(Var_TF(2))-sqrt(Var_TF(3))) )
              WC_TF(3)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                      ( (sqrt(Var_TF(3))-sqrt(Var_TF(1)))*(sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
              WC_TM(1) = (tm_bnd   -Var_TM(8))*(tm_bnd   -Var_TM(9)) / &
                       ( (Var_TM(7)-Var_TM(8))*(Var_TM(7)-Var_TM(9)) )
              WC_TM(2) = (tm_bnd   -Var_TM(7))*(tm_bnd   -Var_TM(9)) / &
                       ( (Var_TM(8)-Var_TM(7))*(Var_TM(8)-Var_TM(9)) )
              WC_TM(3) = (tm_bnd   -Var_TM(7))*(tm_bnd   -Var_TM(8)) / &
                       ( (Var_TM(9)-Var_TM(7))*(Var_TM(9)-Var_TM(8)) )
           endif
         else
            write(*,*) 'the cross section of this range is extrapolated! XSsetFB: TMO= ',tm_bnd, ' TFU= ', tf_bnd
            write(*,*) 'at node position (Ixy, Iz) ', Ixy, Iz
            write(*,*) 'i_type:   ', i_type
            write(*,*) 'if_state: ', if_state
            write(*,*) 'tm=       ', var_tm(:)
            write(*,*) 'tf=       ', var_tf(:)
            stop
         endif
      elseif ((i_type>=3) .and. (if_branch==4)) then
         if (tm_bnd>=Var_ref_TM(1)) then
           !Hot state
           if_state = 1
           W_TM(4) = (tm_bnd       -Var_ref_TM(3)) / &
                   ( (Var_ref_TM(1)-Var_ref_TM(3)) )
           W_TM(5) = 1d0 - W_TM(4)
         elseif (tm_bnd<Var_ref_TM(1)) then
           !Cold state
           if_state = 2
           ! calculate the weight factors for TM
           W_TM(4) = (tm_bnd       -Var_ref_TM(5))*(tm_bnd       -Var_ref_TM(7))*(tm_bnd       -Var_ref_TM(1)) / &
                   ( (Var_ref_TM(4)-Var_ref_TM(5))*(Var_ref_TM(4)-Var_ref_TM(7))*(Var_ref_TM(4)-Var_ref_TM(1)) )
           W_TM(5) = (tm_bnd       -Var_ref_TM(4))*(tm_bnd       -Var_ref_TM(7))*(tm_bnd       -Var_ref_TM(1)) / &
                   ( (Var_ref_TM(5)-Var_ref_TM(4))*(Var_ref_TM(5)-Var_ref_TM(7))*(Var_ref_TM(5)-Var_ref_TM(1)) )
           W_TM(6) = (tm_bnd       -Var_ref_TM(4))*(tm_bnd       -Var_ref_TM(5))*(tm_bnd       -Var_ref_TM(1)) / &
                   ( (Var_ref_TM(7)-Var_ref_TM(4))*(Var_ref_TM(7)-Var_ref_TM(5))*(Var_ref_TM(7)-Var_ref_TM(1)) )
           W_TM(7) = (tm_bnd       -Var_ref_TM(4))*(tm_bnd       -Var_ref_TM(5))*(tm_bnd       -Var_ref_TM(7)) / &
                   ( (Var_ref_TM(1)-Var_ref_TM(4))*(Var_ref_TM(1)-Var_ref_TM(5))*(Var_ref_TM(1)-Var_ref_TM(7)) )
         endif
      else
         if_state=1
         ! calculate the weight factors for TF
         W_TF(1) = (sqrt(tf_bnd)   -sqrt(Var_TF(2)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                 ( (sqrt(Var_TF(1))-sqrt(Var_TF(2)))*(sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
         W_TF(2) = (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                 ( (sqrt(Var_TF(2))-sqrt(Var_TF(1)))*(sqrt(Var_TF(2))-sqrt(Var_TF(3))) )
         W_TF(3) = (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                 ( (sqrt(Var_TF(3))-sqrt(Var_TF(1)))*(sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
         W_TF(4) = (sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                 ( (sqrt(Var_TF(2))-sqrt(Var_TF(3))) )
         W_TF(5) = 1d0 - W_TF(4)
         ! calculate the weight factors for TM
         W_TM(1) = (tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(1)-Var_TM(2))*(Var_TM(1)-Var_TM(3)) )
         W_TM(2) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(2)-Var_TM(1))*(Var_TM(2)-Var_TM(3)) )
         W_TM(3) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(2)) / &
                 ( (Var_TM(3)-Var_TM(1))*(Var_TM(3)-Var_TM(2)) )
         W_TM(4) = (tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(1)-Var_TM(3)) )
         W_TM(5) = 1d0 - W_TM(4)
      endif

      if (i_type<3) then
         if (if_state==1 .and. if_branch==4) then
            I_Reg = 1
            W_SP_1(17) = W_TF(1)*W_TM(1)*W_PPM(1)
            W_SP_1(18) = W_TF(1)*W_TM(1)*W_PPM(2)
            W_SP_1(19) = W_TF(1)*W_TM(1)*W_PPM(3)
            W_SP_1(20) = W_TF(1)*W_TM(1)*W_PPM(4)
            W_SP_1(14) = W_TF(2)*W_TM(1)*W_PPM(1)
            W_SP_1(13) = W_TF(2)*W_TM(1)*W_PPM(2)
            W_SP_1(15) = W_TF(2)*W_TM(1)*W_PPM(3)
            W_SP_1(16) = W_TF(2)*W_TM(1)*W_PPM(4)
            W_SP_1(21) = W_TF(3)*W_TM(1)*W_PPM(1)
            W_SP_1(22) = W_TF(3)*W_TM(1)*W_PPM(2)
            W_SP_1(23) = W_TF(3)*W_TM(1)*W_PPM(3)
            W_SP_1(24) = W_TF(3)*W_TM(1)*W_PPM(4)
            W_SP_1( 5) = W_TF(1)*W_TM(2)*W_PPM(1)
            W_SP_1( 6) = W_TF(1)*W_TM(2)*W_PPM(2)
            W_SP_1( 7) = W_TF(1)*W_TM(2)*W_PPM(3)
            W_SP_1( 8) = W_TF(1)*W_TM(2)*W_PPM(4)
            W_SP_1( 2) = W_TF(2)*W_TM(2)*W_PPM(1)
            W_SP_1( 1) = W_TF(2)*W_TM(2)*W_PPM(2) -1d0
            W_SP_1( 3) = W_TF(2)*W_TM(2)*W_PPM(3)
            W_SP_1( 4) = W_TF(2)*W_TM(2)*W_PPM(4)
            W_SP_1( 9) = W_TF(3)*W_TM(2)*W_PPM(1)
            W_SP_1(10) = W_TF(3)*W_TM(2)*W_PPM(2)
            W_SP_1(11) = W_TF(3)*W_TM(2)*W_PPM(3)
            W_SP_1(12) = W_TF(3)*W_TM(2)*W_PPM(4)
            W_SP_1(29) = W_TF(1)*W_TM(3)*W_PPM(1)
            W_SP_1(30) = W_TF(1)*W_TM(3)*W_PPM(2)
            W_SP_1(31) = W_TF(1)*W_TM(3)*W_PPM(3)
            W_SP_1(32) = W_TF(1)*W_TM(3)*W_PPM(4)
            W_SP_1(26) = W_TF(2)*W_TM(3)*W_PPM(1)
            W_SP_1(25) = W_TF(2)*W_TM(3)*W_PPM(2)
            W_SP_1(27) = W_TF(2)*W_TM(3)*W_PPM(3)
            W_SP_1(28) = W_TF(2)*W_TM(3)*W_PPM(4)
            W_SP_1(33) = W_TF(3)*W_TM(3)*W_PPM(1)
            W_SP_1(34) = W_TF(3)*W_TM(3)*W_PPM(2)
            W_SP_1(35) = W_TF(3)*W_TM(3)*W_PPM(3)
            W_SP_1(36) = W_TF(3)*W_TM(3)*W_PPM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 2
            W_SP_1(37) = W_TM(2)
            W_SP_1(38) = W_TM(1)
            W_SP_1(39) = W_TM(3)
            W_SP_1( 1) = W_SP_1( 1) - W_TM(2)
            W_SP_1(13) = W_SP_1(13) - W_TM(1)
            W_SP_1(25) = W_SP_1(25) - W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 3
            W_SP_1(37) = 0d0
            W_SP_1(38) = 0d0
            W_SP_1(39) = 0d0
            W_SP_1(40) = W_TM(2)
            W_SP_1(41) = W_TM(1)
            W_SP_1(42) = W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 4
            W_SP_1(40) = 0d0
            W_SP_1(41) = 0d0
            W_SP_1(42) = 0d0
            W_SP_1(43) = W_TM(2)
            W_SP_1(44) = W_TM(1)
            W_SP_1(45) = W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 5
            W_SP_1(43) = 0d0
            W_SP_1(44) = 0d0
            W_SP_1(45) = 0d0
            W_SP_1(46) = W_TM(2)
            W_SP_1(47) = W_TM(1)
            W_SP_1(48) = W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
         elseif(if_state==1 .and. if_branch/=4) then
            I_Reg = 1
            W_SP_1(1)  = W_TF (4) * W_TM(2) + W_PPM(2) * W_TM(2) - W_TM(2) - 1d0
            W_SP_1(2)  = W_PPM(1) * W_TM(2)
            W_SP_1(3)  = W_PPM(3) * W_TM(2)
            W_SP_1(4)  = W_PPM(4) * W_TM(2)
            W_SP_1(5)  = W_TF (2) * W_TM(1) + W_PPM(2) * W_TM(1) - W_TM(1)
            W_SP_1(6)  = W_TF (4) * W_TM(3) + W_PPM(2) * W_TM(3) - W_TM(3)
            W_SP_1(7)  = W_PPM(1) * W_TM(1)
            W_SP_1(8)  = W_PPM(3) * W_TM(1)
            W_SP_1(9)  = W_PPM(4) * W_TM(1)
            W_SP_1(10) = W_TF (1) * W_TM(1)
            W_SP_1(11) = W_TF (3) * W_TM(1)
            W_SP_1(12) = W_TF (5) * W_TM(2)
            W_SP_1(13) = W_PPM(1) * W_TM(3)
            W_SP_1(14) = W_PPM(3) * W_TM(3)
            W_SP_1(15) = W_PPM(4) * W_TM(3)
            W_SP_1(16) = W_TF (5) * W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 2
            W_SP_1(17) = W_TM(2)
            W_SP_1(18) = W_TM(1)
            W_SP_1(19) = W_TM(3)
            W_SP_1(1) = W_SP_1(1) - W_TM(2)
            W_SP_1(5) = W_SP_1(5) - W_TM(1)
            W_SP_1(6) = W_SP_1(6) - W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 3
            W_SP_1(17) = 0d0
            W_SP_1(18) = 0d0
            W_SP_1(19) = 0d0
            W_SP_1(20) = W_TM(2)
            W_SP_1(21) = W_TM(1)
            W_SP_1(22) = W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 4
            W_SP_1(20) = 0d0
            W_SP_1(21) = 0d0
            W_SP_1(22) = 0d0
            W_SP_1(23) = W_TM(2)
            W_SP_1(24) = W_TM(1)
            W_SP_1(25) = W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 5
            W_SP_1(23) = 0d0
            W_SP_1(24) = 0d0
            W_SP_1(25) = 0d0
            W_SP_1(26) = W_TM(2)
            W_SP_1(27) = W_TM(1)
            W_SP_1(28) = W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
         elseif(if_state==2) then
            I_Reg = 1
            W_SP_1( 53) = WC_TF(1)*WC_TM(1)*W_PPM(1)
            W_SP_1( 54) = WC_TF(1)*WC_TM(1)*W_PPM(2)
            W_SP_1( 55) = WC_TF(1)*WC_TM(1)*W_PPM(3)
            W_SP_1( 56) = WC_TF(1)*WC_TM(1)*W_PPM(4)
            W_SP_1( 50) = WC_TF(2)*WC_TM(1)*W_PPM(1)
            W_SP_1( 49) = WC_TF(2)*WC_TM(1)*W_PPM(2)
            W_SP_1( 51) = WC_TF(2)*WC_TM(1)*W_PPM(3)
            W_SP_1( 52) = WC_TF(2)*WC_TM(1)*W_PPM(4)
            W_SP_1( 61) = WC_TF(1)*WC_TM(2)*W_PPM(1)
            W_SP_1( 62) = WC_TF(1)*WC_TM(2)*W_PPM(2)
            W_SP_1( 63) = WC_TF(1)*WC_TM(2)*W_PPM(3)
            W_SP_1( 64) = WC_TF(1)*WC_TM(2)*W_PPM(4)
            W_SP_1( 58) = WC_TF(2)*WC_TM(2)*W_PPM(1)
            W_SP_1( 57) = WC_TF(2)*WC_TM(2)*W_PPM(2)
            W_SP_1( 59) = WC_TF(2)*WC_TM(2)*W_PPM(3)
            W_SP_1( 60) = WC_TF(2)*WC_TM(2)*W_PPM(4)
            W_SP_1( 69) = WC_TF(1)*WC_TM(3)*W_PPM(1)
            W_SP_1( 70) = WC_TF(1)*WC_TM(3)*W_PPM(2)
            W_SP_1( 71) = WC_TF(1)*WC_TM(3)*W_PPM(3)
            W_SP_1( 72) = WC_TF(1)*WC_TM(3)*W_PPM(4)
            W_SP_1( 66) = WC_TF(2)*WC_TM(3)*W_PPM(1)
            W_SP_1( 65) = WC_TF(2)*WC_TM(3)*W_PPM(2)
            W_SP_1( 67) = WC_TF(2)*WC_TM(3)*W_PPM(3)
            W_SP_1( 68) = WC_TF(2)*WC_TM(3)*W_PPM(4)
            W_SP_1( 77) = WC_TF(1)*WC_TM(4)*W_PPM(1)
            W_SP_1( 78) = WC_TF(1)*WC_TM(4)*W_PPM(2)
            W_SP_1( 79) = WC_TF(1)*WC_TM(4)*W_PPM(3)
            W_SP_1( 80) = WC_TF(1)*WC_TM(4)*W_PPM(4)
            W_SP_1( 74) = WC_TF(2)*WC_TM(4)*W_PPM(1)
            W_SP_1( 73) = WC_TF(2)*WC_TM(4)*W_PPM(2)
            W_SP_1( 75) = WC_TF(2)*WC_TM(4)*W_PPM(3)
            W_SP_1( 76) = WC_TF(2)*WC_TM(4)*W_PPM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 2
            W_SP_1(81) = WC_TM(1)
            W_SP_1(82) = WC_TM(2)
            W_SP_1(83) = WC_TM(3)
            W_SP_1(84) = WC_TM(4)
            W_SP_1(49) = W_SP_1(49) - WC_TM(1)
            W_SP_1(57) = W_SP_1(57) - WC_TM(2)
            W_SP_1(65) = W_SP_1(65) - WC_TM(3)
            W_SP_1(73) = W_SP_1(73) - WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 3
            W_SP_1(81) = 0d0
            W_SP_1(82) = 0d0
            W_SP_1(83) = 0d0
            W_SP_1(84) = 0d0
            W_SP_1(85) = WC_TM(1)
            W_SP_1(86) = WC_TM(2)
            W_SP_1(87) = WC_TM(3)
            W_SP_1(88) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 4
            W_SP_1(85) = 0d0
            W_SP_1(86) = 0d0
            W_SP_1(87) = 0d0
            W_SP_1(88) = 0d0
            W_SP_1(89) = WC_TM(1)
            W_SP_1(90) = WC_TM(2)
            W_SP_1(91) = WC_TM(3)
            W_SP_1(92) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 5
            W_SP_1(89) = 0d0
            W_SP_1(90) = 0d0
            W_SP_1(91) = 0d0
            W_SP_1(92) = 0d0
            W_SP_1(93) = WC_TM(1)
            W_SP_1(94) = WC_TM(2)
            W_SP_1(95) = WC_TM(3)
            W_SP_1(96) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

         elseif(if_state==3) then
            I_Reg = 1
            W_SP_1( 77) = WC_TF(1)*WC_TM(1)*W_PPM(1)
            W_SP_1( 78) = WC_TF(1)*WC_TM(1)*W_PPM(2)
            W_SP_1( 79) = WC_TF(1)*WC_TM(1)*W_PPM(3)
            W_SP_1( 80) = WC_TF(1)*WC_TM(1)*W_PPM(4)
            W_SP_1( 74) = WC_TF(2)*WC_TM(1)*W_PPM(1)
            W_SP_1( 73) = WC_TF(2)*WC_TM(1)*W_PPM(2)
            W_SP_1( 75) = WC_TF(2)*WC_TM(1)*W_PPM(3)
            W_SP_1( 76) = WC_TF(2)*WC_TM(1)*W_PPM(4)
            W_SP_1( 97) = WC_TF(3)*WC_TM(1)*W_PPM(1)
            W_SP_1( 98) = WC_TF(3)*WC_TM(1)*W_PPM(2)
            W_SP_1( 99) = WC_TF(3)*WC_TM(1)*W_PPM(3)
            W_SP_1(100) = WC_TF(3)*WC_TM(1)*W_PPM(4)
            W_SP_1(102) = WC_TF(1)*WC_TM(2)*W_PPM(1)
            W_SP_1(101) = WC_TF(1)*WC_TM(2)*W_PPM(2)
            W_SP_1(103) = WC_TF(1)*WC_TM(2)*W_PPM(3)
            W_SP_1(104) = WC_TF(1)*WC_TM(2)*W_PPM(4)
            W_SP_1(105) = WC_TF(2)*WC_TM(2)*W_PPM(1)
            W_SP_1(106) = WC_TF(2)*WC_TM(2)*W_PPM(2)
            W_SP_1(107) = WC_TF(2)*WC_TM(2)*W_PPM(3)
            W_SP_1(108) = WC_TF(2)*WC_TM(2)*W_PPM(4)
            W_SP_1(109) = WC_TF(3)*WC_TM(2)*W_PPM(1)
            W_SP_1(110) = WC_TF(3)*WC_TM(2)*W_PPM(2)
            W_SP_1(111) = WC_TF(3)*WC_TM(2)*W_PPM(3)
            W_SP_1(112) = WC_TF(3)*WC_TM(2)*W_PPM(4)
            W_SP_1(114) = WC_TF(1)*WC_TM(3)*W_PPM(1)
            W_SP_1(113) = WC_TF(1)*WC_TM(3)*W_PPM(2)
            W_SP_1(115) = WC_TF(1)*WC_TM(3)*W_PPM(3)
            W_SP_1(116) = WC_TF(1)*WC_TM(3)*W_PPM(4)
            W_SP_1(117) = WC_TF(2)*WC_TM(3)*W_PPM(1)
            W_SP_1(118) = WC_TF(2)*WC_TM(3)*W_PPM(2)
            W_SP_1(119) = WC_TF(2)*WC_TM(3)*W_PPM(3)
            W_SP_1(120) = WC_TF(2)*WC_TM(3)*W_PPM(4)
            W_SP_1(121) = WC_TF(3)*WC_TM(3)*W_PPM(1)
            W_SP_1(122) = WC_TF(3)*WC_TM(3)*W_PPM(2)
            W_SP_1(123) = WC_TF(3)*WC_TM(3)*W_PPM(3)
            W_SP_1(124) = WC_TF(3)*WC_TM(3)*W_PPM(4)
            W_SP_1(126) = WC_TF(1)*WC_TM(4)*W_PPM(1)
            W_SP_1(125) = WC_TF(1)*WC_TM(4)*W_PPM(2)
            W_SP_1(127) = WC_TF(1)*WC_TM(4)*W_PPM(3)
            W_SP_1(128) = WC_TF(1)*WC_TM(4)*W_PPM(4)
            W_SP_1(129) = WC_TF(2)*WC_TM(4)*W_PPM(1)
            W_SP_1(130) = WC_TF(2)*WC_TM(4)*W_PPM(2)
            W_SP_1(131) = WC_TF(2)*WC_TM(4)*W_PPM(3)
            W_SP_1(132) = WC_TF(2)*WC_TM(4)*W_PPM(4)
            W_SP_1(133) = WC_TF(3)*WC_TM(4)*W_PPM(1)
            W_SP_1(134) = WC_TF(3)*WC_TM(4)*W_PPM(2)
            W_SP_1(135) = WC_TF(3)*WC_TM(4)*W_PPM(3)
            W_SP_1(136) = WC_TF(3)*WC_TM(4)*W_PPM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 2
            W_SP_1( 84) = WC_TM(1)
            W_SP_1(137) = WC_TM(2)
            W_SP_1(138) = WC_TM(3)
            W_SP_1(139) = WC_TM(4)
            W_SP_1( 73) = W_SP_1( 73) - WC_TM(1)
            W_SP_1(106) = W_SP_1(106) - WC_TM(2)
            W_SP_1(118) = W_SP_1(118) - WC_TM(3)
            W_SP_1(130) = W_SP_1(130) - WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 3
            W_SP_1( 84) = 0d0
            W_SP_1(137) = 0d0
            W_SP_1(138) = 0d0
            W_SP_1(139) = 0d0
            W_SP_1( 88) = WC_TM(1)
            W_SP_1(140) = WC_TM(2)
            W_SP_1(141) = WC_TM(3)
            W_SP_1(142) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 4
            W_SP_1( 88) = 0d0
            W_SP_1(140) = 0d0
            W_SP_1(141) = 0d0
            W_SP_1(142) = 0d0
            W_SP_1( 92) = WC_TM(1)
            W_SP_1(143) = WC_TM(2)
            W_SP_1(144) = WC_TM(3)
            W_SP_1(145) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 5
            W_SP_1( 92) = 0d0
            W_SP_1(143) = 0d0
            W_SP_1(144) = 0d0
            W_SP_1(145) = 0d0
            W_SP_1( 96) = WC_TM(1)
            W_SP_1(146) = WC_TM(2)
            W_SP_1(147) = WC_TM(3)
            W_SP_1(148) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

         elseif(if_state==4) then
            I_Reg = 1
            W_SP_1(126) = WC_TF(1)*WC_TM(1)*W_PPM(1)
            W_SP_1(125) = WC_TF(1)*WC_TM(1)*W_PPM(2)
            W_SP_1(127) = WC_TF(1)*WC_TM(1)*W_PPM(3)
            W_SP_1(128) = WC_TF(1)*WC_TM(1)*W_PPM(4)
            W_SP_1(129) = WC_TF(2)*WC_TM(1)*W_PPM(1)
            W_SP_1(130) = WC_TF(2)*WC_TM(1)*W_PPM(2)
            W_SP_1(131) = WC_TF(2)*WC_TM(1)*W_PPM(3)
            W_SP_1(132) = WC_TF(2)*WC_TM(1)*W_PPM(4)
            W_SP_1(133) = WC_TF(3)*WC_TM(1)*W_PPM(1)
            W_SP_1(134) = WC_TF(3)*WC_TM(1)*W_PPM(2)
            W_SP_1(135) = WC_TF(3)*WC_TM(1)*W_PPM(3)
            W_SP_1(136) = WC_TF(3)*WC_TM(1)*W_PPM(4)
            W_SP_1(150) = WC_TF(1)*WC_TM(2)*W_PPM(1)
            W_SP_1(149) = WC_TF(1)*WC_TM(2)*W_PPM(2)
            W_SP_1(151) = WC_TF(1)*WC_TM(2)*W_PPM(3)
            W_SP_1(152) = WC_TF(1)*WC_TM(2)*W_PPM(4)
            W_SP_1(153) = WC_TF(2)*WC_TM(2)*W_PPM(1)
            W_SP_1(154) = WC_TF(2)*WC_TM(2)*W_PPM(2)
            W_SP_1(155) = WC_TF(2)*WC_TM(2)*W_PPM(3)
            W_SP_1(156) = WC_TF(2)*WC_TM(2)*W_PPM(4)
            W_SP_1(157) = WC_TF(3)*WC_TM(2)*W_PPM(1)
            W_SP_1(158) = WC_TF(3)*WC_TM(2)*W_PPM(2)
            W_SP_1(159) = WC_TF(3)*WC_TM(2)*W_PPM(3)
            W_SP_1(160) = WC_TF(3)*WC_TM(2)*W_PPM(4)
            W_SP_1(161) = WC_TF(1)*WC_TM(3)*W_PPM(1)
            W_SP_1(162) = WC_TF(1)*WC_TM(3)*W_PPM(2)
            W_SP_1(163) = WC_TF(1)*WC_TM(3)*W_PPM(3)
            W_SP_1(164) = WC_TF(1)*WC_TM(3)*W_PPM(4)
            W_SP_1( 17) = WC_TF(2)*WC_TM(3)*W_PPM(1)
            W_SP_1( 18) = WC_TF(2)*WC_TM(3)*W_PPM(2)
            W_SP_1( 19) = WC_TF(2)*WC_TM(3)*W_PPM(3)
            W_SP_1( 20) = WC_TF(2)*WC_TM(3)*W_PPM(4)
            W_SP_1( 14) = WC_TF(3)*WC_TM(3)*W_PPM(1)
            W_SP_1( 13) = WC_TF(3)*WC_TM(3)*W_PPM(2)
            W_SP_1( 15) = WC_TF(3)*WC_TM(3)*W_PPM(3)
            W_SP_1( 16) = WC_TF(3)*WC_TM(3)*W_PPM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 2
            W_SP_1(139) = WC_TM(1)
            W_SP_1(165) = WC_TM(2)
            W_SP_1(166) = WC_TM(3)
            W_SP_1(130) = W_SP_1(130) - WC_TM(1)
            W_SP_1(154) = W_SP_1(154) - WC_TM(2)
            W_SP_1( 18) = W_SP_1( 18) - WC_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 3
            W_SP_1(139) = 0d0
            W_SP_1(165) = 0d0
            W_SP_1(166) = 0d0
            W_SP_1(142) = WC_TM(1)
            W_SP_1(167) = WC_TM(2)
            W_SP_1(168) = WC_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 4
            W_SP_1(142) = 0d0
            W_SP_1(167) = 0d0
            W_SP_1(168) = 0d0
            W_SP_1(145) = WC_TM(1)
            W_SP_1(169) = WC_TM(2)
            W_SP_1(170) = WC_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
            I_Reg = 5
            W_SP_1(145) = 0d0
            W_SP_1(169) = 0d0
            W_SP_1(170) = 0d0
            W_SP_1(148) = WC_TM(1)
            W_SP_1(171) = WC_TM(2)
            W_SP_1(172) = WC_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
         else
            write(*,*) 'abnormal interpolation!! XSsetFB'
            stop
         endif

         do i=1,N_SP_FA
            if ((abs(W_SP(i,Ixy,Iz))<1d-8).and.(abs(W_SP(i,Ixy,Iz))>1d-30)) then
               W_SP(i,Ixy,Iz) = 0d0
            endif
         enddo

         if (if_state==1 .and. if_branch==4) then
            MacroXSset(:) = macro_xs_table(ixy,iz)%xs(:, 1)  + &
                            macro_xs_table(ixy,iz)%xs(:, 1) * W_SP( 1, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 2) * W_SP( 2, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 3) * W_SP( 3, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 4) * W_SP( 4, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 5) * W_SP( 5, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 6) * W_SP( 6, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 7) * W_SP( 7, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 8) * W_SP( 8, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 9) * W_SP( 9, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,10) * W_SP(10, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,11) * W_SP(11, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,12) * W_SP(12, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,13) * W_SP(13, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,14) * W_SP(14, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,15) * W_SP(15, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,16) * W_SP(16, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,17) * W_SP(17, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,18) * W_SP(18, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,19) * W_SP(19, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,20) * W_SP(20, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,21) * W_SP(21, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,22) * W_SP(22, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,23) * W_SP(23, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,24) * W_SP(24, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,25) * W_SP(25, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,26) * W_SP(26, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,27) * W_SP(27, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,28) * W_SP(28, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,29) * W_SP(29, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,30) * W_SP(30, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,31) * W_SP(31, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,32) * W_SP(32, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,33) * W_SP(33, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,34) * W_SP(34, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,35) * W_SP(35, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,36) * W_SP(36, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,37) * W_SP(37, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,38) * W_SP(38, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,39) * W_SP(39, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,40) * W_SP(40, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,41) * W_SP(41, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,42) * W_SP(42, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,43) * W_SP(43, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,44) * W_SP(44, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,45) * W_SP(45, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,46) * W_SP(46, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,47) * W_SP(47, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,48) * W_SP(48, Ixy, Iz)
         elseif (if_state==1 .and. if_branch/=4) then
            MacroXSset(:) = macro_xs_table(ixy,iz)%xs(:, 1)  + &
                            macro_xs_table(ixy,iz)%xs(:, 1) * W_SP( 1, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 2) * W_SP( 2, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 3) * W_SP( 3, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 4) * W_SP( 4, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 5) * W_SP( 5, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 6) * W_SP( 6, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 7) * W_SP( 7, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 8) * W_SP( 8, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 9) * W_SP( 9, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,10) * W_SP(10, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,11) * W_SP(11, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,12) * W_SP(12, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,13) * W_SP(13, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,14) * W_SP(14, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,15) * W_SP(15, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,16) * W_SP(16, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,17) * W_SP(17, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,18) * W_SP(18, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,19) * W_SP(19, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,20) * W_SP(20, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,21) * W_SP(21, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,22) * W_SP(22, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,23) * W_SP(23, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,24) * W_SP(24, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,25) * W_SP(25, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,26) * W_SP(26, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,27) * W_SP(27, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,28) * W_SP(28, Ixy, Iz)
         elseif (if_state==2) then
            MacroXSset(:) = macro_xs_table(ixy,iz)%xs(:, 1) - &
                            macro_xs_table(ixy,iz)%xs(:, 1) + &
                            macro_xs_table(ixy,iz)%xs(:,53) * W_SP(53, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,54) * W_SP(54, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,55) * W_SP(55, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,56) * W_SP(56, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,50) * W_SP(50, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,49) * W_SP(49, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,51) * W_SP(51, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,52) * W_SP(52, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,61) * W_SP(61, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,62) * W_SP(62, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,63) * W_SP(63, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,64) * W_SP(64, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,58) * W_SP(58, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,57) * W_SP(57, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,59) * W_SP(59, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,60) * W_SP(60, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,69) * W_SP(69, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,70) * W_SP(70, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,71) * W_SP(71, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,72) * W_SP(72, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,66) * W_SP(66, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,65) * W_SP(65, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,67) * W_SP(67, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,68) * W_SP(68, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,77) * W_SP(77, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,78) * W_SP(78, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,79) * W_SP(79, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,80) * W_SP(80, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,74) * W_SP(74, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,73) * W_SP(73, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,75) * W_SP(75, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,76) * W_SP(76, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,81) * W_SP(81, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,82) * W_SP(82, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,83) * W_SP(83, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,84) * W_SP(84, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,85) * W_SP(85, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,86) * W_SP(86, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,87) * W_SP(87, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,88) * W_SP(88, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,89) * W_SP(89, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,90) * W_SP(90, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,91) * W_SP(91, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,92) * W_SP(92, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,93) * W_SP(93, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,94) * W_SP(94, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,95) * W_SP(95, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,96) * W_SP(96, Ixy, Iz)
         elseif (if_state==3) then
            MacroXSset(:) = macro_xs_table(ixy,iz)%xs(:, 1) - &
                            macro_xs_table(ixy,iz)%xs(:, 1) + &
                            macro_xs_table(ixy,iz)%xs(:, 77) * W_SP( 77, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 78) * W_SP( 78, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 79) * W_SP( 79, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 80) * W_SP( 80, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 74) * W_SP( 74, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 73) * W_SP( 73, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 75) * W_SP( 75, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 76) * W_SP( 76, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 97) * W_SP( 97, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 98) * W_SP( 98, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 99) * W_SP( 99, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,100) * W_SP(100, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,102) * W_SP(102, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,101) * W_SP(101, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,103) * W_SP(103, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,104) * W_SP(104, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,105) * W_SP(105, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,106) * W_SP(106, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,107) * W_SP(107, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,108) * W_SP(108, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,109) * W_SP(109, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,110) * W_SP(110, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,111) * W_SP(111, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,112) * W_SP(112, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,114) * W_SP(114, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,113) * W_SP(113, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,115) * W_SP(115, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,116) * W_SP(116, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,117) * W_SP(117, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,118) * W_SP(118, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,119) * W_SP(119, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,120) * W_SP(120, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,121) * W_SP(121, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,122) * W_SP(122, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,123) * W_SP(123, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,124) * W_SP(124, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,126) * W_SP(126, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,125) * W_SP(125, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,127) * W_SP(127, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,128) * W_SP(128, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,129) * W_SP(129, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,130) * W_SP(130, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,131) * W_SP(131, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,132) * W_SP(132, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,133) * W_SP(133, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,134) * W_SP(134, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,135) * W_SP(135, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,136) * W_SP(136, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 84) * W_SP( 84, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,137) * W_SP(137, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,138) * W_SP(138, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,139) * W_SP(139, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 88) * W_SP( 88, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,140) * W_SP(140, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,141) * W_SP(141, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,142) * W_SP(142, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 92) * W_SP( 92, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,143) * W_SP(143, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,144) * W_SP(144, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,145) * W_SP(145, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 96) * W_SP( 96, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,146) * W_SP(146, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,147) * W_SP(147, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,148) * W_SP(148, Ixy, Iz)
         elseif (if_state==4) then
            MacroXSset(:) = macro_xs_table(ixy,iz)%xs(:, 1) - &
                            macro_xs_table(ixy,iz)%xs(:, 1) + &
                            macro_xs_table(ixy,iz)%xs(:,126) * W_SP(126, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,125) * W_SP(125, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,127) * W_SP(127, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,128) * W_SP(128, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,129) * W_SP(129, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,130) * W_SP(130, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,131) * W_SP(131, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,132) * W_SP(132, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,133) * W_SP(133, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,134) * W_SP(134, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,135) * W_SP(135, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,136) * W_SP(136, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,150) * W_SP(150, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,149) * W_SP(149, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,151) * W_SP(151, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,152) * W_SP(152, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,153) * W_SP(153, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,154) * W_SP(154, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,155) * W_SP(155, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,156) * W_SP(156, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,157) * W_SP(157, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,158) * W_SP(158, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,159) * W_SP(159, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,160) * W_SP(160, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,161) * W_SP(161, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,162) * W_SP(162, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,163) * W_SP(163, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,164) * W_SP(164, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 17) * W_SP( 17, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 18) * W_SP( 18, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 19) * W_SP( 19, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 20) * W_SP( 20, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 14) * W_SP( 14, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 13) * W_SP( 13, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 15) * W_SP( 15, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 16) * W_SP( 16, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,139) * W_SP(139, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,165) * W_SP(165, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,166) * W_SP(166, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,142) * W_SP(142, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,167) * W_SP(167, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,168) * W_SP(168, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,145) * W_SP(145, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,169) * W_SP(169, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,170) * W_SP(170, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,148) * W_SP(148, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,171) * W_SP(171, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,172) * W_SP(172, Ixy, Iz)
         endif
      else ! reflector XS
         if (if_state==1) then
            W_SP(1, Ixy, Iz) = W_TM(4) * W_PPM(1)
            W_SP(2, Ixy, Iz) = W_TM(4) * W_PPM(2)
            W_SP(3, Ixy, Iz) = W_TM(4) * W_PPM(3)
            W_SP(4, Ixy, Iz) = W_TM(4) * W_PPM(4)
            W_SP(5, Ixy, Iz) = W_TM(5) * W_PPM(1)
            W_SP(6, Ixy, Iz) = W_TM(5) * W_PPM(2)
            W_SP(7, Ixy, Iz) = W_TM(5) * W_PPM(3)
            W_SP(8, Ixy, Iz) = W_TM(5) * W_PPM(4)
            macroXSset( : ) = macro_xs_table(ixy,iz)%xs(:,1) * W_SP(1, Ixy, Iz) + &
                              macro_xs_table(ixy,iz)%xs(:,2) * W_SP(2, Ixy, Iz) + &
                              macro_xs_table(ixy,iz)%xs(:,3) * W_SP(3, Ixy, Iz) + &
                              macro_xs_table(ixy,iz)%xs(:,4) * W_SP(4, Ixy, Iz) + &
                              macro_xs_table(ixy,iz)%xs(:,5) * W_SP(5, Ixy, Iz) + &
                              macro_xs_table(ixy,iz)%xs(:,6) * W_SP(6, Ixy, Iz) + &
                              macro_xs_table(ixy,iz)%xs(:,7) * W_SP(7, Ixy, Iz) + &
                              macro_xs_table(ixy,iz)%xs(:,8) * W_SP(8, Ixy, Iz)
         elseif (if_state>1) then
            W_SP( 9, Ixy, Iz) = W_TM(4) * W_PPM(1)
            W_SP(10, Ixy, Iz) = W_TM(4) * W_PPM(2)
            W_SP(11, Ixy, Iz) = W_TM(4) * W_PPM(3)
            W_SP(12, Ixy, Iz) = W_TM(4) * W_PPM(4)
            W_SP(13, Ixy, Iz) = W_TM(5) * W_PPM(1)
            W_SP(14, Ixy, Iz) = W_TM(5) * W_PPM(2)
            W_SP(15, Ixy, Iz) = W_TM(5) * W_PPM(3)
            W_SP(16, Ixy, Iz) = W_TM(5) * W_PPM(4)
            W_SP(17, Ixy, Iz) = W_TM(6) * W_PPM(1)
            W_SP(18, Ixy, Iz) = W_TM(6) * W_PPM(2)
            W_SP(19, Ixy, Iz) = W_TM(6) * W_PPM(3)
            W_SP(20, Ixy, Iz) = W_TM(6) * W_PPM(4)
            W_SP(21, Ixy, Iz) = W_TM(7) * W_PPM(1)
            W_SP(22, Ixy, Iz) = W_TM(7) * W_PPM(2)
            W_SP(23, Ixy, Iz) = W_TM(7) * W_PPM(3)
            W_SP(24, Ixy, Iz) = W_TM(7) * W_PPM(4)
            XSset( : ) = macro_xs_table(ixy,iz)%xs(:, 9) * W_SP( 9, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,10) * W_SP(10, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,11) * W_SP(11, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,12) * W_SP(12, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,13) * W_SP(13, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,14) * W_SP(14, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,15) * W_SP(15, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,16) * W_SP(16, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,17) * W_SP(17, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,18) * W_SP(18, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,19) * W_SP(19, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,20) * W_SP(20, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,21) * W_SP(21, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,22) * W_SP(22, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,23) * W_SP(23, Ixy, Iz) + &
                         macro_xs_table(ixy,iz)%xs(:,24) * W_SP(24, Ixy, Iz)
         endif
      endif

      return
      end subroutine MacroXSsetFB


      subroutine MacroXSsetFB_1D(Ixy, Iz)
      use Inc_3D, only: T_Fuel, T_Mod
      use Inc_CR, only: Reg_VF
      use Inc_RP, only: I_LP_1N, AxialComp
      use Inc_TH, only: PPM, TM_In
      use inc_th, only: flag_th_chanwise, i_chan_1n, chanwise_t_inlet
      use Inc_Detector
      use Inc_INP, only: flag_force_tm, force_tm, force_dtm
      use Inc_INP, only: flag_force_tf, force_tf, force_dtf
      use Inc_INP, only: flag_force_tdop, force_tdop, force_dtdop
      use Inc_INP, only: flag_force_ftc, force_ftc, flag_force_mtc, force_mtc
      use Mod_GetNode, only: new_asym_itab
      implicit none
      integer, intent(in) :: Ixy, Iz
      integer :: Ixy_1N
      integer(4) :: i
      real(8) :: tm_in_bak, chg_Tm, chg_tf, tmp_tdop
      real(8) :: tf_bnd, tm_bnd


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [MacroXSsetFB_1D] in Mod_LiteVer'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (flag_th_chanwise) then
         tm_in_bak=chanwise_t_inlet(i_chan_1n(i_4nto1n(ixy)))+degtok
      else
         tm_in_bak=Tm_in+DegToK
      endif

      if (flag_force_tm) then
         chg_Tm=tm_in_bak+(T_Mod(Ixy,Iz)-tm_in_bak)*force_tm+force_dtm
      else
         chg_Tm=T_Mod(Ixy,Iz)
      endif
      if (flag_force_tf) then
         chg_Tf=t_fuel(ixy,iz)*force_tf+force_dtf
      elseif (flag_force_tdop) then
         tmp_tdop=sqrt(t_fuel(ixy,iz))
         tmp_tdop=tmp_tdop*force_tdop+force_dtdop
         chg_Tf=tmp_tdop*tmp_tdop
      else
         chg_Tf=t_fuel(ixy,iz)
      endif
      if (flag_force_ftc) chg_tf=chg_tf+force_ftc
      if (flag_force_mtc) chg_tm=chg_tm+force_mtc

      tf_bnd=min(tf_max,chg_tf)
      tm_bnd=max(tm_min,min(tm_max,chg_Tm))
      W_SP(:, Ixy, Iz) = 0d0
      W_SP_1  = 0d0

      Ixy_1N = I_4Nto1N(Ixy)
      I_Tab = new_asym_itab(axialcomp(i_lp_1n(ixy_1n),iz),ixy)
      I_Type = Type_Tab(I_Tab)

      if_state=0
      if(if_branch==4) WC_TF=0d0
      if(if_branch==4) WC_TM=0d0
      W_PPM=0d0
      W_TM =0d0
      W_TF =0d0

      ! calculate the weight factors for PPM
      W_PPM(1) = (PPM       -Var_PPM(2))*(PPM       -Var_PPM(3))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(1)-Var_PPM(2))*(Var_PPM(1)-Var_PPM(3))*(Var_PPM(1)-Var_PPM(4)) )
      W_PPM(2) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(3))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(2)-Var_PPM(1))*(Var_PPM(2)-Var_PPM(3))*(Var_PPM(2)-Var_PPM(4)) )
      W_PPM(3) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(2))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(3)-Var_PPM(1))*(Var_PPM(3)-Var_PPM(2))*(Var_PPM(3)-Var_PPM(4)) )
      W_PPM(4) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(2))*(PPM       -Var_PPM(3)) / &
               ( (Var_PPM(4)-Var_PPM(1))*(Var_PPM(4)-Var_PPM(2))*(Var_PPM(4)-Var_PPM(3)) )

      if ((i_type<3) .and. (if_branch==4)) then
         if ((tm_bnd>=Var_ref_TM(1)).and.(tf_bnd>=Var_ref_TF(1))) then
            if_state=1
            W_TF(1) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(2)))*(sqrt(tf_bnd)       -sqrt(Var_ref_TF(3))) / &
                    ( (sqrt(Var_ref_TF(1))-sqrt(Var_ref_TF(2)))*(sqrt(Var_ref_TF(1))-sqrt(Var_ref_TF(3))) )
            W_TF(2) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(1)))*(sqrt(tf_bnd)       -sqrt(Var_ref_TF(3))) / &
                    ( (sqrt(Var_ref_TF(2))-sqrt(Var_ref_TF(1)))*(sqrt(Var_ref_TF(2))-sqrt(Var_ref_TF(3))) )
            W_TF(3) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(1)))*(sqrt(tf_bnd)       -sqrt(Var_ref_TF(2))) / &
                    ( (sqrt(Var_ref_TF(3))-sqrt(Var_ref_TF(1)))*(sqrt(Var_ref_TF(3))-sqrt(Var_ref_TF(2))) )
            W_TF(4) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(3))) / &
                    ( (sqrt(Var_ref_TF(2))-sqrt(Var_ref_TF(3))) )
            W_TF(5) = 1d0 - W_TF(4)
            ! calculate the weight factors for TM
            W_TM(1) = (tm_bnd       -Var_ref_TM(2))*(tm_bnd       -Var_ref_TM(3)) / &
                    ( (Var_ref_TM(1)-Var_ref_TM(2))*(Var_ref_TM(1)-Var_ref_TM(3)) )
            W_TM(2) = (tm_bnd       -Var_ref_TM(1))*(tm_bnd       -Var_ref_TM(3)) / &
                    ( (Var_ref_TM(2)-Var_ref_TM(1))*(Var_ref_TM(2)-Var_ref_TM(3)) )
            W_TM(3) = (tm_bnd       -Var_ref_TM(1))*(tm_bnd       -Var_ref_TM(2)) / &
                    ( (Var_ref_TM(3)-Var_ref_TM(1))*(Var_ref_TM(3)-Var_ref_TM(2)) )
            W_TM(4) = (tm_bnd       -Var_ref_TM(3)) / &
                    ( (Var_ref_TM(1)-Var_ref_TM(3)) )
            W_TM(5) = 1d0 - W_TM(4)
         elseif (((tm_bnd<=Var_TM(9)).and.(tf_bnd<=Var_TF(2))).and.((tm_bnd>=Var_TM(1)-5d0).and.(tf_bnd>=Var_TF(1)-5d0))) then
            if(tm_bnd<=400d0) then
               if_state=2
               WC_TF(1)= (sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                       ( (sqrt(Var_TF(1))-sqrt(Var_TF(2))) )
               WC_TF(2)= (sqrt(tf_bnd)   -sqrt(Var_TF(1))) / &
                       ( (sqrt(Var_TF(2))-sqrt(Var_TF(1))) )
               WC_TM(1) = (tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(3))*(tm_bnd   -Var_TM(4)) / &
                        ( (Var_TM(1)-Var_TM(2))*(Var_TM(1)-Var_TM(3))*(Var_TM(1)-Var_TM(4)) )
               WC_TM(2) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(3))*(tm_bnd   -Var_TM(4)) / &
                        ( (Var_TM(2)-Var_TM(1))*(Var_TM(2)-Var_TM(3))*(Var_TM(2)-Var_TM(4)) )
               WC_TM(3) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(4)) / &
                        ( (Var_TM(3)-Var_TM(1))*(Var_TM(3)-Var_TM(2))*(Var_TM(3)-Var_TM(4)) )
               WC_TM(4) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(3)) / &
                        ( (Var_TM(4)-Var_TM(1))*(Var_TM(4)-Var_TM(2))*(Var_TM(4)-Var_TM(3)) )
            elseif(tm_bnd>400d0 .and. tm_bnd<=500d0) then
               if_state=3
               WC_TF(1)= (sqrt(tf_bnd)   -sqrt(Var_TF(2)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(1))-sqrt(Var_TF(2)))*(sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
               WC_TF(2)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(2))-sqrt(Var_TF(1)))*(sqrt(Var_TF(2))-sqrt(Var_TF(3))) )
               WC_TF(3)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                       ( (sqrt(Var_TF(3))-sqrt(Var_TF(1)))*(sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
               WC_TM(1) = (tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(6))*(tm_bnd   -Var_TM(7)) / &
                        ( (Var_TM(4)-Var_TM(5))*(Var_TM(4)-Var_TM(6))*(Var_TM(4)-Var_TM(7)) )
               WC_TM(2) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(6))*(tm_bnd   -Var_TM(7)) / &
                        ( (Var_TM(5)-Var_TM(4))*(Var_TM(5)-Var_TM(6))*(Var_TM(5)-Var_TM(7)) )
               WC_TM(3) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(7)) / &
                        ( (Var_TM(6)-Var_TM(4))*(Var_TM(6)-Var_TM(5))*(Var_TM(6)-Var_TM(7)) )
               WC_TM(4) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(6)) / &
                        ( (Var_TM(7)-Var_TM(4))*(Var_TM(7)-Var_TM(5))*(Var_TM(7)-Var_TM(6)) )
            elseif(tm_bnd>500d0) then
               if_state=4
               WC_TF(1)= (sqrt(tf_bnd)   -sqrt(Var_TF(2)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(1))-sqrt(Var_TF(2)))*(sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
               WC_TF(2)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(2))-sqrt(Var_TF(1)))*(sqrt(Var_TF(2))-sqrt(Var_TF(3))) )
               WC_TF(3)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                       ( (sqrt(Var_TF(3))-sqrt(Var_TF(1)))*(sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
               WC_TM(1) = (tm_bnd   -Var_TM(8))*(tm_bnd   -Var_TM(9)) / &
                        ( (Var_TM(7)-Var_TM(8))*(Var_TM(7)-Var_TM(9)) )
               WC_TM(2) = (tm_bnd   -Var_TM(7))*(tm_bnd   -Var_TM(9)) / &
                        ( (Var_TM(8)-Var_TM(7))*(Var_TM(8)-Var_TM(9)) )
               WC_TM(3) = (tm_bnd   -Var_TM(7))*(tm_bnd   -Var_TM(8)) / &
                        ( (Var_TM(9)-Var_TM(7))*(Var_TM(9)-Var_TM(8)) )
            endif
         else
            write(*,*) 'the cross section of this range is extrapolated! XSsetFB: TMO= ',tm_bnd, ' TFU= ', tf_bnd
            write(*,*) 'at node position (Ixy, Iz) ', Ixy, Iz
            write(*,*) 'i_type:   ', i_type
            write(*,*) 'if_state: ', if_state
            write(*,*) 'tm=       ', var_tm(:)
            write(*,*) 'tf=       ', var_tf(:)
            stop
         endif
      elseif ((i_type>=3) .and. (if_branch==4)) then
         if (tm_bnd>Var_ref_TM(1)) then
           !Hot state
           if_state = 1
           W_TM(4) = (tm_bnd   -Var_ref_TM(3)) / &
                   ( (Var_ref_TM(1)-Var_ref_TM(3)) )
           W_TM(5) = 1d0 - W_TM(4)
         elseif (tm_bnd<=Var_ref_TM(1)) then
           !Cold state
           if_state = 2
           ! calculate the weight factors for TM
           W_TM(4) = (tm_bnd       -Var_ref_TM(5))*(tm_bnd       -Var_ref_TM(7))*(tm_bnd       -Var_ref_TM(1)) / &
                   ( (Var_ref_TM(4)-Var_ref_TM(5))*(Var_ref_TM(4)-Var_ref_TM(7))*(Var_ref_TM(4)-Var_ref_TM(1)) )
           W_TM(5) = (tm_bnd       -Var_ref_TM(4))*(tm_bnd       -Var_ref_TM(7))*(tm_bnd       -Var_ref_TM(1)) / &
                   ( (Var_ref_TM(5)-Var_ref_TM(4))*(Var_ref_TM(5)-Var_ref_TM(7))*(Var_ref_TM(5)-Var_ref_TM(1)) )
           W_TM(6) = (tm_bnd       -Var_ref_TM(4))*(tm_bnd       -Var_ref_TM(5))*(tm_bnd       -Var_ref_TM(1)) / &
                   ( (Var_ref_TM(7)-Var_ref_TM(4))*(Var_ref_TM(7)-Var_ref_TM(5))*(Var_ref_TM(7)-Var_ref_TM(1)) )
           W_TM(7) = (tm_bnd       -Var_ref_TM(4))*(tm_bnd       -Var_ref_TM(5))*(tm_bnd       -Var_ref_TM(7)) / &
                   ( (Var_ref_TM(1)-Var_ref_TM(4))*(Var_ref_TM(1)-Var_ref_TM(5))*(Var_ref_TM(1)-Var_ref_TM(7)) )
         endif
      else
         if_state=1
         ! calculate the weight factors for TF
         W_TF(3) = (sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                 ( (sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
         W_TF(2) = (sqrt(Var_TF(3))-sqrt(tf_bnd)   ) / &
                 ( (sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
         ! calculate the weight factors for TM
         W_TM(1) = (tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(1)-Var_TM(2))*(Var_TM(1)-Var_TM(3)) )
         W_TM(2) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(2)-Var_TM(1))*(Var_TM(2)-Var_TM(3)) )
         W_TM(3) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(2)) / &
                 ( (Var_TM(3)-Var_TM(1))*(Var_TM(3)-Var_TM(2)) )
         W_TM(4) = (tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(1)-Var_TM(3)) )
         W_TM(5) = 1d0 - W_TM(4)
      endif

      if ( I_Type < 3 ) then
         if(if_state==1 .and. if_branch==4) then
            I_Reg = 1
            W_SP_1(17) = W_TF(1)*W_TM(1)*W_PPM(1)
            W_SP_1(18) = W_TF(1)*W_TM(1)*W_PPM(2)
            W_SP_1(19) = W_TF(1)*W_TM(1)*W_PPM(3)
            W_SP_1(20) = W_TF(1)*W_TM(1)*W_PPM(4)
            W_SP_1(14) = W_TF(2)*W_TM(1)*W_PPM(1)
            W_SP_1(13) = W_TF(2)*W_TM(1)*W_PPM(2)
            W_SP_1(15) = W_TF(2)*W_TM(1)*W_PPM(3)
            W_SP_1(16) = W_TF(2)*W_TM(1)*W_PPM(4)
            W_SP_1(21) = W_TF(3)*W_TM(1)*W_PPM(1)
            W_SP_1(22) = W_TF(3)*W_TM(1)*W_PPM(2)
            W_SP_1(23) = W_TF(3)*W_TM(1)*W_PPM(3)
            W_SP_1(24) = W_TF(3)*W_TM(1)*W_PPM(4)
            W_SP_1( 5) = W_TF(1)*W_TM(2)*W_PPM(1)
            W_SP_1( 6) = W_TF(1)*W_TM(2)*W_PPM(2)
            W_SP_1( 7) = W_TF(1)*W_TM(2)*W_PPM(3)
            W_SP_1( 8) = W_TF(1)*W_TM(2)*W_PPM(4)
            W_SP_1( 2) = W_TF(2)*W_TM(2)*W_PPM(1)
            W_SP_1( 1) = W_TF(2)*W_TM(2)*W_PPM(2) -1d0
            W_SP_1( 3) = W_TF(2)*W_TM(2)*W_PPM(3)
            W_SP_1( 4) = W_TF(2)*W_TM(2)*W_PPM(4)
            W_SP_1( 9) = W_TF(3)*W_TM(2)*W_PPM(1)
            W_SP_1(10) = W_TF(3)*W_TM(2)*W_PPM(2)
            W_SP_1(11) = W_TF(3)*W_TM(2)*W_PPM(3)
            W_SP_1(12) = W_TF(3)*W_TM(2)*W_PPM(4)
            W_SP_1(29) = W_TF(1)*W_TM(3)*W_PPM(1)
            W_SP_1(30) = W_TF(1)*W_TM(3)*W_PPM(2)
            W_SP_1(31) = W_TF(1)*W_TM(3)*W_PPM(3)
            W_SP_1(32) = W_TF(1)*W_TM(3)*W_PPM(4)
            W_SP_1(26) = W_TF(2)*W_TM(3)*W_PPM(1)
            W_SP_1(25) = W_TF(2)*W_TM(3)*W_PPM(2)
            W_SP_1(27) = W_TF(2)*W_TM(3)*W_PPM(3)
            W_SP_1(28) = W_TF(2)*W_TM(3)*W_PPM(4)
            W_SP_1(33) = W_TF(3)*W_TM(3)*W_PPM(1)
            W_SP_1(34) = W_TF(3)*W_TM(3)*W_PPM(2)
            W_SP_1(35) = W_TF(3)*W_TM(3)*W_PPM(3)
            W_SP_1(36) = W_TF(3)*W_TM(3)*W_PPM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 2
            W_SP_1(37) = W_TM(2)
            W_SP_1(38) = W_TM(1)
            W_SP_1(39) = W_TM(3)
            W_SP_1( 1) = W_SP_1( 1) - W_TM(2)
            W_SP_1(13) = W_SP_1(13) - W_TM(1)
            W_SP_1(25) = W_SP_1(25) - W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 3
            W_SP_1(37) = 0d0
            W_SP_1(38) = 0d0
            W_SP_1(39) = 0d0
            W_SP_1(40) = W_TM(2)
            W_SP_1(41) = W_TM(1)
            W_SP_1(42) = W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 4
            W_SP_1(40) = 0d0
            W_SP_1(41) = 0d0
            W_SP_1(42) = 0d0
            W_SP_1(43) = W_TM(2)
            W_SP_1(44) = W_TM(1)
            W_SP_1(45) = W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 5
            W_SP_1(43) = 0d0
            W_SP_1(44) = 0d0
            W_SP_1(45) = 0d0
            W_SP_1(46) = W_TM(2)
            W_SP_1(47) = W_TM(1)
            W_SP_1(48) = W_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

         elseif(if_state==1 .and. if_branch/=4) then
            I_Reg = 1
            W_SP_1(1)  = W_PPM(2) + W_TF(2) + W_TM(2) - 3d0
            W_SP_1(2)  = W_PPM(1)
            W_SP_1(3)  = W_PPM(3)
            W_SP_1(4)  = W_PPM(4)
            W_SP_1(5)  = W_TM(1)
            W_SP_1(6)  = W_TM(3)
            W_SP_1(12) = W_TF(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 2
            W_SP_1(17) = 1d0
            W_SP_1(1) = W_SP_1(1) - 1d0
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 3
            W_SP_1(17) = 0d0
            W_SP_1(20) = 1d0
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 4
            W_SP_1(20) = 0d0
            W_SP_1(23) = 1d0
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 5
            W_SP_1(23) = 0d0
            W_SP_1(26) = 1d0
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
         elseif(if_state==2) then
            I_Reg = 1
            W_SP_1( 53) = WC_TF(1)*WC_TM(1)*W_PPM(1)
            W_SP_1( 54) = WC_TF(1)*WC_TM(1)*W_PPM(2)
            W_SP_1( 55) = WC_TF(1)*WC_TM(1)*W_PPM(3)
            W_SP_1( 56) = WC_TF(1)*WC_TM(1)*W_PPM(4)
            W_SP_1( 50) = WC_TF(2)*WC_TM(1)*W_PPM(1)
            W_SP_1( 49) = WC_TF(2)*WC_TM(1)*W_PPM(2)
            W_SP_1( 51) = WC_TF(2)*WC_TM(1)*W_PPM(3)
            W_SP_1( 52) = WC_TF(2)*WC_TM(1)*W_PPM(4)
            W_SP_1( 61) = WC_TF(1)*WC_TM(2)*W_PPM(1)
            W_SP_1( 62) = WC_TF(1)*WC_TM(2)*W_PPM(2)
            W_SP_1( 63) = WC_TF(1)*WC_TM(2)*W_PPM(3)
            W_SP_1( 64) = WC_TF(1)*WC_TM(2)*W_PPM(4)
            W_SP_1( 58) = WC_TF(2)*WC_TM(2)*W_PPM(1)
            W_SP_1( 57) = WC_TF(2)*WC_TM(2)*W_PPM(2)
            W_SP_1( 59) = WC_TF(2)*WC_TM(2)*W_PPM(3)
            W_SP_1( 60) = WC_TF(2)*WC_TM(2)*W_PPM(4)
            W_SP_1( 69) = WC_TF(1)*WC_TM(3)*W_PPM(1)
            W_SP_1( 70) = WC_TF(1)*WC_TM(3)*W_PPM(2)
            W_SP_1( 71) = WC_TF(1)*WC_TM(3)*W_PPM(3)
            W_SP_1( 72) = WC_TF(1)*WC_TM(3)*W_PPM(4)
            W_SP_1( 66) = WC_TF(2)*WC_TM(3)*W_PPM(1)
            W_SP_1( 65) = WC_TF(2)*WC_TM(3)*W_PPM(2)
            W_SP_1( 67) = WC_TF(2)*WC_TM(3)*W_PPM(3)
            W_SP_1( 68) = WC_TF(2)*WC_TM(3)*W_PPM(4)
            W_SP_1( 77) = WC_TF(1)*WC_TM(4)*W_PPM(1)
            W_SP_1( 78) = WC_TF(1)*WC_TM(4)*W_PPM(2)
            W_SP_1( 79) = WC_TF(1)*WC_TM(4)*W_PPM(3)
            W_SP_1( 80) = WC_TF(1)*WC_TM(4)*W_PPM(4)
            W_SP_1( 74) = WC_TF(2)*WC_TM(4)*W_PPM(1)
            W_SP_1( 73) = WC_TF(2)*WC_TM(4)*W_PPM(2)
            W_SP_1( 75) = WC_TF(2)*WC_TM(4)*W_PPM(3)
            W_SP_1( 76) = WC_TF(2)*WC_TM(4)*W_PPM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 2
            W_SP_1(81) = WC_TM(1)
            W_SP_1(82) = WC_TM(2)
            W_SP_1(83) = WC_TM(3)
            W_SP_1(84) = WC_TM(4)
            W_SP_1(49) = W_SP_1(49) - WC_TM(1)
            W_SP_1(57) = W_SP_1(57) - WC_TM(2)
            W_SP_1(65) = W_SP_1(65) - WC_TM(3)
            W_SP_1(73) = W_SP_1(73) - WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 3
            W_SP_1(81) = 0d0
            W_SP_1(82) = 0d0
            W_SP_1(83) = 0d0
            W_SP_1(84) = 0d0
            W_SP_1(85) = WC_TM(1)
            W_SP_1(86) = WC_TM(2)
            W_SP_1(87) = WC_TM(3)
            W_SP_1(88) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 4
            W_SP_1(85) = 0d0
            W_SP_1(86) = 0d0
            W_SP_1(87) = 0d0
            W_SP_1(88) = 0d0
            W_SP_1(89) = WC_TM(1)
            W_SP_1(90) = WC_TM(2)
            W_SP_1(91) = WC_TM(3)
            W_SP_1(92) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 5
            W_SP_1(89) = 0d0
            W_SP_1(90) = 0d0
            W_SP_1(91) = 0d0
            W_SP_1(92) = 0d0
            W_SP_1(93) = WC_TM(1)
            W_SP_1(94) = WC_TM(2)
            W_SP_1(95) = WC_TM(3)
            W_SP_1(96) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
         elseif(if_state==3) then
            I_Reg = 1
            W_SP_1( 77) = WC_TF(1)*WC_TM(1)*W_PPM(1)
            W_SP_1( 78) = WC_TF(1)*WC_TM(1)*W_PPM(2)
            W_SP_1( 79) = WC_TF(1)*WC_TM(1)*W_PPM(3)
            W_SP_1( 80) = WC_TF(1)*WC_TM(1)*W_PPM(4)
            W_SP_1( 74) = WC_TF(2)*WC_TM(1)*W_PPM(1)
            W_SP_1( 73) = WC_TF(2)*WC_TM(1)*W_PPM(2)
            W_SP_1( 75) = WC_TF(2)*WC_TM(1)*W_PPM(3)
            W_SP_1( 76) = WC_TF(2)*WC_TM(1)*W_PPM(4)
            W_SP_1( 97) = WC_TF(3)*WC_TM(1)*W_PPM(1)
            W_SP_1( 98) = WC_TF(3)*WC_TM(1)*W_PPM(2)
            W_SP_1( 99) = WC_TF(3)*WC_TM(1)*W_PPM(3)
            W_SP_1(100) = WC_TF(3)*WC_TM(1)*W_PPM(4)
            W_SP_1(102) = WC_TF(1)*WC_TM(2)*W_PPM(1)
            W_SP_1(101) = WC_TF(1)*WC_TM(2)*W_PPM(2)
            W_SP_1(103) = WC_TF(1)*WC_TM(2)*W_PPM(3)
            W_SP_1(104) = WC_TF(1)*WC_TM(2)*W_PPM(4)
            W_SP_1(105) = WC_TF(2)*WC_TM(2)*W_PPM(1)
            W_SP_1(106) = WC_TF(2)*WC_TM(2)*W_PPM(2)
            W_SP_1(107) = WC_TF(2)*WC_TM(2)*W_PPM(3)
            W_SP_1(108) = WC_TF(2)*WC_TM(2)*W_PPM(4)
            W_SP_1(109) = WC_TF(3)*WC_TM(2)*W_PPM(1)
            W_SP_1(110) = WC_TF(3)*WC_TM(2)*W_PPM(2)
            W_SP_1(111) = WC_TF(3)*WC_TM(2)*W_PPM(3)
            W_SP_1(112) = WC_TF(3)*WC_TM(2)*W_PPM(4)
            W_SP_1(114) = WC_TF(1)*WC_TM(3)*W_PPM(1)
            W_SP_1(113) = WC_TF(1)*WC_TM(3)*W_PPM(2)
            W_SP_1(115) = WC_TF(1)*WC_TM(3)*W_PPM(3)
            W_SP_1(116) = WC_TF(1)*WC_TM(3)*W_PPM(4)
            W_SP_1(117) = WC_TF(2)*WC_TM(3)*W_PPM(1)
            W_SP_1(118) = WC_TF(2)*WC_TM(3)*W_PPM(2)
            W_SP_1(119) = WC_TF(2)*WC_TM(3)*W_PPM(3)
            W_SP_1(120) = WC_TF(2)*WC_TM(3)*W_PPM(4)
            W_SP_1(121) = WC_TF(3)*WC_TM(3)*W_PPM(1)
            W_SP_1(122) = WC_TF(3)*WC_TM(3)*W_PPM(2)
            W_SP_1(123) = WC_TF(3)*WC_TM(3)*W_PPM(3)
            W_SP_1(124) = WC_TF(3)*WC_TM(3)*W_PPM(4)
            W_SP_1(126) = WC_TF(1)*WC_TM(4)*W_PPM(1)
            W_SP_1(125) = WC_TF(1)*WC_TM(4)*W_PPM(2)
            W_SP_1(127) = WC_TF(1)*WC_TM(4)*W_PPM(3)
            W_SP_1(128) = WC_TF(1)*WC_TM(4)*W_PPM(4)
            W_SP_1(129) = WC_TF(2)*WC_TM(4)*W_PPM(1)
            W_SP_1(130) = WC_TF(2)*WC_TM(4)*W_PPM(2)
            W_SP_1(131) = WC_TF(2)*WC_TM(4)*W_PPM(3)
            W_SP_1(132) = WC_TF(2)*WC_TM(4)*W_PPM(4)
            W_SP_1(133) = WC_TF(3)*WC_TM(4)*W_PPM(1)
            W_SP_1(134) = WC_TF(3)*WC_TM(4)*W_PPM(2)
            W_SP_1(135) = WC_TF(3)*WC_TM(4)*W_PPM(3)
            W_SP_1(136) = WC_TF(3)*WC_TM(4)*W_PPM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 2
            W_SP_1( 84) = WC_TM(1)
            W_SP_1(137) = WC_TM(2)
            W_SP_1(138) = WC_TM(3)
            W_SP_1(139) = WC_TM(4)
            W_SP_1( 73) = W_SP_1( 73) - WC_TM(1)
            W_SP_1(106) = W_SP_1(106) - WC_TM(2)
            W_SP_1(118) = W_SP_1(118) - WC_TM(3)
            W_SP_1(130) = W_SP_1(130) - WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 3
            W_SP_1( 84) = 0d0
            W_SP_1(137) = 0d0
            W_SP_1(138) = 0d0
            W_SP_1(139) = 0d0
            W_SP_1( 88) = WC_TM(1)
            W_SP_1(140) = WC_TM(2)
            W_SP_1(141) = WC_TM(3)
            W_SP_1(142) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 4
            W_SP_1( 88) = 0d0
            W_SP_1(140) = 0d0
            W_SP_1(141) = 0d0
            W_SP_1(142) = 0d0
            W_SP_1( 92) = WC_TM(1)
            W_SP_1(143) = WC_TM(2)
            W_SP_1(144) = WC_TM(3)
            W_SP_1(145) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 5
            W_SP_1( 92) = 0d0
            W_SP_1(143) = 0d0
            W_SP_1(144) = 0d0
            W_SP_1(145) = 0d0
            W_SP_1( 96) = WC_TM(1)
            W_SP_1(146) = WC_TM(2)
            W_SP_1(147) = WC_TM(3)
            W_SP_1(148) = WC_TM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

         elseif(if_state==4) then
            I_Reg = 1
            W_SP_1(126) = WC_TF(1)*WC_TM(1)*W_PPM(1)
            W_SP_1(125) = WC_TF(1)*WC_TM(1)*W_PPM(2)
            W_SP_1(127) = WC_TF(1)*WC_TM(1)*W_PPM(3)
            W_SP_1(128) = WC_TF(1)*WC_TM(1)*W_PPM(4)
            W_SP_1(129) = WC_TF(2)*WC_TM(1)*W_PPM(1)
            W_SP_1(130) = WC_TF(2)*WC_TM(1)*W_PPM(2)
            W_SP_1(131) = WC_TF(2)*WC_TM(1)*W_PPM(3)
            W_SP_1(132) = WC_TF(2)*WC_TM(1)*W_PPM(4)
            W_SP_1(133) = WC_TF(3)*WC_TM(1)*W_PPM(1)
            W_SP_1(134) = WC_TF(3)*WC_TM(1)*W_PPM(2)
            W_SP_1(135) = WC_TF(3)*WC_TM(1)*W_PPM(3)
            W_SP_1(136) = WC_TF(3)*WC_TM(1)*W_PPM(4)
            W_SP_1(150) = WC_TF(1)*WC_TM(2)*W_PPM(1)
            W_SP_1(149) = WC_TF(1)*WC_TM(2)*W_PPM(2)
            W_SP_1(151) = WC_TF(1)*WC_TM(2)*W_PPM(3)
            W_SP_1(152) = WC_TF(1)*WC_TM(2)*W_PPM(4)
            W_SP_1(153) = WC_TF(2)*WC_TM(2)*W_PPM(1)
            W_SP_1(154) = WC_TF(2)*WC_TM(2)*W_PPM(2)
            W_SP_1(155) = WC_TF(2)*WC_TM(2)*W_PPM(3)
            W_SP_1(156) = WC_TF(2)*WC_TM(2)*W_PPM(4)
            W_SP_1(157) = WC_TF(3)*WC_TM(2)*W_PPM(1)
            W_SP_1(158) = WC_TF(3)*WC_TM(2)*W_PPM(2)
            W_SP_1(159) = WC_TF(3)*WC_TM(2)*W_PPM(3)
            W_SP_1(160) = WC_TF(3)*WC_TM(2)*W_PPM(4)
            W_SP_1(161) = WC_TF(1)*WC_TM(3)*W_PPM(1)
            W_SP_1(162) = WC_TF(1)*WC_TM(3)*W_PPM(2)
            W_SP_1(163) = WC_TF(1)*WC_TM(3)*W_PPM(3)
            W_SP_1(164) = WC_TF(1)*WC_TM(3)*W_PPM(4)
            W_SP_1( 17) = WC_TF(2)*WC_TM(3)*W_PPM(1)
            W_SP_1( 18) = WC_TF(2)*WC_TM(3)*W_PPM(2)
            W_SP_1( 19) = WC_TF(2)*WC_TM(3)*W_PPM(3)
            W_SP_1( 20) = WC_TF(2)*WC_TM(3)*W_PPM(4)
            W_SP_1( 14) = WC_TF(3)*WC_TM(3)*W_PPM(1)
            W_SP_1( 13) = WC_TF(3)*WC_TM(3)*W_PPM(2)
            W_SP_1( 15) = WC_TF(3)*WC_TM(3)*W_PPM(3)
            W_SP_1( 16) = WC_TF(3)*WC_TM(3)*W_PPM(4)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 2
            W_SP_1(139) = WC_TM(1)
            W_SP_1(165) = WC_TM(2)
            W_SP_1(166) = WC_TM(3)
            W_SP_1(130) = W_SP_1(130) - WC_TM(1)
            W_SP_1(154) = W_SP_1(154) - WC_TM(2)
            W_SP_1( 18) = W_SP_1( 18) - WC_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 3
            W_SP_1(139) = 0d0
            W_SP_1(165) = 0d0
            W_SP_1(166) = 0d0
            W_SP_1(142) = WC_TM(1)
            W_SP_1(167) = WC_TM(2)
            W_SP_1(168) = WC_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 4
            W_SP_1(142) = 0d0
            W_SP_1(167) = 0d0
            W_SP_1(168) = 0d0
            W_SP_1(145) = WC_TM(1)
            W_SP_1(169) = WC_TM(2)
            W_SP_1(170) = WC_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif

            I_Reg = 5
            W_SP_1(145) = 0d0
            W_SP_1(169) = 0d0
            W_SP_1(170) = 0d0
            W_SP_1(148) = WC_TM(1)
            W_SP_1(171) = WC_TM(2)
            W_SP_1(172) = WC_TM(3)
            if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
               W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
            endif
         else
            write(*,*) 'abnormal interpolation!!'
            stop
         endif

         do i=1,N_SP_FA
            if ((abs(W_SP(i,Ixy,Iz))<1d-8).and.(abs(W_SP(i,Ixy,Iz))>1d-100)) then
                W_SP(i,Ixy,Iz) = 0d0
            endif
         enddo

         if (if_state==1 .and. if_branch==4) then
            MacroXSset(:) = macro_xs_table(ixy,iz)%xs(:, 1)  + &
                            macro_xs_table(ixy,iz)%xs(:, 17) * W_SP( 17, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 18) * W_SP( 18, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 19) * W_SP( 19, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 20) * W_SP( 20, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 14) * W_SP( 14, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 13) * W_SP( 13, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 15) * W_SP( 15, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 16) * W_SP( 16, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 21) * W_SP( 21, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 22) * W_SP( 22, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 23) * W_SP( 23, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 24) * W_SP( 24, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  5) * W_SP(  5, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  6) * W_SP(  6, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  7) * W_SP(  7, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  8) * W_SP(  8, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  2) * W_SP(  2, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  1) * W_SP(  1, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  3) * W_SP(  3, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  4) * W_SP(  4, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  9) * W_SP(  9, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 10) * W_SP( 10, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 11) * W_SP( 11, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 12) * W_SP( 12, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 29) * W_SP( 29, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 30) * W_SP( 30, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 31) * W_SP( 31, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 32) * W_SP( 32, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 26) * W_SP( 26, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 25) * W_SP( 25, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 27) * W_SP( 27, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 28) * W_SP( 28, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 33) * W_SP( 33, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 34) * W_SP( 34, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 35) * W_SP( 35, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 36) * W_SP( 36, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 37) * W_SP( 37, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 38) * W_SP( 38, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 39) * W_SP( 39, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 40) * W_SP( 40, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 41) * W_SP( 41, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 42) * W_SP( 42, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 43) * W_SP( 43, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 44) * W_SP( 44, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 45) * W_SP( 45, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 46) * W_SP( 46, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 47) * W_SP( 47, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 48) * W_SP( 48, Ixy, Iz)
         elseif (if_state==1 .and. if_branch/=4) then
            MacroXSset(:) = macro_xs_table(ixy,iz)%xs(:,  1)  + &
                            macro_xs_table(ixy,iz)%xs(:,  1) * W_SP(  1, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  2) * W_SP(  2, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  3) * W_SP(  3, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  4) * W_SP(  4, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  5) * W_SP(  5, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,  6) * W_SP(  6, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 12) * W_SP( 12, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 17) * W_SP( 17, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 20) * W_SP( 20, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 23) * W_SP( 23, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 26) * W_SP( 26, Ixy, Iz)
         elseif (if_state==2) then
            MacroXSset(:) = macro_xs_table(ixy,iz)%xs(:,  1) - &
                            macro_xs_table(ixy,iz)%xs(:,  1) + &
                            macro_xs_table(ixy,iz)%xs(:, 53) * W_SP( 53, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 54) * W_SP( 54, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 55) * W_SP( 55, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 56) * W_SP( 56, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 50) * W_SP( 50, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 49) * W_SP( 49, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 51) * W_SP( 51, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 52) * W_SP( 52, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 61) * W_SP( 61, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 62) * W_SP( 62, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 63) * W_SP( 63, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 64) * W_SP( 64, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 58) * W_SP( 58, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 57) * W_SP( 57, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 59) * W_SP( 59, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 60) * W_SP( 60, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 69) * W_SP( 69, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 70) * W_SP( 70, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 71) * W_SP( 71, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 72) * W_SP( 72, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 66) * W_SP( 66, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 65) * W_SP( 65, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 67) * W_SP( 67, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 68) * W_SP( 68, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 77) * W_SP( 77, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 78) * W_SP( 78, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 79) * W_SP( 79, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 80) * W_SP( 80, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 74) * W_SP( 74, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 73) * W_SP( 73, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 75) * W_SP( 75, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 76) * W_SP( 76, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 81) * W_SP( 81, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 82) * W_SP( 82, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 83) * W_SP( 83, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 84) * W_SP( 84, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 85) * W_SP( 85, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 86) * W_SP( 86, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 87) * W_SP( 87, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 88) * W_SP( 88, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 89) * W_SP( 89, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 90) * W_SP( 90, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 91) * W_SP( 91, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 92) * W_SP( 92, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 93) * W_SP( 93, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 94) * W_SP( 94, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 95) * W_SP( 95, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 96) * W_SP( 96, Ixy, Iz)
         elseif (if_state==3) then
            MacroXSset(:) = macro_xs_table(ixy,iz)%xs(:,  1) - &
                            macro_xs_table(ixy,iz)%xs(:,  1) + &
                            macro_xs_table(ixy,iz)%xs(:, 77) * W_SP( 77, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 78) * W_SP( 78, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 79) * W_SP( 79, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 80) * W_SP( 80, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 74) * W_SP( 74, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 73) * W_SP( 73, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 75) * W_SP( 75, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 76) * W_SP( 76, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 97) * W_SP( 97, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 98) * W_SP( 98, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 99) * W_SP( 99, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,100) * W_SP(100, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,102) * W_SP(102, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,101) * W_SP(101, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,103) * W_SP(103, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,104) * W_SP(104, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,105) * W_SP(105, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,106) * W_SP(106, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,107) * W_SP(107, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,108) * W_SP(108, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,109) * W_SP(109, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,110) * W_SP(110, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,111) * W_SP(111, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,112) * W_SP(112, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,114) * W_SP(114, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,113) * W_SP(113, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,115) * W_SP(115, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,116) * W_SP(116, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,117) * W_SP(117, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,118) * W_SP(118, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,119) * W_SP(119, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,120) * W_SP(120, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,121) * W_SP(121, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,122) * W_SP(122, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,123) * W_SP(123, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,124) * W_SP(124, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,126) * W_SP(126, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,125) * W_SP(125, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,127) * W_SP(127, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,128) * W_SP(128, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,129) * W_SP(129, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,130) * W_SP(130, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,131) * W_SP(131, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,132) * W_SP(132, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,133) * W_SP(133, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,134) * W_SP(134, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,135) * W_SP(135, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,136) * W_SP(136, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 84) * W_SP( 84, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,137) * W_SP(137, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,138) * W_SP(138, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,139) * W_SP(139, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 88) * W_SP( 88, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,140) * W_SP(140, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,141) * W_SP(141, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,142) * W_SP(142, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 92) * W_SP( 92, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,143) * W_SP(143, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,144) * W_SP(144, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,145) * W_SP(145, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 96) * W_SP( 96, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,146) * W_SP(146, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,147) * W_SP(147, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,148) * W_SP(148, Ixy, Iz)
         elseif (if_state==4) then
            MacroXSset(:) = macro_xs_table(ixy,iz)%xs(:,  1) - &
                            macro_xs_table(ixy,iz)%xs(:,  1) + &
                            macro_xs_table(ixy,iz)%xs(:,126) * W_SP(126, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,125) * W_SP(125, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,127) * W_SP(127, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,128) * W_SP(128, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,129) * W_SP(129, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,130) * W_SP(130, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,131) * W_SP(131, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,132) * W_SP(132, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,133) * W_SP(133, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,134) * W_SP(134, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,135) * W_SP(135, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,136) * W_SP(136, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,150) * W_SP(150, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,149) * W_SP(149, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,151) * W_SP(151, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,152) * W_SP(152, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,153) * W_SP(153, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,154) * W_SP(154, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,155) * W_SP(155, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,156) * W_SP(156, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,157) * W_SP(157, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,158) * W_SP(158, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,159) * W_SP(159, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,160) * W_SP(160, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,161) * W_SP(161, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,162) * W_SP(162, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,163) * W_SP(163, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,164) * W_SP(164, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 17) * W_SP( 17, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 18) * W_SP( 18, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 19) * W_SP( 19, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 20) * W_SP( 20, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 14) * W_SP( 14, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 13) * W_SP( 13, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 15) * W_SP( 15, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:, 16) * W_SP( 16, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,139) * W_SP(139, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,165) * W_SP(165, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,166) * W_SP(166, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,142) * W_SP(142, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,167) * W_SP(167, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,168) * W_SP(168, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,145) * W_SP(145, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,169) * W_SP(169, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,170) * W_SP(170, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,148) * W_SP(148, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,171) * W_SP(171, Ixy, Iz) + &
                            macro_xs_table(ixy,iz)%xs(:,172) * W_SP(172, Ixy, Iz)
         endif
      else ! reflector XS
         if (if_state==1) then
            W_SP(1, Ixy, Iz) = W_TM(4) * W_PPM(1)
            W_SP(2, Ixy, Iz) = W_TM(4) * W_PPM(2)
            W_SP(3, Ixy, Iz) = W_TM(4) * W_PPM(3)
            W_SP(4, Ixy, Iz) = W_TM(4) * W_PPM(4)
            W_SP(5, Ixy, Iz) = W_TM(5) * W_PPM(1)
            W_SP(6, Ixy, Iz) = W_TM(5) * W_PPM(2)
            W_SP(7, Ixy, Iz) = W_TM(5) * W_PPM(3)
            W_SP(8, Ixy, Iz) = W_TM(5) * W_PPM(4)
            MacroXSset(:) =  macro_xs_table(ixy,iz)%xs(:, 1) * W_SP( 1, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:, 2) * W_SP( 2, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:, 3) * W_SP( 3, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:, 4) * W_SP( 4, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:, 5) * W_SP( 5, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:, 6) * W_SP( 6, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:, 7) * W_SP( 7, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:, 8) * W_SP( 8, Ixy, Iz)
         elseif (if_state>1) then
            W_SP( 9, Ixy, Iz) = W_TM(4) * W_PPM(1)
            W_SP(10, Ixy, Iz) = W_TM(4) * W_PPM(2)
            W_SP(11, Ixy, Iz) = W_TM(4) * W_PPM(3)
            W_SP(12, Ixy, Iz) = W_TM(4) * W_PPM(4)
            W_SP(13, Ixy, Iz) = W_TM(5) * W_PPM(1)
            W_SP(14, Ixy, Iz) = W_TM(5) * W_PPM(2)
            W_SP(15, Ixy, Iz) = W_TM(5) * W_PPM(3)
            W_SP(16, Ixy, Iz) = W_TM(5) * W_PPM(4)
            W_SP(17, Ixy, Iz) = W_TM(6) * W_PPM(1)
            W_SP(18, Ixy, Iz) = W_TM(6) * W_PPM(2)
            W_SP(19, Ixy, Iz) = W_TM(6) * W_PPM(3)
            W_SP(20, Ixy, Iz) = W_TM(6) * W_PPM(4)
            W_SP(21, Ixy, Iz) = W_TM(7) * W_PPM(1)
            W_SP(22, Ixy, Iz) = W_TM(7) * W_PPM(2)
            W_SP(23, Ixy, Iz) = W_TM(7) * W_PPM(3)
            W_SP(24, Ixy, Iz) = W_TM(7) * W_PPM(4)
            MacroXSset(:) =  macro_xs_table(ixy,iz)%xs(:, 9) * W_SP( 9, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,10) * W_SP(10, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,11) * W_SP(11, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,12) * W_SP(12, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,13) * W_SP(13, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,14) * W_SP(14, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,15) * W_SP(15, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,16) * W_SP(16, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,17) * W_SP(17, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,18) * W_SP(18, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,19) * W_SP(19, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,20) * W_SP(20, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,21) * W_SP(21, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,22) * W_SP(22, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,23) * W_SP(23, Ixy, Iz) + &
                             macro_xs_table(ixy,iz)%xs(:,24) * W_SP(24, Ixy, Iz)
         endif
      endif

      return
      end subroutine MacroXSsetFB_1D


      subroutine MacroHFFsetFB(Ixy, Iz)
      use Inc_3D, only: T_Fuel, T_Mod
      use Inc_CR, only: Reg_VF
      use Inc_RP, only: I_LP_1N, AxialComp
      use Inc_TH, only: PPM, TM_In
      use inc_th, only: flag_th_chanwise, i_chan_1n, chanwise_t_inlet
      use Inc_XS_File
      use Inc_INP, only: flag_force_tm, force_tm, force_dtm
      use Inc_INP, only: flag_force_tf, force_tf, force_dtf
      use Inc_INP, only: flag_force_tdop, force_tdop, force_dtdop
      use Inc_INP, only: flag_force_ftc, force_ftc, flag_force_mtc, force_mtc
      use Mod_GetNode, only: new_asym_itab
      implicit none
      integer, intent(in) :: Ixy, Iz
      integer :: Ixy_1N, I_PinX
      integer(4) :: i
      real(8) :: tm_in_bak, chg_Tm, chg_tf, tmp_tdop
      real(8) :: tf_bnd, tm_bnd


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [MacroHFFsetFB] in Mod_LiteVer'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (flag_th_chanwise) then
         tm_in_bak=chanwise_t_inlet(i_chan_1n(i_4nto1n(ixy)))+degtok
      else
         tm_in_bak=Tm_in+DegToK
      endif

      if (flag_force_tm) then
         chg_Tm=tm_in_bak+(T_Mod(Ixy,Iz)-tm_in_bak)* force_tm + force_dtm
      else
         chg_Tm=T_Mod(Ixy,Iz)
      endif
      if (flag_force_tf) then
         chg_Tf=t_fuel(ixy,iz)*force_tf+force_dtf
      elseif (flag_force_tdop) then
         tmp_tdop=sqrt(t_fuel(ixy,iz))
         tmp_tdop=tmp_tdop*force_tdop+force_dtdop
         chg_Tf=tmp_tdop*tmp_tdop
      else
         chg_Tf=t_fuel(ixy,iz)
      endif
      if (flag_force_ftc) chg_tf=chg_tf+force_ftc
      if (flag_force_mtc) chg_tm=chg_tm+force_mtc

      tf_bnd=min(tf_max,chg_tf)
      tm_bnd=max(tm_min,min(tm_max,chg_Tm))

      W_SP(:, Ixy, Iz) = 0d0
      W_SP_1 = 0d0

      ! calculate the weight factors for PPM
      W_PPM(1) = (PPM       -Var_PPM(2))*(PPM       -Var_PPM(3))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(1)-Var_PPM(2))*(Var_PPM(1)-Var_PPM(3))*(Var_PPM(1)-Var_PPM(4)) )
      W_PPM(2) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(3))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(2)-Var_PPM(1))*(Var_PPM(2)-Var_PPM(3))*(Var_PPM(2)-Var_PPM(4)) )
      W_PPM(3) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(2))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(3)-Var_PPM(1))*(Var_PPM(3)-Var_PPM(2))*(Var_PPM(3)-Var_PPM(4)) )
      W_PPM(4) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(2))*(PPM       -Var_PPM(3)) / &
               ( (Var_PPM(4)-Var_PPM(1))*(Var_PPM(4)-Var_PPM(2))*(Var_PPM(4)-Var_PPM(3)) )

      if (if_branch == 4) then
         if (T_Mod(Ixy,Iz)>=(Var_TM(1)-10d0).and.(T_Fuel(Ixy,Iz)>=Var_TF(1)-10d0)) then
            !Hot state
            if_state = 1
            ! calculate the weight factors for TF
            W_TF(1) = (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                    ( (sqrt(Var_TF(4))-sqrt(Var_TF(1)))*(sqrt(Var_TF(4))-sqrt(Var_TF(3))) )
            W_TF(2) = (sqrt(tf_bnd)   -sqrt(Var_TF(4)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                    ( (sqrt(Var_TF(1))-sqrt(Var_TF(4)))*(sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
            W_TF(3) = (sqrt(tf_bnd)   -sqrt(Var_TF(4)))*(sqrt(tf_bnd)   -sqrt(Var_TF(1))) / &
                    ( (sqrt(Var_TF(3))-sqrt(Var_TF(4)))*(sqrt(Var_TF(3))-sqrt(Var_TF(1))) )
            W_TF(4) = (sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                    ( (sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
            W_TF(5) = D1 - W_TF(4)
            ! calculate the weight factors for TM
            W_TM(1) = (tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(3)) / &
                    ( (Var_TM(1)-Var_TM(2))*(Var_TM(1)-Var_TM(3)) )
            W_TM(2) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(3)) / &
                    ( (Var_TM(2)-Var_TM(1))*(Var_TM(2)-Var_TM(3)) )
            W_TM(3) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(2)) / &
                    ( (Var_TM(3)-Var_TM(1))*(Var_TM(3)-Var_TM(2)) )
            W_TM(4) = (tm_bnd   -Var_TM(3)) / &
                    ( (Var_TM(1)-Var_TM(3)) )
            W_TM(5) = D1 - W_TM(4)
         else
            if ((T_Fuel(Ixy, Iz)-T_Mod(Ixy, Iz))>1d-1) then
               !Cold and Transient State
               if_state = 2
               W_TF(1) = (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(4))-sqrt(Var_TF(1)))*(sqrt(Var_TF(4))-sqrt(Var_TF(3))) )
               W_TF(2) = (sqrt(tf_bnd)   -sqrt(Var_TF(4)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(1))-sqrt(Var_TF(4)))*(sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
               W_TF(3) = (sqrt(tf_bnd)   -sqrt(Var_TF(4)))*(sqrt(tf_bnd)   -sqrt(Var_TF(1))) / &
                       ( (sqrt(Var_TF(3))-sqrt(Var_TF(4)))*(sqrt(Var_TF(3))-sqrt(Var_TF(1))) )
               W_TF(4) = (sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
               W_TF(5) = D1 - W_TF(4)
               ! calculate the weight factors for TM
               W_TM(1) = (tm_bnd   -Var_TM(6))*(tm_bnd   -Var_TM(8)) / &
                       ( (Var_TM(4)-Var_TM(6))*(Var_TM(4)-Var_TM(8)) )
               W_TM(2) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(8)) / &
                       ( (Var_TM(6)-Var_TM(4))*(Var_TM(6)-Var_TM(8)) )
               W_TM(3) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(6)) / &
                       ( (Var_TM(8)-Var_TM(4))*(Var_TM(8)-Var_TM(6)) )
               W_TM(4) = (tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(7))*(tm_bnd   -Var_TM(1)) / &
                       ( (Var_TM(4)-Var_TM(5))*(Var_TM(4)-Var_TM(7))*(Var_TM(4)-Var_TM(1)) )
               W_TM(5) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(7))*(tm_bnd   -Var_TM(1)) / &
                       ( (Var_TM(5)-Var_TM(4))*(Var_TM(5)-Var_TM(7))*(Var_TM(5)-Var_TM(1)) )
               W_TM(6) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(1)) / &
                       ( (Var_TM(7)-Var_TM(4))*(Var_TM(7)-Var_TM(5))*(Var_TM(7)-Var_TM(1)) )
               W_TM(7) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(7)) / &
                       ( (Var_TM(1)-Var_TM(4))*(Var_TM(1)-Var_TM(5))*(Var_TM(1)-Var_TM(7)) )
            else
               !Cold and Steady State
               if_state = 3
               ! calculate the weight factors for TM
               W_TM(1) = (tm_bnd   -Var_TM(7))*(tm_bnd   -Var_TM(1)) / &
                       ( (Var_TM(4)-Var_TM(7))*(Var_TM(4)-Var_TM(1)) )
               W_TM(2) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(1)) / &
                       ( (Var_TM(7)-Var_TM(4))*(Var_TM(7)-Var_TM(1)) )
               W_TM(3) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(7)) / &
                       ( (Var_TM(1)-Var_TM(4))*(Var_TM(1)-Var_TM(7)) )
               W_TM(4) = (tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(7))*(tm_bnd   -Var_TM(1)) / &
                       ( (Var_TM(4)-Var_TM(5))*(Var_TM(4)-Var_TM(7))*(Var_TM(4)-Var_TM(1)) )
               W_TM(5) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(7))*(tm_bnd   -Var_TM(1)) / &
                       ( (Var_TM(5)-Var_TM(4))*(Var_TM(5)-Var_TM(7))*(Var_TM(5)-Var_TM(1)) )
               W_TM(6) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(1)) / &
                       ( (Var_TM(7)-Var_TM(4))*(Var_TM(7)-Var_TM(5))*(Var_TM(7)-Var_TM(1)) )
               W_TM(7) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(7)) / &
                       ( (Var_TM(1)-Var_TM(4))*(Var_TM(1)-Var_TM(5))*(Var_TM(1)-Var_TM(7)) )
            endif
         endif
      else
         ! calculate the weight factors for TF
         W_TF(1) = (sqrt(tf_bnd)   -sqrt(Var_TF(2)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                 ( (sqrt(Var_TF(1))-sqrt(Var_TF(2)))*(sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
         W_TF(2) = (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                 ( (sqrt(Var_TF(2))-sqrt(Var_TF(1)))*(sqrt(Var_TF(2))-sqrt(Var_TF(3))) )
         W_TF(3) = (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                 ( (sqrt(Var_TF(3))-sqrt(Var_TF(1)))*(sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
         W_TF(4) = (sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                 ( (sqrt(Var_TF(2))-sqrt(Var_TF(3))) )
         W_TF(5) = D1 - W_TF(4)
         ! calculate the weight factors for TM
         W_TM(1) = (tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(1)-Var_TM(2))*(Var_TM(1)-Var_TM(3)) )
         W_TM(2) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(2)-Var_TM(1))*(Var_TM(2)-Var_TM(3)) )
         W_TM(3) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(2)) / &
                 ( (Var_TM(3)-Var_TM(1))*(Var_TM(3)-Var_TM(2)) )
         W_TM(4) = (tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(1)-Var_TM(3)) )
         W_TM(5) = D1 - W_TM(4)
      endif

      Ixy_1N = I_4Nto1N(Ixy)
      I_Tab  = AxialComp( I_LP_1N( Ixy_1N ), Iz )
      I_Type = Type_Tab(I_Tab)
      I_Tab = new_asym_itab(I_Tab,Ixy)

      if (if_state==1) then
         I_Reg = 1
         W_SP_1(1) = W_TF(4) * W_TM(2) + W_PPM(2) * W_TM(2) - W_TM(2) - 1d0
         W_SP_1(2) = W_PPM(1) * W_TM(2)
         W_SP_1(3) = W_PPM(3) * W_TM(2)
         W_SP_1(4) = W_PPM(4) * W_TM(2)
         W_SP_1(5) = W_TF(2) * W_TM(1) + W_PPM(2) * W_TM(1) - W_TM(1)
         W_SP_1(6) = W_TF(4) * W_TM(3) + W_PPM(2) * W_TM(3) - W_TM(3)
         W_SP_1(7) = W_PPM(1) * W_TM(1)
         W_SP_1(8) = W_PPM(3) * W_TM(1)
         W_SP_1(9) = W_PPM(4) * W_TM(1)
         W_SP_1(10) = W_TF(1) * W_TM(1)
         W_SP_1(11) = W_TF(3) * W_TM(1)
         W_SP_1(12) = W_TF(5) * W_TM(2)
         W_SP_1(13) = W_PPM(1) * W_TM(3)
         W_SP_1(14) = W_PPM(3) * W_TM(3)
         W_SP_1(15) = W_PPM(4) * W_TM(3)
         W_SP_1(16) = W_TF(5) * W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 2
         W_SP_1(17) = W_TM(2)
         W_SP_1(18) = W_TM(1)
         W_SP_1(19) = W_TM(3)
         W_SP_1(1) = W_SP_1(1) - W_TM(2)
         W_SP_1(5) = W_SP_1(5) - W_TM(1)
         W_SP_1(6) = W_SP_1(6) - W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 3
         W_SP_1(17) = 0d0
         W_SP_1(18) = 0d0
         W_SP_1(19) = 0d0
         W_SP_1(20) = W_TM(2)
         W_SP_1(21) = W_TM(1)
         W_SP_1(22) = W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 4
         W_SP_1(20) = 0d0
         W_SP_1(21) = 0D0
         W_SP_1(22) = 0D0
         W_SP_1(23) = W_TM(2)
         W_SP_1(24) = W_TM(1)
         W_SP_1(25) = W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 5
         W_SP_1(23) = 0d0
         W_SP_1(24) = 0D0
         W_SP_1(25) = 0D0
         W_SP_1(26) = W_TM(2)
         W_SP_1(27) = W_TM(1)
         W_SP_1(28) = W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
      elseif (if_branch==4 .and. if_state==2) then
         I_Reg = 1
         W_SP_1(29) = W_PPM(1) * W_TM(1)
         W_SP_1(30) = W_PPM(2) * W_TM(1) - W_TM(1)
         W_SP_1(31) = W_PPM(3) * W_TM(1)
         W_SP_1(32) = W_PPM(4) * W_TM(1)
         W_SP_1(33) = W_PPM(1) * W_TM(2)
         W_SP_1(34) = W_PPM(2) * W_TM(2) - W_TM(2)
         W_SP_1(35) = W_PPM(3) * W_TM(2)
         W_SP_1(36) = W_PPM(4) * W_TM(2)
         W_SP_1(37) = W_PPM(1) * W_TM(3)
         W_SP_1(38) = W_PPM(2) * W_TM(3) - W_TM(3)
         W_SP_1(39) = W_PPM(3) * W_TM(3)
         W_SP_1(40) = W_PPM(4) * W_TM(3)
         W_SP_1(41) = W_TF(1)  * W_TM(1)
         W_SP_1(42) = W_TF(2)  * W_TM(1)
         W_SP_1(43) = W_TF(3)  * W_TM(1)
         W_SP_1(44) = W_TF(4)  * W_TM(2)
         W_SP_1(45) = W_TF(5)  * W_TM(2)
         W_SP_1(46) = W_TF(4)  * W_TM(3)
         W_SP_1(47) = W_TF(5)  * W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 2
         W_SP_1(60) = W_TM(1)
         W_SP_1(61) = W_TM(2)
         W_SP_1(62) = W_TM(3)
         W_SP_1(30) = W_SP_1(30) - W_TM(1)
         W_SP_1(34) = W_SP_1(34) - W_TM(2)
         W_SP_1(38) = W_SP_1(38) - W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 3
         W_SP_1(60) = 0d0
         W_SP_1(61) = 0d0
         W_SP_1(62) = 0d0
         W_SP_1(66) = W_TM(1)
         W_SP_1(67) = W_TM(2)
         W_SP_1(68) = W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 4
         W_SP_1(66) = 0d0
         W_SP_1(67) = 0d0
         W_SP_1(68) = 0d0
         W_SP_1(72) = W_TM(1)
         W_SP_1(73) = W_TM(2)
         W_SP_1(74) = W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 5
         W_SP_1(72) = 0d0
         W_SP_1(73) = 0d0
         W_SP_1(74) = 0d0
         W_SP_1(78) = W_TM(1)
         W_SP_1(79) = W_TM(2)
         W_SP_1(80) = W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
      elseif (if_branch==4 .and. if_state==3) then
         I_Reg = 1
         W_SP_1(48) = W_PPM(1) * W_TM(1)
         W_SP_1(49) = W_PPM(2) * W_TM(1)
         W_SP_1(50) = W_PPM(3) * W_TM(1)
         W_SP_1(51) = W_PPM(4) * W_TM(1)
         W_SP_1(52) = W_PPM(1) * W_TM(2)
         W_SP_1(53) = W_PPM(2) * W_TM(2)
         W_SP_1(54) = W_PPM(3) * W_TM(2)
         W_SP_1(55) = W_PPM(4) * W_TM(2)
         W_SP_1(56) = W_PPM(1) * W_TM(3)
         W_SP_1(57) = W_PPM(2) * W_TM(3)
         W_SP_1(58) = W_PPM(3) * W_TM(3)
         W_SP_1(59) = W_PPM(4) * W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 2
         W_SP_1(63) = W_PPM(3) * W_TM(1)
         W_SP_1(64) = W_PPM(3) * W_TM(2)
         W_SP_1(65) = W_PPM(3) * W_TM(3)
         W_SP_1(50) = W_SP_1(50) - W_PPM(3) * W_TM(1)
         W_SP_1(54) = W_SP_1(54) - W_PPM(3) * W_TM(2)
         W_SP_1(58) = W_SP_1(58) - W_PPM(3) * W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 3
         W_SP_1(63) = 0d0
         W_SP_1(64) = 0d0
         W_SP_1(65) = 0d0
         W_SP_1(69) = W_PPM(3) * W_TM(1)
         W_SP_1(70) = W_PPM(3) * W_TM(2)
         W_SP_1(71) = W_PPM(3) * W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 4
         W_SP_1(69) = 0d0
         W_SP_1(70) = 0d0
         W_SP_1(71) = 0d0
         W_SP_1(75) = W_PPM(3) * W_TM(1)
         W_SP_1(76) = W_PPM(3) * W_TM(2)
         W_SP_1(77) = W_PPM(3) * W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
         I_Reg = 5
         W_SP_1(75) = 0d0
         W_SP_1(76) = 0d0
         W_SP_1(77) = 0d0
         W_SP_1(71) = W_PPM(3) * W_TM(1)
         W_SP_1(82) = W_PPM(3) * W_TM(2)
         W_SP_1(83) = W_PPM(3) * W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         end if
      else
         write(*,*) 'abnormal interpolation!!'
         stop
      endif

      do i=1,N_SP_FA
         if ((abs(W_SP(i,Ixy,Iz))<1d-8).and.(abs(W_SP(i,Ixy,Iz))>1d-30)) then
            W_SP(i,Ixy,Iz) = 0d0
         end if
      enddo
      do I_PinX = 1, XS_File_NbyN
         if (if_state==1) then
            HFFset(I_PinX, :) = macro_xs_table(ixy,iz)%HFF(I_PinX, :, 1) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 1) * W_SP( 1, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 2) * W_SP( 2, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 3) * W_SP( 3, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 4) * W_SP( 4, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 5) * W_SP( 5, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 6) * W_SP( 6, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 7) * W_SP( 7, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 8) * W_SP( 8, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 9) * W_SP( 9, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,10) * W_SP(10, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,11) * W_SP(11, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,12) * W_SP(12, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,13) * W_SP(13, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,14) * W_SP(14, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,15) * W_SP(15, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,16) * W_SP(16, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,17) * W_SP(17, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,18) * W_SP(18, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,19) * W_SP(19, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,20) * W_SP(20, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,21) * W_SP(21, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,22) * W_SP(22, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,23) * W_SP(23, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,24) * W_SP(24, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,25) * W_SP(25, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,26) * W_SP(26, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,27) * W_SP(27, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,28) * W_SP(28, Ixy, Iz)
         elseif (if_branch==4 .and. if_state==2) then
            HFFset(I_PinX, :) = macro_xs_table(ixy,iz)%HFF(I_PinX, :, 1) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,29) * W_SP(29, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,30) * W_SP(30, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,31) * W_SP(31, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,32) * W_SP(32, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,33) * W_SP(33, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,34) * W_SP(34, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,35) * W_SP(35, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,36) * W_SP(36, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,37) * W_SP(37, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,38) * W_SP(38, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,39) * W_SP(39, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,40) * W_SP(40, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,41) * W_SP(41, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,42) * W_SP(42, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,43) * W_SP(43, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,44) * W_SP(44, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,45) * W_SP(45, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,46) * W_SP(46, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,47) * W_SP(47, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,60) * W_SP(60, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,61) * W_SP(61, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,62) * W_SP(62, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,66) * W_SP(66, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,67) * W_SP(67, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,68) * W_SP(68, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,72) * W_SP(72, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,73) * W_SP(73, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,74) * W_SP(74, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,78) * W_SP(78, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,79) * W_SP(79, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,80) * W_SP(80, Ixy, Iz)
         elseif (if_branch==4 .and. if_state==3) then
            HFFset(I_PinX, :) = macro_xs_table(ixy,iz)%HFF(I_PinX, :, 1) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,48) * W_SP(48, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,49) * W_SP(49, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,50) * W_SP(50, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,51) * W_SP(51, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,52) * W_SP(52, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,53) * W_SP(53, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,54) * W_SP(54, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,55) * W_SP(55, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,56) * W_SP(56, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,57) * W_SP(57, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,58) * W_SP(58, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,59) * W_SP(59, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,63) * W_SP(63, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,64) * W_SP(64, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,65) * W_SP(65, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,69) * W_SP(69, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,70) * W_SP(70, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,71) * W_SP(71, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,75) * W_SP(75, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,76) * W_SP(76, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,77) * W_SP(77, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,81) * W_SP(81, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,82) * W_SP(82, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,83) * W_SP(83, Ixy, Iz)
         endif
      enddo

      return
      end subroutine MacroHFFsetFB


      subroutine MacroHFFsetFB_1D(Ixy, Iz)
      use Inc_3D, only: T_Fuel, T_Mod
      use Inc_CR, only: Reg_VF
      use Inc_RP, only: I_LP_1N, AxialComp
      use Inc_TH, only: PPM, TM_In
      use inc_th, only: flag_th_chanwise, i_chan_1n, chanwise_t_inlet
      use Inc_XS_File
      use Inc_INP, only: flag_force_tm, force_tm, force_dtm
      use Inc_INP, only: flag_force_tf, force_tf, force_dtf
      use Inc_INP, only: flag_force_tdop, force_tdop, force_dtdop
      use Inc_INP, only: flag_force_ftc, force_ftc, flag_force_mtc, force_mtc
      use Mod_GetNode, only: new_asym_itab
      implicit none
      integer, intent(in) :: Ixy, Iz
      integer :: Ixy_1N, I_PinX
      integer(4) :: i
      real(8) :: tm_in_bak, chg_Tm, chg_tf, tmp_tdop
      real(8) :: tf_bnd, tm_bnd


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [MacroHFFsetFB_1D] in Mod_LiteVer'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (flag_th_chanwise) then
         tm_in_bak=chanwise_t_inlet(i_chan_1n(i_4nto1n(ixy)))+degtok
      else
         tm_in_bak=Tm_in+DegToK
      endif

      if (flag_force_tm) then
         chg_Tm=tm_in_bak+(T_Mod(Ixy,Iz)-tm_in_bak)* force_tm + force_dtm
      else
         chg_Tm=T_Mod(Ixy,Iz)
      endif
      if (flag_force_tf) then
         chg_Tf=t_fuel(ixy,iz)*force_tf+force_dtf
      elseif (flag_force_tdop) then
         tmp_tdop=sqrt(t_fuel(ixy,iz))
         tmp_tdop=tmp_tdop*force_tdop+force_dtdop
         chg_Tf=tmp_tdop*tmp_tdop
      else
         chg_Tf=t_fuel(ixy,iz)
      endif
      if (flag_force_ftc) chg_tf=chg_tf+force_ftc
      if (flag_force_mtc) chg_tm=chg_tm+force_mtc

      tf_bnd=min(tf_max,chg_tf)
      tm_bnd=max(tm_min,min(tm_max,chg_Tm))
      W_SP(:, Ixy, Iz) = 0d0
      W_SP_1 = 0d0

      Ixy_1N = I_4Nto1N(Ixy)
      I_Tab = new_asym_itab(axialcomp(i_lp_1n(ixy_1n),iz),ixy)
      I_Type = Type_Tab(I_Tab)

      if_state=0
      if(if_branch==4) WC_TF=0d0
      if(if_branch==4) WC_TM=0d0
      W_PPM=0d0
      W_TM =0d0
      W_TF =0d0

      ! calculate the weight factors for PPM
      W_PPM(1) = (PPM       -Var_PPM(2))*(PPM       -Var_PPM(3))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(1)-Var_PPM(2))*(Var_PPM(1)-Var_PPM(3))*(Var_PPM(1)-Var_PPM(4)) )
      W_PPM(2) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(3))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(2)-Var_PPM(1))*(Var_PPM(2)-Var_PPM(3))*(Var_PPM(2)-Var_PPM(4)) )
      W_PPM(3) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(2))*(PPM       -Var_PPM(4)) / &
               ( (Var_PPM(3)-Var_PPM(1))*(Var_PPM(3)-Var_PPM(2))*(Var_PPM(3)-Var_PPM(4)) )
      W_PPM(4) = (PPM       -Var_PPM(1))*(PPM       -Var_PPM(2))*(PPM       -Var_PPM(3)) / &
               ( (Var_PPM(4)-Var_PPM(1))*(Var_PPM(4)-Var_PPM(2))*(Var_PPM(4)-Var_PPM(3)) )

      if ((i_type<3) .and. (if_branch==4)) then
         if ((tm_bnd>=Var_ref_TM(1)).and.(tf_bnd>=Var_ref_TF(1))) then
            if_state=1
            W_TF(1) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(2)))*(sqrt(tf_bnd)       -sqrt(Var_ref_TF(3))) / &
                    ( (sqrt(Var_ref_TF(1))-sqrt(Var_ref_TF(2)))*(sqrt(Var_ref_TF(1))-sqrt(Var_ref_TF(3))) )
            W_TF(2) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(1)))*(sqrt(tf_bnd)       -sqrt(Var_ref_TF(3))) / &
                    ( (sqrt(Var_ref_TF(2))-sqrt(Var_ref_TF(1)))*(sqrt(Var_ref_TF(2))-sqrt(Var_ref_TF(3))) )
            W_TF(3) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(1)))*(sqrt(tf_bnd)       -sqrt(Var_ref_TF(2))) / &
                    ( (sqrt(Var_ref_TF(3))-sqrt(Var_ref_TF(1)))*(sqrt(Var_ref_TF(3))-sqrt(Var_ref_TF(2))) )
            W_TF(4) = (sqrt(tf_bnd)       -sqrt(Var_ref_TF(3))) / &
                    ( (sqrt(Var_ref_TF(2))-sqrt(Var_ref_TF(3))) )
            W_TF(5) = 1d0 - W_TF(4)
            ! calculate the weight factors for TM
            W_TM(1) = (tm_bnd       -Var_ref_TM(2))*(tm_bnd       -Var_ref_TM(3)) / &
                    ( (Var_ref_TM(1)-Var_ref_TM(2))*(Var_ref_TM(1)-Var_ref_TM(3)) )
            W_TM(2) = (tm_bnd       -Var_ref_TM(1))*(tm_bnd       -Var_ref_TM(3)) / &
                    ( (Var_ref_TM(2)-Var_ref_TM(1))*(Var_ref_TM(2)-Var_ref_TM(3)) )
            W_TM(3) = (tm_bnd       -Var_ref_TM(1))*(tm_bnd       -Var_ref_TM(2)) / &
                    ( (Var_ref_TM(3)-Var_ref_TM(1))*(Var_ref_TM(3)-Var_ref_TM(2)) )
            W_TM(4) = (tm_bnd       -Var_ref_TM(3)) / &
                    ( (Var_ref_TM(1)-Var_ref_TM(3)) )
            W_TM(5) = 1d0 - W_TM(4)
         elseif (((tm_bnd<=Var_TM(9)).and.(tf_bnd<=Var_TF(2))).and.((tm_bnd>=Var_TM(1)-5d0).and.(tf_bnd>=Var_TF(1)-5d0))) then
            if(tm_bnd<=400d0) then
               if_state=2
               WC_TF(1)= (sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                       ( (sqrt(Var_TF(1))-sqrt(Var_TF(2))) )
               WC_TF(2)= (sqrt(tf_bnd)   -sqrt(Var_TF(1))) / &
                       ( (sqrt(Var_TF(2))-sqrt(Var_TF(1))) )
               WC_TM(1) = (tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(3))*(tm_bnd   -Var_TM(4)) / &
                        ( (Var_TM(1)-Var_TM(2))*(Var_TM(1)-Var_TM(3))*(Var_TM(1)-Var_TM(4)) )
               WC_TM(2) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(3))*(tm_bnd   -Var_TM(4)) / &
                        ( (Var_TM(2)-Var_TM(1))*(Var_TM(2)-Var_TM(3))*(Var_TM(2)-Var_TM(4)) )
               WC_TM(3) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(4)) / &
                        ( (Var_TM(3)-Var_TM(1))*(Var_TM(3)-Var_TM(2))*(Var_TM(3)-Var_TM(4)) )
               WC_TM(4) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(3)) / &
                        ( (Var_TM(4)-Var_TM(1))*(Var_TM(4)-Var_TM(2))*(Var_TM(4)-Var_TM(3)) )
            elseif(tm_bnd>400d0 .and. tm_bnd<=500d0) then
               if_state=3
               WC_TF(1)= (sqrt(tf_bnd)   -sqrt(Var_TF(2)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(1))-sqrt(Var_TF(2)))*(sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
               WC_TF(2)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(2))-sqrt(Var_TF(1)))*(sqrt(Var_TF(2))-sqrt(Var_TF(3))) )
               WC_TF(3)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                       ( (sqrt(Var_TF(3))-sqrt(Var_TF(1)))*(sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
               WC_TM(1) = (tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(6))*(tm_bnd   -Var_TM(7)) / &
                        ( (Var_TM(4)-Var_TM(5))*(Var_TM(4)-Var_TM(6))*(Var_TM(4)-Var_TM(7)) )
               WC_TM(2) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(6))*(tm_bnd   -Var_TM(7)) / &
                        ( (Var_TM(5)-Var_TM(4))*(Var_TM(5)-Var_TM(6))*(Var_TM(5)-Var_TM(7)) )
               WC_TM(3) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(7)) / &
                        ( (Var_TM(6)-Var_TM(4))*(Var_TM(6)-Var_TM(5))*(Var_TM(6)-Var_TM(7)) )
               WC_TM(4) = (tm_bnd   -Var_TM(4))*(tm_bnd   -Var_TM(5))*(tm_bnd   -Var_TM(6)) / &
                        ( (Var_TM(7)-Var_TM(4))*(Var_TM(7)-Var_TM(5))*(Var_TM(7)-Var_TM(6)) )
            elseif(tm_bnd>500d0) then
               if_state=4
               WC_TF(1)= (sqrt(tf_bnd)   -sqrt(Var_TF(2)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(1))-sqrt(Var_TF(2)))*(sqrt(Var_TF(1))-sqrt(Var_TF(3))) )
               WC_TF(2)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(3))) / &
                       ( (sqrt(Var_TF(2))-sqrt(Var_TF(1)))*(sqrt(Var_TF(2))-sqrt(Var_TF(3))) )
               WC_TF(3)= (sqrt(tf_bnd)   -sqrt(Var_TF(1)))*(sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                       ( (sqrt(Var_TF(3))-sqrt(Var_TF(1)))*(sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
               WC_TM(1) = (tm_bnd   -Var_TM(8))*(tm_bnd   -Var_TM(9)) / &
                        ( (Var_TM(7)-Var_TM(8))*(Var_TM(7)-Var_TM(9)) )
               WC_TM(2) = (tm_bnd   -Var_TM(7))*(tm_bnd   -Var_TM(9)) / &
                        ( (Var_TM(8)-Var_TM(7))*(Var_TM(8)-Var_TM(9)) )
               WC_TM(3) = (tm_bnd   -Var_TM(7))*(tm_bnd   -Var_TM(8)) / &
                        ( (Var_TM(9)-Var_TM(7))*(Var_TM(9)-Var_TM(8)) )
            endif
         else
            write(*,*) 'the cross section of this range is extrapolated! XSsetFB: TMO= ',tm_bnd, ' TFU= ', tf_bnd
            write(*,*) 'at node position (Ixy, Iz) ', Ixy, Iz
            write(*,*) 'i_type:   ', i_type
            write(*,*) 'if_state: ', if_state
            write(*,*) 'tm=       ', var_tm(:)
            write(*,*) 'tf=       ', var_tf(:)
            stop
         endif
      elseif ((i_type>=3) .and. (if_branch==4)) then
         if (tm_bnd>Var_ref_TM(1)) then
           !Hot state
           if_state = 1
           W_TM(4) = (tm_bnd   -Var_ref_TM(3)) / &
                   ( (Var_ref_TM(1)-Var_ref_TM(3)) )
           W_TM(5) = 1d0 - W_TM(4)
         elseif (tm_bnd<=Var_ref_TM(1)) then
           !Cold state
           if_state = 2
           ! calculate the weight factors for TM
           W_TM(4) = (tm_bnd       -Var_ref_TM(5))*(tm_bnd       -Var_ref_TM(7))*(tm_bnd       -Var_ref_TM(1)) / &
                   ( (Var_ref_TM(4)-Var_ref_TM(5))*(Var_ref_TM(4)-Var_ref_TM(7))*(Var_ref_TM(4)-Var_ref_TM(1)) )
           W_TM(5) = (tm_bnd       -Var_ref_TM(4))*(tm_bnd       -Var_ref_TM(7))*(tm_bnd       -Var_ref_TM(1)) / &
                   ( (Var_ref_TM(5)-Var_ref_TM(4))*(Var_ref_TM(5)-Var_ref_TM(7))*(Var_ref_TM(5)-Var_ref_TM(1)) )
           W_TM(6) = (tm_bnd       -Var_ref_TM(4))*(tm_bnd       -Var_ref_TM(5))*(tm_bnd       -Var_ref_TM(1)) / &
                   ( (Var_ref_TM(7)-Var_ref_TM(4))*(Var_ref_TM(7)-Var_ref_TM(5))*(Var_ref_TM(7)-Var_ref_TM(1)) )
           W_TM(7) = (tm_bnd       -Var_ref_TM(4))*(tm_bnd       -Var_ref_TM(5))*(tm_bnd       -Var_ref_TM(7)) / &
                   ( (Var_ref_TM(1)-Var_ref_TM(4))*(Var_ref_TM(1)-Var_ref_TM(5))*(Var_ref_TM(1)-Var_ref_TM(7)) )
         endif
      else
         if_state=1
         ! calculate the weight factors for TF
         W_TF(3) = (sqrt(tf_bnd)   -sqrt(Var_TF(2))) / &
                 ( (sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
         W_TF(2) = (sqrt(Var_TF(3))-sqrt(tf_bnd)   ) / &
                 ( (sqrt(Var_TF(3))-sqrt(Var_TF(2))) )
         ! calculate the weight factors for TM
         W_TM(1) = (tm_bnd   -Var_TM(2))*(tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(1)-Var_TM(2))*(Var_TM(1)-Var_TM(3)) )
         W_TM(2) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(2)-Var_TM(1))*(Var_TM(2)-Var_TM(3)) )
         W_TM(3) = (tm_bnd   -Var_TM(1))*(tm_bnd   -Var_TM(2)) / &
                 ( (Var_TM(3)-Var_TM(1))*(Var_TM(3)-Var_TM(2)) )
         W_TM(4) = (tm_bnd   -Var_TM(3)) / &
                 ( (Var_TM(1)-Var_TM(3)) )
         W_TM(5) = 1d0 - W_TM(4)
      endif

      if(if_state==1 .and. if_branch==4) then
         I_Reg = 1
         W_SP_1(17) = W_TF(1)*W_TM(1)*W_PPM(1)
         W_SP_1(18) = W_TF(1)*W_TM(1)*W_PPM(2)
         W_SP_1(19) = W_TF(1)*W_TM(1)*W_PPM(3)
         W_SP_1(20) = W_TF(1)*W_TM(1)*W_PPM(4)
         W_SP_1(14) = W_TF(2)*W_TM(1)*W_PPM(1)
         W_SP_1(13) = W_TF(2)*W_TM(1)*W_PPM(2)
         W_SP_1(15) = W_TF(2)*W_TM(1)*W_PPM(3)
         W_SP_1(16) = W_TF(2)*W_TM(1)*W_PPM(4)
         W_SP_1(21) = W_TF(3)*W_TM(1)*W_PPM(1)
         W_SP_1(22) = W_TF(3)*W_TM(1)*W_PPM(2)
         W_SP_1(23) = W_TF(3)*W_TM(1)*W_PPM(3)
         W_SP_1(24) = W_TF(3)*W_TM(1)*W_PPM(4)
         W_SP_1( 5) = W_TF(1)*W_TM(2)*W_PPM(1)
         W_SP_1( 6) = W_TF(1)*W_TM(2)*W_PPM(2)
         W_SP_1( 7) = W_TF(1)*W_TM(2)*W_PPM(3)
         W_SP_1( 8) = W_TF(1)*W_TM(2)*W_PPM(4)
         W_SP_1( 2) = W_TF(2)*W_TM(2)*W_PPM(1)
         W_SP_1( 1) = W_TF(2)*W_TM(2)*W_PPM(2) -1d0
         W_SP_1( 3) = W_TF(2)*W_TM(2)*W_PPM(3)
         W_SP_1( 4) = W_TF(2)*W_TM(2)*W_PPM(4)
         W_SP_1( 9) = W_TF(3)*W_TM(2)*W_PPM(1)
         W_SP_1(10) = W_TF(3)*W_TM(2)*W_PPM(2)
         W_SP_1(11) = W_TF(3)*W_TM(2)*W_PPM(3)
         W_SP_1(12) = W_TF(3)*W_TM(2)*W_PPM(4)
         W_SP_1(29) = W_TF(1)*W_TM(3)*W_PPM(1)
         W_SP_1(30) = W_TF(1)*W_TM(3)*W_PPM(2)
         W_SP_1(31) = W_TF(1)*W_TM(3)*W_PPM(3)
         W_SP_1(32) = W_TF(1)*W_TM(3)*W_PPM(4)
         W_SP_1(26) = W_TF(2)*W_TM(3)*W_PPM(1)
         W_SP_1(25) = W_TF(2)*W_TM(3)*W_PPM(2)
         W_SP_1(27) = W_TF(2)*W_TM(3)*W_PPM(3)
         W_SP_1(28) = W_TF(2)*W_TM(3)*W_PPM(4)
         W_SP_1(33) = W_TF(3)*W_TM(3)*W_PPM(1)
         W_SP_1(34) = W_TF(3)*W_TM(3)*W_PPM(2)
         W_SP_1(35) = W_TF(3)*W_TM(3)*W_PPM(3)
         W_SP_1(36) = W_TF(3)*W_TM(3)*W_PPM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 2
         W_SP_1(37) = W_TM(2)
         W_SP_1(38) = W_TM(1)
         W_SP_1(39) = W_TM(3)
         W_SP_1( 1) = W_SP_1( 1) - W_TM(2)
         W_SP_1(13) = W_SP_1(13) - W_TM(1)
         W_SP_1(25) = W_SP_1(25) - W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 3
         W_SP_1(37) = 0d0
         W_SP_1(38) = 0d0
         W_SP_1(39) = 0d0
         W_SP_1(40) = W_TM(2)
         W_SP_1(41) = W_TM(1)
         W_SP_1(42) = W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 4
         W_SP_1(40) = 0d0
         W_SP_1(41) = 0d0
         W_SP_1(42) = 0d0
         W_SP_1(43) = W_TM(2)
         W_SP_1(44) = W_TM(1)
         W_SP_1(45) = W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 5
         W_SP_1(43) = 0d0
         W_SP_1(44) = 0d0
         W_SP_1(45) = 0d0
         W_SP_1(46) = W_TM(2)
         W_SP_1(47) = W_TM(1)
         W_SP_1(48) = W_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

      elseif(if_state==1 .and. if_branch/=4) then
         I_Reg = 1
         W_SP_1(1)  = W_PPM(2) + W_TF(2) + W_TM(2) - 3d0
         W_SP_1(2)  = W_PPM(1)
         W_SP_1(3)  = W_PPM(3)
         W_SP_1(4)  = W_PPM(4)
         W_SP_1(5)  = W_TM(1)
         W_SP_1(6)  = W_TM(3)
         W_SP_1(12) = W_TF(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 2
         W_SP_1(17) = 1d0
         W_SP_1(1) = W_SP_1(1) - 1d0
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 3
         W_SP_1(17) = 0d0
         W_SP_1(20) = 1d0
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 4
         W_SP_1(20) = 0d0
         W_SP_1(23) = 1d0
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 5
         W_SP_1(23) = 0d0
         W_SP_1(26) = 1d0
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif
      elseif(if_state==2) then
         I_Reg = 1
         W_SP_1( 53) = WC_TF(1)*WC_TM(1)*W_PPM(1)
         W_SP_1( 54) = WC_TF(1)*WC_TM(1)*W_PPM(2)
         W_SP_1( 55) = WC_TF(1)*WC_TM(1)*W_PPM(3)
         W_SP_1( 56) = WC_TF(1)*WC_TM(1)*W_PPM(4)
         W_SP_1( 50) = WC_TF(2)*WC_TM(1)*W_PPM(1)
         W_SP_1( 49) = WC_TF(2)*WC_TM(1)*W_PPM(2)
         W_SP_1( 51) = WC_TF(2)*WC_TM(1)*W_PPM(3)
         W_SP_1( 52) = WC_TF(2)*WC_TM(1)*W_PPM(4)
         W_SP_1( 61) = WC_TF(1)*WC_TM(2)*W_PPM(1)
         W_SP_1( 62) = WC_TF(1)*WC_TM(2)*W_PPM(2)
         W_SP_1( 63) = WC_TF(1)*WC_TM(2)*W_PPM(3)
         W_SP_1( 64) = WC_TF(1)*WC_TM(2)*W_PPM(4)
         W_SP_1( 58) = WC_TF(2)*WC_TM(2)*W_PPM(1)
         W_SP_1( 57) = WC_TF(2)*WC_TM(2)*W_PPM(2)
         W_SP_1( 59) = WC_TF(2)*WC_TM(2)*W_PPM(3)
         W_SP_1( 60) = WC_TF(2)*WC_TM(2)*W_PPM(4)
         W_SP_1( 69) = WC_TF(1)*WC_TM(3)*W_PPM(1)
         W_SP_1( 70) = WC_TF(1)*WC_TM(3)*W_PPM(2)
         W_SP_1( 71) = WC_TF(1)*WC_TM(3)*W_PPM(3)
         W_SP_1( 72) = WC_TF(1)*WC_TM(3)*W_PPM(4)
         W_SP_1( 66) = WC_TF(2)*WC_TM(3)*W_PPM(1)
         W_SP_1( 65) = WC_TF(2)*WC_TM(3)*W_PPM(2)
         W_SP_1( 67) = WC_TF(2)*WC_TM(3)*W_PPM(3)
         W_SP_1( 68) = WC_TF(2)*WC_TM(3)*W_PPM(4)
         W_SP_1( 77) = WC_TF(1)*WC_TM(4)*W_PPM(1)
         W_SP_1( 78) = WC_TF(1)*WC_TM(4)*W_PPM(2)
         W_SP_1( 79) = WC_TF(1)*WC_TM(4)*W_PPM(3)
         W_SP_1( 80) = WC_TF(1)*WC_TM(4)*W_PPM(4)
         W_SP_1( 74) = WC_TF(2)*WC_TM(4)*W_PPM(1)
         W_SP_1( 73) = WC_TF(2)*WC_TM(4)*W_PPM(2)
         W_SP_1( 75) = WC_TF(2)*WC_TM(4)*W_PPM(3)
         W_SP_1( 76) = WC_TF(2)*WC_TM(4)*W_PPM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 2
         W_SP_1(81) = WC_TM(1)
         W_SP_1(82) = WC_TM(2)
         W_SP_1(83) = WC_TM(3)
         W_SP_1(84) = WC_TM(4)
         W_SP_1(49) = W_SP_1(49) - WC_TM(1)
         W_SP_1(57) = W_SP_1(57) - WC_TM(2)
         W_SP_1(65) = W_SP_1(65) - WC_TM(3)
         W_SP_1(73) = W_SP_1(73) - WC_TM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 3
         W_SP_1(81) = 0d0
         W_SP_1(82) = 0d0
         W_SP_1(83) = 0d0
         W_SP_1(84) = 0d0
         W_SP_1(85) = WC_TM(1)
         W_SP_1(86) = WC_TM(2)
         W_SP_1(87) = WC_TM(3)
         W_SP_1(88) = WC_TM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 4
         W_SP_1(85) = 0d0
         W_SP_1(86) = 0d0
         W_SP_1(87) = 0d0
         W_SP_1(88) = 0d0
         W_SP_1(89) = WC_TM(1)
         W_SP_1(90) = WC_TM(2)
         W_SP_1(91) = WC_TM(3)
         W_SP_1(92) = WC_TM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 5
         W_SP_1(89) = 0d0
         W_SP_1(90) = 0d0
         W_SP_1(91) = 0d0
         W_SP_1(92) = 0d0
         W_SP_1(93) = WC_TM(1)
         W_SP_1(94) = WC_TM(2)
         W_SP_1(95) = WC_TM(3)
         W_SP_1(96) = WC_TM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif
      elseif(if_state==3) then
         I_Reg = 1
         W_SP_1( 77) = WC_TF(1)*WC_TM(1)*W_PPM(1)
         W_SP_1( 78) = WC_TF(1)*WC_TM(1)*W_PPM(2)
         W_SP_1( 79) = WC_TF(1)*WC_TM(1)*W_PPM(3)
         W_SP_1( 80) = WC_TF(1)*WC_TM(1)*W_PPM(4)
         W_SP_1( 74) = WC_TF(2)*WC_TM(1)*W_PPM(1)
         W_SP_1( 73) = WC_TF(2)*WC_TM(1)*W_PPM(2)
         W_SP_1( 75) = WC_TF(2)*WC_TM(1)*W_PPM(3)
         W_SP_1( 76) = WC_TF(2)*WC_TM(1)*W_PPM(4)
         W_SP_1( 97) = WC_TF(3)*WC_TM(1)*W_PPM(1)
         W_SP_1( 98) = WC_TF(3)*WC_TM(1)*W_PPM(2)
         W_SP_1( 99) = WC_TF(3)*WC_TM(1)*W_PPM(3)
         W_SP_1(100) = WC_TF(3)*WC_TM(1)*W_PPM(4)
         W_SP_1(102) = WC_TF(1)*WC_TM(2)*W_PPM(1)
         W_SP_1(101) = WC_TF(1)*WC_TM(2)*W_PPM(2)
         W_SP_1(103) = WC_TF(1)*WC_TM(2)*W_PPM(3)
         W_SP_1(104) = WC_TF(1)*WC_TM(2)*W_PPM(4)
         W_SP_1(105) = WC_TF(2)*WC_TM(2)*W_PPM(1)
         W_SP_1(106) = WC_TF(2)*WC_TM(2)*W_PPM(2)
         W_SP_1(107) = WC_TF(2)*WC_TM(2)*W_PPM(3)
         W_SP_1(108) = WC_TF(2)*WC_TM(2)*W_PPM(4)
         W_SP_1(109) = WC_TF(3)*WC_TM(2)*W_PPM(1)
         W_SP_1(110) = WC_TF(3)*WC_TM(2)*W_PPM(2)
         W_SP_1(111) = WC_TF(3)*WC_TM(2)*W_PPM(3)
         W_SP_1(112) = WC_TF(3)*WC_TM(2)*W_PPM(4)
         W_SP_1(114) = WC_TF(1)*WC_TM(3)*W_PPM(1)
         W_SP_1(113) = WC_TF(1)*WC_TM(3)*W_PPM(2)
         W_SP_1(115) = WC_TF(1)*WC_TM(3)*W_PPM(3)
         W_SP_1(116) = WC_TF(1)*WC_TM(3)*W_PPM(4)
         W_SP_1(117) = WC_TF(2)*WC_TM(3)*W_PPM(1)
         W_SP_1(118) = WC_TF(2)*WC_TM(3)*W_PPM(2)
         W_SP_1(119) = WC_TF(2)*WC_TM(3)*W_PPM(3)
         W_SP_1(120) = WC_TF(2)*WC_TM(3)*W_PPM(4)
         W_SP_1(121) = WC_TF(3)*WC_TM(3)*W_PPM(1)
         W_SP_1(122) = WC_TF(3)*WC_TM(3)*W_PPM(2)
         W_SP_1(123) = WC_TF(3)*WC_TM(3)*W_PPM(3)
         W_SP_1(124) = WC_TF(3)*WC_TM(3)*W_PPM(4)
         W_SP_1(126) = WC_TF(1)*WC_TM(4)*W_PPM(1)
         W_SP_1(125) = WC_TF(1)*WC_TM(4)*W_PPM(2)
         W_SP_1(127) = WC_TF(1)*WC_TM(4)*W_PPM(3)
         W_SP_1(128) = WC_TF(1)*WC_TM(4)*W_PPM(4)
         W_SP_1(129) = WC_TF(2)*WC_TM(4)*W_PPM(1)
         W_SP_1(130) = WC_TF(2)*WC_TM(4)*W_PPM(2)
         W_SP_1(131) = WC_TF(2)*WC_TM(4)*W_PPM(3)
         W_SP_1(132) = WC_TF(2)*WC_TM(4)*W_PPM(4)
         W_SP_1(133) = WC_TF(3)*WC_TM(4)*W_PPM(1)
         W_SP_1(134) = WC_TF(3)*WC_TM(4)*W_PPM(2)
         W_SP_1(135) = WC_TF(3)*WC_TM(4)*W_PPM(3)
         W_SP_1(136) = WC_TF(3)*WC_TM(4)*W_PPM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 2
         W_SP_1( 84) = WC_TM(1)
         W_SP_1(137) = WC_TM(2)
         W_SP_1(138) = WC_TM(3)
         W_SP_1(139) = WC_TM(4)
         W_SP_1( 73) = W_SP_1( 73) - WC_TM(1)
         W_SP_1(106) = W_SP_1(106) - WC_TM(2)
         W_SP_1(118) = W_SP_1(118) - WC_TM(3)
         W_SP_1(130) = W_SP_1(130) - WC_TM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 3
         W_SP_1( 84) = 0d0
         W_SP_1(137) = 0d0
         W_SP_1(138) = 0d0
         W_SP_1(139) = 0d0
         W_SP_1( 88) = WC_TM(1)
         W_SP_1(140) = WC_TM(2)
         W_SP_1(141) = WC_TM(3)
         W_SP_1(142) = WC_TM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 4
         W_SP_1( 88) = 0d0
         W_SP_1(140) = 0d0
         W_SP_1(141) = 0d0
         W_SP_1(142) = 0d0
         W_SP_1( 92) = WC_TM(1)
         W_SP_1(143) = WC_TM(2)
         W_SP_1(144) = WC_TM(3)
         W_SP_1(145) = WC_TM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 5
         W_SP_1( 92) = 0d0
         W_SP_1(143) = 0d0
         W_SP_1(144) = 0d0
         W_SP_1(145) = 0d0
         W_SP_1( 96) = WC_TM(1)
         W_SP_1(146) = WC_TM(2)
         W_SP_1(147) = WC_TM(3)
         W_SP_1(148) = WC_TM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

      elseif(if_state==4) then
         I_Reg = 1
         W_SP_1(126) = WC_TF(1)*WC_TM(1)*W_PPM(1)
         W_SP_1(125) = WC_TF(1)*WC_TM(1)*W_PPM(2)
         W_SP_1(127) = WC_TF(1)*WC_TM(1)*W_PPM(3)
         W_SP_1(128) = WC_TF(1)*WC_TM(1)*W_PPM(4)
         W_SP_1(129) = WC_TF(2)*WC_TM(1)*W_PPM(1)
         W_SP_1(130) = WC_TF(2)*WC_TM(1)*W_PPM(2)
         W_SP_1(131) = WC_TF(2)*WC_TM(1)*W_PPM(3)
         W_SP_1(132) = WC_TF(2)*WC_TM(1)*W_PPM(4)
         W_SP_1(133) = WC_TF(3)*WC_TM(1)*W_PPM(1)
         W_SP_1(134) = WC_TF(3)*WC_TM(1)*W_PPM(2)
         W_SP_1(135) = WC_TF(3)*WC_TM(1)*W_PPM(3)
         W_SP_1(136) = WC_TF(3)*WC_TM(1)*W_PPM(4)
         W_SP_1(150) = WC_TF(1)*WC_TM(2)*W_PPM(1)
         W_SP_1(149) = WC_TF(1)*WC_TM(2)*W_PPM(2)
         W_SP_1(151) = WC_TF(1)*WC_TM(2)*W_PPM(3)
         W_SP_1(152) = WC_TF(1)*WC_TM(2)*W_PPM(4)
         W_SP_1(153) = WC_TF(2)*WC_TM(2)*W_PPM(1)
         W_SP_1(154) = WC_TF(2)*WC_TM(2)*W_PPM(2)
         W_SP_1(155) = WC_TF(2)*WC_TM(2)*W_PPM(3)
         W_SP_1(156) = WC_TF(2)*WC_TM(2)*W_PPM(4)
         W_SP_1(157) = WC_TF(3)*WC_TM(2)*W_PPM(1)
         W_SP_1(158) = WC_TF(3)*WC_TM(2)*W_PPM(2)
         W_SP_1(159) = WC_TF(3)*WC_TM(2)*W_PPM(3)
         W_SP_1(160) = WC_TF(3)*WC_TM(2)*W_PPM(4)
         W_SP_1(161) = WC_TF(1)*WC_TM(3)*W_PPM(1)
         W_SP_1(162) = WC_TF(1)*WC_TM(3)*W_PPM(2)
         W_SP_1(163) = WC_TF(1)*WC_TM(3)*W_PPM(3)
         W_SP_1(164) = WC_TF(1)*WC_TM(3)*W_PPM(4)
         W_SP_1( 17) = WC_TF(2)*WC_TM(3)*W_PPM(1)
         W_SP_1( 18) = WC_TF(2)*WC_TM(3)*W_PPM(2)
         W_SP_1( 19) = WC_TF(2)*WC_TM(3)*W_PPM(3)
         W_SP_1( 20) = WC_TF(2)*WC_TM(3)*W_PPM(4)
         W_SP_1( 14) = WC_TF(3)*WC_TM(3)*W_PPM(1)
         W_SP_1( 13) = WC_TF(3)*WC_TM(3)*W_PPM(2)
         W_SP_1( 15) = WC_TF(3)*WC_TM(3)*W_PPM(3)
         W_SP_1( 16) = WC_TF(3)*WC_TM(3)*W_PPM(4)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 2
         W_SP_1(139) = WC_TM(1)
         W_SP_1(165) = WC_TM(2)
         W_SP_1(166) = WC_TM(3)
         W_SP_1(130) = W_SP_1(130) - WC_TM(1)
         W_SP_1(154) = W_SP_1(154) - WC_TM(2)
         W_SP_1( 18) = W_SP_1( 18) - WC_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 3
         W_SP_1(139) = 0d0
         W_SP_1(165) = 0d0
         W_SP_1(166) = 0d0
         W_SP_1(142) = WC_TM(1)
         W_SP_1(167) = WC_TM(2)
         W_SP_1(168) = WC_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 4
         W_SP_1(142) = 0d0
         W_SP_1(167) = 0d0
         W_SP_1(168) = 0d0
         W_SP_1(145) = WC_TM(1)
         W_SP_1(169) = WC_TM(2)
         W_SP_1(170) = WC_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif

         I_Reg = 5
         W_SP_1(145) = 0d0
         W_SP_1(169) = 0d0
         W_SP_1(170) = 0d0
         W_SP_1(148) = WC_TM(1)
         W_SP_1(171) = WC_TM(2)
         W_SP_1(172) = WC_TM(3)
         if ( abs(Reg_VF(I_Reg)) > 1d-30 ) then
            W_SP(:, Ixy, Iz) = W_SP(:, Ixy, Iz) + W_SP_1 * Reg_VF(I_Reg)
         endif
      else
         write(*,*) 'abnormal interpolation!!'
         stop
      endif

      do i=1,N_SP_FA
         if ((abs(W_SP(i,Ixy,Iz))<1d-8).and.(abs(W_SP(i,Ixy,Iz))>1d-100)) then
             W_SP(i,Ixy,Iz) = 0d0
         endif
      enddo

      do I_PinX = 1, XS_File_NbyN
         if (if_state==1 .and. if_branch==4) then
            HFFset(I_PinX, :) = macro_xs_table(ixy,iz)%HFF(I_PinX, :,  1) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 17) * W_SP( 17, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 18) * W_SP( 18, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 19) * W_SP( 19, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 20) * W_SP( 20, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 14) * W_SP( 14, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 13) * W_SP( 13, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 15) * W_SP( 15, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 16) * W_SP( 16, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 21) * W_SP( 21, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 22) * W_SP( 22, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 23) * W_SP( 23, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 24) * W_SP( 24, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  5) * W_SP(  5, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  6) * W_SP(  6, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  7) * W_SP(  7, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  8) * W_SP(  8, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  2) * W_SP(  2, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  1) * W_SP(  1, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  3) * W_SP(  3, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  4) * W_SP(  4, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  9) * W_SP(  9, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 10) * W_SP( 10, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 11) * W_SP( 11, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 12) * W_SP( 12, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 29) * W_SP( 29, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 30) * W_SP( 30, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 31) * W_SP( 31, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 32) * W_SP( 32, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 26) * W_SP( 26, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 25) * W_SP( 25, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 27) * W_SP( 27, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 28) * W_SP( 28, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 33) * W_SP( 33, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 34) * W_SP( 34, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 35) * W_SP( 35, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 36) * W_SP( 36, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 37) * W_SP( 37, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 38) * W_SP( 38, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 39) * W_SP( 39, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 40) * W_SP( 40, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 41) * W_SP( 41, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 42) * W_SP( 42, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 43) * W_SP( 43, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 44) * W_SP( 44, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 45) * W_SP( 45, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 46) * W_SP( 46, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 47) * W_SP( 47, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 48) * W_SP( 48, Ixy, Iz)
         elseif (if_state==1 .and. if_branch/=4) then
            HFFset(I_PinX, :) = macro_xs_table(ixy,iz)%HFF(I_PinX, :,  1) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  1) * W_SP(  1, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  2) * W_SP(  2, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  3) * W_SP(  3, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  4) * W_SP(  4, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  5) * W_SP(  5, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  6) * W_SP(  6, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 12) * W_SP( 12, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 17) * W_SP( 17, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 20) * W_SP( 20, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 23) * W_SP( 23, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 26) * W_SP( 26, Ixy, Iz)
         elseif (if_state==2) then
            HFFset(I_PinX, :) = macro_xs_table(ixy,iz)%HFF(I_PinX, :,  1) - &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  1) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 53) * W_SP( 53, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 54) * W_SP( 54, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 55) * W_SP( 55, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 56) * W_SP( 56, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 50) * W_SP( 50, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 49) * W_SP( 49, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 51) * W_SP( 51, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 52) * W_SP( 52, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 61) * W_SP( 61, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 62) * W_SP( 62, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 63) * W_SP( 63, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 64) * W_SP( 64, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 58) * W_SP( 58, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 57) * W_SP( 57, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 59) * W_SP( 59, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 60) * W_SP( 60, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 69) * W_SP( 69, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 70) * W_SP( 70, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 71) * W_SP( 71, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 72) * W_SP( 72, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 66) * W_SP( 66, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 65) * W_SP( 65, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 67) * W_SP( 67, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 68) * W_SP( 68, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 77) * W_SP( 77, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 78) * W_SP( 78, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 79) * W_SP( 79, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 80) * W_SP( 80, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 74) * W_SP( 74, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 73) * W_SP( 73, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 75) * W_SP( 75, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 76) * W_SP( 76, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 81) * W_SP( 81, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 82) * W_SP( 82, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 83) * W_SP( 83, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 84) * W_SP( 84, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 85) * W_SP( 85, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 86) * W_SP( 86, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 87) * W_SP( 87, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 88) * W_SP( 88, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 89) * W_SP( 89, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 90) * W_SP( 90, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 91) * W_SP( 91, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 92) * W_SP( 92, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 93) * W_SP( 93, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 94) * W_SP( 94, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 95) * W_SP( 95, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 96) * W_SP( 96, Ixy, Iz)
         elseif (if_state==3) then
            HFFset(I_PinX, :) = macro_xs_table(ixy,iz)%HFF(I_PinX, :,  1) - &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  1) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 77) * W_SP( 77, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 78) * W_SP( 78, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 79) * W_SP( 79, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 80) * W_SP( 80, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 74) * W_SP( 74, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 73) * W_SP( 73, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 75) * W_SP( 75, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 76) * W_SP( 76, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 97) * W_SP( 97, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 98) * W_SP( 98, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 99) * W_SP( 99, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,100) * W_SP(100, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,102) * W_SP(102, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,101) * W_SP(101, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,103) * W_SP(103, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,104) * W_SP(104, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,105) * W_SP(105, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,106) * W_SP(106, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,107) * W_SP(107, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,108) * W_SP(108, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,109) * W_SP(109, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,110) * W_SP(110, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,111) * W_SP(111, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,112) * W_SP(112, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,114) * W_SP(114, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,113) * W_SP(113, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,115) * W_SP(115, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,116) * W_SP(116, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,117) * W_SP(117, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,118) * W_SP(118, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,119) * W_SP(119, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,120) * W_SP(120, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,121) * W_SP(121, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,122) * W_SP(122, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,123) * W_SP(123, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,124) * W_SP(124, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,126) * W_SP(126, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,125) * W_SP(125, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,127) * W_SP(127, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,128) * W_SP(128, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,129) * W_SP(129, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,130) * W_SP(130, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,131) * W_SP(131, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,132) * W_SP(132, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,133) * W_SP(133, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,134) * W_SP(134, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,135) * W_SP(135, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,136) * W_SP(136, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 84) * W_SP( 84, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,137) * W_SP(137, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,138) * W_SP(138, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,139) * W_SP(139, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 88) * W_SP( 88, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,140) * W_SP(140, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,141) * W_SP(141, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,142) * W_SP(142, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 92) * W_SP( 92, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,143) * W_SP(143, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,144) * W_SP(144, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,145) * W_SP(145, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 96) * W_SP( 96, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,146) * W_SP(146, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,147) * W_SP(147, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,148) * W_SP(148, Ixy, Iz)
         elseif (if_state==4) then
            HFFset(I_PinX, :) = macro_xs_table(ixy,iz)%HFF(I_PinX, :,  1) - &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,  1) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,126) * W_SP(126, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,125) * W_SP(125, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,127) * W_SP(127, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,128) * W_SP(128, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,129) * W_SP(129, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,130) * W_SP(130, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,131) * W_SP(131, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,132) * W_SP(132, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,133) * W_SP(133, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,134) * W_SP(134, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,135) * W_SP(135, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,136) * W_SP(136, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,150) * W_SP(150, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,149) * W_SP(149, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,151) * W_SP(151, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,152) * W_SP(152, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,153) * W_SP(153, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,154) * W_SP(154, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,155) * W_SP(155, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,156) * W_SP(156, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,157) * W_SP(157, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,158) * W_SP(158, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,159) * W_SP(159, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,160) * W_SP(160, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,161) * W_SP(161, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,162) * W_SP(162, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,163) * W_SP(163, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,164) * W_SP(164, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 17) * W_SP( 17, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 18) * W_SP( 18, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 19) * W_SP( 19, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 20) * W_SP( 20, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 14) * W_SP( 14, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 13) * W_SP( 13, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 15) * W_SP( 15, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :, 16) * W_SP( 16, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,139) * W_SP(139, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,165) * W_SP(165, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,166) * W_SP(166, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,142) * W_SP(142, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,167) * W_SP(167, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,168) * W_SP(168, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,145) * W_SP(145, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,169) * W_SP(169, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,170) * W_SP(170, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,148) * W_SP(148, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,171) * W_SP(171, Ixy, Iz) + &
                                macro_xs_table(ixy,iz)%HFF(I_PinX, :,172) * W_SP(172, Ixy, Iz)
         endif
      enddo

      return
      end subroutine MACROHFFsetFB_1D

      end module Mod_LiteVer


#endif
