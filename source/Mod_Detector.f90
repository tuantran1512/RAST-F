
      MODULE Mod_Detector

      USE Inc_Constant
      USE Inc_Flag
      USE Inc_Geometry
      USE Inc_RP

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE ExCoreDET
      USE Inc_3D
      USE Inc_Detector
      USE Inc_maXS
      USE Inc_Option, only: N_Group
      use mod_alloc
      implicit none
      real(8) :: Sum_Real1, Sum_Real2, nu_per_kap
      integer :: line, ii, ig
      integer :: ix, iy, iz, ixy

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [ExCoreDET] in Mod_Detector'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [ExCoreDET] in Mod_Detector'
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

      if (.not.flag_exdet) return

      if (.not. allocated(drf_core)) then
         call alloc0(drf_core,1,Nxy,1,Nz,1,3)
         call alloc0(dr_core,1,3)
         drf_core = 0d0
         dr_core = 0d0
         if (opt_core == 4) then
            if (opt_drf == 1) then
               line = size(raw_drf_opr_skn1,2)
               do ii = 1,line
                  ix = int(raw_drf_opr_skn1(1,ii),4)
                  iy = int(raw_drf_opr_skn1(2,ii),4)
                  iz = int(raw_drf_opr_skn1(3,ii),4)
                  Ixy = IxIyToIxy(ix,iy)
                  drf_core(Ixy,iz+1,1) = raw_drf_opr_skn1(4,ii)
                  drf_core(Ixy,iz+1,2) = raw_drf_opr_skn1(5,ii)
                  drf_core(Ixy,iz+1,3) = raw_drf_opr_skn1(6,ii)
               enddo
            elseif (opt_drf == 2) then
               line = size(raw_drf_opr_ygn3,2)
               do ii = 1,line
                  ix = int(raw_drf_opr_ygn3(1,ii),4)
                  iy = int(raw_drf_opr_ygn3(2,ii),4)
                  iz = int(raw_drf_opr_ygn3(3,ii),4)
                  Ixy = IxIyToIxy(ix,iy)
                  drf_core(Ixy,iz+1,1) = raw_drf_opr_ygn3(4,ii)
                  drf_core(Ixy,iz+1,3) = raw_drf_opr_ygn3(5,ii)
               enddo
            elseif (opt_drf == 3) then
               line = size(raw_drf_wh_ucn12,2)
               do ii = 1,line
                  ix = int(raw_drf_wh_ucn12(1,ii),4)
                  iy = int(raw_drf_wh_ucn12(2,ii),4)
                  iz = int(raw_drf_wh_ucn12(3,ii),4)
                  Ixy = IxIyToIxy(ix,iy)
                  drf_core(Ixy,iz+1,1) = raw_drf_wh_ucn12(4,ii)
                  drf_core(Ixy,iz+1,3) = raw_drf_wh_ucn12(5,ii)
               enddo
            endif
         elseif (opt_core == 1) then
            if (opt_drf == 1) then
               line = size(raw_drf_opr_skn1,2)
               do ii = 1,line
                  ix = int(raw_drf_opr_skn1(1,ii),4)
                  iy = int(raw_drf_opr_skn1(2,ii),4)
                  iz = int(raw_drf_opr_skn1(3,ii),4)
                  Ixy = ixiytoixy(17+ix,17+iy)
                  drf_core(ixy,iz+1,1) = raw_drf_opr_skn1(4,ii) / 4d0
                  drf_core(ixy,iz+1,2) = raw_drf_opr_skn1(5,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_opr_skn1(6,ii) / 4d0
                  Ixy = ixiytoixy(18-ix,17+iy)
                  drf_core(ixy,iz+1,1) = raw_drf_opr_skn1(4,ii) / 4d0
                  drf_core(ixy,iz+1,2) = raw_drf_opr_skn1(5,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_opr_skn1(6,ii) / 4d0
                  Ixy = ixiytoixy(18-ix,18-iy)
                  drf_core(ixy,iz+1,1) = raw_drf_opr_skn1(4,ii) / 4d0
                  drf_core(ixy,iz+1,2) = raw_drf_opr_skn1(5,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_opr_skn1(6,ii) / 4d0
                  Ixy = ixiytoixy(17+ix,18-iy)
                  drf_core(ixy,iz+1,1) = raw_drf_opr_skn1(4,ii) / 4d0
                  drf_core(ixy,iz+1,2) = raw_drf_opr_skn1(5,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_opr_skn1(6,ii) / 4d0
               enddo
            elseif (opt_drf == 2) then
               line = size(raw_drf_opr_ygn3,2)
               do ii = 1,line
                  ix = int(raw_drf_opr_ygn3(1,ii),4)
                  iy = int(raw_drf_opr_ygn3(2,ii),4)
                  iz = int(raw_drf_opr_ygn3(3,ii),4)
                  Ixy = ixiytoixy(17+ix,17+iy)
                  drf_core(ixy,iz+1,1) = raw_drf_opr_ygn3(4,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_opr_ygn3(5,ii) / 4d0
                  Ixy = ixiytoixy(18-ix,17+iy)
                  drf_core(ixy,iz+1,1) = raw_drf_opr_ygn3(4,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_opr_ygn3(5,ii) / 4d0
                  Ixy = ixiytoixy(18-ix,18-iy)
                  drf_core(ixy,iz+1,1) = raw_drf_opr_ygn3(4,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_opr_ygn3(5,ii) / 4d0
                  Ixy = ixiytoixy(17+ix,18-iy)
                  drf_core(ixy,iz+1,1) = raw_drf_opr_ygn3(4,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_opr_ygn3(5,ii) / 4d0
               enddo
            elseif (opt_drf == 3) then
               line = size(raw_drf_wh_ucn12,2)
               do ii = 1,line
                  ix = int(raw_drf_wh_ucn12(1,ii),4)
                  iy = int(raw_drf_wh_ucn12(2,ii),4)
                  iz = int(raw_drf_wh_ucn12(3,ii),4)
                  Ixy = ixiytoixy(17+ix,17+iy)
                  drf_core(ixy,iz+1,1) = raw_drf_wh_ucn12(4,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_wh_ucn12(5,ii) / 4d0
                  Ixy = ixiytoixy(18-ix,17+iy)
                  drf_core(ixy,iz+1,1) = raw_drf_wh_ucn12(4,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_wh_ucn12(5,ii) / 4d0
                  Ixy = ixiytoixy(18-ix,18-iy)
                  drf_core(ixy,iz+1,1) = raw_drf_wh_ucn12(4,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_wh_ucn12(5,ii) / 4d0
                  Ixy = ixiytoixy(17+ix,18-iy)
                  drf_core(ixy,iz+1,1) = raw_drf_wh_ucn12(4,ii) / 4d0
                  drf_core(ixy,iz+1,3) = raw_drf_wh_ucn12(5,ii) / 4d0
               enddo
            endif
         endif
      endif

      Sum_Real1 = 0d0
      Sum_Real2 = 0d0
      do ig = 1, N_Group
         do iz = izfuelbot, izfueltop
            do ixy = 1, nxy
               if (kap_maxs_f_3d(ixy,iz,ig)<1d-30) cycle
               Sum_Real1 = sum_real1 + nu_maxs_f_3d(ixy,iz,ig) / kap_maxs_f_3d(ixy,iz,ig)  &
                         * maxs_f_3d(ixy,iz,ig) * flux(ixy,iz,ig)
               Sum_Real2 = sum_real2 + maxs_f_3d(ixy,iz,ig) * flux(ixy,iz,ig)
            enddo
         enddo
      enddo
      if (sum_real2<1d-30) then
         nu_per_kap = 1d0
      else
         nu_per_kap = sum_real1 / sum_real2
      endif

      do ig = 1, 3
         Sum_Real1 = 0d0
         do iz = IzFuelbot, izfueltop
            do ixy = 1,nxy
               Sum_Real1 = sum_real1 + drf_core(ixy,iz,ig) * power(ixy,iz)
            enddo
         enddo
         !js+dr_core(ig) = sum_real1
         dr_core(ig) = sum_real1 * nu_per_kap
      enddo

      RETURN
      END SUBROUTINE Excoredet


#ifdef siarhei_delete 
      SUBROUTINE Init_Det_sig
      USE Inc_3D
      USE Inc_Detector
      USE Inc_maXS
      USE Inc_Geometry
      USE Inc_Option , only: n_group
      use mod_alloc
      IMPLICIT NONE
      REAL(8) :: Sum_Real
      INTEGER :: Sum_Int
      INTEGER :: Ixy_1N, m, iz
      REAL(8) :: z_index(0:nz)

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Init_Det_sig] in Mod_Detector'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      Sum_Int = 0
      do ixy_1n = 1, Nxy_1n
         if (DET_1N(ixy_1n)==1) then
            Sum_Int = Sum_int + 1
            DET_XY_LIST(sum_int) = ixy_1n
         endif
      enddo

      if (.not.allocated(det_z_start_end)) call alloc(det_z_start_end,det_n_z,2)
      if (.not.allocated(det_h_start_end)) call alloc(det_h_start_end,det_n_z,2)

      if (opt_det_z_grid) then
         do m = 1, det_n_z
            det_z_start_end(m,1) = det_z_beg(m)
            det_z_start_end(m,2) = det_z_end(m)
            det_h_start_end(m,1) = h_z(det_z_beg(m))
            det_h_start_end(m,2) = h_z(det_z_end(m))
         enddo
      else
         Sum_Real = 0d0
         do Iz = 1, IzFuelbot-1
            Sum_Real = Sum_real + h_z(iz)
         enddo
         do m = 1, det_n_z
            det_y(m) = det_y(m) + sum_real
         enddo
         z_index = 0d0
         do Iz = 1, Nz
            z_index(Iz) = z_index(iz-1) + h_z(iz)
         enddo
         do m = 1, det_n_z
            do iz = 1, nz
               if ( det_y(m) > z_index(iz) ) then
                  det_z_start_end(m,1) = iz + 1
                  det_h_start_end(m,1) = z_index(iz+1) - det_y(m)
               endif
               if ( det_y(m) + det_h > z_index(iz) ) then
                  det_z_start_end(m,2) = iz + 1
                  det_h_start_end(m,2) = det_y(m) + det_h - z_index(iz)
               endif
            enddo
         enddo
      endif

      prt_det=0
      select case (opt_rr_det)
      case (0)
         prt_det=0
      case (1)
         do i_det=1,n_detector
            if (DET_ZA(i_det)==45103) then
               prt_det=i_det
            endif
         enddo
         if (prt_det==0) then
            write(*,*) ' error: there is no detector mixs in xs file ...'
            write(*,*) ' - only flux is printed'
         endif
      case (2)
         do i_det=1,n_detector
            if (DET_ZA(i_det)==23051) then
               prt_det=i_det
            endif
         enddo
         if (prt_det==0) then
            write(*,*) ' error: there is no detector mixs in xs file ...'
            write(*,*) ' - only flux is printed'
         endif
      case (3)
         do i_det=1,n_detector
            if (DET_ZA(i_det)==27059) then
               prt_det=i_det
            endif
         enddo
         if (prt_det==0) then
            write(*,*) ' error: there is no detector mixs in xs file ...'
            write(*,*) ' - only flux is printed'
         endif
      case (4)
         do i_det=1,n_detector
            if (DET_ZA(i_det)==47107) then
               prt_det=i_det
            endif
         enddo
         if (prt_det==0) then
            write(*,*) ' error: there is no detector mixs in xs file ...'
            write(*,*) ' - only flux is printed'
         endif
      case (5)
         do i_det=1,n_detector
            if (DET_ZA(i_det)==92235) then
               prt_det=i_det
            endif
         enddo
         if (prt_det==0) then
            write(*,*) ' error: there is no detector mixs in xs file ...'
            write(*,*) ' - only flux is printed'
         endif
      end select

      if(.not.allocated(det_sig_rr)) call alloc(det_sig_rr,det_n_xy,det_n_z)
      if(.not.allocated(det_rr_xyz)) call alloc(det_rr_xyz,nxy_1n,nz,n_group)
      if (opt_prt_wprime.or.opt_prt_detpow) then
         if(.not.allocated(det_power_xyz)) call alloc( det_power_xyz,nxy_1n,nz)
         if(.not.allocated(det_detpow))    call alloc( det_detpow,det_n_xy,det_n_z)
         if(.not.allocated(det_wprime))    call alloc( det_wprime,det_n_xy,det_n_z)
      endif

      return
      END SUBROUTINE Init_Det_sig
#endif 


      SUBROUTINE Cal_Det_sig
      use Inc_3D
      use Inc_Detector
      use Inc_maXS
      use Inc_Geometry
      use Inc_Option, only: n_group
     ! use mod_xsfb, only: DetXSFB
      use inc_option, only: opt_mode
      implicit none
      integer :: ixy_1n, m, ixy, iz, ig, n
      real(8) :: sum_real, sum_real_z
      integer :: sum_int

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Cal_Det_sig] in Mod_Detector'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [Cal_Det_sig] in Mod_Detector'
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


      if (.not.flag_det_sig.or.opt_mode==3) return

!      if (.not.flag_init_det_sig) call Init_Det_Sig
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif !      call DetXSFB
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
      do ig = 1, n_group
         do iz = 1, nz
            do ixy_1n = 1, nxy_1n
               sum_real = 0d0
               sum_int = 0
               do m = 1, 4
                  ixy = i_1nto4n(ixy_1n, m)
                  if (ixy /= 0) then
                     sum_int = sum_int + 1
                     if (opt_prt_det==10.or.opt_prt_det==11.or.opt_prt_det==12) then
                        sum_real = sum_real + flux(ixy,iz,ig)
                     elseif (opt_prt_det==20.or.opt_prt_det==21.or.opt_prt_det==22) then
                        sum_real = sum_real + flux(ixy,iz,ig)*det_fr(ixy,iz,ig)
                     elseif (opt_prt_det==30.or.opt_prt_det==31.or.opt_prt_det==32) then
                        sum_real = sum_real + flux(ixy,iz,ig)*det_mixs(ixy,iz,ig,prt_det)
                     elseif (opt_prt_det==40.or.opt_prt_det==41.or.opt_prt_det==42) then
                        sum_real = sum_real + flux(ixy,iz,ig)*det_fr(ixy,iz,ig)*det_mixs(ixy,iz,ig,prt_det)
                     endif
                  endif
               enddo
               det_rr_xyz(ixy_1n,iz,ig) = sum_real / real(sum_int,8)
            enddo
         enddo
      enddo

      do n = 1, DET_N_XY
         ixy_1n = DET_XY_list(n)
         do m = 1, DET_n_z
            sum_real = 0d0
            sum_real_z = 0d0
            do ig = det_g_beg, det_g_end
               do iz = det_z_start_end(m,1), det_z_start_end(m,2)
                  if (iz == det_z_start_end(m,1)) then
                     sum_real = sum_real &
                          + det_rr_xyz(ixy_1n,iz,ig)*det_h_start_end(m,1)
                     sum_real_z = sum_real_z + det_h_start_end(m,1)
                  elseif (iz == det_z_start_end(m,2)) then
                     sum_real = sum_real &
                          + det_rr_xyz(ixy_1n,iz,ig)*det_h_start_end(m,2)
                     sum_real_z = sum_real_z + det_h_start_end(m,2)
                  else
                     sum_real = sum_real &
                          + det_rr_xyz(ixy_1n,iz,ig)*h_z(iz)
                     sum_real_z = sum_real_z + h_z(iz)
                  endif
               enddo
            enddo
            det_sig_rr(n,m) = sum_real / sum_real_z
         enddo
      enddo

      if (opt_prt_wprime.or.opt_prt_detpow) then
         do ixy_1n = 1, nxy_1n
            do iz = 1, nz
               sum_real = 0d0
               sum_int = 0
               do m = 1, 4
                  ixy = i_1nto4n(ixy_1n, m)
                  if (ixy /= 0) then
                     sum_int = sum_int + 1
                     sum_real = sum_real + power(ixy,iz) * nodevolume(ixy,iz) ! W
                  endif
               enddo
               det_power_xyz(ixy_1n,iz) = sum_real / real(sum_int,8) * 4d0 / 1d+6  ! MW
            enddo
         enddo

         do n = 1, DET_N_XY
            ixy_1n = DET_XY_list(n)
            do m = 1, DET_n_z
               sum_real = 0d0
               sum_real_z = 0d0
               do iz = det_z_start_end(m,1), det_z_start_end(m,2)
                  if (iz == det_z_start_end(m,1)) then
                     sum_real = sum_real + det_power_xyz(ixy_1n,iz)*det_h_start_end(m,1)/h_z(iz)
                  elseif (iz == det_z_start_end(m,2)) then
                     sum_real = sum_real + det_power_xyz(ixy_1n,iz)*det_h_start_end(m,2)/h_z(iz)
                  else
                     sum_real = sum_real + det_power_xyz(ixy_1n,iz)*h_z(iz)/h_z(iz)
                  endif
               enddo
               det_detpow(n,m) = sum_real ! MW
               det_wprime(n,m) = det_detpow(n,m) / det_sig_rr(n,m) / signal_to_power * 1d+18
            enddo
         enddo
      endif

      return
      END SUBROUTINE Cal_Det_sig


#ifdef siarhei_delete 
      SUBROUTINE INCORE_det(nn_time)

      USE Inc_Detector, only: &
          kpp, kdd, &
          DET_b_A1, DET_b_a2, det_b_a2m, det_b_a3, det_b_b2, det_b_b3, &
          DET_b_V52, DEt_d_v51, det_d_v52, &
          INIT_COND_DET, i_bu_old, &
          det_n_z, det_n_xy, det_sig_phi_g, dET_SIG_PHI_G_Old, &
          DET_Flux_XYZ_g, det_flux_xyz_g_old, &
          DET_XY_LIST, det_z_start_end, det_h_start_end, det_i_fa
      USE Inc_3D, ONLY: opt_det, NN_TIME2, OPT_DET_Power, I_ITER_JR
      USE Inc_Geometry
      USE Inc_RP
      USE Inc_Depletion
      USE Inc_miXS
      USE Inc_Nuclide
      USE Inc_Option, Only: n_group
      USE Inc_Time, ONLy: dt
!     ! USE Mod_Operator, only: get_inva, get_invc
      USE Inc_XS_File, only: i_bu
      USE Inc_History
      USE Inc_PinVar
      use inc_pinpow, only: flag_pinpow
      use mod_alloc
      IMPLICIT NONE
      COMPLEX(8), DIMENsion(6, 6) :: sum_cram_mat
      COMPLEX(8), DIMENsion(2, 2) :: sum_cram_mat_v
      COMPLEX(8), DIMENsion(4, 4) :: sum_cram_mat_si
      REAL(8), DIMENSIOn(6, 6) :: a, eat
      REAL(8), DIMENSIOn(2, 2) :: a_v, eat_v
      REAL(8), DIMENSIOn(4, 4) :: a_si, eat_si
      REAL(8), DIMENSIOn(6) :: vec_n, vec_b
      REAL(8), DIMENSIOn(2) :: vec_n_v, vec_b_v
      REAL(8), DIMENSIOn(4) :: vec_n_si, vec_b_si
      INTEGER :: Ixy, Iz, ig, i
      INTEGER :: NN_TIMe
      REAL(8) :: DET_b_a4
      REAL(8) :: DET_IT_a1, det_it_a2, det_it_a2m, det_it_a3, det_it_b2, det_it_b3
      REAL(8) :: DET_IT_a4
      REAL(8) :: DET_d_a1, det_d_a2, det_d_a2m, det_d_a3, det_d_b2, det_d_b3
      REAL(8) :: DET_d_a4
      REAL(8) :: DET_c_a1  , det_a_a1  , det_f_a1
      REAL(8) :: DET_c_a2  , det_a_a2  , det_f_a2
      REAL(8) :: DET_c_a2m , det_a_a2m , det_f_a2m
      REAL(8), allocatable, save:: a_a1_jr(:), a_a2_jr(:), a_a2m_jr(:)
      INTEGER :: DET_N_fa
      INTEGER :: DET_N_rf
      REAL(8) :: Sum_Real, sum_real2, sum_real_old, sum_real2_old
      INTEGER :: Ixy_1N, m, n
      INTEGER :: Ixy_1N_before
      REAL(8) :: DET_c_a3, det_a_a3, det_f_a3
      REAL(8) :: DET_c_b2, det_a_b2, det_f_b2
      REAL(8) :: DET_c_b3, det_a_b3, det_f_b3
      REAL(8) :: DET_a_a4, det_f_a4
      REAL(8) :: DET_c_a4
      REAL(8) :: DET_r_a1
      REAL(8) :: DET_r_a2
      REAL(8) :: DET_r_a2m
      REAL(8) :: DET_r_a3
      REAL(8) :: DET_r_a4
      REAL(8) :: DET_r_b2
      REAL(8) :: DET_r_b3
      REAL(8) :: DET_w_a1
      REAL(8) :: DET_w_a2
      REAL(8) :: DET_w_a2m
      REAL(8) :: DET_w_a3
      REAL(8) :: DET_w_b2
      REAL(8) :: DET_w_b3
      REAL(8) :: NAA
      REAL(8) :: DET_f_v1, det_f_v2
      REAL(8) :: DET_a_v1, det_a_v2
      REAL(8) :: DET_c_v1, det_c_v2
      REAL(8) :: DET_r_v1, det_r_v2
      INTEGER :: flag_lfa
      INTEGER :: det_p_1, det_p_2
      INTEGER :: mid_det_p, mid_mid_det_p_test
      REAL(8) :: abundance

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [INCORE_det] in Mod_Detector'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.flag_pinpow) then
         stop "flag_pinpow should be turned on"
      endif

      if (.not.allocated(det_sig_phi_g))     call alloc(det_sig_phi_g    ,det_n_xy, det_n_z, 2)
      if (.not.allocated(det_sig_phi_g_old)) call alloc(det_sig_phi_g_old,det_n_xy, det_n_z, 2)

      NN_TIME2 = NN_TIMe

      !OPEN(2020, FILE='detector_pinfluxes.res',position="append")
      !OPEN(2021, FILE='detector_fafluxes.res',position="append")
      !OPEN(3030, FILE='test_absorption',position="append")
      !OPEN(3031, FILE='rast-k flux',position="append")

      Kdd = 1-Kpp
      NAA = 6.02252D23
      flag_lfa=1
      if (.not.allocated(a_A1_JR )) call alloc (a_A1_JR ,izfueltop-izfuelbot+1)
      if (.not.allocated(a_A2_JR )) call alloc (a_A2_JR ,izfueltop-izfuelbot+1)
      if (.not.allocated(a_A2m_JR)) call alloc (a_A2m_JR,izfueltop-izfuelbot+1)
      !write(3030,"(A, i5)") "ixy is", i_iter_jr

      I_ITER_JR = I_ITEr_jr  + 1
      IF (INIT_COND_DET==1 .and. opt_det_power==1) then
         DET_N_FA = 0
         DET_N_RF = 0
         Ixy_1N_before = 0
         mid_mid_det_p_test = npin*1
         mid_det_p = mod(mid_mid_det_p_test,2)
         IF (mid_det_p == 0d0) then
            det_p_1 = npin/2
            det_p_2 = npin/2 + 1
         ELSE
            det_p_1 = npin/2 + 1
            det_p_2 = npin/2 + 1
         END IF

         IF (NN_TIME.EQ.2) then
            DO Ixy = 1, nxy
               Ixy_1N = i_4nto1n(ixy)
               IF ( I_Farf_1n( ixy_1n ) == 0 ) cycle
               IF ( I_Farf_1n( ixy_1n ) == 2 ) then
                  IF (Ixy_1n_before < ixy_1n ) then
                     DEt_n_fa = det_n_fa + 1
                     DEt_i_fa(ixy_1n) = det_n_fa
                  ELSE if (ixy_1n_before > ixy_1n ) then
                     DEt_n_fa = det_i_fa(ixy_1n)
                  END If
               ELSE IF ( i_farf_1n( ixy_1n )  == 1 ) then
                  IF (Ixy_1n_before < ixy_1n ) then
                     DEt_i_fa(ixy_1n) = det_n_fa
                  END If
               END IF
               Ixy_1N_before = ixy_1n
            END DO
         END IF

         DO Ixy = 1, Nxy
            Ixy_1N = I_4nto1n(ixy)
            IF ( I_FARF_1n( ixy_1n ) == 0 ) cycle
            DO Iz = IzFuelbot, izfueltop
               IF ( I_Farf_1n( ixy_1n ) == 2 ) then
                  DET_Flux_xyz_g_old(ixy_1n,iz,:) &
                     & = ( phihom(det_p_1,det_p_1,det_I_FA(Ixy_1N),Iz,:) &
                     &   + phihom(det_p_1,det_p_2,DET_I_FA(Ixy_1N),Iz,:) &
                     &   + phihom(det_p_2,det_p_1,det_i_fA(Ixy_1N),Iz,:) &
                     &   + phihom(det_p_2,det_p_2,DET_I_FA(Ixy_1N),Iz,:) ) / 4
                  DET_Flux_xyz_g(ixy_1n,iz,:) &
                     & = ( phihom(det_p_1,det_p_1,det_I_FA(Ixy_1N),Iz,:) &
                     &   + phihom(det_p_1,det_p_2,DET_I_FA(Ixy_1N),Iz,:) &
                     &   + phihom(det_p_2,det_p_1,det_i_fA(Ixy_1N),Iz,:) &
                     &   + phihom(det_p_2,det_p_2,DET_I_FA(Ixy_1N),Iz,:) ) / 4
               END IF
            END DO
         END DO

         Do n = 1, DET_n_xy
            Ixy_1N = DET_xy_list(n)
            Do m = 1, Det_n_z
               Sum_Real = 0d0
               Sum_Real2 = 0d0
               Sum_Real_old = 0d0
               Sum_Real2_old = 0d0
               Do Ig = 1, n_group
                  Do Iz = det_z_start_end(m,1), det_z_start_end(m,2)
                     IF(iz == det_z_start_end(m,1)) then
                        sum_real = sum_real + det_flux_xyz_g(ixy_1n,iz,ig) * det_h_start_end(m,1)
                        sum_real2 = sum_real2 + det_h_start_end(m,1)
                        sum_real_old = sum_real_old + det_flux_xyz_g_old(ixy_1n,iz,Ig) * det_h_start_end(m,1)
                        sum_real2_old = sum_real2_old + det_h_start_end(m,1)
                     ELseif (iz == det_z_start_end(m,2)) then
                        sum_real = sum_real + det_flux_xyz_g(ixy_1n,iz,ig) * det_h_start_end(m,2)
                        sum_real2 = sum_real2 + det_h_start_end(m,2)
                        sum_real_old = sum_real_old + det_flux_xyz_g_old(ixy_1n,iz,Ig) * det_h_start_end(m,2)
                        sum_real2_old = sum_real2_old + det_h_start_end(m,2)
                     ELse
                        sum_real = sum_real + det_flux_xyz_g(ixy_1n,iz,ig) * h_z(Iz)
                        sum_real2 = sum_real2 + h_z(iz)
                        sum_real_old = sum_real_old + det_flux_xyz_g_old(ixy_1n,iz,Ig) * h_z(Iz)
                        sum_real2_old = sum_real2_old + h_z(iz)
                     ENdif
                  ENDDO
                  DET_Sig_phi_g(n,m,ig) = sum_real / sum_real2
                  DET_Sig_phi_g_old(n,m,ig) = sum_real_old / sum_real2_old
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      DO n = 1, DET_N_Xy
         Do m = 1, DET_n_z
            IF ( OPT_Det == 1 ) then
               DET_miXS_a_a1  (:,:,1) = 0.d0
               DET_miXS_a_a2  (:,:,1) = 0.d0
               DET_miXS_a_a2m (:,:,1) = 0.d0
               DET_miXS_a_a3  (:,:,1) = 0.d0
               DET_miXS_a_b2  (:,:,1) = 0.d0
               DET_miXS_a_b3  (:,:,1) = 0.d0

               DET_miXS_a_a1  (:,:,2) =   145.160d0 * barn
               DET_miXS_a_a2  (:,:,2) =    39.375d0 * barn
               DET_miXS_a_a2m (:,:,2) =   805.320d0 * barn
               DET_miXS_a_a3  (:,:,2) = 12475.000d0 * barn
               DET_miXS_a_b2  (:,:,2) =     0.600d0 * barn
               DET_miXS_a_b3  (:,:,2) =    20.520d0 * barn

               DET_b_A1  = 0.d0
               DET_b_A2  = 1.63865d-2
               DET_b_A2m = 2.66186d-3
               DET_b_A3  = 5.44517d-6
               DET_b_B2  = 0.d0
               DET_b_B3  = 0.d0

               DET_IT_A1  = 0.d0
               DET_IT_A2  = 0.d0
               DET_IT_A2m = 2.66186d-3
               DET_IT_A3  = 0.d0
               DET_IT_B2  = 0.d0
               DET_IT_B3  = 0.d0

               DET_d_A1  = 0.d0
               DET_d_A2  = det_b_a2
               DET_d_A2m = det_it_a2m
               DET_d_A3  = det_b_a3
               DET_d_B2  = 0.d0
               DET_d_B3  = 0.d0

               DET_w_A1  = 0.d0
               DET_w_A2  = 0.923d0
               DET_w_A2m = 0.077d0
               DET_w_A3  = 0.d0
               DET_w_B2  = 0.d0
               DET_w_B3  = 0.d0

               IF (NN_Time.eq.2) then
                  IF (Init_cond_det==0) then
                    DET_n_a1_old  (:,:) = (12.41 / 102.90550 )* naa
                    DET_n_a2_old  (:,:) = 0.d0
                    DET_n_a2m_old (:,:) = 0.d0
                    DET_n_a3_old  (:,:) = 0.d0
                    DET_n_b2_old  (:,:) = 0.d0
                    DET_n_b3_old  (:,:) = 0.d0
                  ELSE if (init_cond_det==1) then
                    DET_n_a1_old (:,:) = (12.41 / 102.90550 )*naa
                    DET_a_a1  = dot_product( det_mixs_a_a1 (n,m,:), det_sig_phi_g_Old(n,m,:) )
                    DET_a_a2  = dot_product( det_mixs_a_a2 (n,m,:), det_sig_phi_g_Old(n,m,:) )
                    DET_a_a2m = dot_product( det_mixs_a_a2m(n,m,:), det_sig_phi_g_Old(n,m,:) )
                    DET_f_a1  = 0.d0
                    DET_c_a1  = det_a_a1  - det_f_a1
                    DET_n_a2m_old(:,:) = det_w_a2m*(det_c_a1*det_n_a1_old(:,:))/(DET_a_A2m+DET_b_A2m)
                    DET_n_a2_old (:,:) = (det_w_a2*det_c_a1*det_n_a1_old(:,:)+det_b_A2m*DET_N_A2m_Old(:,:))/(DET_a_A2+DET_b_A2)
                    DET_n_a3_old (:,:) = 0.d0
                    DET_n_b2_old (:,:) = 0.d0
                    DET_n_b3_old (:,:) = 0.d0
                  END If
               END IF

            ELSE IF (OPt_det == 2) then
               DET_miXS_a_v1(:,:,1) = 0.d0
               DET_miXS_a_v1(:,:,2) = 5.0 * barn
               DET_b_V52= 0.003086415444652d0
               DET_d_V52=det_b_v52
               DET_d_V51=0d0
               IF (NN_Time.eq.2) then
                  IF (Init_cond_det==0) then
                     DEt_n_v1_old  (:,:) = (6.1d0/50.9439595d0)*naa
                     DEt_n_v2_old  (:,:) = 0.d0
                     DEt_n_cr52_old(:,:) = 0.d0
                  ELSE if (init_cond_det==1) then
                     DEt_n_v1_old  (:,:) = (6.1d0/50.9439595d0)*naa
                     DEt_n_v2_old  (:,:) = n_v1_old(:,:)*det_sig_phi_g_old(:,:,2)*DET_miXS_a_V1(:,:,2)/DET_b_V52
                     DEt_n_cr52_old(:,:) = 0.d0
                  END If
               END IF

            ELSE IF (OPt_det == 3) then
               DET_miXS_a_a1  (:,:,1) = 0.d0
               DET_miXS_a_a2  (:,:,1) = 0.d0
               DET_miXS_a_a2m (:,:,1) = 0.d0
               DET_miXS_a_a3  (:,:,1) = 0.d0
               DET_miXS_a_b2  (:,:,1) = 0.d0
               DET_miXS_a_b3  (:,:,1) = 0.d0
               DET_miXS_a_a1  (:,:,2) = 37.00d0 * barn
               DET_miXS_a_a2  (:,:,2) = 2.0d0 * barn
               DET_miXS_a_a2m (:,:,2) = 0.0d0 * barn
               DET_miXS_a_a3  (:,:,2) = 0.0d0 * barn
               DET_miXS_a_b2  (:,:,2) = 2.400d0 * barn
               DET_miXS_a_b3  (:,:,2) = 2.5094 * barn

               DET_b_A1  = 0.d0
               DET_b_A2  = 4.170693152d-9
               DET_b_A2m = 0.001103702399d0
               DET_b_A3  = 1.16691444539e-04
               DET_b_B2  = 0.d0
               DET_b_B3  = 0.d0
               DET_IT_A1  = 0.d0
               DET_IT_A2  = 0.d0
               DET_IT_A2m = 0.001103702399
               DET_IT_A3  = 0.d0
               DET_IT_B2  = 0.d0
               DET_IT_B3  = 0.d0
               DET_d_A1  = 0.d0
               DET_d_A2  = det_b_a2
               DET_d_A2m = det_it_a2m
               DET_d_A3  = det_b_a3
               DET_d_B2  = 0.d0
               DET_d_B3  = 0.d0
               DET_w_A1  = 0.d0
               DET_w_A2  = 0.444d0
               DET_w_A2m = 1- det_w_a2
               DET_w_A3  = 0.d0
               DET_w_B2  = 0.d0
               DET_w_B3  = 0.d0

               IF (NN_Time.eq.2) then
                  IF (Init_cond_det==0) then
                     DEt_n_a1_old (:,:) = (8.9d0/58.93d0)*naa
                     DEt_n_a2_old (:,:) = 0.d0
                     DEt_n_a2m_old(:,:) = 0.d0
                     DEt_n_a3_old (:,:) = 0.d0
                     DEt_n_b2_old (:,:) = 0.d0
                     DEt_n_b3_old (:,:) = 0.d0
                  ELSE if (init_cond_det==1) then
                     DEt_n_a1_old (:, :) = (8.9d0/58.93d0)*naa
                     DEt_a_a1  = dot_product( det_mixs_a_a1 (n,m,:), det_sig_phi_G_Old(n,m,:) )
                     DEt_a_a2  = dot_product( det_mixs_a_a2 (n,m,:), det_sig_phi_G_Old(n,m,:) )
                     DEt_a_a2m = dot_product( det_mixs_a_a2m(n,m,:), det_sig_phi_G_Old(n,m,:) )
                     DEt_f_a1  = 0.d0
                     DEt_c_a1  = det_a_a1 - det_f_a1
                     DEt_n_a2m_old(:,:) = det_w_a2m*(det_c_a1*det_n_a1_old(:,:))/(DET_a_A2m+DET_b_A2m)
                     DEt_n_a2_old (:,:) =(det_w_a2*det_c_a1*det_n_a1_old(:,:)+det_b_A2m*DET_N_A2m_Old(:,:))/(DET_a_A2+DET_b_A2)
                     DEt_n_a3_old (:,:) = 0.d0
                     DEt_n_b2_old (:,:) = 0.d0
                     DEt_n_b3_old (:,:) = 0.d0
                  END If
               END IF

            ELSE IF (OPt_det == 4) then
               DET_miXS_a_a1(:,:,1) = 0.d0
               DET_miXS_a_a2(:,:,1) = 0.d0
               DET_miXS_a_a3(:,:,1) = 0.d0
               DET_miXS_a_a4(:,:,1) = 0.d0
               DET_miXS_a_a1(:,:,2) = 23.80d0 * barn
               DET_miXS_a_a2(:,:,2) = 0.0d0 * barn
               DET_miXS_a_a3(:,:,2) = 58.80d0 * barn
               DET_miXS_a_a4(:,:,2) = 0.0d0 * barn

               DET_b_A1 = 0.d0
               DET_b_A2 = 0.0048d0
               DET_b_A3 = 0.d0
               DET_b_A4 = 0.0282d0
               DET_IT_A1 = 0.d0
               DET_IT_A2 = 0.d0
               DET_IT_A3 = 0.d0
               DET_IT_A4 = 0.d0
               DET_d_A1 = det_b_a1
               DET_d_A2 = det_b_a2
               DET_d_A3 = det_b_a3
               DET_d_A4 = det_b_a4

               IF (NN_Time.eq.2) then
                  abundance=0.51839
                  IF (Init_cond_det==0) then
                     DEt_n_a1_old(:,:) = abundance*(10.50d0/107.8682d0)*naa
                     DEt_n_a2_old(:,:) = 0.d0
                     DEt_n_a3_old(:,:) = (1-abundance)*(10.50d0/107.8682d0)*naa
                     DEt_n_a4_old(:,:) = 0.d0
                  ELSE if (init_cond_det==1) then
                     DEt_n_a1_old(:,:) = abundance*(10.50d0/107.8682d0)*naa
                     DEt_n_a3_old(:,:) = (1-abundance)*(10.50d0/107.8682d0)*naa
                     DEt_a_a1 = dot_product(det_mixs_a_a1(n,m,:), det_sig_phi_g_old(n,m,:))
                     DEt_a_a3 = dot_product(det_mixs_a_a3(n,m,:), det_sig_phi_g_old(n,m,:))
                     DEt_f_a1 = 0.d0
                     DEt_c_a1 = det_a_a1 - det_f_a1
                     DEt_f_a3 = 0.d0
                     DEt_c_a3 = det_a_a3 - det_f_a3
                     DEt_n_a2_old (:,:) = det_c_a1*n_a1_old (:,:) / det_b_a2
                     DEt_n_a4_old (:,:) = det_c_a3*n_a3_old (:,:) / det_b_a4
                  END If
               END IF
            END IF

            IF ( (OPT_Det==1).or.(opt_det==3) ) then
               DET_a_A1  = dot_product(det_mixs_a_a1 (n,m,:),det_sig_phi_g_old(n,m,:))
               DET_a_A2  = dot_product(det_mixs_a_a2 (n,m,:),det_sig_phi_g_old(n,m,:))
               DET_a_A2m = dot_product(det_mixs_a_a2m(n,m,:),det_sig_phi_g_old(n,m,:))
               DET_a_A3  = dot_product(det_mixs_a_a3 (n,m,:),det_sig_phi_g_old(n,m,:))
               DET_a_B2  = dot_product(det_mixs_a_b2 (n,m,:),det_sig_phi_g_old(n,m,:))
               IF ( abs(det_a_a1 -0d0)<1d-10 ) det_a_a1  = dot_product(det_mixs_a_A1 (n,m,:),DET_SIG_PHI_G(n,m,:))
               IF ( abs(det_a_a2 -0d0)<1d-10 ) det_a_a2  = dot_product(det_mixs_a_A2 (n,m,:),DET_SIG_PHI_G(n,m,:))
               IF ( abs(det_a_a2m-0d0)<1d-10 ) det_a_a2m = dot_product(det_mixs_a_A2m(n,m,:),DET_SIG_PHI_G(n,m,:))
               IF ( abs(det_a_a3 -0d0)<1d-10 ) det_a_a3  = dot_product(det_mixs_a_A3 (n,m,:),DET_SIG_PHI_G(n,m,:))
               IF ( abs(det_a_b2 -0d0)<1d-10 ) det_a_b2  = dot_product(det_mixs_a_B2 (n,m,:),DET_SIG_PHI_G(n,m,:))
               IF ( abs(det_a_b3 -0d0)<1d-10 ) det_a_b3  = dot_product(det_mixs_a_B3 (n,m,:),DET_SIG_PHI_G(n,m,:))

               DET_f_A1  = 0.d0
               DET_f_A2  = 0.d0
               DET_f_A2m = 0.d0
               DET_f_A3  = 0.d0
               DET_f_B2  = 0.d0
               DET_f_B3  = 0.d0

               DET_c_A1  = det_a_a1  - det_f_a1
               DET_c_A2  = det_a_a2  - det_f_a2
               DET_c_A2m = det_a_a2m - det_f_a2m
               DET_c_A3  = det_a_a3  - det_f_a3
               DET_c_B2  = det_a_b2  - det_f_b2
               DET_c_B3  = det_a_b3  - det_f_b3

               DET_r_A1  = det_a_a1  + det_d_a1
               DET_r_A2  = det_a_a2  + det_d_a2
               DET_r_A2m = det_a_a2m + det_d_a2m
               DET_r_A3  = det_a_a3  + det_d_a3
               DET_r_B2  = det_a_b2  + det_d_b2
               DET_r_B3  = det_a_b3  + det_d_b3

               A = 0.D0;
               A(1,1) = -det_r_a1
               A(2,2) = -det_r_a2
               A(2,1) =  det_w_a2*det_c_a1
               A(2,3) =  det_d_a2m
               A(3,3) = -det_r_a2m
               A(3,1) =  det_w_a2m*det_c_a1
               A(4,4) = -det_r_a3
               A(4,2) =  det_c_a2
               A(4,3) =  det_c_a2m
               A(5,5) = -det_r_b2
               A(5,2) =  det_d_a2
               A(5,3) =  det_d_a2m
               A(6,6) = -det_r_b3
               A(6,4) =  det_d_a3
               A(6,5) =  det_c_b2

               Vec_N(1) = det_n_a1_old  (n,m)
               Vec_N(2) = det_n_a2_old  (n,m)
               Vec_N(3) = det_n_a2m_old (n,m)
               Vec_N(4) = det_n_a3_old  (n,m)
               Vec_N(5) = det_n_b2_old  (n,m)
               Vec_N(6) = det_n_b3_old  (n,m)

               Vec_b(:) = 0d0

            ELSE IF (OPt_det==2) then
               DET_a_V1 = dot_product( det_mixs_a_v1(n,m,:),det_sig_phi_g_old(n,m,:) )
               DET_a_V2=0.d0
               IF ( abs(det_a_v1-0d0)<1d-10 ) det_a_v1 = dot_product( det_mixs_a_V1(n,m,:),DET_SIG_PHI_G(n,m,:) )

               DET_f_V1 = 0.d0
               DET_f_V2 = 0.d0

               DET_c_V1 = det_a_v1 - det_f_v1
               DET_c_V2 = det_a_v2 - det_f_v2
               DET_r_V1 = det_a_v1 + det_d_v51
               DET_r_V2 = det_a_v2 + det_d_v52

               A_V = 0.d0
               A_V(1,1) = -det_r_v1
               A_V(2,1) = +det_c_v1
               A_V(2,2) = -det_d_v52

               Vec_N_V(1) = det_n_v1_old(n,m)
               Vec_N_V(2) = det_n_v2_old(n,m)

               Vec_b_V(:) = 0d0

            ELSE IF (OPt_det==4) then
               DET_a_A1 = dot_product( det_mixs_a_a1(n,m,:),det_sig_phi_g_old(n,m,:) )
               DET_a_A2 = dot_product( det_mixs_a_a2(n,m,:),det_sig_phi_g_old(n,m,:) )
               DET_a_A3 = dot_product( det_mixs_a_a3(n,m,:),det_sig_phi_g_old(n,m,:) )
               DET_a_A4 = dot_product( det_mixs_a_a4(n,m,:),det_sig_phi_g_old(n,m,:) )
               IF ( abs(det_a_a1-0d0)<1d-10 ) det_a_a1 = dot_product( det_mixs_a_A1(n,m,:),DET_SIG_PHI_G(n,m,:) )
               IF ( abs(det_a_a2-0d0)<1d-10 ) det_a_a2 = dot_product( det_mixs_a_A2(n,m,:),DET_SIG_PHI_G(n,m,:) )
               IF ( abs(det_a_a3-0d0)<1d-10 ) det_a_a3 = dot_product( det_mixs_a_A3(n,m,:),DET_SIG_PHI_G(n,m,:) )
               IF ( abs(det_a_a4-0d0)<1d-10 ) det_a_a4 = dot_product( det_mixs_a_A4(n,m,:),DET_SIG_PHI_G(n,m,:) )

               DET_f_A1 = 0.d0
               DET_f_A2 = 0.d0
               DET_f_A3 = 0.d0
               DET_f_A4 = 0.d0

               DET_c_A1 = det_a_a1 - det_f_a1
               DET_c_A2 = det_a_a2 - det_f_a2
               DET_c_A3 = det_a_a3 - det_f_a3
               DET_c_A4 = det_a_a4 - det_f_a4

               DET_r_A1 = det_a_a1 + det_d_a1
               DET_r_A2 = det_a_a2 + det_d_a2
               DET_r_A3 = det_a_a3 + det_d_a3
               DET_r_A4 = det_a_a4 + det_d_a4

               A_Si=0d0
               A_Si(1,1) = -det_r_a1
               A_Si(2,2) = -det_d_a2
               A_Si(2,1) =  det_c_a1
               A_Si(3,3) = -det_r_a3
               A_Si(4,3) =  det_c_a3
               A_Si(4,4) = -det_d_a4

               Vec_N_Si(1) = det_n_a1_old(n,m)
               Vec_N_Si(2) = det_n_a2_old(n,m)
               Vec_N_Si(3) = det_n_a3_old(n,m)
               Vec_N_Si(4) = det_n_a4_old(n,m)

               Vec_b_Si(:) = 0d0

            END IF

            IF ( (OPT_Det==1).or.(opt_det==3) ) then
               Sum_CRAM_mat = (d0,d0)
              ! DO i = 1, order_cram/2
              !    Sum_Cram_mat = sum_cram_mat + alpha(i)*get_invc(dt*a-theta(i)*mat_I_FP)
              ! END DO
               eAt = D0
               eAt = alpha_0*REAL(Mat_I_FP,8) + D2*REAL(Sum_CRAM_Mat,8)
              ! Vec_N = MATMUL(eAt,Vec_N)+MATMUL((eAt-REAL(Mat_I_FP,8)),MATMUL(Get_invA(A),Vec_b))
               DET_N_A1 (n,m) = Vec_N(1)
               DET_N_A2 (n,m) = Vec_N(2)
               DET_N_A2m(n,m) = Vec_N(3)
               DET_N_A3 (n,m) = Vec_N(4)
               DET_N_B2 (n,m) = Vec_N(5)
               DET_N_B3 (n,m) = Vec_N(6)

               DET_I_beta_A2 (n,m) = DET_d_A2 *DET_N_A2 (n,m)
               DET_I_beta_A2m(n,m) = DET_d_A2m*DET_N_A2m(n,m)

               IF(OPT_Det==3) then
                   I_beta_A2m(Ixy,Iz) = DET_I_beta_A2m(n,m)*0.0024_8
               ENDIF
               DET_I_beta_A3(n,m) = DET_d_A3*DET_N_A3(n,m)
               DET_I_beta_B2(n,m) = 0d0
               DET_I_beta_b3(n,m) = 0d0

               IF (OPT_Det==1) then
                  DET_I_beta_rh_d(n,m) = Kdd*DET_I_beta_A2(n,m)
                  DET_I_beta_rh_p(n,m) = Kpp*(DET_a_A2+DET_a_A2m)*DET_N_A1(n,m)
               ELSE
                  DET_I_beta_rh_d(n,m) = (DET_I_beta_A2(n,m)+DET_I_beta_A2m(n,m) &
                     & +DET_I_beta_A3(n,m)+DET_I_beta_B2(n,m)+DET_I_beta_B3(n,m))
                  DET_I_beta_rh_p(n,m) = Kpp*DET_a_A1*DET_N_A1(n,m)
               ENDIF
               DET_I_beta_rh(n,m) = DET_I_beta_rh_p(n,m) + DET_I_beta_rh_d(n,m)

            ELSE IF (OPT_Det==2) THEN
                 Sum_CRAM_Mat_V = (D0,D0)
              ! DO i = 1, Order_CRAM/2
              !    Sum_CRAM_Mat_V = Sum_CRAM_Mat_V+alpha(i)*Get_invC(dT*A_V-theta(i)*Mat_I_FP)
              ! END DO
               eAt_V = D0
               eAt_V = alpha_0*REAL(Mat_I_FP,8) + D2*REAL(Sum_CRAM_Mat_V,8)
              ! Vec_N_V = MATMUL(eAt_V,Vec_N_V)+MATMUL((eAt_V-REAL(Mat_I_FP,8)),MATMUL(Get_invA(A_V),Vec_b_V))
               DET_N_V1 (n,m) = Vec_N_V(1)
               DET_N_V2 (n,m) = Vec_N_V(2)
               DET_I_beta_V2(n,m) = DET_d_V52*DET_N_V2(n,m)
            ELSE IF (OPT_Det==4) THEN
               Sum_CRAM_Mat_Si = ( D0, D0 )
              ! DO i = 1, Order_CRAM/2
              !    Sum_CRAM_Mat_Si = Sum_CRAM_Mat_Si+alpha(i)*Get_invC(dT*A_Si-theta(i)*Mat_I_FP)
              ! END DO
               eAt_Si = D0
               eAt_Si = alpha_0*REAL(Mat_I_FP,8) + D2*REAL(Sum_CRAM_Mat_Si, 8);
               Vec_N_Si = MATMUL(eAt_Si,Vec_N_Si)!+MATMUL((eAt_Si-Mat_I_FP),MATMUL(Get_invA(A_Si),Vec_b_Si))
               DET_N_A1(n,m) = Vec_N_Si(1)
               DET_N_A2(n,m) = Vec_N_Si(2)
               DET_N_A3(n,m) = Vec_N_Si(3)
               DET_N_A4(n,m) = Vec_N_Si(4)
               DET_I_beta_si(n,m) = DET_d_A2*DET_N_A2(n,m)+DET_d_A4*DET_N_A4(n,m)
            END IF
         END DO
      END DO

      I_BU_OLD=I_BU

      !close(2020)
      !close(2021)
      !close(3030)

      RETURN
      END SUBROUTINE INCORE_det
#endif 


#ifdef siarhei_delete 
      SUBROUTINE Compensate_Detector (NN_TIME)

      USE Inc_Constant
      USE Inc_Detector
      USE Inc_miXS
      USE Inc_Nuclide
      USE Inc_Geometry
      USE Inc_3D, ONLY: OPT_Det
      USE Inc_Depletion
!     ! USE Mod_Operator, ONLY: Get_InvA, Get_InvC
      USE Inc_History
      USE Inc_Pinvar
      use mod_alloc
      IMPLICIT NONE
      INTEGER :: NN_TIME, Iz, m
      REAL(8) :: ta, a
      REAL(8) :: Sv
      REAL(8) :: DET_N_V1_miXs
      REAL(8) :: zm,Kp,Kd, R
      REAL(8), DIMENSION (2,1):: X2m_temp ,x_temp
      REAL(8), DIMENSION (1,2):: CT, A1, A4
      REAL(8), DIMENSION (1,3):: CT_rh, A1_rh, A4_rh
      REAL(8), DIMENSION (2,1):: CT_trans, A3, KN_temp
      REAL(8), DIMENSION (3,1):: CT_trans_rh, A3_rh, KN_temp_rh
      REAL(8), DIMENSION (1,1):: y_FKH_temp,  A2,CTPMCT , kesi_temp
      REAL(8), DIMENSION (2,2):: Ak,Ak_trans,I,Q,Pm_temp,P_temp, P1_temp,P2_temp, A5
      REAL(8), DIMENSION (3,3):: Ak_rh,Ak_rh_trans, I_rh, Q_rh,Pm_temp_rh,P_temp_rh, P1_temp_rh,P2_temp_rh, A5_rh
      REAL(8), DIMENSION (3,1):: X2m_temp_rh, x_temp_rh
      REAL(8) :: randn
      REAL(8) :: g, pp1,qp,rp
      REAL(8) :: alfa1,alfa2
      REAL(8) :: zg,zgm
      INTEGER :: D_Matrix

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Compensate_Detector] in Mod_Detector'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      D_Matrix = 1000 + N_BU
      if (.not.allocated(temp_flux       )) call alloc ( temp_flux        ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(temp_flux_r     )) call alloc ( temp_flux_r      ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(I_beta_av_temp  )) call alloc ( I_beta_av_temp   ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(DET_flux_av_temp)) call alloc ( DET_flux_av_temp ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(x_kal           )) call alloc ( x_kal            ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(xg_kal          )) call alloc ( xg_kal           ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(xgm_kal         )) call alloc ( xgm_kal          ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(dam2            )) call alloc ( dam2             ,DET_N_XY, DET_N_Z)
      if (.not.allocated(xm_FKH          )) call alloc ( xm_FKH           ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(P               )) call alloc ( P                ,DET_N_XY, DET_N_Z, 2,2,D_Matrix)
      if (.not.allocated(Pm              )) call alloc ( Pm               ,DET_N_XY, DET_N_Z, 2,2,D_Matrix)
      if (.not.allocated(y_FKH           )) call alloc ( y_FKH            ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(KN              )) call alloc ( KN               ,DET_N_XY, DET_N_Z, 2,1,D_Matrix)
      if (.not.allocated(kesi_FKH        )) call alloc ( kesi_FKH         ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(X2m             )) call alloc ( X2m              ,DET_N_XY, DET_N_Z, 2,1,D_Matrix)
      if (.not.allocated(X2m_rh          )) call alloc ( X2m_rh           ,DET_N_XY, DET_N_Z, 3,1,D_Matrix)
      if (.not.allocated(KN_rh           )) call alloc ( KN_rh            ,DET_N_XY, DET_N_Z, 3,1,D_Matrix)
      if (.not.allocated(Pm_rh           )) call alloc ( Pm_rh            ,DET_N_XY, DET_N_Z, 3,3,D_Matrix)
      if (.not.allocated(xg              )) call alloc ( xg               ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(xgm             )) call alloc ( xgm              ,DET_N_XY, DET_N_Z, D_Matrix)
      if (.not.allocated(P_rh            )) call alloc ( P_rh             ,DET_N_XY, DET_N_Z, 3,3,D_Matrix)

      IF (OPT_Det==2) THEN
         ta=1d0/DET_b_V52
         a=1d0/ta
         Do m = 1, DET_N_XY
            Do Iz = 1, DET_N_Z
               DET_N_V1_miXs = DET_miXS_a_V1  (m, Iz, 2) * DET_N_V1  (m, Iz)
               Sv=DET_N_V1_miXs
               zm=exp(-a*Tn_kal)
               Kp=Kpp
               Kd=Kdd
               Ak(1,1)= zm
               Ak(1,2)= 1d0
               Ak(2,1)= 0d0
               Ak(2,2)= 1d0
               CT(1,1)=Kd*(1d0-zm)
               CT(1,2)=Kp
               I(1,1)=1d0
               I(1,2)=0d0
               I(2,1)=0d0
               I(2,2)=1d0
               Q(1,1)=0.0001d0
               Q(1,2)=0d0
               Q(2,1)=0d0
               Q(2,2)=15d0
               P(m,Iz,1,1,1)=1.5d0
               P(m,Iz,1,2,1)=0.0d0
               P(m,Iz,2,1,1)=0.0d0
               P(m,Iz,2,2,1)=1.5d0
               R=8d0

!               randn=random_normal()
            dummy_filler = 1 ! @$^ siarhei_plot 
               temp_flux(m,Iz,1)       = DET_SIG_PHI_G_Old(m,Iz,2)
               temp_flux(m,Iz,NN_TIME) = DET_SIG_PHI_G(m,Iz,2)
               dam2(m,Iz) = Noise_percent*temp_flux(m,Iz,NN_TIME)*randn
               temp_flux_r(m,Iz,1)       = (DET_SIG_PHI_G_Old(m,Iz,2))+dam2(m,Iz)
               temp_flux_r(m,Iz,NN_TIME) = temp_flux(m,Iz,NN_TIME)+dam2(m,Iz)

               x_kal(m,Iz,1)=0d0
               xm_FKH(m,Iz,NN_TIME)     = zm*x_kal(m,Iz,(NN_TIME-1))+temp_flux_r(m,Iz,(NN_TIME-1))
               X2m (m,Iz,1,1,NN_TIME)   = xm_FKH(m,Iz,NN_TIME)
               X2m (m,Iz,2,1,NN_TIME)   = temp_flux_r(m,Iz,(NN_TIME-1))
               X2m_temp(1,1) = X2m (m,Iz,1,1,NN_TIME)
               X2m_temp(2,1) = X2m (m,Iz,2,1,NN_TIME)

               y_FKH_temp          = MATMUL(CT,X2m_temp)
               y_FKH(m,Iz,NN_TIME) = Sv*y_FKH_temp(1,1)+dam2(m,Iz)*randn

               P_temp(1,1) = P(m,Iz,1,1,(NN_TIME-1))
               P_temp(1,2) = P(m,Iz,1,2,(NN_TIME-1))
               P_temp(2,1) = P(m,Iz,2,1,(NN_TIME-1))
               P_temp(2,2) = P(m,Iz,2,2,(NN_TIME-1))

               P_temp   = MATMUL(Ak,P_temp)
               Ak_trans = TRANSPOSE(Ak)
               P1_temp  = MATMUL(P_temp, Ak_trans)
               P2_temp  = P1_temp+Q

               Pm(m,Iz,1,1,NN_TIME) = P2_temp(1,1)
               Pm(m,Iz,1,2,NN_TIME) = P2_temp(1,2)
               Pm(m,Iz,2,1,NN_TIME) = P2_temp(2,1)
               Pm(m,Iz,2,2,NN_TIME) = P2_temp(2,2)

               CT_trans=TRANSPOSE(CT)
               Pm_temp(1,1) = Pm(m,Iz,1,1,NN_TIME)
               Pm_temp(1,2) = Pm(m,Iz,1,2,NN_TIME)
               Pm_temp(2,1) = Pm(m,Iz,2,1,NN_TIME)
               Pm_temp(2,2) = Pm(m,Iz,2,2,NN_TIME)

               A1      = MATMUL(CT,Pm_temp)
               A2      = MATMUL(A1,CT_trans)
               CTPMCT  = 1d0/(A2+R)
               A3      = MATMUL(Pm_temp,CT_trans)
               KN_temp = MATMUL(A3,CTPMCT)
               KN(m,Iz,1,1,NN_TIME) = KN_temp(1,1)
               KN(m,Iz,2,1,NN_TIME) = KN_temp(2,1)
               kesi_temp = y_FKH(m,Iz,NN_TIME)
               kesi_temp = kesi_temp-MATMUL(CT,X2m_temp)
               kesi_FKH(m,Iz,NN_TIME) = kesi_temp(1,1)
               x_temp = X2m_temp+MATMUL(KN_temp,kesi_temp)
               x_kal(m,Iz,NN_TIME) = x_temp(1,1)
               A4 = MATMUL(CT,Pm_temp)
               A5 = MATMUL(KN_temp,A4)

               P_temp = Pm_temp-A5
               P(m,Iz,1,1,NN_TIME) = P_temp(1,1)
               P(m,Iz,1,2,NN_TIME) = P_temp(1,2)
               P(m,Iz,2,1,NN_TIME) = P_temp(2,1)
               P(m,Iz,2,2,NN_TIME) = P_temp(2,2)
            ENDDO
         ENDDO
      ELSE IF (OPT_Det==1) THEN
         g     = DET_b_A2*DET_b_A2m/(DET_b_A2-DET_b_A2m)
         zg    = exp(-DET_b_A2*Tn_kal)
         zgm   = exp(-DET_b_A2m*Tn_kal)
         pp1   = Kpp
         qp    = 0.923d0*(1d0-pp1)
         rp    = 0.077d0*(1d0-pp1)
         alfa1 = (1d0/pp1)*(qp-(rp*g/DET_b_A2))*(1d0-zg)
         alfa2 = (1d0/pp1)*(rp*g/DET_b_A2)*(1d0-zgm)
         Do m = 1, DET_N_XY
            Do Iz = 1, DET_N_Z
               Ak_rh(1,1) = zg
               Ak_rh(1,2) = 0d0
               Ak_rh(1,3) = alfa1
               Ak_rh(2,1) = 0d0
               Ak_rh(2,2) = zgm
               Ak_rh(2,3) = alfa2
               Ak_rh(3,1) = 0d0
               Ak_rh(3,2) = 0d0
               Ak_rh(3,3) = 1d0

               CT_rh(1,1) = Kpp
               CT_rh(1,2) = Kpp
               CT_rh(1,3) = Kpp

               I_rh(1,1) = 1d0
               I_rh(1,2) = 0d0
               I_rh(1,3) = 0d0
               I_rh(2,1) = 0d0
               I_rh(2,2) = 1d0
               I_rh(2,3) = 0d0
               I_rh(3,1) = 0d0
               I_rh(3,2) = 0d0
               I_rh(3,3) = 1d0

               Q_rh(1,1) = 0.0082d0
               Q_rh(1,2) = 0d0
               Q_rh(1,3) = 0d0
               Q_rh(2,1) = 0d0
               Q_rh(2,2) = 0.0082d0
               Q_rh(2,3) = 0d0
               Q_rh(3,1) = 0d0
               Q_rh(3,2) = 0d0
               Q_rh(3,3) = 0.05d0

               P_rh(m,Iz,1,1,1)=2d0
               P_rh(m,Iz,1,2,1)=0d0
               P_rh(m,Iz,1,3,1)=0d0
               P_rh(m,Iz,2,1,1)=0d0
               P_rh(m,Iz,2,2,1)=2d0
               P_rh(m,Iz,2,3,1)=0d0
               P_rh(m,Iz,3,1,1)=0d0
               P_rh(m,Iz,3,2,1)=0d0
               P_rh(m,Iz,3,3,1)=2d0

               R=8.0
!               randn=random_normal()
            dummy_filler = 1 ! @$^ siarhei_plot 
               temp_flux(m,Iz,1)         = DET_SIG_PHI_G(m,Iz,2)
               temp_flux(m,Iz,NN_TIME)   = DET_SIG_PHI_G(m,Iz,2)
               dam2(m,Iz)                = Noise_percent*temp_flux(m,Iz,NN_TIME)*randn
               temp_flux_r(m,Iz,1)       = temp_flux(m,Iz,1)
               temp_flux_r(m,Iz,NN_TIME) = temp_flux(m,Iz,NN_TIME)+dam2(m,Iz)

               xg (m,Iz,1)       = 0d0
               xgm(m,Iz,1)       = 0d0
               xg(m,Iz,NN_TIME)  = zg*xg(m,Iz,(NN_TIME-1))+alfa1*temp_flux_r(m,Iz,(NN_TIME-1))
               xgm(m,Iz,NN_TIME) = zgm*xgm(m,Iz,(NN_TIME-1))+alfa2*temp_flux_r(m,Iz,(NN_TIME-1))

               X2m_rh(m,Iz,1,1,NN_TIME) = xg(m,Iz,NN_TIME)
               X2m_rh(m,Iz,2,1,NN_TIME) = xgm(m,Iz,NN_TIME)
               X2m_rh(m,Iz,3,1,NN_TIME) = temp_flux_r(m,Iz,(NN_TIME-1))
               X2m_temp_rh(1,1) = X2m_rh(m,Iz,1,1,NN_TIME)
               X2m_temp_rh(2,1) = X2m_rh(m,Iz,2,1,NN_TIME)
               X2m_temp_rh(3,1) = X2m_rh(m,Iz,3,1,NN_TIME)

               y_FKH_temp           = MATMUL(CT_rh,X2m_temp_rh)
               y_FKH (m,Iz,NN_TIME) = y_FKH_temp(1,1)+dam2(m,Iz)*randn

               P_temp_rh(1,1) = P_rh(m,Iz,1,1,(NN_TIME-1))
               P_temp_rh(1,2) = P_rh(m,Iz,1,2,(NN_TIME-1))
               P_temp_rh(1,3) = P_rh(m,Iz,1,3,(NN_TIME-1))
               P_temp_rh(2,1) = P_rh(m,Iz,2,1,(NN_TIME-1))
               P_temp_rh(2,2) = P_rh(m,Iz,2,2,(NN_TIME-1))
               P_temp_rh(2,3) = P_rh(m,Iz,2,3,(NN_TIME-1))
               P_temp_rh(3,1) = P_rh(m,Iz,3,1,(NN_TIME-1))
               P_temp_rh(3,2) = P_rh(m,Iz,3,2,(NN_TIME-1))
               P_temp_rh(3,3) = P_rh(m,Iz,3,3,(NN_TIME-1))

               P_temp_rh   = MATMUL(Ak_rh,P_temp_rh)
               Ak_rh_trans = TRANSPOSE(Ak_rh)
               P1_temp_rh  = MATMUL(P_temp_rh,Ak_rh_trans)
               P2_temp_rh  = P1_temp_rh+Q_rh

               Pm_rh(m,Iz,1,1,NN_TIME) = P2_temp_rh(1,1)
               Pm_rh(m,Iz,1,2,NN_TIME) = P2_temp_rh(1,2)
               Pm_rh(m,Iz,1,3,NN_TIME) = P2_temp_rh(1,3)
               Pm_rh(m,Iz,2,1,NN_TIME) = P2_temp_rh(2,1)
               Pm_rh(m,Iz,2,2,NN_TIME) = P2_temp_rh(2,2)
               Pm_rh(m,Iz,2,3,NN_TIME) = P2_temp_rh(2,3)
               Pm_rh(m,Iz,3,1,NN_TIME) = P2_temp_rh(3,1)
               Pm_rh(m,Iz,3,2,NN_TIME) = P2_temp_rh(3,2)
               Pm_rh(m,Iz,3,3,NN_TIME) = P2_temp_rh(3,3)

               CT_trans_rh = TRANSPOSE(CT_rh)

               Pm_temp_rh(1,1) = Pm_rh(m,Iz,1,1,NN_TIME)
               Pm_temp_rh(1,2) = Pm_rh(m,Iz,1,2,NN_TIME)
               Pm_temp_rh(1,3) = Pm_rh(m,Iz,1,3,NN_TIME)
               Pm_temp_rh(2,1) = Pm_rh(m,Iz,2,1,NN_TIME)
               Pm_temp_rh(2,2) = Pm_rh(m,Iz,2,2,NN_TIME)
               Pm_temp_rh(2,3) = Pm_rh(m,Iz,2,3,NN_TIME)
               Pm_temp_rh(3,1) = Pm_rh(m,Iz,3,1,NN_TIME)
               Pm_temp_rh(3,2) = Pm_rh(m,Iz,3,2,NN_TIME)
               Pm_temp_rh(3,3) = Pm_rh(m,Iz,3,3,NN_TIME)

               A1_rh  = MATMUL(CT_rh,Pm_temp_rh)
               A2     = MATMUL(A1_rh,CT_trans_rh)
               CTPMCT = 1/(A2+R)
               A3_rh  = MATMUL(Pm_temp_rh,CT_trans_rh)
               KN_temp_rh              = MATMUL(A3_rh,CTPMCT)
               KN_rh(m,Iz,1,1,NN_TIME) = KN_temp_rh(1,1)
               KN_rh(m,Iz,2,1,NN_TIME) = KN_temp_rh(2,1)
               kesi_temp              = y_FKH(m,Iz,NN_TIME)
               kesi_temp              = kesi_temp- MATMUL(CT_rh,X2m_temp_rh)
               kesi_FKH(m,Iz,NN_TIME) = kesi_temp(1,1)
               x_temp_rh           = X2m_temp_rh+MATMUL(KN_temp_rh,kesi_temp)
               x_kal(m,Iz,NN_TIME) = x_temp_rh(1,1)
               A4_rh = MATMUL(CT_rh,Pm_temp_rh)
               A5_rh = MATMUL(KN_temp_rh,A4_rh)

               P_temp_rh = P_temp_rh-A5_rh
               P_rh(m,Iz,1,1,NN_TIME) = P_temp_rh(1,1)
               P_rh(m,Iz,1,2,NN_TIME) = P_temp_rh(1,2)
               P_rh(m,Iz,1,3,NN_TIME) = P_temp_rh(1,3)
               P_rh(m,Iz,2,1,NN_TIME) = P_temp_rh(2,1)
               P_rh(m,Iz,2,2,NN_TIME) = P_temp_rh(2,2)
               P_rh(m,Iz,2,3,NN_TIME) = P_temp_rh(2,3)
               P_rh(m,Iz,3,1,NN_TIME) = P_temp_rh(3,1)
               P_rh(m,Iz,3,2,NN_TIME) = P_temp_rh(3,2)
               P_rh(m,Iz,3,3,NN_TIME) = P_temp_rh(3,3)
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE Compensate_Detector
#endif 


#ifdef siarhei_delete 
      FUNCTION random_normal()

      IMPLICIT NONE

      REAL(8) :: random_normal
      REAL(8) :: s  =  0.449871d0
      REAL(8) :: t  = -0.386595d0
      REAL(8) :: a  =  0.196000d0
      REAL(8) :: b  =  0.254720d0
      REAL(8) :: r1 =  0.275970d0
      REAL(8) :: r2 =  0.278460d0
      REAL(8) :: u, v, x, y, q
      REAL(8) :: half = 0.5d0

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [random_normal] in Mod_Detector'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      DO WHILE (.TRUE.)
         CALL RANDOM_NUMBER(u)
         CALL RANDOM_NUMBER(v)
         v = 1.7156 * (v - half)
         x = u - s
         y = ABS(v) - t
         q = x**2 + y*(a*y - b*x)
         IF (q < r1) EXIT
         IF (q > r2) CYCLE
         IF (v**2 < -4.0*LOG(u)*u**2) EXIT
      END DO

!      random_normal = v/u
            dummy_filler = 1 ! @$^ siarhei_plot 

      RETURN
      END FUNCTION random_normal
#endif 


#ifdef siarhei_delete 
      subroutine cal_det_pow
!      USE Mod_SSnodal, ONLY: NodalSolver_St
!      USE Mod_Pinpow, ONLY: PinPOW
      USE Inc_Control
      USE Inc_TH
      USE Inc_Depletion
      USE Inc_File
      USE Inc_3D
      USE Inc_PinPOW
      USE Inc_Time
      USE Inc_Utile
      USE Inc_XS_File, ONLY: I_BU
      USE Mod_Save
!      USE Mod_Depletion
      USE Mod_GetSome
      USE Inc_Detector
      USE Read_File, ONLY: Read_RI, read_mri
      USE Write_O
      USE Mod_XSFB
      use inc_crud, only:opt_crud,ind_bu,ind_bu_predictor
      use Inc_maXS, only: kap_maXS_f_3D, maXS_f_3D
      use Inc_Detector
      implicit none
      real(8) :: keff_Old_XSFB
      integer :: Ixy, Ixy_FA, Iz
      integer :: I_BRCH, I_XSFB
      real(8) :: mid_err
      INTEGER(4) ::  i, j, k, m, tmp1, tmp2, tmp_iz, tmp_ixy
      real(8)    :: tmp_num, tmp_mid, tmp_sum, tmp_sum_1, tmp_sum_2, tmp_sum_vol, tmp_vol
      real(8),allocatable,dimension(:,:,:) :: DET_NPOW_POW2
      real(8) :: tmp_det_vol

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [cal_det_pow] in Mod_Detector'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.flag_det_pow) return

      if (flag_incore_pow) then
         ! rr --> power
         tmp_sum   = 0d0
         tmp_num   = 0d0
         tmp_sum_1 = 0d0
         tmp_sum_2 = 0d0
         ! average maXS_f_3D
         Do j = 1, DET_NL
            Do i = 1, DET_NPOW
               Do m = 1,4
                  if (I_1Nto4N(DET_NPOW_Ixy_L(i,1),m) == 0) cycle
                  ixy = I_1Nto4N(DET_NPOW_Ixy_L(i,1),m)
                  Do k = DET_NPOW_Iz_L(j,1),DET_NPOW_Iz_L(j,2)
                     tmp_sum_1 = tmp_sum_1 + maXS_f_3D(ixy,k,1)
                     tmp_sum_2 = tmp_sum_2 + maXS_f_3D(ixy,k,2)
                     tmp_num = tmp_num + 1d0
                  Enddo
               Enddo
            Enddo
         Enddo
         tmp_sum_1 = tmp_sum_1 / tmp_num
         tmp_sum_2 = tmp_sum_2 / tmp_num
         ! power
         If (allocated(DET_NPOW_Ixy)) deallocate(DET_NPOW_Ixy)
         If (allocated(DET_NPOW_Iz))  deallocate(DET_NPOW_Iz)
         If (allocated(DET_NPOW_POW)) deallocate(DET_NPOW_POW)
         allocate(DET_NPOW_Ixy(int(tmp_num)))
         allocate(DET_NPOW_Iz (int(tmp_num)))
         allocate(DET_NPOW_POW(int(tmp_num)))
         DET_NPOW_Ixy=0
         DET_NPOW_Iz =0
         DET_NPOW_POW=0d0
         tmp1 = 0
         Do j = 1, DET_NL
            Do i = 1, DET_NPOW
               Do m = 1, 4
                  if (DET_NPOW_Ixy_L(i,1+m) .EQ. 0) cycle
                  Do k = DET_NPOW_Iz_L(j,1),DET_NPOW_Iz_L(j,2)
                     tmp1 = tmp1+1
                     DET_NPOW_Iz (tmp1)=k
                     DET_NPOW_Ixy(tmp1)=DET_NPOW_Ixy_L(i,1+m)
                     DET_NPOW_POW(tmp1)=DET_NPOW_POW_L(i,j)*tmp_sum_1*kap_maXS_f_3D(DET_NPOW_Ixy(tmp1),k,1)/maXS_f_3D(DET_NPOW_Ixy(tmp1),k,1)*flux(DET_NPOW_Ixy(tmp1),k,1)
                     DET_NPOW_POW(tmp1)=DET_NPOW_POW(tmp1)+DET_NPOW_POW_L(i,j)*tmp_sum_2*kap_maXS_f_3D(DET_NPOW_Ixy(tmp1),k,2)/maXS_f_3D(DET_NPOW_Ixy(tmp1),k,2)*flux(DET_NPOW_Ixy(tmp1),k,2)
                  Enddo
               Enddo
            Enddo
         Enddo
         DET_NPOW     = tmp1
         flag_det_pow_sol_once = .FALSE.
      endif
      flag_det_pow_sol = .true.
      if (flag_det_pow_sol_once) then
         if (.not.flag_det_pow_opt) then
            DET_NPOW = DET_NPOW_Old
            If (allocated(DET_NPOW_RATIO))    deallocate(DET_NPOW_RATIO)
            !! SET RATIO MATRIX
            tmp1 = 0
            tmp2 = 0
            Do i = 1, DET_NL
               tmp1 = DET_NPOW_Iz_L(i,2)-DET_NPOW_Iz_L(i,1)+1
               if (tmp2 < tmp1) tmp2 = tmp1
            Enddo
            tmp_det_vol = ((meshsize_x(1)**2d0)*40d0)*4d0
            Do j = 1, DET_NL
               Do i = 1, DET_NPOW
                  DET_NPOW_POW_L(i,j) = DET_NPOW_POW_L(i,j) / tmp_det_vol * 1d+6
               Enddo
            Enddo
            allocate(DET_NPOW_RATIO(DET_NPOW*4,DET_NL,tmp2))
            !! CALCULATE POWER and CALCULATE RATIO
            If (allocated(DET_NPOW_POW2)) deallocate(DET_NPOW_POW2)
            allocate(DET_NPOW_POW2(DET_NPOW*4,DET_NL,tmp2))
            tmp_sum = 0d0
            tmp1 = 0
            Do j = 1, DET_NL
               Do i = 1, DET_NPOW
                  tmp_sum = 0d0
                  tmp_num = 0d0
                  Do m = 1,4
                     if (I_1Nto4N(DET_NPOW_Ixy_L(i,1),m) == 0) cycle
                     tmp_num = tmp_num+1d0
                     DET_NPOW_Ixy_L(i,1+m) = I_1Nto4N(DET_NPOW_Ixy_L(i,1),m)
                     Do k = DET_NPOW_Iz_L(j,1),DET_NPOW_Iz_L(j,2)
                        DET_NPOW_RATIO(4*(i-1)+m,j,k-DET_NPOW_Iz_L(j,1)+1) = Power(DET_NPOW_Ixy_L(i,1+m),k)
                        tmp1 = tmp1+1
                        if (Power(DET_NPOW_Ixy_L(i,1+m),k)>0) then
                           tmp_sum = tmp_sum+Power(DET_NPOW_Ixy_L(i,1+m),k)*meshsize_z(k)
                        endif
                     Enddo
                  Enddo
                  tmp_mid = tmp_sum/tmp_num/(sum(meshsize_z(DET_NPOW_Iz_L(j,1):DET_NPOW_Iz_L(j,2))))
                  Do m = 1,4
                     if (I_1Nto4N(DET_NPOW_Ixy_L(i,1),m) == 0) cycle
                     Do k = DET_NPOW_Iz_L(j,1),DET_NPOW_Iz_L(j,2)
                        DET_NPOW_RATIO(4*(i-1)+m,j,k-DET_NPOW_Iz_L(j,1)+1) = Power(DET_NPOW_Ixy_L(i,1+m),k) / tmp_mid
                        DET_NPOW_POW2(4*(i-1)+m,j,k-DET_NPOW_Iz_L(j,1)+1) = DET_NPOW_RATIO(4*(i-1)+m,j,k-DET_NPOW_Iz_L(j,1)+1)*DET_NPOW_POW_L(i,j)
                     Enddo
                  Enddo
               Enddo
            Enddo
            If (allocated(DET_NPOW_Ixy)) deallocate(DET_NPOW_Ixy)
            If (allocated(DET_NPOW_Iz))  deallocate(DET_NPOW_Iz)
            If (allocated(DET_NPOW_POW)) deallocate(DET_NPOW_POW)
            allocate(DET_NPOW_Ixy(tmp1))
            allocate(DET_NPOW_Iz (tmp1))
            allocate(DET_NPOW_POW(tmp1))
            DET_NPOW_Ixy=0
            DET_NPOW_Iz =0
            DET_NPOW_POW=0d0
            tmp1 = 0
            Do j = 1, DET_NL
               Do i = 1, DET_NPOW
                  Do m = 1, 4
                     if (DET_NPOW_Ixy_L(i,1+m) .EQ. 0) cycle
                     Do k = DET_NPOW_Iz_L(j,1),DET_NPOW_Iz_L(j,2)
                        tmp1 = tmp1+1
                        DET_NPOW_Iz (tmp1)=k
                        DET_NPOW_Ixy(tmp1)=DET_NPOW_Ixy_L(i,1+m)
                        DET_NPOW_POW(tmp1)=DET_NPOW_POW2(4*(i-1)+m,j,k-DET_NPOW_Iz_L(j,1)+1)
                     Enddo
                  Enddo
               Enddo
            Enddo
            DET_NPOW     = tmp1
            deallocate(DET_NPOW_POW2)
            flag_det_pow_sol_once = .FALSE.
         else
            DET_NPOW = DET_NPOW_Old
            If (allocated(DET_NPOW_RATIO))    deallocate(DET_NPOW_RATIO)
            !! SET RATIO MATRIX
            tmp1 = 0
            tmp2 = 0
            Do i = 1, DET_NL
               tmp1 = DET_NPOW_Iz_L(i,2)-DET_NPOW_Iz_L(i,1)+1
               if (tmp2 < tmp1) tmp2 = tmp1
            Enddo
            tmp_det_vol = ((meshsize_x(1)**2d0)*40d0)*4d0
            Do j = 1, DET_NL
               Do i = 1, DET_NPOW
                  DET_NPOW_POW_L(i,j) = DET_NPOW_POW_L(i,j) / tmp_det_vol * 1d+6
               Enddo
            Enddo
            allocate(DET_NPOW_RATIO(DET_NPOW*4,DET_NL,tmp2))
            !! CALCULATE POWER and CALCULATE RATIO
            If (allocated(DET_NPOW_POW2)) deallocate(DET_NPOW_POW2)
            allocate(DET_NPOW_POW2(DET_NPOW*4,DET_NL,tmp2))
            tmp_sum = 0d0
            tmp1 = 0
            Do j = 1, DET_NL
               Do i = 1, DET_NPOW
                  tmp_sum = 0d0
                  tmp_num = 0d0
                  Do m = 1,4
                     if (I_1Nto4N(DET_NPOW_Ixy_L(i,1),m) == 0) cycle
                     tmp_num = tmp_num+1d0
                     DET_NPOW_Ixy_L(i,1+m) = I_1Nto4N(DET_NPOW_Ixy_L(i,1),m)
                     Do k = DET_NPOW_Iz_L(j,1),DET_NPOW_Iz_L(j,2)
                        DET_NPOW_RATIO(4*(i-1)+m,j,k-DET_NPOW_Iz_L(j,1)+1) = Power(DET_NPOW_Ixy_L(i,1+m),k)
                        tmp1 = tmp1+1
                        if (Power(DET_NPOW_Ixy_L(i,1+m),k)>0) then
                           tmp_sum = tmp_sum+Power(DET_NPOW_Ixy_L(i,1+m),k)*meshsize_z(k)
                        endif
                     Enddo
                  Enddo
                  tmp_mid = tmp_sum/tmp_num/(sum(meshsize_z(DET_NPOW_Iz_L(j,1):DET_NPOW_Iz_L(j,2))))
                  Do m = 1,4
                     if (I_1Nto4N(DET_NPOW_Ixy_L(i,1),m) == 0) cycle
                     Do k = DET_NPOW_Iz_L(j,1),DET_NPOW_Iz_L(j,2)
                        DET_NPOW_RATIO(4*(i-1)+m,j,k-DET_NPOW_Iz_L(j,1)+1) = Power(DET_NPOW_Ixy_L(i,1+m),k) / tmp_mid
                        DET_NPOW_POW2(4*(i-1)+m,j,k-DET_NPOW_Iz_L(j,1)+1) = DET_NPOW_POW_L(i,j)
                     Enddo
                  Enddo
               Enddo
            Enddo
            If (allocated(DET_NPOW_Ixy)) deallocate(DET_NPOW_Ixy)
            If (allocated(DET_NPOW_Iz))  deallocate(DET_NPOW_Iz)
            If (allocated(DET_NPOW_POW)) deallocate(DET_NPOW_POW)
            allocate(DET_NPOW_Ixy(tmp1))
            allocate(DET_NPOW_Iz (tmp1))
            allocate(DET_NPOW_POW(tmp1))
            DET_NPOW_Ixy=0
            DET_NPOW_Iz =0
            DET_NPOW_POW=0d0
            tmp1 = 0
            Do j = 1, DET_NL
               Do i = 1, DET_NPOW
                  Do m = 1, 4
                     if (DET_NPOW_Ixy_L(i,1+m) .EQ. 0) cycle
                     Do k = DET_NPOW_Iz_L(j,1),DET_NPOW_Iz_L(j,2)
                        tmp1 = tmp1+1
                        DET_NPOW_Iz (tmp1)=k
                        DET_NPOW_Ixy(tmp1)=DET_NPOW_Ixy_L(i,1+m)
                        DET_NPOW_POW(tmp1)=DET_NPOW_POW2(4*(i-1)+m,j,k-DET_NPOW_Iz_L(j,1)+1)
                     Enddo
                  Enddo
               Enddo
            Enddo
            DET_NPOW     = tmp1
            deallocate(DET_NPOW_POW2)
            flag_det_pow_sol_once = .FALSE.
         endif
      endif
      I_BRCH = 0
      BoronSearch_det: do while (.TRUE.)
         I_BRCH = I_BRCH + 1
         I_XSFB = 0
         Flag_Conv_Critical = .TRUE.
         Flag_Conv_PowLev   = .TRUE.
         keff_OOld = keff
         keff_Old  = keff
         XSFeedback_det: do while (.TRUE.)
            I_XSFB = I_XSFB + 1
            keff_Old_XSFB = Keff

            call XSFB
!            call NodalSolver_St
            dummy_filler = 1 ! @$^ siarhei_plot 
            DET_P_ERR = 10d-30

            Do i = 1,DET_NPOW
               DO Ixy = 1,Nxy
                  DO Iz = 1,NZ
                     if ((Ixy == DET_NPOW_Ixy(i)) .and. (Iz==DET_NPOW_Iz(i))) then
                        mid_err = (power(Ixy,Iz)-DET_NPOW_POW(i))/DET_NPOW_POW(i)*100
                        if (abs(DET_P_ERR) < abs(mid_err)) DET_P_ERR = abs(mid_err)
                     endif
                  ENDDO
               ENDDO
            ENDDO

            tmp1 = 0
            Do j = 1, DET_NL
               Do i = 1, DET_NPOW_Old
                  tmp_sum = 0d0
                  tmp_sum_vol = 0d0
                  Do m = 1, 4
                     if (DET_NPOW_Ixy_L(i,1+m) .EQ. 0) cycle
                     Do k = DET_NPOW_Iz_L(j,1),DET_NPOW_Iz_L(j,2)
                        tmp_vol=0d0
                        tmp1 = tmp1+1
                        tmp_iz  = int(DET_NPOW_Iz(tmp1) )
                        tmp_ixy = int(DET_NPOW_Ixy(tmp1))
                        tmp_vol = (meshsize_x(1)**2d0)*meshsize_z(tmp_iz)/(1d+6)
                        tmp_sum_vol = tmp_sum_vol+tmp_vol
                        tmp_sum = tmp_sum+power(tmp_ixy,tmp_iz)*tmp_vol
                     Enddo
                  Enddo
                  tmp_sum=tmp_sum/tmp_sum_vol*(meshsize_x(1)**2d0)*4d0*40d0*(1d+6)
               Enddo
            Enddo

            if (DET_P_ERR < 0.03d0) EXIT XSFeedback_det

            if (I_XSFB>1) then
               stop 'RAST-K does not come here 1'
            endif
            if (I_XSFB==I_XSFB_Max) then
               EXIT ! Reached to Maximum XS Feedback Iteration Number
            endif
            stop 'RAST-K does not come here 2'
         enddo XSFeedback_det

         call Get_Avg(Avg_Power, Power, 0)
         do Ixy_FA=1,Nxy_FA
            Ixy=I_FA(Ixy_FA)
            do Iz=IzFuelBot,IzFuelTop
               Normal_Power(Ixy,Iz)=Power(Ixy,Iz)/Avg_Power
            enddo
         enddo

         if ((I_BU>1).or.(OPT_RST==2)) then
            do Ixy = 1, Nxy
               if (I_FARF_1N(I_4Nto1N(Ixy))==1) CYCLE
               do Iz=IzFuelBot,IzFuelTop
                  BU(Ixy,Iz)=BU_Predictor(Ixy,Iz)+Normal_Power(Ixy,Iz)*(Half*dBU)

               enddo
            enddo
         endif

         if (DABS(LevFactor-D1)>EPS_PowLev) Flag_Conv_PowLev=.FALSE.
         if (DABS(keff-k_target)>EPS_Critical) Flag_Conv_Critical=.FALSE.
         write(*,*) 'keff-old',keff_Old, 'keff-k',keff, 'k_target',k_target
         if (I_BRCH .EQ. 3) stop
         if ((Flag_Conv_Critical).and.(Flag_Conv_PowLev)) EXIT BoronSearch_det
         if ((I_BRCH==I_CriSearch_Max).or.(OPT_Mode==4)) EXIT BoronSearch_det
         ! == check the power condition == !!
         stop 'RAST-K does not come here 3'
      enddo BoronSearch_det
      flag_det_pow_sol = .false.

      RETURN
      end subroutine cal_det_pow
#endif 


      END MODULE Mod_Detector

