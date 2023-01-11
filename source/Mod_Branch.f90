
      MODULE Mod_Branch

      USE Inc_Constant
      USE Inc_File
      USE Inc_Flag
      USE Inc_Geometry
      USE Inc_RP
      USE Inc_Option
!      use mod_litever

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

      CONTAINS

      subroutine cal_branch(i_buu)
      use inc_3d
      use inc_branch
      use inc_option
      use inc_flag
      use inc_th
      use Mod_THdriver
      use inc_cr
      use inc_depletion
      use write_o !, only: write_simulator, write_simulator_axial
      use Mod_GetSome
      use inc_branch, only: flag_coef
      implicit none
      integer(4),intent(in) :: i_buu
      integer(4) :: ic, i_c, i_p
      integer(4) :: i_step1=0
      integer(4) :: i_step2=0
      integer(4) :: i_step3=0
      integer(4) :: i_step4=0
      integer(4) :: i_step5=0
      integer(4) :: i_s6=0
      real(8) :: keff_save, scram_worth, tpd1

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [cal_branch] in Mod_Branch'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [cal_branch] in Mod_Branch'
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

      if (.not.flag_branch) return

      flag_doing_branch=.true.
!      call backup_data
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif 
      do ic=1,n_case
         if(.not.branch_case(ic)%if_step(i_buu)) cycle

         if(branch_case(ic)%if_refresh) then
!            call reset_data
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif
         endif

         !write(*,'(a,i0," ",a)') 'Case# =',ic,trim(branch_case(ic)%case_name)

         !c_coef=branch_case(ic)%coeftype
         !c_case=ic

         if (branch_case(ic)%if_macro_xs) then
!            call MacroXs_Table_alloc
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif
!            call Micro_Xstable_2_Macro_Xstable
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif
         endif

         if(branch_case(ic)%if_flag_tffb  ) flag_tffb   = branch_case(ic)%flag_tffb
         if(branch_case(ic)%if_flag_tmfb  ) flag_tmfb   = branch_case(ic)%flag_tmfb
         if(branch_case(ic)%if_flag_xsfb  ) flag_xsfb   = branch_case(ic)%flag_xsfb
         if(branch_case(ic)%if_flag_thfb  ) flag_thfb   = branch_case(ic)%flag_thfb
         if(branch_case(ic)%if_opt_mode   ) opt_mode    = branch_case(ic)%opt_mode
         if(branch_case(ic)%if_opt_xe     ) opt_xe      = branch_case(ic)%opt_xe
         if(branch_case(ic)%if_opt_sm     ) opt_sm      = branch_case(ic)%opt_sm
         if(branch_case(ic)%if_opt_gd     ) opt_gd      = branch_case(ic)%opt_gd
         if(branch_case(ic)%if_ppm        ) ppm         = branch_case(ic)%ppm
         if(branch_case(ic)%if_ppower     ) ppower      = branch_case(ic)%ppower
         if(branch_case(ic)%if_opt_findtavg) opt_findtavg = branch_case(ic)%opt_findtavg

         if(branch_case(ic)%if_core_massflow) then
            core_massflow = branch_case(ic)%core_massflow
            fa_massflow   = core_massflow * chanvf_save / Core_N_FA
         endif

         if(branch_case(ic)%if_tm_in      ) then
            tm_in       = branch_case(ic)%tm_in
            t_mod(:,:)  = tm_in + DegToK
!            call nested_cal_dmod
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif
!            call restore_tm_in
#ifdef siarhei_plot
            dummy_filler = 1 ! @$^ siarhei_plot 
#endif
         endif

         flag_coef  = branch_case(ic)%if_czp
         if(branch_case(ic)%if_czp        ) then
            t_fuel(:,:)       =  t_mod(:,:)
         endif

         !if(branch_case(ic)%if_hzp    .and. branch_case(ic)%if_ppower) then
         !    t_fuel(:,:)      =  (1- branch_case(ic)%ppower )*t_mod(:,:)  + branch_case(ic)%ppower  * t_fuel(:,:)
         !endif

         if(branch_case(ic)%if_tf_in      ) then
            tf_in       = branch_case(ic)%tf_in
         endif
         if(branch_case(ic)%if_d_TF      ) then
            t_fuel(:,:) = t_fuel(:,:) + branch_case(ic)%d_TF
         endif

         if(branch_case(ic)%if_d_TM       ) then
            t_mod(:,:)  = t_mod(:,:)  + branch_case(ic)%d_TM
            tm_in       = tm_in       + branch_case(ic)%d_TM
!            call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
         endif

         if(branch_case(ic)%if_cr_bot     ) then
            cr_bot(:)= branch_case(ic)%cr_bot(:)
            CR_top(:)= CR_bot(:)+length_CR(:)
         endif

         write(*,'(a)') '* Branch case summary'
         write(*,'(a,f12.3)')                                         '- Cycle burnup :', cycle_bu(i_buu)

         if(branch_case(ic)%if_flag_xsfb  )   write(*,'(a,L12)')      '- flag_xsfb    :', branch_case(ic)%flag_xsfb
         if(branch_case(ic)%if_flag_thfb  )   write(*,'(a,L12)')      '- flag_thfb    :', branch_case(ic)%flag_thfb
         if(branch_case(ic)%if_flag_tffb  )   write(*,'(a,L12)')      '- flag_tffb    :', branch_case(ic)%flag_tffb
         if(branch_case(ic)%if_flag_tmfb  )   write(*,'(a,L12)')      '- flag_tmfb    :', branch_case(ic)%flag_tmfb
         if(branch_case(ic)%if_opt_mode   )   write(*,'(a,I12)')      '- opt_mode     :', branch_case(ic)%opt_mode
         if(branch_case(ic)%if_opt_xe     )   write(*,'(a,I12)')      '- opt_xe       :', branch_case(ic)%opt_xe
         if(branch_case(ic)%if_opt_sm     )   write(*,'(a,I12)')      '- opt_sm       :', branch_case(ic)%opt_sm
         if(branch_case(ic)%if_opt_gd     )   write(*,'(a,I12)')      '- opt_gd       :', branch_case(ic)%opt_gd
         if(branch_case(ic)%if_ppm        )   write(*,'(a,f12.3)')    '- ppm          :', branch_case(ic)%ppm
         if(branch_case(ic)%if_ppower     )   write(*,'(a,f12.3)')    '- ppower       :', branch_case(ic)%ppower
         if(branch_case(ic)%if_core_massflow) write(*,'(a,es12.5)')   '- core_massflow:', branch_case(ic)%core_massflow
         if(branch_case(ic)%if_tm_in      )   write(*,'(a,f12.3)')    '- tm_in        :', branch_case(ic)%tm_in
         if(branch_case(ic)%if_tf_in      )   then
            write(*,'(a,f12.3)')    '- tf_in        :', branch_case(ic)%tf_in
         elseif(branch_case(ic)%if_czp      ) then
            branch_case(ic)%tf_in=branch_case(ic)%tm_in
            write(*,'(a,f12.3)')    '- tf_in = tm_in:', branch_case(ic)%tf_in
         end if
         if(branch_case(ic)%if_d_TF       )   write(*,'(a,f12.3)')    '- d_TF         :', branch_case(ic)%d_TF
         if(branch_case(ic)%if_d_TM       )   write(*,'(a,f12.3)')    '- d_TM         :', branch_case(ic)%d_TM
         if(branch_case(ic)%if_cr_bot     )   write(*,'(a,100f12.3)') '- cr_bot       :', branch_case(ic)%cr_bot

         if(branch_case(ic)%if_ndr == 1 ) then
            call alloc(branch_case(ic)%ndr2 , branch_case(ic)%n_ppm, branch_case(ic)%n_mt)
            do i_p = 1, branch_case(ic)%n_ppm
               ppm = branch_case(ic)%ndr_ppm(i_p)
               do i_c = 1, branch_case(ic)%n_mt
                  flag_thfb = .false.
                  flag_tffb = .false.
                  flag_tmfb = .false.
                  ppower      = branch_case(ic)%ndr_ppower(i_c)
                  tm_in       = branch_case(ic)%ndr_mt(i_c)
                  t_mod(:,:)  = tm_in + DegToK
                  if(branch_case(ic)%if_czp        ) then
                     t_fuel(:,:) = t_mod(:,:)
                  endif
                  t_mod(:,:)  = t_mod(:,:)  +  branch_case(ic)%d_TM1(1)
                  tm_in       = tm_in       +  branch_case(ic)%d_TM1(1)
!                  call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
                  write(*, '(A6, f12.3 , A6, f12.3)') 'tm_in = ' , tm_in ,'ppm = ' , ppm

!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  keff_save = keff

                  t_mod(:,:)  = t_mod(:,:)  -  branch_case(ic)%d_TM1(1)+ branch_case(ic)%d_TM1(2)
                  tm_in       = tm_in       -  branch_case(ic)%d_TM1(1)+ branch_case(ic)%d_TM1(2)
!                  call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
                  write(*, '(A6, f12.3 , A6, f12.3)') 'tm_in = ' , tm_in ,'ppm = ' , ppm
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 

                  branch_case(ic)%ndr2(i_p,i_c) = 100000 * LOG(keff_save/keff) / (branch_case(ic)%d_TM1(1) - branch_case(ic)%d_TM1(2))
                  t_mod(:,:)  = t_mod(:,:)  - branch_case(ic)%d_TM1(2)
                  tm_in       = tm_in       - branch_case(ic)%d_TM1(2)
!                  call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 

               enddo
            enddo

         elseif(branch_case(ic)%if_ndr == 2 ) then
            i_step2 = i_step2 + 1
            do i_p = 1, branch_case(ic)%n_ppm
               flag_thfb = .false.
               flag_tffb = .false.
               flag_tmfb = .false.
               if (branch_case(ic)%opt_mtc==1) then
                  if ( allocated(cr_bot) ) then
                     cr_bot(:)= branch_case(ic)%cr_bot(:)
                     CR_top(:)= CR_bot(:)+length_CR(:)
                  endif
                  ppm = branch_case(ic)%ndr_ppm(i_p)
                  t_mod(:,:)  = t_mod(:,:)  + branch_case(ic)%d_TM
                  tm_in       = tm_in       + branch_case(ic)%d_TM
!                  call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
                  write(*, '(A6, f12.3)') 'tm_in = ' , tm_in
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  keff_save = keff
                  t_mod(:,:)  = t_mod(:,:)  - 2*branch_case(ic)%d_TM
                  tm_in       = tm_in       - 2*branch_case(ic)%d_TM
!                  call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
                  write(*, '(A6, f12.3)') 'tm_in = ' , tm_in
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  branch_case(ic)%ndr2(i_p ,i_step2 ) = 100000 * LOG(keff_save/keff) / (2*branch_case(ic)%d_TM) ! IL MTC_PPM
               endif
               if ( allocated(cr_bot) ) then
                  cr_bot(:) = length_cr_tip(:)+length_cr(:)
                  cr_top(:) = cr_bot(:)+length_cr(:)
               endif
               t_mod(:,:)  = t_mod(:,:)  + 2*branch_case(ic)%d_TM
               tm_in       = tm_in       + 2*branch_case(ic)%d_TM
!               call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
               write(*, '(A6, f12.3)') 'tm_in = ' , tm_in
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               keff_save = keff
               t_mod(:,:)  = t_mod(:,:)  - 2*branch_case(ic)%d_TM
               tm_in       = tm_in       - 2*branch_case(ic)%d_TM
!               call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
               write(*, '(A6, f12.3)') 'tm_in = ' , tm_in
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               branch_case(ic)%ndr2_1(i_p ,i_step2 ) = 100000 * LOG(keff_save/keff) / (2*branch_case(ic)%d_TM) ! ARO MTC_PPM
               t_mod(:,:)  = t_mod(:,:)  + branch_case(ic)%d_TM
               tm_in       = tm_in       + branch_case(ic)%d_TM
!               call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 

            enddo

         elseif(branch_case(ic)%if_ndr == 3 ) then
            flag_thfb = .true.
            flag_tffb = .true.
            flag_tmfb = .true.
            t_mod(:,:)  = t_mod(:,:)  + branch_case(ic)%d_TM
            tm_in       = tm_in       + branch_case(ic)%d_TM
!            call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
            write(*, '(A6, f12.3)') 'tm_in = ' , tm_in
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            keff_save = keff

            flag_thfb = .false.
            flag_tffb = .false.
            flag_tmfb = .false.
            t_mod(:,:)  = t_mod(:,:)  - 2*branch_case(ic)%d_TM
            tm_in       = tm_in       - 2*branch_case(ic)%d_TM
!            call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
            write(*, '(A6, f12.3)') 'tm_in = ' , tm_in
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            branch_case(ic)%ndr1(i_buu) = 100000 * LOG(keff_save/keff) / (2*branch_case(ic)%d_TM)
            t_mod(:,:)  = t_mod(:,:)  + branch_case(ic)%d_TM
            tm_in       = tm_in       + branch_case(ic)%d_TM
!            call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 

         elseif(branch_case(ic)%if_ndr == 4 ) then
            i_step1 = i_step1 + 1
            flag_thfb = .true.
            flag_tffb = .true.
            flag_tmfb = .true.
            t_fuel(:,:) = t_mod(:,:)
            ppower = 1e-4
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            do i_p = 1 , branch_case(ic)%n_ppower
               ppower = branch_case(ic)%ndr_ppower(i_p)
               flag_thfb = .true.
               flag_tffb = .true.
               flag_tmfb = .true.
               t_fuel(:,:)  = t_fuel(:,:)  + branch_case(ic)%d_TM1(1)
               call Get_Avg( branch_case(ic)%TF_Avg(i_buu), T_Fuel, 0 )
               write(*, '(A6, f12.3)') 'tf_in = ' , branch_case(ic)%TF_Avg(i_buu) + branch_case(ic)%d_TF
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               flag_thfb = .false.
               flag_tffb = .true.
               flag_tmfb = .false.
               keff_save = keff
               t_fuel(:,:)  = t_fuel(:,:) - branch_case(ic)%d_TM1(1) + branch_case(ic)%d_TM1(2)
               write(*, '(A6, f12.3)') 'tf_in = ' , branch_case(ic)%TF_Avg(i_buu) - branch_case(ic)%d_TF
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               branch_case(ic)%ndr2(i_p , i_step1) = 100000 * LOG(keff_save/keff) / (branch_case(ic)%d_TM1(1) - branch_case(ic)%d_TM1(2))
               write(*,*) branch_case(ic)%ndr2(i_p , i_step1) ,ppower
               t_fuel(:,:)  = t_fuel(:,:)  - branch_case(ic)%d_TM1(2)
            enddo
         elseif(branch_case(ic)%if_ndr == 5 ) then
            call alloc(branch_case(ic)%ndr2 , branch_case(ic)%n_ppm,branch_case(ic)%n_mt)
            call alloc(branch_case(ic)%ndr_ft , branch_case(ic)%n_mt)
            flag_thfb = .false.
            flag_tffb = .false.
            flag_tmfb = .false.
            do i_p = 1, branch_case(ic)%n_ppm
               ppm = branch_case(ic)%ndr_ppm(i_p)
               do i_c = 1, branch_case(ic)%n_mt
                  ppower      = branch_case(ic)%ndr_ppower(i_c)
                  tm_in       = branch_case(ic)%ndr_mt(i_c)
                  t_mod(:,:)  = tm_in + DegToK
                  if(branch_case(ic)%if_czp        ) then
                     t_fuel(:,:)       =  t_mod(:,:)

                  endif
                  t_mod(:,:)  = t_mod(:,:)  + branch_case(ic)%d_TM1(1)
                  t_fuel(:,:) = t_fuel(:,:) + branch_case(ic)%d_TM1(1)
!                  call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  call Get_Avg( branch_case(ic)%TF_Avg(i_buu), T_Fuel, 0 )
                  branch_case(ic)%NDR_FT(i_c) = branch_case(ic)%TF_Avg(i_buu)
                  write(*, '(A10, f12.3, A10, f12.3 , A10, f12.3)') 'tm_in = ' , tm_in  ,'tf_avg = ' ,branch_case(ic)%TF_Avg(i_buu)  , ' ppm  = '  , ppm
                  keff_save = keff
                  t_mod(:,:)  = t_mod(:,:)  - branch_case(ic)%d_TM1(1) + branch_case(ic)%d_TM1(2)
                  t_fuel(:,:) = t_fuel(:,:) - branch_case(ic)%d_TM1(1) + branch_case(ic)%d_TM1(2)
!                  call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  call Get_Avg( branch_case(ic)%TF_Avg(i_buu), T_Fuel, 0 )
                  write(*, '(A10, f12.3, A10, f12.3 , A10, f12.3)') 'tm_in = ' , tm_in  ,'tf_avg = ' ,branch_case(ic)%TF_Avg(i_buu)  , ' ppm  = '  , ppm
                  branch_case(ic)%ndr2(i_p,i_c) = 100000 * LOG(keff_save/keff) / (branch_case(ic)%d_TM1(1) - branch_case(ic)%d_TM1(2))
                  t_mod(:,:)  = t_mod(:,:)  - branch_case(ic)%d_TM1(2)
                  t_fuel(:,:) = t_fuel(:,:) - branch_case(ic)%d_TM1(2)
!                  call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 

               enddo
            enddo

         elseif(branch_case(ic)%if_ndr == 6 ) then
            opt_mode = 1
            call alloc(branch_case(ic)%ndr2 , branch_case(ic)%n_ppm,branch_case(ic)%n_mt)
            call alloc(branch_case(ic)%ndr_ft , branch_case(ic)%n_mt)
            do i_p = 1, branch_case(ic)%n_ppm
               ppm = branch_case(ic)%ndr_ppm(i_p)
               do i_c = 1, branch_case(ic)%n_mt
                  flag_thfb = .true.
                  flag_tffb = .true.
                  flag_tmfb = .true.
                  ppower      = branch_case(ic)%ndr_ppower(i_c)
                  tm_in       = branch_case(ic)%ndr_mt(i_c)
                  t_mod(:,:)  = tm_in + DegToK
                  if(branch_case(ic)%if_czp        ) then
                     t_fuel(:,:)       =  t_mod(:,:)
                  endif
                  ppm = ppm + branch_case(ic)%d_TM1(1)
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  call Get_Avg( branch_case(ic)%TF_Avg(i_buu), T_Fuel, 0 )

                  keff_save = keff

                  ppm = ppm - branch_case(ic)%d_TM1(1) + branch_case(ic)%d_TM1(2)

!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  call Get_Avg( branch_case(ic)%TF_Avg(i_buu), T_Fuel, 0 )
                  branch_case(ic)%ndr2(i_p,i_c) = 100000 * LOG(keff_save/keff) / (branch_case(ic)%d_TM1(1) - branch_case(ic)%d_TM1(2))

                  ppm = ppm - branch_case(ic)%d_TM1(2)
              enddo
           enddo
         elseif(branch_case(ic)%if_ndr == 7 ) then
            flag_thfb = .true.
            flag_tffb = .true.
            flag_tmfb = .true.
            opt_mode = 1
            ppm = ppm + 10
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            call Get_Avg( branch_case(ic)%TF_Avg(i_buu), T_Fuel, 0 )
            write(*, '(A10, f12.3, A10, f12.3 , A10, f12.3)') 'tm_in = ' , tm_in  ,'tf_avg = ' ,branch_case(ic)%TF_Avg(i_buu)  , ' ppm  = '  , ppm
            keff_save = keff
            ppm = ppm - 10
            flag_thfb = .false.
            flag_tffb = .false.
            flag_tmfb = .false.
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            call Get_Avg( branch_case(ic)%TF_Avg(i_buu), T_Fuel, 0 )
            write(*, '(A10, f12.3, A10, f12.3 , A10, f12.3)') 'tm_in = ' , tm_in  ,'tf_avg = ' ,branch_case(ic)%TF_Avg(i_buu)  , ' ppm  = '  , ppm
            branch_case(ic)%ndr1(i_buu) = 100000 * LOG(keff_save/keff) /10

         elseif(branch_case(ic)%if_ndr == 8 ) then
            opt_mode = 2
            flag_thfb = .true.
            flag_tffb = .true.
            flag_tmfb = .true.
            call alloc(branch_case(ic)%ndr1 , branch_case(ic)%n_ppower)
            do i_p = 1 , branch_case(ic)%n_ppower
               ppower = branch_case(ic)%NDR_PPOWER(i_p)
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               call Get_Avg( branch_case(ic)%TF_Avg(i_buu), T_Fuel, 0 )
               branch_case(ic)%ndr1(i_p) =  branch_case(ic)%TF_Avg(i_buu)
            enddo

        ! elseif(branch_case(ic)%if_ndr == 9 ) then
        !    opt_mode = 1
        !    call alloc(branch_case(ic)%ndr2 , branch_case(ic)%n_ppm, branch_case(ic)%n_ppower)
        !    do i_p = 1, branch_case(ic)%n_ppm
        !       ppm = branch_case(ic)%ndr_ppm(i_p)
        !       ppower = 1.0
        !       flag_thfb = .false.
        !       flag_tffb = .true.
        !       flag_tmfb = .true.
!        !       call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
        !       ppower = 0.5
        !       flag_thfb = .false.
        !       flag_tffb = .true.
        !       flag_tmfb = .true.
!        !       call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
        !       ppower = 1e-4
        !       flag_thfb = .false.
        !       flag_tffb = .true.
        !       flag_tmfb = .true.
!        !       call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
        !       flag_thfb = .true.
        !       flag_tffb = .true.
        !       flag_tmfb = .true.
!        !       call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
        !       keff_save = keff
        !       ppower      = branch_case(ic)%ndr_ppower(1)
!        !       call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
        !       branch_case(ic)%ndr2(i_p,1) = 0.0
        !       do i_c = 2, branch_case(ic)%n_ppower
        !          ppower      = branch_case(ic)%ndr_ppower(i_c)
!        !          call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
        !          branch_case(ic)%ndr2(i_p,i_c) = 100000 * LOG(keff/keff_save)
        !       enddo
        !    enddo

         elseif(branch_case(ic)%if_ndr == 9 ) then
            opt_mode = 1
            call alloc(branch_case(ic)%ndr2 , branch_case(ic)%n_ppm, branch_case(ic)%n_ppower)
            do i_p = 1, branch_case(ic)%n_ppm
               ppm = branch_case(ic)%ndr_ppm(i_p)
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               keff_save = keff
               flag_thfb = .true.
               flag_tffb = .true.
               flag_tmfb = .true.
               ppower = 1.0
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               tpd1 = 100000 * LOG(keff_save/keff)
               branch_case(ic)%ndr2(i_p,1) = 0
               do i_c = 2, branch_case(ic)%n_ppower
                  ppower      = branch_case(ic)%ndr_ppower(i_c)
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  branch_case(ic)%ndr2(i_p,i_c) = 100000 * LOG(keff_save/keff) - tpd1
               enddo
            enddo

         elseif(branch_case(ic)%if_ndr == 10 ) then !total power coefficient
            opt_mode = 1
            call alloc(branch_case(ic)%ndr2 , branch_case(ic)%n_ppm, branch_case(ic)%n_ppower)
            do i_p = 1, branch_case(ic)%n_ppm
               ppm = branch_case(ic)%ndr_ppm(i_p)
               ppower = 1.0
               flag_thfb = .true.
               flag_tffb = .true.
               flag_tmfb = .true.
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               ppower = 0.5
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               ppower = 1e-4
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  do i_c = 1, branch_case(ic)%n_ppower
                  ppower      = branch_case(ic)%ndr_ppower(i_c)
                  flag_thfb = .true.
                  flag_tffb = .true.
                  flag_tmfb = .true.
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  keff_save = keff
                  ppower = ppower + 0.01  !USNRC HRTD Rev 1208  2.1-15
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  branch_case(ic)%ndr2(i_p,i_c) = 100000 * LOG(keff/keff_save)
               enddo
            enddo

         elseif(branch_case(ic)%if_ndr == 11 ) then !doppler power defect
            opt_mode = 1
            i_step4 =  i_step4 + 1
            do i_p = 1, branch_case(ic)%n_ppower
               if (i_p == 1 ) then
                  ppower = 1.0
                  flag_thfb = .false.
                  flag_tffb = .true.
                  flag_tmfb = .true.
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  ppower = 0.5
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  ppower = 1e-4
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  flag_thfb = .true.
                  flag_tffb = .true.
                  flag_tmfb = .true.
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  keff_save = keff
                  branch_case(ic)%ndr2(i_p, i_step4) = 0.0
               else
                  flag_tmfb = .false.
                  ppower      = branch_case(ic)%ndr_ppower(i_p)
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  branch_case(ic)%ndr2(i_p , i_step4) = 100000 * LOG(keff/keff_save)
               endif
            enddo

         elseif(branch_case(ic)%if_ndr == 12 ) then !
            i_step5 =  i_step5 + 1
            opt_mode = 1
            ppower = 1.0
            flag_thfb = .false.
            flag_tffb = .true.
            flag_tmfb = .true.
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            ppower = 0.5
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            ppower = 1e-4
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            do i_p = 1, branch_case(ic)%n_ppower
               ppower      = branch_case(ic)%ndr_ppower(i_p)
               flag_thfb = .true.
               flag_tffb = .true.
               flag_tmfb = .true.
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               keff_save = keff
               ppower = ppower+ 0.01
               flag_tmfb = .false.
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               branch_case(ic)%ndr2(i_p , i_step5) = 100000 * LOG(keff/keff_save)
            enddo

         elseif(branch_case(ic)%if_ndr == 13 ) then ! Control rod worth
            opt_mode = 1
            !branch_case(ic)%i_substep = branch_case(ic)%i_substep + 1
            if (.not.allocated(branch_case(ic)%ndr2)) then
               call alloc(branch_case(ic)%ndr2 , branch_case(ic)%n_crw, crw_cases)
            endif
            i_s6 = i_s6 + 1
            write(*,*)crw_cases , i_s6
            if (branch_case(ic)%ppower  < 0.5) then
               flag_thfb = .true.
               flag_tffb = .true.
               flag_tmfb = .true.
               opt_xe   = 0
            else
               flag_thfb = .false.
               flag_tffb = .false.
               flag_tmfb = .false.
               opt_xe   = 1
            endif
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            keff_save = keff
            flag_thfb = .false.
            flag_tffb = .false.
            flag_tmfb = .false.
            do i_c = 1 ,  branch_case(ic)%n_crw
               cr_bot(branch_case(ic)%i_crw(i_c))= Length_CR_Tip( branch_case(ic)%i_crw(i_c))
               CR_top(branch_case(ic)%i_crw(i_c))= CR_bot(branch_case(ic)%i_crw(i_c))+length_CR(branch_case(ic)%i_crw(i_c))
               write(*,'(a)') 'cr_bot'
               write(*,'(100f7.2)') cr_bot(:)
!               call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
               branch_case(ic)%ndr2(i_c ,i_s6) = 100000 * LOG(keff_save/keff)
            enddo
            cr_bot(:) = length_cr_tip(:)+length_cr(:)
            cr_top(:) = cr_bot(:)+length_cr(:)
         elseif(branch_case(ic)%if_ndr == 14 ) then ! shutdown margin
            opt_mode = 1
            i_step3 = i_step3 + 1
            opt_xe   = 1
            flag_thfb = .true.
            flag_tffb = .true.
            flag_tmfb = .true.
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            keff_save = keff
            flag_thfb = .false.
            flag_tffb = .false.
            flag_tmfb = .false.
            do i_c = 1 ,  branch_case(ic)%n_crw
               cr_bot(branch_case(ic)%i_crw(i_c))= Length_CR_Tip( branch_case(ic)%i_crw(i_c))
               CR_top(branch_case(ic)%i_crw(i_c))= CR_bot(branch_case(ic)%i_crw(i_c))+length_CR(branch_case(ic)%i_crw(i_c))
            enddo
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            scram_worth = 100000 * LOG(keff_save/keff) * (1- branch_case(ic)%unc_sdm*0.01 )
            ppower = 1e-4
            opt_xe   = 0
            flag_thfb = .true.
            flag_tffb = .true.
            flag_tmfb = .true.
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            keff_save = keff
            ppower = 1.0
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            branch_case(ic)%ndr2(1,i_step3) = scram_worth
            branch_case(ic)%ndr2(2,i_step3) =  100000 * LOG(keff_save/keff)
            branch_case(ic)%ndr2(3,i_step3) = scram_worth - 100000 * LOG(keff_save/keff)

         elseif(branch_case(ic)%if_ndr == 15 ) then
            opt_mode = 1
            call alloc(branch_case(ic)%ndr2 , branch_case(ic)%n_ppm, branch_case(ic)%n_mt)
            do i_p = 1, branch_case(ic)%n_ppm
               ppm = branch_case(ic)%ndr_ppm(i_p)
               do i_c = 1, branch_case(ic)%n_mt
                  ppower      = branch_case(ic)%ndr_ppower(i_c)
                  tm_in       = branch_case(ic)%ndr_mt(i_c)
                  t_mod(:,:)  = tm_in + DegToK
                  if(branch_case(ic)%if_czp        ) then
                     t_fuel(:,:)       =  t_mod(:,:)
                  endif
                  if (i_c == 1 ) then
                     flag_thfb = .false.
                     flag_tffb = .false.
                     flag_tmfb = .false.
!                     call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                     keff_save = keff
                     branch_case(ic)%ndr2(i_p,i_c) = 0.0
                  else
                     flag_thfb = .false.
                     flag_tffb = .false.
                     flag_tmfb = .false.
                     write(*, '(A6, f12.3 , A6, f12.3)') 'tm_in = ' , tm_in ,'ppm = ' , ppm
!                     call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                     branch_case(ic)%ndr2(i_p,i_c) = 100000 * LOG(keff/keff_save)
                  endif
               enddo
            enddo

         elseif(branch_case(ic)%if_ndr == 16 ) then
            opt_mode = 1
            call alloc(branch_case(ic)%ndr1 , branch_case(ic)%n_crw)
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            keff_save = keff
            do i_p = 1, branch_case(ic)%n_crw
                CR_bot(branch_case(ic)%dcrm)= branch_case(ic)%cr_bot(i_p)
                CR_top(branch_case(ic)%dcrm)= CR_bot(branch_case(ic)%dcrm)+length_CR(branch_case(ic)%dcrm)
!                call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                branch_case(ic)%ndr1(i_p) = 100000* LOG(keff/keff_save)
            enddo

         elseif(branch_case(ic)%if_ndr == 17 ) then
            call alloc(branch_case(ic)%ndr2 , branch_case(ic)%n_ppm,branch_case(ic)%n_mt)
            call alloc(branch_case(ic)%ndr_ft , branch_case(ic)%n_mt)
            do i_p = 1, branch_case(ic)%n_ppm
               ppm = branch_case(ic)%ndr_ppm(i_p)
               do i_c = 1, branch_case(ic)%n_mt
                  flag_thfb = .false.
                  flag_tffb = .false.
                  flag_tmfb = .false.
                  ppower      = branch_case(ic)%ndr_ppower(i_c)
                  tm_in       = branch_case(ic)%ndr_mt(i_c)
                  t_mod(:,:)  = tm_in + DegToK
                  if(branch_case(ic)%if_czp        ) then
                     t_fuel(:,:)       =  t_mod(:,:)
                  endif
                  t_fuel(:,:) = t_fuel(:,:) + branch_case(ic)%d_TM1(1)
!                  call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  call Get_Avg( branch_case(ic)%TF_Avg(i_buu), T_Fuel, 0 )
                  branch_case(ic)%NDR_FT(i_c) = branch_case(ic)%TF_Avg(i_buu)
                  write(*, '(A10, f12.3, A10, f12.3 , A10, f12.3)') 'tm_in = ' , tm_in  ,'tf_avg = ' ,branch_case(ic)%TF_Avg(i_buu)  , ' ppm  = '  , ppm
                  keff_save = keff
                  t_fuel(:,:) = t_fuel(:,:) - branch_case(ic)%d_TM1(1) + branch_case(ic)%d_TM1(2)
!                  call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
!                  call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
                  call Get_Avg( branch_case(ic)%TF_Avg(i_buu), T_Fuel, 0 )
                  write(*, '(A10, f12.3, A10, f12.3 , A10, f12.3)') 'tm_in = ' , tm_in  ,'tf_avg = ' ,branch_case(ic)%TF_Avg(i_buu)  , ' ppm  = '  , ppm
                  branch_case(ic)%ndr2(i_p,i_c) = 100000 * LOG(keff_save/keff) / (branch_case(ic)%d_TM1(1) - branch_case(ic)%d_TM1(2))
                  t_fuel(:,:) = t_fuel(:,:) - branch_case(ic)%d_TM1(2)
!                  call nested_cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 

               enddo
            enddo

         endif

         flag_macroxs=branch_case(ic)%if_macro_xs

         if(branch_case(ic)%if_ndr == 0 ) then
!            call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
            if (no_axial == 0) then
!               call write_simulator(branch_case(ic)%case_name,i_buu)
            dummy_filler = 1 ! @$^ siarhei_plot 
            else
!               call write_simulator_axial(branch_case(ic)%case_name,i_buu)
            dummy_filler = 1 ! @$^ siarhei_plot 
            endif

            branch_case(ic)%keff(i_buu)=keff
            call Get_Avg( branch_case(ic)%TF_Avg(i_buu), T_Fuel, 0 )
            call Get_Avg( branch_case(ic)%TM_Avg(i_buu), T_Mod , 0 )
            call Get_Avg( branch_case(ic)%DM_Avg(i_buu), D_Mod , 0 )
            branch_case(ic)%ppm_out(i_buu)=ppm

            if(branch_case(ic)%if_d_TM) then
               tm_in = tm_in - branch_case(ic)%d_TM
!               call nested_Cal_dmod
            dummy_filler = 1 ! @$^ siarhei_plot 
            endif
         endif

         if (branch_case(ic)%if_macro_xs) then
!            call MacroXS_Table_dealloc
            dummy_filler = 1 ! @$^ siarhei_plot 
         endif

!         call reset_data
            dummy_filler = 1 ! @$^ siarhei_plot 
         if(if_ndr_logic ) then
!           call Search_EigVal_branch
            dummy_filler = 1 ! @$^ siarhei_plot 
         endif

      enddo

      flag_doing_branch=.false.

      flag_macroxs=.false.

!      call reset_data
            dummy_filler = 1 ! @$^ siarhei_plot 
      return

      end subroutine cal_branch

#ifdef siarhei_delete 
      subroutine nested_cal_dmod
      use Mod_THfunc, only: fdens
      use inc_rp
      use inc_geometry
      use inc_3d
!      USE MOD_H2O,  ONLY: H2O_T2D
      USE Inc_XS_File, ONLY: if_th

      implicit none
      integer(4) :: ixy,iz,ixy_fa

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [nested_cal_dmod] in Mod_Branch'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      do ixy_fa=1,nxy_fa
         ixy=i_Fa(ixy_fa)
         do iz=1,nz
            if(if_th==4)then
               d_mod(ixy,iz)=1 !H2O_T2D(t_mod(ixy,iz))
            else
!               d_mod(ixy,iz)=fdens(t_mod(ixy,iz)-DegToK)*1d-3
            dummy_filler = 1 ! @$^ siarhei_plot 
            endif
         enddo
      enddo
      return
      end subroutine nested_cal_dmod
#endif 


#ifdef siarhei_delete 
      subroutine backup_data
      use inc_3d
      use inc_branch
      use inc_option
      use inc_flag
      use inc_th
      use inc_cr
      use inc_nuclide
      use inc_geometry
      use mod_alloc
      use inc_3d, only: keff
      implicit none

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [backup_data] in Mod_Branch'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (allocated(CR_BOT)) then
         CALL Alloc( CR_bot_back , N_CR )
      endif

      if (.not. allocated(T_Mod_back)) then
         CALL Alloc( T_Mod_back  , Nxy , Nz )
         CALL Alloc( D_Mod_back  , Nxy , Nz )
         CALL Alloc( T_Fuel_back , Nxy , Nz )

         CALL Alloc( Flux_back, Nxy , Nz , N_group )
         CALL Alloc( Flux_old_back, Nxy , Nz , N_group )

         CALL Alloc( FisSrc_back , Nxy , Nz )
         CALL Alloc( Power_back , Nxy , Nz )
         CALL Alloc( Normal_Power_back , Nxy , Nz )

         CALL Alloc( N_U34_back  , Nxy , Nz )
         CALL Alloc( N_U35_back  , Nxy , Nz )
         CALL Alloc( N_U36_back  , Nxy , Nz )
         CALL Alloc( N_U37_back  , Nxy , Nz )
         CALL Alloc( N_U38_back  , Nxy , Nz )
         CALL Alloc( N_Np37_back , Nxy , Nz )
         CALL Alloc( N_Np38_back , Nxy , Nz )
         CALL Alloc( N_Np39_back , Nxy , Nz )
         CALL Alloc( N_Pu38_back , Nxy , Nz )
         CALL Alloc( N_Pu39_back , Nxy , Nz )
         CALL Alloc( N_Pu40_back , Nxy , Nz )
         CALL Alloc( N_Pu41_back , Nxy , Nz )
         CALL Alloc( N_Pu42_back , Nxy , Nz )
         CALL Alloc( N_Pu43_back , Nxy , Nz )
         CALL Alloc( N_Am41_back , Nxy , Nz )
         CALL Alloc( N_As42_back , Nxy , Nz )
         CALL Alloc( N_Am42_back , Nxy , Nz )
         CALL Alloc( N_Am43_back , Nxy , Nz )
         CALL Alloc( N_Am44_back , Nxy , Nz )
         CALL Alloc( N_Cm42_back , Nxy , Nz )
         CALL Alloc( N_Cm43_back , Nxy , Nz )
         CALL Alloc( N_Cm44_back , Nxy , Nz )
         CALL Alloc( N_I35_back  , Nxy , Nz )
         CALL Alloc( N_Xe35_back , Nxy , Nz )
         CALL Alloc( N_Nd47_back , Nxy , Nz )
         CALL Alloc( N_Nd48_back , Nxy , Nz )
         CALL Alloc( N_Nd49_back , Nxy , Nz )
         CALL Alloc( N_Pm47_back , Nxy , Nz )
         CALL Alloc( N_Ps48_back , Nxy , Nz )
         CALL Alloc( N_Pm48_back , Nxy , Nz )
         CALL Alloc( N_Pm49_back , Nxy , Nz )
         CALL Alloc( N_Sm47_back , Nxy , Nz )
         CALL Alloc( N_Sm48_back , Nxy , Nz )
         CALL Alloc( N_Sm49_back , Nxy , Nz )
         CALL Alloc( N_Gd52_back , Nxy , Nz )
         CALL Alloc( N_Gd54_back , Nxy , Nz )
         CALL Alloc( N_Gd55_back , Nxy , Nz )
         CALL Alloc( N_Gd56_back , Nxy , Nz )
         CALL Alloc( N_Gd57_back , Nxy , Nz )
         CALL Alloc( N_Gd58_back , Nxy , Nz )
         CALL Alloc( N_Gd60_back , Nxy , Nz )

         CALL Alloc( N_U34_old_back  , Nxy , Nz )
         CALL Alloc( N_U35_old_back  , Nxy , Nz )
         CALL Alloc( N_U36_old_back  , Nxy , Nz )
         CALL Alloc( N_U37_old_back  , Nxy , Nz )
         CALL Alloc( N_U38_old_back  , Nxy , Nz )
         CALL Alloc( N_Np37_old_back , Nxy , Nz )
         CALL Alloc( N_Np38_old_back , Nxy , Nz )
         CALL Alloc( N_Np39_old_back , Nxy , Nz )
         CALL Alloc( N_Pu38_old_back , Nxy , Nz )
         CALL Alloc( N_Pu39_old_back , Nxy , Nz )
         CALL Alloc( N_Pu40_old_back , Nxy , Nz )
         CALL Alloc( N_Pu41_old_back , Nxy , Nz )
         CALL Alloc( N_Pu42_old_back , Nxy , Nz )
         CALL Alloc( N_Pu43_old_back , Nxy , Nz )
         CALL Alloc( N_Am41_old_back , Nxy , Nz )
         CALL Alloc( N_As42_old_back , Nxy , Nz )
         CALL Alloc( N_Am42_old_back , Nxy , Nz )
         CALL Alloc( N_Am43_old_back , Nxy , Nz )
         CALL Alloc( N_Am44_old_back , Nxy , Nz )
         CALL Alloc( N_Cm42_old_back , Nxy , Nz )
         CALL Alloc( N_Cm43_old_back , Nxy , Nz )
         CALL Alloc( N_Cm44_old_back , Nxy , Nz )
         CALL Alloc( N_I35_old_back  , Nxy , Nz )
         CALL Alloc( N_Xe35_old_back , Nxy , Nz )
         CALL Alloc( N_Nd47_old_back , Nxy , Nz )
         CALL Alloc( N_Nd48_old_back , Nxy , Nz )
         CALL Alloc( N_Nd49_old_back , Nxy , Nz )
         CALL Alloc( N_Pm47_old_back , Nxy , Nz )
         CALL Alloc( N_Ps48_old_back , Nxy , Nz )
         CALL Alloc( N_Pm48_old_back , Nxy , Nz )
         CALL Alloc( N_Pm49_old_back , Nxy , Nz )
         CALL Alloc( N_Sm47_old_back , Nxy , Nz )
         CALL Alloc( N_Sm48_old_back , Nxy , Nz )
         CALL Alloc( N_Sm49_old_back , Nxy , Nz )
         CALL Alloc( N_Gd52_old_back , Nxy , Nz )
         CALL Alloc( N_Gd54_old_back , Nxy , Nz )
         CALL Alloc( N_Gd55_old_back , Nxy , Nz )
         CALL Alloc( N_Gd56_old_back , Nxy , Nz )
         CALL Alloc( N_Gd57_old_back , Nxy , Nz )
         CALL Alloc( N_Gd58_old_back , Nxy , Nz )
         CALL Alloc( N_Gd60_old_back , Nxy , Nz )

         CALL Alloc( N_B0_back , Nxy , Nz )
         CALL Alloc( BU_back , Nxy , Nz )
      endif

      flag_tmfb_back   = flag_tmfb
      flag_tffb_back   = flag_tffb

      opt_mode_back      = opt_mode
      opt_xe_back        = opt_xe
      opt_sm_back        = opt_sm
      opt_gd_back        = opt_gd
      OPT_findtavg_back  = OPT_findtavg
      flag_xsfb_back     = flag_xsfb
      flag_thfb_back     = flag_thfb
      ppm_back           = ppm
      ppower_back        = ppower
      core_massflow_back = core_massflow
      tm_in_back         = tm_in
      keff_back          = keff
      if (allocated(CR_BOT)) then
         cr_bot_back        = cr_bot
      endif

      flux_back          = flux
      flux_old_back      = flux_old

      FisSrc_back        = FisSrc
      Power_back         = Power
      Normal_Power_back  = Normal_Power

      t_mod_back         = t_mod
      d_mod_back         = d_mod
      t_fuel_back        = t_fuel
      tf_in_back         = tf_in

      N_U34_back         = N_U34
      N_U35_back         = N_U35
      N_U36_back         = N_U36
      N_U37_back         = N_U37
      N_U38_back         = N_U38
      N_Np37_back        = N_Np37
      N_Np38_back        = N_Np38
      N_Np39_back        = N_Np39
      N_Pu38_back        = N_Pu38
      N_Pu39_back        = N_Pu39
      N_Pu40_back        = N_Pu40
      N_Pu41_back        = N_Pu41
      N_Pu42_back        = N_Pu42
      N_Pu43_back        = N_Pu43
      N_Am41_back        = N_Am41
      N_As42_back        = N_As42
      N_Am42_back        = N_Am42
      N_Am43_back        = N_Am43
      N_Am44_back        = N_Am44
      N_Cm42_back        = N_Cm42
      N_Cm43_back        = N_Cm43
      N_Cm44_back        = N_Cm44
      N_U34_old_back     = N_U34_old
      N_U35_old_back     = N_U35_old
      N_U36_old_back     = N_U36_old
      N_U37_old_back     = N_U37_old
      N_U38_old_back     = N_U38_old
      N_Np37_old_back    = N_Np37_old
      N_Np38_old_back    = N_Np38_old
      N_Np39_old_back    = N_Np39_old
      N_Pu38_old_back    = N_Pu38_old
      N_Pu39_old_back    = N_Pu39_old
      N_Pu40_old_back    = N_Pu40_old
      N_Pu41_old_back    = N_Pu41_old
      N_Pu42_old_back    = N_Pu42_old
      N_Pu43_old_back    = N_Pu43_old
      N_Am41_old_back    = N_Am41_old
      N_As42_old_back    = N_As42_old
      N_Am42_old_back    = N_Am42_old
      N_Am43_old_back    = N_Am43_old
      N_Am44_old_back    = N_Am44_old
      N_Cm42_old_back    = N_Cm42_old
      N_Cm43_old_back    = N_Cm43_old
      N_Cm44_old_back    = N_Cm44_old

      N_I35_back         = N_I35
      N_Xe35_back        = N_Xe35
      N_Nd47_back        = N_Nd47
      N_Nd48_back        = N_Nd48
      N_Nd49_back        = N_Nd49
      N_Pm47_back        = N_Pm47
      N_Ps48_back        = N_Ps48
      N_Pm48_back        = N_Pm48
      N_Pm49_back        = N_Pm49
      N_Sm47_back        = N_Sm47
      N_Sm48_back        = N_Sm48
      N_Sm49_back        = N_Sm49
      N_I35_old_back     = N_I35_old
      N_Xe35_old_back    = N_Xe35_old
      N_Nd47_old_back    = N_Nd47_old
      N_Nd48_old_back    = N_Nd48_old
      N_Nd49_old_back    = N_Nd49_old
      N_Pm47_old_back    = N_Pm47_old
      N_Ps48_old_back    = N_Ps48_old
      N_Pm48_old_back    = N_Pm48_old
      N_Pm49_old_back    = N_Pm49_old
      N_Sm47_old_back    = N_Sm47_old
      N_Sm48_old_back    = N_Sm48_old
      N_Sm49_old_back    = N_Sm49_old

      N_Gd52_back        = N_Gd52
      N_Gd54_back        = N_Gd54
      N_Gd55_back        = N_Gd55
      N_Gd56_back        = N_Gd56
      N_Gd57_back        = N_Gd57
      N_Gd58_back        = N_Gd58
      N_Gd60_back        = N_Gd60
      N_Gd52_old_back    = N_Gd52_old
      N_Gd54_old_back    = N_Gd54_old
      N_Gd55_old_back    = N_Gd55_old
      N_Gd56_old_back    = N_Gd56_old
      N_Gd57_old_back    = N_Gd57_old
      N_Gd58_old_back    = N_Gd58_old
      N_Gd60_old_back    = N_Gd60_old

      N_B0_back          = N_B0
      BU_back            = BU

      return
      end subroutine backup_data
#endif 

#ifdef siarhei_delete 
      subroutine reset_data
      use inc_3d
      use inc_branch
      use inc_option
      use inc_flag
      use inc_th
      use inc_cr
      use inc_nuclide
      use inc_geometry
      use Mod_THdriver
      use inc_3d, only: keff
      implicit none

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [reset_data] in Mod_Branch'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      flag_tmfb   = flag_tmfb_back
      flag_tffb   = flag_tffb_back

      opt_mode      = opt_mode_back
      opt_xe        = opt_xe_back
      opt_sm        = opt_sm_back
      opt_gd        = opt_gd_back
      OPT_findtavg  = OPT_findtavg_back
      flag_xsfb     = flag_xsfb_back
      flag_thfb     = flag_thfb_back
      ppm           = ppm_back
      ppower        = ppower_back
      core_massflow = core_massflow_back
      tm_in         = tm_in_back
      keff          = keff_back
      if (allocated(CR_BOT)) then
         cr_bot        = cr_bot_back
         CR_top(:)=CR_bot(:)+length_CR(:)
      endif

      flux          = flux_back
      flux_old      = flux_old_back

      FisSrc        = FisSrc_back
      Power         = Power_back
      Normal_Power  = Normal_Power_back

      t_mod         = t_mod_back
      d_mod         = d_mod_back
      t_fuel        = t_fuel_back

      N_U34         = N_U34_back
      N_U35         = N_U35_back
      N_U36         = N_U36_back
      N_U37         = N_U37_back
      N_U38         = N_U38_back
      N_Np37        = N_Np37_back
      N_Np38        = N_Np38_back
      N_Np39        = N_Np39_back
      N_Pu38        = N_Pu38_back
      N_Pu39        = N_Pu39_back
      N_Pu40        = N_Pu40_back
      N_Pu41        = N_Pu41_back
      N_Pu42        = N_Pu42_back
      N_Pu43        = N_Pu43_back
      N_Am41        = N_Am41_back
      N_As42        = N_As42_back
      N_Am42        = N_Am42_back
      N_Am43        = N_Am43_back
      N_Am44        = N_Am44_back
      N_Cm42        = N_Cm42_back
      N_Cm43        = N_Cm43_back
      N_Cm44        = N_Cm44_back
      N_U34_old     = N_U34_old_back
      N_U35_old     = N_U35_old_back
      N_U36_old     = N_U36_old_back
      N_U37_old     = N_U37_old_back
      N_U38_old     = N_U38_old_back
      N_Np37_old    = N_Np37_old_back
      N_Np38_old    = N_Np38_old_back
      N_Np39_old    = N_Np39_old_back
      N_Pu38_old    = N_Pu38_old_back
      N_Pu39_old    = N_Pu39_old_back
      N_Pu40_old    = N_Pu40_old_back
      N_Pu41_old    = N_Pu41_old_back
      N_Pu42_old    = N_Pu42_old_back
      N_Pu43_old    = N_Pu43_old_back
      N_Am41_old    = N_Am41_old_back
      N_As42_old    = N_As42_old_back
      N_Am42_old    = N_Am42_old_back
      N_Am43_old    = N_Am43_old_back
      N_Am44_old    = N_Am44_old_back
      N_Cm42_old    = N_Cm42_old_back
      N_Cm43_old    = N_Cm43_old_back
      N_Cm44_old    = N_Cm44_old_back

      N_I35         = N_I35_back
      N_Xe35        = N_Xe35_back
      N_Nd47        = N_Nd47_back
      N_Nd48        = N_Nd48_back
      N_Nd49        = N_Nd49_back
      N_Pm47        = N_Pm47_back
      N_Ps48        = N_Ps48_back
      N_Pm48        = N_Pm48_back
      N_Pm49        = N_Pm49_back
      N_Sm47        = N_Sm47_back
      N_Sm48        = N_Sm48_back
      N_Sm49        = N_Sm49_back
      N_I35_old     = N_I35_old_back
      N_Xe35_old    = N_Xe35_old_back
      N_Nd47_old    = N_Nd47_old_back
      N_Nd48_old    = N_Nd48_old_back
      N_Nd49_old    = N_Nd49_old_back
      N_Pm47_old    = N_Pm47_old_back
      N_Ps48_old    = N_Ps48_old_back
      N_Pm48_old    = N_Pm48_old_back
      N_Pm49_old    = N_Pm49_old_back
      N_Sm47_old    = N_Sm47_old_back
      N_Sm48_old    = N_Sm48_old_back
      N_Sm49_old    = N_Sm49_old_back

      N_Gd52        = N_Gd52_back
      N_Gd54        = N_Gd54_back
      N_Gd55        = N_Gd55_back
      N_Gd56        = N_Gd56_back
      N_Gd57        = N_Gd57_back
      N_Gd58        = N_Gd58_back
      N_Gd60        = N_Gd60_back
      N_Gd52_old    = N_Gd52_old_back
      N_Gd54_old    = N_Gd54_old_back
      N_Gd55_old    = N_Gd55_old_back
      N_Gd56_old    = N_Gd56_old_back
      N_Gd57_old    = N_Gd57_old_back
      N_Gd58_old    = N_Gd58_old_back
      N_Gd60_old    = N_Gd60_old_back

      N_B0          = N_B0_back
      BU            = BU_back

      fa_massflow   = core_massflow * chanvf_save / Core_N_FA
!      if(flag_thfb) call restore_massflow
            dummy_filler = 1 ! @$^ siarhei_plot 

      c_coef=''
      c_case=0

      return
      end subroutine reset_data
#endif 


#ifdef siarhei_delete 
      SUBROUTINE Search_EigVal_BRANCH
!      USE Mod_SSnodal, ONLY: NodalSolver_St
!     ! USE Mod_Pinpow, only:pinpow_f
      USE Inc_3D, ONLY: Power, Normal_Power, avg_power
      USE Inc_File
      USE Inc_INP
      USE Inc_3D, ONLY: keff
      USE Inc_RST
      USE Inc_TH
      USE Mod_GetSome
      USE Read_File, ONLY: Read_RI
      USE Mod_Save, ONLY: Save_History, XYZ_Indexing
      USE Write_O
      USE Mod_XSFB
      use inc_branch, only: flag_stack
      use mod_getsome, only: levelbalancing
      IMPLICIT NONE
      INTEGER :: Ixy, Ixy_FA, Iz

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [Search_EigVal_BRANCH] in Mod_Branch'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      Core_Power     = Core_Power_100 * PPower
      Core_Power_0   = Core_Power
      PPower_0       = PPower
      FA_Power     = Core_Power / Core_N_FA
      !!!!!FA_Power_0   = FA_Power
      Avg_CorePower = Core_Power / Tot_FuelVol
      !!!!!Avg_CorePower0 = Core_Power / Tot_FuelVol

      if (.not. flag_stack) then
         WRITE(*, '(A)') "==      EIGENVALUE CALCULATION START      =="
      endif

      ! XS Feedback
      CALL XSFB

      ! Nodal Solver
!      CALL NodalSolver_St
            dummy_filler = 1 ! @$^ siarhei_plot 

      ! Update Power
      CALL Get_POW

      ! Flux Level Balancing for Design Power
      CALL Get_Avg(Avg_Power, Power, 0)
      LevFactor = Avg_CorePower/Avg_Power
      CALL LevelBalancing

      ! Power Normalization
      CALL Get_Avg(Avg_Power, Power, 0)

      DO Ixy_FA = 1, Nxy_FA
         Ixy = I_FA(Ixy_FA)

         DO Iz = IzFuelBot, IzFuelTop
            Normal_Power(Ixy, Iz) = Power(Ixy, Iz) / Avg_Power
         END DO
      END DO

      ! Get XYZ Indexing (1N/1FA)
!      CALL XYZ_Indexing
            dummy_filler = 1 ! @$^ siarhei_plot 

      ! Pin Power Reconstruction
     ! call pinpow_f

      if (.not. flag_stack) then
         WRITE(*, "(A, F8.6, /)") "*** keff = ", keff
      endif
      RETURN
      END SUBROUTINE Search_EigVal_BRANCH
#endif 

      END MODULE Mod_Branch

