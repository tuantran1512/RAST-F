   module Mod_PinPow_Hex
#ifdef siarhei_ppr
      use Inc_PinPow_Hex
      use Inc_TPEN
#endif

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub
#endif

#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none
      contains

#ifdef siarhei_ppr
      subroutine open_file_hex(namee,file_id)
         implicit none
         character(200)    :: namee
         character(200)    :: filename
         integer           :: file_id

#ifdef siarhei_plot
         ! @$^ siarhei_fr

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [open_file_hex] in Mod_PinPow_Hex'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         write(678,'(A)') '^#^ Entered [open_file_hex] in Mod_PinPow_Hex'
         ! =-=-=-=-=-=-=-=
#endif
         filename = trim(namee)//'_pin_power_output.dat'
         open(file_id,file=filename,status='replace',&
         access='sequential',form='formatted',action='write')

      end subroutine open_file_hex

      subroutine advance_yes_file_plot(file_id)
         implicit none
         integer           :: file_id

#ifdef siarhei_plot
         ! @$^ siarhei_fr

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [advance_yes_file_plot] in Mod_PinPow_Hex'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         current_sub = '^#^ Entered [advance_yes_file_plot] in Mod_PinPow_Hex'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         write(file_id,'(A)',advance='yes')

      end subroutine advance_yes_file_plot

      subroutine close_file_hex(file_id)
         implicit none
         integer           :: file_id


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [close_file_hex] in Mod_PinPow_Hex'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         close(file_id)

      end subroutine close_file_hex


      subroutine allocate_pinpow_hex
         use Inc_PinPow_Hex
         use Inc_TPEN
         use Inc_Option,   only: n_group
         use Inc_Geometry, only: Nz
        ! use Inc_TPEN, only: nassy
         implicit none

         integer                    :: actual_xy(2) ! size of FA

#ifdef siarhei_plot
         ! @$^ siarhei_fr

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [allocate_pinpow_hex] in Mod_PinPow_Hex'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         current_sub = '^#^ Entered [allocate_pinpow_hex] in Mod_PinPow_Hex'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif


         if (.not.allocated(inp_HFF_hex)) then
            select case(N_FA_parts)
               case(6) ! triangle
                  actual_xy(1) = 2*(radNum-1) + 1
                  actual_xy(2) = radNum
                  allocate(inp_HFF_hex(actual_xy(1),actual_xy(2),nassy))
               case(2) ! half-hex
                  actual_xy(1) = 4*(radNum-1) + 1
                  actual_xy(2) = radNum
                  allocate(inp_HFF_hex(actual_xy(1),actual_xy(2),nassy))
               case(1) ! full hex
                  actual_xy(1) = 4*(radNum-1) + 1
                  actual_xy(2) = 2*(radNum-1) + 1
                  allocate(inp_HFF_hex(actual_xy(1),actual_xy(2),nassy))
               case default ! triangle
                  actual_xy(1) = 2*(radNum-1) + 1
                  actual_xy(2) = radNum
                  allocate(inp_HFF_hex(actual_xy(1),actual_xy(2),nassy))
            end select
            inp_HFF_hex = -7.777 ! irrelevant surrounding of a FA
         end if

         actual_xy(1) = 4*(radNum-1) + 1
         actual_xy(2) = 2*(radNum-1) + 1
         if (.not.allocated(HFF_hex)) then
            allocate(HFF_hex(actual_xy(1),actual_xy(2),nassy))
            allocate(zero_HFF_hex(actual_xy(2),actual_xy(1)))
            allocate(x_coord_fa_hex(actual_xy(2),actual_xy(1)))
            allocate(y_coord_fa_hex(actual_xy(2),actual_xy(1)))
            allocate(p_coord_fa_hex(actual_xy(2),actual_xy(1)))
            allocate(u_coord_fa_hex(actual_xy(2),actual_xy(1)))
            allocate(rad_coord_fa_hex(actual_xy(2),actual_xy(1)))
            allocate(sin_coord_fa_hex(actual_xy(2),actual_xy(1)))
            allocate(ang_coord_fa_hex(actual_xy(2),actual_xy(1)))
            allocate(x_y_u_p_coord_fa_hex(actual_xy(2),actual_xy(1),4))
            allocate(location_in_hex(actual_xy(2),actual_xy(1)))

            HFF_hex = -7.777 ! irrelevant surrounding of a FA
            zero_HFF_hex = -7.777
            x_coord_fa_hex = 0d0
            y_coord_fa_hex = 0d0
            p_coord_fa_hex = 0d0
            u_coord_fa_hex = 0d0
            rad_coord_fa_hex = 0d0
            sin_coord_fa_hex = 0d0
            ang_coord_fa_hex = 0d0
            x_y_u_p_coord_fa_hex = 0d0
         end if

         ! allocate fluxes as (n_triangles,n_group, n_FAs)
         if (.not.allocated(Fl_momx       )) allocate(Fl_momx(6, n_group, nassy, Nz))
         if (.not.allocated(Fl_momy       )) allocate(Fl_momy(6, n_group, nassy, Nz))
         if (.not.allocated(Fl_avg        )) allocate(Fl_avg(6, n_group, nassy, Nz))

         if (.not.allocated(corner_flux_h )) allocate(corner_flux_h(nxpnt, n_group, Nz))
         if (.not.allocated(central_flux_h)) allocate(central_flux_h(n_group, nassy, Nz))
         if (.not.allocated(inner_b_flx_h )) allocate(inner_b_flx_h(6, n_group, nassy, Nz))
         if (.not.allocated(outer_b_flx_h )) allocate(outer_b_flx_h(6, n_group, nassy, Nz))

         if (.not.allocated(hom_pow_hex )) &
            allocate(hom_pow_hex(actual_xy(2),actual_xy(1),nassy,Nz))
         if (.not.allocated(xy_pow_hex )) &
            allocate(xy_pow_hex(actual_xy(2),actual_xy(1),nassy,Nz))
         hom_pow_hex = 0d0
         xy_pow_hex = 0d0

         Fl_momx = 0d0
         Fl_momy = 0d0
         Fl_avg  = 0d0

      end subroutine allocate_pinpow_hex


      subroutine read_hff_hex(filename)
         use Inc_PinPow_Hex
         use Inc_TPEN, only: nassy
         implicit none
         character(300)                 :: oneline, long_line
         character(5)                   :: text_gap = '     ' ! 5 spaces
         integer                        :: i,j,k,rad
         character(15)                  :: tmp_inp
         logical                        :: red_flag, green_flag
         integer                        :: curr_step, x_count, y_count,real_len
         character(20)                  :: read_style
         integer                        :: actual_xy(2),actual_center(2) !size FA
         real(8)                        :: tmp_real
         character(len=:),allocatable   :: actual_oneline
         character(200)                 :: filename

#ifdef siarhei_plot
         ! @$^ siarhei_fr

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [read_hff_hex] in Mod_PinPow_Hex'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         current_sub = '^#^ Entered [read_hff_hex] in Mod_PinPow_Hex'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         open(11,file=filename,status='old',&
         access='sequential',form='formatted',action='read')
         red_flag = .false.  ! stop reading HFF
         green_flag = .false. ! start reading HFF
         curr_step = -4444
         x_count = -4444
         y_count = -4444

         do while(.not.red_flag)
            oneline = ''
            read(11,'(a)') oneline
            ! check for empty lines
            if (len(trim(oneline)).eq.0) then
               if (green_flag) then
                  green_flag = .false.
                  write(*,*) len(trim(oneline)), oneline
                  goto 112
               else
                  goto 112
               end if
            end if
            ! look for subcards
            do i = 1,5
               if (oneline(i:i+1).eq.'!#') then
                  red_flag = .true.
                  green_flag = .false.
                  goto 112
               else if (oneline(i:i+2).eq.'HFF') then
                  green_flag = .true.
                  curr_step = 1
                  x_count = 1
                  y_count = 1
                  goto 112
               else if (oneline(i:i+5).eq.'HEXDIM') then
                  read(oneline,*) tmp_inp, N_FA_parts, radNum, radPtch
                 ! write(*,"(a6,x,i,x,i,x,f5.3)") tmp_inp, N_FA_parts, radNum, radPtch

                  call allocate_pinpow_hex

                  write(*,*) '[subroutine read_hff_hex] Arrays allocated'
                  goto 112
               end if
            end do

            ! read HFF based on the input options
            if (green_flag) then
               actual_oneline = trim(oneline(6:))
               real_len = int(len(actual_oneline)/5)
               write(*,*) len(actual_oneline)
               write(tmp_inp,*) real_len
               tmp_inp = '('//trim(adjustl(tmp_inp))//'f5.3)'
               write(*,*) 'Chosen format of read:',tmp_inp
               read(actual_oneline,trim(tmp_inp),end=911,err=911) &
                  inp_HFF_hex(y_count,1:real_len,1)
               deallocate(actual_oneline)
               do i = 1,real_len
                  if(inp_HFF_hex(y_count,i,1).eq.0.0) &
                     inp_HFF_hex(y_count,i,1) = -7.777 !
                     ! irrelevant surrounding of a FA
               end do
               y_count = y_count + 1
            end if

      ! 914         read(oneline,'(A)',advance='yes')
911         continue
112         cycle
         end do

         close(11) ! close HFF file

         ! choose printing layout (geometry)
         select case(N_FA_parts)
            case(6) ! triangle
               actual_xy(1) = 2*(radNum-1) + 1
               actual_xy(2) = radNum
            case(2) ! half-hex
               actual_xy(1) = 4*(radNum-1) + 1
               actual_xy(2) = radNum
            case(1) ! full hex
               actual_xy(1) = 4*(radNum-1) + 1
               actual_xy(2) = 2*(radNum-1) + 1
            case default ! triangle
               actual_xy(1) = 2*(radNum-1) + 1
               actual_xy(2) = radNum
         end select
         ! print input array to the console screen
         do i = 1, actual_xy(1)
            long_line = ''
            do j = 1, actual_xy(2)
               if (inp_HFF_hex(i,j,1).gt.9.0) then
                  write(*,'(f5.3)',advance='no') 0.000
               else if (inp_HFF_hex(i,j,1).lt.0.0) then
                  write(*,'(a5)',advance='no') text_gap
               else
                  write(*,'(f5.3)',advance='no') inp_HFF_hex(i,j,1)
               end if
            end do
            write(*,'(A)',advance='yes')
         end do

         ! prepare to convert input to full FA
         actual_xy(1) = 4*(radNum-1) + 1
         actual_xy(2) = 2*(radNum-1) + 1
         actual_center(1) = 2*(radNum-1) + 1
         actual_center(2) = radNum


         if (N_FA_parts.eq.6) then ! 1/6 FA
            do i = actual_center(2), 1, -1
               do j = actual_center(1), 1, -1
                  ! copy the 1/6 slice to the full FA layout
                  HFF_hex(j+(radNum-1),i,1) = inp_HFF_hex(j,i,1)
               end do
            end do
            do i = actual_center(2), 1, -1
               do rad = 1,i-1
                  ! Add NW slice (1/6)
                  HFF_hex(actual_center(1)-rad*2-(actual_center(2)-i),i,1) = &
                     HFF_hex(actual_center(1)+(actual_center(2)-i)-rad,i-rad,1)
                  write(*,'(A5,I3,A1,I3,A4,I3,A1,I3)') 'From ',actual_center(1)+(actual_center(2)-i)-rad,&
                     ',',i-rad,' to ',actual_center(1)-rad*2-(actual_center(2)-i), &
                     ',',i
                 ! HFF_hex(actual_center(1)-rad*2-(actual_center(2)-i),i,1) = &
                 !    HFF_hex(actual_center(1)+(actual_center(2)-i)-rad,i-rad,1)
                  ! Add SW slice (1/6)
                  HFF_hex(actual_center(1)+rad*2+(actual_center(2)-i),i,1) = &
                     HFF_hex(actual_center(1)-(actual_center(2)-i)+rad,i-rad,1)
               end do
            end do
            ! reflect and turn around as if the input was 1/2 FA
            do i = actual_center(2)-1,1,-1
               do j = actual_xy(1), 1, -1
                 ! HFF_hex(actual_xy(1)-j+1,actual_xy(2)-i+1,1) = HFF_hex(j,i,1)
                  HFF_hex(j,actual_xy(2)-i+1,1) = HFF_hex(actual_xy(1)-j+1,i,1)
               end do
            end do

         else if (N_FA_parts.eq.2) then ! half-FA
            do i = actual_center(2)-1,1,-1
               do j = actual_xy(1), 1, -1
                  ! copy from input
                  HFF_hex(j,i,1) = inp_HFF_hex(j,i,1)
                  ! reflect the copied part to the opposite part
                  HFF_hex(actual_xy(1)-j+1,actual_xy(2)-i+1,1) = HFF_hex(j,i,1)
               end do
            end do
            ! copy the central column
            HFF_hex(:,actual_center(2),1) = inp_HFF_hex(:,actual_center(2),1)

         else if (N_FA_parts.eq.1) then ! full FA
            ! copy from input
            HFF_hex(:,:,:) = inp_HFF_hex(:,:,:)

         else
            write(*,*) "[subroutine read_hff_hex] Error! Unsupported HFF geometry"
            pause
         end if
         do concurrent (i=1:actual_xy(2),j=1:actual_xy(1))
        ! do i = 1,actual_xy(2)
        !    do j = 1,actual_xy(1)
            if (HFF_hex(j,i,1).gt.9.0) HFF_hex(j,i,1) = 0.000
              !$ create map in here < - - subroutine
         !   end do
         end do
         ! temporary until STN HFF implemented
         do concurrent (i=2:nassy)
        ! do i = 2, nassy
            HFF_hex(:,:,i) = HFF_hex(:,:,1)
         end do
         ! print full FA array to the console screen
         ! and fill out the empty FA with zero data
         do i = 1, actual_xy(1)
            long_line = ''
            do j = 1, actual_xy(2)
               if (HFF_hex(i,j,1).lt.0.0) then
                  write(*,'(a5)',advance='no') text_gap
                !  write(*,'(f4.1)',advance='no') HFF_hex(i,j,1)
                 ! write(*,'(a2)',advance='no') ' '
               else
                  write(*,'(f5.3)',advance='no') HFF_hex(i,j,1)
                  zero_HFF_hex(j,i) = 0.000
                 ! write(*,'(a2)',advance='no') ' '
               end if
            end do
            write(*,'(A)',advance='yes')
         end do

         do i = 1, actual_xy(1)
            long_line = ''
            do j = 1, actual_xy(2)
               if (zero_HFF_hex(j,i).lt.0.0) then
                  write(*,'(a5)',advance='no') text_gap
               else
                  write(*,'(f5.3)',advance='no') zero_HFF_hex(j,i)
               end if
            end do
            write(*,'(A)',advance='yes')
         end do

         if (allocated(inp_HFF_hex)) deallocate(inp_HFF_hex)

      end subroutine read_hff_hex


      subroutine get_plotting_boundaries_core(xmin,xmax,ymin,ymax)
         use Inc_TPEN, only:nx_hex,ny_hex
         use Inc_PinPow_Hex, only: saved_nodel_hex
         implicit none
         integer,intent(inout)     :: xmin,xmax,ymin,ymax
         logical                   :: empty_col,empty_row
         integer                   :: i,j

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [get_plotting_boundaries_core] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_plotting_boundaries_core] in Mod_PinPow_Hex'
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

         xmin = -4444
         empty_col = .true.
! modified by Tuan
         do j = -3, nx_hex+4
            do i = -1, ny_hex+2
!         do j = 1, nx_hex+8
!            do i = 1, ny_hex+4
               if (saved_nodel_hex(j,i).ne.0) then
                  empty_col = .false.
                !  write(*,*)j,i,saved_nodel_hex(j,i)
                  if (xmin.lt.0) xmin = j
                  exit
               end if
            end do
            if (.not.empty_col) exit
         end do
         xmax = -4444
         empty_col = .true.
! modified by Tuan
         do j = nx_hex+4 ,xmin, -1
            do i = -1, ny_hex+2
!         do j = nx_hex+8 ,xmin, -1
!            do i = 1, ny_hex+4
               if (saved_nodel_hex(j,i).ne.0) then
                  empty_col = .false.
                !  write(*,*)j,i,saved_nodel_hex(j,i)
                  if (xmax.lt.0) xmax = j
                  exit
               end if
            end do
            if (.not.empty_col) exit
         end do
         ! row-wise
         ymin = -4444
         empty_row = .true.
! modified by Tuan
         do i = -1, ny_hex+2
            do j = -3, nx_hex+4
!         do i = 1, ny_hex+4
!            do j = 1, nx_hex+8
               if (saved_nodel_hex(j,i).ne.0) then
                  empty_row = .false.
                !  write(*,*)j,i,saved_nodel_hex(j,i)
                  if (ymin.lt.0) ymin = i
                  exit
               end if
            end do
            if (.not.empty_row) exit
         end do
         ymax = -4444
         empty_row = .true.
         do i =  ny_hex+4,ymin,-1
            do j = 1, nx_hex+8
               if (saved_nodel_hex(j,i).ne.0) then
                  empty_row = .false.
                !  write(*,*)j,i,saved_nodel_hex(j,i)
                  if (ymax.lt.0) ymax = i
                  exit
               end if
            end do
            if (.not.empty_row) exit
         end do
         write(*,*) "[Mod_PinPow_Hex -- get_plotting_boundaries_core] ",&
            "Non-zero borders of saved_nodel_hex: (",xmin,&
            ',',ymin,')  (',xmax,',',ymax,')'

      end subroutine get_plotting_boundaries_core

      subroutine calculate_pin_power_hex
         use Inc_PinPow_Hex
         use Inc_TPEN
         use Inc_Option,   only: n_group
         use Inc_Geometry, only: Nz
         use Inc_File,     only: Name_INP, Len_INP

         implicit none
         integer                                   :: i,j,k,gr,z_ax
         integer                                   :: xmin,xmax,ymin,ymax
         integer                                   :: actual_xy(2)

         character(5)                              :: z_axis_string
         character(200)                            :: filename,print_name_tmp
         character(50)                             :: prefix
         character(10)                             :: point_flx
         integer                                   :: ci,cj,ck
         character(10000),allocatable,dimension(:) :: long_line
         character(10000)                          :: empty_line
         character(5)                   :: text_gap = '     ' ! 5 spaces
         real(8)                                   :: dummy_large = 9999999
         integer                                   :: dummy_int
         logical                                   :: not_balanced
         integer,allocatable,dimension(:)          :: lft_cnt,rgt_cnt

         real(8)                                   :: avg_pin_mlinpow,avg_pin_integpow
         integer                                   :: max_Fq_hex_loc(4),max_FdH_hex_loc(3)
         integer                                   :: max_Fxy_hex_loc(4),max_Fz_hex_loc(1)
         real(8)                                   :: max_Fq_hex,max_FdH_hex
         real(8)                                   :: max_Fxy_hex,max_Fz_hex

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [calculate_pin_power_hex] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [calculate_pin_power_hex] in Mod_PinPow_Hex'
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

         ! new array for Ф(x,y) values
         actual_xy(1) = 4*(radNum-1) + 1
         actual_xy(2) = 2*(radNum-1) + 1
         allocate(xy_Flux_hex(actual_xy(2),actual_xy(1),nassy,n_group,Nz))
         allocate(lin_pow_hex(actual_xy(2),actual_xy(1),nassy,Nz))
         allocate(xy_Fq_hex(actual_xy(2),actual_xy(1),nassy,Nz))
         allocate(xy_FdH_hex(actual_xy(2),actual_xy(1),nassy))
         allocate(xy_Fxy_hex(actual_xy(2),actual_xy(1),nassy,Nz))
         allocate(xy_Fz_hex(Nz))
         allocate(integral_pin_pow_hex(actual_xy(2),actual_xy(1),nassy))
         xy_Flux_hex = -4444
         xy_Fq_hex = 0d0
         xy_FdH_hex = 0d0
         xy_Fxy_hex = 0d0
         xy_Fz_hex = 0d0

         ! check for negative fluxes (temporary)
         do z_ax = 1,Nz
            do i = 1,nassy
               do j = 1,n_group
                  do k = 1,6
                     if (inner_b_flx_h(k,j,i,z_ax).lt.0.0) then
                        write(*,*) 'Inner',inner_b_flx_h(k,j,i,z_ax),k,j,i,z_ax
                     end if
                     if (outer_b_flx_h(k,j,i,z_ax).lt.0.0) then
                        write(*,*) 'Outer',outer_b_flx_h(k,j,i,z_ax),k,j,i,z_ax
                     end if
                  end do
               end do
            end do
            do i = 1,n_group
               do j = 1,nxpnt
                  if (corner_flux_h(j,i,z_ax).lt.0.0) then
                     write(*,*) 'Corner',corner_flux_h(j,i,z_ax),j,i,z_ax
                     ! Siarhei <----- CHECK this situation!!!
                     ! corner_flux_h(j,i,z_ax) = 0d0
                  end if
               end do
               do j = 1,nassy
                  if (central_flux_h(i,j,z_ax).lt.0.0) then
                     write(*,*) 'Central',central_flux_h(i,j,z_ax),i,j,z_ax
                  end if
               end do
            end do
         end do


         ! fill out the arrays for each specific coordinate (x,y,p,u)
         call set_coordinates_fa_hex(actual_xy(2),actual_xy(1))


         do i = 1,actual_xy(1) ! y
            do j = 1,actual_xy(2) ! x
               ! find required coordinates (x,y,u,p) inside one triangle
               ! (x,y,angle,radius,height,center_x)
               x_y_u_p_coord_fa_hex(j,i,1:4) = get_x_y_u_p_hex( &
                                             x_coord_fa_hex(j,i),  &
                                             y_coord_fa_hex(j,i),  &
                                             ang_coord_fa_hex(j,i),&
                                             rad_coord_fa_hex(j,i),&
                                             triangle_height_hex,  &
                                             center_x_triangle_hex)
               location_in_hex(j,i) = define_location_hex(ang_coord_fa_hex(j,i))
            end do
         end do

         ! needed for setting up actual core boundaries
         ! finding the actual LP borders (as a rectangle for better output printing)
         ! column-wise
         call get_plotting_boundaries_core(xmin,xmax,ymin,ymax)

         ! obtain Ф(x,y) for each 'non-zero' node
         ! function determine_2d_flux_hex(x_y_u_p,          &
         !                     Fl_momx,Fl_momy,             &
         !                     Fl_avg,Fl_c_x,Fl_c_p,        &
         !                     Fl_c_u,Fl_s_x,Fl_s_p,Fl_s_u, &
         !                     radpitch,radnumber)

         do z_ax = 1,Nz
            write(*,'(A12,I2)') 'Axial layer ',z_ax
            do gr = 1,n_group
               do cj = ymin, ymax
                  do ci = xmin, xmax
                     if (saved_nodel_hex(ci,cj).gt.0) then
                        do i = 1, actual_xy(1)
                           do j = 1, actual_xy(2)
                              select case(location_in_hex(j,i))
                              case(1) ! triangle 1
                                 xy_Flux_hex(j,i,saved_nodel_hex(ci,cj),gr,z_ax) =      &
                                       determine_2d_flux_hex(                      &
                                          x_y_u_p_coord_fa_hex(j,i,:),             &
                                          Fl_momx(1,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_momy(1,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_avg(1,gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          central_flux_h(gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          corner_flux_h(neigpt(2,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          corner_flux_h(neigpt(1,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          outer_b_flx_h(1,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(6,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(1,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          radPtch,radNum)
                              case(2) ! triangle 2
                                 xy_Flux_hex(j,i,saved_nodel_hex(ci,cj),gr,z_ax) =      &
                                       determine_2d_flux_hex(                      &
                                          x_y_u_p_coord_fa_hex(j,i,:),             &
                                          Fl_momx(2,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_momy(2,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_avg(2,gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          central_flux_h(gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          corner_flux_h(neigpt(3,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          corner_flux_h(neigpt(2,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          outer_b_flx_h(2,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(1,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(2,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          radPtch,radNum)
                              case(3) ! triangle 3
                                 xy_Flux_hex(j,i,saved_nodel_hex(ci,cj),gr,z_ax) =      &
                                       determine_2d_flux_hex(                      &
                                          x_y_u_p_coord_fa_hex(j,i,:),             &
                                          Fl_momx(3,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_momy(3,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_avg(3,gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          central_flux_h(gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          corner_flux_h(neigpt(4,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          corner_flux_h(neigpt(3,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          outer_b_flx_h(3,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(2,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(3,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          radPtch,radNum)
                              case(4) ! triangle 4
                                 xy_Flux_hex(j,i,saved_nodel_hex(ci,cj),gr,z_ax) =      &
                                       determine_2d_flux_hex(                      &
                                          x_y_u_p_coord_fa_hex(j,i,:),             &
                                          Fl_momx(4,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_momy(4,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_avg(4,gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          central_flux_h(gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          corner_flux_h(neigpt(5,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          corner_flux_h(neigpt(4,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          outer_b_flx_h(4,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(3,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(4,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          radPtch,radNum)
                              case(5) ! triangle 5
                                 xy_Flux_hex(j,i,saved_nodel_hex(ci,cj),gr,z_ax) =      &
                                       determine_2d_flux_hex(                      &
                                          x_y_u_p_coord_fa_hex(j,i,:),             &
                                          Fl_momx(5,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_momy(5,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_avg(5,gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          central_flux_h(gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          corner_flux_h(neigpt(6,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          corner_flux_h(neigpt(5,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          outer_b_flx_h(5,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(4,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(5,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          radPtch,radNum)
                              case(6) ! triangle 6
                                 xy_Flux_hex(j,i,saved_nodel_hex(ci,cj),gr,z_ax) =      &
                                       determine_2d_flux_hex(                      &
                                          x_y_u_p_coord_fa_hex(j,i,:),             &
                                          Fl_momx(6,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_momy(6,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          Fl_avg(6,gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          central_flux_h(gr,saved_nodel_hex(ci,cj),z_ax),  &
                                          corner_flux_h(neigpt(1,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          corner_flux_h(neigpt(6,saved_nodel_hex(ci,cj)),gr,z_ax), &
                                          outer_b_flx_h(6,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(5,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          inner_b_flx_h(6,gr,saved_nodel_hex(ci,cj),z_ax), &
                                          radPtch,radNum)
                              end select
                           end do
                        end do
                     end if
                  end do
               end do
            end do
         end do
         write(*,*) '[subroutine calculate_pin_power_hex] Finished calculating Ф(x,y)'

         write(*,*) '[subroutine calculate_pin_power_hex] Start restoring flux #/cm3'
       !  call restore_flux_hex
         write(*,*) '[subroutine calculate_pin_power_hex] Flux restored to #/cm3'
         write(*,*) '[subroutine calculate_pin_power_hex] run calculate_power_in_pins_hex'
         do z_ax=1,Nz
            call calculate_power_in_pins_hex(nassy,actual_xy(2), &
                                             actual_xy(1),n_group,z_ax)
         end do
         write(*,*) '[subroutine calculate_pin_power_hex] run multiply_pow_and_hff_hex'
         do z_ax=1,Nz
            call multiply_pow_and_hff_hex(nassy,actual_xy(2),actual_xy(1),z_ax)
         end do

         ! level balancing with Core_Power
         call restore_power_hex

         ! calculate Fq,Fxy,Fz
         call get_Fq_Fxy_Fz_values_hex &
               (xy_pow_hex, lin_pow_hex, &
                  xy_Fq_hex,xy_Fxy_hex,xy_Fz_hex,avg_pin_mlinpow)
         ! calculate FdH
         call get_FdH_values_hex &
               (xy_pow_hex, &
                  integral_pin_pow_hex,xy_FdH_hex,avg_pin_integpow)


         !!!! move to a separate subroutine
         write(*,*) '============================================================='
         write(*,*) '============================================================='
         max_Fq_hex = maxval(xy_Fq_hex)
         max_Fq_hex_loc = maxloc(xy_Fq_hex)
         write(*,*) 'Max Fq value = ',max_Fq_hex
         write(*,*) 'Max Fq location = ',max_Fq_hex_loc
         write(*,*) 'Average pin mesh linear power = ',avg_pin_mlinpow
         write(*,*) 'Max pin mesh linear power = ', &
            lin_pow_hex(max_Fq_hex_loc(1),max_Fq_hex_loc(2), &
               max_Fq_hex_loc(3),max_Fq_hex_loc(4))
         write(*,*) 'Actual Max pin mesh linear power = ',maxval(lin_pow_hex)
         write(*,*) 'Actual Location of Max pin mesh linear power = ', &
                           maxloc(lin_pow_hex)
         write(*,*) '============================================================='
         write(*,*) '============================================================='

         max_Fxy_hex = maxval(xy_Fxy_hex)
         max_Fxy_hex_loc = maxloc(xy_Fxy_hex)
         write(*,*) 'Max Fxy value = ',max_Fxy_hex
         write(*,*) 'Max Fxy location = ',max_Fxy_hex_loc
         write(*,*) '============================================================='
         write(*,*) '============================================================='

         max_Fz_hex = maxval(xy_Fz_hex)
         max_Fz_hex_loc = maxloc(xy_Fz_hex)
         write(*,*) 'Max Fz value = ',max_Fz_hex
         write(*,*) 'Max Fz location = ',max_Fz_hex_loc
         write(*,*) '============================================================='
         write(*,*) '============================================================='

         max_FdH_hex = maxval(xy_FdH_hex)
         max_FdH_hex_loc = maxloc(xy_FdH_hex)
         write(*,*) 'Max FdH value = ',max_FdH_hex
         write(*,*) 'Max FdH location = ',max_FdH_hex_loc
         write(*,*) 'Average pin rod linear power = ',avg_pin_integpow
         write(*,*) 'Max pin rod linear power = ', &
            integral_pin_pow_hex(max_FdH_hex_loc(1), &
               max_FdH_hex_loc(2),max_FdH_hex_loc(3))
         write(*,*) 'Actual Max pin rod linear power = ',maxval(integral_pin_pow_hex)
         write(*,*) 'Actual Location of Max pin rod linear power = ', &
                           maxloc(integral_pin_pow_hex)
         write(*,*) '============================================================='
         write(*,*) '============================================================='


         ! plot axial layers of the full core
         ! xy_pow_hex(j,i,saved_nodel_hex(ci,cj),z_ax)
         ! xy_Fq_hex(j,i,saved_nodel_hex(ci,cj),z_ax)
         ! xy_FdH_hex(j,i,saved_nodel_hex(ci,cj))
         ! integral_pin_pow_hex(j,i,saved_nodel_hex(ci,cj))
         if(.not.allocated(lft_cnt)) allocate(lft_cnt(size(saved_nodel_hex,2)))
         if(.not.allocated(rgt_cnt)) allocate(rgt_cnt(size(saved_nodel_hex,2)))
         call adjust_fa_locations_for_plot & !align the core to the center(x,y)
                           (saved_nodel_hex,xmin,xmax,ymin,ymax,lft_cnt,rgt_cnt)
         prefix = 'Pin_Power_Meshes'
         do z_ax=1,Nz
            call plot_core_axial_layers &
            (z_ax,prefix,xy_pow_hex,xmin,xmax, &
               ymin,ymax,actual_xy,lft_cnt,rgt_cnt)
         end do
         prefix = 'Table_of_FAs'
         do z_ax=1,Nz
            call plot_individual_fas &
               (z_ax,prefix,xy_pow_hex,xmin,xmax, &
               ymin,ymax,actual_xy)
         end do

         if (allocated(lft_cnt             )) deallocate(lft_cnt             )
         if (allocated(rgt_cnt             )) deallocate(rgt_cnt             )
         if (allocated(xy_Fq_hex           )) deallocate(xy_Fq_hex           )
         if (allocated(xy_FdH_hex          )) deallocate(xy_FdH_hex          )
         if (allocated(lin_pow_hex         )) deallocate(lin_pow_hex         )
         if (allocated(integral_pin_pow_hex)) deallocate(integral_pin_pow_hex)
         if (allocated(xy_Fxy_hex          )) deallocate(xy_Fxy_hex          )
         if (allocated(xy_Fz_hex           )) deallocate(xy_Fz_hex           )
         call deallocate_pin_pow_hex
       !  deallocate(saved_nodel_hex_adj)

      end subroutine calculate_pin_power_hex

      subroutine adjust_fa_locations_for_plot &
            (old_svd_ndl_hx,xmin,xmax,ymin,ymax,lft_cnt,rgt_cnt)

         implicit none
         integer,allocatable,dimension(:,:)     :: old_svd_ndl_hx
         integer                                :: ci,cj,xmin,xmax,ymin,ymax
         integer                                :: i,j,k
         integer,allocatable,dimension(:)       :: lft_cnt,rgt_cnt

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [adjust_fa_locations_for_plot] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [adjust_fa_locations_for_plot] in Mod_PinPow_Hex'
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

         do cj = ymin, ymax
            do ci = xmin, xmax
               if (old_svd_ndl_hx(ci,cj).eq.0) then
                  if (old_svd_ndl_hx(ci+1,cj).ne.0) then
                     do i = ci+1,xmax+1
                        old_svd_ndl_hx(i-1,cj) = old_svd_ndl_hx(i,cj)
                     end do
                  end if
               end if
            end do
            ! aligning the FAs around the middle of the X axis (with an uncertainty of 2)
            do while(.true.)
               lft_cnt(cj) = 0 ; rgt_cnt(cj) = 0
               do ci = xmin,xmax
                  if (old_svd_ndl_hx(ci,cj).eq.0) then
                     lft_cnt(cj) = lft_cnt(cj) + 1
                  else
                     exit
                  end if
               end do
               do ci = xmax,xmin,-1
                  if (old_svd_ndl_hx(ci,cj).eq.0) then
                     rgt_cnt(cj) = rgt_cnt(cj) + 1
                  else
                     exit
                  end if
               end do
               if ((lft_cnt(cj)-rgt_cnt(cj)).ge.1) then
                  do ci = xmin,xmax+1
                     old_svd_ndl_hx(ci,cj) = old_svd_ndl_hx(ci+1,cj)
                  end do
                  cycle
               else if ((lft_cnt(cj)-rgt_cnt(cj)).le.-2) then
                  do ci = xmax+1,xmin,-1
                     old_svd_ndl_hx(ci,cj) = old_svd_ndl_hx(ci-1,cj)
                  end do
                  cycle
               else
                  exit
               end if
            end do
         end do

      end subroutine adjust_fa_locations_for_plot

      subroutine plot_core_axial_layers &
            (z_ax,prefix,what_to_print,xmin,xmax,ymin,ymax,actual_xy,lft_cnt,rgt_cnt)
         use Inc_File,     only: Name_INP, Len_INP
         use Inc_PinPow_Hex, only: HFF_hex,saved_nodel_hex
         implicit none
         integer                    :: z_ax,i,j,k
         integer                    :: ci,cj,xmin,xmax,ymin,ymax
         integer                    :: actual_xy(2)
         character(200)             :: filename,print_name_tmp
         character(50)              :: prefix
         character(5)               :: z_axis_string
         real(8),allocatable,dimension(:,:,:,:) :: what_to_print
         character(10000),allocatable,dimension(:) :: long_line
         character(10000)                          :: empty_line
         character(10)                             :: point_flx
         integer,allocatable,dimension(:)       :: lft_cnt,rgt_cnt
         real(8)                    :: adj_factor
         ! just print out the result

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [plot_core_axial_layers] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [plot_core_axial_layers] in Mod_PinPow_Hex'
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

         adj_factor = 1d0
         if (maxval(what_to_print(:,:,:,z_ax)).ge.1000000) then
            adj_factor = 1000000d0
         else if (maxval(what_to_print(:,:,:,z_ax)).ge.100000) then
            adj_factor = 100000d0
         else if (maxval(what_to_print(:,:,:,z_ax)).ge.10000) then
            adj_factor = 10000d0
         else if (maxval(what_to_print(:,:,:,z_ax)).ge.1000) then
            adj_factor = 1000d0
         else
            adj_factor = 1d0
         end if

         filename = Name_INP(1:Len_INP)//prefix//'.dat'
         print_name_tmp = ''
         write(z_axis_string,'(I3)') z_ax
         print_name_tmp = Name_INP(1:Len_INP)//'_'//trim(prefix)//'_'//trim(z_axis_string)
         call open_file_hex(print_name_tmp,44)
         empty_line = ''
         do ci = xmin, xmax
            do j = 1, actual_xy(2)
               empty_line = trim(empty_line)//',0000000'
            end do
         end do
         if(.not.allocated(long_line)) allocate(long_line(actual_xy(1)))
         do cj = ymin, ymax
            long_line = ''
            ! adding half of zero FA to align corners
            if ((lft_cnt(cj).lt.rgt_cnt(cj))) then
               do i = 1, actual_xy(1)
                  do j = 1,actual_xy(2)/2 + actual_xy(2)/8 + 1
                     long_line(i) = trim(long_line(i))//','//'0000000'
                  end do
               end do
            end if
            ! print out the value of adjusted saved_nodel_hex(ci,cj)
            do ci = xmin, xmax
               ! printing out the 2D fluxes to a file for Python draw
               do i = 1, actual_xy(1)
                  do j = 1, actual_xy(2)
                     if (saved_nodel_hex(ci,cj).eq.0) then
                        write(point_flx,'(a7)') '0000000'
                     else
                        if (HFF_hex(i,j,saved_nodel_hex(ci,cj)).lt.0.0) then
                           write(point_flx,'(a7)') '0000000'
                        else
                           write(point_flx,'(f7.3)') &
                              (what_to_print(j,i,saved_nodel_hex(ci,cj),z_ax)/adj_factor)
                        end if
                     end if
                     long_line(i) = trim(long_line(i))//','//trim(point_flx)
                     if(j.eq.actual_xy(2)) then ! separate FAs vertical line (trail)
                        do k = 1,actual_xy(2)/4
                           long_line(i) = trim(long_line(i))//','//'0000000'
                        end do
                     end if
                  end do
               end do
            end do
            do i = 1, actual_xy(1)
               write(44,'(A)') trim(long_line(i))
               empty_line = trim(empty_line)//',0000000' ! separate FAs bottom line
            end do
            write(44,'(A)') trim(empty_line)
         end do
         call close_file_hex(44)

         if(allocated(long_line)) deallocate(long_line)

      end subroutine plot_core_axial_layers

      subroutine plot_individual_fas &
         (z_ax,prefix,what_to_print,xmin,xmax,ymin,ymax,actual_xy)
         use Inc_File,     only: Name_INP, Len_INP
         use Inc_PinPow_Hex, only: HFF_hex,saved_nodel_hex,&
                              xy_Fq_hex,xy_Fxy_hex,xy_Fz_hex,xy_FdH_hex
         use Inc_TPEN, only: nassy
         implicit none
         integer                    :: z_ax,i,j,faa
         integer                    :: ci,cj,xmin,xmax,ymin,ymax
         real(8),allocatable,dimension(:,:,:,:) :: what_to_print
         integer                    :: actual_xy(2)
         character(10)              :: point_flx
         character(10000)           :: empty_line
         character(50)              :: prefix
         character(200)             :: filename,print_name_tmp
         character(5)               :: z_axis_string
         real(8)                    :: adj_factor

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [plot_individual_fas] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [plot_individual_fas] in Mod_PinPow_Hex'
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

         adj_factor = 1d0
         if (maxval(what_to_print(:,:,:,z_ax)).ge.1000000) then
            adj_factor = 1000000d0
         else if (maxval(what_to_print(:,:,:,z_ax)).ge.100000) then
            adj_factor = 100000d0
         else if (maxval(what_to_print(:,:,:,z_ax)).ge.10000) then
            adj_factor = 10000d0
         else if (maxval(what_to_print(:,:,:,z_ax)).ge.1000) then
            adj_factor = 1000d0
         else
            adj_factor = 1d0
         end if

        ! what_to_print(:,:,:,z_ax) = what_to_print(:,:,:,z_ax)/adj_factor

         filename = Name_INP(1:Len_INP)//prefix//'.dat'
         print_name_tmp = ''
         write(z_axis_string,'(I3)') z_ax
         print_name_tmp = Name_INP(1:Len_INP)//'_'//trim(prefix)//'_'//trim(z_axis_string)
         call open_file_hex(print_name_tmp,44)
         do faa = 1,nassy
            call print_max_fa_peaking_factors(44,faa,z_ax)
            write(44,'(A30,(1p,E10.3))') ' Power Values are multiple of:',adj_factor
            call advance_yes_file_plot(44)
            write(44,'(A7)',advance='no') '       '
            do j = 1, actual_xy(2)
              ! write(point_flx,'(I2)') j
               write(44,'(I2)',advance='no') j
               write(44,'(A6)',advance='no') '      '
            end do
            call advance_yes_file_plot(44)
            do i = 1, actual_xy(1)
               do j = 1, actual_xy(2)
                  if (HFF_hex(i,j,faa).lt.0.0) then
                     write(point_flx,'(a7)') '-------'!'0000000'
                  else
                     write(point_flx,'(f7.3)') &
                        (what_to_print(j,i,faa,z_ax)/adj_factor)
                  end if
                  if (j.eq.1) then
                     write(44,'(I2)',advance='no') i
                     write(44,'(A5)',advance='no') '     '
                  end if
                  if(.not.j.eq.actual_xy(2)) then
                     write(44,'(A)',advance='no') trim(point_flx)//'|'!','
                  else
                     write(44,'(A)',advance='no') trim(point_flx)
                  end if
               end do
               call advance_yes_file_plot(44)
            end do

            call advance_yes_file_plot(44)
         end do
         call close_file_hex(44)

      end subroutine plot_individual_fas

      ! @$^ Siarhei - this subroutine for adjusting Write_O to plot hex maps
      subroutine preprocess_hex_2d_plot&
            (file_id,input_2d_array,inp_name,inp_units,bu_step,bu_val,acc,zero_el)
         use Inc_RP, only: I_FARF_1N_2D
         implicit none
         real(8),allocatable,dimension(:,:) :: input_2d_array
         integer                       :: i,j,ci,cj,bu_step,file_id
         integer                       :: dim_x,dim_y,old_dim_x,old_dim_y
         integer                       :: count_x,count_max
         integer,allocatable,dimension(:)   :: count_rows,count_cols
         character(200)    :: nameee
         character(50)     :: inp_name,inp_units
         logical           :: start_count,stop_count,zero_el
         real(8)           :: bu_val
         character(1)      :: acc
         dim_x = 0 ; dim_y = 0

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [preprocess_hex_2d_plot] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [preprocess_hex_2d_plot&] in Mod_PinPow_Hex'
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

         old_dim_y = size(input_2d_array,1)-1
         old_dim_x = size(input_2d_array,2)-1
         dim_y = size(I_FARF_1N_2D,1)
         dim_x = size(I_FARF_1N_2D,2)
         if(.not.allocated(count_rows))      allocate(count_rows(old_dim_y))
         if(.not.allocated(count_cols))      allocate(count_cols(old_dim_x))
         write(*,*) 'Plot Dimensions (x,y) - old:',old_dim_x,old_dim_y
         write(*,*) 'I_FARF_1N_2D (x,y):',dim_x,dim_y
         count_rows(:) = 0 ; count_cols(:) = 0

         ! preprocess number of elements (for chess order)
         count_x = 0 ; count_max = 0
         start_count = .false. ; stop_count = .true.
         do i = 1,old_dim_y
            count_x = 0
            do j = 1,old_dim_x
              ! if (abs(input_2d_array(i,j)).gt.1D-10) start_count = .true.
               if (I_FARF_1N_2D(i,j).gt.0) start_count = .true.
               if (start_count) then
                 ! if (abs(input_2d_array(i,j)).lt.1D-10.and.&
                 !       abs(input_2d_array(i,j+1)).lt.1D-10) then
                  if (I_FARF_1N_2D(i,j).eq.0) then
                           start_count = .false.
                  else
                        count_rows(i) = count_rows(i) + 1
                        count_cols(j) = count_cols(j) + 1
                  end if
               end if
            end do
            if (count_rows(i).gt.count_max) count_max = count_rows(i)
            write(*,*) '(preprocess_hex_2d_plot) Line ',&
               i,'# elements:',count_x,inp_name
         end do
         write(*,*) '(preprocess_hex_2d_plot) Max_count:',count_max
        ! nameee = 'TMP_Plot_T_FUEL.out'
        ! call open_file_hex(nameee,file_id)
         call advance_yes_file_plot(file_id)
         count_x = 1
         write(file_id,'(A14,I2,A16,F7.'//acc//',A8)') ' Burnup Step: ',bu_step,&
                                         ', Burnup Value: ',bu_val,' GWD/MTU'
         write(file_id,'(A12,A50)') ' Parameter: ',inp_name
         if (zero_el) then
            write(file_id,'(A16,F7.'//acc//',A1)',advance='no') ' Maximum Value: ',&
                                                      0d0,' '
         else
            if (maxval(input_2d_array).lt.1D-10) then
               write(file_id,'(A16,F7.'//acc//',A1)',advance='no') ' Minimum Value: ',&
                                                   minval(input_2d_array),' '
            else
               write(file_id,'(A16,F7.'//acc//',A1)',advance='no') ' Maximum Value: ',&
                                                      maxval(input_2d_array),' '
            end if
         end if
         write(file_id,'(A)') trim(inp_units)
         write(file_id,'(A)') ' - - - - - - - - - - - - - - - - - - - - - - -'
         call advance_yes_file_plot(file_id)
         do j = 0,(2*count_max-1)
            if (j.eq.0) then
               write(file_id,'(A8)',advance='no') 'Row/Col '
            else
               write(file_id,'(I4)',advance='no') j
               write(file_id,'(A4)',advance='no') '    '
            end if
         end do
         call advance_yes_file_plot(file_id)
         do j = 0,(2*count_max-1)
            write(file_id,'(A8)',advance='no') '--------'
         end do
         call advance_yes_file_plot(file_id)
         do i = 1,old_dim_y
            if(count_rows(i).gt.1D-10) then
               write(file_id,'(I4)',advance='no') count_x
               write(file_id,'(A4)',advance='no') '    '
               count_x = count_x + 1
               do j = 1,ceiling(real(count_max-count_rows(i))/1)
                  write(file_id,'(A8)',advance='no') '        '
               end do

               do j = 1,old_dim_x
                  if(count_cols(j).gt.1D-10) then
                    ! if (abs(input_2d_array(i,j)).gt.1D-10) start_count = .true.
                     if (I_FARF_1N_2D(i,j).gt.0) start_count = .true.
                     if (start_count) then
                       ! if (abs(input_2d_array(i,j)).lt.1D-10.and.&
                       !       abs(input_2d_array(i,j+1)).lt.1D-10) then
                        if (I_FARF_1N_2D(i,j).eq.0) then
                                 start_count = .false.
                        else
                          ! if (abs(input_2d_array(i,j)).gt.1D-10) then
                           if (I_FARF_1N_2D(i,j).gt.0) then
                              if (zero_el) then
                                 write(file_id,'(F7.'//acc//')',advance='no') 0d0
                              else
                                 write(file_id,'(F7.'//acc//')',advance='no') input_2d_array(i,j)
                              end if
                              write(file_id,'(A1)',advance='no') ' '
                              write(file_id,'(A8)',advance='no') '        '
                           else
                              write(file_id,'(A8)',advance='no') '--||||--'
                              write(file_id,'(A8)',advance='no') '        '
                           end if
                        end if
                     end if
                  end if
               end do
               call advance_yes_file_plot(file_id)
            end if
         end do
         do j = 0,(2*count_max-1)
            write(file_id,'(A8)',advance='no') '--------'
         end do
         call advance_yes_file_plot(file_id)
         do j = 0,(2*count_max-1)
            write(file_id,'(A8)',advance='no') '--------'
         end do
         call advance_yes_file_plot(file_id)
        ! call close_file_hex(file_id)

         if(allocated(count_rows))      deallocate(count_rows)
         if(allocated(count_cols))      deallocate(count_cols)

      end subroutine preprocess_hex_2d_plot

      subroutine print_max_fa_peaking_factors(file_id,fa_numb,z_ax)
         use Inc_PinPow_Hex, only: xy_Fq_hex,xy_Fxy_hex,xy_FdH_hex
         implicit none
         integer           :: file_id,fa_numb,z_ax
         integer           :: Fq_loc(2),Fxy_loc(2),FdH_loc(2)

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [print_max_fa_peaking_factors] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [print_max_fa_peaking_factors] in Mod_PinPow_Hex'
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


         Fq_loc  = maxloc(xy_Fq_hex(:,:,fa_numb,z_ax))
         Fxy_loc = maxloc(xy_Fxy_hex(:,:,fa_numb,z_ax))
         FdH_loc = maxloc(xy_FdH_hex(:,:,fa_numb))

         write(file_id,*) '- - - - - - - - - - - - - - - - - -'
         write(file_id,'(A15,I3,A15,I3)') &
               ' Fuel assembly:',fa_numb,' | axial layer:',z_ax

         write(file_id,'(A9,F7.3,A21)',advance='no') &
               ' Max  Fq:',maxval(xy_Fq_hex(:,:,fa_numb,z_ax)),'   | Location (x,y): '
         write(file_id,*) Fq_loc
         write(file_id,'(A9,F7.3,A21)',advance='no') &
               ' Max Fxy:',maxval(xy_Fxy_hex(:,:,fa_numb,z_ax)),'   | Location (x,y): '
         write(file_id,*) Fxy_loc
         write(file_id,'(A9,F7.3,A21)',advance='no') &
               ' Max FdH:',maxval(xy_FdH_hex(:,:,fa_numb)),'   | Location (x,y): '
         write(file_id,*) FdH_loc
       !  write(file_id,'(A)',advance='no') '  '


      end subroutine print_max_fa_peaking_factors

      subroutine restore_flux_hex
         use Mod_GetSome
         use Inc_3D, only: Power,Avg_Power
         use Inc_PinPow_Hex, only: xy_Flux_hex
         use Inc_TH!, only: LevFactor, Avg_CorePower
         implicit none
         logical                    :: Flag_Conv_PowLev

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [restore_flux_hex] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [restore_flux_hex] in Mod_PinPow_Hex'
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

         Flag_Conv_PowLev = .false.
         write(*,*) '[restore_flux_hex] Core_Power = ',Core_Power
         do while(.not.Flag_Conv_PowLev)
            call Get_POW
            call Get_LinPOW
            call Get_Avg(Avg_Power, Power, 0)

            LevFactor = Avg_CorePower/Avg_Power

            call LevelBalancing
            xy_Flux_hex = xy_Flux_hex * LevFactor

            if (abs(LevFactor - 1d0).gt.1.d-3) then
               Flag_Conv_PowLev = .false.
            else
               Flag_Conv_PowLev = .true.
               exit
            end if

            call Get_POW
            call Get_LinPOW
            call Get_Avg(Avg_Power, Power, 0)
            call get_normal_power
         end do


      end subroutine restore_flux_hex


      subroutine restore_power_hex
         use Mod_GetSome
         use Inc_PinPow_Hex, only: radNum,xy_pow_hex
         use Inc_3D, only: Power
         use Inc_TH!, only: LevFactor, Avg_CorePower
         use Inc_TPEN,     only: nassy
         use Inc_Geometry, only: Nz
         implicit none
         logical                    :: Flag_Conv_PowLev
         real(8)                    :: init_nodal_power(nassy,Nz)
         integer                    :: actual_xy(2)
         real(8)                    :: local_Avg_Power,local_Tot_Power

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [restore_power_hex] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [restore_power_hex] in Mod_PinPow_Hex'
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

         actual_xy(1) = 4*(radNum-1) + 1 ! i - x_ax
         actual_xy(2) = 2*(radNum-1) + 1 ! j - y_ax
         local_Avg_Power = 0d0

         Flag_Conv_PowLev = .false.
         write(*,*) '[restore_flux_hex] Core_Power = ',Core_Power
         do while(.not.Flag_Conv_PowLev)
           ! call Get_POW
           ! call Get_LinPOW
            init_nodal_power = 0d0
            local_Tot_Power = 0d0
            local_Tot_Power = update_nodal_power( &
           ! init_nodal_power = update_nodal_power( &
                        xy_pow_hex,actual_xy(1),actual_xy(2),Nz)
           ! call Get_Avg(local_Avg_Power, init_nodal_power, 0)

            LevFactor = Core_Power/local_Tot_Power
           ! LevFactor = Avg_CorePower/local_Avg_Power

           ! call LevelBalancing
            xy_pow_hex = xy_pow_hex * LevFactor

            if (abs(LevFactor - 1d0).gt.1.d-3) then
               Flag_Conv_PowLev = .false.
            else
               Flag_Conv_PowLev = .true.
               exit
            end if

           ! call Get_POW
           ! call Get_LinPOW
           ! call Get_Avg(Avg_Power, Power, 0)
           ! call get_normal_power
         end do


      end subroutine restore_power_hex

      function update_nodal_power(tmp_pin_power,x_ax,y_ax,z_ax) &
               result(hex_core_power) !updated_node_power)
         use Inc_PinPow_Hex, only: HFF_hex
         use Inc_TPEN,     only: nassy
         use Inc_Geometry, only: Nz
         implicit none
         real(8),dimension(:,:,:,:)           :: tmp_pin_power
         integer                              :: i,j,n_fa,n_z
         integer                              :: x_ax,y_ax,z_ax
         real(8)                              :: updated_node_power(nassy,Nz)
         real(8)                              :: hex_core_power

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [update_nodal_power] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [update_nodal_power] in Mod_PinPow_Hex'
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

         updated_node_power = 0d0
         hex_core_power = 0d0
         do concurrent (n_z=1:z_ax,n_fa=1:nassy,i = 1:x_ax,j = 1:y_ax)
      !   do n_z=1,z_ax
      !      do n_fa=1,nassy
      !         do i = 1, x_ax
      !            do j = 1, y_ax
                    ! if (tmp_pin_power(j,i,n_fa,n_z).gt.0.0) then
            if (HFF_hex(i,j,n_fa).ge.0.0) then
               hex_core_power = hex_core_power + &
                                tmp_pin_power(j,i,n_fa,n_z)

            end if
      !            end do
      !         end do
      !      end do
         end do


      end function update_nodal_power

      ! get Fq and Fxy for all pin meshes
      subroutine get_Fq_Fxy_Fz_values_hex &
               (init_pin_mesh_power, pin_mesh_linear_power, &
                  Fq_mesh_pin,Fxy_mesh_pin,Fz_layer,avg_pin_mesh_lin_pow)

         use Inc_PinPow_Hex, only: radNum,HFF_hex
         use Inc_3D, only: Power
         use Inc_TH!, only: LevFactor, Avg_CorePower
         use Inc_TPEN,     only: nassy
         use Inc_Geometry, only: Nz,GridSize_z!(z_ax)
         implicit none
         logical                    :: Flag_Conv_PowLev
         real(8),dimension(:,:,:,:) :: init_pin_mesh_power
         integer                    :: actual_xy(2)
         integer                    :: i,j,n_fa,n_z,count_mesh_pins,count_mesh
         integer                    :: x_ax,y_ax,z_ax
         real(8),dimension(:)      ,intent(inout) :: Fz_layer
         real(8),dimension(:,:,:,:),intent(inout) :: pin_mesh_linear_power
         real(8)                   ,intent(inout) :: avg_pin_mesh_lin_pow
         real(8)                                  :: avg_z_pin_mesh_lin_pow(Nz)
         real(8)                                  :: z_layer_pow(Nz)
         real(8)                                  :: total_z_layer_pow
         real(8),dimension(:,:,:,:),intent(inout) :: Fq_mesh_pin,Fxy_mesh_pin

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [get_Fq_Fxy_Fz_values_hex] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_Fq_Fxy_Fz_values_hex] in Mod_PinPow_Hex'
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

         actual_xy(1) = 4*(radNum-1) + 1 ! i - x_ax
         actual_xy(2) = 2*(radNum-1) + 1 ! j - y_ax
         x_ax = actual_xy(1) ; y_ax = actual_xy(2) ; z_ax = Nz
         Fq_mesh_pin = 0d0 ; Fxy_mesh_pin = 0d0 ; Fz_layer = 0d0
         count_mesh_pins = 0
         count_mesh = 0
         pin_mesh_linear_power = 0d0
         avg_pin_mesh_lin_pow = 0d0
         avg_z_pin_mesh_lin_pow = 0d0
         z_layer_pow = 0d0

         total_z_layer_pow = 0d0

         do n_z=1,z_ax
            do concurrent (n_fa=1:nassy,i = 1:x_ax,j = 1:y_ax)
      !      do n_fa=1,nassy
      !         do i = 1, x_ax
      !            do j = 1, y_ax
                    ! if (init_pin_mesh_power(j,i,n_fa,n_z).gt.0.0) then
               if (HFF_hex(i,j,n_fa).ge.0.0.and.&
                  init_pin_mesh_power(j,i,n_fa,n_z).gt.0.0) then
                     count_mesh_pins = count_mesh_pins + 1
                     count_mesh = count_mesh + 1
                     pin_mesh_linear_power(j,i,n_fa,n_z) = &
                           init_pin_mesh_power(j,i,n_fa,n_z)/GridSize_z(n_z)
                     avg_pin_mesh_lin_pow = avg_pin_mesh_lin_pow + &
                           pin_mesh_linear_power(j,i,n_fa,n_z)
                     avg_z_pin_mesh_lin_pow(n_z) = avg_z_pin_mesh_lin_pow(n_z) + &
                        pin_mesh_linear_power(j,i,n_fa,n_z)
                     z_layer_pow(n_z) = z_layer_pow(n_z) + &
                           init_pin_mesh_power(j,i,n_fa,n_z)
                     if (pin_mesh_linear_power(j,i,n_fa,n_z).le.0.0) &
                        write(*,*) "Fq_pin_mesh",j,i,n_fa,n_z, &
                           pin_mesh_linear_power(j,i,n_fa,n_z)

               end if
      !            end do
      !         end do
            end do
            total_z_layer_pow = total_z_layer_pow + z_layer_pow(n_z)
            avg_z_pin_mesh_lin_pow(n_z) = &
               avg_z_pin_mesh_lin_pow(n_z)/count_mesh
            count_mesh = 0
         end do
         avg_pin_mesh_lin_pow = avg_pin_mesh_lin_pow/count_mesh_pins

         Fq_mesh_pin = pin_mesh_linear_power/avg_pin_mesh_lin_pow
         do concurrent (n_z=1:z_ax)
        ! do n_z=1,z_ax
            Fxy_mesh_pin(:,:,:,n_z) = &
               pin_mesh_linear_power(:,:,:,n_z)/avg_z_pin_mesh_lin_pow(n_z)
            Fz_layer(n_z) = z_layer_pow(n_z)/total_z_layer_pow * z_ax ! multiply for normalization
         end do


      end subroutine get_Fq_Fxy_Fz_values_hex


      subroutine get_FdH_values_hex &
               (init_pin_mesh_power, &
                  pin_integral_linear_power,FdH_pin,avg_pin_pow)

         use Inc_PinPow_Hex, only: radNum,HFF_hex,radPtch
         use Inc_3D, only: Power
         use Inc_TH!, only: LevFactor, Avg_CorePower
         use Inc_TPEN,     only: nassy
         use Inc_Geometry, only: Nz,GridSize_z!(z_ax)
         implicit none
         logical                    :: Flag_Conv_PowLev
         real(8),dimension(:,:,:,:) :: init_pin_mesh_power
         integer                    :: actual_xy(2)
         integer                    :: i,j,n_fa,n_z,count_pins
         integer                    :: x_ax,y_ax,z_ax
         real(8),dimension(:,:,:),intent(inout)   :: pin_integral_linear_power
         real(8)                 ,intent(inout)   :: avg_pin_pow
         real(8),dimension(:,:,:),intent(inout)   :: FdH_pin

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [get_FdH_values_hex] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_FdH_values_hex] in Mod_PinPow_Hex'
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

         actual_xy(1) = 4*(radNum-1) + 1 ! i - x_ax
         actual_xy(2) = 2*(radNum-1) + 1 ! j - y_ax
         x_ax = actual_xy(1) ; y_ax = actual_xy(2) ; z_ax = Nz
         count_pins = 0
         pin_integral_linear_power = 0d0
         avg_pin_pow = 0d0
         FdH_pin = 0d0
         ! remove z-axis mesh
         do concurrent (n_z=1:z_ax)
        ! do n_z=1,z_ax
            pin_integral_linear_power = pin_integral_linear_power(:,:,:) + &
               init_pin_mesh_power(:,:,:,n_z)/1000.0 ! &
                 ! *(3**1.5)*0.5*radPtch*radPtch*GridSize_z(n_z)/1000.0
         end do
         ! calculate average and FdH
         do concurrent (n_fa=1:nassy,i = 1:x_ax,j = 1:y_ax)
      !   do n_fa=1,nassy
      !      do i = 1, x_ax
      !         do j = 1, y_ax
                 ! if (pin_integral_linear_power(j,i,n_fa).gt.0.0) then
            if (HFF_hex(i,j,n_fa).ge.0.0.and.&
               pin_integral_linear_power(j,i,n_fa).gt.0.0) then
                  count_pins = count_pins + 1
                  avg_pin_pow = avg_pin_pow + &
                        pin_integral_linear_power(j,i,n_fa)
                 ! if (pin_integral_linear_power(j,i,n_fa).le.0.0) &
                 !       write(*,*) "FdH_pin_mesh",j,i,n_fa, &
                 !          pin_integral_linear_power(j,i,n_fa)
            end if
      !         end do
      !      end do
         end do

         avg_pin_pow = avg_pin_pow/count_pins
         FdH_pin = pin_integral_linear_power/avg_pin_pow

      end subroutine get_FdH_values_hex


      subroutine calculate_power_in_pins_hex(fa_count,x,y,energy_groups,z_axis)
         use Inc_PinPow_Hex
         use Inc_maXS, only: kap_maXS_f_3D ! (ixy,iz,ig)
         implicit none
         integer                    :: fa_count,x,y,energy_groups,z_axis
         integer                    :: i,j,k,gr

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [calculate_power_in_pins_hex] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [calculate_power_in_pins_hex] in Mod_PinPow_Hex'
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

         do concurrent (k=1:fa_count,i = 1:y,j = 1:x,gr = 1:energy_groups)
      !   do k = 1,fa_count
      !      do i = 1,y
      !         do j = 1,x
      !            do gr = 1,energy_groups
            hom_pow_hex(j,i,k,z_axis) = hom_pow_hex(j,i,k,z_axis) + &
               kap_maXS_f_3D(k,z_axis,gr) *        &
               xy_Flux_hex(j,i,k,gr,z_axis)
      !            end do
      !         end do
      !      end do
         end do

      end subroutine calculate_power_in_pins_hex


      subroutine multiply_pow_and_hff_hex(fa_count,x,y,z_axis)
         use Inc_PinPow_Hex
         USE Inc_TPEN, only: kap_maXS_f_3D_MG ! (ixy,iz,ig)
         implicit none
         integer                    :: fa_count,x,y,z_axis
         integer                    :: i,j,k

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [multiply_pow_and_hff_hex] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [multiply_pow_and_hff_hex] in Mod_PinPow_Hex'
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

         do concurrent (k=1:fa_count,i = 1:y,j = 1:x)
      !   do k = 1,fa_count
      !      do i = 1,y
      !         do j = 1,x
            if(HFF_hex(i,j,k).ge.0) then
               xy_pow_hex(j,i,k,z_axis) = &
                  hom_pow_hex(j,i,k,z_axis) * HFF_hex(i,j,k)
            else
               xy_pow_hex(j,i,k,z_axis) = -7.777
                  end if
      !         end do
      !      end do
         end do


      end subroutine multiply_pow_and_hff_hex

      subroutine set_coordinates_fa_hex(x,y)
         use Inc_PinPow_Hex
         implicit none
         integer                    :: i,j,k
         integer                    :: x,y
         real(8)                    :: cosine_dir,sin_dir
         real(8)                    :: a, R, h, h_minus_R, sqrt3
         real(8),parameter          :: pi = 3.1415926535897932
         ! FA map with xy for each triangle center
         real(8)                    :: central_cell(2)


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [set_coordinates_fa_hex] in Mod_PinPow_Hex'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

         sqrt3 = 3**0.5
         a = radPtch*(radNum-1)!+0.5*radPtch
         R = a/sqrt3
         h = a*sqrt3/2
         h_minus_R = h - R

         sin_dir = 0.5*radPtch
         cosine_dir = sin_dir*sqrt3

         triangle_height_hex = h ! save triangle height to global
         center_x_triangle_hex = h_minus_R ! save the center coordinate of x (triangle)
         central_cell(1) = h ! x
         central_cell(2) = a ! y

         ! assign global scale in cm (real distance units)
         do concurrent (i = 1:y,j = 1:x)
      !   do i = 1,y
      !      do j = 1,x
            x_coord_fa_hex(j,i) = (j-1)*cosine_dir+0.5*radPtch
            y_coord_fa_hex(j,i) = (i-1)*sin_dir+0.5*radPtch
            ! find radial distance from the center
            rad_coord_fa_hex(j,i) = (central_cell(1)-x_coord_fa_hex(j,i))**2
            rad_coord_fa_hex(j,i) = rad_coord_fa_hex(j,i) + &
               (central_cell(2)-y_coord_fa_hex(j,i))**2
            rad_coord_fa_hex(j,i) = rad_coord_fa_hex(j,i)**0.5
            ! find sin(angle) to the horizontal line y = a
            if(rad_coord_fa_hex(j,i).eq.0) then
               sin_coord_fa_hex(j,i) = -4444
               ang_coord_fa_hex(j,i) = -4444
            else
               sin_coord_fa_hex(j,i) = central_cell(2) - y_coord_fa_hex(j,i)
               sin_coord_fa_hex(j,i) = sin_coord_fa_hex(j,i)/rad_coord_fa_hex(j,i)
               ang_coord_fa_hex(j,i) = asin(sin_coord_fa_hex(j,i))
               if (x_coord_fa_hex(j,i).gt.central_cell(1)) then
                  ! arrange angle from -pi/2 to 3pi/2
                  ang_coord_fa_hex(j,i) = pi - ang_coord_fa_hex(j,i)
               end if
            end if

      !      end do
         end do
         write(*,*) "[Mod_PinPow_Hex -- set_coordinates_fa_hex] ",&
            'Finished setting up the fuel pin coordinates (x,y,Radius,Angle)'

      end subroutine set_coordinates_fa_hex

      ! seems to be not used in the module !$ <--- check if still needed
      function find_distance_2d_points(x1,y1,x2,y2) &
               result(dist)
         implicit none
         real(8)              :: x1,y1,x2,y2
         real(8)              :: dist

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [find_distance_2d_points] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [find_distance_2d_points] in Mod_PinPow_Hex'
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

         dist = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)
         dist = dist**0.5

      end function find_distance_2d_points


      function define_location_hex(angle) &
                                 result(location)
         implicit none
         real(8)              :: angle
         real(8),parameter    :: pi = 3.1415926535897932
         integer              :: location

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [define_location_hex] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [define_location_hex] in Mod_PinPow_Hex'
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

         if(angle.gt.(-pi/6).and.angle.le.(pi/6)) then
            location = 1          ! triangle 1
         else if (angle.gt.(pi/6).and.angle.le.(pi/2)) then
            location = 2          ! triangle 2
         else if (angle.gt.(pi/2).and.angle.le.(5*pi/6)) then
            location = 3          ! triangle 3
         else if (angle.gt.(5*pi/6).and.angle.le.(7*pi/6)) then
            location = 4          ! triangle 4
         else if (angle.gt.(7*pi/6).and.angle.lt.(9*pi/6)) then
            location = 5          ! triangle 5
         else
            location = 6          ! triangle 6
         end if

      end function define_location_hex

      ! convert any point in full-size FA to a single triangle point (x,y,u,p)
      function get_x_y_u_p_hex(x,y,angle,radius,height,center_x) &
               result(x_y_u_p)

         implicit none
         real(8)              :: x_y_u_p(4) ! order as given in the name
         real(8)              :: x,y,angle,radius,new_angle,height
         real(8)              :: center_x
         real(8),parameter    :: pi = 3.1415926535897932

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [get_x_y_u_p_hex] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [get_x_y_u_p_hex] in Mod_PinPow_Hex'
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

         if(angle.gt.(-pi/6).and.angle.le.(pi/6)) then
            new_angle = angle          ! triangle 1
         else if (angle.gt.(pi/6).and.angle.le.(pi/2)) then
            new_angle = angle - pi/3   ! triangle 2
         else if (angle.gt.(pi/2).and.angle.le.(5*pi/6)) then
            new_angle = angle - 2*pi/3 ! triangle 3
         else if (angle.gt.(5*pi/6).and.angle.le.(7*pi/6)) then
            new_angle = angle - 3*pi/3 ! triangle 4
         else if (angle.gt.(7*pi/6).and.angle.lt.(9*pi/6)) then
            new_angle = angle - 4*pi/3 ! triangle 5
         else
            new_angle = angle + pi/3   ! triangle 6
         end if
         x_y_u_p(1) = center_x - (height - radius*cos(new_angle))       ! new x
         x_y_u_p(2) = -radius*sin(new_angle)                            ! new y
         x_y_u_p(3) = 0.5*(height-center_x)-radius*sin(pi/6-new_angle)  ! new u
         x_y_u_p(4) = 0.5*(height-center_x)+radius*sin(-pi/6-new_angle) ! new p

      end function get_x_y_u_p_hex


      ! determine expansion coeffs for all FAs and all nodes
      function determine_2d_flux_hex(x_y_u_p,              &
                              Fl_momx,Fl_momy,             &
                              Fl_avg,Fl_c_x,Fl_c_p,        &
                              Fl_c_u,Fl_s_x,Fl_s_p,Fl_s_u, &
                              radpitch,radnumber) &
               result(flux_2d_hex)
        ! use Inc_PinPow_Hex, only: radPtch,radNum
        ! use Inc_TPEN
        ! use Inc_Option,   only: n_group
         implicit none
         real(8)                    :: flux_2d_hex
         real(8)                    :: x_y_u_p(4)
         real(8)                    :: Fl_momx,Fl_momy,Fl_avg
         real(8)                    :: Fl_c_x,Fl_c_p,Fl_c_u
         real(8)                    :: Fl_s_x,Fl_s_p,Fl_s_u
         real(8)                    :: C_g0,A_gx,A_gy,B_gx
         real(8)                    :: B_gu,B_gp,C_gx,C_gu,C_gp

         integer                    :: i,j,k,radnumber
         real(8)                    :: a, R, a_minus_R, sqrt3,radpitch

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [determine_2d_flux_hex] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [determine_2d_flux_hex] in Mod_PinPow_Hex'
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

         ! a - side of a triangle, R - radius of circumscribed circle
         ! a_minus_R - distance from center to the boundary plane
         sqrt3 = 3**0.5
         a = radpitch*(radnumber-1)
         R = a/sqrt3
         a_minus_R = a - R

        ! smooth 2D flux function from Cho et al
        ! Ф(x,y) = C_g0 + A_gx*x + A_gy*y + B_gx*x*x + &
        !          B_gu*u*u + B_gp*p*p + C_gx*x*x*x + &
        !          C_gu*u*u*u + C_gp*p*p*p

         ! coefficients for the flux function in (x,u,p) coordinates
         C_g0 = 1d0/27d0 * (Fl_c_x+Fl_c_u+Fl_c_p - &
            12*(Fl_s_x+Fl_s_u+Fl_s_p) + 60*Fl_avg) !
         A_gx = 2*sqrt3/(9*a) * ((2*Fl_c_x-Fl_c_u-Fl_c_p) + &
            90*Fl_momx) !
         A_gy = 2/(3*a) * (-(Fl_c_u-Fl_c_p)+30*Fl_momy) !
         B_gx = -4/(3*a*a) * ((Fl_c_u+Fl_c_p)-12*Fl_s_x + &
               10*Fl_avg+60*Fl_momx) !
         B_gu = -4/(3*a*a) * ((Fl_c_x+Fl_c_p)-12*Fl_s_u + &
               10*Fl_avg-30*(Fl_momx+Fl_momy)) !
         B_gp = -4/(3*a*a) * ((Fl_c_x+Fl_c_u)-12*Fl_s_p + &
               10*Fl_avg-30*(Fl_momx-Fl_momy)) !
         C_gx = -16*sqrt3/(9*a*a*a) * (-3*Fl_c_x-(Fl_c_u+Fl_c_p) + &
               9*Fl_s_x+3*(Fl_s_u+Fl_s_p) - 10*Fl_avg-60*Fl_momx) !
         C_gu = -16*sqrt3/(9*a*a*a) * (-3*Fl_c_u-(Fl_c_p+Fl_c_x) + &
               9*Fl_s_u+3*(Fl_s_p+Fl_s_x) - 10*Fl_avg+30*(Fl_momx+Fl_momy)) !
         C_gp = -16*sqrt3/(9*a*a*a) * (-3*Fl_c_p-(Fl_c_x+Fl_c_u) + &
               9*Fl_s_p+3*(Fl_s_x+Fl_s_u) - 10*Fl_avg+30*(Fl_momx-Fl_momy)) !

         ! finally, find the value of 2d flux (for one group)
         flux_2d_hex =  C_g0 + &
                        A_gx*x_y_u_p(1) + &
                        A_gy*x_y_u_p(2) + &
                        B_gx*x_y_u_p(1)*x_y_u_p(1) + &
                        B_gu*x_y_u_p(3)*x_y_u_p(3) + &
                        B_gp*x_y_u_p(4)*x_y_u_p(4) + &
                        C_gx*x_y_u_p(1)*x_y_u_p(1)*x_y_u_p(1) + &
                        C_gu*x_y_u_p(3)*x_y_u_p(3)*x_y_u_p(3) + &
                        C_gp*x_y_u_p(4)*x_y_u_p(4)*x_y_u_p(4)

      end function determine_2d_flux_hex
      !!!!!$ check how to treat Z direction (combine, or calculate each axial level?)

      subroutine deallocate_pin_pow_hex
         use Inc_PinPow_Hex
         implicit none

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '^#^ Entered [deallocate_pin_pow_hex] in Mod_PinPow_Hex'

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [deallocate_pin_pow_hex] in Mod_PinPow_Hex'
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

         if (allocated(HFF_hex))       deallocate(HFF_hex)
         if (allocated(saved_nodel_hex))deallocate(saved_nodel_hex)

        ! if (allocated(C_g0          )) deallocate(C_g0)
        ! if (allocated(A_gx          )) deallocate(A_gx)
        ! if (allocated(A_gy          )) deallocate(A_gy)
        ! if (allocated(B_gx          )) deallocate(B_gx)
        ! if (allocated(B_gu          )) deallocate(B_gu)
        ! if (allocated(B_gp          )) deallocate(B_gp)
        ! if (allocated(C_gx          )) deallocate(C_gx)
        ! if (allocated(C_gu          )) deallocate(C_gu)
        ! if (allocated(C_gp          )) deallocate(C_gp)


       !  if (allocated(Fl_c_x        )) deallocate(Fl_c_x )
       !  if (allocated(Fl_c_p        )) deallocate(Fl_c_p )
       !  if (allocated(Fl_c_u        )) deallocate(Fl_c_u )
       !  if (allocated(Fl_s_x        )) deallocate(Fl_s_x )
       !  if (allocated(Fl_s_p        )) deallocate(Fl_s_p )
       !  if (allocated(Fl_s_u        )) deallocate(Fl_s_u )
         if (allocated(Fl_momx       )) deallocate(Fl_momx)
         if (allocated(Fl_momy       )) deallocate(Fl_momy)
         if (allocated(Fl_avg        )) deallocate(Fl_avg )

         if (allocated(corner_flux_h )) deallocate(corner_flux_h )
         if (allocated(central_flux_h)) deallocate(central_flux_h)
         if (allocated(inner_b_flx_h )) deallocate(inner_b_flx_h)
         if (allocated(outer_b_flx_h )) deallocate(outer_b_flx_h)

         if (allocated(Fl_momx_av    )) deallocate(Fl_momx_av)
         if (allocated(Fl_momy_av    )) deallocate(Fl_momy_av)
         if (allocated(Fl_avg_av     )) deallocate(Fl_avg_av)

         if (allocated(corner_flux_h_av)) deallocate(corner_flux_h_av)
         if (allocated(central_flux_h_av)) deallocate(central_flux_h_av)
         if (allocated(inner_b_flx_h_av)) deallocate(inner_b_flx_h_av)
         if (allocated(outer_b_flx_h_av)) deallocate(outer_b_flx_h_av)

         if (allocated(xy_Flux_hex)) deallocate(xy_Flux_hex)
         if (allocated(x_coord_fa_hex)) deallocate(x_coord_fa_hex)
         if (allocated(y_coord_fa_hex)) deallocate(y_coord_fa_hex)
         if (allocated(p_coord_fa_hex)) deallocate(p_coord_fa_hex)
         if (allocated(u_coord_fa_hex)) deallocate(u_coord_fa_hex)
         if (allocated(rad_coord_fa_hex)) deallocate(rad_coord_fa_hex)
         if (allocated(sin_coord_fa_hex)) deallocate(sin_coord_fa_hex)
         if (allocated(ang_coord_fa_hex)) deallocate(ang_coord_fa_hex)
         if (allocated(x_y_u_p_coord_fa_hex)) deallocate(x_y_u_p_coord_fa_hex)
         if (allocated(location_in_hex)) deallocate(location_in_hex)

         if (allocated(hom_pow_hex )) deallocate(hom_pow_hex)
         if (allocated(xy_pow_hex )) deallocate(xy_pow_hex)


         write(*,*) '[subroutine deallocate_pin_pow_hex] All arrays deallocated'

      end subroutine deallocate_pin_pow_hex
#endif
   end module Mod_PinPow_Hex

 !  program Pinpow_test
 !     use Mod_PinPow_Hex
 !     use Inc_PinPow_Hex
 !     implicit none
 !     character(50)              :: filename
 !     filename = 'sample_triange_one_sixth_input.inp'
 !    ! filename = 'sample_hexagonal_input.inp'
 !    ! filename = 'sample_half_hex_input.inp'
 !     call read_hff_hex(filename)
 !     call deallocate_pin_pow_hex
 !  end program Pinpow_test
