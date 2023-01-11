   ! module for Depletion in Fast Reactor
!//    !
!//    !
!// #ifdef tuan_fr
!//    module FastDepletion
!// 
!// 
!// 
!// #ifdef siarhei_plot
!//       use Inc_Control,only:last_saved,current_sub,dummy_filler
!// #endif
!// 
!//       implicit none
!//    contains
!// 
!//       ! fully read the depletion file and save data
!//       ! to 'all_names' array (Inc_FastDepletion)
!//       subroutine fr_read_depletion_file(name_of_file)
!//          use Inc_FastDepletion
!// 
!//          implicit none
!// 
!//          ! block #1 variables
!// !         character(1)               :: spce ! first space col
!//          character(10)              :: tmpname ! nuclide name
!//          integer                    :: tmpfiflag ! fis branch
!//          integer                    :: dumid    ! dummy value
!//          integer                    :: dumflg   ! dummy value
!//          real(8)                    :: dummass  ! dummy value
!//          real(8)                    :: dumeng   ! dummy value
!//          real(8)                    :: dumceng  ! dummy value
!//          real(8)                    :: dumdheng ! dummy value
!// 
!//          ! block #2 variables
!//          character(5)               :: arrow ! type of reaction
!//          character(10)              :: prod_name ! product name
!//          character(10)              :: react_name ! reaction name
!//          integer                    :: tmpdecaynumb ! number of decay reactions
!//          integer                    :: tmpntrnnumb ! number of neutron-induced reactions
!//          real(8)                    :: tmpdecconst ! decay constant of nuclide
!//          character(7)               :: dumunit ! dummy /second
!//          real(8)                    :: tmpprob ! reaction probability
!// 
!//          !block #3 variables
!//          real(8),allocatable        :: tmpbranchdata(:)
!//          logical                    :: switch_br = .false. ! find 2nd part of the table
!//          integer                    :: point1 = 1
!// 
!//          ! general variables
!//          character(3),allocatable   :: buff_char3(:)
!//          character(100)             :: oneline
!//          character(100)             :: tmpline
!//          character(1)               :: buff_char
!//          character(18)              :: name_of_file
!//          integer                    :: inp_f = 15
!//          integer                    :: i,j,k
!//          integer                    :: cnt1,cnt2,cnt3,cnt4,cnt5 ! counts
!//          logical                    :: end_of_file = .false.
!//          logical                    :: eob1 = .false.
!//          logical                    :: eob2 = .false.
!//          logical                    :: eob3 = .false.
!// 
!// 
!// 
!// #ifdef siarhei_plot
!//          ! @$^ siarhei_fr
!//          current_sub = &
!//             '+#+ Entered [fr_read_depletion_file] in FastDepletion'
!//          if(trim(last_saved).ne.trim(current_sub)) then
!//             write(678,'(A)') trim(current_sub)
!//             last_saved = current_sub
!//          end if
!//          ! =-=-=-=-=-=-=-=
!// #endif
!// 
!// 109      format (7X,I4,7X,I4,5X,I4,7X,I4,6X,I4,7X,I4,4X,I4,7X,I4)
!// 110      format (1X,A10,I9,F10.3,I10,F7.3,I10,F10.3,E10.3)
!// 120      format (1X,A10,I5,I5,ES15.6,3X,A7)
!// 121      format (A5,A10,A10,ES16.8)
!// 130      format (1X,A9,5I12)
!// 131      format (1X,A9,5ES12.5)
!// 
!// 
!//           write(*,*) "Start initial setting for depletion"
!//          open(unit=inp_f,file=name_of_file,status='old',   &
!//             access='sequential',form='formatted',action='read')
!// !            access='sequential',form='formatted',action='read')
!// ! write(*,*) 'aaa1'
!//          cnt1 = 1
!//          cnt2 = 0
!//          cnt5 = 0
!// !         fission_flag = 0D0
!// 
!//          if (allocated(all_names)) deallocate(all_names)
!// 
!//          ! find important headline data
!//          do while (.not. headline)
!//             read(inp_f,'(a)') oneline
!//             if (oneline(5:11).ne.'numNuk') then
!//                cycle
!//             else
!//                headline = .true.
!//                read(inp_f,'(a)') oneline
!//                read(oneline,109) numNuk,numHM,numFP,maxDPerN,  &
!//                   numDpath,maxNPerN,numNpath,numFPYset
!// !               write(*,*) "Headline recorded"
!// !               write(*,*) numNuk,numHM,numFP,maxDPerN,  &
!// !                  numDpath,maxNPerN,numNpath,numFPYset
!//                write(*,*) "# Nuclide | # Heavy Nuclide  | # Fission Product "
!//                write(*,*) numNuk,numHM,numFP
!// !               write(*,*) '|',numNuk,'|',numHM,'|',numFP,'|',(numNuk+345)
!//                if (.not.allocated(all_names)) allocate(all_names(1:numNuk))
!// !               if (.not.allocated(all_names%neutron)) allocate(all_names%neutron(1:numNuk))
!// !               if (.not.allocated(decay)) allocate(decay(1:numNuk))
!// ! need to remove when kappa fission is correct
!//                if (.not.allocated(kappaf)) allocate(kappaf(1:numNuk))
!//                exit
!//             end if
!//          end do
!// 
!//          ! dealing with block #1 data
!//          do while (.not.eob1.and..not.end_of_file)
!//             read(inp_f,'(a)') oneline
!//             ! find the beginning of a block
!//             if (oneline(3:5).eq."+++") then
!//                do i=1,(len(oneline)-8) ! find the block number
!//                   if (oneline(i:i+7)=='Block #1') then
!//                      block1 = .true.
!//                      block2 = .false.
!//                      block3 = .false.
!// !                     write(*,*) "Block 1 started"
!//                      exit
!//                   end if
!//                end do
!//             else if (oneline(3:5).eq."===") then ! end of block
!//                if (block1) then
!// !                  write(*,*) "Block 1 finished"
!//                   block1 = .false.
!//                   eob1 = .true.
!//                end if
!//             end if
!// 
!//             ! saving block #1 data
!//             if (block1.and.oneline(1:1).ne.'%') then
!//                read(oneline,110) tmpname,dumid,dummass, &
!//                   tmpfiflag,dumeng,dumflg,dumceng,dumdheng
!// !                  tmpfiflag,kappaf,dumflg,dumceng,dumdheng
!//                if (cnt1.le.numHM) then
!//                   all_names(cnt1) % name = tmpname
!//                   all_names(cnt1) % ID = dumid
!// 
!//                   all_names(cnt1) % fission_branch = tmpfiflag
!//                else
!//                   all_names(cnt1) % name = tmpname
!//               all_names(cnt1) % ID = dumid
!// 
!//                   all_names(cnt1) % fission_branch = -4444
!//                end if
!//                cnt1 = cnt1 + 1
!//             end if
!//          end do
!// 
!//          ! dealing with block #2 data
!//          do while (.not.eob2.and..not.end_of_file)
!//             read(inp_f,'(a)') oneline
!//             ! find the beginning of a block
!//             if (oneline(3:5).eq."+++") then
!//                do i=1,(len(oneline)-8)
!//                   if (oneline(i:i+7)=='Block #2') then
!//                      block1 = .false.
!//                      block2 = .true.
!//                      block3 = .false.
!// !                     write(*,*) "Block 2 started"
!//                      exit
!//                   end if
!//                end do
!//             else if (oneline(3:5).eq."===") then ! end of block
!//                if (block2) then
!// !                  write(*,*) "Block 2 finished"
!//                   block2 = .false.
!//                   eob2 = .true.
!//                end if
!//             end if
!// 
!//             ! saving block #2 data
!//             if (block2.and.oneline(1:1).ne.'%') then
!// 
!//                if (oneline(2:4).eq.'|=>') then ! neutron reactions
!//                   read(oneline,121) arrow,prod_name,react_name,tmpprob
!//                   all_names(cnt2) % neutron(cnt3) % product_name = prod_name
!//                   all_names(cnt2) % neutron(cnt3) % probability = tmpprob
!//                   all_names(cnt2) % neutron(cnt3) % what_reaction = react_name
!//                   do k = 1,numNuk
!//                      if (trim(all_names(k) % name).eq.prod_name) &
!//                         all_names(cnt2) % neutron(cnt3) % coord = k
!//                   end do
!//                   cnt3 = cnt3 + 1
!// 
!//                else if (oneline(2:4).eq.'|->') then ! decay reactions
!//                   read(oneline,121) arrow,prod_name,react_name,tmpprob
!//                   all_names(cnt2) % decay(cnt4) % product_name = prod_name
!//                   all_names(cnt2) % decay(cnt4) % probability = tmpprob
!//                   all_names(cnt2) % decay(cnt4) % what_reaction = react_name
!//                   do k = 1,numNuk
!//                      if (trim(all_names(k) % name).eq.prod_name) &
!//                         all_names(cnt2) % decay(cnt4) % coord = k
!//                   end do
!//                   cnt4 = cnt4 + 1
!// 
!//                else ! read parent nuclide data
!//                   read(oneline,120) tmpname,tmpdecaynumb,   &
!//                      tmpntrnnumb,tmpdecconst,dumunit
!//                   cnt2 = cnt2 + 1
!//                   if (all_names(cnt2) % name.eq.tmpname) then
!//                      all_names(cnt2) % len_neutron = tmpntrnnumb
!//                      all_names(cnt2) % len_decay = tmpdecaynumb
!//                      if (.not.allocated(all_names(cnt2) % neutron)) &
!//                         allocate(all_names(cnt2) % neutron(tmpntrnnumb))
!//                      if (.not.allocated(all_names(cnt2) % decay)) &
!//                         allocate(all_names(cnt2) % decay(tmpdecaynumb))
!//                      all_names(cnt2) % decay_const = tmpdecconst
!//                      cnt3 = 1
!//                      cnt4 = 1
!//                   else
!//                      write(*,*) "Error! Name doesn't match"
!//                      write(*,*) 'all_names|',all_names(cnt2) % name,'|'
!//                      write(*,*) 'tmpname|',tmpname,'|'
!//                   end if
!//                end if
!//             end if
!//          end do
!// 
!//          ! allocate fission yield probabilities for all FP
!//          if (eob2) then
!//             if (.not.allocated(tmpbranchdata)) &
!//                allocate(tmpbranchdata(1:maxNPerN))
!//             do i = 1,numNuk
!//                if (all_names(i) % fission_branch<0) then
!//                   if (.not.allocated(all_names(i) % fission_yield)) &
!//                      allocate(all_names(i) % fission_yield(1:numFPYset))
!//                end if
!//             end do
!//          end if
!// 
!//          ! dealing with block #3 data
!//          do while (.not.eob3.and..not.end_of_file)
!//             read(inp_f,'(a)') oneline
!//             ! find the beginning of a block
!//             if (oneline(3:5).eq."+++") then
!//                do i=1,(len(oneline)-8) ! determine the activation of a block
!//                   if (oneline(i:i+7)=='Block #3') then
!//                      block1 = .false.
!//                      block2 = .false.
!//                      block3 = .true.
!// !                     write(*,*) "Block 3 started"
!//                      exit
!//                   end if
!//                end do
!//             else if (oneline(3:5).eq."===") then ! end of block and file
!//                if (block3) then
!//                   close(inp_f)
!// !                  write(*,*) "Block 3 finished"
!//                   end_of_file = .true. ! end of file
!//                   block3 = .false.
!//                   eob2 = .true.
!// !                  write(*,*) "End of file!!! "
!// 
!//                end if
!//             end if
!// 
!//             if (block3.and.oneline(1:1).ne.'%') then
!//                switch_br = .true.
!//               ! write(*,*) point1,oneline
!//                read(oneline,131) tmpname,tmpbranchdata(1:5)
!// 
!//                if (all_names(numHM+cnt5) % name.eq.tmpname) then
!//                   all_names(numHM+cnt5) &
!//                      % fission_yield(point1:point1+4) = tmpbranchdata(1:5)
!//                   cnt5 = cnt5 + 1
!//                else
!//                   do i = (numHM+cnt5),numNuk
!//                      if (all_names(i) % name.eq.tmpname) then
!//                         cnt5 = i - numHM + 1
!//                         all_names(i) &
!//                            % fission_yield(point1:point1+4) = tmpbranchdata(1:5)
!//                         exit
!//                      end if
!//                   end do
!//                end if
!//             end if
!// 
!//             if (block3.and.switch_br.and.oneline(1:1).eq.'%') then
!//                point1 = 6
!//                cnt5 = 0
!//             end if
!// 
!//          end do
!//                   close(inp_f)
!// 
!// 
!// 
!//          if (allocated(tmpbranchdata)) deallocate(tmpbranchdata)
!// 
!//          write(*,*) "Finished reading depletion file"
!// 
!//       end subroutine fr_read_depletion_file
!// 
!//       ! set coeffs for CRAM (can be updated in future with new coeffs/order)
!//       subroutine fr_set_cram_coefficients
!//          use Inc_FastDepletion
!//          use Inc_Geometry, ONLY: Nxy, Nz
!//          use Inc_Option,   only: n_group
!//          use Inc_miXS,     only: miXS_Hex, miXS_n2n
!//          use Inc_File,     ONLY: NuNum
!//          use Inc_Depletion, ONLY: N_BU
!//          USE Inc_RP, ONLY: AxialComp, I_LP_1N
!// 
!// 
!//          implicit none
!//          integer                             :: totbu
!// 
!//          real(8)                             :: temp_read
!//          integer                             :: i,j
!//          integer                             :: Ixy,Iz
!// 
!// 
!// #ifdef siarhei_fr
!// 
!// #ifdef siarhei_plot
!//          ! @$^ siarhei_fr
!//          current_sub = &
!//             '+#+ Entered [fr_set_cram_coefficients] in FastDepletion'
!//          if(trim(last_saved).ne.trim(current_sub)) then
!//             write(678,'(A)') trim(current_sub)
!//             last_saved = current_sub
!//          end if
!//          ! =-=-=-=-=-=-=-=
!// #endif
!// 
!//          if (NuNum.ne.numNuk) numNuk = NuNum !@$^ Siarhei_FR
!// #endif
!// 
!//          totbu = N_BU
!//          cram_alpha_0 = 2.124853710495237488D-16
!//          order_cram = 16
!//          if (.not.allocated(all_concentrations)) &
!//             allocate(all_concentrations(1:Nxy, 1:Nz, 1:numNuk,0:totbu))
!//          if (allocated(A_matrix)) deallocate(A_matrix)
!//          if (allocated(cram_alpha)) deallocate(cram_alpha)
!//          if (allocated(cram_theta)) deallocate(cram_theta)
!//          if (allocated(vec_N_old)) deallocate(vec_N_old)
!// 
!//          if (.not.allocated(xs_a)) allocate(xs_a(Nxy,Nz,numNuk,N_group))
!//          if (.not.allocated(xs_f)) allocate(xs_f(Nxy,Nz,numNuk,N_group))
!//          if (.not.allocated(xs_n2n)) allocate(xs_n2n(Nxy,Nz,numNuk,N_group))
!//          if (.not.allocated(cram_alpha)) allocate(cram_alpha(order_cram/2))
!//          if (.not.allocated(cram_theta)) allocate(cram_theta(order_cram/2))
!//          if (.not.allocated(vec_N_old)) allocate(vec_N_old(numNuk))
!//          if (.not.allocated(mat_I_all)) allocate(mat_I_all(numNuk,numNuk))
!//          if (.not.allocated(sum_CRAM_vec)) allocate(sum_CRAM_vec(numNuk))
!// 
!//          cram_alpha = ( 0.0D0 , 0.0D0 )
!//          cram_theta = ( 0.0D0 , 0.0D0 )
!// 
!//          cram_alpha(1) = ( - 5.0901521865224915650D-7 , - 2.4220017652852287970D-5 )
!//          cram_alpha(2) = ( + 2.1151742182466030907D-4 , + 4.3892969647380673918D-3 )
!//          cram_alpha(3) = ( + 1.1339775178483930527D+2 , + 1.0194721704215856450D+2 )
!//          cram_alpha(4) = ( + 1.5059585270023467528D+1 , - 5.7514052776421819979D+0 )
!//          cram_alpha(5) = ( - 6.4500878025539646595D+1 , - 2.2459440762652096056D+2 )
!//          cram_alpha(6) = ( - 1.4793007113557999718D+0 , + 1.7686588323782937906D+0 )
!//          cram_alpha(7) = ( - 6.2518392463207918892D+1 , - 1.1190391094283228480D+1 )
!//          cram_alpha(8) = ( + 4.1023136835410021273D-2 , - 1.5743466173455468191D-1 )
!// 
!//          cram_theta(1) = ( - 1.0843917078696988026D+1 , + 1.9277446167181652284D+1 )
!//          cram_theta(2) = ( - 5.2649713434426468895D+0 , + 1.6220221473167927305D+1 )
!//          cram_theta(3) = ( + 5.9481522689511774808D+0 , + 3.5874573620183222829D+0 )
!//          cram_theta(4) = ( + 3.5091036084149180974D+0 , + 8.4361989858843750826D+0 )
!//          cram_theta(5) = ( + 6.4161776990994341923D+0 , + 1.1941223933701386874D+0 )
!//          cram_theta(6) = ( + 1.4193758971856659786D+0 , + 1.0925363484496722585D+1 )
!//          cram_theta(7) = ( + 4.9931747377179963991D+0 , + 5.9968817136039422260D+0 )
!//          cram_theta(8) = ( - 1.4139284624888862114D+0 , + 1.3497725698892745389D+1 )
!//          write(*,*) "Should be dimensions: ", numNuk,totbu
!//          if (.not.allocated(A_matrix)) allocate(A_matrix(1:numNuk,1:numNuk))
!//          A_matrix = 0d0
!//          mat_I_all = (0d0,0d0)
!// 
!//          do i = 1, numNuk
!//             mat_I_all(i,i) = (1d0,0d0)
!// 
!//          end do
!//             xs_f(:,:,:,:) = miXS_Hex(:,:,:,3,:) *1d-24 ! convert from barn to cm^2
!// 
!//             DO Ixy = 1,Nxy
!//                 DO Iz = 1,Nz
!//                 END DO
!//             END DO
!// 
!//             xs_a(:,:,:,:) = miXS_Hex(:,:,:,2,:)  *1d-24
!//       end subroutine fr_set_cram_coefficients
!// 
!//       ! fill out the A matrix
!//       subroutine fr_update_xs_data(a,b)
!//          use Inc_FastDepletion
!//          use Inc_TPEN, only: fluxf
!//          implicit none
!//          integer                             :: a,b
!//          integer                             :: i,j
!//          integer                             :: brch
!//          integer                             :: cd
!//          real(8)                             :: lamb
!// 
!// 
!// #ifdef siarhei_plot
!//          ! @$^ siarhei_fr
!//          current_sub = &
!//             '+#+ Entered [fr_update_xs_data] in FastDepletion'
!//          if(trim(last_saved).ne.trim(current_sub)) then
!//             write(678,'(A)') trim(current_sub)
!//             last_saved = current_sub
!//          end if
!//          ! =-=-=-=-=-=-=-=
!// #endif
!//           
!//         write(19999,*) 'Ixy, Iz:  ',a,b
!//         write(19999,*) 'abs XS 1: ',xs_a(a,b,1,:)
!//         write(19999,*) 'abs XS U234: ',xs_a(a,b,6,:)
!//         write(19999,*) 'abs XS U235: ',xs_a(a,b,7,:)
!//         write(19999,*) 'abs XS U238: ',xs_a(a,b,10,:)
!//         write(19999,*) 'fission XS 1: ',xs_f(a,b,1,:)
!//         write(19999,*) 'fission XS U234: ',xs_f(a,b,6,:)
!//         write(19999,*) 'fission XS U235: ',xs_f(a,b,7,:)
!//         write(19999,*) 'fission XS U238: ',xs_f(a,b,10,:)
!//         write(19999,*) 'fluxf: ',fluxf(a,b,:)
!//         write(19999,*) '   '
!// 
!// 
!//          A_matrix = 0d0
!//          do i = 1, numNuk
!// 
!//     lamb = all_names(i) % decay_const
!// 
!//     A_matrix(i,i) = A_matrix(i,i) - &
!//        lamb - dot_product(xs_a(a,b,i,:),fluxf(a,b,:))
!//  ! deal with neutron reactions
!//     do j = 1, all_names(i) % len_neutron
!//        if (trim(all_names(i) % neutron(j) % what_reaction).eq.&
!//           'CAPTURE') then
!//           cd = all_names(i) % neutron(j) % coord
!// 
!// 
!//           A_matrix(cd,i) = A_matrix(cd,i) + dot_product(xs_a(a,b,i,:),fluxf(a,b,:))*&
!//              all_names(i) % neutron(j) % probability  - dot_product(xs_f(a,b,i,:),fluxf(a,b,:))*&
!//              all_names(i) % neutron(j) % probability  !-n2n_U35
!// 
!// 
!//        else if (trim(all_names(i) % neutron(j) % what_reaction).eq.&
!//           'N2N') then
!// !           cd = all_names(i) % neutron(j) % coord
!// !            cd = all_names(i) % neutron(j) % coord
!// !            A_matrix(i,i) = A_matrix(i,i) - dot_product(xs_n2n(a,b,i,:),fluxf(a,b,:))*&
!// !               all_names(i) % neutron(j) % probability
!// !
!// !            A_matrix(cd,i) = A_matrix(cd,i) + dot_product(xs_n2n(a,b,i,:),fluxf(a,b,:))*&
!// !               all_names(i) % neutron(j) % probability    !  write(*,*) "N2n"
!//           ! A_matrix(i,i) = A_matrix(i,i) - dot_product(xs_n2n(a,b,i,:),flux(a,b,:))*&
!//           !    all_names(i) % neutron(j) % probability
!// !           if (i==2) then
!// !           A_matrix(cd,i) = A_matrix(cd,i) + n2n_U35*&
!// !              all_names(i) % neutron(j) % probability
!// !          else if (i==5) then
!// !           A_matrix(cd,i) = A_matrix(cd,i) + n2n_U38*&
!// !              all_names(i) % neutron(j) % probability         !  write(*,*) "N2n"
!// !          else
!// !           A_matrix(cd,i) = A_matrix(cd,i) + n2n_Pu39*&
!// !              all_names(i) % neutron(j) % probability         !  write(*,*) "N2n"
!// !
!// !          end if
!// 
!//        else
!//        end if
!//     end do
!// 
!// 
!// 
!// ! deal with decay reactions
!// 
!//     do j = 1, all_names(i) % len_decay
!//        cd = all_names(i) % decay(j) % coord
!// 
!//        A_matrix(cd,i) = A_matrix(cd,i) + &
!//           all_names(i) % decay(j) % probability * &
!//           lamb
!// 
!//     end do
!// 
!//             ! deal with fission
!//             if (i.le.numHM) then
!//                brch = all_names(i) % fission_branch
!//                do j = numHM+1,numNuk
!//                   A_matrix(j,i) = A_matrix(j,i) + dot_product(xs_f(a,b,i,:),fluxf(a,b,:))*&
!//                      all_names(j) % fission_yield(brch)
!//                end do
!//             end if
!//          end do
!// ! ////////////////////////////////////////////
!// 
!// !#endif
!// 
!//       end subroutine fr_update_xs_data
!// 
!// 
!// 
!// 
!// 
!//       ! function for finding time step in EFPD
!//       real(8) function fr_find_depletion_time_step(bu1,bu2,core_pow)
!//          USE Inc_Depletion , ONLY: Tot_MTU
!// 
!//         ! use Inc_3D, only: Core_Power
!// 
!//          implicit none
!//          real(8)                    :: core_pow ! GWth
!//          real(8)                    :: bu1,bu2
!// 
!//          fr_find_depletion_time_step = abs(bu2-bu1)*&
!//             Tot_MTU/core_pow
!// 
!//       end function fr_find_depletion_time_step
!// 
!// 
!//       ! aux subroutine for writing the cram matrix
!//       ! to a new file (due to the size)
!//       subroutine fr_print_cram_matrix
!//          use Inc_FastDepletion
!//          implicit none
!//          integer           :: out_f=16
!//          integer           :: i,j
!//          character(3000)   :: temp_line
!//          character(15)     :: temp_elem
!//          character(3)      :: temp_index
!// 
!// 
!// #ifdef siarhei_plot
!//          ! @$^ siarhei_fr
!//          current_sub = &
!//             '+#+ Entered [fr_print_cram_matrix] in FastDepletion'
!//          if(trim(last_saved).ne.trim(current_sub)) then
!//             write(678,'(A)') trim(current_sub)
!//             last_saved = current_sub
!//          end if
!//          ! =-=-=-=-=-=-=-=
!// #endif
!// 
!//          write(*,*) "Start writing to a file - CRAM_matrix.out"
!//          open(unit=out_f,file='CRAM_matrix.out',status='replace',action='write')
!// 
!//         ! write(*,*) "Num Nuk = ", numNuk
!//          temp_line = ''
!//          do j = 1,numNuk
!//             temp_index = ''
!// 
!//             if ((j/10).ge.10) then
!//                write(temp_index,'(I3)') j
!//                temp_line = trim(adjustl(temp_line))//'         '//trim(adjustl(temp_index))
!//             else if ((j/10).ge.1) then
!//                write(temp_index,'(I2)') j
!//                temp_line = trim(adjustl(temp_line))//'          '//trim(adjustl(temp_index))
!//             else
!//                write(temp_index,'(I1)') j
!//                temp_line = trim(adjustl(temp_line))//'           '//trim(adjustl(temp_index))
!//             end if
!//            ! write(*,*) temp_line
!//          end do
!//          write(out_f,'(A)') '            '//trim(adjustl(temp_line))
!//          write(*,*) "CRAM head written", numNuk
!// 
!//          temp_line = ''
!//          do i = 1,numNuk
!//                                   write(*,*) "CRAM head written1 "
!// 
!//             if (len(trim(all_names(i) % name)).le.5) then
!//                temp_line = trim(all_names(i) % name)//'_'
!//                                   write(*,*) "CRAM head written2 "
!// 
!//             else
!//                temp_line = trim(all_names(i) % name)
!//                                   write(*,*) "CRAM head written3 "
!// 
!//             end if
!//                      write(*,*) "CRAM head written1 "
!// 
!//             do j = 1,numNuk
!//                temp_elem = ''
!//                if (A_matrix(i,j)<0) then
!//                   write(temp_elem,'(ES12.4)') A_matrix(i,j)
!//                else
!//                   write(temp_elem,'(ES12.5)') A_matrix(i,j)
!//                end if
!//                temp_line = trim(adjustl(temp_line))//trim(temp_elem)
!//             end do
!//                      write(*,*) "CRAM head written2"
!// 
!//             write(out_f,'(A)') temp_line
!//            ! write(*,'(A)') temp_line
!//          end do
!//          close(out_f)
!//          write(*,*) "Finished writing to a file - CRAM_matrix.out"
!// 
!//       end subroutine fr_print_cram_matrix
!// !
!//       subroutine fr_print_concentrations(Ixy,Iz)
!//          use Inc_FastDepletion
!//          use Inc_miXS, only: N_FR
!//         
!//          use Inc_Depletion, ONLY: N_BU,Cycle_BU
!//          implicit none
!//          integer           :: out_f=16
!//          integer           :: i,j
!//          character(500)    :: temp_line
!//          character(20)     :: temp_elem
!//          integer           :: totbu
!//          integer           :: Ixy, Iz
!// 
!//          totbu = N_BU
!// 
!// !         open(unit=out_f,file='nuclide_concentrations.out', &
!// !            status='old')
!// !            status='replace',action='write')
!// !         write(*,*) Ixy, Iz
!// 
!//          temp_line = ''
!//          do j = 1,totbu
!//             temp_elem = ''
!//             write(temp_elem,'(F5.2)') Cycle_BU(j)
!// !            write(*,'(F5.2)') Cycle_BU(j)
!//             if ((Cycle_BU(j)/10).ge.1) then
!//                temp_line = &
!//                trim(adjustl(temp_line))//'         '//trim(adjustl(temp_elem))
!//             else
!//                temp_line = &
!//                trim(adjustl(temp_line))//'          '//trim(adjustl(temp_elem))
!//             end if
!//          end do
!//          write(1994,'(A)') 'Nuclide      '//trim(adjustl(temp_line))
!// !         write(*,'(A)') '             '//trim(adjustl(temp_line))
!// 
!//          do i = 1,numNuk
!//             temp_line = ''
!//             if (len(trim(all_names(i) % name)).le.5) then
!//                temp_line = trim(all_names(i) % name)//'_'
!//             else
!//                temp_line = trim(all_names(i) % name)
!//             end if
!//             do j = 1, totbu
!//                temp_elem = ''
!//                write(temp_elem,'(ES12.5)') N_FR(Ixy,Iz,i)/1D-24
!//                temp_line = trim(adjustl(temp_line))//'   '//trim(adjustl(temp_elem))
!//             end do
!//             write(1994,'(A)') trim(adjustl(temp_line))
!// !            write(*,'(A)') trim(adjustl(temp_line))
!//             temp_line = ''
!//          end do
!// !         close(out_f)
!// 
!//       end subroutine fr_print_concentrations
!// 
!// 
!//       ! function Get_InvC copied from Mod_Operator and slightly modified
!// !      function fr_get_invMatx(C) result(invC)
!// !         use Inc_FastDepletion
!// !      implicit none
!// !
!// !      complex(8), dimension(:, :), intent(IN) :: C
!// !      complex(8), dimension(:, :), allocatable :: invC, X
!// !      complex(8) :: Quotient
!// !      integer :: N
!// !      integer :: i, j, m
!// !
!// !
!// !#ifdef siarhei_plot
!// !         ! @$^ siarhei_fr
!// !         current_sub = &
!// !            '+#+ Entered [fr_get_invMatx] in FastDepletion'
!// !         if(trim(last_saved).ne.trim(current_sub)) then
!// !            write(678,'(A)') trim(current_sub)
!// !            last_saved = current_sub
!// !         end if
!// !         ! =-=-=-=-=-=-=-=
!// !#endif
!// !
!// !      i    = 0
!// !      j    = 0
!// !      m    = 0
!// !
!// !      N = size(C,2) ! should be 221 (numNuk)
!// !
!// !      allocate(invC(N, N))
!// !      allocate(X(N, 2*N))  ! X = [C | eye]
!// !      invC = (0.0d0, 0.0d0)
!// !      X    = (0.0d0, 0.0d0)
!// !
!// !      do i = 1, N ! copy C to X (to the left hand side)
!// !         do j = 1, N
!// !            X(i, j) = C(i, j)
!// !         end do
!// !      end do
!// !
!// !      do i = 1, N
!// !         X(i, (i + N)) = (1.0d0, 0.0d0) ! fill in the 'eye' part of the matrix
!// !      end do                            ! right hand side part
!// !
!// !      do i = 1, N
!// !         Quotient = X(i, i) ! diagonal from the C matrix (should be non-zero!!!)
!// !         do j = i, 2*N ! divide the right upper side of the matrix by diag elem
!// !            X(i, j) = X(i, j)/Quotient ! < - here is why should be non-zero
!// !         end do
!// !         if (i == N) exit
!// !         do m = (i + 1), N
!// !            Quotient = X(m, i) ! manipulate with the 'eye' matrix
!// !            do j = i, 2*N
!// !               X(m, j) = X(m, j) - Quotient*X(i, j) ! X was already modified by previous
!// !            end do                  ! Quotient, we return it back, by than adjusting the
!// !         end do                     ! 'eye' part
!// !      end do
!// !
!// !      ! This Procedure Makes Matrix X_Left => Reduced Echelon Form,
!// !      ! Then X_Right(Identity Matrix) => inv(C)
!// !      do i = N, 2, - 1
!// !         do m = (i - 1), 1, - 1
!// !            Quotient = X(m, i)
!// !            do j = 1, 2*N
!// !               X(m, j) = X(m, j) - Quotient*X(i, j)
!// !            end do
!// !         end do
!// !      end do
!// !      do i = 1, N
!// !         do j = 1, N
!// !            invC(i, j) = X(i, (j + N))
!// !         end do
!// !      end do
!// !
!// !      deallocate(X)
!// !
!// !      return
!// !
!// !      end function fr_get_invMatx
!// 
!// 
!// 
!//       subroutine fr_CRAM_solver
!//          use Inc_FastDepletion
!//          use Inc_Time, only: dT
!//          use Inc_RP
!//          use Inc_Geometry, only: izFuelTop, izFuelBot
!// !        ! USE Mod_Operator, ONLY: Get_InvC
!//          use Inc_File,     ONLY: NuNum
!// ! remove in future
!//          USE Inc_XS_File, ONLY: I_BU
!// !
!// !         USE Read_File, ONLY:  update_maXS
!// 
!//          use Inc_miXS, ONLY: N_FR, N_FR_old
!// 
!//          implicit none
!//          integer                             :: Ixy_FA,Ixy, Iz
!//          integer                             :: i,j
!// 
!//          INTEGER(4),ALLOCATABLE  :: SPS_n_non0(:)                        !                             @@@ DO NOT END
!//          INTEGER(4),ALLOCATABLE  :: SPS_is2J(:)                          !                             @@@ DO NOT END
!//          INTEGER(4)              :: SPS_SUM_non0=0                 !                             @@@ DO NOT END
!//          INTEGER(4),PARAMETER    :: n_CRAM_k=8
!// 
!//      ! local
!//          INTEGER(4)              :: ii
!//          INTEGER(4)              :: jj
!//          INTEGER(4)              :: is
!//          INTEGER(4)              :: k
!//          INTEGER(4)              :: tmp_non0_idx
!//          REAL(8)                 :: diag_matA(NuNum)
!//          COMPLEX(8)              :: div(NuNum)
!//      !    REAL(8)                 :: vec_N_old(NuNum)
!// 
!//          COMPLEX(8)              :: complex_x(NuNum)
!//          COMPLEX(8)              :: complex_b(NuNum)
!//          REAL(8)                 :: complex_xsum(NuNum)
!//          REAL(8), ALLOCATABLE    :: SPS_matA(:)
!//     REAL(8)                                     :: alpha_real(0:n_CRAM_k)
!//     REAL(8)                                     :: alpha_img(1:n_CRAM_k)
!//     REAL(8)                                     :: theta_real(1:n_CRAM_k)
!//     REAL(8)                                     :: theta_img(1:n_CRAM_k)
!//     data alpha_real(0:n_CRAM_k)   &
!//     / +2.1248537104952237488e-16 &
!//     , -5.0901521865224915650e-7  &
!//     , +2.1151742182466030907e-4  &
!//     , +1.1339775178483930527e+2  &
!//     , +1.5059585270023467528e+1  &
!//     , -6.4500878025539646595e+1  &
!//     , -1.4793007113557999718e+0  &
!//     , -6.2518392463207918892e+1  &
!//     , +4.1023136835410021273e-2  /
!//     data alpha_img(1:n_CRAM_k)    &
!//     / -2.4220017652852287970e-5  &
!//     , +4.3892969647380673918e-3  &
!//     , +1.0194721704215856450e+2  &
!//     , -5.7514052776421819979e+0  &
!//     , -2.2459440762652096056e+2  &
!//     , +1.7686588323782937906e+0  &
!//     , -1.1190391094283228480e+1  &
!//     , -1.5743466173455468191e-1  /
!//     data theta_real(1:n_CRAM_k)   &
!//     / -1.0843917078696988026e+1  &
!//     , -5.2649713434426468895e+0  &
!//     , +5.9481522689511774808e+0  &
!//     , +3.5091036084149180974e+0  &
!//     , +6.4161776990994341923e+0  &
!//     , +1.4193758971856659786e+0  &
!//     , +4.9931747377179963991e+0  &
!//     , -1.4139284624888862114e+0  /
!//     data theta_img(1:n_CRAM_k)    &
!//     / +1.9277446167181652284e+1  &
!//     , +1.6220221473167927305e+1  &
!//     , +3.5874573620183222829e+0  &
!//     , +8.4361989858843750826e+0  &
!//     , +1.1941223933701386874e+0  &
!//     , +1.0925363484496722585e+1  &
!//     , +5.9968817136039422260e+0  &
!//     , +1.3497725698892745389e+1  /
!// 
!// 
!// #ifdef siarhei_plot
!//          ! @$^ siarhei_fr
!//          current_sub = &
!//             '+#+ Entered [fr_CRAM_solver] in FastDepletion'
!//          if(trim(last_saved).ne.trim(current_sub)) then
!//             write(678,'(A)') trim(current_sub)
!//             last_saved = current_sub
!//          end if
!//          ! =-=-=-=-=-=-=-=
!// #endif
!//       
!//       write(20000,*) 'I_BU: ', I_BU
!// 
!//       DO Iz = IzFuelBot, IzFuelTop
!// #ifdef siarhei_plot
!//             dummy_filler = 1 ! @$^ siarhei_plot 
!// #endif 
!//           DO Ixy_FA = 1, Nxy_FA
!//               Ixy = I_FA(Ixy_FA)
!// 
!//                  do jj=1,NuNum
!//                      vec_N_old(jj)=N_FR_old(Ixy,Iz,jj)*1d24  ! #/(barn-cm)
!// write(20000,*) 'N_FR_old(Ixy,Iz,jj) : ', Ixy,IZ, jj, N_FR_old(Ixy,Iz,jj)
!//                  enddo
!//                  ! Set full matrix A
!//                  call fr_update_xs_data(Ixy,Iz)
!//              !    NuNum = NuNum
!//                  if(allocated(SPS_n_non0)) deallocate(SPS_n_non0)
!// 
!//                  ALLOCATE(SPS_n_non0(NuNum))
!//                  SPS_n_non0(:) = 0
!//                  do ii=1,NuNum
!//                      tmp_non0_idx = 0
!//                      do jj=1,NuNum
!//                          if(ii == jj) cycle
!//                          if(abs(A_matrix(ii,jj))>0.0E0) then
!//                              tmp_non0_idx = tmp_non0_idx + 1
!//                          endif
!//                      enddo
!//                      SPS_n_non0(ii) = tmp_non0_idx
!//                  enddo
!// 
!//                  SPS_SUM_non0=sum(SPS_n_non0(:))
!//                  if(allocated(SPS_is2J)) deallocate(SPS_is2J)
!//                  allocate(SPS_is2J(SPS_SUM_non0))
!//                  SPS_is2J = 0
!//                  is=0
!//                  do ii=1,NuNum
!//                      !tmp_non0_idx = 0
!//                      do jj=1, NuNum
!//                          if(ii == jj) cycle
!//                          if(abs(A_matrix(ii,jj))>0) then
!//                              !tmp_non0_idx = tmp_non0_idx + 1
!//                              is=is+1
!//                              SPS_is2J(is) = jj
!//                          endif
!//                      enddo
!//                  enddo
!// 
!// 
!//                  if(allocated(SPS_matA)) deallocate(SPS_matA)
!//                  allocate(SPS_matA(SPS_sum_non0))
!// 
!//                     ! Make matrix At
!//                  SPS_matA(:)=0d0
!//                  is=0
!//                  do ii=1,NuNum
!//                      do jj=1,SPS_n_non0(ii)
!//                          is=is+1
!//                          SPS_matA(is)=A_matrix(ii,SPS_is2j(is))*dT
!// 
!//                      enddo
!//                      diag_matA(ii)=A_matrix(ii,ii)*dT
!//                  enddo
!// write(20000,*) 'dT:  ',dT
!// 
!//                  ! Solve n=n0*exp(At) with CRAM
!//                  complex_xsum(:)=0d0
!//                  complex_x(:)   =1d0
!//                  do k=1,n_CRAM_k
!//                      do ii=1,NuNum
!//                          div(ii)=1d0/(diag_matA(ii)-theta_real(k)-theta_img(k)*(0.0E0,1.0E0))
!//                          complex_b(ii)=(alpha_real(k)+alpha_img(k)*(0.0E0,1.0E0))*vec_N_old(ii)
!//                      enddo
!// 
!//                      call NESTED_SGS_AXB
!// 
!//                      do ii=1,NuNum
!//                          complex_xsum(ii)=complex_xsum(ii)+real(complex_x(ii))
!//                      enddo
!//                  enddo
!//                  
!//                  write(20000,*) 'Ixy: ', Ixy, '  Iz: ', Iz
!// 
!//                  do ii=1,NuNum
!//                      vec_N_old(ii)=vec_N_old(ii)*alpha_real(0)+2d0*complex_xsum(ii)
!//                      vec_N_old(ii)=max(0.0E0,vec_N_old(ii))
!//                      if(vec_N_old(ii)<1d-30) vec_N_old(ii)=0.0E0
!//                      N_FR(Ixy,Iz,ii)=vec_N_old(ii)*1D-24
!// !                     all_concentrations(Ixy, Iz,ii,I_BU) = vec_N_old(ii)*1d-24
!//                  write(20000,*) N_FR(Ixy,Iz,ii)
!// 
!//                  enddo
!// 
!// 
!//              END DO
!//          END DO
!// 
!// 
!//     return
!//     contains
!// 
!// 
!//     SUBROUTINE NESTED_SGS_AXB
!//     implicit none
!//     INTEGER(4)                                  :: i
!//     INTEGER(4)                                  :: m
!//     INTEGER(4)                                  :: is
!//     INTEGER(4)                                  :: iter
!//     INTEGER(4),PARAMETER                        :: max_iter=1000
!//     REAL(8)                                     :: max_err
!//     COMPLEX(8)                                  :: Xo(NuNum)
!//     COMPLEX(8)                                  :: tmp1
!// 
!// 
!// #ifdef siarhei_plot
!//          ! @$^ siarhei_fr
!//          current_sub = &
!//             '+#+ Entered [NESTED_SGS_AXB] in FastDepletion'
!//          if(trim(last_saved).ne.trim(current_sub)) then
!//             write(678,'(A)') trim(current_sub)
!//             last_saved = current_sub
!//          end if
!//          ! =-=-=-=-=-=-=-=
!// #endif
!// 
!//     do iter=1,max_iter
!//         Xo(:)=complex_x(:)
!//         is=0
!//         do i=1,NuNum
!//             tmp1=0.0
!//             do m=1,SPS_n_non0(i)
!//                 is  =is+1
!//                 tmp1=tmp1+sps_matA(is)*complex_x( SPS_is2j(is) ) ! do such as GS
!//             enddo
!//             complex_x(i)=(complex_b(i)-tmp1)*div(i)
!// #ifdef RFORM
!//             !if(abs(real(complex_x(i)))<1.0E-35_r8k) complex_x(i)=1.0E-60_r8k
!//             if(abs((complex_x(i)))<1.0d-35) complex_x(i)=1.0d-60
!// #else
!//             !if(abs(real(complex_x(i)))<1.0E-100_r8k) complex_x(i)=1.0E-99_r8k
!//             if(abs((complex_x(i)))<1.0d-100) complex_x(i)=1.0d-99
!// #endif
!//         enddo
!//         max_err=0d0
!//         do i=1,NuNum
!// #ifdef RFORM
!//             !max_Err=max(max_err,abs((complex_x(i)-xo(i))/(1.0E-30_r8k+complex_x(i))))
!//             max_Err=max(max_err,abs((complex_x(i)-xo(i))/max(1.0d-25,abs(xo(i)),abs(complex_x(i)))))
!// #else
!//             max_Err=max(max_err,abs((complex_x(i)-xo(i))/(1.0d-70+complex_x(i))))
!// #endif
!//         enddo
!// 
!//         if(max_err<1.0d-06 .and. iter>=3) then
!//             exit
!//         elseif(iter==max_iter) then
!//             write(*,*) 'SGS_AXB solver diverges in SGS_AXB routine. max_err=',max_err
!//             do i=1,NuNum
!// #ifdef RFORM
!//                 max_Err=max(max_err,abs((complex_x(i)-xo(i))/max(1.0d-25,abs(xo(i)),abs(complex_x(i)))))
!// #else
!//                 max_Err=max(max_err,abs((complex_x(i)-xo(i))/(1.0d-70+complex_x(i))))
!// #endif
!//                 write(*,*) i,complex_x(i),xo(i),max_err,abs((complex_x(i)-xo(i))/max(1.0d-25,abs(xo(i)),abs(complex_x(i))))
!//             enddo
!//             !call end_stream
!//         endif
!//     enddo
!//     return
!//     END SUBROUTINE NESTED_SGS_AXB
!// 
!//       end subroutine fr_CRAM_solver
!// 
!//       ! deallocation of previously allocated memory for FastDepletion
!//       subroutine fr_deallocate_arrays
!//          use Inc_FastDepletion
!//          implicit none
!// 
!// 
!// #ifdef siarhei_plot
!//          ! @$^ siarhei_fr
!//          current_sub = &
!//             '+#+ Entered [fr_deallocate_arrays] in FastDepletion'
!//          if(trim(last_saved).ne.trim(current_sub)) then
!//             write(678,'(A)') trim(current_sub)
!//             last_saved = current_sub
!//          end if
!//          ! =-=-=-=-=-=-=-=
!// #endif
!// 
!//          if (allocated(all_names)) deallocate(all_names)
!// !         if (allocated(all_concentrations)) deallocate(all_concentrations)
!//          if (allocated(A_matrix)) deallocate(A_matrix)
!//          if (allocated(cram_alpha)) deallocate(cram_alpha)
!//          if (allocated(cram_theta)) deallocate(cram_theta)
!//          if (allocated(vec_N_old)) deallocate(vec_N_old)
!//          if (allocated(mat_I_all)) deallocate(mat_I_all)
!//          if (allocated(xs_a)) deallocate(xs_a)
!//          if (allocated(xs_f)) deallocate(xs_f)
!//          if (allocated(xs_n2n)) deallocate(xs_n2n)
!//          if (allocated(sum_CRAM_vec)) deallocate(sum_CRAM_vec)
!// 
!//       end subroutine fr_deallocate_arrays
!// 
!// 
!//       ! will be replaced with some RAST-K data
!//       subroutine fr_load_initial_fuel(a, b, fuel_opt)
!//          use Inc_FastDepletion
!//          use Inc_miXS, ONLY: N_FR
!//          USE Inc_RP, ONLY: AxialComp, I_LP_1N
!//          implicit none
!//          ! fuel_opt:
!//          ! 1 - UO2 fuel initial load
!//          ! 2 - MOX fuel initial load
!//          ! 3 - ThU fuel initial load
!//          ! 4 - ThPu fuel initial load
!//          integer                    :: fuel_opt
!//          integer                    :: coolant_opt
!//          integer                    :: clad_opt
!//          integer                    :: a,b
!//          integer                    :: i,j
!//          real(8)                    :: multiplier = 3.201D21
!// 
!//          all_concentrations(a,b,:,0) = 0.0
!// !         write(*,*) shape(all_concentrations),numNuk
!//          select case(fuel_opt)
!// 
!//          case (1) ! U
!//              do i = 1,numNuk
!//                if (trim(all_names(i) % name).eq.'U-234') &
!//                   all_concentrations(a,b,i,0) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'U-235') &
!//                   all_concentrations(a,b,i,0) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'U-238') &                 
!//                    all_concentrations(a,b,i,0) = N_FR(a,b,i)
!// !               if (trim(all_names(i) % name).eq.'ZR095') &
!// !                  all_concentrations(a,b,i,1) = 1.00000E-06
!// 
!//              end do
!//          case (2) ! MOX
!//             do i = 1,numNuk
!//                if (trim(all_names(i) % name).eq.'U-235') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)             
!//                if (trim(all_names(i) % name).eq.'U-238') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)       
!//                if (trim(all_names(i) % name).eq.'PU238') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'PU239') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'PU240') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'PU241') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'PU242') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'AM241') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//             end do   
!// 
!//          case (3) ! ThU
!//             do i = 1,numNuk
!//                if (trim(all_names(i) % name).eq.'U-235') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'U-238') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'U-234') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'TH232') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//             end do
!// 
!//          case (4) ! ThPu
!//             do i = 1,numNuk
!//                if (trim(all_names(i) % name).eq.'TH232') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'PU238') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'PU239') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'PU240') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'PU241') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'PU242') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//                if (trim(all_names(i) % name).eq.'AM241') &
!//                   all_concentrations(a,b,i,1) = N_FR(a,b,i)
!//             end do
!//             
!//          case default ! just a copy of the UO2 branch
!//             do i = 1,numNuk
!//                if (trim(all_names(i) % name).eq.'U-235') &
!//                   all_concentrations(a,b,i,1) = 2.715D20
!//                if (trim(all_names(i) % name).eq.'U-238') &
!//                   all_concentrations(a,b,i,1) = 6.201D21
!//                if (trim(all_names(i) % name).eq.'U-234') &
!//                   all_concentrations(a,b,i,1) = 1.206D18
!//             end do
!// 
!//          end select
!// 
!// 
!//       end subroutine fr_load_initial_fuel
!// 
!//    end module FastDepletion
!// #endif
   ! module for Depletion in Fast Reactor
   !
   !
#ifdef tuan_fr
   module FastDepletion

     ! use Module1
     ! use Module2
     ! use Module3

     ! use Inc_Constant
     ! use Inc_Geometry
     ! use Inc_RP

      implicit none
   contains
      
      ! fully read the depletion file and save data 
      ! to 'all_names' array (Inc_FastDepletion) 
      subroutine fr_read_depletion_file(name_of_file)
         use Inc_FastDepletion

         implicit none
         
         ! block #1 variables
         character(1)               :: spce ! first space col
         character(10)              :: tmpname ! nuclide name         
         integer                    :: tmpfiflag ! fis branch
         integer                    :: dumid    ! dummy value
         integer                    :: dumflg   ! dummy value
         real(8)                    :: dummass  ! dummy value
         real(8)                    :: dumeng   ! dummy value
         real(8)                    :: dumceng  ! dummy value
         real(8)                    :: dumdheng ! dummy value
         
         ! block #2 variables
         character(5)               :: arrow ! type of reaction  
         character(10)              :: prod_name ! product name        
         character(10)              :: react_name ! reaction name
         integer                    :: tmpdecaynumb ! number of decay reactions
         integer                    :: tmpntrnnumb ! number of neutron-induced reactions
         real(8)                    :: tmpdecconst ! decay constant of nuclide
         character(7)               :: dumunit ! dummy /second
         real(8)                    :: tmpprob ! reaction probability

         !block #3 variables
         real(8),allocatable        :: tmpbranchdata(:)
         logical                    :: switch_br = .false. ! find 2nd part of the table
         integer                    :: point1 = 1

         ! general variables
         character(3),allocatable   :: buff_char3(:)
         character(100)             :: oneline
         character(100)             :: tmpline
         character(1)               :: buff_char
         character(50)              :: name_of_file
         integer                    :: inp_f = 15
         integer                    :: i,j,k
         integer                    :: cnt1,cnt2,cnt3,cnt4,cnt5 ! counts
         logical                    :: end_of_file = .false.
         logical                    :: eob1 = .false.
         logical                    :: eob2 = .false.
         logical                    :: eob3 = .false.


109      format (7X,I4,7X,I4,5X,I4,7X,I4,6X,I4,7X,I4,4X,I4,7X,I4)
110      format (1X,A10,I9,F10.3,I10,F7.3,I10,F10.3,E10.3)
120      format (1X,A10,I5,I5,ES15.6,3X,A7)
121      format (A5,A10,A10,ES16.8)
130      format (1X,A9,5I12)
131      format (1X,A9,5ES12.5)


          write(*,*) "Start initial setting for depletion"
         open(unit=inp_f,file=name_of_file,status='old',   &
            access='sequential',form='formatted',action='read')

         cnt1 = 1  
         cnt2 = 0    
         cnt5 = 0
!         fission_flag = 0D0

         if (allocated(all_names)) deallocate(all_names)

         ! find important headline data
         do while (.not.headline)
            read(inp_f,'(a)') oneline                        
            if (oneline(5:11).ne.'numNuk') then
               cycle
            else 
               headline = .true.
               read(inp_f,'(a)') oneline
               read(oneline,109) numNuk,numHM,numFP,maxDPerN,  &
                  numDpath,maxNPerN,numNpath,numFPYset
 !              write(*,*) "Headline recorded"
 !              write(*,*) numNuk,numHM,numFP,maxDPerN,  &
 !                 numDpath,maxNPerN,numNpath,numFPYset
 !              write(*,*) '|',numNuk,'|',numHM,'|',numFP,'|',(numNuk+345)
               if (.not.allocated(all_names)) allocate(all_names(1:numNuk))                  
! need to remove when kappa fission is correct
               if (.not.allocated(kappaf)) allocate(kappaf(1:numNuk))                  
               exit
            end if 
         end do                 
                           
         ! dealing with block #1 data              
         do while (.not.eob1.and..not.end_of_file)   
            read(inp_f,'(a)') oneline
            ! find the beginning of a block
            if (oneline(3:5).eq."+++") then
               do i=1,(len(oneline)-8) ! find the block number
                  if (oneline(i:i+7)=='Block #1') then 
                     block1 = .true.
                     block2 = .false.
                     block3 = .false.
!                     write(*,*) "Block 1 started"
                     exit                  
                  end if 
               end do
            else if (oneline(3:5).eq."===") then ! end of block
               if (block1) then 
                  write(*,*) "Block 1 finished"
                  block1 = .false.
                  eob1 = .true.
               end if
            end if

            ! saving block #1 data
            if (block1.and.oneline(1:1).ne.'%') then
               read(oneline,110) tmpname,dumid,dummass, &
                  tmpfiflag,dumeng,dumflg,dumceng,dumdheng            
!                  tmpfiflag,kappaf,dumflg,dumceng,dumdheng            
               if (cnt1.le.numHM) then               
                  all_names(cnt1) % name = tmpname     
                  all_names(cnt1) % ID = dumid   
				  
                  all_names(cnt1) % fission_branch = tmpfiflag                   
               else                 
                  all_names(cnt1) % name = tmpname
				  all_names(cnt1) % ID = dumid

                  all_names(cnt1) % fission_branch = -4444              
               end if
               cnt1 = cnt1 + 1            
            end if
         end do

         ! dealing with block #2 data
         do while (.not.eob2.and..not.end_of_file)   
            read(inp_f,'(a)') oneline
            ! find the beginning of a block
            if (oneline(3:5).eq."+++") then
               do i=1,(len(oneline)-8) 
                  if (oneline(i:i+7)=='Block #2') then
                     block1 = .false.
                     block2 = .true.
                     block3 = .false.
!                     write(*,*) "Block 2 started"
                     exit
                  end if
               end do
            else if (oneline(3:5).eq."===") then ! end of block
               if (block2) then 
!                  write(*,*) "Block 2 finished"
                  block2 = .false.
                  eob2 = .true.
               end if
            end if

            ! saving block #2 data
            if (block2.and.oneline(1:1).ne.'%') then

               if (oneline(2:4).eq.'|=>') then ! neutron reactions
                  read(oneline,121) arrow,prod_name,react_name,tmpprob                  
                  all_names(cnt2) % neutron(cnt3) % product_name = prod_name
                  all_names(cnt2) % neutron(cnt3) % probability = tmpprob    
                  all_names(cnt2) % neutron(cnt3) % what_reaction = react_name
                  do k = 1,numNuk
                     if (trim(all_names(k) % name).eq.prod_name) &
                        all_names(cnt2) % neutron(cnt3) % coord = k   
                  end do     
                  cnt3 = cnt3 + 1

               else if (oneline(2:4).eq.'|->') then ! decay reactions
                  read(oneline,121) arrow,prod_name,react_name,tmpprob                  
                  all_names(cnt2) % decay(cnt4) % product_name = prod_name
                  all_names(cnt2) % decay(cnt4) % probability = tmpprob   
                  all_names(cnt2) % decay(cnt4) % what_reaction = react_name
                  do k = 1,numNuk
                     if (trim(all_names(k) % name).eq.prod_name) &
                        all_names(cnt2) % decay(cnt4) % coord = k   
                  end do            
                  cnt4 = cnt4 + 1
                  
               else ! read parent nuclide data
                  read(oneline,120) tmpname,tmpdecaynumb,   &
                     tmpntrnnumb,tmpdecconst,dumunit                
                  cnt2 = cnt2 + 1
                  if (all_names(cnt2) % name.eq.tmpname) then
                     all_names(cnt2) % len_neutron = tmpntrnnumb
                     all_names(cnt2) % len_decay = tmpdecaynumb
                     if (.not.allocated(all_names(cnt2) % neutron)) & 
                        allocate(all_names(cnt2) % neutron(tmpntrnnumb))               
                     if (.not.allocated(all_names(cnt2) % decay)) & 
                        allocate(all_names(cnt2) % decay(tmpdecaynumb))                    
                     all_names(cnt2) % decay_const = tmpdecconst                
                     cnt3 = 1
                     cnt4 = 1                     
                  else
!                     write(*,*) "Error! Name doesn't match"
!                     write(*,*) 'all_names|',all_names(cnt2) % name,'|'
!                     write(*,*) 'tmpname|',tmpname,'|'
                  end if
               end if
            end if
         end do

         ! allocate fission yield probabilities for all FP
         if (eob2) then
            if (.not.allocated(tmpbranchdata)) &
               allocate(tmpbranchdata(1:maxNPerN))            
            do i = 1,numNuk
               if (all_names(i) % fission_branch<0) then
                  if (.not.allocated(all_names(i) % fission_yield)) &
                     allocate(all_names(i) % fission_yield(1:numFPYset)) 
               end if
            end do
         end if

         ! dealing with block #3 data
         do while (.not.eob3.and..not.end_of_file)   
            read(inp_f,'(a)') oneline
            ! find the beginning of a block
            if (oneline(3:5).eq."+++") then
               do i=1,(len(oneline)-8) ! determine the activation of a block
                  if (oneline(i:i+7)=='Block #3') then
                     block1 = .false.
                     block2 = .false.
                     block3 = .true.
!                     write(*,*) "Block 3 started"
                     exit
                  end if
               end do
            else if (oneline(3:5).eq."===") then ! end of block and file
               if (block3) then 
                  close(inp_f)
!                  write(*,*) "Block 3 finished"
                  end_of_file = .true. ! end of file
                  block3 = .false.
                  eob2 = .true.
!                  write(*,*) "End of file!!! "

               end if
            end if
     
            if (block3.and.oneline(1:1).ne.'%') then
               switch_br = .true.
              ! write(*,*) point1,oneline
               read(oneline,131) tmpname,tmpbranchdata(1:5)
               
               if (all_names(numHM+cnt5) % name.eq.tmpname) then
                  all_names(numHM+cnt5) &
                     % fission_yield(point1:point1+4) = tmpbranchdata(1:5)
                  cnt5 = cnt5 + 1
               else
                  do i = (numHM+cnt5),numNuk
                     if (all_names(i) % name.eq.tmpname) then
                        cnt5 = i - numHM + 1
                        all_names(i) &
                           % fission_yield(point1:point1+4) = tmpbranchdata(1:5)
                        exit
                     end if
                  end do
               end if
            end if

            if (block3.and.switch_br.and.oneline(1:1).eq.'%') then
               point1 = 6
               cnt5 = 0
            end if

         end do      
         


         if (allocated(tmpbranchdata)) deallocate(tmpbranchdata)      

         write(*,*) "Finished reading depletion file"
         
      end subroutine fr_read_depletion_file

      ! set coeffs for CRAM (can be updated in future with new coeffs/order)
      subroutine fr_set_cram_coefficients
         use Inc_FastDepletion
         use Inc_Geometry, ONLY: Nxy, Nz
         use Inc_Option,   only: n_group
         use Inc_miXS,     only: miXS_Hex, miXS_n2n
         use Inc_File,     ONLY: NuNum
         use Inc_Depletion, ONLY: N_BU
         USE Inc_RP, ONLY: AxialComp, I_LP_1N


         implicit none
         integer                             :: totbu
         
!         integer                             :: xs_f_   = 11
!         integer                             :: xs_n2n_ = 12
!         integer                             :: xs_a_   = 13
         real(8)                             :: temp_read
         integer                             :: i,j
         integer                             :: Ixy,Iz
         



!         totbu = size(burnups)
         totbu = N_BU
         cram_alpha_0 = 2.124853710495237488D-16
         order_cram = 16

         if (allocated(all_concentrations)) deallocate(all_concentrations)
         if (allocated(A_matrix)) deallocate(A_matrix)
         if (allocated(cram_alpha)) deallocate(cram_alpha)
         if (allocated(cram_theta)) deallocate(cram_theta)
         if (allocated(vec_N_old)) deallocate(vec_N_old)
   
         if (.not.allocated(xs_a)) allocate(xs_a(Nxy,Nz,numNuk,N_group))
         if (.not.allocated(xs_f)) allocate(xs_f(Nxy,Nz,numNuk,N_group))
         if (.not.allocated(xs_n2n)) allocate(xs_n2n(Nxy,Nz,numNuk,N_group))
         if (.not.allocated(cram_alpha)) allocate(cram_alpha(order_cram/2))
         if (.not.allocated(cram_theta)) allocate(cram_theta(order_cram/2)) 
         if (.not.allocated(vec_N_old)) allocate(vec_N_old(numNuk))
         if (.not.allocated(mat_I_all)) allocate(mat_I_all(numNuk,numNuk))
         if (.not.allocated(sum_CRAM_vec)) allocate(sum_CRAM_vec(numNuk))

         cram_alpha = ( 0.0D0 , 0.0D0 )
         cram_theta = ( 0.0D0 , 0.0D0 )
            
         cram_alpha(1) = ( - 5.0901521865224915650D-7 , - 2.4220017652852287970D-5 )
         cram_alpha(2) = ( + 2.1151742182466030907D-4 , + 4.3892969647380673918D-3 )
         cram_alpha(3) = ( + 1.1339775178483930527D+2 , + 1.0194721704215856450D+2 )
         cram_alpha(4) = ( + 1.5059585270023467528D+1 , - 5.7514052776421819979D+0 )
         cram_alpha(5) = ( - 6.4500878025539646595D+1 , - 2.2459440762652096056D+2 )
         cram_alpha(6) = ( - 1.4793007113557999718D+0 , + 1.7686588323782937906D+0 )
         cram_alpha(7) = ( - 6.2518392463207918892D+1 , - 1.1190391094283228480D+1 )
         cram_alpha(8) = ( + 4.1023136835410021273D-2 , - 1.5743466173455468191D-1 )
                                                                              
         cram_theta(1) = ( - 1.0843917078696988026D+1 , + 1.9277446167181652284D+1 )
         cram_theta(2) = ( - 5.2649713434426468895D+0 , + 1.6220221473167927305D+1 )
         cram_theta(3) = ( + 5.9481522689511774808D+0 , + 3.5874573620183222829D+0 )
         cram_theta(4) = ( + 3.5091036084149180974D+0 , + 8.4361989858843750826D+0 )
         cram_theta(5) = ( + 6.4161776990994341923D+0 , + 1.1941223933701386874D+0 )
         cram_theta(6) = ( + 1.4193758971856659786D+0 , + 1.0925363484496722585D+1 )
         cram_theta(7) = ( + 4.9931747377179963991D+0 , + 5.9968817136039422260D+0 )
         cram_theta(8) = ( - 1.4139284624888862114D+0 , + 1.3497725698892745389D+1 ) 
         write(*,*) "Should be dimensions: ", numNuk,totbu
         if (.not.allocated(all_concentrations)) &
            allocate(all_concentrations(1:Nxy, 1:Nz, 1:numNuk,0:totbu))
         if (.not.allocated(A_matrix)) allocate(A_matrix(1:numNuk,1:numNuk))         
!         all_concentrations = 0d0
         A_matrix = 0d0
         mat_I_all = (0d0,0d0)

         do i = 1, numNuk
            mat_I_all(i,i) = (1d0,0d0)

         end do
            xs_f(:,:,:,:) = miXS_Hex(:,:,:,3,:) *1d-24 ! convert from barn to cm^2

            DO Ixy = 1,Nxy
                DO Iz = 1,Nz
! for 221 nu         
!                    xs_n2n(Ixy,Iz,:,:) = miXS_n2n(AxialComp(I_LP_1N( Ixy ), Iz),1, :,:) *1d-24
                END DO
            END DO

            xs_a(:,:,:,:) = miXS_Hex(:,:,:,2,:)  *1d-24
      end subroutine fr_set_cram_coefficients

      ! fill out the A matrix
      subroutine fr_update_xs_data(a,b)
         use Inc_FastDepletion
         use Inc_TPEN, only: fluxf
        ! use ifport ! for rand()
         implicit none
         integer                             :: a,b
         integer                             :: i,j
         integer                             :: brch
         integer                             :: cd
!         real(8),intent(in)                  :: dt
         real(8)                             :: lamb

         
         A_matrix = 0d0
!         flux(a,b,:)=1d15
         write(1994,*) 'flux(a,b,:):', a,b
         write(1994,'(24ES13.6)') fluxf(a,b,:)
!         write(*,'(24ES13.6)') fluxf(a,b,:)
         do i = 1, numNuk

    lamb = all_names(i) % decay_const  

    A_matrix(i,i) = A_matrix(i,i) - & 
       lamb - dot_product(xs_a(a,b,i,:),fluxf(a,b,:))
 ! deal with neutron reactions
    do j = 1, all_names(i) % len_neutron 
       if (trim(all_names(i) % neutron(j) % what_reaction).eq.&
          'CAPTURE') then
          cd = all_names(i) % neutron(j) % coord

!           A_matrix(i,i) = A_matrix(i,i) - ab(i)*&
!              all_names(i) % neutron(j) % probability

          A_matrix(cd,i) = A_matrix(cd,i) + dot_product(xs_a(a,b,i,:),fluxf(a,b,:))*&
             all_names(i) % neutron(j) % probability  - dot_product(xs_f(a,b,i,:),fluxf(a,b,:))*&
             all_names(i) % neutron(j) % probability  !-n2n_U35  
          
      
       else if (trim(all_names(i) % neutron(j) % what_reaction).eq.&
          'N2N') then
!           cd = all_names(i) % neutron(j) % coord
!            cd = all_names(i) % neutron(j) % coord
!            A_matrix(i,i) = A_matrix(i,i) - dot_product(xs_n2n(a,b,i,:),fluxf(a,b,:))*&
!               all_names(i) % neutron(j) % probability 
! 
!            A_matrix(cd,i) = A_matrix(cd,i) + dot_product(xs_n2n(a,b,i,:),fluxf(a,b,:))*&
!               all_names(i) % neutron(j) % probability    !  write(*,*) "N2n"
          ! A_matrix(i,i) = A_matrix(i,i) - dot_product(xs_n2n(a,b,i,:),flux(a,b,:))*&
          !    all_names(i) % neutron(j) % probability 
!           if (i==2) then
!           A_matrix(cd,i) = A_matrix(cd,i) + n2n_U35*&
!              all_names(i) % neutron(j) % probability 
!          else if (i==5) then
!           A_matrix(cd,i) = A_matrix(cd,i) + n2n_U38*&
!              all_names(i) % neutron(j) % probability         !  write(*,*) "N2n"
!          else
!           A_matrix(cd,i) = A_matrix(cd,i) + n2n_Pu39*&
!              all_names(i) % neutron(j) % probability         !  write(*,*) "N2n"
!              
!          end if
          
       else
        !  write(*,*) "Reaction not recognized"
       end if
    end do



! deal with decay reactions   

    do j = 1, all_names(i) % len_decay                
       cd = all_names(i) % decay(j) % coord               

       A_matrix(cd,i) = A_matrix(cd,i) + & 
          all_names(i) % decay(j) % probability * &  
          lamb

    end do

            ! deal with fission 
            if (i.le.numHM) then
               brch = all_names(i) % fission_branch
               do j = numHM+1,numNuk
                  A_matrix(j,i) = A_matrix(j,i) + dot_product(xs_f(a,b,i,:),fluxf(a,b,:))*& 
                     all_names(j) % fission_yield(brch)
               end do
            end if      
         end do
! ////////////////////////////////////////////         


             
      end subroutine fr_update_xs_data

      



      ! function for finding time step in EFPD
      real(8) function fr_find_depletion_time_step(bu1,bu2,core_pow) 
         USE Inc_Depletion , ONLY: Tot_MTU

        ! use Inc_3D, only: Core_Power

         implicit none
         real(8)                    :: core_pow ! GWth
         real(8)                    :: bu1,bu2

         fr_find_depletion_time_step = abs(bu2-bu1)*&
            Tot_MTU/core_pow

      end function fr_find_depletion_time_step

      ! will be replaced with some RAST-K data
      subroutine fr_load_initial_fuel(a, b, fuel_opt)
         use Inc_FastDepletion
         use Inc_miXS, ONLY: N_FR
         USE Inc_RP, ONLY: AxialComp, I_LP_1N
         implicit none
         ! fuel_opt:
         ! 1 - UO2 fuel initial load
         ! 2 - MOX fuel initial load
         ! 3 - ThU fuel initial load
         ! 4 - ThPu fuel initial load
         integer                    :: fuel_opt
         integer                    :: coolant_opt
         integer                    :: clad_opt
         integer                    :: a,b
         integer                    :: i,j
         real(8)                    :: multiplier = 3.201D21

         all_concentrations(a,b,:,0) = 0.0
!         write(*,*) shape(all_concentrations),numNuk
         select case(fuel_opt)

         case (1) ! U
             do i = 1,numNuk
               if (trim(all_names(i) % name).eq.'U-234') &
                  all_concentrations(a,b,i,0) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'U-235') &
                  all_concentrations(a,b,i,0) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'U-238') &                 
                   all_concentrations(a,b,i,0) = N_FR(a,b,i)
!               if (trim(all_names(i) % name).eq.'ZR095') &
!                  all_concentrations(a,b,i,1) = 1.00000E-06

             end do
         case (2) ! MOX
            do i = 1,numNuk
               if (trim(all_names(i) % name).eq.'U-235') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)             
               if (trim(all_names(i) % name).eq.'U-238') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)       
               if (trim(all_names(i) % name).eq.'PU238') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'PU239') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'PU240') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'PU241') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'PU242') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'AM241') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
            end do   

         case (3) ! ThU
            do i = 1,numNuk
               if (trim(all_names(i) % name).eq.'U-235') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'U-238') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'U-234') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'TH232') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
            end do

         case (4) ! ThPu
            do i = 1,numNuk
               if (trim(all_names(i) % name).eq.'TH232') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'PU238') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'PU239') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'PU240') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'PU241') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'PU242') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
               if (trim(all_names(i) % name).eq.'AM241') &
                  all_concentrations(a,b,i,1) = N_FR(a,b,i)
            end do
            
         case default ! just a copy of the UO2 branch
            do i = 1,numNuk
               if (trim(all_names(i) % name).eq.'U-235') &
                  all_concentrations(a,b,i,1) = 2.715D20
               if (trim(all_names(i) % name).eq.'U-238') &
                  all_concentrations(a,b,i,1) = 6.201D21
               if (trim(all_names(i) % name).eq.'U-234') &
                  all_concentrations(a,b,i,1) = 1.206D18
            end do

         end select


      end subroutine fr_load_initial_fuel


      ! aux subroutine for writing the cram matrix 
      ! to a new file (due to the size)
      subroutine fr_print_cram_matrix
         use Inc_FastDepletion
         implicit none
         integer           :: out_f=16
         integer           :: i,j
         character(3000)   :: temp_line
         character(15)     :: temp_elem
         character(3)      :: temp_index

         write(*,*) "Start writing to a file - CRAM_matrix.out"
         open(unit=out_f,file='CRAM_matrix.out',status='replace',action='write')

        ! write(*,*) "Num Nuk = ", numNuk
         temp_line = ''
         do j = 1,numNuk
            temp_index = ''
            
            if ((j/10).ge.10) then 
               write(temp_index,'(I3)') j
               temp_line = trim(adjustl(temp_line))//'         '//trim(adjustl(temp_index))
            else if ((j/10).ge.1) then 
               write(temp_index,'(I2)') j
               temp_line = trim(adjustl(temp_line))//'          '//trim(adjustl(temp_index))
            else 
               write(temp_index,'(I1)') j
               temp_line = trim(adjustl(temp_line))//'           '//trim(adjustl(temp_index))
            end if
           ! write(*,*) temp_line
         end do
         write(out_f,'(A)') '            '//trim(adjustl(temp_line))
         write(*,*) "CRAM head written", numNuk

         temp_line = ''         
         do i = 1,numNuk
                                  write(*,*) "CRAM head written1 "

            if (len(trim(all_names(i) % name)).le.5) then 
               temp_line = trim(all_names(i) % name)//'_'
                                  write(*,*) "CRAM head written2 "

            else
               temp_line = trim(all_names(i) % name)
                                  write(*,*) "CRAM head written3 "
               
            end if
                     write(*,*) "CRAM head written1 "

            do j = 1,numNuk
               temp_elem = ''
               if (A_matrix(i,j)<0) then
                  write(temp_elem,'(ES12.4)') A_matrix(i,j)
               else
                  write(temp_elem,'(ES12.5)') A_matrix(i,j)
               end if
               temp_line = trim(adjustl(temp_line))//trim(temp_elem)
            end do
                     write(*,*) "CRAM head written2"

            write(out_f,'(A)') temp_line   
           ! write(*,'(A)') temp_line      
         end do
         close(out_f)
         write(*,*) "Finished writing to a file - CRAM_matrix.out"    
    
      end subroutine fr_print_cram_matrix
!
      subroutine fr_print_concentrations(Ixy,Iz)
         use Inc_FastDepletion
         use Inc_Depletion, ONLY: N_BU,Cycle_BU
         implicit none
         integer           :: out_f=16
         integer           :: i,j
         character(500)    :: temp_line
         character(20)     :: temp_elem
         integer           :: totbu
         integer           :: Ixy, Iz

         totbu = N_BU

!         open(unit=out_f,file='nuclide_concentrations.out', &
!            status='old')
!            status='replace',action='write')
!         write(*,*) Ixy, Iz

         temp_line = ''
         do j = 1,totbu
            temp_elem = ''
            write(temp_elem,'(F5.2)') Cycle_BU(j)
!            write(*,'(F5.2)') Cycle_BU(j)
            if ((Cycle_BU(j)/10).ge.1) then 
               temp_line = &            
               trim(adjustl(temp_line))//'         '//trim(adjustl(temp_elem))
            else 
               temp_line = &            
               trim(adjustl(temp_line))//'          '//trim(adjustl(temp_elem))
            end if
         end do
         write(1994,'(A)') 'Nuclide      '//trim(adjustl(temp_line))
!         write(*,'(A)') '             '//trim(adjustl(temp_line))

         do i = 1,numNuk
            temp_line = ''        
            if (len(trim(all_names(i) % name)).le.5) then 
               temp_line = trim(all_names(i) % name)//'_'
            else
               temp_line = trim(all_names(i) % name)
            end if
            do j = 1, totbu
               temp_elem = ''
               write(temp_elem,'(ES12.5)') all_concentrations(Ixy,Iz,i,j)/1D-24
               temp_line = trim(adjustl(temp_line))//'   '//trim(adjustl(temp_elem))
            end do
            write(1994,'(A)') trim(adjustl(temp_line))
!            write(*,'(A)') trim(adjustl(temp_line))
            temp_line = ''
         end do
!         close(out_f)

      end subroutine fr_print_concentrations


      ! function Get_InvC copied from Mod_Operator and slightly modified
      function fr_get_invMatx(C) result(invC)
         use Inc_FastDepletion
      implicit none

      complex(8), dimension(:, :), intent(IN) :: C              
      complex(8), dimension(:, :), allocatable :: invC, X 
      complex(8) :: Quotient 
      integer :: N             
      integer :: i, j, m
      
      i    = 0
      j    = 0
      m    = 0
      
      N = size(C,2) ! should be 221 (numNuk)
     
      allocate(invC(N, N))
      allocate(X(N, 2*N))  ! X = [C | eye]
      invC = (0.0d0, 0.0d0)
      X    = (0.0d0, 0.0d0)               
      
      do i = 1, N ! copy C to X (to the left hand side)
         do j = 1, N
            X(i, j) = C(i, j) 
         end do
      end do

      do i = 1, N 
         X(i, (i + N)) = (1.0d0, 0.0d0) ! fill in the 'eye' part of the matrix 
      end do                            ! right hand side part

      do i = 1, N
         Quotient = X(i, i) ! diagonal from the C matrix (should be non-zero!!!)
         do j = i, 2*N ! divide the right upper side of the matrix by diag elem
            X(i, j) = X(i, j)/Quotient ! < - here is why should be non-zero
         end do
         if (i == N) exit   
         do m = (i + 1), N 
            Quotient = X(m, i) ! manipulate with the 'eye' matrix
            do j = i, 2*N
               X(m, j) = X(m, j) - Quotient*X(i, j) ! X was already modified by previous 
            end do                  ! Quotient, we return it back, by than adjusting the 
         end do                     ! 'eye' part
      end do
      
      ! This Procedure Makes Matrix X_Left => Reduced Echelon Form,
      ! Then X_Right(Identity Matrix) => inv(C)
      do i = N, 2, - 1         
         do m = (i - 1), 1, - 1
            Quotient = X(m, i)
            do j = 1, 2*N
               X(m, j) = X(m, j) - Quotient*X(i, j)
            end do
         end do
      end do
      do i = 1, N
         do j = 1, N
            invC(i, j) = X(i, (j + N))
         end do
      end do
      
      deallocate(X)
         
      return

      end function fr_get_invMatx         


      ! CRAM matrix solving subroutine 
!!    !  Heavy nuclide+fission product
!      subroutine fr_CRAM_solver
!         use Inc_FastDepletion
!         use Inc_Time, only: dT
!         use Inc_RP
!         use Inc_Geometry, only: izFuelTop, izFuelBot
!         USE Mod_Operator, ONLY: Get_InvC
!! remove in future
!         USE Inc_XS_File, ONLY: I_BU
!!
!         USE Read_File, ONLY:  update_maXS
!
!         use Inc_miXS, ONLY: N_FR, N_FR_old
!		 
!         implicit none
!         integer                             :: Ixy_FA,Ixy, Iz
!         integer                             :: i,j
!
!! need to replace in future
!!                 chosen_time_step = &
!!                fr_find_depletion_time_step(Cycle_BU(I_BU-1),Cycle_BU(I_BU),Core_Power*1D-9)
!!                dt=chosen_time_step*24*3600
!         DO Ixy_FA = 1, Nxy_FA
!             Ixy = I_FA(Ixy_FA)
!             DO Iz = IzFuelBot, IzFuelTop
!
!             vec_N_old(:) = N_FR_old(Ixy,Iz,:)*1d24
!                 call fr_update_xs_data(Ixy,Iz)
!!                 sum_CRAM_vec = (0d0,0d0)
!!                 do i = 1,(order_cram/2)
!!                    sum_CRAM_vec(:) = sum_CRAM_vec(:) + cram_alpha(i)*&
!!                       matmul(Get_invC(dT*A_matrix - cram_theta(i)*&
!!                       mat_I_all),vec_N_old)
!!                 end do
!!                 vec_N_old = cram_alpha_0*vec_N_old+2d0*real(sum_CRAM_vec,8)
!!                 do i = 1,numNuk
!!                    if (vec_N_old(i)<1d-30) vec_N_old(i) = 0d0
!!                 end do
!
!!write(1994,*) vec_N_old
!
!
!                     all_concentrations(Ixy, Iz,:,I_BU) = vec_N_old(:)*1d-24
!                     N_FR( Ixy , Iz,:) = vec_N_old(:) *1d-24
!
!!                     CALL miXS_func(Ixy,Iz)
!!                      call update_maXS(Ixy,Iz)
!             END DO
!         END DO
!
!      end subroutine fr_CRAM_solver

      subroutine fr_CRAM_solver
         use Inc_FastDepletion
         use Inc_Time, only: dT
         use Inc_RP
         use Inc_Geometry, only: izFuelTop, izFuelBot
!         USE Mod_Operator, ONLY: Get_InvC
         use Inc_File,     ONLY: NuNum
! remove in future
         USE Inc_XS_File, ONLY: I_BU
!
!         USE Read_File, ONLY:  update_maXS

         use Inc_miXS, ONLY: N_FR, N_FR_old

         implicit none
         integer                             :: Ixy_FA,Ixy, Iz
         integer                             :: i,j

         INTEGER(4),ALLOCATABLE  :: SPS_n_non0(:)                        !                             @@@ DO NOT END
         INTEGER(4),ALLOCATABLE  :: SPS_is2J(:)                          !                             @@@ DO NOT END
         INTEGER(4)              :: SPS_SUM_non0=0                 !                             @@@ DO NOT END
         INTEGER(4),PARAMETER    :: n_CRAM_k=8
     
     ! local
         INTEGER(4)              :: ii
         INTEGER(4)              :: jj
         INTEGER(4)              :: is
         INTEGER(4)              :: k
         INTEGER(4)              :: tmp_non0_idx
         REAL(8)                 :: diag_matA(NuNum)
         COMPLEX(8)              :: div(NuNum)
     !    REAL(8)                 :: vec_N_old(NuNum)
     
         COMPLEX(8)              :: complex_x(NuNum)
         COMPLEX(8)              :: complex_b(NuNum)
         REAL(8)                 :: complex_xsum(NuNum)
         REAL(8), ALLOCATABLE    :: SPS_matA(:)
    REAL(8)                                     :: alpha_real(0:n_CRAM_k)
    REAL(8)                                     :: alpha_img(1:n_CRAM_k)
    REAL(8)                                     :: theta_real(1:n_CRAM_k)
    REAL(8)                                     :: theta_img(1:n_CRAM_k)
    data alpha_real(0:n_CRAM_k)   &
    / +2.1248537104952237488e-16 & 
    , -5.0901521865224915650e-7  & 
    , +2.1151742182466030907e-4  & 
    , +1.1339775178483930527e+2  & 
    , +1.5059585270023467528e+1  & 
    , -6.4500878025539646595e+1  & 
    , -1.4793007113557999718e+0  & 
    , -6.2518392463207918892e+1  & 
    , +4.1023136835410021273e-2  / 
    data alpha_img(1:n_CRAM_k)    &
    / -2.4220017652852287970e-5  &   
    , +4.3892969647380673918e-3  &
    , +1.0194721704215856450e+2  &
    , -5.7514052776421819979e+0  &
    , -2.2459440762652096056e+2  &
    , +1.7686588323782937906e+0  &
    , -1.1190391094283228480e+1  &
    , -1.5743466173455468191e-1  /
    data theta_real(1:n_CRAM_k)   &
    / -1.0843917078696988026e+1  &
    , -5.2649713434426468895e+0  &
    , +5.9481522689511774808e+0  &
    , +3.5091036084149180974e+0  &
    , +6.4161776990994341923e+0  &
    , +1.4193758971856659786e+0  &
    , +4.9931747377179963991e+0  &
    , -1.4139284624888862114e+0  /
    data theta_img(1:n_CRAM_k)    &
    / +1.9277446167181652284e+1  &
    , +1.6220221473167927305e+1  &
    , +3.5874573620183222829e+0  &
    , +8.4361989858843750826e+0  &
    , +1.1941223933701386874e+0  &
    , +1.0925363484496722585e+1  &
    , +5.9968817136039422260e+0  &
    , +1.3497725698892745389e+1  /
     
      DO Ixy_FA = 1, Nxy_FA
             Ixy = I_FA(Ixy_FA)
             DO Iz = IzFuelBot, IzFuelTop

                 do jj=1,NuNum
                     vec_N_old(jj)=N_FR_old(Ixy,Iz,jj)*1d24  ! #/(barn-cm)
                 enddo

                 ! Set full matrix A
                 call fr_update_xs_data(Ixy,Iz)
             !    NuNum = NuNum
                 if(allocated(SPS_n_non0)) deallocate(SPS_n_non0)
                 allocate(SPS_n_non0(NuNum))
                 SPS_n_non0(:) = 0
                 do ii=1,NuNum
                     tmp_non0_idx = 0
                     do jj=1,NuNum
                         if(ii == jj) cycle
                         if(abs(A_matrix(ii,jj))>0.0E0) then
                             tmp_non0_idx = tmp_non0_idx + 1
                         endif
                     enddo
                     SPS_n_non0(ii) = tmp_non0_idx
                 enddo

                 SPS_SUM_non0=sum(SPS_n_non0(:))
                 if(allocated(SPS_is2J)) deallocate(SPS_is2J)
                 allocate(SPS_is2J(SPS_SUM_non0))
                 SPS_is2J = 0
                 is=0
                 do ii=1,NuNum
                     !tmp_non0_idx = 0
                     do jj=1, NuNum
                         if(ii == jj) cycle
                         if(abs(A_matrix(ii,jj))>0) then
                             !tmp_non0_idx = tmp_non0_idx + 1
                             is=is+1
                             SPS_is2J(is) = jj
                         endif
                     enddo
                 enddo


                 if(allocated(SPS_matA)) deallocate(SPS_matA)
                 allocate(SPS_matA(SPS_sum_non0))

                    ! Make matrix At
                 SPS_matA(:)=0d0
                 is=0
                 do ii=1,NuNum
                     do jj=1,SPS_n_non0(ii)
                         is=is+1
                         SPS_matA(is)=A_matrix(ii,SPS_is2j(is))*dT

                     enddo
                     diag_matA(ii)=A_matrix(ii,ii)*dT
                 enddo
Write(1994,*)  'dT:   ', dT

                 ! Solve n=n0*exp(At) with CRAM
                 complex_xsum(:)=0d0
                 complex_x(:)   =1d0
                 do k=1,n_CRAM_k
                     do ii=1,NuNum
                         div(ii)=1d0/(diag_matA(ii)-theta_real(k)-theta_img(k)*(0.0E0,1.0E0))
                         complex_b(ii)=(alpha_real(k)+alpha_img(k)*(0.0E0,1.0E0))*vec_N_old(ii)
                     enddo

                     call NESTED_SGS_AXB

                     do ii=1,NuNum
                         complex_xsum(ii)=complex_xsum(ii)+real(complex_x(ii))
                     enddo
                 enddo

                 do ii=1,NuNum
                     vec_N_old(ii)=vec_N_old(ii)*alpha_real(0)+2d0*complex_xsum(ii)
                     vec_N_old(ii)=max(0.0E0,vec_N_old(ii))
                     if(vec_N_old(ii)<1d-30) vec_N_old(ii)=0.0E0
                     N_FR(Ixy,Iz,ii)=vec_N_old(ii)*1D-24
                     all_concentrations(Ixy, Iz,ii,I_BU) = vec_N_old(ii)*1d-24
                 enddo



             END DO
         END DO


    return
    contains


    SUBROUTINE NESTED_SGS_AXB
    implicit none
    INTEGER(4)                                  :: i
    INTEGER(4)                                  :: m
    INTEGER(4)                                  :: is
    INTEGER(4)                                  :: iter
    INTEGER(4),PARAMETER                        :: max_iter=1000
    REAL(8)                                     :: max_err
    COMPLEX(8)                                  :: Xo(NuNum)
    COMPLEX(8)                                  :: tmp1

    do iter=1,max_iter
        Xo(:)=complex_x(:)
        is=0
        do i=1,NuNum
            tmp1=0.0
            do m=1,SPS_n_non0(i)
                is  =is+1
                tmp1=tmp1+sps_matA(is)*complex_x( SPS_is2j(is) ) ! do such as GS
            enddo
            complex_x(i)=(complex_b(i)-tmp1)*div(i)
#ifdef RFORM
            !if(abs(real(complex_x(i)))<1.0E-35_r8k) complex_x(i)=1.0E-60_r8k
            if(abs((complex_x(i)))<1.0d-35) complex_x(i)=1.0d-60
#else
            !if(abs(real(complex_x(i)))<1.0E-100_r8k) complex_x(i)=1.0E-99_r8k
            if(abs((complex_x(i)))<1.0d-100) complex_x(i)=1.0d-99
#endif
        enddo
        max_err=0d0
        do i=1,NuNum
#ifdef RFORM
            !max_Err=max(max_err,abs((complex_x(i)-xo(i))/(1.0E-30_r8k+complex_x(i))))
            max_Err=max(max_err,abs((complex_x(i)-xo(i))/max(1.0d-25,abs(xo(i)),abs(complex_x(i)))))
#else
            max_Err=max(max_err,abs((complex_x(i)-xo(i))/(1.0d-70+complex_x(i))))
#endif
        enddo

        if(max_err<1.0d-06 .and. iter>=3) then
            exit
        elseif(iter==max_iter) then
            write(*,*) 'SGS_AXB solver diverges in SGS_AXB routine. max_err=',max_err
            do i=1,NuNum
#ifdef RFORM
                max_Err=max(max_err,abs((complex_x(i)-xo(i))/max(1.0d-25,abs(xo(i)),abs(complex_x(i)))))
#else
                max_Err=max(max_err,abs((complex_x(i)-xo(i))/(1.0d-70+complex_x(i))))
#endif
                write(*,*) i,complex_x(i),xo(i),max_err,abs((complex_x(i)-xo(i))/max(1.0d-25,abs(xo(i)),abs(complex_x(i))))
            enddo
            !call end_stream
        endif
    enddo
    return
    END SUBROUTINE NESTED_SGS_AXB

      end subroutine fr_CRAM_solver


			



! previous

!         DO Ixy_FA = 1, Nxy_FA
!             Ixy = I_FA(Ixy_FA)
!             DO Iz = IzFuelBot, IzFuelTop
!
!                 do jj=1,NuNum
!                     vec_N_old(jj)=N_FR_old(Ixy,Iz,jj)*1d24  ! #/(barn-cm)         
!                 enddo
!             
!                 ! Set full matrix A
!                 call fr_update_xs_data(Ixy,Iz)
!             !    NuNum = NuNum
!                 if(allocated(SPS_n_non0)) deallocate(SPS_n_non0)
!                 allocate(SPS_n_non0(NuNum))
!                 SPS_n_non0(:) = 0
!                 do ii=1,NuNum
!                     tmp_non0_idx = 0
!                     do jj=1,NuNum
!                         if(ii == jj) cycle
!                         if(abs(A_matrix(jj,ii))>0.0E0) then
!                             tmp_non0_idx = tmp_non0_idx + 1
!                         endif
!                     enddo
!                     SPS_n_non0(ii) = tmp_non0_idx
!                 enddo
!             
!                 SPS_SUM_non0=sum(SPS_n_non0(:))
!                 if(allocated(SPS_is2J)) deallocate(SPS_is2J)
!                 allocate(SPS_is2J(SPS_SUM_non0))
!                 SPS_is2J = 0
!                 is=0
!                 do ii=1,NuNum
!                     !tmp_non0_idx = 0
!                     do jj=1, NuNum
!                         if(ii == jj) cycle
!                         if(abs(A_matrix(jj,ii))>0) then
!                             !tmp_non0_idx = tmp_non0_idx + 1
!                             is=is+1
!                             SPS_is2J(is) = jj
!                         endif
!                     enddo
!                 enddo
!             
!             
!                 if(allocated(SPS_matA)) deallocate(SPS_matA)
!                 allocate(SPS_matA(SPS_sum_non0))
!             	
!             	    ! Make matrix At
!                 SPS_matA(:)=0d0
!                 is=0
!                 do ii=1,NuNum
!                     do jj=1,SPS_n_non0(ii)
!                         is=is+1
!                         SPS_matA(is)=A_matrix(SPS_is2j(is),ii)*dT
!        
!                     enddo
!                     diag_matA(ii)=A_matrix(ii,ii)*dT
!                 enddo
!             	
!
!                 ! Solve n=n0*exp(At) with CRAM
!                 complex_xsum(:)=0d0 
!                 complex_x(:)   =1d0
!                 do k=1,n_CRAM_k
!                     do ii=1,NuNum
!                         div(ii)=1d0/(diag_matA(ii)-theta_real(k)-theta_img(k)*(0.0E0,1.0E0))
!                         complex_b(ii)=(alpha_real(k)+alpha_img(k)*(0.0E0,1.0E0))*vec_N_old(ii)
!                     enddo
!             
!                     call NESTED_SGS_AXB
!             
!                     do ii=1,NuNum
!                         complex_xsum(ii)=complex_xsum(ii)+real(complex_x(ii))
!                     enddo
!                 enddo
!             
!                 do ii=1,NuNum
!                     vec_N_old(ii)=vec_N_old(ii)*alpha_real(0)+2d0*complex_xsum(ii)
!                     vec_N_old(ii)=max(0.0E0,vec_N_old(ii))
!                     if(vec_N_old(ii)<1d-30) vec_N_old(ii)=0.0E0
!                     N_FR(Ixy,Iz,ii)=vec_N_old(ii)*1D-24
!                     all_concentrations(Ixy, Iz,ii,I_BU) = vec_N_old(ii)*1d-24
!                 enddo
!             
!
!
!             END DO
!         END DO
!
!
!    return
!    contains
!
!
!    SUBROUTINE NESTED_SGS_AXB
!    implicit none
!    INTEGER(4)                                  :: i
!    INTEGER(4)                                  :: m
!    INTEGER(4)                                  :: is
!    INTEGER(4)                                  :: iter
!    INTEGER(4),PARAMETER                        :: max_iter=1000
!    REAL(8)                                     :: max_err
!    COMPLEX(8)                                  :: Xo(NuNum)
!    COMPLEX(8)                                  :: tmp1
!
!    do iter=1,max_iter
!        Xo(:)=complex_x(:)
!        is=0
!        do i=1,NuNum
!            tmp1=0.0
!            do m=1,SPS_n_non0(i)
!                is  =is+1
!                tmp1=tmp1+sps_matA(is)*complex_x( SPS_is2j(is) ) ! do such as GS 
!            enddo
!            complex_x(i)=(complex_b(i)-tmp1)*div(i)
!#ifdef RFORM
!            !if(abs(real(complex_x(i)))<1.0E-35_r8k) complex_x(i)=1.0E-60_r8k
!            if(abs((complex_x(i)))<1.0d-35) complex_x(i)=1.0d-60
!#else
!            !if(abs(real(complex_x(i)))<1.0E-100_r8k) complex_x(i)=1.0E-99_r8k
!            if(abs((complex_x(i)))<1.0d-100) complex_x(i)=1.0d-99
!#endif
!        enddo
!        max_err=0d0
!        do i=1,NuNum
!#ifdef RFORM
!            !max_Err=max(max_err,abs((complex_x(i)-xo(i))/(1.0E-30_r8k+complex_x(i))))
!            max_Err=max(max_err,abs((complex_x(i)-xo(i))/max(1.0d-25,abs(xo(i)),abs(complex_x(i)))))
!#else
!            max_Err=max(max_err,abs((complex_x(i)-xo(i))/(1.0d-70+complex_x(i))))
!#endif
!        enddo
!
!        if(max_err<1.0d-06 .and. iter>=3) then
!            exit
!        elseif(iter==max_iter) then
!            write(*,*) 'SGS_AXB solver diverges in SGS_AXB routine. max_err=',max_err
!            do i=1,NuNum
!#ifdef RFORM
!                max_Err=max(max_err,abs((complex_x(i)-xo(i))/max(1.0d-25,abs(xo(i)),abs(complex_x(i)))))
!#else
!                max_Err=max(max_err,abs((complex_x(i)-xo(i))/(1.0d-70+complex_x(i))))
!#endif
!                write(*,*) i,complex_x(i),xo(i),max_err,abs((complex_x(i)-xo(i))/max(1.0d-25,abs(xo(i)),abs(complex_x(i))))
!            enddo
!            !call end_stream
!        endif
!    enddo
!    return
!    END SUBROUTINE NESTED_SGS_AXB
!
!      end subroutine fr_CRAM_solver
	  


      ! deallocation of previously allocated memory for FastDepletion
      subroutine fr_deallocate_arrays
         use Inc_FastDepletion
         implicit none

         if (allocated(all_names)) deallocate(all_names)
         if (allocated(all_concentrations)) deallocate(all_concentrations)
         if (allocated(A_matrix)) deallocate(A_matrix)         
         if (allocated(cram_alpha)) deallocate(cram_alpha)
         if (allocated(cram_theta)) deallocate(cram_theta)
         if (allocated(vec_N_old)) deallocate(vec_N_old)
         if (allocated(mat_I_all)) deallocate(mat_I_all)
         if (allocated(xs_a)) deallocate(xs_a)
         if (allocated(xs_f)) deallocate(xs_f)
         if (allocated(xs_n2n)) deallocate(xs_n2n)
         if (allocated(sum_CRAM_vec)) deallocate(sum_CRAM_vec)

      end subroutine fr_deallocate_arrays
   

   end module FastDepletion
#endif
