#ifdef siarhei_delete


      MODULE Mod_Sol2N



#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

         implicit none

      type UNMNodal
         real(8) :: width
         real(8), allocatable  :: fiss_spec  (:)
         real(8), allocatable  :: flux_avg   (:)
         real(8), allocatable  :: source0    (:)
         real(8), allocatable  :: source1    (:)
         real(8), allocatable  :: source2    (:)
         real(8), allocatable  :: adf        (:)
         real(8), allocatable  :: diff_coeff (:)
         real(8), allocatable  :: sigma_r    (:)
         real(8), allocatable  :: nu_sigma_f (:)
         real(8), allocatable  :: sigma_s    (:,:)
      end type UNMNodal
      integer :: n_groups
      real(8), allocatable :: face_flux (:)
      real(8), allocatable :: face_curr (:)
      real(8), allocatable :: corr_fact (:,:)
      real(8), allocatable :: C0        (:,:)
      real(8), allocatable :: C1        (:,:)
      real(8), allocatable :: C2        (:,:)
      real(8), allocatable :: C3        (:,:)
      real(8), allocatable :: C4        (:,:)
      real(8), allocatable :: B         (:,:)
      real(8), allocatable :: D         (:,:)
      real(8), allocatable :: M0        (:,:,:)
      real(8), allocatable :: M1        (:,:,:)
      real(8), allocatable :: M2        (:,:,:)
      real(8), allocatable :: M3        (:,:,:)
      real(8), allocatable :: M4        (:,:,:)
      real(8), allocatable :: MN        (:,:)
      real(8), allocatable :: MNN       (:,:)
      real(8), allocatable :: MR        (:,:)
      real(8), allocatable :: ME        (:,:)
      real(8), allocatable :: ME22      (:,:)
      real(8), allocatable :: nPinv     (:)
      real(8), allocatable :: QUNM      (:)
      real(8), allocatable :: Q0        (:)
      real(8), allocatable :: Q1        (:)
      real(8), allocatable :: Q2        (:)
      real(8), allocatable :: FlFr      (:)
      real(8), allocatable :: BlBr      (:)
      real(8), allocatable :: FlS1      (:)
      real(8), allocatable :: BlS2      (:)
      real(8), allocatable :: QQ        (:)
      real(8), allocatable :: sol       (:)
      integer :: OPT_NEM, OPT_ANM, LEFT, RGHT
      parameter( OPT_NEM=0, OPT_ANM=1 )
      parameter( LEFT=1, RGHT=2 )

      contains

      subroutine allocate_one_dim_node(node,ngrp)
      implicit none
      type(UNMNodal) :: node
      integer :: ngrp
      allocate( node%fiss_spec(ngrp) )

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [allocate_one_dim_node] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      allocate( node%flux_avg(ngrp) )
      allocate( node%source0(ngrp) )
      allocate( node%source1(ngrp) )
      allocate( node%source2(ngrp) )
      allocate( node%adf(ngrp) )
      allocate( node%diff_coeff(ngrp) )
      allocate( node%sigma_r(ngrp) )
      allocate( node%nu_sigma_f(ngrp) )
      allocate( node%sigma_s(ngrp,ngrp) )
      node%fiss_spec(:)  = 0.0d0
      node%flux_avg(:)   = 0.0d0
      node%source0(:)    = 0.0d0
      node%source1(:)    = 0.0d0
      node%source2(:)    = 0.0d0
      node%adf(:)        = 1.0d0
      node%diff_coeff(:) = 0.0d0
      node%sigma_r(:)    = 0.0d0
      node%nu_sigma_f(:) = 0.0d0
      node%sigma_s(:,:)  = 0.0d0
      return
      end subroutine allocate_one_dim_node

      subroutine deallocate_one_dim_node(node)
      implicit none
      type(UNMNodal) :: node
      deallocate( node%fiss_spec )

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [deallocate_one_dim_node] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      deallocate( node%flux_avg )
      deallocate( node%source0 )
      deallocate( node%source1 )
      deallocate( node%source2 )
      deallocate( node%adf )
      deallocate( node%diff_coeff )
      deallocate( node%sigma_r )
      deallocate( node%nu_sigma_f )
      deallocate( node%sigma_s )
      return
      end subroutine deallocate_one_dim_node

      subroutine initialize_two_node_unm_solver(ngrp)
      implicit none
      integer :: ngrp
      n_groups = ngrp

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [initialize_two_node_unm_solver] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      allocate( face_flux(ngrp) )
      allocate( face_curr(ngrp) )
      allocate( corr_fact(ngrp,2) )
      allocate( C0(ngrp,2) )
      allocate( C1(ngrp,2) )
      allocate( C2(ngrp,2) )
      allocate( C3(ngrp,2) )
      allocate( C4(ngrp,2) )
      allocate(  B(ngrp,2) )
      allocate(  D(ngrp,2) )
      allocate( M0(ngrp,ngrp,2) )
      allocate( M1(ngrp,ngrp,2) )
      allocate( M2(ngrp,ngrp,2) )
      allocate( M3(ngrp,ngrp,2) )
      allocate( M4(ngrp,ngrp,2) )
      allocate( MN(ngrp,ngrp) )
      allocate( MNN(ngrp,ngrp) )
      allocate( MR(ngrp,ngrp) )
      allocate( ME(2*ngrp,2*ngrp) )
      allocate( ME22(ngrp,ngrp) )
      allocate( nPinv(ngrp) )
      allocate( QUNM(ngrp) )
      allocate( Q0(ngrp) )
      allocate( Q1(ngrp) )
      allocate( Q2(ngrp) )
      allocate( FlFr(ngrp) )
      allocate( BlBr(ngrp) )
      allocate( FlS1(ngrp) )
      allocate( BlS2(ngrp) )
      allocate( QQ(2*ngrp) )
      allocate( sol(2*ngrp) )
      return
      end subroutine initialize_two_node_unm_solver

      subroutine finalize_two_node_unm_solver
      implicit none

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [finalize_two_node_unm_solver] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      deallocate( face_flux )
      deallocate( face_curr )
      deallocate( corr_fact )
      deallocate( C0 )
      deallocate( C1 )
      deallocate( C2 )
      deallocate( C3 )
      deallocate( C4 )
      deallocate( B )
      deallocate( D )
      deallocate( M0 )
      deallocate( M1 )
      deallocate( M2 )
      deallocate( M3 )
      deallocate( M4 )
      deallocate( MN )
      deallocate( MNN )
      deallocate( MR )
      deallocate( ME )
      deallocate( ME22 )
      deallocate( nPinv )
      deallocate( QUNM )
      deallocate( Q0 )
      deallocate( Q1 )
      deallocate( Q2 )
      deallocate( FlFr )
      deallocate( BlBr )
      deallocate( FlS1 )
      deallocate( BlS2 )
      deallocate( QQ )
      deallocate( sol )
      return
      end subroutine finalize_two_node_unm_solver


      subroutine setup_resp_mtrx_NEM(node,side,k_eff)
      implicit none
      type(UNMNodal) :: node
      integer :: side
      real(8) :: k_eff
      real(8) :: node_width, delta
      integer :: g, gg


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [setup_resp_mtrx_NEM] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      node_width = node%width;
      M0(:,:,side) = 0.0
      do g=1,n_groups
         B(g,side) = node%diff_coeff(g)/node_width
         D(g,side) = node%diff_coeff(g)
         M0(g,g,side) = node%sigma_r(g)
         if (abs(node%fiss_spec(g))>1d-10) then
            do gg=1,n_groups
               M0(g,gg,side) = M0(g,gg,side) - node%fiss_spec(g)*node%nu_sigma_f(gg)/k_eff
            enddo
         endif
         do gg=1,n_groups
            M0(g,gg,side) = M0(g,gg,side) - node%sigma_s(g,gg)
            M1(g,gg,side) = -M0(g,gg,side)
            M2(g,gg,side) = -M0(g,gg,side)
            M3(g,gg,side) = 1.000000000000000000e-01*M0(g,gg,side)
            M4(g,gg,side) = 7.142857142857142857e-02*M0(g,gg,side)
         enddo
         delta = B(g,side)/node_width
         M3(g,g,side) = M3(g,g,side) +  6.0*delta
         M4(g,g,side) = M4(g,g,side) + 10.0*delta
      enddo

      return
      end subroutine setup_resp_mtrx_NEM


      subroutine setup_resp_mtrx_ANM( node, side, k_eff )
      implicit none
      type(UNMNodal) :: node
      integer :: side
      real(8) :: k_eff
      real(8) :: node_width, h2, h4, Tg, sum_elm, temp, temp11, temp12, temp21, temp22, one_over_d
      integer :: g, gg, p, q
      integer :: k


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [setup_resp_mtrx_ANM] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      node_width = node%width
      h2 = node_width*node_width
      h4 = h2*h2
      M0(:,:,side) = 0.0
      do g=1,n_groups
         Tg = 1.0/node%diff_coeff(g)
         B(g,side) = node%diff_coeff(g)/node_width
         D(g,side) = node%diff_coeff(g)
         M0(g,g,side) = node%sigma_r(g)
         if (abs(node%fiss_spec(g))>1d-10) then
            do gg=1,n_groups
               M0(g,gg,side) = M0(g,gg,side) - node%fiss_spec(g)*node%nu_sigma_f(gg)/k_eff
            enddo
         endif
         do gg=1,n_groups
            M0(g,gg,side) = M0(g,gg,side) - node%sigma_s(g,gg)
            M1(g,gg,side) = -M0(g,gg,side)
            M2(g,gg,side) = -M0(g,gg,side)
            MN(g,gg)      =  M0(g,gg,side)*Tg
         enddo
      enddo
      do g=1,n_groups
         do gg=1,n_groups
            sum_elm = 0.0
            do k=1,n_groups
               sum_elm = sum_elm + MN(g,k)*MN(k,gg)
            enddo
            MNN(g,gg) = sum_elm
         enddo
      enddo
      temp11 = 1.8181818181818182e-01*h2
      temp12 = 7.5757575757575758e-04*h4
      temp21 = 1.3636363636363636e-02*h2
      temp22 = 1.8037518037518038e-05*h4
      do g=1,n_groups
         do gg=1,n_groups
            M3(g,gg,side) = temp11*MN(g,gg) + temp12*MNN(g,gg);
            MR(g,gg)      = temp21*MN(g,gg) + temp22*MNN(g,gg);
         enddo
         M3(g,g,side) = M3(g,g,side) + 6.0
         MR(g,g)      = MR(g,g)      + 1.0
      enddo
      do g=1,n_groups
         one_over_d = 1.0/MR(g,g)
         do gg=g+1,n_groups
            MR(g,gg) = MR(g,gg)*one_over_d
         enddo
         do gg=1,n_groups
            M3(g,gg,side) = M3(g,gg,side)*one_over_d
         enddo
         do p=g+1,n_groups
            do q=g+1,n_groups
               MR(p,q) = MR(p,q) - MR(p,g)*MR(g,q)
            enddo
            do q=1,n_groups
               M3(p,q,side) = M3(p,q,side) - MR(p,g)*M3(g,q,side)
            enddo
         enddo
      enddo
      do g=n_groups,1,-1
         do p=g-1,1,-1
            do q=1,n_groups
               M3(p,q,side) = M3(p,q,side) - MR(p,g)*M3(g,q,side)
            enddo
         enddo
         temp = node%diff_coeff(g)/h2
         do q=1,n_groups
            M3(g,q,side) = M3(g,q,side)*temp
         enddo
      enddo
      temp11 = 1.5384615384615385e-01*h2
      temp12 = 3.7462537462537463e-04*h4
      temp21 = 8.2417582417582418e-03*h2
      temp22 = 6.9375069375069375e-06*h4
      do g=1,n_groups
         do gg=1,n_groups
            M4(g,gg,side) = temp11*MN(g,gg) + temp12*MNN(g,gg)
            MR(g,gg)      = temp21*MN(g,gg) + temp22*MNN(g,gg)
         enddo
         M4(g,g,side) = M4(g,g,side) + 10.0
         MR(g,g)      = MR(g,g)      + 1.0
      enddo
      do g=1,n_groups
         one_over_d = 1.0/MR(g,g)
         do gg=g+1,n_groups
            MR(g,gg) = MR(g,gg)*one_over_d
         enddo
         do gg=1,n_groups
            M4(g,gg,side) = M4(g,gg,side)*one_over_d
         enddo
         do p=g+1,n_groups
            do q=g+1,n_groups
               MR(p,q) = MR(p,q) - MR(p,g)*MR(g,q)
            enddo
            do q=1,n_groups
               M4(p,q,side) = M4(p,q,side) - MR(p,g)*M4(g,q,side)
            enddo
         enddo
      enddo
      do g=n_groups,1,-1
         do p=g-1,1,-1
            do q=1,n_groups
               M4(p,q,side) = M4(p,q,side) - MR(p,g)*M4(g,q,side)
            enddo
         enddo
         temp = node%diff_coeff(g)/h2
         do q=1,n_groups
            M4(g,q,side) = M4(g,q,side)*temp
         enddo
      enddo

      return
      end subroutine setup_resp_mtrx_ANM

      subroutine solve_even_modes( node, side )
      implicit none
      type(UNMNodal) :: node
      integer :: side
      real(8) :: width, sum
      integer ::g,gg,i,j


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [solve_even_modes] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      width = node%width
      do g=1,n_groups
         sum = 0.0
         do gg=1,n_groups
            sum = sum + M0(g,gg,side)*node%flux_avg(gg)
         enddo
         Q0(g) = (1.0/12.0)*width*(sum-node%source0(g))/B(g,side)
      enddo
      do g=1,n_groups
         sum = node%source2(g)
         do gg=1,n_groups
            MN(g,gg) = M4(g,gg,side) - (1.0/6.0)*M2(g,gg,side)
            sum = sum + M2(g,gg,side)*Q0(gg)
         enddo
         QUNM(g) = -sum
      enddo
      do g=1,n_groups
         do i=g+1,n_groups
            do j=g+1,n_groups
               MN(i,j) = MN(i,j) - MN(i,g)/MN(g,g)*MN(g,j)
            enddo
            QUNM(i) = QUNM(i) - MN(i,g)/MN(g,g)*QUNM(g)
         enddo
      enddo
      do g=n_groups,1,-1
         sum = QUNM(g)
         do i=n_groups,g+1,-1
            sum = sum - MN(g,i)*C4(i,side)
         enddo
         C4(g,side) = sum / MN(g,g)
         C2(g,side) = Q0(g) - (1.0/6.0)*C4(g,side)
         C0(g,side) = node%flux_avg(g)
      enddo

      return
      end subroutine solve_even_modes

      subroutine solve_left_bndry_odd_modes( ib, node_r )
      implicit none
      type(UNMNodal) :: node_r
      integer :: ib
      integer :: g, gg, r, s
      real(8) :: adf, temp


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [solve_left_bndry_odd_modes] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if ( ib == 3 )then
         do g=1,n_groups
            C1(g,RGHT) = C0(g,RGHT)+C2(g,RGHT)
            Q2(g) = -node_r%source1(g)
            do gg=1,n_groups
               ME22(g,gg) = M3(g,gg,RGHT)
            enddo
         enddo
         do g=1,n_groups                          !   Elemenate the matrix M1r
            do r=1,n_groups
               Q2(r) = Q2(r) - M1(r,g,RGHT)*C1(g,RGHT)
            enddo
         enddo
      else
         do g=1,n_groups
            adf = node_r%adf(g)
            nPinv(g) = B(g,RGHT)/(2.0*B(g,RGHT)+0.5)
            Q1(g) = (0.5*adf*(C0(g,RGHT)+C2(g,RGHT)))/B(g,RGHT) + 6.0*C2(g,RGHT) + C4(g,RGHT)
            Q2(g) = -node_r%source1(g)
            do gg=1,n_groups
               ME22(g,gg) = M3(g,gg,RGHT)
            enddo
         enddo
         do g=1,n_groups                          !   Elemenate the matrix M1r
            do r=1,n_groups
               temp = M1(r,g,RGHT)*nPinv(g)
               ME22(r,g) = ME22(r,g) - temp
               Q2(r) = Q2(r) - temp*Q1(g)
            enddo
         enddo
      endif

      do g=1,n_groups                          !   Elemenate the Lower Triangular Part of E
         temp = 1.0/ME22(g,g)
         do r=g+1,n_groups
            do s=g+1,n_groups
               ME22(r,s) = ME22(r,s) - ME22(r,g)*temp*ME22(g,s)
            enddo
            Q2(r) = Q2(r) - ME22(r,g)*temp*Q2(g)
         enddo
      enddo
      do r=n_groups,1,-1                       !   Solve for C3r
         temp = Q2(r)
         do s=n_groups,r+1,-1
            temp = temp - ME22(r,s)*C3(s,RGHT)
         enddo
         C3(r,RGHT) = temp/ME22(r,r)
      enddo

      IF ( ib == 2 )then
          do g=1,n_groups                          !   Solve for C1r
             C1(g,RGHT) = nPinv(g)*(Q1(g)-C3(g,RGHT))
          enddo
      endif

      return
      end subroutine solve_left_bndry_odd_modes


      subroutine solve_odd_modes( node_l, node_r )
      implicit none
      type(UNMNodal) :: node_l, node_r
      real(8) :: sum
      integer :: g, gg, i, j
      integer :: k


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [solve_odd_modes] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      do g=1,n_groups
         FLFr(g) = node_r%adf(g)/node_l%adf(g)
         BlBr(g) = B(g,RGHT)/B(g,LEFT)
         FLS1(g) = FLFr(g)*(C0(g,RGHT)+C2(g,RGHT))              &
              &                 - (C0(g,LEFT)+C2(g,LEFT))
         BlS2(g) = - BlBr(g)*(6.0*C2(g,RGHT)+C4(g,RGHT))        &
              &                 - (6.0*C2(g,LEFT)+C4(g,LEFT))
      enddo
      do g=1,n_groups
         sum = 0.0
         do gg=1,n_groups
            ME(g,gg) = 2*M3(g,gg,LEFT)*(BlBr(gg)+FlFr(gg))      &
                 &                     - M1(g,gg,LEFT)*FlFr(gg)
            ME(g,n_groups+gg) = M3(g,gg,LEFT)*BlBr(gg)
            ME(n_groups+g,gg) = M1(g,gg,RGHT)
            ME(n_groups+g,n_groups+gg) = M3(g,gg,RGHT)
            sum = sum + M1(g,gg,LEFT)*FlS1(gg)                  &
                 &              + M3(g,gg,LEFT)*(BlS2(gg)-2*FlS1(gg))
         enddo
         QQ(g) = -node_l%source1(g) - sum
         QQ(n_groups+g) = -node_r%source1(g)
      enddo

      do i=1,2*n_groups
         do j=i+1,2*n_groups
            do k=i+1,2*n_groups
               ME(j,k) = ME(j,k) - ME(j,i)/ME(i,i)*ME(i,k)
            enddo
            QQ(j) = QQ(j) - ME(j,i)/ME(i,i)*QQ(i)
         enddo
      enddo
      do j=2*n_groups,1,-1
         sum = QQ(j)
         do k=2*n_groups,j+1,-1
            sum = sum - ME(j,k)*sol(k)
         enddo
         sol(j) = sum/ME(j,j)
      enddo
      do g=1,n_groups
         C1(g,RGHT) = sol(g)
         C3(g,RGHT) = sol(n_groups+g)
         C1(g,LEFT) = FlS1(g) - FlFr(g)*C1(g,RGHT)
         C3(g,LEFT) = BlS2(g) - 2*C1(g,LEFT) + BlBr(g)*(2*C1(g,RGHT)+C3(g,RGHT))
      enddo
      end subroutine solve_odd_modes


      subroutine solve_right_bndry_odd_modes( ib, node_l )
      implicit none
      type(UNMNodal) :: node_l
      integer :: ib
      integer :: g, gg, r, s
      real(8) :: adf, temp


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [solve_right_bndry_odd_modes] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      IF ( ib == 3 ) then
         do g=1,n_groups
            adf = node_l%adf(g)
            Q1(g) = -adf*(C0(g,LEFT)+C2(g,LEFT))
            Q2(g) = -node_l%source1(g)
            do gg=1,n_groups
               ME22(g,gg) = M3(g,gg,LEFT)
            enddo
         enddo
         do g=1,n_groups                          !   Elemenate the matrix M1r
            do r=1,n_groups
               temp = M1(r,g,LEFT)
               Q2(r) = Q2(r) - temp*Q1(g)
            enddo
         enddo
      else
         do g=1,n_groups
            adf = node_l%adf(g)
            nPinv(g) = B(g,LEFT)/(2.0*B(g,LEFT)+0.5)
            Q1(g) = -(0.5*adf*(C0(g,LEFT)+C2(g,LEFT)))/B(g,LEFT) - 6.0*C2(g,LEFT) - C4(g,LEFT)
            Q2(g) = -node_l%source1(g)
            do gg=1,n_groups
               ME22(g,gg) = M3(g,gg,LEFT)
            enddo
         enddo
         do g=1,n_groups                          !   Elemenate the matrix M1r
            do r=1,n_groups
               temp = M1(r,g,LEFT)*nPinv(g)
               ME22(r,g) = ME22(r,g) - temp
               Q2(r) = Q2(r) - temp*Q1(g)
            enddo
         enddo
      endif

      do g=1,n_groups                          !   Elemenate the Lower Triangular Part of E
         temp = 1.0/ME22(g,g)
         do r=g+1,n_groups
            do s=g+1,n_groups
               ME22(r,s) = ME22(r,s) - ME22(r,g)*temp*ME22(g,s)
            enddo
            Q2(r) = Q2(r) - ME22(r,g)*temp*Q2(g)
         enddo
      enddo
      do r=n_groups,1,-1                       !   Solve for C3l
         temp = Q2(r)
         do s=n_groups,r+1,-1
            temp = temp - ME22(r,s)*C3(s,LEFT)
         enddo
         C3(r,LEFT) = temp/ME22(r,r)
      enddo

      IF ( ib == 3 )then
         do g=1,n_groups                          !   Solve for C1r
            C1(g,LEFT) = -(C0(g,LEFT)+C2(g,LEFT))
         enddo
      else
         do g=1,n_groups                          !   Solve for C1l
            C1(g,LEFT) = nPinv(g)*(Q1(g)-C3(g,LEFT))
         enddo
      endif

      return
      end subroutine solve_right_bndry_odd_modes


      subroutine solve_left_bndry_node_unm( ib, node_r, k_eff, opt )
      implicit none
      type(UNMNodal) :: node_r
      real(8) :: k_eff
      integer :: opt, ib
      real(8) :: flux_r
      integer :: g


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [solve_left_bndry_node_unm] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if ( opt .EQ. OPT_NEM ) then
         call setup_resp_mtrx_NEM( node_r, RGHT, k_eff )
      else
         call setup_resp_mtrx_ANM( node_r, RGHT, k_eff )
      endif
      call solve_even_modes( node_r, RGHT )
      call solve_left_bndry_odd_modes( ib, node_r )
      do g=1,n_groups
         flux_r  = C0(g,RGHT) - C1(g,RGHT) + C2(g,RGHT)
         face_flux(g) = flux_r*node_r%adf(g)
         face_curr(g) = -B(g,RGHT)*(2.0*C1(g,RGHT)-6.0*C2(g,RGHT)+C3(g,RGHT)-C4(g,RGHT))
         corr_fact(g,RGHT) = -(2*B(g,RGHT)*(flux_r-C0(g,RGHT))-face_curr(g))/(flux_r+C0(g,RGHT))
      enddo

      return
      end subroutine solve_left_bndry_node_unm


      subroutine solve_two_node_unm( node_l, node_r, k_eff, opt )
      implicit none
      type(UNMNodal) :: node_l, node_r
      real(8) :: k_eff
      integer :: opt
      real(8) :: flux_l, flux_r
      integer :: g


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [solve_two_node_unm] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if ( opt .EQ. OPT_NEM ) then
         call setup_resp_mtrx_NEM( node_l, LEFT, k_eff )
         call setup_resp_mtrx_NEM( node_r, RGHT, k_eff )
      else
         call setup_resp_mtrx_ANM( node_l, LEFT, k_eff )
         call setup_resp_mtrx_ANM( node_r, RGHT, k_eff )
      endif
      call solve_even_modes( node_l, LEFT )
      call solve_even_modes( node_r, RGHT )
      call solve_odd_modes( node_l, node_r )
      do g=1,n_groups
         flux_l  = C0(g,LEFT) + C1(g,LEFT) + C2(g,LEFT)
         face_flux(g) = flux_l*node_l%adf(g)
         face_curr(g) = -B(g,LEFT)*(2.0*C1(g,LEFT)+6.0*C2(g,LEFT)+C3(g,LEFT)+C4(g,LEFT))
         flux_r  = face_flux(g)/node_r%adf(g)
         corr_fact(g,LEFT) = -(2*B(g,LEFT)*(flux_l-C0(g,LEFT))+face_curr(g))/(flux_l+C0(g,LEFT))
         corr_fact(g,RGHT) = -(2*B(g,RGHT)*(flux_r-C0(g,RGHT))-face_curr(g))/(flux_r+C0(g,RGHT))
      enddo

      return
      end subroutine solve_two_node_unm


      subroutine solve_right_bndry_node_unm( ib, node_l, k_eff, opt )
      implicit none
      type(UNMNodal) :: node_l
      real(8) :: k_eff
      integer :: opt, ib
      real(8) :: flux_l
      integer :: g


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [solve_right_bndry_node_unm] in Mod_Sol2N'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if ( opt .EQ. OPT_NEM ) then
         call setup_resp_mtrx_NEM( node_l, LEFT, k_eff )
      else
         call setup_resp_mtrx_ANM( node_l, LEFT, k_eff )
      endif
      call solve_even_modes( node_l, LEFT )
      call solve_right_bndry_odd_modes( ib, node_l )
      do g=1,n_groups
         flux_l  = C0(g,LEFT) + C1(g,LEFT) + C2(g,LEFT)
         face_flux(g) = flux_l*node_l%adf(g)
         face_curr(g) = -B(g,LEFT)*(2.0*C1(g,LEFT)+6.0*C2(g,LEFT)+C3(g,LEFT)+C4(g,LEFT))
         corr_fact(g,LEFT) = -(2*B(g,LEFT)*(flux_l-C0(g,LEFT))+face_curr(g))/(flux_l+C0(g,LEFT))
      enddo

      return
      end subroutine solve_right_bndry_node_unm


      end module Mod_Sol2N


#endif
