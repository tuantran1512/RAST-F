#ifdef tuan_tr_test
!
      MODULE Mod_Soldhat
      USE Inc_maXS
      USE Inc_Control
      USE Inc_FluxVar
      USE Inc_Lscoef
!     ! USE Mod_Sol2N
      USE Inc_Geometry
      USE Inc_3D
      USE Inc_DF
      USE Inc_XS_File, only: Leakage_Table, Flag_Leakage, leak_ratio, Flag_LC_mean
      USE Inc_RP, only: AxialComp, I_LP_1N
      USE Inc_Option
      USE Inc_Transient, ONLY: Flag_Transient
      use Mod_GetNode, only: new_asym_itab
      use mod_charedit, only: print_msg


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

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

      CONTAINS

#ifdef siarhei_delete
      subroutine updatedhat
      implicit none
      type(UNMNodal) node_l, node_r
      real(8) :: rl(n_group)
      real(8) :: rlp(n_group)
      real(8) :: rm1(n_group)
      real(8) :: rmp1(n_group)
      real(8) :: rm2(n_group)
      real(8) :: rmp2(n_group)
      real(8) :: df(n_group)
      real(8) :: dfp(n_group)
      integer :: g, gg, kp, kl, klp, lm, lp
      integer :: l, k, m, i, j
      integer :: ibeg, jbeg, isub, kbot, ktop, kbeg
      integer :: Ixy, Iz, Ixy_1N, I_Tab


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [updatedhat] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call initialize_two_node_unm_solver(n_group)
      call Allocate_one_dim_node(node_l,n_group)
      call Allocate_one_dim_node(node_r,n_group)

      call updcurnxy
      call updtlxy

      ! sweep along nodes in radial plane.
      isub=1
      kbot=1
      ktop=nz
      if (isub==1) then
         kbeg=1
      else
         kbeg=kbot-1
      endif

      do k=kbot,ktop
         ! sweep along x-direction (west-east).
         do j=1,Ny
            ! solve one-node problem for western b.c. of row j if j(edge)<>0 and
            ! update dnw(edge) using dfw(edge) which exists.
            if (idomx==1) then
               i=Ix_Start_y(j)
               l=nodel(i,j)
               select case (BC_Lx)
               case (1)
                  do m=1,n_group
                     dnw(m,l,k)=0d0
                  enddo
               case (2,3)
                  do m=1,n_group
                     rl(m) =L_0y_3D(m,l,k)+L_0z_3D(m,l,k)-src_tr(m,l,k)
                     rm1(m)=L_1x_3D(m,l,k)
                     rm2(m)=L_2x_3D(m,l,k)
                     rm1(m)= rm1(m)*0.5d0
                     rm2(m)= rm2(m)*0.5d0
                  enddo
                  node_r%width=MeshSize_x(i)
                  do g=1,n_group
                     node_r%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                     node_r%source0(g)    =-rl(g)
                     node_r%source1(g)    =-rm1(g)
                     node_r%source2(g)    =-rm2(g)
                     node_r%adf(g)        = ADF_Lx(l,k,g)
                     node_r%flux_avg(g)   = Flux(l,k,g)
                     node_r%diff_coeff(g) = D_3D(l,k,g)
                     node_r%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                     node_r%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                     do gg=1,n_group
                        node_r%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                     enddo
                  enddo
                  select case (opt_nodal)
                  case (1)
                     call solve_left_bndry_node_unm(BC_Lx,node_r,keff,OPT_NEM)
                  case (2)
                     call solve_left_bndry_node_unm(BC_Lx,node_r,keff,OPT_ANM)
                  case (3)
                     call solve_left_bndry_node_unm(BC_Lx,node_r,keff,OPT_ANM)
                  case default
                     call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                     stop
                  end select
                  do m=1,n_group
                     dnw(m,l,k)=-dfw(m,l,k)-face_curr(m)/Flux(l,k,m)
                  enddo
               end select
            endif

            ! solve for two-node problems in the x-direction.
            ibeg=Ix_Start_y(j)
            do i=ibeg,Ix_End_y(j)-1
               l=nodel(i,j)
               lm=nodel(i+1,j)
               do m=1,n_group
                  rl(m)   = L_0y_3D(m,l,k)+L_0z_3D(m,l,k)-src_tr(m,l,k)
                  rlp(m)  = L_0y_3D(m,lm,k)+L_0z_3D(m,lm,k)-src_tr(m,lm,k)
                  rm1(m)  = L_1x_3D(m,l,k)
                  rmp1(m) = L_1x_3D(m,lm,k)
                  rm2(m)  = L_2x_3D(m,l,k)
                  rmp2(m) = L_2x_3D(m,lm,k)
                  df(m)   = ADF_Rx(l,k,m)
                  dfp(m)  = ADF_Lx(lm,k,m)
                  rm1(m)  = rm1(m)*0.5d0
                  rm2(m)  = rm2(m)*0.5d0
                  rmp1(m) = rmp1(m)*0.5d0
                  rmp2(m) = rmp2(m)*0.5d0
               enddo
               node_l%width = MeshSize_x(i)
               node_r%width = MeshSize_x(i+1)
               do g=1,n_group
                  node_l%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                  node_r%fiss_spec(g)  = maxs_chi_3d(lm,k,g)
                  node_l%source0(g)    =-rl(g)
                  node_r%source0(g)    =-rlp(g)
                  node_l%source1(g)    =-rm1(g)
                  node_r%source1(g)    =-rmp1(g)
                  node_l%source2(g)    =-rm2(g)
                  node_r%source2(g)    =-rmp2(g)
                  node_l%adf(g)        = df(g)
                  node_r%adf(g)        = dfp(g)
                  node_l%flux_avg(g)   = Flux(l,k,g)
                  node_r%flux_avg(g)   = Flux(lm,k,g)
                  node_l%diff_coeff(g) = D_3D(l,k,g)
                  node_r%diff_coeff(g) = D_3D(lm,k,g)
                  node_l%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                  node_r%sigma_r(g)    = maXs_r_3D(lm,k,g)+maXS_scat_3D(g,g,lm,k)
                  node_l%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                  node_r%nu_sigma_f(g) = nu_maXs_f_3D(lm,k,g)
                  do gg=1,n_group
                     node_l%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                     node_r%sigma_s(g,gg) = maXS_scat_3D(gg,g,lm,k)
                  enddo
               enddo
               select case (opt_nodal)
               case (1)
                  call solve_two_node_unm(node_l,node_r,keff,OPT_NEM)
               case (2)
                  call solve_two_node_unm(node_l,node_r,keff,OPT_ANM)
               case (3)
                  call solve_two_node_unm(node_l,node_r,keff,OPT_ANM)
               case default
                  call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                  stop
               end select
               do m=1,n_group
                  dnw(m,lm,k)=-(dfw(m,lm,k)*(Flux(lm,k,m)-Flux(l,k,m))+face_curr(m)) &
                             &             /(Flux(lm,k,m)+Flux(l,k,m))
               enddo
            enddo

            ! solve one-node problem for eastern b.c. of row j if j(edge)<>0 and
            ! update dnw(edge) using dfw(edge) which exists.
            if (idomx==ndomx) then
               i=Ix_End_y(j)
               l=nodel(i,j)
               lp=nodel(i+1,j)
               select case (BC_Rx)
               case (1)
                  dnw(1,lp,k)=0d0
                  dnw(2,lp,k)=0d0
               case (2,3)
                  do m=1,n_group
                     rl(m)=L_0y_3D(m,l,k)+L_0z_3D(m,l,k)-src_tr(m,l,k)
                     rm1(m)=L_1x_3D(m,l,k)
                     rm2(m)=L_2x_3D(m,l,k)
                     rm1(m)=rm1(m)*0.5d0
                     rm2(m)=rm2(m)*0.5d0
                  enddo
                  node_l%width = MeshSize_x(i)
                  do g=1,n_group
                     node_l%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                     node_l%source0(g)    =-rl(g)
                     node_l%source1(g)    =-rm1(g)
                     node_l%source2(g)    =-rm2(g)
                     node_l%adf(g)        = ADF_Rx(l,k,g)
                     node_l%flux_avg(g)   = Flux(l,k,g)
                     node_l%diff_coeff(g) = D_3D(l,k,g)
                     node_l%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                     node_l%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                     do gg=1,n_group
                        node_l%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                     enddo
                  enddo
                  select case (opt_nodal)
                  case (1)
                     call solve_right_bndry_node_unm(BC_Rx,node_l,keff,OPT_NEM)
                  case (2)
                     call solve_right_bndry_node_unm(BC_Rx,node_l,keff,OPT_ANM)
                  case (3)
                     call solve_right_bndry_node_unm(BC_Rx,node_l,keff,OPT_ANM)
                  case default
                     call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                     stop
                  end select
                  DO m=1,n_group
                     dnw(m,lp,k)=-face_curr(m)/Flux(l,k,m)+dfw(m,lp,k)
                  END DO
               end select
            endif
         enddo

         ! sweep along y-direction (north-south).
         do i=1,nx
            if (idomy==1) then
               j=Iy_Start_x(i)
               l=nodel(i,j)
               ! solve one-node problem for northern b.c. of row i if j(edge)<>0 and
               ! update dnn(edge) using dfn(edge) which exists.
               select case (BC_Ly)
               case (1)
                  do m=1,n_group
                     dnn(m,l,k)=0d0
                  enddo
               case (2,3)
                  do m=1,n_group
                     rl(m)=L_0x_3D(m,l,k)+L_0z_3D(m,l,k)-src_tr(m,l,k)
                     rm1(m)=L_1y_3D(m,l,k)
                     rm2(m)=L_2y_3D(m,l,k)
                     rm1(m)=rm1(m)*0.5d0
                     rm2(m)=rm2(m)*0.5d0
                  enddo
                  node_r%width = MeshSize_y(j)
                  do g=1,n_group
                     node_r%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                     node_r%source0(g)    =-rl(g)
                     node_r%source1(g)    =-rm1(g)
                     node_r%source2(g)    =-rm2(g)
                     node_r%adf(g)        = ADF_Ly(l,k,g)
                     node_r%flux_avg(g)   = Flux(l,k,g)
                     node_r%diff_coeff(g) = D_3D(l,k,g)
                     node_r%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                     node_r%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                     do gg=1,n_group
                        node_r%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                     enddo
                  enddo
                  select case (opt_nodal)
                  case (1)
                     call solve_left_bndry_node_unm(BC_Ly,node_r,keff,OPT_NEM)
                  case (2)
                     call solve_left_bndry_node_unm(BC_Ly,node_r,keff,OPT_ANM)
                  case (3)
                     call solve_left_bndry_node_unm(BC_Ly,node_r,keff,OPT_ANM)
                  case default
                     call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                     stop
                  end select
                  do m=1,n_group
                     dnn(m,l,k)=-dfn(m,l,k)-face_curr(m)/Flux(l,k,m)
                  enddo
               end select
            endif

            ! solve for two-node problems in the y-direction.
            jbeg=Iy_Start_x(i)
            do j=jbeg,Iy_End_x(i)-1
               l=nodel(i,j)
               lm=nodel(i,j+1)
               do m=1,n_group
                  rl(m)   = L_0x_3D(m,l,k)+L_0z_3D(m,l,k)-src_tr(m,l,k)
                  rlp(m)  = L_0x_3D(m,lm,k)+L_0z_3D(m,lm,k)-src_tr(m,lm,k)
                  rm1(m)  = L_1y_3D(m,l,k)
                  rmp1(m) = L_1y_3D(m,lm,k)
                  rm2(m)  = L_2y_3D(m,l,k)
                  rmp2(m) = L_2y_3D(m,lm,k)
                  rm1(m)  = rm1(m)*0.5d0
                  rm2(m)  = rm2(m)*0.5d0
                  rmp1(m) = rmp1(m)*0.5d0
                  rmp2(m) = rmp2(m)*0.5d0
                  df(m)   = ADF_Ry(l,k,m)
                  dfp(m)  = ADF_Ly(lm,k,m)
               enddo
               node_l%width = MeshSize_y(j)
               node_r%width = MeshSize_y(j+1)
               do g=1,n_group
                  node_l%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                  node_r%fiss_spec(g)  = maxs_chi_3d(lm,k,g)
                  node_l%source0(g)    =-rl(g)
                  node_r%source0(g)    =-rlp(g)
                  node_l%source1(g)    =-rm1(g)
                  node_r%source1(g)    =-rmp1(g)
                  node_l%source2(g)    =-rm2(g)
                  node_r%source2(g)    =-rmp2(g)
                  node_l%adf(g)        = df(g)
                  node_r%adf(g)        = dfp(g)
                  node_l%flux_avg(g)   = Flux(l,k,g)
                  node_r%flux_avg(g)   = Flux(lm,k,g)
                  node_l%diff_coeff(g) = D_3D(l,k,g)
                  node_r%diff_coeff(g) = D_3D(lm,k,g)
                  node_l%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                  node_r%sigma_r(g)    = maXs_r_3D(lm,k,g)+maXS_scat_3D(g,g,lm,k)
                  node_l%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                  node_r%nu_sigma_f(g) = nu_maXs_f_3D(lm,k,g)
                  do gg=1,n_group
                     node_l%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                     node_r%sigma_s(g,gg) = maXS_scat_3D(gg,g,lm,k)
                  enddo
               enddo
               select case (opt_nodal)
               case (1)
                  call solve_two_node_unm(node_l,node_r,keff,OPT_NEM)
               case (2)
                  call solve_two_node_unm(node_l,node_r,keff,OPT_ANM)
               case (3)
                  call solve_two_node_unm(node_l,node_r,keff,OPT_ANM)
               case default
                  call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                  stop
               end select
               do m=1,n_group
                  dnn(m,lm,k)=-(dfn(m,lm,k)*(Flux(lm,k,m)-Flux(l,k,m))+face_curr(m)) &
                             &             /(Flux(lm,k,m)+Flux(l, k, m))
               enddo
            enddo

            ! solve one-node problem for southern b.c. of row j if j(edge)<>0 and
            ! update dnw(edge) using dfw(edge) which exists.
            if (idomy==ndomy) then
               j=Iy_End_x(i)
               l=nodel(i,j)
               lp=nodel(i,j+1)
               select case (BC_Ry)
               case (1)
                  dnn(1,lp,k)=0d0
                  dnn(2,lp,k)=0d0
               case (2,3)
                  do m=1,n_group
                     rl(m)=L_0x_3D(m,l,k)+L_0z_3D(m,l,k)-src_tr(m,l,k)
                     rm1(m)=L_1y_3D(m,l,k)
                     rm2(m)=L_2y_3D(m,l,k)
                     rm1(m)=rm1(m)*0.5d0
                     rm2(m)=rm2(m)*0.5d0
                  enddo
                  node_l%width = MeshSize_y(j)
                  do g=1,n_group
                     node_l%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                     node_l%source0(g)    =-rl(g)
                     node_l%source1(g)    =-rm1(g)
                     node_l%source2(g)    =-rm2(g)
                     node_l%adf(g)        = ADF_Ry(l,k,g)
                     node_l%flux_avg(g)   = Flux(l,k,g)
                     node_l%diff_coeff(g) = D_3D(l,k,g)
                     node_l%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                     node_l%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                     do gg=1,n_group
                        node_l%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                     enddo
                  enddo
                  select case (opt_nodal)
                  case (1)
                     call solve_right_bndry_node_unm(BC_Ry,node_l,keff,OPT_NEM)
                  case (2)
                     call solve_right_bndry_node_unm(BC_Ry,node_l,keff,OPT_ANM)
                  case (3)
                     call solve_right_bndry_node_unm(BC_Ry,node_l,keff,OPT_ANM)
                  case default
                     call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                     stop
                  end select
                  do m=1,n_group
                     dnn(m,lp,k)=-face_curr(m)/Flux(l,k,m)+dfn(m,lp,k)
                  enddo
               end select
            endif
         enddo
      enddo
#ifdef tuan_tr_test

#else
      call update_axialdf
#endif      
      ! sweep along nodes in axial plane.

      ! solve one-node problem for bottom b.c. (k=1) if j(edge)<>0
      ! and update dll(edge) using dfb(edge) which exists.
      do l=1,Nxy
         if (idomz==1) then
            k=1
            select case (BC_Lz)
            case (1)
               do m=1,n_group
                  dnb(m,l,k)=0d0
                  if (flag_axialdf) phisurf_inpz(m,l,k)=flux(l,k,m)
               enddo

            case (2,3)
               do m=1,n_group
                  rl(m)=L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
                  rm1(m)=L_1z_3D(m,l,k)
                  rm2(m)=L_2z_3D(m,l,k)
                  rm1(m)=rm1(m)*0.5d0
                  rm2(m)=rm2(m)*0.5d0
               enddo
               node_r%width = MeshSize_z(k)
               do g=1,n_group
                  node_r%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                  node_r%source0(g)    =-rl(g)
                  node_r%source1(g)    =-rm1(g)
                  node_r%source2(g)    =-rm2(g)
                  node_r%adf(g)        = ADF_Lz(l,k,g)*corr_ADF_Lz(l,k,g)
                  node_r%flux_avg(g)   = Flux(l,k,g)
                  node_r%diff_coeff(g) = D_3D(l,k,g)
                  node_r%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                  node_r%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                  do gg=1,n_group
                     node_r%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                  enddo
               enddo
               select case (opt_nodal)
               case (1)
                  call solve_left_bndry_node_unm(BC_Lz,node_r,keff,OPT_NEM)
               case (2)
                  call solve_left_bndry_node_unm(BC_Lz,node_r,keff,OPT_ANM)
               case (3)
                  call solve_left_bndry_node_unm(BC_Lz,node_r,keff,OPT_ANM)
               case default
                  call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                  stop
               end select
               do m=1,n_group
                  dnb(m,l,k)=-dfb(m,l,k)-face_curr(m)/Flux(l,k,m)
                  if (flag_axialdf) phisurf_inpz(m,l,k)=face_flux(m)
               enddo
            end select
         endif

         ! now we solve for the (Nz-1) two-node problems in the z-direction
         do k=kbeg,ktop-1
            kl=k
            klp=k+1
            do m=1,n_group
               rl(m)   = L_0x_3D(m,l,kl)+L_0y_3D(m,l,kl)-src_tr(m,l,kl)
               rlp(m)  = L_0x_3D(m,l,klp)+L_0y_3D(m,l,klp)-src_tr(m,l,klp)
               rm1(m)  = L_1z_3D(m,l,kl)
               rmp1(m) = L_1z_3D(m,l,klp)
               rm2(m)  = L_2z_3D(m,l,kl)
               rmp2(m) = L_2z_3D(m,l,klp)
               df(m)   = ADF_Rz(l,kl,m)*corr_ADF_Rz(l,kl,m)
               dfp(m)  = ADF_Lz(l,klp,m)*corr_ADF_Lz(l,klp,m)
               rm1(m)  = rm1(m)*0.5d0
               rm2(m)  = rm2(m)*0.5d0
               rmp1(m) = rmp1(m)*0.5d0
               rmp2(m) = rmp2(m)*0.5d0
            enddo
            node_l%width = MeshSize_z(kl)
            node_r%width = MeshSize_z(klp)
            do g=1,n_group
               node_l%fiss_spec(g)  = maxs_chi_3d(l,kl,g)
               node_r%fiss_spec(g)  = maxs_chi_3d(l,klp,g)
               node_l%source0(g)    =-rl(g)
               node_r%source0(g)    =-rlp(g)
               node_l%source1(g)    =-rm1(g)
               node_r%source1(g)    =-rmp1(g)
               node_l%source2(g)    =-rm2(g)
               node_r%source2(g)    =-rmp2(g)
               node_l%adf(g)        = df(g)
               node_r%adf(g)        = dfp(g)
               node_l%flux_avg(g)   = Flux(l,kl,g)
               node_r%flux_avg(g)   = Flux(l,klp,g)
               node_l%diff_coeff(g) = D_3D(l,kl,g)
               node_r%diff_coeff(g) = D_3D(l,klp,g)
               node_l%sigma_r(g)    = maXs_r_3D(l,kl,g)+maXS_scat_3D(g,g,l,kl)
               node_r%sigma_r(g)    = maXs_r_3D(l,klp,g)+maXS_scat_3D(g,g,l,klp)
               node_l%nu_sigma_f(g) = nu_maXs_f_3D(l,kl,g)
               node_r%nu_sigma_f(g) = nu_maXs_f_3D(l,klp,g)
               do gg=1,n_group
                  node_l%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,kl)
                  node_r%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,klp)
               enddo
            enddo
            select case (opt_nodal)
            case (1)
               call solve_two_node_unm(node_l,node_r,keff,OPT_NEM)
            case (2)
               call solve_two_node_unm(node_l,node_r,keff,OPT_ANM)
            case (3)
               call solve_two_node_unm(node_l,node_r,keff,OPT_ANM)
            case default
               call print_msg(3,"updatedhat: opt_nodal is not assgiend")
               stop
            end select
            do m=1,n_group
               dnb(m,l,klp)=-(dfb(m,l,klp)*(Flux(l,klp,m)-Flux(l,kl,m))+face_curr(m)) &
                           &              /(Flux(l,klp,m)+Flux(l,kl,m))
               if (flag_axialdf) phisurf_inpz(m,l,klp)=face_flux(m)
            enddo
         enddo

         ! here we solve for the one-node problem for k=Nz only if j(edge)<>0
         ! and we update dlu(edge) using dcu(edge) which exists
         if (idomz==ndomz) then
            k=Nz
            kp=Nz+1
            select case (BC_Rz)
            case (1)
               do m=1,n_group
                  dnb(m,l,kp)=0d0
                  if (flag_axialdf) phisurf_inpz(m,l,kp)=flux(l,k,m)
               enddo
            case (2,3)
               do m=1,n_group
                  rl(m) =L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
                  rm1(m)=L_1z_3D(m,l,k)
                  rm2(m)=L_2z_3D(m,l,k)
                  rm1(m)=rm1(m)*0.5d0
                  rm2(m)=rm2(m)*0.5d0
               enddo
               node_l%width = MeshSize_z(k)
               do g=1,n_group
                  node_l%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                  node_l%source0(g)    =-rl(g)
                  node_l%source1(g)    =-rm1(g)
                  node_l%source2(g)    =-rm2(g)
                  node_l%adf(g)        = ADF_Rz(l,k,g)*corr_ADF_Rz(l,k,g)
                  node_l%flux_avg(g)   = Flux(l,k,g)
                  node_l%diff_coeff(g) = D_3D(l,k,g)
                  node_l%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                  node_l%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                  do gg=1,n_group
                     node_l%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                  enddo
               enddo
               select case (opt_nodal)
               case (1)
                  call solve_right_bndry_node_unm(BC_Rz,node_l,keff,OPT_NEM)
               case (2)
                  call solve_right_bndry_node_unm(BC_Rz,node_l,keff,OPT_ANM)
               case (3)
                  call solve_right_bndry_node_unm(BC_Rz,node_l,keff,OPT_ANM)
               case default
                  call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                  stop
               end select
               do m=1,n_group
                  dnb(m,l,kp)=-face_curr(m)/Flux(l,k,m)+dfb(m,l,kp)
                  if (flag_axialdf) phisurf_inpz(m,l,kp)=face_flux(m)
               enddo
            end select
         endif
      enddo

      if (Flag_Leakage) then
         call updcurnxy
         buck2_lc=0d0
         do k=1,Nz
            kp=k+1
            do i=1,Nx
               do j=1,Ny
                  l=nodel(i,j)
                  lp=nodel(i+1,j)
                  lm=nodel(i,j+1)
                  do m=1,n_group
                     buck2_lc(l,k,m) = buck2_lc(l,k,m) &
                       + (j_net_x_3D(m,lp,k)-j_net_x_3D(m,l,k)) / MeshSize_x(i) &
                       + (j_net_y_3D(m,lm,k)-j_net_y_3D(m,l,k)) / MeshSize_y(j) &
                       + (j_net_z_3D(m,l,kp)-j_net_z_3D(m,l,k)) / MeshSize_z(k)
                  enddo
               enddo
            enddo
         enddo
         do Ixy=1,Nxy
            do Iz=1,Nz
               if (Flag_LC_mean) then
                  I_Tab=0
               else
                  Ixy_1N=I_4Nto1N(Ixy)
                  I_Tab=AxialComp( I_LP_1N(Ixy_1N), Iz)
                  I_Tab=new_asym_itab(I_Tab,Ixy)
               endif
               do m=1,n_group
                  dbuck2_lc(Ixy,Iz,m) = buck2_lc(Ixy,Iz,m) / Flux (Ixy,Iz,m)
               enddo
               leak_ratio(Ixy,Iz,1) = (1 + Leakage_Table(1,I_Tab) * dbuck2_lc(Ixy,Iz,1) &
                          / (maXS_a_3D(Ixy,Iz,1)+maXS_s_3D(Ixy,Iz,1)) )
               leak_ratio(Ixy,Iz,2) = (1 + Leakage_Table(3,I_Tab) * dbuck2_lc(Ixy,Iz,1) &
                          / (maXS_a_3D(Ixy,Iz,1)+maXS_s_3D(Ixy,Iz,1)) )
               leak_ratio(Ixy,Iz,3) = (1 + Leakage_Table(2,I_Tab) * dbuck2_lc(Ixy,Iz,1) &
                          / (maXS_a_3D(Ixy,Iz,1)+maXS_s_3D(Ixy,Iz,1)) )
               leak_ratio(Ixy,Iz,4) = (1 + Leakage_Table(4,I_Tab) * dbuck2_lc(Ixy,Iz,2) &
                          / (maXS_a_3D(Ixy,Iz,2)+maXS_s_3D(Ixy,Iz,2)) )
            enddo
         enddo
      endif

      call deallocate_one_dim_node(node_r)
      call deallocate_one_dim_node(node_l)
      call finalize_two_node_unm_solver()

      return
      end subroutine updatedhat
#endif

      subroutine updcurnxy
      implicit none
      integer :: i,j,k,l,m,ktop,kbot
      integer :: lm,lp,kp,km


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [updcurnxy] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      kbot=1
      ktop=nz

      do k=kbot,ktop
         do j=1,ny
            do i=ix_start_y(j)+1,ix_end_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               do m=1,n_group
                  J_net_x_3D(m,l,k)=-dfw(m,l,k)*(flux(l,k,m)-flux(lm,k,m)) &
                                  & -dnw(m,l,k)*(flux(l,k,m)+flux(lm,k,m))
               enddo
            enddo
            if (idomx==1) then
               i=ix_start_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               if (BC_Lx==1) then
                  do m=1,n_group
                     j_net_x_3D(m,l,k)=0d0
                  enddo
               elseif (BC_Lx== 3.or.BC_Lx==2) then
                  do m=1,n_group
                     j_net_x_3D(m,l,k)=-(dfw(m,l,k)+dnw(m,l,k))*flux(l,k,m)
                  enddo
               endif
            endif
            if (idomx==ndomx) then
               i=ix_end_y(j)
               l=nodel(i,j)
               lp=nodel(i+1,j)
               if (BC_Rx==1) then
                  do m=1,n_group
                     j_net_x_3D(m,lp,k)=0
                  enddo
               elseif (BC_Rx==3.or.BC_Rx==2) then
                  do m=1,n_group
                     j_net_x_3D(m,lp,k)=(dfw(m,lp,k)-dnw(m,lp,k))*flux(l,k,m)
                  enddo
               endif
            endif
         enddo

         do i=1,nx
            do j=iy_start_x(i)+1,iy_end_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               do m=1,n_group
                  j_net_y_3D(m,l,k)=-dfn(m,l,k)*(flux(l,k,m)-flux(lm,k,m)) &
                                  & -dnn(m,l,k)*(flux(l,k,m)+flux(lm,k,m))
               enddo
            enddo
            if (idomy==1) then
               j=iy_start_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               if (BC_Ly==1) then
                  do m=1,n_group
                     j_net_y_3D(m,l,k)=0d0
                  enddo
               elseif (BC_Ly==3.or.BC_Ly==2) then
                  do m=1,n_group
                     j_net_y_3D(m,l,k)=-(dfn(m,l,k)+dnn(m,l,k))*flux(l,k,m)
                  enddo
               endif
            endif
            if (idomy==ndomy) then
               j=iy_end_x(i)
               l=nodel(i,j)
               lp=nodel(i,j+1)
               if (BC_Ry==1) then
                  do m=1,n_group
                     j_net_y_3D(m,lp,k)=0d0
                  enddo
               elseif (BC_Ry==3.or.BC_Ry==2) then
                  do m=1,n_group
                     j_net_y_3D(m,lp,k)=(dfn(m,lp,k)-dnn(m,lp,k))*flux(l,k,m)
                  enddo
               endif
            endif

         enddo
      enddo

      do l=1,nxy
         if (idomz==1) then
            k=1
            if (BC_Lz==1) then
               do m=1,n_group
                  j_net_z_3D(m,l,k)=0d0
               enddo
            elseif (BC_Lz==3.or.BC_Lz==2) then
               do m=1,n_group
                  j_net_z_3D(m,l,k)=-(dfb(m,l,k)+dnb(m,l,k))*flux(l,k,m)
               enddo
            endif
         endif
         if (idomz==ndomz) then
            k=nz
            kp=k+1
            if (BC_Rz==1) then
               do m=1,n_group
                  j_net_z_3D(m,l,kp)=0d0
               enddo
            elseif (BC_Rz==3.or.BC_Rz==2) then
               do m=1,n_group
                  j_net_z_3D(m,l,kp)=(dfb(m,l,kp)-dnb(m,l,kp))*flux(l,k,m)
               enddo
            endif
         endif
         do k=2,nz
            km=k-1
            do m=1,n_group
               j_net_z_3D(m,l,k)=-dfb(m,l,k)*(flux(l,k,m)-flux(l,km,m)) &
                               & -dnb(m,l,k)*(flux(l,k,m)+flux(l,km,m))
            enddo
         enddo
      enddo

      return
      end subroutine updcurnxy


      subroutine updtlxy
      use inc_solver, only: rvdelt, betap
      use inc_extsrc, only: flag_extsrc, extsrc
#ifdef js_sp3
      use inc_sp3, only: l2_2x_3d
      use inc_sp3, only: l2_2y_3d
      use inc_sp3, only: l2_2z_3d
      use inc_sp3, only: j2_net_x_3d
      use inc_sp3, only: j2_net_y_3d
      use inc_sp3, only: j2_net_z_3d
#endif
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: ixy2iproc, ixyip2ixy
      use inc_parallel, only: nodelmpi
#endif
      implicit none
      integer :: i,j,k,l,m,le,ls,lm,lp
      integer :: km1,kp1,jp1,isub,ktop,kbot
      real(8) :: w,wp,wm,rhx,rhy,rhz,rvol
      real(8) :: tl(n_group),tlp(n_group),tlm(n_group)
      integer :: l0, l00, le0, ls0, lm0, lp0


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [updtlxy] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      isub=1
      kbot=1
      ktop=nz

      do k=kbot,ktop
         do j=1,ny
            jp1=j+1
            rhy=1d0/MeshSize_y(j)
            do i=Ix_Start_y(j),Ix_End_y(j)
               l=nodel(i,j)
               rvol=1d0/NodeVolume(l,k)
               rhx=1d0/MeshSize_x(i)
               le=nodel(i+1,j)
               ls=nodel(i,jp1)
               l0=l
               l00=l
               le0=le
               ls0=ls
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  le0=nodelmpi(i+1,j)
                  ls0=nodelmpi(i,jp1)
               endif
#endif
               do m=1,n_group
                  if (flag_transient) then
                     if (n_group>2) then
                        src_tr(m,l00,k)=src(m,l,k)*rvol-rvdelt(m,l0,k)*flux(l,k,m) &
                        & -maxs_chi_3d(l,k,m)*(1d0-betap(l0,k))*FisSrc(l,k)*rvol
                     else
                        if (m==1) then
                           src_tr(1,l00,k)=src(1,l00,k)*rvol-rvdelt(1,l0,k)*flux(l,k,1) &
                           & -(1d0-betap(l0,k))*FisSrc(l,k)*rvol
                        else
                           src_tr(2,l00,k)=src(2,l,k)*rvol-rvdelt(2,l0,k)*flux(l,k,2)
                        endif
                     endif
                  else
                     src_tr(m,l00,k)=0d0
                     if (m==1) then
                        if (flag_extsrc) then
                           src_tr(1,l00,k)=extsrc(l,k)*rvol
                        endif
                     endif
                  endif
                  L_0x_3D(m,l00,k)=(j_net_x_3D(m,le0,k)-j_net_x_3D(m,l00,k))*rhx
                  L_0y_3D(m,l00,k)=(j_net_y_3D(m,ls0,k)-j_net_y_3D(m,l00,k))*rhy
               enddo
            enddo
         enddo
      enddo

#ifdef js_sp3
      if (opt_nodal==3) then
         do k=kbot,ktop
            do j=1,ny
               jp1=j+1
               rhy=1d0/MeshSize_y(j)
               do i=Ix_Start_y(j),Ix_End_y(j)
                  l=nodel(i,j)
                  rvol=1d0/NodeVolume(l,k)
                  rhx=1d0/MeshSize_x(i)
                  le=nodel(i+1,j)
                  ls=nodel(i,jp1)
                  l0=l
                  l00=l
                  le0=le
                  ls0=ls
#ifdef js_mpi
                  if (comm%usempi) then
                     if (iproc/=ixy2iproc(l)) cycle
                     l0=ixyip2ixy(l,iproc)
                     l00=nodelmpi(i,j)
                     le0=nodelmpi(i+1,j)
                     ls0=nodelmpi(i,jp1)
                  endif
#endif
                  do m=1,n_group
                     L2_2x_3D(m,l00,k)=(j2_net_x_3D(m,le0,k)-j2_net_x_3D(m,l00,k))*rhx
                     L2_2y_3D(m,l00,k)=(j2_net_y_3D(m,ls0,k)-j2_net_y_3D(m,l00,k))*rhy
                  enddo
               enddo
            enddo
         enddo
         do k=kbot,ktop
            rhz=1d0/MeshSize_z(k)
            kp1=k+1
            do l=1,nxy
               l0=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
               endif
#endif
               do m=1,n_group
                  L2_2z_3D(m,l0,k)=(j2_net_z_3D(m,l0,kp1)-j2_net_z_3D(m,l0,k))*rhz
               enddo
            enddo
         enddo
         ! update source term due to the leakage of 2nd moment
         do k=1,nz
            kp1=k+1
            km1=k-1
            if (kp1>nz) kp1=nz
            if (km1<1) km1=1
            do l=1,nxy
               l0=l
               l00=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(ltox(l),ltoy(l))
               endif
#endif
               do m=1,n_group
                  src_tr(m,l00,k)=src_tr(m,l00,k) &
                  & -(L2_2x_3D(m,l00,k)+L2_2y_3D(m,l00,k)+L2_2z_3D(m,l0,k))
               enddo
            enddo
         enddo
      endif ! opt_nodal==3
#endif

      if (nz==1) return !return if 2d

      do k=kbot,ktop
         rhz=1d0/MeshSize_z(k)
         kp1=k+1
         do l=1,nxy
            l0=l
            l00=l
#ifdef js_mpi
            if (comm%usempi) then
               if (iproc/=ixy2iproc(l)) cycle
               l0=ixyip2ixy(l,iproc)
               l00=nodelmpi(ltox(l),ltoy(l))
            endif
#endif
            do m=1,n_group
               L_0z_3D(m,l00,k)=(j_net_z_3D(m,l0,kp1)-j_net_z_3D(m,l0,k))*rhz
            enddo
         enddo
      enddo

      do k=kbot,ktop
         do j=1,ny
            do i=Ix_Start_y(j),Ix_End_y(j)
               w=MeshSize_x(i)
               l=nodel(i,j)
               l0=l
               l00=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
               endif
#endif
               do m=1,n_group
                  tl(m)=L_0y_3D(m,l00,k)+L_0z_3D(m,l00,k)-src_tr(m,l00,k)
               enddo
               if (i>Ix_Start_y(j)) then
                  wm=MeshSize_x(i-1)
                  lm=nodel(i-1,j)
                  lm0=lm
#ifdef js_mpi
                  if (comm%usempi) then
                     lm0=nodelmpi(i-1,j)
                  endif
#endif
                  do m=1,n_group
                     tlm(m)=L_0y_3D(m,lm0,k)+L_0z_3D(m,lm0,k)-src_tr(m,lm0,k)
                  enddo
               elseif (BC_Lx==1) then
                  wm=w
                  do m=1,n_group
                     tlm(m)=tl(m)
                  enddo
               else
                  wm=0d0
                  do m=1,n_group
                     tlm(m)=0d0
                  enddo
               endif
               if (i<Ix_End_y(j)) then
                  wp=MeshSize_x(i+1)
                  lp=nodel(i+1,j)
                  lp0=lp
#ifdef js_mpi
                  if (comm%usempi) then
                     lp0=nodelmpi(i+1,j)
                  endif
#endif
                  do m=1,n_group
                     tlp(m)=L_0y_3D(m,lp0,k)+L_0z_3D(m,lp0,k)-src_tr(m,lp0,k)
                  enddo
               elseif (BC_Rx==1) then
                  wp=w
                  do m=1,n_group
                     tlp(m)=tl(m)
                  enddo
               else
                  wp=0d0
                  do m=1,n_group
                     tlp(m)=0d0
                  enddo
               endif
               do m=1,n_group
                  L_1x_3D(m,l00,k)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
                  L_2x_3D(m,l00,k)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
               enddo
            enddo
         enddo

         do i=1,nx
            do j=Iy_Start_x(i),Iy_End_x(i)
               w=MeshSize_y(j)
               l=nodel(i,j)
               l0=l
               l00=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
               endif
#endif
               do m=1,n_group
                  tl(m)=L_0x_3D(m,l00,k)+L_0z_3D(m,l00,k)-src_tr(m,l00,k)
               enddo
               if (j>Iy_Start_x(i)) then
                  wm=MeshSize_y(j-1)
                  lm=nodel(i,j-1)
                  lm0=lm
#ifdef js_mpi
                  if (comm%usempi) then
                     lm0=nodelmpi(i,j-1)
                  endif
#endif
                  do m=1,n_group
                     tlm(m)=L_0x_3D(m,lm0,k)+L_0z_3D(m,lm0,k)-src_tr(m,lm0,k)
                  enddo
               elseif (BC_Ly==1) then
                  wm=w
                  do m=1,n_group
                     tlm(m)=tl(m)
                  enddo
               else
                  wm=0d0
                  do m=1,n_group
                     tlm(m)=0d0
                  enddo
               endif
               if (j<Iy_End_x(i)) then
                  wp=MeshSize_y(j+1)
                  lp=nodel(i,j+1)
                  lp0=lp
#ifdef js_mpi
                  if (comm%usempi) then
                     lp0=nodelmpi(i,j+1)
                  endif
#endif
                  do m=1,n_group
                    tlp(m)=L_0x_3D(m,lp0,k)+L_0z_3D(m,lp0,k)-src_tr(m,lp0,k)
                  enddo
               elseif (BC_Ry==1) then
                  wp=w
                  do m=1,n_group
                    tlp(m)=tl(m)
                  enddo
               else
                  wp=0d0
                  do m=1,n_group
                    tlp(m)=0d0
                  enddo
               endif
               do m=1,n_group
                  L_1y_3D(m,l00,k)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
                  L_2y_3D(m,l00,k)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
               enddo
            enddo
         enddo
      enddo

      if (isub==1) then
         k=1
         kp1=k+1
         w=MeshSize_z(k)
         wp=MeshSize_z(kp1)
         do l=1,nxy
            l0=l
            l00=l
#ifdef js_mpi
            if (comm%usempi) then
               if (iproc/=ixy2iproc(l)) cycle
               l0=ixyip2ixy(l,iproc)
               l00=nodelmpi(ltox(l),ltoy(l))
            endif
#endif
            do m=1,n_group
               tl(m)=L_0x_3D(m,l00,k)+L_0y_3D(m,l00,k)-src_tr(m,l00,k)
               tlp(m)=L_0x_3D(m,l00,kp1)+L_0y_3D(m,l00,kp1)-src_tr(m,l00,kp1)
               if (BC_Lz==1) then
                  wm=w
                  tlm(m)=tl(m)
               else
                  wm=0d0
                  tlm(m)=0d0
               endif
               L_1z_3D(m,l0,k)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2z_3D(m,l0,k)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      endif

      if (isub==ndomz) then
         k=nz
         km1=k-1
         w=MeshSize_z(k)
         wm=MeshSize_z(km1)
         do l=1,nxy
            l0=l
            l00=l
#ifdef js_mpi
            if (comm%usempi) then
               if (iproc/=ixy2iproc(l)) cycle
               l0=ixyip2ixy(l,iproc)
               l00=nodelmpi(ltox(l),ltoy(l))
            endif
#endif
            do m=1,n_group
               tl(m)=L_0x_3D(m,l00,k)+L_0y_3D(m,l00,k)-src_tr(m,l00,k)
               tlm(m)=L_0x_3D(m,l00,km1)+L_0y_3D(m,l00,km1)-src_tr(m,l00,km1)
               if (BC_Rz==1) then
                  wp=w
                  tlp(m)=tl(m)
               else
                  wp=0d0
                  tlp(m)=0d0
               endif
               L_1z_3D(m,l0,k)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2z_3D(m,l0,k)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      endif

      do k=2,nz-1
         kp1=k+1
         km1=k-1
         w=MeshSize_z(k)
         wm=MeshSize_z(km1)
         wp=MeshSize_z(kp1)
         do l=1,nxy
            l0=l
            l00=l
#ifdef js_mpi
            if (comm%usempi) then
               if (iproc/=ixy2iproc(l)) cycle
               l0=ixyip2ixy(l,iproc)
               l00=nodelmpi(ltox(l),ltoy(l))
            endif
#endif
            do m=1,n_group
               tl(m)=L_0x_3D(m,l00,k)+L_0y_3D(m,l00,k)-src_tr(m,l00,k)
               tlp(m)=L_0x_3D(m,l00,kp1)+L_0y_3D(m,l00,kp1)-src_tr(m,l00,kp1)
               tlm(m)=L_0x_3D(m,l00,km1)+L_0y_3D(m,l00,km1)-src_tr(m,l00,km1)
               L_1z_3D(m,l0,k)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2z_3D(m,l0,k)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      enddo

      return
      end subroutine updtlxy


      FUNCTION pol1 (w,wp,wm,tl,tlp,tlm)
      IMPLICIT NONE
      REAL(8) :: pol1, w, wp, wm, tl, tlp, tlm, dd
      dd=one/((w+wp)*(w+wm)*(w+wm+wp))

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [pol1] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      pol1=dd*w*((tlp-tl)*(w+2*wm)*(w+wm)+(tl-tlm)*(w+2*wp)*(w+wp))
      RETURN
      END FUNCTION pol1

      FUNCTION pol2 (w,wp,wm,tl,tlp,tlm)
      IMPLICIT NONE
      REAL(8) :: pol2,w, wp, wm, tl, tlp, tlm, dd
      dd=one/((w+wp)*(w+wm)*(w+wm+wp))

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [pol2] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      pol2=dd*w*w*((tlp-tl)*(w+wm)+(tlm-tl)*(w+wp))
      RETURN
      END FUNCTION pol2


      subroutine init_subaxial_nodal
      implicit none
      real(8) :: hzori
      integer :: iz, izs, ize


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [init_subaxial_nodal] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.allocated(corr_ADF_Lz)) allocate(corr_ADF_Lz(nxy,nz,n_group))
      if (.not.allocated(corr_ADF_Rz)) allocate(corr_ADF_Rz(nxy,nz,n_group))
      corr_ADF_Lz=1d0
      corr_ADF_Rz=1d0

      if (.not.flag_axialdf) return
      if (nz==1) then
         flag_axialdf=.false.
         return
      endif

      ! geometry for sub-axial mesh (around 1cm)
      if (.not.allocated(nz_subz)) allocate(nz_subz(1:nz))
      nz_subz=1
      do iz=1,nz
         hzori=meshsize_z(iz)
         nz_subz(iz)=max(1,int(hzori))
      enddo
      nsubz=sum(nz_subz)
      if (.not.allocated(inp2sub_beg)) allocate(inp2sub_beg(1:nz))
      if (.not.allocated(inp2sub_end)) allocate(inp2sub_end(1:nz))
      inp2sub_beg=0
      inp2sub_end=0
      inp2sub_beg(1)=1
      inp2sub_end(1)=nz_subz(1)
      do iz=2,nz
         inp2sub_beg(iz)=inp2sub_end(iz-1)+1
         inp2sub_end(iz)=inp2sub_end(iz-1)+nz_subz(iz)
      enddo
      if (.not.allocated(sub2inp)) allocate(sub2inp(1:nsubz))
      sub2inp=0
      izs=1
      ize=1
      do iz=1,nz
         izs=inp2sub_beg(iz)
         ize=inp2sub_end(iz)
         sub2inp(izs:ize)=iz
      enddo
      if (.not.allocated(hsubz)) allocate(hsubz(1:nsubz))
      hsubz=0d0
      izs=1
      ize=1
      do iz=1,nz
         izs=inp2sub_beg(iz)
         ize=inp2sub_end(iz)
         hsubz(izs:ize)=meshsize_z(iz)/real(nz_subz(iz),8)
      enddo

      if (.not.allocated(L_1subz_3D)) allocate(L_1subz_3D(1:n_group,1:nxy,0:nsubz))
      if (.not.allocated(L_2subz_3D)) allocate(L_2subz_3D(1:n_group,1:nxy,0:nsubz))
      L_1subz_3D=0d0
      L_2subz_3D=0d0
      if (.not.allocated(phisurf_inpz)) allocate(phisurf_inpz(n_group,nxy,nz+1))
      if (.not.allocated(phisurf_subz)) allocate(phisurf_subz(n_group,nxy,nsubz+1))
      phisurf_inpz=0d0
      phisurf_subz=0d0

      return
      end subroutine init_subaxial_nodal

#ifdef siarhei_delete

      subroutine update_axialdf
      implicit none
      type(UNMNodal) node_l, node_r
      real(8) :: rl(n_group)
      real(8) :: rlp(n_group)
      real(8) :: rm1(n_group)
      real(8) :: rmp1(n_group)
      real(8) :: rm2(n_group)
      real(8) :: rmp2(n_group)
      real(8) :: df(n_group)
      real(8) :: dfp(n_group)
      real(8) :: tl(n_group)
      real(8) :: tlp(n_group)
      real(8) :: tlm(n_group)
      integer :: g, gg
      integer :: l, m
      integer :: izs, ize
      integer :: k, kp, km
      integer :: subk, subkp, subkm
      real(8) :: w, wp, wm


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [update_axialdf] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.flag_axialdf) return
      if (sum(phisurf_inpz(:,:,:))<1d-10) return

      ! update sub-axial mesh transverse leakage
      if (idomz==1) then
         subk=1
         subkp=subk+1
         w=hsubz(subk)
         wp=hsubz(subkp)
         k=sub2inp(subk)
         kp=sub2inp(subkp)
         do l=1,nxy
            do m=1,n_group
               tl(m)=L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
               tlp(m)=L_0x_3D(m,l,kp)+L_0y_3D(m,l,kp)-src_tr(m,l,kp)
               if (BC_Lz==1) then
                  wm=w
                  tlm(m)=tl(m)
               else
                  wm=0d0
                  tlm(m)=0d0
               endif
               L_1subz_3D(m,l,subk)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2subz_3D(m,l,subk)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      endif
      do subk=2,nsubz-1
         subkp=subk+1
         subkm=subk-1
         k=sub2inp(subk)
         kp=sub2inp(subkp)
         km=sub2inp(subkm)
         w=hsubz(subk)
         wm=hsubz(subkm)
         wp=hsubz(subkp)
         do l=1,nxy
            do m=1,n_group
               tl(m)=L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
               tlp(m)=L_0x_3D(m,l,kp)+L_0y_3D(m,l,kp)-src_tr(m,l,kp)
               tlm(m)=L_0x_3D(m,l,km)+L_0y_3D(m,l,km)-src_tr(m,l,km)
               L_1subz_3D(m,l,subk)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2subz_3D(m,l,subk)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      enddo
      if (idomz==ndomz) then
         subk=nsubz
         subkm=subk-1
         k=sub2inp(subk)
         km=sub2inp(subkm)
         w=hsubz(subk)
         wm=hsubz(subkm)
         do l=1,nxy
            do m=1,n_group
               tl(m)=L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
               tlm(m)=L_0x_3D(m,l,km)+L_0y_3D(m,l,km)-src_tr(m,l,km)
               if (BC_Rz==1) then
                  wp=w
                  tlp(m)=tl(m)
               else
                  wp=0d0
                  tlp(m)=0d0
               endif
               L_1subz_3D(m,l,subk)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2subz_3D(m,l,subk)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      endif

      ! solve 1d nodal for sub-axial mesh
      call Allocate_one_dim_node(node_l,n_group)
      call Allocate_one_dim_node(node_r,n_group)
      do l=1,Nxy
         if (idomz==1) then
            subk=1
            k=sub2inp(subk)
            select case (BC_Lz)
            case (1)
               phisurf_subz(:,l,subk)=flux(l,k,:)
            case (2,3)
               do m=1,n_group
                  rl(m)=L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
                  rm1(m)=L_1subz_3D(m,l,subk) !re-calculate line 1503~1589
                  rm2(m)=L_2subz_3D(m,l,subk) !same... search MeshSize_z
                  rm1(m)=rm1(m)*0.5d0
                  rm2(m)=rm2(m)*0.5d0
               enddo
               node_r%width = hsubz(subk)
               do g=1,n_group
                  node_r%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                  node_r%source0(g)    =-rl(g)
                  node_r%source1(g)    =-rm1(g)
                  node_r%source2(g)    =-rm2(g)
                  node_r%adf(g)        = 1d0
                  node_r%flux_avg(g)   = Flux(l,k,g)
                  node_r%diff_coeff(g) = D_3D(l,k,g)
                  node_r%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                  node_r%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                  do gg=1,n_group
                     node_r%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                  enddo
               enddo
               select case (opt_nodal)
               case (1)
                  call solve_left_bndry_node_unm(BC_Lz,node_r,keff,OPT_NEM)
               case (2)
                  call solve_left_bndry_node_unm(BC_Lz,node_r,keff,OPT_ANM)
               case (3)
                  call solve_left_bndry_node_unm(BC_Lz,node_r,keff,OPT_ANM)
               case default
                  call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                  stop
               end select
               do m=1,n_group
                  phisurf_subz(m,l,subk)=face_flux(m)
               enddo
            end select
         endif

         do subk=1,nsubz-1
            subkp=subk+1
            k=sub2inp(subk)
            kp=sub2inp(subkp)
            do m=1,n_group
               rl(m)   = L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
               rlp(m)  = L_0x_3D(m,l,kp)+L_0y_3D(m,l,kp)-src_tr(m,l,kp)
               rm1(m)  = L_1subz_3D(m,l,subk)
               rmp1(m) = L_1subz_3D(m,l,subkp)
               rm2(m)  = L_2subz_3D(m,l,subk)
               rmp2(m) = L_2subz_3D(m,l,subkp)
               df(m)   = ADF_Rz(l,k,m)
               dfp(m)  = ADF_Lz(l,kp,m)
               rm1(m)  = rm1(m)*0.5d0
               rm2(m)  = rm2(m)*0.5d0
               rmp1(m) = rmp1(m)*0.5d0
               rmp2(m) = rmp2(m)*0.5d0
            enddo
            node_l%width = hsubz(subk)
            node_r%width = hsubz(subkp)
            do g=1,n_group
               node_l%fiss_spec(g)  = maxs_chi_3d(l,k,g)
               node_r%fiss_spec(g)  = maxs_chi_3d(l,kp,g)
               node_l%source0(g)    =-rl(g)
               node_r%source0(g)    =-rlp(g)
               node_l%source1(g)    =-rm1(g)
               node_r%source1(g)    =-rmp1(g)
               node_l%source2(g)    =-rm2(g)
               node_r%source2(g)    =-rmp2(g)
               node_l%adf(g)        = df(g)
               node_r%adf(g)        = dfp(g)
               node_l%flux_avg(g)   = Flux(l,k,g)
               node_r%flux_avg(g)   = Flux(l,kp,g)
               node_l%diff_coeff(g) = D_3D(l,k,g)
               node_r%diff_coeff(g) = D_3D(l,kp,g)
               node_l%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
               node_r%sigma_r(g)    = maXs_r_3D(l,kp,g)+maXS_scat_3D(g,g,l,kp)
               node_l%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
               node_r%nu_sigma_f(g) = nu_maXs_f_3D(l,kp,g)
               do gg=1,n_group
                  node_l%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                  node_r%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,kp)
               enddo
            enddo
            select case (opt_nodal)
            case (1)
               call solve_two_node_unm(node_l,node_r,keff,OPT_NEM)
            case (2)
               call solve_two_node_unm(node_l,node_r,keff,OPT_ANM)
            case (3)
               call solve_two_node_unm(node_l,node_r,keff,OPT_ANM)
            case default
               call print_msg(3,"updatedhat: opt_nodal is not assgiend")
               stop
            end select
            do m=1,n_group
               phisurf_subz(m,l,subkp)=face_flux(m)
            enddo
         enddo

         if (idomz==ndomz) then
            subk=nsubz
            subkp=nsubz+1
            k=sub2inp(subk)
            kp=k+1
            select case (BC_Rz)
            case (1)
               phisurf_subz(:,l,subkp)=flux(l,k,:)
            case (2,3)
               do m=1,n_group
                  rl(m) =L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
                  rm1(m)=L_1subz_3D(m,l,subk)
                  rm2(m)=L_2subz_3D(m,l,subk)
                  rm1(m)=rm1(m)*0.5d0
                  rm2(m)=rm2(m)*0.5d0
               enddo
               node_l%width = hsubz(subk)
               do g=1,n_group
                  node_l%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                  node_l%source0(g)    =-rl(g)
                  node_l%source1(g)    =-rm1(g)
                  node_l%source2(g)    =-rm2(g)
                  node_l%adf(g)        = ADF_Rz(l,k,g)
                  node_l%flux_avg(g)   = Flux(l,k,g)
                  node_l%diff_coeff(g) = D_3D(l,k,g)
                  node_l%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                  node_l%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                  do gg=1,n_group
                     node_l%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                  enddo
               enddo
               select case (opt_nodal)
               case (1)
                  call solve_right_bndry_node_unm(BC_Rz,node_l,keff,OPT_NEM)
               case (2)
                  call solve_right_bndry_node_unm(BC_Rz,node_l,keff,OPT_ANM)
               case (3)
                  call solve_right_bndry_node_unm(BC_Rz,node_l,keff,OPT_ANM)
               case default
                  call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                  stop
               end select
               do m=1,n_group
                  phisurf_subz(m,l,subkp)=face_flux(m)
               enddo
            end select
         endif
      enddo

      call deallocate_one_dim_node(node_r)
      call deallocate_one_dim_node(node_l)

      ! calculate correction factor of  axial DF
      do l=1,nxy
         do k=1,nz
            izs=inp2sub_beg(k)
            ize=inp2sub_end(k)
            do m=1,n_group
               corr_ADF_Lz(l,k,m)=phisurf_inpz(m,l,k)/phisurf_subz(m,l,izs)
               if (corr_ADF_Lz(l,k,m)>2d0.or.corr_ADF_Lz(l,k,m)<0d0) then
                  corr_ADF_Lz(l,k,m)=1d0
               else
                  corr_ADF_Lz(l,k,m)=corr_ADF_Lz(l,k,m)/ADF_Lz(l,k,m)
               endif
            enddo
         enddo
      enddo
      do l=1,nxy
         do k=1,nz
            izs=inp2sub_beg(k)
            ize=inp2sub_end(k)
            do m=1,n_group
               corr_ADF_Rz(l,k,m)=phisurf_inpz(m,l,k+1)/phisurf_subz(m,l,ize+1)
               if (corr_ADF_Rz(l,k,m)>2d0.or.corr_ADF_Rz(l,k,m)<0d0) then
                  corr_ADF_Rz(l,k,m)=1d0
               else
                  corr_ADF_Rz(l,k,m)=corr_ADF_Rz(l,k,m)/ADF_rz(l,k,m)
               endif
            enddo
         enddo
      enddo

      return
      end subroutine update_axialdf
#endif
      
#ifdef js_sp3
      subroutine updatedhat2d1d
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: ixy2iproc, ixyip2ixy
      use inc_parallel, only: nodelmpi
      use mod_parallel, only: barrier
#endif
      implicit none
      type(UNMNodal) node_l, node_r
      real(8) :: rl(n_group), rlp(n_group)
      real(8) :: rm1(n_group), rmp1(n_group)
      real(8) :: rm2(n_group), rmp2(n_group)
      real(8) :: df(n_group), dfp(n_group)
      integer :: g, gg, isub, kp, kl, klp
      integer :: kbot, ktop, l, k, m
      integer :: l0, l00


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [updatedhat2d1d] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call initialize_two_node_unm_solver( n_group )
      call allocate_one_dim_node( node_l, n_group )
      call allocate_one_dim_node( node_r, n_group )
      ! call updtlxy2d1d => updcurnxy, updtlxy, updtlcoeff2
      call updcurnxy_2d1d
      call updtlxy

      isub=1
      kbot=1
      ktop=nz

      ! sweep along nodes in axial plane.

      ! update dnb(edge) using dfb(edge) which exists.
      do l=1,nxy
         l0=l
         l00=l
#ifdef js_mpi
         if (comm%usempi) then
            if (iproc/=ixy2iproc(l)) cycle
            l0=ixyip2ixy(l,iproc)
            l00=nodelmpi(ltox(l),ltoy(l))
         endif
#endif
         ! solve one-node problem for bottom b.c. (k=1)
         k=kbot
         if (BC_Lz==1) then ! reflective
            do m=1,n_group
               dnb(m,l0,k)=0d0
            enddo
         elseif (BC_Lz==2.or.BC_Lz==3) then ! Black / Fluxzero
            do m=1,n_group
               rl(m) =L_0x_3D(m,l00,k)+L_0y_3D(m,l00,k)-src_tr(m,l00,k)
               rm1(m)=L_1z_3D(m,l0,k)
               rm2(m)=L_2z_3D(m,l0,k)
               rm1(m)=rm1(m)*0.5d0
               rm2(m)=rm2(m)*0.5d0
            enddo
            node_r%width=meshsize_z(k)
            do g=1,n_group
               node_r%fiss_spec(g)  = maxs_chi_3d(l,k,g)
               node_r%source0(g)    =-rl(g)
               node_r%source1(g)    =-rm1(g)
               node_r%source2(g)    =-rm2(g)
               node_r%adf(g)        = 1d0 !js+ADF_Lz(l,k,g)
               node_r%flux_avg(g)   = flux(l,k,g)
               node_r%diff_coeff(g) = d_3d(l,k,g)
               node_r%sigma_r(g)    = maxs_r_3d(l,k,g)+maxs_scat_3d(g,g,l,k)
               node_r%nu_sigma_f(g) = nu_maxs_f_3d(l,k,g)
               do gg=1,n_group
                  node_r%sigma_s(g,gg) = maxs_scat_3d(gg,g,l,k)
               enddo
            enddo
            call solve_left_bndry_node_unm( BC_Lz, node_r, keff, OPT_ANM )
            do m=1,n_group
               dnb(m,l0,k)=-dfb(m,l0,k)-face_curr(m)/flux(l,k,m)
            enddo
         endif

         ! solve for the (nz-1) two-node problems in the z-direction
         do k=kbot,ktop-1
            kl=k
            klp=k+1
            do m=1,n_group
               rl(m)   = L_0x_3D(m,l00,kl) +L_0y_3D(m,l00,kl) -src_tr(m,l00,kl)
               rlp(m)  = L_0x_3D(m,l00,klp)+L_0y_3D(m,l00,klp)-src_tr(m,l00,klp)
               rm1(m)  = L_1z_3D(m,l0,kl)
               rmp1(m) = L_1z_3D(m,l0,klp)
               rm2(m)  = L_2z_3D(m,l0,kl)
               rmp2(m) = L_2z_3D(m,l0,klp)
               df(m)   = 1d0 !js+ADF_Rz(l,kl,m)
               dfp(m)  = 1d0 !js+ADF_Lz(l,klp,m)
               rm1(m)  = rm1(m)*0.5d0
               rmp1(m) = rmp1(m)*0.5d0
               rm2(m)  = rm2(m)*0.5d0
               rmp2(m) = rmp2(m)*0.5d0
            enddo
            node_l%width = meshsize_z(kl)
            node_r%width = meshsize_z(klp)
            do g=1,n_group
               node_l%fiss_spec(g)  = maxs_chi_3d(l,kl,g)
               node_r%fiss_spec(g)  = maxs_chi_3d(l,klp,g)
               node_l%source0(g)    = -rl(g)
               node_r%source0(g)    = -rlp(g)
               node_l%source1(g)    = -rm1(g)
               node_r%source1(g)    = -rmp1(g)
               node_l%source2(g)    = -rm2(g)
               node_r%source2(g)    = -rmp2(g)
               node_l%adf(g)        = df(g)
               node_r%adf(g)        = dfp(g)
               node_l%flux_avg(g)   = flux(l,kl,g)
               node_r%flux_avg(g)   = flux(l,klp,g)
               node_l%diff_coeff(g) = d_3d(l,kl,g)
               node_r%diff_coeff(g) = d_3d(l,klp,g)
               node_l%sigma_r(g)    = maxs_r_3d(l,kl,g) +maxs_scat_3d(g,g,l,kl)
               node_r%sigma_r(g)    = maxs_r_3d(l,klp,g)+maxs_scat_3d(g,g,l,klp)
               node_l%nu_sigma_f(g) = nu_maxs_f_3d(l,kl,g)
               node_r%nu_sigma_f(g) = nu_maxs_f_3d(l,klp,g)
               do gg=1,n_group
                  node_l%sigma_s(g,gg)=maxs_scat_3d(gg,g,l,kl)
                  node_r%sigma_s(g,gg)=maxs_scat_3d(gg,g,l,klp)
               enddo
            enddo
            call solve_two_node_unm( node_l, node_r, keff, OPT_ANM )
            do m=1,n_group
               dnb(m,l0,klp)=-(dfb(m,l0,klp)*(flux(l,klp,m)-flux(l,kl,m))+face_curr(m)) &
                              /(flux(l,klp,m)+flux(l,kl,m))
            enddo
         enddo

         ! solve for the one-node problem for k=nz only if j(edge)<>0
         k=nz
         kp=nz+1
         if (BC_Rz==1) then ! Refelctive
            do m=1,n_group
               dnb(m,l0,kp)=0d0
            enddo
         elseif (BC_Rz==2.or.BC_Rz==3) then ! Black / Flux zero
            do m=1,n_group
               rl(m) =L_0x_3D(m,l00,k)+L_0y_3D(m,l00,k)-src_tr(m,l00,k)
               rm1(m)=L_1z_3D(m,l0,k)
               rm2(m)=L_2z_3D(m,l0,k)
               rm1(m)=rm1(m)*0.5d0
               rm2(m)=rm2(m)*0.5d0
            enddo
            node_l%width=meshsize_z(k)
            do g=1,n_group
               node_l%fiss_spec(g)  = maxs_chi_3d(l,k,g)
               node_l%source0(g)    = -rl(g)
               node_l%source1(g)    = -rm1(g)
               node_l%source2(g)    = -rm2(g)
               node_l%adf(g)        = 1d0 !js+ADF_Rz(l,k,g)
               node_l%flux_avg(g)   = flux(l,k,g)
               node_l%diff_coeff(g) = d_3d(l,k,g)
               node_l%sigma_r(g)    = maxs_r_3d(l,k,g)+maxs_scat_3d(g,g,l,k)
               node_l%nu_sigma_f(g) = nu_maxs_f_3d(l,k,g)
               do gg=1,n_group
                  node_l%sigma_s(g,gg) = maxs_scat_3d(gg,g,l,k)
               enddo
            enddo
            call solve_right_bndry_node_unm ( BC_Rz, node_l, keff, OPT_ANM )
            do m=1,n_group
               dnb(m,l0,kp)=-face_curr(m)/flux(l,k,m)+dfb(m,l0,kp)
            enddo
         endif
      enddo

      call deallocate_one_dim_node( node_r )
      call deallocate_one_dim_node( node_l )
      call finalize_two_node_unm_solver( )

      !js+if (opt_nodal==3) then
      !js+   call upddhat2d1d2 ! not yet
      !js+endif

#ifdef js_mpi
      if (comm%usempi) call barrier
#endif

      return
      end subroutine updatedhat2d1d
#endif
#ifdef js_sp3
      subroutine updcurnxy_2d1d
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: ixy2iproc, ixyip2ixy
      use inc_parallel, only: nodelmpi
#endif
      implicit none
      integer :: i,j,k,l,m,ktop,kbot
      integer :: lm,lp,kp,km
      integer :: l0, l00, lm0, lp0


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [updcurnxy_2d1d] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      kbot=1
      ktop=nz

      do k=kbot,ktop
         do j=1,ny
            do i=ix_start_y(j)+1,ix_end_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i-1,j)
                  lp=nodel(i+1,j)
                  lp0=nodelmpi(i+1,j)
               endif
#endif
               do m=1,n_group
                  j_net_x_3D(m,l00,k)=-dfw(m,l00,k)*(flux(l,k,m)-flux(lm,k,m))
                  !J_net_x_3D(m,l00,k)=-dfw(m,l00,k)*(flux(l,k,m)-flux(lm,k,m)) &
                  !                  & -dnw(m,l00,k)*(flux(l,k,m)+flux(lm,k,m))
#ifdef js_mpi
                  if (comm%usempi) then
                     j_net_x_3D(m,lp0,k)=-dfw(m,lp0,k)*(flux(lp,k,m)-flux(l,k,m))
                  endif
#endif
               enddo
            enddo
            if (idomx==1) then
               i=ix_start_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 901
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i-1,j)
               endif
#endif
               if (BC_Lx==1) then
                  do m=1,n_group
                     j_net_x_3D(m,l00,k)=0d0
                  enddo
               elseif (BC_Lx== 3.or.BC_Lx==2) then
                  do m=1,n_group
                     j_net_x_3D(m,l00,k)=-dfw(m,l00,k)*flux(l,k,m)
                     !j_net_x_3D(m,l00,k)=-(dfw(m,l00,k)+dnw(m,l00,k))*flux(l,k,m)
                  enddo
               endif
#ifdef js_mpi
901            continue
#endif
            endif
            if (idomx==ndomx) then
               i=ix_end_y(j)
               l=nodel(i,j)
               lp=nodel(i+1,j)
               l0=l
               l00=l
               lp0=lp
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 902
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lp0=nodelmpi(i+1,j)
               endif
#endif
               if (BC_Rx==1) then
                  do m=1,n_group
                     j_net_x_3D(m,lp0,k)=0
                  enddo
               elseif (BC_Rx==3.or.BC_Rx==2) then
                  do m=1,n_group
                     j_net_x_3D(m,lp0,k)=dfw(m,lp0,k)*flux(l,k,m)
                     !j_net_x_3D(m,lp0,k)=(dfw(m,lp0,k)-dnw(m,lp0,k))*flux(l,k,m)
                  enddo
               endif
#ifdef js_mpi
902            continue
#endif
            endif
         enddo

         do i=1,nx
            do j=iy_start_x(i)+1,iy_end_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i,j-1)
                  lp=nodel(i,j+1)
                  lp0=nodelmpi(i,j+1)
               endif
#endif
               do m=1,n_group
                  j_net_y_3D(m,l00,k)=-dfn(m,l00,k)*(flux(l,k,m)-flux(lm,k,m))
                  !j_net_y_3D(m,l00,k)=-dfn(m,l00,k)*(flux(l,k,m)-flux(lm,k,m)) &
                  !                  & -dnn(m,l00,k)*(flux(l,k,m)+flux(lm,k,m))
#ifdef js_mpi
                  if (comm%usempi) then
                     j_net_y_3D(m,lp0,k)=-dfn(m,lp0,k)*(flux(lp,k,m)-flux(l,k,m))
                  endif
#endif
               enddo
            enddo
            if (idomy==1) then
               j=iy_start_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 903
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i,j-1)
               endif
#endif
               if (BC_Ly==1) then
                  do m=1,n_group
                     j_net_y_3D(m,l00,k)=0d0
                  enddo
               elseif (BC_Ly==3.or.BC_Ly==2) then
                  do m=1,n_group
                     j_net_y_3D(m,l00,k)=-dfn(m,l00,k)*flux(l,k,m)
                     !j_net_y_3D(m,l00,k)=-(dfn(m,l00,k)+dnn(m,l00,k))*flux(l,k,m)
                  enddo
               endif
#ifdef js_mpi
903            continue
#endif
            endif
            if (idomy==ndomy) then
               j=iy_end_x(i)
               l=nodel(i,j)
               lp=nodel(i,j+1)
               l0=l
               l00=l
               lp0=lp
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 904
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lp0=nodelmpi(i,j+1)
               endif
#endif
               if (BC_Ry==1) then
                  do m=1,n_group
                     j_net_y_3D(m,lp0,k)=0d0
                  enddo
               elseif (BC_Ry==3.or.BC_Ry==2) then
                  do m=1,n_group
                     j_net_y_3D(m,lp0,k)=dfn(m,lp0,k)*flux(l,k,m)
                     !j_net_y_3D(m,lp0,k)=(dfn(m,lp0,k)-dnn(m,lp0,k))*flux(l,k,m)
                  enddo
               endif
#ifdef js_mpi
904            continue
#endif
            endif

         enddo
      enddo

      do l=1,nxy
         l0=l
#ifdef js_mpi
         if (comm%usempi) then
            if (iproc/=ixy2iproc(l)) cycle
            l0=ixyip2ixy(l,iproc)
         endif
#endif
         if (idomz==1) then
            k=1
            if (BC_Lz==1) then
               do m=1,n_group
                  j_net_z_3D(m,l0,k)=0d0
               enddo
            elseif (BC_Lz==3.or.BC_Lz==2) then
               do m=1,n_group
                  j_net_z_3D(m,l0,k)=-(dfb(m,l0,k)+dnb(m,l0,k))*flux(l,k,m)
               enddo
            endif
         endif
         if (idomz==ndomz) then
            k=nz
            kp=k+1
            if (BC_Rz==1) then
               do m=1,n_group
                  j_net_z_3D(m,l0,kp)=0d0
               enddo
            elseif (BC_Rz==3.or.BC_Rz==2) then
               do m=1,n_group
                  j_net_z_3D(m,l0,kp)=(dfb(m,l0,kp)-dnb(m,l0,kp))*flux(l,k,m)
               enddo
            endif
         endif
         do k=2,nz
            km=k-1
            do m=1,n_group
               j_net_z_3D(m,l0,k)=-dfb(m,l0,k)*(flux(l,k,m)-flux(l,km,m)) &
                                & -dnb(m,l0,k)*(flux(l,k,m)+flux(l,km,m))
            enddo
         enddo
      enddo

#ifdef js_sp3
      if (opt_nodal==3) then
         call updcur2
      endif
#endif

      return
      end subroutine updcurnxy_2d1d
#endif
#ifdef js_sp3
      subroutine updcur2
      use inc_sp3, only: j2_net_x_3d, j2_net_y_3d, j2_net_z_3d
      !use inc_sp3, only: dnb2
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: ixy2iproc, ixyip2ixy
      use inc_parallel, only: nodelmpi
#endif
      implicit none
      integer :: i, j, k, m, l, lm, lp, kp, km1
      integer :: isub
      integer :: kend, kbeg, kbot, ktop
      integer :: l0, l00, lm0, lp0


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [updcur2] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      isub=1
      kbot=1
      ktop=nz
      do k=kbot,ktop
         do j=1,ny
            do i=Ix_Start_y(j)+1,Ix_End_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i-1,j)
                  lp=nodel(i+1,j)
                  lp0=nodelmpi(i+1,j)
               endif
#endif
               do m=1,n_group
                  j2_net_x_3D(m,l00,k)=-2d0*dfw(m,l00,k)*(flux2(l,k,m)-flux2(lm,k,m))
#ifdef js_mpi
                  if (comm%usempi) then
                     j2_net_x_3D(m,lp0,k)=-2d0*dfw(m,lp0,k)*(flux2(lp,k,m)-flux2(l,k,m))
                  endif
#endif
               enddo
            enddo
            if (idomx==1) then
               i=Ix_Start_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 901
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i-1,j)
               endif
#endif
               if (BC_Lx==1) then
                  do m=1,n_group
                     j2_net_x_3D(m,l00,k)=0d0
                  enddo
               elseif (BC_Lx==3.or.BC_Lx==2) then
                  do m=1,n_group
                     j2_net_x_3D(m,l00,k)=-2d0*dfw(m,l00,k)*flux2(l,k,m)
                  enddo
               endif
#ifdef js_mpi
901            continue
#endif
            endif
            if (idomx==ndomx) then
               i=Ix_End_y(j)
               l=nodel(i,j)
               lp=nodel(i+1,j)
               l0=l
               l00=l
               lp0=lp
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 902
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lp0=nodelmpi(i+1,j)
               endif
#endif
               if (BC_Rx==1) then
                  do m=1,n_group
                     j2_net_x_3D(m,lp0,k)=0d0
                  enddo
               elseif (BC_Rx==3.or.BC_Rx==2) then
                  do m=1,n_group
                     j2_net_x_3D(m,lp0,k)=2d0*dfw(m,lp0,k)*flux2(l,k,m)
                  enddo
               endif
#ifdef js_mpi
902            continue
#endif
            endif

         enddo

         do i=1,nx
            do j=Iy_Start_x(i)+1,Iy_End_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i,j-1)
                  lp=nodel(i,j+1)
                  lp0=nodelmpi(i,j+1)
               endif
#endif
               do m=1,n_group
                  j2_net_y_3D(m,l00,k)=-2d0*dfn(m,l00,k)*(flux2(l,k,m)-flux2(lm,k,m))
#ifdef js_mpi
                  if (comm%usempi) then
                     j2_net_y_3D(m,lp0,k)=-2d0*dfn(m,lp0,k)*(flux2(lp,k,m)-flux2(l,k,m))
                  endif
#endif
               enddo
            enddo
            if (idomy==1) then
               j=Iy_Start_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 903
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i,j-1)
               endif
#endif
               if (BC_Ly==1) then
                  do m=1,n_group
                     j2_net_y_3D(m,l00,k)=0d0
                  enddo
               elseif (BC_Ly==3.or.BC_Ly==2) then
                  do m=1,n_group
                     j2_net_y_3D(m,l00,k)=-2d0*dfn(m,l00,k)*flux2(l,k,m)
                  enddo
               endif
#ifdef js_mpi
903            continue
#endif
            endif
            if (idomy==ndomy) then
               j=Iy_End_x(i)
               l=nodel(i,j)
               lp=nodel(i,j+1)
               l0=l
               l00=l
               lp0=lp
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 904
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lp0=nodelmpi(i,j+1)
               endif
#endif
               if (BC_Ry==1) then
                  do m=1,n_group
                     j2_net_y_3D(m,lp0,k)=0d0
                  enddo
               elseif (BC_Ry==3.or.BC_Ry==2) then
                  do m=1,n_group
                     j2_net_y_3D(m,lp0,k)=2d0*dfn(m,lp0,k)*flux2(l,k,m)
                  enddo
               endif
#ifdef js_mpi
904            continue
#endif
            endif
         enddo
      enddo

      kbeg=kbot
      kend=ktop+1
      if (isub==1) then
         k=1
         if (BC_Lz==1) then
            do l=1,nxy
               l0=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
               endif
#endif
               do m=1,n_group
                  j2_net_z_3D(m,l0,k)=0d0
               enddo
            enddo
         elseif (BC_Lz==3.or.BC_Lz==2) then
            do l=1,nxy
               l0=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
               endif
#endif
               do m=1,n_group
                  j2_net_z_3D(m,l0,k)=-2d0*dfb(m,l0,k)*flux2(l,k,m)
                  !j2_net_z_3D(m,l0,k)=-(2d0*dfb(m,l0,k)+dnb2(m,l0,k))*flux2(l,k,m)
               enddo
            enddo
         endif
         kbeg=2
      endif
      if (isub==ndomz) then
         k=nz
         kp=k+1
         if (BC_Rz==1) then
            do l=1,nxy
               l0=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
               endif
#endif
               do m=1,n_group
                  j2_net_z_3D(m,l0,kp)=0d0
               enddo
            enddo
         elseif (BC_Rz==3.or.BC_Rz==2) then
            do l=1,nxy
               l0=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
               endif
#endif
               do m=1,n_group
                  j2_net_z_3D(m,l0,kp)=2d0*dfb(m,l0,kp)*flux2(l,k,m)
                  !j2_net_z_3D(m,l0,kp)=(2d0*dfb(m,l0,kp)-dnb2(m,l0,kp))*flux2(l,k,m)
               enddo
            enddo
         endif
         kend=nz
      endif
      do k=kbeg,kend
         km1=k-1
         do l=1,nxy
            l0=l
#ifdef js_mpi
            if (comm%usempi) then
               if (iproc/=ixy2iproc(l)) cycle
               l0=ixyip2ixy(l,iproc)
            endif
#endif
            do m=1,n_group
               j2_net_z_3D(m,l0,k)=-2d0*dfb(m,l0,k)*(flux2(l,k,m)-flux2(l,km1,m))
               !j2_net_z_3D(m,l0,k)=-2d0*dfb(m,l0,k)*(flux2(l,k,m)-flux2(l,km1,m)) &
               !                       -dnb2(m,l0,k)*(flux2(l,k,m)+flux2(l,km1,m))
            enddo
         enddo
      enddo

      return
      end subroutine updcur2
#endif

      END MODULE Mod_Soldhat

#elif siarhei_delete
      MODULE Mod_Soldhat
      USE Inc_maXS
      USE Inc_Control
      USE Inc_FluxVar
      USE Inc_Lscoef
!     ! USE Mod_Sol2N
      USE Inc_Geometry
      USE Inc_3D
      USE Inc_DF
      USE Inc_XS_File, only: Leakage_Table, Flag_Leakage, leak_ratio, Flag_LC_mean
      USE Inc_RP, only: AxialComp, I_LP_1N
      USE Inc_Option
      USE Inc_Transient, ONLY: Flag_Transient
      use Mod_GetNode, only: new_asym_itab
      use mod_charedit, only: print_msg


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      IMPLICIT NONE

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

      CONTAINS



      subroutine updcurnxy
      implicit none
      integer :: i,j,k,l,m,ktop,kbot
      integer :: lm,lp,kp,km


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [updcurnxy] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      kbot=1
      ktop=nz

      do k=kbot,ktop
         do j=1,ny
            do i=ix_start_y(j)+1,ix_end_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               do m=1,n_group
                  J_net_x_3D(m,l,k)=-dfw(m,l,k)*(flux(l,k,m)-flux(lm,k,m)) &
                                  & -dnw(m,l,k)*(flux(l,k,m)+flux(lm,k,m))
               enddo
            enddo
            if (idomx==1) then
               i=ix_start_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               if (BC_Lx==1) then
                  do m=1,n_group
                     j_net_x_3D(m,l,k)=0d0
                  enddo
               elseif (BC_Lx== 3.or.BC_Lx==2) then
                  do m=1,n_group
                     j_net_x_3D(m,l,k)=-(dfw(m,l,k)+dnw(m,l,k))*flux(l,k,m)
                  enddo
               endif
            endif
            if (idomx==ndomx) then
               i=ix_end_y(j)
               l=nodel(i,j)
               lp=nodel(i+1,j)
               if (BC_Rx==1) then
                  do m=1,n_group
                     j_net_x_3D(m,lp,k)=0
                  enddo
               elseif (BC_Rx==3.or.BC_Rx==2) then
                  do m=1,n_group
                     j_net_x_3D(m,lp,k)=(dfw(m,lp,k)-dnw(m,lp,k))*flux(l,k,m)
                  enddo
               endif
            endif
         enddo

         do i=1,nx
            do j=iy_start_x(i)+1,iy_end_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               do m=1,n_group
                  j_net_y_3D(m,l,k)=-dfn(m,l,k)*(flux(l,k,m)-flux(lm,k,m)) &
                                  & -dnn(m,l,k)*(flux(l,k,m)+flux(lm,k,m))
               enddo
            enddo
            if (idomy==1) then
               j=iy_start_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               if (BC_Ly==1) then
                  do m=1,n_group
                     j_net_y_3D(m,l,k)=0d0
                  enddo
               elseif (BC_Ly==3.or.BC_Ly==2) then
                  do m=1,n_group
                     j_net_y_3D(m,l,k)=-(dfn(m,l,k)+dnn(m,l,k))*flux(l,k,m)
                  enddo
               endif
            endif
            if (idomy==ndomy) then
               j=iy_end_x(i)
               l=nodel(i,j)
               lp=nodel(i,j+1)
               if (BC_Ry==1) then
                  do m=1,n_group
                     j_net_y_3D(m,lp,k)=0d0
                  enddo
               elseif (BC_Ry==3.or.BC_Ry==2) then
                  do m=1,n_group
                     j_net_y_3D(m,lp,k)=(dfn(m,lp,k)-dnn(m,lp,k))*flux(l,k,m)
                  enddo
               endif
            endif

         enddo
      enddo

      do l=1,nxy
         if (idomz==1) then
            k=1
            if (BC_Lz==1) then
               do m=1,n_group
                  j_net_z_3D(m,l,k)=0d0
               enddo
            elseif (BC_Lz==3.or.BC_Lz==2) then
               do m=1,n_group
                  j_net_z_3D(m,l,k)=-(dfb(m,l,k)+dnb(m,l,k))*flux(l,k,m)
               enddo
            endif
         endif
         if (idomz==ndomz) then
            k=nz
            kp=k+1
            if (BC_Rz==1) then
               do m=1,n_group
                  j_net_z_3D(m,l,kp)=0d0
               enddo
            elseif (BC_Rz==3.or.BC_Rz==2) then
               do m=1,n_group
                  j_net_z_3D(m,l,kp)=(dfb(m,l,kp)-dnb(m,l,kp))*flux(l,k,m)
               enddo
            endif
         endif
         do k=2,nz
            km=k-1
            do m=1,n_group
               j_net_z_3D(m,l,k)=-dfb(m,l,k)*(flux(l,k,m)-flux(l,km,m)) &
                               & -dnb(m,l,k)*(flux(l,k,m)+flux(l,km,m))
            enddo
         enddo
      enddo

      return
      end subroutine updcurnxy


      subroutine updtlxy
      use inc_solver, only: rvdelt, betap
      use inc_extsrc, only: flag_extsrc, extsrc
#ifdef js_sp3
      use inc_sp3, only: l2_2x_3d
      use inc_sp3, only: l2_2y_3d
      use inc_sp3, only: l2_2z_3d
      use inc_sp3, only: j2_net_x_3d
      use inc_sp3, only: j2_net_y_3d
      use inc_sp3, only: j2_net_z_3d
#endif
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: ixy2iproc, ixyip2ixy
      use inc_parallel, only: nodelmpi
#endif
      implicit none
      integer :: i,j,k,l,m,le,ls,lm,lp
      integer :: km1,kp1,jp1,isub,ktop,kbot
      real(8) :: w,wp,wm,rhx,rhy,rhz,rvol
      real(8) :: tl(n_group),tlp(n_group),tlm(n_group)
      integer :: l0, l00, le0, ls0, lm0, lp0


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [updtlxy] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      isub=1
      kbot=1
      ktop=nz

      do k=kbot,ktop
         do j=1,ny
            jp1=j+1
            rhy=1d0/MeshSize_y(j)
            do i=Ix_Start_y(j),Ix_End_y(j)
               l=nodel(i,j)
               rvol=1d0/NodeVolume(l,k)
               rhx=1d0/MeshSize_x(i)
               le=nodel(i+1,j)
               ls=nodel(i,jp1)
               l0=l
               l00=l
               le0=le
               ls0=ls
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  le0=nodelmpi(i+1,j)
                  ls0=nodelmpi(i,jp1)
               endif
#endif
               do m=1,n_group
                  if (flag_transient) then
                     if (n_group>2) then
                        src_tr(m,l00,k)=src(m,l,k)*rvol-rvdelt(m,l0,k)*flux(l,k,m) &
                        & -maxs_chi_3d(l,k,m)*(1d0-betap(l0,k))*FisSrc(l,k)*rvol
                     else
                        if (m==1) then
                           src_tr(1,l00,k)=src(1,l00,k)*rvol-rvdelt(1,l0,k)*flux(l,k,1) &
                           & -(1d0-betap(l0,k))*FisSrc(l,k)*rvol
                        else
                           src_tr(2,l00,k)=src(2,l,k)*rvol-rvdelt(2,l0,k)*flux(l,k,2)
                        endif
                     endif
                  else
                     src_tr(m,l00,k)=0d0
                     if (m==1) then
                        if (flag_extsrc) then
                           src_tr(1,l00,k)=extsrc(l,k)*rvol
                        endif
                     endif
                  endif
                  L_0x_3D(m,l00,k)=(j_net_x_3D(m,le0,k)-j_net_x_3D(m,l00,k))*rhx
                  L_0y_3D(m,l00,k)=(j_net_y_3D(m,ls0,k)-j_net_y_3D(m,l00,k))*rhy
               enddo
            enddo
         enddo
      enddo

#ifdef js_sp3
      if (opt_nodal==3) then
         do k=kbot,ktop
            do j=1,ny
               jp1=j+1
               rhy=1d0/MeshSize_y(j)
               do i=Ix_Start_y(j),Ix_End_y(j)
                  l=nodel(i,j)
                  rvol=1d0/NodeVolume(l,k)
                  rhx=1d0/MeshSize_x(i)
                  le=nodel(i+1,j)
                  ls=nodel(i,jp1)
                  l0=l
                  l00=l
                  le0=le
                  ls0=ls
#ifdef js_mpi
                  if (comm%usempi) then
                     if (iproc/=ixy2iproc(l)) cycle
                     l0=ixyip2ixy(l,iproc)
                     l00=nodelmpi(i,j)
                     le0=nodelmpi(i+1,j)
                     ls0=nodelmpi(i,jp1)
                  endif
#endif
                  do m=1,n_group
                     L2_2x_3D(m,l00,k)=(j2_net_x_3D(m,le0,k)-j2_net_x_3D(m,l00,k))*rhx
                     L2_2y_3D(m,l00,k)=(j2_net_y_3D(m,ls0,k)-j2_net_y_3D(m,l00,k))*rhy
                  enddo
               enddo
            enddo
         enddo
         do k=kbot,ktop
            rhz=1d0/MeshSize_z(k)
            kp1=k+1
            do l=1,nxy
               l0=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
               endif
#endif
               do m=1,n_group
                  L2_2z_3D(m,l0,k)=(j2_net_z_3D(m,l0,kp1)-j2_net_z_3D(m,l0,k))*rhz
               enddo
            enddo
         enddo
         ! update source term due to the leakage of 2nd moment
         do k=1,nz
            kp1=k+1
            km1=k-1
            if (kp1>nz) kp1=nz
            if (km1<1) km1=1
            do l=1,nxy
               l0=l
               l00=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(ltox(l),ltoy(l))
               endif
#endif
               do m=1,n_group
                  src_tr(m,l00,k)=src_tr(m,l00,k) &
                  & -(L2_2x_3D(m,l00,k)+L2_2y_3D(m,l00,k)+L2_2z_3D(m,l0,k))
               enddo
            enddo
         enddo
      endif ! opt_nodal==3
#endif

      if (nz==1) return !return if 2d

      do k=kbot,ktop
         rhz=1d0/MeshSize_z(k)
         kp1=k+1
         do l=1,nxy
            l0=l
            l00=l
#ifdef js_mpi
            if (comm%usempi) then
               if (iproc/=ixy2iproc(l)) cycle
               l0=ixyip2ixy(l,iproc)
               l00=nodelmpi(ltox(l),ltoy(l))
            endif
#endif
            do m=1,n_group
               L_0z_3D(m,l00,k)=(j_net_z_3D(m,l0,kp1)-j_net_z_3D(m,l0,k))*rhz
            enddo
         enddo
      enddo

      do k=kbot,ktop
         do j=1,ny
            do i=Ix_Start_y(j),Ix_End_y(j)
               w=MeshSize_x(i)
               l=nodel(i,j)
               l0=l
               l00=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
               endif
#endif
               do m=1,n_group
                  tl(m)=L_0y_3D(m,l00,k)+L_0z_3D(m,l00,k)-src_tr(m,l00,k)
               enddo
               if (i>Ix_Start_y(j)) then
                  wm=MeshSize_x(i-1)
                  lm=nodel(i-1,j)
                  lm0=lm
#ifdef js_mpi
                  if (comm%usempi) then
                     lm0=nodelmpi(i-1,j)
                  endif
#endif
                  do m=1,n_group
                     tlm(m)=L_0y_3D(m,lm0,k)+L_0z_3D(m,lm0,k)-src_tr(m,lm0,k)
                  enddo
               elseif (BC_Lx==1) then
                  wm=w
                  do m=1,n_group
                     tlm(m)=tl(m)
                  enddo
               else
                  wm=0d0
                  do m=1,n_group
                     tlm(m)=0d0
                  enddo
               endif
               if (i<Ix_End_y(j)) then
                  wp=MeshSize_x(i+1)
                  lp=nodel(i+1,j)
                  lp0=lp
#ifdef js_mpi
                  if (comm%usempi) then
                     lp0=nodelmpi(i+1,j)
                  endif
#endif
                  do m=1,n_group
                     tlp(m)=L_0y_3D(m,lp0,k)+L_0z_3D(m,lp0,k)-src_tr(m,lp0,k)
                  enddo
               elseif (BC_Rx==1) then
                  wp=w
                  do m=1,n_group
                     tlp(m)=tl(m)
                  enddo
               else
                  wp=0d0
                  do m=1,n_group
                     tlp(m)=0d0
                  enddo
               endif
               do m=1,n_group
                  L_1x_3D(m,l00,k)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
                  L_2x_3D(m,l00,k)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
               enddo
            enddo
         enddo

         do i=1,nx
            do j=Iy_Start_x(i),Iy_End_x(i)
               w=MeshSize_y(j)
               l=nodel(i,j)
               l0=l
               l00=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
               endif
#endif
               do m=1,n_group
                  tl(m)=L_0x_3D(m,l00,k)+L_0z_3D(m,l00,k)-src_tr(m,l00,k)
               enddo
               if (j>Iy_Start_x(i)) then
                  wm=MeshSize_y(j-1)
                  lm=nodel(i,j-1)
                  lm0=lm
#ifdef js_mpi
                  if (comm%usempi) then
                     lm0=nodelmpi(i,j-1)
                  endif
#endif
                  do m=1,n_group
                     tlm(m)=L_0x_3D(m,lm0,k)+L_0z_3D(m,lm0,k)-src_tr(m,lm0,k)
                  enddo
               elseif (BC_Ly==1) then
                  wm=w
                  do m=1,n_group
                     tlm(m)=tl(m)
                  enddo
               else
                  wm=0d0
                  do m=1,n_group
                     tlm(m)=0d0
                  enddo
               endif
               if (j<Iy_End_x(i)) then
                  wp=MeshSize_y(j+1)
                  lp=nodel(i,j+1)
                  lp0=lp
#ifdef js_mpi
                  if (comm%usempi) then
                     lp0=nodelmpi(i,j+1)
                  endif
#endif
                  do m=1,n_group
                    tlp(m)=L_0x_3D(m,lp0,k)+L_0z_3D(m,lp0,k)-src_tr(m,lp0,k)
                  enddo
               elseif (BC_Ry==1) then
                  wp=w
                  do m=1,n_group
                    tlp(m)=tl(m)
                  enddo
               else
                  wp=0d0
                  do m=1,n_group
                    tlp(m)=0d0
                  enddo
               endif
               do m=1,n_group
                  L_1y_3D(m,l00,k)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
                  L_2y_3D(m,l00,k)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
               enddo
            enddo
         enddo
      enddo

      if (isub==1) then
         k=1
         kp1=k+1
         w=MeshSize_z(k)
         wp=MeshSize_z(kp1)
         do l=1,nxy
            l0=l
            l00=l
#ifdef js_mpi
            if (comm%usempi) then
               if (iproc/=ixy2iproc(l)) cycle
               l0=ixyip2ixy(l,iproc)
               l00=nodelmpi(ltox(l),ltoy(l))
            endif
#endif
            do m=1,n_group
               tl(m)=L_0x_3D(m,l00,k)+L_0y_3D(m,l00,k)-src_tr(m,l00,k)
               tlp(m)=L_0x_3D(m,l00,kp1)+L_0y_3D(m,l00,kp1)-src_tr(m,l00,kp1)
               if (BC_Lz==1) then
                  wm=w
                  tlm(m)=tl(m)
               else
                  wm=0d0
                  tlm(m)=0d0
               endif
               L_1z_3D(m,l0,k)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2z_3D(m,l0,k)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      endif

      if (isub==ndomz) then
         k=nz
         km1=k-1
         w=MeshSize_z(k)
         wm=MeshSize_z(km1)
         do l=1,nxy
            l0=l
            l00=l
#ifdef js_mpi
            if (comm%usempi) then
               if (iproc/=ixy2iproc(l)) cycle
               l0=ixyip2ixy(l,iproc)
               l00=nodelmpi(ltox(l),ltoy(l))
            endif
#endif
            do m=1,n_group
               tl(m)=L_0x_3D(m,l00,k)+L_0y_3D(m,l00,k)-src_tr(m,l00,k)
               tlm(m)=L_0x_3D(m,l00,km1)+L_0y_3D(m,l00,km1)-src_tr(m,l00,km1)
               if (BC_Rz==1) then
                  wp=w
                  tlp(m)=tl(m)
               else
                  wp=0d0
                  tlp(m)=0d0
               endif
               L_1z_3D(m,l0,k)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2z_3D(m,l0,k)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      endif

      do k=2,nz-1
         kp1=k+1
         km1=k-1
         w=MeshSize_z(k)
         wm=MeshSize_z(km1)
         wp=MeshSize_z(kp1)
         do l=1,nxy
            l0=l
            l00=l
#ifdef js_mpi
            if (comm%usempi) then
               if (iproc/=ixy2iproc(l)) cycle
               l0=ixyip2ixy(l,iproc)
               l00=nodelmpi(ltox(l),ltoy(l))
            endif
#endif
            do m=1,n_group
               tl(m)=L_0x_3D(m,l00,k)+L_0y_3D(m,l00,k)-src_tr(m,l00,k)
               tlp(m)=L_0x_3D(m,l00,kp1)+L_0y_3D(m,l00,kp1)-src_tr(m,l00,kp1)
               tlm(m)=L_0x_3D(m,l00,km1)+L_0y_3D(m,l00,km1)-src_tr(m,l00,km1)
               L_1z_3D(m,l0,k)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2z_3D(m,l0,k)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      enddo

      return
      end subroutine updtlxy


      FUNCTION pol1 (w,wp,wm,tl,tlp,tlm)
      IMPLICIT NONE
      REAL(8) :: pol1, w, wp, wm, tl, tlp, tlm, dd
      dd=one/((w+wp)*(w+wm)*(w+wm+wp))

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [pol1] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      pol1=dd*w*((tlp-tl)*(w+2*wm)*(w+wm)+(tl-tlm)*(w+2*wp)*(w+wp))
      RETURN
      END FUNCTION pol1

      FUNCTION pol2 (w,wp,wm,tl,tlp,tlm)
      IMPLICIT NONE
      REAL(8) :: pol2,w, wp, wm, tl, tlp, tlm, dd
      dd=one/((w+wp)*(w+wm)*(w+wm+wp))

#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [pol2] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      pol2=dd*w*w*((tlp-tl)*(w+wm)+(tlm-tl)*(w+wp))
      RETURN
      END FUNCTION pol2


      subroutine init_subaxial_nodal
      implicit none
      real(8) :: hzori
      integer :: iz, izs, ize


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [init_subaxial_nodal] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.allocated(corr_ADF_Lz)) allocate(corr_ADF_Lz(nxy,nz,n_group))
      if (.not.allocated(corr_ADF_Rz)) allocate(corr_ADF_Rz(nxy,nz,n_group))
      corr_ADF_Lz=1d0
      corr_ADF_Rz=1d0

      if (.not.flag_axialdf) return
      if (nz==1) then
         flag_axialdf=.false.
         return
      endif

      ! geometry for sub-axial mesh (around 1cm)
      if (.not.allocated(nz_subz)) allocate(nz_subz(1:nz))
      nz_subz=1
      do iz=1,nz
         hzori=meshsize_z(iz)
         nz_subz(iz)=max(1,int(hzori))
      enddo
      nsubz=sum(nz_subz)
      if (.not.allocated(inp2sub_beg)) allocate(inp2sub_beg(1:nz))
      if (.not.allocated(inp2sub_end)) allocate(inp2sub_end(1:nz))
      inp2sub_beg=0
      inp2sub_end=0
      inp2sub_beg(1)=1
      inp2sub_end(1)=nz_subz(1)
      do iz=2,nz
         inp2sub_beg(iz)=inp2sub_end(iz-1)+1
         inp2sub_end(iz)=inp2sub_end(iz-1)+nz_subz(iz)
      enddo
      if (.not.allocated(sub2inp)) allocate(sub2inp(1:nsubz))
      sub2inp=0
      izs=1
      ize=1
      do iz=1,nz
         izs=inp2sub_beg(iz)
         ize=inp2sub_end(iz)
         sub2inp(izs:ize)=iz
      enddo
      if (.not.allocated(hsubz)) allocate(hsubz(1:nsubz))
      hsubz=0d0
      izs=1
      ize=1
      do iz=1,nz
         izs=inp2sub_beg(iz)
         ize=inp2sub_end(iz)
         hsubz(izs:ize)=meshsize_z(iz)/real(nz_subz(iz),8)
      enddo

      if (.not.allocated(L_1subz_3D)) allocate(L_1subz_3D(1:n_group,1:nxy,0:nsubz))
      if (.not.allocated(L_2subz_3D)) allocate(L_2subz_3D(1:n_group,1:nxy,0:nsubz))
      L_1subz_3D=0d0
      L_2subz_3D=0d0
      if (.not.allocated(phisurf_inpz)) allocate(phisurf_inpz(n_group,nxy,nz+1))
      if (.not.allocated(phisurf_subz)) allocate(phisurf_subz(n_group,nxy,nsubz+1))
      phisurf_inpz=0d0
      phisurf_subz=0d0

      return
      end subroutine init_subaxial_nodal

#ifdef siarhei_delete
      subroutine update_axialdf
      implicit none
      type(UNMNodal) node_l, node_r
      real(8) :: rl(n_group)
      real(8) :: rlp(n_group)
      real(8) :: rm1(n_group)
      real(8) :: rmp1(n_group)
      real(8) :: rm2(n_group)
      real(8) :: rmp2(n_group)
      real(8) :: df(n_group)
      real(8) :: dfp(n_group)
      real(8) :: tl(n_group)
      real(8) :: tlp(n_group)
      real(8) :: tlm(n_group)
      integer :: g, gg
      integer :: l, m
      integer :: izs, ize
      integer :: k, kp, km
      integer :: subk, subkp, subkm
      real(8) :: w, wp, wm


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [update_axialdf] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      if (.not.flag_axialdf) return
      if (sum(phisurf_inpz(:,:,:))<1d-10) return

      ! update sub-axial mesh transverse leakage
      if (idomz==1) then
         subk=1
         subkp=subk+1
         w=hsubz(subk)
         wp=hsubz(subkp)
         k=sub2inp(subk)
         kp=sub2inp(subkp)
         do l=1,nxy
            do m=1,n_group
               tl(m)=L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
               tlp(m)=L_0x_3D(m,l,kp)+L_0y_3D(m,l,kp)-src_tr(m,l,kp)
               if (BC_Lz==1) then
                  wm=w
                  tlm(m)=tl(m)
               else
                  wm=0d0
                  tlm(m)=0d0
               endif
               L_1subz_3D(m,l,subk)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2subz_3D(m,l,subk)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      endif
      do subk=2,nsubz-1
         subkp=subk+1
         subkm=subk-1
         k=sub2inp(subk)
         kp=sub2inp(subkp)
         km=sub2inp(subkm)
         w=hsubz(subk)
         wm=hsubz(subkm)
         wp=hsubz(subkp)
         do l=1,nxy
            do m=1,n_group
               tl(m)=L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
               tlp(m)=L_0x_3D(m,l,kp)+L_0y_3D(m,l,kp)-src_tr(m,l,kp)
               tlm(m)=L_0x_3D(m,l,km)+L_0y_3D(m,l,km)-src_tr(m,l,km)
               L_1subz_3D(m,l,subk)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2subz_3D(m,l,subk)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      enddo
      if (idomz==ndomz) then
         subk=nsubz
         subkm=subk-1
         k=sub2inp(subk)
         km=sub2inp(subkm)
         w=hsubz(subk)
         wm=hsubz(subkm)
         do l=1,nxy
            do m=1,n_group
               tl(m)=L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
               tlm(m)=L_0x_3D(m,l,km)+L_0y_3D(m,l,km)-src_tr(m,l,km)
               if (BC_Rz==1) then
                  wp=w
                  tlp(m)=tl(m)
               else
                  wp=0d0
                  tlp(m)=0d0
               endif
               L_1subz_3D(m,l,subk)=pol1(w,wp,wm,tl(m),tlp(m),tlm(m))
               L_2subz_3D(m,l,subk)=pol2(w,wp,wm,tl(m),tlp(m),tlm(m))
            enddo
         enddo
      endif

      ! solve 1d nodal for sub-axial mesh
      call Allocate_one_dim_node(node_l,n_group)
      call Allocate_one_dim_node(node_r,n_group)
      do l=1,Nxy
         if (idomz==1) then
            subk=1
            k=sub2inp(subk)
            select case (BC_Lz)
            case (1)
               phisurf_subz(:,l,subk)=flux(l,k,:)
            case (2,3)
               do m=1,n_group
                  rl(m)=L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
                  rm1(m)=L_1subz_3D(m,l,subk) !re-calculate line 1503~1589
                  rm2(m)=L_2subz_3D(m,l,subk) !same... search MeshSize_z
                  rm1(m)=rm1(m)*0.5d0
                  rm2(m)=rm2(m)*0.5d0
               enddo
               node_r%width = hsubz(subk)
               do g=1,n_group
                  node_r%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                  node_r%source0(g)    =-rl(g)
                  node_r%source1(g)    =-rm1(g)
                  node_r%source2(g)    =-rm2(g)
                  node_r%adf(g)        = 1d0
                  node_r%flux_avg(g)   = Flux(l,k,g)
                  node_r%diff_coeff(g) = D_3D(l,k,g)
                  node_r%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                  node_r%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                  do gg=1,n_group
                     node_r%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                  enddo
               enddo
               select case (opt_nodal)
               case (1)
                  call solve_left_bndry_node_unm(BC_Lz,node_r,keff,OPT_NEM)
               case (2)
                  call solve_left_bndry_node_unm(BC_Lz,node_r,keff,OPT_ANM)
               case (3)
                  call solve_left_bndry_node_unm(BC_Lz,node_r,keff,OPT_ANM)
               case default
                  call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                  stop
               end select
               do m=1,n_group
                  phisurf_subz(m,l,subk)=face_flux(m)
               enddo
            end select
         endif

         do subk=1,nsubz-1
            subkp=subk+1
            k=sub2inp(subk)
            kp=sub2inp(subkp)
            do m=1,n_group
               rl(m)   = L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
               rlp(m)  = L_0x_3D(m,l,kp)+L_0y_3D(m,l,kp)-src_tr(m,l,kp)
               rm1(m)  = L_1subz_3D(m,l,subk)
               rmp1(m) = L_1subz_3D(m,l,subkp)
               rm2(m)  = L_2subz_3D(m,l,subk)
               rmp2(m) = L_2subz_3D(m,l,subkp)
               df(m)   = ADF_Rz(l,k,m)
               dfp(m)  = ADF_Lz(l,kp,m)
               rm1(m)  = rm1(m)*0.5d0
               rm2(m)  = rm2(m)*0.5d0
               rmp1(m) = rmp1(m)*0.5d0
               rmp2(m) = rmp2(m)*0.5d0
            enddo
            node_l%width = hsubz(subk)
            node_r%width = hsubz(subkp)
            do g=1,n_group
               node_l%fiss_spec(g)  = maxs_chi_3d(l,k,g)
               node_r%fiss_spec(g)  = maxs_chi_3d(l,kp,g)
               node_l%source0(g)    =-rl(g)
               node_r%source0(g)    =-rlp(g)
               node_l%source1(g)    =-rm1(g)
               node_r%source1(g)    =-rmp1(g)
               node_l%source2(g)    =-rm2(g)
               node_r%source2(g)    =-rmp2(g)
               node_l%adf(g)        = df(g)
               node_r%adf(g)        = dfp(g)
               node_l%flux_avg(g)   = Flux(l,k,g)
               node_r%flux_avg(g)   = Flux(l,kp,g)
               node_l%diff_coeff(g) = D_3D(l,k,g)
               node_r%diff_coeff(g) = D_3D(l,kp,g)
               node_l%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
               node_r%sigma_r(g)    = maXs_r_3D(l,kp,g)+maXS_scat_3D(g,g,l,kp)
               node_l%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
               node_r%nu_sigma_f(g) = nu_maXs_f_3D(l,kp,g)
               do gg=1,n_group
                  node_l%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                  node_r%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,kp)
               enddo
            enddo
            select case (opt_nodal)
            case (1)
               call solve_two_node_unm(node_l,node_r,keff,OPT_NEM)
            case (2)
               call solve_two_node_unm(node_l,node_r,keff,OPT_ANM)
            case (3)
               call solve_two_node_unm(node_l,node_r,keff,OPT_ANM)
            case default
               call print_msg(3,"updatedhat: opt_nodal is not assgiend")
               stop
            end select
            do m=1,n_group
               phisurf_subz(m,l,subkp)=face_flux(m)
            enddo
         enddo

         if (idomz==ndomz) then
            subk=nsubz
            subkp=nsubz+1
            k=sub2inp(subk)
            kp=k+1
            select case (BC_Rz)
            case (1)
               phisurf_subz(:,l,subkp)=flux(l,k,:)
            case (2,3)
               do m=1,n_group
                  rl(m) =L_0x_3D(m,l,k)+L_0y_3D(m,l,k)-src_tr(m,l,k)
                  rm1(m)=L_1subz_3D(m,l,subk)
                  rm2(m)=L_2subz_3D(m,l,subk)
                  rm1(m)=rm1(m)*0.5d0
                  rm2(m)=rm2(m)*0.5d0
               enddo
               node_l%width = hsubz(subk)
               do g=1,n_group
                  node_l%fiss_spec(g)  = maxs_chi_3d(l,k,g)
                  node_l%source0(g)    =-rl(g)
                  node_l%source1(g)    =-rm1(g)
                  node_l%source2(g)    =-rm2(g)
                  node_l%adf(g)        = ADF_Rz(l,k,g)
                  node_l%flux_avg(g)   = Flux(l,k,g)
                  node_l%diff_coeff(g) = D_3D(l,k,g)
                  node_l%sigma_r(g)    = maXs_r_3D(l,k,g)+maXS_scat_3D(g,g,l,k)
                  node_l%nu_sigma_f(g) = nu_maXs_f_3D(l,k,g)
                  do gg=1,n_group
                     node_l%sigma_s(g,gg) = maXS_scat_3D(gg,g,l,k)
                  enddo
               enddo
               select case (opt_nodal)
               case (1)
                  call solve_right_bndry_node_unm(BC_Rz,node_l,keff,OPT_NEM)
               case (2)
                  call solve_right_bndry_node_unm(BC_Rz,node_l,keff,OPT_ANM)
               case (3)
                  call solve_right_bndry_node_unm(BC_Rz,node_l,keff,OPT_ANM)
               case default
                  call print_msg(3,"updatedhat: opt_nodal is not assgiend")
                  stop
               end select
               do m=1,n_group
                  phisurf_subz(m,l,subkp)=face_flux(m)
               enddo
            end select
         endif
      enddo

      call deallocate_one_dim_node(node_r)
      call deallocate_one_dim_node(node_l)

      ! calculate correction factor of  axial DF
      do l=1,nxy
         do k=1,nz
            izs=inp2sub_beg(k)
            ize=inp2sub_end(k)
            do m=1,n_group
               corr_ADF_Lz(l,k,m)=phisurf_inpz(m,l,k)/phisurf_subz(m,l,izs)
               if (corr_ADF_Lz(l,k,m)>2d0.or.corr_ADF_Lz(l,k,m)<0d0) then
                  corr_ADF_Lz(l,k,m)=1d0
               else
                  corr_ADF_Lz(l,k,m)=corr_ADF_Lz(l,k,m)/ADF_Lz(l,k,m)
               endif
            enddo
         enddo
      enddo
      do l=1,nxy
         do k=1,nz
            izs=inp2sub_beg(k)
            ize=inp2sub_end(k)
            do m=1,n_group
               corr_ADF_Rz(l,k,m)=phisurf_inpz(m,l,k+1)/phisurf_subz(m,l,ize+1)
               if (corr_ADF_Rz(l,k,m)>2d0.or.corr_ADF_Rz(l,k,m)<0d0) then
                  corr_ADF_Rz(l,k,m)=1d0
               else
                  corr_ADF_Rz(l,k,m)=corr_ADF_Rz(l,k,m)/ADF_rz(l,k,m)
               endif
            enddo
         enddo
      enddo

      return
      end subroutine update_axialdf
#endif
      
#ifdef js_sp3
      subroutine updatedhat2d1d
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: ixy2iproc, ixyip2ixy
      use inc_parallel, only: nodelmpi
      use mod_parallel, only: barrier
#endif
      implicit none
      type(UNMNodal) node_l, node_r
      real(8) :: rl(n_group), rlp(n_group)
      real(8) :: rm1(n_group), rmp1(n_group)
      real(8) :: rm2(n_group), rmp2(n_group)
      real(8) :: df(n_group), dfp(n_group)
      integer :: g, gg, isub, kp, kl, klp
      integer :: kbot, ktop, l, k, m
      integer :: l0, l00


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [updatedhat2d1d] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      call initialize_two_node_unm_solver( n_group )
      call allocate_one_dim_node( node_l, n_group )
      call allocate_one_dim_node( node_r, n_group )
      ! call updtlxy2d1d => updcurnxy, updtlxy, updtlcoeff2
      call updcurnxy_2d1d
      call updtlxy

      isub=1
      kbot=1
      ktop=nz

      ! sweep along nodes in axial plane.

      ! update dnb(edge) using dfb(edge) which exists.
      do l=1,nxy
         l0=l
         l00=l
#ifdef js_mpi
         if (comm%usempi) then
            if (iproc/=ixy2iproc(l)) cycle
            l0=ixyip2ixy(l,iproc)
            l00=nodelmpi(ltox(l),ltoy(l))
         endif
#endif
         ! solve one-node problem for bottom b.c. (k=1)
         k=kbot
         if (BC_Lz==1) then ! reflective
            do m=1,n_group
               dnb(m,l0,k)=0d0
            enddo
         elseif (BC_Lz==2.or.BC_Lz==3) then ! Black / Fluxzero
            do m=1,n_group
               rl(m) =L_0x_3D(m,l00,k)+L_0y_3D(m,l00,k)-src_tr(m,l00,k)
               rm1(m)=L_1z_3D(m,l0,k)
               rm2(m)=L_2z_3D(m,l0,k)
               rm1(m)=rm1(m)*0.5d0
               rm2(m)=rm2(m)*0.5d0
            enddo
            node_r%width=meshsize_z(k)
            do g=1,n_group
               node_r%fiss_spec(g)  = maxs_chi_3d(l,k,g)
               node_r%source0(g)    =-rl(g)
               node_r%source1(g)    =-rm1(g)
               node_r%source2(g)    =-rm2(g)
               node_r%adf(g)        = 1d0 !js+ADF_Lz(l,k,g)
               node_r%flux_avg(g)   = flux(l,k,g)
               node_r%diff_coeff(g) = d_3d(l,k,g)
               node_r%sigma_r(g)    = maxs_r_3d(l,k,g)+maxs_scat_3d(g,g,l,k)
               node_r%nu_sigma_f(g) = nu_maxs_f_3d(l,k,g)
               do gg=1,n_group
                  node_r%sigma_s(g,gg) = maxs_scat_3d(gg,g,l,k)
               enddo
            enddo
            call solve_left_bndry_node_unm( BC_Lz, node_r, keff, OPT_ANM )
            do m=1,n_group
               dnb(m,l0,k)=-dfb(m,l0,k)-face_curr(m)/flux(l,k,m)
            enddo
         endif

         ! solve for the (nz-1) two-node problems in the z-direction
         do k=kbot,ktop-1
            kl=k
            klp=k+1
            do m=1,n_group
               rl(m)   = L_0x_3D(m,l00,kl) +L_0y_3D(m,l00,kl) -src_tr(m,l00,kl)
               rlp(m)  = L_0x_3D(m,l00,klp)+L_0y_3D(m,l00,klp)-src_tr(m,l00,klp)
               rm1(m)  = L_1z_3D(m,l0,kl)
               rmp1(m) = L_1z_3D(m,l0,klp)
               rm2(m)  = L_2z_3D(m,l0,kl)
               rmp2(m) = L_2z_3D(m,l0,klp)
               df(m)   = 1d0 !js+ADF_Rz(l,kl,m)
               dfp(m)  = 1d0 !js+ADF_Lz(l,klp,m)
               rm1(m)  = rm1(m)*0.5d0
               rmp1(m) = rmp1(m)*0.5d0
               rm2(m)  = rm2(m)*0.5d0
               rmp2(m) = rmp2(m)*0.5d0
            enddo
            node_l%width = meshsize_z(kl)
            node_r%width = meshsize_z(klp)
            do g=1,n_group
               node_l%fiss_spec(g)  = maxs_chi_3d(l,kl,g)
               node_r%fiss_spec(g)  = maxs_chi_3d(l,klp,g)
               node_l%source0(g)    = -rl(g)
               node_r%source0(g)    = -rlp(g)
               node_l%source1(g)    = -rm1(g)
               node_r%source1(g)    = -rmp1(g)
               node_l%source2(g)    = -rm2(g)
               node_r%source2(g)    = -rmp2(g)
               node_l%adf(g)        = df(g)
               node_r%adf(g)        = dfp(g)
               node_l%flux_avg(g)   = flux(l,kl,g)
               node_r%flux_avg(g)   = flux(l,klp,g)
               node_l%diff_coeff(g) = d_3d(l,kl,g)
               node_r%diff_coeff(g) = d_3d(l,klp,g)
               node_l%sigma_r(g)    = maxs_r_3d(l,kl,g) +maxs_scat_3d(g,g,l,kl)
               node_r%sigma_r(g)    = maxs_r_3d(l,klp,g)+maxs_scat_3d(g,g,l,klp)
               node_l%nu_sigma_f(g) = nu_maxs_f_3d(l,kl,g)
               node_r%nu_sigma_f(g) = nu_maxs_f_3d(l,klp,g)
               do gg=1,n_group
                  node_l%sigma_s(g,gg)=maxs_scat_3d(gg,g,l,kl)
                  node_r%sigma_s(g,gg)=maxs_scat_3d(gg,g,l,klp)
               enddo
            enddo
            call solve_two_node_unm( node_l, node_r, keff, OPT_ANM )
            do m=1,n_group
               dnb(m,l0,klp)=-(dfb(m,l0,klp)*(flux(l,klp,m)-flux(l,kl,m))+face_curr(m)) &
                              /(flux(l,klp,m)+flux(l,kl,m))
            enddo
         enddo

         ! solve for the one-node problem for k=nz only if j(edge)<>0
         k=nz
         kp=nz+1
         if (BC_Rz==1) then ! Refelctive
            do m=1,n_group
               dnb(m,l0,kp)=0d0
            enddo
         elseif (BC_Rz==2.or.BC_Rz==3) then ! Black / Flux zero
            do m=1,n_group
               rl(m) =L_0x_3D(m,l00,k)+L_0y_3D(m,l00,k)-src_tr(m,l00,k)
               rm1(m)=L_1z_3D(m,l0,k)
               rm2(m)=L_2z_3D(m,l0,k)
               rm1(m)=rm1(m)*0.5d0
               rm2(m)=rm2(m)*0.5d0
            enddo
            node_l%width=meshsize_z(k)
            do g=1,n_group
               node_l%fiss_spec(g)  = maxs_chi_3d(l,k,g)
               node_l%source0(g)    = -rl(g)
               node_l%source1(g)    = -rm1(g)
               node_l%source2(g)    = -rm2(g)
               node_l%adf(g)        = 1d0 !js+ADF_Rz(l,k,g)
               node_l%flux_avg(g)   = flux(l,k,g)
               node_l%diff_coeff(g) = d_3d(l,k,g)
               node_l%sigma_r(g)    = maxs_r_3d(l,k,g)+maxs_scat_3d(g,g,l,k)
               node_l%nu_sigma_f(g) = nu_maxs_f_3d(l,k,g)
               do gg=1,n_group
                  node_l%sigma_s(g,gg) = maxs_scat_3d(gg,g,l,k)
               enddo
            enddo
            call solve_right_bndry_node_unm ( BC_Rz, node_l, keff, OPT_ANM )
            do m=1,n_group
               dnb(m,l0,kp)=-face_curr(m)/flux(l,k,m)+dfb(m,l0,kp)
            enddo
         endif
      enddo

      call deallocate_one_dim_node( node_r )
      call deallocate_one_dim_node( node_l )
      call finalize_two_node_unm_solver( )

      !js+if (opt_nodal==3) then
      !js+   call upddhat2d1d2 ! not yet
      !js+endif

#ifdef js_mpi
      if (comm%usempi) call barrier
#endif

      return
      end subroutine updatedhat2d1d
#endif
#ifdef js_sp3
      subroutine updcurnxy_2d1d
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: ixy2iproc, ixyip2ixy
      use inc_parallel, only: nodelmpi
#endif
      implicit none
      integer :: i,j,k,l,m,ktop,kbot
      integer :: lm,lp,kp,km
      integer :: l0, l00, lm0, lp0


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [updcurnxy_2d1d] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      kbot=1
      ktop=nz

      do k=kbot,ktop
         do j=1,ny
            do i=ix_start_y(j)+1,ix_end_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i-1,j)
                  lp=nodel(i+1,j)
                  lp0=nodelmpi(i+1,j)
               endif
#endif
               do m=1,n_group
                  j_net_x_3D(m,l00,k)=-dfw(m,l00,k)*(flux(l,k,m)-flux(lm,k,m))
                  !J_net_x_3D(m,l00,k)=-dfw(m,l00,k)*(flux(l,k,m)-flux(lm,k,m)) &
                  !                  & -dnw(m,l00,k)*(flux(l,k,m)+flux(lm,k,m))
#ifdef js_mpi
                  if (comm%usempi) then
                     j_net_x_3D(m,lp0,k)=-dfw(m,lp0,k)*(flux(lp,k,m)-flux(l,k,m))
                  endif
#endif
               enddo
            enddo
            if (idomx==1) then
               i=ix_start_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 901
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i-1,j)
               endif
#endif
               if (BC_Lx==1) then
                  do m=1,n_group
                     j_net_x_3D(m,l00,k)=0d0
                  enddo
               elseif (BC_Lx== 3.or.BC_Lx==2) then
                  do m=1,n_group
                     j_net_x_3D(m,l00,k)=-dfw(m,l00,k)*flux(l,k,m)
                     !j_net_x_3D(m,l00,k)=-(dfw(m,l00,k)+dnw(m,l00,k))*flux(l,k,m)
                  enddo
               endif
#ifdef js_mpi
901            continue
#endif
            endif
            if (idomx==ndomx) then
               i=ix_end_y(j)
               l=nodel(i,j)
               lp=nodel(i+1,j)
               l0=l
               l00=l
               lp0=lp
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 902
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lp0=nodelmpi(i+1,j)
               endif
#endif
               if (BC_Rx==1) then
                  do m=1,n_group
                     j_net_x_3D(m,lp0,k)=0
                  enddo
               elseif (BC_Rx==3.or.BC_Rx==2) then
                  do m=1,n_group
                     j_net_x_3D(m,lp0,k)=dfw(m,lp0,k)*flux(l,k,m)
                     !j_net_x_3D(m,lp0,k)=(dfw(m,lp0,k)-dnw(m,lp0,k))*flux(l,k,m)
                  enddo
               endif
#ifdef js_mpi
902            continue
#endif
            endif
         enddo

         do i=1,nx
            do j=iy_start_x(i)+1,iy_end_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i,j-1)
                  lp=nodel(i,j+1)
                  lp0=nodelmpi(i,j+1)
               endif
#endif
               do m=1,n_group
                  j_net_y_3D(m,l00,k)=-dfn(m,l00,k)*(flux(l,k,m)-flux(lm,k,m))
                  !j_net_y_3D(m,l00,k)=-dfn(m,l00,k)*(flux(l,k,m)-flux(lm,k,m)) &
                  !                  & -dnn(m,l00,k)*(flux(l,k,m)+flux(lm,k,m))
#ifdef js_mpi
                  if (comm%usempi) then
                     j_net_y_3D(m,lp0,k)=-dfn(m,lp0,k)*(flux(lp,k,m)-flux(l,k,m))
                  endif
#endif
               enddo
            enddo
            if (idomy==1) then
               j=iy_start_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 903
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i,j-1)
               endif
#endif
               if (BC_Ly==1) then
                  do m=1,n_group
                     j_net_y_3D(m,l00,k)=0d0
                  enddo
               elseif (BC_Ly==3.or.BC_Ly==2) then
                  do m=1,n_group
                     j_net_y_3D(m,l00,k)=-dfn(m,l00,k)*flux(l,k,m)
                     !j_net_y_3D(m,l00,k)=-(dfn(m,l00,k)+dnn(m,l00,k))*flux(l,k,m)
                  enddo
               endif
#ifdef js_mpi
903            continue
#endif
            endif
            if (idomy==ndomy) then
               j=iy_end_x(i)
               l=nodel(i,j)
               lp=nodel(i,j+1)
               l0=l
               l00=l
               lp0=lp
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 904
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lp0=nodelmpi(i,j+1)
               endif
#endif
               if (BC_Ry==1) then
                  do m=1,n_group
                     j_net_y_3D(m,lp0,k)=0d0
                  enddo
               elseif (BC_Ry==3.or.BC_Ry==2) then
                  do m=1,n_group
                     j_net_y_3D(m,lp0,k)=dfn(m,lp0,k)*flux(l,k,m)
                     !j_net_y_3D(m,lp0,k)=(dfn(m,lp0,k)-dnn(m,lp0,k))*flux(l,k,m)
                  enddo
               endif
#ifdef js_mpi
904            continue
#endif
            endif

         enddo
      enddo

      do l=1,nxy
         l0=l
#ifdef js_mpi
         if (comm%usempi) then
            if (iproc/=ixy2iproc(l)) cycle
            l0=ixyip2ixy(l,iproc)
         endif
#endif
         if (idomz==1) then
            k=1
            if (BC_Lz==1) then
               do m=1,n_group
                  j_net_z_3D(m,l0,k)=0d0
               enddo
            elseif (BC_Lz==3.or.BC_Lz==2) then
               do m=1,n_group
                  j_net_z_3D(m,l0,k)=-(dfb(m,l0,k)+dnb(m,l0,k))*flux(l,k,m)
               enddo
            endif
         endif
         if (idomz==ndomz) then
            k=nz
            kp=k+1
            if (BC_Rz==1) then
               do m=1,n_group
                  j_net_z_3D(m,l0,kp)=0d0
               enddo
            elseif (BC_Rz==3.or.BC_Rz==2) then
               do m=1,n_group
                  j_net_z_3D(m,l0,kp)=(dfb(m,l0,kp)-dnb(m,l0,kp))*flux(l,k,m)
               enddo
            endif
         endif
         do k=2,nz
            km=k-1
            do m=1,n_group
               j_net_z_3D(m,l0,k)=-dfb(m,l0,k)*(flux(l,k,m)-flux(l,km,m)) &
                                & -dnb(m,l0,k)*(flux(l,k,m)+flux(l,km,m))
            enddo
         enddo
      enddo

#ifdef js_sp3
      if (opt_nodal==3) then
         call updcur2
      endif
#endif

      return
      end subroutine updcurnxy_2d1d
#endif
#ifdef js_sp3
      subroutine updcur2
      use inc_sp3, only: j2_net_x_3d, j2_net_y_3d, j2_net_z_3d
      !use inc_sp3, only: dnb2
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: ixy2iproc, ixyip2ixy
      use inc_parallel, only: nodelmpi
#endif
      implicit none
      integer :: i, j, k, m, l, lm, lp, kp, km1
      integer :: isub
      integer :: kend, kbeg, kbot, ktop
      integer :: l0, l00, lm0, lp0


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [updcur2] in Mod_Soldhat'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif

      isub=1
      kbot=1
      ktop=nz
      do k=kbot,ktop
         do j=1,ny
            do i=Ix_Start_y(j)+1,Ix_End_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i-1,j)
                  lp=nodel(i+1,j)
                  lp0=nodelmpi(i+1,j)
               endif
#endif
               do m=1,n_group
                  j2_net_x_3D(m,l00,k)=-2d0*dfw(m,l00,k)*(flux2(l,k,m)-flux2(lm,k,m))
#ifdef js_mpi
                  if (comm%usempi) then
                     j2_net_x_3D(m,lp0,k)=-2d0*dfw(m,lp0,k)*(flux2(lp,k,m)-flux2(l,k,m))
                  endif
#endif
               enddo
            enddo
            if (idomx==1) then
               i=Ix_Start_y(j)
               l=nodel(i,j)
               lm=nodel(i-1,j)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 901
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i-1,j)
               endif
#endif
               if (BC_Lx==1) then
                  do m=1,n_group
                     j2_net_x_3D(m,l00,k)=0d0
                  enddo
               elseif (BC_Lx==3.or.BC_Lx==2) then
                  do m=1,n_group
                     j2_net_x_3D(m,l00,k)=-2d0*dfw(m,l00,k)*flux2(l,k,m)
                  enddo
               endif
#ifdef js_mpi
901            continue
#endif
            endif
            if (idomx==ndomx) then
               i=Ix_End_y(j)
               l=nodel(i,j)
               lp=nodel(i+1,j)
               l0=l
               l00=l
               lp0=lp
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 902
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lp0=nodelmpi(i+1,j)
               endif
#endif
               if (BC_Rx==1) then
                  do m=1,n_group
                     j2_net_x_3D(m,lp0,k)=0d0
                  enddo
               elseif (BC_Rx==3.or.BC_Rx==2) then
                  do m=1,n_group
                     j2_net_x_3D(m,lp0,k)=2d0*dfw(m,lp0,k)*flux2(l,k,m)
                  enddo
               endif
#ifdef js_mpi
902            continue
#endif
            endif

         enddo

         do i=1,nx
            do j=Iy_Start_x(i)+1,Iy_End_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i,j-1)
                  lp=nodel(i,j+1)
                  lp0=nodelmpi(i,j+1)
               endif
#endif
               do m=1,n_group
                  j2_net_y_3D(m,l00,k)=-2d0*dfn(m,l00,k)*(flux2(l,k,m)-flux2(lm,k,m))
#ifdef js_mpi
                  if (comm%usempi) then
                     j2_net_y_3D(m,lp0,k)=-2d0*dfn(m,lp0,k)*(flux2(lp,k,m)-flux2(l,k,m))
                  endif
#endif
               enddo
            enddo
            if (idomy==1) then
               j=Iy_Start_x(i)
               l=nodel(i,j)
               lm=nodel(i,j-1)
               l0=l
               l00=l
               lm0=lm
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 903
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lm0=nodelmpi(i,j-1)
               endif
#endif
               if (BC_Ly==1) then
                  do m=1,n_group
                     j2_net_y_3D(m,l00,k)=0d0
                  enddo
               elseif (BC_Ly==3.or.BC_Ly==2) then
                  do m=1,n_group
                     j2_net_y_3D(m,l00,k)=-2d0*dfn(m,l00,k)*flux2(l,k,m)
                  enddo
               endif
#ifdef js_mpi
903            continue
#endif
            endif
            if (idomy==ndomy) then
               j=Iy_End_x(i)
               l=nodel(i,j)
               lp=nodel(i,j+1)
               l0=l
               l00=l
               lp0=lp
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) goto 904
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  lp0=nodelmpi(i,j+1)
               endif
#endif
               if (BC_Ry==1) then
                  do m=1,n_group
                     j2_net_y_3D(m,lp0,k)=0d0
                  enddo
               elseif (BC_Ry==3.or.BC_Ry==2) then
                  do m=1,n_group
                     j2_net_y_3D(m,lp0,k)=2d0*dfn(m,lp0,k)*flux2(l,k,m)
                  enddo
               endif
#ifdef js_mpi
904            continue
#endif
            endif
         enddo
      enddo

      kbeg=kbot
      kend=ktop+1
      if (isub==1) then
         k=1
         if (BC_Lz==1) then
            do l=1,nxy
               l0=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
               endif
#endif
               do m=1,n_group
                  j2_net_z_3D(m,l0,k)=0d0
               enddo
            enddo
         elseif (BC_Lz==3.or.BC_Lz==2) then
            do l=1,nxy
               l0=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
               endif
#endif
               do m=1,n_group
                  j2_net_z_3D(m,l0,k)=-2d0*dfb(m,l0,k)*flux2(l,k,m)
                  !j2_net_z_3D(m,l0,k)=-(2d0*dfb(m,l0,k)+dnb2(m,l0,k))*flux2(l,k,m)
               enddo
            enddo
         endif
         kbeg=2
      endif
      if (isub==ndomz) then
         k=nz
         kp=k+1
         if (BC_Rz==1) then
            do l=1,nxy
               l0=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
               endif
#endif
               do m=1,n_group
                  j2_net_z_3D(m,l0,kp)=0d0
               enddo
            enddo
         elseif (BC_Rz==3.or.BC_Rz==2) then
            do l=1,nxy
               l0=l
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
               endif
#endif
               do m=1,n_group
                  j2_net_z_3D(m,l0,kp)=2d0*dfb(m,l0,kp)*flux2(l,k,m)
                  !j2_net_z_3D(m,l0,kp)=(2d0*dfb(m,l0,kp)-dnb2(m,l0,kp))*flux2(l,k,m)
               enddo
            enddo
         endif
         kend=nz
      endif
      do k=kbeg,kend
         km1=k-1
         do l=1,nxy
            l0=l
#ifdef js_mpi
            if (comm%usempi) then
               if (iproc/=ixy2iproc(l)) cycle
               l0=ixyip2ixy(l,iproc)
            endif
#endif
            do m=1,n_group
               j2_net_z_3D(m,l0,k)=-2d0*dfb(m,l0,k)*(flux2(l,k,m)-flux2(l,km1,m))
               !j2_net_z_3D(m,l0,k)=-2d0*dfb(m,l0,k)*(flux2(l,k,m)-flux2(l,km1,m)) &
               !                       -dnb2(m,l0,k)*(flux2(l,k,m)+flux2(l,km1,m))
            enddo
         enddo
      enddo

      return
      end subroutine updcur2
#endif

      END MODULE Mod_Soldhat

#endif
!
!