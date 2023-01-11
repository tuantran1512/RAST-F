
      module Mod_SolFDM
      use Inc_Constant
      use Inc_Geometry
      use Inc_Option, only: N_Group


#ifdef siarhei_plot
      use Inc_Control,only:last_saved,current_sub,dummy_filler
#endif

      implicit none

      contains

      subroutine fdmcoupl
      use inc_lscoef, only: dfw, dfn, dfb
      use inc_maxs, only: d_3d
#ifdef js_mpi
      use inc_parallel, only: comm, iproc
      use inc_parallel, only: ixy2iproc, ixyip2ixy
      use inc_parallel, only: nodelmpi
#endif
#ifdef tuan_tr_test
      use Mod_SolTPEN, only: DtilHex
!      use Inc_TPEN, only: if_hexgeometry 
#endif
      implicit none
      real(8) :: buff_real
      real(8) :: hzk, hzkb, hyj, hyjn, hxi, hxiw, albf
      real(8) :: difl, diflw, difln, diflb
      integer :: k, kb, j, i, l, lw, ln, m, ib, ie, le, js, ls, kt
      integer :: l0, l00, le0, ls0
#ifdef js_mpi
      real(8) :: difle, difls, hxie, hyjs
#endif


#ifdef siarhei_plot
         ! @$^ siarhei_fr
         current_sub = &
            '+#+ Entered [fdmcoupl] in Mod_SolFDM'
         if(trim(last_saved).ne.trim(current_sub)) then
            write(678,'(A)') trim(current_sub)
            last_saved = current_sub
         end if
         ! =-=-=-=-=-=-=-=
#endif
#ifdef tuan_tr_test
!      if (if_hexgeometry) then
         call DtilHex
         return
!      endif
#endif
      

      do k=1,nz
         kb=k-1
         hzk=MeshSize_z(k)
         hzkb=MeshSize_z(kb)
         if (k==1.and.BC_Lz==1.and.hzk/=hzkb) hzkb=hzk
         do j=1,ny
            hyj=MeshSize_y(j)
            hyjn=MeshSize_y(j-1)
            if (j==1.and.BC_Ly==1.and.hyj/=hyjn) hyjn=hyj
            do i=ix_start_y(j),ix_end_y(j)
               hxi=MeshSize_x(i)
               hxiw=MeshSize_x(i-1)
               if (i==ix_start_y(j).and.BC_Lx==1.and.hxi/=hxiw) hxiw=hxi
               l=nodel(i,j)
               lw=nodel(i-1,j)
               ln=nodel(i,j-1)
               l0=l
               l00=l
#ifdef js_mpi
               le=nodel(i+1,j)
               ls=nodel(i,j+1)
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  l0=ixyip2ixy(l,iproc)
                  l00=nodelmpi(i,j)
                  le0=nodelmpi(i+1,j)
                  ls0=nodelmpi(i,j+1)
                  hxie=meshsize_x(i+1)
                  hyjs=meshsize_y(j+1)
               endif
#endif
               do m=1,n_group
                  difl=d_3d(l,k,m)
                  diflw=d_3d(lw,k,m)
                  difln=d_3d(ln,k,m)
                  diflb=d_3d(l,kb,m)
#ifdef js_mpi
                  if (comm%usempi) then
                     difle=d_3d(le,k,m)
                     difls=d_3d(ls,k,m)
                  endif
#endif
                  dfw(m,l00,k)=2d0*difl*diflw/(difl*hxiw+diflw*hxi)
                  dfn(m,l00,k)=2d0*difl*difln/(difl*hyjn+difln*hyj)
                  dfb(m,l0,k)=2d0*difl*diflb/(difl*hzkb+diflb*hzk)
#ifdef js_mpi
                  if (comm%usempi) then
                     dfw(m,le0,k)=2d0*difle*difl/(difle*hxi+difl*hxie)
                     dfn(m,ls0,k)=2d0*difls*difl/(difls*hyj+difl*hyjs)
                  endif
#endif
               enddo
            enddo
         enddo
      enddo
      ! x-direct
      do ib=1,2
         do k=1,nz
            do j=1,ny
               if (ib==1) then
                  if (BC_Lx==1) then
                     albf=0d0
                  elseif (BC_Lx==2) then
                     albf=half
                  elseif (BC_Lx==3) then
                     albf=big
                  endif
                  i=Ix_Start_y(j)
                  ie=Ix_Start_y(j)
               else !ib==2
                  if (BC_Rx==1) then
                     albf=0d0
                  elseif (BC_Rx==2) then
                     albf=half
                  elseif (BC_Rx==3) then
                     albf=big
                  endif
                  i=Ix_End_y(j)
                  ie=i+1
               endif
               hxi=MeshSize_x(i)
               l=nodel(i,j)
               le=nodel(ie,j)
               le0=le
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  le0=nodelmpi(ie,j)
               endif
#endif
               do m=1,n_group
                  buff_real=d_3d(l,k,m)
                  dfw(m,le0,k)=bcf(buff_real,albf,hxi) ! check later siarhei_rev
!            dummy_filler = 1 ! @$^ siarhei_plot 
               enddo
            enddo
         enddo
      enddo
      ! y-direct
      do ib=1,2
         do k=1,nz
            do i=1,nx
               if (ib==1) then
                  if (BC_Ly==1) then
                     albf=0d0
                  elseif (BC_Ly==2) then
                     albf=half
                  elseif (BC_Ly==3) then
                     albf=big
                  endif
                  j=Iy_Start_x(i)
                  js=Iy_Start_x(i)
               else !ib==2
                  if (BC_Ry==1) then
                     albf=0d0
                  elseif (BC_Ry==2) then
                     albf=half
                  elseif (BC_Ry==3) then
                     albf=big
                  endif
                  j=Iy_End_x(i)
                  js=j+1
               endif
               hyj=MeshSize_y(j)
               l=nodel(i,j)
               ls=nodel(i,js)
               ls0=ls
#ifdef js_mpi
               if (comm%usempi) then
                  if (iproc/=ixy2iproc(l)) cycle
                  ls0=nodelmpi(i,js)
               endif
#endif
               do m=1,n_group
                  buff_real=d_3d(l,k,m)
                  dfn(m,ls0,k)=bcf(buff_real,albf,hyj) ! check later siarhei_rev
         !   dummy_filler = 1 ! @$^ siarhei_plot 
               enddo
            enddo
         enddo
      enddo
      ! - z-direct
      do ib=1,2
         if (ib==1) then
            if (BC_Lz==1) then
               albf=0d0
            elseif (BC_Lz==2) then
               albf=half
            elseif (BC_Lz==3) then
               albf=big
            endif
            k=1
            kt=1
         else !ib==2
            if (BC_Rz==1) then
               albf=0d0
            elseif (BC_Rz==2) then
               albf=half
            elseif (BC_Rz==3) then
               albf=big
            endif
            k=nz
            kt=k+1
         endif
         hzk=MeshSize_z(k)
         do l=1,nxy
            l0=l
#ifdef js_mpi
            if (comm%usempi) then
               if (iproc/=ixy2iproc(l)) cycle
               l0=ixyip2ixy(l,iproc)
            endif
#endif
            do m=1,n_group
               buff_real=d_3d(l,k,m)
               dfb(m,l0,kt)=bcf(buff_real,albf,hzk) ! check later siarhei_rev
          !  dummy_filler = 1 ! @$^ siarhei_plot 
            enddo
         enddo
      enddo

      return
      end subroutine fdmcoupl


!!!#ifdef siarhei_delete  ! check later siarhei_rev
      function bcf(dif,albf,h)
      ! function to determine the contribution of boudary condition to the diagonal
      ! elements for the boundary nodes
      ! albf - albedo (jnet/phi)
      !  = 0     for zero net current
      !  = inf   for zero flux
      !  = 0.5   for zero incoming current
      implicit none
      real(8) :: bcf, dif, albf, h, alb

      alb=albf/dif
      bcf=2*dif*alb/(2+alb*h)
!            dummy_filler = 1 ! @$^ siarhei_plot 

      return
      end function bcf
!!!#endif 

      end module Mod_SolFDM
