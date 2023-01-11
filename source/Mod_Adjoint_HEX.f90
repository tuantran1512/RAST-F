#ifdef tuan_tr_test
      MODULE Mod_Adjoint_HEX
  
      implicit none
      
      contains
 
      subroutine Adjoint_hex
!      use Mod_SolLU,    only: ilufac2d
      use Mod_SStpen,   only: IluFacHex
      use Mod_SolLS,    only: setls!, SolveLinearSystem
      use Inc_File,     only: name_adj, W_ADJ
      use Inc_Control,  only: iout_Max, Period_CMFD, EPS_keff_St
      use Inc_Control,  only: eigvd, reigv, reigvs, reigvsd, eshift, eshift0
      use Inc_Option,   only: n_group
      use Inc_Geometry, only: nz, nxy, nxy_1N
      use Inc_Geometry, only: izFuelTop, izFuelBot
      use Inc_Geometry, only: nodeVolume, tot_vol
      use Inc_Geometry, only: LtoLFA
      use Inc_3D,       only: keff
      use Inc_3D,       only: flux, flux_old, FisSrc, flux_adj
      use Inc_RP,       only: I_FA, nxy_FA, nxy_FA_1N
      use Inc_maXS,     only: maxs_chi_3d
      use Inc_Lscoef,   only: af
      use Inc_FluxVar,  only: src
      use Inc_Geometry, only: meshSize_z
      use Inc_Geometry, only: ny_1n
      use Inc_XYZ,      only: fFlux_adj_XYZ, tFlux_adj_XYZ
      use Inc_XYZ,      only: fFlux_adj_XY, tFlux_adj_XY
      use Inc_Geometry, only: ix_startFA_y_1N
      use Inc_Geometry, only: ix_endFA_y_1N
      use Inc_Utile,    only: blank_line
      use Inc_Utile,    only: div_line2, div_line
      use Mod_Alloc
      use Inc_TPEN
      use Inc_maXS
      use Mod_Soltpen,   only: setlshex, solvelinearsystem_hex
      implicit none
      real(8), allocatable :: psia(:,:,:) 
      real(8), allocatable :: psiad(:,:,:)
      real(8), allocatable :: psiri(:,:)  
      real(8), allocatable :: volasy(:)   
      real(8), allocatable :: volax(:)    
      integer :: ixy_fa, ixy
      integer :: M2D(2)
      real(8) :: eigvf, eigvs, vol, reigvdel, fast
      real(8) :: fsadj, gamman, gammad, gamma, sum_adj, fnorm, ther
      integer :: k, l, i, m
      integer :: iz, iin, lfa, iout
      real(8) :: tmp1, tmp2, tmp3, tmp4
      integer :: iy_1n
      character(200) :: FMT
      real(8) :: valt
      integer(4) :: ig, ii, it, km1, ih1, ih2
      
      open(W_ADJ, FILE = name_adj)
      
      ! store forward flux
      if (.not.allocated(psia  )) call alloc (psia  ,n_group,nxy,nz)
      if (.not.allocated(psiad )) call alloc (psiad ,n_group,nxy,nz)
      if (.not.allocated(psiri )) call alloc (psiri ,nxy,nz)
      if (.not.allocated(volasy)) call alloc (volasy,nxy_1N)
      if (.not.allocated(volax )) call alloc (volax ,nxy)
      write(W_ADJ,*) Div_Line2(1:70)
      write(W_ADJ,*) "Adjoint Calculation Data"
      write(W_ADJ,*) ""
      write(W_ADJ,*) "[ Adjoint iteration table ]"
      write(W_ADJ,*) "[ Axially Integrated Adjoint Flux Distribution ]"
      write(W_ADJ,*) "[ Radially Integrated Adjoint Flux Distribution ]"
      write(W_ADJ,*) "[ 3-D Adjoint Flux Distribution ]"
      write(W_ADJ,*) ""
      
      eigvf=keff
      do m=1,n_group
         do k=1,nz
            do l=1,nxy
               flux_old(l,k,m)=flux(l,k,m)
            enddo
         enddo
      enddo

      ! reset the forward linear system
      reigvs = 0d0
      CALL setlshex(1)
      
      keff  =1d0
      eigvs =1d0+eshift0
      reigv =1d0/keff
      reigvs=1d0/eigvs
      !! For hexageometry problem
      DO ih1=1,nxy
         DO it=1,ipntr(0,ih1)
            ih2=ipntr(it,ih1)
            IF(ih2.GT.ih1) THEN
               DO ii=1,ipntr(0,ih2)
                  IF(ipntr(ii,ih2).EQ.ih1) EXIT
               ENDDO
               DO iz=1,nz
                  DO ig=1,n_group
                     valt=cmat(ig,it,ih1,iz)
                     cmat(ig,it,ih1,iz)=cmat(ig,ii,ih2,iz)
                     cmat(ig,ii,ih2,iz)=valt
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      
      DO k=2,nz
         km1=k-1
         DO l=1,nxy
            DO m=1,n_group
               valt=cmat(m,8,l,km1)
               cmat(m,8,l,km1)=cmat(m,7,l,k)
               cmat(m,7,l,k)=valt
            ENDDO
         ENDDO
      ENDDO
      DO k=1,nz
         DO l=1,nxy
            dcmat(1,l,k)=maXs_r_3D(imap(l),k,1)-nu_maXs_f_3D(imap(l),k,1)*reigvs &
                 +dsum(1,l,k)
            !dcmat(2,l,k)=-maXS_scat_3D(1, 2, imap(l), k)
            dcmat(2,l,k)=-maXs_s_3D(imap(l),k,1)
            dcmat(3,l,k)=-nu_maXs_f_3D(imap(l),k,2)*reigvs
            dcmat(4,l,k)=maXs_r_3D(imap(l),k,2)+dsum(2,l,k)
            vol=nodevolume(imap(l),k)
            DO i=1,4
               dcmat(i,l,k)=dcmat(i,l,k)*vol
            ENDDO
         ENDDO
      ENDDO
 
      ! build new preconditioner
      call IluFacHex 
  
      ! outer iteration (ioutp)
      write(W_ADJ,*) Div_Line2(1:70)
      WRITE(W_ADJ,*) "[ Adjoint iteration table ]"
      WRITE(W_ADJ,*) ""
      WRITE(W_ADJ,'(a,4x,2(a,7x),a)') " Iteration Eigenvalue","Gamma_a",  &
           "Gamma_d","Gamma_n"
      
      !! === ITERATION === !!
      do iout=1,iout_max
         ! source generation
         reigvdel=reigv-reigvs
         do k=1,nz
            do l=1,nxy
               !! === indexing ====== !!
               !! af   : hex index    !!
               !! src  : square index !!
               !! psia : square index !!
               !! =================== !!
               fsadj=maxs_chi_3d(l,k,1)*flux(l,k,1) &
                    +maxs_chi_3d(l,k,2)*flux(l,k,2)
               psia(1,l,k)=af(1,imapsol(l),k)*fsadj
               psia(2,l,k)=af(2,imapsol(l),k)*fsadj
               src(1,l,k)=psia(1,l,k)*reigvdel
               src(2,l,k)=psia(2,l,k)*reigvdel
            enddo
         enddo

         call SolveLinearSystem_hex(iin,.FALSE.) 
         
         ! update adjoint fission source
         gamman = 0d0
         gammad = 0d0
         do ixy_fa=1,nxy_fa
            ixy=i_fa(ixy_fa)
            do iz=izFuelBot,IzFuelTop
               psiad(1,ixy,iz)=psia(1,ixy,iz)
               psiad(2,ixy,iz)=psia(2,ixy,iz)
               fsadj=maxs_chi_3d(ixy,iz,1)*flux(ixy,iz,1) &
                  & +maxs_chi_3d(ixy,iz,2)*flux(ixy,iz,2)
               psia(1,ixy,iz)=af(1,imapsol(ixy),iz)*fsadj
               psia(2,ixy,iz)=af(2,imapsol(ixy),iz)*fsadj
               gamman=gamman+psia(1,ixy,iz)*psia(1,ixy,iz) &
                     &      +psia(2,ixy,iz)*psia(2,ixy,iz)
               gammad=gammad+psiad(1,ixy,iz)*psia(1,ixy,iz) &
                     &      +psiad(2,ixy,iz)*psia(2,ixy,iz)
            enddo
         enddo
         ! update adjoint eigenvalue
         eigvd=keff
         gamma=gammad/gamman
         keff=1/(reigv*gamma+(1-gamma)*reigvs)
         if (iout<Period_CMFD) THEN
            eigvs=keff+eshift0
         else
            eigvs=keff+eshift
         endif
         reigv=1d0/keff
         reigvsd=reigvs
         reigvs=1d0/eigvs
         
         ! update linear system
         reigvdel=reigvs-reigvsd
         !! l = hex index 
         !!@@ for lattice
         !!do k=1,nz
         !!   do l=1,nxy
         !!      !! === index ====================== !!
         !!      !! am : hex index                   !!
         !!      !! af : square index                !!
         !!      !! scat: no need in hex calculation !!
         !!      !! ================================ !!
         !!      am(1,l,k)=am(1,l,k)-af(1,l,k)*reigvdel
         !!      scat(l,k)=af(2,l,k)*reigvs
         !!   enddo
         !!enddo
         DO k=1,nz
            DO l=1,nxy
               vol=nodevolume(imap(l),k)
               dcmat(1,l,k)=(maXs_r_3D(imap(l),k,1)+dsum(1,l,k) &
                    -nu_maXs_f_3D(imap(l),k,1)*reigvs*maXs_chi_3d(imap(l),k,1))*vol
         
               dcmat(2,l,k)=(-maXs_s_3D(imap(l),k,1) &
                    -nu_maXs_f_3D(imap(l),k,1)*reigvs*maXs_chi_3D(imap(l),k,2))*vol
               dcmat(3,l,k)=-nu_maXs_f_3D(imap(l),k,2)*reigvs*maXs_chi_3D(imap(l),k,1)*vol
               dcmat(4,l,k)=(maXs_r_3D(imap(l),k,2)+dsum(2,l,k) &
                    -nu_maXs_f_3D(imap(l),k,2)*reigvs*maXs_chi_3D(imap(l),k,2))*vol
            ENDDO
         ENDDO
         WRITE(W_ADJ,'(i4,f15.7,1p,4e15.6)') iout, keff, gamma, gammad, gamman
         if (abs(keff-eigvd)<=EPS_keff_St*0.1d0) exit
      enddo
      !! === END ITERATION === !!
      WRITE(W_ADJ,'(a,1x,f8.6,a14,1x,i4,a12)') &
           & " Adjoint k-eff= ",keff, " Obtained in ",iout,"  Iterations"
      WRITE(W_ADJ,*) ""
      ! store adjoint flux
      sum_adj=0d0
      do k=1,nz
         do l = 1,nxy
            flux_adj(l,k,1)=flux(l,k,1)
            flux_adj(l,k,2)=flux(l,k,2)
            flux(l,k,1)=flux_old(l,k,1)
            flux(l,k,2)=flux_old(l,k,2)
            sum_adj=sum_adj+FisSrc(l,k)*flux_adj(l,k,1)
         enddo
      enddo
      keff=eigvf
      reigv=1d0/keff
      
      fnorm=Tot_Vol/sum_adj
      do lfa=1,nxy_FA_1N
         psiri(lfa,3)=0d0
         psiri(lfa,4)=0d0
         volasy(lfa)=0d0
      enddo
      volax=0d0
  
      do ixy_fa=1,nxy_fa
         ixy=I_FA(ixy_fa)
         lfa=ltolfa(ixy)
         psiri(ixy,1)=0d0
         psiri(ixy,2)=0d0
         do iz=izFuelBot,izFuelTop    
            vol=nodeVolume(ixy,iz)
            psiri(ixy,1)=psiri(ixy,1)+flux_adj(ixy,iz,1)*vol
            psiri(ixy,2)=psiri(ixy,2)+flux_adj(ixy,iz,2)*vol
            psiri(lfa,3)=psiri(lfa,3)+flux_adj(ixy,iz,1)*vol
            psiri(lfa,4)=psiri(lfa,4)+flux_adj(ixy,iz,2)*vol
            volasy(lfa)=volasy(lfa)+vol
            volax(ixy)=volax(ixy)+vol
         enddo
      enddo
      
      ! assembly-wise power
      do lfa=1,nxy_FA_1N
         psiri(lfa,3)=psiri(lfa,3)/volasy(lfa)*fnorm
         psiri(lfa,4)=psiri(lfa,4)/volasy(lfa)*fnorm
      enddo
      
      write(W_ADJ,*) Div_Line2(1:70)
      write(W_ADJ,*) "[ Axially Integrated Adjoint Flux Distribution ]"
      write(W_ADJ,*) ""
      call ADJ_XYZ_Indexing
      
      ! Print 2D Radial Adjoint Flux Distribution
      M2D = MAXLOC( fFlux_adj_XY(:,:) )
      write(W_ADJ,*) "2-D Fast Adjoint Flux Distribution"
      write(W_ADJ,*) ""
      write(W_ADJ,'(A,ES15.5, A9 )') " Peak Fast Flux = ", MAXVAL( fFlux_adj_XY(:,:) ) , "[#/cm3]"         
      write(W_ADJ,'(A,I11, A2, I2, A/)') " Peak Position  = ", M2D(2), ", ", M2D(1), "  (Row, Column)"
  
      do Iy_1N=1,Ny_1N
         if (Ix_StartFA_y_1N(Iy_1N)>0) then
            if (Ix_StartFA_y_1N(Iy_1N)==1) then
               write(FMT,*) "( ", Ix_EndFA_y_1N(Iy_1N) - Ix_StartFA_y_1N(Iy_1N) + 1, " Es12.5 )"
               write(W_ADJ,FMT) fFlux_adj_XY(Ix_StartFA_y_1N(Iy_1N):Ix_EndFA_y_1N(Iy_1N),Iy_1N)
            else
               write(FMT,*) "( ", Ix_StartFA_y_1N(Iy_1N) - 1, " A12, ", &
                  & Ix_EndFA_y_1N(Iy_1N) - Ix_StartFA_y_1N(Iy_1N) + 1, " Es12.5 )"
               write(W_ADJ,FMT) Blank_Line( 1 : Ix_StartFA_y_1N(Iy_1N) - 1 ),  &
                 fFlux_adj_XY( Ix_StartFA_y_1N(Iy_1N):Ix_EndFA_y_1N(Iy_1N), Iy_1N )
            endif
         endif
      enddo
      WRITE(W_ADJ,*) ""
      write(W_ADJ,*) Div_Line(1:70)
      
      M2D = MAXLOC( tFlux_adj_XY(:,:) )
      write(W_ADJ,*) "2-D Thermal Adjoint Flux Distribution"
      write(W_ADJ,*) ""
      write(W_ADJ,'(A,ES15.5, A9 )') " Peak Thermal Flux = ", MAXVAL( tFlux_adj_XY(:,:) ) , "[#/cm3]"         
      write(W_ADJ,'(A,I11, A2, I2, A/)') " Peak Position     = ", M2D(2), ", ", M2D(1), "  (Row, Column)"
      
      do Iy_1N=1,Ny_1N
         if (Ix_StartFA_y_1N(Iy_1N)>0) then
            if (Ix_StartFA_y_1N(Iy_1N)==1) then
               write(FMT,*) "( ", Ix_EndFA_y_1N(Iy_1N) - Ix_StartFA_y_1N(Iy_1N) + 1, " Es12.5 )"
               write(W_ADJ,FMT) tFlux_adj_XY(Ix_StartFA_y_1N(Iy_1N):Ix_EndFA_y_1N(Iy_1N),Iy_1N)
            else
               write(FMT,*) "( ", Ix_StartFA_y_1N(Iy_1N) - 1, " A12, ", &
                  & Ix_EndFA_y_1N(Iy_1N) - Ix_StartFA_y_1N(Iy_1N) + 1, " Es12.5 )"
               write(W_ADJ,FMT) Blank_Line( 1 : Ix_StartFA_y_1N(Iy_1N) - 1 ),  &
                 tFlux_adj_XY( Ix_StartFA_y_1N(Iy_1N):Ix_EndFA_y_1N(Iy_1N), Iy_1N )
            endif
         endif
      enddo
      write(W_ADJ,*) ""
      
      ! Print 1D Axial Adjoint Flux Distribution
      write(W_ADJ,*) Div_Line2(1:70)
      write(W_ADJ,*) "[ Radially Integrated Adjoint Flux Distribution ]"
      write(W_ADJ,*) ""
      write(W_ADJ,*) " Iz  Height[cm]  Height[%]   Fast Flux    Thermal Flux"
      tmp3=sum(meshSize_z(1:nz))
      tmp1=0d0
      do k=1,nz
         tmp1=tmp1+meshSize_z(k)
         tmp2=meshSize_z(k)/2.d0
         tmp4=(tmp1-tmp2)/tmp3*100.d0
         
         fast=0.d0
         ther=0.d0
         do l=1,nxy
            fast=fast+flux_adj(l,k,1)
            ther=ther+flux_adj(l,k,2) 
         enddo
         fast=fast/nxy
         ther=ther/nxy
         WRITE(W_ADJ,'(i3, F11.3, F11.2, 2Es15.5)') k, tmp1-tmp2, tmp4, fast, ther
      enddo
      write(W_ADJ,*) ""
      
      ! Print radial adjoint flux distirubiton
      write(W_ADJ,*) Div_Line2(1:70)
      write(W_ADJ,*) "[ 3-D Adjoint Flux Distribution ]"
      write(W_ADJ,*) ""
      do iz=IzFuelBot,IzFuelTop
         M2D = MAXLOC( fFlux_adj_XYZ(:,:,iz) )
         write(W_ADJ,*) "3-D Fast Adjoint Flux Distribution"
         write(W_ADJ,*) ""
         write(W_ADJ,'(A,F15.2, A4 )')  " Z bottom       = ", sum(MeshSize_z(1:iz-1)), "cm"   
         write(W_ADJ,'(A,F15.2, A4 )')  " Height         = ", MeshSize_z(iz), "cm"
         write(W_ADJ,'(A,ES15.5, A9 )') " Peak Fast Flux = ", MAXVAL( fFlux_adj_XYZ(:,:,iz) ) , "[#/cm3]"         
         write(W_ADJ,'(A,I11, A2, I2, A/)') " Peak Position  = ", M2D(2), ", ", M2D(1), "  (Row, Column)"
         
         do Iy_1N=1,Ny_1N
            if (Ix_StartFA_y_1N(Iy_1N)>0) then
               if (Ix_StartFA_y_1N(Iy_1N)==1) then
                  write(FMT,*) "( ", Ix_EndFA_y_1N(Iy_1N) - Ix_StartFA_y_1N(Iy_1N) + 1, " Es12.5 )"
                  write(W_ADJ,FMT) fFlux_adj_XYZ( Ix_StartFA_y_1N(Iy_1N):Ix_EndFA_y_1N(Iy_1N), Iy_1N, Iz )
               else
                  write(FMT,*) "( ", Ix_StartFA_y_1N(Iy_1N) - 1, " A12, ", &
                     & Ix_EndFA_y_1N(Iy_1N) - Ix_StartFA_y_1N(Iy_1N) + 1, " Es12.5 )"
                  write(W_ADJ,FMT) Blank_Line( 1 : Ix_StartFA_y_1N(Iy_1N) - 1 ),  &
                     fFlux_adj_XYZ( Ix_StartFA_y_1N(Iy_1N):Ix_EndFA_y_1N(Iy_1N), Iy_1N, Iz )
               endif
            endif
         enddo
         write(W_ADJ,*) ""
         write(W_ADJ,*) Div_Line(1:70)
      enddo
      write(W_ADJ,*) Div_Line2(1:70)
      
      do iz=IzFuelBot,IzFuelTop
         M2D = MAXLOC( tFlux_adj_XYZ(:,:,iz) )    
         write(W_ADJ,*) "3-D Thermal Adjoint Flux Distribution"
         write(W_ADJ,*) ""
         write(W_ADJ,'(A,F15.2, A4 )')  " Z bottom       = ", sum(MeshSize_z(1:iz-1)), "cm"   
         write(W_ADJ,'(A,F15.2, A4 )')  " Height         = ", MeshSize_z(iz), "cm"
         write(W_ADJ,'(A,ES15.5, A9 )') " Peak Fast Flux = ", MAXVAL( fFlux_adj_XYZ(:,:,iz) ) , "[#/cm3]"         
         write(W_ADJ,'(A,I11, A2, I2, A/)') " Peak Position  = ", M2D(2), ", ", M2D(1), "  (Row, Column)"
         
         do Iy_1N=1,Ny_1N
            if (Ix_StartFA_y_1N(Iy_1N)>0) then
               if (Ix_StartFA_y_1N(Iy_1N)==1) then
                  write(FMT,*) "( ", Ix_EndFA_y_1N(Iy_1N) - Ix_StartFA_y_1N(Iy_1N) + 1, " Es12.5 )"
                  write(W_ADJ,FMT) tFlux_adj_XYZ( Ix_StartFA_y_1N(Iy_1N):Ix_EndFA_y_1N(Iy_1N), Iy_1N, Iz )
               else
                  write(FMT,*) "( ", Ix_StartFA_y_1N(Iy_1N) - 1, " A12, ", &
                     & Ix_EndFA_y_1N(Iy_1N) - Ix_StartFA_y_1N(Iy_1N) + 1, " Es12.5 )"
                  write(W_ADJ,FMT) Blank_Line( 1 : Ix_StartFA_y_1N(Iy_1N) - 1 ),  &
                     tFlux_adj_XYZ( Ix_StartFA_y_1N(Iy_1N):Ix_EndFA_y_1N(Iy_1N), Iy_1N, Iz )
               endif
            endif
         enddo
         write(W_ADJ,*) ""
         write(W_ADJ,*) Div_Line(1:70)
      enddo
      
      close(W_ADJ)
      write(*,*) "*** Adjoint Calculation ... OK"
      return 
      end subroutine Adjoint_hex

      subroutine ADJ_XYZ_Indexing
      use Inc_3D
      use Inc_File
      use Inc_INP
      use Inc_Option
      use Inc_XYZ
      use Inc_RP
      use Inc_Geometry
      implicit none
      real(8) :: FoldFactor
      real(8) :: FoldFactor_1N
      real(8) :: sum_real(2)
      integer :: ixy_1N, iz, ix_1N, iy_1N
      integer :: FoldCount, m

      FoldCount=0
      FoldFactor=0d0
      FoldFactor_1N=0d0
      do iy_1N=1,Ny_1N
         do ix_1N=Ix_Start_y_1N(iy_1N),Ix_End_y_1N(iy_1N)
            ixy_1N=ixiy_1NToixy_1N(ix_1N,iy_1N)
            do iz=1,Nz
               Sum_Real=0d0
               FoldCount=0
               do m=1,4
                  if (i_1Nto4N(ixy_1N,m)==0) then
                      FoldCount=FoldCount+1
                      cycle
                  else
                  sum_real(1)=sum_real(1)+flux_adj(i_1Nto4N(ixy_1n,m),iz,1)
                  sum_real(2)=sum_real(2)+flux_adj(i_1Nto4N(ixy_1n,m),iz,2)
                  endif
               enddo
               FoldFactor=4d0/(4d0-REAL(FoldCount,8))
               fFlux_adj_XYZ(ix_1N,iy_1N,iz)=FoldFactor*(Sum_Real(1)/4d0)
               tFlux_adj_XYZ(ix_1N,iy_1N,iz)=FoldFactor*(Sum_Real(2)/4d0)
            ENDDO
         ENDDO
      ENDDO

      do iy_1N=1,Ny_1N
         do ix_1N=Ix_Start_y_1N(iy_1N),Ix_End_y_1N(iy_1N)
             fFlux_adj_XY(ix_1N,iy_1N)=sum(fFlux_adj_XYZ(ix_1N,iy_1N,1:Nz))
             tFlux_adj_XY(ix_1N,iy_1N)=sum(tFlux_adj_XYZ(ix_1N,iy_1N,1:Nz))
         enddo
      enddo

      return
      end subroutine ADJ_XYZ_Indexing

      END MODULE Mod_Adjoint_HEX
#endif
 
