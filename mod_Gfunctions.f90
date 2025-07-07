module GreensFunctions
  use DefineHamiltonian
  use OMP_LIB
  implicit none
  integer :: INFO

  real*8, parameter :: epsilon = 1e-9
  real*8, allocatable, dimension(:) :: Hub, omega, Ev
  real*8 :: pulay, dw, up
  
  complex*16, allocatable, dimension(:,:) :: GammaL, GammaR, Eigenvec, G_nil
!  complex*16, allocatable, dimension(:,:) ::  SigmaL, Sigma1, SigmaR
  complex*16, allocatable, dimension(:,:) :: work1, work2, work3, work4
  
  logical :: verb
  type :: GF
     complex*16, allocatable, dimension(:,:,:) :: r, a, L, G
  end type GF
  type(GF) :: GF0
  
  type :: GF_full
     complex*16, allocatable, dimension(:,:,:) :: R, A, L, G
  end type GF_full
  type(GF_full) :: GFf

contains 
  
   !.....fermi-dirac distribution
  real*8 function fermi_dist(w, V)
    implicit none
    real*8 ::  arg, w, V
    arg = (w - mu + V/hbar)*beta
    fermi_dist =  1.d0/(exp(hbar*arg) +1.d0)
  end function fermi_dist

!====================================================
!========== Self consistency field calculations =====
!====================================================

!.............need G0f%L to calculate GL_of_0
!.............Calculates GFs at every omega and Voltage simultaneously 
subroutine SCF_GFs(Volt, first)
  implicit none
  integer :: iw, iteration,i, Vname
  real*8 :: Volt, err, diff, st, et
  character(len=30) :: fn1
  logical :: first

  iteration = 0

  write(22,*) '........SCF Calculations at Voltage:', Volt, '..........'
  if (first) then 
     st =0.d0; et= 0.d0
     call CPU_TIME(st)
     call G0_R_A()
     call CPU_TIME(et)
     write(22,'(A,F10.8,A,A,A,I4)') 'G0_R_A runtime:', (et-st), 'seconds', '   ', 'Iteration:', iteration
     
     st =0.d0; et= 0.d0
     call CPU_TIME(st)
     call G0_L_G(Volt)
     call CPU_TIME(et)
     write(22,'(A,F10.8,A,A,A,I4)') 'G0_L_G runtime:', (et-st), 'seconds', '   ', 'Iteration:', iteration
  end if
  print *, '>>>>>>>>>>VOLTAGE:', Volt

  if(verb) then
     call print_3matrix(0,GF0%R,'GFR')
     call print_3matrix(0,GF0%L,'GFL')
  end if
  
!........ printing the spectral function at 0 iteration (embedding only) 
!            and calculating spec(1,iw)

  call print_sf(iteration)
  
  Vname = abs(Volt)
  write(fn1,'(i0)') Vname
  if (Volt .ge. 0) then 
     open(17,file='err_V_'//trim(fn1)//'.dat',status='unknown')
  else
     open(17,file='err_V_n'//trim(fn1)//'.dat',status='unknown')
  end if

  DO
     iteration = iteration + 1
     write(*,*) '.... ITERATION = ',iteration,' ....'

!.......real variable interactions turns off the Interaction component of the sigmas 
!.................full Gr and Ga, Eq. (5) and (6)

    ! write(3,'(/a,i3,i5/)') 'iteration = ',iteration

     st =0.d0; et= 0.d0
     call CPU_TIME(st)
     call GL_of_0()
     call CPU_TIME(et)
     write(22,'(A,F10.8,A,A,A,I4)') 'GL_of_0 runtime:', (et-st), 'seconds', '   ','Iteration:', iteration
     
     st = 0.d0; et = 0.d0
     st = OMP_GET_WTIME()
     
     !$OMP PARALLEL DO &
     !$OMP& PRIVATE(iw, INFO)

     do iw = 1, N_of_w
        call G_full(iw, Volt)
     end do
     !$OMP END PARALLEL DO
     
     
     et = OMP_GET_WTIME()
     write(22,'(A,F10.8,A,A,A,I4)') 'G_full w-loop runtime:', (et-st), 'seconds', '   ', 'Iteration:', iteration
     
     
     if(verb) then
        call print_3matrix(iteration,GFf%R,'GFR')
        call print_3matrix(iteration,GFf%L,'GFL')
     end if
     
     !..... calculation of the error

     !$OMP PARALLEL DO PRIVATE(iw, i, diff) REDUCTION(+:err)
     do iw = 1, N_of_w
        do i=1,Natoms
           diff=2.d0*hbar*(AIMAG(GFf%R(i,i,iw))-AIMAG(GF0%R(i,i,iw)))
           err=err +diff*diff
        end do
     end do
     !$OMP END PARALLEL DO
     
     write(*,*) 'err = ',sqrt(err)
     !     write(*,*) iteration, sqrt(err)

     !$OMP CRITICAL
     GF0%R = pulay*GFf%R + (1.0d0-pulay)*GF0%R
     GF0%A = pulay*GFf%A + (1.0d0-pulay)*GF0%A
     GF0%L = pulay*GFf%L + (1.0d0-pulay)*GF0%L
     GF0%G = pulay*GFf%G + (1.0d0-pulay)*GF0%G !GF0%L + GF0%R - GF0%A
     !$OMP END CRITICAL

!... printing the spectral function

     call print_sf(iteration)  
     if (sqrt(err) .lt. epsilon .or. order .eq. 0) then
        write(*,*)'... REACHED REQUIRED ACCURACY ...'
        exit
     end if
  END DO
  close(17)
end subroutine SCF_GFs


!=====================================================
!================== Full GFs =========================
!===================================================== 

subroutine G_full(iw, Volt) !... Full Greens function, leaves Retarded and Advanced in the work arrays, application of Eq. (16) and (17), but with the full Sigmas, Eq. (3), (7) and (8) in CHE
  implicit none
  integer :: i, j, iw, sp, sp1, ii, jj
  real*8 :: Volt, w 
  complex*16 :: OmR, SigL, SigG
  complex*16, dimension(Natoms,Natoms) ::  SigmaL, Sigma1, SigmaR, SigmaG,  w1, w2
  
  !............full SigmaR due to interactions Eq. (7)
  SigmaR = (0.d0, 0.d0); SigmaL = (0.d0, 0.d0); SigmaG = (0.d0, 0.d0)
  
  if (order .eq. 1) then
     call first_order_sigma(Sigma1)
     SigmaR = Sigma1
  else if (order .eq. 2) then
     call first_order_sigma(Sigma1)

     do i = 1, Natoms, 2
        do j = 1, Natoms, 2
           
           do sp = 0, 1
              do sp1 = 0, 1
                 ii = i +sp; jj= j + sp1
                 
                 OmR = Omega_r(i,j,sp,sp1,iw)
                 
                 SigmaR(ii,jj) = SigmaR(ii,jj)+ Sigma1(ii,jj)  + (Hub(j)*Hub(i)*OmR)*hbar**2 
              end do
           end do
           
        end do
     end do

     !..............full SigmaL, Eq. (3) and (4)     
     !.....Interaction contribution of both Sigmas     
     do i = 1, Natoms, 2
        do j = 1, Natoms, 2

           do sp = 0, 1
              do sp1 = 0, 1
                 ii = i +sp; jj= j + sp1
                 
                 SigL = int_SigL(i,j, sp, sp1, iw)
                 SigmaL(ii,jj) = SigmaL(ii,jj) + Hub(j)*Hub(i)*SigL*hbar**2

                 SigG = int_SigG(i,j, sp, sp1, iw)
                 SigmaG(ii,jj) = SigmaG(ii,jj) + Hub(j)*Hub(i)*SigG*hbar**2
              end do
           end do
          
        end do
     end do
  end if
  
  !............full Gr and Ga, Eq. (5) and (6)
  
  w = omega(iw)
  w1 = -H + 0.5d0*(im/hbar)*(GammaL + GammaR) - SigmaR    

  do i = 1 , Natoms
     w1(i,i) = w1(i,i) + hbar*(w +im*0.01)
  end do

  call Inverse_complex(Natoms, w1, info)
  call Hermitian_Conjg(w1, Natoms, w2) 

  GFf%R(:,:,iw) = w1; GFf%A(:,:,iw) = w2 
  
  !.....Embedding contribution of Sigma
  SigmaL = SigmaL + im*(fermi_dist(w, Volt)*GammaL + fermi_dist(w, 0.d0)*GammaR)/hbar
  SigmaG = SigmaG + im*((fermi_dist(w, Volt)-1.d0)*GammaL + (fermi_dist(w, 0.d0)-1.d0)*GammaR)/hbar
  
  !.............full GL and GG, Eq. (16) and (17)
  GFf%L(:,:,iw) = matmul(matmul(w1, SigmaL), w2) !.. GL = Gr * SigmaL * Ga
 ! GFf%G(:,:,iw) = matmul(matmul(w1, SigmaG), w2) !GFf%L(:,:,iw) +  GFf%R(:,:,iw) -  GFf%A(:,:,iw)

  GFf%G(:,:,iw) = GFf%L(:,:,iw) + GFf%R(:,:,iw) - GFf%A(:,:,iw)
  
  write(3,*) '==================with Identity================'
  write(3,*) GF0%L(1,1,iw) + GF0%R(1,1,iw) - GF0%A(1,1,iw), GF0%L(2,2,iw) + GF0%R(2,2,iw) - GF0%A(2,2,iw), &
       GF0%L(3,3,iw) + GF0%R(3,3,iw) - GF0%A(3,3,iw)
  
  
  write(3,*) GFf%L(1,1,iw) + GFf%R(1,1,iw) - GFf%A(1,1,iw), GFf%L(2,2,iw) + GFf%R(2,2,iw) - GFf%A(2,2,iw), &
       GFf%L(3,3,iw) + GFf%R(3,3,iw) - GFf%A(3,3,iw) 
  
  write(3,*) '-----------GF Retarded ------------'
  do i = 1, Natoms
     write(3,*) (GFf%R(i,j,1), j= 1, Natoms)
  end do
  write(3,*) '-----------------------------------------'
  write(3,*) '-----------GF Lesser ------------'
  do i = 1, Natoms
     write(3,*) (GFf%R(i,j,1), j= 1, Natoms)
  end do
  write(3,*) '-----------------------------------------'

STOP
end subroutine G_full

 
!=====================================================
!======== Non-interacting GFs ========================
!=====================================================  
  
!......Calculates lesser Greens function at time = 0
subroutine GL_of_0()
  !.....Lesser Greens function at time = 0
  !.....ei - eigenvalues of the Hamiltonian 
  implicit none
  integer :: i, j, k1
  complex*16 :: s
  real*8 :: pp
  
  pp=delta/(2.d0*pi)
  do i = 1, Natoms 
     s =(0.d0, 0.d0)
     do k1 = 1, N_of_w
        s = s + GF0%L(i,i,k1)
     end do
     G_nil(i,i)=s*pp
  end do
end subroutine GL_of_0

subroutine G0_R_A()
  !............non-interacting Greens functions: GR and Ga,  Eq. (5) and Eq. (6) in 'Current_Hubbard_Equations' document (CHE)
  implicit none
  integer :: j, i
  real*8 :: w
  
    do j = 1, N_of_w
       work1 = -H + 0.5d0*(im/hbar) * (GammaL + GammaR) !LK <========= must be +
       w = omega(j)
       do i = 1 , Natoms
          work1(i,i) = work1(i,i) + hbar*(w +im*0.01)
        end do
       
       call Inverse_complex(Natoms, work1, info)
       call Hermitian_Conjg(work1, Natoms, work2)
       
       GF0%r(:,:,j) = work1
       GF0%a(:,:,j) = work2
    end do

end subroutine G0_R_A

subroutine G0_L_G(Volt)
  !............non-interacting Greens functions: G> and G< for all omega on the grid, Eq. (16) and (17) in CHE
  implicit none
  real*8 :: Volt, w
  integer :: j, i
  
    work1 = (0.d0, 0.d0) ; work2 =(0.d0, 0.d0) ; work3 = (0.d0, 0.d0);  work4 = (0.d0, 0.d0)
    do j = 1 , N_of_w
       w = omega(j)
       work1 = GF0%r(:,:,j) 
       work2 = GF0%a(:,:,j) 
       work3 = matmul(matmul(work1, im*(fermi_dist(w, Volt)*GammaL + fermi_dist(w, 0.d0)*GammaR)), work2) 
       !  work4 = matmul(matmul(work1, im*((fermi_dist(w, Volt)-1.d0)*GammaL + (fermi_dist(w, 0.d0)-1.d0)*GammaR)), work2)
       GF0%L(:,:,j) = work3
       !   GF0%G(:,:,j) = work4
       GF0%G(:,:,j) =  GF0%L(:,:,j) + GF0%R(:,:,j) - GF0%A(:,:,j)
    end do

end subroutine G0_L_G

!=====================================================
!========Calcualtions needed for full GFs=============
!===================================================== 

subroutine first_order_sigma(Sigma1)
  implicit none
  integer :: i, s, s1 
  complex*16 :: Hartree
  complex*16, dimension(:,:) :: Sigma1(Natoms, Natoms) 
  
  Sigma1 = (0.d0, 0.d0); Hartree = (0.d0, 0.d0)
  
  do i = 1, Natoms, 2 !..orbitals 

     Hartree = G_nil(i,i)+G_nil(i+1,i+1)
     !.. spins 
     do s = 0, 1
        do s1 = 0, 1
           if (s .eq. s1) Sigma1(i+s, i+s1) = -im*hbar*Hub(i)*Hartree 
           Sigma1(i+s, i+s1) = Sigma1(i+s, i+s1) + im*hbar*Hub(i)*G_nil(i+s, i+s1)
        end do
     end do
 
  end do
  
!  write(3,*) '-----------Hartree-Fock------------'
!  write(3,*) Sigma1(1,1), Sigma1(2,2), Sigma1(3,3)
!  write(3,*) '-----------------------------------'
  
end subroutine first_order_sigma

!......................Calculation of Omega terms for the self-energies, Eq. (9) in CHE
complex*16 function Omega_r(i, j, sp, sp1, iw)
  implicit none
  integer :: i, j, iw, k_1, k_2, k_3, m, n, s, s1, sp, sp1, ii, jj
  complex*16 :: Omr
  real*8 :: pp
  
  ii = i +sp; jj= j + sp1
  
  Omr = (0.d0, 0.d0)
  do k_1 = 1, N_of_w
     do k_2 = 1, N_of_w
        k_3 = iw- k_1 +k_2

        if (k_3 .ge. 1 .and. k_3 .le. N_of_w) then       
           do s = 0, 1 !...Sum over orbitals
              do s1 = 0, 1
                 m = i+s
                 n = j+s1
                
                    !.....both second order diagram contributions 
                 Omr = Omr - GF0%r(ii,m,k_1)*GF0%L(m,n,k_2)*GF0%L(n,jj,k_3) &
                      - GF0%L(ii,m,k_1)*GF0%a(m,n,k_2)*GF0%L(n,jj,k_3) &
                      - GF0%L(ii,m,k_1)*GF0%G(m,n,k_2)*GF0%a(n,jj,k_3) & !..Eq. (24) in CHE
                      
                      + GF0%r(m,n,k_1)*GF0%L(n,m,k_2)*GF0%L(ii,jj,k_3) & 
                      + GF0%L(m,n,k_1)*GF0%a(n,m,k_2)*GF0%L(ii,jj,k_3) & !... Eq. (21) in CHe
                      + GF0%L(m,n,k_1)*GF0%G(n,m,k_2)*GF0%a(ii,jj,k_3)
              end do
           end do
        end if
        
     end do
  end do
  
  pp = (delta/2.d0*pi)
  Omega_R = Omr*pp*pp
end function Omega_r

complex*16 function int_SigL(i,j,sp,sp1,iw) !... interaction contributions of Eq. (23) + Eq. (26) 
  implicit none
  integer :: i, j, k1, k2, k3, iw, m, n, sp, sp1, s, s1, ii, jj
  complex*16 ::  SigL
  real*8 :: pp

  ii = i +sp; jj= j + sp1

  SigL=(0.0d0, 0.0d0)
  do k1 = 1, N_of_w
     do k2 = 1, N_of_w
        k3 = iw- k1 +k2
        
        if (k3 .ge. 1 .and. k3 .le. N_of_w) then
           do s = 0, 1 !...Sum over spin
              do s1 = 0, 1
                 m = i+s
                 n = j+s1
                 
                 SigL = SigL + GF0%L(m,n,k1)*GF0%G(n,m,k2)*GF0%L(ii,jj,k3) &
                      - GF0%L(ii,m,k1)*GF0%G(m,n,k2)*GF0%L(n,jj,k3) 
                 
              end do
           end do
        end if
        
     end do
  end do
  pp = (delta/2.d0*pi)
  int_SigL = SigL*pp
end function int_SigL

complex*16 function int_SigG(i,j,sp,sp1,iw) !... interaction contributions of Eq. (23) + Eq. (26) 
  implicit none
  integer :: i, j, k1, k2, k3, iw, m, n, sp, sp1, s, s1, ii, jj
  complex*16 ::  SigG
  real*8 :: pp

  ii = i +sp; jj= j + sp1

  SigG=(0.0d0, 0.0d0)
  do k1 = 1, N_of_w
     do k2 = 1, N_of_w
        k3 = iw- k1 +k2
        
        if (k3 .ge. 1 .and. k3 .le. N_of_w) then
           do s = 0, 1 !...Sum over spin
              do s1 = 0, 1
                 m = i+s
                 n = j+s1
                 
                 SigG = SigG + GF0%G(m,n,k1)*GF0%L(n,m,k2)*GF0%G(ii,jj,k3) &
                      - GF0%G(ii,m,k1)*GF0%L(m,n,k2)*GF0%G(n,jj,k3) 
                 
              end do
           end do
        end if
     end do
  end do
  pp = (delta/2.d0*pi)

  int_SigG = SigG*pp
end function int_SigG


subroutine print_sf(iteration)
  integer :: iteration,n,iw,j
  character :: fn*3 
  
  n=Natoms ; if(n.gt.10) n=10
  write(fn,'(i0)') iteration
  open(7,file='sf_'//trim(fn)//'.dat',status='unknown')
  do iw=1,N_of_w
     write(7,'(f10.5,x,10(e12.5,x))') omega(iw),(- 2.d0*hbar*AIMAG(GF0%R(j,j,iw)),j=1,n)
  end do
  close(7)
  write(*,*) '... written sf_'//trim(fn)//'.dat'
end subroutine print_sf

subroutine print_3matrix(iteration,X,name)
  integer :: iteration,n,iw,j,i
  character :: fn*3,name*3
  complex*16 :: X(Natoms,Natoms,N_of_w)
   
  n=Natoms ; if(n.gt.10) n=10
  write(fn,'(i0)') iteration
  open(7,file=name//'_'//trim(fn)//'.dat',status='unknown')
  do iw=1,N_of_w
     write(7,'(/i3,x,f10.5)') iw,omega(iw)
     do i=1,n
        write(7,*) i,j,(X(i,j,iw),j=1,n)
     end do
  end do
  close(7)
  write(*,*) '... written '//name//'_'//trim(fn)//'.dat'
end subroutine print_3matrix

end module GreensFunctions
