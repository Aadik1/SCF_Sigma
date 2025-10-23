module GreensFunctions
  use DefineHamiltonian
  use OMP_LIB
  implicit none
  integer :: INFO

  real*8, parameter :: epsilon = 1e-9
  real*8, allocatable, dimension(:) :: Hub, omega, Ev
  real*8 :: pulay, dw, up
  
  complex*16, allocatable, dimension(:,:) :: GammaL, GammaR, Eigenvec, G_nil
  complex*16, allocatable, dimension(:,:,:) ::  SigmaL, SigmaR, SigmaG
  complex*16, allocatable, dimension(:,:) :: work1, work2, work3, work4
  
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
    real*8 :: arg, w, V

    arg = (w - mu + V/hbar)*beta
    fermi_dist =  1.d0/(exp(hbar*arg) +1.d0)
    
  end function fermi_dist

!====================================================
!========== Self consistency field calculations =====
!====================================================

  subroutine SCF_calc(Volt)
    implicit none
    integer :: i, iw, iteration
    real*8 :: Volt, w, err, diff
    complex*16, allocatable, dimension(:,:,:) ::  SigmaRi, SigmaLi
    complex*16, allocatable, dimension(:,:) :: Sigma1
    
    
    allocate(Sigma1(Natoms, Natoms))
    allocate(SigmaRi(Natoms, Natoms, N_of_W)); allocate(SigmaLi(Natoms, Natoms, N_of_w))
    
    !............full SigmaR due to interactions Eq. (7)
    SigmaR = (0.d0, 0.d0); SigmaL= (0.d0, 0.d0); SigmaG = (0.d0, 0.d0)
    
    iteration = 0
    
    write(*,*) '........SCF Calculations at Voltage:', Volt, '..........'
    
    call G0_R_A()
    call G0_L_G(Volt)
    
    DO
       iteration = iteration + 1
       write(*,*) '.... ITERATION = ',iteration,' ....'
       
       call GL_of_0()
       
       call all_sigmas(SigmaRi, SigmaLi, Sigma1)

        !...Simple Pulay mixing scheme for Sigmas     
       SigmaR = pulay*SigmaRi + (1.d0 - pulay)*SigmaR
       SigmaL = pulay*SigmaLi + (1.d0 - pulay)*SigmaL
       
       
       !$OMP PARALLEL DO &
       !$OMP& PRIVATE(iw, INFO)
       do iw = 1, N_of_w 
          call G_full(iw, Volt)
       end do
       !$OMP END PARALLEL DO
       
       err=0.0d0
       do iw = 1, N_of_w
          do i=1,Natoms
             diff=(SigmaR(i,i,iw))-(SigmaRi(i,i,iw))
             err=err +diff*diff
          end do
       end do
       write(*,*) 'err = ',sqrt(err)

       if (err .le. epsilon .or. order .eq. 0) then
       write(*,*)'... REACHED REQUIRED ACCURACY ...'
       exit
    end if
    
 end DO
 
 deallocate(SigmaLi, SigmaRi, Sigma1)
end subroutine SCF_Calc

subroutine all_sigmas(SigmaRi, SigmaLi, Sigma1)
  implicit none
  integer :: i, j, iw, sp, sp1, ii, jj
  complex*16 :: SigL, SigG, Omr
  complex*16 ::  SigmaRi(:,:,:), SigmaLi(:,:,:)
  complex*16 :: Sigma1(:,:)

  SigmaRi = (0.d0, 0.d0); SigmaLi = (0.d0, 0.d0)
  
  !$OMP PARALLEL DO &
  !$OMP& PRIVATE(i,j,sp,sp1,ii,jj,Omr,SigL,SigG,Sigma1) 
  
     !.....Calculate Sigmas      
       do iw = 1, N_of_w
          
          if (order .eq. 1) then
             call first_order_sigma(Sigma1)
             SigmaRi(:,:,iw) = Sigma1
          else if (order .eq. 2) then
             call first_order_sigma(Sigma1)
             
             do i = 1, Natoms, 2
                do sp = 0, 1
                   ii = i +sp
                   
                   do j = 1, Natoms, 2
                      do sp1 = 0, 1
                         jj= j + sp1
                         
                         Omr = Omega_r(i,j,sp,sp1,iw)
                         SigmaRi(ii,jj,iw) =  Sigma1(ii,jj) + (Hub(j)*Hub(i)*OmR)*hbar**2
                         
                      end do
                   end do
                   
                end do
             end do
             
             !..............full SigmaL, Eq. (3) and (4)     
             !.....Interaction contribution of both Sigmas     
             
             do i = 1, Natoms, 2
                do sp = 0, 1
                   ii = i +sp
                   
                   do j = 1, Natoms, 2
                      do sp1 = 0, 1
                         jj= j + sp1                 
                         
                         call int_SigLnG(i,j, sp, sp1, iw, SigL, SigG)
                         SigmaLi(ii,jj,iw) = Hub(j)*Hub(i)*SigL*hbar**2
                         
                      end do
                   end do
                   
                end do
             end do
            
       end if
    end do
    !$OMP END PARALLEL DO   
  end subroutine all_sigmas

!=====================================================
!================== Full GFs =========================
!=====================================================
  
!.............need G0f%L to calculate GL_of_0
!.............Calculates GFs at every omega and Voltage simultaneously
subroutine G_full(iw, Volt) !... Full Greens function, leaves Retarded and Advanced in the work arrays, application of Eq. (16) and (17), but with the full Sigmas, Eq. (3), (7) and (8) in CHE
  implicit none
  integer :: i, j, iw, sp, sp1, ii, jj, N
  real*8 :: Volt, w 
  complex*16, allocatable, dimension(:,:) :: work_1, work_2
  !............full Gr and Ga, Eq. (5) and (6)

  allocate(work_1(Natoms, Natoms), work_2(Natoms, Natoms))

  w = omega(iw)
  work_1 = -H + 0.5d0*(im/hbar)*(GammaL + GammaR) - SigmaR(:,:,iw) ! LK <== must be + for emb and minus for interaction sigma
  do i = 1 , Natoms
     work_1(i,i) = work_1(i,i) + hbar*w !(w+im*0.01)
  end do
  
  call Inverse_complex(Natoms, work_1, info)
  call Hermitian_Conjg(work_1, Natoms, work_2)
  
  GF0%r(:,:,iw) = work_1; GF0%a(:,:,iw) = work_2
  !.....Embedding contribution of both Sigmas
  
  !.............full GL and GG, Eq. (16) and (17)
  GF0%l(:,:,iw) = matmul(matmul(GF0%r(:,:,iw), (im*(fermi_dist(w, Volt)*GammaL + fermi_dist(w, 0.d0)*GammaR)/hbar) +SigmaL(:,:,iw)), GF0%a(:,:,iw)) !.. GL = Gr * SigmaL * Ga     
  GF0%g(:,:,iw) = GF0%l(:,:,iw) + GF0%r(:,:,iw) - GF0%a(:,:,iw)
  
  deallocate(work_1, work_2)   
end subroutine G_full

!=====================================================
!======== Non-interacting GFs ========================
!=====================================================  
  
!......Calculates lesser Greens function at time = 0
subroutine GL_of_0()
  !.....Lesser Greens function at time = 0
  !.....ei - eigenvalues of the Hamiltonian 
  implicit none
  integer :: i, j, k1, s1, s2
  complex*16 :: s
  real*8 :: pp
  
  pp=delta/(2.d0*pi)
  G_nil = (0.d0, 0.d0)
  do i = 1, Natoms,2
     do s1=0,1
        do s2=0,1
           
           s =(0.d0, 0.d0)
           do k1 = 1, N_of_w
              s = s + GF0%L(i+s1,i+s2,k1)
           end do
           G_nil(i+s1,i+s2)=s*pp
           
        end do
     end do
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
          work1(i,i) = work1(i,i) + hbar*w!(w +im*0.01)
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
    
    work1 = (0.d0, 0.d0) ; work2 =(0.d0, 0.d0) ; work3 = (0.d0, 0.d0); work4 = (0.d0, 0.d0)
    do j = 1 , N_of_w
       w = omega(j)
       work1 = GF0%r(:,:,j) 
       work2 = GF0%a(:,:,j)

       work3 = matmul(matmul(work1, im*(fermi_dist(w, Volt)*GammaL + fermi_dist(w, 0.d0)*GammaR)), work2) 
      ! work4 = matmul(matmul(work1, im*((fermi_dist(w, Volt)-1.d0)*GammaL + (fermi_dist(w, 0.d0)-1.d0)*GammaR)), work2)
       GF0%L(:,:,j) = work3
       !GF0%G(:,:,j) = work4
       GF0%G(:,:,j) =  GF0%L(:,:,j) + GF0%R(:,:,j) - GF0%A(:,:,j)
    end do
  end subroutine G0_L_G 

!=====================================================
!========Calcualtions needed for full GFs=========== ==
!=====================================================

  subroutine first_order_sigma(Sigma1)
    implicit none
    integer :: i, s, s1, N 
    complex*16 :: Hartree
    complex*16 :: Sigma1(:,:)
    
    Sigma1 = (0.d0, 0.d0); Hartree = (0.d0, 0.d0)
    
    do i = 1, Natoms, 2 !..orbitals 
       
       Hartree = G_nil(i,i)+G_nil(i+1,i+1)
       !.. spins 
       do s = 0, 1
          do s1 = 0, 1
             if (s .eq. s1) then 
                Sigma1(i+s, i+s1) = im*hbar*Hub(i)*(G_nil(i+s, i+s1) - Hartree)
             else
                Sigma1(i+s, i+s1) = im*hbar*Hub(i)*G_nil(i+s, i+s1)
             end if
          end do
       end do
       
    end do
    
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
                 m = i+s1; n= j+s
                 !.....both second order diagram contributions
                 !_________ 3rd diagram
                 Omr = Omr + GF0%r(ii,jj,k_1)*GF0%L(n,m,k_2)*GF0%G(m,n,k_3) & 
                      + GF0%L(ii,jj,k_1)*GF0%L(n,m,k_2)*GF0%r(m,n,k_3) & 
                      + GF0%L(ii,jj,k_1)*GF0%a(n,m,k_2)*GF0%L(m,n,k_3) &
                      
                      - GF0%r(ii,n,k_1)*GF0%L(n,m,k_2)*GF0%G(m,jj,k_3) &
                      - GF0%L(ii,n,k_1)*GF0%L(n,m,k_2)*GF0%r(m,jj,k_3) &
                      - GF0%L(ii,n,k_1)*GF0%a(n,m,k_2)*GF0%L(m,jj,k_3) 
                      !_________ 4th diagram
                   
              end do
           end do
        end if
        
     end do
  end do
  
  pp = delta/(2.d0*pi)
  Omega_R = Omr*pp*pp
end function Omega_r

subroutine int_SigLnG(i,j,sp,sp1,iw,SigL,SigG) !... interaction contributions of Eq. (23) + Eq. (26) 
  implicit none
  integer :: i, j, k1, k2, k3, iw, m, n, sp, sp1, s, s1, ii, jj
  complex*16 ::  SigL, SigG
  real*8 :: pp

  ii = i +sp; jj= j + sp1

  SigL=(0.0d0, 0.0d0) ! SigG=(0.0d0, 0.0d0)
  
  do k1 = 1, N_of_w
     do k2 = 1, N_of_w
        k3 = iw- k1 +k2

        if (k3 .ge. 1 .and. k3 .le. N_of_w) then   
           do s = 0, 1 !...Sum over spin
              do s1 = 0, 1
                 m = i+s1; n= j+s
                 
                 SigL = SigL + GF0%L(ii,jj,k1)*GF0%G(n,m,k2)*GF0%L(m,n,k3) &
                      - GF0%L(ii,n,k1)*GF0%G(n,m,k2)*GF0%L(m,jj,k3) 
                      
                
                 SigG = SigG + GF0%G(m,n,k1)*GF0%L(n,m,k2)*GF0%G(ii,jj,k3) &
                      - GF0%G(ii,m,k1)*GF0%L(m,n,k2)*GF0%G(n,jj,k3) 
                 
              end do
           end do
        end if
        
     end do
  end do
  pp = delta/(2.d0*pi)
  SigL = SigL*pp*pp ; SigG = SigG*pp*pp
end subroutine int_SigLnG

end module GreensFunctions
