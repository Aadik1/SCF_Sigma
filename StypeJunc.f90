program StypeJunction_Spin
  use DefineHamiltonian
  use GreensFunctions
  use OMP_Lib
  implicit none
  real*8 :: V1, J_up, J_down, total_time
  integer ::  k, i
  character(len=30) :: vfn
  logical :: first
  
  !.........................Deifnes Hamiltonian
  call input_SOC()
  
  allocate(H(Natoms,Natoms))
  allocate(GammaL(Natoms, Natoms)); allocate(GammaR(Natoms, Natoms))
  call SOC_Hamiltonian()
   
  GammaL = (0.d0, 0.d0); GammaL(1,1) = GammaL_up;  GammaL(2,2) = GammaL_dw
  
  GammaR = (0.d0, 0.d0); GammaR(Natoms-1, Natoms-1) = GammaR_up;  GammaR(Natoms, Natoms) = GammaR_dw
  
  !......................Defines the level width funcitons for L,R-leads to central region, i.e. the respective couplings
  
  allocate(Hub(Natoms)) !...Hubbard only on orbital indices, not spin indices
 
  Hub = 0.d0 
  do i = 1, Natoms, 2
     Hub(i) = Hubbard !>...read(11, *) Hub(i)
  end do
  close(11)
  
  !.......................Finds the eigenvalues of the Central Hamiltonian and uses it to define the time-independent G(0)
  
  allocate(Eigenvec(Natoms, Natoms))  
  allocate(Ev(Natoms))
  
  Eigenvec = H
  call complex_eigen_symm_martix(Eigenvec, Natoms, Ev)
  
  w_init = Ev(1)-dw; w_fin = Ev(Natoms)+up
  N_of_w = (w_fin - w_init)/delta ; write(*,*) 'N_of_w:', N_of_w, 'w_fin:', w_fin, 'w_init:', w_init
  
  write(3,*) 'Temp:', T, 'delta:', delta, 'Hamiltonian Dimension:', Natoms, '# of omegas:', N_of_w
  
  allocate(omega(N_of_w))
  
  do i = 1, N_of_w 
     omega(i) = w_init + i*delta
  end do

  call PrintFunctions()
  deallocate(Eigenvec, Ev)
  
!.......................Allocates Greens Functions for the integrals within the Sigmas across all omegas and organises them into user defined type GF0 and GFf in 'mod_Gfunctions.f90'
  allocate(GF0%r(Natoms, Natoms, N_of_w));  allocate(GF0%a(Natoms, Natoms, N_of_w))
  allocate(GF0%L(Natoms, Natoms, N_of_w));  allocate(GF0%G(Natoms, Natoms, N_of_w))

  allocate(GFf%L(Natoms, Natoms, N_of_w));  allocate(GFf%G(Natoms, Natoms, N_of_w))
  allocate(GFf%R(Natoms, Natoms, N_of_w));  allocate(GFf%A(Natoms, Natoms, N_of_w))

  allocate(G_nil(Natoms, Natoms))
!.......................Allocates all the full self-energies and full Greens functions needed in the current
!  allocate(SigmaL(Natoms, Natoms)); allocate(Sigma1(Natoms, Natoms))
!  allocate(SigmaR(Natoms, Natoms))
  
!...............calculate GR and GA for all voltages on the omega grid
  allocate(work1(Natoms, Natoms)); allocate(work2(Natoms, Natoms)); allocate(work3(Natoms, Natoms)); allocate(work4(Natoms, Natoms))
  
  !.......................Calculates and plots Voltage vs Current curve
  !...For now, sticking to a single voltage to run pulay

  open(3, file='Print.dat', status='unknown')
  
  write(vfn,'(i0)') order
  open(30, file='Volt_Current_'//trim(vfn)//'.dat', status='unknown')
  
  print *, 'Pre-voltage' 
  first=.true.
  do k = 0, Volt_range
     V1 = V + k*0.05
     
     call SCF_GFs(V1, first)
     
     if(first) first=.false.
     
     call Current(V1, J_up, J_down)
     write(30, *) V1, J_up, J_down
     print *, 'Progress:', k/(Volt_range*0.01), '%', J_up, J_down
  end do
  
  close(3)
  close(30)

  deallocate(GFf%L, GFf%G,GFf%R, GFf%A)
  deallocate(GF0%r, GF0%a, GF0%L, GF0%G) 
  deallocate(work1, work2, work3, work4)
  deallocate(H, Hub, omega)
  !deallocate(SigmaL, Sigma1, SigmaR);
  deallocate(GammaL, GammaR, G_nil)
end program StypeJunction_Spin


