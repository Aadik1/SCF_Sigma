subroutine input_SOC()
  use GreensFunctions 
  implicit none

  open(22, file='inputSOC.dat', status='old')
  read(22, *) T, V, mu, Volt_range
  read(22,*) order
  read(22,*) dw,up,delta
  read(22,*) pulay
  read(22,*) hel_radius, hel_length, N_ions, N_turns, hand
  read(22,*) E_CC_up, E_CC_down, t_hop, lamb, Hubbard
  read(22,*) GammaL_up, GammaL_dw, GammaR_up, GammaR_dw
  close(22) 

  beta = 1.d0/(kb*T)
  Natoms = 2*N_ions*N_turns !...# of ions per turn of the helix  and the # of turns multiplied by 2
  
  write(*,*) 'T:', T, 'V:', V, 'mu:', mu, 'Volt_range:', Volt_range, 'Hubbard:', Hubbard
  write(*,*) 'Order:', order, 'Natoms:', Natoms
  write(*,*) 'delta:', delta
  write(*,*) 'pulay:', pulay
end subroutine input_SOC

subroutine PrintFunctions()
  use GreensFunctions
  implicit none
  integer :: i, j
  complex*16 :: diff

  open(12, file='info_Hamiltonian.dat', status='replace')
  
  write(12,'(/a)') '... Hamiltonian:'
  do i = 1, Natoms
     write(12, '(10(a,2f10.5,a))') (' [',H(i,j),'] ', j=1,Natoms)
  end do
  
  !.... checking if Hermitian  
  
  write(12,'(/a)') '... Checking the Hamiltonian is Hermitian:'
  do i=1,Natoms
     do j=i,Natoms
        diff=H(i,j)-conjg(H(j,i))
        write(12,*) i,j,diff
     end do
  end do
  
  write(12,'(/a)') '... Eigenvalues:'
  do i = 1, Natoms
     write(12, *) i,'.', Ev(i)
  end do
  
  write(12,'(/a)') '... Eigenvectors:'
  do j = 1, Natoms
     write(12, '(i3,10(a,2f10.5,a))') j,(' [',Eigenvec(i,j),'] ', i = 1, Natoms)
  end do

  close(12)
end subroutine PrintFunctions

subroutine trans(iw, Volt, trans_up, trans_down) !....square bracket terms of Eq. (2) in CHE
  use GreensFunctions
  implicit none
  integer :: iw,i,ii
  complex*16 :: trace1, trace2
  real*8 :: Volt, w, trans_up, trans_down

  w = omega(iw)
  
  work1 = GFf%L(:,:,iw)
  work2 = GFf%G(:,:,iw)
  
  work3 = im*matmul(GammaL, (fermi_dist(w, Volt)-1.d0)*work1 - fermi_dist(w, Volt)*work2)

  trace1 = (0.d0,0.d0);   trace2 = (0.d0,0.d0)
  do i = 1, Natoms, 2
     ii = i + 1 

     trace1 = trace1 + work3(i,i)
     trace2 = trace2 + work3(ii,ii)
     
  end do

  trans_up = real(trace1)/(2.d0*pi)
  trans_down = real(trace2)/(2.d0*pi)
end subroutine trans

subroutine Current(Volt, J_up, J_down)
  use GreensFunctions
  implicit none
  real*8 :: Volt, J_up, J_down, trans_up, trans_down
  integer :: iw
  
  J_up = 0.d0; J_down = 0.d0
  do iw = 1, N_of_w
     call trans(iw, Volt, trans_up, trans_down)
     J_up = J_up + trans_up
     J_down = J_down + trans_down
  end do
  
  J_up = J_up*(delta/hbar)
  J_down = J_down*(delta/hbar)
  
end subroutine Current
