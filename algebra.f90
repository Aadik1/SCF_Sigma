  subroutine BvecA_to_C(A,B,n) !..A is a square matrix, B is a vector, kills original A
    implicit none 
    integer :: n,i,k 
    complex*16, dimension(:,:) :: A(n,n),B(n),C(n,n)
    complex*16 :: s
    
    do i=1,n
       do k=1,n
          s=(0.d0,0.d0)
          s=s+B(k)*A(i,k)
          C(i,k)=s
       end do
    end do
    A = C
  end subroutine BvecA_to_C
  
  subroutine Hermitian_conjg(A,N,B)
    implicit none 
    complex*16 :: A(N,N),B(N,N)
    integer :: i,j,N
    
    do i=1,N
       do j=1,N
          B(i,j)=conjg(A(j,i))
       end do
    end do
  end subroutine Hermitian_Conjg
  
  subroutine Trace_of_A(A,n,trace) !... complex trace
    implicit none 
    integer :: n,i
    complex*16, dimension(:,:) :: A(n,n)
    complex*16 :: s, trace
    
    s = 0.d0
    do i=1,n
       s = s + A(i,i)
    end do
    trace=s 
  end subroutine Trace_of_A

  subroutine vec_cross_product(a,b,c, N)
    implicit none
    integer :: N
    real*8, dimension(N) :: a, b, c
    
    c = 0.d0
    
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - b(2)*a(1)
  end subroutine vec_cross_product

 real*8  function mod_vec(r)
    implicit none
    real*8, dimension(3) :: r
    real*8 :: modr

    modr = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
    mod_vec = modr
  end function mod_vec
  
  subroutine AB_to_C(A,B,C,n)
    implicit none 
    integer :: n,i,j,k
    complex*16, dimension(:,:) :: A(n,n),B(n,n),C(n,n)
    complex*16 :: s
    
    do i=1,n
       do j=1,n
          s=(0.d0,0.d0)
          do k=1,n
             s=s+A(i,k)*B(k,j)
          end do
          C(i,j)=s
       end do
    end do
  end subroutine AB_to_C
  
  subroutine Trace_of_AB(A,B,n,trace)
    implicit none 
    integer :: n,i,j
    complex*16, dimension(:,:) :: A(n,n), B(n,n)
    complex*16 :: s
    real*8 :: trace
    
    s = 0.d0
    do i=1,n
       do j=1,n
          s = s + A(i,j)*B(j,i)
       end do
    end do
    trace=s
  end subroutine Trace_of_AB
  
  subroutine Filelen(N, filename)
    implicit none
    integer :: N, iostat
    character(len=300) :: line
    character(*), intent(in) :: filename
    
    open(unit=11, file=filename, status='old', iostat=iostat)
    if (iostat /= 0) then
       write(*,*) 'Error reading file', iostat
       stop
    end if
    N = 0
    do
       read(11, '(A)', iostat=iostat) line
       if (iostat /= 0) exit
       N = N + 1
    end do
    close(11)
    
  end subroutine Filelen
  
  !....................defines the rotation vector about the z-axis
  subroutine VecAssign(Ri, wi)
    implicit none
    real*8, dimension(:) :: Ri(3), wi(3)
    
    wi = 0.d0
    wi(1) = Ri(2)
    wi(2) = -Ri(1)
  end subroutine VecAssign
