!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! U S I N G   L A P A C K !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gen_eigen_symm_matrix(A,S,N,E)
!...............................................
!  Generalised eigenporoblem: A*X=E*S*X 
!...............................................
! INPUT
! A,S - symmetric matrices;
!  INFO = 0 - succeeded  
!  INFO =/= 0 - not succeeded
!...............................................
! OUTPUT
!  (i) eigenvectors of A are in columns of A
!      (the original A is destroyed)
!  (ii) eigenvalues of A are in E
!...............................................
  
  implicit none
  real*8 :: A(N,N),E(N),S(N,N)
  integer :: N,INFO,LWORK
  real*8,dimension(:),allocatable ::  WORK

!_______ determine first the dimension of WORK and then compute


  allocate(WORK(1))
  call dsygv(1,'V','U',N,A,N,S,N,E,WORK,-1,INFO)
  LWORK=WORK(1)
!  write(*,*) LWORK
  deallocate(WORK)
  allocate(WORK(LWORK))
  call dsygv(1,'V','U',N,A,N,S,N,E,WORK,LWORK,INFO)
  deallocate(WORK)
  if(INFO.eq.0) then
     write(9,*)'... eigenproblem successful'
  else
     write(*,*)'... eigenproblem UNsuccessful' ; stop
  end if
end subroutine gen_eigen_symm_matrix

subroutine complex_eigen_symm_martix(A,N,E)
!...............................................
!  Complex eigenporoblem: A*X=E*X 
!...............................................
! INPUT
! A,S - Hermitian matrices;
!  INFO = 0 - succeeded  
!  INFO =/= 0 - not succeeded
!...............................................
! OUTPUT
!  (i) eigenvectors of A are in columns of A
!      (the original A is destroyed)
!  (ii) eigenvalues of A are in E
!...............................................
  
  implicit none
  integer :: N,INFO,LWORK
  complex*16,dimension(:),allocatable ::  WORK
  real*8, dimension(:),allocatable :: RWORK
  complex*16 :: A(N,N)
  real*8 :: E(N)
  
!_______ determine first the dimension of WORK and then compute

  allocate(WORK(1)); allocate(RWORK(3*N-2))  
  call ZHEEV('V','U',N,A,N,E,WORK,-1,RWORK,INFO)
  LWORK=int(WORK(1))
!  write(*,*) LWORK
  deallocate(WORK)
  allocate(WORK(LWORK))
  call ZHEEV('V','U',N,A,N,E,WORK,LWORK,RWORK,INFO)
  deallocate(WORK)
  deallocate(RWORK)
  if(INFO.eq.0) then
     write(9,*)'... eigenproblem successful'
  else
     write(*,*)'... eigenproblem UNsuccessful' ; stop
  end if
end subroutine complex_eigen_symm_martix


subroutine eigen_symm_matrix(A,N,E)
!...............................................
! INPUT
! A - symmetric matrix; 
!  INFO = 0 - succeeded  
!  INFO =/= 0 - not succeeded
!...............................................
! OUTPUT
!  (i) eigenvectors of A are in columns of A
!      (the original A is destroyed)
!  (ii) eigenvalues of A are in E
!...............................................
  implicit none
  real*8 :: A(N,N),E(N)
  integer :: N,N6, info
  real*8,dimension(:),allocatable ::  WORK
  logical :: verb
  
!_______ determine first the dimension of WORK and then compute

  allocate(WORK(1))
  call dsyev('V','U',N,A,N,E,WORK,-1,INFO)
  N6=WORK(1) 
  deallocate(WORK) ; allocate(WORK(N6))
  call dsyev('V','U',N,A,N,E,WORK,N6,INFO)
  deallocate(WORK)
  if(INFO.eq.0) then
     if(verb) write(9,*)'... eigenproblem successful'
  else
     write(*,*)'... eigenproblem UNsuccessful' ; stop
  end if
end subroutine eigen_symm_matrix

subroutine inverse_root(N,A,INFO,case)
!...............................................................
! INPUT
! A - symmetirc matrix; hence its A^(-1/2) or A^(1/2) is also
!     symmetric
!  INFO = 0 - succeeded  
!  INFO = 1 - not succeeded as diagonalisation routine failed
!  INFO =-1 - not succeeded as A has negative eigenvalues
!  INFO = 2 - not succedded as A is not symmetric
!...............................................................
! case =  1 - calculates root of A ( = sqrt(A) )
! case = -1 - calculates the inverse root of A ( = 1/sqrt(A) )
!...............................................................
! OUTPUT
! A = 1/sqrt(A) [the original A is destroyed]  
!...............................................................

  implicit none
  real*8,dimension(:,:) ::  A(N,N)
  real*8,dimension(:,:),allocatable :: B 
  real*8,dimension(:) ::  E(N)
  real*8,dimension(:),allocatable :: WORK
  integer :: INFO,i,j,N,k,LWORK,case
  real*8 :: s,small=1.0e-6
  logical :: verb

      INFO=0
!......... check if A is symmetric

      do i=1,N
         do j=i,N
            if(A(i,j).ne.A(j,i)) then
               INFO=2; return
            end if
         end do
      end do
      
!......... run eigenproblem routine

!_______ determine first the dimension of WORK

      allocate(WORK(1))
      call dsyev('V','U',N,A,N,E,WORK,-1,INFO)
      LWORK=WORK(1)
      deallocate(WORK) 
      if(verb) write(9,*) 'dim = ',LWORK
      
      allocate(WORK(LWORK))
      call dsyev('V','U',N,A,N,E,WORK,LWORK,INFO)
      write(9,*) 'INFO=',INFO
      deallocate(WORK)
      if(INFO.ne.0) then
         write(9,*)'... DSYEV: inverse root NOT successful'
         INFO=1 ; return
      end if
      
!......... check if A is positive definite

!      write(*,*)"Eigenvalues:"
!      write(*,'(a,i5,a,e12.6)') ('E(',k,') = ',E(k),k=1,N)
      do k=1,N
         if(E(k).le.small) then
            write(*,'(a,i5,a,e12.6)')'ERROR: negative (small) ',k,'-th eigenvalue: ',E(k)
            INFO=-1 ; return
         end if
      end do

!........ construct the inverse root of the matrix

      if(verb) write(9,*) '... attempting inverse root'

      allocate(B(N,N))
      do i=1,N
         do j=i,N
            s=0.0d0
            do k=1,N
               if(case.eq.1) then
                  s=s+sqrt(E(k)) * A(i,k)*A(j,k)
               else
                  s=s+1.0d0/sqrt(E(k)) * A(i,k)*A(j,k)
               end if
            end do
            B(i,j)=s ; B(j,i)=B(i,j)
         end do
      end do
      A=B
      deallocate(B)

      if(verb) write(9,*) '... inverse root succeeded'
      
end subroutine inverse_root

subroutine Inverse_complex(N,A,INFO)
!........................................................................
! Calculates an inverse of a general complex matrix A
!
! INPUT
!  A -  N x N matrix 
!  N - diumension of A
!........................................................................
!  OUTPUT
!  A   - inverse of A [i.e,, the original A will be destroyed]
! INFO =  0 - success
!      = -i - the i-th argument of ZGETRI is illegal
!      =  i - the matrix is singular, no iverse possible
!........................................................................
  implicit none
  complex*16,dimension(:,:) ::  A(N,N)
  complex*16,dimension(:),allocatable ::  WORK
  integer,dimension(:),allocatable :: IPIV
  integer :: INFO
  integer :: N,LWORK !i,j,k


      allocate(IPIV(N),WORK(1))
  
!_______  perform LU factorisation

      call ZGETRF(N,N,A,N,IPIV,info)
      
!_______ determine the dimension of WORK

      call zgetri(N,A,N,IPIV,WORK,-1,INFO)
      LWORK=WORK(1) 
      deallocate(WORK)

      allocate(WORK(LWORK))
      call zgetri(N,A,N,IPIV,WORK,LWORK,INFO)
      deallocate(WORK,IPIV)
      if(INFO.ne.0) then
         write(9,*) 'ZGENRI: INFO = ',INFO ; stop
      end if

end subroutine Inverse_complex
    
   

