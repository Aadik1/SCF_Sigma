program parallel_test_1
use OMP_LIB 
implicit none

integer :: i, partial, total

!$OMP PARALLEL PRIVATE(partial) SHARED(total)
partial = 0.d0; total = 0.d0

!$OMP DO
do i = 0,1000
   partial = partial + i
end do
!$OMP END DO

!$OMP CRITICAL
total = total+partial 
!$OMP END CRITICAL

!$OMP END PARALLEL

print *, 'Total Sum of numbers from 1 to 1000:', total
end program parallel_test_1
