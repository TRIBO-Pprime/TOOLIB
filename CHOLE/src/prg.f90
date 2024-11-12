program test_utils
use data_arch, only : I4, R8
use cholesky,  only : choldc, cholsl
implicit none
integer(kind = I4), parameter :: tai = 100
real(kind=R8), dimension(tai, tai) :: a, aa
real(kind=R8), dimension(tai)      :: b, x, p
integer(kind = I4) :: i, info

   call random_number(a(1:tai,1:tai))
   call random_number(b(1:tai))

   a = matmul(a, transpose(a))

   forall(i = 1:tai) a(i, i) = 100*abs(a(i, i))

   aa(1:tai, 1:tai) = a(1:tai, 1:tai)

   call choldc( a    = a(1:tai,1:tai), &  !
                n    = tai,            &  !
                np   = tai,            &  !
                p    = p(1:tai),       &  !
                info = info )             !

   call cholsl( a    = a(1:tai,1:tai), &  !
                n    = tai,            &  !
                np   = tai,            &  !
                p    = p(1:tai),       &  !
                b    = b(1:tai),       &  !
                x    = x(1:tai),       &  !
                info = info )             !

   write(*,*) sum( abs(matmul(aa, x)-b) ), info

stop
endprogram test_utils
