module cholesky
use data_arch, only: I4, R8
implicit none

private

public :: choldc, cholsl

contains

   subroutine choldc (a, n, np, p, info)
   !================================================================================================
   !< @note Given a positive definite symmetric matrix a(1:n,1:n), with
   !        physical dimensions np, this routine constructs its Cholesky
   !        decomposition, A=L L^T. On input, only the upper triangle of
   !        a need to be given; it is not modified. The Cholesky factor L
   !        is returned in the lower triangle of a, except for its diagonal
   !        elements which are returned in p(1:n).
   !        (c) Numerical recipes
   !
   !  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in   )                    :: n     !! *system size*
   integer(kind=I4), intent(in   )                    :: np    !! *matrix size*
   integer(kind=I4), intent(  out)                    :: info  !! *information ouput*
   real   (kind=R8), intent(inout), dimension(np, np) :: a     !! *system matrix*
   real   (kind=R8), intent(  out), dimension(np)     :: p     !! *diagonal elements*

      integer(kind=I4) :: i, j, k
      real   (kind=R8) :: ssum

      info = 0
o:    do i = 1, n
         do j = i, n
            ssum = a(i,j)
            do k = i-1, 1, -1
               ssum = ssum -a(i,k)*a(j,k)
            enddo
            if(i==j) then
               if (ssum<=0.) then
                  info = 1
                  p(np) = 0._R8
                  p(1) = -ssum
                  exit o
               endif
               p(i) = sqrt(ssum)
            else
               a(j, i) = ssum/p(i)
            endif
         enddo
      enddo o
   return
   endsubroutine choldc


   subroutine cholsl (a, n, np, p, b, x, info)
   !================================================================================================
   !< @note Solves the set of linear equations A x = b, where A is a positive-
   !        definite symmetric matrix with physical dimensions np. A and P are
   !        are input as the output from choldc. Only the lower triangle of A
   !        is accessed. B(1:n) is inout as the right-hand side vector. The
   !        solution vector is returned in X(1:n). A, n, np, and P are not
   !        modified and can be left in place for successive calls with different
   !        right-hand sides B. B is not modified unless you identify B and X in
   !        the calling sequence, which is allowed.
   !        (c) after Numerical recipes
   !
   !  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in   )                    :: n     !! *system size*
   integer(kind=I4), intent(in   )                    :: np    !! *matrix size*
   integer(kind=I4), intent(inout)                    :: info  !! *information ouput*
   real   (kind=R8), intent(in   ), dimension(np, np) :: a     !! *system matrix*
   real   (kind=R8), intent(in   ), dimension(np)     :: p     !! *diagonal elements*
   real   (kind=R8), intent(inout), dimension(np)     :: b     !! *rhs*
   real   (kind=R8), intent(inout), dimension(np)     :: x     !! *solution vector*

      integer(kind=I4) :: i, k
      real   (kind=R8) :: ssum

      if (info==1) return
      do i = 1, n
         ssum = b(i)
         do k = i-1, 1, -1
            ssum = ssum -a(i,k)*x(k)
         enddo
         x(i) = ssum/p(i)
      enddo
      do i = n, 1, -1
         ssum = x(i)
         do k = i+1,n
            ssum = ssum -a(k,i)*x(k)
         enddo
         x(i) = ssum/p(i)
      enddo
   return
   endsubroutine cholsl

endmodule cholesky

!~    subroutine solve_tridiag(a, b, c, d, x, n)
!~    implicit none
!~    integer(kind=I4), intent(in )                 :: n !   n - number of equations
!~    real   (kind=R8), intent(in ), dimension(1:n) :: a !   a - sub-diagonal (means it is the diagonal below the main diagonal)
!~    real   (kind=R8), intent(in ), dimension(1:n) :: b !   b - the main diagonal
!~    real   (kind=R8), intent(in ), dimension(1:n) :: c !   c - sup-diagonal (means it is the diagonal above the main diagonal)
!~    real   (kind=R8), intent(in ), dimension(1:n) :: d !   d - right part
!~    real   (kind=R8), intent(out), dimension(1:n) :: x !   x - the answer

!~       real(kind=R8), dimension(1:n) :: cp, dp

!~       real   (kind=R8) :: m
!~       integer(kind=I4) :: i

!~       ! initialize c-prime and d-prime
!~       cp(1) = c(1)/b(1)
!~       dp(1) = d(1)/b(1)
!~       ! solve for vectors c-prime and d-prime
!~       do i = 2, n
!~          m = b(i) -cp(i-1)*a(i)
!~          cp(i) = c(i)/m
!~          dp(i) = (d(i)-dp(i-1)*a(i))/m
!~       enddo
!~       ! initialize x
!~       x(n) = dp(n)
!~       ! solve for x from the vectors c-prime and d-prime
!~       do i = n -1, 1, -1
!~          x(i) = dp(i) -cp(i)*x(i+1)
!~       enddo

!~    return
!~    endsubroutine solve_tridiag
