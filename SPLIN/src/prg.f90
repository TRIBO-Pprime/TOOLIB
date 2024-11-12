program test_bspline
use bspline
use data_arch, only : I4, R8, PI_R8
implicit none

   integer(kind=I4), parameter ::  nx = 128,  ny =  256
   integer(kind=I4), parameter :: nnx = 700, nny = 1400

   integer(kind=I4) :: deg

   real(kind=R8) :: mad

   deg = 2
   call interp_surf(nx = nx, ny = ny, nnx = nnx, nny = nny, deg = deg, mad = mad)
   write(*, *) 'Deg, Max. abs. diff. :', deg, mad

   deg = 3
   call interp_surf(nx = nx, ny = ny, nnx = nnx, nny = nny, deg = deg, mad = mad)
   write(*, *) 'Deg, Max. abs. diff. :', deg, mad

   deg = 4
   call interp_surf(nx = nx, ny = ny, nnx = nnx, nny = nny, deg = deg, mad = mad)
   write(*, *) 'Deg, Max. abs. diff. :', deg, mad

   deg = 5
   call interp_surf(nx = nx, ny = ny, nnx = nnx, nny = nny, deg = deg, mad = mad)
   write(*, *) 'Deg, Max. abs. diff. :', deg, mad


stop
contains

   subroutine interp_surf(nx, ny, nnx, nny, deg, mad)
   implicit none
   integer(kind=I4), intent(in ) :: nx, ny, nnx, nny, deg
   real(kind=R8),    intent(out) :: mad

      integer(kind=I4) :: i, j, inbvx, inbvy, iloy, iflag
      real(kind=R8)    :: val

      real(kind=R8), dimension(1: nx, 1: ny) :: coeff, tab
      real(kind=R8), dimension(1:nnx, 1:nny) :: tab_ref, tab_int

      real(kind=R8), dimension(1:(nx + deg)) :: tx
      real(kind=R8), dimension(1:(ny + deg)) :: ty

      real(kind=R8), dimension(1: nx) :: x
      real(kind=R8), dimension(1: ny) :: y

      real(kind=R8), dimension(1:nnx) :: xx
      real(kind=R8), dimension(1:nny) :: yy

      x(1:nx) = [(-1. + (i - 1) * 2. / (nx - 1), i = 1, nx)]
      y(1:ny) = [(-1. + (j - 1) * 2. / (ny - 1), j = 1, ny)]

      do j = 1, ny
      do i = 1, nx
         tab(i, j) = func( x(i), y(j) )
      enddo
      enddo

      iflag = 0
      call db2ink(   x = x(1:nx),                     &  ! Array of x abcissae. Must be strictly increasing.
                    nx = nx,                          &  ! Number of x abcissae
                     y = y(1:ny),                     &  ! Array of y abcissae. Must be strictly increasing.
                    ny = ny,                          &  ! Number of y abcissae
                   fcn = tab(1:nx, 1:ny),             &  ! Array of function values to interpolate. fcn(i,j) should
                                                         !    contain the function value at the point (x(i),y(j))
                    kx = deg,                         &  ! The order of spline pieces in x (>= 2, < nx). (order = polynomial degree + 1)
                    ky = deg,                         &  ! The order of spline pieces in y (>= 2, < ny). (order = polynomial degree + 1)
                    tx = tx(1:(nx + deg)),            &  ! The knots in the x direction for the spline interpolant.
                                                         !    If iflag=0 these are chosen by [[db2ink]].
                                                         !    If iflag=1 these are specified by the user.
                                                         !    Must be non-decreasing.
                    ty = ty(1:(ny + deg)),            &  ! The knots in the y direction for the spline interpolant.
                                                         !    If iflag=0 these are chosen by [[db2ink]].
                                                         !    If iflag=1 these are specified by the user.
                                                         !    Must be non-decreasing.
                 bcoef = coeff(1:nx, 1:ny),           &  ! Array of coefficients of the b-spline interpolant.
                 iflag = iflag)                          ! **on input:**  0 = knot sequence chosen by [[db2ink]].
                                                         !                1 = knot sequence chosen by user.
                                                         ! **on output:** 1 = successful execution.
                                                         !                2 = iflag out of range.
                                                         !                3 = nx out of range.
                                                         !                4 = kx out of range.
                                                         !                5 = x not strictly increasing.
                                                         !                6 = tx not non-decreasing.
                                                         !                7 = ny out of range.
                                                         !                8 = ky out of range.
                                                         !                9 = y not strictly increasing.
                                                         !               10 = ty not non-decreasing.

      if (iflag/=1) error stop 'error calling db2ink'

      xx(1:nnx) = [(-1. + (i - 1) * 2. / (nnx - 1), i = 1, nnx)]
      yy(1:nny) = [(-1. + (j - 1) * 2. / (nny - 1), j = 1, nny)]

      inbvx = 1
      inbvy = 1
      iloy  = 1

      do j = 1, nny

         do i = 1, nnx

            call db2val(xval = xx(i),              &  ! xval     !! x coordinate of evaluation point.
                        yval = yy(j),              &  ! yval     !! y coordinate of evaluation point.
                         idx = 0,                  &  ! idx      !! x derivative of piecewise polynomial to evaluate.
                         idy = 0,                  &  ! idy      !! y derivative of piecewise polynomial to evaluate.
                          tx = tx(1:(nx + deg)),   &  ! tx       !! sequence of knots defining the piecewise polynomial in the x direction. (same as in last call to [[db2ink]])
                          ty = ty(1:(ny + deg)),   &  ! ty       !! sequence of knots defining the piecewise polynomial in the y direction. (same as in last call to [[db2ink]])
                          nx = nx,                 &  ! nx       !! the number of interpolation points in x. (same as in last call to [[db2ink]])
                          ny = ny,                 &  ! ny       !! the number of interpolation points in y. (same as in last call to [[db2ink]])
                          kx = deg,                &  ! kx       !! order of polynomial pieces in x. (same as in last call to [[db2ink]])
                          ky = deg,                &  ! ky       !! order of polynomial pieces in y. (same as in last call to [[db2ink]])
                       bcoef = coeff(1:nx, 1:ny),  &  ! bcoef    !! the b-spline coefficients computed by [[db2ink]].
                           f = val,                &  ! f        !! interpolated value &
                       iflag = iflag,              &  ! iflag    !! status flag: 0 : no errors, /=0 : error
                       inbvx = inbvx,              &  ! inbvx    !! initialization parameter which must be set to 1 the first time this routine is called, and must not be changed by the user.
                       inbvy = inbvy,              &  ! inbvy    !! initialization parameter which must be set to 1 the first time this routine is called, and must not be changed by the user.
                        iloy = iloy)                  ! iloy     !! initialization parameter which must be set to 1 the first time this routine is called, and must not be changed by the user.

            tab_int(i, j) = val

         enddo
      enddo

      do j = 1, nny
      do i = 1, nnx
         tab_ref(i, j) = func( xx(i), yy(j) )
      enddo
      enddo

      mad = maxval( abs(tab_ref - tab_int ) )

   return
   endsubroutine interp_surf

   pure elemental real(kind=R8) function func(xi, yj)
   implicit none
     real(kind=R8), intent(in) :: xi, yj

     func = sin(10 * PI_R8 * xi) * sin(10 * PI_R8 * yj)

   return
   endfunction func

endprogram test_bspline
