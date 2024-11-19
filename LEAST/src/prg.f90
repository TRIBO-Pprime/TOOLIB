!< author: Arthur Francisco
!< version: 1.0
!< date: 15 mai 2012
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!< **Least squares, linear and non linear. Example of use.**
!< </span>
program test_least
use data_arch, only : I4, R8, PI_R8
use least_squares, only : moindres_carres, moindres_carres_lineaire
implicit none

integer(kind=I4) :: n1, n2, info

real(kind=R8), dimension(:), allocatable :: hij, beta

real(kind=R8), dimension(:,:), allocatable :: vec_xy, Jf

   n1 = 128
   n2 = 128

   allocate( hij(1:n1 * n2) )
   allocate( vec_xy(1:n1 * n2, 1:2) )

   allocate( beta(1:3) )

   call calc_imp_acf( long    = n1,                            &  !
                      larg    = n2,                            &  !
                      tau1    =  0.3_R8,                       &  !
                      tau2    =  0.1_R8,                       &  !
                      alpha   = -0.2_R8,                       &  !
                      ang     = PI_R8 / 6.,                    &  !
                      vec_xy  = vec_xy(1:n1 * n2, 1:2),        &  !
                      vec_acf = hij(1:n1 * n2) )                  !

   beta(1:3) = [0.1, 0.1, 0.1]

   call moindres_carres( nb_var     = 3,                       &  !
                         nb_pts     = n1 * n2,                 &  !
                         hij        = hij(1:n1 * n2),          &  !
                         vec_xy     = vec_xy(1:n1 * n2, 1:2),  &  !
                         beta       = beta(1:3),               &  !
                         f          = f,                       &  !
                         df         = df,                      &  !
                         typ        = 'no_type',               &  !
                         eps        = 0.001_R8,                &  !
                         relax      = 0.5_R8,                  &  !
                         nb_var_der = 3,                       &  !
                         info       = info )                      !

   write(*,*) 'Non linear approach to a non linear problem. Solution is 0.3 0.1 +/-0.5'
   write(*,*) beta(1), beta(2), sin(mod(beta(3), 2*PI_R8))

   !=======================================================

   allocate( Jf(1:n1 * n2, 1:3) )

   beta(1:3) = [0.1, 0.1, 0.1]

   call calc_Jf( long   = n1,                                  &  !
                 larg   = n2,                                  &  !
                 nb_var = 3,                                   &  !
                 vec_xy = vec_xy(1:n1 * n2, 1:2),              &  !
                 var    = beta(1:3),                           &  !
                 Jf     = Jf(1:n1 * n2, 1:3) )                    !

   call moindres_carres_lineaire( nb_var = 3,                  &  !
                                  nb_pts = n1 * n2,            &  !
                                  hij    = hij(1:n1 * n2),     &  !
                                  beta   = beta(1:3),          &  !
                                  Jf     = Jf(1:n1 * n2, 1:3) )   !

   write(*,*) 'Linear approach to a non linear problem. Solution is 0.3 0.1 +/-0.5'
   write(*,*) beta(1), beta(2), sin(mod(beta(3), 2*PI_R8))

   !=======================================================

   deallocate( hij, vec_xy, beta, Jf )

contains

   real(kind=R8) function f(xi, yi, var, nb_var, typ)
   !! Kind of particular 2D autocorrelation function
   use data_arch, only : I4, R8
   implicit none
   integer(kind=I4), intent(in   )                      :: nb_var !! *number of parameters*
   real   (kind=R8), intent(in   )                      :: xi     !! *x coordinates*
   real   (kind=R8), intent(in   )                      :: yi     !! *y coordinates*
   real   (kind=R8), intent(inout), dimension(1:nb_var) :: var    !! *parameter vector*
   character(len=*), intent(in)                         :: typ    !! *not used here*

      f = autocov_impo(xi=xi, xj=yi, tau1=var(1), tau2=var(2), alpha=-0.2_R8, ang=var(3))

   return
   endfunction f


   real(kind=R8) function df(xi, yi, var, nb_var, ivar, typ)
   !! Kind of particular 2D autocorrelation function. Numerical derivatives.
   use data_arch, only : I4, R8
   implicit none
   integer(kind=I4), intent(in   )                      :: nb_var !! *number of parameters*
   integer(kind=I4), intent(in   )                      :: ivar   !! *ith parameter*
   real   (kind=R8), intent(in   )                      :: xi     !! *x coordinates*
   real   (kind=R8), intent(in   )                      :: yi     !! *y coordinates*
   real   (kind=R8), intent(inout), dimension(1:nb_var) :: var    !! *parameter vector*
   character(len=*), intent(in)                         :: typ    !! *not used here*

      real(kind=R8) :: v1, v2, v3, dv1, dv2, dv3

      v1 = var(1) ; dv1 = 0.001 !v1 / 100
      v2 = var(2) ; dv2 = 0.001 !v2 / 100
      v3 = var(3) ; dv3 = 0.100 !v3 / 100

      select case(ivar)

         case(1)

            df = ( autocov_impo(xi=xi, xj=yi, tau1=v1 + 0.5 * dv1, tau2=v2, alpha=-0.2_R8, ang=v3) -        &  !
                   autocov_impo(xi=xi, xj=yi, tau1=v1 - 0.5 * dv1, tau2=v2, alpha=-0.2_R8, ang=v3) ) / dv1     !

         case(2)

            df = ( autocov_impo(xi=xi, xj=yi, tau1=v1, tau2=v2 + 0.5 * dv2, alpha=-0.2_R8, ang=v3) -        &  !
                   autocov_impo(xi=xi, xj=yi, tau1=v1, tau2=v2 - 0.5 * dv2, alpha=-0.2_R8, ang=v3) ) / dv2     !

         case(3)

            df = ( autocov_impo(xi=xi, xj=yi, tau1=v1, tau2=v2, alpha=-0.2_R8, ang=v3 + 0.5 * dv3) -        &  !
                   autocov_impo(xi=xi, xj=yi, tau1=v1, tau2=v2, alpha=-0.2_R8, ang=v3 - 0.5 * dv3) ) / dv3     !

      endselect

   return
   endfunction df


   subroutine calc_Jf(long, larg, nb_var, vec_xy, var, Jf)
   !! Determine the Jacobian matrix of f
   implicit none
   integer(kind=I4), intent(in) :: long     !! *number of points along x*
   integer(kind=I4), intent(in) :: larg     !! *number of points along y*
   integer(kind=I4), intent(in) :: nb_var   !! *number of parameters to be determined*
   real   (kind=R8), intent(in   ), dimension(1:long * larg, 1:2)      :: vec_xy   !! *x and y coordinates of evaluation points*
   real   (kind=R8), intent(out  ), dimension(1:long * larg, 1:nb_var) :: Jf       !! *Jacobian matrix*
   real   (kind=R8), intent(inout), dimension(1:nb_var)                :: var      !! *parameters vector*

      integer(kind=I4) :: i, j, k, ivar

      k = 0
      do j = 1, larg
      do i = 1, long

         k = k + 1

         do ivar = 1, nb_var

            Jf(k, ivar) = df( xi     = vec_xy(k, 1),  &  !
                              yi     = vec_xy(k, 2),  &  !
                              var    = var(1:3),      &  !
                              nb_var = 3,             &  !
                              ivar   = ivar,          &  !
                              typ    = "no_type")        !
         enddo

      enddo
      enddo

   return
   endsubroutine calc_Jf


   real(kind=R8) function autocov_impo(xi, xj, tau1, tau2, alpha, ang)
   !================================================================================================
   !<@note Function that returns \( \exp \left(\alpha \sqrt{\left(\frac{x}{\tau_1}\right)^2+
   !<                                                        \left(\frac{y}{\tau_2}\right)^2}
   !<                                      \right) \)
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   real(kind=R8), intent(in) :: tau1   !! *correlation length along \(x\)*
   real(kind=R8), intent(in) :: tau2   !! *correlation length along \(y\)*
   real(kind=R8), intent(in) :: alpha  !! *log(z)* where *z* is often 0.2
   real(kind=R8), intent(in) :: xi     !! *\(x\) coordinate*
   real(kind=R8), intent(in) :: xj     !! *\(y\) coordinate*
   real(kind=R8), intent(in) :: ang    !! *angle* (rad)

      real(kind=R8) :: x, y

      x = +cos(ang) * xi + sin(ang) * xj
      y = -sin(ang) * xi + cos(ang) * xj

      autocov_impo = exp( alpha * sqrt( (x / tau1)**2 + (y / tau2)**2 ) )

   return
   endfunction autocov_impo


   subroutine calc_imp_acf(long, larg, tau1, tau2, alpha, ang, vec_acf, vec_xy)
   !================================================================================================
   !<@note Function that returns the theoretical autocorrelation function in an array.<br/>
   !< The autocorrelation function is supposed to be obtained from a real surface which must be periodic
   !< or nearly periodic (because of the use of FFTs).
   !< In addition, the surface is supposed to be 0 mean and normalized (\(\sigma = 1 \)),
   !< therefore *acf* is zero-mean and normalized so that its max value is 1.<br/>
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in) :: long   !! *surface acf width*
   integer(kind=I4), intent(in) :: larg   !! *surface acf height*
   real   (kind=R8), intent(in) :: tau1   !! *first correlation length*
   real   (kind=R8), intent(in) :: tau2   !! *surface second correlation length*
   real   (kind=R8), intent(in) :: alpha  !! *parameter that controls the expondential decrease*
   real   (kind=R8), intent(in) :: ang    !! *acf ellipsis angle*
   real   (kind=R8), dimension(1:long * larg),      intent(out) :: vec_acf  !! *resulting acf*
   real   (kind=R8), dimension(1:long * larg, 1:2), intent(out) :: vec_xy   !! *points coordinates*

      integer(kind=I4) :: i, j, k, long2, larg2
      real   (kind=R8) :: xi, xj, r

      long2 = long / 2 + 1
      larg2 = larg / 2 + 1

      k = 0
      do j = 1, larg
      do i = 1, long

         k = k + 1

         xi = real(i - long2, kind=R8) / long
         xj = real(j - larg2, kind=R8) / larg

         vec_xy(k, 1) = xi
         vec_xy(k, 2) = xj

         call random_number(r)

         r = 1. + 0.05 * (2 * (0.5 - r))

         vec_acf(k) = r * autocov_impo( xi    = xi,       &  ! IN
                                        xj    = xj,       &  ! IN
                                        tau1  = tau1,     &  ! IN
                                        tau2  = tau2,     &  ! IN
                                        alpha = alpha,    &  ! IN
                                        ang   = ang )        ! IN

      enddo
      enddo

   return
   endsubroutine calc_imp_acf


endprogram test_least
