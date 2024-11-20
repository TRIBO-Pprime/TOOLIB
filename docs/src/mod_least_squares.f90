!< author: Arthur Francisco
!< version: 1.0.0
!< date: july, 23 2018
!<
!< <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<    **Least squares, linear and non linear**
!< </span>

module least_squares
use data_arch, only : I4, R8
use cholesky, only : choldc, cholsl
implicit none

private

logical(kind=I4), parameter :: verbose = .false.

public :: moindres_carres, moindres_carres_lineaire

contains

   subroutine moindres_carres(nb_var, nb_pts, hij, vec_xy, beta, f, df, typ, eps, relax, nb_var_der, info)
   !================================================================================================
   !< @note Function that returns the parameters of a function that approximates a data set. The parameters
   !<        determination is achieved by non linear least squares approximation.
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in   )                             :: nb_var      !! *number of parameters to be determined*
   integer(kind=I4), intent(in   )                             :: nb_pts      !! *number of points for function evaluation*
   integer(kind=I4), intent(in   )                             :: nb_var_der  !! *number of derivatives*
   real   (kind=R8), intent(in   ), dimension(1:nb_pts)        :: hij         !! *vector of evaluation points*
   real   (kind=R8), intent(in   ), dimension(1:nb_pts, 1:2)   :: vec_xy      !! *x and y coordinates of evaluation points*
   real   (kind=R8), intent(inout), dimension(1:nb_var)        :: beta        !! *parameters vector*
   real   (kind=R8), intent(in   )                             :: eps         !! *stop criterion*
   real   (kind=R8), intent(in   )                             :: relax       !! *relaxation parameter*
   character(len=*), intent(in   )                             :: typ         !! *kind of function used*
   integer(kind=I4), intent(out  )                             :: info        !! *information from Cholesky resolution*

   interface
      real(kind=R8) function f(xi, yi, var, nb_var, typ)
         use data_arch, only : I4, R8
         implicit none
         integer(kind=I4), intent(in   )                      :: nb_var
         real   (kind=R8), intent(in   )                      :: xi
         real   (kind=R8), intent(in   )                      :: yi
         real   (kind=R8), intent(inout), dimension(1:nb_var) :: var
         character(len=*), intent(in)                         :: typ
      endfunction f
      real(kind=R8) function df(xi, yi, var, nb_var, ivar, typ)
         use data_arch, only : I4, R8
         implicit none
         integer(kind=I4), intent(in   )                      :: nb_var
         integer(kind=I4), intent(in   )                      :: ivar
         real   (kind=R8), intent(in   )                      :: xi
         real   (kind=R8), intent(in   )                      :: yi
         real   (kind=R8), intent(inout), dimension(1:nb_var) :: var
         character(len=*), intent(in)                         :: typ
      endfunction df
   endinterface

      integer(kind=I4) :: i, ii, j, compteur
      integer(kind=I4) :: n, pseudo_newton
      real   (kind=R8) :: xi, yi, m_abs_scd

      real   (kind=R8),              dimension(1:nb_var_der)               :: scd
      real   (kind=R8),              dimension(1:nb_var_der, 1:nb_var_der) :: mat, sav_mat
      real   (kind=R8), allocatable, dimension(:,:)                        :: Jf, JfT
      real   (kind=R8), allocatable, dimension(:)                          :: vec_f, r, pt

      n = nint( sqrt(real(nb_pts, kind=R8)) )

      allocate(  Jf(1:nb_pts, 1:nb_var_der) )
      allocate( JfT(1:nb_var_der, 1:nb_pts) )
      allocate( vec_f(1:nb_pts) )
      allocate(     r(1:nb_pts) )
      allocate( pt(1:nb_var_der))

      compteur      = 0
      pseudo_newton = 10

w:    do
         do ii = 1, nb_pts
            xi = vec_xy(ii, 1)
            yi = vec_xy(ii, 2)
            vec_f(ii) = f(xi = xi,              & !
                          yi = yi,              & !
                         var = beta(1:nb_var),  & !
                      nb_var = nb_var,          & !
                         typ = typ)               !
         enddo

         if (mod(compteur, pseudo_newton)==0) then
            ! http://en.wikipedia.org/wiki/Least_squares
            do j = 1, nb_var_der
               do ii = 1, nb_pts
                  xi = vec_xy(ii, 1)
                  yi = vec_xy(ii, 2)
                  Jf(ii, j) = df(xi = xi,             & !
                                 yi = yi,             & !
                                var = beta(1:nb_var), & !
                             nb_var = nb_var,         & !
                               ivar = j,              & !
                                typ = typ)              !
               enddo
            enddo
            JfT(1:nb_var_der, 1:nb_pts) = transpose( Jf(1:nb_pts, 1:nb_var_der) )
            sav_mat(1:nb_var_der, 1:nb_var_der) = matmul( JfT(1:nb_var_der, 1:nb_pts), Jf(1:nb_pts, 1:nb_var_der) )
         endif
         r(1:nb_pts) = hij(1:nb_pts) -vec_f(1:nb_pts)
         scd(1:nb_var_der) = matmul( JfT(1:nb_var_der, 1:nb_pts), r(1:nb_pts) )

         mat(1:nb_var_der, 1:nb_var_der) = sav_mat(1:nb_var_der, 1:nb_var_der)
         forall(i=1:nb_var_der) mat(i,i) = 1.5*mat(i,i) ! type levenberg-marquardt

         call choldc(a = mat(1:nb_var_der, 1:nb_var_der),   & !
                     n = nb_var_der,                        & !
                    np = nb_var_der,                        & !
                     p = pt(1:nb_var_der),                  & !
                  info = info)                                !

         call cholsl(a = mat(1:nb_var_der, 1:nb_var_der),   & !
                     n = nb_var_der,                        & !
                    np = nb_var_der,                        & !
                     p = pt(1:nb_var_der),                  & !
                     b = scd(1:nb_var_der),                 & !
                     x = scd(1:nb_var_der),                 & !
                  info = info)                                !

         beta(1:nb_var_der) = beta(1:nb_var_der) +relax*scd(1:nb_var_der)
         m_abs_scd = 100*maxval( abs(scd(1:nb_var_der)/abs(beta(1:nb_var_der))) )

         if ((m_abs_scd) < eps) exit w
         compteur = compteur +1
         if (mod(compteur,1000) == 0) then
            if (verbose) write(*,*) compteur, m_abs_scd
         endif

      enddo w

      deallocate( Jf, JfT, vec_f, r, pt )

   return
   endsubroutine moindres_carres


   subroutine moindres_carres_lineaire(nb_var, nb_pts, hij, beta, Jf)
   !================================================================================================
   !< @note Function that returns the parameters of a function that approximates a data set. The parameters
   !<        determination is achieved by linear least squares approximation.
   !<
   !<  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                                  :: nb_var   !! *number of parameters to be determined*
   integer(kind=I4), intent(in )                                  :: nb_pts   !! *number of points for function evaluation*
   real   (kind=R8), intent(in ), dimension(1:nb_pts)             :: hij      !! *vector of evaluation points*
   real   (kind=R8), intent(in ), dimension(1:nb_pts, 1:nb_var)   :: Jf       !! *Jacobian*
   real   (kind=R8), intent(out), dimension(1:nb_var)             :: beta     !! *parameters vector*

      integer(kind=I4) :: info

      real(kind=R8),              dimension(1:nb_var)           :: scd, pt
      real(kind=R8),              dimension(1:nb_var, 1:nb_var) :: mat
      real(kind=R8), allocatable, dimension(:,:)                :: JfT

      allocate( JfT(1:nb_var, 1:nb_pts) )

      JfT(1:nb_var, 1:nb_pts) = transpose(  Jf(1:nb_pts, 1:nb_var) )
      mat(1:nb_var, 1:nb_var)  =   matmul( JfT(1:nb_var, 1:nb_pts),  Jf(1:nb_pts, 1:nb_var) )
      scd(1:nb_var)            =   matmul( JfT(1:nb_var, 1:nb_pts), hij(1:nb_pts) )

      call choldc(a = mat(1:nb_var, 1:nb_var),  & !
                  n = nb_var,                   & !
                 np = nb_var,                   & !
                  p = pt(1:nb_var),             & !
               info = info)                       !

      call cholsl(a = mat(1:nb_var, 1:nb_var),  & !
                  n = nb_var,                   & !
                 np = nb_var,                   & !
                  p =  pt(1:nb_var),            & !
                  b = scd(1:nb_var),            & !
                  x = scd(1:nb_var),            & !
                  info = info)                    !

      beta(1:nb_var) = scd(1:nb_var)

      deallocate( JfT )

   return
   endsubroutine moindres_carres_lineaire

endmodule least_squares
