program test_tchebychev
use data_arch,     only : I4, R8, HIG_R8
use miscellaneous, only : get_unit
use tchebychev,    only : least_squares_tcheby, tcheby, coeff_tcheby_vers_monome

implicit none

integer(kind = I4) :: i, iu, nbvar

integer(kind = I4), parameter :: n1 = 1024, n2 = 512, i1 = 9, i2 = 11

real(kind=R8), allocatable, dimension(:,:) :: tab_height, tab_result
real(kind=R8), allocatable, dimension(:)   :: tab_coef1, tab_coef2

character(len=128) :: str

   nbvar = (i1 + 1) * (i2 + 1)

   allocate( tab_height(1:n1, 1:n2) ) ; tab_height(1:n1, 1:n2) = HIG_R8
   allocate( tab_result(1:n1, 1:n2) ) ; tab_result(1:n1, 1:n2) = HIG_R8

   allocate( tab_coef1(1:nbvar) ) ; tab_coef1(1:nbvar) = HIG_R8
   allocate( tab_coef2(1:nbvar) ) ; tab_coef2(1:nbvar) = HIG_R8

   call genere_surf_poly( long1    = n1,                        &  !
                          long2    = n2,                        &  !
                          deg1     = i1,                        &  !
                          deg2     = i2,                        &  !
                          tab_out  = tab_height(1:n1, 1:n2),    &  !
                          tab_coef = tab_coef1(1:nbvar) )          !

   call least_squares_tcheby( tab_in  = tab_height(1:n1, 1:n2),   &  !
                              tab_out = tab_result(1:n1, 1:n2),   &  !
                              long1   = n1,                       &  !
                              long2   = n2,                       &  !
                              nvarx   = i1,                       &  !
                              nvary   = i2,                       &  !
                              verif   = .true.,                   &  !
                         multi_thread = .true.)                      !

   write(*, *)
   write(*, *) 'If the surface is well approximated, the differences are very small.'
   write(*, *) 'Indeed, the approximated surface is a polynomial of the same degree along x, y as the approximating one.'
   write(*, *) 'Result :', sum( abs(tab_result - tab_height) )

   write(*, *)
   write(*, *) 'If the translation from Tchebychev to ordinary polynomial is well carried out,'
   write(*, *) 'there should be no big difference between the two surfaces.'

   call get_unit(iu)
   open(iu, file = 'verif_tcheby_vers_monome.txt')
      write(str, '(a,i3.3,a)') '(', nbvar,'E18.8,a)'
      read(iu, trim(str)) ( tab_coef2(i), i = 1, nbvar )
   close(iu)

   call convert_to_poly ( long1    = n1,                        &  !
                          long2    = n2,                        &  !
                          deg1     = i1,                        &  !
                          deg2     = i2,                        &  !
                          tab_coef = tab_coef2(1:nbvar),        &  !
                          tab_out  = tab_height(1:n1, 1:n2) )      !

   write(*, *) 'Max. abs. diff. :', maxval( abs((tab_result - tab_height)/(1._R8 + tab_height) ) )

   write(*, *)
   write(*, *) 'Given the fact that the reference surface is made of the product of two polynomials,'
   write(*, *) 'the coefficients found must be the same as the imposed ones.'
   write(*, *) 'Max. abs. diff. :', maxval( abs((tab_coef1 - tab_coef2)/tab_coef1) )

   deallocate( tab_height, tab_result, tab_coef1, tab_coef2 )

contains

   !==================================================================================================
   !> @brief         Génération d'une surface polynômiale pour vérification des procédures d'approximation
   !==================================================================================================
   subroutine convert_to_poly(long1, long2, deg1, deg2, tab_coef, tab_out)
   implicit none
   integer(kind=I4), intent(in )                                     :: long1    !! taille x
   integer(kind=I4), intent(in )                                     :: long2    !! taille y
   integer(kind=I4), intent(in )                                     :: deg1     !!
   integer(kind=I4), intent(in )                                     :: deg2     !!
   real   (kind=R8), intent(in ), dimension(1:(deg1+1) * (deg2+1))   :: tab_coef !! tableau résultant : surface
   real   (kind=R8), intent(out), dimension(1:long1, 1:long2)        :: tab_out  !! tableau résultant : surface

      real(kind=R8)    :: xi, xj
      integer(kind=I4) :: i, j, k1, k2, k1k2

      tab_out = 0._R8
      do j = 1, long2
         xj = -1. + (j - 1) * 2. / (long2 - 1)

         do i = 1, long1
            xi = -1. + (i - 1) * 2. / (long1 - 1)

            k1k2 = 0
            do k2 = 0, deg2
            do k1 = 0, deg1

               k1k2 = k1k2 + 1
               tab_out(i, j) = tab_out(i, j) + tab_coef(k1k2) * (xi ** k1) * (xj ** k2)

            enddo
            enddo
         enddo

      enddo
   return
   endsubroutine convert_to_poly


   subroutine genere_surf_poly(long1, long2, deg1, deg2, tab_out, tab_coef)
   implicit none
   integer(kind=I4), intent(in )                                     :: long1    !! taille x
   integer(kind=I4), intent(in )                                     :: long2    !! taille y
   integer(kind=I4), intent(in )                                     :: deg1     !!
   integer(kind=I4), intent(in )                                     :: deg2     !!
   real   (kind=R8), intent(out), dimension(1:long1, 1:long2)        :: tab_out  !! tableau résultant : surface
   real   (kind=R8), intent(out), dimension(1:(deg1+1) * (deg2+1))   :: tab_coef !! tableau résultant : surface

      real(kind=R8)    :: xi, xj
      integer(kind=I4) :: i, j, ij, k1, k2, k1k2

      real(kind=R8), dimension(0:deg1) :: ai
      real(kind=R8), dimension(0:deg2) :: aj

      real(kind=R8), dimension(0:deg1) :: ai_m

      real(kind=R8), dimension(1:long1) :: val_xi_t
      real(kind=R8), dimension(1:long1) :: val_xi_m

      call random_number( ai(0:deg1) )
      call random_number( aj(0:deg2) )

      ai = +9 * ai - 4._R8
      aj = -5 * aj + 8._R8
      ! =========================== SIMPLE CHECK =====================================================

      ! A combination of Tchebychev polynomials must
      ! yield the same results as a classical polynomial
      !--------------------------------------------
      val_xi_t = 0._R8
      do i = 1, long1
         xi = -1. + (i - 1) * 2. / (long1 - 1)

         do k1 = 0, deg1
            val_xi_t(i) = val_xi_t(i) + ai(k1) * tcheby(n = k1, x = xi)
         enddo
      enddo
      !--------------------------------------------
      ! Towards equivalent classical polynomial
      call coeff_tcheby_vers_monome(coeff_t = ai(0:deg1), coeff_m = ai_m(0:deg1), deg = deg1)

      val_xi_m = 0._R8
      do i = 1, long1
         xi = -1. + (i - 1) * 2. / (long1 - 1)

         do k1 = 0, deg1
            val_xi_m(i) = val_xi_m(i) + ai_m(k1) * (xi**k1)
         enddo
      enddo
      !--------------------------------------------
      write(*,*) 'Equivalence of Tchebychev and classical polynomials'
      write(*,*) 'Difference must be negligible'
      write(*,*) ' Diff = ', maxval( abs( val_xi_m - val_xi_t ) )
      !--------------------------------------------
      ! =========================== END SIMPLE CHECK ==================================================

      tab_coef = 0._R8
      ij = 0
      do j = 0, deg2
         do i = 0, deg1
            ij = ij + 1
            tab_coef(ij) = ai(i) * aj(j)
         enddo
      enddo

      tab_out = 0._R8
      do j = 1, long2
         xj = -1. + (j - 1)* 2. / (long2 - 1)

         do i = 1, long1
            xi = -1. + (i - 1) * 2. / (long1 - 1)

            k1k2 = 0
            do k2 = 0, deg2
            do k1 = 0, deg1

               k1k2 = k1k2 + 1
               tab_out(i, j) = tab_out(i, j) + tab_coef(k1k2) * (xi ** k1) * (xj ** k2)

            enddo
            enddo
         enddo

      enddo
   return
   endsubroutine genere_surf_poly

endprogram test_tchebychev
