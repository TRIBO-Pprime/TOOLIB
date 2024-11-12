program test_intpl
use data_arch, only : I4, R8, PI_R8
use intpl
implicit none

integer(kind=I4), parameter :: n1 = 1024, n2 = 2*n1

   write(*,*) ' ORDER 1' ; call test_interp_pond( ordre = 1 )
   write(*,*) ' ORDER 3' ; call test_interp_pond( ordre = 3 )
   write(*,*) ' ORDER 5' ; call test_interp_pond( ordre = 5 )
   write(*,*) ' ORDER 7' ; call test_interp_pond( ordre = 7 )

stop
contains

   !-!< --------------------------------------------------------------------     \n
   !-!< subroutine de test de l'interpolation/restriction                        \n
   !-!< --------------------------------------------------------------------
   subroutine test_interp_pond( ordre )
   implicit none
   integer(kind=I4), intent(in) :: ordre

      integer(kind=I4) :: i, ii, j, jj, istat1, istat2

      real(kind=R8), allocatable, dimension(:)   :: tableau1D, tab_f1D, tab_g1D, erreur1D
      real(kind=R8), allocatable, dimension(:,:) :: tableau2D, tab_f2D, tab_g2D, erreur2D

      type(tborne) :: bfin, bgro

      allocate( tableau1D(0:n2),  tab_f1D(0:n2),   &  !
                  tab_g1D(0:n1), erreur1D(0:n2),   &  !
                  stat = istat1 )                     !

      allocate( tableau2D(0:n2, 0:n2),  tab_f2D(0:n2, 0:n2),   &  !
                  tab_g2D(0:n1, 0:n1), erreur2D(0:n2, 0:n2),   &  !
                  stat = istat2 )                                 !

      if ( istat1 + istat2 /=0 ) stop '"test_interp_pond" problem of allocation'

      !-!<.............................................................
      do ii = 0, n1
         i = 2*ii
         tableau1D(ii) = def_tab1D(i) !-!< initialisation à une fonction ressemblant à UN noyau élastique
      enddo

      call interp1D( tabgros = tableau1D(0:n1),    &  !
                     lb_gros = 0,                  &  !
                     tabfin  = tab_f1D(0:n2),      &  !
                     lb_fin  = 0,                  &  !
                     ub_gros = n1,                 &  !
                     ordre   = ordre )                !-!< interpolation de ce tableau

      do i = 0, n2
         erreur1D(i) = 200*abs( ( def_tab1D(i) - tab_f1D(i) ) / &  !
                                ( def_tab1D(i) + tab_f1D(i) ) )    !-!< comparaison théorique/interpolé
      enddo
      write(*,*) 'interp 1D - max err, mean err: :         ', maxval( erreur1D(0:n2) ), sum( erreur1D(0:n2) ) / (n2 + 1)

      !-!<.............................................................

      do i = 0, n2
         tableau1D(i) = def_tab1D(i) !-!< initialisation à une fonction ressemblant à UN noyau élastique
      enddo

      call restrict1D( tabfin  = tableau1D(0:n2),     &  !
                       lb_fin  = 0,                   &  !
                       tabgros = tab_g1D(0:n1),       &  !
                       lb_gros = 0,                   &  !
                       ub_gros = n1,                  &  !
                       ordre   = ordre )                 !-!< restriction de ce tableau

      do ii = 0, n1
         erreur1D(ii) = 200*abs( ( def_tab1D(2*ii) - tab_g1D(ii) ) / &  !
                                 ( def_tab1D(2*ii) + tab_g1D(ii) ) )    !-!< comparaison théorique/interpolé
      enddo
      write(*,*) 'weight 1D - max err, mean err:           ',  maxval( erreur1D(0:n1) ), sum( erreur1D(0:n1) ) / (n1 + 1)

      !-!<.............................................................

      do i = 0, n2
         tableau1D(i) = def_tab1D(i) !-!< initialisation à une fonction ressemblant à UN noyau élastique
      enddo

      call restrict1D( tabfin  = tableau1D(0:n2),  &  !
                       lb_fin  = 0,                &  !
                       tabgros = tab_g1D(0:n1),    &  !
                       lb_gros = 0,                &  !
                       ub_gros = n1,               &  !
                       ordre   = ordre )              !-!< restriction de ce tableau

      call interp1D( tabgros = tab_g1D(0:n1),      &  !
                     lb_gros = 0,                  &  !
                     tabfin  = tab_f1D(0:n2),      &  !
                     lb_fin  = 0,                  &  !
                     ub_gros = n1,                 &  !
                     ordre   = ordre )                !-!< interpolation de ce tableau

      do i = 0, n2
         erreur1D(i) = 200*abs( ( tableau1D(i) - tab_f1D(i) ) / &  !
                                ( tableau1D(i) + tab_f1D(i) ) )    !-!< comparaison avant/après
      enddo
      write(*,*) 'weight + interp 1D - max err, mean err:  ', maxval( erreur1D(0:n2) ), sum( erreur1D(0:n2) ) / (n2 + 1)

      !-!<.............................................................
      !-!< Idem 2D
      !-!<.............................................................

      do jj = 0, n1
      do ii = 0, n1
         i = 2*ii
         j = 2*jj
         tableau2D(ii, jj) = def_tab2D(i, j)
      enddo
      enddo

      bfin = tborne( lb1 = 0, ub1 = n2, lb2 = 0, ub2 = n2 )
      bgro = tborne( lb1 = 0, ub1 = n1, lb2 = 0, ub2 = n1 )

      call interp2D( tabgro = tableau2D(0:n1, 0:n1),  &  !
                     bgro   = bgro,                   &  !
                     tabfin = tab_f2D(0:n2, 0:n2),    &  !
                     bfin   = bfin,                   &  !
                     ordre  = ordre )                    !

      do j = 0, n2
      do i = 0, n2
         erreur2D(i, j) = 200*abs( ( def_tab2D(i,j) - tab_f2D(i,j) ) / &  !
                                   ( def_tab2D(i,j) + tab_f2D(i,j) ) )    !
      enddo
      enddo
      write(*,*) 'interp 2D - max err, mean err:           ',  maxval( erreur2D(0:n2, 0:n2) ), sum( erreur2D(0:n2, 0:n2) ) / (n2 + 1)**2

      !-!<.............................................................

      do j = 0, n2
      do i = 0, n2
         tableau2D(i, j) = def_tab2D(i, j)
      enddo
      enddo

      bfin = tborne( lb1 = 0, ub1 = n2, lb2 = 0, ub2 = n2 )
      bgro = tborne( lb1 = 0, ub1 = n1, lb2 = 0, ub2 = n1 )

      call restrict2D( tabfin  = tableau2D(0:n2, 0:n2),  &  !
                       bfin    = bfin,                   &  !
                       tabgros = tab_g2D(0:n1, 0:n1),    &  !
                       bgros   = bgro,                   &  !
                       ordre   = ordre )                    !

      do jj = 0, n1
      do ii = 0, n1
         erreur2D(ii, jj) = 200*abs( ( def_tab2D(2*ii, 2*jj) - tab_g2D(ii, jj) ) / &  !
                                     ( def_tab2D(2*ii, 2*jj) + tab_g2D(ii, jj) ) )    !
      enddo
      enddo

      write(*,*) 'weight 2D - max err, mean err:           ', maxval( erreur2D(0:n1, 0:n1) ), sum( erreur2D(0:n1, 0:n1) ) / (n1 + 1)**2

      !-!<.............................................................

      do j = 0, n2
      do i = 0, n2
         tableau2D(i, j) = def_tab2D(i, j)
      enddo
      enddo

      bfin = tborne( lb1 = 0, ub1 = n2, lb2 = 0, ub2 = n2 )
      bgro = tborne( lb1 = 0, ub1 = n1, lb2 = 0, ub2 = n1 )

      call restrict2D( tabfin  = tableau2D(0:n2, 0:n2),     &  !
                       bfin    = bfin,                      &  !
                       tabgros = tab_g2D(0:n1, 0:n1),       &  !
                       bgros   = bgro,                      &  !
                       ordre   = ordre )                       !

      call interp2D( tabgro = tab_g2D(0:n1, 0:n1),       &  !
                     bgro   = bgro,                      &  !
                     tabfin = tab_f2D(0:n2, 0:n2),       &  !
                     bfin   = bfin,                      &  !
                     ordre  = ordre )                       !

      do j = 0, n2
      do i = 0, n2
         erreur2D(i,j) = 200*abs( ( tableau2D(i, j) - tab_f2D(i, j) ) / &  !
                                  ( tableau2D(i, j) + tab_f2D(i, j ) ) )   !
      enddo
      enddo
      write(*,*) 'weight + interp 2D - max err, mean err:  ',  maxval( erreur2D(0:n2, 0:n2) ), sum( erreur2D(0:n2, 0:n2) ) / (n2 + 1)**2

      deallocate( tableau1D, tab_f1D, tab_g1D, erreur1D )
      deallocate( tableau2D, tab_f2D, tab_g2D, erreur2D )

   return
   endsubroutine test_interp_pond

   real(kind=R8) function def_tab1D(i)
   implicit none
   integer(kind=I4), intent(in) :: i

!~       def_tab1D = 1._R8 / ( 1._R8 + (i/64._R8)**2._R8 )
      def_tab1D = sin( 2 * i * PI_R8 / n1 ) + 1.1_R8

   return
   endfunction def_tab1D

   real(kind=R8) function def_tab2D(i, j)
   implicit none
   integer(kind=I4), intent(in) :: i, j

!~       def_tab2D = 1._R8 / (1._R8 + ( (i-j)/64._R8 )**2._R8 )
      def_tab2D = sin( 2 * i * PI_R8 / n1 ) * cos( 2 * j * PI_R8 / n1 ) + 1.1_R8

   return
   endfunction def_tab2D


endprogram test_intpl
