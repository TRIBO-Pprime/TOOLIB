!------------------------------------------------------------------------------------------
!-!> @file       mod_interp_pond.f90
!-!> @authors    Arthur Francisco
!-!> @version    1.0
!-!> @date       15 mai 2012
!-!> @brief      Module d'interpolation et restriction
!-!> @details    Ce module implémente des interpolations/pondérations sur des points régulièrement espacés                         \n
!-!>             La technique utilisées est celle de Lagrange qui présente l'inconvénient d'être instable pour les ordres élevés   \n
!-!>             (instabilité de Runge) Des polynômes de Tchebychev auraient été plus efficaces mais ils nécessitent des points    \n
!-!>             localisés différemment.
!------------------------------------------------------------------------------------------
module intpl
use data_arch, only : I4, R8
implicit none

private

!-!< coefficient d'interpolation des ordres 1, 3, 5, 7 obtenus par les polynômes de Lagrange
!-!< "c" pour coefficient, "i" interpolation, "1" "3" "5" "7" pour l'ordre
real(kind=R8), parameter ::   ci1_00 = +1.0_R8/2.0_R8
real(kind=R8), parameter ::   ci1_01 = +1.0_R8/2.0_R8

real(kind=R8), parameter ::   ci3_00 = -1.0_R8/16.0_R8
real(kind=R8), parameter ::   ci3_01 = +9.0_R8/16.0_R8
real(kind=R8), parameter ::   ci3_02 = +9.0_R8/16.0_R8
real(kind=R8), parameter ::   ci3_03 = -1.0_R8/16.0_R8

real(kind=R8), parameter ::   ci5_00 = +  3.0_R8/256.0_R8
real(kind=R8), parameter ::   ci5_01 = - 25.0_R8/256.0_R8
real(kind=R8), parameter ::   ci5_02 = +150.0_R8/256.0_R8
real(kind=R8), parameter ::   ci5_03 = +150.0_R8/256.0_R8
real(kind=R8), parameter ::   ci5_04 = - 25.0_R8/256.0_R8
real(kind=R8), parameter ::   ci5_05 = +  3.0_R8/256.0_R8

real(kind=R8), parameter ::   ci7_00 = -   5.0_R8/2048.0_R8
real(kind=R8), parameter ::   ci7_01 = +  49.0_R8/2048.0_R8
real(kind=R8), parameter ::   ci7_02 = - 245.0_R8/2048.0_R8
real(kind=R8), parameter ::   ci7_03 = +1225.0_R8/2048.0_R8
real(kind=R8), parameter ::   ci7_04 = +1225.0_R8/2048.0_R8
real(kind=R8), parameter ::   ci7_05 = - 245.0_R8/2048.0_R8
real(kind=R8), parameter ::   ci7_06 = +  49.0_R8/2048.0_R8
real(kind=R8), parameter ::   ci7_07 = -   5.0_R8/2048.0_R8

!-!< coefficient de pondération des ordres 0, 1, 3, 5, 7 obtenus par les polynômes de Lagrange
!-!< -> transposée des coeff précédents, divisée par 2

real(kind=R8), parameter ::   cp0_00 = +1.0_R8

real(kind=R8), parameter ::   cp1_00 = +1.0_R8/4.0_R8
real(kind=R8), parameter ::   cp1_01 = +2.0_R8/4.0_R8
real(kind=R8), parameter ::   cp1_02 = +1.0_R8/4.0_R8

real(kind=R8), parameter ::   cp3_00 = - 1.0_R8/32.0_R8
real(kind=R8), parameter ::   cp3_01 = + 0.0_R8/32.0_R8
real(kind=R8), parameter ::   cp3_02 = + 9.0_R8/32.0_R8
real(kind=R8), parameter ::   cp3_03 = +16.0_R8/32.0_R8
real(kind=R8), parameter ::   cp3_04 = + 9.0_R8/32.0_R8
real(kind=R8), parameter ::   cp3_05 = + 0.0_R8/32.0_R8
real(kind=R8), parameter ::   cp3_06 = - 1.0_R8/32.0_R8

real(kind=R8), parameter ::   cp5_00 = +  3.0_R8/512.0_R8
real(kind=R8), parameter ::   cp5_01 = +  0.0_R8/512.0_R8
real(kind=R8), parameter ::   cp5_02 = - 25.0_R8/512.0_R8
real(kind=R8), parameter ::   cp5_03 = +  0.0_R8/512.0_R8
real(kind=R8), parameter ::   cp5_04 = +150.0_R8/512.0_R8
real(kind=R8), parameter ::   cp5_05 = +256.0_R8/512.0_R8
real(kind=R8), parameter ::   cp5_06 = +150.0_R8/512.0_R8
real(kind=R8), parameter ::   cp5_07 = +  0.0_R8/512.0_R8
real(kind=R8), parameter ::   cp5_08 = - 25.0_R8/512.0_R8
real(kind=R8), parameter ::   cp5_09 = +  0.0_R8/512.0_R8
real(kind=R8), parameter ::   cp5_10 = +  3.0_R8/512.0_R8

real(kind=R8), parameter ::   cp7_00 = -   5.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_01 = +   0.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_02 = +  49.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_03 = +   0.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_04 = - 245.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_05 = +   0.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_06 = +1225.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_07 = +2048.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_08 = +1225.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_09 = +   0.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_10 = - 245.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_11 = +   0.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_12 = +  49.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_13 = +   0.0_R8/4096.0_R8
real(kind=R8), parameter ::   cp7_14 = -   5.0_R8/4096.0_R8

type tborne
   integer(kind=I4) :: lb1 !! lower bound 1
   integer(kind=I4) :: ub1 !! upper bound 1
   integer(kind=I4) :: lb2 !! lower bound 2
   integer(kind=I4) :: ub2 !! upper bound 2
endtype tborne

public :: interp1D, restrict1D, interp2D, restrict2D, tborne

contains

!-!< -------------------------------------------------------------   \n
!-!< Fonction interpolant des points régulièrement espacés 1D        \n
!-!< -------------------------------------------------------------
function interp(tab, lb, ind, ordre)
implicit none
real(kind=R8)                                :: interp    !-!< valeur particulière interpolée
integer(kind=I4), intent(in)                 :: lb        !-!< borne inférieure
integer(kind=4),  intent(in)                 :: ind       !-!< position de l'élément "milieu"
integer(kind=4),  intent(in)                 :: ordre     !-!< ordre de l'interp 1, 3, 5 ou 7
real(kind=R8),    intent(in), dimension(lb:) :: tab       !-!< tableau 1D à interpoler

   select case (ordre)

      case(1)
         interp = ci1_00*tab(ind  ) +  &  !
                  ci1_01*tab(ind+1)       !

      case(3)
         interp = ci3_00*tab(ind-1) +  &  !
                  ci3_01*tab(ind  ) +  &  !
                  ci3_02*tab(ind+1) +  &  !
                  ci3_03*tab(ind+2)       !

      case(5)
         interp = ci5_00*tab(ind-2) +  &  !
                  ci5_01*tab(ind-1) +  &  !
                  ci5_02*tab(ind  ) +  &  !
                  ci5_03*tab(ind+1) +  &  !
                  ci5_04*tab(ind+2) +  &  !
                  ci5_05*tab(ind+3)       !

      case(7)
         interp = ci7_00*tab(ind-3) +  &  !
                  ci7_01*tab(ind-2) +  &  !
                  ci7_02*tab(ind-1) +  &  !
                  ci7_03*tab(ind  ) +  &  !
                  ci7_04*tab(ind+1) +  &  !
                  ci7_05*tab(ind+2) +  &  !
                  ci7_06*tab(ind+3) +  &  !
                  ci7_07*tab(ind+4)       !

      case default
         stop 'Bad choice in function "interp"'

   endselect

return
endfunction interp


!-!< ---------------------------------------------------------------------------------    \n
!-!< Interpolation des points régulièrement espacés 1D avec prise en compte des bords     \n
!-!< ---------------------------------------------------------------------------------
subroutine interp1D(tabgros, lb_gros, tabfin, lb_fin, ub_gros, ordre)
implicit none
integer(kind=I4), intent(in)                       :: lb_gros  !-!< indice inférieur
integer(kind=I4), intent(in)                       :: lb_fin   !-!< indice inférieur de tab_fin
integer(kind=I4), intent(in)                       :: ub_gros  !-!< taille de tabgros
integer(kind=I4), intent(in)                       :: ordre    !-!< ordre de l'interpolation
real(kind=R8),    intent(in), dimension(lb_gros:)  :: tabgros  !-!< tableau grossier à interpoler
real(kind=R8),   intent(out), dimension(lb_fin :)  :: tabfin   !-!< tableau résultant, 2 fois plus fin

   integer(kind=I4) :: l_inf, l_sup, i, ii
   real(kind=R8)    :: tmp0, dtmp

   real(kind=R8), dimension(       -ordre/2:        ordre  ) :: tab_inf
   real(kind=R8), dimension(ub_gros-ordre  :ub_gros+ordre/2) :: tab_sup

   !-!< bornes pour déterminer les limites d'utilisation de la fonction interp
   l_inf = ordre/2
   l_sup = ub_gros -l_inf

   !-!< extension du tableau par prolongement de la dérivée
   tab_inf(0:ordre) = tabgros(0:ordre)
   tmp0 = tab_inf(0)
   dtmp = tab_inf(0)-tab_inf(1)

   do i = 1, l_inf
      tab_inf(-i) = tmp0 + i*dtmp
   enddo

   do ii = 0, l_inf
      i = 2*ii
      tabfin(i  ) = tabgros(ii)
      tabfin(i+1) = interp( tab   = tab_inf,    &  !
                            lb    = -ordre/2,   &  !
                            ind   = ii,         &  !
                            ordre = ordre )        !
   enddo

   !-!< utilisation d'interp dans les limites normales
   do ii = l_inf+1, l_sup-1
      i = 2*ii
      tabfin(i  ) = tabgros(ii)
      tabfin(i+1) = interp( tab   = tabgros,    &  !
                            lb    = lb_gros,    &  !
                            ind   = ii,         &  !
                            ordre = ordre )        !
   enddo

   !-!< extension du tableau par prolongement de la dérivée
   tab_sup(ub_gros-ordre:ub_gros) = tabgros(ub_gros-ordre:ub_gros)
   tmp0 = tab_sup(ub_gros)
   dtmp = tab_sup(ub_gros)-tab_sup(ub_gros-1)

   do i = 1, l_inf
      tab_sup(i+ub_gros) = tmp0 + i*dtmp
   enddo

   do ii = l_sup, ub_gros-1
      i = 2*ii
      tabfin(i  ) = tabgros(ii)
      tabfin(i+1) = interp( tab   = tab_sup,          &  !
                            lb    = ub_gros - ordre,  &  !
                            ind   = ii,               &  !
                            ordre = ordre )              !
   enddo
   tabfin(2*ub_gros) = tabgros(ub_gros)

return
endsubroutine interp1D


!-!<----------------------------------------------------------------------------       \n
!-!< Routine d'interpolation sur grille 2D utilisant "interp1D" suivant les lignes
!-!< et les colonnes.                                                                  \n
!-!<----------------------------------------------------------------------------
subroutine interp2D(tabgro, bgro, tabfin, bfin, ordre)
implicit none
type(tborne),     intent(in )                                                  :: bfin, bgro  !-!< indices des tableaux
integer(kind=I4), intent(in )                                                  :: ordre       !-!< ordre de l'interpolation
real(kind=R8),    intent(in ), dimension(bgro%lb1:bgro%ub1, bgro%lb2:bgro%ub2) :: tabgro      !-!< tableau grossier départ
real(kind=R8),    intent(out), dimension(bfin%lb1:bfin%ub1, bfin%lb2:bfin%ub2) :: tabfin      !-!< tableau résultant fin

   integer(kind=I4) :: ii, j

   real(kind=R8), dimension(bgro%lb1:bgro%ub1, bfin%lb2:2*bgro%ub2) :: tab_tmp

   do ii = bgro%lb1, bgro%ub1

      call interp1D( tabgros = tabgro(ii,:),    &  !
                     lb_gros = bgro%lb2,        &  !
                     tabfin  = tab_tmp(ii,:),   &  !
                     lb_fin  = bfin%lb2,        &  !
                     ub_gros = bgro%ub2,        &  !
                     ordre   = ordre)              !
   enddo

   do j = bfin%lb2, 2*bgro%ub2

      call interp1D( tabgros = tab_tmp(:,j),    &  !
                     lb_gros = bgro%lb1,        &  !
                     tabfin  = tabfin(:,j),     &  !
                     lb_fin  = bfin%lb1,        &  !
                     ub_gros = bgro%ub1,        &  !
                     ordre   = ordre )             !
   enddo

return
endsubroutine interp2D


!-!<------------------------------------------------------------------------------------------------------------
!-!<                          RESTRICTION
!-!<------------------------------------------------------------------------------------------------------------
!-!< ----------------------------------------------------------------------------      \n
!-!< Fonction retournant la valeur pondérée en UN endroit donné d'UN tableau 1D        \n
!-!< ----------------------------------------------------------------------------
function restrict(tab, lb, ind, ordre)
implicit none
real(kind=R8)                               :: restrict  !-!< valeur particulière pondérée
integer(kind=4), intent(in)                 :: lb        !-!< borne inférieure
integer(kind=4), intent(in)                 :: ind       !-!< position de l'élément "milieu"
integer(kind=4), intent(in)                 :: ordre     !-!< ordre de la restriction 1, 3, 5 ou 7
real(kind=R8),   intent(in), dimension(lb:) :: tab       !-!< tableau 1D à réduire

   select case (ordre)

      case(0)
         restrict =  cp0_00*tab(ind)

      case(1)
         restrict =  cp1_00*tab(ind-1) +  &  !
                     cp1_01*tab(ind  ) +  &  !
                     cp1_02*tab(ind+1)       !

      case(3)
         restrict =  cp3_00*tab(ind-3) +  &  !
                     cp3_01*tab(ind-2) +  &  !
                     cp3_02*tab(ind-1) +  &  !
                     cp3_03*tab(ind  ) +  &  !
                     cp3_04*tab(ind+1) +  &  !
                     cp3_05*tab(ind+2) +  &  !
                     cp3_06*tab(ind+3)       !

      case(5)
         restrict =  cp5_00*tab(ind-5) +  &  !
                     cp5_01*tab(ind-4) +  &  !
                     cp5_02*tab(ind-3) +  &  !
                     cp5_03*tab(ind-2) +  &  !
                     cp5_04*tab(ind-1) +  &  !
                     cp5_05*tab(ind  ) +  &  !
                     cp5_06*tab(ind+1) +  &  !
                     cp5_07*tab(ind+2) +  &  !
                     cp5_08*tab(ind+3) +  &  !
                     cp5_09*tab(ind+4) +  &  !
                     cp5_10*tab(ind+5)       !

      case(7)
         restrict =  cp7_00*tab(ind-7) +  &  !
                     cp7_01*tab(ind-6) +  &  !
                     cp7_02*tab(ind-5) +  &  !
                     cp7_03*tab(ind-4) +  &  !
                     cp7_04*tab(ind-3) +  &  !
                     cp7_05*tab(ind-2) +  &  !
                     cp7_06*tab(ind-1) +  &  !
                     cp7_07*tab(ind  ) +  &  !
                     cp7_08*tab(ind+1) +  &  !
                     cp7_09*tab(ind+2) +  &  !
                     cp7_10*tab(ind+3) +  &  !
                     cp7_11*tab(ind+4) +  &  !
                     cp7_12*tab(ind+5) +  &  !
                     cp7_13*tab(ind+6) +  &  !
                     cp7_14*tab(ind+7)       !

      case default
         stop 'Bad choice in function "restrict"'

   endselect

return
endfunction restrict


!-!<----------------------------------------------------------------------------    \n
!-!< subroutine pondérant UN tableau 1D entier en prenant les bords en compte       \n
!-!<----------------------------------------------------------------------------
subroutine restrict1D(tabfin, lb_fin, tabgros, lb_gros, ub_gros, ordre)
implicit none
integer(kind=I4), intent(in )                      :: lb_fin      !-!< indice inférieur de tab_fin
integer(kind=I4), intent(in )                      :: lb_gros     !-!< indice inférieur
integer(kind=I4), intent(in )                      :: ub_gros     !-!< taille de tabgros
integer(kind=I4), intent(in )                      :: ordre       !-!< ordre de la restriction
real(kind=R8),    intent(in ), dimension(lb_fin:)  :: tabfin      !-!< tableau de départ
real(kind=R8),    intent(out), dimension(lb_gros:) :: tabgros     !-!< tableau grossier résultant

   integer(kind=I4) :: l_inf, l_sup, i, ii
   real(kind=R8)    :: tmp0, dtmp

   real(kind=R8), dimension(          -   ordre:          2*ordre) :: tab_inf
   real(kind=R8), dimension(2*ub_gros - 2*ordre:2*ub_gros + ordre) :: tab_sup

   !-!< bornes pour déterminer les limites d'utilisation de la fonction restrict
   l_inf = ordre/2
   l_sup = ub_gros -l_inf

   !-!< extension du tableau par prolongement de la dérivée
   tab_inf(0:2*ordre) = tabfin(0:2*ordre)
   tmp0 = tab_inf(0)
   dtmp = tab_inf(0)-tab_inf(1)

   do i = 1, ordre
      tab_inf(-i) = tmp0 + i*dtmp
   enddo

   do ii = 0, l_inf
      i = 2*ii
      tabgros(ii) = restrict( tab   = tab_inf,  &  !
                              lb    = -ordre,   &  !
                              ind   = i,        &  !
                              ordre = ordre )      !
   enddo

   !-!< utilisation d'interp dans les limites normales
   do ii = l_inf+1, l_sup-1
      i = 2*ii
      tabgros(ii) = restrict( tab   = tabfin,   &  !
                              lb    = lb_fin,   &  !
                              ind   = i,        &  !
                              ordre = ordre )      !
   enddo

   !-!< extension du tableau par prolongement de la dérivée
   tab_sup(2*ub_gros-2*ordre:2*ub_gros) = tabfin(2*ub_gros-2*ordre:2*ub_gros)
   tmp0 = tab_sup(2*ub_gros)
   dtmp = tab_sup(2*ub_gros)-tab_sup(2*ub_gros-1)

   do i = 1, ordre
      tab_sup(2*ub_gros+i) = tmp0 + i*dtmp
   enddo

   do ii = l_sup, ub_gros
      i = 2*ii
      tabgros(ii) = restrict( tab   = tab_sup,              &  !
                              lb    = 2*ub_gros - 2*ordre,  &  !
                              ind   = i,                    &  !
                              ordre = ordre )                  !
   enddo

return
endsubroutine restrict1D


!-!<----------------------------------------------------------------------------       \n
!-!< Routine de restriction sur grille 2D utilisant "restrict1D" suivant les lignes
!-!< et les colonnes.                                                                  \n
!-!<----------------------------------------------------------------------------
subroutine restrict2D(tabfin, bfin, tabgros, bgros, ordre)
implicit none
type(tborne),     intent(in )                                                      :: bfin, bgros     !-!< indices des tableaux
integer(kind=I4), intent(in )                                                      :: ordre           !-!< ordre de l'interpolation
real(kind=R8),    intent(in ), dimension( bfin%lb1: bfin%ub1,  bfin%lb2: bfin%ub2) :: tabfin          !-!< tableau de départ fin
real(kind=R8),    intent(out), dimension(bgros%lb1:bgros%ub1, bgros%lb2:bgros%ub2) :: tabgros         !-!< tableau grossier résultant

   integer(kind=I4) :: ii, j

   real(kind=R8), dimension(bgros%lb1:bgros%ub1, bfin%lb2:2*bgros%ub2) :: tab_tmp

   do j = bfin%lb2, 2*bgros%ub2

      call restrict1D( tabfin  = tabfin(:,j),   &  !
                       lb_fin  = bfin%lb1,      &  !
                       tabgros = tab_tmp(:,j),  &  !
                       lb_gros = bgros%lb1,     &  !
                       ub_gros = bgros%ub1,     &  !
                       ordre   = ordre )           !
   enddo

   do ii = bgros%lb1, bgros%ub1

      call restrict1D( tabfin  = tab_tmp(ii,:), &  !
                       lb_fin  = bfin%lb2,      &  !
                       tabgros = tabgros(ii,:), &  !
                       lb_gros = bgros%lb2,     &  !
                       ub_gros = bgros%ub2,     &  !
                       ordre   = ordre )           !
   enddo

return
endsubroutine restrict2D


!-!< --------------------------------------------------------------------     \n
!-!< subroutine générant les coefficients pour l'interpolation d'ordre k      \n
!-!< --------------------------------------------------------------------
subroutine genere_coeff_lagrange()
implicit none

   integer(kind=I4) :: i, j, k, n, c
   real(kind=R8)    :: coeff

   do

      write(*,*) 'n' ; read(*,*) n ; if (n==0) exit
      write(*,*) 'k' ; read(*,*) k
      write(*,*) 'c' ; read(*,*) c

      do i = 0, n

         coeff = 1.0d0
         do j = 0, n
            if (j==i) cycle
            coeff = coeff * ( 0.5_R8 * ( 2 * k - 2 * j + 1 ) )/( i - j )
         enddo
         write(*,*) coeff*c

      enddo

   enddo

return
endsubroutine genere_coeff_lagrange

endmodule intpl


