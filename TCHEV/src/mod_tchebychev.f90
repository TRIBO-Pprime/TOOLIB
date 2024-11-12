
!< author: Arthur Francisco
!  version: 1.0.0
!  date: february, 27 2023
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!     **Routines to subtract a least square polynomial from a surface**
!  </span>

module tchebychev
use data_arch,     only : I4, R8
use miscellaneous, only : get_unit
use least_squares, only : moindres_carres_lineaire
!$ use omp_lib

implicit none

private

public :: tcheby, tab_tcheby, tab_poly_tcheby, tab_Jf_tcheby, coeff_tcheby_vers_monome, coeff_tcheby_xy_vers_monome, & !
          coeff_poly_tcheby_xy_vers_poly_monome, least_squares_tcheby

contains

   !==================================================================================================
   !> @brief         Valeur en x du polynôme de Tchebichev de degré n
   !> @details       @f{align*}{
   !> P_n(x)=\sum_{i=1}^{n}\left[ 2^{i-1} t_i x^i +\sum_{k=1}^{i/2} \frac{i(i-k-1)!}{k!(i-2k)!} (-1)^k 2^{i-2k-1} t_i x^{i-2k} \right]
   !>                @f}
   !>                Site : http://fr.wikipedia.org/wiki/Polynômes_de_Tchebychev
   !==================================================================================================
   function tcheby(n, x)
   implicit none
   real(kind=R8) :: tcheby
   integer(kind=I4), intent(in) :: n !! degré du polynôme
   real(kind=R8), intent(in)    :: x !! variable

      integer(kind=I4) :: k, kk
      real(kind=R8)    :: y, tmp

      if (n==0) then
         tcheby = 1._R8
         return
      endif
      tmp = 0.
      do k = 1, n/2
         y = 1._R8
         do kk = 1, k-1
            y = y*(real((n-k -kk),kind = R8)/kk)
         enddo
         tmp = tmp + y*((-1)**k)*((2*x)**(n -2*k))/k
      enddo
      tcheby = ((2*x)**n)/2 +n*tmp/2.

   return
   endfunction tcheby
   !==================================================================================================
   !> @brief         Valeurs tabulées de polynômes de Tchebychev
   !> @details       On se donne des points xi le long d'UN axe et on calcule les valeurs d'UN      \n
   !>                polynôme de Tchebychev de degré j en ces points.
   !==================================================================================================
   subroutine tab_tcheby(deg, nx, vec_x, tab_tche)
   implicit none
   integer(kind=I4), intent(in ) :: deg !! degré max du polynôme
   integer(kind=I4), intent(in ) :: nx  !! nbre de points de calcul
   real(kind=R8),    intent(in ), dimension(1:nx) :: vec_x !! vecteur des points de calcul
   real(kind=R8),    intent(out), allocatable, dimension(:,:)   :: tab_tche !! tableau des valeurs calculées

      integer(kind=I4) :: i, j
      real(kind=R8)    :: xi

      allocate( tab_tche(1:nx, 1:deg+1) )

      do i = 1, nx
         xi = vec_x(i) ! points de l'axe
         do j = 0, deg ! pour UN degré de polynôme donné, valeurs en ces points
            tab_tche(i, j+1) = tcheby(j, xi)
         enddo
      enddo

   return
   endsubroutine tab_tcheby
   !==================================================================================================
   !> @brief         Surface définie par UN produit de polynômes de Tchebychev en x et y
   !> @details       Le domaine étant discrétisé et UN ensemble de coefficients donnés (provenant  \n
   !>                 d'une approximation par moindres carrés) on a la valeur de la fonction       \n
   !>                 surface en chaque point.
   !==================================================================================================
   subroutine tab_poly_tcheby(nx1, nx2, nvarx, nvary, nb_var, tab_tche1, tab_tche2, var, tab_poly_tche, multi_thread)
   implicit none
   integer(kind=I4), intent(in ) :: nb_var   !!
   integer(kind=I4), intent(in ) :: nx1      !! nbre de points de calcul selon x
   integer(kind=I4), intent(in ) :: nx2      !! nbre de points de calcul selon y
   integer(kind=I4), intent(in ) :: nvarx    !! degré max de Tchebychev utilisé selon x
   integer(kind=I4), intent(in ) :: nvary    !! degré max de Tchebychev utilisé selon y
   real(kind=R8),    intent(in ), dimension(1:nx1, 1:nvarx+1) :: tab_tche1       !! tableau des valeurs calculées
   real(kind=R8),    intent(in ), dimension(1:nx2, 1:nvary+1) :: tab_tche2       !! tableau des valeurs calculées
   real(kind=R8),    intent(in ), dimension(1:nb_var)         :: var             !! vecteur des coefficients
   real(kind=R8),    intent(out), dimension(1:nx1, 1:nx2)     :: tab_poly_tche   !! tableau résultant : surface
   logical(kind=I4), intent(in ), optional :: multi_thread

      real(kind=R8)    :: tmp1, tmp2
      integer(kind=I4) :: ivar, jvar, ij, ipt, jpt, ibatch, nb_threads
      logical(kind=I4) :: mlth

      mlth = .false.
      if ( present( multi_thread ) ) mlth = multi_thread

      ! surface d'UN produit de polynômes de tchebytchev

      nb_threads = omp_get_num_procs()
      ibatch = max( nx2/nb_threads, 1 )

      !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nb_threads) IF(mlth)
      !$OMP DO SCHEDULE (STATIC,ibatch) PRIVATE(ipt, tmp1, ij, jvar, tmp2, ivar)

      do jpt = 1, nx2
         do ipt = 1, nx1            ! Un dommaine discrétisé (nx+1)*(nx+1) est balayé

            tmp1 = 0._R8
            ij   = 0                ! chaque point est unique
            do jvar = 0, nvary      ! En chaque point, la fonction calculée est le produit des sommes
               tmp2 = 0._R8         !  de polynômes de Tchebychev. En effet l'approximation est faite
               do ivar = 0, nvarx   !  avec UN polynôme à variables séparées.
                  ij = ij +1
                  tmp2 = tmp2 +var(ij)*tab_tche1(ipt, ivar +1)
               enddo
               tmp1 = tmp1 +tmp2*tab_tche2(jpt, jvar +1)
            enddo
            tab_poly_tche(ipt, jpt) = tmp1

         enddo
      enddo

      !$OMP END DO
      !$OMP END PARALLEL


   return
   endsubroutine tab_poly_tcheby
   !==================================================================================================
   !> @brief         Tableau des dérivées par rapport aux coefficients de tab_tche
   !==================================================================================================
   subroutine tab_Jf_tcheby(nx1, nx2, nb_pts, nvarx, nvary, nb_var, tab_tche1, tab_tche2, tab_Jf, imask, multi_thread)
   implicit none
   integer(kind=I4), intent(in ) :: nb_var   !! nombre de fonctions de base utilisées
   integer(kind=I4), intent(in ) :: nx1      !! nbre de points de calcul selon x
   integer(kind=I4), intent(in ) :: nx2      !! nbre de points de calcul selon y
   integer(kind=I4), intent(in ) :: nb_pts   !!
   integer(kind=I4), intent(in ) :: nvarx    !! degré max de Tchebychev utilisé selon x
   integer(kind=I4), intent(in ) :: nvary    !! degré max de Tchebychev utilisé selon y
   real(kind=R8),    intent(in ), dimension(1:nx1, 1:nvarx+1) :: tab_tche1    !! tableau des valeurs calculées
   real(kind=R8),    intent(in ), dimension(1:nx2, 1:nvary+1) :: tab_tche2    !! tableau des valeurs calculées
   integer(kind=I4), intent(in ), dimension(1:nx1, 1:nx2), optional :: imask  !! masque
   real(kind=R8),    intent(out), allocatable, dimension(:,:) :: tab_Jf       !! tableau des dérivées
   logical(kind=I4), intent(in ), optional :: multi_thread

      integer(kind=I4) :: ivar, jvar, ij, ipt, jpt, ijpt, ibatch, nb_threads
      logical(kind=I4) :: lmask
      logical(kind=I4) :: mlth

      mlth = .false.
      if ( present( multi_thread ) ) mlth = multi_thread

      nb_threads = omp_get_num_procs()
      ibatch = max( nx2/nb_threads, 1 )

      lmask = present(imask)

      allocate( tab_Jf(1:nb_pts, 1:nb_var) )
      ! table des dérivées des poly par rapport aux coeff, c'est donc UN produit de polynômes

      !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nb_threads) IF(mlth)
      !$OMP DO SCHEDULE (STATIC,ibatch) PRIVATE(ijpt, ipt, jpt, ij, ivar, jvar)

      do jpt = 1, nx2
         do ipt = 1, nx1

            if ( lmask ) then
               if ( imask(ipt, jpt)==0 ) cycle
            endif

            ijpt = (jpt-1)*nx1 + ipt

            ij = 0
            do jvar = 0, nvary
               do ivar = 0, nvarx
                  ij = ij +1
                  tab_Jf(ijpt, ij) = tab_tche1(ipt, ivar +1)*tab_tche2(jpt, jvar +1)
               enddo
            enddo
         enddo
      enddo

      !$OMP END DO
      !$OMP END PARALLEL

   return
   endsubroutine tab_Jf_tcheby
   !==================================================================================================
   ! MODIFIE !!!
   !> @brief         Transformation des coefficients de Tchebychev en coefficients de monômes
   !> @details       Un monôme de degré p reçoit de Ti(x) (polynôme de Tcheby de degré i) :
   !>                @f{align*}{
   !> 2^{i-1} t_i \text{ si i=p}
   !> m(p)=m(p)+\frac{i(i-k-1)!}{k!(i-2k)!} (-1)^k 2^{i-2k-1} t(i) \text{ si } i=2k+p}
   !>                @f}
   !==================================================================================================
   subroutine coeff_tcheby_vers_monome(coeff_t, coeff_m, deg)
   implicit none
   integer(kind=I4), intent(in) :: deg !! degré du polynôme
   real(kind=R8), intent(in ), dimension(0:deg) :: coeff_t  !! coefficients de la CL de polynômes de Tchebychev
   real(kind=R8), intent(out), dimension(0:deg) :: coeff_m  !! coefficients de la CL de monômes

      integer(kind=I4) :: k, n, p, q
      real(kind=R8) :: tkm2q, tmp

      n = deg
      coeff_m(0:deg) = 0._R8

      select case (n)

         case (0)
            coeff_m(0)   =   coeff_t(0)
         case (1)
            coeff_m(0:1) = [ coeff_t(0)             , coeff_t(1) ]
         case (2)
            coeff_m(0:2) = [ coeff_t(0) - coeff_t(2), coeff_t(1)                 , 2 * coeff_t(2) ]
         case (3)
            coeff_m(0:3) = [ coeff_t(0) - coeff_t(2), coeff_t(1) - 3 * coeff_t(3), 2 * coeff_t(2), 4 * coeff_t(3) ]
         case (4:)
            coeff_m(0:3) = [ coeff_t(0) - coeff_t(2), coeff_t(1) - 3 * coeff_t(3), 2 * coeff_t(2), 4 * coeff_t(3) ]

            do k = 4, n

               q = 0
               coeff_m(k-2*q) = coeff_m(k-2*q) + coeff_t(k) * (2**(k-1))

               q = 1
               coeff_m(k-2*q) = coeff_m(k-2*q) + coeff_t(k) * (-k/2._R8) * (2**(k-2))

               do q = 2, k/2

                  tmp = 1._R8
                  do p = 1, q - 1
                     tmp = tmp * ( (k - q - p) / (1._R8*p) )
                  enddo

                  tkm2q = (k/2._R8) * ((-1)**q) * (2**(k-2*q)) * (1._R8/q) * tmp

                  coeff_m(k-2*q) = coeff_m(k-2*q) + coeff_t(k) * tkm2q

               enddo

            enddo
      endselect

   return
   endsubroutine coeff_tcheby_vers_monome
   !==================================================================================================
   !> @brief         Transformation du produit ti(x)*tj(y) de Tchebychev en coefficients de monômes x^i.y^j
   !==================================================================================================
   subroutine coeff_tcheby_xy_vers_monome(tab_coeff_m, deg_x, deg_y)
   implicit none
   integer(kind=I4), intent(in ) :: deg_x !! degré du polynôme en x
   integer(kind=I4), intent(in ) :: deg_y !! degré du polynôme en y
   real(kind=R8),    intent(out), dimension(0:deg_x, 0:deg_y) :: tab_coeff_m  !! coefficients de la CL de monômes x^i * y^j

      integer(kind=I4) :: i, j
      real(kind=R8), dimension(0:deg_x) :: coeff_tx, coeff_mx
      real(kind=R8), dimension(0:deg_y) :: coeff_ty, coeff_my

      ! single Tchebychev polynomial, deg_x
      coeff_tx(0:deg_x) = 0._R8
      coeff_tx(  deg_x) = 1._R8
      call coeff_tcheby_vers_monome(coeff_tx(0:deg_x), coeff_mx(0:deg_x), deg_x)

      ! single Tchebychev polynomial, deg_y
      coeff_ty(0:deg_y) = 0._R8
      coeff_ty(  deg_y) = 1._R8
      call coeff_tcheby_vers_monome(coeff_ty(0:deg_y), coeff_my(0:deg_y), deg_y)

      tab_coeff_m(0:deg_x, 0:deg_y) = 0._R8
      do j = 0, deg_y
      do i = 0, deg_x
         tab_coeff_m(i, j) = coeff_mx(i)*coeff_my(j)
      enddo
      enddo

   return
   endsubroutine coeff_tcheby_xy_vers_monome
   !==================================================================================================
   !> @brief         Transformation d'une CL de produits de polynômes de Tchebychev en x et y en polynôme classique
   !==================================================================================================
   subroutine coeff_poly_tcheby_xy_vers_poly_monome(var, coeff_m, deg_x, deg_y)
   implicit none
   integer(kind=I4), intent(in ) :: deg_x !! degré du polynôme en x
   integer(kind=I4), intent(in ) :: deg_y !! degré du polynôme en y
   real(kind=R8),    intent(in ), dimension(1:(deg_x+1)*(deg_y+1)) :: var     !! coefficients du produits de polynômes de Tchebychev
   real(kind=R8),    intent(out), dimension(1:(deg_x+1)*(deg_y+1)) :: coeff_m !! coefficients du polynôme classique en x et y

      integer(kind=I4) :: i, j, ij
      real(kind=R8), dimension(0:deg_x, 0:deg_y) :: tab_poly__m, tab_coeff_m

      coeff_m(1:(deg_x+1)*(deg_y+1)) = 0._R8
      tab_poly__m(0:deg_x, 0:deg_y)  = 0._R8
      ij = 0
      do j = 0, deg_y
      do i = 0, deg_x
         ij = ij +1
         call coeff_tcheby_xy_vers_monome(tab_coeff_m(0:i, 0:j), deg_x=i, deg_y=j)
         tab_coeff_m(0:i, 0:j) = var(ij)*tab_coeff_m(0:i, 0:j)
         tab_poly__m(0:i, 0:j) = tab_poly__m(0:i, 0:j) +tab_coeff_m(0:i, 0:j)
      enddo
      enddo

      ij = 0
      do j = 0, deg_y
      do i = 0, deg_x
         ij = ij +1
         coeff_m(ij) = tab_poly__m(i, j)
      enddo
      enddo

   return
   endsubroutine coeff_poly_tcheby_xy_vers_poly_monome
!~    !==================================================================================================
!~    !> @brief         Même principe que Tchebychev avec des monômes (plus simple, mais moindres     \n
!~    !> @return        polynome résultat : fonction en (xi, yi)
!~    !==================================================================================================
!~    real(kind=R8) function polynome(xi, yi, var, nb_var)
!~    implicit none
!~    integer(kind=I4), intent(in) :: nb_var    !! nombre de fonctions de base utilisées
!~    real(kind=R8),    intent(in) :: xi        !! abscisse d'un points
!~    real(kind=R8),    intent(in) :: yi        !! ordonnée d'un points
!~    real(kind=R8),    intent(in), dimension(1:nb_var) :: var    !! vecteur des coefficients

!~       real(kind=R8)    :: tmp1, tmp2
!~       integer(kind=I4) :: ivar, jvar, nvar, ij

!~       nvar = nint(sqrt(real(nb_var, kind = R8))) -1
!~       tmp1 = 0._R8
!~       ij   = 0
!~       do ivar = 0, nvar
!~          tmp2 = 0._R8
!~          do jvar = 0, nvar
!~             ij = ij +1
!~             tmp2 = tmp2 +var(ij)*(yi**jvar)
!~          enddo
!~          tmp1 = tmp1 +tmp2*(xi**ivar)
!~       enddo
!~       polynome = tmp1

!~    return
!~    endfunction polynome
!~    !==================================================================================================
!~    !> @brief         Même principe que Tchebychev
!~    !> @return        dpolynome polynôme dérivé
!~    !==================================================================================================
!~    real(kind=R8) function dpolynome(xi, yi, nb_var, ind)
!~    implicit none
!~    integer(kind=I4), intent(in) :: nb_var !! nombre de fonctions de base utilisées
!~    integer(kind=I4), intent(in) :: ind    !! indice par rapport auquel on dérive
!~    real(kind=R8),    intent(in) :: xi     !! abscisse d'un points
!~    real(kind=R8),    intent(in) :: yi     !! ordonnée d'un points

!~       integer(kind=I4) :: ivar, jvar, nvar, ij

!~       nvar = nint(sqrt(real(nb_var, kind = R8))) -1
!~       ij   = 0
!~       dpolynome = 0
!~       do ivar = 0, nvar
!~          do jvar = 0, nvar
!~             ij = ij +1
!~             if (ij==ind) then
!~                dpolynome = (xi**ivar)*(yi**jvar)
!~                return
!~             endif
!~          enddo
!~       enddo

!~    return
!~    endfunction dpolynome
!~    !==================================================================================================
!~    !> @brief         Vérification des approximations par moindres carrés avec des polynômes
!~    !==================================================================================================
!~    subroutine verif_pol(tab1, tab2, var, long, larg, nb_var)
!~    implicit none
!~    integer(kind=I4), intent(in   ) :: long      !! taille x de la grille
!~    integer(kind=I4), intent(in   ) :: larg      !! taille y de la grille
!~    integer(kind=I4), intent(in   ) :: nb_var    !! nombre de fonctions de base utilisées
!~    real(kind=R8),    intent(inout), dimension(1:long, 1:larg) :: tab1   !! surface dont on a fait une approximation
!~    real(kind=R8),    intent(  out), dimension(1:long, 1:larg) :: tab2   !! approximation par moindres carrés
!~    real(kind=R8),    intent(in   ), dimension(1:nb_var)       :: var    !! vecteur des coefficients

!~       real(kind=R8)    :: xi, yi, x
!~       integer(kind=I4) :: i, j

!~       do j = 1, larg
!~       do i = 1, long
!~          xi = real(i-1, kind=R8)/(long -1)
!~          yi = real(j-1, kind=R8)/(larg -1)
!~          x = polynome(xi, yi, var, nb_var)
!~          tab1(i, j) = tab1(i, j) -x
!~          tab2(i, j) = x
!~       enddo
!~       enddo

!~    return
!~    endsubroutine verif_pol
!~    !==================================================================================================
!~    !> @brief         Vérification des approximations par moindres carrés avec des polynômes de Tcheby
!~    !==================================================================================================
!~    subroutine verif_tch(tab1, tab2, long, larg, tab_poly_tche)
!~    implicit none
!~    integer(kind=I4), intent(in   ) :: long      !! taille x de la grille
!~    integer(kind=I4), intent(in   ) :: larg      !! taille y de la grille
!~    real(kind=R8),    intent(inout), dimension(1:long, 1:larg)     :: tab1           !! surface dont on a fait une approximation qui lui est soustraite
!~    real(kind=R8),    intent(  out), dimension(1:long, 1:larg)     :: tab2           !! approximation par moindres carrés
!~    real(kind=R8),    intent(in   ), dimension(1:long+1, 1:larg+1) :: tab_poly_tche  !! surface dont on a fait une approximation

!~       real(kind=R8)    :: x
!~       integer(kind=I4) :: i, j

!~       do j = 1, larg
!~       do i = 1, long
!~          x = tab_poly_tche(i, j)
!~          tab1(i, j) = tab1(i, j) -x
!~          tab2(i, j) = x
!~       enddo
!~       enddo

!~    return
!~    endsubroutine verif_tch
!~    !==================================================================================================


   subroutine least_squares_tcheby(tab_in, tab_out, long1, long2, nvarx, nvary, imask, verif, multi_thread)
   !================================================================================================
   !< @note Function that returns
   !
   !  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                                        :: long1    !! taille x
   integer(kind=I4), intent(in )                                        :: long2    !! taille y
   integer(kind=I4), intent(in )                                        :: nvarx    !! degré du polynôme en x
   integer(kind=I4), intent(in )                                        :: nvary    !! degré du polynôme en y
   real   (kind=R8), intent(in ), dimension(1:long1, 1:long2)           :: tab_in   !! surface dont on a fait une approximation qui lui est soustraite
   real   (kind=R8), intent(out), dimension(1:long1, 1:long2)           :: tab_out  !! tableau résultant : surface
   integer(kind=I4), intent(in ), dimension(1:long1, 1:long2), optional :: imask    !! masque
   logical(kind=I4), intent(in )                             , optional :: verif    !! dump
   logical(kind=I4), intent(in )                             , optional :: multi_thread

      real(kind=R8), allocatable, dimension(:)   :: var1, vec_x1, vec_x2, vec_hij
      real(kind=R8), allocatable, dimension(:,:) :: tab_tche1, tab_tche2, Jf

      real(kind=R8), dimension(1:(nvarx +1)*(nvary +1)) :: coeff_mxy

      integer  (kind=I4) :: i, iu, nbvar, nbpts
      character(len=128) :: string

      if (nvarx==0 .and. nvary==0) then
         tab_out(1:long1, 1:long2) = 0._R8
         return
      endif

      nbvar = (nvarx +1)*(nvary +1)
      if (present(imask)) then
         nbpts = sum(imask)
      else
         nbpts = long1 * long2
      endif

      allocate( var1(1:nbvar) )
      allocate( vec_x1(1:long1), vec_x2(1:long2) )

         var1 = 0
         do i = 1, long1
            vec_x1(i) = -1. + (i - 1) * 2. / (long1 - 1)!real(i-long1/2, kind=R8)/(long1/2)
         enddo
         do i = 1, long2
            vec_x2(i) = -1. + (i - 1) * 2. / (long2 - 1)!real(i-long2/2, kind=R8)/(long2/2)
         enddo

         call tab_tcheby(deg = nvarx, nx = long1, vec_x = vec_x1(1:long1), tab_tche = tab_tche1)
         call tab_tcheby(deg = nvary, nx = long2, vec_x = vec_x2(1:long2), tab_tche = tab_tche2)

         call tab_Jf_tcheby(nx1 = long1,                          &  !
                            nx2 = long2,                          &  !
                         nb_pts = nbpts,                          &  !
                          nvarx = nvarx,                          &  !
                          nvary = nvary,                          &  !
                         nb_var = nbvar,                          &  !
                      tab_tche1 = tab_tche1(1:long1, 1:nvarx+1),  &  !
                      tab_tche2 = tab_tche2(1:long2, 1:nvary+1),  &  !
                         tab_Jf = Jf,                             &  !
                          imask = imask(1:long1, 1:long2),        &  !
                   multi_thread = multi_thread)                      !

         allocate(vec_hij(1:nbpts))

         if (present(imask)) then
            vec_hij(1:nbpts) = pack(tab_in(1:long1, 1:long2), mask = (imask(1:long1, 1:long2)==1))
         else
            vec_hij(1:nbpts) = reshape(tab_in(1:long1, 1:long2), shape = [nbpts])
         endif

         call moindres_carres_lineaire(nb_var = nbvar,                  &  !
                                       nb_pts = nbpts,                  &  !
                                          hij = vec_hij(1:nbpts),       &  !
                                         beta = var1(1:nbvar),          &  !
                                           Jf = Jf(1:nbpts, 1:nbvar))      !
         deallocate(vec_hij)


         if (present(verif) .and. verif) then

            call coeff_poly_tcheby_xy_vers_poly_monome(var = var1(1:nbvar),         &  !
                                                   coeff_m = coeff_mxy(1:nbvar),    &  !
                                                     deg_x = nvarx,                 &  !
                                                     deg_y = nvary)                    !

            write(string, '(a,i3.3,a)') '(', nbvar,'E18.8,a)'

            call get_unit(iu)
            open(iu, file = 'verif_tcheby_vers_monome.txt')
               write(iu, trim(string)) (coeff_mxy(i), i=1, nbvar), " coeff tchebychev"

               ! à ce stade on a les coefficients var1(:) tels que :
               ! surface_approchante(x,y) = var1(1).t0(x)t0(y) +var1(2).t1(x)t0(y) +var1(3).t2(x)t0(y) +...+var1(nvarx+1 +1).t0(x)t1(y) +...+var1(nbvar).t_{nvarx+1}(x)t_{nvary+1}(y)
               write(iu,'(a)') 'pour vérifier avec gwyddion, il faut préalablement décocher "surface carrée", donner à l''image un facteur d''échelle x et y identiques, '// &  !
                               'et mettre un décalage pour que x soit compris entre -1 et 1 tout comme y. '//                                                               &  !
                               'coeff_mxy donne les coeffs de 1 x x^2 x^3 ... x^i y xy x^2.y ...x^i.y^j'
            close(iu)
         endif

         call tab_poly_tcheby(nx1 = long1,                           &  !
                              nx2 = long2,                           &  !
                            nvarx = nvarx,                           &  !
                            nvary = nvary,                           &  !
                           nb_var = nbvar,                           &  !
                        tab_tche1 = tab_tche1(1:long1, 1:nvarx+1),   &  !
                        tab_tche2 = tab_tche2(1:long2, 1:nvary+1),   &  !
                              var = var1(1:nbvar),                   &  !
                    tab_poly_tche = tab_out(1:long1, 1:long2),       &  !
                     multi_thread = multi_thread)                       !

      deallocate( Jf    )
      deallocate( vec_x1, vec_x2 )
      deallocate( var1  )
      deallocate( tab_tche1, tab_tche2 )

   return
   endsubroutine least_squares_tcheby



endmodule tchebychev
