!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: november, 01 2024
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<     **PIKAIA oop example of use**
!<  </span>
!<
!< @note
!<  The example chosen is the generation of a surface with given statistical moments.
!<  The implementation is part of the analytical developments presented in:<br/>
!<  A. Francisco and N. Brunetière, "A hybrid method for fast and efficient rough surface generation",
!<  Proc IMechE Part J:J Engineering Tribology, 2016, Vol. 230(7) 747–768, DOI: 10.1177/1350650115612116
!<
!<  The surface heights are generated thanks the *tangent* function. This function is indeed very near the common material
!<  curves. Thus, providing the right lower and upper bounds makes it possible to match the desired statistical moments.
!< @endnote
!<
!< @note
!<  The main difficulty to cope with is the high non linearity and the convergence only possible near the solution.
!<  Classical optimization methods fail in reaching the true solution: the lower and upper limits ```x(1:2```.<br/>
!<  The **PIKAIA** program, which implements a genetic algorithm, is utilized to found a "good" solution, then a
!<  very classical optimization routine is used to find the nearly exact solution.
!< @endnote
!<
!< @warning
!<  **PIKAIA** cost function must be chosen so as to be maximized.<br/>
!<  For clarity reasons, the optimization function must be minimized because it is the deviation to the solution.
!<  @warning
!<
!< @note
!<  First each thread runs the **PIKAIA** program to find a rather good starting point ```x```. The parameters
!<  (as defined by ```CTRL```) are the same. The only difference is the starting population. Thus, several concurrent
!<  populations evolve during ```CTRL(02)``` generations.<br/>
!<  Then, the optimization subroutine starts from the best population.
!< @endnote
!<
!< @note
!<  ```make all ```    <br/>
!<  ```./prg```        <br/>
!<  The surface is in ascii mode in the "/out" directory
!< @endnote

program test_algen
use omp_lib
use sort_arrays,   only : sort_array2
use data_arch,     only : I4, R8, PI_R8
use miscellaneous, only : get_unit
use pikaia_oop,    only : pikaia_class

implicit none

integer(kind=I4), parameter :: nn  = 512        !! *Surface side dimension*
real(kind=R8),    parameter :: SSK = -01.0_R8   !! *Imposed surface height skewness*
real(kind=R8),    parameter :: SKU = +20.0_R8   !! *Imposed surface height kurtosis*

integer(kind=I4) :: status                      !! **PIKAIA** *status*
integer(kind=I4) :: i, j, k                     !! *loop indices*
integer(kind=I4) :: iu                          !! *i/o unit*
integer(kind=I4) :: nb_th                       !! *nb threads*

real(kind=R8), dimension(1:2) :: xx             !! *chromosom for* **PIKAIA**
real(kind=R8), dimension(1:2) :: xl, xu         !! *lower and upper bonds of xx*

real(kind=R8)    :: f_xx                        !! *cost function or optimization function at* ```xx```
real(kind=R8)    :: tab_sur(1:nn*nn)            !! *output surface heights*

type moment_stat
!! <span style="color:green">Statistical moment type</span>
   real(kind=R8) :: mu !! *mean*
   real(kind=R8) :: va !! *variance*
   real(kind=R8) :: si !! *standard deviation*
   real(kind=R8) :: sk !! *skewness*
   real(kind=R8) :: ku !! *kurtosis*
endtype moment_stat

type(moment_stat) :: mom   !! *a statistical moment variable*

type(pikaia_class)   :: p              !! **PIKAIA** *class instanciation*
real(kind=R8)        :: f              !! *best cost*
real(kind=R8)        :: tstart,tend    !! *time variables*


   !$ real(kind=R8)    :: ostart,oend    !! *openmp time variables*
   !$ integer(kind=I4) :: tid            !! *active thread*

   nb_th = 1

   ! get the number of available threads
   !$OMP PARALLEL PRIVATE(nb_th, TID)
   !$
   !$ tid = omp_get_thread_num()
   !$
   !$ if (tid == 0) then
   !$    nb_th = omp_get_num_threads()
   !$    write(*,'(A)') '--------------'
   !$    write(*,'(A,1X,I5)') 'number of OMP threads: ', nb_th
   !$    write(*,'(A)') '--------------'
   !$ endif
   !$OMP END PARALLEL

   xx(1:2) = 0.0_R8
   xl(1:2) = 0.0_R8
   xu(1:2) = 1.0_R8

   !initialize the class:
   call p%init(           n = 2,                   &  ! IN           ; the parameter space dimension, i.e., the number of adjustable parameters (size of the x vector).
                         xl = xl,                  &  ! IN, DIM(n)   ;  vector of lower bounds for x
                         xu = xu,                  &  ! IN, DIM(n)   ;  vector of upper bounds for x
                          f = cost,                &  !              ; user-supplied scalar function of n variables, which must have the pikaia_func procedure interface.
                     status = status,              &  ! OUT          ; status output flag (0 if there were no errors)
                    !iter_f = report_iteration,    &  !     OPT      ; user-supplied subroutine that will report the best solution for each generation. It must have the iter_func procedure interface.
                         np = 100,                 &  ! IN, OPT      ; number of individuals in a population (default is 100)
                       ngen = 1000,                &  ! IN, OPT      ; maximum number of iterations
                         nd = 9,                   &  ! IN           ; number of significant digits (i.e., number of genes) retained in chromosomal encoding
                     pcross = 0.85_R8,             &  ! IN, OPT      ; crossover probability; must be <= 1.0 (default is 0.85). If crossover takes place, either one or two splicing points are used, with equal probabilities
                     pmutmn = 0.0005_R8,           &  ! IN, OPT      ; minimum mutation rate; must be >= 0.0 (default is 0.0005)
                     pmutmx = 0.25_R8,             &  ! IN, OPT      ; maximum mutation rate; must be <= 1.0 (default is 0.25)
                       pmut = 0.005_R8,            &  ! IN, OPT      ; initial mutation rate; should be small (default is 0.005) (Note: the mutation rate is the probability that any one gene locus will mutate in any one generation.)
                       imut = 2,                   &  ! IN, OPT      ; mutation mode; 1/2/3/4/5 (default is 2).
                                                      !              1=one-point mutation, fixed rate.
                                                      !              2=one-point, adjustable rate based on fitness.
                                                      !              3=one-point, adjustable rate based on distance.
                                                      !              4=one-point+creep, fixed rate.
                                                      !              5=one-point+creep, adjustable rate based on fitness.
                                                      !              6=one-point+creep, adjustable rate based on distance.
                       fdif = 1._R8,               &  ! IN, OPT      ; relative fitness differential; range from 0 (none) to 1 (maximum). (default is 1.0)
                       irep = 3,                   &  ! IN, OPT      ; reproduction plan; 1/2/3=Full generational replacement/Steady-state-replace-random/Steady- state-replace-worst (default is 3)
                     ielite = 0,                   &  ! IN, OPT      ; elitism flag; 0/1=off/on (default is 0) (Applies only to reproduction plans 1 and 2)
                       ivrb = 0,                   &  ! IN, OPT      ; printed output 0/1/2=None/Minimal/Verbose
            convergence_tol = 1.0e-6_R8,           &  ! IN, OPT      ; convergence tolerance; must be > 0.0 (default is 0.0001)
         convergence_window = 200,                 &  ! IN, OPT      ; convergence window; must be >= 0 This is the number of consecutive solutions within the tolerance for convergence to be declared (default is 20)
         initial_guess_frac = 0.1_R8,              &  ! IN, OPT      ; raction of the initial population to set equal to the initial guess. Range from 0 (none) to 1.0 (all). (default is 0.1 or 10%).
                      iseed = 999)                    ! IN, OPT      ; random seed value; must be > 0 (default is 999)

   !Now call pikaia:
   call cpu_time(tstart)
   !$ ostart = omp_get_wtime()

   call p%solve(      x = xx(1:2),     &  ! INOUT, DIM(*) ;
                      f = f,           &  !   OUT         ;
                 status = status) !,   &  !   OUT         ;
!~                     omp = .false. )    ! IN, OPTIONAL

   !$ oend = omp_get_wtime()
   call cpu_time(tend)

   write(*, *) 'Total time spent: ', tend - tstart
   write(*, *) 'Real time spent:  ', oend - ostart

   write(*,*) 'Absolute diff after GA: ', abs_diff_sk_ku(chrom = xx(1:2))

   ! optimization function ran with this "good" solution
   call newton_raphson_downhill( x     = xx(1:2),           & ! starting/ending point
                                 fvec  = f_xx,              & ! best deviation to the wanted solution
                                 eps   = 1.e-6_R8,          & ! dx
                                 ndir  = 32*nb_th,          & ! number of directions to explore
                                 rel   = 0.9_R8)              ! relaxation parameter

   write(*,*) 'Absolute diff after refinement: ', abs_diff_sk_ku(chrom = xx(1:2))


   write(*,*) 'A (nn x nn) SSK, SKU random surface will be now generated'

   ! generation of nn*nn heights with the tangent function with xx(1:2) as bound parameters.
   ! mu = 0.
   ! va = si = 1.
   ! sk = SSK
   ! ku = SKU
   call profil_theo_trie_1D(tab = tab_sur(1:nn*nn),   &  ! resulting vector of heights
                             lg = nn*nn,              &  ! vector length
                              x = xx(1:2),            &  ! best variable couple
                             mx = mom)                   ! resulting moment

   ! height shuffle
   call melange(tab = tab_sur(1:nn*nn),   &  !
                 lg = nn*nn)                 !

   ! random surface output in ascii format
   call get_unit(iunit = iu)
   open(unit = iu, file = 'out/surf.dat')
      k = 1
      do i = 1, nn
      do j = 1, nn
         write(iu, *) i, j, tab_sur(k)
         k = k +1
      enddo
      enddo
   close(unit=iu)

   ! statistical moments check
   write(*,*) 'statistical moments:'
   write(*,*)  mom%mu, mom%si, mom%sk, mom%ku

   ! precision
   write(*,'(a, e8.2)') 'maximum difference (%) ', max( 100*abs((mom%sk -SSK)/SSK), & !
                                                        100*abs((mom%ku -SKU)/SKU) )

stop

contains

   subroutine cost(me, x, f)
   implicit none
   class(pikaia_class), intent(inout)               :: me
   real(kind=R8)      , intent(in   ), dimension(:) :: x
   real(kind=R8)      , intent(  out)               :: f

      f = 1./(1. + abs_diff_sk_ku(chrom = x(1:2)))

   return
   endsubroutine cost


   subroutine newton_raphson_downhill(x, fvec, eps, ndir, rel)
   implicit none
   integer(kind=I4), intent(in   )                 :: ndir
   real(kind=R8),    intent(in   )                 :: eps
   real(kind=R8),    intent(in   )                 :: rel
   real(kind=R8),    intent(inout), dimension(1:2) :: x
   real(kind=R8),    intent(out  )                 :: fvec

      real(kind=R8), dimension(1:2)    :: x0
      real(kind=R8), dimension(1:ndir) :: ct, st, ff, x1, x2

      real(kind=R8)    :: f0, fd1, fd2, ti, dfdx1, dfdx2, dx1, dx2, dt
      integer(kind=I4) :: i

      do i = 1, ndir
         ti = PI_R8*( -1. +(i -1)*2./ndir )
         ct(i) = cos(ti)
         st(i) = sin(ti)
      enddo
      x0 = x

      do
         f0  = abs_diff_sk_ku(chrom = (/x0(1)     , x0(2)     /))
         fd1 = abs_diff_sk_ku(chrom = (/x0(1) +eps, x0(2)     /))
         fd2 = abs_diff_sk_ku(chrom = (/x0(1)     , x0(2) +eps/))

         dfdx1 = (fd1 -f0)/eps
         dfdx2 = (fd2 -f0)/eps

         !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nb_th)
         !$OMP DO PRIVATE(i, dt, dx1, dx2)
            do i = 1, ndir
               dt  = -f0/( eps*(dfdx1*ct(i) +dfdx2*st(i)) )
               dx1 = eps*ct(i)*dt
               dx2 = eps*st(i)*dt

               x1(i) = x0(1) +rel*dx1
               x2(i) = x0(2) +rel*dx2

               ff(i) = abs_diff_sk_ku(chrom = (/x1(i), x2(i)/))
            enddo
         !$OMP END DO
         !$OMP END PARALLEL

         i  = minloc(ff, 1)
         x0 = (/x1(i), x2(i)/)

         if (ff(i)<1.e-8) exit
         write(*,*) x0(1), x0(2), ff(i)
      enddo

      x    = x0
      fvec = ff(i)
   return
   endsubroutine newton_raphson_downhill


   real(kind=R8) function abs_diff_sk_ku(chrom)
   implicit none
   real(kind=R8), intent(in), dimension(1:2):: chrom

      real(kind=R8) :: sk, ku

      call sk_ku(xx = chrom(1:2), sk = sk, ku = ku)
      abs_diff_sk_ku = abs(sk -SSK) +abs(ku -SKU)

   return
   endfunction abs_diff_sk_ku


   real(kind=R8) function cost_func(chrom)
   implicit none
   real(kind=R8), intent(in), dimension(1:2) :: chrom

      cost_func = 1./(1. + abs_diff_sk_ku(chrom(1:2)))

   return
   endfunction cost_func


   subroutine sk_ku(xx, sk, ku)
   implicit none
   real(kind=R8), intent(in ), dimension(1:2) :: xx
   real(kind=R8), intent(out)                 :: sk
   real(kind=R8), intent(out)                 :: ku

      real(kind=R8)    :: xa, xb, mu, si, a, b, un
      real(kind=R8)    :: h, hh, b1, b2, alp, bet
      integer(kind=I4) :: ia, ib, deb, fin, npts, long, i, k

      real(kind=R8), dimension(1:2) :: x
      long  = nn

      do k = 1, 2
        x(k) = max( xx(k), 1.e-4_R8 )
      enddo

      ia = long
      ib = long

      npts = long*long

      deb = 1    +ia
      fin = npts -ib

      a = x(1)
      b = x(2)

      hh = (2._R8 -a -b)/(npts -1)
      h  = (pi_R8/2)*hh

      xa = a +ia*hh
      xb = b +ib*hh

      b1 = -pi_R8/2 *(1.-a)
      b2 = +pi_R8/2 *(1.-b)

      alp = -(b2-npts*b1)/(b2-b1)
      bet =      (npts-1)/(b2-b1)

      un = 1._R8
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      mu = log(1._R8/sin((pi_R8*xb)/2.)*sin((pi_R8*xa)/2.))
      mu = (un/h)*mu +add_tang(1, deb, fin, alp, bet, mu=0._R8, si=1._R8)
      do i = 1, ia-1
         mu = mu +tang(i*un, 1, alp, bet, mu=0._R8, si=1._R8)
      enddo
      do i = npts, npts -(ib -2), -1
         mu = mu +tang(i*un, 1, alp, bet, mu=0._R8, si=1._R8)
      enddo
      mu = mu/npts
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      si = (pi_R8*(-2 + xa + xb))/2. - (mu**2*pi_R8*(-2 + xa + xb))/2. + 1._R8/tan((pi_R8*xa)/2.) + 1._R8/tan((pi_R8*xb)/2.) - 2*mu*Log(1._R8/sin((pi_R8*xb)/2.)*Sin((pi_R8*xa)/2.))
      si = (un/h)*si +add_tang(2, deb, fin, alp, bet, mu, si=1._R8)
      do i = 1, ia-1
         si = si+ tang(i*un, 2, alp, bet, mu, si=1._R8)
      enddo
      do i = npts, npts -(ib -2), -1
         si = si +tang(i*un, 2, alp, bet, mu, si=1._R8)
      enddo
      si = si/npts
      si = sqrt(si)
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      sk = (6*mu*pi_R8 - 2*mu**3*pi_R8 - 3*mu*pi_R8*xa + mu**3*pi_R8*xa - 3*mu*pi_R8*xb + mu**3*pi_R8*xb - 6*mu*1._R8/tan((pi_R8*xa)/2.) - 6*mu*1._R8/tan((pi_R8*xb)/2.) - 1._R8/sin((pi_R8*xa)/2.)**2 + 1._R8/sin((pi_R8*xb)/2.)**2 - 2*Log(Sin((pi_R8*xa)/2.)) + 6*mu**2*Log(Sin((pi_R8*xa)/2.)) + 2*Log(Sin((pi_R8*xb)/2.)) - 6*mu**2*Log(Sin((pi_R8*xb)/2.)))/(2.*si**3)
      sk = (un/h)*sk +add_tang(3, deb, fin, alp, bet, mu, si)
      do i = 1, ia-1
         sk = sk +tang(i*un, 3, alp, bet, mu, si)
      enddo
      do i = npts, npts -(ib -2), -1
         sk = sk +tang(i*un, 3, alp, bet, mu, si)
      enddo
      sk = sk/npts
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      ku = (6*pi_R8 - 36*mu**2*pi_R8 + 6*mu**4*pi_R8 - 3*pi_R8*xa + 18*mu**2*pi_R8*xa - 3*mu**4*pi_R8*xa - 3*pi_R8*xb + 18*mu**2*pi_R8*xb - 3*mu**4*pi_R8*xb + 4*(-2 + 9*mu**2)*1._R8/tan((pi_R8*xa)/2.) + 4*(-2 + 9*mu**2)*1._R8/tan((pi_R8*xb)/2.) + 12*mu*1._R8/sin((pi_R8*xa)/2.)**2 - 12*mu*1._R8/sin((pi_R8*xb)/2.)**2 + 24*mu*Log(Sin((pi_R8*xa)/2.)) - 24*mu**3*Log(Sin((pi_R8*xa)/2.)) - 24*mu*Log(Sin((pi_R8*xb)/2.)) + 24*mu**3*Log(Sin((pi_R8*xb)/2.)) + 1._R8/sin((pi_R8*xa)/2.)**4*Sin(pi_R8*xa) + 1._R8/sin((pi_R8*xb)/2.)**4*Sin(pi_R8*xb))/(6.*si**4)
      ku = (un/h)*ku +add_tang(4, deb, fin, alp, bet, mu, si)
      do i = 1, ia-1
         ku = ku +tang(i*un, 4, alp, bet, mu, si)
      enddo
      do i = npts, npts -(ib -2), -1
         ku = ku +tang(i*un, 4, alp, bet, mu, si)
      enddo
      ku = ku/npts
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .

   return
   endsubroutine sk_ku


   real(kind=R8) function add_tang(n, deb, fin, alp, bet, mu, si)
   implicit none
   real(kind=R8),    intent(in) :: alp, bet, mu, si
   integer(kind=i4), intent(in) :: n, fin, deb

      real(kind=R8) :: xdeb, xfin

      xdeb = deb
      xfin = fin
      add_tang = (1./12)*( +9*(tang(xdeb +0.0_R8, n, alp, bet, mu, si)+tang(xfin -0.0_R8, n, alp, bet, mu, si)) &
                           +1*(tang(xdeb +1.0_R8, n, alp, bet, mu, si)+tang(xfin -1.0_R8, n, alp, bet, mu, si)) &
                           -4*(tang(xdeb +0.5_R8, n, alp, bet, mu, si)+tang(xfin -0.5_R8, n, alp, bet, mu, si)) )
   return
   endfunction add_tang


   real(kind=R8) function tang(xi, n, alp, bet, mu, si)
   implicit none
   real(kind=R8),    intent(in) :: xi, alp, bet, mu, si
   integer(kind=i4), intent(in) :: n

      real(kind=R8) :: tmp

      tmp = (xi +alp)/bet
      tang = ( (tan(tmp) -mu)/si )**n

   return
   endfunction tang


   subroutine calc_moments_1D(tab, mx, nb_mom, lg)
   implicit none
   integer(kind=I4) , intent(in )                  :: lg
   integer(kind=I4) , intent(in )                  :: nb_mom
   real(kind=R8)    , intent(in ), dimension(1:lg) :: tab
   type(moment_stat), intent(out)                  :: mx

      integer(kind=I4) :: i
      real(kind=R8)    :: tmp

      mx%mu = 0
      mx%si = 0
      mx%va = 0
      mx%Sk = 0
      mx%Ku = 0

      do i = 1, lg
         mx%mu = mx%mu +tab(i)/lg
      enddo
      if (nb_mom==1) return

      do i = 1, lg
         mx%va = mx%va +((tab(i) -mx%mu)**2)/lg
      enddo
      mx%si = sqrt( mx%va )
      if (nb_mom==2) return
      if (mx%si < 1.e-15_R8) then
         mx%Sk = 0
         mx%Ku = 0
      else
         do i = 1, lg
            tmp = (tab(i) -mx%mu)/mx%si
            mx%Sk = mx%Sk +(tmp**3)/lg
            mx%Ku = mx%Ku +(tmp**4)/lg
         enddo
      endif

   return
   endsubroutine calc_moments_1D


   subroutine profil_theo_trie_1D(tab, lg, x, mx)
   implicit none
   integer(kind=I4) , intent(in )                  :: lg
   real(kind=R8)    , intent(out), dimension(1:lg) :: tab
   real(kind=R8)    , intent(in ), dimension(1:2)  :: x
   type(moment_stat), intent(out)                  :: mx

      real(kind=R8)    :: b1, b2, alp, bet
      integer(kind=I4) :: i

      b1 = -PI_R8/2 *(1. -x(1))
      b2 = +PI_R8/2 *(1. -x(2))
      alp = -(b2- lg*b1)/(b2 -b1)
      bet =      (lg- 1)/(b2 -b1)
      do i = 1, lg
         tab(i) = tan( (i +alp)/bet )
      enddo
      call calc_moments_1D(tab, mx, nb_mom=4, lg=lg)
      tab(1:lg) = (tab(1:lg) -mx%mu)/mx%si

      mx%mu = 0._R8
      mx%si = 1._R8

   return
   endsubroutine profil_theo_trie_1D


   subroutine melange(tab, lg)
   implicit none
   integer(kind=i4), intent(in   )                  :: lg
   real(kind=r8)   , intent(inout), dimension(1:lg) :: tab

      real(kind=r8), dimension(1:lg) :: tmp
      integer(kind=i4) :: i

      call random_number(harvest=tmp)

      call sort_array2(tab_inout = tmp(1:lg),            &  !
                            tab1 = tab(1:lg), n = lg)       !

   return
   endsubroutine melange

endprogram test_algen
