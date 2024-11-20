!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: November, 18 2024
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<     **Routines to work with FFTs. Example of use**
!<  </span>

program test_fftw3
!$ use omp_lib
use data_arch, only : I4, R8, R4
use miscellaneous, only : get_unit
use fftw3
implicit none


integer(kind=I4) :: i, k
integer(kind=I4) :: ni
integer(kind=I4) :: nj
integer(kind=I4) :: nb_iter
real(kind=R8)    :: error
real(kind=R4)    :: t1, t2

real(kind=R8),    dimension(:,:), allocatable :: tab     !! real array containing the information to process
real(kind=R8),    dimension(:,:), allocatable :: tabr1   !! real array containing the information to process
real(kind=R8),    dimension(:,:), allocatable :: tabr2   !! real array containing the information to process
complex(kind=R8), dimension(:,:), allocatable :: tab1    !! input array ```FORWARD```,  output array ```BACKWARD```
complex(kind=R8), dimension(:,:), allocatable :: tab2    !! input array ```BACKWARD```, output array ```FORWARD```

   ! COMPLEX -> COMPLEX -> COMPLEX

   !----------------------------------------------------------------------------------------------------------------
   ! multithread activation
   !----------------------------------------------------------------------------------------------------------------
   NB_THREADS_FFT = omp_get_max_threads()
   call fftw_plan_with_nthreads(nthreads = NB_THREADS_FFT)

   ! first case: random array forward, then backward transformed
   ! The difference ```error``` between original and processed data is calculated
   ni = 4096
   nj = 4096

   allocate( tab( 1:ni, 1:nj) )
   allocate( tab1(1:ni, 1:nj),  tab2(1:ni, 1:nj) )

   call init_fftw3( long = ni, larg = nj )

      call random_number( harvest = tab(1:ni, 1:nj) )
      tab1(1:ni, 1:nj) = cmplx( tab(1:ni, 1:nj), 0._R8, kind = R8 )

      call calc_fftw3( sens = FORWARD,  tab_in = tab1, tab_ou = tab2, long = ni, larg = nj )
      call calc_fftw3( sens = BACKWARD, tab_in = tab2, tab_ou = tab1, long = ni, larg = nj )

      error = 100*maxval( abs( (tab(1:ni, 1:nj) - real(tab1(1:ni, 1:nj), kind = R8)) / tab(1:ni, 1:nj) ) )

      write(*, *) 'C-C-C error = ', error

   call end_fftw3()

   deallocate( tab, tab1, tab2 )

   ! second case: height array forward, then backward transformed
   ! The difference between original and processed data can be assessed with the resulting "600x300_FB.dat"
   ni = 600
   nj = 300

   call read_surf( nom_fic = "./sur/600x300.dat", tab_s = tab, nx = ni, ny = nj )

   allocate( tab1(1:ni, 1:nj),  tab2(1:ni, 1:nj) ) ; tab1(1:ni, 1:nj) = cmplx( tab(1:ni, 1:nj), 0._R8, kind=R8)

   call cpu_time(t1)

   call init_fftw3( long = ni, larg = nj )

      call calc_fftw3( sens = FORWARD,  tab_in = tab1, tab_ou = tab2, long = ni, larg = nj )
      call calc_fftw3( sens = BACKWARD, tab_in = tab2, tab_ou = tab1, long = ni, larg = nj )

   call end_fftw3()

   call cpu_time(t2)

   write(*,*) t2-t1

   tab(1:ni, 1:nj) = real( tab1(1:ni, 1:nj), kind = R8 )
   call save_surf( nom_fic = "./sur/600x300_FB_comp.dat", tab_s = tab, nx = ni, ny = nj )

   deallocate( tab, tab1, tab2 )

   !----------------------------------------------------------------------------------------------------------------
   ! multithread deactivation
   ! ```NB_THREADS_MAX``` are simultaneously computed on random 512x512 arrays
   ! The difference ```error``` between original and processed data is calculated
   call fftw_plan_with_nthreads( nthreads = 1 )

   ni = 0512
   nj = 0512
   nb_iter = 100

   allocate( tab( 1:ni, 1:nj) )
   allocate( tab1(1:ni, 1:nj),  tab2(1:ni, 1:nj) )
   call get_unit(k)

   open( unit = k, file = "out/error_comp.txt" )

   NB_THREADS_FFT = omp_get_max_threads()

   call tab_init_fftw3(long = ni, larg = nj, plan_flag = FFTW_MEASURE)

   !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(NB_THREADS_FFT)

   !$OMP DO SCHEDULE (STATIC, NB_ITER/NB_THREADS_FFT) PRIVATE(tab, tab1, tab2)
   do i = 1, nb_iter

      call random_number( harvest = tab(1:ni, 1:nj) )

      tab(1:ni, 1:nj) = tab(1:ni, 1:nj) + 1._R8

      tab1(1:ni, 1:nj) = cmplx( tab(1:ni, 1:nj), 0._R8, kind = R8 )
      call tab_calc_fftw3( sens = FORWARD,  tab_in = tab1, tab_ou = tab2, long = ni, larg = nj )
      call tab_calc_fftw3( sens = BACKWARD, tab_in = tab2, tab_ou = tab1, long = ni, larg = nj )
      error = 100*maxval( abs( (tab(1:ni, 1:nj) - real(tab1(1:ni, 1:nj), kind = R8)) / tab(1:ni, 1:nj) ) )
      write(k,*) error

   enddo
   !$OMP END DO

   !$OMP END PARALLEL

   call tab_end_fftw3()

   close(k)
   deallocate( tab, tab1, tab2 )


   ! REAL -> COMPLEX -> REAL

   !----------------------------------------------------------------------------------------------------------------
   ! multithread activation
   !----------------------------------------------------------------------------------------------------------------
   NB_THREADS_FFT = omp_get_max_threads()
   call fftw_plan_with_nthreads(nthreads = NB_THREADS_FFT)

   ! first case: random array forward, then backward transformed
   ! The difference ```error``` between original and processed data is calculated
   ni = 4096
   nj = 4096

   allocate( tab( 1:ni, 1:nj) )
   allocate( tabr1(1:ni, 1:nj), tabr2(1:ni, 1:nj) )
   allocate( tab2(1:ni, 1:nj) )

   call init_fftw3_real( long = ni, larg = nj, plan_flag = FFTW_ESTIMATE )

      call random_number( harvest = tab(1:ni, 1:nj) )

      tabr1(1:ni, 1:nj) = tab(1:ni, 1:nj)

      call calc_fftw3_real_fwd( tab_in = tabr1, tab_ou = tab2, long = ni, larg = nj )
      call calc_fftw3_real_bwd( tab_in = tab2, tab_ou = tabr2, long = ni, larg = nj )

      error = 100*maxval( abs( (tab(1:ni, 1:nj) - real(tabr2(1:ni, 1:nj), kind = R8)) / tab(1:ni, 1:nj) ) )

      write(*, *) 'R-C-R: error = ', error

   call end_fftw3()

   deallocate( tab, tabr1, tabr2, tab2 )


   ! second case: height array forward, then backward transformed
   ! The difference between original and processed data can be assessed with the resulting "600x300_FB.dat"
   ni = 600
   nj = 300

   call read_surf( nom_fic = "./sur/600x300.dat", tab_s = tab, nx = ni, ny = nj )

   allocate( tabr1(1:ni, 1:nj),  tab2(1:ni, 1:nj) ) ; tabr1(1:ni, 1:nj) = tab(1:ni, 1:nj)

   call cpu_time(t1)

   call init_fftw3_real( long = ni, larg = nj, plan_flag = FFTW_ESTIMATE )

      call calc_fftw3_real_fwd( tab_in = tabr1, tab_ou = tab2,  long = ni, larg = nj )
      call calc_fftw3_real_bwd( tab_in = tab2,  tab_ou = tabr1, long = ni, larg = nj )

   call end_fftw3()

   call cpu_time(t2)

   write(*,*) t2-t1

   tab(1:ni, 1:nj) = tabr1(1:ni, 1:nj)

   call save_surf( nom_fic = "./sur/600x300_FB_real.dat", tab_s = tab, nx = ni, ny = nj )

   deallocate( tab, tabr1, tab2 )

   !----------------------------------------------------------------------------------------------------------------
   ! multithread deactivation
   ! ```NB_THREADS_MAX``` are simultaneously computed on random 512x512 arrays
   ! The difference ```error``` between original and processed data is calculated

   call fftw_plan_with_nthreads( nthreads = 1 )

   ni = 0512
   nj = 0512
   nb_iter = 100

   allocate( tab( 1:ni, 1:nj) )
   allocate( tab2( 1:ni, 1:nj) )
   allocate( tabr1(1:ni, 1:nj),  tabr2(1:ni, 1:nj) )

   call get_unit(k)

   open( unit = k, file = "out/error_real.txt" )

   NB_THREADS_FFT = omp_get_max_threads()

   call tab_init_fftw3_real( long = ni, larg = nj, plan_flag = FFTW_MEASURE )

   !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(NB_THREADS_FFT)

   !$OMP DO SCHEDULE (STATIC, NB_ITER/NB_THREADS_FFT) PRIVATE(tab, tabr1, tabr2, tab2)
   do i = 1, nb_iter

      call random_number( harvest = tab(1:ni, 1:nj) )

      tab(1:ni, 1:nj) = tab(1:ni, 1:nj) + 1._R8

      tabr1(1:ni, 1:nj) = tab(1:ni, 1:nj)

      call tab_calc_fftw3_real_fwd( tab_in = tabr1, tab_ou = tab2, long = ni, larg = nj )
      call tab_calc_fftw3_real_bwd( tab_in = tab2, tab_ou = tabr2, long = ni, larg = nj )
      error = 100*maxval( abs( (tab(1:ni, 1:nj) - tabr2(1:ni, 1:nj)) / tab(1:ni, 1:nj) ) )

      write(k,*) error

   enddo
   !$OMP END DO

   !$OMP END PARALLEL

   call tab_end_fftw3_real()

   close(k)
   deallocate( tab, tabr1, tabr2, tab2 )

stop
contains

   !=========================================================================================
   subroutine read_surf(nom_fic, tab_s, nx, ny)
   !! Subroutine that opens a surface file ```.dat```
   implicit none
   character(len=*), intent(in )                               :: nom_fic  !! *file name*
   integer(kind=I4), intent(in )                               :: nx       !! *number of pixels along x*
   integer(kind=I4), intent(in )                               :: ny       !! *number of pixels along y*
   real   (kind=R8), intent(out), dimension(:,:), allocatable  :: tab_s    !! *height array*

      integer(kind=I4) :: i, j, k
      real(kind=R8)    :: x, y

      allocate( tab_s(1:nx, 1:ny) )

      call get_unit(k)

      open( unit = k, file = trim(nom_fic), status = 'old')

         do i = 1, nx
            do j = 1, ny
               read(k, *) x, y, tab_s(i, j)
            enddo
         enddo

      close(k)

   return
   endsubroutine read_surf


   !=========================================================================================
   subroutine save_surf(nom_fic, tab_s, nx, ny)
   !! Subroutine that saves a surface file ```.dat```
   implicit none
   character(len=*), intent(in)                          :: nom_fic  !! *file name*
   integer(kind=I4), intent(in)                          :: nx       !! *number of pixels along x*
   integer(kind=I4), intent(in)                          :: ny       !! *number of pixels along y*
   real   (kind=R8), intent(in), dimension(1:nx, 1:ny)   :: tab_s    !! *height array*

      integer(kind=I4) :: i, j, k

      call get_unit(k)

      open( unit = k, file = trim(nom_fic), status = 'unknown')

         do i = 1, nx
            do j = 1, ny
               write(k, *) i, j, tab_s(i, j)
            enddo
         enddo

      close(k)

   return
   endsubroutine save_surf

endprogram test_fftw3

