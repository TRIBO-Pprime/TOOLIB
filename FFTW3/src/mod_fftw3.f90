!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: april, 9 2023
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!< **A fortran api to FFTW3**
!< </span>
!<
!< + FFT distributed on **multiple threads**
!< + **multiple FFT** simultaneously computed.
!<


module fftw3
use, intrinsic :: iso_c_binding
use data_arch, only: I4, R8, PI_R8, UN
!$ use omp_lib
implicit none

private

   integer(kind=I4) :: NB_THREADS_FFT = 4

   integer(kind=I4), parameter :: FORWARD  = +1       !! *just as suggested, it means  forward transformation required*
   integer(kind=I4), parameter :: BACKWARD = -1       !! *just as suggested, it means backward transformation required*

   integer(kind=I4), dimension(2) :: FFT_DIM = [0, 0] !! *Sizes currently allocated*

   real(kind=R8), parameter :: PAD_FFT = 1.50_R8      !! *dimension multiplier for 0-padding*

   logical(kind=I4) :: MULTI_FFTW_ALLOCATED = .false. !! *the fftw arrays are allocated and the plans are defined*
   logical(kind=I4) :: SINGL_FFTW_ALLOCATED = .false. !! *the fftw arrays are allocated and the plans are defined*

   !------------------      classical way of using FFT   -----------------------------------------------------------------------------------------
   !------------------ 1 FFT computed on several threads -----------------------------------------------------------------------------------------
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: cmp_f_i  !! *memory address of the input  array for a ```FORWARD``` transformation*
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: cmp_f_o  !! *memory address of the output array for a ```FORWARD``` transformation*

   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: cmp_b_i  !! *memory address of the input  array for a ```BACKWARD``` transformation*
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: cmp_b_o  !! *memory address of the output array for a ```BACKWARD``` transformation*

   real(C_DOUBLE), dimension(:,:), pointer :: rea_f_i             !! *memory address of the input  array for a ```FORWARD``` transformation*
   real(C_DOUBLE), dimension(:,:), pointer :: rea_b_o             !! *memory address of the output array for a ```BACKWARD``` transformation*

   type(C_PTR) :: p_f_i    !! *memory address for a plan ```FORWARD``` in  input*
   type(C_PTR) :: p_f_o    !! *memory address for a plan ```FORWARD``` in output*

   type(C_PTR) :: p_b_i    !! *memory address for a plan ```BACKWARD``` in  input*
   type(C_PTR) :: p_b_o    !! *memory address for a plan ```BACKWARD``` in output*

   type(C_PTR) :: plan_f   !! *plan  ```FORWARD```*
   type(C_PTR) :: plan_b   !! *plan ```BACKWARD```*


  !------------------      parallel FFTs   -----------------------------------------------------------------------------------------------------
  !------------------ several FFT treated on 1 thread each -------------------------------------------------------------------------------------
   !< @note
   !<   Because FFTW3 is built so that it works on the same memory zone, for concurrent executions,
   !<   a zone per thread is created.
   !< @endnote
   type tab_fftw
      complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: tab
   endtype tab_fftw

   !< @note
   !<   Because FFTW3 is built so that it works on the same memory zone, for concurrent executions,
   !<   a zone per thread is created.
   !< @endnote
   type tab_fftw_real
      real(C_DOUBLE), dimension(:,:), pointer :: tab
   endtype tab_fftw_real

   type(tab_fftw), dimension(:), allocatable :: tab_cmp_f_i       !! *array of memory addresses of the  input arrays for a ```FORWARD``` transformation*
   type(tab_fftw), dimension(:), allocatable :: tab_cmp_f_o       !! *array of memory addresses of the output arrays for a ```FORWARD``` transformation*

   type(tab_fftw), dimension(:), allocatable :: tab_cmp_b_i       !! *array of memory addresses of the  input arrays for a ```BACKWARD``` transformation*
   type(tab_fftw), dimension(:), allocatable :: tab_cmp_b_o       !! *array of memory addresses of the output arrays for a ```BACKWARD``` transformation*

   type(tab_fftw_real), dimension(:), allocatable :: tab_rea_f_i !! *array of memory addresses of the  input arrays for a ```FORWARD``` transformation*
   type(tab_fftw_real), dimension(:), allocatable :: tab_rea_f_o !! *array of memory addresses of the output arrays for a ```FORWARD``` transformation*

   type(tab_fftw_real), dimension(:), allocatable :: tab_rea_b_i !! *array of memory addresses of the  input arrays for a ```BACKWARD``` transformation*
   type(tab_fftw_real), dimension(:), allocatable :: tab_rea_b_o !! *array of memory addresses of the output arrays for a ```BACKWARD``` transformation*

   type(C_PTR), dimension(:), allocatable :: tab_p_f_i      !! *array of memory addresses for a plan ```FORWARD``` in  input*
   type(C_PTR), dimension(:), allocatable :: tab_p_f_o      !! *array of memory addresses for a plan ```FORWARD``` in output*

   type(C_PTR), dimension(:), allocatable :: tab_p_b_i      !! *array of memory addresses for a plan ```BACKWARD``` in  input*
   type(C_PTR), dimension(:), allocatable :: tab_p_b_o      !! *array of memory addresses for a plan ```BACKWARD``` in output*

   type(C_PTR), dimension(:), allocatable :: tab_plan_f     !! *plan  ```FORWARD```*
   type(C_PTR), dimension(:), allocatable :: tab_plan_b     !! *plan ```BACKWARD```*

   include "fftw3.f03"

public :: tab_init_fftw3, tab_calc_fftw3, tab_end_fftw3, fftw_plan_with_nthreads, FORWARD, BACKWARD,  &  !
          tab_init_fftw3_real, tab_calc_fftw3_real_bwd, tab_calc_fftw3_real_fwd, tab_end_fftw3_real,  &  !
              init_fftw3,     calc_fftw3,     end_fftw3,                                              &  !
              init_fftw3_real, calc_fftw3_real_fwd, calc_fftw3_real_bwd,                              &  !
              MULTI_FFTW_ALLOCATED, SINGL_FFTW_ALLOCATED, NB_THREADS_FFT, FFT_DIM,                    &  !
              apod, PAD_FFT, extend, FFTW_ESTIMATE, FFTW_MEASURE, FFTW_EXHAUSTIVE                        !

contains

   !=========================================================================================
   subroutine init_fftw3(long, larg)
   !! Subroutine to initialize the FFTW3 process *1 FFT distributed on several threads*.
   !! Complex case.
   implicit none
   integer(kind=I4), intent(in) :: long  !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg  !! *second 2D array dimension*

      call alloc_fftw3(long, larg)
      call make_plan_fftw3(long, larg)

      SINGL_FFTW_ALLOCATED = .true.
      FFT_DIM(1:2) = [long, larg]

   return
   endsubroutine init_fftw3


   !=========================================================================================
   subroutine init_fftw3_real(long, larg, plan_flag)
   !! Subroutine to initialize the FFTW3 process *1 FFT distributed on several threads*
   !! Real case.
   implicit none
   integer(kind=I4), intent(in) :: long      !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg      !! *second 2D array dimension*
   integer(kind=I4), intent(in) :: plan_flag !! *planning option, [[fftw3(module):FFTW_ESTIMATE]] for example*

      call alloc_fftw3_real(long, larg)
      call make_plan_fftw3_real(long, larg, plan_flag)

      SINGL_FFTW_ALLOCATED = .true.
      FFT_DIM(1:2) = [long, larg]

   return
   endsubroutine init_fftw3_real


   !=========================================================================================
   subroutine tab_init_fftw3(long, larg, plan_flag)
   !! Subroutine to initialize the FFTW3 process *several FFT on single thread each*
   !! Complex case.
   implicit none
   integer(kind=I4), intent(in) :: long      !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg      !! *second 2D array dimension*
   integer(kind=I4), intent(in) :: plan_flag !! *planning option, [[fftw3(module):FFTW_ESTIMATE]] for example*

      allocate( tab_cmp_f_i( 0:NB_THREADS_FFT -1) )
      allocate( tab_cmp_f_o( 0:NB_THREADS_FFT -1) )

      allocate( tab_cmp_b_i( 0:NB_THREADS_FFT -1) )
      allocate( tab_cmp_b_o( 0:NB_THREADS_FFT -1) )

      allocate( tab_p_f_i( 0:NB_THREADS_FFT -1) )
      allocate( tab_p_f_o( 0:NB_THREADS_FFT -1) )

      allocate( tab_p_b_i( 0:NB_THREADS_FFT -1) )
      allocate( tab_p_b_o( 0:NB_THREADS_FFT -1) )

      allocate( tab_plan_f(0:NB_THREADS_FFT -1) )
      allocate( tab_plan_b(0:NB_THREADS_FFT -1) )

      call tab_alloc_fftw3(long, larg)
      call tab_make_plan_fftw3(long, larg, plan_flag)

      MULTI_FFTW_ALLOCATED = .true.

   return
   endsubroutine tab_init_fftw3


   !=========================================================================================
   subroutine tab_init_fftw3_real(long, larg, plan_flag)
   !! Subroutine to initialize the FFTW3 process *several FFT on single thread each*
   !! Real case.
   implicit none
   integer(kind=I4), intent(in) :: long      !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg      !! *second 2D array dimension*
   integer(kind=I4), intent(in) :: plan_flag !! *planning option, [[fftw3(module):FFTW_ESTIMATE]] for example*

      allocate( tab_rea_f_i( 0:NB_THREADS_FFT -1) )
      allocate( tab_cmp_f_o( 0:NB_THREADS_FFT -1) )

      allocate( tab_cmp_b_i( 0:NB_THREADS_FFT -1) )
      allocate( tab_rea_b_o( 0:NB_THREADS_FFT -1) )

      allocate( tab_p_f_i( 0:NB_THREADS_FFT -1) )
      allocate( tab_p_f_o( 0:NB_THREADS_FFT -1) )

      allocate( tab_p_b_i( 0:NB_THREADS_FFT -1) )
      allocate( tab_p_b_o( 0:NB_THREADS_FFT -1) )

      allocate( tab_plan_f(0:NB_THREADS_FFT -1) )
      allocate( tab_plan_b(0:NB_THREADS_FFT -1) )

      call tab_alloc_fftw3_real(long, larg)
      call tab_make_plan_fftw3_real(long, larg, plan_flag)

      MULTI_FFTW_ALLOCATED = .true.

   return
   endsubroutine tab_init_fftw3_real


   !=========================================================================================
   !< @note
   !<   Subroutine that transforms forward or backward a double complex array. For speed reasons
   !<   FFTW will always work on the same memory area, until the plans are destroyed of course.
   !<    *1 FFT distributed on several threads*
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine calc_fftw3(sens, tab_in, tab_ou, long, larg)
   implicit none
   integer(kind=I4), intent(in )                            :: sens     !! *```=FORWARD``` or ```=BACKWARD```*
   integer(kind=I4), intent(in )                            :: long     !! *first  2D array dimension*
   integer(kind=I4), intent(in )                            :: larg     !! *second 2D array dimension*
   complex(kind=R8), dimension(1:long, 1:larg), intent(in ) :: tab_in   !! *array to transform*
   complex(kind=R8), dimension(1:long, 1:larg), intent(out) :: tab_ou   !! *transformed array*

      if ( any( FFT_DIM(1:2) /= [long, larg] ) ) then

         if ( sum(FFT_DIM(1:2)) /= 0 ) call end_fftw3()

         call fftw_plan_with_nthreads(nthreads = omp_get_num_procs())
         call init_fftw3(long = long, larg = larg)

      endif

      select case(sens)

         case(FORWARD)
            cmp_f_i(1:long, 1:larg) = tab_in(1:long, 1:larg)
            call fftw_execute_dft(plan_f, cmp_f_i(1:long, 1:larg), cmp_f_o(1:long, 1:larg))
            tab_ou(1:long, 1:larg) = cmp_f_o(1:long, 1:larg)/sqrt(real(long*larg, kind=r8))

         case(BACKWARD)
            cmp_f_i(1:long, 1:larg) = tab_in(1:long, 1:larg)
            call fftw_execute_dft(plan_b, cmp_f_i(1:long, 1:larg), cmp_f_o(1:long, 1:larg))
            tab_ou(1:long, 1:larg) = cmp_f_o(1:long, 1:larg)/sqrt(real(long*larg, kind=r8))

      endselect

   return
   endsubroutine calc_fftw3


   !=========================================================================================
   !< @note
   !<   Subroutine that transforms forward a double real array. For speed reasons
   !<   FFTW will always work on the same memory area, until the plans are destroyed of course.
   !<    *1 FFT distributed on several threads*
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine calc_fftw3_real_fwd(tab_in, tab_ou, long, larg, planner_flag)
   implicit none
   integer(kind=I4), intent(in )                            :: long           !! *first  2D array dimension*
   integer(kind=I4), intent(in )                            :: larg           !! *second 2D array dimension*
   real   (kind=R8), dimension(1:long, 1:larg), intent(in ) :: tab_in         !! *array to transform*
   complex(kind=R8), dimension(1:long, 1:larg), intent(out) :: tab_ou         !! *transformed array*
   integer(kind=I4), intent(in),                optional    :: planner_flag   !! *planning option, [[fftw3(module):FFTW_ESTIMATE]] for example*

      integer(kind=I4) :: plan_flag

      if ( .not. present(planner_flag) ) then

         plan_flag = FFTW_ESTIMATE

      else

         plan_flag = planner_flag

      endif

      if ( any( FFT_DIM(1:2) /= [long, larg] ) ) then

         if ( sum(FFT_DIM(1:2)) /= 0 ) call end_fftw3()

         call fftw_plan_with_nthreads(nthreads = omp_get_num_procs())
         call init_fftw3_real(long = long, larg = larg, plan_flag = plan_flag)

      endif

      rea_f_i(1:long, 1:larg) = tab_in(1:long, 1:larg)

      call fftw_execute_dft_r2c(plan_f, rea_f_i(1:long, 1:larg), cmp_f_o(1:long, 1:larg))

      tab_ou(1:long, 1:larg) = cmp_f_o(1:long, 1:larg) / sqrt(real(long*larg, kind=r8))

   return
   endsubroutine calc_fftw3_real_fwd


   !=========================================================================================
   !< @note
   !<   Subroutine that transforms backward a double real array. For speed reasons
   !<   FFTW will always work on the same memory area, until the plans are destroyed of course.
   !<    *1 FFT distributed on several threads*
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine calc_fftw3_real_bwd(tab_in, tab_ou, long, larg, planner_flag)
   implicit none
   integer(kind=I4), intent(in )                            :: long           !! *first  2D array dimension*
   integer(kind=I4), intent(in )                            :: larg           !! *second 2D array dimension*
   complex(kind=R8), dimension(1:long, 1:larg), intent(in ) :: tab_in         !! *array to transform*
   real   (kind=R8), dimension(1:long, 1:larg), intent(out) :: tab_ou         !! *transformed array*
   integer(kind=I4), intent(in),                optional    :: planner_flag   !! *planning option, [[fftw3(module):FFTW_ESTIMATE]] for example*

      integer(kind=I4) :: plan_flag

      if ( .not. present(planner_flag) ) then

         plan_flag = FFTW_ESTIMATE

      else

         plan_flag = planner_flag

      endif

      if ( any( FFT_DIM(1:2) /= [long, larg] ) ) then

         if ( sum(FFT_DIM(1:2)) /= 0 ) call end_fftw3()

         call fftw_plan_with_nthreads(nthreads = omp_get_num_procs())
         call init_fftw3_real(long = long, larg = larg, plan_flag = plan_flag)

      endif

      cmp_b_i(1:long, 1:larg) = tab_in(1:long, 1:larg)

      call fftw_execute_dft_c2r(plan_b, cmp_b_i(1:long, 1:larg), rea_b_o(1:long, 1:larg))

      tab_ou(1:long, 1:larg) = rea_b_o(1:long, 1:larg) / sqrt(real(long*larg, kind=r8))

   return
   endsubroutine calc_fftw3_real_bwd


   !=========================================================================================
   !< @note
   !<   Subroutine that transforms forward or backward a double complex array. For speed reasons
   !<   FFTW will always work on the same memory area, until the plans are destroyed of course.
   !<    *several FFT on single thread each*
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine tab_calc_fftw3(sens, tab_in, tab_ou, long, larg)
   implicit none
   integer(kind=I4), intent(in ) :: sens                                !! *```=FORWARD``` or ```=BACKWARD```*
   integer(kind=I4), intent(in ) :: long                                !! *first  2D array dimension*
   integer(kind=I4), intent(in ) :: larg                                !! *second 2D array dimension*
   complex(kind=R8), dimension(1:long, 1:larg), intent(in ) :: tab_in   !! *array to transform*
   complex(kind=R8), dimension(1:long, 1:larg), intent(out) :: tab_ou   !! *transformed array*

      integer(kind=I4) :: ithread

      ithread = omp_get_thread_num()
      select case(sens)

         case(FORWARD)
            tab_cmp_f_i(ithread)%tab(1:long, 1:larg) = tab_in(1:long, 1:larg)
            call fftw_execute_dft(tab_plan_f(ithread), tab_cmp_f_i(ithread)%tab(1:long, 1:larg),   &  !
                                                       tab_cmp_f_o(ithread)%tab(1:long, 1:larg))      !

            tab_ou(1:long, 1:larg) = tab_cmp_f_o(ithread)%tab(1:long, 1:larg)/sqrt(real(long*larg, kind=r8))

         case(BACKWARD)
            tab_cmp_b_i(ithread)%tab(1:long, 1:larg) = tab_in(1:long, 1:larg)
            call fftw_execute_dft(tab_plan_b(ithread), tab_cmp_b_i(ithread)%tab(1:long, 1:larg),   &  !
                                                       tab_cmp_b_o(ithread)%tab(1:long, 1:larg))      !

            tab_ou(1:long, 1:larg) = tab_cmp_b_o(ithread)%tab(1:long, 1:larg)/sqrt(real(long*larg, kind=r8))

      endselect

   return
   endsubroutine tab_calc_fftw3


   !=========================================================================================
   !< @note
   !<   Subroutine that transforms forward a real array. For speed reasons
   !<   FFTW will always work on the same memory area, until the plans are destroyed of course.
   !<    *several FFT on single thread each*
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine tab_calc_fftw3_real_fwd(tab_in, tab_ou, long, larg)
   implicit none
   integer(kind=I4), intent(in ) :: long                                !! *first  2D array dimension*
   integer(kind=I4), intent(in ) :: larg                                !! *second 2D array dimension*
   real   (kind=R8), dimension(1:long, 1:larg), intent(in ) :: tab_in   !! *array to transform*
   complex(kind=R8), dimension(1:long, 1:larg), intent(out) :: tab_ou   !! *transformed array*

      integer(kind=I4) :: ithread

      ithread = omp_get_thread_num()

      tab_rea_f_i(ithread)%tab(1:long, 1:larg) = tab_in(1:long, 1:larg)

      call fftw_execute_dft_r2c(tab_plan_f(ithread), tab_rea_f_i(ithread)%tab(1:long, 1:larg),   &  !
                                                     tab_cmp_f_o(ithread)%tab(1:long, 1:larg))      !

      tab_ou(1:long, 1:larg) = tab_cmp_f_o(ithread)%tab(1:long, 1:larg) / sqrt(real(long*larg, kind=r8))

   return
   endsubroutine tab_calc_fftw3_real_fwd


   !=========================================================================================
   !< @note
   !<   Subroutine that transforms backward a real array. For speed reasons
   !<   FFTW will always work on the same memory area, until the plans are destroyed of course.
   !<    *several FFT on single thread each*
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine tab_calc_fftw3_real_bwd(tab_in, tab_ou, long, larg)
   implicit none
   integer(kind=I4), intent(in ) :: long                                !! *first  2D array dimension*
   integer(kind=I4), intent(in ) :: larg                                !! *second 2D array dimension*
   complex(kind=R8), dimension(1:long, 1:larg), intent(in ) :: tab_in   !! *array to transform*
   real   (kind=R8), dimension(1:long, 1:larg), intent(out) :: tab_ou   !! *transformed array*

      integer(kind=I4) :: ithread

      ithread = omp_get_thread_num()

      tab_cmp_b_i(ithread)%tab(1:long, 1:larg) = tab_in(1:long, 1:larg)

      call fftw_execute_dft_c2r(tab_plan_b(ithread), tab_cmp_b_i(ithread)%tab(1:long, 1:larg),   &  !
                                                     tab_rea_b_o(ithread)%tab(1:long, 1:larg))      !

      tab_ou(1:long, 1:larg) = tab_rea_b_o(ithread)%tab(1:long, 1:larg) / sqrt(real(long*larg, kind=r8))

   return
   endsubroutine tab_calc_fftw3_real_bwd


   !=========================================================================================
   subroutine end_fftw3()
   !! FFTW3 is no more useful from here. *1 FFT distributed on several threads*
   implicit none

      if ( SINGL_FFTW_ALLOCATED ) then

         call destroy_plan_fftw3()
         call desalloc_fftw3()

      endif

      SINGL_FFTW_ALLOCATED = .false.

      FFT_DIM(1:2) = [0, 0]

   return
   endsubroutine end_fftw3


   !=========================================================================================
   subroutine tab_end_fftw3()
   !! FFTW3 is no more useful from here. *several FFT on single thread each*
   implicit none

      if ( MULTI_FFTW_ALLOCATED ) then

         call tab_destroy_plan_fftw3()
         call tab_desalloc_fftw3()

         deallocate(   tab_p_f_i,   tab_p_f_o,   tab_p_b_i,   tab_p_b_o )
         deallocate( tab_cmp_f_i, tab_cmp_f_o, tab_cmp_b_i, tab_cmp_b_o )
         deallocate( tab_plan_f, tab_plan_b )

      endif

      MULTI_FFTW_ALLOCATED = .false.

      FFT_DIM(1:2) = [0, 0]

   return
   endsubroutine tab_end_fftw3


   !=========================================================================================
   subroutine tab_end_fftw3_real()
   !! FFTW3 is no more useful from here. *several FFT on single thread each*
   implicit none

      call tab_destroy_plan_fftw3()
      call tab_desalloc_fftw3()

      deallocate(   tab_p_f_i,   tab_p_f_o,   tab_p_b_i,   tab_p_b_o )
      deallocate( tab_rea_f_i, tab_cmp_f_o, tab_cmp_b_i, tab_rea_b_o )
      deallocate( tab_plan_f, tab_plan_b )

      MULTI_FFTW_ALLOCATED = .false.

      FFT_DIM(1:2) = [0, 0]

   return
   endsubroutine tab_end_fftw3_real


   !=========================================================================================
   !< @note
   !<   Allocation of the memory needed by the transformations, forward and backward.
   !<    *1 FFT distributed on several threads*
   !<
   !<   The space remains allocated as long as transformations are needed.
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine alloc_fftw3(long, larg)
   implicit none
   integer(kind=I4), intent(in) :: long   !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg   !! *second 2D array dimension*

      ! forward
      p_f_i = fftw_alloc_complex(int(long*larg, C_SIZE_T))
      p_f_o = fftw_alloc_complex(int(long*larg, C_SIZE_T))
      call c_f_pointer(p_f_i, cmp_f_i, (/long,larg/))
      call c_f_pointer(p_f_o, cmp_f_o, (/long,larg/))

      ! backward
      p_b_i = fftw_alloc_complex(int(long*larg, C_SIZE_T))
      p_b_o = fftw_alloc_complex(int(long*larg, C_SIZE_T))
      call c_f_pointer(p_b_i, cmp_b_i, (/long,larg/))
      call c_f_pointer(p_b_o, cmp_b_o, (/long,larg/))

   return
   endsubroutine alloc_fftw3


   !=========================================================================================
   !< @note
   !<   Allocation of the memory needed by the transformations, forward and backward, for the
   !<   real case. *1 FFT distributed on several threads*
   !<
   !<   The space remains allocated as long as transformations are needed.
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine alloc_fftw3_real(long, larg)
   implicit none
   integer(kind=I4), intent(in) :: long   !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg   !! *second 2D array dimension*

      ! forward
      p_f_i = fftw_alloc_real(int(long*larg, C_SIZE_T))
      p_f_o = fftw_alloc_complex(int(long*larg, C_SIZE_T))
      call c_f_pointer(p_f_i, rea_f_i, (/long,larg/))
      call c_f_pointer(p_f_o, cmp_f_o, (/long,larg/))

      ! backward
      p_b_i = fftw_alloc_complex(int(long*larg, C_SIZE_T))
      p_b_o = fftw_alloc_real(int(long*larg, C_SIZE_T))
      call c_f_pointer(p_b_i, cmp_b_i, (/long,larg/))
      call c_f_pointer(p_b_o, rea_b_o, (/long,larg/))

   return
   endsubroutine alloc_fftw3_real


   !=========================================================================================
   !< @note
   !<   Allocation of the memory needed by the transformations, forward and backward.
   !<    *several FFT on single thread each*
   !<
   !<   The space remains allocated as long as transformations are needed.
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine tab_alloc_fftw3(long, larg)
   implicit none
   integer(kind=I4), intent(in) :: long   !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg   !! *second 2D array dimension*

      integer(kind=I4) :: ithread

      ! forward
      do ithread = 0, NB_THREADS_FFT -1

         tab_p_f_i(ithread) = fftw_alloc_complex(int(long*larg, C_SIZE_T))
         tab_p_f_o(ithread) = fftw_alloc_complex(int(long*larg, C_SIZE_T))
         call c_f_pointer(tab_p_f_i(ithread), tab_cmp_f_i(ithread)%tab, (/long,larg/))
         call c_f_pointer(tab_p_f_o(ithread), tab_cmp_f_o(ithread)%tab, (/long,larg/))

      enddo

      ! backward
      do ithread = 0, NB_THREADS_FFT -1

         tab_p_b_i(ithread) = fftw_alloc_complex(int(long*larg, C_SIZE_T))
         tab_p_b_o(ithread) = fftw_alloc_complex(int(long*larg, C_SIZE_T))
         call c_f_pointer(tab_p_b_i(ithread), tab_cmp_b_i(ithread)%tab, (/long,larg/))
         call c_f_pointer(tab_p_b_o(ithread), tab_cmp_b_o(ithread)%tab, (/long,larg/))

      enddo

   return
   endsubroutine tab_alloc_fftw3


   !=========================================================================================
   !< @note
   !<   Allocation of the memory needed by the transformations, forward and backward, for the
   !<   real case. *several FFT on single thread each*
   !<
   !<   The space remains allocated as long as transformations are needed.
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine tab_alloc_fftw3_real(long, larg)
   implicit none
   integer(kind=I4), intent(in) :: long   !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg   !! *second 2D array dimension*

      integer(kind=I4) :: ithread

      ! forward
      do ithread = 0, NB_THREADS_FFT -1

         tab_p_f_i(ithread) = fftw_alloc_real(int(long*larg, C_SIZE_T))
         tab_p_f_o(ithread) = fftw_alloc_complex(int(long*larg, C_SIZE_T))
         call c_f_pointer(tab_p_f_i(ithread), tab_rea_f_i(ithread)%tab, (/long,larg/))
         call c_f_pointer(tab_p_f_o(ithread), tab_cmp_f_o(ithread)%tab, (/long,larg/))

      enddo

      ! backward
      do ithread = 0, NB_THREADS_FFT -1

         tab_p_b_i(ithread) = fftw_alloc_complex(int(long*larg, C_SIZE_T))
         tab_p_b_o(ithread) = fftw_alloc_real(int(long*larg, C_SIZE_T))
         call c_f_pointer(tab_p_b_i(ithread), tab_cmp_b_i(ithread)%tab, (/long,larg/))
         call c_f_pointer(tab_p_b_o(ithread), tab_rea_b_o(ithread)%tab, (/long,larg/))

      enddo

   return
   endsubroutine tab_alloc_fftw3_real


   !=========================================================================================
   !< @note
   !<   When no more transformation is needed, the memory is released.
   !<    *1 FFT distributed on several threads*
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine desalloc_fftw3()
   implicit none

      ! forward
      call fftw_free(p_f_i) ; p_f_i = C_NULL_PTR
      call fftw_free(p_f_o) ; p_f_o = C_NULL_PTR

      ! backward
      call fftw_free(p_b_i) ; p_b_i = C_NULL_PTR
      call fftw_free(p_b_o) ; p_b_o = C_NULL_PTR

   return
   endsubroutine desalloc_fftw3


   !=========================================================================================
   !< @note
   !<   When no more transformation is needed, the memory is released.
   !<     *several FFT on single thread each*
   !< @endnote
   !-----------------------------------------------------------------------------------------
   subroutine tab_desalloc_fftw3()
   implicit none

      integer(kind=I4) :: ithread

      do ithread = 0, NB_THREADS_FFT -1

        ! forward
         call fftw_free(tab_p_f_i(ithread)) ; tab_p_f_i(ithread) = C_NULL_PTR
         call fftw_free(tab_p_f_o(ithread)) ; tab_p_f_o(ithread) = C_NULL_PTR
         ! backward
         call fftw_free(tab_p_b_i(ithread)) ; tab_p_b_i(ithread) = C_NULL_PTR
         call fftw_free(tab_p_b_o(ithread)) ; tab_p_b_o(ithread) = C_NULL_PTR

      enddo

   return
   endsubroutine tab_desalloc_fftw3


   !=========================================================================================
   !< @note
   !<   Creates forward and backward plans. *1 FFT distributed on several threads*
   !<
   !<   Until no more transformation is needed, the plans remain as they are.
   !< @endnote
   !<
   !< @warning
   !<    In C, the order line/column is reversed, so the 2nd dimension ```larg``` of the array
   !<    is first provided in ```fftw_plan_dft_2d```
   !<
   !<    [calling from fortran](http://www.fftw.org/doc/Calling-FFTW-from-Modern-Fortran.html#Calling-FFTW-from-Modern-Fortran)
   !< @endwarning
   !-----------------------------------------------------------------------------------------
   subroutine make_plan_fftw3(long, larg)
   implicit none
   integer(kind=I4), intent(in) :: long   !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg   !! *second 2D array dimension*

      ! forward
      plan_f = fftw_plan_dft_2d(larg, long, cmp_f_i, cmp_f_o, FFTW_FORWARD,  flags=FFTW_ESTIMATE)

      ! backward
      plan_b = fftw_plan_dft_2d(larg, long, cmp_b_i, cmp_b_o, FFTW_BACKWARD, flags=FFTW_ESTIMATE)

   return
   endsubroutine make_plan_fftw3


   !=========================================================================================
   !< @note
   !<   Creates forward and backward plans. *1 FFT distributed on several threads*
   !<
   !<   Until no more transformation is needed, the plans remain as they are.
   !< @endnote
   !<
   !< @warning
   !<    In C, the order line/column is reversed, so the 2nd dimension ```larg``` of the array
   !<    is first provided in ```fftw_plan_dft_2d```
   !<
   !<    [calling from fortran](http://www.fftw.org/doc/Calling-FFTW-from-Modern-Fortran.html#Calling-FFTW-from-Modern-Fortran)
   !< @endwarning
   !-----------------------------------------------------------------------------------------
   subroutine make_plan_fftw3_real(long, larg, plan_flag)
   implicit none
   integer(kind=I4), intent(in) :: long         !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg         !! *second 2D array dimension*
   integer(kind=I4), intent(in) :: plan_flag    !! *planning option, [[fftw3(module):FFTW_ESTIMATE]] for example*

      ! forward
      plan_f = fftw_plan_dft_r2c_2d(n0 = larg, n1 = long, in = rea_f_i, out = cmp_f_o, flags = plan_flag)

      ! backward
      plan_b = fftw_plan_dft_c2r_2d(n0 = larg, n1 = long, in = cmp_b_i, out = rea_b_o, flags = plan_flag)

   return
   endsubroutine make_plan_fftw3_real


   !=========================================================================================
   !< @note
   !<   Creates forward and backward plans. *several FFT on single thread each*
   !<
   !<   Until no more transformation is needed, the plans remain as they are.
   !< @endnote
   !<
   !< @warning
   !<    In C, the order line/column is reversed, so the 2nd dimension ```larg``` of the array
   !<    is first provided in ```fftw_plan_dft_2d```
   !<
   !<    [calling from fortran](http://www.fftw.org/doc/Calling-FFTW-from-Modern-Fortran.html#Calling-FFTW-from-Modern-Fortran)
   !< @endwarning
   !-----------------------------------------------------------------------------------------
   subroutine tab_make_plan_fftw3(long, larg, plan_flag)
   implicit none
   integer(kind=I4), intent(in) :: long         !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg         !! *second 2D array dimension*
   integer(kind=I4), intent(in) :: plan_flag    !! *planning option, [[fftw3(module):FFTW_ESTIMATE]] for example*

      integer(kind=I4) :: ithread

      do ithread = 0, NB_THREADS_FFT -1

         ! forward
         tab_plan_f(ithread) = fftw_plan_dft_2d(larg, long, tab_cmp_f_i(ithread)%tab,  &  !
                                                            tab_cmp_f_o(ithread)%tab,  &  !
                                                            FORWARD, plan_flag)       !

         ! backward
         tab_plan_b(ithread) = fftw_plan_dft_2d(larg, long, tab_cmp_b_i(ithread)%tab,  &  !
                                                            tab_cmp_b_o(ithread)%tab,  &  !
                                                            BACKWARD, plan_flag)      !

      enddo

   return
   endsubroutine tab_make_plan_fftw3


   !=========================================================================================
   !< @note
   !<   Creates forward and backward plans. *several FFT on single thread each*
   !<
   !<   Until no more transformation is needed, the plans remain as they are.
   !< @endnote
   !<
   !< @warning
   !<    In C, the order line/column is reversed, so the 2nd dimension ```larg``` of the array
   !<    is first provided in ```fftw_plan_dft_2d```
   !<
   !<    [calling from fortran](http://www.fftw.org/doc/Calling-FFTW-from-Modern-Fortran.html#Calling-FFTW-from-Modern-Fortran)
   !< @endwarning
   !-----------------------------------------------------------------------------------------
   subroutine tab_make_plan_fftw3_real(long, larg, plan_flag)
   implicit none
   integer(kind=I4), intent(in) :: long         !! *first  2D array dimension*
   integer(kind=I4), intent(in) :: larg         !! *second 2D array dimension*
   integer(kind=I4), intent(in) :: plan_flag    !! *planning option, [[fftw3(module):FFTW_ESTIMATE]] for example*

      integer(kind=I4) :: ithread

      do ithread = 0, NB_THREADS_FFT -1

         ! forward
         tab_plan_f(ithread) = fftw_plan_dft_r2c_2d(larg, long, tab_rea_f_i(ithread)%tab,  &  !
                                                                tab_cmp_f_o(ithread)%tab,  &  !
                                                                plan_flag)                    !

         ! backward
         tab_plan_b(ithread) = fftw_plan_dft_c2r_2d(larg, long, tab_cmp_b_i(ithread)%tab,  &  !
                                                                tab_rea_b_o(ithread)%tab,  &  !
                                                                plan_flag)                    !

      enddo

   return
   endsubroutine tab_make_plan_fftw3_real


   !=========================================================================================
   subroutine destroy_plan_fftw3()
   !! Plans are no more needed as no additional transformation will occur. *1 FFT distributed on several threads*
   implicit none

      ! forward
      call fftw_destroy_plan(plan_f)

      ! backward
      call fftw_destroy_plan(plan_b)

   return
   endsubroutine destroy_plan_fftw3


   !=========================================================================================
   subroutine tab_destroy_plan_fftw3()
   !! Plans are no more needed as no additional transformation will occur. *several FFT on single thread each*
   implicit none

      integer(kind=I4) :: ithread

      do ithread = 0, NB_THREADS_FFT -1

         ! forward
         call fftw_destroy_plan(tab_plan_f(ithread))

         ! backward
         call fftw_destroy_plan(tab_plan_b(ithread))

      enddo

   return
   endsubroutine tab_destroy_plan_fftw3


   !=========================================================================================
   !< @note Function that extends an array for FFT processing.
   !<
   !< + nx2 = 2 * ( nint(PAD_FFT_FILTER * nx)/2 )
   !< + ny2 = 2 * ( nint(PAD_FFT_FILTER * ny)/2 )
   !<
   !<  @endnote
   !----------------------------------------------------------------------------------------
   subroutine extend(tab_in, tab_out, nx, ny, nx2, ny2, ext, type_apo)
   implicit none
   integer(kind=I4), intent(in )                          :: nx       !! *2D input array length*
   integer(kind=I4), intent(in )                          :: ny       !! *2D input array width*
   integer(kind=I4), intent(in )                          :: nx2      !! *2D output array length*
   integer(kind=I4), intent(in )                          :: ny2      !! *2D output array width*
   real   (kind=R8), intent(in ), dimension(1:nx,  1:ny ) :: tab_in   !! *input array*
   real   (kind=R8), intent(out), dimension(1:nx2, 1:ny2) :: tab_out  !! *apodized array*
   character(len=*), intent(in )                          :: ext      !! *extension*
   character(len=*), intent(in ), optional                :: type_apo !! *apodization type*

      integer(kind=I4) :: i, j, ibx, iby, iex, iey

      real(kind=R8), dimension(:,:), allocatable :: tab_tmp

      ibx = ceiling( (nx2 - nx)/2. ) ; iex = ibx + nx - 1
      iby = ceiling( (ny2 - ny)/2. ) ; iey = iby + ny - 1

      allocate( tab_tmp(1:nx2, 1:ny2) )

      tab_tmp(1:nx2, 1:ny2) = 0

      tab_tmp(ibx:iex, iby:iey) = tab_in(1:nx, 1:ny)

      select case ( ext )

         case( 'symmetry' )

            do i = 1, ibx - 1
              tab_tmp(ibx - i, iby:iey) = tab_tmp(ibx + i, iby:iey)
            enddo

            do i = iex + 1, nx2
              tab_tmp(i, iby:iey) = tab_tmp(iex - (i - iex), iby:iey)
            enddo

            do j = 1, iby - 1
              tab_tmp(1:nx2, iby - j) = tab_tmp(1:nx2, iby + j)
            enddo

            do j = iey + 1, ny2
              tab_tmp(1:nx2, j) = tab_tmp(1:nx2, iey - (j - iey))
            enddo

         case( 'constant' )

            do i = 1, ibx - 1
              tab_tmp(i, iby:iey) = tab_tmp(ibx, iby:iey)
            enddo

            do i = iex + 1, nx2
              tab_tmp(i, iby:iey) = tab_tmp(iex, iby:iey)
            enddo

            do j = 1, iby - 1
              tab_tmp(1:nx2, j) = tab_tmp(1:nx2, iby)
            enddo

            do j = iey + 1, ny2
              tab_tmp(1:nx2, j) = tab_tmp(1:nx2, iey)
            enddo

         case( 'zero' )

      endselect

      if ( present(type_apo) ) then

         call apod(  tab_in = tab_tmp(1:nx2, 1:ny2),  &  !
                    tab_out = tab_out(1:nx2, 1:ny2),  &  !
                       long = nx2,                    &  !
                       larg = ny2,                    &  !
                   type_apo = type_apo )                 !

      else

         tab_out(1:nx2, 1:ny2) = tab_tmp(1:nx2, 1:ny2)

      endif

      deallocate( tab_tmp )

   return
   endsubroutine extend


   !=========================================================================================
   !< @note Function that returns an apodized array.
   !<
   !<   To prevent gaps from appearing after FFT (because of non periodic waves), the surface must
   !<   be transformed, but not too much ...
   !<
   !<  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine apod(tab_in, tab_out, long, larg, type_apo, param)
   implicit none
   integer(kind=I4), intent(in )                            :: long     !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg     !! *2D array width*
   character(len=*), intent(in )                            :: type_apo !! *apodization type*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab_in   !! *input array*
   real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: tab_out  !! *apodized array*
   real   (kind=R8), intent(in ), optional                  :: param    !! *apodized array*

      real   (kind=R8) :: a0, a1, a2, ro, u, v, fct, pun, mun, eps, mi, mj, li, lj, ri, rj, fcti
      integer(kind=I4) :: i, j, dis, ii, jj, lo2, la2

      select case ( type_apo(1:6) )

         case('blackm')
            a0 = 21./50.
            a1 = 01./02.
            a2 = 02./25.
            ri = real(long)  ; rj = real(larg)
            mi = (ri + UN)/2 ; mj = (rj + UN)/2
            li = (ri - UN)/2 ; lj = (rj - UN)/2
            do i = 1, long
            u = (i - mi)/li
               do j = 1, larg
                  v = (j - mj)/lj
                  fct = ( a0 + a1 * cos(PI_R8*u) + a2 * cos(2*PI_R8*u) ) *    &  !
                        ( a0 + a1 * cos(PI_R8*v) + a2 * cos(2*PI_R8*v) )         !
                  tab_out(i, j) = tab_in(i, j)*fct
               enddo
            enddo

         case('tuckey') !http://en.wikipedia.org/wiki/Window_function#Tukey_window
            mun = -UN
            pun = +UN
            eps = 0.25_R8

            if ( present(param) ) eps = param

            dis = nint(eps*long/2)
            tab_out(1:long, 1:larg) = UN

            do i = 0, dis
               ro = 2. * i / ( eps*(long-1) ) - UN
               tab_out(i+1, 1:larg) = 0.5_R8*( UN + cos(PI_R8*ro) )
            enddo
            do i = long - 1 - dis, long - 1
               ro = 2.*i/( eps*(long-1) ) + UN - 2./eps
               tab_out(i+1, 1:larg) = 0.5_R8*( UN + cos(PI_R8*ro) )
            enddo

            dis = nint(eps*larg/2)
            do j = 0, dis
               ro = 2. * j / ( eps*(larg-1) ) - UN
               tab_out(1:long, j+1) = tab_out(1:long, j+1) * 0.5_R8 * ( UN + cos(PI_R8*ro) )
            enddo
            do j = larg-1 -dis, larg-1
               ro = 2. * j / ( eps*(larg-1) ) + UN - 2. / eps
               tab_out(1:long, j+1) = tab_out(1:long, j+1) * 0.5_R8 * ( UN + cos(PI_R8*ro) )
            enddo
            tab_out(1:long, 1:larg) = tab_out(1:long, 1:larg) * tab_in(1:long, 1:larg)

         case('hann__')
            do i = 1, long
               fcti = 0.5 * (1.0 - cos(2.0 * PI_R8 * (i - 1) / (long - 1)))
               do j = 1, larg
                  fct = fcti * 0.5 * (1.0 - cos(2.0 * PI_R8 * (j - 1) / (larg - 1)))
                  tab_out(i, j) = tab_in(i, j)*fct
               enddo
            enddo

         case('welch_')
            ri  = long / 2.000001_R8
            rj  = larg / 2.000001_R8
            lo2 = ceiling( ri )
            la2 = ceiling( rj )
            do ii = lo2 - long, lo2 - 1
               u = (ii / ri)**2
               i = max(ii + lo2, 1)
               do jj = la2 - larg, la2 - 1
                  v = (jj / rj)**2
                  j = max(jj + la2, 1)
                  if ( u + v  > UN ) then
                     tab_out(i, j) = 0
                     cycle
                  endif
                  tab_out(i, j) = tab_in(i, j) * ( 1._R8 - (u + v) )
               enddo
            enddo

         case('no_apo')
            tab_out(1:long, 1:larg) = tab_in(1:long, 1:larg)

         case default
            stop 'apod, apodization type bad choice'

      endselect

   return
   endsubroutine apod

!~    subroutine sample_grid(w, h, wf, hf, tab_in, tab_ou)
!~    implicit none
!~    integer(kind=I4), intent(in) :: w   !! *input surface width*
!~    integer(kind=I4), intent(in) :: h   !! *input surface height*
!~    integer(kind=I4), intent(in) :: wf  !! *output surface width*
!~    integer(kind=I4), intent(in) :: hf  !! *output surface height*
!~    real   (kind=R8), dimension(1:wf, 1:hf), intent( in) :: tab_in  !! *input surface*
!~    real   (kind=R8), dimension(1:w , 1:hf), intent(out) :: tab_ou  !! *output surface*

!~       integer(kind=I4) :: iw, ih, iwf, ihf
!~       real   (kind=R8) :: dw, udw, dh, udh, h1, h2, h3, h4, hh

!~       do iw = 1, w

!~          iwf = floor( real( wf * (iw-1), kind = R8) / w ) + 1
!~          dw  =        real( wf * (iw-1), kind = R8) / w   + 1 - iwf
!~          udw = 1._R8 - dw

!~          do ih = 1, h

!~             ihf = floor( real( hf * (ih-1), kind = R8) / h ) + 1
!~             dh  =        real( hf * (ih-1), kind = R8) / h   + 1 - ihf
!~             udh = 1._R8 - dh

!~             h1 = tab_in(iwf    , ihf    )
!~             h2 = tab_in(iwf + 1, ihf    )
!~             h3 = tab_in(iwf + 1, ihf + 1)
!~             h4 = tab_in(iwf    , ihf + 1)

!~             hh = h1 * udw * udh + &  !
!~                  h2 *  dw * udh + &  !
!~                  h3 *  dw *  dh + &  !
!~                  h4 * udw *  dh

!~             tab_ou(iw, ih) = hh

!~          enddo

!~       enddo

!~    return
!~    endsubroutine sample_grid

endmodule fftw3
