!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: march, 07 2023
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<     **Various subroutines. Example of use.**
!<  </span>

program test_data_arch
use data_arch
use miscellaneous, only : get_unit, trans_center2corner, trans_corner2center, progress_bar_terminal
implicit none

   integer(kind=I4) :: uu, i, j
   integer(kind=I4), parameter :: nx = 11, ny = 11
   character(len=3) :: snx

   real(kind=R8), dimension(nx, ny) :: array_in, array_out

   call get_unit(uu)

   open( unit = uu, file = 'out/GLOB_VAR_VALUES.txt' )
      write(uu, '(I9,T40,a)') I4,               'I4'

      write(uu, '(I9,T40,a)') I8,               'I8'
      write(uu, '(I9,T40,a)') R4,               'R4'
      write(uu, '(I9,T40,a)') R8,               'R8'
      write(uu, '(I12,T40,a)') HIG_I4,          'HIG_I4'

      write(uu, '(I9,T40,a)') OPU,              'OPU'
      write(uu, '(I9,T40,a)') IPU,              'IPU'
      write(uu, '(I9,T40,a)') ERU,              'ERU'

      write(uu, '(E20.12,T40,a)') UN,           'UN'

      write(uu, '(E20.12,T40,a)') PI_R4,        'PI_R4'
      write(uu, '(E20.12,T40,a)') PI_R8,        'PI_R8'
      write(uu, '(E20.12,T40,a)') EPS_R4,       'EPS_R4'
      write(uu, '(E20.12,T40,a)') EPS_R8,       'EPS_R8'
      write(uu, '(E20.12,T40,a)') HIG_R8,       'HIG_R8'
      write(uu, '(E20.12,T40,a)') HIG_E8,       'HIG_E8'
      write(uu, '(E20.12,T40,a)') EPS_E8,       'EPS_E8'

      write(uu, '(I9,T40,a)') EXPO_MAX,         'EXPO_MAX'

   close( uu )

   write(snx, '(I3.3)') nx

   write(*, *) '============== INITIAL ARRAY =================='
   array_in = reshape( [ ( 0., i = 1, nx * ny ) ], [nx, ny] )
   array_in( nx/2:nx/2 + 2, ny/2:ny/2 + 2 ) = 1.

   do i = 1, nx
      write(*, '('//snx//'10I1)') ( int(array_in(i, j)), j = 1, ny )
   enddo

   write(*, *) '============ CENTER => CORNER ================='
   call trans_center2corner( tab_in = array_in, tab_out = array_out, long = nx, larg = ny)
   do i = 1, nx
      write(*, '('//snx//'10I1)') ( int(array_out(i, j)), j = 1, ny )
   enddo

   write(*, *) '============== INITIAL ARRAY =================='
   array_in = array_out
   do i = 1, nx
      write(*, '('//snx//'10I1)') ( int(array_in(i, j)), j = 1, ny )
   enddo

   write(*, *) '============ CORNER => CENTER ================='
   call trans_corner2center( tab_in = array_in, tab_out = array_out, long = nx, larg = ny)
   do i = 1, nx
      write(*, '('//snx//'10I1)') ( int(array_out(i, j)), j = 1, ny )
   enddo

   call progress_bar_terminal(val = 0, max_val = 100, init = .true.)
   do i = 10, 100, 10
      call sleep(1)
      call progress_bar_terminal(val = i, max_val = 100, init = .false.)
   enddo

stop
endprogram test_data_arch
