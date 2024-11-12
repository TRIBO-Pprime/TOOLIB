!< author: Arthur Francisco
!  version: 1.1.0
!  date: april, 6 2023
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!     **General type definition**
!  </span>
module miscellaneous
use data_arch, only : I4, R8, OPU, IPU, ERU
implicit none

private

public :: get_unit, trans_center2corner, trans_corner2center, progress_bar_terminal

contains

   !================================================================================================
   !< @note from [John Burkardt website](https://people.sc.fsu.edu/~jburkardt/f_src)
   !     Provide for a free unit.
   !  @endnote
   !------------------------------------------------------------------------------------------------
   subroutine get_unit(iunit)
   implicit none
   integer(kind=I4), intent(out) :: iunit

      integer(kind=I4) :: i
      integer(kind=I4) :: ios
      logical(kind=I4) :: lopen

      iunit = 0
      do i = 10, 99

         if (i /= OPU .and. i /= IPU .and. i /= ERU) then
            inquire (unit = i, opened = lopen, iostat = ios)
            if (ios == 0) then
               if ( .not. lopen ) then
                  iunit = i
                  return
               endif
            endif
         endif

      enddo

   return
   endsubroutine get_unit

   subroutine trans_center2corner(tab_in, tab_out, long, larg)
   implicit none
   integer(kind=I4), intent(in )                    :: long     !! *2D array length*
   integer(kind=I4), intent(in )                    :: larg     !! *2D array width*
   class(*), intent(in ), dimension(1:long, 1:larg) :: tab_in   !! *2D array to transform*
   class(*), intent(out), dimension(1:long, 1:larg) :: tab_out  !! *transformed 2D array*

      select type(tab_in)

         type is( real(kind=R8) )
            select type(tab_out)
               type is( real(kind=R8) )
                  call trans_center2corner_real(tab_in(1:long, 1:larg), tab_out(1:long, 1:larg), long = long, larg = larg)
            endselect

         type is( complex(kind=R8) )
            select type(tab_out)
               type is( complex(kind=R8) )
                  call trans_center2corner_cmpl(tab_in(1:long, 1:larg), tab_out(1:long, 1:larg), long = long, larg = larg)
            endselect

      endselect

   return
   endsubroutine trans_center2corner


   subroutine trans_corner2center(tab_in, tab_out, long, larg)
   implicit none
   integer(kind=I4), intent(in )                    :: long     !! *2D array length*
   integer(kind=I4), intent(in )                    :: larg     !! *2D array width*
   class(*), intent(in ), dimension(1:long, 1:larg) :: tab_in   !! *2D array to transform*
   class(*), intent(out), dimension(1:long, 1:larg) :: tab_out  !! *transformed 2D array*

      select type(tab_in)

         type is( real(kind=R8) )
            select type(tab_out)
               type is( real(kind=R8) )
                  call trans_corner2center_real(tab_in(1:long, 1:larg), tab_out(1:long, 1:larg), long = long, larg = larg)
            endselect

         type is( complex(kind=R8) )
            select type(tab_out)
               type is( complex(kind=R8) )
                  call trans_corner2center_cmpl(tab_in(1:long, 1:larg), tab_out(1:long, 1:larg), long = long, larg = larg)
            endselect

      endselect

   return
   endsubroutine trans_corner2center


   subroutine trans_center2corner_real(tab_in, tab_out, long, larg)
   !================================================================================================
   !< @note Subroutine to transform an array so that the center is in the corners
   !
   !  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long     !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg     !! *2D array width*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab_in   !! *2D array to transform*
   real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: tab_out  !! *transformed 2D array*

      integer(kind=I4) :: i, j, ii, jj

      ii = 0
      jj = 0

      if ( long == 2 * (long/2) ) ii = 1
      if ( larg == 2 * (larg/2) ) jj = 1

      do j = 1, larg
         do i = 1, long
            tab_out(i, j) = tab_in( mod( i + long/2 - ii, long ) + 1, &  !
                                    mod( j + larg/2 - jj, larg ) + 1 )
         enddo
      enddo

   return
   endsubroutine trans_center2corner_real


   subroutine trans_center2corner_cmpl(tab_in, tab_out, long, larg)
   !================================================================================================
   !< @note Function to transform
   !
   !  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long     !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg     !! *2D array width*
   complex(kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab_in   !! *2D array to transform*
   complex(kind=R8), intent(out), dimension(1:long, 1:larg) :: tab_out  !! *transformed 2D array*

      integer(kind=I4) :: i, j, ii, jj

      ii = 0
      jj = 0

      if ( long == 2 * (long/2) ) ii = 1
      if ( larg == 2 * (larg/2) ) jj = 1

      do j = 1, larg
         do i = 1, long
            tab_out(i, j) = tab_in( mod( i + long/2 - ii, long ) + 1, &  !
                                    mod( j + larg/2 - jj, larg ) + 1 )
         enddo
      enddo

   return
   endsubroutine trans_center2corner_cmpl

   subroutine trans_corner2center_real(tab_in, tab_out, long, larg)
   !================================================================================================
   !< @note Function to transform an acf array so that the acf maximum is in the center
   !
   !  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long     !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg     !! *2D array width*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab_in   !! *2D array to transform*
   real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: tab_out  !! *transformed 2D array*

      integer(kind=I4) :: i, j

      do j = 1, larg
         do i = 1, long
            tab_out(i, j) = tab_in( mod( i + long/2 - 1, long ) + 1, &  !
                                    mod( j + larg/2 - 1, larg ) + 1 )
         enddo
      enddo

   return
   endsubroutine trans_corner2center_real


   subroutine trans_corner2center_cmpl(tab_in, tab_out, long, larg)
   !================================================================================================
   !< @note Function to transform an acf array so that the acf maximum is in the center
   !
   !  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long     !! *2D array length*
   integer(kind=I4), intent(in )                            :: larg     !! *2D array width*
   complex(kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab_in   !! *2D array to transform*
   complex(kind=R8), intent(out), dimension(1:long, 1:larg) :: tab_out  !! *transformed 2D array*

      integer(kind=I4) :: i, j

      do j = 1, larg
         do i = 1, long
            tab_out(i, j) = tab_in( mod( i + long/2 - 1, long ) + 1, &  !
                                    mod( j + larg/2 - 1, larg ) + 1 )
         enddo
      enddo

   return
   endsubroutine trans_corner2center_cmpl

   subroutine progress_bar_terminal(val, max_val, init)
   implicit none
   integer(kind=I4), intent(in) :: val
   integer(kind=I4), intent(in) :: max_val
   logical(kind=I4), intent(in) :: init

      character(len=102) :: bar
      integer(kind=I4)   :: ival

      if ( init ) then

         write(*, *)

         write(bar, '(a)') '[' // repeat('.', 100) // ']'

         write(*, '(a)', advance = 'no') bar

         return

      endif

      ival = nint( 99.99 * ( real(val, kind = R8) / max_val ) )

      write(bar, '(a)') '[' // repeat('*', ival) // repeat('.', 100 - ival) // ']'

      write(*, '(a)', advance = 'no') repeat(achar(8), 102) // bar

      if ( val == max_val ) then

         write(*, *) ' ... done'

         write(*, *)

      endif

   return
   endsubroutine progress_bar_terminal

endmodule miscellaneous

