!< author: Arthur Francisco
!<  version: 1.0.1
!<  date: feb, 24 2023
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<     **Various routines to sort real/integer arrays**
!<  </span>

module sort_arrays
use data_arch, only : I4, R8
implicit none

private

public :: init_order, sort_array2

contains

   subroutine sort_array2(tab_inout, tab0, tab1, tab2, tab3, n)
   !! Sort 1D arrays, real or integer, according the first one
   implicit none
   integer(kind=I4), intent(in)                            :: n            !! *size of the arrays*
   class(*)        , intent(inout), dimension(n)           :: tab_inout    !! *reference array to sort*
   integer(kind=I4), intent(inout), dimension(n), optional :: tab0         !! *second array to sort according the order of the first one*
   class(*)        , intent(inout), dimension(n), optional :: tab1         !! *third array to sort according the order of the first one*
   class(*)        , intent(inout), dimension(n), optional :: tab2         !! *4th array to sort according the order of the first one*
   class(*)        , intent(inout), dimension(n), optional :: tab3         !! *5th array to sort according the order of the first one*

      integer(kind=I4), allocatable, dimension(:) :: tab_order

      allocate( tab_order(1:n) )

      if ( .not.present(tab0) ) then

         call init_order(order = tab_order(1:n), n = n)

      else

         tab_order(1:n) = tab0(1:n)

      endif

      select type(tab_inout)

         type is( integer(kind=I4) )
            call sort_array_integer_with_order(g = 1, d = n, itabref = tab_inout(1:n), order = tab_order(1:n))

         type is( real(kind=R8) )
            call sort_array_real_with_order   (g = 1, d = n, rtabref = tab_inout(1:n), order = tab_order(1:n))

      endselect

      if ( present(tab1) ) call change_array_order(tab_inout = tab1(1:n), order = tab_order(1:n), n = n)
      if ( present(tab2) ) call change_array_order(tab_inout = tab2(1:n), order = tab_order(1:n), n = n)
      if ( present(tab3) ) call change_array_order(tab_inout = tab3(1:n), order = tab_order(1:n), n = n)

      if ( present(tab0) ) then

          tab0(1:n) = tab_order(1:n)

      endif

      deallocate( tab_order )

   return
   endsubroutine sort_array2


   subroutine change_array_order(tab_inout, order, n)
   !! Given an order vector, sort a real or integer vector
   implicit none
   integer(kind=I4), intent(in)                  :: n          !! *size of the arrays*
   class(*)        , intent(inout), dimension(n) :: tab_inout  !! *array to sort*
   integer(kind=I4), intent(inout), dimension(n) :: order      !! *order vector*

      integer(kind=I4) :: i

      integer(kind=I4), allocatable, dimension(:) :: tab_int
      real(kind=R8),    allocatable, dimension(:) :: tab_real


      select type(tab_inout)

         type is( integer(kind=I4) )

            allocate( tab_int(1:n) )

            tab_int(1:n) = tab_inout(1:n)

            do i = 1, n

               tab_inout( i ) = tab_int( order(i) )

            enddo

            deallocate( tab_int )

         type is( real(kind=R8) )

            allocate( tab_real(1:n) )

            tab_real(1:n) = tab_inout(1:n)

            do i = 1, n

               tab_inout( i ) = tab_real( order(i) )

            enddo

            deallocate( tab_real )

      endselect

   return
   endsubroutine change_array_order


   subroutine init_order(order, n)
   !! Vector initialization: 1 ... n
   implicit none
   integer(kind=I4), intent(in)                :: n      !! *size of the vector*
   integer(kind=I4), dimension(n), intent(out) :: order  !! *order vector*

      integer(kind=I4) :: i

      order = [ integer(kind=I4) :: (i, i = 1, n) ]

   return
   endsubroutine init_order

   !=========================================================================================
   recursive subroutine sort_array_integer_with_order(g, d, itabref, order)
   !! Sort a vector of integers and store the order
   implicit none
   integer(kind=I4), intent(in   )               :: g          !! *left index*
   integer(kind=I4), intent(in   )               :: d          !! *right index*
   integer(kind=I4), intent(inout), dimension(:) :: itabref    !! *vector to sort*
   integer(kind=I4), intent(inout), dimension(:) :: order      !! *sort order*

      integer(kind=I4) :: i, j, mil, itmp
      integer(kind=I4) :: tmp, cle

      i = g
      j = d
      mil = (g+d)/2
      cle = itabref(mil)

      if (g>=d) return

      do while (i<=j)
         do while (itabref(i)<cle)
            i = i + 1
         enddo
         do while (itabref(j)>cle)
            j = j - 1
         enddo
         if (i<=j) then
            ! échange des éléments du tableau
            tmp = itabref(i)
            itabref(i) = itabref(j)
            itabref(j) = tmp

            ! échange des éléments du tableau
            itmp     = order(i)
            order(i) = order(j)
            order(j) = itmp

            ! échange des éléments du vecteur position
            i = i + 1
            j = j - 1
         endif
      enddo

      if (g<j) call sort_array_integer_with_order(g, j, itabref, order)
      if (d>i) call sort_array_integer_with_order(i, d, itabref, order)

   return
   endsubroutine sort_array_integer_with_order

   !=========================================================================================
   recursive subroutine sort_array_real_with_order(g, d, rtabref, order)
   !! Sort a vector of reals and store the order
   implicit none
   integer(kind=I4), intent(in)                  :: g          !! *left index*
   integer(kind=I4), intent(in)                  :: d          !! *right index*
   real(kind=R8),    intent(inout), dimension(:) :: rtabref    !! *vector to sort*
   integer(kind=I4), intent(inout), dimension(:) :: order      !! *sort order*

      integer(kind=I4) :: i, j, mil, itmp
      real(kind=R8)    :: tmp, cle

      i = g
      j = d
      mil = (g+d)/2
      cle = rtabref(mil)

      if (g>=d) return

      do while (i<=j)
         do while (rtabref(i)<cle)
            i = i + 1
         enddo
         do while (rtabref(j)>cle)
            j = j - 1
         enddo
         if (i<=j) then
            ! échange des éléments du tableau
            tmp = rtabref(i)
            rtabref(i) = rtabref(j)
            rtabref(j) = tmp

            ! échange des éléments du tableau
            itmp     = order(i)
            order(i) = order(j)
            order(j) = itmp

            ! échange des éléments du vecteur position
            i = i + 1
            j = j - 1
         endif
      enddo

      if (g<j) call sort_array_real_with_order(g, j, rtabref, order)
      if (d>i) call sort_array_real_with_order(i, d, rtabref, order)

   return
   endsubroutine sort_array_real_with_order


   !=========================================================================================
   recursive subroutine sort_array_integer(g, d, itabref)
   !! Sort a vector of integers
   implicit none
   integer(kind=I4), intent(in   )               :: g          !! *left index*
   integer(kind=I4), intent(in   )               :: d          !! *right index*
   integer(kind=I4), intent(inout), dimension(:) :: itabref    !! *vector to sort*

      integer(kind=I4) :: i, j, mil
      integer(kind=I4) :: tmp, cle

      i = g
      j = d
      mil = (g+d)/2
      cle = itabref(mil)

      if (g>=d) return

      do while (i<=j)
         do while (itabref(i)<cle)
            i = i + 1
         enddo
         do while (itabref(j)>cle)
            j = j - 1
         enddo
         if (i<=j) then
            ! échange des éléments du tableau
            tmp = itabref(i)
            itabref(i) = itabref(j)
            itabref(j) = tmp

            ! échange des éléments du vecteur position
            i = i + 1
            j = j - 1
         endif
      enddo

      if (g<j) call sort_array_integer(g, j, itabref)
      if (d>i) call sort_array_integer(i, d, itabref)

   return
   endsubroutine sort_array_integer

   !=========================================================================================
   recursive subroutine sort_array_real(g, d, rtabref)
   !! Sort a vector of reals
   implicit none
   integer(kind=I4), intent(in)                  :: g          !! *left index*
   integer(kind=I4), intent(in)                  :: d          !! *right index*
   real(kind=R8),    intent(inout), dimension(:) :: rtabref    !! *vector to sort*

      integer(kind=I4) :: i, j, mil
      real(kind=R8)    :: tmp, cle

      i = g
      j = d
      mil = (g+d)/2
      cle = rtabref(mil)

      if (g>=d) return

      do while (i<=j)
         do while (rtabref(i)<cle)
            i = i + 1
         enddo
         do while (rtabref(j)>cle)
            j = j - 1
         enddo
         if (i<=j) then
            ! échange des éléments du tableau
            tmp = rtabref(i)
            rtabref(i) = rtabref(j)
            rtabref(j) = tmp

            ! échange des éléments du vecteur position
            i = i + 1
            j = j - 1
         endif
      enddo

      if (g<j) call sort_array_real(g, j, rtabref)
      if (d>i) call sort_array_real(i, d, rtabref)

   return
   endsubroutine sort_array_real

!~    !-----------------------------------------------------------------------------------------
!~    recursive subroutine sort_int_1int_1real(g, d, itabref, itab1, rtab2)
!~    implicit none
!~    integer(kind=I4), intent(in) :: g, d
!~    integer(kind=I4), dimension(:), intent(inout) :: itabref
!~    integer(kind=I4), dimension(:), intent(inout) :: itab1
!~    real(kind=R8), dimension(:), intent(inout)    :: rtab2
!~       integer(kind=I4) :: i, j, mil, cle, itmp
!~       real(kind=R8)    :: rtmp
!~       i = g
!~       j = d
!~       mil = (g+d)/2
!~       cle = itabref(mil)

!~       if (g>=d) return

!~       do while (i<=j)
!~          do while (itabref(i)<cle)
!~             i = i + 1
!~          enddo
!~          do while (itabref(j)>cle)
!~             j = j - 1
!~          enddo
!~          if (i<=j) then
!~             ! échange des éléments du tableau
!~             itmp       = itabref(i)
!~             itabref(i) = itabref(j)
!~             itabref(j) = itmp
!~             ! échange des éléments du vecteur 2
!~             itmp     = itab1(i)
!~             itab1(i) = itab1(j)
!~             itab1(j) = itmp
!~             ! échange des éléments du vecteur 3
!~             rtmp     = rtab2(i)
!~             rtab2(i) = rtab2(j)
!~             rtab2(j) = rtmp

!~             i = i + 1
!~             j = j - 1
!~          endif
!~       enddo

!~       if (g<j) call sort_int_1int_1real(g, j, itabref, itab1, rtab2)
!~       if (d>i) call sort_int_1int_1real(i, d, itabref, itab1, rtab2)

!~    return
!~    endsubroutine sort_int_1int_1real

!~    recursive subroutine sort_int_1real(g, d, itabref, rtab1)
!~    implicit none
!~    integer(kind=I4), intent(in) :: g, d
!~    integer(kind=I4), dimension(:), intent(inout) :: itabref
!~    real(kind=R8), dimension(:), intent(inout)    :: rtab1
!~       integer(kind=I4) :: i, j, mil, cle, itmp
!~       real(kind=R8)    :: rtmp
!~       i = g
!~       j = d
!~       mil = (g+d)/2
!~       cle = itabref(mil)

!~       if (g>=d) return

!~       do while (i<=j)
!~          do while (itabref(i)<cle)
!~             i = i + 1
!~          enddo
!~          do while (itabref(j)>cle)
!~             j = j - 1
!~          enddo
!~          if (i<=j) then
!~             ! échange des éléments du tableau
!~             itmp       = itabref(i)
!~             itabref(i) = itabref(j)
!~             itabref(j) = itmp
!~             ! échange des éléments du vecteur 3
!~             rtmp     = rtab1(i)
!~             rtab1(i) = rtab1(j)
!~             rtab1(j) = rtmp

!~             i = i + 1
!~             j = j - 1
!~          endif
!~       enddo

!~       if (g<j) call sort_int_1real(g, j, itabref, rtab1)
!~       if (d>i) call sort_int_1real(i, d, itabref, rtab1)

!~    return
!~    endsubroutine sort_int_1real

endmodule sort_arrays
