!< author: Arthur Francisco
!<  version: 1.0.1
!<  date: feb, 24 2023
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<     **Various routines to sort real/integer arrays. Example of use.**
!<  </span>

program main
use data_arch,     only : I4, R8
use miscellaneous, only : get_unit
use sort_arrays,   only : sort_array2, init_order
implicit none

integer(kind=I4), parameter :: n = 100

real(kind=R8),    dimension(n) :: vec_real1,     vec_real2
real(kind=R8),    dimension(n) :: vec_sol_real1, vec_sol_real2
integer(kind=I4), dimension(n) :: vec_int1,      vec_int2,     vec_order
integer(kind=I4), dimension(n) :: vec_sol_int1,  vec_sol_int2, vec_sol_order

integer(kind=I4) :: i, k, u10, u20
real(kind=R8)    :: r

!=============== Sort an integer vector, and other vectors with the same order

! ------------- get non ordered vectors
call get_unit(u10) ; open( unit = u10, file = 'try/non_ordered_vec.csv' )

   read(u10, *) ! column name
   do i = 1, n
      read(u10, *) vec_int1(i), vec_real1(i), &  !
                   vec_int2(i), vec_real2(i), k  !
   enddo

close(u10)

! ------------- sort with respect to the first column (integer values)
call sort_array2(tab_inout = vec_int1(1:n),           &  !
                      tab1 = vec_real1(1:n),          &  !
                      tab2 = vec_int2(1:n),           &  !
                      tab3 = vec_real2(1:n), n = n)      !

! ------------- get solution
call get_unit(u10) ; open( unit = u10, file = 'try/ordered_vec_col1_int.csv' )

   read(u10, *) ! column name
   do i = 1, n
      read(u10, *) vec_sol_int1(i), vec_sol_real1(i), &  !
                   vec_sol_int2(i), vec_sol_real2(i), k  !
   enddo

close(u10)

! ------------- compare
call get_unit(u20) ; open( unit = u20, file = 'out/res_ordered_vec_col1_int.csv' )

   do i = 1, n
      write(u20,*) vec_int1(i) - vec_sol_int1(i), vec_real1(i) - vec_sol_real1(i), &  !
                   vec_int2(i) - vec_sol_int2(i), vec_real2(i) - vec_sol_real2(i)
   enddo
   !----------------------------- REMARK --------------------------------------
   ! As the sorted vector is integers, some values may be identical, leading to
   ! several possibilities for the subsequently sorted vectors. Look elements
   ! 9 and 10 for example.
   !---------------------------------------------------------------------------

close(u20)

!=============== Sort a real vector, and other vectors with the same order

! ------------- get non ordered vectors
call get_unit(u10) ; open( unit = u10, file = 'try/non_ordered_vec.csv' )

   read(u10, *) ! column name
   do i = 1, n
      read(u10, *) vec_int1(i), vec_real1(i), &  !
                   vec_int2(i), vec_real2(i), k  !
   enddo

close(u10)


! ------------- sort with respect to the second column (real values)
call sort_array2(tab_inout = vec_real1(1:n),          &  !
                      tab1 = vec_int1(1:n),           &  !
                      tab2 = vec_int2(1:n),           &  !
                      tab3 = vec_real2(1:n), n = n)      !


! ------------- get solution
call get_unit(u10) ; open( unit = u10, file = 'try/ordered_vec_col2_real.csv' )

   read(u10, *) ! column name
   do i = 1, n
      read(u10, *) vec_sol_int1(i), vec_sol_real1(i), &  !
                   vec_sol_int2(i), vec_sol_real2(i), k  !
   enddo

close(u10)

! ------------- compare
call get_unit(u20) ; open( unit = u20, file = 'out/res_ordered_vec_col2_real.csv' )

   do i = 1, n
      write(u20,*) vec_int1(i) - vec_sol_int1(i), vec_real1(i) - vec_sol_real1(i), &  !
                   vec_int2(i) - vec_sol_int2(i), vec_real2(i) - vec_sol_real2(i)     !
   enddo

close(u20)

!=============== Sort a real vector and return the order vector

! ------------- get non ordered vectors
call get_unit(u10) ; open( unit = u10, file = 'try/non_ordered_vec.csv' )

   read(u10, *) ! column name
   do i = 1, n
      read(u10, *) k, vec_real1(i), k, r, k
   enddo

close(u10)


! ------------- sort with respect to the second column (real values, same case as above)
call init_order(order = vec_order(1:n), n = n)

call sort_array2(tab_inout = vec_real1(1:n),          &  !
                      tab0 = vec_order(1:n), n = n)      !

! ------------- get solution of order
call get_unit(u10) ; open( unit = u10, file = 'try/ordered_vec_col2_real.csv' )

   read(u10, *) ! column name
   do i = 1, n
      read(u10, *) k, r, k, r, vec_sol_order(i)
   enddo

close(u10)

! ------------- compare
call get_unit(u20) ; open( unit = u20, file = 'out/res_order_vec.csv' )

   do i = 1, n
      write(u20,*) vec_order(i) - vec_sol_order(i)
   enddo

close(u20)

stop

endprogram main
