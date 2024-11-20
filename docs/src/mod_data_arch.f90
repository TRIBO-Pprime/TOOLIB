!< author: Arthur Francisco
!<  version: 1.1.0
!<  date: april, 6 2023
!
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<     **General parameters definition**
!<  </span>
module data_arch
use, intrinsic :: iso_fortran_env, only : output_unit, input_unit, error_unit,   &  !
                                          int32, int64,                          &  !
                                          real32, real64                            !
implicit none
public ! don't keep global iso_fortran_env parameters

   integer(kind=int32), parameter :: I4 = int32

   integer(kind=I4), parameter :: I8 = int64
   integer(kind=I4), parameter :: R4 = real32
   integer(kind=I4), parameter :: R8 = real64
   integer(kind=I4), parameter :: HIG_I4 = huge(1)

   integer(kind=I4), parameter :: OPU = output_unit   !! *Output unit*
   integer(kind=I4), parameter :: IPU = input_unit    !! *Input unit*
   integer(kind=I4), parameter :: ERU = error_unit    !! *Error unit*

   real(kind=R8), parameter :: UN = 1.0_R8

   real(kind=R4), parameter :: PI_R4  = acos(-1._R4)
   real(kind=R8), parameter :: PI_R8  = acos(-1._R8)
   real(kind=R4), parameter :: EPS_R4 = tiny(1._R4)
   real(kind=R8), parameter :: EPS_R8 = tiny(1._R8)
   real(kind=R8), parameter :: HIG_R8 = huge(1._R8)
   real(kind=R8), parameter :: HIG_E8 = log(HIG_R8)
   real(kind=R8), parameter :: EPS_E8 = log(EPS_R8)

   integer(kind=I4), parameter :: EXPO_MAX = exponent(HIG_R8)

!~    integer(kind=I4) :: NB_THREADS_MAX

private :: output_unit, input_unit, error_unit, int32, int64, real32, real64

endmodule data_arch

