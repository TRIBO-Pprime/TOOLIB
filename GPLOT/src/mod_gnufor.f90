module gnufor
implicit none

private

public :: write_xyy_data, write_xyy_plots, write_xy_data, write_xy_plot, run_gnuplot,  &  !
          test01, test02, test03, test04, test05, test06

contains

subroutine get_unit ( iunit )
!
!*******************************************************************************
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none
!
  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
endsubroutine get_unit
function pi ( )
!
!*******************************************************************************
!
!! PI returns the value of pi.
!
!
!  Modified:
!
!    04 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real PI, the value of pi.
!
  implicit none
!
  real pi
!
  pi = 3.14159265358979323846264338327950288419716939937510E+00

  return
endfunction pi
subroutine run_gnuplot ( command_file_name )
!
!*******************************************************************************
!
!! RUN_GNUPLOT runs GNUPLOT with a given command file.
!
!
!  Discussion:
!
!    The GNUPLOT program, version 3.7, must be available.  To check whether
!    this is so, try typing
!
!      which gnuplot
!
!    If the response is
!
!      gnuplot: command not found
!
!    then you're going to have to make GNUPLOT available.
!
!    At ISU, this may require that you issue the command
!
!      add gnu
!
!    You may need to set the environment variable GNUTERM:
!
!      setenv GNUTERM x11
!
!    so that GNUPLOT automatically displays to your X window terminal.
!
!
!    This routine expects that there is a text file containing the appropriate
!    commands to GNUPLOT to display your picture.  There are a number of
!    routines in this package that will do this for simple plotting tasks.
!    Most of them require that you also set up a file of data to be plotted.
!
!    Once this routine invokes GNUPLOT, a graphics window should open
!    up, and the FORTRAN program will pause.  Hitting RETURN should advance
!    to the next picture, or terminate the window at the end, allowing the
!    FORTRAN routine to proceed.
!
!
!    You can look at the data and command files created by the routines.
!    Moreover, you can easily modify the command file to change the options
!    used in GNUPLOT, and then run GNUPLOT interactively, as in:
!
!      gnuplot commands
!
!    In particular, if you want a PostScript version of your graphics files,
!    insert the command "set term postscript" at the beginning of the command
!    file and run gnuplot as follows:
!
!      gnuplot commands > mypicture.ps
!
!    You will also have to hit RETURN once for each plot that is made.
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILE_NAME, the name of the
!    command file.
!
  implicit none
!
  character ( len = 512 ) command
  character ( len = * ) command_file_name
  !integer status
!
!  Issue a command to the system that will startup GNUPLOT, using
!  the file we just wrote as input.
!
  write ( command, * ) 'gnuplot ' // '"' // trim ( command_file_name ) // '"'
!~   call execute_command_line(trim ( command ), wait=.false., exitstat=status)
  call system(trim ( command ))

  return
endsubroutine run_gnuplot
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,I4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
endsubroutine timestamp
subroutine write_vector_data ( data_file_name, n, x, y, dx, dy, ierror )
!
!*******************************************************************************
!
!! WRITE_VECTOR_DATA writes vector data to a file, for plotting by GNUPLOT.
!
!
!  Discussion:
!
!    Each vector is described by 4 values, X, Y, dX, dY, indicating that
!    a vector is to be drawn from (X,Y) to (X+dX,Y+dY).
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Input, integer N, the number of vectors.
!
!    Input, real X(N), Y(N), DX(N), DY(N), the vector data.
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  integer n
!
  character ( len = * ) data_file_name
  real dx(n)
  real dy(n)
  integer file_unit
  integer i
  integer ierror
  integer ios
  real x(n)
  real y(n)
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_VECTOR_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_VECTOR_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do i = 1, n
    write ( file_unit, * ) x(i), y(i), dx(i), dy(i)
  end do

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_VECTOR_DATA:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT vector data file "' // &
    trim ( data_file_name ) // '"'

  return
endsubroutine write_vector_data
subroutine write_vector_plot ( command_file_name, data_file_name, &
  ierror )
!
!*******************************************************************************
!
!! WRITE_VECTOR_PLOT writes GNUPLOT commands to plot vectors.
!
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILE_NAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  character ( len = * ) command_file_name
  character ( len = * ) data_file_name
  integer file_unit
  integer ierror
  integer ios
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_VECTOR_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_VECTOR_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a)' ) 'set style arrow'
  write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_file_name )
  write ( file_unit, '(a)' ) 'pause -1  "Hit return to continue"'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_VECTOR_PLOT:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT table plots command file "' // &
    trim ( command_file_name ) // '"'

  return
endsubroutine write_vector_plot
subroutine write_xy_data ( data_file_name, n, x, y, ierror )
!
!*******************************************************************************
!
!! WRITE_XY_DATA writes X(1:N), Y(1:N) data to a file.
!
!
!  Modified:
!
!    23 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Input, integer N, the number of data items.
!
!    Input, real X(N), Y(N), the X and Y data
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  integer n
!
  character ( len = * ) data_file_name
  integer file_unit
  integer i
  integer ierror
  integer ios
  real x(n)
  real y(n)
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do i = 1, n
    write ( file_unit, * ) x(i), y(i)
  end do

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XY_DATA:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XY data file "' // &
    trim ( data_file_name ) // '"'

  return
endsubroutine write_xy_data
subroutine write_xy2_data ( data_file_name, n, x, y, y2, ierror )
!
!*******************************************************************************
  implicit none
!
  integer n
!
  character ( len = * ) data_file_name
  integer file_unit
  integer i
  integer ierror
  integer ios
  real x(n)
  real y(n)
  real y2(n)
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do i = 1, n
    write ( file_unit, * ) x(i), y(i), y2(i)
  end do

  close ( unit = file_unit )

  return
endsubroutine write_xy2_data
subroutine write_xy_plot ( command_file_name, data_file_name, logscale, ierror, nom )
!
!*******************************************************************************
!
!! WRITE_XY_PLOT writes GNUPLOT commands to plot X(1:N), Y(1:N) data.
!
!
!  Modified:
!
!    23 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILE_NAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  character ( len = * ) command_file_name
  character ( len = * ) nom
  character ( len = * ) data_file_name
  integer file_unit
  !integer i
  integer ierror
  integer ios
  logical logscale
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  if (trim(nom)/='*') then
      write ( file_unit, '(a)' ) 'set term pngcairo'
      write ( file_unit, '(a)' ) 'set output "'//trim(nom)//'.png"'
  endif



  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  if (logscale) write ( file_unit, '(a)' ) 'set logscale x'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_file_name ) // '" using 1:2 with lines'
  !write ( file_unit, '(a)' ) 'pause -1  "Hit return to continue"'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XY_PLOT:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XY plot command file "' // trim ( command_file_name ) // '"'

  return
endsubroutine write_xy_plot
subroutine write_xy2_plot ( command_file_name, data_file_name, logscale, ierror, nom )
  implicit none
!
  character ( len = * ) command_file_name
  character ( len = * ) nom
  character ( len = * ) data_file_name
  integer file_unit
  !integer i
  integer ierror
  integer ios
  logical logscale
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  if (trim(nom)/='*') then
      write ( file_unit, '(a)' ) 'set term pngcairo'
      write ( file_unit, '(a)' ) 'set output "'//trim(nom)//'.png"'
  endif



  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  if (logscale) write ( file_unit, '(a)' ) 'set logscale x'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_file_name ) // '" using 1:2 with lines,\'
  write ( file_unit, '(a,i2,a)' )       '"'// trim ( data_file_name ) // '" using 1:3 with lines '
  !write ( file_unit, '(a)' ) 'pause -1  "Hit return to continue"'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  return
endsubroutine write_xy2_plot
subroutine write_y_plot ( command_file_name, data_file_name, ierror )
!
!*******************************************************************************
!
!! WRITE_Y_PLOT writes GNUPLOT commands to plot Y(1:N) data.
!
!
!  Modified:
!
!    23 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILE_NAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  character ( len = * ) command_file_name
  character ( len = * ) data_file_name
  integer file_unit
  !integer i
  integer ierror
  integer ios
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_file_name ) // &
    '" using 1 with lines'
  write ( file_unit, '(a)' ) 'pause -1  "Hit return to continue"'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XY_PLOT:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XY plot command file "' // &
    trim ( command_file_name ) // '"'

  return
endsubroutine write_y_plot
subroutine write_xyy_data ( data_file_name, lda, nrow, ncol, x, ierror )
!
!*******************************************************************************
!
!! WRITE_XYY_DATA writes a table of data to a file, for plotting by GNUPLOT.
!
!
!  Discussion:
!
!    The first column of data is assumed to be the independent variable, X.
!    Separate plots are made of X versus all the other columns of data.
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Input, integer LDA, the leading dimension of X.
!
!    Input, integer NROW, NCOL, the dimensions of X.
!
!    Input, real X(LDA,NCOL), the NROW by NCOL data to be written.
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  integer lda
  integer ncol
!
  character ( len = * ) data_file_name
  integer file_unit
  integer i
  integer ierror
  integer ios
  integer nrow
  real x(lda,ncol)
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYY_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYY_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do i = 1, nrow
    write ( file_unit, * ) x(i,1:ncol)
  end do

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYY_DATA:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XYY data file "' // &
    trim ( data_file_name ) // '"'

  return
endsubroutine write_xyy_data
subroutine write_xyy_plots ( command_file_name, data_file_name, add_lines, title, &
  ncol, ierror )
!
!*******************************************************************************
!
!! WRITE_XYY_PLOTS writes GNUPLOT commands to make multiple (X,Y) plots.
!
!
!  Discussion:
!
!    The first column of data is assumed to be the independent variable, X.
!    Separate plots are made of X versus all the other columns of data.
!
!  Modified:
!
!    23 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILE_NAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Input, integer NCOL, the number of columns of data.
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  character ( len = * ) command_file_name
  character ( len = * ) data_file_name
  character ( len = * ), optional :: add_lines
  integer file_unit
  integer i
  integer ierror
  integer ios
  integer ncol
  character ( len = * ), dimension(1:ncol) :: title

!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYY_PLOTS - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYY_PLOTS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'

  if (present(add_lines)) then
   write ( file_unit, '(a)' ) trim(add_lines)
  endif

  do i = 1, 1
    write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_file_name ) // '" using 1:', i+1, ' with lines title "' // trim(title(i)) // '", \'
  end do
  do i = 2, ncol-2
    write ( file_unit, '(a,i2,a)' )       '"'// trim ( data_file_name ) // '" using 1:', i+1, ' with lines title "' // trim(title(i)) // '", \'
  end do
  do i = ncol-1, ncol-1
    write ( file_unit, '(a,i2,a)' )       '"'// trim ( data_file_name ) // '" using 1:', i+1, ' with lines title "' // trim(title(i)) // '"'
  end do
  write ( file_unit, '(a)' ) 'pause -1  "Hit return to continue"'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYY_PLOTS:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XYY plots command file "' // &
    trim ( command_file_name ) // '"'

  return
endsubroutine write_xyy_plots
subroutine write_xyz_data ( data_file_name, n, x, y, z, ierror )
!
!*******************************************************************************
!
!! WRITE_XYZ_DATA writes X(1:N), Y(1:N), Z(1:N) data to a file.
!
!
!  Modified:
!
!    23 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Input, integer N, the number of data items.
!
!    Input, real X(N), Y(N), Z(N), the X, Y, Z data
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  integer n
!
  character ( len = * ) data_file_name
  integer file_unit
  integer i
  integer ierror
  integer ios
  real x(n)
  real y(n)
  real z(n)
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZ_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZ_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do i = 1, n
    write ( file_unit, * ) x(i), y(i), z(i)
  end do

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYZ_DATA:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XYZ data file "' // &
    trim ( data_file_name ) // '"'

  return
endsubroutine write_xyz_data
subroutine write_xyz_plot ( command_file_name, data_file_name, ierror )
!
!*******************************************************************************
!
!! WRITE_XYZ_PLOT writes commands to plot parametric (X,Y,Z) data.
!
!
!  Discussion:
!
!    This routine tries to write a command file suitable for displaying
!    a 3D arc specified by points (X,Y,Z).  A grid data file, containing
!    values of X, Y and Z, will also be needed.
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILE_NAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  character ( len = * ) command_file_name
  character ( len = * ) data_file_name
  integer file_unit
  !integer i
  integer ierror
  integer ios
  !integer ncol
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZ_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZ_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a)' ) 'set parametric'
  write ( file_unit, '(a)' ) 'splot "' // trim ( data_file_name ) // &
      '" using 1:2:3 with lines'
  write ( file_unit, '(a)' ) 'pause -1  "Hit return to continue"'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYZ_PLOT:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT SPLOT command file "' // &
    trim ( command_file_name ) // '"'

  return
endsubroutine write_xyz_plot
subroutine write_xyzgrid_contour ( command_file_name, data_file_name, ierror )
!
!*******************************************************************************
!
!! WRITE_XYZGRID_CONTOUR writes commands to plot contours of Z(X,Y).
!
!
!  Discussion:
!
!    This routine tries to write a command file suitable for displaying
!    contours of Z(X,Y) gridded data.  A grid data file, containing values
!    of X, Y and Z, will also be needed.
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILE_NAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  character ( len = * ) command_file_name
  character ( len = * ) data_file_name
  integer file_unit
  !integer i
  integer ierror
  integer ios
  !integer ncol
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a)' ) 'set parametric'
  write ( file_unit, '(a)' ) 'set nosurface'
  write ( file_unit, '(a)' ) 'set contour'
  write ( file_unit, '(a)' ) 'set cntrparam levels 10'
  !write ( file_unit, '(a)' ) 'set terminal table'
  !write ( file_unit, '(a)' ) 'set out "table.txt"'
  write ( file_unit, '(a)' ) 'splot "' // trim ( data_file_name ) // &
    '" using 1:2:3 with lines'
  write ( file_unit, '(a)' ) 'set term x11'
  !write ( file_unit, '(a)' ) 'plot "table.txt" with lines'
  write ( file_unit, '(a)' ) 'pause -1  "Hit return to continue"'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XYZGRID contour plot command file "' // &
    trim ( command_file_name ) // '"'

  return
endsubroutine write_xyzgrid_contour
subroutine write_xyzgrid_data ( data_file_name, nx, ny, xyz, ierror )
!
!*******************************************************************************
!
!! WRITE_XYZGRID_DATA writes a file of XYZ grid data.
!
!
!  Discussion:
!
!    It is assumed that values of Z are available on a regular NX by NY grid
!    of (X,Y) points.
!
!    The form of the data file requires that all the data for a given value
!    of Y be listed, followed by a blank line, followed by the data for
!    another value of Y.
!
!  Example:
!
!    Here is a grid data file for a 3 by 3 grid, with Z = X + Y.
!
!    0.0 0.0 0.0
!    1.0 0.0 1.0
!    2.0 0.0 2.0
!
!    0.0 1.0 1.0
!    1.0 1.0 2.0
!    2.0 1.0 3.0
!
!    0.0 2.0 2.0
!    1.0 2.0 3.0
!    2.0 2.0 4.0
!
!  Modified:
!
!    23 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Input, integer NX, NY, the dimensions of the grid.
!
!    Input, real XYZ(3,NX,NY), the XYZ grid data to be written.
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  integer nx
  integer ny
!
  character ( len = * ) data_file_name
  integer file_unit
  integer i
  integer ierror
  integer ios
  integer j
  real xyz(3,nx,ny)
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do j = 1, ny
    do i = 1, nx
      write ( file_unit, * ) xyz(1:3,i,j)
    end do
    write ( file_unit, '(a)' )
  end do

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYZGRID_DATA:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XYZ grid data file "' // &
    trim ( data_file_name ) // '"'

  return
endsubroutine write_xyzgrid_data
subroutine write_xyzgrid_surface ( command_file_name, data_file_name, ierror )
!
!*******************************************************************************
!
!! WRITE_XYZGRID_SURFACE writes a file of GNUPLOT commands to plot a 3D surface.
!
!
!  Discussion:
!
!    This routine tries to write a command file suitable for displaying
!    a surface Z(X,Y).  A grid data file, containing values of X, Y and Z,
!    will also be needed.  The routine WRITE_XYZGRID_DATA can write this file.
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILE_NAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILE_NAME, the name of the data file.
!
!    Output, integer IERROR, nonzero if an error occurred.
!
  implicit none
!
  character ( len = * ) command_file_name
  character ( len = * ) data_file_name
  integer file_unit
  !integer i
  integer ierror
  integer ios
  !integer ncol
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_SURFACE - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_SURFACE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a)' ) 'set parametric'
  write ( file_unit, '(a)' ) 'set pm3d implicit at s'!; set palette gray'
  write ( file_unit, '(a)' ) 'set hidden3d'
  !write ( file_unit, '(a)' ) 'set contour'
!~   write ( file_unit, '(a)' ) 'splot "' // trim ( data_file_name ) // &
!~     '" using 1:2:3 with lines'
  write ( file_unit, '(a)' ) 'splot "' // trim ( data_file_name ) // '"'
  write ( file_unit, '(a)' ) 'pause -1  "Hit return to continue"'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_SURFACE_COMMANDS:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT surface plot command file "' // &
    trim ( command_file_name ) // '"'

  return
endsubroutine write_xyzgrid_surface
!**************************************************************************************************
!**************************************************************************************************
subroutine test01
!
!*******************************************************************************
!
!! TEST01 demonstrates the plotting of Y(X) data.
!
  implicit none
!
  integer, parameter :: n = 101
!
  !real angle
  !real area
  character ( len = 100 ) :: command_file_name = 'tmp/test01_commands.txt'
  character ( len = 100 ) :: data_file_name = 'tmp/test01_data.txt'
  integer i
  integer ierror
  real x(n)
  real, parameter :: xmin = 0.0+00
  real, parameter :: xmax = 20.0E+00
  real y(n)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  To plot a simple set of (X,Y) data,'
  write ( *, '(a)' ) '  WRITE_XY_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_XY_PLOT writes the plot command file.'

  do i = 1, n
    x(i) = ( real ( n - i ) * xmin + real ( i - 1 ) * xmax ) / real ( n - 1 )
    y(i) = sin ( x(i) ) * sin ( 4.0 * x(i) )
  end do

  call write_xy_data ( data_file_name, n, x, y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01'
    write ( *, '(a,i6)' ) '  WRITE_XY_DATA returned IERROR = ', ierror
  end if

  call write_xy_plot ( command_file_name, data_file_name, .true., ierror, '*' )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01'
    write ( *, '(a,i6)' ) '  WRITE_XY_PLOT returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_file_name )

  return
endsubroutine test01
subroutine test02
!
!*******************************************************************************
!
!! TEST02 demonstrates the plotting of a table of data.
!
  implicit none
!
  integer, parameter :: nrow = 101
  integer, parameter :: ncol = 4
!
  integer, parameter :: lda = nrow
!
  real angle
  real area
  character ( len = 100 ) :: command_file_name = 'tmp/test02_commands.txt'
  character ( len = 100 ) :: data_file_name = 'tmp/test02_data.txt'
  character ( len = 100 ), dimension(1:ncol) :: title
  real height
  integer i
  integer ierror
  real, parameter :: r = 50.0E+00
  real width
  real x(lda,ncol)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  To plot X versus multiple sets of Y data,'
  write ( *, '(a)' ) '  WRITE_XYY_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_XYY_PLOT writes the plot command file.'

  do i = 1, nrow

    height = 2.0E+00 * r * real ( i - 1 ) / real ( nrow - 1 )
    width = 2.0E+00 * sqrt ( r**2 - ( r - height )**2 )
    angle = acos ( ( r - height ) / r )
    area = 0.5E+00 * r**2 * 2.0E+00 * acos ( ( r - height ) / r ) &
      - ( r - height ) * sqrt ( height * ( 2.0E+00 * r - height ) )

    x(i,1) = height
    x(i,2) = width
    x(i,3) = angle
    x(i,4) = area

  end do

  title(1) = 'height'
  title(2) = 'width'
  title(3) = 'angle'
  title(4) = 'area'

  call write_xyy_data ( data_file_name, lda, nrow, ncol, x, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02'
    write ( *, '(a,i6)' ) '  WRITE_XYY_DATA returned IERROR = ', ierror
  end if

  call write_xyy_plots ( command_file_name=command_file_name, data_file_name=data_file_name, title=title, ncol=ncol, &
    ierror=ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02'
    write ( *, '(a,i6)' ) '  WRITE_XYY_PLOTS returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_file_name )

  return
endsubroutine test02
subroutine test03
!
!*******************************************************************************
!
!! TEST03 plots parameter (X,Y,Z) data.
!
  implicit none
!
  integer, parameter :: n = 101
!
  character ( len = 100 ) :: command_file_name = 'tmp/test03_commands.txt'
  character ( len = 100 ) :: data_file_name = 'tmp/test03_data.txt'
  integer i
  integer ierror
  integer, parameter :: nturn = 5
  real, parameter :: r = 5.0E+00
  real theta
  real x(n)
  real y(n)
  real z(n)
  real, parameter :: zmax = 10.0E+00
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  To plot a (parametric) set of (X,Y,Z) data,'
  write ( *, '(a)' ) '  WRITE_XYZ_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_XYZ_PLOT writes the plot command file.'

  do i = 1, n

    z(i) = zmax * real ( i - 1 ) / real ( n - 1 )

    theta = ( 2.0E+00 * pi() ) * z(i) * real ( nturn ) / zmax

    x(i) = r * cos ( theta )
    y(i) = r * sin ( theta )

  end do

  call write_xyz_data ( data_file_name, n, x, y, z, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03'
    write ( *, '(a,i6)' ) '  WRITE_XYZ_DATA returned IERROR = ', ierror
  end if

  call write_xyz_plot ( command_file_name, data_file_name, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03'
    write ( *, '(a,i6)' ) '  WRITE_XYZ_PLOT returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_file_name )

  return
endsubroutine test03
subroutine test04
!
!*******************************************************************************
!
!! TEST04 plots vector data.
!
  implicit none
!
  integer, parameter :: nx = 21
  integer, parameter :: ny = 21
  integer, parameter :: n = nx * ny
!
  character ( len = 100 ) :: command_file_name = 'tmp/test04_commands.txt'
  character ( len = 100 ) :: data_file_name = 'tmp/test04_data.txt'
  real dx(n)
  real dy(n)
  integer i
  integer ierror
  integer j
  integer k
  real x(n)
  real, parameter :: xmax = 1.0E+00
  real, parameter :: xmin = -1.0E+00
  real xx
  real y(n)
  real, parameter :: ymax = 1.0E+00
  real, parameter :: ymin = -1.0E+00
  real yy
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  To plot a vector field,'
  write ( *, '(a)' ) '  WRITE_VECTOR_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_VECTOR_PLOT writes the plot command file.'

  k = 0

  do i = 1, nx

    do j = 1, ny

      k = k + 1

      xx = ( real ( nx - i ) * xmin + real ( i - 1 ) * xmax ) / real ( nx - 1 )
      yy = ( real ( ny - j ) * ymin + real ( j - 1 ) * ymax ) / real ( ny - 1 )

      dx(k) = - 0.10E+00 * yy
      dy(k) =   0.10E+00 * xx

      x(k) = xx - 0.5E+00 * dx(k)
      y(k) = yy - 0.5E+00 * dy(k)

    end do

  end do

  call write_vector_data ( data_file_name, n, x, y, dx, dy, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04'
    write ( *, '(a,i6)' ) '  WRITE_VECTOR_DATA returned IERROR = ', ierror
  end if

  call write_vector_plot ( command_file_name, data_file_name, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04'
    write ( *, '(a,i6)' ) '  WRITE_VECTOR_PLOT returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_file_name )

  return
endsubroutine test04
subroutine test05
!
!*******************************************************************************
!
!! TEST05 plots Z(X,Y) grid data as a surface.
!
  implicit none
!
  integer, parameter :: nx = 21
  integer, parameter :: ny = 21

  !integer, parameter :: nrow = nx * ny
!
  character ( len = 100 ) :: command_file_name = 'tmp/test05_commands.txt'
  character ( len = 100 ) :: data_file_name = 'tmp/test05_data.txt'
  integer i
  integer ierror
  integer j
  real x
  real, parameter :: xmax = 1.0E+00
  real, parameter :: xmin = 0.0E+00
  real xyz(3,nx,ny)
  real y
  real, parameter :: ymax = 1.0E+00
  real, parameter :: ymin = 0.0E+00
  real z
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  To plot a gridded set of Z(X,Y) data as a surface,'
  write ( *, '(a)' ) '  WRITE_XYZGRID_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_XYZGRID_SURFACE writes the plot command file.'

  do i = 1, nx

    x = ( real ( nx - i ) * xmin + real ( i - 1 ) * xmax ) / real ( nx - 1 )

    do j = 1, ny

      y = ( real ( ny - j ) * ymin + real ( j - 1 ) * ymax ) / real ( ny - 1 )

      z = sin ( 64.0E+00 * ( x - 0.5E+00 )**2 * ( y - 0.5E+00 )**2 )

      xyz(1,i,j) = x
      xyz(2,i,j) = y
      xyz(3,i,j) = z

    end do

  end do

  call write_xyzgrid_data ( data_file_name, nx, ny, xyz, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05'
    write ( *, '(a,i6)' ) '  WRITE_XYZGRID_DATA returned IERROR = ', ierror
  end if

  call write_xyzgrid_surface ( command_file_name, data_file_name, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05'
    write ( *, '(a,i6)' ) '  WRITE_XYZGRID_SURFACE returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_file_name )

  return
endsubroutine test05
subroutine test06
!
!*******************************************************************************
!
!! TEST06 plots Z(X,Y) grid data as contours.
!
  implicit none
!
  integer, parameter :: nx = 41
  integer, parameter :: ny = 41
!
  character ( len = 100 ) :: command_file_name = 'tmp/test06_commands.txt'
  character ( len = 100 ) :: data_file_name = 'tmp/test06_data.txt'
  integer i
  integer ierror
  integer j
  real x
  real, parameter :: xmax = 1.0E+00
  real, parameter :: xmin = 0.0E+00
  real xyz(3,nx,ny)
  real y
  real, parameter :: ymax = 1.0E+00
  real, parameter :: ymin = 0.0E+00
  real z
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  To plot gridded Z(X,Y) data as contours,'
  write ( *, '(a)' ) '  WRITE_XYZGRID_DATA writes the data file,'
  write ( *, '(a)' ) '  WRITE_XYZGRID_CONTOUR writes the plot command file.'

  do i = 1, nx

    x = ( real ( nx - i ) * xmin + real ( i - 1 ) * xmax ) / real ( nx - 1 )

    do j = 1, ny

      y = ( real ( ny - j ) * ymin + real ( j - 1 ) * ymax ) / real ( ny - 1 )

      z = sin ( 64.0E+00 * ( x - 0.5E+00 )**2 * ( y - 0.5E+00 )**2 )

      xyz(1,i,j) = x
      xyz(2,i,j) = y
      xyz(3,i,j) = z

    end do

  end do

  call write_xyzgrid_data ( data_file_name, nx, ny, xyz, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06'
    write ( *, '(a,i6)' ) '  WRITE_XYZGRID_DATA returned IERROR = ', ierror
  end if

  call write_xyzgrid_contour ( command_file_name, data_file_name, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06'
    write ( *, '(a,i6)' ) '  WRITE_XYZGRID_CONTOUR returned IERROR = ', ierror
  end if

  call run_gnuplot ( command_file_name )

  return
endsubroutine test06













subroutine write_polar_plot ( command_file_name, data_file_name, ierror )
  implicit none
!
  character ( len = * ) command_file_name
  character ( len = * ) data_file_name
  integer file_unit
  integer ierror
  integer ios
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_POLAR_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_POLAR_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set angles degrees'
  write ( file_unit, '(a)' ) 'set polar'
  write ( file_unit, '(a)' ) 'set grid polar'
  write ( file_unit, '(a)' ) 'set xlabel "azimuth"'
  write ( file_unit, '(a)' ) 'set ylabel "r"'
  write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_file_name ) // &
    '" using 1:2 with lines'
  write ( file_unit, '(a)' ) 'pause -1  "Hit return to continue"'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XY_PLOT:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XY plot command file "' // &
    trim ( command_file_name ) // '"'

  return
endsubroutine write_polar_plot


subroutine write_polar_data ( data_file_name, n, x, y, ierror )
  implicit none
!
  integer n
!
  character ( len = * ) data_file_name
  integer file_unit
  integer i
  integer ierror
  integer ios
  real x(n)
  real y(n)
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_POLAR_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_POLAR_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do i = 1, n
    write ( file_unit, * ) x(i), y(i)
  end do

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_POLAR_DATA:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT POLAR data file "' // &
    trim ( data_file_name ) // '"'

  return
endsubroutine write_polar_data



endmodule gnufor
