
!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: february, 27 2023
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<     **Some routines to deal with files**
!<  </span>
!<
!< @todo
!< rm -rf dir-name
!< @endtodo

module files
use data_arch,     only : I4
use miscellaneous, only : get_unit
implicit none

private

public :: is_linux, dir_separator, mkdir, str_replace, list_files, list_dirs, clean_scratch, make_path,  &  !
          path2vec, vec2path, filename, dirname, str_remove_chars

contains

   !================================================================================================
   subroutine clean_scratch()
   !! Subroutine that removes all files with extension .scratch
   implicit none

      call execute_command_line( "find -type f -name ""*.scratch"" | xargs rm > err.scratch" )

   return
   endsubroutine clean_scratch


   !================================================================================================
   subroutine list_files(dir, list, ext)
   !! Subroutine that returns a list of files in a directory
   implicit none
   character(len =   *), intent(in )                            :: dir
   character(len = 512), intent(out), allocatable, dimension(:) :: list
   character(len =   *), intent(in ), optional                  :: ext

      character(len = :), allocatable :: out_file

      character(len = 512) :: string
      integer(kind = I4) :: ui, line, i, max_len, ierr

      if ( present(ext) ) then
         out_file = trim(ext) // ".scratch"
         call execute_command_line( "find " // trim(dir) // " -maxdepth 1 -type f -name ""*." // trim(ext) // """ > " // out_file )
      else
         out_file = "files.scratch"
         call execute_command_line("find " // trim(dir) // " -maxdepth 1 -type f > " // out_file)
      endif

      call get_unit(ui)

      open( unit = ui, file = out_file )

         line = 0
         max_len = 0
         do
            string = ' '
            read( unit = ui, fmt = '(a)', iostat = ierr ) string
            if ( ierr /= 0 ) exit
            max_len = max( max_len, len_trim( string ) )
            line = line + 1
         enddo
         if ( max_len > 512 ) stop 'Error list_files: path too long'

         allocate( list(1:line) )

      rewind( ui )

         do i = 1, line
            read( unit = ui, fmt = '(a)' ) list( i )
         enddo

      close( ui )

      deallocate( out_file )

   return
   endsubroutine list_files


   !================================================================================================
   subroutine list_dirs(str)
   !! Subroutine that returns a list of subdirectories
   implicit none
   character(len = *), intent(in), optional :: str

      if ( present(str) ) then
         call execute_command_line("find . -maxdepth 1 -type d -name ""*" // trim(str) // "*"" > " // trim(str) // ".scratch")
      else
         call execute_command_line("find . -maxdepth 1 -type d > dirs.scratch")
      endif

   return
   endsubroutine list_dirs


   !================================================================================================
   logical(kind = 4) function is_linux()
   !! Function that returns true if the operating system is linux
   implicit none

      character(len = 512) :: path

      call get_environment_variable("PATH", path)
      is_linux = ( path(1:1) == "/" )

   return
   endfunction is_linux


   !================================================================================================
   character(len = 1) function dir_separator()
   !! Function that returns the system directory separator
   implicit none

      logical(kind = I4)  :: os_linux, os_windows

      os_linux   = is_linux()
      os_windows = .not. os_linux

      dir_separator = "/"
      if ( os_windows ) dir_separator = "\"

   return
   endfunction dir_separator


   !================================================================================================
   subroutine mkdir(wkd, directory, sep, exit_status)
   !! Subroutine that creates a directory
   implicit none
   character(len = *), intent(in ) :: wkd
   character(len = *), intent(in ) :: directory
   character(len = 1), intent(in ) :: sep
   integer(kind = I4), intent(out) :: exit_status

      character(len = 512) :: dir

      dir = trim(wkd) // sep // trim(directory)

      call execute_command_line("mkdir -p " // sep // dir, exitstat = exit_status)

   return
   endsubroutine mkdir


   !================================================================================================
   subroutine make_path(wkd, file_path, exit_status)
   !! Subroutine that creates the folders of a file path
   implicit none
   character(len = *), intent(in ) :: wkd
   character(len = *), intent(in ) :: file_path
   integer(kind = I4), intent(out) :: exit_status

      character(len = 512) :: dir
      character(len =   1) :: sep
      integer(kind = I4)   :: isep

      sep = dir_separator()

      isep = index( file_path, sep, back = .true. )

      dir = file_path(1:isep - 1)

      call mkdir(wkd = wkd, directory = trim(dir), sep = sep, exit_status = exit_status)

   return
   endsubroutine make_path


   !================================================================================================
   subroutine path2vec(file_path, vec_path)
   !! Subroutine that creates a vector containing the folders of a file path
   implicit none
   character(len = *),   intent(in ) :: file_path
   character(len = 512), intent(out), dimension(:), allocatable :: vec_path

      character(len = 512) :: dir
      character(len =   1) :: sep
      integer(kind = I4)   :: is1, is2, k

      sep = dir_separator()

      k   = 0
      is2 = len_trim( file_path )
      do

         is2 = index( file_path(1:is2 - 1), sep, back = .true. )

         if ( is2 > 0 ) then

            k = k + 1

         else

            exit

         endif

      enddo

      allocate( vec_path(1:k - 1) )

      is2 = index( file_path, sep, back = .true. )
      dir = file_path(1:is2 - 1)

      k = 0
      do

         is1 = index( file_path(1:is2 - 1), sep, back = .true. )

         if ( is1 > 0 ) then

            k = k + 1
            vec_path(k) = file_path( is1 + 1:is2 - 1 )

            is2 = is1

         else

            exit

         endif

      enddo

      vec_path = vec_path( size(vec_path):1:-1 )

   return
   endsubroutine path2vec


   !================================================================================================
   subroutine vec2path(file_path, vec_path)
   !! Subroutine that creates a path from vector of folders
   implicit none
   character(len = : ),  intent(out), allocatable  :: file_path
   character(len = 512), intent(in ), dimension(:) :: vec_path

      character(len =   1) :: sep
      integer(kind = I4)   :: k

      sep = dir_separator()

      file_path = ''
      do k = 1, size(vec_path)

         file_path = file_path // sep // trim( vec_path(k) )

      enddo

   return
   endsubroutine vec2path


   !================================================================================================
   function filename(file_path)
   !! Subroutine that keeps only the file from a path
   implicit none
   character(len = * ), intent(in) :: file_path
   character(len = :), allocatable :: filename

      character(len = 1) :: sep
      integer(kind = I4) :: ind

      sep = dir_separator()

      ind = index( file_path, sep, back = .true. )

      filename = trim( file_path( ind + 1: ) )

   return
   endfunction filename


   !================================================================================================
   function dirname(file_path)
   !! Subroutine that keeps only the directory from a file path
   implicit none
   character(len = * ), intent(in) :: file_path
   character(len = :), allocatable :: dirname

      character(len = 1) :: sep
      integer(kind = I4) :: ind

      sep = dir_separator()

      ind = index( file_path, sep, back = .true. )

      dirname = trim( file_path( :ind - 1 ) )

   return
   endfunction dirname


   !================================================================================================
   function str_replace(string, old_str, new_str , place)
   !! Function that replaces a string with another string
   implicit none
   character(len = :), allocatable :: str_replace  !! returned string
   character(len = *), intent(in) :: string        !! string to be modified
   character(len = *), intent(in) :: old_str
   character(len = *), intent(in) :: new_str
   integer(kind = I4), intent(in) :: place

      integer(kind = I4) :: ind

      select case ( place )
         case (0, 2)
            ind = index( string, old_str )
         case (1)
            ind = index( string, old_str, back = .true. )
         case default
            stop 'str_replace bad choice'
      endselect

      str_replace = string(1:ind - 1) // new_str // string(ind + len(old_str):len(string))

      if ( place == 2) then
         do
            ind = index( str_replace, old_str )
            if ( ind == 0 ) exit
            str_replace = str_replace(1:ind - 1) // new_str // str_replace(ind + len(old_str):len(str_replace))
         enddo
      endif

   return
   endfunction str_replace


   !================================================================================================
   function str_remove_chars(string, chars)
   !! Function that removes the characters of a string from another string
   implicit none
   character(len = :), allocatable :: str_remove_chars  !! returned string
   character(len = *), intent(in)  :: string            !! string to be modified
   character(len = *), intent(in)  :: chars             !! list of characters to remove

      integer(kind = I4) :: i, j

      str_remove_chars = ''

      ! Loop through each character in the input string
      j = 1
      do i = 1, len_trim(string)

         if ( index( chars, string(i:i) ) == 0 ) then

            str_remove_chars = str_remove_chars // String(i:i)

            j = j + 1

         endif

      enddo

   return
   endfunction str_remove_chars


endmodule files
