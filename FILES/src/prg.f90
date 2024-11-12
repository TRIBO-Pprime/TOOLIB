program test_files
use files

implicit none

logical(kind = 4) :: dir_exists
integer(kind = 4) :: exit_status, i_g

character(len = :), allocatable :: slogan, str

character(len = 512), allocatable, dimension(:) :: vpath

character(len = 512), allocatable, dimension(:) :: list_f90

character(len = 512) :: path, cwd, str_tmp
character(len =   1) :: os_sep

   write(*, *) 'Give a sting from which spaces and hyphen have to be removed'
   read(*, '(a)') str_tmp
   str = str_remove_chars(string = trim(str_tmp), chars = '- ')
   write(*,*) str

   path = "/my_name/is bond/james/bond/  .007"
   call path2vec( file_path = trim(path), vec_path = vpath )

   write(*, *) size(vpath), ' folders in path'
   write(*, '(*(a))') ( trim(vpath(i_g)), ' | ', i_g = 1, ubound(vpath, 1) )

   call vec2path( file_path = str, vec_path = vpath )

   write(*, *) 'complete path: ', trim(str)

   write(*, *) 'file is: "', filename( trim(path) ), '"'

   deallocate( str )

   !===================================
   write(*, *) '----------------------'
   !===================================

   call getcwd( cwd )
   call make_path(wkd = trim(cwd), file_path = "out/a/b/c/d/file.txt", exit_status = exit_status)

   !===================================
   write(*, *) '----------------------'
   !===================================

   call list_files(dir = "div", list = list_f90, ext = "f90")

   do i_g = 1, ubound( list_f90, 1 )
      write(*,*) trim( list_f90(i_g) )
   enddo

   !===================================
   write(*, *) '----------------------'
   !===================================

   slogan = str_replace(  string = 'think different, think Linux',   &  !
                         old_str = 'think',                          &  !
                         new_str = 'be',                             &  !
                           place = 0 )                                  ! first instance only
   write(*, *) slogan

   slogan = str_replace(  string = 'think different, think Linux',   &  !
                         old_str = 'think',                          &  !
                         new_str = 'be',                             &  !
                           place = 1 )                                  ! last instance only
   write(*, *) slogan

   slogan = str_replace(  string = 'think different, think Linux',   &  !
                         old_str = 'think',                          &  !
                         new_str = 'be',                             &  !
                           place = 2 )                                  ! both instances
   write(*, *) slogan

   !===================================
   write(*, *) '----------------------'
   !===================================

   os_sep = dir_separator()
   write(*, *) 'The local directory separator is:', os_sep

   !===================================
   write(*, *) '----------------------'
   !===================================

   call getcwd( path )
   write(*, *) 'The working directory is: ', trim( path )

   !===================================
   write(*, *) '----------------------'
   !===================================

   call mkdir(wkd = trim(path), directory = 'scratch', sep = os_sep, exit_status = exit_status)

   inquire(file = 'scratch', exist = dir_exists)
   if ( dir_exists ) then
      write(*, *) 'Directory ''scratch'' successfuly created'
   else
      write(*, *) 'Directory ''scratch'' not created'
   endif

   !===================================
   write(*, *) '----------------------'
   !===================================

   call chdir( "scratch" )
   call getcwd( path )
   write(*, *) 'The new working directory is: ', trim( path )

   !===================================
   write(*, *) '----------------------'
   !===================================


endprogram test_files
