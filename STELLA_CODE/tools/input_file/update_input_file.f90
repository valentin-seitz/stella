!###############################################################################
!############################## UPDATE INPUT FILE ##############################
!###############################################################################

program update_input_file

   implicit none

   character(500) :: input_file_name
      
   logical :: debug = .true.
   logical :: list

   !call parse_command_line()

   if (debug) write (*, *) '   '
   if (debug) write (*, *) 'update_input_file::read_input_file_name'
   call get_input_file_name(input_file_name)
   write(*,*) 'input_file_name: ', input_file_name
   
   

contains 

!###############################################################################
!################################# FILE UTILS ##################################
!###############################################################################
 
   !============================================================================
   !===================== Read the name of the input file ======================
   !============================================================================  
   subroutine get_input_file_name(input_file_name)

      use command_line, only: cl_getarg, cl_iargc

      implicit none

      character(500), intent(out) :: input_file_name
      integer :: string_length, ierr

      ! Get the first argument from the command line
      ! This should correspond to the input file name
      if (cl_iargc() /= 0) then
         call cl_getarg(1, input_file_name, string_length, ierr)
         if (ierr /= 0) then
            write (*,*) "Error getting input file name."
            write (*,*) "Please specify the name of the input file as the first argument."
         end if
      end if
      
      ! Remove the extension '.in'
      string_length = len_trim(input_file_name)
      if (string_length > 3 .and. input_file_name(string_length - 2:string_length) == ".in") then
         input_file_name = input_file_name(1:string_length - 3)
      end if

   end subroutine get_input_file_name
   
   !============================================================================
   !=========================== COMMAND LINE FLAGS =============================
   !============================================================================  
   ! Parse some basic command line arguments. 
   ! Currently just 'version' and 'help'. 
   !============================================================================ 
   subroutine parse_command_line()
   
      use git_version, only: get_git_version

      implicit none
      
      integer :: arg_count, arg_n
      integer :: arg_length
      character(len=:), allocatable :: argument
      character(len=*), parameter :: endl = new_line('a')

      arg_count = command_argument_count()

      do arg_n = 0, arg_count
         call get_command_argument(1, length=arg_length)
         if (allocated(argument)) deallocate (argument)
         allocate (character(len=arg_length)::argument)
         call get_command_argument(1, argument)

         if ((argument == "--version") .or. (argument == "-v")) then
            write (*, '("stella version ", a)') get_git_version()
            stop
         else if ((argument == "--help") .or. (argument == "-h")) then
            write (*, '(a)') "stella [--version|-v] [--help|-h] [input file]"//endl//endl// &
               "stella is a flux tube gyrokinetic code for micro-stability and turbulence "// &
               "simulations of strongly magnetised plasma"//endl// &
               "For more help, see the documentation at https://stellagk.github.io/stella/"//endl// &
               "or create an issue https://github.com/stellaGK/stella/issues/new"//endl// &
               endl// &
               "  -h, --help     Print this message"//endl// &
               "  -v, --version  Print the stella version"
            stop
         end if
      end do
      
   end subroutine parse_command_line

end program update_input_file
