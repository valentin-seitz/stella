!###############################################################################
!############################## UPDATE INPUT FILE ##############################
!###############################################################################

program update_input_file

   use file_utils, only: init_file_utils, run_name
   use file_utils, only: open_other_input_file, input_unit_exist

   implicit none
   
   ! Fortran uses the unit number to access the file with later read and write statements
   integer :: user_input_unit_number      ! "<input_file_name>.in" (in_unit)
   integer :: input_unit_number           ! ".<input_file_name>.in" (input_unit_no = out_unit)
   integer :: error_unit_number           ! "<input_file_name>.error" (error_unit_no)
   integer :: input_with_defaults_unit_number

   character(500) :: input_file_name
   character(500) :: input_file_with_defaults_name
   
   integer :: num_input_lines
   integer :: unit_number
      
   logical :: debug = .true.
   logical :: list
   logical :: inp, err
   
   ! This is used a lot
   logical :: nml_exist
   
   
   if (debug) write (*, *) 'update_input_file::start'

   !call parse_command_line()

   call init_file_utils(list)

   ! Get a unit number for the <run_name>_with_defaults.in file we will write
   input_file_with_defaults_name = trim(run_name)//"_with_defaults.in"
   call get_unused_unit(input_with_defaults_unit_number)
   open(unit=input_with_defaults_unit_number, file=input_file_with_defaults_name)
      
   ! Write default values to input file
   call write_parameters_physics()

   ! End code 
   close(unit=input_with_defaults_unit_number)
   
   ! For stella, make stella read inputs from the <run_name>_with_defaults.in file 
   call open_other_input_file(input_file_with_defaults_name) 
   
   ! Perform some checks 
   
   
contains 

   subroutine write_parameters_physics()

      use parameters_physics, only: set_default_parameters_physics => set_default_parameters
      use parameters_physics, only: include_parallel_streaming, include_mirror, nonlinear
      use parameters_physics, only: xdriftknob, ydriftknob, wstarknob, adiabatic_option, prp_shear_enabled
      use parameters_physics, only: hammett_flow_shear, include_pressure_variation, include_geometric_variation
      use parameters_physics, only: include_parallel_nonlinearity, suppress_zonal_interaction, full_flux_surface
      use parameters_physics, only: include_apar, include_bpar, radial_variation
      use parameters_physics, only: beta, zeff, tite, nine, rhostar, vnew_ref
      use parameters_physics, only: g_exb, g_exbfac, omprimfac, irhostar
   
      implicit none
      
      integer :: unit_number
      
      namelist /parameters_physics/ include_parallel_streaming, include_mirror, nonlinear, &
        adiabatic_option, prp_shear_enabled, &
        hammett_flow_shear, include_pressure_variation, include_geometric_variation, &
        include_parallel_nonlinearity, suppress_zonal_interaction, full_flux_surface, &
        include_apar, include_bpar, radial_variation, &
        xdriftknob, ydriftknob, wstarknob,  &
        beta, zeff, tite, nine, rhostar, vnew_ref, &
        g_exb, g_exbfac, omprimfac, irhostar
      
      namelist /parameters/ beta, zeff, tite, nine, rhostar, vnew_ref 
      
      
      if (debug) write (*, *) 'update_input_file::write_parameters_physics'
   
      ! Read in the default parameters set in stella
      call set_default_parameters_physics() 
      if (debug) write (*, *) 'default beta = ', beta
      
      ! Now read in the user-specified input parameters in <parameters_physics>
      unit_number = input_unit_exist("parameters_physics", nml_exist) 
      if (nml_exist) read (unit=unit_number, nml=parameters_physics)
      if (debug) write (*, *) 'read beta = ', beta
      
      ! ADD BACKWARDS COMPATIBILITY
      
      ! Write the namelist to <run_name>_with_defaults.in
      write(unit=input_with_defaults_unit_number, nml=parameters_physics)

   end subroutine write_parameters_physics

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
   !============================= Open output file =============================
   !============================================================================
   ! Open an output file named <input_file_name>.<extension>.
   !============================================================================
   subroutine open_output_file(unit, extension, overwrite_in)

      implicit none

      integer, intent(out) :: unit
      logical, intent(in), optional :: overwrite_in
      character(*), intent(in) :: extension

      logical :: overwrite
      character(500) :: temp_string

      ! Initiate the optional argument
      if (present(overwrite_in)) then
         overwrite = overwrite_in
      else
         overwrite = .true.
      end if

      ! Get a unit for the output file that is not currently in use
      call get_unused_unit(unit)

      ! Create the name of the output file
      temp_string = trim(input_file_name)//extension

      ! If overwrite==True: Create a new output file or replace the existing file
      ! If overwrite==False: Append data to the already existing output file
      if (overwrite) then
         open (unit=unit, file=trim(temp_string), status="replace", action="write")
      else
         open (unit=unit, file=trim(temp_string), status="unknown", action="write", position="append")
      end if

   end subroutine open_output_file
   
   !============================================================================
   !===================== Get an unused unit number for I/O ====================
   !============================================================================
   subroutine get_unused_unit(unit)

      implicit none
      
      integer, intent(out) :: unit
      logical :: od
      
      unit = 50
      do
         inquire (unit=unit, opened=od)
         if (.not. od) return
         unit = unit + 1
      end do
      
   end subroutine get_unused_unit
   
   subroutine strip_comments(line)
      implicit none
      character(*), intent(in out) :: line
      logical :: in_single_quotes, in_double_quotes
      integer :: i, length

      length = len_trim(line)
      i = 1
      in_single_quotes = .false.
      in_double_quotes = .false.
      loop: do
         if (in_single_quotes) then
            if (line(i:i) == "'") in_single_quotes = .false.
         else if (in_double_quotes) then
            if (line(i:i) == '"') in_double_quotes = .false.
         else
            select case (line(i:i))
            case ("'")
               in_single_quotes = .true.
            case ('"')
               in_double_quotes = .true.
            case ("!")
               i = i - 1
               exit loop
            end select
         end if
         if (i >= length) exit loop
         i = i + 1
      end do loop
      line = line(1:i)
   end subroutine strip_comments
   
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
