!###############################################################################
!############################## UPDATE INPUT FILE ##############################
!###############################################################################

program update_input_file

   use file_utils, only: init_file_utils, run_name, get_unused_unit
   use file_utils, only: open_other_input_file, input_unit_exist

   implicit none
   
   ! Fortran uses the unit number to access the file with later read and write statements  
   ! This script will write a new file called <run_name>_with_defaults.in
   integer :: input_with_defaults_unit_number 
   character(500) :: input_file_with_defaults_name
   
   ! Turn debug on for now
   logical :: debug = .true.
   
   ! Mandatory argument to use <init_file_utils> even though we do not use it
   logical :: list 
   
   ! This is used a lot
   logical :: nml_exist
   
   !----------------------------------------------------------------------------
   
   if (debug) write (*, *) 'update_input_file::start'

   ! Allow for --version and --help flags for the executable
   call parse_command_line()

   ! Use stella's routine to read and clean the input file
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

!###############################################################################
!################################# PARAMETERS ##################################
!###############################################################################

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
   
   
!###############################################################################
!############################# COMMAND LINE FLAGS ##############################
!###############################################################################
! Parse some basic command line arguments. 
! Currently just 'version' and 'help'. 
    
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
