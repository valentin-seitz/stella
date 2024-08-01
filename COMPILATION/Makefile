####################################################################
#                                                                  #
#       Makefile for the stella gyrokinetic turbulence code        # 
#                                                                  #
####################################################################
#
# Makefile written by Bill Dorland and Ryusuke Numata
# Cleaned-up and modified by Hanne Thienpondt on 07/2024
#
# In this Makefile the flags for the Fortran and C compilers are set.
#
# In the ~/.source.sh or ~/.bashrc file of your computer define, e.g., 
#     export STELLA_SYSTEM = 'marconi'
# which matches stella/Makefiles/Makefile.$STELLA_SYSTEM
#
# * Available Compilers (tested on limited hosts)
#   (must be Fortran 95 Standard compliant)
#
# Intel ifort
# GNU's gfortran and g95
# IBM XL Fortran xlf90
# PathScale Compiler Suite pathf90
# The Portland Group pgf90
# NAGWare f95 (v5.1)
# Lahey/Fujitsu Fortran lf95
# 
# * Frequently Tested Hosts, Systems
#
# Standard Linux
# Standard Mac OS X with MacPorts
# Franklin at NERSC and Jaguar at NCCS (Cray XT4 with PGI)
# Bassi at NERSC (IBM Power 5 with IBM XL Fortran)
# Ranger (Sun Constellation Linux Cluster with Intel)
####################################################################


####################################################################
#                             Switches                             #
####################################################################
# In the brackets, the accepted values are shown, where "undefined" means blank.
# Note that switches with (bin) only check whether they are defined or not,
# hence the value you give it does not matter, e.g., DEBUG=off means DEBUG=on.
# The "?=" sets a variable if it does not already have a value. 
####################################################################

DEBUG ?=                   # Turns on debug mode (bin) 
TEST ?=                    # Turns on test mode (bin) 
OPT ?= on                  # Optimization (on, aggressive, undefined)
STATIC ?=                  # Prevents linking with shared libraries (bin)
DBLE ?= on                 # Promotes precisions of real and complex (bin) 
USE_NETCDF ?= on           # Uses netcdf library (bin)
USE_PARALLEL_NETCDF ?=     # Uses parallel netcdf library
USE_HDF5 ?=                # Uses hdf5 library (bin)
USE_LOCAL_RAN ?=           # Use local random number generator (mt,undefined) (see also README)
USE_LOCAL_SPFUNC ?=        # Use local special functions (bin)
USE_NAGLIB ?=              # Use nag libraray (spfunc,undefined)
USE_SFINCS ?=              # Link to sfincs library at compilation
USE_LAPACK ?= on           # Use LAPACK, needed for test particle collisions
HAS_ISO_C_BINDING ?= on    # Does the compiler support the iso_c_binding features of Fortran? (needed for local parallel LU decomposition) 

# Which FFT library to use (fftw,fftw3,mkl_fftw,undefined)
# We can't have any spaces behind 'fftw3' so we can't put the comment behind it
USE_FFT ?= fftw3

# The following variables can be set in platform-dependent Makefile's
# For example, inside stella/Makefiles/Makefile.marconi

MAKE		 = make
CPP		     = cpp
CPPFLAGS	 = -P -traditional
export FC	 = f90
export MPIFC?= mpif90
H5FC		?= h5fc
H5FC_par	?= h5pfc
F90FLAGS	 =
F90OPTFLAGS	 =
CC		     = cc
MPICC		?= mpicc
H5CC		?= h5cc
H5CC_par	?= h5pcc
CFLAGS		 =
COPTFLAGS 	 =
LD 		     = $(FC)
LDFLAGS 	 = $(F90FLAGS)
ARCH 		 = ar
ARCHFLAGS 	 = cr
RANLIB		 = ranlib
AWK 		 = awk
PERL		 = perl
FORD       	?= ford

MPI_INC ?=
MPI_LIB ?=
FFT_INC ?=
FFT_LIB ?=
NETCDF_INC ?=
NETCDF_LIB ?=
HDF5_INC ?=
HDF5_LIB ?=
NAG_LIB ?=
NAG_PREC ?= dble
SFINCS_LIB ?=
SFINCS_INC ?=
PETSC_LIB ?=
PETSC_INC ?=
LIBSTELL_LIB ?=

####################################################################
#                            DIRECTORIES                           #
####################################################################

# In the make file we can use CURDIR, which is set to the absolute 
# pathname of the current working directory. When GNU make starts 
# (after it has processed any -C options) it sets the variable
# CURDIR to the pathname of the current working directory.
CURRENT_DIRECTORY = $(CURDIR)

# As make reads various makefiles, including any obtained from the 
# MAKEFILES variable, the command line, the default files, or from 
# include directives, their names will be automatically appended to 
# the MAKEFILE_LIST variable. They are added right before make begins
# to parse them. This means that if the first thing a makefile does 
# is examine the last word in this variable, it will be the name of 
# the current makefile. Once the current makefile has used include, 
# however, the last word will be the just-included makefile. 

# Set MAKEFILE_PATH to the location of this Makefile.
# Using 'abspath' will not resolve symbolic links
MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST))) 

# Get the directory in which the Makefile is located, 
# and the directory of the stella code 
PROJECT_DIR := $(realpath $(dir $(MAKEFILE_PATH))/..)
STELLA_DIR := $(PROJECT_DIR)/STELLACODE
EXTERNALS_DIR := $(PROJECT_DIR)/externals
BUILD_DIR := $(PROJECT_DIR)/build_make 
COMPILATION_DIR := $(PROJECT_DIR)/COMPILATION

# Check directories
.PHONY: print-directories
print-directories: 
	@echo "                     " 
	@echo "Directories for stella:" 
	@echo "   - CURRENT_DIRECTORY is $(CURRENT_DIRECTORY)"  
	@echo "   - PROJECT_DIR is $(PROJECT_DIR)" 
	@echo "   - STELLA_DIR is $(STELLA_DIR)" 
	@echo "   - COMPILATION_DIR is $(COMPILATION_DIR)" 
	@echo "   - EXTERNALS_DIR is $(EXTERNALS_DIR)"  
	@echo "                     " 

# Keep the compilated files in seperate folders
export BUILD_F90_DIR_NAME := f90fromfpp
export BUILD_DIR := $(COMPILATION_DIR)/build_make
export BUILD_MODULES_DIR := $(BUILD_DIR)/modules
export BUILD_OBJECTS_DIR := $(BUILD_DIR)/objects
export BUILD_F90_DIR := $(BUILD_DIR)/$(BUILD_F90_DIR_NAME)

# Keep the compilated files in seperate folders
# These commands can be triggered by the target 'create-build-directories'
# or they will be triggered if $(BUILD_DIR) is a prerequisite, but it does not exist yet.
.PHONY: create-build-directories
$(BUILD_DIR) create-build-directories: 
	@echo "                     " 
	@echo "Create the COMPILATION/build_make/* directories"
	@echo "   - $(BUILD_DIR) "
	@echo "   - $(BUILD_MODULES_DIR) "
	@echo "   - $(BUILD_OBJECTS_DIR) "
	@echo "   - $(BUILD_F90_DIR) "
	@mkdir -p $(BUILD_DIR) 
	@mkdir -p $(BUILD_MODULES_DIR)
	@mkdir -p $(BUILD_OBJECTS_DIR)
	@mkdir -p $(BUILD_F90_DIR)   
	@echo "                     "
	 

####################################################################
#                        PLATFORM DEPENDENCE                       #
####################################################################
# Each platform will need slightly different compiler flags, thus
# each platform has its owm 'Makefiles/Makefile.STELLA_SYSTEM' file. 
####################################################################

# A user can include a 'Makefile.local' to change some compiler flags
# Here "sinclude" will include the file, without an error if it doesn't exist
sinclude $(COMPILATION_DIR)/Makefile.local

# Include system-dependent make variables, defined in stella/Makefiles/Makefile.$STELLA_SYSTEM 
# where e.g., export STELLA_SYSTEM='marconi' inside e.g. ~/.bashrc defined this system variable on your computer 
ifndef STELLA_SYSTEM
	ifdef SYSTEM
$(warning WARNING: SYSTEM environment variable is obsolete, use GK_SYSTEM instead)
$(info    )
	GK_SYSTEM = $(SYSTEM)
	else ifdef GK_SYSTEM
$(warning WARNING: GK_SYSTEM environment variable is obsolete, use STELLA_SYSTEM instead)
$(info    )
	STELLA_SYSTEM = $(GK_SYSTEM)
	else
$(info    )
$(info ERROR: STELLA_SYSTEM is not set.)
$(info Please set the following line in your .bashrc file, corresponding to your operating system:)
$(info     >> export STELLA_SYSTEM='gnu_ubuntu')
$(info    )
$(error Aborting stella compilation because STELLA_SYSTEM is not set)
	endif
endif

# Read the make file for this specific operating system, set in STELLA_SYSTEM
# The compiler mode switches (DEBUG, TEST, OPT, STATIC, DBLE) must be set
# before loading Makefile.$(STELLA_SYSTEM) because they may affect compiler options.
include $(COMPILATION_DIR)/Makefiles/Makefile.$(STELLA_SYSTEM) 

# However, Makefile.local may override some options set in Makefile.$(STELLA_SYSTEM),
# thus it is included before and after Makefile.$(_SYSTEM), if it exists
sinclude $(COMPILATION_DIR)/Makefile.local

# Export some variables, so they become available to all child processes
export F90FLAGS
export NETCDF_INC
export NETCDF_LIB

####################################################################
#                        PROCESS SOME FLAGS                        #
####################################################################

# Must invoke full functionality when we do >> make depend
ifeq ($(MAKECMDGOALS),depend)
	MAKE += USE_HDF5=on USE_FFT=fftw3 USE_NETCDF=on \
		USE_LOCAL_BESSEL=on USE_LOCAL_RAN=mt
endif

# Warning that we can't do nonlinear sims without fourier transformations
ifndef USE_FFT
$(warning      )
$(warning USE_FFT is off)
$(warning Be sure that nonlinear run makes no sense)
$(warning      )
endif

ifdef HAS_ISO_C_BINDING
	CPPFLAGS += -DISO_C_BINDING
endif

FC = $(MPIFC)
CC = $(MPICC)
CPPFLAGS += -DMPI

# Flags for <USE_FFT> are (fftw,fftw3,mkl_fftw,undefined)
ifeq ($(USE_FFT),fftw)
	CPPFLAGS += -DFFT=_FFTW_
	ifeq ($(FFT_LIB),)
		FFT_LIB = -lfftw -lrfftw
	endif
endif
ifeq ($(USE_FFT),fftw3)
	CPPFLAGS += -DFFT=_FFTW3_
	ifeq ($(FFT_LIB),)
		FFT_LIB = -lfftw3
	endif
endif
ifeq ($(USE_FFT),mkl_fftw)
	CPPFLAGS += -DFFT=_FFTW_
endif
ifdef USE_NETCDF
	ifeq ($(NETCDF_LIB),)
		NETCDF_LIB = -lnetcdf
    endif
	CPPFLAGS += -DNETCDF
endif
ifdef USE_LAPACK
    LAPACK_LIB ?= -llapack
	CPPFLAGS += -DLAPACK
endif
ifdef USE_HDF5
	FC = $(H5FC_par)
	CC = $(H5CC_par)
	ifdef USE_PARALLEL_NETCDF
		CPPFLAGS += -DNETCDF_PARALLEL
	endif
	CPPFLAGS += -DHDF
endif
ifeq ($(USE_LOCAL_RAN),mt)
	CPPFLAGS += -DRANDOM=_RANMT_
endif
ifdef USE_LOCAL_SPFUNC
	CPPFLAGS += -DSPFUNC=_SPLOCAL_
else
	ifeq ($(findstring spfunc,$(USE_NAGLIB)),spfunc)
		CPPFLAGS += -DSPFUNC=_SPNAG_
	endif
endif
ifdef USE_NAGLIB
	ifeq ($(NAG_PREC),dble)
		ifndef DBLE
$(warning Precision mismatch with NAG libarray)	
		endif
		CPPFLAGS += -DNAG_PREC=_NAGDBLE_
	endif
	ifeq ($(NAG_PREC),sngl)
		ifdef DBLE
$(warning Precision mismatch with NAG libarray)	
		endif
		CPPFLAGS += -DNAG_PREC=_NAGSNGL_
	endif
endif
ifdef USE_SFINCS
	CPPFLAGS += -DUSE_SFINCS
endif

# List of all the used libraries
# Here $LIBSTELL_LIB points to utils/mini_libstell which contains a stellarator resource library
LIBS += $(DEFAULT_LIB) $(MPI_LIB) $(FFT_LIB) $(NETCDF_LIB) $(HDF5_LIB) \
		$(NAG_LIB) $(SFINCS_LIB) $(PETSC_LIB) $(LIBSTELL_LIB) $(LAPACK_LIB)

# Include flags
INC_FLAGS= $(DEFAULT_INC) $(MPI_INC) $(FFT_INC) $(NETCDF_INC) $(HDF5_INC) \
	   $(SFINCS_INC) $(PETSC_INC) $(LAPACK_INC)

# Flags for the fortran compiler and the C compiler
F90FLAGS += $(F90OPTFLAGS)
CFLAGS += $(COPTFLAGS)

####################################################################
#                              RULES                               #
#################################################################### 

# @ supresses the output in the command prompt, remove this to debug	
# -o is a compiler option that specifies the name of the output file
# -c is a compiler option that tells the compiler to generate an object file
# $@ is the name of the target being generated -> 'stella'
# $< is the first prerequisite -> the first file in '$(external_modules_dependencies)'
# $^ are all the prerequisites -> all files in '$(external_modules_dependencies)' and '$(stella_mod)'
# -c $^ -o $@

# We process the following files
.SUFFIXES: .f90 .o 

# For the external libraries, and for stella, we need to build *.o files from *.f90 files
# Build object files '*.o' from '*.f90' files using the Fortran f90 compiler, e.g., ifort = The Intel Fortran compiler
%.o: %.f90
	@$(FC) $(F90FLAGS) $(INC_FLAGS) -J$(BUILD_MODULES_DIR) -o $(BUILD_OBJECTS_DIR)/$@ -c $<
	@echo "    *.f90 to *.o: $(notdir $<)" 
	
####################################################################
#                         EXTERNAL MODULES                         #
####################################################################
# Stella uses the external modules 'git_version' and 'neasyf'. 
# Note that the automated Fortran tests use the external module 
# 'pFunit' is processed in tests/automated_fortran_tests/Makefile
####################################################################

include $(COMPILATION_DIR)/Makefile.externals

############################################################# MORE DIRECTIVES


revision:
	@LANG=C svn info | awk '{if($$1=="Revision:") printf("%20d",$$2) }' > Revision

TAGS:	*.f90 *.fpp */*.f90 */*.fpp
	etags $^


####################################################################
#                              STELLA                              #
####################################################################

include $(COMPILATION_DIR)/Makefile.stella

