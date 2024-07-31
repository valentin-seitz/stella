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
STELLA_PARENT_DIR := $(realpath $(dir $(MAKEFILE_PATH)))
STELLA_CODE_DIR := $(STELLA_PARENT_DIR)/STELLACODE
EXTERNALS_DIR := $(STELLA_PARENT_DIR)/externals
BUILD_DIR := $(STELLA_PARENT_DIR)/build_make

####################################################################
#                        PLATFORM DEPENDENCE                       #
####################################################################
# Each platform will need slightly different compiler flags, thus
# each platform has its owm 'Makefiles/Makefile.STELLA_SYSTEM' file. 
####################################################################

# A user can include a 'Makefile.local' to change some compiler flags
# Here "sinclude" will include the file, without an error if it doesn't exist
sinclude Makefile.local

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
include $(STELLA_PARENT_DIR)/Makefiles/Makefile.$(STELLA_SYSTEM) 

# However, Makefile.local may override some options set in Makefile.$(STELLA_SYSTEM),
# thus it is included before and after Makefile.$(_SYSTEM), if it exists
sinclude Makefile.local

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

.PHONY: create-build-directory
export BUILD_F90_DIR_NAME := f90fromfpp
export BUILD_DIR := $(STELLA_PARENT_DIR)/build_make
export BUILD_MODULES_DIR := $(STELLA_PARENT_DIR)/build_make/modules
export BUILD_OBJECTS_DIR := $(STELLA_PARENT_DIR)/build_make/objects
export BUILD_F90_DIR := $(STELLA_PARENT_DIR)/build_make/$(BUILD_F90_DIR_NAME)
create-build-directory: 
	mkdir -p $(BUILD_DIR) 
	mkdir -p $(BUILD_MODULES_DIR)
	mkdir -p $(BUILD_OBJECTS_DIR)
	mkdir -p $(BUILD_F90_DIR) 
	@echo "Created the build_cmake directies"
	@echo "  "
	
# -o is a compiler option that specifies the name of the output file
# -c is a compiler option that tells the compiler to generate an object file
# $@ is the name of the target being generated -> 'stella'
# $< is the first prerequisite -> the first file in '$(external_modules_dependencies)'
# $^ are all the prerequisites -> all files in '$(external_modules_dependencies)' and '$(stella_mod)'
# -c $^ -o $@

# We process the following files
.SUFFIXES: .fpp .f90 .F90 .c .o .i90

# Build fortran files '*.f90' from '*.fpp' using the C preprocessor (CPP) compiler 
%.f90: %.fpp
	$(CPP) $(CPPFLAGS) -J$(BUILD_MODULES_DIR) $< $(BUILD_F90_DIR)/$@

# Build object files '*.o' from '*.f90' files using the Fortran f90 compiler, e.g., ifort = The Intel Fortran compiler
%.o: %.f90
	$(FC) $(F90FLAGS) $(INC_FLAGS) -J$(BUILD_MODULES_DIR) -o $(BUILD_OBJECTS_DIR)/$@ -c $<

# Build object files '*.o' from '*.c' using the C compiler 
%.o: %.c
	$(CC) $(CFLAGS) -J$(BUILD_MODULES_DIR) -c $(BUILD_F90_DIR)/$<
	
####################################################################
#                         EXTERNAL MODULES                         #
####################################################################
# Stella uses the external modules 'git_version' and 'neasyf'. 
# Note that the automated Fortran tests use the external module 
# 'pFunit' is processed in tests/automated_fortran_tests/Makefile
####################################################################

# --------------------- UTILS AND MINI_LIBSTELL --------------------

# We will compile the static library <mini_libstell> and we link LIBSTELL_LIB to the static library
# In a static library, the modules are bound into the executable file before execution. 
# Static libraries are commonly named libname.a. The .a suffix refers to archive. 
export MINILIBSTELL_DIR=$(EXTERNALS_DIR)/mini_libstell
export UTILS_DIR=$(EXTERNALS_DIR)/utils
LIBSTELL_LIB=$(MINILIBSTELL_DIR)/mini_libstell.a

# stella needs two modules 'utils.a' and 'mini_libstell.a'
# we will build all modules using the 'modules' target 
.PHONY: external_modules 
external_modules: utils.a mini_libstell.a 

# The 'utils.a' static library contains all external python scripts
# Note that in UTIL_OBJ we do not need to include all files in utils,
# because 'stella' will compile all other files which aren't included in 'utils' 
# $@ is the name of the target being generated -> 'stella'
# $< is the first prerequisite -> the first file in '$(external_modules_dependencies)'
# $^ are all the prerequisites -> all files in '$(external_modules_dependencies)' and '$(stella_mod)'
UTIL_OBJ = command_line.o constants.o convert.o fft_work.o file_utils.o gauss_quad.o job_manage.o linear_solve.o mp.o mp_lu_decomposition.o mt19937.o netcdf_utils.o ran.o redistribute.o smooth_step.o sort.o spfunc.o spl.o system_fortran.o text_options.o
#UTIL_OBJ = $(addprefix $(BUILD_OBJECTS_DIR)/, command_line.o constants.o convert.o fft_work.o file_utils.o gauss_quad.o job_manage.o linear_solve.o mp.o mp_lu_decomposition.o mt19937.o netcdf_utils.o ran.o redistribute.o smooth_step.o sort.o spfunc.o spl.o system_fortran.o text_options.o )
utils.a: $(UTIL_OBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(addprefix $(BUILD_OBJECTS_DIR)/, $^)
	$(RANLIB) $@

# The mini_libstell.a is an external library to read VMEC files
# It has no dependencies, so only write the build commands
mini_libstell.a:
	$(MAKE) -C $(MINILIBSTELL_DIR) 
	
# --------------------- GIT_VERSION AND NEASYF ---------------------

# External libraries
GIT_VERSION_PARENT_DIR := $(EXTERNALS_DIR)/git_version
GIT_VERSION_DIR := $(EXTERNALS_DIR)/git_version/src
NEASYF_DIR := $(EXTERNALS_DIR)/neasyf/src

# Files that we know are going to be in the submodule directories.
# This is in case the directories exist but the submodules haven't actually been initialised
GIT_VERSION_SENTINEL := $(GIT_VERSION_PARENT_DIR)/Makefile
NEASYF_SENTINEL := $(NEASYF_DIR)/neasyf.f90

# If make is killed or interrupted during the execution of their recipes, 
# the targets in ".PRECIOUS" are not deleted.
.PRECIOUS: $(GIT_VERSION_SENTINEL) $(NEASYF_SENTINEL)

# Make sure the submodules exist; they all currently share the same recipe
$(GIT_VERSION_SENTINEL) $(NEASYF_SENTINEL) submodules:
	@echo "Downloading submodules"
	git submodule update --init --recursive

# Dump the compilation flags to a file, so we can check if they change between
# invocations of `make`. The `cmp` bit checks if the file contents
# change. Adding a dependency of a file on `.compiler_flags` causes it to be
# rebuilt when the flags change. Taken from
# https://stackoverflow.com/a/3237349/2043465
COMPILER_FLAGS_CONTENTS = "FC = $(FC)\n CPPFLAGS = $(CPPFLAGS)\n F90FLAGS = $(F90FLAGS)\n INC_FLAGS = $(INC_FLAGS)\n CFLAGS = $(CFLAGS)"
COMPILER_FLAGS_CONTENTS += "\n FORTRAN_GIT_DEFS = $(FORTRAN_GIT_DEFS)"
.PHONY: force
.compiler_flags: force
	@echo -e $(COMPILER_FLAGS_CONTENTS) | cmp -s - $@ || echo -e $(COMPILER_FLAGS_CONTENTS) > $@

# Use the make rules from fortran-git-version, which requires ensuring the submodules exist
include $(GIT_VERSION_PARENT_DIR)/Makefile

# We have some extra macros to add when preprocessing this file. We'll
# need the submodule to exist first
git_version_impl.f90: $(GIT_VERSION_DIR)/git_version_impl.F90 | $(GIT_VERSION_SENTINEL)
	$(CPP) $(CPPFLAGS) $(FORTRAN_GIT_DEFS) $< $@

# The fortdep script doesn't know about Fortran submodules, so we need
# to write the dependency ourselves
git_version_impl.o: git_version_impl.f90 git_version.o .compiler_flags
$(GIT_VERSION_DIR)/git_version_impl.F90: .compiler_flags

# There is an undefined reference to `get_git_state` or `get_git_hash` error. 
# We need to manually add the dependency of stella on `$(OBJDIR)/git_version_impl.o`
# by writing 'stella: $(external_modules) $(stella_mod)' when defining the target.
# This file needs to be first in the list of dependencies so that the file that 
# needs the functions can find them. Yes this is unpleasant, no we (probably?) 
# can`t automate this like for the other dependencies. `fortdep` doesn't 
# understand Fortran submodules, so it can't pick this up.
external_modules_dependencies := git_version_impl.o








############################################################### SPECIAL RULES

# comment this out to keep intermediate .f90 files
#.PRECIOUS: $(F90FROMFPP)

.INTERMEDIATE: stella_transforms.f90 stella_io.f90 stella_save.f90 \
		mp.f90 fft_work.f90 response_matrix.f90 sources.f90 \
		fields.f90 mp_lu_decomposition.f90 git_version_impl.f90

############################################################# MORE DIRECTIVES


unlink:
	-rm -f $(F90FROMFPP)

revision:
	@LANG=C svn info | awk '{if($$1=="Revision:") printf("%20d",$$2) }' > Revision

TAGS:	*.f90 *.fpp */*.f90 */*.fpp
	etags $^

############################################################# Documentation

create_namelist_markdown:
	docs/pages/user_manual/namelist_files/combine_namelists.sh

ifneq ("$(wildcard $(shell which $(FORD) 2>/dev/null))","")
check_ford_install:
	@echo "Using ford at $(shell which $(FORD))"
else
check_ford_install:
	@echo "Ford command $(FORD) not in path -- is it installed?\\n\\tConsider installing with 'pip install --user ford' and add ${HOME}/.local/bin to PATH" ; which $(FORD)
endif

doc: docs/stella_docs.md create_namelist_markdown check_ford_install
	$(FORD) $(INC_FLAGS) -r $(GIT_VERSION) $<

cleandoc:
	@echo "FORD docs"
	-rm -rf docs/html


####################################################################
#                              STELLA                              #
####################################################################

include $(STELLA_PARENT_DIR)/Makefile.stella
