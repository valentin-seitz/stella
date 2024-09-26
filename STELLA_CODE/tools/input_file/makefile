
# We set the .DEFAULT_GOAL and for GNU make 3.80 or less, we use the 'all' target 
.DEFAULT_GOAL := update-input-file
all: update-input-file

# Targets which do not create files
.PHONY: all build-stuffs

# Set MAKEFILE_PATH to the location of this Makefile.
# Using 'abspath' will not resolve symbolic links
MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST))) 
PROJECT_DIR := $(realpath $(dir $(MAKEFILE_PATH))/../../../)
COMPILATION_DIR := $(PROJECT_DIR)/COMPILATION
BUILD_DIR := $(COMPILATION_DIR)/build_make
FPP := $(BUILD_DIR)/f90fromfpp
OBJ := $(BUILD_DIR)/objects
MOD := $(BUILD_DIR)/modules
SRC := $(PROJECT_DIR)/STELLA_CODE/tools/input_file
UTILS := $(PROJECT_DIR)/EXTERNALS/utils
GIT_VERSION_DIR := $(PROJECT_DIR)/EXTERNALS/git_version/src

export VPATH := $(VPATH):$(UTILS):$(GIT_VERSION_DIR):$(MOD):$(OBJ):$(FPP):

# Create the directories which will contain the compiled files
# These commands will be triggered if $(BUILD_DIR) is a prerequisite, but it does not exist yet.
$(BUILD_DIR) $(OBJ) $(MOD) $(FPP): 
	@echo "                     " 
	@echo "   Create the COMPILATION/build_make/* directories"
	@echo "      - $(BUILD_DIR) "
	@echo "      - $(MOD) "
	@echo "      - $(OBJ) "
	@echo "      - $(FPP) "
	@mkdir -p $(BUILD_DIR) 
	@mkdir -p $(MOD)
	@mkdir -p $(OBJ)
	@mkdir -p $(FPP)
	@echo "                     "

# Get compilation flags 
include $(COMPILATION_DIR)/Makefile.compilerflags

# Objects needed to compile the program 
objects = $(OBJ)/update_input_file.o
externals = $(OBJ)/file_utils.o $(OBJ)/command_line.o $(OBJ)/git_version.o $(OBJ)/git_version_impl.o

# Dependencies
$(OBJ)/update_input_file.o: $(OBJ)/file_utils.o $(OBJ)/git_version.o $(OBJ)/git_version_impl.o
$(OBJ)/file_utils.o: $(OBJ)/command_line.o $(FPP)/file_utils.f90
$(OBJ)/command_line.o: 
$(OBJ)/git_version.o:

# Libraries
# Here $LIBSTELL_LIB points to utils/mini_libstell which contains a stellarator resource library
LIBS += $(DEFAULT_LIB) $(MPI_LIB) $(FFT_LIB) $(NETCDF_LIB) $(HDF5_LIB) \
		$(NAG_LIB) $(SFINCS_LIB) $(PETSC_LIB) $(LAPACK_LIB) $(LIBSTELL_LIB) 

# Build update-input-file executable
update-input-file: $(objects) $(externals) | $(OBJ)
	@$(LD) $(LDFLAGS) -o $@ $(filter %.o, $^) $(LIBS)
	@echo "       ==> Created update-input-file executable"; echo "  "

#===================================================================
#                           COMPILE F90                            #
#===================================================================
	
# The -J $(MOD) flag makes sure that the compiler writes the *.mod files to the $(MOD) directory
# The -I $(MOD) flag makes sure that the compiler looks for the *.mod files in the $(MOD) directory
MOD_FLAGS_INPUTFILE := -I$(MOD) -J$(MOD)

# Pattern rules to compile the scripts
# Create *.o files from *.f90 files using the $(FC) compiler 
# Note that the *.f90 files can be located in $(SRC) or $(FPP)
# Where the $(FPP) directory contains the prepocessed *.f90 from *.fpp files 
$(OBJ)/%.o: $(SRC)/%.f90 | $(OBJ)
	@echo "   (update-input-file) *.f90 to *.o: $(notdir $<)" 
	@$(FC) $(F90FLAGS) $(INC_FLAGS) $(MOD_FLAGS_INPUTFILE) -o $@ -c $< 
$(OBJ)/%.o: $(UTILS)/%.f90 | $(OBJ)
	@echo "   (update-input-file) *.f90 to *.o: $(notdir $<)" 
	@$(FC) $(F90FLAGS) $(INC_FLAGS) $(MOD_FLAGS_INPUTFILE) -o $@ -c $< 
$(OBJ)/%.o: $(UTILS)/%.f90 | $(OBJ)
	@echo "   (update-input-file) *.f90 to *.o: $(notdir $<)" 
	@$(FC) $(F90FLAGS) $(INC_FLAGS) $(MOD_FLAGS_INPUTFILE) -o $@ -c $< 
$(OBJ)/%.o: $(FPP)/%.f90 | $(OBJ)
	@echo "   (update-input-file) *.f90 to *.o: $(notdir $<)" 
	@$(FC) $(F90FLAGS) $(INC_FLAGS) $(MOD_FLAGS_INPUTFILE) -o $@ -c $< 

	
# FPP compiler
$(FPP)/%.f90: $(FPP)/%.fpp $(FPP)/define.inc | $(FPP)
	@echo "   (update-input-file) *.fpp to *.f90: $(notdir $<)"
	@$(CPP) $(CPPFLAGS) $< $@
$(FPP)/%.fpp: $(UTILS)/%.fpp | $(FPP)
	@echo "   (update-input-file) ln -s STELLA_CODE/$(notdir $<) COMPILATION/$(FPP_NAME)/$(notdir $<)"	
	@ln -sf $< $(FPP)
$(FPP)/define.inc: | $(FPP)
	@echo "   (update-input-file) ln -s EXTERNALS/utils/define.inc COMPILATION/$(FPP_NAME)/define.inc"	
	@ln -sf $(UTILS)/define.inc $(FPP)/define.inc
	
.PRECIOUS: $(OBJ)/%.o $(OBJ_GENERAL)/%.o $(FPP)/%.f90 $(FPP)/%.fpp
	
#===================================================================
#                              CLEAN                               #
#===================================================================

# Remove the compiled files 
clean:
	@$(MAKE) --no-print-directory -C $(COMPILATION_DIR) clean

#===================================================================
#                            GIT VERSION                           #
#===================================================================

include $(COMPILATION_DIR)/Makefile.externals.git_version
