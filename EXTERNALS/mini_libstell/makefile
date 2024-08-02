# The variable LIBSTELL_DIR should either be "mini_libstell", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a directory containing libstell .mod files
# elsewhere on your system.
LIBSTELL_DIR = mini_libstell

# The variable LIBSTELL_FOR_STELLA should either be "mini_libstell/mini_libstell.a", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a libstell.a library elsewhere on your system.
LIBSTELL_FOR_STELLA = mini_libstell/mini_libstell.a

# If F90FLAGS has not been set (e.g. running makefile from geo/vmec_interface directory), then set it to
# 'empty' so we can check for it later
F90FLAGS ?= empty

# Load up the compiler flags from the main Stella makefile. We need to include the NETCDF paths for mini_libstell 
EXTRA_COMPILE_FLAGS=$(F90FLAGS) $(NETCDF_INC)
EXTRA_LINK_FLAGS=$(NETCDF_LIB)

export

# Gather the scripts
LIBSTELL_OBJ_FILES = ezcdf.o ezcdf_inqvar.o handle_err.o stel_constants.o \
	ezcdf_GenGet.o ezcdf_opncls.o read_wout_mod.o safe_open_mod.o stel_kinds.o vsvd0.o \
	parse_extension.o

.PHONY: all clean

# The first target in the file, is the one which will be build
all: mini_libstell.a

%.o: %.f90
	@$(FC) $(EXTRA_COMPILE_FLAGS) -w -DNETCDF -c $<
	@echo "    mini_libstell: *.f90 to *.o:  $(notdir $<)"

%.o: %.f
	@$(FC) $(EXTRA_COMPILE_FLAGS) -DNETCDF -c $<
	@echo "    mini_libstell: *.f to *.o:    $(notdir $<)"

%.o: %.F
	@$(FC) $(EXTRA_COMPILE_FLAGS) -DNETCDF -c $<
	@echo "    mini_libstell: *.F to *.o:    $(notdir $<)"

mini_libstell.a: $(LIBSTELL_OBJ_FILES)
	@ar rcs $@ $(LIBSTELL_OBJ_FILES) 

clean:
	@echo "   - remove compiled files (mini_libstell)"
	@rm -f *.o *.mod *.MOD *~ mini_libstell.a

# Dependencies:
ezcdf.o: ezcdf_GenGet.o ezcdf_opncls.o
ezcdf_inqvar.o:
handle_err.o:
stel_constants.o: stel_kinds.o
ezcdf_GenGet.o: ezcdf_opncls.o ezcdf_inqvar.o
ezcdf_opncls.o:
read_wout_mod.o: stel_kinds.o safe_open_mod.o stel_constants.o vsvd0.o ezcdf.o parse_extension.o
safe_open_mod.o:
stel_kinds.o:
vsvd0.o: stel_kinds.o
parse_extension.o:
