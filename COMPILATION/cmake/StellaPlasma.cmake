############################################################
#  PLASMA linear algebra library (STELLA_ENABLE_PLASMA)
############################################################
# Provides the variable STELLA_PLASMA_LIBRARIES, which lists the
# targets that need to be linked to use PLASMA.
#
# If PLASMA_DIR (cmake variable or environment variable) points
# to the install prefix of PLASMA, that installation is used.
# Otherwise the git submodule EXTERNALS/plasma is compiled as an
# external project in the build directory. PLASMA
# requires the CBLAS and LAPACKE libraries, extra cmake flags for
# PLASMA can be passed through PLASMA_CMAKE_FLAGS (a cmake list).
# Only the PLASMA libraries are compiled, since the PLASMA tester
# needs LAPACKE routines missing from older LAPACK versions, e.g.
# LAPACKE_zlangb in LAPACK 3.10 (Ubuntu 22.04). The libraries are
# used from the PLASMA build directory, since installing PLASMA
# would also install (and therefore require) the tester.
#
# Note that the plasmaConfig.cmake file installed by PLASMA does
# not define its targets, so we define imported targets ourselves.
############################################################

if (NOT STELLA_ENABLE_DOUBLE)
  message(FATAL_ERROR "STELLA_ENABLE_PLASMA requires STELLA_ENABLE_DOUBLE, since stella is bound to the double precision PLASMA routines")
endif()

if (NOT PLASMA_DIR AND DEFINED ENV{PLASMA_DIR})
  set(PLASMA_DIR "$ENV{PLASMA_DIR}")
endif()

if (PLASMA_DIR)

  find_library(_plasma_library plasma PATHS "${PLASMA_DIR}" PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
  find_library(_plasma_core_blas_library plasma_core_blas PATHS "${PLASMA_DIR}" PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
  if (NOT _plasma_library OR NOT _plasma_core_blas_library)
    message(FATAL_ERROR "Could not find the PLASMA libraries in PLASMA_DIR=${PLASMA_DIR}")
  endif()
  message(STATUS "Using PLASMA from ${PLASMA_DIR}")

else()

  get_filename_component(_plasma_source_dir "${CMAKE_CURRENT_LIST_DIR}/../../EXTERNALS/plasma" ABSOLUTE)
  set(_plasma_build_dir "${CMAKE_CURRENT_BINARY_DIR}/EXTERNALS/plasma/build")
  if (NOT EXISTS "${_plasma_source_dir}/CMakeLists.txt")
    message(FATAL_ERROR "The PLASMA submodule is missing, run: git submodule update --init EXTERNALS/plasma")
  endif()
  message(STATUS "Compiling PLASMA from the submodule ${_plasma_source_dir}")

  set(_plasma_library "${_plasma_build_dir}/${CMAKE_SHARED_LIBRARY_PREFIX}plasma${CMAKE_SHARED_LIBRARY_SUFFIX}")
  set(_plasma_core_blas_library "${_plasma_build_dir}/${CMAKE_SHARED_LIBRARY_PREFIX}plasma_core_blas${CMAKE_SHARED_LIBRARY_SUFFIX}")

  include(ExternalProject)
  ExternalProject_Add(plasma_external
    SOURCE_DIR "${_plasma_source_dir}"
    BINARY_DIR "${_plasma_build_dir}"
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=Release
      ${PLASMA_CMAKE_FLAGS}
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target plasma plasma_core_blas
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS "${_plasma_library}" "${_plasma_core_blas_library}"
    )

endif()

# Imported targets for the PLASMA libraries
add_library(stella_plasma_core_blas SHARED IMPORTED GLOBAL)
set_target_properties(stella_plasma_core_blas PROPERTIES IMPORTED_LOCATION "${_plasma_core_blas_library}")
add_library(stella_plasma SHARED IMPORTED GLOBAL)
set_target_properties(stella_plasma PROPERTIES IMPORTED_LOCATION "${_plasma_library}")
if (TARGET plasma_external)
  add_dependencies(stella_plasma_core_blas plasma_external)
  add_dependencies(stella_plasma plasma_external)
endif()
set(STELLA_PLASMA_LIBRARIES stella_plasma stella_plasma_core_blas)
