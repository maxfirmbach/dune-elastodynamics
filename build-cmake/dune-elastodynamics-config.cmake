if(NOT dune-elastodynamics_FOUND)
# Whether this module is installed or not
set(dune-elastodynamics_INSTALLED OFF)

# Settings specific to the module

# Package initialization
# Set prefix to source dir
set(PACKAGE_PREFIX_DIR /home/a11bmafi/Library/dune/dune-elastodynamics)
macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

#report other information
set_and_check(dune-elastodynamics_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(dune-elastodynamics_INCLUDE_DIRS "/home/a11bmafi/Library/dune/dune-elastodynamics")
set(dune-elastodynamics_CXX_FLAGS "-std=c++17 ")
set(dune-elastodynamics_CXX_FLAGS_DEBUG "-g")
set(dune-elastodynamics_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(dune-elastodynamics_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(dune-elastodynamics_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(dune-elastodynamics_DEPENDS "dune-common;dune-geometry;dune-uggrid;dune-grid;dune-localfunctions;dune-istl;dune-functions;dune-foamgrid")
set(dune-elastodynamics_SUGGESTS "")
set(dune-elastodynamics_MODULE_PATH "/home/a11bmafi/Library/dune/dune-elastodynamics/cmake/modules")
set(dune-elastodynamics_LIBRARIES "")

# Lines that are set by the CMake build system via the variable DUNE_CUSTOM_PKG_CONFIG_SECTION


#import the target
if(dune-elastodynamics_LIBRARIES)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/dune-elastodynamics-targets.cmake")
endif()
endif()
