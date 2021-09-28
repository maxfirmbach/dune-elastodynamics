# Install script for directory: /home/a11bmafi/Library/dune/dune-elastodynamics

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  set(CMAKE_MODULE_PATH /home/a11bmafi/Library/dune/dune-elastodynamics/cmake/modules;/home/a11bmafi/Library/dune/dune-functions/cmake/modules;/home/a11bmafi/Library/dune/dune-foamgrid/cmake/modules;/home/a11bmafi/Library/dune/dune-localfunctions/cmake/modules;/home/a11bmafi/Library/dune/dune-grid/cmake/modules;/home/a11bmafi/Library/dune/dune-istl/cmake/modules;/home/a11bmafi/Library/dune/dune-typetree/cmake/modules;/home/a11bmafi/Library/dune/dune-geometry/cmake/modules;/home/a11bmafi/Library/dune/dune-uggrid/cmake/modules;/home/a11bmafi/Library/dune/dune-common/cmake/modules;/usr/lib/cmake/Vc)
              set(DUNE_PYTHON_WHEELHOUSE /usr/local/share/dune/wheelhouse)
              include(DuneExecuteProcess)
              dune_execute_process(COMMAND "/usr/bin/cmake" --build . --target install_python --config $<CONFIG>)
              
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/dunecontrol/dune-elastodynamics" TYPE FILE FILES "/home/a11bmafi/Library/dune/dune-elastodynamics/dune.module")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/dune-elastodynamics" TYPE FILE FILES
    "/home/a11bmafi/Library/dune/dune-elastodynamics/build-cmake/cmake/pkg/dune-elastodynamics-config.cmake"
    "/home/a11bmafi/Library/dune/dune-elastodynamics/build-cmake/dune-elastodynamics-config-version.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/dune-elastodynamics" TYPE FILE FILES "/home/a11bmafi/Library/dune/dune-elastodynamics/config.h.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/a11bmafi/Library/dune/dune-elastodynamics/build-cmake/dune-elastodynamics.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/a11bmafi/Library/dune/dune-elastodynamics/build-cmake/src/cmake_install.cmake")
  include("/home/a11bmafi/Library/dune/dune-elastodynamics/build-cmake/dune/cmake_install.cmake")
  include("/home/a11bmafi/Library/dune/dune-elastodynamics/build-cmake/cmake/modules/cmake_install.cmake")
  include("/home/a11bmafi/Library/dune/dune-elastodynamics/build-cmake/test/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/a11bmafi/Library/dune/dune-elastodynamics/build-cmake/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
