# SPDX-License-Identifier: MIT
# Copyright (C) 2020 Benjamin Fovet

# Find the GMSH library and include directory
# This module defines the following variables :
#    Gmsh_INCLUDE_DIRS: where to find gmsh.
#    Gmsh_LIBRARIES: the gmsh library.
#    Gmsh_FOUND: if false, do not try to use Gmsh.

find_path(Gmsh_INCLUDE_DIR
  NAMES
    gmsh.h
  PATHS
    ${GMSH_DIR}
    ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES
    include)
mark_as_advanced(Gmsh_INCLUDE_DIR)

find_library(Gmsh_LIBRARY
  NAMES
    gmsh libgmsh
  PATHS
    ${GMSH_DIR}
    ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES
    lib lib64)
mark_as_advanced(Gmsh_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmsh
  REQUIRED_VARS Gmsh_LIBRARY Gmsh_INCLUDE_DIR)

if (Gmsh_FOUND)
  set(Gmsh_LIBRARIES "${Gmsh_LIRBARY}")
  set(Gmsh_INCLUDE_DIRS "${Gmsh_INCLUDE_DIR}")
  if (NOT TARGET Gmsh::Gmsh)
      add_library(Gmsh::Gmsh UNKNOWN IMPORTED)
      set_target_properties(Gmsh::Gmsh
        PROPERTIES
          IMPORTED_LOCATION "${Gmsh_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${Gmsh_INCLUDE_DIR}")
  endif()
endif()
