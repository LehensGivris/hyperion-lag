cmake_minimum_required(VERSION 3.10)

## Project initialization ######################################################
project(HyPERION
  VERSION 0.1.0
  DESCRIPTION "Hydrodynamics Platform for Exascale Research, In situ analysis and OptimizatioN"
  HOMEPAGE_URL https://github.com/bfovet/hyperion-lag
  LANGUAGES C CXX)

# CMake modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# Enable In Situ analysis with ParaView Catalyst
option(HYPERION_ENABLE_INSITU "Enable In Situ analysis with ParaView Catalyst" ON)

## System configuration checks #################################################

# Find GMSH
find_package(Gmsh)

# Find ParaView for Catalyst
if (HYPERION_ENABLE_INSITU)
  find_package(ParaView 5.8 REQUIRED COMPONENTS Catalyst PythonCatalyst)
endif()

# Find VTK
find_package(VTK COMPONENTS
  CommonCore
  CommonDataModel
  FiltersGeneral
  FiltersExtraction
  FiltersGeometry
  IOXML
  IOLegacy
  REQUIRED)

# Fetch yaml-cpp
include(${CMAKE_MODULE_PATH}/FetchYamlCpp.cmake)

## Code configuration ##########################################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Directory containing sources
add_subdirectory(src)
