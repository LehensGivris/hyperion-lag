# SPDX-License-Identifier: MIT
# Copyright (C) 2020 Benjamin Fovet

cmake_minimum_required(VERSION 3.10)

## Project initialization ######################################################
project(HyPERION
  VERSION 0.1.0
  DESCRIPTION "Hydrodynamics Platform for Exascale Research, In situ analysis and OptimizatioN"
  HOMEPAGE_URL https://github.com/bfovet/hyperion-lag
  LANGUAGES CXX)

# CMake modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)
enable_testing()

## System configuration checks #################################################

## Code configuration ##########################################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Directory containing sources
add_subdirectory(src)

# Directory containing unit tests
#add_subdirectory(test)
