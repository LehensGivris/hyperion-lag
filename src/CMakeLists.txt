# SPDX-License-Identifier: MIT
# Copyright (C) 2020 Benjamin Fovet

# HypERION main library
add_library(h2p_lib)
set_target_properties(h2p_lib PROPERTIES OUTPUT_NAME h2p)

target_compile_features(h2p_lib PUBLIC cxx_std_17)
target_compile_options(h2p_lib PUBLIC
  $<$<CXX_COMPILER_ID:GNU>:-Wall;-Wpedantic;-Wextra;-fexceptions>)

if (HYPERION_ENABLE_INSITU)
  target_compile_definitions(h2p_lib
    PUBLIC
    ANALYZE_INSITU)
endif()

set(H2P_SOURCES
  HydroVars.cpp
  Hydro.cpp
  HyperionMainDriver.cpp)

if (HYPERION_ENABLE_INSITU)
  list(APPEND H2P_SOURCES
    PRIVATE
    CatalystInsituAdaptor.cpp)
endif()

target_sources(h2p_lib
  PRIVATE
  ${H2P_SOURCES})

target_link_libraries(h2p_lib
  PUBLIC
  Gmsh::Gmsh
  yaml-cpp
  VTK::CommonDataModel
  VTK::CommonCore
  VTK::FiltersGeneral
  VTK::FiltersExtraction
  VTK::FiltersGeometry
  VTK::IOXML
  VTK::IOLegacy)

if (HYPERION_ENABLE_INSITU)
  target_link_libraries(h2p_lib
    PUBLIC
    ParaView::Catalyst
    ParaView::PythonCatalyst)
endif()

# HyPERION executable
add_executable(h2p_exe main.cpp)
set_target_properties(h2p_exe PROPERTIES OUTPUT_NAME h2p)

target_compile_features(h2p_exe PUBLIC cxx_std_17)
target_compile_options(h2p_exe PUBLIC
  $<$<CXX_COMPILER_ID:GNU>:-Wall;-Wpedantic;-Wextra;-fexceptions>)

if (HYPERION_ENABLE_INSITU)
  target_compile_definitions(h2p_exe
    PUBLIC
    ANALYZE_INSITU)
endif()

target_link_libraries(h2p_exe
  PRIVATE
  h2p_lib
  )
