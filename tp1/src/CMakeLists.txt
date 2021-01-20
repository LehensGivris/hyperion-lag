# HypERION main library
add_library(h2p_lib)
set_target_properties(h2p_lib PROPERTIES OUTPUT_NAME h2p)

target_compile_features(h2p_lib PUBLIC cxx_std_17)
target_compile_options(h2p_lib PUBLIC
  $<$<CXX_COMPILER_ID:GNU>:-Wall;-Wpedantic;-Wextra;-fexceptions>)

target_sources(h2p_lib
  PRIVATE
  HydroVars.cpp
  Hydro.cpp
  HyperionMainDriver.cpp)

target_link_libraries(h2p_lib
  PUBLIC
  Gmsh::Gmsh
  yaml-cpp
  ${VTK_LIBRARIES}
  )

# HyPERION executable
add_executable(h2p_exe main.cpp)
set_target_properties(h2p_exe PROPERTIES OUTPUT_NAME h2p)

target_compile_features(h2p_exe PUBLIC cxx_std_17)
target_compile_options(h2p_exe PUBLIC
  $<$<CXX_COMPILER_ID:GNU>:-Wall;-Wpedantic;-Wextra;-fexceptions>)

target_link_libraries(h2p_exe
  PRIVATE
  h2p_lib
  )
