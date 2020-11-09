include(FetchContent)
FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG yaml-cpp-0.6.3
  UPDATE_DISCONNECTED ON
  GIT_SHALLOW ON
)

FetchContent_GetProperties(yaml-cpp)
if(NOT yaml-cpp_POPULATED)
  message(STATUS "Populating yaml-cpp...")
  FetchContent_Populate(yaml-cpp)
  set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "disable yaml tests")
  add_subdirectory(${yaml-cpp_SOURCE_DIR} ${yaml-cpp_BINARY_DIR})
  message(STATUS "Done populating yaml-cpp.")
endif()
