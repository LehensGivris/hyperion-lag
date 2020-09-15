# SPDX-License-Identifier: MIT
# Copyright (C) 2020 Benjamin Fovet

add_executable(h2p_exe main.cpp)
set_target_properties(h2p_exe PROPERTIES OUTPUT_NAME h2p)

target_compile_features(h2p_exe PUBLIC cxx_std_17)
target_compile_options(h2p_exe PUBLIC
  $<$<CXX_COMPILER_ID:GNU>:-Wall;-Wpedantic;-Wextra;-fexceptions>)
