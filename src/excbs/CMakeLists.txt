cmake_minimum_required(VERSION ${CMAKE_MINIMUM_REQUIRED_VERSION})

# excbs: common support library
add_library(excbs STATIC solver.cpp)
target_include_directories(excbs PUBLIC "${EXCBS_INCLUDE_DIR}")
