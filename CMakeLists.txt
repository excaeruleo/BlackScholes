cmake_minimum_required(VERSION 3.21)

project(
  BlackScholes
  VERSION 0.1.0
  DESCRIPTION "Simple Black-Scholes finite-difference example"
  HOMEPAGE_URL https://github.com/excaeruleo/BlackScholes/
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# output build artifacts in top level with per-config subdirectories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

# enable testing
include(CTest)

# compile options
if(MSVC)
  # note: need standards-compliant lambda processor
  add_compile_options(/W4 /Zc:lambda)
else()
  add_compile_options(-Wall)
endif()

# top-level include directory
set(EXCBS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory(src)
