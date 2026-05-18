cmake_minimum_required(VERSION 3.21)

project(
  BlackScholes
  VERSION 0.1.0
  DESCRIPTION "Simple Black-Scholes finite-difference example"
  HOMEPAGE_URL https://github.com/excaeruleo/BlackScholes/
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# compile options
if(MSVC)
  # note: need standards-compliant lambda processor
  add_compile_options(/W4 /Zc:lambda)
else()
  add_compile_options(-Wall)
endif()

# enable testing
include(CTest)

# bs_compare: Black-Scholes Euler forward difference vs. analytical comparison
add_executable(bs_compare black_scholes_solver.cpp)
add_test(NAME bs_compare COMMAND bs_compare)
