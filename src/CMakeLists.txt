cmake_minimum_required(VERSION ${CMAKE_MINIMUM_REQUIRED_VERSION})

# add excbs support library
add_subdirectory(excbs)

# bs_compare: Black-Scholes Euler forward difference vs. analytical comparison
add_executable(bs_compare bs_compare.cpp)
target_link_libraries(bs_compare PRIVATE excbs)
# simple test to run the program
add_test(NAME bs_compare COMMAND bs_compare)
