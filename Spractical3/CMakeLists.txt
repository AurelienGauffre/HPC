cmake_minimum_required(VERSION 2.8.11)

project(String CXX)

# Look for OpenMP in path and add it to compile options
find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_executable(string_omp_v1 string_equation_omp_v1.cxx)
add_executable(string_omp_v2 string_equation_omp_v2.cxx)
add_executable(string_omp_v3 string_equation_omp_v3.cxx)
add_executable(string_seq string_equation.cxx)
