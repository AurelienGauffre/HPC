cmake_minimum_required(VERSION 2.8)

project(PracticalSession1)

# Look for OpenMP in path and add it to compile options
find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# Force compiler to CXX if pass as cmake execute
if (CXX)
  set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

# Compile the various executable
add_executable(hello hello.cxx)
add_executable(pi_Integral pi_Integral.cxx)
add_executable(pi_MonteCarlo pi_MonteCarlo.cxx)
add_executable(mandel mandel.cxx)
add_executable(max max.cxx)
