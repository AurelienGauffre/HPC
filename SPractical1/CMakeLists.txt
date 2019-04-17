cmake_minimum_required(VERSION 3.12)

project(PracticalSession1 CXX)

# Look for OpenMP in path and add it to compile options
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}
${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Force compiler to CXX if pass as cmake execute
if (CXX)
  set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Compile the various executable
add_executable(hello hello_omp.cxx)
add_executable(pi_Inegral pi_Integral_omp.cxx)
add_executable(pi_MonteCarlo pi_MonteCarlo_omp.cxx)
add_executable(mandel mandel_omp.cxx)




