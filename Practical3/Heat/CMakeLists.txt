cmake_minimum_required(VERSION 2.8.11)

project(Heat)

# Look for OpenMP in path and add it to compile options
# find_package(OpenMP)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/*.cxx)
add_executable(heat ${SRC_FILES})
