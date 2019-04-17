cmake_minimum_required(VERSION 2.8)


set(CMAKE_C_COMPILER "/usr/local/opt/llvm/bin/clang" CACHE INTERNAL "")
set(CMAKE_CXX_COMPILER "/usr/local/opt/llvm/bin/clang++" CACHE INTERNAL "")
set(OPENMP_LIBRARIES "/usr/local/opt/llvm/lib")
set(OPENMP_INCLUDES "/usr/local/opt/llvm/include")

project(Practical2)

# Force compiler to CXX if pass as cmake execute
if (CXX)
  set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

# Look for OpenMP in path and add it to compile options
message(STATUS "Checking OpenMP")
if(APPLE)
     if(CMAKE_C_COMPILER_ID MATCHES "Clang")
         set(OpenMP_C "${CMAKE_C_COMPILER}")
         set(OpenMP_C_FLAGS "-fopenmp")
         set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
         set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
         set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
         set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
     endif()
     if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
       set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
       set(OpenMP_CXX_FLAGS "-fopenmp")
       set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
       set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
       set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
       set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
     endif()
 endif()

 find_package(OpenMP REQUIRED)
 if (OPENMP_FOUND)
     include_directories("${OPENMP_INCLUDES}")
     link_directories("${OPENMP_LIBRARIES}")
     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
 endif(OPENMP_FOUND)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})




# Compile the various executable
add_executable(max_subarray max_subarray_omp.cxx)
target_link_libraries(max_subarray Eigen3::Eigen)
# add_executable(matrix_product matrix_product.cxx)
# target_link_libraries(matrix_product Eigen3::Eigen)
