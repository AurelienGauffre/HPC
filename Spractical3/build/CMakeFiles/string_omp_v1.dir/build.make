# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build"

# Include any dependencies generated for this target.
include CMakeFiles/string_omp_v1.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/string_omp_v1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/string_omp_v1.dir/flags.make

CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o: CMakeFiles/string_omp_v1.dir/flags.make
CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o: ../string_equation_omp_v1.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o -c "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/string_equation_omp_v1.cxx"

CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/string_equation_omp_v1.cxx" > CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.i

CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/string_equation_omp_v1.cxx" -o CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.s

CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o.requires:

.PHONY : CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o.requires

CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o.provides: CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o.requires
	$(MAKE) -f CMakeFiles/string_omp_v1.dir/build.make CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o.provides.build
.PHONY : CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o.provides

CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o.provides.build: CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o


# Object files for target string_omp_v1
string_omp_v1_OBJECTS = \
"CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o"

# External object files for target string_omp_v1
string_omp_v1_EXTERNAL_OBJECTS =

string_omp_v1: CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o
string_omp_v1: CMakeFiles/string_omp_v1.dir/build.make
string_omp_v1: CMakeFiles/string_omp_v1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable string_omp_v1"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/string_omp_v1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/string_omp_v1.dir/build: string_omp_v1

.PHONY : CMakeFiles/string_omp_v1.dir/build

CMakeFiles/string_omp_v1.dir/requires: CMakeFiles/string_omp_v1.dir/string_equation_omp_v1.cxx.o.requires

.PHONY : CMakeFiles/string_omp_v1.dir/requires

CMakeFiles/string_omp_v1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/string_omp_v1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/string_omp_v1.dir/clean

CMakeFiles/string_omp_v1.dir/depend:
	cd "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction" "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction" "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build" "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build" "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build/CMakeFiles/string_omp_v1.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/string_omp_v1.dir/depend

