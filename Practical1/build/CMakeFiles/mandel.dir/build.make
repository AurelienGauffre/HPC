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
CMAKE_SOURCE_DIR = /mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1/build

# Include any dependencies generated for this target.
include CMakeFiles/mandel.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mandel.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mandel.dir/flags.make

CMakeFiles/mandel.dir/mandel.cxx.o: CMakeFiles/mandel.dir/flags.make
CMakeFiles/mandel.dir/mandel.cxx.o: ../mandel.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mandel.dir/mandel.cxx.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mandel.dir/mandel.cxx.o -c /mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1/mandel.cxx

CMakeFiles/mandel.dir/mandel.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mandel.dir/mandel.cxx.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1/mandel.cxx > CMakeFiles/mandel.dir/mandel.cxx.i

CMakeFiles/mandel.dir/mandel.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mandel.dir/mandel.cxx.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1/mandel.cxx -o CMakeFiles/mandel.dir/mandel.cxx.s

CMakeFiles/mandel.dir/mandel.cxx.o.requires:

.PHONY : CMakeFiles/mandel.dir/mandel.cxx.o.requires

CMakeFiles/mandel.dir/mandel.cxx.o.provides: CMakeFiles/mandel.dir/mandel.cxx.o.requires
	$(MAKE) -f CMakeFiles/mandel.dir/build.make CMakeFiles/mandel.dir/mandel.cxx.o.provides.build
.PHONY : CMakeFiles/mandel.dir/mandel.cxx.o.provides

CMakeFiles/mandel.dir/mandel.cxx.o.provides.build: CMakeFiles/mandel.dir/mandel.cxx.o


# Object files for target mandel
mandel_OBJECTS = \
"CMakeFiles/mandel.dir/mandel.cxx.o"

# External object files for target mandel
mandel_EXTERNAL_OBJECTS =

mandel: CMakeFiles/mandel.dir/mandel.cxx.o
mandel: CMakeFiles/mandel.dir/build.make
mandel: CMakeFiles/mandel.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mandel"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mandel.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mandel.dir/build: mandel

.PHONY : CMakeFiles/mandel.dir/build

CMakeFiles/mandel.dir/requires: CMakeFiles/mandel.dir/mandel.cxx.o.requires

.PHONY : CMakeFiles/mandel.dir/requires

CMakeFiles/mandel.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mandel.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mandel.dir/clean

CMakeFiles/mandel.dir/depend:
	cd /mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1 /mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1 /mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1/build /mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1/build /mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical1/build/CMakeFiles/mandel.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mandel.dir/depend

