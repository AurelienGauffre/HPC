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
include CMakeFiles/string_seq.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/string_seq.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/string_seq.dir/flags.make

CMakeFiles/string_seq.dir/string_equation.cxx.o: CMakeFiles/string_seq.dir/flags.make
CMakeFiles/string_seq.dir/string_equation.cxx.o: ../string_equation.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/string_seq.dir/string_equation.cxx.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/string_seq.dir/string_equation.cxx.o -c "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/string_equation.cxx"

CMakeFiles/string_seq.dir/string_equation.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/string_seq.dir/string_equation.cxx.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/string_equation.cxx" > CMakeFiles/string_seq.dir/string_equation.cxx.i

CMakeFiles/string_seq.dir/string_equation.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/string_seq.dir/string_equation.cxx.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/string_equation.cxx" -o CMakeFiles/string_seq.dir/string_equation.cxx.s

CMakeFiles/string_seq.dir/string_equation.cxx.o.requires:

.PHONY : CMakeFiles/string_seq.dir/string_equation.cxx.o.requires

CMakeFiles/string_seq.dir/string_equation.cxx.o.provides: CMakeFiles/string_seq.dir/string_equation.cxx.o.requires
	$(MAKE) -f CMakeFiles/string_seq.dir/build.make CMakeFiles/string_seq.dir/string_equation.cxx.o.provides.build
.PHONY : CMakeFiles/string_seq.dir/string_equation.cxx.o.provides

CMakeFiles/string_seq.dir/string_equation.cxx.o.provides.build: CMakeFiles/string_seq.dir/string_equation.cxx.o


# Object files for target string_seq
string_seq_OBJECTS = \
"CMakeFiles/string_seq.dir/string_equation.cxx.o"

# External object files for target string_seq
string_seq_EXTERNAL_OBJECTS =

string_seq: CMakeFiles/string_seq.dir/string_equation.cxx.o
string_seq: CMakeFiles/string_seq.dir/build.make
string_seq: CMakeFiles/string_seq.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable string_seq"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/string_seq.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/string_seq.dir/build: string_seq

.PHONY : CMakeFiles/string_seq.dir/build

CMakeFiles/string_seq.dir/requires: CMakeFiles/string_seq.dir/string_equation.cxx.o.requires

.PHONY : CMakeFiles/string_seq.dir/requires

CMakeFiles/string_seq.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/string_seq.dir/cmake_clean.cmake
.PHONY : CMakeFiles/string_seq.dir/clean

CMakeFiles/string_seq.dir/depend:
	cd "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction" "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction" "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build" "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build" "/mnt/c/Users/Aurélien/Desktop/ENSIMAG/HPC/Practical3/String Correction/build/CMakeFiles/string_seq.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/string_seq.dir/depend
