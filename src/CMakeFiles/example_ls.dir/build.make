# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

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
CMAKE_COMMAND = /home/guillermo/Programs/clion-2017.2.3/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/guillermo/Programs/clion-2017.2.3/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs"

# Include any dependencies generated for this target.
include CMakeFiles/example_ls.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/example_ls.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/example_ls.dir/flags.make

CMakeFiles/example_ls.dir/main.cpp.o: CMakeFiles/example_ls.dir/flags.make
CMakeFiles/example_ls.dir/main.cpp.o: main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/example_ls.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/example_ls.dir/main.cpp.o -c "/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs/main.cpp"

CMakeFiles/example_ls.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/example_ls.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs/main.cpp" > CMakeFiles/example_ls.dir/main.cpp.i

CMakeFiles/example_ls.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/example_ls.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs/main.cpp" -o CMakeFiles/example_ls.dir/main.cpp.s

CMakeFiles/example_ls.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/example_ls.dir/main.cpp.o.requires

CMakeFiles/example_ls.dir/main.cpp.o.provides: CMakeFiles/example_ls.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/example_ls.dir/build.make CMakeFiles/example_ls.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/example_ls.dir/main.cpp.o.provides

CMakeFiles/example_ls.dir/main.cpp.o.provides.build: CMakeFiles/example_ls.dir/main.cpp.o


# Object files for target example_ls
example_ls_OBJECTS = \
"CMakeFiles/example_ls.dir/main.cpp.o"

# External object files for target example_ls
example_ls_EXTERNAL_OBJECTS =

example_ls: CMakeFiles/example_ls.dir/main.cpp.o
example_ls: CMakeFiles/example_ls.dir/build.make
example_ls: liblocalsearch.so
example_ls: CMakeFiles/example_ls.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable example_ls"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example_ls.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/example_ls.dir/build: example_ls

.PHONY : CMakeFiles/example_ls.dir/build

CMakeFiles/example_ls.dir/requires: CMakeFiles/example_ls.dir/main.cpp.o.requires

.PHONY : CMakeFiles/example_ls.dir/requires

CMakeFiles/example_ls.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/example_ls.dir/cmake_clean.cmake
.PHONY : CMakeFiles/example_ls.dir/clean

CMakeFiles/example_ls.dir/depend:
	cd "/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs" "/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs" "/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs" "/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs" "/home/guillermo/UGR/MH/TRABAJO FINAL/INFO/localsearchs/CMakeFiles/example_ls.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/example_ls.dir/depend

