# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.18.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.18.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/tomoleary1/Documents/projects/quEST/QuEST

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/tomoleary1/Documents/projects/quEST/QuEST

# Include any dependencies generated for this target.
include CMakeFiles/iam.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/iam.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/iam.dir/flags.make

CMakeFiles/iam.dir/vqe/iam.c.o: CMakeFiles/iam.dir/flags.make
CMakeFiles/iam.dir/vqe/iam.c.o: vqe/iam.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/tomoleary1/Documents/projects/quEST/QuEST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/iam.dir/vqe/iam.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/iam.dir/vqe/iam.c.o -c /Users/tomoleary1/Documents/projects/quEST/QuEST/vqe/iam.c

CMakeFiles/iam.dir/vqe/iam.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/iam.dir/vqe/iam.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/tomoleary1/Documents/projects/quEST/QuEST/vqe/iam.c > CMakeFiles/iam.dir/vqe/iam.c.i

CMakeFiles/iam.dir/vqe/iam.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/iam.dir/vqe/iam.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/tomoleary1/Documents/projects/quEST/QuEST/vqe/iam.c -o CMakeFiles/iam.dir/vqe/iam.c.s

# Object files for target iam
iam_OBJECTS = \
"CMakeFiles/iam.dir/vqe/iam.c.o"

# External object files for target iam
iam_EXTERNAL_OBJECTS =

iam: CMakeFiles/iam.dir/vqe/iam.c.o
iam: CMakeFiles/iam.dir/build.make
iam: QuEST/libQuEST.dylib
iam: CMakeFiles/iam.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/tomoleary1/Documents/projects/quEST/QuEST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable iam"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/iam.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/iam.dir/build: iam

.PHONY : CMakeFiles/iam.dir/build

CMakeFiles/iam.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/iam.dir/cmake_clean.cmake
.PHONY : CMakeFiles/iam.dir/clean

CMakeFiles/iam.dir/depend:
	cd /Users/tomoleary1/Documents/projects/quEST/QuEST && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/tomoleary1/Documents/projects/quEST/QuEST /Users/tomoleary1/Documents/projects/quEST/QuEST /Users/tomoleary1/Documents/projects/quEST/QuEST /Users/tomoleary1/Documents/projects/quEST/QuEST /Users/tomoleary1/Documents/projects/quEST/QuEST/CMakeFiles/iam.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/iam.dir/depend

