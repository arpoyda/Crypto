# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

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
CMAKE_COMMAND = /home/arpoyda/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/211.7628.27/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/arpoyda/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/211.7628.27/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/arpoyda/Documents/mceliece

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/arpoyda/Documents/mceliece

# Include any dependencies generated for this target.
include CMakeFiles/McElieceCrypto.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/McElieceCrypto.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/McElieceCrypto.dir/flags.make

CMakeFiles/McElieceCrypto.dir/main.cpp.o: CMakeFiles/McElieceCrypto.dir/flags.make
CMakeFiles/McElieceCrypto.dir/main.cpp.o: main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/arpoyda/Documents/mceliece/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/McElieceCrypto.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/McElieceCrypto.dir/main.cpp.o -c /home/arpoyda/Documents/mceliece/main.cpp

CMakeFiles/McElieceCrypto.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/McElieceCrypto.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/arpoyda/Documents/mceliece/main.cpp > CMakeFiles/McElieceCrypto.dir/main.cpp.i

CMakeFiles/McElieceCrypto.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/McElieceCrypto.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/arpoyda/Documents/mceliece/main.cpp -o CMakeFiles/McElieceCrypto.dir/main.cpp.s

# Object files for target McElieceCrypto
McElieceCrypto_OBJECTS = \
"CMakeFiles/McElieceCrypto.dir/main.cpp.o"

# External object files for target McElieceCrypto
McElieceCrypto_EXTERNAL_OBJECTS =

McElieceCrypto: CMakeFiles/McElieceCrypto.dir/main.cpp.o
McElieceCrypto: CMakeFiles/McElieceCrypto.dir/build.make
McElieceCrypto: libMcElieceCryptoLib.so
McElieceCrypto: CMakeFiles/McElieceCrypto.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/arpoyda/Documents/mceliece/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable McElieceCrypto"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/McElieceCrypto.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/McElieceCrypto.dir/build: McElieceCrypto

.PHONY : CMakeFiles/McElieceCrypto.dir/build

CMakeFiles/McElieceCrypto.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/McElieceCrypto.dir/cmake_clean.cmake
.PHONY : CMakeFiles/McElieceCrypto.dir/clean

CMakeFiles/McElieceCrypto.dir/depend:
	cd /home/arpoyda/Documents/mceliece && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/arpoyda/Documents/mceliece /home/arpoyda/Documents/mceliece /home/arpoyda/Documents/mceliece /home/arpoyda/Documents/mceliece /home/arpoyda/Documents/mceliece/CMakeFiles/McElieceCrypto.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/McElieceCrypto.dir/depend

