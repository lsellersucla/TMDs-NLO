# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 4.0

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
CMAKE_COMMAND = /opt/homebrew/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/luke/Desktop/TMDs-NLO/lib/apfel

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/luke/desktop/TMDs-NLO/lib/apfel

# Include any dependencies generated for this target.
include pywrap/CMakeFiles/_apfel.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include pywrap/CMakeFiles/_apfel.dir/compiler_depend.make

# Include the progress variables for this target.
include pywrap/CMakeFiles/_apfel.dir/progress.make

# Include the compile flags for this target's objects.
include pywrap/CMakeFiles/_apfel.dir/flags.make

pywrap/CMakeFiles/_apfel.dir/codegen:
.PHONY : pywrap/CMakeFiles/_apfel.dir/codegen

pywrap/CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.o: pywrap/CMakeFiles/_apfel.dir/flags.make
pywrap/CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.o: pywrap/apfelPYTHON_wrap.cxx
pywrap/CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.o: pywrap/CMakeFiles/_apfel.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/luke/desktop/TMDs-NLO/lib/apfel/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object pywrap/CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.o"
	cd /Users/luke/desktop/TMDs-NLO/lib/apfel/pywrap && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT pywrap/CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.o -MF CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.o.d -o CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.o -c /Users/luke/desktop/TMDs-NLO/lib/apfel/pywrap/apfelPYTHON_wrap.cxx

pywrap/CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.i"
	cd /Users/luke/desktop/TMDs-NLO/lib/apfel/pywrap && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/luke/desktop/TMDs-NLO/lib/apfel/pywrap/apfelPYTHON_wrap.cxx > CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.i

pywrap/CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.s"
	cd /Users/luke/desktop/TMDs-NLO/lib/apfel/pywrap && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/luke/desktop/TMDs-NLO/lib/apfel/pywrap/apfelPYTHON_wrap.cxx -o CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.s

# Object files for target _apfel
_apfel_OBJECTS = \
"CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.o"

# External object files for target _apfel
_apfel_EXTERNAL_OBJECTS =

pywrap/_apfel.so: pywrap/CMakeFiles/_apfel.dir/apfelPYTHON_wrap.cxx.o
pywrap/_apfel.so: pywrap/CMakeFiles/_apfel.dir/build.make
pywrap/_apfel.so: lib/libAPFEL.0.0.0.dylib
pywrap/_apfel.so: /opt/anaconda3/envs/ml/lib/libLHAPDF.dylib
pywrap/_apfel.so: pywrap/CMakeFiles/_apfel.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/luke/desktop/TMDs-NLO/lib/apfel/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module _apfel.so"
	cd /Users/luke/desktop/TMDs-NLO/lib/apfel/pywrap && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/_apfel.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
pywrap/CMakeFiles/_apfel.dir/build: pywrap/_apfel.so
.PHONY : pywrap/CMakeFiles/_apfel.dir/build

pywrap/CMakeFiles/_apfel.dir/clean:
	cd /Users/luke/desktop/TMDs-NLO/lib/apfel/pywrap && $(CMAKE_COMMAND) -P CMakeFiles/_apfel.dir/cmake_clean.cmake
.PHONY : pywrap/CMakeFiles/_apfel.dir/clean

pywrap/CMakeFiles/_apfel.dir/depend:
	cd /Users/luke/desktop/TMDs-NLO/lib/apfel && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/luke/Desktop/TMDs-NLO/lib/apfel /Users/luke/Desktop/TMDs-NLO/lib/apfel/pywrap /Users/luke/desktop/TMDs-NLO/lib/apfel /Users/luke/desktop/TMDs-NLO/lib/apfel/pywrap /Users/luke/desktop/TMDs-NLO/lib/apfel/pywrap/CMakeFiles/_apfel.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : pywrap/CMakeFiles/_apfel.dir/depend

