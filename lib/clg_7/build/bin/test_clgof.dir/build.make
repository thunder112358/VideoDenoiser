# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/ivan/yasser/VideoDenoiser/lib/clg_7

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ivan/yasser/VideoDenoiser/lib/clg_7/build

# Include any dependencies generated for this target.
include CMakeFiles/../bin/test_clgof.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/../bin/test_clgof.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/../bin/test_clgof.dir/flags.make

CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.o: CMakeFiles/../bin/test_clgof.dir/flags.make
CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.o: ../source/clg_of_main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/lib/clg_7/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.o   -c /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/clg_of_main.c

CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/clg_of_main.c > CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.i

CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/clg_of_main.c -o CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.s

CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.o: CMakeFiles/../bin/test_clgof.dir/flags.make
CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.o: ../source/clg_of.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/lib/clg_7/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.o   -c /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/clg_of.c

CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/clg_of.c > CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.i

CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/clg_of.c -o CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.s

CMakeFiles/../bin/test_clgof.dir/source/iio.c.o: CMakeFiles/../bin/test_clgof.dir/flags.make
CMakeFiles/../bin/test_clgof.dir/source/iio.c.o: ../source/iio.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/lib/clg_7/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/../bin/test_clgof.dir/source/iio.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/../bin/test_clgof.dir/source/iio.c.o   -c /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/iio.c

CMakeFiles/../bin/test_clgof.dir/source/iio.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/../bin/test_clgof.dir/source/iio.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/iio.c > CMakeFiles/../bin/test_clgof.dir/source/iio.c.i

CMakeFiles/../bin/test_clgof.dir/source/iio.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/../bin/test_clgof.dir/source/iio.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/iio.c -o CMakeFiles/../bin/test_clgof.dir/source/iio.c.s

CMakeFiles/../bin/test_clgof.dir/source/util.c.o: CMakeFiles/../bin/test_clgof.dir/flags.make
CMakeFiles/../bin/test_clgof.dir/source/util.c.o: ../source/util.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/lib/clg_7/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/../bin/test_clgof.dir/source/util.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/../bin/test_clgof.dir/source/util.c.o   -c /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/util.c

CMakeFiles/../bin/test_clgof.dir/source/util.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/../bin/test_clgof.dir/source/util.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/util.c > CMakeFiles/../bin/test_clgof.dir/source/util.c.i

CMakeFiles/../bin/test_clgof.dir/source/util.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/../bin/test_clgof.dir/source/util.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/util.c -o CMakeFiles/../bin/test_clgof.dir/source/util.c.s

CMakeFiles/../bin/test_clgof.dir/source/mask.c.o: CMakeFiles/../bin/test_clgof.dir/flags.make
CMakeFiles/../bin/test_clgof.dir/source/mask.c.o: ../source/mask.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/lib/clg_7/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/../bin/test_clgof.dir/source/mask.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/../bin/test_clgof.dir/source/mask.c.o   -c /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/mask.c

CMakeFiles/../bin/test_clgof.dir/source/mask.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/../bin/test_clgof.dir/source/mask.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/mask.c > CMakeFiles/../bin/test_clgof.dir/source/mask.c.i

CMakeFiles/../bin/test_clgof.dir/source/mask.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/../bin/test_clgof.dir/source/mask.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/mask.c -o CMakeFiles/../bin/test_clgof.dir/source/mask.c.s

CMakeFiles/../bin/test_clgof.dir/source/zoom.c.o: CMakeFiles/../bin/test_clgof.dir/flags.make
CMakeFiles/../bin/test_clgof.dir/source/zoom.c.o: ../source/zoom.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/lib/clg_7/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object CMakeFiles/../bin/test_clgof.dir/source/zoom.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/../bin/test_clgof.dir/source/zoom.c.o   -c /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/zoom.c

CMakeFiles/../bin/test_clgof.dir/source/zoom.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/../bin/test_clgof.dir/source/zoom.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/zoom.c > CMakeFiles/../bin/test_clgof.dir/source/zoom.c.i

CMakeFiles/../bin/test_clgof.dir/source/zoom.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/../bin/test_clgof.dir/source/zoom.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/zoom.c -o CMakeFiles/../bin/test_clgof.dir/source/zoom.c.s

CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.o: CMakeFiles/../bin/test_clgof.dir/flags.make
CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.o: ../source/bicubic_interpolation.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/lib/clg_7/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.o   -c /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/bicubic_interpolation.c

CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/bicubic_interpolation.c > CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.i

CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ivan/yasser/VideoDenoiser/lib/clg_7/source/bicubic_interpolation.c -o CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.s

# Object files for target ../bin/test_clgof
__/bin/test_clgof_OBJECTS = \
"CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.o" \
"CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.o" \
"CMakeFiles/../bin/test_clgof.dir/source/iio.c.o" \
"CMakeFiles/../bin/test_clgof.dir/source/util.c.o" \
"CMakeFiles/../bin/test_clgof.dir/source/mask.c.o" \
"CMakeFiles/../bin/test_clgof.dir/source/zoom.c.o" \
"CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.o"

# External object files for target ../bin/test_clgof
__/bin/test_clgof_EXTERNAL_OBJECTS =

../bin/test_clgof: CMakeFiles/../bin/test_clgof.dir/source/clg_of_main.c.o
../bin/test_clgof: CMakeFiles/../bin/test_clgof.dir/source/clg_of.c.o
../bin/test_clgof: CMakeFiles/../bin/test_clgof.dir/source/iio.c.o
../bin/test_clgof: CMakeFiles/../bin/test_clgof.dir/source/util.c.o
../bin/test_clgof: CMakeFiles/../bin/test_clgof.dir/source/mask.c.o
../bin/test_clgof: CMakeFiles/../bin/test_clgof.dir/source/zoom.c.o
../bin/test_clgof: CMakeFiles/../bin/test_clgof.dir/source/bicubic_interpolation.c.o
../bin/test_clgof: CMakeFiles/../bin/test_clgof.dir/build.make
../bin/test_clgof: CMakeFiles/../bin/test_clgof.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ivan/yasser/VideoDenoiser/lib/clg_7/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking C executable ../bin/test_clgof"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/../bin/test_clgof.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/../bin/test_clgof.dir/build: ../bin/test_clgof

.PHONY : CMakeFiles/../bin/test_clgof.dir/build

CMakeFiles/../bin/test_clgof.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/../bin/test_clgof.dir/cmake_clean.cmake
.PHONY : CMakeFiles/../bin/test_clgof.dir/clean

CMakeFiles/../bin/test_clgof.dir/depend:
	cd /home/ivan/yasser/VideoDenoiser/lib/clg_7/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ivan/yasser/VideoDenoiser/lib/clg_7 /home/ivan/yasser/VideoDenoiser/lib/clg_7 /home/ivan/yasser/VideoDenoiser/lib/clg_7/build /home/ivan/yasser/VideoDenoiser/lib/clg_7/build /home/ivan/yasser/VideoDenoiser/lib/clg_7/build/bin/test_clgof.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/../bin/test_clgof.dir/depend

