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
CMAKE_SOURCE_DIR = /home/ivan/OpticalFlow/MeshFlow_Video_Denoising

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build

# Include any dependencies generated for this target.
include source/CMakeFiles/SRC_MODULE_dynamic.dir/depend.make

# Include the progress variables for this target.
include source/CMakeFiles/SRC_MODULE_dynamic.dir/progress.make

# Include the compile flags for this target's objects.
include source/CMakeFiles/SRC_MODULE_dynamic.dir/flags.make

source/CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.o: source/CMakeFiles/SRC_MODULE_dynamic.dir/flags.make
source/CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.o: ../source/Fast_klt.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object source/CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.o"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.o -c /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/Fast_klt.cpp

source/CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.i"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/Fast_klt.cpp > CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.i

source/CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.s"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/Fast_klt.cpp -o CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.s

source/CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.o: source/CMakeFiles/SRC_MODULE_dynamic.dir/flags.make
source/CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.o: ../source/Mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object source/CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.o"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.o -c /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/Mesh.cpp

source/CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.i"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/Mesh.cpp > CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.i

source/CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.s"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/Mesh.cpp -o CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.s

source/CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.o: source/CMakeFiles/SRC_MODULE_dynamic.dir/flags.make
source/CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.o: ../source/MeshFlow.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object source/CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.o"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.o -c /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/MeshFlow.cpp

source/CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.i"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/MeshFlow.cpp > CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.i

source/CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.s"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/MeshFlow.cpp -o CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.s

source/CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.o: source/CMakeFiles/SRC_MODULE_dynamic.dir/flags.make
source/CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.o: ../source/MotionDenoiser.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object source/CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.o"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.o -c /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/MotionDenoiser.cpp

source/CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.i"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/MotionDenoiser.cpp > CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.i

source/CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.s"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/MotionDenoiser.cpp -o CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.s

source/CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.o: source/CMakeFiles/SRC_MODULE_dynamic.dir/flags.make
source/CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.o: ../source/VideoIO.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object source/CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.o"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.o -c /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/VideoIO.cpp

source/CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.i"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/VideoIO.cpp > CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.i

source/CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.s"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source/VideoIO.cpp -o CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.s

# Object files for target SRC_MODULE_dynamic
SRC_MODULE_dynamic_OBJECTS = \
"CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.o" \
"CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.o" \
"CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.o" \
"CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.o" \
"CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.o"

# External object files for target SRC_MODULE_dynamic
SRC_MODULE_dynamic_EXTERNAL_OBJECTS =

source/libSRC_MODULE.so.1.2: source/CMakeFiles/SRC_MODULE_dynamic.dir/Fast_klt.cpp.o
source/libSRC_MODULE.so.1.2: source/CMakeFiles/SRC_MODULE_dynamic.dir/Mesh.cpp.o
source/libSRC_MODULE.so.1.2: source/CMakeFiles/SRC_MODULE_dynamic.dir/MeshFlow.cpp.o
source/libSRC_MODULE.so.1.2: source/CMakeFiles/SRC_MODULE_dynamic.dir/MotionDenoiser.cpp.o
source/libSRC_MODULE.so.1.2: source/CMakeFiles/SRC_MODULE_dynamic.dir/VideoIO.cpp.o
source/libSRC_MODULE.so.1.2: source/CMakeFiles/SRC_MODULE_dynamic.dir/build.make
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_stitching.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_aruco.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_bgsegm.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_bioinspired.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_ccalib.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_dnn_objdetect.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_dnn_superres.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_dpm.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_face.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_freetype.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_fuzzy.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_hdf.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_hfs.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_img_hash.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_line_descriptor.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_quality.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_reg.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_rgbd.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_saliency.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_shape.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_stereo.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_structured_light.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_superres.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_surface_matching.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_tracking.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_videostab.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_viz.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_xobjdetect.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_xphoto.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_highgui.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_datasets.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_plot.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_text.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_dnn.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_ml.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_phase_unwrapping.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_optflow.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_ximgproc.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_video.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_videoio.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_imgcodecs.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_objdetect.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_calib3d.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_features2d.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_flann.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_photo.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_imgproc.so.4.2.0
source/libSRC_MODULE.so.1.2: /usr/lib/x86_64-linux-gnu/libopencv_core.so.4.2.0
source/libSRC_MODULE.so.1.2: source/CMakeFiles/SRC_MODULE_dynamic.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX shared library libSRC_MODULE.so"
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SRC_MODULE_dynamic.dir/link.txt --verbose=$(VERBOSE)
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && $(CMAKE_COMMAND) -E cmake_symlink_library libSRC_MODULE.so.1.2 libSRC_MODULE.so.1 libSRC_MODULE.so

source/libSRC_MODULE.so.1: source/libSRC_MODULE.so.1.2
	@$(CMAKE_COMMAND) -E touch_nocreate source/libSRC_MODULE.so.1

source/libSRC_MODULE.so: source/libSRC_MODULE.so.1.2
	@$(CMAKE_COMMAND) -E touch_nocreate source/libSRC_MODULE.so

# Rule to build all files generated by this target.
source/CMakeFiles/SRC_MODULE_dynamic.dir/build: source/libSRC_MODULE.so

.PHONY : source/CMakeFiles/SRC_MODULE_dynamic.dir/build

source/CMakeFiles/SRC_MODULE_dynamic.dir/clean:
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source && $(CMAKE_COMMAND) -P CMakeFiles/SRC_MODULE_dynamic.dir/cmake_clean.cmake
.PHONY : source/CMakeFiles/SRC_MODULE_dynamic.dir/clean

source/CMakeFiles/SRC_MODULE_dynamic.dir/depend:
	cd /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ivan/OpticalFlow/MeshFlow_Video_Denoising /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/source /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/build/source/CMakeFiles/SRC_MODULE_dynamic.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : source/CMakeFiles/SRC_MODULE_dynamic.dir/depend

