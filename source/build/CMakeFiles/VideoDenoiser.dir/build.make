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
CMAKE_SOURCE_DIR = /home/ivan/yasser/VideoDenoiser/source

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ivan/yasser/VideoDenoiser/source/build

# Include any dependencies generated for this target.
include CMakeFiles/VideoDenoiser.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/VideoDenoiser.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/VideoDenoiser.dir/flags.make

CMakeFiles/VideoDenoiser.dir/main.cpp.o: CMakeFiles/VideoDenoiser.dir/flags.make
CMakeFiles/VideoDenoiser.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/VideoDenoiser.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VideoDenoiser.dir/main.cpp.o -c /home/ivan/yasser/VideoDenoiser/source/main.cpp

CMakeFiles/VideoDenoiser.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VideoDenoiser.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/yasser/VideoDenoiser/source/main.cpp > CMakeFiles/VideoDenoiser.dir/main.cpp.i

CMakeFiles/VideoDenoiser.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VideoDenoiser.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/yasser/VideoDenoiser/source/main.cpp -o CMakeFiles/VideoDenoiser.dir/main.cpp.s

CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.o: CMakeFiles/VideoDenoiser.dir/flags.make
CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.o: ../DenseOF.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.o -c /home/ivan/yasser/VideoDenoiser/source/DenseOF.cpp

CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/yasser/VideoDenoiser/source/DenseOF.cpp > CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.i

CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/yasser/VideoDenoiser/source/DenseOF.cpp -o CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.s

CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.o: CMakeFiles/VideoDenoiser.dir/flags.make
CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.o: ../Fast_klt.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.o -c /home/ivan/yasser/VideoDenoiser/source/Fast_klt.cpp

CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/yasser/VideoDenoiser/source/Fast_klt.cpp > CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.i

CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/yasser/VideoDenoiser/source/Fast_klt.cpp -o CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.s

CMakeFiles/VideoDenoiser.dir/Mesh.cpp.o: CMakeFiles/VideoDenoiser.dir/flags.make
CMakeFiles/VideoDenoiser.dir/Mesh.cpp.o: ../Mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/VideoDenoiser.dir/Mesh.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VideoDenoiser.dir/Mesh.cpp.o -c /home/ivan/yasser/VideoDenoiser/source/Mesh.cpp

CMakeFiles/VideoDenoiser.dir/Mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VideoDenoiser.dir/Mesh.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/yasser/VideoDenoiser/source/Mesh.cpp > CMakeFiles/VideoDenoiser.dir/Mesh.cpp.i

CMakeFiles/VideoDenoiser.dir/Mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VideoDenoiser.dir/Mesh.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/yasser/VideoDenoiser/source/Mesh.cpp -o CMakeFiles/VideoDenoiser.dir/Mesh.cpp.s

CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.o: CMakeFiles/VideoDenoiser.dir/flags.make
CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.o: ../MeshFlow.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.o -c /home/ivan/yasser/VideoDenoiser/source/MeshFlow.cpp

CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/yasser/VideoDenoiser/source/MeshFlow.cpp > CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.i

CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/yasser/VideoDenoiser/source/MeshFlow.cpp -o CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.s

CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.o: CMakeFiles/VideoDenoiser.dir/flags.make
CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.o: ../MotionDenoiser.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.o -c /home/ivan/yasser/VideoDenoiser/source/MotionDenoiser.cpp

CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/yasser/VideoDenoiser/source/MotionDenoiser.cpp > CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.i

CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/yasser/VideoDenoiser/source/MotionDenoiser.cpp -o CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.s

CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.o: CMakeFiles/VideoDenoiser.dir/flags.make
CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.o: ../VideoIO.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ivan/yasser/VideoDenoiser/source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.o -c /home/ivan/yasser/VideoDenoiser/source/VideoIO.cpp

CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ivan/yasser/VideoDenoiser/source/VideoIO.cpp > CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.i

CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ivan/yasser/VideoDenoiser/source/VideoIO.cpp -o CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.s

# Object files for target VideoDenoiser
VideoDenoiser_OBJECTS = \
"CMakeFiles/VideoDenoiser.dir/main.cpp.o" \
"CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.o" \
"CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.o" \
"CMakeFiles/VideoDenoiser.dir/Mesh.cpp.o" \
"CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.o" \
"CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.o" \
"CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.o"

# External object files for target VideoDenoiser
VideoDenoiser_EXTERNAL_OBJECTS =

VideoDenoiser: CMakeFiles/VideoDenoiser.dir/main.cpp.o
VideoDenoiser: CMakeFiles/VideoDenoiser.dir/DenseOF.cpp.o
VideoDenoiser: CMakeFiles/VideoDenoiser.dir/Fast_klt.cpp.o
VideoDenoiser: CMakeFiles/VideoDenoiser.dir/Mesh.cpp.o
VideoDenoiser: CMakeFiles/VideoDenoiser.dir/MeshFlow.cpp.o
VideoDenoiser: CMakeFiles/VideoDenoiser.dir/MotionDenoiser.cpp.o
VideoDenoiser: CMakeFiles/VideoDenoiser.dir/VideoIO.cpp.o
VideoDenoiser: CMakeFiles/VideoDenoiser.dir/build.make
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_calib3d.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_core.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_dnn.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_features2d.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_flann.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_gapi.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_highgui.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_imgcodecs.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_imgproc.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_ml.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_objdetect.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_photo.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_stitching.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_video.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_videoio.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_alphamat.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_aruco.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_bgsegm.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_bioinspired.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_ccalib.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_datasets.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_dnn_objdetect.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_dnn_superres.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_dpm.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_face.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_freetype.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_fuzzy.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_hdf.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_hfs.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_img_hash.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_intensity_transform.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_line_descriptor.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_mcc.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_optflow.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_phase_unwrapping.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_plot.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_quality.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_rapid.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_reg.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_rgbd.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_saliency.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_sfm.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_shape.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_signal.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_stereo.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_structured_light.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_superres.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_surface_matching.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_text.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_tracking.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_videostab.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_viz.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_wechat_qrcode.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_xfeatures2d.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_ximgproc.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_xobjdetect.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_xphoto.a
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libfreetype.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libz.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libexpat.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkDomainsChemistryOpenGL2-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libjpeg.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libpng.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libtiff.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersGeneric-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersHyperTree-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelDIY2-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelFlowPaths-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelGeometry-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelImaging-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelMPI-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelStatistics-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersPoints-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersProgrammable-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersPython-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libpython3.8.so
VideoDenoiser: /usr/lib/libvtkWrappingTools-7.1.a
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersReebGraph-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersSMP-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersSelection-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersTexture-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersVerdict-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkverdict-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkGUISupportQt-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkGUISupportQtSQL-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOAMR-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libsz.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libdl.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libm.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5_hl.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOEnSight-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libnetcdf_c++.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libnetcdf.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOExport-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingGL2PSOpenGL2-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgl2ps.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOFFMPEG-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOMovie-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libtheoraenc.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libtheoradec.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libogg.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOGDAL-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOGeoJSON-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libjsoncpp.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOImport-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOInfovis-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libxml2.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOMINC-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOMPIImage-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOMPIParallel-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOParallel-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIONetCDF-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOMySQL-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOODBC-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOPLY-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOParallelExodus-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOParallelLSDyna-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOParallelNetCDF-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOParallelXML-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOPostgreSQL-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOTecplotTable-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOVPIC-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkVPIC-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOVideo-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOXdmf2-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkxdmf2-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkImagingMorphological-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkImagingStatistics-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkImagingStencil-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkInfovisBoostGraphAlgorithms-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkInteractionImage-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkLocalExample-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkParallelMPI4Py-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingContextOpenGL2-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingExternal-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingFreeTypeFontConfig-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingImage-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingLOD-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingMatplotlib-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingParallel-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingParallelLIC-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingSceneGraph-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingVolumeAMR-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingVolumeOpenGL2-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkTestingGenericBridge-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkTestingIOSQL-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkTestingRendering-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkViewsContext2D-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkViewsGeovis-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkViewsInfovis-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkWrappingJava-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgflags.so.2.2.2
VideoDenoiser: libtvl1flow.a
VideoDenoiser: libclg7.a
VideoDenoiser: libsms_of.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/3rdparty/lib/libade.a
VideoDenoiser: /home/ivan/anaconda3/lib/libhdf5.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/librt.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libpthread.so
VideoDenoiser: /home/ivan/anaconda3/lib/libz.so
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_shape.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv.sfm.correspondence.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv.sfm.multiview.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv.sfm.numeric.a
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libglog.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgflags.so.2.2.2
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_highgui.a
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgtk-3.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgdk-3.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libpangocairo-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libpango-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libharfbuzz.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libatk-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libcairo-gobject.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libcairo.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgdk_pixbuf-2.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgio-2.0.so
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_datasets.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_plot.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_text.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_ml.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_phase_unwrapping.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_optflow.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_ximgproc.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_video.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_videoio.a
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstbase-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstreamer-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstapp-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstriff-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstpbutils-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstaudio-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstvideo-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstbase-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstreamer-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstapp-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstriff-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstpbutils-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstaudio-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgstvideo-1.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libgobject-2.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libglib-2.0.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libavcodec.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libavformat.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libavutil.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libswscale.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libavresample.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libc.so
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_imgcodecs.a
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libwebp.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libImath.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libIlmImf.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libIex.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libHalf.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libIlmThread.so
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_objdetect.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_calib3d.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_dnn.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/3rdparty/lib/liblibprotobuf.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_features2d.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_flann.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_photo.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_imgproc.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/lib/libopencv_core.a
VideoDenoiser: /usr/lib/x86_64-linux-gnu/liblapack.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libcblas.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libatlas.so
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/3rdparty/lib/libittnotify.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/3rdparty/lib/libippiw.a
VideoDenoiser: /home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/3rdparty/ippicv/ippicv_lnx/icv/lib/intel64/libippicv.a
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libjpeg.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libpng.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libtiff.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkDomainsChemistry-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersFlowPaths-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.12.8
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.12.8
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.12.8
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOGeometry-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOExodus-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkexoIIc-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libnetcdf_c++.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libnetcdf.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOLSDyna-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libsz.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libdl.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libm.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5_hl.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libxml2.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkWrappingPython38Core-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkPythonInterpreter-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libpython3.8.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersParallel-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkParallelMPI-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingLICOpenGL2-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersAMR-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkParallelCore-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOLegacy-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingOpenGL2-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libGLEW.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libSM.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libICE.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libX11.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libXext.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libXt.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkImagingMath-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOSQL-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkGeovisCore-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkproj4-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkChartsCore-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingContext2D-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersImaging-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkInfovisLayout-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkInfovisCore-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkViewsCore-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkInteractionWidgets-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersHybrid-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkImagingGeneral-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkImagingSources-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersModeling-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkInteractionStyle-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersExtraction-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersStatistics-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkImagingFourier-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkalglib-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkImagingHybrid-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOImage-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkDICOMParser-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkmetaio-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libz.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingAnnotation-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkImagingColor-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingVolume-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkImagingCore-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOXML-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOXMLParser-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkIOCore-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingLabel-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingFreeType-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkRenderingCore-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkCommonColor-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersGeometry-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersSources-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersGeneral-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkCommonComputationalGeometry-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkFiltersCore-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkCommonExecutionModel-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkCommonDataModel-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkCommonTransforms-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkCommonMisc-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkCommonMath-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkCommonSystem-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libfreetype.so
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtkCommonCore-7.1.so.7.1p.1
VideoDenoiser: /usr/lib/x86_64-linux-gnu/libvtksys-7.1.so.7.1p.1
VideoDenoiser: CMakeFiles/VideoDenoiser.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ivan/yasser/VideoDenoiser/source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable VideoDenoiser"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/VideoDenoiser.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/VideoDenoiser.dir/build: VideoDenoiser

.PHONY : CMakeFiles/VideoDenoiser.dir/build

CMakeFiles/VideoDenoiser.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/VideoDenoiser.dir/cmake_clean.cmake
.PHONY : CMakeFiles/VideoDenoiser.dir/clean

CMakeFiles/VideoDenoiser.dir/depend:
	cd /home/ivan/yasser/VideoDenoiser/source/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ivan/yasser/VideoDenoiser/source /home/ivan/yasser/VideoDenoiser/source /home/ivan/yasser/VideoDenoiser/source/build /home/ivan/yasser/VideoDenoiser/source/build /home/ivan/yasser/VideoDenoiser/source/build/CMakeFiles/VideoDenoiser.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/VideoDenoiser.dir/depend

