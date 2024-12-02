cmake_minimum_required(VERSION 2.8)
project(VideoDenoiser)
set(CMAKE_CXX_STANDARD 11)  # Specify the C++ standard
cmake_policy(SET CMP0028 NEW)

# Set the path to your OpenCV build and installation directories
set(OpenCV_DIR "/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build")  # Adjust this path
set(OpenCV_INCLUDE_DIRS "/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/install/include")  # Adjust this path
set(OpenCV_LIB_DIR "/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/install/lib")  # Adjust this path

# Find OpenCV package
find_package(OpenCV REQUIRED PATHS ${OpenCV_DIR} NO_DEFAULT_PATH)
# Find Eigen3 package
find_package(Eigen3 REQUIRED)
find_package(Iconv REQUIRED)
find_package(VTK REQUIRED)

option(WITH_VTK "" OFF)
set(WITH_VTK OFF CACHE BOOL "" FORCE)
# Add gflags package
find_package(gflags REQUIRED)

# Add TVL1 library
add_library(tvl1flow STATIC
    ../lib/tvl1flow/tvl1flow_lib.c
    ../lib/tvl1flow/mask.c
    ../lib/tvl1flow/bicubic_interpolation.c
    ../lib/tvl1flow/zoom.c
)

# Add CLG library
add_library(clg7 STATIC
    ../lib/clg_7/source/bicubic_interpolation.c
    ../lib/clg_7/source/clg_of.c
    ../lib/clg_7/source/zoom.c
    ../lib/clg_7/source/clg_of_main.c
    ../lib/clg_7/source/iio.c
    ../lib/clg_7/source/mask.c
    ../lib/clg_7/source/util.c
     
    )

# Add SMS optical flow library
add_library(sms_of STATIC
    ../lib/sms_optic_flow_2.0/spatial/brox_optic_flow.h
    ../lib/sms_optic_flow_2.0/spatial/mask.h
    ../lib/sms_optic_flow_2.0/spatial/bicubic_interpolation.h
    ../lib/sms_optic_flow_2.0/spatial/zoom.h
    ../lib/sms_optic_flow_2.0/spatial/gaussian.h
    ../lib/sms_optic_flow_2.0/temporal/brox_optic_flow.h
    ../lib/sms_optic_flow_2.0/temporal/mask.h
    ../lib/sms_optic_flow_2.0/temporal/bicubic_interpolation.h
    ../lib/sms_optic_flow_2.0/temporal/zoom.h
    ../lib/sms_optic_flow_2.0/temporal/gaussian.h
)

# Set SMS library properties
set_target_properties(sms_of PROPERTIES 
    LINKER_LANGUAGE CXX
    CXX_STANDARD 11
)

# Add compile definitions for SMS
target_compile_definitions(sms_of PRIVATE
    SMS_SPATIAL_IMPLEMENTATION
    SMS_TEMPORAL_IMPLEMENTATION
)

# Source files
set(SOURCES
    main.cpp
    DenseOF.cpp
    Fast_klt.cpp
    Mesh.cpp
    MeshFlow.cpp
    MotionDenoiser.cpp
    VideoIO.cpp
)

# Header files
set(HEADERS
    DenseOF.h
    Fast_klt.h
    Mesh.h
    MeshFlow.h
    MotionDenoiser.h
    VideoIO.h
)

# Include OpenCV directories
include_directories(${OpenCV_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS} ../lib/tvl1flow ../lib/clg_7/source ../lib/sms_optic_flow_2.0/spatial ../lib/sms_optic_flow_2.0/temporal)

# Link OpenCV libraries
link_directories(${OpenCV_LIB_DIR})
add_executable(VideoDenoiser ${SOURCES} ${HEADERS})
target_link_libraries(VideoDenoiser ${OpenCV_LIBS}   Eigen3::Eigen Iconv::Iconv ${VTK_LIBRARIES} gflags tvl1flow clg7 sms_of) 
