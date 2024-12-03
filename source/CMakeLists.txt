cmake_minimum_required(VERSION 2.8)
project(VideoDenoiser)
set(CMAKE_CXX_STANDARD 11)  # Specify the C++ standard
cmake_policy(SET CMP0028 NEW)

# Set the path to your OpenCV build and installation directories
set(OpenCV_DIR "/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build")
set(OpenCV_INCLUDE_DIRS "/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/install/include")
set(OpenCV_LIB_DIR "/home/ivan/OpticalFlow/MeshFlow_Video_Denoising/opencv/build/install/lib")

# Find required packages
find_package(OpenCV REQUIRED PATHS ${OpenCV_DIR} NO_DEFAULT_PATH)
find_package(Eigen3 REQUIRED)
find_package(Iconv REQUIRED)
find_package(VTK REQUIRED)
find_package(gflags REQUIRED)

# FFTW3 Configuration
set(FFTW3_INCLUDE_DIRS "/usr/include")
set(FFTW3_LIBRARIES -lfftw3 -lfftw3f)

# Add sptwoIPOL library
add_library(sptwoIPOL STATIC
    ../lib/sptwoIPOL/src/library/libBasic.cpp
    ../lib/sptwoIPOL/src/library/libImage.cpp
    ../lib/sptwoIPOL/src/libNoiseVideo/libNLPCA.cpp
    ../lib/sptwoIPOL/src/libNoiseVideo/libFlow.cpp
    ../lib/sptwoIPOL/src/library/libImageFormatPM.cpp
    ../lib/sptwoIPOL/src/library/libImageFormats.c
    ../lib/sptwoIPOL/src/libNoiseVideo/libDenoisingVideo.cpp
)

target_include_directories(sptwoIPOL PUBLIC
    ../lib/sptwoIPOL/src
    ../lib/sptwoIPOL/src/library
    ../lib/sptwoIPOL/src/libNoiseVideo
    ${FFTW3_INCLUDE_DIRS}
)

# Add IIO library
add_library(iio STATIC
    ../lib/iio/iio.c
)

target_link_libraries(sptwoIPOL
    iio
    ${FFTW3_LIBRARIES}
    -ltiff
    -lpng
    -ljpeg
)

# Add sptwo processor library
add_library(sptwo_processor
    sptwo_processor.cpp
)

target_link_libraries(sptwo_processor
    ${OpenCV_LIBS}
    sptwoIPOL
)

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
    sptwo_processor.cpp
)

# Header files
set(HEADERS
    DenseOF.h
    Fast_klt.h
    Mesh.h
    MeshFlow.h
    MotionDenoiser.h
    VideoIO.h
    sptwo_processor.h
)

# Include directories
include_directories(
    ${OpenCV_INCLUDE_DIRS} 
    ${VTK_INCLUDE_DIRS} 
    ../lib/tvl1flow 
    ../lib/clg_7/source 
    ../lib/sms_optic_flow_2.0/spatial 
    ../lib/sms_optic_flow_2.0/temporal
    ../lib/sptwoIPOL/src
    ../lib/sptwoIPOL/src/library
    ../lib/sptwoIPOL/src/libNoiseVideo
)

# Link OpenCV libraries
link_directories(${OpenCV_LIB_DIR})

# Create executable and link libraries
add_executable(VideoDenoiser ${SOURCES} ${HEADERS})
target_link_libraries(VideoDenoiser 
    ${OpenCV_LIBS}
    Eigen3::Eigen 
    Iconv::Iconv 
    ${VTK_LIBRARIES} 
    gflags 
    tvl1flow 
    clg7 
    sms_of
    sptwoIPOL
    sptwo_processor
)
