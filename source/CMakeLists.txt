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
include_directories(${OpenCV_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})

# Link OpenCV libraries
link_directories(${OpenCV_LIB_DIR})
add_executable(VideoDenoiser ${SOURCES} ${HEADERS})
target_link_libraries(VideoDenoiser ${OpenCV_LIBS}   Eigen3::Eigen Iconv::Iconv ${VTK_LIBRARIES} gflags) 
