cmake_minimum_required(VERSION 2.8.3)
project(mono_tracker)

set(CMAKE_CXX_COMPILER "/usr/bin/g++-5")
set(CMAKE_C_COMPILER "/usr/bin/gcc-5")

# set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set macro variables from *.cmake file
find_package(catkin REQUIRED COMPONENTS
    roscpp
    std_msgs
    dvs_msgs
    sensor_msgs
    cv_bridge
    message_filters
    )
find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


# catkin-provided CMake macro
catkin_package()

# add include directories
include_directories(
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
    )
    
# add executable files
add_executable(mono_tracker
    src/mtracker_node.cpp
    src/mtracker.cpp
    ../msckf_estimator/src/utils/Attitude.cpp
    ../msckf_estimator/src/utils/utils.cpp
    )

# Link libraries
target_link_libraries(mono_tracker
    ${catkin_LIBRARIES}
    ${OpenCV_LIBS}
    ${OpenMP_CXX_FLAGS} # Link OpenMP library
)

