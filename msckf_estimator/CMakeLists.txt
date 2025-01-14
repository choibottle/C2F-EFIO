cmake_minimum_required(VERSION 2.8.3)
project(smsckf_nesl)

set(CMAKE_CXX_COMPILER "/usr/bin/g++-5")
set(CMAKE_C_COMPILER "/usr/bin/gcc-5")

# set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  sensor_msgs
  genmsg
  nav_msgs
  tf
)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenCV REQUIRED)

catkin_package()

include_directories(
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
)

add_executable(estimator_node
    src/estimator_node.cpp
    src/core/msckf.cpp
    src/core/mapmanager.cpp
    src/utils/Attitude.cpp
    src/utils/utils.cpp
    src/utils/NavEarth.cpp)

target_link_libraries(estimator_node
    ${catkin_LIBRARIES}
    ${OpenCV_LIBRARIES})
