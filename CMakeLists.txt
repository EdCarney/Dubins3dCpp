cmake_minimum_required(VERSION 3.21)

project(Dubins3D VERSION 1.0)

# specify C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# configure a header file to pass in some config settings
configure_file(src/Dubins3d.h.in Dubins3d.h)

# add libraries
add_library(Dubins3d src/DubinsManeuver2d.cpp src/DubinsManeuver3d.cpp)

list(APPEND EXTRA_LIBS Dubins3d)

# add the executables
add_executable(test src/test.cpp)

target_link_libraries(test PUBLIC ${EXTRA_LIBS})

target_include_directories(test PUBLIC "${PROJECT_BINARY_DIR}")

# set include directories
set(Dubins3D_build_include_dirs
  "${Dubins3D_SOURCE_DIR}/include"
  "${Dubins3D_SOURCE_DIR}")
include_directories(${Dubins3D_build_include_dirs})