cmake_minimum_required(VERSION 3.21)

project(Dubins3D VERSION 1.0)

# specify C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# configure a header file to pass in some config settings
configure_file(Dubins3d.h.in Dubins3d.h)

# add libraries
add_library(Dubins3d DubinsManeuver2d.cpp DubinsManeuver3d.cpp)

list(APPEND EXTRA_LIBS Dubins3d)

# add the executables
add_executable(test test.cpp)

target_link_libraries(test PUBLIC ${EXTRA_LIBS})

target_include_directories(test PUBLIC "${PROJECT_BINARY_DIR}")