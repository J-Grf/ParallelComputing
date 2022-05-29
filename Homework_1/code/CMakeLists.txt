cmake_minimum_required(VERSION 2.8)
PROJECT(read-mesh)

set(VTK_DIR "/usr/lib64/vtk")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
add_compile_options(-g -std=c++17)

include_directories(.)
file(GLOB SOURCES "*.cpp")

add_executable(read-mesh ${SOURCES})

target_link_libraries(read-mesh ${VTK_LIBRARIES})

