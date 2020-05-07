cmake_minimum_required(VERSION 2.8)
 
PROJECT(prog03)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(prog03 prog03.cpp vtkEasyTransfer.hpp)
 
target_link_libraries(prog03 ${VTK_LIBRARIES})