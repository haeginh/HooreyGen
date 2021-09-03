cmake_minimum_required(VERSION 3.9.0)

project(fileGen)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 99)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(LIBIGL REQUIRED QUIET)
set(LIBIGL_DIR ${LIBIGL_INCLUDE_DIR}/..)
include_directories(${LIBIGL_DIR})

add_executable(fileGen fileGen.cc)
target_link_libraries(fileGen Eigen3::Eigen igl::core)
