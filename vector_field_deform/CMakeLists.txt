cmake_minimum_required(VERSION 2.8)

enable_language(Fortran)

set( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};$ENV{HOME}/usr/share/cmake/Modules/")

set(CMAKE_CXX_FLAGS "-std=c++11 -fPIC ${CMAKE_CXX_FLAGS} -fopenmp")

find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

set(SRC_DIR "{${CMAKE_CURRENT_SOURCE_DIR},$ENV{HOME}/usr/share/maxima}")

add_custom_command(OUTPUT quad_scalar_field.f90
  COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b quad_scalar_field.mac
  DEPENDS quad_scalar_field.mac)

#add_library(quad_scalar_field SHARED quad_scalar_field.f90)

file(GLOB src *.cpp)
file(GLOB headers *.h)

add_library(vector_fields SHARED ${src} ${headers} quad_scalar_field.f90)
add_subdirectory(test)