cmake_minimum_required(VERSION 2.6)

set( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};$ENV{HOME}/usr/share/cmake/Modules/")



include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

include_geo_sim_sdk()
link_geo_sim_sdk()

#include_directories(CMAKE_CURRENT_SOURCE_DIR/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -fPIC")

project(test)

###
##FIND_PACKAGE
###

find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

find_package(OpenMesh REQUIRED)
include_directories(${OPENMESH_INCLUDE_DIRS})


add_executable(test_sampling test_sampling.cpp)
add_executable(test_gen_surf test_surf_genator.cpp)

target_link_libraries(test_sampling isotropic_approximat ${OPENMESH_LIBRARIES})
target_link_libraries(test_gen_surf isotropic_approximat ${OPENMESH_LIBRARIES})