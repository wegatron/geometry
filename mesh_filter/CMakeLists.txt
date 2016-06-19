cmake_minimum_required(VERSION 2.6)

project(mesh_filter)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};$ENV{HOME}/usr/share/cmake/Modules/")

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
#include($ENV{HOME}/usr/share/cmake/Modules/zsw.cmake)

#include_directories(CMAKE_CURRENT_SOURCE_DIR/include)

include_geo_sim_sdk()
link_geo_sim_sdk()

#include_beta()
#link_beta()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

###
##FIND_PACKAGE
###

#set(FLANN_ROOT E:/3rd_install/flann)
#find_package(Flann REQUIRED)
#include_directories(${FLANN_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

#for openmesh
# include_directories(E:/usr/Windows/64/vc2010)
# set(OPENMESH_LIBRARIES E:/usr/Windows/64/vc2010/OpenMesh/Core/bin/OpenMeshCore.lib)
find_package(OpenMesh REQUIRED)
include_directories(${OPENMESH_INCLUDE_DIRS})


# add_subdirectory(src)

file(GLOB headers *.h)
file(GLOB src *.cpp)

#add_library(mesh_filter_basic SHARED ${src} ${headers})

#target_link_libraries(mesh_filter_basic zswlib ${OPENMESH_LIBRARIES})

#add_subdirectory(lhf)
#add_subdirectory(lf)
#add_subdirectory(bnf)
add_subdirectory(app)
#add_subdirectory(test)
