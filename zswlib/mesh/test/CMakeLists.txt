cmake_minimum_required(VERSION 2.6)

project(zswlib)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};$ENV{HOME}/usr/share/cmake/Modules/")

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
#include($ENV{HOME}/usr/share/cmake/Modules/zsw.cmake)

#include_directories(CMAKE_CURRENT_SOURCE_DIR/include)

include_geo_sim_sdk()
link_geo_sim_sdk()

# include_beta()
# link_beta()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

###
##FIND_PACKAGE
###

# set(FLANN_ROOT E:/3rd_install/flann)
# find_package(Flann REQUIRED)
# include_directories(${FLANN_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)


#for openmesh
# include_directories(E:/usr/Windows/64/vc2010)
# set(OPENMESH_LIBRARIES E:/usr/Windows/64/vc2010/OpenMesh/Core/bin/OpenMeshCore.lib)
find_package(OpenMesh REQUIRED)
include_directories($OPENMESH_INCLUDE_DIRS$)

# add_subdirectory(mesh)
# add_subdirectory(log)

add_executable(test_rRingVertex ../mesh_op.cpp test_rRingVertex.cpp)
target_link_libraries(test_rRingVertex ${OPENMESH_LIBRARIES})
