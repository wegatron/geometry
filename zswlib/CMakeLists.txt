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
#find_package(OpenMesh REQUIRED)
#include_directories($OPENMESH_INCLUDE_DIRS$)

# add_subdirectory(mesh)

file(GLOB common_headers *.h)
file(GLOB mesh_headers mesh/*.h io/*.h)
file(GLOB mesh_src mesh/*.cpp)

add_library(zswlib SHARED ${mesh_src} log/zsw_log.cpp)
#target_link_libraries(zswlib ${OPENMESH_LIBRARIES})

# aux_source_directory_with_headers(. sources headers)

install(FILES ${common_headers} log/zsw_log.h
  DESTINATION include/zswlib)

install(FILES ${mesh_headers}
  DESTINATION include/zswlib/mesh)

geo_sim_sdk_install_cxx(TARGETS zswlib)

#add_subdirectory(test)
