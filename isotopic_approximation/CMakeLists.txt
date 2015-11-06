cmake_minimum_required(VERSION 2.6)

set( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};$ENV{HOME}/usr/share/cmake/Modules/")

project(isotopic_approximation)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

include_geo_sim_sdk()
link_geo_sim_sdk()
link_directories($ENV{HOME}/usr/lib)
#include_directories(CMAKE_CURRENT_SOURCE_DIR/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -fPIC")

###
##FIND_PACKAGE
###

find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

find_package(OpenMesh REQUIRED)
include_directories(${OPENMESH_INCLUDE_DIRS})

find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories(${BOOST_INCLUDE_DIRES})

# add_subdirectory(src)

file(GLOB headers *.h)
file(GLOB src *.cpp)

add_library(isotropic_approximate SHARED ${headers} ${src})

target_link_libraries(isotropic_approximate zswlib ${OPENMESH_LIBRARIES}
  ipopt coinhsl coinmumps coinmetis blas lapack dl mpi)

add_subdirectory(test)
add_subdirectory(app)

# add_subdirectory()

# aux_source_directory_with_headers(. src headers)

# add_library_from_dir(name mode dir)
# add_executable_from_dir(name dir)

# default_install_cxx(TARGET name)
# default_install_c(TARGET name)
