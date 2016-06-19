cmake_minimum_required(VERSION 2.6)

set( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};$ENV{HOME}/usr/share/cmake/Modules/")

project(isotopic_approximation)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

include_geo_sim_sdk()
link_geo_sim_sdk()
link_directories($ENV{HOME}/usr/lib)
#include_directories(CMAKE_CURRENT_SOURCE_DIR/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -fPIC")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS} -std=c++11 -fopenmp -fPIC -O3 -ggdb")

###
##FIND_PACKAGE
###

find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})


add_executable(testt basic_data_structure.cpp main.cpp)

target_link_libraries(testt
  ${CGAL_LIBRARY}
  ${GMP_LIBRARIES})

# add_subdirectory()

# aux_source_directory_with_headers(. src headers)

# add_library_from_dir(name mode dir)
# add_executable_from_dir(name dir)

# default_install_cxx(TARGET name)
# default_install_c(TARGET name)
