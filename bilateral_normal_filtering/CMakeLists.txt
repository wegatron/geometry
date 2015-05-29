cmake_minimum_required(VERSION 2.6)

project(bilter_filter)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

#include_directories(CMAKE_CURRENT_SOURCE_DIR/include)

include_geo_sim_sdk()
link_geo_sim_sdk()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

###
##FIND_PACKAGE
###

# add_subdirectory(src)

file(GLOB headers *.h)
file(GLOB src *.cpp)

add_library(bilateral_filter SHARED ${src} ${headers})
add_subdirectory(test)

# aux_source_directory_with_headers(. sources headers)

# add_library_from_dir(name mode dir)
# add_executable_from_dir(name dir)

# default_install_cxx(TARGET name)
# default_install_c(TARGET name)
