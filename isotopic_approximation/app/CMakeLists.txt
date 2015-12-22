cmake_minimum_required(VERSION 2.6)

project(approximate)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

#include_directories(CMAKE_CURRENT_SOURCE_DIR/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -fPIC")

###
##FIND_PACKAGE
###

# add_subdirectory(src)

add_executable(approximate main.cpp)
add_executable(gen_noise noise_maker.cpp)

target_link_libraries(approximate
  isotropic_approximate
  zswlib
  ${OPENMESH_LIBRARIES}
  ${CGAL_LIBRARY}
  ${GMP_LIBRARIES}
  )

target_link_libraries(gen_noise zswlib ${OPENMESH_LIBRARIES})
# add_library(lib_name SHARED headers src)
# add_subdirectory()

# aux_source_directory_with_headers(. src headers)

# add_library_from_dir(name mode dir)
# add_executable_from_dir(name dir)

# default_install_cxx(TARGET name)
# default_install_c(TARGET name)
