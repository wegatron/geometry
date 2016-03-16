cmake_minimum_required(VERSION 2.6)

project(approximate)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

#include_directories(CMAKE_CURRENT_SOURCE_DIR/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -fPIC")

###
##FIND_PACKAGE
###

# add_subdirectory(src)

# add_executable(test_tri test_triangulation.cpp)
# add_executable(debug_normal_cond debug_normal_condition.cpp)

# target_link_libraries(test_tri
#   iso_app
#   zswlib
#   ${OPENMESH_LIBRARIES}
#   ${CGAL_LIBRARY}
#   ${GMP_LIBRARIES}
#   )

add_executable(test_mls_basic test_mls_basic.cpp)
target_link_libraries(test_mls_basic
  iso_app
  zswlib
  ${OPENMESH_LIBRARIES}
  ${CGAL_LIBRARY}
  ${GMP_LIBRARIES}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

# target_link_libraries(debug_normal_cond
#   iso_app
#   zswlib
#   ${OPENMESH_LIBRARIES}
#   ${CGAL_LIBRARY}
#   ${GMP_LIBRARIES}
#   )
