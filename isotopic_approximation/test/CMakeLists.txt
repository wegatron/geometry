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

find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})

# add_executable(test_sampling2 test_sampling2.cpp)
# add_executable(test_gen_surf test_surf_genator.cpp)
# add_executable(test_cgal_triangulation test_triangle.cpp)
# add_executable(test_kernelRegionJudge test_KernelRegionJudge.cpp)
add_executable(test_triangulation test_triangulation.cpp)
add_executable(test_linkcondition test_linkcondition.cpp)
add_executable(test_simptolerance test_simptolerance.cpp)
add_executable(debug_simp_jpts debug_simp_jpts.cpp)
add_executable(test_tetquality test_tetquality.cpp)
add_executable(test_point2tridis test_point2tri_sqdis.cpp)
add_executable(test_mutual_tessellation test_mutual_tessellation.cpp)
add_executable(test_collapse_bc test_collapse_bc.cpp)

# target_link_libraries(test_sampling2 isotropic_approximate ${OPENMESH_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES})
# target_link_libraries(test_gen_surf isotropic_approximate ${OPENMESH_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES})
# target_link_libraries(test_cgal_triangulation ${CGAL_LIBRARY} ${GMP_LIBRARIES})
# target_link_libraries(test_kernelRegionJudge isotropic_approximate ${CGAL_LIBRARY} ${GMP_LIBRARIES}
#   ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
target_link_libraries(test_triangulation isotropic_approximate ${OPENMESH_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES} zswlib)
target_link_libraries(test_linkcondition isotropic_approximate ${OPENMESH_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES}
${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

target_link_libraries(test_simptolerance isotropic_approximate ${OPENMESH_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES}
${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

target_link_libraries(debug_simp_jpts isotropic_approximate ${OPENMESH_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES}
${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

target_link_libraries(test_tetquality isotropic_approximate ${OPENMESH_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES}
${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

target_link_libraries(test_point2tridis isotropic_approximate ${OPENMESH_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES}
${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

target_link_libraries(test_mutual_tessellation isotropic_approximate ${OPENMESH_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

target_link_libraries(test_collapse_bc isotropic_approximate ${OPENMESH_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES})
### old unused
# add_executable(test_optimizer test_optimizer.cpp)
# target_link_libraries(test_optimizer ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
#   isotropic_approximate ${CGAL_LIBRARY} ${GMP_LIBRARIES}
#   ipopt coinhsl coinmumps coinmetis blas lapack dl mpi)
# add_executable(test_collapseEdge test_collapseEdge.cpp)
# target_link_libraries(test_collapseEdge isotropic_approximate ${OPENMESH_LIBRARIES} ${CGAL_LIBRARY} ${GMP_LIBRARIES}
#   ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
