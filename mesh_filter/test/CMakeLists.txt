cmake_minimum_required(VERSION 2.6)

# add_executable(bnf test_bilateral_normal_filter.cpp)
# target_link_libraries(bnf bilateral_filter)

#add_executable(diffv test_analysis.cpp)
#target_link_libraries(diffv bilateral_filter jtf-mesh2)

#add_executable(analysis analysis.cpp)
#target_link_libraries(analysis bilateral_filter jtf-mesh2)

#add_executable(sort_mesh_v sort_mesh_vertex.cpp)
#target_link_libraries(sort_mesh_v bilateral_filter jtf-mesh2)

#add_executable(test_flann test_flann.cpp)
#target_link_libraries(test_flann jtf-mesh2 ${FLANN_LIBRARIES})

# add_executable(test_arguement test_boost_arguement.cpp)
# target_link_libraries(test_arguement libboost_program_options-vc100-mt-1_58.lib)

# add_executable(test_processNVRing test_processNVRing.cpp)
# target_link_libraries(test_processNVRing bilateral_filter)

add_executable(test_rRingFacesByEdge test_rRingFacesByEdge.cpp)
target_link_libraries(test_rRingFacesByEdge bnf mesh_filter_basic zswlib ${OPENMESH_LIBRARIES})
