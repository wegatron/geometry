cmake_minimum_required(VERSION 2.6)

add_executable(test_mesh test_mesh.cpp)
target_link_libraries(test_mesh jtf-mesh)

add_executable(bnf test_bilateral_normal_filter.cpp)
target_link_libraries(bnf bilateral_filter jtf-mesh)