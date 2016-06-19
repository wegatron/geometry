cmake_minimum_required(VERSION 2.6)

# add_executable(filter main.cpp)
# target_link_libraries(filter bnf)

add_executable(flip flip_edge.cpp ../bnf/EdgeFlipMean.cpp)

target_link_libraries(flip ${OPENMESH_LIBRARIES})
