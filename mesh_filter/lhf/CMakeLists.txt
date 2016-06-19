cmake_minimum_required(VERSION 2.6)

file(GLOB src *.cpp)
file(GLOB headers *.h)

add_library(lhf SHARED ${src} ${headers})

target_link_libraries(lhf mesh_filter_basic zswlib ${OPENMESH_LIBRARIES})
