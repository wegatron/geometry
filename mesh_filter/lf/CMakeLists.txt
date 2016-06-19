cmake_minimum_required(VERSION 2.6)

file(GLOB src *.cpp)
file(GLOB headers *.h)

add_library(lf SHARED ${src} ${headers})

target_link_libraries(lf mesh_filter_basic zswlib ${OPENMESH_LIBRARIES})
