cmake_minimum_required(VERSION 2.6)


file(GLOB src *.cpp)
file(GLOB headers *.h)
add_library(bnf SHARED ${src} ${headers})

target_link_libraries(bnf zswlib ${OPENMESH_LIBRARIES})
