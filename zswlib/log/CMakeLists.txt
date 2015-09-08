cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH "$ENV{HOME}/usr/share/cmake/Modules/;${CMAKE_MODULE_PATH}")

project(zsw_log)

set(APPNAME zsw_log)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

include_geo_sim_sdk()
link_geo_sim_sdk()

set(CMAKE_C_FLAGS "$ENV{CFLAGS} -O3")
set(CMAKE_CXX_FLAGS "-std=c++11 -fPIC -fopenmp ${CMAKE_CXX_FLAGS}")

###
##FIND_PACKAGE
###

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
  message("CAN NOT FIND BOOST")
endif(Boost_FOUND)

# add_subdirectory(src)

# file(GLOB source *.cpp)

add_library(zswlog zsw_log.h zsw_log.cpp)

add_executable(test test.cpp)

target_link_libraries( test
  zswlog
  )


install(DIRECTORY include/
  DESTINATION include/
  FILES_MATCHING PATTERN "**.h" REGEX ".svn.*" EXCLUDE)

# add_library_from_dir(name mode dir)
# add_executable_from_dir(name dir)

## install header
# install(DIRECTORY include/
#   DESTINATION include/
#   FILES_MATCHING PATTERN "**.h" REGEX ".svn.*" EXCLUDE)

#install bin or library in /os/bit/compiler/
geo_sim_sdk_install_cxx(TARGETS zswlog)

## install bin or library in /os/bit/
# geo_sim_sdk_install_c(TARGET name)
