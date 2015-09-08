cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH "$ENV{HOME}/usr/share/cmake/Modules/;${CMAKE_MODULE_PATH}")

project(matrix_conv)

set(APPNAME test_matrix_conv)

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

include_geo_sim_sdk()
link_geo_sim_sdk()

set(CMAKE_C_FLAGS "$ENV{CFLAGS} -O3")
set(CMAKE_CXX_FLAGS "-std=c++0x -fPIC -fopenmp ${CMAKE_CXX_FLAGS}")

###
##FIND_PACKAGE
###

find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
  message("CAN NOT FIND BOOST")
endif(Boost_FOUND)


aux_source_directory_with_headers(. sources headers)

add_executable(${APPNAME} ${sources} ${headers})

# target_link_libraries( ${APPNAME}

#   )

# add_library_from_dir(name mode dir)
# add_executable_from_dir(name dir)

## install header
# install(DIRECTORY include/
#   DESTINATION include/
#   FILES_MATCHING PATTERN "**.h" REGEX ".svn.*" EXCLUDE)

## install bin or library in /os/bit/compiler/
# geo_sim_sdk_install_cxx(TARGET name)

## install bin or library in /os/bit/
# geo_sim_sdk_install_c(TARGET name)
