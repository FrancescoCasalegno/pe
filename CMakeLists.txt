cmake_minimum_required(VERSION 3.0)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)
project(BBP_PE VERSION 0.0.1 LANGUAGES CXX C)

find_package(CORE REQUIRED)
find_package(PETSc REQUIRED)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()



set(SOURCES
app.cc
assemble.cc
integrate.cc
linsys.cc
post.cc
pre.cc
utils.cc
)


set(HEADERS
app.h
integrate.h
linsys.h
utils.h
)


#set(SCOREC_INSTALL_DIR "/home/casalegn/Documents/core_install")
#include_directories(${SCOREC_INSTALL_DIR})
include_directories(${CORE_INCLUDE_DIR})
include_directories(${PETSC_INCLUDE_DIR})

add_library(pe ${SOURCES})

#target_include_directories(pe PUBLIC
#    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#    $<INSTALL_INTERFACE:include>
#)

add_executable(pe_exec main.cc)
target_link_libraries(pe_exec pe ${PETSC_LIBRARIES} ${CORE_LIBRARIES})

#bob_export_target(pe_exec)
#bob_end_subdir()

