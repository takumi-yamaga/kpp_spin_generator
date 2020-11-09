#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_FLAGS "-O2 -std=c++11 -Wall")
project(generator)
#
#----------------------------------------------------------------------------
# Setup include directory for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories($ENV{ROOTSYS}/include)
#
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# include Root
#
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/cmake)
find_package(ROOT)
include_directories(${ROOT_INCLUDE_DIR})
#
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(generator kpp_spin_generator.cc ${sources} ${headers})
target_link_libraries(generator ${ROOT_LIBRARIES})
set(CMAKE_MACOSX_RPATH 1)

