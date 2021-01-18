cmake_minimum_required(VERSION 3.13)
project(MRTPhaseSimulation)

set(CMAKE_CXX_STANDARD 20)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(HDF5 REQUIRED COMPONENTS C CXX)
include_directories(${HDF5_INCLUDE_DIRS})

ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)
find_package(Boost COMPONENTS system filesystem thread log log_setup unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(VerifyIsochrones)
add_subdirectory(SOscFile)
add_subdirectory(SOscModels)
add_subdirectory(test)
add_subdirectory(libSPhaseFile)
add_subdirectory(libSDEToolbox)
add_subdirectory(libMPIFunctions)
add_subdirectory(Experiments)