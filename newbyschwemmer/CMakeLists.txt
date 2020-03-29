cmake_minimum_required(VERSION 3.15)
project(newbyschwemmer)

set(CMAKE_CXX_STANDARD 20)

add_executable(newbyschwemmer main.cpp)

find_package(GSL REQUIRED)
target_link_libraries(newbyschwemmer GSL::gsl
                                     GSL::gslcblas)

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
target_link_libraries(newbyschwemmer ${HDF5_LIBRARIES})