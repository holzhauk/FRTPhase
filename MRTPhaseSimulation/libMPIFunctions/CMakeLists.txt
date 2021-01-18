file(GLOB_RECURSE source_list ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)
file(GLOB_RECURSE include_list ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
add_library(libMPIFunctions SHARED ${source_list} ${include_list})
target_include_directories(libMPIFunctions PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
set_target_properties(libMPIFunctions PROPERTIES PUBLIC_HEADER include/libMPIFunctions/libMPIFunctions.h)
target_link_libraries(libMPIFunctions
        libSPhaseFile
        stdc++fs
        ${MPI_LIBRARIES})

add_executable(libMPIFunctions.test ${CMAKE_CURRENT_SOURCE_DIR}/tests/libMPIFunctions.test.cpp)
target_link_libraries(libMPIFunctions.test
        Boost::unit_test_framework
        stdc++fs
        libSPhaseFile
        libMPIFunctions
        ${MPI_LIBRARIES})