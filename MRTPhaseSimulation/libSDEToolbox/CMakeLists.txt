file(GLOB_RECURSE source_list ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)
file(GLOB_RECURSE include_list ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
add_library(libSDEToolbox SHARED ${source_list} ${include_list})
target_include_directories(libSDEToolbox PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
set_target_properties(libSDEToolbox PROPERTIES PUBLIC_HEADER include/libSDEToolbox/libSDEToolbox.h)
target_link_libraries(libSDEToolbox
        libSPhaseFile
        libMPIFunctions
        ${MPI_LIBRARIES}
        stdc++fs)

add_executable(libSDEToolbox.test ${CMAKE_CURRENT_SOURCE_DIR}/tests/libSDEToolbox.test.cpp)
target_link_libraries(libSDEToolbox.test
        Boost::unit_test_framework
        stdc++fs
        libSPhaseFile
        libMPIFunctions
        ${MPI_LIBRARIES}
        libSDEToolbox)