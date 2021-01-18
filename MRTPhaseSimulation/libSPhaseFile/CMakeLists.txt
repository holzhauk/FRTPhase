file(GLOB_RECURSE source_list ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)
file(GLOB_RECURSE include_list ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
add_library(libSPhaseFile SHARED ${source_list} ${include_list})
target_include_directories( libSPhaseFile PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
set_target_properties(libSPhaseFile PROPERTIES PUBLIC_HEADER include/libSPhaseFile/libSPhaseFile.h)
target_link_libraries(libSPhaseFile ${HDF5_C_LIBRARIES}
        ${HDF5_CXX_LIBRARIES}
        Boost::filesystem
        stdc++fs)

add_executable(libSPhaseFile.test ${CMAKE_CURRENT_SOURCE_DIR}/tests/libSPhaseFile.test.cpp)
#target_include_directories(libSPhaseFile.test ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(libSPhaseFile.test
        Boost::unit_test_framework
        ${HDF5_C_LIBRARIES}
        ${HDF5_CXX_LIBRARIES}
        stdc++fs
        libSPhaseFile)