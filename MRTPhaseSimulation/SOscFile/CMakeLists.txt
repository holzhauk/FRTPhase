file(GLOB_RECURSE source_list "*.cpp" "*.h" "*.hpp")
add_library( SOscFile SHARED ${source_list} )
target_include_directories( SOscFile PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
                                    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src )
target_link_libraries( SOscFile ${HDF5_C_LIBRARIES}
                                ${HDF5_CXX_LIBRARIES}
                                stdc++fs )