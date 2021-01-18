#file(GLOB_RECURSE source_list ".cpp" ".h" ".hpp")
add_executable(test_SOscFile test_SOscFile.cpp)
target_link_libraries(test_SOscFile
                        Boost::unit_test_framework)