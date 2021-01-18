file(GLOB_RECURSE source_list "*.cpp" "*.h" "*.hpp")
add_library( SOscModels SHARED ${source_list} include/SOscModels/SOscModels.h src/SOscModels/Solver.cpp src/SOscModels/Solver.h)
target_include_directories( SOscModels PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src )