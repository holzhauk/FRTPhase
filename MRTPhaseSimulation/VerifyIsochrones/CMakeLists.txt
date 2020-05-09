file( GLOB_RECURSE source_list "*.cpp" "*.h" "*.hpp" )
add_executable( VerifyIsochrones ${source_list} )
target_include_directories( VerifyIsochrones PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/VerifyIsochrones )
target_link_libraries( VerifyIsochrones PRIVATE SOscModels PRIVATE SOscFile
        ${MPI_LIBRARIES}
        Boost::filesystem Boost::log_setup Boost::log Boost::system Boost::thread
        stdc++fs )
