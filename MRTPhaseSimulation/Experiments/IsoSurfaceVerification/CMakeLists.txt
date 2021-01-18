add_executable(NewbySchwemmer NewbySchwemmer.cpp)
target_link_libraries(NewbySchwemmer std++fs
        ${MPI_LIBRARIES}
        libSPhaseFile
        libMPIFunctions
        libSDEToolbox)
