add_executable(IsoSurfaceVerification IsoSurfaceVerification.cpp)
target_link_libraries(IsoSurfaceVerification stdc++fs
        ${MPI_LIBRARIES}
        libSPhaseFile
        libMPIFunctions
        libSDEToolbox)
