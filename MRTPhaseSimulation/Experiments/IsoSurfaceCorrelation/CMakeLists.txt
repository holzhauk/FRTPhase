add_executable(SerialCorrelation SerialCorrelation.cpp)
target_link_libraries(SerialCorrelation stdc++fs
        ${MPI_LIBRARIES}
        ${HDF5_C_LIBRARIES}
        ${HDF5_CXX_LIBRARIES}
        libSPhaseFile
        libMPIFunctions
        libSDEToolbox)