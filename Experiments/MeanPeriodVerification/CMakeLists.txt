add_executable(MeanPeriodVerification MeanPeriodVerification.cpp)
target_link_libraries(MeanPeriodVerification stdc++fs
        ${MPI_LIBRARIES}
        libSPhaseFile
        libMPIFunctions
        libSDEToolbox)