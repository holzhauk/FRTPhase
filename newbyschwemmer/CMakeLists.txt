cmake_minimum_required(VERSION 3.15)
project(newbyschwemmer)

set(CMAKE_CXX_STANDARD 20)

add_executable(newbyschwemmer main.cpp
        Isochrone.h Isochrone.cpp
        IsochroneSet.h IsochroneSet.cpp
        IsoPlanarSOsc.cpp IsoPlanarSOsc.h
        NewSOsc.cpp NewSOsc.h
        Domain.cpp Domain.h
        MFPTs.cpp MFPTs.h
        MFPTSet.cpp MFPTSet.h
        Trajectory.cpp Trajectory.h
        TrajectorySet.cpp TrajectorySet.h
        NormMoments.cpp NormMoments.h
        ModelFactory.cpp ModelFactory.h
        SimConfig_t.cpp SimConfig_t.h
        Adapters.cpp Adapters.h)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(newbyschwemmer ${MPI_LIBRARIES})

find_package(HDF5 REQUIRED COMPONENTS C CXX)
include_directories(${HDF5_INCLUDE_DIRS})
target_link_libraries(newbyschwemmer ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})

ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)
find_package(Boost COMPONENTS system filesystem thread log log_setup REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(newbyschwemmer Boost::filesystem
                                     Boost::log_setup
                                     Boost::log
                                     Boost::system
                                     Boost::thread)

target_link_libraries(newbyschwemmer stdc++fs)
