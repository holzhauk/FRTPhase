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
        NormMoments.cpp NormMoments.h ModelFactory.cpp ModelFactory.h SimConfig_t.cpp SimConfig_t.h Adapters.cpp Adapters.h)

find_package(MPI REQUIRED)
target_link_libraries(newbyschwemmer ${MPI_LIBRARIES})

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set(_hdf5_libs hdf5 hdf5_cpp)
target_link_libraries(newbyschwemmer ${HDF5_LIBRARIES} ${_hdf5_libs})

find_package(Boost REQUIRED COMPONENTS filesystem)
target_link_libraries(newbyschwemmer Boost::filesystem)
