//
// Created by konstantin on 1/13/21.
//

#ifndef MRTPHASESIMULATION_MPIISOSURFACEFILE_H
#define MRTPHASESIMULATION_MPIISOSURFACEFILE_H

#include <mpi.h>
#include <libSPhaseFile/libSPhaseFile.h>

void MPI_Share(int world_rank, IsoSurfaceFile &isoSurfaceFile);

#endif //MRTPHASESIMULATION_MPIISOSURFACEFILE_H
