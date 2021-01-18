//
// Created by konstantin on 5/16/20.
//

#include "Solver.h"

bool Solver::in_time() {
    return (x[2] >= config->x0[2]) && (x[2] < config->T);
}

std::array<double, 3> Solver::run() {
    x[0] = config->x0[0];
    x[1] = config->x0[1];
    x[2] = config->x0[2];
    while(this->in_time()) {
        x = this->evolve();
    }
    return x;
}

MPI_Datatype IsoPlanarSOsc::config_t::mpiType() {
    MPI_Datatype mpi_config_t;
    int blengths[3] = {1, 1, 3};
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[3];
    offsets[0] = offsetof(IsoPlanarSOsc::config_t, dt);
    offsets[1] = offsetof(IsoPlanarSOsc::config_t, T);
    offsets[2] = offsetof(IsoPlanarSOsc::config_t, x0);
    MPI_Type_create_struct(3, blengths, offsets, types, &mpi_config_t);
    MPI_Type_commit(&mpi_config_t);
    return mpi_config_t;
}