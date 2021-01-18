//
// Created by konstantin on 1/16/21.
//
#include "MPISimConfig.h"

void MPI_Share(int world_rank, Config::Simulation& simConfig) {
    SimConfig_t config(simConfig);
    MPI_Datatype MPI_SimConfig_t = config.get_mpi_type();
    MPI_Bcast(&config, 1, MPI_SimConfig_t, 0, MPI_COMM_WORLD);
    MPI_Type_free(&MPI_SimConfig_t);
    simConfig = ConfigAdapter(config);
}

unsigned int SubEnsembleSize(int world_rank, int world_size, const int& ensemble_size){
    // calculate the individual ensemble sizes
    div_t int_div_res = div(static_cast<int>(ensemble_size), world_size); // int_dev_res = eSize / worldSize
    vector<int> SubEnsembleSizes = vector<int>(world_size, int_div_res.quot);
    for (int i = 0; i < int_div_res.rem; i++){
        SubEnsembleSizes[(world_size-1)-i] += 1;
    }
    return static_cast<unsigned int>(SubEnsembleSizes[world_rank]);
}



SimConfig_t::SimConfig_t(Config::Simulation& simConfig) {
    dt = simConfig.dt;
    t0 = simConfig.t0;
    T = simConfig.T;
    ensembleSize = simConfig.EnsembleSize;
    sampleSize = simConfig.SampleSize;
}

MPI_Datatype SimConfig_t::get_mpi_type() {
    // registration of mpi datatype for communication
    // of extension of data to be shared across the nodes
    MPI_Datatype MPI_SimConfig_t;
    int blengths[5] = {1, 1, 1, 1, 1};
    MPI_Datatype types[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UNSIGNED, MPI_UNSIGNED};
    MPI_Aint offsets[5];
    offsets[0] = offsetof(SimConfig_t, dt);
    offsets[1] = offsetof(SimConfig_t, t0);
    offsets[2] = offsetof(SimConfig_t, T);
    offsets[3] = offsetof(SimConfig_t, ensembleSize);
    offsets[4] = offsetof(SimConfig_t, sampleSize);
    MPI_Type_create_struct(5, blengths, offsets, types, &MPI_SimConfig_t);
    MPI_Type_commit(&MPI_SimConfig_t);
    // return the composit type registered with mpi
    return MPI_SimConfig_t;
}

ConfigAdapter::ConfigAdapter(SimConfig_t& simConfig) {
    this->t0 = simConfig.t0;
    this->dt = simConfig.dt;
    this->T = simConfig.T;
    this->SampleSize = simConfig.sampleSize;
    this->EnsembleSize = simConfig.ensembleSize;
}