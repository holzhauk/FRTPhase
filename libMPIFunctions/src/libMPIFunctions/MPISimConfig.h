//
// Created by konstantin on 1/16/21.
//

#ifndef MRTPHASESIMULATION_MPISIMCONFIG_H
#define MRTPHASESIMULATION_MPISIMCONFIG_H

#include <cstdlib>
#include <string>
#include <mpi.h>
#include <libSPhaseFile/libSPhaseFile.h>

using namespace std;

void MPI_Share(int world_rank, Config::Simulation& simConfig);
void MPI_Share(int world_rank, string& modelName);
unsigned int SubEnsembleSize(int world_rank, int world_size, const int& ensemble_size);

struct SimConfig_t {
    double dt;
    double t0;
    double T;
    unsigned int ensembleSize;
    unsigned int sampleSize;
    SimConfig_t(Config::Simulation& simConfig);
    MPI_Datatype get_mpi_type();
};

struct ConfigAdapter: public Config::Simulation {
public:
    ConfigAdapter(SimConfig_t& simConfig);
};

#endif //MRTPHASESIMULATION_MPISIMCONFIG_H
