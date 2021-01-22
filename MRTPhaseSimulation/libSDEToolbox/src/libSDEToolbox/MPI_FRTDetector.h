//
// Created by konstantin on 1/17/21.
//

#ifndef MRTPHASESIMULATION_MPI_FRTDETECTOR_H
#define MRTPHASESIMULATION_MPI_FRTDETECTOR_H

#include <array>
#include <mpi.h>
#include <libSPhaseFile/libSPhaseFile.h>
#include <libMPIFunctions/libMPIFunctions.h>
#include "Sampler.h"
#include "SDEIntegrator.h"

using State_t = array<double, 2>;

namespace MPI {

    class FRTDetector {
    private:
        bool has_initiated_mpi;
        int world_rank, world_size;
    public:
        FRTDetector();
        FRTDetector(int world_rank, int world_size);
        ~FRTDetector();
        bool is_master();
        int get_worldSize();
        FRTData run(Config::Simulation& config,
                        InterpolatedCurve& curve,
                        Sampler* sampler_ptr,
                        SDEIntegrator* integrator_ptr);
    };

};



class SimConfig2IntegratorConfig: public SDEIntegrator::config_t{
public:
    SimConfig2IntegratorConfig(const Config::Simulation& config);
    SimConfig2IntegratorConfig& operator = (const SimConfig2IntegratorConfig& other);
};

class Pen {
private:
    FRTData& dataSet;
public:
    Pen(FRTData& dataSet): dataSet(dataSet) {};
    void write_down(const State_t& x0);
    void write_down(const StationaryStats& sStats);
    void write_down(const FRTStats& frtStats);
};
#endif //MRTPHASESIMULATION_MPI_FRTDETECTOR_H
