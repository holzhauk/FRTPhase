//
// Created by konstantin on 1/17/21.
//
#include <cmath>
#include "MPI_FRTDetector.h"

MPI::FRTDetector::FRTDetector() {
    // configure and init mpi
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    this->has_initiated_mpi = true;
}

MPI::FRTDetector::FRTDetector(int world_rank, int world_size) {
    this->world_size = world_size;
    this->world_rank = world_rank;
    this->has_initiated_mpi = false;
}

MPI::FRTDetector::~FRTDetector() {
    if (this->has_initiated_mpi == true)
        MPI_Finalize();
}

bool MPI::FRTDetector::is_master() {
    return (world_rank == 0);
}

int MPI::FRTDetector::get_worldSize() {
    return int(world_size);
}

FRTData MPI::FRTDetector::run(Config::Simulation& config,
                              InterpolatedCurve& curve,
                              Sampler* sampler_ptr,
                              SDEIntegrator* integrator_ptr) {

    FRTData dataSet(curve.get_name());
    Pen pen(dataSet);
    vector<Pos_t> Samples = sampler_ptr->get_samples(config.SampleSize, curve);
    SDEIntegrator::config_t integratorConfig = SimConfig2IntegratorConfig(config);
    integrator_ptr->configure(integratorConfig);

    double omegaBar = curve.get_omegaBar();
    if (omegaBar == 0.0)
        throw FRTDetectorNoSenseOfRotation();
    bool is_pos_sense_of_rotation = (omegaBar > 0.0);

    for (auto x0: Samples){
        pen.write_down(x0);
        unsigned int subE = SubEnsembleSize(world_rank, world_size, config.EnsembleSize);
        FRTStats frtStats;

        Data<FRTStats> frtData(config);
        for (unsigned int e = 0; e < subE; e++){
            integrator_ptr->reset(x0, config.t0);
            array<double, 2> x = x0;
            double T = config.t0;
            while(!curve.is_first_return_event(x, is_pos_sense_of_rotation) && integrator_ptr->is_in_time()){
                tie(x, T) = integrator_ptr->evolve();
            }
            frtData.add_data_point(T);
        }
        frtStats = MPI_Dist_StatsFromData(frtData, config.EnsembleSize);
        pen.write_down(frtStats);
    }

    return dataSet;
}

SimConfig2IntegratorConfig::SimConfig2IntegratorConfig(const Config::Simulation& config) {
    this->x0 = {0.0, 0.0};
    this->t0 = config.t0;
    this->dt = config.dt;
    this->T = config.T;
}

SimConfig2IntegratorConfig& SimConfig2IntegratorConfig::operator = (const SimConfig2IntegratorConfig& other) {
    this->x0 = other.x0;
    this->t0 = other.t0;
    this->dt = other.dt;
    this->T = other.T;
    return *this;
}

void Pen::write_down(const Pos_t& x0) {
    dataSet.x0[0].push_back(x0[0]);
    dataSet.x0[1].push_back(x0[1]);
}

void Pen::write_down(const FRTStats& frtStats) {
    dataSet.mFRT.push_back(frtStats.mFRT);
    dataSet.varFRT.push_back(frtStats.varFRT);
}