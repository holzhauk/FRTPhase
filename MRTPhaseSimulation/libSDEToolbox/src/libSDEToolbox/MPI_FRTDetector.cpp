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
    world_size = world_size;
    world_rank = world_rank;
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
                              Sampler& sampler,
                              SDEIntegrator& integrator) {

    FRTData dataSet(curve.get_name());
    Pen pen(dataSet);
    vector<State_t> Samples = sampler.get_samples(config.SampleSize, curve);
    SDEIntegrator::config_t integratorConfig = SimConfig2IntegratorConfig(config);
    integrator.configure(integratorConfig);

    for (auto x0: Samples){

        pen.write_down(x0);
        Data<StationaryStats> sData(config);
        unsigned int subE = SubEnsembleSize(world_rank, world_size, config.EnsembleSize);

        for (unsigned int e = 0; e < subE; e++){
            integrator.reset(x0, config.t0);
            auto [xT, T] = integrator.integrate();
            sData.add_data_point(xT[1]);
        }

        StationaryStats sStats = MPI_Dist_StatsFromData(sData, config.EnsembleSize);
        pen.write_down(sStats);
        FRTStats frtStats;
        if ((sStats.mPhiT >= x0[1] + 2*M_PI) || (sStats.mPhiT <= x0[1] - 2*M_PI)){
            bool pos_sense_of_rotation = false;
            if (sStats.mPhiT >= x0[1] + 2*M_PI)
                pos_sense_of_rotation = true;

            Data<FRTStats> frtData(config);
            for (unsigned int e = 0; e < subE; e++){
                integrator.reset(x0, config.t0);
                array<double, 2> x = x0;
                double T = config.t0;
                while(!curve.is_first_return_event(x, pos_sense_of_rotation) && integrator.is_in_time()){
                    tie(x, T) = integrator.evolve();
                }
                frtData.add_data_point(T);
            }
            frtStats = MPI_Dist_StatsFromData(frtData, config.EnsembleSize);
        } else {
            frtStats.mFRT = config.T;
            frtStats.varFRT = 0.0;
        }
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

void Pen::write_down(const State_t& x0) {
    dataSet.x0[0].push_back(x0[0]);
    dataSet.x0[1].push_back(x0[1]);
}

void Pen::write_down(const StationaryStats& sStats) {
    dataSet.mPhiT.push_back(sStats.mPhiT);
    dataSet.varPhiT.push_back(sStats.varPhiT);
    dataSet.mT.push_back(sStats.mT);
    dataSet.varT.push_back(sStats.varT);
}

void Pen::write_down(const FRTStats& frtStats) {
    dataSet.mFRT.push_back(frtStats.mFRT);
    dataSet.varFRT.push_back(frtStats.varFRT);
}