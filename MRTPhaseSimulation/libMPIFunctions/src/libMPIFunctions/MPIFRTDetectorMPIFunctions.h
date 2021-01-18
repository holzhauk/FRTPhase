//
// Created by konstantin on 1/17/21.
//

#ifndef MRTPHASESIMULATION_MPIFRTDETECTORMPIFUNCTIONS_H
#define MRTPHASESIMULATION_MPIFRTDETECTORMPIFUNCTIONS_H

#include <mpi.h>
#include <vector>
#include <libSPhaseFile/libSPhaseFile.h>

using namespace std;

struct StationaryStats{
    // asymptotic angle PhiT = phi(T)
    double mPhiT = 0.0;
    double varPhiT = 0.0;
    // stationary period
    double mT = 0.0;
    double varT = 0.0;
    StationaryStats() = default;
    MPI_Datatype get_mpi_type();
    StationaryStats& operator = (const StationaryStats& other);
};

struct FRTStats{
    // first return time FRT
    double mFRT = 0.0;
    double varFRT = 0.0;
    FRTStats() = default;
    MPI_Datatype get_mpi_type();
    FRTStats& operator = (const FRTStats& other);
};

template<class X>
struct Data {
    using XSums = X;
    template<class Y>
    friend Y MPI_Dist_StatsFromData(const Data<Y>& data, const unsigned int& ensemble_size);
private:
    Config::Simulation config;
    vector<double> data;
    X get_sums() const;
public:
    Data(const Config::Simulation& config): config(config) {};
    void add_data_point(double& dp);
};

template<class X>
void Data<X>::add_data_point(double& dp) {
    data.push_back(dp);
}

#endif //MRTPHASESIMULATION_MPIFRTDETECTORMPIFUNCTIONS_H
