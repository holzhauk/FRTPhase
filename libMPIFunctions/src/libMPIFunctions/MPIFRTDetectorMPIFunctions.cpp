//
// Created by konstantin on 1/17/21.
//

#include <cmath>
#include "MPIFRTDetectorMPIFunctions.h"

void MPI_StationaryStats_SUM(StationaryStats* in, StationaryStats* inout, int* len, MPI_Datatype *dtptr);
void MPI_FRTStats_SUM(FRTStats* in, FRTStats* inout, int* len, MPI_Datatype *dtptr);

MPI_Datatype StationaryStats::get_mpi_type() {
    // registration of mpi datatype for communication
    // of extension of data to be shared across the nodes
    MPI_Datatype MPI_StationaryStats_t;
    int blengths[4] = {1, 1, 1, 1};
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[4];
    offsets[0] = offsetof(StationaryStats, mPhiT);
    offsets[1] = offsetof(StationaryStats, varPhiT);
    offsets[2] = offsetof(StationaryStats, mT);
    offsets[3] = offsetof(StationaryStats, varT);
    MPI_Type_create_struct(4, blengths, offsets, types, &MPI_StationaryStats_t);
    MPI_Type_commit(&MPI_StationaryStats_t);
    // return the composit type registered with mpi
    return MPI_StationaryStats_t;
}

StationaryStats& StationaryStats::operator=(const StationaryStats& other) {
    this->mPhiT = other.mPhiT;
    this->varPhiT = other.varPhiT;
    this->mT = other.mT;
    this->varT = other.varT;
    return *this;
}

void MPI_StationaryStats_SUM(StationaryStats* in, StationaryStats* inout, int* len, MPI_Datatype *dtptr){
    StationaryStats result;
    for (int i = 0; i < *len; i++){
        result.mPhiT = in->mPhiT + inout->mPhiT;
        result.varPhiT = in->varPhiT + inout->varPhiT;
        result.mT = in->mT + inout->mT;
        result.varT = in->varT + inout->varT;
        *inout = result;
        in++;
        inout++;
    }
}

MPI_Datatype FRTStats::get_mpi_type() {
    // registration of mpi datatype for communication
    // of extension of data to be shared across the nodes
    MPI_Datatype MPI_FRTStats_t;
    int blengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[2];
    offsets[0] = offsetof(FRTStats, mFRT);
    offsets[1] = offsetof(FRTStats, varFRT);
    MPI_Type_create_struct(2, blengths, offsets, types, &MPI_FRTStats_t);
    MPI_Type_commit(&MPI_FRTStats_t);
    // return the composit type registered with mpi
    return MPI_FRTStats_t;
}

FRTStats& FRTStats::operator=(const FRTStats& other) {
    this->mFRT = other.mFRT;
    this->varFRT = other.varFRT;
    return *this;
}

void MPI_FRTStats_SUM(FRTStats* in, FRTStats* inout, int* len, MPI_Datatype *dtptr){
    FRTStats result;
    for (int i = 0; i < *len; i++){
        result.mFRT = in->mFRT + inout->mFRT;
        result.varFRT = in->varFRT + inout->varFRT;
        *inout = result;
        in++;
        inout++;
    }
}

template<>
Data<StationaryStats>::XSums Data<StationaryStats>::get_sums() const{
    StationaryStats sums;
    double Tbar;
    for (auto d: data){
        sums.mPhiT += d;
        sums.varPhiT += pow(d, 2);
        Tbar = 2*M_PI*config.T / d;
        sums.mT += Tbar;
        sums.varT += pow(Tbar, 2);
    }
    return sums;
}

template<>
StationaryStats MPI_Dist_StatsFromData(const Data<StationaryStats>& data, const unsigned int& ensemble_size) {
    StationaryStats stationaryStats;
    StationaryStats sums = data.get_sums();
    MPI_Op MPI_SStats_SUM;
    MPI_Datatype MPI_SStats_t = stationaryStats.get_mpi_type();
    MPI_Op_create((MPI_User_function*)MPI_StationaryStats_SUM, 1, &MPI_SStats_SUM);
    MPI_Allreduce(&sums, &stationaryStats, 1, MPI_SStats_t, MPI_SStats_SUM, MPI_COMM_WORLD);
    MPI_Type_free(&MPI_SStats_t);
    MPI_Op_free(&MPI_SStats_SUM);

    auto E = static_cast<double>(ensemble_size);
    stationaryStats.mPhiT /= E;
    stationaryStats.varPhiT /= E;
    stationaryStats.varPhiT -= pow(stationaryStats.mPhiT, 2);
    stationaryStats.mT /= E;
    stationaryStats.varT /= E;
    stationaryStats.varT -= pow(stationaryStats.mT, 2);

    return stationaryStats;
}

template<>
Data<FRTStats>::XSums Data<FRTStats>::get_sums() const {
    FRTStats sums;
    for (auto d: data){
        sums.mFRT += d;
        sums.varFRT += pow(d, 2);
    }
    return sums;
}

template<>
FRTStats MPI_Dist_StatsFromData(const Data<FRTStats>& data, const unsigned int& ensemble_size){
    FRTStats frtStats;
    FRTStats sums = data.get_sums();
    MPI_Op MPI_FRTStats_SUM_Op_t;
    MPI_Datatype MPI_FRTStats_t = frtStats.get_mpi_type();
    MPI_Op_create((MPI_User_function*)MPI_FRTStats_SUM, 1, &MPI_FRTStats_SUM_Op_t);
    MPI_Allreduce(&sums, &frtStats, 1, MPI_FRTStats_t, MPI_FRTStats_SUM_Op_t, MPI_COMM_WORLD);
    MPI_Type_free(&MPI_FRTStats_t);
    MPI_Op_free(&MPI_FRTStats_SUM_Op_t);

    auto E = static_cast<double>(ensemble_size);
    frtStats.mFRT /= E;
    frtStats.varFRT /= E;
    frtStats.varFRT -= pow(frtStats.mFRT, 2);

    return frtStats;
}