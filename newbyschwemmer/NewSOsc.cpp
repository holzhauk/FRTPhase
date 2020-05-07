//
// Created by konstantin on 4/7/20.
//

#include "NewSOsc.h"

NewSOsc::NewSOsc(Domain& domain, NewSOsc::config_t& config, IsoPlanarSOsc::pSet_t* pSet_ptr)
    :IsoPlanarSOsc(domain, config, pSet_ptr) {
    auto NewS_pSet_ptr = static_cast<NewSOsc::pSet_t*>(pSet_ptr);
    omega = NewS_pSet_ptr->omega;
    gamma = NewS_pSet_ptr->gamma;
    c = NewS_pSet_ptr->c;
}

void NewSOsc::configure(Domain& domain, NewSOsc::config_t& config, IsoPlanarSOsc::pSet_t* pSet_ptr) {
    auto NewS_pSet_ptr = static_cast<NewSOsc::pSet_t*>(pSet_ptr);
    omega = NewS_pSet_ptr->omega;
    gamma = NewS_pSet_ptr->gamma;
    c = NewS_pSet_ptr->c;
    IsoPlanarSOsc::configure(domain, config, pSet_ptr);
}

/*
 * Define configuration interface
 */
void NewSOsc::pSet_t::load(std::map<std::string, double> &pMap) {
    D = pMap["D"];
    omega = pMap["omega"];
    gamma = pMap["gamma"];
    c = pMap["c"];
}

/*
 * Define pure virtual methods describing deterministic
 * dynamics of the system
 */
double NewSOsc::g(double Rho) { // radial dynamics
    return -gamma * Rho * (pow(Rho, 2.0) - 1.0);
}

double NewSOsc::f(double Rho) {// angular dynamics
    return omega * (1.0 + gamma * c * pow(1.0 - Rho, 2.0));
}

MPI_Datatype NewSOsc::pSet_t::mpiType() {
    MPI_Datatype mpi_pSet_t;
    int blengths[4] = {1, 1, 1, 1};
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[4];
    offsets[0] = offsetof(NewSOsc::pSet_t, D);
    offsets[1] = offsetof(NewSOsc::pSet_t, omega);
    offsets[2] = offsetof(NewSOsc::pSet_t, gamma);
    offsets[3] = offsetof(NewSOsc::pSet_t, c);
    MPI_Type_create_struct(4, blengths, offsets, types, &mpi_pSet_t);
    MPI_Type_commit(&mpi_pSet_t);
    return mpi_pSet_t;
}
