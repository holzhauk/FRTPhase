//
// Created by konstantin on 4/6/20.
//

#include "IsoPlanarSOsc.h"

IsoPlanarSOsc::IsoPlanarSOsc(Domain& domain, config_t& config, IsoPlanarSOsc::pSet_t* pSet_ptr) {
    this->configure(domain, config, pSet_ptr);
}

void IsoPlanarSOsc::configure(Domain& domain, config_t& config, IsoPlanarSOsc::pSet_t* pSet_ptr) {
    Disk = domain;
    D = pSet_ptr->D;
    dt = config.dt;
    T = config.T;
    std::copy(std::begin(config.x0), std::end(config.x0), std::begin(x0));
}

std::array<double, 3> IsoPlanarSOsc::evolve() {
    /*
     * The time evolution of the system is determined by
     * the solution of the SDE.
     * Here, the integration method of Heun is implemented
     *      -> Stratonovich interpretation
     */
    auto [Rho_0, Phi_0, t0] = x0;

    // draw Wiener increments
    double dW_Rho = sqrt(dt)*ND_Rho(rd_Rho);
    double dW_Phi = sqrt(dt)*ND_Phi(rd_Phi);

    /*
     * Heun predictor step.
     * Don't need prediction for Phi, because
     * the dynamic system is isotropic.
     */

    double Rho_p;
    Rho_p = Rho_0 + (this->g(Rho_0) + D / Rho_0)*dt + sqrt(2*D)*dW_Rho;

    // Heun corrector step
    double Rho_n, Phi_n;
    Rho_n = Rho_0 + ((this->g(Rho_0) + D / Rho_0) + (this->g(Rho_p) + D / Rho_p))*dt / 2
                + sqrt(2*D)*dW_Rho;
    Phi_n = Phi_0 + (this->f(Rho_0) + this->f(Rho_p))*dt / 2
                + sqrt(2*D)*(1 / Rho_0 + 1 / Rho_p)*dW_Phi / 2;

    /*
     * apply reflecting boundary condition
     */
    std::array<double, 2> pos = {Rho_n, Phi_n};
    pos = Disk.apply(pos);

    // new flow vector
    x0 = {pos[0], pos[1], t0 + dt};
    return x0;
}

bool IsoPlanarSOsc::in_time() const {
    return (x0[2] < T);
}

std::array<double, 3> IsoPlanarSOsc::get_state() const {
    return x0;
}

/*
 * Iterator operator overloads
 */
bool IsoPlanarSOsc::IsoPlanarSOscIt::operator!=(IsoPlanarSOscIt& iterator) {
    return (this->t < iterator.t + iterator.model_ref.dt);
}

IsoPlanarSOsc& IsoPlanarSOsc::IsoPlanarSOscIt::operator++() {
    t = model_ref.evolve()[2];
    return model_ref;
}

IsoPlanarSOsc& IsoPlanarSOsc::IsoPlanarSOscIt::operator*() const {
    return this->model_ref;
}

/*
 * Configuration interface
 */
void IsoPlanarSOsc::pSet_t::load(std::map<std::string, double> &pMap) {
    D = pMap["D"];
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