//
// Created by konstantin on 4/6/20.
//

#include <iostream>

#include "IsoPlanarSOsc.h"

IsoPlanarSOsc::IsoPlanarSOsc(Domain& d){
    Disk = d;
    /*
     * assign default values
     */
    D = 0.198;
    x0 = {0.1, M_PI / 2.0, 0.0};
    dt = 0.01;
    T = 10.0;
}

IsoPlanarSOsc::IsoPlanarSOsc(Domain& d, sim_config_t& sim_config) {
    Disk = d;
    // default noise intensity
    D = 0.198;
    // values from configuration structure
    x0 = sim_config.x0;
    dt = sim_config.dt;
    T = sim_config.T;
}

IsoPlanarSOsc::IsoPlanarSOsc(Domain& domain, sim_config_t& sim_config, double noise_intensity) {
    Disk = domain;
    D = noise_intensity;
    x0 = sim_config.x0;
    dt = sim_config.dt;
    T = sim_config.T;
}

IsoPlanarSOsc::IsoPlanarSOsc(Domain& domain, sim_config_C_t& sim_config) {
    Disk = domain;
    D = 0.198;
    x0[0] = sim_config.x0[0]; // Rho0
    x0[1] = sim_config.x0[1]; // Phi0
    x0[2] = sim_config.x0[2]; // t0
    dt = sim_config.dt;
    T = sim_config.T;
}

IsoPlanarSOsc::IsoPlanarSOsc(Domain& domain, sim_config_C_t& sim_config, double noise_intensity) {
    Disk = domain;
    D = noise_intensity;
    x0[0] = sim_config.x0[0]; // Rho0
    x0[1] = sim_config.x0[1]; // Phi0
    x0[2] = sim_config.x0[2]; // t0
    dt = sim_config.dt;
    T = sim_config.T;
}

void IsoPlanarSOsc::configure(Domain& domain, sim_config_t& config, double noise_intensity) {
    x0 = config.x0;
    dt = config.dt;
    T = config.T;
    D = noise_intensity;
    Disk = domain;
}

void IsoPlanarSOsc::configure(Domain& domain, sim_config_C_t& config, double noise_intensity) {
    for (int i = 0; i < x0.size(); i++){
        x0[i] = config.x0[i];
    }
    dt = config.dt;
    T = config.T;
    D = noise_intensity;
    Disk = domain;
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