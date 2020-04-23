//
// Created by konstantin on 4/7/20.
//

#include "NewSOsc.h"

NewSOsc::NewSOsc(Domain& domain):IsoPlanarSOsc(domain) {
    // default values = antirotating case
    omega = 1.0;
    gamma = 15.0;
    c = -15.0;
}

NewSOsc::NewSOsc(Domain& domain, sim_config_t& sim_config):IsoPlanarSOsc(domain, sim_config) {
    // default values = antirotating case
    omega = 1.0;
    gamma = 15.0;
    c = -15.0;
}

NewSOsc::NewSOsc(Domain& domain, pset_NewS_t& pset_News, sim_config_t& sim_config)
    :IsoPlanarSOsc(domain, sim_config, pset_News.D) {
    omega = pset_News.omega;
    gamma = pset_News.gamma;
    c = pset_News.c;
}

NewSOsc::NewSOsc(Domain& domain, sim_config_C_t& sim_config)
    :IsoPlanarSOsc(domain, sim_config) {
    // default values = antirotating case
    omega = 1.0;
    gamma = 15.0;
    c = -15.0;
}

NewSOsc::NewSOsc(Domain& domain, pset_NewS_t& pset_NewS, sim_config_C_t& sim_config)
    :IsoPlanarSOsc(domain, sim_config, pset_NewS.D) {
    omega = pset_NewS.omega;
    gamma = pset_NewS.omega;
    c = pset_NewS.c;

}

void NewSOsc::configure(Domain& domain, pset_NewS_t& pset_NewS, sim_config_t& config) {
    omega = pset_NewS.omega;
    gamma = pset_NewS.gamma;
    c = pset_NewS.c;
    IsoPlanarSOsc::configure(domain, config, pset_NewS.D);
}

void NewSOsc::configure(Domain& domain, pset_NewS_t& pset_NewS, sim_config_C_t& config) {
    omega = pset_NewS.omega;
    gamma = pset_NewS.gamma;
    c = pset_NewS.c;
    IsoPlanarSOsc::configure(domain, config, pset_NewS.D);
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

