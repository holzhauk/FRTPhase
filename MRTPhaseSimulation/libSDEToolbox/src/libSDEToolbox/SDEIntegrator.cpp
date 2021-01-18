//
// Created by konstantin on 1/10/21.
//

#include "SDEIntegrator.h"

void SDEIntegrator::configure(const config_t& config){
    x = config.x0;
    t = config.t0;
    dt = config.dt;
    T = config.T;
}

void SDEIntegrator::reset(array<double, 2>& x0, double& t0) {
    x = x0;
    t = t0;
}

state_t SDEIntegrator::integrate() {
    while(t < T)
        tie(x, t) = this->evolve();
    return tie(x, t);
};

bool SDEIntegrator::is_in_time() {
    return (t < T);
}
