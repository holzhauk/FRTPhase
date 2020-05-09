//
// Created by konstantin on 5/4/20.
//

#include "Adapters.h"

using namespace Adapters;

Config2SimConfig::Config2SimConfig(SimConfig_t &simConfig) {
    dt = simConfig.Simulation.dt;
    T = simConfig.Simulation.T;
    x0[0] = 1.5;
    x0[1] = 0.0;
    x0[2] = simConfig.Simulation.t0;
}