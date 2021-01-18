//
// Created by konstantin on 1/8/21.
//

#ifndef MRTPHASESIMULATION_SDEINTEGRATOR_H
#define MRTPHASESIMULATION_SDEINTEGRATOR_H

#include <array>
#include <memory>
#include <tuple>

#include "Domain.h"
#include "IsotropicPlanarSSDE.h"

using namespace std;
using state_t = tuple<array<double, 2>, double>;

class SDEIntegrator {
public:
    struct config_t {
        array<double, 2> x0;
        double t0;
        double dt;
        double T;
    };
protected:
    array<double, 2> x;
    double t;
    double dt;
    double T;
public:
    SDEIntegrator(const config_t& config){
       this->configure(config);
    };
    virtual ~SDEIntegrator() {};
    virtual state_t evolve() = 0;
    void configure(const config_t& config);
    void reset(array<double, 2>& x0, double& t0);
    state_t integrate();
    bool is_in_time();
};

#endif //MRTPHASESIMULATION_SDEINTEGRATOR_H
