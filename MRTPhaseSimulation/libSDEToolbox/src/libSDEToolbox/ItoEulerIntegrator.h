//
// Created by konstantin on 1/10/21.
//

#ifndef MRTPHASESIMULATION_ITOEULERINTEGRATOR_H
#define MRTPHASESIMULATION_ITOEULERINTEGRATOR_H

#include <random>
#include "SDEIntegrator.h"

class ItoEulerIntegrator: public SDEIntegrator {
private:
    random_device rdn_dev {};
    mt19937_64 rn_gen {rdn_dev()};
    normal_distribution<double> norm_dist {0.0, 1.0};
    Domain* domain_ptr;
    IsotropicPlanarSSDE* sde_ptr;
public:
    ItoEulerIntegrator(Domain* domain_ptr,
                        IsotropicPlanarSSDE* sde_ptr);
    ItoEulerIntegrator(const config_t& config,
                       Domain* domain_ptr,
                       IsotropicPlanarSSDE* sde_ptr);
    state_t evolve();
};


#endif //MRTPHASESIMULATION_ITOEULERINTEGRATOR_H
