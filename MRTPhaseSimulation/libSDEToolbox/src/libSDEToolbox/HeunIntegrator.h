//
// Created by konstantin on 1/10/21.
//

#ifndef MRTPHASESIMULATION_HEUNINTEGRATOR_H
#define MRTPHASESIMULATION_HEUNINTEGRATOR_H

#include <random>
#include "SDEIntegrator.h"
#include "Domain.h"
#include "IsotropicPlanarSSDEAdditiveNoise.h"

class HeunIntegrator: public SDEIntegrator {
private:
    random_device rdn_dev {};
    mt19937_64 rn_gen {rdn_dev()};
    normal_distribution<double> norm_dist {0.0, 1.0};
    Domain* domain_ptr;
    IsotropicPlanarSSDEAdditiveNoise* sde_ptr;
public:
    HeunIntegrator(Domain* domain_ptr,
                   IsotropicPlanarSSDEAdditiveNoise* sde_ptr);
    HeunIntegrator(const config_t& config,
                   Domain* domain_ptr,
                   IsotropicPlanarSSDEAdditiveNoise* sde_ptr);
    state_t evolve();
};


#endif //MRTPHASESIMULATION_HEUNINTEGRATOR_H
