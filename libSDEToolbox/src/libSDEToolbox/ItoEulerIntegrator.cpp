//
// Created by konstantin on 1/10/21.
//

#include <cmath>
#include "ItoEulerIntegrator.h"

ItoEulerIntegrator::ItoEulerIntegrator(Domain* domain_ptr,
                                       IsotropicPlanarSSDE* sde_ptr){
    this->domain_ptr = domain_ptr;
    this->sde_ptr = sde_ptr;
}

ItoEulerIntegrator::ItoEulerIntegrator(const config_t& config,
                                       Domain* domain_ptr,
                                       IsotropicPlanarSSDE* sde_ptr): SDEIntegrator(config){
    this->domain_ptr = domain_ptr;
    this->sde_ptr = sde_ptr;
}

state_t ItoEulerIntegrator::evolve(){
    array<double, 2> x0 = x;
    //rho
    double dW = sqrt(dt)*norm_dist(rn_gen);
    x[0] = x0[0] + sde_ptr->g(x0[0])*dt + sde_ptr->q_rho(x0[0])*dW;
    //phi
    dW = sqrt(dt)*norm_dist(rn_gen);
    x[1] = x0[1] + sde_ptr->f(x0[0])*dt + sde_ptr->q_phi(x0[0])*dW;

    x = domain_ptr->apply_boundary_conditions(x0, x);

    //t
    t += dt;

    return tie(x, t);
}