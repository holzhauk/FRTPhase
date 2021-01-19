//
// Created by konstantin on 1/10/21.
//

#include <cmath>
#include "HeunIntegrator.h"

HeunIntegrator::HeunIntegrator(unique_ptr<Domain>& domain_ptr,
                               unique_ptr<IsotropicPlanarSSDEAdditiveNoise>& sde_ptr) {
    this->domain_ptr = move(domain_ptr);
    this->sde_ptr = move(sde_ptr);
}

HeunIntegrator::HeunIntegrator(const config_t& config,
                               unique_ptr<Domain>& domain_ptr,
                               unique_ptr<IsotropicPlanarSSDEAdditiveNoise>& sde_ptr): SDEIntegrator(config) {
    this->domain_ptr = move(domain_ptr);
    this->sde_ptr = move(sde_ptr);
}

state_t HeunIntegrator::evolve() {

    array<double, 2> x0 = x;

    /*
     * Heun predictor step.
     * Don't need prediction for Phi, because
     * the dynamic system is isotropic.
     */
    double Rho_p;
    double dW = sqrt(dt)*norm_dist(rn_gen);
    double D_1 = pow(sde_ptr->q_rho(x0[0]), 2.0) / 2.0;
    Rho_p = x0[0] + ( sde_ptr->g(x0[0]) + D_1 / x0[0] )*dt +
                            sde_ptr->q_rho(x0[0])*dW;

    // Heun corrector step
    x[0] = x0[0] + ((sde_ptr->g(x0[0]) + D_1 / x0[0]) + (sde_ptr->g(Rho_p) + D_1 / Rho_p))*dt / 2
            + sde_ptr->q_rho(x0[0])*dW;

    // phi iteration
    dW = sqrt(dt)*norm_dist(rn_gen);
    x[1] = x[1] + (sde_ptr->f(x0[0]) + sde_ptr->f(Rho_p))*dt / 2
            + sde_ptr->q_phi(x0[0])*(1 / x0[0] + 1 / Rho_p)*dW / 2;

    x = domain_ptr->apply_boundary_conditions(x0, x);

    t += dt;

    return tie(x, t);
}