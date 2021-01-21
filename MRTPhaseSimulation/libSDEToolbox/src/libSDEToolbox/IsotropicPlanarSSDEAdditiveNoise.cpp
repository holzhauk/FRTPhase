//
// Created by konstantin on 1/10/21.
//

#include <cmath>
#include <iostream>
#include "IsotropicPlanarSSDEAdditiveNoise.h"

IsotropicPlanarSSDEAdditiveNoise::IsotropicPlanarSSDEAdditiveNoise(double D_rho,
                                                                   double D_phi) {
    this->D_rho = D_rho;
    this->D_phi = D_phi;
}

double IsotropicPlanarSSDEAdditiveNoise::q_rho(double& rho) {
    return sqrt(2*D_rho);
}

double IsotropicPlanarSSDEAdditiveNoise::q_phi(double& rho) {
    return sqrt(2*D_phi) / rho;
}