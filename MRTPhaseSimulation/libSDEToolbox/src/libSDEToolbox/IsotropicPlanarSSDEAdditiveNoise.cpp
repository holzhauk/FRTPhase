//
// Created by konstantin on 1/10/21.
//

#include <cmath>
#include "IsotropicPlanarSSDEAdditiveNoise.h"

IsotropicPlanarSSDEAdditiveNoise::IsotropicPlanarSSDEAdditiveNoise(double D_rho,
                                                                   double D_phi) {
    pSet["D_rho"] = D_rho;
    pSet["D_phi"] = D_phi;
}

double IsotropicPlanarSSDEAdditiveNoise::q_rho(double& rho) {
    return sqrt(2*pSet["D_rho"]);
}

double IsotropicPlanarSSDEAdditiveNoise::q_phi(double& phi) {
    return sqrt(2*pSet["D_phi"]);
}