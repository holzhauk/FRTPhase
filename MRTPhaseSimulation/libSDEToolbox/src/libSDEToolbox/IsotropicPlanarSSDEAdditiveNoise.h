//
// Created by konstantin on 1/10/21.
//

#ifndef MRTPHASESIMULATION_ISOTROPICPLANARSSDEADDITIVENOISE_H
#define MRTPHASESIMULATION_ISOTROPICPLANARSSDEADDITIVENOISE_H

#include "IsotropicPlanarSSDE.h"

class IsotropicPlanarSSDEAdditiveNoise: public IsotropicPlanarSSDE {
public:
    IsotropicPlanarSSDEAdditiveNoise(double D_rho, double D_phi);
    double q_rho(double& rho);
    double q_phi(double& phi);
};

#endif //MRTPHASESIMULATION_ISOTROPICPLANARSSDEADDITIVENOISE_H
