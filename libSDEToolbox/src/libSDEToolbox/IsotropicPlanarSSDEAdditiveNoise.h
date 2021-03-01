//
// Created by konstantin on 1/10/21.
//

#ifndef MRTPHASESIMULATION_ISOTROPICPLANARSSDEADDITIVENOISE_H
#define MRTPHASESIMULATION_ISOTROPICPLANARSSDEADDITIVENOISE_H

#include "IsotropicPlanarSSDE.h"

class IsotropicPlanarSSDEAdditiveNoise: public IsotropicPlanarSSDE {
private:
    double D_rho, D_phi;
public:
    IsotropicPlanarSSDEAdditiveNoise(double D_rho, double D_phi);
    double q_rho(double& rho) override;
    double q_phi(double& rho) override;
};

#endif //MRTPHASESIMULATION_ISOTROPICPLANARSSDEADDITIVENOISE_H
