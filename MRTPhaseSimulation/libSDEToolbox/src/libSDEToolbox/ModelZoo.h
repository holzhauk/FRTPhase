//
// Created by konstantin on 2/1/21.
//

#ifndef MRTPHASESIMULATION_MODELZOO_H
#define MRTPHASESIMULATION_MODELZOO_H

#include "IsotropicPlanarSSDE.h"
#include "IsotropicPlanarSSDEAdditiveNoise.h"

namespace ModelZoo {

    class NewbySchwemmer : public IsotropicPlanarSSDEAdditiveNoise {
    public:
        NewbySchwemmer(ParameterSet pSet) :
                IsotropicPlanarSSDEAdditiveNoise(pSet["D"], pSet["D"]) {
            this->pSet = pSet;
        };

        double g(double &rho) override {
            return -pSet["gamma"] * rho * (pow(rho, 2.0) - 1.0) + pSet["D"] / rho;
        };

        double f(double &rho) override {
            return pSet["omega"] * (1 + pSet["gamma"] * pSet["c"] * pow((1.0 - rho), 2.0));
        };
    };

}

#endif //MRTPHASESIMULATION_MODELZOO_H
