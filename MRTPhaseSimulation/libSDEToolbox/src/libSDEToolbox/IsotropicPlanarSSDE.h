//
// Created by konstantin on 1/8/21.
//

#ifndef MRTPHASESIMULATION_ISOTROPICPLANARSSDE_H
#define MRTPHASESIMULATION_ISOTROPICPLANARSSDE_H

#include <libSPhaseFile/libSPhaseFile.h>

class IsotropicPlanarSSDE {
protected:
    ParameterSet pSet;
public:
    //IsotropicPlanarSSDE(ParameterSet& pSet): pSet(pSet){};
    virtual ~IsotropicPlanarSSDE() {};

    virtual double g(double& rho) = 0;
    virtual double q_rho(double& rho) = 0;
    virtual double f(double& rho) = 0;
    virtual double q_phi(double& rho) = 0;
};


#endif //MRTPHASESIMULATION_ISOTROPICPLANARSSDE_H
