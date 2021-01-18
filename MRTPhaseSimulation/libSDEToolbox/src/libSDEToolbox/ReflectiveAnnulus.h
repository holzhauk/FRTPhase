//
// Created by konstantin on 1/10/21.
//

#ifndef MRTPHASESIMULATION_ANNULUS_H
#define MRTPHASESIMULATION_ANNULUS_H

#include "Domain.h"

class ReflectiveAnnulus : public Domain{
private:
    double rho_min; // inner boundary radius
    double rho_max; // outer boundary radius
    double width;
public:
    ReflectiveAnnulus(double rho_min, double rho_max);
    array<double, 2> apply_boundary_conditions(array<double, 2>& x_init, array<double, 2>& x_final);
};


#endif //MRTPHASESIMULATION_ANNULUS_H
